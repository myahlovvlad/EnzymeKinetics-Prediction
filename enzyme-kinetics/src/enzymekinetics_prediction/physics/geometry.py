"""Canonical equal-weight structural metrics and compact-assembly PBC geometry.

Coordinates, lengths and row box vectors are in nm. Only axis-aligned positive
orthorhombic cells are supported. Reconstruction is spatial, not temporal
unwrapping: it assumes bonded molecules and a compact nearest-image assembly.
"""

import numpy as np

STRUCTURAL_METHOD = "orthorhombic-bonded-assembly-kabsch-v1"


def coordinates(value):
    xyz = np.asarray(value, dtype=float)
    if (
        xyz.ndim != 3
        or xyz.shape[0] < 1
        or xyz.shape[1] < 1
        or xyz.shape[2] != 3
        or not np.isfinite(xyz).all()
    ):
        raise ValueError("Expected finite [frames, atoms, 3] coordinates in nm")
    return xyz


def orthorhombic_lengths(box_vectors_nm):
    box = np.asarray(box_vectors_nm, dtype=float)
    if box.shape[-2:] != (3, 3) or not np.isfinite(box).all():
        raise ValueError("Expected finite (..., 3, 3) box vectors in nm")
    lengths = np.diagonal(box, axis1=-2, axis2=-1)
    diagonal = lengths[..., :, None] * np.eye(3)
    if np.any(lengths <= 0) or not np.allclose(box, diagonal, rtol=0, atol=1e-6):
        raise ValueError("Only positive axis-aligned orthorhombic cells are supported")
    return lengths.copy()


def minimum_image(displacement, lengths_nm):
    delta = np.asarray(displacement, dtype=float)
    lengths = np.asarray(lengths_nm, dtype=float)
    if lengths.shape != (3,) or not np.isfinite(lengths).all() or np.any(lengths <= 0):
        raise ValueError("Expected positive finite orthorhombic lengths (3,)")
    if delta.ndim < 1 or delta.shape[-1] != 3 or not np.isfinite(delta).all():
        raise ValueError("Expected finite displacements (..., 3)")
    # Half-box ties follow NumPy ties-to-even; these images are nonunique.
    return delta - lengths * np.rint(delta / lengths)


def _unique_image(displacement, lengths_nm):
    wrapped = minimum_image(displacement, lengths_nm)
    if np.any(
        np.isclose(np.abs(wrapped), np.asarray(lengths_nm) / 2, rtol=0, atol=1e-10)
    ):
        raise ValueError("Ambiguous half-box image in molecular reconstruction")
    return wrapped


def atom_indices(value, count):
    raw = np.asarray(value)
    if (
        raw.ndim != 1
        or raw.size == 0
        or not np.issubdtype(raw.dtype, np.integer)
        or np.any(raw < 0)
        or np.any(raw >= count)
        or len(np.unique(raw)) != len(raw)
    ):
        raise ValueError("Expected nonempty unique integer atom indices in range")
    return raw.astype(int)


def reconstruct_assembly(coordinates_nm, box_lengths_nm, groups):
    """Place already-whole groups nearest the first group's geometric center.

    Returns a copy and per-frame/group translations applied in nm. All atoms
    must belong to exactly one group. The first group is the explicit anchor.
    This assumption is appropriate for a compact complex, not diffusing solutes.
    """
    xyz = coordinates(coordinates_nm).copy()
    lengths = np.asarray(box_lengths_nm, dtype=float)
    if lengths.shape != (len(xyz), 3):
        raise ValueError("Expected orthorhombic lengths for every frame")
    ids = [atom_indices(g, xyz.shape[1]) for g in groups]
    flat = np.concatenate(ids) if ids else np.array([], dtype=int)
    if sorted(flat.tolist()) != list(range(xyz.shape[1])):
        raise ValueError("Groups must partition all atoms without overlap")
    shifts = np.zeros((len(xyz), len(ids), 3))
    for frame, box, moved in zip(xyz, lengths, shifts, strict=True):
        minimum_image(np.zeros(3), box)  # Validate even a single-group system.
        anchor = frame[ids[0]].mean(axis=0)
        for i, group in enumerate(ids[1:], 1):
            delta = frame[group].mean(axis=0) - anchor
            moved[i] = _unique_image(delta, box) - delta
            frame[group] += moved[i]
    return xyz, shifts


def make_whole(coordinates_nm, box_vectors_nm, bonds):
    """Reconstruct connected bond graphs, then the compact assembly.

    Bond cycles must close consistently. Missing connectivity is not inferred
    from distances; callers must validate molecular identity in their topology.
    """
    xyz = coordinates(coordinates_nm).copy()
    lengths = orthorhombic_lengths(box_vectors_nm)
    if lengths.shape != (len(xyz), 3):
        raise ValueError("Expected one box vector matrix per frame")
    graph = [set() for _ in range(xyz.shape[1])]
    for bond in bonds:
        a, b = atom_indices(bond, xyz.shape[1])
        graph[a].add(b)
        graph[b].add(a)
    groups, edges, seen = [], [], set()
    for root in range(len(graph)):
        if root in seen:
            continue
        group = [root]
        seen.add(root)
        for a in group:
            for b in sorted(graph[a]):
                if b not in seen:
                    seen.add(b)
                    group.append(b)
                    edges.append((a, b))
        groups.append(np.asarray(group))
    for frame, original, box in zip(
        xyz, coordinates(coordinates_nm), lengths, strict=True
    ):
        for a, b in edges:
            frame[b] = frame[a] + _unique_image(original[b] - original[a], box)
        for a, neighbors in enumerate(graph):
            for b in neighbors:
                if b > a and not np.allclose(
                    frame[b] - frame[a],
                    _unique_image(original[b] - original[a], box),
                    rtol=0,
                    atol=1e-6,
                ):
                    raise ValueError("Bond cycle has inconsistent periodic images")
    whole, shifts = reconstruct_assembly(xyz, lengths, groups)
    return whole, groups, shifts


def structural_metrics(coordinates_nm):
    """Metrics on whole, coherently imaged coordinates; no hidden PBC guessing.

    Fit all supplied atoms equally to frame 0 with a proper Kabsch rotation.
    RMSF/DCCM use the time mean after that fit, population normalization 1/T.
    Rg is geometric (unweighted). Zero-variance DCCM entries are zero by convention.
    """
    xyz = coordinates(coordinates_nm)
    reference = xyz[0] - xyz[0].mean(axis=0)
    aligned = np.empty_like(xyz)
    for i, frame in enumerate(xyz):
        centered = frame - frame.mean(axis=0)
        u, _, vt = np.linalg.svd(centered.T @ reference)
        correction = np.eye(3)
        correction[-1, -1] = np.linalg.det(u @ vt)
        aligned[i] = centered @ u @ correction @ vt
    fluct = aligned - aligned.mean(axis=0)
    covariance = np.einsum("tix,tjx->ij", fluct, fluct) / len(xyz)
    normalization = np.sqrt(np.diag(covariance)[:, None] * np.diag(covariance)[None, :])
    return {
        "coordinates_nm": xyz.copy(),
        "aligned_nm": aligned,
        "rmsd_nm": np.sqrt(np.mean(np.sum((aligned - reference) ** 2, axis=2), axis=1)),
        "rmsf_nm": np.sqrt(np.mean(np.sum(fluct**2, axis=2), axis=0)),
        "radius_gyration_nm": np.sqrt(
            np.mean(
                np.sum((xyz - xyz.mean(axis=1, keepdims=True)) ** 2, axis=2), axis=1
            )
        ),
        "dccm": np.divide(
            covariance,
            normalization,
            out=np.zeros_like(covariance),
            where=normalization > 1e-15,
        ),
    }
