"""Euclidean force projection for fixed distance constraints.

The supplied displacement must match the constraint geometry used by the
system. OpenMM distance constraints here use unwrapped Cartesian positions.
This module does not perform minimization or decide its stop reason.
"""

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import lsmr


def project_forces(positions, forces, pairs, *, atol=1e-11, maxiter=None):
    """Return orthogonal projection onto null(J), plus diagnostic residual.

    J rows contain unit bond vectors at the two affected atoms. Solve the
    least-squares problem J.T lambda ~= F using LSMR; this is also valid for
    dependent constraints. Check the normal-equation residual explicitly and
    refuse results that do not satisfy the requested tangent tolerance.
    """
    x = np.asarray(positions, dtype=float)
    f = np.asarray(forces, dtype=float)
    p = np.asarray(pairs, dtype=int)
    if x.ndim != 2 or x.shape[1] != 3 or f.shape != x.shape:
        raise ValueError("positions and forces must have matching (N,3) shape")
    if not np.isfinite(x).all() or not np.isfinite(f).all():
        raise ValueError("positions and forces must be finite")
    if p.size == 0:
        return f.copy(), {"tangent_residual_max": 0.0, "iterations": 0}
    if (
        p.ndim != 2
        or p.shape[1] != 2
        or (p < 0).any()
        or (p >= len(x)).any()
        or (p[:, 0] == p[:, 1]).any()
    ):
        raise ValueError("invalid constraint indices")
    d = x[p[:, 0]] - x[p[:, 1]]
    length = np.linalg.norm(d, axis=1)
    if np.any(length < 1e-12):
        raise ValueError("degenerate distance constraint")
    u = d / length[:, None]
    rows = np.repeat(np.arange(len(p)), 6)
    cols = np.concatenate(
        (3 * p[:, 0, None] + np.arange(3), 3 * p[:, 1, None] + np.arange(3)), axis=1
    ).ravel()
    values = np.concatenate((u, -u), axis=1).ravel()
    j = coo_matrix((values, (rows, cols)), shape=(len(p), 3 * len(x))).tocsr()
    result = lsmr(j.T, f.ravel(), atol=atol, btol=atol, maxiter=maxiter)
    projected = f.ravel() - j.T @ result[0]
    residual = np.asarray(j @ projected)
    residual_max = float(np.max(np.abs(residual)))
    threshold = max(1e-8, 1e-8 * float(np.max(np.abs(f))))
    if not np.isfinite(projected).all() or residual_max > threshold:
        raise RuntimeError(
            f"constraint projection did not converge: residual={residual_max}, threshold={threshold}, lsmr_stop={result[1]}"
        )
    return projected.reshape(x.shape), {
        "tangent_residual_max": residual_max,
        "iterations": int(result[2]),
        "lsmr_stop": int(result[1]),
    }


def convergence_confirmed(blocks):
    """Apply the preregistered strict thresholds to the last three blocks."""
    if len(blocks) < 3:
        return False
    return all(
        b.get("iterations") == 100
        and b.get("stop_reason") == "converged"
        and np.isfinite(
            [
                b["energy_change_kJ_mol"],
                b["projected_force_max"],
                b["max_displacement_angstrom"],
            ]
        ).all()
        and abs(b["energy_change_kJ_mol"]) < 0.1
        and b["projected_force_max"] < 10.0
        and b["max_displacement_angstrom"] < 0.01
        for b in blocks[-3:]
    )
