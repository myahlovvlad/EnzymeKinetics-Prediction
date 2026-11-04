"""Periodic side-chain and backbone geometry for a GLU residue.

The high-force backbone C is not the side-chain carboxyl CD. Coordinates and
orthorhombic box lengths are in nm; angles are in degrees.
"""

import numpy as np

from enzymekinetics_prediction.physics.geometry import minimum_image


def image(v, box):
    """Compatibility wrapper for the canonical orthorhombic minimum image."""
    return minimum_image(v, box)


def angle(a, b):
    denominator = np.linalg.norm(a) * np.linalg.norm(b)
    if denominator < 1e-15:
        raise ValueError("degenerate angle")
    return float(np.degrees(np.arccos(np.clip(np.dot(a, b) / denominator, -1, 1))))


def geometry(atoms, box):
    """Return backbone C–O and carboxyl CD–OE1/OE2 distances, angle and improper.

    Improper is the signed dihedral CG–CD–OE1–OE2; planarity is checked
    relative to 0 or 180 degrees. All bond vectors use minimum image.
    """
    required = {"C", "O", "CG", "CD", "OE1", "OE2"}
    if not required <= atoms.keys():
        raise ValueError(f"missing atoms: {sorted(required - atoms.keys())}")
    x = {name: np.asarray(atoms[name], dtype=float) for name in required}
    if any(value.shape != (3,) for value in x.values()):
        raise ValueError("coordinates must have shape (3,)")
    cg = image(x["CG"] - x["CD"], box)
    o1 = image(x["OE1"] - x["CD"], box)
    o2 = image(x["OE2"] - x["CD"], box)
    axis = o1 / np.linalg.norm(o1) if np.linalg.norm(o1) > 1e-15 else None
    if axis is None:
        raise ValueError("degenerate carboxyl bond")
    left = cg - np.dot(cg, axis) * axis
    right = (o2 - o1) - np.dot(o2 - o1, axis) * axis
    if np.linalg.norm(left) < 1e-15 or np.linalg.norm(right) < 1e-15:
        raise ValueError("degenerate improper")
    improper = float(
        np.degrees(np.arctan2(np.dot(np.cross(axis, left), right), np.dot(left, right)))
    )
    return {
        "backbone_C_O_nm": float(np.linalg.norm(image(x["O"] - x["C"], box))),
        "sidechain_CD_OE1_nm": float(np.linalg.norm(o1)),
        "sidechain_CD_OE2_nm": float(np.linalg.norm(o2)),
        "sidechain_OE1_CD_OE2_deg": angle(o1, o2),
        "sidechain_CG_CD_OE1_OE2_improper_deg": improper,
    }
