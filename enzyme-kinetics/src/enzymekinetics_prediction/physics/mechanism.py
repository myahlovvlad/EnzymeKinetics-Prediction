"""Explicit two-chemical-step KSI mechanism and stationary-point validation.

A reaction graph states what must be calculated. Missing transition-state
coordinates and barriers remain missing, rather than being inferred from kcat.
"""

import numpy as np

KSI_MECHANISM = {
    "system": "Comamonas testosteroni P00947; author residue numbering",
    "states": [
        {
            "id": "E_S",
            "kind": "separated",
            "description": "Free enzyme and Delta5 steroid",
        },
        {
            "id": "ES",
            "kind": "minimum",
            "description": "Bound Delta5 steroid; Asp38 base, Tyr14/Asp99 oxyanion donors",
        },
        {
            "id": "TS1",
            "kind": "saddle_index_1",
            "description": "C4-beta H abstraction by Asp38; C4-H breaking and Asp38 O-H formation",
        },
        {
            "id": "EI",
            "kind": "minimum",
            "description": "Dienolate intermediate; Asp38 protonated; charge stabilized by Tyr14/Asp99",
        },
        {
            "id": "TS2",
            "kind": "saddle_index_1",
            "description": "Asp38-H donation to C6; conjugated Delta4 product formation",
        },
        {
            "id": "EP",
            "kind": "minimum",
            "description": "Bound Delta4 steroid; Asp38 returns to base state",
        },
        {
            "id": "E_P",
            "kind": "separated",
            "description": "Product released, enzyme regenerated",
        },
    ],
    "edges": [
        ["E_S", "ES"],
        ["ES", "TS1"],
        ["TS1", "EI"],
        ["EI", "TS2"],
        ["TS2", "EP"],
        ["EP", "E_P"],
    ],
    "coordinates": {
        "q1": "r(C4,H_transfer)-r(O_Asp38,H_transfer)",
        "q2": "r(O_Asp38,H_transfer)-r(C6,H_transfer)",
        "additional": "O3--Tyr14/Asp99 distances, closure, solvent polarization and protonation microstates",
    },
    "required_checks": [
        "Same atoms, total charge and spin along each path",
        "Converged gradient for full unfrozen degrees of freedom",
        "Exactly one relevant imaginary mode at each full TS; minima have no negative internal modes",
        "Forward/backward path connectivity to specified minima, not Hessian index alone",
        "Thermal/solvent/ensemble free-energy corrections before rate predictions",
    ],
    "source": "https://pubs.acs.org/doi/10.1021/jacs.6b04796",
    "computed_full_TS": False,
}


def classify_stationary(
    gradient, hessian, gradient_tolerance=1e-3, curvature_tolerance=1e-4
):
    g = np.asarray(gradient, float).reshape(-1)
    h = np.asarray(hessian, float)
    if (
        h.shape != (len(g), len(g))
        or not np.isfinite(g).all()
        or not np.isfinite(h).all()
    ):
        raise ValueError("Invalid gradient/Hessian")
    if not np.allclose(h, h.T, atol=1e-6):
        raise ValueError("Hessian must be symmetric")
    eig = np.linalg.eigvalsh(h)
    negative = int((eig < -curvature_tolerance).sum())
    stationary = bool(np.max(np.abs(g)) <= gradient_tolerance)
    classification = (
        (
            "minimum"
            if negative == 0
            else "index_1_saddle"
            if negative == 1
            else "higher_index_saddle"
        )
        if stationary
        else "nonstationary"
    )
    return {
        "classification": classification,
        "stationary": stationary,
        "negative_modes": negative,
        "eigenvalues": eig.tolist(),
        "gradient_max_abs": float(np.max(np.abs(g))),
        "limitation": "Classification applies only to supplied coordinates/subspace; TS path connectivity requires separate verification",
    }


def numerical_hessian(gradient_function, coordinates, step=1e-3):
    x = np.asarray(coordinates, float)
    flat = x.reshape(-1)
    if not np.isfinite(x).all() or step <= 0:
        raise ValueError("Invalid finite difference geometry or step")
    columns = []
    for i in range(len(flat)):
        delta = np.zeros_like(flat)
        delta[i] = step
        a = np.asarray(gradient_function((flat + delta).reshape(x.shape))).reshape(-1)
        b = np.asarray(gradient_function((flat - delta).reshape(x.shape))).reshape(-1)
        columns.append((a - b) / (2 * step))
    h = np.column_stack(columns)
    return (h + h.T) / 2
