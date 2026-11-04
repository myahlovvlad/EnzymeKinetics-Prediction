"""Sequence pI/charge estimates and structure-based pKa output with explicit scope."""

from functools import lru_cache
from pathlib import Path
import re

from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint
import numpy as np
from scipy.optimize import brentq

from .artifacts import artifact, digest_file


@lru_cache(maxsize=32768)
def sequence_titration(sequence, ph=7.0):
    if not sequence or any(a not in "ACDEFGHIKLMNPQRSTVWY" for a in sequence):
        return (float("nan"), float("nan"), float("nan"))
    if not np.isfinite(ph) or not 0 <= ph <= 14:
        return (
            float(brentq(IsoelectricPoint(sequence).charge_at_pH, 0, 14)),
            float("nan"),
            float("nan"),
        )
    model = IsoelectricPoint(sequence)
    pi = float(brentq(model.charge_at_pH, 0, 14))
    return pi, float(model.charge_at_pH(ph)), float(ph - pi)


def propka_summary(log):
    """Parse summary, not determinant pKa shifts; require unambiguous chain/residue."""
    text = Path(log).read_text()
    # PROPKA places the complete summary after this heading in PDB2PQR logs.
    section = text.split("SUMMARY OF THIS PREDICTION")
    if len(section) < 2:
        raise ValueError("PROPKA summary not found")
    records = []
    for line in section[-1].splitlines():
        match = re.match(
            r"^\s*(ASP|GLU|HIS|CYS|TYR|LYS|ARG|N\+|C-)\s+(\d+)\s+(\S)\s+([-\d.]+)\s+([-\d.]+)",
            line,
        )
        if match:
            records.append(
                {
                    "resname": match[1],
                    "resid": int(match[2]),
                    "chain": match[3],
                    "pka_predicted": float(match[4]),
                    "model_pka": float(match[5]),
                }
            )
    pi = re.search(r"The pI is\s+([\d.]+) \(folded\) and\s+([\d.]+)", text)
    if not records:
        raise ValueError("No parseable PROPKA residues")
    return artifact(
        "PROPKA-log-summary",
        {"log_sha256": digest_file(log)},
        {},
        {
            "residue_pka": records,
            "folded_pi": float(pi[1]) if pi else None,
            "unfolded_pi": float(pi[2]) if pi else None,
        },
        {"residue_pka": "pH units", "folded_pi": "pH units", "unfolded_pi": "pH units"},
        limitations=[
            "Empirical structure-based predictions, not experimental microscopic pKa",
            "pI is a global net-charge property; it does not replace catalytic-site pKa",
            "No temperature-dependent pKa shift is inferred without ionization enthalpies",
        ],
    )


def ionized_fraction(ph, pka, acid=True):
    if not np.isfinite(ph) or not np.isfinite(pka):
        raise ValueError("Finite pH and pKa required")
    return float(
        1 / (1 + 10.0 ** np.clip((pka - ph) if acid else (ph - pka), -300, 300))
    )


def water_dielectric(temperature_K):
    """Malmberg/Maryott 1956 pure-water fit, 0..100 C; not buffer calibration."""
    celsius = temperature_K - 273.15
    if not np.isfinite(celsius) or not 0 <= celsius <= 100:
        raise ValueError("Pure-water dielectric fit is limited to 273.15..373.15 K")
    return float(
        87.740 - 0.40008 * celsius + 9.398e-4 * celsius**2 - 1.410e-6 * celsius**3
    )
