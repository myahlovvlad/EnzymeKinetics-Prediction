"""TST and an explicit E+S <-> ES <-> EP -> E+P initial-rate mechanism."""

from dataclasses import dataclass
import math

import numpy as np
from scipy.constants import R, Boltzmann, Planck, speed_of_light


def positive(value, name, allow_zero=False):
    if not math.isfinite(value) or (value < 0 if allow_zero else value <= 0):
        raise ValueError(f"Invalid {name}: {value}")
    return float(value)


def eyring(delta_g_kj_mol, temperature_K, transmission=1.0):
    positive(temperature_K, "temperature_K")
    positive(transmission, "transmission")
    if not math.isfinite(delta_g_kj_mol):
        raise ValueError("Nonfinite activation free energy")
    return (
        transmission
        * Boltzmann
        * temperature_K
        / Planck
        * math.exp(-1000 * delta_g_kj_mol / (R * temperature_K))
    )


def tunneling(imaginary_wavenumber_cm, temperature_K, method="wigner"):
    """High-T approximations; frequency must come from the relevant TS mode.

    Bell here denotes the inverted-parabolic-barrier thermal factor x/sin(x),
    valid above its crossover temperature. It is not a fitted universal factor.
    """
    positive(imaginary_wavenumber_cm, "absolute imaginary frequency")
    positive(temperature_K, "temperature_K")
    x = (
        Planck
        * speed_of_light
        * 100
        * imaginary_wavenumber_cm
        / (2 * Boltzmann * temperature_K)
    )
    if method == "wigner":
        if x >= math.pi:
            raise ValueError("Wigner approximation below crossover is not supported")
        return 1 + x * x / 6
    if method == "bell-parabolic":
        if x >= math.pi:
            raise ValueError("Parabolic Bell factor requires T above crossover")
        return x / math.sin(x)
    raise ValueError(f"Unknown tunneling method: {method}")


@dataclass(frozen=True)
class CycleRates:
    kon_M_inv_s: float
    koff_s: float
    chemistry_forward_s: float
    chemistry_reverse_s: float
    release_s: float

    def __post_init__(self):
        for name, value in vars(self).items():
            positive(value, name, allow_zero=name in {"koff_s", "chemistry_reverse_s"})

    def parameters(self):
        a, b, c, d, e = (
            self.kon_M_inv_s,
            self.koff_s,
            self.chemistry_forward_s,
            self.chemistry_reverse_s,
            self.release_s,
        )
        kcat = c * e / (c + d + e)
        km = (b * (d + e) + c * e) / (a * (c + d + e))
        return {"kcat_s_inv": kcat, "Km_M": km, "kcat_over_Km_M_inv_s": kcat / km}

    def stationary(self, substrate_M, enzyme_total_M=1e-6):
        positive(substrate_M, "substrate_M", allow_zero=True)
        positive(enzyme_total_M, "enzyme_total_M")
        a = self.kon_M_inv_s * substrate_M
        b, c, d, e = (
            self.koff_s,
            self.chemistry_forward_s,
            self.chemistry_reverse_s,
            self.release_s,
        )
        generator = np.array([[-a, a, 0], [b, -b - c, c], [e, d, -d - e]], dtype=float)
        system = generator.T.copy()
        system[-1] = 1
        probabilities = np.linalg.solve(system, [0, 0, 1])
        if np.any(probabilities < -1e-10):
            raise ValueError("Invalid stationary enzyme probabilities")
        return {
            "E_ES_EP_fractions": probabilities.tolist(),
            "velocity_M_per_s": float(enzyme_total_M * e * probabilities[2]),
        }
