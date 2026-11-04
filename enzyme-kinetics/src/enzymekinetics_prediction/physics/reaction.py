"""Mass-action E+S <-> ES <-> EP <-> E+P with explicit conservation laws.

Rate constants are inputs, not inferred from molecular descriptors. Activities
are approximated by ideal dilute concentrations. This is a progress-curve
simulation of a stated kinetic mechanism, not an atomistic bond-breaking model.
"""

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.constants import R

from .artifacts import artifact, write_json
from .kinetics import CycleRates, positive

SPECIES = ["E", "ES", "EP", "S", "P"]
# Columns are substrate binding, chemistry, product release.
STOICHIOMETRY = np.array(
    [[-1, 0, 1], [1, -1, 0], [0, 1, -1], [-1, 0, 0], [0, 0, 1]], float
)


@dataclass(frozen=True)
class ReversibleCycle:
    kon_substrate_M_inv_s: float
    koff_substrate_s: float
    chemistry_forward_s: float
    chemistry_reverse_s: float
    kon_product_M_inv_s: float
    release_product_s: float

    def __post_init__(self):
        for key, value in asdict(self).items():
            positive(value, key)

    def fluxes(self, concentrations_M):
        e, es, ep, s, p = concentrations_M
        return np.array(
            [
                self.kon_substrate_M_inv_s * e * s - self.koff_substrate_s * es,
                self.chemistry_forward_s * es - self.chemistry_reverse_s * ep,
                self.release_product_s * ep - self.kon_product_M_inv_s * e * p,
            ]
        )

    def initial_rate_parameters(self):
        return CycleRates(
            self.kon_substrate_M_inv_s,
            self.koff_substrate_s,
            self.chemistry_forward_s,
            self.chemistry_reverse_s,
            self.release_product_s,
        ).parameters()

    def equilibrium_constant(self):
        return (
            self.kon_substrate_M_inv_s
            * self.chemistry_forward_s
            * self.release_product_s
            / (
                self.koff_substrate_s
                * self.chemistry_reverse_s
                * self.kon_product_M_inv_s
            )
        )

    def check_thermodynamics(
        self, standard_reaction_delta_g_kj_mol, temperature_K, tolerance_kj_mol=0.01
    ):
        positive(temperature_K, "temperature_K")
        positive(tolerance_kj_mol, "tolerance_kj_mol")
        implied = -R * temperature_K * np.log(self.equilibrium_constant()) / 1000
        if (
            not np.isfinite(standard_reaction_delta_g_kj_mol)
            or abs(implied - standard_reaction_delta_g_kj_mol) > tolerance_kj_mol
        ):
            raise ValueError(
                "Rate constants violate the supplied equilibrium free energy/Haldane relation"
            )
        return float(implied)

    def simulate(self, initial_M, times_s, rtol=1e-8, atol_M=1e-14):
        y = np.asarray(initial_M, float)
        times = np.asarray(times_s, float)
        if y.shape != (5,) or not np.isfinite(y).all() or (y < 0).any():
            raise ValueError(
                "Provide five nonnegative finite initial concentrations E,ES,EP,S,P in M"
            )
        if (
            times.ndim != 1
            or len(times) < 2
            or not np.isfinite(times).all()
            or times[0] != 0
            or (np.diff(times) <= 0).any()
        ):
            raise ValueError("Times must be strictly increasing from zero, in seconds")
        positive(rtol, "rtol")
        positive(atol_M, "atol_M")
        solution = solve_ivp(
            lambda t, c: STOICHIOMETRY @ self.fluxes(c),
            (0, times[-1]),
            y,
            method="BDF",
            t_eval=times,
            rtol=rtol,
            atol=atol_M,
        )
        if not solution.success or not np.isfinite(solution.y).all():
            raise RuntimeError(f"Kinetic integration failed: {solution.message}")
        if solution.y.min() < -10 * atol_M:
            raise RuntimeError("Negative concentrations exceed solver tolerance")
        enzyme = solution.y[:3].sum(axis=0)
        substrate = solution.y[[1, 2, 3, 4]].sum(axis=0)
        for actual, total in [
            (enzyme, y[:3].sum()),
            (substrate, y[[1, 2, 3, 4]].sum()),
        ]:
            if not np.allclose(actual, total, rtol=rtol * 10, atol=atol_M * 10):
                raise RuntimeError("Mass conservation failed")
        frame = pd.DataFrame(solution.y.T, columns=[s + "_M" for s in SPECIES])
        frame.insert(0, "time_s", times)
        frame["product_flux_M_per_s"] = [self.fluxes(c)[2] for c in solution.y.T]
        evidence = artifact(
            "mass-action-reversible-cycle",
            {"rates": asdict(self), "initial_M": y.tolist()},
            {"rtol": rtol, "atol_M": atol_M, "times_s": times.tolist()},
            {
                "equilibrium_P_over_S": self.equilibrium_constant(),
                "enzyme_conservation_max_error_M": float(
                    np.max(np.abs(enzyme - enzyme[0]))
                ),
                "substrate_conservation_max_error_M": float(
                    np.max(np.abs(substrate - substrate[0]))
                ),
                **self.initial_rate_parameters(),
            },
            {
                "equilibrium_P_over_S": "dimensionless",
                "enzyme_conservation_max_error_M": "M",
                "substrate_conservation_max_error_M": "M",
                "kcat_s_inv": "s^-1",
                "Km_M": "M",
                "kcat_over_Km_M_inv_s": "M^-1 s^-1",
            },
            limitations=[
                "Ideal dilute mass-action concentrations; no spatial transport, changing pH or crowding",
                "Input rate constants require independent experimental or free-energy justification",
                "Reversible rates imply an equilibrium constant; independent thermodynamic data can reject inconsistent rates",
            ],
        )
        return frame, evidence


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    config = json.loads(args.config.read_text())
    allowed = {
        "rates",
        "initial_M",
        "duration_s",
        "frames",
        "provenance",
        "standard_reaction_delta_g_kj_mol",
        "temperature_K",
    }
    if set(config) - allowed or not config.get("provenance"):
        raise ValueError("Provide recognized fields and explicit rate provenance")
    if args.output.exists() and any(args.output.iterdir()):
        raise ValueError("Choose a new output directory")
    positive(config["duration_s"], "duration_s")
    if not isinstance(config["frames"], int) or not 2 <= config["frames"] <= 100000:
        raise ValueError("frames must be an integer from 2 to 100000")
    rates = ReversibleCycle(**config["rates"])
    if "standard_reaction_delta_g_kj_mol" in config:
        rates.check_thermodynamics(
            config["standard_reaction_delta_g_kj_mol"], config["temperature_K"]
        )
    frame, evidence = rates.simulate(
        config["initial_M"], np.linspace(0, config["duration_s"], config["frames"])
    )
    evidence["inputs"]["provenance"] = config["provenance"]
    # Provenance must participate in the immutable artifact identity.
    evidence = artifact(
        evidence["backend"],
        evidence["inputs"],
        evidence["settings"],
        evidence["values"],
        evidence["units"],
        limitations=evidence["limitations"],
    )
    args.output.mkdir(parents=True, exist_ok=True)
    frame.to_csv(args.output / "progress.csv", index=False)
    write_json(args.output / "reaction.json", evidence)
    print(json.dumps(evidence["values"]))


if __name__ == "__main__":
    main()
