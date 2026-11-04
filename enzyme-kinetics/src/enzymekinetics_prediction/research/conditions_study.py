"""Registered static active-site/pH/temperature/salt study on WT tKSI assembly."""

import argparse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import itertools
import json

import numpy as np
import pandas as pd

from enzymekinetics_prediction.physics.artifacts import write_json, digest_file
from enzymekinetics_prediction.physics.conditions import (
    propka_summary,
    sequence_titration,
    ionized_fraction,
)
from enzymekinetics_prediction.physics.electrostatics import (
    prepare_pqr,
    run_apbs,
    read_pqr,
)
from enzymekinetics_prediction.physics.sites import flatten_assembly, site_report
from enzymekinetics_prediction.paths import DATA_DIR, OUTPUT_DIR
from .journal import Journal


def run(directory, workers=6):
    directory = Path(directory)
    journal = Journal(directory)
    registration = json.loads((directory / "registration.json").read_text())
    if not registration["registered_before_scientific_runs"]:
        raise ValueError("Register hypotheses first")
    source = DATA_DIR / "research/p00947"
    normalized = directory / "8CHO-wt-dimer.pdb"
    with journal.run(
        "H1-static-sites", {}, {"8CHO": digest_file(source / "8CHO-assembly.pdb")}
    ):
        mapping = flatten_assembly(source / "8CHO-assembly.pdb", normalized)
        write_json(directory / "assembly.json", mapping)
        structures = {
            "WT_experimental": site_report(normalized, "A"),
            "AF_monomer": site_report(DATA_DIR / "alphafold_pdb/P00947.pdb", "A"),
            "D38N_equilenin_analog": site_report(
                source / "1QJG-assembly.pdb", "C", "EQU"
            ),
            "Y14F_D38N_5alpha_androstanedione_analog": site_report(
                source / "1OHS-asymmetric.pdb", "A", "5SD"
            ),
        }
        journal.result("H1-static-sites", structures)
    sequence = "MNTPEHMTAVVQRYVAALNAGDLDGIVALFADDATVEDPVGSEPRSGTAAIREFYANSLKLPLAVELTQEVRAVANEAAFAFIVSFEYQGRKTVVAPIDHFRFNGAGKVVSMRALFGEKNIHAGA"
    # Catalytic oxyanion-hole probe is a geometric reference point, not a bound ligand atom.
    from Bio.PDB import PDBParser

    chain = PDBParser(QUIET=True).get_structure("wt", normalized)[0]["A"]
    probe = (chain[14]["OH"].coord + chain[99]["OD2"].coord) / 2
    with journal.run(
        "H2-conditions-factorial",
        {
            "ph": [5, 7, 9],
            "temperature_K": [288.15, 298.15, 310.15],
            "ionic_strength_M": [0, 0.15, 0.3],
            "probe_A": probe.tolist(),
            "workers": workers,
        },
    ):
        pqrs = {}
        for ph in [5, 7, 9]:
            pqr = prepare_pqr(normalized, directory / f"ph-{ph}/protein.pqr", ph)
            pqrs[ph] = pqr
            write_json(
                directory / f"ph-{ph}/propka.json",
                propka_summary(pqr.with_suffix(".log")),
            )

        def calculate(condition):
            ph, temp, salt = condition
            location = directory / f"factorial/ph{ph}-T{temp}-I{salt}"
            result = run_apbs(
                pqrs[ph],
                location,
                probe,
                temperature_K=temp,
                ionic_strength_M=salt,
                threads=1,
            )
            write_json(location / "result.json", result)
            _, charge, _ = read_pqr(pqrs[ph])
            return {
                "ph": ph,
                "temperature_K": temp,
                "ionic_strength_M": salt,
                "potential_V": result["values"]["potential_V"][0],
                "field_norm_V_per_A": float(
                    np.linalg.norm(result["values"]["field_V_per_A"][0])
                ),
                "protein_net_charge_e": float(charge.sum()),
                "sequence_pi": sequence_titration(sequence, ph)[0],
                **{
                    f"field_{axis}_V_per_A": result["values"]["field_V_per_A"][0][i]
                    for i, axis in enumerate("xyz")
                },
            }

        results = list(
            ThreadPoolExecutor(workers).map(
                calculate,
                itertools.product([5, 7, 9], [288.15, 298.15, 310.15], [0, 0.15, 0.3]),
            )
        )
        pd.DataFrame(results).to_csv(
            directory / "conditions-factorial.csv", index=False
        )
        grids = []
        for points in [65, 97, 129]:
            r = run_apbs(
                pqrs[7],
                directory / f"grid-{points}",
                probe,
                grid_points=points,
                threads=1,
            )
            write_json(directory / f"grid-{points}/result.json", r)
            grids.append(
                {
                    "grid": points,
                    "potential_V": r["values"]["potential_V"][0],
                    "field_V_per_A": r["values"]["field_V_per_A"][0],
                }
            )
        final = np.array(grids[-1]["field_V_per_A"])
        previous = np.array(grids[-2]["field_V_per_A"])
        relative = float(
            np.linalg.norm(final - previous) / max(np.linalg.norm(final), 1e-12)
        )
        pkas = json.loads((directory / "ph-7/propka.json").read_text())
        sites = [
            {**r, "deprotonated_fraction_ph7": ionized_fraction(7, r["pka_predicted"])}
            for r in pkas["values"]["residue_pka"]
            if r["resid"] in [14, 38, 99] and r["chain"] == "A"
        ]
        journal.result(
            "H2-conditions-factorial",
            {
                "runs": len(results),
                "probe_A": probe.tolist(),
                "grid_comparison": grids,
                "relative_field_difference_97_to_129": relative,
                "grid_converged_at_10percent": relative <= 0.1,
                "site_pka_predictions": sites,
                "folded_pi": pkas["values"]["folded_pi"],
                "limitations": [
                    "Deterministic sensitivity, not statistical replicates; no p-values",
                    "No ligand at probe and no inferred activation barrier",
                    "pKa not temperature corrected; static protonation at nominal pH, not constant-pH ensemble",
                    "Dielectric solvent fixed to 78.5; temperature-dependent solvent response not fitted",
                ],
            },
        )
    journal.event("study", "static-stages-complete")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path, default=OUTPUT_DIR / "experiments/active-site-v3"
    )
    parser.add_argument("--workers", type=int, default=6)
    args = parser.parse_args()
    if not 1 <= args.workers <= 8:
        raise ValueError("Use 1..8 workers to bound APBS memory")
    run(args.output, args.workers)


if __name__ == "__main__":
    main()
