"""Refined factorial with validated recovery of completed native APBS files."""

from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import itertools
import json
import numpy as np
import pandas as pd
from enzymekinetics_prediction.physics.electrostatics import run_apbs, sample_dx
from enzymekinetics_prediction.physics.artifacts import (
    artifact,
    write_json,
    digest_file,
)
from enzymekinetics_prediction.research.journal import Journal


def run_refined(directory, workers=6, threads=4):
    p = Path(directory)
    j = Journal(p)
    x = json.loads((p / "H2-conditions-factorial.json").read_text())
    probe = x["probe_A"]
    with j.run(
        "H2e-refined-factorial-resume",
        {
            "grid": 129,
            "focus_length_A": 20,
            "workers": workers,
            "threads_per_job": threads,
            "parser": "vectorized numpy.fromstring",
        },
    ):

        def calculate(parameters):
            ph, temp, salt = parameters
            location = p / f"factorial-focused/ph{ph}-T{temp}-I{salt}"
            pqr = p / f"ph-{ph}/protein.pqr"
            if (location / "result.json").exists():
                r = json.loads((location / "result.json").read_text())
                settings = r["settings"]
                if (
                    settings.get("temperature_K") != temp
                    or settings.get("ionic_strength_M") != salt
                    or settings.get("grid") != [129] * 3
                    or not np.allclose(
                        np.atleast_2d(settings.get("query_A", [])), np.atleast_2d(probe)
                    )
                ):
                    raise ValueError("Cached APBS conditions or probe changed")
                if r["inputs"]["pqr_sha256"] != digest_file(pqr) or r["inputs"][
                    "dx_sha256"
                ] != digest_file(location / "potential.dx"):
                    raise ValueError("Cached APBS inputs changed")
            elif (location / "potential.dx").exists() and "Thanks for using APBS!" in (
                location / "apbs.log"
            ).read_text():
                script = (location / "apbs.in").read_text()
                required = [
                    "dime 129 129 129",
                    f"temp {temp}",
                    f"ion charge 1 conc {salt}",
                    "fglen 20.0 20.0 20.0",
                    "sdie 78.5",
                ]
                if not all(value in script for value in required) or digest_file(
                    location / "input.pqr"
                ) != digest_file(pqr):
                    raise ValueError("Cannot recover APBS with inconsistent inputs")
                phi, field = sample_dx(location / "potential.dx", probe, temp)
                r = artifact(
                    "APBS-recovered-completed-DX",
                    {
                        "pqr_sha256": digest_file(pqr),
                        "dx_sha256": digest_file(location / "potential.dx"),
                        "apbs_input_sha256": digest_file(location / "apbs.in"),
                        "apbs_log_sha256": digest_file(location / "apbs.log"),
                    },
                    {
                        "temperature_K": temp,
                        "ionic_strength_M": salt,
                        "grid": [129] * 3,
                        "focus_length_A": 20,
                        "query_A": probe,
                        "recovery": "Complete APBS log and expected input parameters verified; native solver not repeated",
                    },
                    {"potential_V": phi.tolist(), "field_V_per_A": field.tolist()},
                    {"potential_V": "Volt", "field_V_per_A": "Volt/Angstrom"},
                )
                write_json(location / "result.json", r)
            else:
                r = run_apbs(
                    pqr,
                    location,
                    probe,
                    temperature_K=temp,
                    ionic_strength_M=salt,
                    threads=threads,
                    grid_points=129,
                    focus_length_A=20,
                    timeout=300,
                )
                write_json(location / "result.json", r)
            f = np.array(r["values"]["field_V_per_A"][0])
            return {
                "ph": ph,
                "temperature_K": temp,
                "ionic_strength_M": salt,
                "potential_V": r["values"]["potential_V"][0],
                "field_norm_V_per_A": float(np.linalg.norm(f)),
                **{f"field_{a}_V_per_A": float(f[i]) for i, a in enumerate("xyz")},
            }

        with ThreadPoolExecutor(workers) as pool:
            rows = list(
                pool.map(
                    calculate,
                    itertools.product(
                        [5, 7, 9], [288.15, 298.15, 310.15], [0, 0.15, 0.3]
                    ),
                )
            )
        pd.DataFrame(rows).to_csv(p / "conditions-factorial-focused.csv", index=False)
        checks = []
        for ph, temp, salt in [(5, 288.15, 0), (9, 310.15, 0.3)]:
            r = run_apbs(
                p / f"ph-{ph}/protein.pqr",
                p / f"focused-check-ph{ph}",
                probe,
                temperature_K=temp,
                ionic_strength_M=salt,
                threads=threads,
                grid_points=97,
                focus_length_A=20,
            )
            write_json(p / f"focused-check-ph{ph}/result.json", r)
            fine = next(
                row
                for row in rows
                if row["ph"] == ph
                and row["temperature_K"] == temp
                and row["ionic_strength_M"] == salt
            )
            f = np.array([fine[f"field_{a}_V_per_A"] for a in "xyz"])
            c = np.array(r["values"]["field_V_per_A"][0])
            change = float(np.linalg.norm(f - c) / np.linalg.norm(f))
            checks.append(
                {
                    "ph": ph,
                    "temperature_K": temp,
                    "ionic_strength_M": salt,
                    "relative_change_97_to_129": change,
                    "passes_10percent": change <= 0.1,
                }
            )
        j.result(
            "H2e-refined-factorial",
            {
                "runs": 27,
                "checks_at_extremes": checks,
                "all_checked_conditions_pass": all(
                    c["passes_10percent"] for c in checks
                ),
                "limitations": [
                    "Grid checks at three conditions, not experimental validation of pKa or dielectric model",
                    "Empirical fixed protonation and apo structure, not ES or constant-pH sampling",
                ],
            },
        )
