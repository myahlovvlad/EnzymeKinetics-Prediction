"""Run registered active-site experiments with explicit stages and a resource budget."""

import argparse
import json
from pathlib import Path

from enzymekinetics_prediction.paths import PROJECT_DIR, DATA_DIR
from .journal import Journal


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--registration",
        type=Path,
        default=PROJECT_DIR / "configs/active-site-hypotheses.json",
    )
    parser.add_argument(
        "--stage",
        choices=[
            "static",
            "refine",
            "prepare",
            "replicas",
            "torsions",
            "surface",
            "reanalyze",
        ],
        required=True,
    )
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--equilibration-steps", type=int, default=5000)
    parser.add_argument("--production-steps", type=int, default=5000)
    args = parser.parse_args()
    journal = Journal(args.output)
    journal.register(json.loads(args.registration.read_text()))
    if not 1 <= args.workers <= 8:
        raise ValueError("Use 1..8 APBS workers to bound memory")
    if args.stage == "static":
        if (args.output / "H2-conditions-factorial.json").exists():
            raise ValueError("Use a fresh study directory to repeat static experiments")
        from .conditions_study import run

        run(args.output, args.workers)
    elif args.stage == "refine":
        import os

        if args.workers * args.threads > len(os.sched_getaffinity(0)):
            raise ValueError("APBS jobs exceed the CPU affinity budget")
        from .refined_conditions import run_refined

        run_refined(args.output, args.workers, args.threads)
    elif args.stage == "prepare":
        from enzymekinetics_prediction.physics.dynamics import MDSettings, run_md

        result = run_md(
            args.output / "8CHO-wt-dimer.pdb",
            args.output / "prepared-wt-dimer",
            MDSettings(
                threads=args.threads,
                equilibration_steps=0,
                production_steps=4,
                report_interval=2,
                protonation_variants={"A:99": "ASH", "B:99": "ASH"},
            ),
        )
        journal.result("WT-dimer-preparation", result)
    elif args.stage == "replicas":
        import os

        if args.threads * 3 > len(os.sched_getaffinity(0)):
            raise ValueError("Three replicas exceed CPU affinity budget")
        from .ensemble import run_replicas

        run_replicas(
            args.output / "prepared-wt-dimer",
            args.output,
            args.threads,
            args.equilibration_steps,
            args.production_steps,
        )
    elif args.stage == "torsions":
        from .torsions import run

        run(args.output / "prepared-wt-dimer", args.output, args.threads)
    elif args.stage == "surface":
        from .proton_surface import run

        run(
            args.output / "8CHO-wt-dimer.pdb",
            DATA_DIR / "research/p00947/3NHX-asymmetric.pdb",
            args.output / "product-template-surface",
            grid=9,
        )
    elif args.stage == "reanalyze":
        from .ensemble import analyze_replica

        for folder in sorted(args.output.glob("replica-*")):
            result = analyze_replica(folder)
            journal.result(folder.name + "-analysis", result)


if __name__ == "__main__":
    main()
