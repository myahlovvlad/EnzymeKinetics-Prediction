"""Benchmark one 24-thread WT trajectory against three 8-thread replicas on this PC."""

import os
import sys
import subprocess
import time
import json
from pathlib import Path
from enzymekinetics_prediction.research.performance import benchmark
from enzymekinetics_prediction.physics.artifacts import write_json
from enzymekinetics_prediction.research.journal import Journal

import argparse
from enzymekinetics_prediction.paths import OUTPUT_DIR


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--study", type=Path, default=OUTPUT_DIR / "experiments/active-site-v3"
    )
    args = parser.parse_args()
    if len(os.sched_getaffinity(0)) < 24:
        raise ValueError(
            "This registered comparison requires 24 available logical CPUs; choose smaller configurations in research.performance"
        )
    root = args.study
    p = root / "prepared-wt-dimer"
    j = Journal(root)
    with j.run(
        "performance-WT-throughput",
        {"single_threads": [8, 12, 24], "concurrent": "3x8", "steps": 20, "repeats": 2},
    ):
        single = benchmark(
            p / "system.xml", p / "topology.pdb", [8, 12, 24], steps=20, repeats=2
        )
        output = args.study / "performance"
        output.mkdir(exist_ok=True)
        write_json(output / "wt-single.json", single)
        command = [
            sys.executable,
            "-m",
            "enzymekinetics_prediction.research.performance",
            "--system",
            str(p / "system.xml"),
            "--topology",
            str(p / "topology.pdb"),
            "--threads",
            "8",
            "--steps",
            "20",
            "--repeats",
            "2",
        ]
        start = time.perf_counter()
        children = []
        for i in range(3):
            children.append(
                subprocess.Popen(
                    command + ["--output", str(output / f"wt-concurrent-{i}.json")],
                    env={
                        **os.environ,
                        "OMP_NUM_THREADS": "1",
                        "OPENBLAS_NUM_THREADS": "1",
                    },
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                )
            )
        for child in children:
            if child.wait() != 0:
                raise RuntimeError("Concurrent benchmark failed")
        concurrent = [
            json.loads((output / f"wt-concurrent-{i}.json").read_text())
            for i in range(3)
        ]
        aggregate = sum(d["records"][0]["steps_per_s"] for d in concurrent)
        best = max(r["steps_per_s"] for r in single["records"])
        j.result(
            "performance-WT-throughput",
            {
                "single": single,
                "concurrent_3x8": concurrent,
                "aggregate_3x8_steps_per_s": aggregate,
                "best_single_steps_per_s": best,
                "concurrent_wall_seconds_including_startup": time.perf_counter()
                - start,
                "recommended_schedule": "3x8 concurrent"
                if aggregate > best
                else f"sequential replicas at {single['best_threads']} threads",
                "limitation": "Short throughput pilot; no claim of long-time scientific convergence",
            },
        )


if __name__ == "__main__":
    main()
