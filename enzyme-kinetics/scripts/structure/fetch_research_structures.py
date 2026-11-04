"""Fetch primary RCSB inputs for the registered P00947 study; no implicit downloads."""

import argparse
from pathlib import Path

import requests

from enzymekinetics_prediction.paths import DATA_DIR
from enzymekinetics_prediction.physics.artifacts import digest_file, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DATA_DIR / "research/p00947")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    sources = []
    for code, kind, suffix in [
        ("8CHO", "assembly", "pdb1"),
        ("1QJG", "assembly", "pdb1"),
        ("1OHS", "asymmetric", "pdb"),
        ("3NHX", "asymmetric", "pdb"),
    ]:
        url = f"https://files.rcsb.org/download/{code}.{suffix}"
        path = args.output / f"{code}-{kind}.pdb"
        if not path.exists():
            response = requests.get(url, timeout=60)
            response.raise_for_status()
            if b"ATOM" not in response.content:
                raise ValueError("RCSB response does not contain coordinates")
            path.write_bytes(response.content)
        sources.append({"url": url, "file": path.name, "sha256": digest_file(path)})
    write_json(
        args.output / "required-inputs.json",
        {
            "sources": sources,
            "scope": "8CHO WT assembly; 1QJG D38N equilenin; 1OHS Y14F/D38N saturated analog, not native substrate",
        },
    )
    print(args.output)


if __name__ == "__main__":
    main()
