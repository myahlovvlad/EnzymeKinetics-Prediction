"""Check Python imports and report optional external executables."""

import argparse
import importlib
import importlib.metadata
import platform
import shutil
import sys

PROFILES = {
    "base": [
        ("numpy", "numpy"),
        ("scipy", "scipy"),
        ("pandas", "pandas"),
        ("matplotlib", "matplotlib"),
        ("sklearn", "scikit-learn"),
        ("yaml", "PyYAML"),
        ("PIL", "Pillow"),
    ],
    "chemistry": [
        ("rdkit", "rdkit"),
        ("Bio", "biopython"),
        ("requests", "requests"),
        ("tqdm", "tqdm"),
        ("sympy", "sympy"),
    ],
    "quantum": [("rdkit", "rdkit"), ("pyscf", "pyscf")],
    "dashboard": [
        ("streamlit", "streamlit"),
        ("plotly", "plotly"),
        ("pyarrow", "pyarrow"),
    ],
    "ml": [("torch", "torch")],
    "notebooks": [("jupyterlab", "jupyterlab"), ("ipykernel", "ipykernel")],
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", action="append", choices=PROFILES, default=[])
    args = parser.parse_args()
    print(f"Python {platform.python_version()}: {sys.executable}")
    failures = 0
    checked = set()
    for profile in ["base", *args.profile]:
        for module, distribution in PROFILES[profile]:
            if module in checked:
                continue
            checked.add(module)
            try:
                importlib.import_module(module)
                print(f"OK {distribution}=={importlib.metadata.version(distribution)}")
            except Exception as exc:  # noqa: BLE001 - report native import failures too
                print(f"FAIL {distribution}: {exc}")
                failures += 1
    for executable in ["java", "xtb", "apbs", "pdb2pqr", "vina"]:
        print(
            f"External {executable}: {shutil.which(executable) or 'not found (optional)'}"
        )
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
