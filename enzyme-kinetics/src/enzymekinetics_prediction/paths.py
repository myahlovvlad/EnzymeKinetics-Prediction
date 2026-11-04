"""Shared paths for the local research checkout, independent of the current directory.

Set ENZYME_KINETICS_PROJECT_DIR when using an installed wheel outside the checkout.
"""

import os
from pathlib import Path

PROJECT_DIR = (
    Path(
        os.environ.get(
            "ENZYME_KINETICS_PROJECT_DIR", Path(__file__).resolve().parents[2]
        )
    )
    .expanduser()
    .resolve()
)
REPOSITORY_DIR = PROJECT_DIR.parent
DATA_DIR = PROJECT_DIR / "data"
OUTPUT_DIR = PROJECT_DIR / "outputs"
MANIFEST_DIR = DATA_DIR / "manifests"
VENDOR_DIR = PROJECT_DIR / "vendor"
