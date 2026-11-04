"""Compatibility entry point: Ridge and PLS with grouped validation and saved artifacts."""

import sys
from enzymekinetics_prediction.learning.training import main

if __name__ == "__main__":
    if "--models" not in sys.argv:
        sys.argv.extend(["--models", "mean", "ridge", "pls"])
    main()
