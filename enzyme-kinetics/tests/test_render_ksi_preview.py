import importlib.util
from pathlib import Path

import pytest

pytest.importorskip("imageio")
spec = importlib.util.spec_from_file_location("render_ksi_preview", Path(__file__).resolve().parents[1] / "scripts/render_ksi_preview.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_missing_complex_is_rejected():
    with pytest.raises((KeyError, ValueError)):
        module.collect({0: {}})
