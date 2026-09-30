"""
Authors: Andre L. Belem (https://github.com/andrebelem) and F.R.I.D.A.Y. (https://observatoriooceanografico.org/people/friday-bot)

Shared fixtures. Tests marked ``large`` read the original files in data/ read-only
and skip cleanly when those files are absent. No test materializes the EDS cube.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

import inspect_bruker

ROOT = Path(__file__).resolve().parent.parent
# Child processes (``python -m microxrf_to_netcdf``) must find the package without an installation.
os.environ["PYTHONPATH"] = os.pathsep.join(filter(None, [str(ROOT / "src"), os.environ.get("PYTHONPATH")]))

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
STEM = "GRF17A_9-29cm_slab3_Elemental_map"
BCF_PATH = DATA_DIR / f"{STEM}.bcf"
RTX_PATH = DATA_DIR / f"{STEM}.rtx"


def _require(path: Path) -> Path:
    if not path.is_file():
        pytest.skip(f"original acquisition file not available: {path.name}")
    return path


@pytest.fixture(scope="session")
def bcf_path() -> Path:
    return _require(BCF_PATH)


@pytest.fixture(scope="session")
def rtx_path() -> Path:
    return _require(RTX_PATH)


@pytest.fixture(scope="session")
def rtx_report(rtx_path: Path) -> dict:
    """Incremental (1 MiB chunks) XML analysis of the RTX; payload text is never retained."""
    return inspect_bruker.inspect_rtx(rtx_path)


@pytest.fixture(scope="session")
def bcf_result(bcf_path: Path) -> dict:
    """Lazy RosettaSciIO read requested exactly as inspect_bruker.py does."""
    result = inspect_bruker.inspect_bcf(bcf_path)
    assert result["error"] is None, result["error"]
    return result


@pytest.fixture(scope="session")
def video_signal(bcf_result: dict) -> dict:
    return next(s for s in bcf_result["signals"] if s["metadata"]["General"]["title"] == "Video")


@pytest.fixture(scope="session")
def edx_signal(bcf_result: dict) -> dict:
    return next(s for s in bcf_result["signals"] if s["metadata"]["General"]["title"] == "EDX")


@pytest.fixture(scope="session")
def rtx_decoded(rtx_path: Path) -> dict:
    """Streaming decode of the RTX payload (rtx_payload) keeping only what the tests compare.

    Kept in memory: the 15 planes of the ``Mapdaten`` image (about 26 MB) and the 3 planes of the
    second ``Video Mosaic`` (about 17 MB). The rest of the 63 MB payload is never held.
    """
    import numpy as np

    from microxrf_to_netcdf import rtx_payload

    kept: dict[tuple[int, int], object] = {}

    def sink(image: int, plane: int, raw: bytes) -> None:
        if image == 2:
            kept[(image, plane)] = np.frombuffer(raw, "<u2").reshape(240, 1800).copy()
        elif image == 1:
            kept[(image, plane)] = np.frombuffer(raw, "u1").reshape(948, 6000).copy()

    report = rtx_payload.inspect_rtx_payload(rtx_path, sink)
    assert report["error"] is None, report["error"]
    return {"report": report, "planes": kept}
