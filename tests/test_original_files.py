"""
Authors: Andre L. Belem (https://github.com/andrebelem) and F.R.I.D.A.Y. (https://observatoriooceanografico.org/people/friday-bot)

Regression tests for the values recorded in FINDINGS.md. They read the original
acquisition files in data/ (read-only) and are skipped when the files are absent.
The (240, 1800, 4096) EDS cube is only inspected through its lazy Dask metadata
(shape, dtype, chunks); no test computes or converts its values.
"""

from __future__ import annotations

import dask.array as da
import pytest

from conftest import DATA_DIR
from inspect_bruker import inspect_png, normalize_date

pytestmark = pytest.mark.large


# ---------------------------------------------------------------- BCF ----

def test_bcf_file_size_and_signature(bcf_path):
    assert bcf_path.stat().st_size == 617_320_728
    with bcf_path.open("rb") as source:
        assert source.read(8) == b"AAMVHFSS"


def test_bcf_returns_video_and_edx_datasets_only(bcf_result):
    """Only Video and EDX are returned: no explicit elemental-map dataset."""
    titles = sorted(s["metadata"]["General"]["title"] for s in bcf_result["signals"])
    assert titles == ["EDX", "Video"]


def test_video_dataset_shape_and_dtype(video_signal):
    assert video_signal["data"].shape == (240, 1800)
    assert video_signal["data"].dtype == "uint16"


def test_edx_dataset_is_lazy_with_expected_shape_and_dtype(edx_signal):
    data = edx_signal["data"]
    assert isinstance(data, da.Array)
    assert data.shape == (240, 1800, 4096)
    assert data.dtype == "uint8"


def test_edx_lazy_array_is_a_single_chunk(edx_signal):
    """
    Observed RosettaSciIO 0.14.0 behavior (rsciio.bruker parse_hypermap): the lazy cube is
    ``da.from_delayed`` of one parse call, i.e. ONE chunk of 1,769,472,000 bytes. It is lazy
    until first access, but it is NOT chunked, so any slice would decode the whole cube.
    If a RosettaSciIO upgrade changes this, this test fails and FINDINGS.md must be updated.
    """
    data = edx_signal["data"]
    assert data.numblocks == (1, 1, 1)
    assert data.chunksize == (240, 1800, 4096)
    assert data.nbytes == 1_769_472_000
    assert len(data.dask.layers) == 2  # the delayed parse call and from_delayed; nothing computed


def test_spatial_calibration(video_signal, edx_signal):
    for signal in (video_signal, edx_signal):
        axes = {a["name"]: a for a in signal["axes"]}
        assert axes["height"]["size"] == 240 and axes["width"]["size"] == 1800
        assert axes["height"]["scale"] == pytest.approx(100.0148121593, rel=1e-9)
        assert axes["width"]["scale"] == pytest.approx(100.0148121593, rel=1e-9)
        assert axes["height"]["offset"] == 0 and axes["width"]["offset"] == 0
        assert axes["height"]["units"] == "µm" and axes["width"]["units"] == "µm"
    axes = {a["name"]: a for a in edx_signal["axes"]}
    assert axes["height"]["scale"] == axes["width"]["scale"]


def test_spectral_calibration(edx_signal):
    energy = next(a for a in edx_signal["axes"] if a["name"] == "Energy")
    assert energy["size"] == 4096
    assert energy["units"] == "keV"
    assert energy["scale"] == pytest.approx(0.010001, rel=1e-9)
    assert energy["offset"] == pytest.approx(-0.96079607, rel=1e-9)
    assert energy["offset"] + 4095 * energy["scale"] == pytest.approx(39.99, abs=0.01)


def test_edx_axis_order_matches_data_shape(edx_signal):
    assert [a["name"] for a in edx_signal["axes"]] == ["height", "width", "Energy"]
    assert [a["size"] for a in edx_signal["axes"]] == list(edx_signal["data"].shape)


def test_edx_instrument_and_acquisition_metadata(edx_signal):
    meta = edx_signal["metadata"]
    sem = meta["Acquisition_instrument"]["SEM"]
    assert sem["Detector"]["EDS"]["detector_type"] == "XFlash 430"
    assert sem["Detector"]["EDS"]["elevation_angle"] == 50.0
    assert sem["Detector"]["EDS"]["real_time"] == 12960.0
    # Units of beam_energy are unverified (FINDINGS.md): only the raw value is pinned.
    assert sem["beam_energy"] == 50
    assert meta["Sample"]["name"] == "Mapdaten"
    assert meta["General"]["date"] == "2026-07-30"
    assert meta["General"]["time"] == "09:55:52"
    assert meta["Signal"]["signal_type"] == "EDS_SEM"


def test_bcf_original_metadata_dimensions(video_signal):
    dsp = video_signal["original_metadata"]["DSP Configuration"]
    assert (dsp["ImageWidth"], dsp["ImageHeight"]) == (1800, 240)


# ---------------------------------------------------------------- RTX ----

def test_rtx_file_size_and_encoding_declaration(rtx_path):
    assert rtx_path.stat().st_size == 43_723_950
    with rtx_path.open("rb") as source:
        head = source.read(64)
    assert head.startswith(b'<?xml version="1.0" encoding="WI')
    assert b"1252" in head


def test_rtx_is_well_formed_with_expected_structure(rtx_report):
    assert rtx_report["well_formed"] is True, rtx_report["error"]
    assert rtx_report["root"] == "TRTProject"
    assert rtx_report["major_sections"] == ["RTHeader", "RTData"]
    assert dict(rtx_report["element_counts"]) == {
        "TRTProject": 1, "RTHeader": 1, "ProjectHeader": 1, "Date": 1, "Time": 1,
        "Creator": 1, "Comment": 1, "RTCompression": 1, "RTData": 1,
    }
    assert set(rtx_report["hierarchy"]) == {
        ("TRTProject", "RTHeader"), ("TRTProject", "RTData"),
        ("RTHeader", "ProjectHeader"), ("RTHeader", "RTCompression"),
        ("ProjectHeader", "Date"), ("ProjectHeader", "Time"),
        ("ProjectHeader", "Creator"), ("ProjectHeader", "Comment"),
    }


def test_rtx_declares_base64_and_zlib(rtx_report):
    compression = [a for t, _, a in rtx_report["scalar_values"] if t == "RTCompression"]
    assert compression == [{"compressor": "zlib", "encoder": "base64"}]


def test_rtx_payload_is_a_single_base64_like_text_element(rtx_report):
    assert dict(rtx_report["base64_like_elements"]) == {"RTData": 1}
    assert rtx_report["nul_bytes"] == 0
    assert rtx_report["bytes_scanned"] == 43_723_950


def test_rtx_header_date_and_time(rtx_report):
    values = {t: v for t, v, _ in rtx_report["scalar_values"] if t in ("Date", "Time")}
    assert values == {"Date": "30.7.2026", "Time": "15:58:47"}


def test_rtx_and_bcf_share_calendar_date(rtx_report, edx_signal):
    rtx_date = next(v for t, v, _ in rtx_report["scalar_values"] if t == "Date")
    assert normalize_date(rtx_date) == edx_signal["metadata"]["General"]["date"]


# ------------------------------------------------------ companion images ----

@pytest.mark.parametrize("name, size", [("CaFe.png", (1800, 240)), ("Video Mosaic.png", (948, 6000))])
def test_companion_png_dimensions(name, size):
    path = DATA_DIR / name
    if not path.is_file():
        pytest.skip(f"{name} not available")
    details = inspect_png(path)
    assert (details["width"], details["height"]) == size
