"""
Authors: Andre L. Belem (https://github.com/andrebelem) and F.R.I.D.A.Y. (https://observatoriooceanografico.org/people/friday-bot)

Regression tests for the RTX values recorded in FINDINGS.md section 8. They read the original RTX
(and, for the cross-checks, a bounded part of the BCF and the PNG files) read-only and are skipped
when the files are absent. The EDS cube is never computed: only the first BCF scan lines are streamed.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from PIL import Image

from conftest import DATA_DIR, RTX_PATH
from microxrf_to_netcdf import rtx_payload as rtx

pytestmark = pytest.mark.large

TOOL = Path(__file__).resolve().parent.parent / "tools" / "diagnostics" / "inspect_rtx.py"

ELEMENTS = ["Ca-KA", "K-KA", "S-KA", "Si-K", "Ti-KA", "Cr-KA", "Mn-KA", "Fe-KA", "Ni-KA", "Zn-KA",
            "Sr-KA", "Pd-KA", "Rh-KA", "Al-K"]


@pytest.fixture(scope="module")
def report(rtx_decoded):
    return rtx_decoded["report"]


@pytest.fixture(scope="module")
def images(report):
    return dict(enumerate(report["payload"]["images"]))


def _pearson(a, b) -> float:
    return float(np.corrcoef(np.asarray(a, float).ravel(), np.asarray(b, float).ravel())[0, 1])


def _map_rectangles(image) -> list[dict]:
    """Non-degenerate rectangle overlays named 'Map' (the all-zero one is an unset placeholder)."""
    return [o["rect"] for o in image["overlays"] if o["name"] == "Map" and o["rect"]["Right"]]


# ---- decoding ---------------------------------------------------------------------------------

def test_decoding_succeeds_with_recorded_sizes(report):
    assert report["file_bytes"] == 43_723_950
    assert report["outer"]["attributes"] == {"compressor": "zlib", "encoder": "base64"}
    decode = report["decode"]
    assert decode["base64_chars"] == 43_723_624
    assert decode["decoded_bytes"] == 62_915_687
    assert decode["zlib_finished"] is True and decode["zlib_unused_bytes"] == 0
    assert decode["outer_tail"].strip() == "</TRTProject>"
    assert decode["sha256_decoded"] == "832f6f4c090e5413033621335127f7095ceb9ff7f9f7c3604bca2ffba00d4b61"


def test_decompressed_payload_is_xml_text_not_binary(report):
    payload = report["payload"]
    assert payload["signature_text"].startswith('<CompData><ClassInstance Type="TRTProject" Name="Bruker project"')
    assert report["decode"]["first_bytes_hex"].startswith("3c 43 6f 6d 70 44 61 74 61 3e")
    assert (payload["root"], payload["elements"], payload["max_depth"]) == ("CompData", 5305, 14)
    assert payload["class_instance_types"] == {
        "TRTProject": 1, "TRTImageData": 3, "TRTImageOverlay": 7, "TRTImageConfigurationData": 3,
        "TRTRectangleOverlayElement": 6, "TRTLineOverlayElement": 6, "TRTTextOverlayElement": 7,
        "TRTMapConfigurationData": 1, "TRTImagePalette": 1, "TRTMapViewSettings": 1}


def test_input_file_unchanged_by_decoding(rtx_path):
    before = rtx_path.stat()
    rtx.inspect_rtx_payload(rtx_path)
    after = rtx_path.stat()
    assert (before.st_size, before.st_mtime_ns) == (after.st_size, after.st_mtime_ns)


def test_peak_memory_of_decoding_is_far_below_payload_size(rtx_path):
    """Decode in a fresh process. Payload: 62.9 MB decoded, 43.7 MB encoded."""
    done = subprocess.run([sys.executable, str(TOOL), str(rtx_path)], capture_output=True, text=True,
                          encoding="utf-8", check=True)
    line = next(text for text in done.stdout.splitlines() if text.startswith("- time:"))
    peak = line.rsplit("peak working set:", 1)[1].strip().removesuffix("MiB").strip()
    if peak == "None":
        pytest.skip("peak working set is only measured on Windows")
    assert float(peak) < 150.0   # measured about 56 MiB including Python, NumPy and expat


# ---- structure ---------------------------------------------------------------------------------

def test_three_images_with_recorded_geometry(images):
    assert [i["name"] for i in images.values()] == ["Video Mosaic", "Video Mosaic", "Mapdaten"]
    geometry = [tuple(i["scalars"][k] for k in ("Width", "Height", "ItemSize", "PlaneCount", "MultiImage"))
                for i in images.values()]
    assert geometry == [("6000", "948", "1", "3", "0"), ("6000", "948", "1", "3", "0"), ("1800", "240", "2", "15", "1")]
    assert all(i["geometry_consistent"] for i in images.values())
    assert [i["scalars"]["Time"] for i in images.values()] == ["9:38:26", "9:55:51", "9:55:53"]
    assert {i["scalars"]["Date"] for i in images.values()} == {"30.7.2026"}


def test_mapdaten_plane_descriptions(images):
    descriptions = [images[2]["planes"][k]["Description"] for k in range(15)]
    assert descriptions == ["Video 1"] + ELEMENTS
    assert all("Description" not in images[i]["planes"][k] for i in (0, 1) for k in range(3))


def test_recorded_plane_statistics(images):
    planes = images[2]["planes"]
    assert (planes[0]["min"], planes[0]["max"], planes[0]["sum"]) == (3263, 65280, 11_806_561_069)
    assert (planes[1]["max"], planes[1]["sum"]) == (1607, 233_524_498)     # Ca-KA
    assert (planes[8]["max"], planes[8]["sum"]) == (2159, 44_031_298)      # Fe-KA
    assert all(planes[k]["dtype"] == "<u2" for k in range(15))
    assert images[0]["planes"][0]["sum"] == images[1]["planes"][0]["sum"] == 560_840_716


def test_first_two_mosaics_have_identical_planes_but_different_annotations(images):
    for k in range(3):
        assert images[0]["planes"][k]["sha256"] == images[1]["planes"][k]["sha256"]
    assert _map_rectangles(images[0]) == [{"Left": 740, "Top": 89, "Right": 5181, "Bottom": 621}]
    assert _map_rectangles(images[1]) == [{"Left": 192, "Top": 89, "Right": 5941, "Bottom": 855}]


def test_annotation_texts_and_display_settings(images):
    assert [o["text"] for o in images[2]["overlays"] if "text" in o] == ["40000 \xb5m", "Mapdaten", "HV: 50,0 kV", "Ca", "Fe"]
    assert [o["text"] for o in images[1]["overlays"] if "text" in o] == ["40000 \xb5m", "Mosaik"]
    used = {k: v["MapColor"] for k, v in images[2]["map_display"].items() if v.get("MapUsed") == "1"}
    assert used == {"MapImage1": "16711680", "MapImage8": "255"}


def test_no_spectral_detector_stage_or_quantification_elements(report):
    """Negative finding for this file: the payload vocabulary has no such elements."""
    vocabulary = " ".join(report["payload"]["tag_vocabulary"]).lower()
    for word in ("spectr", "detector", "energy", "roi", "region", "stage", "quant", "peak", "sample", "instrument",
                 "software", "version", "creator", "filename", "path", "guid", "uuid"):
        assert word not in vocabulary, word
    assert len(report["payload"]["tag_vocabulary"]) == 86


# ---- relationship to the BCF (internal evidence) ------------------------------------------------

def test_rtx_calibration_equals_bcf_axis_scale(images, video_signal):
    axes = {a["name"]: a for a in video_signal["axes"]}
    scalars = images[2]["scalars"]
    assert rtx.parse_decimal(scalars["XCalibration"]) == pytest.approx(axes["width"]["scale"], rel=1e-12)
    assert rtx.parse_decimal(scalars["YCalibration"]) == pytest.approx(axes["height"]["scale"], rel=1e-12)
    assert (int(scalars["Width"]), int(scalars["Height"])) == (axes["width"]["size"], axes["height"]["size"])


def test_rtx_video_plane_is_bit_identical_to_bcf_video(rtx_decoded, video_signal):
    bcf_video = np.asarray(video_signal["data"])       # 240 x 1800 uint16, 864,000 bytes
    plane = rtx_decoded["planes"][(2, 0)]
    assert bcf_video.dtype == plane.dtype == np.dtype("<u2") and bcf_video.shape == plane.shape == (240, 1800)
    assert np.array_equal(bcf_video, plane)


def test_rtx_mosaic_map_rectangle_matches_bcf_footprint_within_one_mosaic_pixel(images, video_signal):
    """Derived: rectangle size in mosaic pixels x mosaic calibration vs BCF size x BCF scale (both um)."""
    mosaic = images[1]
    scale = rtx.parse_decimal(mosaic["scalars"]["XCalibration"])
    (rect,) = _map_rectangles(mosaic)
    axes = {a["name"]: a for a in video_signal["axes"]}
    for extent, axis in ((rect["Right"] - rect["Left"], "width"), (rect["Bottom"] - rect["Top"], "height")):
        bcf_um = axes[axis]["size"] * axes[axis]["scale"]
        assert abs((extent + 1) * scale - bcf_um) < scale, axis    # inclusive rectangle, tolerance one mosaic pixel


def test_scale_bar_length_gives_mosaic_units_of_micrometres(images):
    """Derived: the 40000 um scale bar spans about 1278 mosaic pixels at 31.309 per pixel (units um)."""
    scale = rtx.parse_decimal(images[1]["scalars"]["XCalibration"])
    assert 1278 * scale == pytest.approx(40000, rel=0.002)


def test_rtx_element_planes_track_bcf_energy_windows(rtx_decoded, edx_signal):
    """First 40 scan lines only: streamed one line per band through the unmodified decoder (bounded memory).

    Windows are +-0.15 keV around tabulated K-alpha energies converted with the BCF energy axis. The RTX
    maps are NOT equal to these raw window sums (they are a different, processed product); Pearson r only.
    """
    import probe_bcf_streaming as probe

    energy = next(a for a in edx_signal["axes"] if a["name"] == "Energy")
    lines = 40
    kev = {"Ca-KA": 3.692, "Mn-KA": 5.899, "Fe-KA": 6.404, "Zn-KA": 8.638}
    channel = {n: int(round((e - energy["offset"]) / energy["scale"])) for n, e in kev.items()}
    window = {n: np.zeros((lines, 1800), np.int64) for n in kev}
    _, item = probe.open_stream(DATA_DIR / "GRF17A_9-29cm_slab3_Elemental_map.bcf")
    for first, band in probe.stream_bands(item, 1, 4096, np.uint16):
        for n, c in channel.items():
            window[n][first] = band[0, :, c - 15:c + 16].sum(axis=1, dtype=np.int64)
        if first == lines - 1:
            break
    for n in kev:
        plane = rtx_decoded["planes"][(2, ELEMENTS.index(n) + 1)][:lines]
        assert _pearson(plane, window[n]) > 0.99, n
    # An unrelated plane must not correlate as strongly (sanity check of the test itself).
    assert _pearson(rtx_decoded["planes"][(2, ELEMENTS.index("Al-K") + 1)][:lines], window["Fe-KA"]) < 0.5


# ---- relationship to the PNG files -------------------------------------------------------------

def test_cafe_png_is_the_ca_fe_composite_of_the_rtx_planes(rtx_decoded):
    picture = np.asarray(Image.open(DATA_DIR / "CaFe.png").convert("RGB"))
    assert picture.shape == (240, 1800, 3)
    ca, fe = rtx_decoded["planes"][(2, 1)], rtx_decoded["planes"][(2, 8)]
    assert _pearson(picture[..., 2], ca) > 0.97      # blue <-> Ca (MapColor 16711680 = 0xFF0000, BGR order)
    assert _pearson(picture[..., 0], fe) > 0.90      # red  <-> Fe (MapColor 255)
    assert abs(_pearson(picture[..., 1], ca)) < 0.1 and abs(_pearson(picture[..., 1], fe)) < 0.1   # green unused
    assert _pearson(picture[..., 2], fe) < 0.2 and _pearson(picture[..., 0], ca) < 0.2


def test_video_mosaic_png_is_the_rtx_mosaic_rotated_clockwise(rtx_decoded):
    picture = np.asarray(Image.open(DATA_DIR / "Video Mosaic.png").convert("RGB"))
    assert picture.shape == (6000, 948, 3)
    for channel in range(3):
        plane = rtx_decoded["planes"][(1, channel)]
        assert _pearson(picture[..., channel], np.rot90(plane, -1)) > 0.94   # about 0.95 to 0.97 measured; not pixel-exact
        assert _pearson(picture[..., channel], np.rot90(plane, 1)) < 0.7
        assert _pearson(picture[..., channel], plane.T) < 0.7
