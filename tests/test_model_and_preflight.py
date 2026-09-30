"""
Authors: Andre L. Belem (https://github.com/andrebelem) and F.R.I.D.A.Y. (https://observatoriooceanografico.org/people/friday-bot)

Fast unit tests: the RTX scan/residual reader, the data-model rules, registration helpers, the dtype policy and
the resource preflight (including a hypothetical dataset far larger than RAM).
"""

from __future__ import annotations

import collections
from dataclasses import replace

import numpy as np
import pytest

from microxrf_to_netcdf import preflight, registration
from microxrf_to_netcdf.config import ConversionConfig
from microxrf_to_netcdf.errors import GridMismatchError, PreflightError, RTXFormatError
from microxrf_to_netcdf.model import build_layout, iso_from_rtx, split_label, video_sha256
from microxrf_to_netcdf.preflight import BCFScan, check_destination, choose_counts_dtype, clip, make_plan
from microxrf_to_netcdf.rtx import classify_images, scan_rtx, stream_rtx_planes

import synthetic_data as sd
import test_rtx_synthetic as rs


@pytest.fixture()
def pair(tmp_path):
    counts = sd.random_counts(10, 13, 64, seed=2)
    source = sd.synthetic_source(counts)
    rtx_path = tmp_path / "pair.rtx"
    extra = sd.write_pair_rtx(rtx_path, source.info.video)
    return source, rtx_path, extra


def test_scan_reports_geometry_and_keeps_all_non_pixel_content(pair):
    source, rtx_path, extra = pair
    scan = scan_rtx(rtx_path)
    assert [(i.name, i.width, i.height, i.itemsize, i.plane_count) for i in scan.images] == [
        ("Video Mosaic", 72, 50, 1, 3), ("Video Mosaic", 72, 50, 1, 3), ("Mapdaten", 13, 10, 2, 4)]
    assert scan.images[2].x_calibration == 100.0 and scan.images[0].x_calibration == 25.0
    assert [p.description for p in scan.images[2].planes] == ["Video 1", "Ca-KA", "Fe-KA", "Al-K"]
    assert scan.outer_attributes == {"compressor": "zlib", "encoder": "base64"} and scan.outer_header_fields["Date"] == "1.2.2026"
    assert "<Plane1><Description>Ca-KA</Description>" in scan.residual_xml
    assert "<Valid>0</Valid>" in scan.residual_xml and "LineCounter" in scan.residual_xml
    assert scan.residual_xml.count("<Data></Data>") == 4 + 3 + 3          # every plane's Data text is removed
    assert len(scan.residual_xml) < 20_000


def test_residual_xml_of_a_larger_plane_stays_small(tmp_path):
    video = np.zeros((200, 300), dtype=np.uint16)
    rtx_path = tmp_path / "big.rtx"
    sd.write_pair_rtx(rtx_path, video, mosaics="none")
    scan = scan_rtx(rtx_path)
    assert scan.payload_bytes > 100_000 and len(scan.residual_xml) < 5_000


def test_stream_delivers_every_plane_once_and_matches_the_scan_hashes(pair):
    _, rtx_path, _ = pair
    import hashlib
    scan = scan_rtx(rtx_path)
    seen = {}
    stream_rtx_planes(rtx_path, lambda image, plane, raw: seen.__setitem__((image, plane), hashlib.sha256(raw).hexdigest()))
    expected = {(i.index, p.index): p.sha256 for i in scan.images for p in i.planes}
    assert seen == expected


def test_layout_identifies_video_element_maps_and_deduplicates_mosaics(pair):
    source, rtx_path, _ = pair
    layout = build_layout(source.info, scan_rtx(rtx_path))
    assert layout.video_plane.description == "Video 1" and layout.video_plane.sha256 == video_sha256(source.info.video)
    assert [p.description for p in layout.element_planes] == ["Ca-KA", "Fe-KA", "Al-K"]
    assert [(m.pixels_variable, m.writes_pixels, m.duplicate_of) for m in layout.mosaics] == [
        ("pixels", True, None), ("pixels", False, "rtx:0")]
    assert all(c["status"] == "passed" for c in layout.checks)


def test_classify_refuses_an_image_it_would_otherwise_drop(tmp_path):
    video = np.arange(6 * 8, dtype=np.uint16).reshape(6, 8)
    rtx_path = tmp_path / "odd.rtx"
    odd = rs.image_xml("Odd", 4, 3, 2, planes=[np.arange(12)] * 2)
    grid = sd.image_xml("Mapdaten", 8, 6, 2, [("Video 1", video), ("Ca-KA", video)], "100,0", "1.1.2026", "1:00:00")
    rs.write_rtx(rtx_path, rs.payload_bytes(odd, grid))
    with pytest.raises(RTXFormatError, match="neither the map image nor a"):
        classify_images(scan_rtx(rtx_path), 6, 8)


def test_classify_requires_exactly_one_image_on_the_bcf_grid(tmp_path):
    rtx_path = tmp_path / "none.rtx"
    grid = sd.image_xml("Mapdaten", 8, 6, 2, [("Video 1", np.zeros((6, 8)))], "100,0", "1.1.2026", "1:00:00")
    rs.write_rtx(rtx_path, rs.payload_bytes(grid))
    with pytest.raises(GridMismatchError, match="exactly one"):
        classify_images(scan_rtx(rtx_path), 6, 9)


def test_unsupported_item_size_is_refused(tmp_path):
    rtx_path = tmp_path / "item4.rtx"
    body = rs.image_xml("Mapdaten", 2, 2, 2).replace("<ItemSize>2</ItemSize>", "<ItemSize>4</ItemSize>")
    rs.write_rtx(rtx_path, rs.payload_bytes(body))
    with pytest.raises(RTXFormatError):
        scan_rtx(rtx_path)


def test_label_split_and_timestamp_parsing():
    assert split_label("Ca-KA") == ("Ca", "KA") and split_label("Al-K") == ("Al", "K")
    assert split_label("Video 1") == ("Video 1", "") and split_label("Fe") == ("Fe", "")
    assert iso_from_rtx("30.7.2026", "9:38:26") == "2026-07-30T09:38:26"
    assert iso_from_rtx("garbage", "9:38:26") == "" and iso_from_rtx(None, None) == ""


def test_map_rectangle_ignores_the_degenerate_record_and_needs_exactly_one_real_one():
    degenerate = {"type": "TRTRectangleOverlayElement", "name": "Map", "rect": {"Left": 0, "Top": 0, "Right": 0, "Bottom": 0}}
    real = {"type": "TRTRectangleOverlayElement", "name": "Map", "rect": {"Left": 192, "Top": 89, "Right": 5941, "Bottom": 855}}
    assert registration.find_map_rectangle([degenerate, real]) == {"Left": 192, "Top": 89, "Right": 5941, "Bottom": 855}
    assert registration.find_map_rectangle([degenerate]) is None
    assert registration.find_map_rectangle([real, real]) is None
    assert registration.find_map_rectangle([{"type": "TRTTextOverlayElement", "name": "Map"}]) is None


def test_footprint_size_is_compared_with_the_acquisition_extent():
    rect = {"Left": 192, "Top": 89, "Right": 5941, "Bottom": 855}
    summary = registration.footprint_summary(rect, (31.308984675955, 31.3089846759549), (240, 1800),
                                             (100.0148121593, 100.014812159301))
    assert summary["rect_width_px_inclusive"] == 5750 and summary["rect_height_px_inclusive"] == 767
    assert summary["size_matches"] and abs(summary["rect_width_um"] - 180026.66) < 1
    wrong = registration.footprint_summary({"Left": 740, "Top": 89, "Right": 5181, "Bottom": 621},
                                           (31.31, 31.31), (240, 1800), (100.01, 100.01))
    assert not wrong["size_matches"]


def test_counts_dtype_policy():
    assert choose_counts_dtype(0, "auto")[0] == "uint8" and choose_counts_dtype(255, "auto")[0] == "uint8"
    assert choose_counts_dtype(256, "auto")[0] == "uint16" and choose_counts_dtype(65536, "auto")[0] == "uint32"
    assert choose_counts_dtype(10, "uint16")[0] == "uint16"
    with pytest.raises(PreflightError, match="truncate or wrap"):
        choose_counts_dtype(256, "uint8")
    with pytest.raises(PreflightError, match="exceeds uint32"):
        choose_counts_dtype(2**32, "auto")


def test_chunks_are_clipped_to_the_dimensions():
    assert clip((4, 60, 4096), (2, 1800, 100)) == (2, 60, 100)


def test_config_validation():
    for bad in (dict(band_lines=0), dict(decode_dtype="uint8"), dict(counts_dtype="int8"), dict(complevel=10),
                dict(counts_chunks=(0, 1, 1)), dict(memory_fraction=0)):
        with pytest.raises(ValueError):
            ConversionConfig(**bad).validate()
    ConversionConfig().validate()


Usage = collections.namedtuple("Usage", "total used free")


def test_a_dataset_far_larger_than_ram_is_planned_by_streaming_and_not_rejected_for_its_size(pair, monkeypatch, tmp_path):
    """37 TiB of logical counts on a machine with 8 GiB of free RAM: fine with one-line bands; the size is not the issue."""
    source, rtx_path, _ = pair
    scan = scan_rtx(rtx_path)
    layout = build_layout(source.info, scan)
    monkeypatch.setattr(preflight, "available_memory_bytes", lambda: 8 * 2**30)
    monkeypatch.setattr(preflight.shutil, "disk_usage", lambda path: Usage(2**60, 0, 2**60))
    huge = replace(source.info, height=100_000, width=100_000, channels=4096)
    stats = BCFScan(maximum=100, largest_line_bytes=200_000_000)
    config = ConversionConfig(band_lines=1, counts_chunks=(1, 60, 4096), decode_dtype="uint16")
    huge.header_bytes = 1_000_000
    plan = make_plan(huge, stats, layout, scan, config, tmp_path / "o" / "x.nc")
    assert plan.estimated_logical_total > 37 * 2**40 > 1000 * 8 * 2**30
    assert plan.memory_ok and plan.memory_estimated_peak_mib < 0.8 * 8 * 1024
    assert plan.estimated_logical_bytes["acquisition/counts"] == 100_000 * 100_000 * 4096
    with pytest.raises(PreflightError, match="reduce band_lines"):      # only the decoding band can be too large
        make_plan(huge, stats, layout, scan, ConversionConfig(band_lines=4, counts_chunks=(4, 60, 4096)),
                  tmp_path / "o" / "x.nc")


def test_disk_requirement_uses_the_uncompressed_worst_case_plus_a_margin(pair, monkeypatch, tmp_path):
    source, rtx_path, _ = pair
    scan = scan_rtx(rtx_path)
    layout = build_layout(source.info, scan)
    monkeypatch.setattr(preflight.shutil, "disk_usage", lambda path: Usage(10**12, 0, 10**12))
    plan = make_plan(source.info, BCFScan(maximum=5, largest_line_bytes=1000), layout, scan,
                     ConversionConfig(band_lines=4, counts_chunks=(4, 5, 64)), tmp_path / "o" / "x.nc")
    assert plan.disk_worst_case_bytes > plan.estimated_logical_total
    assert plan.disk_required_bytes >= plan.disk_worst_case_bytes + 512 * 2**20
    assert plan.as_dict()["status"] == "ESTIMATES (not measured)"


def test_destination_checks(tmp_path):
    inputs = [tmp_path / "a.rtx"]
    inputs[0].write_bytes(b"x")
    check_destination(tmp_path / "out" / "ok.nc", inputs, False)
    with pytest.raises(PreflightError, match="next to the original"):
        check_destination(tmp_path / "ok.nc", inputs, False)
    existing = tmp_path / "out" / "there.nc"
    existing.parent.mkdir()
    existing.write_bytes(b"x")
    with pytest.raises(PreflightError, match="already exists"):
        check_destination(existing, inputs, False)
    check_destination(existing, inputs, True)
    with pytest.raises(PreflightError, match="'data'"):
        check_destination(tmp_path / "Data" / "x.nc", inputs, False)
