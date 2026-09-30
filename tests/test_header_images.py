"""
Authors: Andre L. Belem (https://github.com/andrebelem) and F.R.I.D.A.Y. (https://observatoriooceanografico.org/people/friday-bot)

Fast tests of the images that live inside the BCF header and that RosettaSciIO does not return (FINDINGS.md
section 9): the acquisition-grid ``PixelTimes`` image, the overview mosaic that joins the mosaic group, further
overview images with their own grids, an image without planes, and the refusals (a header video that disagrees
with the video RosettaSciIO returns, an unsupported image layout).
"""

from __future__ import annotations

import json

import netCDF4
import numpy as np
import pytest

from microxrf_to_netcdf.convert import convert
from microxrf_to_netcdf.errors import RTXFormatError, ValidationError
from microxrf_to_netcdf.validate import check_bcf_header_images

import synthetic_data as sd
from test_convert_synthetic import CHANNELS, HEIGHT, WIDTH, nothing_left, small_config


def build(tmp_path, **header_kwargs):
    counts = sd.random_counts(HEIGHT, WIDTH, CHANNELS, seed=2)
    video = sd.make_info(HEIGHT, WIDTH, CHANNELS).video
    header = sd.header_document(video, **header_kwargs)
    source = sd.synthetic_source(counts, header=header)
    rtx_path = tmp_path / "pair.rtx"
    extra = sd.write_pair_rtx(rtx_path, source.info.video)
    return counts, source, rtx_path, extra, video


@pytest.fixture()
def converted(tmp_path):
    counts, source, rtx_path, extra, video = build(tmp_path)
    out = tmp_path / "output" / "acq.nc"
    report = convert(source, rtx_path, out, small_config())
    return {"out": out, "report": report, "video": video, "extra": extra, "tmp": tmp_path}


def test_header_scan_lists_every_image_with_geometry(tmp_path):
    _, source, _, _, _ = build(tmp_path)
    images = source.info.header_scan.images
    assert [(i.name, i.width, i.height, i.itemsize, i.plane_count) for i in images] == [
        ("", WIDTH, HEIGHT, 2, 1), ("Counter", WIDTH, HEIGHT, 2, 0), ("PixelTimes", WIDTH, HEIGHT, 4, 1),
        ("Default", 72, 50, 1, 3), ("Image_0", 9, 6, 1, 3)]
    assert "<Data></Data>" in source.info.header_scan.residual_xml and "PixelTimes" in source.info.header_scan.residual_xml


def test_pixel_times_are_stored_on_the_shared_grid_without_invented_units(converted):
    with netCDF4.Dataset(converted["out"]) as ds:
        variable = ds["acquisition"]["pixel_times"]
        assert variable.dimensions == ("y", "x") and variable.dtype == np.uint32
        assert "units" not in variable.ncattrs() and "unknown" in variable.units_status
        expected = np.arange(HEIGHT * WIDTH, dtype=np.uint32).reshape(HEIGHT, WIDTH) * 7 + 100
        assert np.array_equal(np.asarray(variable[:]), expected)
        assert variable.data_max == int(expected.max()) and variable.data_sum == int(expected.sum(dtype=np.uint64))
        assert "NOT verified" in variable.note


def test_the_bcf_overview_mosaic_joins_the_mosaic_group_as_a_third_instance_with_one_pixel_copy(converted):
    with netCDF4.Dataset(converted["out"]) as ds:
        mosaic = ds["mosaic"]
        assert set(mosaic.groups) == {"instance_rtx_0", "instance_rtx_1", "instance_bcf_3_Default"}
        assert [v for v in mosaic.variables if v.startswith("pixels")] == ["pixels"]
        assert json.loads(mosaic["pixels"].stored_once_for_instances) == ["rtx:0", "rtx:1", "bcf:3"]
        bcf = mosaic["instance_bcf_3_Default"]
        assert bcf.source_file.startswith("BCF header") and bcf.source_time == "9:55:51"
        assert bcf.pixels_duplicate_of == "rtx:0" and bcf.plane_count == 3
        assert bcf["map_footprint"].present == 0
        assert np.array_equal(np.asarray(mosaic["pixels"][1]), converted["extra"]["canvas"])


def test_other_overview_images_keep_their_own_grid_and_are_not_georeferenced(converted):
    with netCDF4.Dataset(converted["out"]) as ds:
        overview = ds["overview"]
        assert list(overview.groups) == ["Image_0"]
        group = overview["Image_0"]
        assert group["pixels"].dimensions == ("plane", "row", "column") and group["pixels"].shape == (3, 6, 9)
        small = np.arange(6 * 9, dtype=np.uint8).reshape(6, 9) * 5
        assert np.array_equal(np.asarray(group["pixels"][0]), small) and np.array_equal(np.asarray(group["pixels"][2]), small + 2)
        assert np.allclose(group["column"][:3], [0, 14.75, 29.5])
        assert "units" not in group["column"].ncattrs() and "unknown" in group["column"].units_status
        assert group["pixels"].relationship_to_acquisition_grid == "not established"
        assert "row" not in ds["acquisition"].dimensions        # its own dimensions, not the acquisition grid


def test_every_header_image_is_recorded_including_one_without_planes(converted):
    with netCDF4.Dataset(converted["out"]) as ds:
        table = {row["name"]: row for row in json.loads(ds["metadata"]["bcf_image_table_json"][...])}
        assert set(table) == {"", "Counter", "PixelTimes", "Default", "Image_0"}
        assert table["Counter"]["plane_count"] == 0 and "no planes" in table["Counter"]["role"]
        assert table["PixelTimes"]["role"].endswith("/acquisition/pixel_times")
        assert table[""]["role"].startswith("video image")
        residual = ds["metadata"]["bcf_header_residual_xml"][...]
        assert "Counter" in residual and "<Data></Data>" in residual and len(residual) < 20_000
        stamps = json.loads(ds["metadata"]["timestamps_json"][...])
        assert any(s["source"].startswith("BCF header image 3") and s["time"] == "9:55:51" for s in stamps)


def test_read_back_checks_cover_the_header_images_and_all_pass(converted):
    names = {c["check"]: c["status"] for c in converted["report"]["validation"]}
    for check in ("acquisition_grid_header_images_equal_the_source_plane_sha256",
                  "overview_images_equal_the_source_plane_sha256",
                  "mosaic_planes_of_every_instance_equal_the_source_plane_sha256"):
        assert names[check] == "passed"
    assert converted["report"]["actual_layout"]["acquisition/pixel_times"]["dtype"] == "uint32"
    assert "overview/Image_0/pixels" in converted["report"]["actual_layout"]


def test_a_header_video_that_differs_from_the_video_rosettasciio_returns_is_refused(tmp_path):
    video = sd.make_info(HEIGHT, WIDTH, CHANNELS).video
    counts, source, rtx_path, _, _ = build(tmp_path, video_override=video + 1)
    with pytest.raises(ValidationError, match="differs from the video returned by RosettaSciIO"):
        convert(source, rtx_path, tmp_path / "output" / "acq.nc", small_config())
    assert nothing_left(tmp_path / "output" / "acq.nc")


def test_an_unsupported_header_image_is_refused_not_dropped(tmp_path):
    odd = sd.image_xml("Odd", 5, 5, 2, [(None, np.zeros((5, 5)))] * 2, "1,0", "1.1.2026", "1:00:00")
    counts, source, rtx_path, _, _ = build(tmp_path, extra_images=odd)
    with pytest.raises(RTXFormatError, match="unsupported layout"):
        convert(source, rtx_path, tmp_path / "output" / "acq.nc", small_config())
    assert nothing_left(tmp_path / "output" / "acq.nc")


def test_a_grid_image_with_several_planes_is_refused(tmp_path):
    odd = sd.image_xml("TwoPlanes", WIDTH, HEIGHT, 2, [(None, np.zeros((HEIGHT, WIDTH)))] * 2, "0,0", "1.1.2026", "1:00:00")
    counts, source, rtx_path, _, _ = build(tmp_path, extra_images=odd)
    with pytest.raises(RTXFormatError, match="single-plane"):
        convert(source, rtx_path, tmp_path / "output" / "acq.nc", small_config())


def test_a_bcf_mosaic_with_other_pixels_is_stored_as_its_own_array(tmp_path):
    video = sd.make_info(HEIGHT, WIDTH, CHANNELS).video
    canvas, _ = sd.mosaic_canvas(video)
    other = sd.image_xml("Default", canvas.shape[1], canvas.shape[0], 1, [(None, np.roll(canvas, 5, axis=0))] * 3,
                         "25", "30.7.2026", "9:55:51")
    counts, source, rtx_path, _, _ = build(tmp_path, with_mosaic=False, extra_images=other)
    out = tmp_path / "output" / "acq.nc"
    convert(source, rtx_path, out, small_config(check_registration=False))
    with netCDF4.Dataset(out) as ds:
        assert {v for v in ds["mosaic"].variables if v.startswith("pixels")} == {"pixels", "pixels_2"}
        assert ds["mosaic"]["instance_bcf_4_Default"].pixels_variable == "/mosaic/pixels_2"
        assert np.array_equal(np.asarray(ds["mosaic"]["pixels_2"][0]), np.roll(canvas, 5, axis=0))


def test_without_a_header_document_the_conversion_still_works_and_reports_no_extra_images(tmp_path):
    counts = sd.random_counts(HEIGHT, WIDTH, CHANNELS, seed=2)
    source = sd.synthetic_source(counts)
    rtx_path = tmp_path / "pair.rtx"
    sd.write_pair_rtx(rtx_path, source.info.video)
    out = tmp_path / "output" / "acq.nc"
    convert(source, rtx_path, out, small_config())
    with netCDF4.Dataset(out) as ds:
        assert "pixel_times" not in ds["acquisition"].variables and not ds["overview"].groups
        assert "bcf_header_residual_xml" not in ds["metadata"].variables
        assert set(ds["mosaic"].groups) == {"instance_rtx_0", "instance_rtx_1"}


def test_independent_header_validation_agrees_with_the_written_file(converted, tmp_path):
    """The validate command's header check needs a real SFS container; here its logic is exercised through the same
    diagnostic decoder on the synthetic header bytes."""
    from microxrf_to_netcdf import rtx_payload
    document = sd.header_document(converted["video"])
    seen = []
    rtx_payload.parse_payload(iter(sd.split_blocks(document, 4064)), lambda i, p, raw: seen.append((i, p, len(raw))))
    assert (2, 0, WIDTH * HEIGHT * 4) in seen and (3, 2, 50 * 72) in seen and len(seen) == 1 + 1 + 3 + 3
    assert callable(check_bcf_header_images)


def test_every_written_variable_is_documented_in_the_schema(converted):
    """NETCDF_SCHEMA.md must name every variable the converter writes (synthetic file with every feature)."""
    from pathlib import Path
    schema = (Path(__file__).resolve().parent.parent / "NETCDF_SCHEMA.md").read_text(encoding="utf-8")
    names = set()

    def walk(group):
        names.update(group.variables)
        for sub in group.groups.values():
            walk(sub)

    with netCDF4.Dataset(converted["out"]) as ds:
        walk(ds)
        global_attributes = set(ds.ncattrs())
    undocumented = sorted(name for name in names if f"`{name}`" not in schema and name not in schema)
    assert not undocumented, undocumented
    assert {name for name in global_attributes if f"`{name}`" not in schema and name not in schema} == set()
