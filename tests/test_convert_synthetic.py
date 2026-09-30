"""
Authors: Andre L. Belem (https://github.com/andrebelem) and F.R.I.D.A.Y. (https://observatoriooceanografico.org/people/friday-bot)

Fast end-to-end tests of the converter on synthetic BCF/RTX pairs: sequential processing, several band sizes,
overflow protection, valid zeros, lazy reopening, chunk boundaries, grid and video mismatches, truncated and
malformed inputs, insufficient disk, interruption and cleanup, destination safety and provenance.
"""

from __future__ import annotations

import collections
import hashlib
import json
import shutil
from pathlib import Path

import netCDF4
import numpy as np
import pytest
import xarray as xr

from microxrf_to_netcdf import preflight
from microxrf_to_netcdf.config import ConversionConfig
from microxrf_to_netcdf.convert import convert
from microxrf_to_netcdf.errors import (BCFStreamError, ConversionError, GridMismatchError, PreflightError, RTXFormatError,
                               ValidationError)
from microxrf_to_netcdf.preflight import partial_path
from microxrf_to_netcdf.validate import check_lazy_structure

import synthetic_data as sd

HEIGHT, WIDTH, CHANNELS = 10, 13, 64


def small_config(**overrides) -> ConversionConfig:
    values = dict(band_lines=4, counts_chunks=(4, 5, 64), map_chunks=(1, 5, 7), video_chunks=(5, 7),
                  mosaic_chunks=(1, 16, 20), complevel=2)
    values.update(overrides)
    return ConversionConfig(**values)


def build(tmp_path: Path, big=None, high=6, mosaics="consistent", **rtx_kwargs):
    counts = sd.random_counts(HEIGHT, WIDTH, CHANNELS, seed=2, high=high, big=big)
    source = sd.synthetic_source(counts)
    rtx_path = tmp_path / "pair.rtx"
    extra = sd.write_pair_rtx(rtx_path, source.info.video, mosaics=mosaics, **rtx_kwargs)
    return counts, source, rtx_path, extra


@pytest.fixture()
def out(tmp_path):
    return tmp_path / "output" / "acq.nc"


@pytest.fixture()
def converted(tmp_path, out):
    counts, source, rtx_path, extra = build(tmp_path)
    report = convert(source, rtx_path, out, small_config())
    return {"counts": counts, "source": source, "rtx": rtx_path, "extra": extra, "report": report, "out": out}


def nothing_left(out: Path) -> bool:
    return not out.exists() and not partial_path(out).exists() and not (out.parent.exists() and any(out.parent.iterdir()))


# ---------------------------------------------------------------------------------------------------
# Content and lazy reopening
# ---------------------------------------------------------------------------------------------------


def test_conversion_produces_a_complete_validated_file(converted):
    out, report = converted["out"], converted["report"]
    assert out.is_file() and not partial_path(out).exists()
    assert report["output_size_bytes"] == out.stat().st_size
    assert hashlib.sha256(out.read_bytes()).hexdigest() == report["output_sha256"]
    assert (out.parent / (out.name + ".sha256")).read_text().split()[0] == report["output_sha256"]
    assert (out.parent / (out.name + ".report.json")).is_file()
    assert all(c["status"] == "passed" for c in report["validation"])
    with netCDF4.Dataset(out) as ds:
        assert ds.conversion_status == "complete" and ds.microxrf_to_netcdf_schema_version == "1.0.0"
        assert ds.source_bcf_filename == "synthetic.bcf" and ds.source_rtx_filename == "pair.rtx"


def test_lazy_reopen_with_xarray_and_dask_preserves_integer_counts_and_zeros(converted):
    counts = converted["counts"]
    with xr.open_dataset(converted["out"], group="acquisition", chunks={}) as acquisition:
        array = acquisition["counts"]
        assert hasattr(array.data, "dask") and array.dims == ("y", "x", "energy")
        assert array.dtype == np.uint8 and array.data.chunksize == (4, 5, 64)
        assert "_FillValue" not in array.encoding and "_FillValue" not in array.attrs
        values = array.values
    assert values.dtype == np.uint8 and np.array_equal(values, counts)
    assert int((values == 0).sum()) == int((counts == 0).sum()) > 0  # valid zeros stay integer zeros, never NaN


def test_chunk_boundaries_and_partial_reads(converted):
    counts = converted["counts"]
    with xr.open_dataset(converted["out"], group="acquisition", chunks={}) as acquisition:
        array = acquisition["counts"]
        window = array.isel(y=slice(2, 9), x=slice(3, 12), energy=slice(10, 50)).values  # crosses y and x chunks
        assert np.array_equal(window, counts[2:9, 3:12, 10:50])
        for y, x in [(0, 0), (3, 4), (4, 5), (9, 12), (7, 9)]:  # chunk corners and the ragged edge
            assert np.array_equal(array.isel(y=y, x=x).values, counts[y, x])
        assert np.array_equal(array.isel(energy=17).values, counts[:, :, 17])


@pytest.mark.parametrize("band_lines", [4, 8, 12])
def test_band_size_does_not_change_the_output(tmp_path, band_lines):
    counts, source, rtx_path, _ = build(tmp_path)
    out = tmp_path / "o" / f"b{band_lines}.nc"
    convert(source, rtx_path, out, small_config(band_lines=band_lines))
    with netCDF4.Dataset(out) as ds:
        stored = np.asarray(ds["acquisition"]["counts"][:])
        assert ds["acquisition"]["counts"].chunking() == [4, 5, 64]
    assert np.array_equal(stored, counts)


def test_band_lines_must_align_with_the_chunk_extent(tmp_path, out):
    _, source, rtx_path, _ = build(tmp_path)
    with pytest.raises(PreflightError, match="multiple"):
        convert(source, rtx_path, out, small_config(band_lines=3))
    assert nothing_left(out)


def test_coordinates_calibration_and_provenance(converted):
    source = converted["source"]
    with netCDF4.Dataset(converted["out"]) as ds:
        acq = ds["acquisition"]
        assert np.array_equal(acq["energy"][:], source.info.calib_abs + source.info.calib_lin * np.arange(CHANNELS))
        assert acq["energy"].units == "keV" and acq["y"].units == "um" and acq["x"].units == "um"
        assert np.array_equal(acq["y"][:], np.arange(HEIGHT) * 100.0)
        calibration = acq["energy_calibration"]
        assert calibration.calib_abs_raw == source.info.calib_abs and calibration.calib_lin_raw == source.info.calib_lin
        assert np.array_equal(acq["video"][:], source.info.video)
        provenance = json.loads(ds["metadata"]["provenance_json"][...])
        assert provenance["parameters"]["band_lines"] == 4 and provenance["software"]["microxrf_to_netcdf"]
        assert provenance["bcf_scan"]["total_counts"] == int(converted["counts"].sum())
        assert provenance["preflight_estimates"]["status"].startswith("ESTIMATES")
        assert provenance["actual_layout"]["acquisition/counts"]["chunks"] == [4, 5, 64]
        assert json.loads(ds["metadata"]["unresolved_properties_json"][...])
        assert ds["metadata"]["bcf_original_metadata_json"][...].startswith("{")
        assert "TRTProject" in ds["metadata"]["rtx_payload_residual_xml"][...]


def test_element_maps_keep_identifiers_and_carry_no_invented_units(converted):
    maps = converted["extra"]["maps"]
    with xr.open_dataset(converted["out"], group="acquisition", chunks={}) as acquisition:
        element_maps = acquisition["element_maps"]
        assert list(acquisition["element"].values) == ["Ca", "Fe", "Al"]
        assert list(acquisition["element_line"].values) == ["KA", "KA", "K"]
        assert list(acquisition["element_label"].values) == ["Ca-KA", "Fe-KA", "Al-K"]
        assert element_maps.dtype == np.uint16 and "units" not in element_maps.attrs
        assert element_maps.attrs["processing_method"] == "unknown" and "unknown" in element_maps.attrs["units_status"]
        for k, (_, plane) in enumerate(maps):
            assert np.array_equal(element_maps.isel(element=k).values, plane)
        assert acquisition["video"].dims == element_maps.dims[1:]  # one shared grid


def test_mosaic_pixels_are_stored_once_but_each_instance_keeps_its_own_metadata(converted):
    with netCDF4.Dataset(converted["out"]) as ds:
        mosaic = ds["mosaic"]
        assert [v for v in mosaic.variables if v.startswith("pixels")] == ["pixels"]
        assert np.array_equal(mosaic["pixels"][0], converted["extra"]["canvas"])
        assert set(mosaic.groups) == {"instance_rtx_0", "instance_rtx_1"}
        first, second = mosaic["instance_rtx_0"], mosaic["instance_rtx_1"]
        assert first.source_time == "9:38:26" and second.source_time == "9:55:51"
        assert first.timestamp_iso_naive == "2026-07-30T09:38:26"
        assert second.pixels_duplicate_of == "rtx:0"
        assert json.loads(first["annotations_json"][...]) != json.loads(second["annotations_json"][...])
        assert "assumed scientifically interchangeable" in mosaic["pixels"].equal_pixels_note
        assert "mosaic_x" in mosaic.dimensions and "x" not in mosaic.dimensions  # its own grid


def test_different_mosaic_pixels_are_stored_separately(tmp_path, out):
    _, source, rtx_path, _ = build(tmp_path, mosaics="different")
    convert(source, rtx_path, out, small_config(check_registration=False))
    with netCDF4.Dataset(out) as ds:
        assert {v for v in ds["mosaic"].variables if v.startswith("pixels")} == {"pixels", "pixels_2"}


def test_registration_is_verified_for_a_consistent_mosaic_and_recorded_as_unverified_otherwise(tmp_path, out):
    _, source, rtx_path, _ = build(tmp_path)
    report = convert(source, rtx_path, out, small_config())
    with netCDF4.Dataset(out) as ds:
        good = ds["mosaic"]["instance_rtx_1"]["map_footprint"]
        assert good.registration_status == "verified" and good.size_matches == 1
        assert "verified" in good.relationship_to_acquisition_grid
        assert json.loads(good.registration_json)["correlation_by_orientation"]["identity"] > 0.99
        other = ds["mosaic"]["instance_rtx_0"]["map_footprint"]  # its Map rectangle has another size
        assert other.registration_status == "not_verified"
    assert any("registration" in c["check"] for c in report["validation"])


def test_a_mirrored_mosaic_is_recorded_as_not_verified(tmp_path, out):
    _, source, rtx_path, _ = build(tmp_path, mosaics="flipped")
    convert(source, rtx_path, out, small_config())
    with netCDF4.Dataset(out) as ds:
        footprint = ds["mosaic"]["instance_rtx_1"]["map_footprint"]
        result = json.loads(footprint.registration_json)
        status, relationship = footprint.registration_status, footprint.relationship_to_acquisition_grid
    assert status == "not_verified"
    by_orientation = result["correlation_by_orientation"]
    assert by_orientation["flip_vertical"] > 0.99 > by_orientation["identity"]
    assert "not verified" in relationship


# ---------------------------------------------------------------------------------------------------
# dtype policy and overflow
# ---------------------------------------------------------------------------------------------------


def test_dtype_is_chosen_from_the_observed_maximum(tmp_path):
    for big, expected in (({(3, 4, 5): 255}, "uint8"), ({(3, 4, 5): 256}, "uint16"), ({(3, 4, 5): 70000}, "uint32")):
        case = tmp_path / expected
        case.mkdir()
        counts, source, rtx_path, _ = build(case, big=big)
        out = case / "o" / "a.nc"
        report = convert(source, rtx_path, out, small_config())
        assert report["counts_dtype"] == expected
        with netCDF4.Dataset(out) as ds:
            stored = ds["acquisition"]["counts"]
            assert stored.dtype == np.dtype(expected) and int(np.max(stored[3])) == list(big.values())[0]
            assert np.array_equal(np.asarray(stored[:]), counts)


def test_an_explicit_dtype_that_would_truncate_is_refused_before_writing(tmp_path, out):
    _, source, rtx_path, _ = build(tmp_path, big={(1, 1, 1): 300})
    with pytest.raises(PreflightError, match="cannot hold the observed maximum"):
        convert(source, rtx_path, out, small_config(counts_dtype="uint8"))
    assert nothing_left(out)


def test_a_count_that_wraps_in_uint16_is_still_stored_exactly(tmp_path, out):
    """70000 wraps to 4464 in a uint16 decode. The exact uint32 scan pass sees it, so the write pass decodes in
    uint32 and the output is exact; the report says so."""
    counts, source, rtx_path, _ = build(tmp_path, big={(3, 4, 5): 70000, (8, 2, 9): 40000})
    report = convert(source, rtx_path, out, small_config())
    assert report["counts_dtype"] == "uint32"
    with netCDF4.Dataset(out) as ds:
        stored = np.asarray(ds["acquisition"]["counts"][:])
        assert stored[3, 4, 5] == 70000 and stored[8, 2, 9] == 40000 and np.array_equal(stored, counts)
        decoding = json.loads(ds["metadata"]["provenance_json"][...])["decoding"]
        assert decoding["scan_pass"].startswith("uint32") and decoding["write_and_readback_passes"].startswith("uint32")


def test_uint16_write_pass_is_used_when_the_exact_maximum_fits_and_is_cross_checked(tmp_path, out):
    counts, source, rtx_path, _ = build(tmp_path, big={(3, 4, 5): 40000})
    report = convert(source, rtx_path, out, small_config())
    assert report["counts_dtype"] == "uint16"
    checks = {c["check"]: c for c in report["validation"]}
    assert checks["write_pass_decode_equals_exact_uint32_scan_decode"]["status"] == "passed"
    with netCDF4.Dataset(out) as ds:
        decoding = json.loads(ds["metadata"]["provenance_json"][...])["decoding"]
        assert decoding["write_and_readback_passes"].startswith("uint16") and decoding["widened_bands_in_write_pass"] >= 1


def test_uint32_decoding_gives_the_same_file(tmp_path, out):
    counts, source, rtx_path, _ = build(tmp_path, big={(1, 1, 1): 300})
    convert(source, rtx_path, out, small_config(decode_dtype="uint32"))
    with netCDF4.Dataset(out) as ds:
        assert np.array_equal(np.asarray(ds["acquisition"]["counts"][:]), counts)


# ---------------------------------------------------------------------------------------------------
# Grid, video and input faults: clear failure, no output
# ---------------------------------------------------------------------------------------------------


def test_mismatched_spatial_grid_is_refused(tmp_path, out):
    counts = sd.random_counts(HEIGHT, WIDTH, CHANNELS, seed=2)
    source = sd.synthetic_source(counts)
    rtx_path = tmp_path / "pair.rtx"
    sd.write_pair_rtx(rtx_path, source.info.video, grid_size=(HEIGHT, WIDTH + 1),
                      video_override=np.zeros((HEIGHT, WIDTH + 1), dtype=np.uint16))
    with pytest.raises(GridMismatchError, match="exactly one"):
        convert(source, rtx_path, out, small_config())
    assert nothing_left(out)


def test_mismatched_calibration_is_refused(tmp_path, out):
    _, source, _, _ = build(tmp_path)
    rtx_path = tmp_path / "other.rtx"
    sd.write_pair_rtx(rtx_path, source.info.video, calibration="50,0")
    with pytest.raises(GridMismatchError, match="calibration"):
        convert(source, rtx_path, out, small_config())
    assert nothing_left(out)


def test_video_plane_that_differs_from_the_bcf_video_is_a_validation_failure(tmp_path, out):
    _, source, _, _ = build(tmp_path)
    rtx_path = tmp_path / "other.rtx"
    sd.write_pair_rtx(rtx_path, source.info.video, video_override=source.info.video + 1)
    with pytest.raises(ValidationError, match="bit-identical|refusing to merge"):
        convert(source, rtx_path, out, small_config())
    assert nothing_left(out)


def test_truncated_rtx_fails_and_leaves_no_output(tmp_path, out):
    _, source, rtx_path, _ = build(tmp_path)
    data = rtx_path.read_bytes()
    rtx_path.write_bytes(data[: len(data) // 2])
    with pytest.raises(RTXFormatError):
        convert(source, rtx_path, out, small_config())
    assert nothing_left(out)


def test_truncated_bcf_stream_fails_before_anything_is_written(tmp_path, out):
    counts = sd.random_counts(HEIGHT, WIDTH, CHANNELS, seed=2)
    stream = sd.build_stream(counts)
    source = sd.SyntheticBCF(sd.synthetic_source(counts).info, lambda: iter(sd.split_blocks(stream[:-40], 64)))
    rtx_path = tmp_path / "pair.rtx"
    sd.write_pair_rtx(rtx_path, source.info.video)
    with pytest.raises(BCFStreamError, match="truncated"):
        convert(source, rtx_path, out, small_config())
    assert nothing_left(out)


def test_decoded_pixels_above_the_header_sum_spectrum_are_refused(tmp_path, out):
    counts = sd.random_counts(HEIGHT, WIDTH, CHANNELS, seed=2)
    source = sd.synthetic_source(counts, sum_spectrum=np.zeros(CHANNELS, dtype=np.uint64))
    rtx_path = tmp_path / "pair.rtx"
    sd.write_pair_rtx(rtx_path, source.info.video)
    with pytest.raises(ValidationError, match="exceed the BCF header sum spectrum"):
        convert(source, rtx_path, out, small_config())
    assert nothing_left(out)


def test_rtx_without_element_maps_is_refused(tmp_path, out):
    _, source, _, _ = build(tmp_path)
    rtx_path = tmp_path / "nomaps.rtx"
    sd.write_pair_rtx(rtx_path, source.info.video, element_names=())
    with pytest.raises(RTXFormatError, match="no element-map plane"):
        convert(source, rtx_path, out, small_config())
    assert nothing_left(out)


# ---------------------------------------------------------------------------------------------------
# Resources: disk and memory
# ---------------------------------------------------------------------------------------------------

Usage = collections.namedtuple("Usage", "total used free")


def test_insufficient_disk_space_fails_safely_before_creating_anything(tmp_path, out, monkeypatch):
    _, source, rtx_path, _ = build(tmp_path)
    monkeypatch.setattr(preflight.shutil, "disk_usage", lambda path: Usage(10**9, 10**9 - 1000, 1000))
    with pytest.raises(PreflightError, match="insufficient free disk space"):
        convert(source, rtx_path, out, small_config())
    assert nothing_left(out)


def test_low_disk_can_be_overridden_explicitly_and_is_recorded(tmp_path, out, monkeypatch):
    _, source, rtx_path, _ = build(tmp_path)
    monkeypatch.setattr(preflight.shutil, "disk_usage", lambda path: Usage(10**9, 10**9 - 1000, 1000))
    report = convert(source, rtx_path, out, small_config(allow_low_disk=True))
    assert out.is_file() and any("OVERRIDDEN" in w for w in report["plan"]["warnings"])


def test_unknown_free_disk_space_is_a_failure_not_a_guess(tmp_path, out, monkeypatch):
    _, source, rtx_path, _ = build(tmp_path)

    def broken(path):
        raise OSError("no such volume")

    monkeypatch.setattr(preflight.shutil, "disk_usage", broken)
    with pytest.raises(PreflightError, match="cannot determine free disk space"):
        convert(source, rtx_path, out, small_config())
    assert nothing_left(out)


def test_insufficient_memory_fails_before_writing_unless_overridden(tmp_path, out, monkeypatch):
    _, source, rtx_path, _ = build(tmp_path)
    monkeypatch.setattr(preflight, "available_memory_bytes", lambda: 50 * 2**20)
    with pytest.raises(PreflightError, match="estimated peak memory"):
        convert(source, rtx_path, out, small_config())
    assert nothing_left(out)
    report = convert(source, rtx_path, out, small_config(allow_low_memory=True))
    assert report["plan"]["memory_mib"]["sufficient"] is False


def test_logical_size_above_ram_is_not_a_reason_to_refuse(tmp_path, out, monkeypatch):
    """Only the estimated PEAK memory is compared with free RAM; the dataset size may exceed it freely."""
    _, source, rtx_path, _ = build(tmp_path)
    monkeypatch.setattr(preflight, "available_memory_bytes", lambda: 2 * 2**30)
    report = convert(source, rtx_path, out, small_config())
    assert report["estimated_logical_bytes"] > 0 and out.is_file()


def test_dry_run_reports_the_plan_and_writes_nothing(tmp_path, out):
    _, source, rtx_path, _ = build(tmp_path, big={(1, 1, 1): 300})
    plan = convert(source, rtx_path, out, small_config(), dry_run=True)
    assert plan["dry_run"] and plan["plan"]["counts_dtype"] == "uint16"
    sizes = plan["plan"]["estimated_logical_bytes"]
    assert sizes["acquisition/counts"] == HEIGHT * WIDTH * CHANNELS * 2
    assert plan["plan"]["chunks"]["counts"] == [4, 5, 64]
    assert plan["plan"]["disk"]["sufficient"] and plan["plan"]["memory_mib"]["estimated_peak"] > 0
    assert not out.parent.exists() or not any(out.parent.iterdir())


# ---------------------------------------------------------------------------------------------------
# Interruption, cleanup and destination safety
# ---------------------------------------------------------------------------------------------------


def test_failure_in_the_middle_leaves_no_partial_or_final_file(tmp_path, out):
    _, source, rtx_path, _ = build(tmp_path)

    def fault(stage, index):
        if stage == "band" and index == 1:
            raise RuntimeError("disk removed")

    with pytest.raises(ConversionError, match="partial file was removed"):
        convert(source, rtx_path, out, small_config(), _fault=fault)
    assert nothing_left(out)
    convert(source, rtx_path, out, small_config())  # a later run is unaffected
    assert out.is_file()


@pytest.mark.parametrize("stage", ["after_create", "after_write"])
def test_interruption_at_other_stages_is_cleaned_up(tmp_path, out, stage):
    _, source, rtx_path, _ = build(tmp_path)

    def fault(current, index):
        if current == stage:
            raise KeyboardInterrupt

    with pytest.raises(KeyboardInterrupt):
        convert(source, rtx_path, out, small_config(), _fault=fault)
    assert nothing_left(out)


def test_a_partial_from_a_crashed_run_is_never_mistaken_for_a_product(tmp_path, out):
    _, source, rtx_path, _ = build(tmp_path)
    out.parent.mkdir(parents=True)
    partial_path(out).write_bytes(b"half a file")
    with pytest.raises(PreflightError, match="partial file"):
        convert(source, rtx_path, out, small_config())
    assert not out.exists() and partial_path(out).read_bytes() == b"half a file"


def test_existing_output_is_not_overwritten_unless_requested(converted, tmp_path):
    out = converted["out"]
    before = out.read_bytes()
    with pytest.raises(PreflightError, match="already exists"):
        convert(converted["source"], converted["rtx"], out, small_config())
    assert out.read_bytes() == before
    convert(converted["source"], converted["rtx"], out, small_config(overwrite=True))
    assert out.is_file()


def test_unsafe_destinations_are_refused(tmp_path):
    _, source, rtx_path, _ = build(tmp_path)
    with pytest.raises(PreflightError, match="next to the original"):
        convert(source, rtx_path, tmp_path / "next_to_input.nc", small_config())
    with pytest.raises(PreflightError, match="'data'"):
        convert(source, rtx_path, tmp_path / "data" / "sub" / "x.nc", small_config())
    with pytest.raises(PreflightError, match="should end in .nc"):
        convert(source, rtx_path, tmp_path / "o" / "x.txt", small_config())
    with pytest.raises(PreflightError, match="one of the input files"):
        convert(source, rtx_path, rtx_path, small_config())
    assert rtx_path.exists()


def test_inputs_are_not_modified(tmp_path, out):
    _, source, rtx_path, _ = build(tmp_path)
    before = (hashlib.sha256(rtx_path.read_bytes()).hexdigest(), rtx_path.stat().st_mtime_ns)
    convert(source, rtx_path, out, small_config())
    assert (hashlib.sha256(rtx_path.read_bytes()).hexdigest(), rtx_path.stat().st_mtime_ns) == before


def test_lazy_structure_check_flags_a_missing_variable(converted, tmp_path):
    broken = tmp_path / "broken.nc"
    shutil.copy(converted["out"], broken)
    with netCDF4.Dataset(broken, "a") as ds:
        ds["acquisition"].renameVariable("line_counter", "line_counter_renamed")
    failed = [c for c in check_lazy_structure(broken, "uint8", 3) if c["status"] == "FAILED"]
    assert failed and "line_counter" in failed[0]["detail"]
