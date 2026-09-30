"""
Authors: Andre L. Belem (https://github.com/andrebelem) and F.R.I.D.A.Y. (https://observatoriooceanografico.org/people/friday-bot)

Large-file tests: convert the ORIGINAL BCF + RTX in data/ (read-only) into a temporary NetCDF-4 file by running
the command line in a child process, then validate the result lazily against independent references (the public
RosettaSciIO path, the prototype framer of tools/benchmarks/probe_bcf_streaming.py, the RTX payload decoder and the
companion PNG files). No test loads the complete EDS cube: reads are bands, lines, windows or single pixels.
Skipped cleanly when the original files are absent.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import netCDF4
import numpy as np
import pytest
import xarray as xr

from conftest import BCF_PATH, DATA_DIR, RTX_PATH

pytestmark = pytest.mark.large
ROOT = Path(__file__).resolve().parent.parent
CUBE_BYTES = 240 * 1800 * 4096


@pytest.fixture(scope="module")
def real(tmp_path_factory, bcf_path, rtx_path):
    before = {p: (p.stat().st_size, p.stat().st_mtime_ns) for p in (bcf_path, rtx_path)}
    out = tmp_path_factory.mktemp("real_conversion") / "GRF17A_9-29cm_slab3.nc"
    result = subprocess.run(
        [sys.executable, "-m", "microxrf_to_netcdf", "convert", "--bcf", str(bcf_path), "--rtx", str(rtx_path),
         "--out", str(out), "--quiet"], cwd=ROOT, capture_output=True, text=True, timeout=900)
    assert result.returncode == 0, result.stderr[-3000:]
    summary = json.loads(result.stdout)
    report = json.loads(out.with_name(out.name + ".report.json").read_text(encoding="utf-8"))
    return {"out": out, "summary": summary, "report": report, "before": before}


def test_converter_ran_with_bounded_memory_far_below_the_cube_size(real):
    peak_mib = real["summary"]["peak_working_set_mib"]
    assert peak_mib is not None and peak_mib < 800, peak_mib     # measured about 380 MiB for 4-line bands
    assert peak_mib * 2**20 < 0.5 * CUBE_BYTES                   # the cube (1.77 GB as uint8) was never materialized


def test_output_is_complete_and_recorded_correctly(real):
    out, report = real["out"], real["report"]
    assert out.is_file() and not out.with_name(f".{out.name}.partial").exists()
    assert report["output_size_bytes"] == out.stat().st_size < CUBE_BYTES // 4   # compressed; measured about 212 MiB
    assert all(c["status"] == "passed" for c in report["validation"])
    assert report["counts_dtype"] == "uint8"
    counts = report["actual_layout"]["acquisition/counts"]
    assert counts["chunks"] == [2, 30, 4096] and counts["filters"]["zlib"] and counts["filters"]["complevel"] == 4
    assert report["plan"]["status"].startswith("ESTIMATES")


def test_inputs_are_unchanged(real):
    for path, stamp in real["before"].items():
        assert (path.stat().st_size, path.stat().st_mtime_ns) == stamp


def test_dimensions_variables_and_lazy_reopen(real):
    with xr.open_dataset(real["out"], group="acquisition", chunks={}) as acquisition:
        counts = acquisition["counts"]
        assert counts.dims == ("y", "x", "energy") and counts.shape == (240, 1800, 4096)
        assert hasattr(counts.data, "dask") and counts.dtype == np.uint8 and counts.data.chunksize == (2, 30, 4096)
        maps = acquisition["element_maps"]
        assert maps.shape == (14, 240, 1800) and maps.dtype == np.uint16 and "units" not in maps.attrs
        assert list(acquisition["element_label"].values) == [
            "Ca-KA", "K-KA", "S-KA", "Si-K", "Ti-KA", "Cr-KA", "Mn-KA", "Fe-KA", "Ni-KA", "Zn-KA", "Sr-KA", "Pd-KA",
            "Rh-KA", "Al-K"]
        assert acquisition["video"].shape == (240, 1800)
    with xr.open_dataset(real["out"], group="mosaic", chunks={}) as mosaic:
        assert mosaic["pixels"].shape == (3, 948, 6000) and set(mosaic.dims) >= {"mosaic_y", "mosaic_x", "plane"}
    with netCDF4.Dataset(real["out"]) as ds:
        assert set(ds.groups) == {"acquisition", "mosaic", "overview", "metadata"}
        assert set(ds["mosaic"].groups) == {"instance_rtx_0", "instance_rtx_1", "instance_bcf_3_Default"}
        assert set(ds["overview"].groups) == {"Image_0", "Image_1", "Image_2"}


def test_spectral_and_spatial_calibration_equal_the_public_reader(real, bcf_path):
    from rsciio.bruker import file_reader
    edx = next(s for s in file_reader(str(bcf_path), lazy=True, select_type="spectrum_image"))
    axes = {a["name"]: a for a in edx["axes"]}
    with netCDF4.Dataset(real["out"]) as ds:
        acq = ds["acquisition"]
        calibration = acq["energy_calibration"]
        assert calibration.calib_abs_raw == axes["Energy"]["offset"] and calibration.calib_lin_raw == axes["Energy"]["scale"]
        assert np.allclose(acq["energy"][:], axes["Energy"]["offset"] + axes["Energy"]["scale"] * np.arange(4096),
                           rtol=0, atol=1e-12)
        assert acq["energy"].units == "keV"
        assert acq["y"][1] == pytest.approx(axes["height"]["scale"], rel=1e-12)
        assert acq["x"][1] == pytest.approx(axes["width"]["scale"], rel=1e-12)


def test_every_pixel_spectrum_equals_the_diagnostic_prototype_decoder(real, bcf_path):
    """All 240 lines, all 4096 channels, against probe_bcf_streaming.py (an independent framer, same decoder)."""
    import probe_bcf_streaming as probe
    _, item = probe.open_stream(bcf_path)
    with netCDF4.Dataset(real["out"]) as ds:
        counts = ds["acquisition"]["counts"]
        counts.set_auto_maskandscale(False)
        lines = 0
        for first, band in probe.stream_bands(item, 4, 4096, np.uint16):
            stored = np.asarray(counts[first:first + band.shape[0]])
            assert stored.dtype == np.uint8 and np.array_equal(stored, band.astype(np.uint8)), first
            lines += band.shape[0]
        assert lines == 240


def test_valid_zero_counts_and_total_counts(real):
    with xr.open_dataset(real["out"], group="acquisition", chunks={}) as acquisition:
        band = acquisition["counts"].isel(y=slice(100, 104)).values
        assert band.dtype == np.uint8 and int((band == 0).sum()) > 0.5 * band.size   # zero counts are stored as 0
        total = int(acquisition["sum_spectrum_from_pixels"].sum())
        assert total == 677_896_870                       # FINDINGS.md 7.2, streamed decode of the real file
        assert int(acquisition["sum_spectrum_header"].sum()) == 680_903_120   # recorded in the BCF header
        assert int(acquisition["counts"].isel(y=slice(0, 8)).max().values) <= 149


def test_representative_pixel_spectrum_window_and_channel_image(real):
    import probe_bcf_streaming as probe
    _, item = probe.open_stream(BCF_PATH)
    reference = {}
    for first, band in probe.stream_bands(item, 1, 4096, np.uint16):
        if first in (50, 100):
            reference[first] = band[0].copy()
        if first == 100:
            break
    with xr.open_dataset(real["out"], group="acquisition", chunks={}) as acquisition:
        array = acquisition["counts"]
        assert np.array_equal(array.isel(y=100, x=900).values, reference[100][900])
        assert np.array_equal(array.isel(y=slice(50, 51), x=slice(1000, 1100), energy=slice(0, 400)).values[0],
                              reference[50][1000:1100, :400])
        assert np.array_equal(array.isel(y=100, energy=468).values, reference[100][:, 468])


def test_element_maps_video_and_mosaic_equal_the_rtx_decoder(real, rtx_decoded, video_signal):
    with netCDF4.Dataset(real["out"]) as ds:
        acq = ds["acquisition"]
        for k in range(14):
            assert np.array_equal(np.asarray(acq["element_maps"][k]), rtx_decoded["planes"][(2, k + 1)]), k
        video = np.asarray(acq["video"][:])
        assert np.array_equal(video, np.asarray(video_signal["data"]))       # BCF Video dataset
        assert np.array_equal(video, rtx_decoded["planes"][(2, 0)])          # RTX Mapdaten plane 0
        for plane in range(3):
            assert np.array_equal(np.asarray(ds["mosaic"]["pixels"][plane]), rtx_decoded["planes"][(1, plane)])


def test_mosaic_orientation_and_element_composite_match_the_png_files(real):
    from microxrf_to_netcdf.validate import check_companion_images
    checks = check_companion_images(real["out"], DATA_DIR)
    assert len(checks) == 2 and all(c["status"] == "passed" for c in checks), checks


def test_mosaic_footprint_relationship_is_recorded(real):
    with netCDF4.Dataset(real["out"]) as ds:
        verified = ds["mosaic"]["instance_rtx_1"]["map_footprint"]
        assert (verified.rect_left, verified.rect_top, verified.rect_right, verified.rect_bottom) == (192, 89, 5941, 855)
        assert verified.registration_status == "verified"
        assert json.loads(verified.registration_json)["correlation_by_orientation"]["identity"] > 0.99
        other = ds["mosaic"]["instance_rtx_0"]["map_footprint"]
        assert other.registration_status == "not_verified" and other.size_matches == 0


def test_metadata_and_provenance(real):
    with netCDF4.Dataset(real["out"]) as ds:
        assert ds.conversion_status == "complete" and ds.microxrf_to_netcdf_schema_version == "1.0.0"
        assert ds.source_bcf_size_bytes == 617_320_728 and ds.source_rtx_size_bytes == 43_723_950
        assert len(ds.source_bcf_sha256) == 64 and len(ds.source_rtx_sha256) == 64
        instrument = ds["acquisition"]["instrument"]
        assert instrument.detector_type == "XFlash 430" and instrument.sample_name == "Mapdaten"
        assert instrument.beam_energy_recorded == 50.0 and "NOT stored" in instrument.beam_energy_units_status
        assert "units" not in instrument.ncattrs()
        metadata = ds["metadata"]
        original = json.loads(metadata["bcf_original_metadata_json"][...])
        assert original["Spectrum"]["CalibLin"] == 0.010001 and original["Detector"]["Type"] == "XFlash 430"
        stamps = json.loads(metadata["timestamps_json"][...])
        assert {s["time"] for s in stamps} >= {"09:55:52", "15:58:47", "9:38:26", "9:55:51", "9:55:53"}
        residual = metadata["rtx_payload_residual_xml"][...]
        assert "HV: 50,0 kV" in residual and "LineCounter" in residual and "40000 µm" in residual
        assert len(residual) < 1_000_000                                # pixel text was removed
        assert metadata["rtx_payload_residual_xml"].payload_sha256 == \
            "832f6f4c090e5413033621335127f7095ceb9ff7f9f7c3604bca2ffba00d4b61"
        provenance = json.loads(metadata["provenance_json"][...])
        assert provenance["software"]["rosettasciio"] == "0.14.0" and provenance["bcf_scan"]["maximum_count"] == 149


def test_independent_source_validation_passes(real, bcf_path, rtx_path):
    """The validate command's checks: public RosettaSciIO path (binned 8 x 8, all channels; first 400 channels of
    every pixel), the diagnostic RTX decoder, source hashes and the PNG files."""
    from microxrf_to_netcdf.validate import validate_against_sources
    checks = validate_against_sources(real["out"], bcf_path, rtx_path, DATA_DIR)
    failed = [c for c in checks if c["status"] != "passed"]
    assert not failed, failed
    names = {c["check"] for c in checks}
    assert "counts_binned_8x8_all_channels_equal_public_rosettasciio_path" in names
    assert "counts_first_channels_of_every_pixel_equal_public_decoder" in names
    assert "every_rtx_plane_equals_the_diagnostic_decoder_output" in names


def test_the_lazy_opening_example_of_the_schema_document_runs_on_the_real_output(real):
    """Execute the python block of NETCDF_SCHEMA.md section 13 against the converted file."""
    import re
    text = (ROOT / "NETCDF_SCHEMA.md").read_text(encoding="utf-8")
    section = text[text.index("## 13. Opening the file lazily"):]
    code = re.search(r"```python\n(.*?)```", section, re.S).group(1)
    code = code.replace('"output/GRF17A_9-29cm_slab3.nc"', repr(str(real["out"])))
    namespace: dict = {}
    exec(compile(code, "NETCDF_SCHEMA.md", "exec"), namespace)
    assert namespace["counts"].dtype == np.uint8 and namespace["spectrum"].shape == (4096,)
    assert namespace["window"].shape == (8, 100, 400) and namespace["ca_map"].shape == (240, 1800)
    assert namespace["roi"].dims == ("y", "x") and namespace["mosaic"].shape == (3, 948, 6000)
    assert 0 < float(namespace["ca_map"].max()) <= 1607          # RTX Ca-KA maximum (FINDINGS 8.3)


def test_bcf_header_images_that_rosettasciio_does_not_return_are_in_the_file(real, rtx_decoded):
    """FINDINGS.md section 9: the BCF header holds 7 images; RosettaSciIO returns only the video."""
    with netCDF4.Dataset(real["out"]) as ds:
        pixel_times = ds["acquisition"]["pixel_times"]
        assert pixel_times.dtype == np.uint32 and pixel_times.shape == (240, 1800) and "units" not in pixel_times.ncattrs()
        values = np.asarray(pixel_times[:])
        assert (int(values.min()), int(values.max()), int(values.sum(dtype=np.uint64))) == (14500, 31400, 12_955_261_600)
        mosaic = ds["mosaic"]
        bcf = mosaic["instance_bcf_3_Default"]
        assert bcf.pixels_duplicate_of == "rtx:0" and bcf.source_time == "9:55:51"
        assert [v for v in mosaic.variables if v.startswith("pixels")] == ["pixels"]     # one copy for three instances
        assert json.loads(mosaic["pixels"].stored_once_for_instances) == ["rtx:0", "rtx:1", "bcf:3"]
        sizes = {name: ds["overview"][name]["pixels"].shape for name in ds["overview"].groups}
        assert sizes == {"Image_0": (3, 768, 1024), "Image_1": (3, 768, 1024), "Image_2": (3, 480, 752)}
        table = {row["name"]: row for row in json.loads(ds["metadata"]["bcf_image_table_json"][...])}
        assert len(table) == 7 and table["Counter"]["plane_count"] == 0
        assert ds["metadata"]["bcf_header_residual_xml"].header_decoded_bytes == 34_006_136
        assert len(ds["metadata"]["bcf_header_residual_xml"][...]) < 100_000
