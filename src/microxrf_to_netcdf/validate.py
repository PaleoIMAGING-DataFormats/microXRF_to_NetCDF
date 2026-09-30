"""
Authors: Andre L. Belem (https://github.com/andrebelem) and F.R.I.D.A.Y. (https://observatoriooceanografico.org/people/friday-bot)

Validation of a converted file.

``check_lazy_structure`` is run by the converter on every output: it reopens the file lazily with xarray and
Dask (no data is loaded) and checks dimensions, dtypes, the absence of ``_FillValue`` on counts and the presence
of every schema variable.

``validate_against_sources`` is the independent, source-facing validation used by ``microxrf_to_netcdf validate``:

* the public RosettaSciIO path (``BCF_reader.parse_hypermap(downsample=8)``) against the counts binned 8 x 8
  over ALL channels, and the unmodified decoder on the real stream truncated to the first channels against the
  counts of EVERY pixel;
* the RTX payload decoder ``rtx_payload`` against every stored element map, video and mosaic plane;
* the source hashes recorded in the file against the input files;
* optionally the companion PNG images (element composite and rotated mosaic);
* optionally (``strict``) uint16 versus uint32 decoding of every band, which proves that no count wrapped.

Memory: one band of the counts (a few scan lines), one public-decoder result (about 55 MiB binned, about
350 MiB for the first channels of all pixels), one image plane. The cube is never materialized.
"""

from __future__ import annotations

import hashlib
import json
import struct
from pathlib import Path
from typing import Any, Callable

import netCDF4
import numpy as np

Log = Callable[[str], None]
REQUIRED_ACQUISITION = ("y", "x", "energy", "channel", "counts", "video", "element_maps", "element",
                        "element_line", "element_label", "energy_calibration", "spatial_calibration",
                        "sum_spectrum_header", "sum_spectrum_from_pixels", "line_counter", "instrument")
LOW_CHANNELS = 400
BIN = 8


def _check(name: str, passed: bool, detail: str) -> dict[str, Any]:
    return {"check": name, "status": "passed" if passed else "FAILED", "detail": detail}


def check_lazy_structure(path: Path, counts_dtype: str, n_elements: int) -> list[dict[str, Any]]:
    """Reopen ``path`` lazily with xarray and Dask and check the schema. Loads no data."""
    import xarray as xr
    checks: list[dict[str, Any]] = []
    with xr.open_dataset(path, group="acquisition", chunks={}) as acquisition:
        counts = acquisition["counts"]
        checks.append(_check("lazy_open_counts_is_dask_with_dims_y_x_energy",
                             hasattr(counts.data, "dask") and counts.dims == ("y", "x", "energy"),
                             f"{type(counts.data).__name__}, dims {counts.dims}, shape {counts.shape}, "
                             f"chunks {counts.data.chunksize}"))
        checks.append(_check("lazy_open_counts_dtype_is_unsigned_integer_as_planned",
                             str(counts.dtype) == counts_dtype and "_FillValue" not in counts.encoding
                             and "_FillValue" not in counts.attrs,
                             f"xarray dtype {counts.dtype} (planned {counts_dtype}); no _FillValue attribute"))
        missing = [name for name in REQUIRED_ACQUISITION if name not in acquisition.variables]
        checks.append(_check("acquisition_group_has_every_schema_variable", not missing, f"missing: {missing}"))
        maps = acquisition["element_maps"]
        checks.append(_check("lazy_open_element_maps_uint16_shared_grid",
                             str(maps.dtype) == "uint16" and maps.dims == ("element", "y", "x")
                             and maps.shape[0] == n_elements and maps.shape[1:] == counts.shape[:2]
                             and "units" not in maps.attrs,
                             f"{maps.dtype}, dims {maps.dims}, shape {maps.shape}; no units attribute (unknown)"))
    with xr.open_dataset(path, group="mosaic", chunks={}) as mosaic:
        pixels = mosaic["pixels"] if "pixels" in mosaic.variables else None
        checks.append(_check("mosaic_group_has_pixels_on_its_own_dimensions",
                             pixels is not None and pixels.dims == ("plane", "mosaic_y", "mosaic_x")
                             and str(pixels.dtype) == "uint8", f"{None if pixels is None else (pixels.dims, pixels.shape)}"))
    with netCDF4.Dataset(path) as ds:
        overview_groups = list(ds["overview"].groups)
    for name in overview_groups:
        with xr.open_dataset(path, group=f"overview/{name}", chunks={}) as overview:
            pixels = overview["pixels"]
            checks.append(_check(f"overview_{name}_lazy_open_uint8_own_grid",
                                 pixels.dims == ("plane", "row", "column") and str(pixels.dtype) == "uint8"
                                 and hasattr(pixels.data, "dask"), f"{pixels.dims} {pixels.shape}"))
    with xr.open_dataset(path, group="metadata") as metadata:
        needed = ("bcf_original_metadata_json", "rtx_payload_residual_xml", "provenance_json", "validation_json",
                  "timestamps_json", "rtx_image_table_json")
        missing = [name for name in needed if name not in metadata.variables]
        checks.append(_check("metadata_group_has_every_schema_variable", not missing, f"missing: {missing}"))
    return checks


def _reference_binned(bcf_path: Path) -> np.ndarray:
    from rsciio.bruker._api import BCF_reader
    return BCF_reader(str(bcf_path)).parse_hypermap(index=0, downsample=BIN, lazy=False).astype(np.int64)


def _reference_first_channels(bcf_path: Path, shape: tuple[int, int], channels: int) -> np.ndarray:
    from rsciio.bruker import unbcf_fast
    from rsciio.bruker._api import SFS_reader
    item = SFS_reader(str(bcf_path)).get_file("EDSDatabase/SpectrumData0")
    return unbcf_fast.parse_to_numpy(item, (*shape, channels), np.uint16)


def check_counts_against_public_decoder(nc_path: Path, bcf_path: Path, log: Log) -> list[dict[str, Any]]:
    checks = []
    with netCDF4.Dataset(nc_path) as ds:
        counts = ds["acquisition"]["counts"]
        counts.set_auto_maskandscale(False)
        height, width, channels = counts.shape
        log("validate: public RosettaSciIO decode (downsample=8, all channels)")
        if height % BIN == 0 and width % BIN == 0:
            reference = _reference_binned(bcf_path)
            binned = np.zeros(reference.shape, np.int64)
            sums = np.zeros(channels, np.int64)
            for y0 in range(0, height, BIN):
                block = np.asarray(counts[y0:y0 + BIN])
                sums += block.sum(axis=(0, 1), dtype=np.int64)
                binned[y0 // BIN] = block.reshape(BIN, width // BIN, BIN, channels).sum(axis=(0, 2), dtype=np.int64)
            checks.append(_check("counts_binned_8x8_all_channels_equal_public_rosettasciio_path",
                                 np.array_equal(binned, reference),
                                 f"shape {reference.shape}, total {int(reference.sum())}; every one of the {channels} "
                                 "channels enters the comparison"))
            stored = np.asarray(ds["acquisition"]["sum_spectrum_from_pixels"][:])
            checks.append(_check("stored_sum_spectrum_from_pixels_equals_sum_of_counts_read_back",
                                 np.array_equal(stored, sums), f"total {int(sums.sum())}"))
            del reference
        else:
            checks.append(_check("counts_binned_8x8_all_channels_equal_public_rosettasciio_path", False,
                                 f"grid {height} x {width} not divisible by {BIN}; check not run"))
        limit = min(LOW_CHANNELS, channels)
        log(f"validate: unmodified decoder on the real stream, first {limit} channels of every pixel")
        low = _reference_first_channels(bcf_path, (height, width), limit)
        equal = True
        for y0 in range(0, height, 4):
            if not np.array_equal(np.asarray(counts[y0:y0 + 4, :, :limit]).astype(np.uint16), low[y0:y0 + 4]):
                equal = False
                break
        checks.append(_check("counts_first_channels_of_every_pixel_equal_public_decoder", equal,
                             f"{height * width} pixels x {limit} channels compared exactly"))
        del low
    return checks


def check_rtx_planes(nc_path: Path, rtx_path: Path, log: Log) -> list[dict[str, Any]]:
    """Every RTX plane, decoded by ``rtx_payload``, against what is stored in the file."""
    from . import rtx_payload
    mismatches: list[str] = []
    compared = [0]
    with netCDF4.Dataset(nc_path) as ds:
        acquisition = ds["acquisition"]
        source_plane = list(np.asarray(acquisition["element_map_source_plane"][:]))
        mosaic_vars: dict[int, str] = {}
        for name, group in ds["mosaic"].groups.items():
            if group.source_file == "RTX payload":
                mosaic_vars[int(group.source_image_index)] = group.pixels_variable.split("/")[-1]

        def sink(image: int, plane: int, raw: bytes) -> None:
            compared[0] += 1
            if image in mosaic_vars:
                stored = np.asarray(ds["mosaic"][mosaic_vars[image]][plane])
                ok = np.array_equal(stored.reshape(-1), np.frombuffer(raw, "u1"))
            else:
                data = np.frombuffer(raw, "<u2")
                if plane in source_plane:
                    stored = np.asarray(acquisition["element_maps"][source_plane.index(plane)])
                else:
                    stored = np.asarray(acquisition["video"][:])
                ok = np.array_equal(stored.reshape(-1), data)
            if not ok:
                mismatches.append(f"image {image} plane {plane}")

        report = rtx_payload.inspect_rtx_payload(rtx_path, sink)
    return [_check("every_rtx_plane_equals_the_diagnostic_decoder_output", not mismatches and report["error"] is None,
                   f"{compared[0]} planes compared; error {report['error']}; differing {mismatches}")]


def check_bcf_header_images(nc_path: Path, bcf_path: Path, log: Log) -> list[dict[str, Any]]:
    """Every image in the BCF header (RosettaSciIO returns one), decoded again by the payload decoder
    ``rtx_payload`` straight from the SFS blocks, against where the file stored it."""
    import re

    from . import rtx_payload
    from rsciio.bruker._api import SFS_reader

    mismatches: list[str] = []
    compared = [0]
    with netCDF4.Dataset(nc_path) as ds:
        table = json.loads(ds["metadata"]["bcf_image_table_json"][...])
        roles = {row["header_image_index"]: row["role"] for row in table}

        def stored_plane(index: int, plane: int):
            role = roles[index]
            found = re.search(r"/acquisition/(\w+)$", role)
            if found:
                return np.asarray(ds["acquisition"][found.group(1)][:])
            found = re.search(r"/mosaic/(\w+)$", role)
            if found:
                variable = ds["mosaic"][found.group(1)].getncattr("pixels_variable").split("/")[-1]
                return np.asarray(ds["mosaic"][variable][plane])
            found = re.search(r"/overview/(\w+)$", role)
            if found:
                return np.asarray(ds["overview"][found.group(1)]["pixels"][plane])
            if "see /acquisition/video" in role:
                return np.asarray(ds["acquisition"]["video"][:])
            return None

        def sink(image: int, plane: int, raw: bytes) -> None:
            compared[0] += 1
            stored = stored_plane(image, plane)
            itemsize = {1: "u1", 2: "<u2", 4: "<u4"}
            if stored is None or not np.array_equal(
                    stored.reshape(-1).astype(itemsize[stored.dtype.itemsize]), np.frombuffer(
                        raw, itemsize[stored.dtype.itemsize])):
                mismatches.append(f"header image {image} plane {plane}")

        item = SFS_reader(str(bcf_path)).get_file("EDSDatabase/HeaderData")
        payload = rtx_payload.parse_payload(item.get_iter_and_properties()[0], sink)
        expected_planes = sum(len(im["planes"]) for im in payload["images"])
    log(f"validate: {compared[0]} BCF header planes compared")
    return [_check("every_bcf_header_image_plane_equals_the_diagnostic_decoder_output",
                   not mismatches and compared[0] == expected_planes > 0,
                   f"{compared[0]} of {expected_planes} planes in {len(payload['images'])} header images compared "
                   f"(RosettaSciIO returns only one image); differing: {mismatches}")]


def check_provenance(nc_path: Path, bcf_path: Path, rtx_path: Path) -> list[dict[str, Any]]:
    def digest(path: Path) -> str:
        h = hashlib.sha256()
        with path.open("rb") as source:
            for block in iter(lambda: source.read(1 << 20), b""):
                h.update(block)
        return h.hexdigest()

    checks = []
    with netCDF4.Dataset(nc_path) as ds:
        checks.append(_check("file_is_marked_complete", ds.conversion_status == "complete", str(ds.conversion_status)))
        checks.append(_check("recorded_source_hashes_equal_input_files",
                             ds.source_bcf_sha256 == digest(bcf_path) and ds.source_rtx_sha256 == digest(rtx_path),
                             f"bcf {ds.source_bcf_sha256[:16]}..., rtx {ds.source_rtx_sha256[:16]}..."))
        checks.append(_check("recorded_source_sizes_equal_input_files",
                             int(ds.source_bcf_size_bytes) == bcf_path.stat().st_size
                             and int(ds.source_rtx_size_bytes) == rtx_path.stat().st_size, "sizes match"))
        provenance = json.loads(ds["metadata"]["provenance_json"][...])
        checks.append(_check("provenance_lists_software_and_parameters",
                             {"microxrf_to_netcdf", "rosettasciio", "netCDF4", "numpy"} <= set(provenance["software"])
                             and "band_lines" in provenance["parameters"], json.dumps(provenance["software"])[:200]))
    return checks


def _pearson(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.corrcoef(a.ravel().astype(np.float64), b.ravel().astype(np.float64))[0, 1])


def check_companion_images(nc_path: Path, data_dir: Path) -> list[dict[str, Any]]:
    """Mosaic orientation (rotation) and the Ca/Fe composite against the supplied PNG files."""
    from PIL import Image
    checks = []
    mosaic_png, cafe_png = data_dir / "Video Mosaic.png", data_dir / "CaFe.png"
    with netCDF4.Dataset(nc_path) as ds:
        if mosaic_png.is_file():
            png = np.asarray(Image.open(mosaic_png).convert("RGB"))
            best = []
            for name in ds["mosaic"].variables:
                if name.startswith("pixels"):
                    planes = np.asarray(ds["mosaic"][name][:])
                    rotated = [np.rot90(planes[c], k=-1) for c in range(3)]
                    if rotated[0].shape != png.shape[:2]:
                        best.append((name, "shape mismatch"))
                        continue
                    best.append((name, [round(_pearson(rotated[c], png[:, :, c]), 4) for c in range(3)]))
            ok = bool(best) and all(isinstance(r, list) and min(r) > 0.9 for _, r in best)
            checks.append(_check("mosaic_pixels_rotated_90_clockwise_correlate_with_video_mosaic_png", ok,
                                 f"Pearson r per RGB channel: {best}"))
        if cafe_png.is_file():
            png = np.asarray(Image.open(cafe_png).convert("RGB"))
            labels = list(ds["acquisition"]["element_label"][:])
            maps = ds["acquisition"]["element_maps"]
            if "Ca-KA" in labels and "Fe-KA" in labels:
                ca = np.asarray(maps[labels.index("Ca-KA")])
                fe = np.asarray(maps[labels.index("Fe-KA")])
                r_ca, r_fe = _pearson(ca, png[:, :, 2]), _pearson(fe, png[:, :, 0])
                checks.append(_check("element_maps_ca_fe_correlate_with_cafe_png_blue_and_red", r_ca > 0.9 and r_fe > 0.85,
                                     f"Ca-KA vs blue r {r_ca:.4f}; Fe-KA vs red r {r_fe:.4f}"))
    return checks


def check_decode_dtype_agreement(bcf_path: Path, band_lines: int, log: Log) -> list[dict[str, Any]]:
    """Decode every band as uint16 and as uint32 and compare: proves that no uint16 count wrapped."""
    from .bcf import decode_band, frame_lines, open_bcf
    source = open_bcf(bcf_path)
    info = source.info
    framer = frame_lines(source.blocks(), info.height, info.width)
    header = next(framer)
    pending, bands, maximum = [], 0, 0
    ok = True
    for record in framer:
        pending.append(record.data)
        if len(pending) == band_lines or record.index == info.height - 1:
            small = decode_band(header, pending, info.width, info.channels, np.uint16)
            large = decode_band(header, pending, info.width, info.channels, np.uint32)
            ok &= bool(np.array_equal(small, large))
            maximum = max(maximum, int(large.max()))
            bands += 1
            pending = []
    log(f"validate: uint16 vs uint32 decode of {bands} bands")
    return [_check("uint16_and_uint32_decodes_identical_no_wraparound", ok,
                   f"{bands} bands, maximum count {maximum}")]


def validate_against_sources(nc_path: Path, bcf_path: Path, rtx_path: Path, data_dir: Path | None = None,
                             strict: bool = False, log: Log = lambda _: None) -> list[dict[str, Any]]:
    nc_path, bcf_path, rtx_path = Path(nc_path), Path(bcf_path), Path(rtx_path)
    with netCDF4.Dataset(nc_path) as ds:
        counts_dtype = str(ds["acquisition"]["counts"].dtype)
        n_elements = len(ds["acquisition"].dimensions["element"])
        band_lines = int(json.loads(ds["metadata"]["provenance_json"][...])["parameters"]["band_lines"])
    checks = check_lazy_structure(nc_path, counts_dtype, n_elements)
    checks += check_provenance(nc_path, bcf_path, rtx_path)
    checks += check_counts_against_public_decoder(nc_path, bcf_path, log)
    checks += check_rtx_planes(nc_path, rtx_path, log)
    checks += check_bcf_header_images(nc_path, bcf_path, log)
    if data_dir is not None:
        checks += check_companion_images(nc_path, Path(data_dir))
    if strict:
        checks += check_decode_dtype_agreement(bcf_path, band_lines, log)
    return checks
