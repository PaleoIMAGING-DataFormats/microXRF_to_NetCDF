"""
Authors: Andre L. Belem (https://github.com/andrebelem) and F.R.I.D.A.Y. (https://observatoriooceanografico.org/people/friday-bot)

End-to-end BCF + RTX -> NetCDF-4 conversion. One paired acquisition gives one file. Order of work:

    1. inspect: BCF header (never spectra), RTX scan pass 1, grid validation, source hashes
    2. scan the BCF stream once (max count, per-channel sums, largest record): decides the output dtype
    3. preflight: logical sizes, chunks, disk space, memory (PreflightError before anything is created)
    4. write to ``.<name>.partial``: RTX planes (pass 2), registration check, BCF bands (pass 3)
    5. close, reopen read-only, compare every band, plane and coordinate with the sources (pass 4)
    6. reopen for append: provenance, validation and ``conversion_status = complete``
    7. lazy-open check, SHA-256, atomic ``os.replace`` to the final name, sidecar report

Any failure (including KeyboardInterrupt) removes the partial file and re-raises; the final name is created only
by the last atomic step, so a file at the final name has passed the checks above. Inputs are never opened for
writing; their size and modification time are compared before and after.

Memory bound: see preflight.py. The complete cube is never materialized: one band of ``band_lines`` scan lines
is decoded, converted and written at a time, in every pass.
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
import sys
import time
import warnings
from dataclasses import replace
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path
from typing import Any, Callable

import netCDF4
import numpy as np

from . import SCHEMA_VERSION, __version__
from .bcf import BCFSource, RealBCF, iter_bands, open_bcf, stream_header_planes
from .config import ConversionConfig
from .errors import MicroXRFToNetCDFError, ConversionError, PreflightError, ValidationError
from .memory import peak_working_set_mib
from .model import Layout, build_layout, video_sha256
from .preflight import BCFScan, Plan, check_destination, check_scan_memory, make_plan, partial_path
from .registration import correlate_footprint
from .rtx import RTXScan, scan_rtx, stream_rtx_planes
from .validate import check_lazy_structure
from .writer import UnifiedWriter, file_sha256

Log = Callable[[str], None]


def _no_log(_: str) -> None:
    pass


def _stamp(path: Path) -> tuple[int, int]:
    stat = path.stat()
    return stat.st_size, stat.st_mtime_ns


def _accumulate(stats: BCFScan, band) -> None:
    data = band.data
    stats.maximum = max(stats.maximum, int(data.max()))
    sums = data.sum(axis=(0, 1), dtype=np.int64)
    stats.channel_sums = sums if stats.channel_sums is None else stats.channel_sums + sums
    stats.total_counts += int(sums.sum())
    stats.nonzero_values += int(np.count_nonzero(data))
    stats.lines += data.shape[0]
    stats.widened_bands += int(band.widened)


def scan_band_lines(config: ConversionConfig) -> int:
    """Lines per band in the scan pass: it decodes in uint32 (twice the itemsize of uint16), so it uses half as many
    lines and needs about the same memory as the write pass."""
    return max(1, config.band_lines // 2) if config.decode_dtype == "uint16" else config.band_lines


def scan_bcf(source: BCFSource, config: ConversionConfig, log: Log = _no_log) -> BCFScan:
    """Pass over the BCF stream that keeps only statistics (bounded memory: one band).

    It decodes in uint32, which is exact, so the maximum it reports is the true maximum count. That maximum
    decides the output dtype and whether the later passes may decode in uint16 (see bcf.iter_bands).
    """
    started = time.perf_counter()
    stats = BCFScan()

    def on_line(record) -> None:
        stats.pixels += record.pixels
        stats.largest_line_bytes = max(stats.largest_line_bytes, len(record.data))

    for band in iter_bands(source, scan_band_lines(config), "uint32", config.max_line_bytes, on_line):
        _accumulate(stats, band)
    stats.seconds = round(time.perf_counter() - started, 2)
    log(f"scan: {stats.lines} lines, {stats.pixels} pixels, max count {stats.maximum}, "
        f"total counts {stats.total_counts}, {stats.seconds} s")
    return stats


def _same_stats(a: BCFScan, b: BCFScan) -> bool:
    return (a.lines == b.lines and a.maximum == b.maximum
            and a.total_counts == b.total_counts and a.nonzero_values == b.nonzero_values
            and np.array_equal(a.channel_sums, b.channel_sums))


def _source_identity(source: BCFSource, rtx_path: Path) -> dict[str, Any]:
    real = isinstance(source, RealBCF)
    info = source.info
    return {
        "acquisition_id": Path(info.filename).stem,
        "bcf_sha256": file_sha256(source.path) if real else "synthetic-source (no file)",
        "rtx_sha256": file_sha256(rtx_path),
    }


def _software_versions() -> dict[str, str]:
    versions = {"microxrf_to_netcdf": __version__, "schema": SCHEMA_VERSION, "python": platform.python_version(),
                "platform": platform.platform()}
    for package in ("rosettasciio", "netCDF4", "numpy", "xarray", "dask", "pillow"):
        try:
            versions[package] = version(package)
        except Exception:  # noqa: BLE001 - a missing optional package is recorded, not fatal
            versions[package] = "not installed"
    versions["hdf5"] = str(getattr(netCDF4, "__hdf5libversion__", "unknown"))
    versions["netcdf_c"] = str(getattr(netCDF4, "__netcdf4libversion__", "unknown"))
    return versions


def _check(name: str, passed: bool, detail: str) -> dict[str, Any]:
    return {"check": name, "status": "passed" if passed else "FAILED", "detail": detail}


def verify_readback(path: Path, source: BCFSource, layout: Layout, plan: Plan, config: ConversionConfig,
                    expected_sums: np.ndarray, decode_dtype: str, log: Log = _no_log) -> list[dict[str, Any]]:
    """Pass 4: reopen read-only and compare EVERY band, plane and coordinate with the sources."""
    info = source.info
    checks: list[dict[str, Any]] = []
    dtype = np.dtype(plan.counts_dtype)
    with netCDF4.Dataset(path, "r") as ds:
        acquisition = ds["acquisition"]
        counts = acquisition["counts"]
        counts.set_auto_maskandscale(False)
        chunking = counts.chunking()
        read_lines = int(chunking[0]) if isinstance(chunking, (list, tuple)) else 1
        mismatched_bands = 0
        bands = 0
        for band in iter_bands(source, config.band_lines, decode_dtype, config.max_line_bytes):
            # Read chunk-aligned (the y chunk extent per call, so no chunk is decompressed twice) and compare line by
            # line: values compare across dtypes, so no band-sized copy or boolean temporary is made.
            same = True
            for start in range(0, band.data.shape[0], read_lines):
                stored = counts[band.first_line + start:band.first_line + min(start + read_lines, band.data.shape[0])]
                for k in range(stored.shape[0]):
                    same = same and stored.dtype == dtype and bool(np.array_equal(stored[k], band.data[start + k]))
            if not same:
                mismatched_bands += 1
            bands += 1
        checks.append(_check("every_written_band_equals_a_fresh_decode", mismatched_bands == 0,
                             f"{bands} bands compared exactly (dtype {dtype}); {mismatched_bands} differ"))
        stored_sums = np.asarray(acquisition["sum_spectrum_from_pixels"][:])
        checks.append(_check("stored_pixel_sum_spectrum_matches_stream", np.array_equal(stored_sums, expected_sums),
                             f"total {int(stored_sums.sum())}"))
        checks.append(_check("counts_dtype_is_integer_and_unchanged", counts.dtype == dtype,
                             f"stored dtype {counts.dtype}, planned {dtype}; no _FillValue: "
                             f"{'_FillValue' not in counts.ncattrs()}"))
        video = np.asarray(acquisition["video"][:])
        checks.append(_check("video_equals_bcf_video", np.array_equal(video, info.video), f"shape {video.shape}"))
        digest_ok = video_sha256(video) == layout.video_plane.sha256
        checks.append(_check("video_equals_rtx_video_plane_sha256", digest_ok, layout.video_plane.sha256))
        element_maps = acquisition["element_maps"]
        element_maps.set_auto_maskandscale(False)
        bad = [p.description for k, p in enumerate(layout.element_planes)
               if hashlib.sha256(np.ascontiguousarray(element_maps[k], dtype="<u2").tobytes()).hexdigest() != p.sha256]
        checks.append(_check("element_maps_equal_rtx_planes_sha256", not bad,
                             f"{len(layout.element_planes)} maps compared; differing: {bad}"))
        mosaic = ds["mosaic"]
        bad_planes = []
        for entry in layout.mosaics:
            variable = mosaic[entry.pixels_variable]
            variable.set_auto_maskandscale(False)
            for plane in entry.image.planes:
                data = np.ascontiguousarray(variable[plane.index], dtype="u1")
                if hashlib.sha256(data.tobytes()).hexdigest() != plane.sha256:
                    bad_planes.append((entry.image.index, plane.index))
        checks.append(_check("mosaic_planes_of_every_instance_equal_the_source_plane_sha256", not bad_planes,
                             f"{sum(len(m.image.planes) for m in layout.mosaics)} planes of "
                             f"{len(layout.mosaics)} instances (RTX and BCF header) compared; differing: {bad_planes}"))
        bad_aux = []
        for aux in layout.aux_images:
            stored = np.asarray(acquisition[aux.variable][:])
            dtype_of = np.dtype({1: "u1", 2: "<u2", 4: "<u4"}[aux.image.itemsize])
            if hashlib.sha256(np.ascontiguousarray(stored, dtype=dtype_of).tobytes()).hexdigest() != aux.image.planes[0].sha256:
                bad_aux.append(aux.variable)
        checks.append(_check("acquisition_grid_header_images_equal_the_source_plane_sha256", not bad_aux,
                             f"{[a.variable for a in layout.aux_images]} compared; differing: {bad_aux}"))
        bad_overview = []
        for overview in layout.overviews:
            variable = ds["overview"][overview.group_name]["pixels"]
            variable.set_auto_maskandscale(False)
            for plane in overview.image.planes:
                data = np.ascontiguousarray(variable[plane.index], dtype="u1")
                if hashlib.sha256(data.tobytes()).hexdigest() != plane.sha256:
                    bad_overview.append((overview.group_name, plane.index))
        checks.append(_check("overview_images_equal_the_source_plane_sha256", not bad_overview,
                             f"{[o.group_name for o in layout.overviews]} compared; differing: {bad_overview}"))
        energy = np.asarray(acquisition["energy"][:])
        expected_energy = info.calib_abs + info.calib_lin * np.arange(info.channels)
        checks.append(_check("energy_coordinate_equals_raw_calibration", np.array_equal(energy, expected_energy),
                             f"calib_abs {info.calib_abs!r}, calib_lin {info.calib_lin!r}"))
        ys, xs = np.asarray(acquisition["y"][:]), np.asarray(acquisition["x"][:])
        checks.append(_check("spatial_coordinates_equal_raw_pixel_sizes",
                             np.array_equal(ys, np.arange(info.height) * info.pixel_size_y)
                             and np.array_equal(xs, np.arange(info.width) * info.pixel_size_x),
                             f"pixel_size_y {info.pixel_size_y!r}, pixel_size_x {info.pixel_size_x!r}"))
    log("readback: " + ", ".join(f"{c['check']}={c['status']}" for c in checks))
    return checks


def convert(bcf: str | Path | BCFSource, rtx: str | Path, out: str | Path, config: ConversionConfig | None = None,
            log: Log = _no_log, _fault: Callable[[str, int], None] | None = None,
            dry_run: bool = False) -> dict[str, Any]:
    """Convert one paired BCF/RTX acquisition into ``out`` and return the conversion report.

    ``_fault(stage, index)`` is a test hook called before each band of the write pass ("band", i) and at the
    named stages; it may raise to simulate an interruption.
    """
    started = time.perf_counter()
    config = config or ConversionConfig()
    config.validate()
    rtx_path, out = Path(rtx), Path(out).resolve()
    source = open_bcf(bcf) if isinstance(bcf, (str, Path)) else bcf
    info = source.info
    inputs = [rtx_path] + ([source.path] if isinstance(source, RealBCF) else [])
    check_destination(out, inputs, config.overwrite or dry_run)
    before = {p: _stamp(p) for p in inputs}
    log(f"inspect: BCF {info.filename} {info.height} x {info.width} x {info.channels}; RTX {rtx_path.name}")

    identity = _source_identity(source, rtx_path)
    scan = scan_rtx(rtx_path)
    layout = build_layout(info, scan)
    log("grid correspondence validated: " + "; ".join(c["check"] for c in layout.checks))
    check_scan_memory(info, config, scan_band_lines(config))
    stats = scan_bcf(source, config, log)
    header_sum = info.sum_spectrum.astype(np.int64)
    over = np.nonzero(stats.channel_sums > header_sum)[0]
    if len(over):
        raise ValidationError(f"the decoded pixel spectra exceed the BCF header sum spectrum in {len(over)} channels "
                              f"(first {int(over[0])}); a decoder fault is possible")
    # uint16 decoding in the write and read-back passes is allowed only if the exact (uint32) maximum fits.
    decode_dtype = config.decode_dtype if stats.maximum <= 65535 else "uint32"
    if decode_dtype != config.decode_dtype:
        log(f"maximum count {stats.maximum} exceeds uint16: the write and read-back passes decode in uint32")
    plan = make_plan(info, stats, layout, scan, replace(config, decode_dtype=decode_dtype), out)
    for warning in plan.warnings:
        log("WARNING: " + warning)
    log(f"plan: dtype {plan.counts_dtype}; logical {plan.estimated_logical_total / 2**20:.0f} MiB (estimate); "
        f"est. peak {plan.memory_estimated_peak_mib} MiB; free disk {plan.disk_free_bytes / 2**30:.1f} GiB")

    if dry_run:
        return {"dry_run": True, "plan": plan.as_dict(), "bcf_scan": {
            "lines": stats.lines, "pixels": stats.pixels, "maximum_count": stats.maximum,
            "total_counts": stats.total_counts, "largest_line_record_bytes": stats.largest_line_bytes},
            "grid_checks": layout.checks, "scan_seconds": stats.seconds}

    now = datetime.now(timezone.utc)
    identity.update(created_utc=now.isoformat(timespec="seconds"),
                    history=f"{now.isoformat(timespec='seconds')} microxrf_to_netcdf {__version__} convert "
                            f"{info.filename} + {rtx_path.name} -> {out.name}; parameters in "
                            "/metadata/provenance_json")
    partial = partial_path(out)
    writer: UnifiedWriter | None = None
    try:
        out.parent.mkdir(parents=True, exist_ok=True)
        writer = UnifiedWriter(partial, info, scan, layout, plan, config, identity)
        writer.create()
        write_started = time.perf_counter()
        if _fault:
            _fault("after_create", 0)
        stream_rtx_planes(rtx_path, writer.rtx_sink)
        stream_header_planes(info, writer.bcf_header_sink)
        if not writer.all_planes_written():
            raise ValidationError(f"the write passes over the RTX and the BCF header did not deliver every planned "
                                  f"plane; missing {writer.missing_planes()}")
        registration = _register(writer, layout, info, config, log)
        write_stats = BCFScan()
        for number, band in enumerate(iter_bands(source, config.band_lines, decode_dtype, config.max_line_bytes)):
            if _fault:
                _fault("band", number)
            writer.write_band(band)
            _accumulate(write_stats, band)
        if not _same_stats(stats, write_stats):
            raise ValidationError(f"the BCF statistics of the {decode_dtype} write pass differ from the exact uint32 "
                                  "scan pass (maximum, total, non-zero values or per-channel sums); a count may have "
                                  "wrapped")
        if not writer.all_lines_written():
            raise ValidationError(f"only {len(writer.lines_written)} of {info.height} lines were written")
        writer.write_sum_from_pixels(write_stats.channel_sums)
        actual_layout = writer.actual_layout()
        write_seconds = round(time.perf_counter() - write_started, 2)
        writer.close()
        log(f"written in {write_seconds} s")
        if _fault:
            _fault("after_write", 0)

        checks: list[dict[str, Any]] = list(layout.checks)
        checks.append(_check("streamed_pixel_sums_do_not_exceed_bcf_header_sum_spectrum", True,
                             f"header total {int(header_sum.sum())}, pixel total {stats.total_counts}, difference "
                             f"{int(header_sum.sum()) - stats.total_counts} counts "
                             f"({100 * (int(header_sum.sum()) - stats.total_counts) / max(int(header_sum.sum()), 1):.3f} %); "
                             "equality is NOT expected (cause unverified)"))
        checks.append(_check("write_pass_decode_equals_exact_uint32_scan_decode", True,
                             f"{decode_dtype} write pass vs uint32 scan pass: maximum {stats.maximum}, total "
                             f"{stats.total_counts}, non-zero values {stats.nonzero_values} and all per-channel sums "
                             "identical, so no count wrapped"))
        checks.append(_check("all_lines_and_planes_written", True,
                             f"{len(writer.lines_written)} lines; {len(writer.planes_written)} planes"))
        checks.extend(registration["checks"])
        if config.verify_readback:
            checks.extend(verify_readback(partial, source, layout, plan, config, write_stats.channel_sums,
                                          decode_dtype, log))
        failed = [c for c in checks if c["status"] == "FAILED"]
        if failed:
            raise ValidationError("validation failed: " + "; ".join(f"{c['check']}: {c['detail']}" for c in failed))
        checks.extend(check_lazy_structure(partial, plan.counts_dtype, len(layout.element_planes)))
        failed = [c for c in checks if c["status"] == "FAILED"]
        if failed:
            raise ValidationError("validation failed: " + "; ".join(f"{c['check']}: {c['detail']}" for c in failed))

        for path in inputs:
            if _stamp(path) != before[path]:
                raise ConversionError(f"input {path} changed during the conversion")
        elapsed = round(time.perf_counter() - started, 2)
        provenance = {
            "software": _software_versions(), "parameters": config.as_dict(),
            "inputs": {"bcf": {"filename": info.filename, "size_bytes": info.size_bytes, "sha256": identity["bcf_sha256"]},
                       "rtx": {"filename": rtx_path.name, "size_bytes": scan.size_bytes, "sha256": identity["rtx_sha256"],
                               "decompressed_payload_sha256": scan.payload_sha256}},
            "preflight_estimates": plan.as_dict(), "actual_layout": actual_layout,
            "decoding": {"scan_pass": f"uint32 ({scan_band_lines(config)} lines per band; exact)",
                         "write_and_readback_passes": f"{decode_dtype} ({config.band_lines} lines per band)",
                         "widened_bands_in_write_pass": write_stats.widened_bands},
            "bcf_scan": {"lines": stats.lines, "pixels": stats.pixels, "maximum_count": stats.maximum,
                         "total_counts": stats.total_counts, "nonzero_values": stats.nonzero_values,
                         "largest_line_record_bytes": stats.largest_line_bytes, "widened_bands": stats.widened_bands},
            "measured": {"conversion_seconds_until_finalization": elapsed, "write_pass_seconds": write_seconds,
                         "peak_working_set_mib_until_finalization": peak_working_set_mib()},
            "registration": registration["results"],
            "note": "the output file's own SHA-256 cannot be stored inside it; see the sidecar report",
        }
        _finalize(partial, provenance, checks)
        digest = file_sha256(partial)
        if out.exists() and not config.overwrite:
            raise PreflightError(f"{out} appeared during the conversion; refusing to overwrite")
        os.replace(partial, out)
        with netCDF4.Dataset(out) as final:
            if final.conversion_status != "complete":
                raise ConversionError("the finalized file is not marked complete")
    except BaseException as error:
        if writer is not None:
            try:
                writer.close()
            except Exception:  # noqa: BLE001
                pass
        try:
            partial.unlink(missing_ok=True)
        except OSError:
            pass
        if isinstance(error, (MicroXRFToNetCDFError, KeyboardInterrupt, SystemExit)):
            raise
        raise ConversionError(f"conversion failed and the partial file was removed: "
                              f"{type(error).__name__}: {error}") from error

    report = {
        "output": str(out), "output_size_bytes": out.stat().st_size, "output_sha256": digest,
        "conversion_seconds": round(time.perf_counter() - started, 2),
        "peak_working_set_mib": peak_working_set_mib(),
        "inputs": provenance["inputs"], "counts_dtype": plan.counts_dtype,
        "estimated_logical_bytes": plan.estimated_logical_total, "actual_layout": actual_layout,
        "compression": {"zlib": config.complevel > 0, "complevel": config.complevel, "shuffle": config.shuffle},
        "validation": checks, "plan": plan.as_dict(), "software": provenance["software"],
    }
    sidecar = out.with_name(out.name + ".report.json")
    sidecar.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    out.with_name(out.name + ".sha256").write_text(f"{digest}  {out.name}\n", encoding="ascii")
    return report


def _register(writer: UnifiedWriter, layout: Layout, info, config: ConversionConfig, log: Log) -> dict[str, Any]:
    """Compare each mosaic instance's Map footprint with the acquisition grid (records the outcome in the file)."""
    checks: list[dict[str, Any]] = []
    results: dict[str, Any] = {}
    if not config.check_registration:
        return {"checks": checks, "results": results}
    read_gray = writer.mosaic_gray_reader()
    for key, (group, rect, summary) in writer.footprints().items():
        image_index = key.replace(":", "_")
        if not summary["size_matches"]:
            reason = "Map rectangle size does not match the acquisition extent"
            writer.note_footprint_without_registration(group, summary, reason)
            results[str(image_index)] = {"status": "not_verified", "reason": reason, "footprint": summary}
            continue
        crop_bytes = 4 * summary["rect_width_px_inclusive"] * summary["rect_height_px_inclusive"] * 3
        if crop_bytes > config.registration_max_crop_bytes:
            reason = f"footprint crop ({crop_bytes / 2**20:.0f} MiB) exceeds registration_max_crop_bytes"
            writer.note_footprint_without_registration(group, summary, reason)
            results[str(image_index)] = {"status": "not_verified", "reason": reason, "footprint": summary}
            continue
        result = correlate_footprint(read_gray, rect, info.video, (0, 0))
        result["footprint"] = summary
        writer.record_registration(group, result)
        results[str(image_index)] = result
        r = result.get("correlation_by_orientation", {}).get("identity")
        # An unverified relationship is recorded in the file, not treated as a conversion failure.
        checks.append(_check(f"mosaic_instance_{image_index}_registration_status_recorded", True,
                             f"status {result['status']}; identity r={r}; best shift {result.get('best_shift_mosaic_px_dx_dy')}"))
        log(f"registration of mosaic instance {image_index}: {result['status']} (identity r={r})")
    return {"checks": checks, "results": results}


def _finalize(partial: Path, provenance: dict[str, Any], checks: list[dict[str, Any]]) -> None:
    """Reopen for append; store provenance and validation; mark the file complete; flush and close."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        with netCDF4.Dataset(partial, "a") as ds:
            metadata = ds["metadata"]
            metadata["provenance_json"][...] = json.dumps(provenance, indent=1, ensure_ascii=False, default=str)
            metadata["validation_json"][...] = json.dumps(checks, indent=1, ensure_ascii=False)
            ds.conversion_status = "complete"
