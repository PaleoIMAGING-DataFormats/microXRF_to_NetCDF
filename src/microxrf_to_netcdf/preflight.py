"""
Authors: Andre L. Belem (https://github.com/andrebelem) and F.R.I.D.A.Y. (https://observatoriooceanografico.org/people/friday-bot)

Resource preflight. Runs BEFORE the output file is created and answers, from the inspected inputs:

1. the logical (uncompressed) size of every major output variable and of the whole dataset;
2. the proposed chunk configuration;
3. the free disk space of the destination filesystem against a conservative worst case;
4. the expected peak process memory from the decoding buffers, output chunks and header parsing.

Sizes here are ESTIMATES. Compression ratios depend on the data, so the compressed size can only be
measured after writing; the worst case used for the disk check assumes no compression at all. Nothing
here limits the output size: a logical size larger than RAM is expected and fine. Only an estimated PEAK
MEMORY above free RAM, or unknown/insufficient free disk space, stops the conversion (unless overridden).
"""

from __future__ import annotations

import math
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from .bcf import BCFInfo
from .config import ConversionConfig
from .errors import PreflightError
from .memory import available_memory_bytes
from .model import Layout
from .rtx import RTXScan

MIB = 2**20
PYTHON_BASELINE_MIB = 110.0        # interpreter + numpy + netCDF4 + RosettaSciIO imported (measured order of magnitude)
HEADER_PEAK_FACTOR = 5.4           # RosettaSciIO header parse peak / HeaderData bytes (182 MiB for 34 MB, FINDINGS 7.1)
DECODE_FACTOR = 2.4                # peak / band size while decoding (FINDINGS 7.3)
RETENTION_MIB = 64.0               # HDF5 chunk cache and allocator retention seen after the band loops (FINDINGS 9)
ESTIMATE_MARGIN = 1.25             # applied to the summed memory estimate
NETCDF_OVERHEAD_FRACTION = 0.01    # allowance for HDF5 metadata and chunk index in the worst case
NETCDF_OVERHEAD_BYTES = 16 * MIB


@dataclass
class BCFScan:
    """Statistics gathered by the first sequential pass over the BCF stream (no data is kept)."""

    lines: int = 0
    pixels: int = 0
    maximum: int = 0
    total_counts: int = 0
    nonzero_values: int = 0
    channel_sums: np.ndarray | None = None
    largest_line_bytes: int = 0
    widened_bands: int = 0
    seconds: float = 0.0


@dataclass
class Plan:
    counts_dtype: str
    counts_dtype_reason: str
    estimated_logical_bytes: dict[str, int]
    estimated_logical_total: int
    chunks: dict[str, tuple[int, ...]]
    band_lines: int
    decode_dtype: str
    disk_free_bytes: int
    disk_worst_case_bytes: int
    disk_required_bytes: int
    disk_ok: bool
    memory_estimate_mib: dict[str, float]
    memory_estimated_peak_mib: float
    memory_available_mib: float | None
    memory_ok: bool | None
    warnings: list[str] = field(default_factory=list)

    def as_dict(self) -> dict[str, Any]:
        return {
            "status": "ESTIMATES (not measured)",
            "counts_dtype": self.counts_dtype, "counts_dtype_reason": self.counts_dtype_reason,
            "estimated_logical_bytes": self.estimated_logical_bytes,
            "estimated_logical_total_bytes": self.estimated_logical_total,
            "chunks": {k: list(v) for k, v in self.chunks.items()},
            "band_lines": self.band_lines, "decode_dtype": self.decode_dtype,
            "disk": {"free_bytes": self.disk_free_bytes, "worst_case_output_bytes": self.disk_worst_case_bytes,
                     "required_free_bytes": self.disk_required_bytes, "sufficient": self.disk_ok},
            "memory_mib": {**self.memory_estimate_mib, "estimated_peak": self.memory_estimated_peak_mib,
                           "available": self.memory_available_mib, "sufficient": self.memory_ok},
            "warnings": self.warnings,
        }


def choose_counts_dtype(maximum: int, requested: str) -> tuple[str, str]:
    """Smallest unsigned integer type holding ``maximum``; an explicit request that is too small is refused."""
    limits = {"uint8": 2**8 - 1, "uint16": 2**16 - 1, "uint32": 2**32 - 1}
    if requested == "auto":
        for name, limit in limits.items():
            if maximum <= limit:
                return name, f"auto: observed maximum count {maximum} fits {name} (max {limit})"
        raise PreflightError(f"observed maximum count {maximum} exceeds uint32")
    if maximum > limits[requested]:
        raise PreflightError(f"requested counts dtype {requested} cannot hold the observed maximum count "
                             f"{maximum}; refusing to truncate or wrap")
    return requested, f"requested {requested}; observed maximum count {maximum} fits (max {limits[requested]})"


def clip(chunks: tuple[int, ...], shape: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(max(1, min(int(c), int(s))) for c, s in zip(chunks, shape))


def nearest_existing_parent(path: Path) -> Path:
    path = path.resolve()
    while not path.exists():
        if path.parent == path:
            break
        path = path.parent
    return path


def make_plan(info: BCFInfo, scan: BCFScan, layout: Layout, rtx: RTXScan, config: ConversionConfig,
              out_path: Path) -> Plan:
    """Estimate sizes, disk and memory. Raises PreflightError when the conversion must not start."""
    config.validate()
    grid, mosaics = layout.grid, [m.image for m in layout.mosaics]
    element_count, unique_mosaic_sets = len(layout.element_planes), layout.unique_mosaic_sets
    header_bytes = info.header_bytes
    header_text = len(info.header_scan.residual_xml.encode("utf-8")) if info.header_scan else 0
    warnings: list[str] = []
    dtype_name, reason = choose_counts_dtype(scan.maximum, config.counts_dtype)
    itemsize = np.dtype(dtype_name).itemsize
    y, x, e = info.height, info.width, info.channels
    counts_chunks = clip(config.counts_chunks, (y, x, e))
    if config.band_lines % counts_chunks[0] != 0 and config.band_lines < y:
        raise PreflightError(f"band_lines ({config.band_lines}) must be a multiple of the counts chunk extent "
                             f"along y ({counts_chunks[0]}) so that every chunk is written once")

    sizes = {
        "acquisition/counts": y * x * e * itemsize,
        "acquisition/video": y * x * 2,
        "acquisition/element_maps": element_count * y * x * 2,
        "acquisition/sum_spectra_and_coordinates": e * (8 + 8 + 8 + 4) + (y + x) * 8 + y * 4,
        "acquisition/header_images_on_grid": sum(a.image.width * a.image.height * a.image.itemsize
                                                 for a in layout.aux_images),
        "mosaic/pixels": unique_mosaic_sets * mosaics[0].width * mosaics[0].height * 3 if mosaics else 0,
        "overview/pixels": sum(o.image.width * o.image.height * 3 for o in layout.overviews),
        "metadata/text": len(rtx.residual_xml.encode("utf-8")) + header_text + 1_000_000,
    }
    total = sum(sizes.values())
    chunks = {
        "counts": counts_chunks,
        "element_maps": clip(config.map_chunks, (element_count, y, x)),
        "video": clip(config.video_chunks, (y, x)),
        "mosaic": clip(config.mosaic_chunks, (3, mosaics[0].height, mosaics[0].width)) if mosaics else (),
    }

    # Disk: worst case = no compression benefit at all. Compressed size is only measured afterwards.
    worst = int(total * (1 + NETCDF_OVERHEAD_FRACTION)) + NETCDF_OVERHEAD_BYTES
    required = worst + max(int(worst * config.disk_margin_fraction), config.disk_margin_bytes)
    parent = nearest_existing_parent(out_path.parent)
    try:
        free = shutil.disk_usage(parent).free
    except OSError as error:
        raise PreflightError(f"cannot determine free disk space of {parent}: {error}; nothing was written "
                             "(use allow_low_disk to override)") from error
    disk_ok = free >= required
    if not disk_ok:
        message = (f"insufficient free disk space on {parent}: {free / MIB:.0f} MiB free, "
                   f"{required / MIB:.0f} MiB required (worst case {worst / MIB:.0f} MiB uncompressed + margin)")
        if not config.allow_low_disk:
            raise PreflightError(message + "; nothing was written (allow_low_disk overrides at your own risk)")
        warnings.append(message + " -- OVERRIDDEN by allow_low_disk")

    # Memory
    decode_itemsize = np.dtype(config.decode_dtype).itemsize
    band_bytes = config.band_lines * x * e * decode_itemsize
    header_images = info.header_scan.images if info.header_scan else []
    max_plane = int(max(im.width * im.height * im.itemsize for im in [grid, *mosaics, *header_images]))
    largest_line = max(scan.largest_line_bytes, 1)
    header_peak = (HEADER_PEAK_FACTOR * header_bytes / MIB) if header_bytes else 0.0
    decode_peak = (DECODE_FACTOR * band_bytes + 3 * config.band_lines * largest_line) / MIB
    write_buffers = float(config.band_lines * x * e * itemsize
                          + math.prod(counts_chunks) * itemsize * math.ceil(x / counts_chunks[1]) * 2) / MIB
    rtx_peak = (3 * max_plane + 3 * (len(rtx.residual_xml) + header_text)) / MIB
    out_band = config.band_lines * x * e * itemsize / MIB
    parts = {
        "python_baseline": PYTHON_BASELINE_MIB,
        "header_parse_transient": round(header_peak, 1),
        "decode_band_and_raw_records": round(decode_peak, 1),
        "output_conversion_and_chunk_cache": round(write_buffers, 1),
        "hdf5_and_allocator_retention": RETENTION_MIB,
        "rtx_plane_and_residual": round(rtx_peak, 1),
    }
    # The stages do not overlap: header parsing, the write pass, the read-back pass and the RTX planes each set a peak.
    peak = float(ESTIMATE_MARGIN * max(PYTHON_BASELINE_MIB + header_peak,
                                       PYTHON_BASELINE_MIB + decode_peak + write_buffers + RETENTION_MIB,
                                       PYTHON_BASELINE_MIB + decode_peak + out_band + RETENTION_MIB,
                                       PYTHON_BASELINE_MIB + rtx_peak))
    available = available_memory_bytes()
    available_mib = None if available is None else round(available / MIB, 1)
    memory_ok = None if available is None else peak <= config.memory_fraction * available / MIB
    if memory_ok is False:
        message = (f"estimated peak memory {peak:.0f} MiB exceeds {config.memory_fraction:.0%} of the "
                   f"{available_mib:.0f} MiB currently free; reduce band_lines")
        if not config.allow_low_memory:
            raise PreflightError(message)
        warnings.append(message + " -- OVERRIDDEN by allow_low_memory")
    if available is None:
        warnings.append("free physical memory could not be determined; no memory decision was made")
    return Plan(
        counts_dtype=dtype_name, counts_dtype_reason=reason, estimated_logical_bytes=sizes,
        estimated_logical_total=total, chunks=chunks, band_lines=config.band_lines,
        decode_dtype=config.decode_dtype, disk_free_bytes=free, disk_worst_case_bytes=worst,
        disk_required_bytes=required, disk_ok=disk_ok, memory_estimate_mib=parts,
        memory_estimated_peak_mib=round(peak, 1), memory_available_mib=available_mib, memory_ok=memory_ok,
        warnings=warnings)


def check_scan_memory(info: BCFInfo, config: ConversionConfig, lines: int) -> float:
    """Estimated peak (MiB) of the scan pass, which decodes in uint32; raises before it starts if it cannot fit."""
    header_peak = (HEADER_PEAK_FACTOR * info.header_bytes / MIB) if info.header_bytes else 0.0
    band_bytes = lines * info.width * info.channels * 4
    decode_peak = DECODE_FACTOR * band_bytes / MIB
    peak = float(ESTIMATE_MARGIN * max(PYTHON_BASELINE_MIB + header_peak, PYTHON_BASELINE_MIB + decode_peak))
    available = available_memory_bytes()
    if available is not None and peak > config.memory_fraction * available / MIB and not config.allow_low_memory:
        raise PreflightError(f"estimated peak memory of the scan pass {peak:.0f} MiB exceeds "
                             f"{config.memory_fraction:.0%} of the {available / MIB:.0f} MiB currently free; "
                             "reduce band_lines")
    return round(peak, 1)


def check_destination(out_path: Path, input_paths: list[Path], overwrite: bool) -> None:
    """Refuse an unsafe destination before anything is read in earnest or written."""
    out_path = out_path.resolve()
    for source in input_paths:
        if out_path == source.resolve():
            raise PreflightError("the output path is one of the input files")
    if out_path.suffix.lower() not in (".nc", ".nc4"):
        raise PreflightError(f"output path {out_path.name} should end in .nc or .nc4")
    for source in input_paths:
        source = source.resolve()
        if out_path.parent == source.parent:
            raise PreflightError(f"the output must not be written next to the original acquisition files "
                                 f"({source.parent}); choose another directory")
    if "data" in [part.lower() for part in out_path.parts[:-1]]:
        raise PreflightError("the output path lies inside a directory named 'data' (reserved for immutable inputs)")
    if out_path.exists() and not overwrite:
        raise PreflightError(f"{out_path} already exists; refusing to overwrite (pass overwrite=True / --overwrite)")
    if out_path.is_dir():
        raise PreflightError(f"{out_path} is a directory")
    partial = partial_path(out_path)
    if partial.exists():
        raise PreflightError(f"a partial file from an earlier run exists: {partial}; remove it after checking it "
                             "is not in use")


def partial_path(out_path: Path) -> Path:
    return out_path.with_name(f".{out_path.name}.partial")
