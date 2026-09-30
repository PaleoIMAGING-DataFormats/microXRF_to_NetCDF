"""
Authors: Andre L. Belem (https://github.com/andrebelem) and F.R.I.D.A.Y. (https://observatoriooceanografico.org/people/friday-bot)

Read-only feasibility probe for sequential (streaming) BCF decoding and incremental NetCDF-4
writing. It is a diagnostic, NOT the converter and NOT a reader: the line framer below is a
throw-away prototype whose only purpose is to produce evidence recorded in FINDINGS.md section 7.

It reuses the installed RosettaSciIO 0.14.0 (never modified): the SFS container reader for the
raw byte stream and the compiled ``unbcf_fast.parse_to_numpy`` decoder for the spectra, fed with a
synthetic stream that holds only a few scan lines at a time.

Subcommands (run from the repository root inside the project environment):

    python tools/benchmarks/probe_bcf_streaming.py structure
    python tools/benchmarks/probe_bcf_streaming.py memory [--band-lines 1]
    python tools/benchmarks/probe_bcf_streaming.py verify [--band-lines 1]
    python tools/benchmarks/probe_bcf_streaming.py netcdf --out <new .nc path outside data/> [--band-lines 4]
                                         [--chunks 4 60 4096] [--complevel 4]

``netcdf`` needs netCDF4 and xarray, which are NOT in requirements.txt yet.
Peak memory is the Windows peak working set of this process (ctypes, no extra dependency).
"""

from __future__ import annotations

import argparse
import collections
import os
import struct
import sys
import time
from pathlib import Path

import numpy as np

DEFAULT_BCF = Path("data/GRF17A_9-29cm_slab3_Elemental_map.bcf")
STREAM_PATH = "EDSDatabase/SpectrumData0"
HEADER_BYTES = 0x1A0          # the decoder seeks here to find the first line (unbcf_fast.pyx)
PIXEL_HEADER = struct.Struct("<IHHIHHHI")  # pixel_x, chan1, chan2, 4 skipped bytes, flag, size1, n_pulses, size2
LOW_CHANNELS = 400            # channels compared per pixel against the public decoder path
BIN = 8                       # spatial binning factor of the whole-cube reference (240 and 1800 divisible by 8)


def peak_memory_mib() -> tuple[float, float]:
    """Return (current, peak) working set of this process in MiB (Windows only)."""
    import ctypes
    from ctypes import wintypes

    class Counters(ctypes.Structure):
        _fields_ = [("cb", wintypes.DWORD), ("faults", wintypes.DWORD),
                    ("peak", ctypes.c_size_t), ("current", ctypes.c_size_t),
                    ("_pad", ctypes.c_size_t * 4),
                    ("pagefile", ctypes.c_size_t), ("peak_pagefile", ctypes.c_size_t)]

    kernel, psapi = ctypes.windll.kernel32, ctypes.windll.psapi
    kernel.GetCurrentProcess.restype = ctypes.c_void_p
    psapi.GetProcessMemoryInfo.argtypes = [ctypes.c_void_p, ctypes.POINTER(Counters), wintypes.DWORD]
    counters = Counters()
    counters.cb = ctypes.sizeof(counters)
    psapi.GetProcessMemoryInfo(kernel.GetCurrentProcess(), ctypes.byref(counters), counters.cb)
    return round(counters.current / 2**20, 1), round(counters.peak / 2**20, 1)


def open_stream(bcf: Path):
    """Return the SFS item holding the spectrum stream (RosettaSciIO, read-only)."""
    from rsciio.bruker._api import SFS_reader
    sfs = SFS_reader(str(bcf))
    return sfs, sfs.get_file(STREAM_PATH)


def frame_lines(blocks):
    """Prototype framer. Yield the header bytes, then (line_index, raw_line_bytes) per scan line.

    It only follows the record framing that ``bin_to_numpy`` follows (line word, 22-byte pixel
    header, payload length); it never decodes spectra. Memory: one line record plus one block.
    """
    buffer, position, source = bytearray(), 0, iter(blocks)

    def need(count):
        while len(buffer) - position < count:
            buffer.extend(next(source))

    need(HEADER_BYTES)
    header = bytes(buffer[:HEADER_BYTES])
    position = HEADER_BYTES
    height, width = struct.unpack_from("<II", header)
    yield header
    for line in range(height):
        del buffer[:position]
        position = 0
        need(4)
        (pixels,) = struct.unpack_from("<I", buffer, position)
        position += 4
        for _ in range(pixels):
            need(PIXEL_HEADER.size)
            pixel_x, _c1, _c2, _skip, flag, _size1, pulses, size2 = PIXEL_HEADER.unpack_from(buffer, position)
            position += PIXEL_HEADER.size
            if pixel_x >= width:  # unbcf_fast is compiled with boundscheck(False): refuse, never write OOB
                raise ValueError(f"pixel_x {pixel_x} >= width {width} on line {line}")
            if flag in (0, 1):
                length = size2
            else:
                length = size2 - 4 + ((4 + 2 * pulses) if pulses > 0 else 4)
            need(length)
            position += length
        yield line, bytes(buffer[:position]), pixels
    yield "trailing", len(buffer) - position


class _OneBlockFile:
    """Duck-typed stand-in for an SFS item: parse_to_numpy only calls get_iter_and_properties()."""

    def __init__(self, block: bytes):
        self.block = block

    def get_iter_and_properties(self):
        return iter([self.block]), len(self.block), 1


def stream_bands(item, band_lines: int, channels: int, dtype):
    """Yield (first_line, array[band, width, channels]) decoded by the unmodified Cython decoder."""
    from rsciio.bruker import unbcf_fast
    framer = frame_lines(item.get_iter_and_properties()[0])
    header = next(framer)
    height, width = struct.unpack_from("<II", header)
    pending, first = [], 0
    for record in framer:
        if record[0] == "trailing":
            if record[1]:
                raise ValueError(f"{record[1]} unexpected trailing bytes after the last line")
            break
        pending.append(record[1])
        if len(pending) == band_lines or record[0] == height - 1:
            patched = bytearray(header)
            struct.pack_into("<I", patched, 0, len(pending))  # decoder loops over this many lines
            block = bytes(patched) + b"".join(pending)
            band = unbcf_fast.parse_to_numpy(_OneBlockFile(block), (len(pending), width, channels), dtype)
            yield first, band
            first += len(pending)
            pending = []


def cmd_structure(args) -> None:
    _, item = open_stream(args.bcf)
    start = time.time()
    counts, lines, pixels, longest, max_x, sizes = collections.Counter(), 0, 0, 0, 0, []
    framer = frame_lines(item.get_iter_and_properties()[0])
    header = next(framer)
    print("stream bytes:", item.size, "| SFS block bytes:", item.get_iter_and_properties()[1],
          "| SFS compression:", item.sfs.compression)
    print("header height, width:", struct.unpack_from("<II", header))
    for record in framer:
        if record[0] == "trailing":
            print("trailing bytes after last line:", record[1])
            break
        lines += 1
        pixels += record[2]
        longest = max(longest, len(record[1]))
        counts[f"pixels_in_line={record[2]}"] += 1
    print("lines:", lines, "| pixels listed:", pixels, "| longest line record bytes:", longest,
          "| counters:", dict(counts))
    print("framing time s:", round(time.time() - start, 1), "| working set / peak MiB:", peak_memory_mib())


def flag_survey(item) -> dict:
    """Count (flag, n_pulses>0) combinations to see which decoder branches the file exercises."""
    survey, buffer, source, position = collections.Counter(), bytearray(), item.get_iter_and_properties()[0], 0

    def need(count):
        while len(buffer) - position < count:
            buffer.extend(next(source))

    need(HEADER_BYTES)
    height, _ = struct.unpack_from("<II", buffer)
    position = HEADER_BYTES
    for _ in range(height):
        del buffer[:position]
        position = 0
        need(4)
        (pixels,) = struct.unpack_from("<I", buffer, position)
        position += 4
        for _ in range(pixels):
            need(22)
            _x, _c1, _c2, _s, flag, _s1, pulses, size2 = PIXEL_HEADER.unpack_from(buffer, position)
            position += 22
            survey[(flag, pulses > 0)] += 1
            length = size2 if flag in (0, 1) else size2 - 4 + ((4 + 2 * pulses) if pulses > 0 else 4)
            need(length)
            position += length
    return {f"flag={f},n_pulses>0={p}": n for (f, p), n in survey.items()}


def reference_binned(bcf: Path) -> np.ndarray:
    """Whole-cube reference through the PUBLIC path, bounded: downsample=BIN sums BIN x BIN pixels."""
    from rsciio.bruker._api import BCF_reader
    return BCF_reader(str(bcf)).parse_hypermap(index=0, downsample=BIN, lazy=False).astype(np.int64)


def reference_low_channels(bcf: Path, channels: int) -> np.ndarray:
    """Per-pixel reference through the unmodified decoder on the REAL stream, truncated to `channels`."""
    from rsciio.bruker import unbcf_fast
    _, item = open_stream(bcf)
    shape = struct.unpack_from("<II", next(frame_lines(item.get_iter_and_properties()[0])))
    return unbcf_fast.parse_to_numpy(item, (*shape, channels), np.uint8)


def cmd_memory(args) -> None:
    """Stream every band (uint16, all 4096 channels) and discard it: isolates the decoder's peak memory."""
    _, item = open_stream(args.bcf)
    before, start, top, bands = peak_memory_mib(), time.time(), 0, 0
    for _, band in stream_bands(item, args.band_lines, 4096, np.uint16):
        top, bands = max(top, int(band.max())), bands + 1
    print(f"band_lines={args.band_lines}: {bands} bands of {args.band_lines * 1800 * 4096 * 2 / 2**20:.1f} MiB (uint16), "
          f"{time.time() - start:.1f} s, max count {top}, working set / peak MiB before {before} after {peak_memory_mib()}")


def cmd_verify(args) -> None:
    _, item = open_stream(args.bcf)
    print("branches exercised by this file:", flag_survey(item))
    start = time.time()
    binned = np.zeros((240 // BIN, 1800 // BIN, 4096), np.int64)
    low = np.zeros((240, 1800, LOW_CHANNELS), np.uint16)
    top, total, band_mib = 0, 0, 0.0
    for first, band in stream_bands(item, args.band_lines, 4096, np.uint16):
        top = max(top, int(band.max()))
        total += int(band.sum(dtype=np.int64))
        band_mib = max(band_mib, band.nbytes / 2**20)
        low[first:first + band.shape[0]] = band[:, :, :LOW_CHANNELS]
        for k in range(band.shape[0]):
            binned[(first + k) // BIN] += band[k].reshape(1800 // BIN, BIN, 4096).sum(1, dtype=np.int64)
    print(f"streamed: band_lines={args.band_lines}, band {band_mib:.1f} MiB (uint16), "
          f"{time.time() - start:.1f} s, max per-pixel channel count {top}, total counts {total}")
    print("working set / peak MiB after streaming (includes the two small references above):", peak_memory_mib())
    print(f"binned x{BIN}, all 4096 channels, equals public-path downsample={BIN}:",
          np.array_equal(binned, reference_binned(args.bcf)))
    print(f"per-pixel, first {LOW_CHANNELS} channels, equals public decoder on the real stream:",
          np.array_equal(low, reference_low_channels(args.bcf, LOW_CHANNELS).astype(np.uint16)))


def cmd_netcdf(args) -> None:
    out = args.out.resolve()
    if "data" in [p.name for p in out.parents] or out.exists():
        sys.exit("refusing: --out must be a NEW file outside data/")
    import netCDF4
    import xarray as xr
    _, item = open_stream(args.bcf)
    start = time.time()
    with netCDF4.Dataset(out, "w", format="NETCDF4") as dataset:
        for name, size in (("y", 240), ("x", 1800), ("energy", 4096)):
            dataset.createDimension(name, size)
        # fill_value=False: NO _FillValue attribute. With _FillValue=0 xarray masks every zero count to NaN.
        counts = dataset.createVariable("counts", "u1", ("y", "x", "energy"), chunksizes=tuple(args.chunks),
                                        zlib=args.complevel > 0, complevel=max(args.complevel, 1),
                                        shuffle=args.complevel > 0, fill_value=False)
        for first, band in stream_bands(item, args.band_lines, 4096, np.uint16):
            if int(band.max()) > 255:
                raise OverflowError("counts exceed uint8; choose a wider dtype")
            counts[first:first + band.shape[0]] = band.astype(np.uint8)
    print(f"written {out.name}: {time.time() - start:.1f} s, {out.stat().st_size / 2**20:.1f} MiB, "
          f"working set / peak MiB {peak_memory_mib()}")
    reopened = xr.open_dataset(out, chunks={})["counts"]
    print("reopened lazily:", type(reopened.data).__name__, reopened.dtype, reopened.shape, "chunksize", reopened.data.chunksize)
    low = reference_low_channels(args.bcf, LOW_CHANNELS)
    spectrum = reopened.isel(y=100, x=900).values
    window = reopened.isel(y=slice(50, 58), x=slice(1000, 1100), energy=slice(0, LOW_CHANNELS)).values
    image = reopened.isel(energy=100).values
    print("pixel spectrum (first channels) equal:", np.array_equal(spectrum[:LOW_CHANNELS], low[100, 900]))
    print("window equal:", np.array_equal(window, low[50:58, 1000:1100]))
    print("energy-channel image equal:", np.array_equal(image, low[:, :, 100]))
    binned = reopened.astype("int64").coarsen(y=BIN, x=BIN).sum().values
    print(f"whole-file binned x{BIN} equals public-path downsample={BIN}:", np.array_equal(binned, reference_binned(args.bcf)))
    print("working set / peak MiB after reopen and checks:", peak_memory_mib())


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bcf", type=Path, default=DEFAULT_BCF)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("structure").set_defaults(func=cmd_structure)
    memory = sub.add_parser("memory")
    memory.add_argument("--band-lines", type=int, default=1)
    memory.set_defaults(func=cmd_memory)
    verify = sub.add_parser("verify")
    verify.add_argument("--band-lines", type=int, default=1)
    verify.set_defaults(func=cmd_verify)
    netcdf = sub.add_parser("netcdf")
    netcdf.add_argument("--out", type=Path, required=True)
    netcdf.add_argument("--band-lines", type=int, default=4)
    netcdf.add_argument("--chunks", type=int, nargs=3, default=(4, 60, 4096))
    netcdf.add_argument("--complevel", type=int, default=4)
    netcdf.set_defaults(func=cmd_netcdf)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
