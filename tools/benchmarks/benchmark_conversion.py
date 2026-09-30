"""
Authors: Andre L. Belem (https://github.com/andrebelem) and F.R.I.D.A.Y. (https://observatoriooceanografico.org/people/friday-bot)

Read-only diagnostic that measures (a) the peak memory of the converter for several band sizes and (b) the size,
write time and read timings of several HDF5 chunk shapes. It converts the real acquisition into a scratch
directory that you name (never ``data/``); each run is a fresh child process so the peak working set is its own.

    python tools/benchmarks/benchmark_conversion.py memory --scratch <dir> [--bands 1 4 8 20]
    python tools/benchmarks/benchmark_conversion.py chunks --scratch <dir>

The peak is the Windows peak working set of the child (whole conversion: scan, write, read-back), from the
report the converter writes. Read timings run on a warm operating-system cache, so they compare decompression
and chunk-index cost, not disk speed; they are indicative and machine-specific.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import netCDF4
import numpy as np

BCF = Path("data/GRF17A_9-29cm_slab3_Elemental_map.bcf")
RTX = Path("data/GRF17A_9-29cm_slab3_Elemental_map.rtx")
ROOT = Path(__file__).resolve().parents[2]
CHILD_ENV = {**os.environ, "PYTHONPATH": os.pathsep.join(filter(None, [str(ROOT / "src"), os.environ.get("PYTHONPATH")]))}
CHUNK_CANDIDATES = {
    "A_default_4x60x4096": (4, 60, 4096),
    "B_spectrum_1x16x4096": (1, 16, 4096),
    "C_energy_split_8x64x256": (8, 64, 256),
    "D_window_2x30x4096": (2, 30, 4096),
}


def run_convert(out: Path, band_lines: int, chunks: tuple[int, int, int], verify: bool) -> dict:
    command = [sys.executable, "-m", "microxrf_to_netcdf", "convert", "--bcf", str(BCF), "--rtx", str(RTX), "--out", str(out),
               "--band-lines", str(band_lines), "--chunks", *map(str, chunks), "--overwrite", "--quiet"]
    if not verify:
        command.append("--no-verify")
    started = time.perf_counter()
    subprocess.run(command, check=True, capture_output=True, text=True, env=CHILD_ENV)
    report = json.loads(out.with_name(out.name + ".report.json").read_text(encoding="utf-8"))
    report["wall_seconds"] = round(time.perf_counter() - started, 1)
    return report


def read_timings(path: Path) -> dict:
    """Time each access on a freshly opened file (cold HDF5 chunk cache; warm operating-system cache)."""
    rng = np.random.default_rng(0)
    pixels = [(int(rng.integers(0, 240)), int(rng.integers(0, 1800))) for _ in range(20)]
    accesses = {
        "pixel_spectrum_ms_mean_of_20": (lambda c: [c[y, x, :] for y, x in pixels], 20),
        "window_16x16_all_energy_ms": (lambda c: c[100:116, 900:916, :], 1),
        "window_8x100x400ch_ms": (lambda c: c[50:58, 1000:1100, 0:400], 1),
        "one_energy_channel_image_ms": (lambda c: c[:, :, 468], 1),
    }
    results = {}
    for name, (function, divisor) in accesses.items():
        times = []
        for _ in range(3):
            with netCDF4.Dataset(path) as ds:
                counts = ds["acquisition"]["counts"]
                counts.set_auto_maskandscale(False)
                start = time.perf_counter()
                function(counts)
                times.append(time.perf_counter() - start)
        results[name] = round(min(times) * 1000 / divisor, 2)
    return results


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    memory = sub.add_parser("memory")
    memory.add_argument("--scratch", type=Path, required=True)
    memory.add_argument("--bands", type=int, nargs="+", default=[1, 4, 8, 20])
    chunks = sub.add_parser("chunks")
    chunks.add_argument("--scratch", type=Path, required=True)
    args = parser.parse_args()
    args.scratch.mkdir(parents=True, exist_ok=True)
    if args.command == "memory":
        for band in args.bands:
            # one-line bands need chunks of one line along y
            chunk = (1, 60, 4096) if band == 1 else (4, 60, 4096)
            report = run_convert(args.scratch / f"mem_band{band}.nc", band, chunk, verify=True)
            print(f"band_lines={band} chunks={chunk}: peak working set {report['peak_working_set_mib']} MiB, "
                  f"{report['conversion_seconds']} s, estimate {report['plan']['memory_mib']['estimated_peak']} MiB "
                  f"(estimate = ESTIMATE, peak = MEASURED)", flush=True)
    else:
        for name, chunk in CHUNK_CANDIDATES.items():
            band = max(4, chunk[0])
            out = args.scratch / f"chunks_{name}.nc"
            report = run_convert(out, band, chunk, verify=False)
            timings = read_timings(out)
            print(json.dumps({"candidate": name, "chunks": chunk, "band_lines": band,
                              "size_MiB": round(report["output_size_bytes"] / 2**20, 1),
                              "conversion_seconds": report["conversion_seconds"], **timings}), flush=True)


if __name__ == "__main__":
    main()
