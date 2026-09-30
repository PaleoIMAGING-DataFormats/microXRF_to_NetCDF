"""
Authors: Andre L. Belem (https://github.com/andrebelem) and F.R.I.D.A.Y. (https://observatoriooceanografico.org/people/friday-bot)

Bounded-memory test on a synthetic dataset whose LOGICAL size exceeds the imposed process-memory limit.

The child process (tests/_bounded_child.py) puts itself in a Windows Job Object with a hard limit on committed
memory before importing anything heavy, then converts a synthetic acquisition streamed on the fly (one scan line
in memory at a time): 48 x 1000 pixels x 16384 channels = 786 MB of uint8 counts against a limit smaller than
that. A control run with a much smaller limit must fail, which shows that the limit is really enforced.
Windows only (skipped elsewhere). Marked ``slow``: it takes tens of seconds.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import netCDF4
import numpy as np
import pytest

import synthetic_data as sd

CHILD = Path(__file__).with_name("_bounded_child.py")
LIMIT_MIB = 500
HEIGHT, WIDTH, CHANNELS, PEAKS, SEED = 48, 1000, 16384, 20, 7

pytestmark = [pytest.mark.slow,
              pytest.mark.skipif(sys.platform != "win32", reason="uses a Windows Job Object memory limit")]


def run_child(limit_mib: int, work: Path):
    return subprocess.run([sys.executable, str(CHILD), str(limit_mib), str(work)], capture_output=True, text=True,
                          timeout=900)


def test_control_run_a_tiny_limit_really_stops_the_process(tmp_path):
    result = run_child(40, tmp_path / "control")
    assert result.returncode != 0
    assert not (tmp_path / "control" / "out" / "big.nc").exists()


def test_conversion_whose_logical_size_exceeds_the_memory_limit_completes_correctly(tmp_path):
    logical = HEIGHT * WIDTH * CHANNELS
    assert logical > LIMIT_MIB * 2**20, "the test is only meaningful if the dataset is larger than the limit"
    result = run_child(LIMIT_MIB, tmp_path / "work")
    assert result.returncode == 0, result.stderr[-3000:]
    info = json.loads(result.stdout.strip().splitlines()[-1])
    assert info["logical_counts_bytes"] == logical > info["limit_mib"] * 2**20
    assert info["counts_dtype"] == "uint8" and info["peak_working_set_mib"] < LIMIT_MIB
    out = Path(info["output"])
    with netCDF4.Dataset(out) as ds:
        assert ds.conversion_status == "complete"
        counts = ds["acquisition"]["counts"]
        assert counts.shape == (HEIGHT, WIDTH, CHANNELS) and counts.dtype == np.uint8
        total = 0
        for y in (0, 17, HEIGHT - 1):                       # sampled lines against the regenerated ground truth
            expected = sd.sparse_expected_line(y, WIDTH, CHANNELS, PEAKS, SEED)
            assert np.array_equal(np.asarray(counts[y]), expected.astype(np.uint8)), y
            total += int(expected.sum())
        assert total > 0
        sums = np.asarray(ds["acquisition"]["sum_spectrum_from_pixels"][:])
        assert int(sums.sum()) == sum(int(sd.sparse_line(y, WIDTH, CHANNELS, PEAKS, SEED)[1].sum()) for y in range(HEIGHT))
    assert out.stat().st_size < logical // 10                # compressed output, far below the logical size
