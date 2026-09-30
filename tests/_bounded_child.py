"""
Authors: Andre L. Belem (https://github.com/andrebelem) and F.R.I.D.A.Y. (https://observatoriooceanografico.org/people/friday-bot)

Child process of tests/test_bounded_memory.py (Windows). It places ITSELF in a Job Object with a hard limit on
committed process memory BEFORE importing numpy, netCDF4 or microxrf_to_netcdf, then converts a synthetic acquisition
whose logical size exceeds that limit. Prints one JSON line. Not a test module.

Usage: python _bounded_child.py <limit_mib> <work_dir> [height width channels peaks]
"""

import ctypes
import json
import sys
from ctypes import wintypes


def limit_process_memory(limit_bytes: int) -> None:
    class IoCounters(ctypes.Structure):
        _fields_ = [(name, ctypes.c_ulonglong) for name in (
            "ReadOperationCount", "WriteOperationCount", "OtherOperationCount", "ReadTransferCount",
            "WriteTransferCount", "OtherTransferCount")]

    class Basic(ctypes.Structure):
        _fields_ = [("PerProcessUserTimeLimit", ctypes.c_longlong), ("PerJobUserTimeLimit", ctypes.c_longlong),
                    ("LimitFlags", wintypes.DWORD), ("MinimumWorkingSetSize", ctypes.c_size_t),
                    ("MaximumWorkingSetSize", ctypes.c_size_t), ("ActiveProcessLimit", wintypes.DWORD),
                    ("Affinity", ctypes.c_size_t), ("PriorityClass", wintypes.DWORD),
                    ("SchedulingClass", wintypes.DWORD)]

    class Extended(ctypes.Structure):
        _fields_ = [("BasicLimitInformation", Basic), ("IoInfo", IoCounters), ("ProcessMemoryLimit", ctypes.c_size_t),
                    ("JobMemoryLimit", ctypes.c_size_t), ("PeakProcessMemoryUsed", ctypes.c_size_t),
                    ("PeakJobMemoryUsed", ctypes.c_size_t)]

    kernel = ctypes.WinDLL("kernel32", use_last_error=True)
    kernel.CreateJobObjectW.restype = ctypes.c_void_p
    kernel.GetCurrentProcess.restype = ctypes.c_void_p
    kernel.SetInformationJobObject.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_void_p, wintypes.DWORD]
    kernel.AssignProcessToJobObject.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
    job = kernel.CreateJobObjectW(None, None)
    info = Extended()
    info.BasicLimitInformation.LimitFlags = 0x100  # JOB_OBJECT_LIMIT_PROCESS_MEMORY
    info.ProcessMemoryLimit = limit_bytes
    if not kernel.SetInformationJobObject(job, 9, ctypes.byref(info), ctypes.sizeof(info)):
        raise OSError(ctypes.get_last_error(), "SetInformationJobObject failed")
    if not kernel.AssignProcessToJobObject(job, kernel.GetCurrentProcess()):
        raise OSError(ctypes.get_last_error(), "AssignProcessToJobObject failed")


def main() -> None:
    import os
    for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):   # before numpy is imported
        os.environ.setdefault(name, "1")
    limit_mib, work = int(sys.argv[1]), sys.argv[2]
    height, width, channels, peaks = (int(v) for v in sys.argv[3:7]) if len(sys.argv) >= 7 else (48, 1000, 16384, 20)
    limit_process_memory(limit_mib * 2**20)

    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
    import numpy as np

    import synthetic_data as sd
    from microxrf_to_netcdf.config import ConversionConfig
    from microxrf_to_netcdf.convert import convert
    from microxrf_to_netcdf.memory import peak_working_set_mib

    seed = 7
    video = (np.arange(height * width, dtype=np.uint16).reshape(height, width) * 3)
    info = sd.make_info(height, width, channels, video=video, sum_spectrum=np.full(channels, 2**40, np.uint64))
    source = sd.SyntheticBCF(info, lambda: sd.sparse_blocks(height, width, channels, peaks, seed))
    work = Path(work)
    work.mkdir(parents=True, exist_ok=True)
    rtx = work / "pair.rtx"
    sd.write_pair_rtx(rtx, video)
    config = ConversionConfig(band_lines=1, counts_chunks=(1, 50, 4096), map_chunks=(1, 24, 250),
                              video_chunks=(24, 250), mosaic_chunks=(1, 64, 512), complevel=1)
    report = convert(source, rtx, work / "out" / "big.nc", config)
    print(json.dumps({"peak_working_set_mib": peak_working_set_mib(), "limit_mib": limit_mib,
                      "logical_counts_bytes": height * width * channels, "output": report["output"],
                      "output_size_bytes": report["output_size_bytes"], "counts_dtype": report["counts_dtype"],
                      "seed": seed, "peaks": peaks}))


if __name__ == "__main__":
    main()
