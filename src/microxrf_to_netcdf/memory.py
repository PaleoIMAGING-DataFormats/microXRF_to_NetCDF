"""
Authors: Andre L. Belem (https://github.com/andrebelem) and F.R.I.D.A.Y. (https://observatoriooceanografico.org/people/friday-bot)

Process-memory helpers. ``peak_working_set_mib`` uses the Windows API (as in the diagnostics) and falls
back to ``resource.getrusage`` elsewhere; both are the peak resident set of THIS process, including the
interpreter and imported libraries. ``available_memory_bytes`` reports free physical memory or None.
"""

from __future__ import annotations

import sys


def _windows_counters():
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
    return counters


def peak_working_set_mib() -> float | None:
    """Peak resident memory of this process in MiB, or None when it cannot be measured."""
    try:
        if sys.platform == "win32":
            return round(_windows_counters().peak / 2**20, 1)
        import resource
        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return round(peak / (2**20 if sys.platform == "darwin" else 2**10), 1)
    except (AttributeError, OSError, ImportError):
        return None


def current_working_set_mib() -> float | None:
    try:
        if sys.platform == "win32":
            return round(_windows_counters().current / 2**20, 1)
    except (AttributeError, OSError):
        pass
    return None


def available_memory_bytes() -> int | None:
    """Free physical memory in bytes, or None when unknown (then no memory preflight decision is made)."""
    try:
        if sys.platform == "win32":
            import ctypes

            class Status(ctypes.Structure):
                _fields_ = [("length", ctypes.c_ulong), ("load", ctypes.c_ulong),
                            ("total_phys", ctypes.c_ulonglong), ("avail_phys", ctypes.c_ulonglong),
                            ("total_page", ctypes.c_ulonglong), ("avail_page", ctypes.c_ulonglong),
                            ("total_virtual", ctypes.c_ulonglong), ("avail_virtual", ctypes.c_ulonglong),
                            ("avail_ext", ctypes.c_ulonglong)]

            status = Status()
            status.length = ctypes.sizeof(status)
            ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(status))
            return int(status.avail_phys)
        with open("/proc/meminfo", encoding="ascii") as meminfo:
            for line in meminfo:
                if line.startswith("MemAvailable:"):
                    return int(line.split()[1]) * 1024
    except (AttributeError, OSError, ValueError):
        pass
    return None
