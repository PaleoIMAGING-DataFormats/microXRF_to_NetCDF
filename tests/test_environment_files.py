"""
Authors: Andre L. Belem (https://github.com/andrebelem) and F.R.I.D.A.Y. (https://observatoriooceanografico.org/people/friday-bot)

Fast tests that requirements.txt, environment.yml and the installed packages agree, and that the packages the
converter needs work together (netCDF4 writes what xarray and Dask read lazily).
"""

from __future__ import annotations

import re
import warnings
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent


def constraints(text: str) -> dict[str, str]:
    found = {}
    for line in text.splitlines():
        line = line.strip().lstrip("- ").strip()
        match = re.fullmatch(r"([A-Za-z0-9_.]+)\s*((?:[<>]=?)[^#]+)", line.split("#")[0].strip())
        if match:
            found[match.group(1).lower()] = match.group(2).replace(" ", "")
    return found


def as_tuple(text: str) -> tuple[int, ...]:
    return tuple(int(part) for part in re.findall(r"\d+", text)[:3])


def test_netcdf_and_xarray_are_listed_with_identical_constraints_in_both_files():
    requirements = constraints((ROOT / "requirements.txt").read_text(encoding="utf-8"))
    environment = constraints((ROOT / "environment.yml").read_text(encoding="utf-8"))
    for package in ("netcdf4", "xarray"):
        assert package in requirements, package
        assert environment[package] == requirements[package]


@pytest.mark.parametrize("package", ["rosettasciio", "dask", "numpy", "netCDF4", "xarray", "pillow", "pytest"])
def test_installed_versions_satisfy_the_requirements(package):
    requirements = constraints((ROOT / "requirements.txt").read_text(encoding="utf-8"))
    spec = requirements[package.lower()]
    lower = re.search(r">=([0-9.]+)", spec)
    upper = re.search(r"<([0-9.]+)", spec)
    installed = as_tuple(version(package))
    assert lower is None or installed >= as_tuple(lower.group(1)), (package, installed, spec)
    assert upper is None or installed < as_tuple(upper.group(1)), (package, installed, spec)


def test_the_stack_works_together_netcdf4_writes_and_xarray_dask_reads_lazily(tmp_path):
    import netCDF4
    import xarray as xr
    path = tmp_path / "stack.nc"
    with netCDF4.Dataset(path, "w", format="NETCDF4") as ds:
        for name, size in (("y", 5), ("x", 6), ("energy", 8)):
            ds.createDimension(name, size)
        counts = ds.createVariable("counts", "u1", ("y", "x", "energy"), chunksizes=(2, 3, 8), zlib=True,
                                   complevel=3, shuffle=True, fill_value=False)
        data = np.zeros((5, 6, 8), dtype=np.uint8)
        data[1, 2, 3] = 7
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            for y0 in range(0, 5, 2):
                counts[y0:y0 + 2] = data[y0:y0 + 2]
    with xr.open_dataset(path, chunks={}) as dataset:
        array = dataset["counts"]
        assert hasattr(array.data, "dask") and array.dtype == np.uint8 and array.data.chunksize == (2, 3, 8)
        assert np.array_equal(array.values, data)          # zeros are zeros, not NaN
        assert int(array.sum().values) == 7
