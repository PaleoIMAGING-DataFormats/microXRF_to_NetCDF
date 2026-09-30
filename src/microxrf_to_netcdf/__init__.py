"""
Authors: Andre L. Belem (https://github.com/andrebelem) and F.R.I.D.A.Y. (https://observatoriooceanografico.org/people/friday-bot)

microXRF to NetCDF: converts a paired Bruker BCF + RTX acquisition into one NetCDF-4 file, sequentially and
with bounded memory. See NETCDF_SCHEMA.md for the output schema and FINDINGS.md for the verified facts.
"""

import os

# The converter does no linear algebra. NumPy's OpenBLAS would nevertheless reserve (commit) hundreds of MiB per
# process on a many-core machine at import time (measured: about 685 MiB versus 17 MiB with one thread, FINDINGS.md
# section 9), which matters under commit-limited environments and job memory limits. Only effective if NumPy has not
# been imported yet; a value set by the user is respected.
for _variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_variable, "1")

__version__ = "0.1.0"
SCHEMA_VERSION = "1.0.0"
