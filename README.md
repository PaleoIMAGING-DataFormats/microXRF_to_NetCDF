# microXRF to NetCDF

`microXRF to NetCDF` converts micro-XRF acquisition data into a single, self-contained, chunked and compressed NetCDF-4 scientific hypercube. Conversion is sequential and uses bounded memory, allowing large spectral-imaging acquisitions to be converted without loading the complete cube into RAM.

The current implementation supports paired **Bruker BCF + RTX** acquisitions and preserves EDS spectra, elemental maps, optical/video images, mosaics, acquisition metadata, and provenance. The resulting NetCDF-4 can be opened lazily with xarray and Dask.

> **Status:** under active development. The current implementation has been validated against the Münster GRF17 Bruker BCF + RTX acquisition. Other acquisitions and format variants still require validation.

## Installation

Python 3.12 is recommended.

```bash
pip install -r requirements.txt
pip install -e . --no-deps
```

## Reference dataset

Development and validation currently use the **GRF17 interlaboratory micro-XRF comparison dataset**.

Reference acquisition:

**GRF_17A_9-29cm_slab 3_Münster_µXRF**

Download the dataset from:

https://drive.proton.me/urls/C2594KQJ90#HQoTFkgIrQuw

Place the required files in `data/`:

```text
data/
├── GRF17A_9-29cm_slab3_Elemental_map.bcf
├── GRF17A_9-29cm_slab3_Elemental_map.rtx
├── CaFe.png
└── Video Mosaic.png
```

Large scientific source files and generated NetCDF files are intentionally excluded from Git.

## Usage

Inspect the acquisition and conversion plan:

```bash
python -m microxrf_to_netcdf plan   --bcf data/GRF17A_9-29cm_slab3_Elemental_map.bcf   --rtx data/GRF17A_9-29cm_slab3_Elemental_map.rtx   --out output/GRF17A_9-29cm_slab3.nc
```

Convert to NetCDF-4:

```bash
python -m microxrf_to_netcdf convert   --bcf data/GRF17A_9-29cm_slab3_Elemental_map.bcf   --rtx data/GRF17A_9-29cm_slab3_Elemental_map.rtx   --out output/GRF17A_9-29cm_slab3.nc
```

Validate the converted dataset:

```bash
python -m microxrf_to_netcdf validate   --nc output/GRF17A_9-29cm_slab3.nc   --bcf data/GRF17A_9-29cm_slab3_Elemental_map.bcf   --rtx data/GRF17A_9-29cm_slab3_Elemental_map.rtx   --data-dir data   --strict
```

## Lazy access

```python
import xarray as xr

path = "output/GRF17A_9-29cm_slab3.nc"

acq = xr.open_dataset(path, group="acquisition", chunks={})
counts = acq["counts"]

spectrum = counts.isel(y=100, x=900)
window = counts.isel(
    y=slice(50, 58),
    x=slice(1000, 1100),
    energy=slice(0, 400),
)

ca_map = acq["element_maps"].sel(element="Ca")
mosaic = xr.open_dataset(path, group="mosaic", chunks={})["pixels"]
```

The arrays remain lazy until values are explicitly computed or loaded.

## Documentation

- [`NETCDF_SCHEMA.md`](NETCDF_SCHEMA.md) — NetCDF-4 data model and conventions.
- [`FINDINGS.md`](FINDINGS.md) — validation results, measurements, and technical findings.
- [`TODO.md`](TODO.md) — remaining scientific and engineering work.
- [`tools/diagnostics/`](tools/diagnostics/) — format inspection utilities.
- [`tools/benchmarks/`](tools/benchmarks/) — streaming, memory, and chunking benchmarks.

## Tests

```bash
python -m pytest
```

Tests requiring the original acquisition are skipped when the reference files are not available under `data/`.


## Authors

Andre L. Belem (https://github.com/andrebelem) and F.R.I.D.A.Y. (https://observatoriooceanografico.org/people/friday-bot)

### Scientific validation disclaimer

The current implementation focuses on the faithful extraction, preservation, and conversion of information found in the source BCF and RTX files.

The scientific meaning, units, calibration, and interpretation of all fields extracted from these proprietary formats have **not yet been independently validated** against Bruker documentation, instrument specifications, or external reference measurements.

Where the meaning of a source field is uncertain, microXRF to NetCDF preserves the source information as directly as possible, retaining the original naming or an explicitly documented mapping into the NetCDF schema. Inclusion of a field in the NetCDF file must therefore **not** be interpreted as confirmation that its physical meaning, units, calibration, or scientific validity have been established.

This distinction is particularly important for derived or insufficiently documented quantities from the BCF and RTX files. See [`FINDINGS.md`](FINDINGS.md) for currently verified observations and unresolved scientific questions.

### Vibe coding disclaimer

This repository was developed using a vibe-coding workflow, with extensive use of AI coding agents for implementation, refactoring, testing, and documentation under human direction and scientific supervision.

The code and scientific outputs should not be assumed to be correct solely because they were generated or validated by automated agents. Reproducibility, independent validation, and review against the original scientific data remain essential.

[`CLAUDE.md`](CLAUDE.md) contains the project architecture, conventions, constraints, and development instructions used during AI-assisted development. It can also be used as project context when starting a new session with the coding agent of your choice.
