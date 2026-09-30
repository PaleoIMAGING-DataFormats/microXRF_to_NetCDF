# CLAUDE.md — microXRF to NetCDF

Permanent instructions for coding agents working in this repository.

This file defines the persistent project context, architecture, conventions, scientific constraints, and development rules for **microXRF to NetCDF**.

It is intended to bootstrap new AI-assisted development sessions. Although the filename `CLAUDE.md` follows the Claude Code convention, its contents are written to be usable as project context with other coding agents as well.

When starting a new coding-agent session, ask the agent to read this file before modifying the repository.


## Project identity and mission

- Project name: **microXRF to NetCDF**. Repository directory: **`microXRF_to_NetCDF`**. Python package and import name: **`microxrf_to_netcdf`**. Conda environment name: **`microxrf_to_netcdf`**.
- **Project scope:** micro-XRF → interoperable NetCDF-4 scientific hypercube. The identity is the scientific modality and the target data model, not a vendor.
- **Current implemented source:** Bruker BCF + RTX (one acquisition pair, validated on one real dataset).
- **Current scientific capabilities:** EDS spectral cube, processed elemental maps, video/images, optical mosaics, acquisition metadata and provenance.
- **Not yet justified:** claiming generic support for all micro-XRF vendors or formats. Do not claim or document support for a format or instrument that has not been tested. Keep vendor-specific reading in vendor-named modules (`bcf.py`, `rtx.py`) so that a future source adapter can be added without weakening the validated Bruker implementation; do not create a speculative adapter abstraction before a second source has been inspected (principle I).
- Bruker, BCF, RTX, EDS, XRF and RosettaSciIO name the vendor, the source formats and the reader used. They are not renamed.
- Mission: read Bruker BCF and RTX acquisition files, investigate and preserve their scientific content, and convert the combined information into NetCDF-4 datasets suitable for lazy access and chunked processing.
- The software must handle acquisitions substantially larger than available RAM.

## Scientific interpretation and validation of source fields

The primary goal of **microXRF to NetCDF** is to extract, preserve, and reorganize information from micro-XRF source files into an open and interoperable NetCDF-4 representation.

### Critical scientific constraint

The scientific meaning, units, calibration, and interpretation of all fields found in the proprietary BCF and RTX formats have **not been independently validated**.

A field being successfully extracted from a BCF or RTX file does **not** mean that its physical meaning, unit, calibration, or scientific interpretation has been established.

The current converter therefore prioritizes **faithful preservation over interpretation**.

When the meaning of a source field is uncertain:

1. Preserve the original value without applying undocumented scientific transformations.
2. Preserve the original source field name whenever practical. If the field is mapped to a different name in the NetCDF schema, the mapping must be explicitly documented.
3. Do not invent, infer, normalize, or assign physical units unless those units have been independently verified.
4. Do not infer physical meaning from field names, data types, numerical ranges, neighboring fields, GUI labels, or apparent correlations with other variables.
5. Do not convert numerical values into alternative physical quantities unless the conversion rule and units have been independently established.
6. Do not silently rename ambiguous fields to scientifically meaningful names.
7. Do not treat values displayed by Bruker software or found in RTX display annotations as authoritative metadata unless their relationship to the corresponding source field has been verified.
8. Derived quantities present in the source files must remain distinguishable from raw measurements. Their processing method must not be inferred when it is undocumented.
9. If a field must be exposed through the NetCDF schema before its meaning is established, document its interpretation status explicitly as unknown, unverified, or source-derived where appropriate.
10. New scientific interpretations must be supported by evidence such as vendor documentation, instrument documentation, controlled experiments, independent reference data, or reproducible comparison with trusted software.

### Extraction validation is not scientific validation

Keep the following distinction explicit throughout the project:

- **Extraction validation** establishes that information was decoded and transferred correctly from the source file.
- **Structural validation** establishes that dimensions, arrays, relationships, and storage structures are represented consistently.
- **Scientific validation** establishes what a quantity physically represents, its units, calibration, processing history, and whether it can be used for a particular scientific interpretation.

Successful extraction or structural validation must never be described as scientific validation.

For example, byte-identical arrays, agreement with an independent decoder, correlation between an RTX elemental map and an EDS-derived spectral window, or successful round-trip storage in NetCDF demonstrate important aspects of data integrity. They do not, by themselves, establish the physical meaning, units, calibration, or processing algorithm of the source quantity.

### NetCDF metadata

NetCDF metadata must reflect the current level of knowledge. Do not add metadata merely to make the dataset appear more complete or more CF-like.

In particular, do not assign or infer physical units, `standard_name`, calibration equations, detector corrections, quantification methods, elemental-map processing methods, coordinate reference interpretations, instrument parameters, or acquisition parameters unless they are supported by verified evidence.

It is preferable to preserve an incomplete but accurate representation of the source metadata than to create a complete-looking but scientifically unsupported representation.

### Documentation of new findings

When the meaning of a previously uncertain field is established:

1. Record the evidence and reasoning in `FINDINGS.md`.
2. Distinguish direct evidence from interpretation.
3. Update `NETCDF_SCHEMA.md` if the interpretation affects the public data model.
4. Add or update tests where the interpretation can be tested programmatically.
5. Only then update the converter or NetCDF metadata.

Unresolved fields and interpretations must remain documented as unresolved until adequate evidence is available.

## Unified acquisition architecture (established 2026-09-29)

These principles govern the converter (`microxrf_to_netcdf`, schema in [NETCDF_SCHEMA.md](NETCDF_SCHEMA.md)):

1. **One paired BCF/RTX acquisition produces one unified NetCDF-4 scientific hypercube.** Different acquisitions remain separate physical files.
2. **Conversion is sequential and memory-bounded for both inputs.** The BCF spectrum stream, the BCF header and the RTX payload are each decoded incrementally (a few scan lines, or one image plane, at a time); nothing decoded is kept beyond the band or plane being written. Passes over an input may be repeated; materializing a whole cube or a whole payload may not.
3. **All scientific content keeps its source identity, calibration, units and provenance.** Raw values and derived values have different names. A unit the source does not state is not assigned (record `units_status`). Content that cannot be interpreted is preserved verbatim (residual XML, raw texts), never dropped, and anything a reader does not return (for example the seven images in the BCF header, of which RosettaSciIO returns one) is carried or reported.
4. **The NetCDF-4 output may contain several groups and several spatial grids.** Data share dimensions only after their grid correspondence is validated. Equal pixels do not make two source images interchangeable: each source instance keeps its own timestamp, annotations and provenance.
5. **No arbitrary maximum output size is imposed.** A logical size larger than RAM is expected. Only an estimated peak memory above free RAM, or missing or insufficient disk space, stops a conversion.
6. **Preflight resource checks and safe failure handling are mandatory.** Before the output is created: inspect dimensions and dtypes, estimate every variable's logical size (marked as an estimate, distinct from measured sizes), report the chunks, check free disk space of the destination filesystem against the uncompressed worst case plus a margin, estimate peak memory. Unknown or insufficient resources fail before writing, or require an explicit override. Output goes to a partial file and is finalized by an atomic rename after verification; failure or interruption removes the partial file; the final name never holds an unverified product.
7. **The file must be suitable for lazy access to spatial and spectral subsets** (xarray and Dask, HDF5 chunking chosen from measured access patterns, documented trade-offs).
8. Preserve all existing findings and completed tasks in [FINDINGS.md](FINDINGS.md) and [TODO.md](TODO.md); correct a superseded statement visibly instead of deleting it.

## Language policy

Everything is written in English: source code, comments, docstrings, tests, configuration, documentation, CLI output, error messages, commit messages and pull requests. Use the exact names "microXRF to NetCDF", `microXRF_to_NetCDF` (repository) and `microxrf_to_netcdf` (package).

## Authorship

Every new Python source file must contain this exact line (as in [inspect_bruker.py](tools/diagnostics/inspect_bruker.py)):

```
Authors: Andre L. Belem (https://github.com/andrebelem) and F.R.I.D.A.Y. (https://observatoriooceanografico.org/people/friday-bot)
```

## Architectural principles

A. RosettaSciIO is the current verified BCF reader. Preserve lazy Dask access whenever possible.
B. The RTX format is not yet fully understood. Do not invent its schema or silently discard its contents.
C. Do not assume BCF and RTX are interchangeable, or that one contains all information present in the other.
D. Design the conversion around an explicit scientific data model, not merely around copying arrays between formats.
E. Preserve spatial and spectral calibration, instrument information, acquisition metadata, units, masks and provenance.
F. NetCDF-4 is the intended output. It must support lazy reopening, suitable chunking and partial access. The converter must use bounded memory and must never materialize a complete spectral cube.
G. Original BCF, RTX and image files are immutable inputs. Conversion always creates new outputs.
H. Keep microXRF to NetCDF independent from Hyper Fusion and HSI Toolbox. Future integration uses documented interfaces, not shared internal code.
I. Do not create speculative readers, conversion modules or storage abstractions before the relevant data structures have been verified.
J. Keep dependencies explicit and tested. Never install packages into the Conda base environment.
K. Every scientific claim and format assumption must rest on actual inspection, documented specifications or reproducible tests.
L. Do not commit, push, rename repositories or modify other projects without explicit authorization.

## Scientific integrity and provenance

- Record what is observed, not what is expected. Separate confirmed facts, derived values, interpretations and open questions, as in [FINDINGS.md](FINDINGS.md).
- Never report a value that was not read from a file or computed by code. Do not infer units; verify them (for example, the beam energy value `50` has unverified units).
- Outputs must record provenance: source file names, sizes, software versions (microXRF to NetCDF, RosettaSciIO, and others), and conversion parameters.

## Repository organization

Current contents:

- [microxrf_to_netcdf/](src/microxrf_to_netcdf/) — the converter package: `bcf.py` (line framer, band decoder, header images), `rtx_payload.py` (Base64 + zlib + expat payload decoder), `rtx.py` (incremental RTX and TRT-document scanner, residual XML), `model.py` (data-model rules), `registration.py` (mosaic footprint check), `preflight.py` (estimates, disk, memory, destination), `config.py`, `writer.py` (incremental NetCDF-4 writer, schema 1.0.0), `convert.py` (orchestration, atomic finalization), `validate.py` (independent checks), `memory.py`, `errors.py`, `cli.py` / `__main__.py` (`python -m microxrf_to_netcdf plan|convert|validate`). `microxrf_to_netcdf.rtx` builds on the verified payload decoder `rtx_payload.py`; do not edit that decoder without re-running its tests. `pyproject.toml` (src layout, dependencies read from `requirements.txt`) makes the package installable.
- [tools/diagnostics/](tools/diagnostics/) — read-only diagnostics: [inspect_bruker.py](tools/diagnostics/inspect_bruker.py) (inventory and report of the files in `data/`), [inspect_rtx.py](tools/diagnostics/inspect_rtx.py) (prints a report of an RTX using the package decoder). [tools/benchmarks/](tools/benchmarks/) — [probe_bcf_streaming.py](tools/benchmarks/probe_bcf_streaming.py) (read-only prototype framer and streaming probe, an independent reference of the tests) and [benchmark_conversion.py](tools/benchmarks/benchmark_conversion.py) (memory and chunk benchmark that converts into a scratch directory you name). Run them from the repository root; they must stay read-only. No script belongs in the repository root.
- [data/](data/) (original acquisitions and images, immutable), `output/` (generated NetCDF-4 files, ignored by git; never inside `data/`), [tests/](tests/) with [pytest.ini](pytest.ini) (markers: `large` reads `data/`; `slow` is synthetic but takes tens of seconds), [reference/](reference/) (saved diagnostic, benchmark, conversion and validation output), [requirements.txt](requirements.txt), [environment.yml](environment.yml), [.gitignore](.gitignore).
- Documents: [FINDINGS.md](FINDINGS.md), [NETCDF_SCHEMA.md](NETCDF_SCHEMA.md), [TODO.md](TODO.md), [README.md](README.md).

Create further packages or directories only when a task in [TODO.md](TODO.md) needs them, and update this section when that happens. Keep `data/` free of generated output.

## BCF and RTX reader responsibilities

- BCF reader (`microxrf_to_netcdf.bcf`): the SFS container and the compiled decoder of the installed RosettaSciIO (0.14.x only, checked at run time), fed a few whole scan lines at a time by the line framer; never RosettaSciIO's single-chunk Dask array. It refuses an unknown `flag`, inconsistent record lengths, `pixel_x >= width`, oversized records, a short stream and trailing bytes. It also reads every image in the BCF header (RosettaSciIO returns only the video) and keeps the header XML without the image text as residual XML. Report anything a reader does not return rather than dropping it silently.
- RTX reader (`microxrf_to_netcdf.rtx`): two sequential passes (scan, then plane-by-plane), Base64 + zlib + expat, one plane in memory, residual XML for everything that is not pixel data (see [FINDINGS.md](FINDINGS.md) section 8). Unsupported layouts (an image that is neither the map image nor a 3-plane 8-bit mosaic, `ItemSize` other than 1, 2 or 4) are refused with a clear error, not dropped. Do not assume the RTX element maps equal raw energy-window sums of the BCF; record their processing and units as unknown.
- A video plane of the RTX that differs from the BCF video is a validation failure, never a silent merge. Mosaic instances with equal pixels are stored once but keep their own metadata.
- Readers only read. They never write to input locations.

## Lazy reading and chunked processing

- Never call `.compute()`, `np.asarray` or similar on a full EDS cube (240 x 1800 x 4096 uint8 in the current dataset, and larger in general).
- Process by chunks whose size is configurable and documented; state the memory bound of every pipeline.
- Diagnostic and test code must also read bounded amounts of data (headers, slices, sampled pixels).

## Memory-bounded NetCDF-4 writing

- Write incrementally, band by band and plane by plane, to a new file. Never overwrite an input, and never write next to the inputs or inside a directory named `data`.
- Chunking and compression are configurable, with documented defaults (see [NETCDF_SCHEMA.md](NETCDF_SCHEMA.md) section 9); the band height must be a multiple of the chunk extent along `y`.
- A written file must reopen lazily and support partial reads with correct values.
- Measure peak RAM; do not assume it. The preflight memory figure is an estimate; the report carries the measured peak.
- Write to `.<name>.partial`, verify, then `os.replace` to the final name. Remove the partial file on any failure, including `KeyboardInterrupt`. The output's own hash goes in a sidecar file, not inside it.

## Dataset and coordinate conventions

Fixed in [NETCDF_SCHEMA.md](NETCDF_SCHEMA.md) (schema 1.0.0). Constraints that must not be relaxed without a schema version change: axis order `(y, x, energy)`; coordinate values `index x raw calibration`; the derived `energy` coordinate is separate from the raw calibration variable; counts stay unsigned integers with **no `_FillValue`** (a valid zero must never become `NaN`); no units are assigned to element maps, `PixelTimes`, overview images, beam energy or elevation angle; the group layout is `/acquisition`, `/mosaic`, `/overview`, `/metadata`. Bump `SCHEMA_VERSION` when a variable is renamed, removed or changes meaning.

## Preservation of original acquisitions

Files in `data/` are read-only inputs. Never modify, move, rename, re-save or delete them. Conversion writes to a new location.

## Testing and validation

- Add automated tests for every behavior and every documented finding.
- Tests requiring the large files must skip cleanly when the files are absent and must not load them entirely.
- Compare converted output against source values on representative samples.
- Mark a task done in [TODO.md](TODO.md) only when a test or recorded observation supports it. Report failing tests as failing.

## Dependency management

- Use the `microxrf_to_netcdf` Conda environment (Python 3.12; install the package with `pip install -e . --no-deps`). Never install into the Conda base environment.
- Keep `requirements.txt` complete and tested in a fresh environment. Add a dependency only when it is used, and record the reason.
- Record verified versions (for example RosettaSciIO 0.14.0); the current table is in section 6 of [FINDINGS.md](FINDINGS.md).
- Run tests with `python -m pytest` inside the environment (`-m "not large and not slow"` for the fast subset). The Python of that environment is `C:\Users\abelem\miniforge3\envs\microxrf_to_netcdf\python.exe` when `python` is not on the shell PATH.
- netCDF4 and xarray are dependencies (netCDF4 1.7.4, xarray 2026.7.0 verified); `requirements.txt` is the single source and `environment.yml` repeats their constraints (`tests/test_environment_files.py` keeps them identical). Do not disable SSL verification to make conda-forge work; the pip route is the verified one.
- RosettaSciIO's lazy EDX array is a single 1.77 GB chunk (FINDINGS.md section 1): do not slice or compute it. The verified bounded-memory route is sequential decoding of a few scan lines at a time with the unmodified compiled decoder (FINDINGS.md sections 7 and 9); Dask rechunking does not provide streaming input. Never set `_FillValue = 0` on count variables (xarray then returns `float32` with `NaN`).
- **Counts precision.** The compiled decoder adds in place in the chosen dtype, so `uint16` wraps silently (70000 becomes 4464, a small value that no "near the ceiling" test flags). The converter scans in `uint32` (exact), picks the output dtype from the exact maximum, decodes the write and read-back passes in `uint16` only if the maximum fits, and requires identical per-channel sums, total, non-zero count and maximum between the passes. Never choose a dtype from RosettaSciIO's `estimate_map_depth`.
- The BCF header sum spectrum is only an upper bound of the pixel sums on the real file (it exceeds them by 0.44 %, cause unknown): never require equality.
- `microxrf_to_netcdf` sets `OPENBLAS_NUM_THREADS`, `OMP_NUM_THREADS` and `MKL_NUM_THREADS` to 1 when they are unset: on this many-core machine `import numpy` otherwise commits about 685 MiB. Import `microxrf_to_netcdf` before NumPy in entry points that must stay within a memory limit.

## Documentation maintenance

Update [FINDINGS.md](FINDINGS.md), [TODO.md](TODO.md), [README.md](README.md) and this file in the same change as the work that affects them. If documents disagree, resolve the contradiction and report it.

## Git safety

- This directory is not currently a Git repository. Do not run `git init` without authorization.
- Do not commit, push, create branches or pull requests, rename repositories or touch other projects (including Hyper Fusion and HSI Toolbox) without explicit authorization.
- Never use destructive Git operations, and never bypass hooks.
- When authorized, commit messages and pull requests are written in English.
