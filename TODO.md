# microXRF to NetCDF — Roadmap

Prioritized scientific development tasks. Each task has an acceptance criterion (AC).
A task is marked `[x]` only when supported by an actual test or recorded observation.
Verified facts live in [FINDINGS.md](FINDINGS.md); the output contract is [NETCDF_SCHEMA.md](NETCDF_SCHEMA.md).

## Architecture of the unified converter (established 2026-09-29)

- One paired BCF/RTX acquisition produces **one** unified NetCDF-4 scientific hypercube. Different acquisitions remain separate physical files.
- Conversion is **sequential and memory-bounded for both inputs** (BCF spectrum stream and header, RTX payload).
- All scientific content retains its **source identity, calibration, units and provenance**; unstated units are not assigned; uninterpretable content is preserved verbatim.
- The NetCDF-4 output may contain **multiple groups and spatial grids** (`/acquisition`, `/mosaic`, `/overview`, `/metadata`).
- **No arbitrary maximum output-file size** is imposed.
- **Preflight resource checks and safe failure handling are mandatory** (estimates versus measurements, disk space against the uncompressed worst case, memory estimate, temporary file, atomic finalization, cleanup on failure).
- The file must be **suitable for lazy access** to spatial and spectral subsets.

Status of the first complete converter (2026-09-29/30; every item below has a test or a recorded measurement, FINDINGS.md section 9):

- [x] Phase 1 architecture principles recorded here and in [CLAUDE.md](CLAUDE.md).
- [x] Phase 2 dependencies: netCDF4 1.7.4 and xarray 2026.7.0 in [requirements.txt](requirements.txt) and [environment.yml](environment.yml); `pip check` clean; suite passes; a fresh virtualenv built from `requirements.txt` alone reproduces the saved diagnostic outputs and passes the fast tests (`tests/test_environment_files.py`).
- [x] Phase 3 versioned schema 1.0.0 documented in [NETCDF_SCHEMA.md](NETCDF_SCHEMA.md) and enforced by a test that every written variable is documented.
- [x] Phase 4 bounded-memory BCF reader (`src/microxrf_to_netcdf/bcf.py`): framer plus the unmodified decoder, configurable band lines (default 4, measured peak 379 MiB for the converter), exact `uint32` scan before any `uint16` decode.
- [x] Phase 5 bounded-memory RTX reader (`src/microxrf_to_netcdf/rtx.py`): two sequential passes, one plane in memory, residual XML, video identity check, mosaics stored once with per-instance metadata.
- [x] Phase 6 preflight (`src/microxrf_to_netcdf/preflight.py`, `python -m microxrf_to_netcdf plan`).
- [x] Phase 7 incremental writer (`src/microxrf_to_netcdf/writer.py`): configurable chunks and zlib, no `_FillValue`, temporary file and atomic finalization.
- [x] Phase 8 end-to-end conversion of the real pair to `output/GRF17A_9-29cm_slab3.nc` (226,985,418 bytes, 51.3 s, 453.4 MiB peak; FINDINGS 9.7).
- [x] Phase 9 scientific validation (29 checks in the file; 21 independent checks by `python -m microxrf_to_netcdf validate --strict`; properties that cannot be validated are listed in FINDINGS 9.8).
- [x] Phase 10 tests and documentation (197 tests: 145 fast, 50 `large`, 2 `slow`).

Remaining work (none of these is done):

- [ ] R1. A **second real acquisition** (other pixel encodings `flag` 0/1, `n_of_pulses` > 0, zlib-compressed SFS container, `FileVersion` 2 with overview images, another RTX layout, another detector). Until then those branches are verified by synthetic streams only.
- [ ] R2. A dataset **larger than the physical RAM** of the machine converted end to end (only a process-memory cap on a synthetic 786 MB cube has been demonstrated).
- [ ] R3. Per-pixel validation of counts in channels 400 to 4,095 against the public RosettaSciIO path without materializing the cube (FINDINGS 9.8 item 1).
- [ ] R4. The open scientific questions of FINDINGS 9.8 (element-map processing and units, `PixelTimes` units, header sum excess, `Valid` and `LineCounter`, timestamps, footprint of RTX mosaic 0, overview-image registration, corner versus centre).
- [ ] R5. A conda-forge environment build from scratch (blocked by an SSL certificate error on the development machine; SSL verification must not be disabled).
- [ ] R6. A documented interoperability interface for Hyper Fusion (P4 last item).

## P0 — Establish reproducibility

- [ ] Verify and document the `microxrf_to_netcdf` Conda environment.
  AC: the environment can be recreated from documented steps, and the diagnostic script runs in it. Nothing is installed in the Conda base environment.
  Status: the environment `microxrf_to_netcdf` exists and the diagnostic scripts run in it; [environment.yml](environment.yml) is written, but `conda env create -f environment.yml` failed with a conda-forge SSL error (re-tried 2026-09-30 with `conda create --dry-run`: same certificate verification error; verification was not disabled), so recreation **through conda** is not verified. The pip route is verified (see the next task).
- [x] Create and maintain a complete, tested `requirements.txt`.
  AC: a fresh environment built only from `requirements.txt` runs [inspect_bruker.py](tools/diagnostics/inspect_bruker.py) successfully.
  Evidence (2026-09-30): a virtualenv created with `python -m venv` and populated with `pip install -r requirements.txt` only (rosettasciio 0.14.0, netCDF4 1.7.4, xarray 2026.7.0, numpy 2.5.3, dask 2026.8.0, ...) passed `pip check`; `inspect_bruker.py --data-dir data` and `inspect_rtx.py` produced output identical to [reference/](reference/) (except volatile time lines) and the 145 fast tests passed in it.
- [x] Preserve the existing working diagnostic script.
  AC: its output on the files in [data/](data/) is saved as a reference and still matches after any change.
  Evidence: [reference/inspect_bruker_output.txt](reference/inspect_bruker_output.txt), 2026-09-29 (file modification times in it are volatile). Script unmodified.
- [x] Add automated regression tests.
  AC: tests run with a single documented command and check the values listed in [FINDINGS.md](FINDINGS.md). Tests that need the large files skip cleanly when they are absent.
  Evidence: `python -m pytest` gives 67 passed (32 fast, 35 `large`) on 2026-09-29 (27 before the RTX work); after the unified converter, **197 passed** (145 fast, 50 `large`, 2 `slow`) on 2026-09-30. Skip behavior when `data/` is absent was exercised on 2026-09-30 on a copy of the repository without `data/`: 147 passed and the 50 `large` tests skipped cleanly (one pre-existing test, `test_peak_memory_of_decoding_is_far_below_payload_size`, failed instead of skipping and was fixed to use the `rtx_path` fixture).

## P1 — Understand the acquisition

**Next priority: the RTX investigation.** It began on 2026-09-29; the results are in FINDINGS.md section 8 and are
limited to the one RTX file in `data/`. Tasks are ticked only where a test or a recorded observation supports them.

- [x] Decode and investigate the RTX payload safely.
  AC: the payload is decoded with bounded memory, the original file is untouched, and each decoding step (Base64, zlib) is documented.
  Evidence: [inspect_rtx.py](tools/diagnostics/inspect_rtx.py) decodes file -> Base64 -> zlib -> expat in about 1.0 s and 54 to 56 MiB peak process working set (payload 62.9 MB); FINDINGS 8.1 and 8.2; tests `test_decoding_succeeds_with_recorded_sizes`, `test_input_file_unchanged_by_decoding`, `test_peak_memory_of_decoding_is_far_below_payload_size`, and the synthetic error-path tests. Output: [reference/inspect_rtx_output.txt](reference/inspect_rtx_output.txt).
- [x] Identify the decoded structure and any internal metadata.
  AC: the structure is described in FINDINGS.md using only observed evidence.
  Evidence: FINDINGS 8.3; tests in `tests/test_rtx_original.py` (`test_decompressed_payload_is_xml_text_not_binary`, `test_three_images_with_recorded_geometry`, `test_mapdaten_plane_descriptions`, `test_no_spectral_detector_stage_or_quantification_elements`).
- [x] Establish the relationship between BCF and RTX.
  AC: FINDINGS.md states, with evidence, whether they are associated and what each contains that the other lacks.
  Evidence: FINDINGS 4 and 8.4 to 8.6: the RTX video plane is bit-identical to the BCF `Video` dataset, calibrations are equal, RTX element maps track BCF spectra pixelwise (`test_rtx_video_plane_is_bit_identical_to_bcf_video`, `test_rtx_calibration_equals_bcf_axis_scale`, `test_rtx_element_planes_track_bcf_energy_windows`). The RTX has no explicit identifier of the BCF, so the link is content-based. Timestamps and filename are recorded separately as weaker evidence.
- [ ] RTX follow-up A: determine how the RTX element maps were computed.
  AC: the method (ROI window, net counts, deconvolution, other), the count units and the meaning of the `KA` / `K` suffixes are established from documentation or from a reproducible reconstruction from the BCF spectra that reproduces the RTX planes, or the question is recorded as unanswerable. Known so far: the maps are not raw +-0.15 keV window sums (FINDINGS 8.4).
- [ ] RTX follow-up B: establish the registration of the optical mosaic with the map grid.
  Progress (FINDINGS 9.3): for RTX mosaic instance 1 the axes are parallel with the same direction and the `Map` rectangle is inclusive (BCF video versus mosaic crop, r 0.9956; flips at most 0.42; best shift (0, 0)); recorded in the output file. Still open: pixel corner versus centre, the rotation of `Video Mosaic.png`, the `Map` rectangle of RTX mosaic 0, and the overview images.
  AC: the axis orientation, the origin, the inclusive/exclusive meaning of the `Rect` values, and the parallelism of mosaic and map axes are verified (for example by matching features of the mosaic and the BCF `Video` image inside the `Map` rectangle), and the rotation of `Video Mosaic.png` is explained. Known so far: the footprint size matches within one mosaic pixel (FINDINGS 8.4).
- [ ] RTX follow-up C: interpret or explicitly leave unresolved the opaque fields `Valid`, `LineCounter`, the two mosaic instances, and the header time 15:58:47.
  AC: each is either explained by evidence or listed in FINDINGS.md as unresolved.
- [ ] RTX follow-up D: test more RTX files and the BCF `EDSDatabase`.
  Progress (FINDINGS 9.1): the BCF `EDSDatabase` was searched: the SFS container holds only `HeaderData` and `SpectrumData0`; the header holds no element maps but **seven images** (RosettaSciIO returns one). No additional RTX file exists yet (R1).
  AC: at least one additional RTX (other `TRT*` types, `ItemSize` values, spectra or ROI tables if any) is inspected without code changes or with documented changes; the BCF `EDSDatabase` is searched for element maps; the findings that hold for one file only are marked as such.
- [ ] Validate spatial and spectral calibration.
  AC: pixel size and energy axis are checked against an independent source (for example known peak positions), and the outcome is documented.
  Partial (not sufficient): the RTX calibration equals the BCF one, and the BCF energy axis places the Ca, Mn, Fe and Zn K-alpha windows where the RTX maps peak, within +-0.15 keV; the mosaic scale bar shows the mosaic unit is micrometres (FINDINGS 8.4). Peak positions have not been fitted, and the BCF pixel size has no external reference.
- [x] Compare extracted images with the supplied PNG files.
  AC: a reproducible comparison states whether `CaFe.png` and `Video Mosaic.png` correspond to BCF data, and how.
  Evidence: FINDINGS 3 and 8.4; `test_cafe_png_is_the_ca_fe_composite_of_the_rtx_planes`, `test_video_mosaic_png_is_the_rtx_mosaic_rotated_clockwise`. `CaFe.png` is the Ca/Fe composite of RTX planes that agree with BCF spectra; `Video Mosaic.png` is the RTX mosaic rotated 90 degrees clockwise (r 0.95 to 0.97, not pixel-exact) and has no counterpart in the BCF. The BCF `Video` image was compared with the RTX only, not with the PNGs.
- [ ] Validate representative spectra and elemental maps.
  AC: spectra from selected pixels and maps for selected elements are inspected, and the results are recorded.
  Partial: RTX maps for Ca, Mn, Fe and Zn were compared with BCF window sums (FINDINGS 8.4). Pixel spectra are compared exactly with independent decoders (FINDINGS 9.7) but no individual spectrum has been inspected or plotted.
- [ ] Verify the original units of the recorded beam energy.
  AC: the units are confirmed from the raw metadata and documented.
  Partial: the RTX annotation `HV: 50,0 kV` supports kilovolts (FINDINGS 8.3), but it is not the BCF raw metadata; this stays open.

## P2 — Design the scientific data model

- [x] Identify all scientific variables and dimensions.
  AC: a table lists each variable, its dimensions, dtype and source (BCF or RTX).
  Evidence: [NETCDF_SCHEMA.md](NETCDF_SCHEMA.md) sections 4 to 7; a test checks that every written variable and global attribute is named there.
- [x] Establish coordinate conventions and units.
  AC: axis orientation, origin, units and spectral axis convention are documented.
  Evidence: NETCDF_SCHEMA.md section 8 (orientation verified for the mosaic, FINDINGS 9.3; corner versus centre stays unverified and is stated).
- [x] Preserve sample identity, instrument metadata and acquisition provenance.
  AC: each such field is mapped to a named attribute or variable, with its source recorded.
  Evidence: `instrument` attributes plus the complete header dictionaries and residual XML of both sources (NETCDF_SCHEMA.md sections 4, 7, 11); tests `test_metadata_and_provenance`, `test_coordinates_calibration_and_provenance`.
- [x] Define how BCF and RTX information will be combined.
  AC: a written rule covers every element of both sources, including anything that cannot be mapped and why.
  Evidence: NETCDF_SCHEMA.md section 11 (a mapping table that includes the BCF-header images RosettaSciIO does not return); unsupported layouts are refused, tested.
- [x] Determine missing-data and mask conventions.
  AC: fill values and mask semantics are documented and tested on real data.
  Evidence: no mask exists in either source; no `_FillValue` is set; valid zeros survive lazy reopening as integer 0 (synthetic and real-file tests); `Valid` is copied verbatim and not interpreted.
- [x] Design an initial NetCDF-4 schema based on verified data.
  AC: the schema is reviewed against the data model and contains no variable lacking a verified source.
  Evidence: schema 1.0.0; every variable lists its source and derived variables are marked as derived. Open interpretations are listed in NETCDF_SCHEMA.md section 12.

## P3 — Implement conversion

Design basis: [FINDINGS.md](FINDINGS.md) section 7 (sequential decoding was verified feasible for the BCF
on one file) and section 8 (the RTX payload structure is now documented for one file; the RTX stage still
waits for the data model of P2). Implement in this order and do not start a task before the
acceptance criteria of its predecessors are met. Create the `microxrf_to_netcdf` package only when the first task
needs it, and update [CLAUDE.md](CLAUDE.md) and [README.md](README.md) in the same change.

- [x] Decide whether sequential BCF decoding is feasible, and how (investigation only, no converter).
  Evidence: FINDINGS.md section 7; [probe_bcf_streaming.py](tools/benchmarks/probe_bcf_streaming.py) and
  [reference/probe_bcf_streaming_output.txt](reference/probe_bcf_streaming_output.txt).

**BCF proof of concept (verified 2026-09-29; a probe, not the converter).** Each item is supported by
[probe_bcf_streaming.py](tools/benchmarks/probe_bcf_streaming.py) output and FINDINGS.md section 7; all are for the one real BCF.

- [x] PoC: sequential scan-line decoding with a small Python framer and the unmodified RosettaSciIO compiled decoder, run on the real BCF (240 lines, 432,000 pixels, framer ends at the last byte). FINDINGS 7.2, 7.3.
- [x] PoC: incremental NetCDF-4 writing with chunks (4, 60, 4096), zlib level 4 and shuffle: about 22 s, 196.8 MiB (197 MiB) output, 260 to 303 MiB peak working set. FINDINGS 7.5 item 5.
- [x] PoC: the output reopened lazily with xarray and Dask (`chunks={}`, `uint8`, chunks (4, 60, 4096)). FINDINGS 7.5 item 5.
- [x] PoC: the tested numerical comparisons matched the original decoder exactly (all 4096 channels binned 8 x 8; every pixel over the first 400 channels; a pixel spectrum, a 3-D window and one channel image after reopening). FINDINGS 7.4, 7.5.
- [x] PoC: the original acquisition files were not modified (the probe only reads; the netCDF output went to a scratch location, never `data/`). The RTX tests additionally assert unchanged size and modification time.

**Still open after the proof of concept** (none of these is done; the ticked items above do not cover them):

- Validation of other BCF pixel encodings (`flag` 0, `flag` 1, `n_of_pulses` > 0) and zlib-compressed SFS containers: P3.2.
- Tests with datasets larger than available RAM (capped-memory and synthetic-stream tests): P3.11, the test plan below, and P4.
- A reproducible environment built from scratch: P0.
- The finalized scientific data model: P2 (FINDINGS 8.7 lists RTX-derived candidates only).
- The production converter and its comprehensive validation: P3.1 to P3.11 and P4.
- [x] P3.1 BCF line stream: production line framer plus band decoder (option A of FINDINGS 7.5). **Done 2026-09-29** (`src/microxrf_to_netcdf/bcf.py`; tests `tests/test_bcf_stream.py`, `tests/test_original_conversion.py`; FINDINGS 9.4, 9.5, 9.7). The converter-level peak for 1, 4, 8 and 20 lines is 312, 379, 643 and 1324 MiB; the decoder-only figures for 1, 8 and 30 lines are in FINDINGS 7.3.
  Public surface: an iterator of `(first_line, band)`, band decoded as `uint16` or wider, configurable
  number of lines per band, SFS and decoder versions reported.
  AC: (a) it only reads and never touches the input; (b) it fails with a clear error, never silently, on an
  unknown `flag`, an inconsistent record length, `pixel_x >= width`, trailing bytes, a short stream or an
  unsupported SFS compression; (c) the RosettaSciIO version is checked at run time against the verified
  range; (d) on the real file the streamed result equals the public decoder exactly (binned 8 x 8 over 4096
  channels and per pixel over the first 400 channels, as in FINDINGS 7.4); (e) the memory formula is
  documented and the peak measured for 1, 8 and 30 lines per band is recorded in FINDINGS.md.
- [ ] P3.2 Synthetic-stream tests for decoder branches the real file does not exercise.
  Status 2026-09-29: the synthetic part is **done** (`tests/test_bcf_stream.py`: flag 0, flag 1, flag 2, flag 2 with extra pulses, several lines, block splits of 1, 7, 100, 4064 bytes and whole, comparison with `unbcf_fast.parse_to_numpy` on the same bytes). The real second acquisition is not obtained (R1), so the task stays open and FINDINGS keeps those branches listed as unverified on real data.
  AC: fast (`not large`) tests build small artificial streams for `flag` 0, `flag` 1, `flag` 2 with
  `n_of_pulses > 0`, several lines, and a stream split at arbitrary block boundaries, and compare framer plus
  decoder with `unbcf_fast.parse_to_numpy` on the same bytes. Separately, a second real acquisition using
  those encodings, or a zlib-compressed SFS container, is obtained and verified; until then FINDINGS.md
  keeps those branches listed as unverified.
- [ ] P3.3 Full-channel, per-pixel validation of the streamed decode.
  Status 2026-09-30: partial (FINDINGS 9.8 item 1). Done: counts binned 8 x 8 over all 4,096 channels equal the public path; the first 400 channels of every pixel equal the public decoder; every pixel over all channels equals the prototype framer (same decoder); uint16 and uint32 decodes are identical. Not done: per-pixel comparison of channels 400 to 4,095 with the public decoder (it cannot return them without materializing about the whole cube).
  AC: the streamed output equals the public decoder pixel by pixel over all 4096 channels, obtained without
  materializing the cube (for example the public decoder on successive channel ranges of at most a few
  hundred channels, compared band by band). Result recorded in FINDINGS.md.
- [x] P3.4 Dtype policy. **Done** (`src/microxrf_to_netcdf/preflight.choose_counts_dtype`, exact uint32 scan in `convert.scan_bcf`; tests `test_dtype_is_chosen_from_the_observed_maximum`, `test_a_count_that_wraps_in_uint16_is_still_stored_exactly`, `test_an_explicit_dtype_that_would_truncate_is_refused_before_writing`; FINDINGS 9.4).
  AC: bands are decoded in a wider integer type; the writer stores `uint8` only if the observed maximum
  count fits, otherwise widens or aborts; the decision and the maximum are stored as provenance; a test
  proves that no wrap-around happens on a synthetic stream whose counts exceed 255.
- [x] P3.5 Add netCDF4 and xarray as tested dependencies (only when the writer task starts). **Done 2026-09-29** (netCDF4 1.7.4, xarray 2026.7.0, NumPy 2.5.3; `pip check` clean; the DeprecationWarning is filtered inside the writer and recorded in requirements.txt and FINDINGS 9.9).
  AC: [requirements.txt](requirements.txt) lists them with the reason and lower-bounded verified versions
  (evaluated so far: netCDF4 1.7.4, xarray 2026.7.0, NumPy 2.5.3); `pip check` passes in `microxrf_to_netcdf`; the
  netCDF4 `Setting the shape on a NumPy array` DeprecationWarning is recorded or resolved; nothing is
  installed in the Conda base environment.
- [x] P3.6 Incremental NetCDF-4 writer with configurable chunking and compression (schema from P2; may be
  prototyped with the probe schema). **Done** (`src/microxrf_to_netcdf/writer.py`; band height validated as a multiple of the chunk extent; no `_FillValue`; completeness checked and recorded; refuses to overwrite or write next to inputs or inside `data/`; tests in `tests/test_convert_synthetic.py`).
  AC: chunk shape, compression level, shuffle and lines per band are parameters with documented defaults;
  the band height is validated as a multiple of the chunk extent along `y`; no `_FillValue` of 0 is set on
  count variables (regression test: reopening through xarray keeps `uint8` and yields no NaN); completion is
  checked (every line written) and recorded, because unwritten chunks do not read as 0; the output goes to
  a new file, and the writer refuses to overwrite an existing file or to write inside `data/`.
- [x] P3.7 Chunk-shape decision from access patterns. **Done** (`benchmark_conversion.py chunks`, four candidates, table in NETCDF_SCHEMA.md section 9 and FINDINGS 9.6; default `(2, 30, 4096)`).
  AC: measured timings on this dataset for (i) one pixel spectrum, (ii) one energy-channel image, (iii) a
  small 3-D window, plus file size and write time, for at least three candidate chunk shapes; the chosen
  default and the reason are recorded in FINDINGS.md.
- [x] P3.8 Coordinates, units and metadata in the output (depends on P2). **Done** (tests `test_coordinates_calibration_and_provenance`, `test_metadata_and_provenance`; what a reader does not return is reported: the BCF-header images are now carried, and the SFS container holds only the two items that are read).
  AC: tests confirm that dimensions, coordinate values, units, calibration, instrument and acquisition
  metadata, and provenance (source file names and sizes, microXRF to NetCDF, RosettaSciIO, netCDF4 and xarray
  versions, conversion parameters) survive the conversion; anything read but not stored is reported.
- [x] P3.9 RTX sequential stage (payload structure documented in FINDINGS 8; waits for P2 and for the RTX follow-ups A, B and D). **Done for the one available file** (`src/microxrf_to_netcdf/rtx.py`; follow-ups A, B and D remain open as scientific questions, each recorded as unresolved in the output).
  AC: Base64 and zlib are decoded incrementally with bounded buffers (feasibility shown by [inspect_rtx.py](tools/diagnostics/inspect_rtx.py): about 55 MiB peak on the real file, no disk-backed intermediate needed); unrecognized elements and fields (`LineCounter`, `Valid`, overlays, display settings) are kept verbatim, not discarded; the element maps, mosaics and annotations are mapped to the P2 model; a file whose RTX declares another compressor, encoder or an `ItemSize` other than 1 or 2 fails with a clear error.
- [x] P3.10 Lazy reopening. **Done** (`test_lazy_reopen_with_xarray_and_dask_preserves_integer_counts_and_zeros`, `test_chunk_boundaries_and_partial_reads`, `test_dimensions_variables_and_lazy_reopen`, `test_representative_pixel_spectrum_window_and_channel_image`).
  AC: the output opens with `xarray.open_dataset(..., chunks={})` without loading data; the variable is a
  Dask array with the documented chunks; partial reads (pixel spectrum, window, channel image) equal the
  source values, checked in tests.
- [x] P3.11 Memory-bound tests, as specified in the plan below. **Done** (`tests/test_original_conversion.py::test_converter_ran_with_bounded_memory_far_below_the_cube_size`, `tests/test_bounded_memory.py` with a Windows Job Object cap on a synthetic cube larger than the cap; FINDINGS 9.5).
  AC: the tests in the plan pass and their measurements are recorded in FINDINGS.md.

### Memory and validation test plan (specification for P3.1 to P3.11; implemented 2026-09-30: items 1 to 6 map to the tests named in P3.11, `tests/test_convert_synthetic.py`, `tests/test_bcf_stream.py` and the version gate test)

1. **Bounded peak memory on the real file (`large`).** Convert the full cube in a child process and sample
   its peak working set (Windows `GetProcessMemoryInfo`, as in the probe; document a portable alternative).
   Assert `peak <= baseline + k x n x W x C x 2 bytes + header allowance`, with `k` and the allowance taken
   from FINDINGS 7.3, and that one-line bands stay far below the 1.77 GB cube.
2. **Larger than available RAM** (the real file cannot show this on a 31.5 GB machine). Two complementary
   tests: (a) run the converter in a child process under a hard memory cap (Windows Job Object
   `JOB_OBJECT_LIMIT_PROCESS_MEMORY`) of, for example, 400 MiB, smaller than the 1.77 GB cube, and require it
   to complete with correct output; (b) a synthetic-stream generator that emits a valid record stream for a
   cube whose decoded size exceeds physical RAM (only framer, decoder and writer see it; nothing is stored
   except the compressed output; low-entropy spectra), run under the same cap after recording the free disk
   space it needs.
3. **Numerical validation.** Streamed and written data equal the public decoder (the FINDINGS 7.4 checks
   and P3.3), sampled pixels and spectra match after lazy reopening, and the total count (677,896,870 for
   the current file) is preserved.
4. **Metadata validation.** Every data-model field appears in the output with its unit and source;
   provenance lists source names, sizes and library versions; a round trip through xarray preserves dtype
   (no float promotion) and attributes.
5. **Robustness.** A truncated stream, a corrupted length field, an unknown `flag`, an existing output path
   and an output path inside `data/` all fail with clear errors and leave the inputs unchanged (size and
   modification time checked before and after).
6. **Regression.** A RosettaSciIO upgrade that changes the record framing or the single-chunk behavior makes
   a test fail and forces FINDINGS.md to be revisited.

## P4 — Validate scientific interoperability

- [x] Compare representative input and output values.
  AC: sampled pixels and spectra are identical between source and output.
  Evidence: FINDINGS 9.7 (every band compared exactly after writing; independent comparisons with the public path and the prototype framer).
- [x] Verify spectral axes, spatial dimensions, orientation and calibration.
  AC: each is checked against the source and documented.
  Evidence: coordinates equal the raw calibration (checked in every conversion and against the public reader); mosaic orientation checked against the BCF video (r 0.9956) and `Video Mosaic.png` (rotation r 0.95 to 0.97). External validation of the calibration itself is still open (FINDINGS 9.8).
- [x] Check metadata and provenance preservation.
  AC: every field in the data model appears in the output, including source file names and sizes.
  Evidence: `check_provenance` in `python -m microxrf_to_netcdf validate` (hashes and sizes recomputed from the inputs) and the schema-documentation test.
- [x] Measure peak RAM usage and conversion performance.
  AC: measurements are recorded with hardware and software versions.
  Evidence: FINDINGS 9.5 to 9.7 (Windows 11, 31.5 GB RAM, Python 3.12.14, versions in the conversion report).
- [ ] Test acquisitions larger than available memory.
  AC: a conversion completes on a dataset larger than RAM while staying within the bound.
  Partial: a synthetic 786 MB (logical) cube converted under a hard 350 MiB and 500 MiB process-memory cap, and failed cleanly at 250 MiB (FINDINGS 9.5). A dataset larger than the physical RAM (31.5 GB) is R2.
- [ ] Assess future interoperability with Hyper Fusion.
  AC: a documented interface (file format and metadata contract) is described, without sharing internal code with Hyper Fusion or HSI Toolbox.
