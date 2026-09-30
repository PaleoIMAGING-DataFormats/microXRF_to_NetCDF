# microXRF to NetCDF — Verified Findings

This document records only what has been verified by the diagnostic scripts
[inspect_bruker.py](tools/diagnostics/inspect_bruker.py), [probe_bcf_streaming.py](tools/benchmarks/probe_bcf_streaming.py) and
[inspect_rtx.py](tools/diagnostics/inspect_rtx.py), and by the tests in [tests/](tests/), on the files in [data/](data/).
Each item is classified as one of:

- **Confirmed** — directly observed and reproducible.
- **Derived** — simple arithmetic on confirmed values, with no new assumption.
- **Interpretation** — plausible reading of confirmed values, not yet verified.
- **Limitation / Unresolved** — known gap or open question.

Do not promote an Interpretation to Confirmed without a reproducible test.

## 1. BCF acquisition

File: `GRF17A_9-29cm_slab3_Elemental_map.bcf`

### Confirmed

| Property | Value |
|---|---|
| File size | 617,320,728 bytes |
| Signature | `AAMVHFSS` |
| Reader | RosettaSciIO 0.14.0 reads the file successfully |
| Video dataset | `uint16`, shape (240, 1800) |
| Spatial sampling | approximately 100.014812 micrometers per pixel |
| EDX dataset | lazy Dask array, `uint8`, shape (240, 1800, 4096) |
| Energy axis | keV, scale 0.010001, offset -0.96079607 |
| Detector | XFlash 430 |
| Recorded beam energy | `50` (original units not yet verified) |
| Elevation angle | 50.0 degrees |
| Real acquisition time | 12960 seconds |
| Sample metadata | `Mapdaten` |
| EDX acquisition timestamp | 2026-07-30 09:55:52 |

The reader returned **no explicit elemental-map dataset**.

Verified on 2026-09-29 in the `microxrf_to_netcdf` environment (see section 6) and pinned by
[tests/test_original_files.py](tests/test_original_files.py):

- The reader returns exactly two datasets, `Video` and `EDX`; no other dataset.
- The height and width axis scales are identical in the reader output (100.0148121593, unit
  `µm`, offset 0). The raw `Microscope.DX` / `DY` values differ only beyond the 12th
  significant digit (100.014812159301 vs 100.0148121593).
- The EDX axis order is (height, width, Energy) and equals the array shape (240, 1800, 4096).
- **The lazy EDX array is a single Dask chunk.** With `lazy=True`, RosettaSciIO 0.14.0 wraps one
  `dask.delayed` parse call in `da.from_delayed`: `numblocks == (1, 1, 1)`, chunk
  240 x 1800 x 4096, 1,769,472,000 bytes. It is lazy until first access but not chunked.
  Importing it and running the whole test suite peaked at about 196 MiB working set (no data
  was computed).
- `inspect_bruker.py` requests `select_type="images"`, which is **not a valid value** for the
  reader (valid: `'image'`, `'spectrum_image'`, `None`). The reader falls through to its
  default branch and returns images and the spectrum image. The diagnostic's "image-only"
  wording is therefore inaccurate; the returned datasets are those of `select_type=None`.

### Derived

- 240 x 1800 pixels at about 100.0148 micrometers per pixel corresponds to about 24.0 mm x 180.0 mm.
  The two axis scales are equal in the reader output (see above).
- Consequence of the single-chunk finding (derived from the Dask graph, **not** executed):
  any slice of the RosettaSciIO lazy EDX array would trigger a full decode of the cube
  (about 1.77 GB as `uint8`, plus decoder working memory). RosettaSciIO's Dask array therefore
  cannot by itself provide bounded-memory chunked access. A bounded-memory route through the
  same decoder was verified afterwards (section 7).
- The energy axis spans about -0.961 to about 40.0 keV over 4096 channels
  (offset + 4095 x 0.010001).

### Limitations / Unresolved

- The units of the recorded beam energy value (`50`) must be checked against the original metadata before use.
- `uint8` per-channel counts: a streamed `uint16` decode of all 4096 channels found a maximum per-pixel,
  per-channel count of 149 (section 7), so `uint8` holds every count of this file without wrapping.
  This is a property of this file only; RosettaSciIO's dtype choice is a heuristic (section 7).
- Spatial and spectral calibration have not been validated against an independent source.

## 2. RTX file

File: `GRF17A_9-29cm_slab3_Elemental_map.rtx`

### Confirmed

| Property | Value |
|---|---|
| File size | 43,723,950 bytes |
| Encoding | valid Windows-1252 XML |
| Root element | `TRTProject` |
| Sections | `RTHeader`, `RTData` |
| `RTData` declaration | Base64 encoding, zlib compression |
| Recorded timestamp | 2026-07-30 15:58:47 (`Date` is written `30.7.2026`) |
| Structure | `TRTProject` > `RTHeader` (`ProjectHeader` > `Date`, `Time`, `Creator`, `Comment`; `RTCompression`) and `RTData`; 9 elements in total |
| `RTCompression` attributes | `compressor="zlib"`, `encoder="base64"` |
| Payload | `RTData` is one Base64-like text element; the file contains 0 NUL bytes |
| Instrument, software, sample, scan, coordinate fields | not present in the XML text (`inspect_bruker.py` reports "not found") |

The payload decoding and its content are recorded in section 8 (2026-09-29).

### Limitations / Unresolved

- ~~The decompressed structure of the `RTData` payload has not been verified.~~ Verified in section 8: the
  payload is a second XML-like document (`CompData`), 62,915,687 bytes.
- ~~Nothing is known about what scientific content the RTX holds.~~ Section 8.3 lists it: three images
  (two optical mosaics and one 15-plane map image) plus annotations. It holds no spectra.
- The `RTHeader` timestamp (15:58:47) is not one of the three image timestamps inside the payload; its meaning
  is unknown (section 8.5).

## 3. Companion images

### Confirmed

| File | Size in pixels |
|---|---|
| `CaFe.png` | 1800 x 240 |
| `Video Mosaic.png` | 948 x 6000 |

### Confirmed by comparison with the RTX payload (section 8.4; tests in `tests/test_rtx_original.py`)

- `CaFe.png` is an RGB composite of the RTX planes `Ca-KA` (blue channel, Pearson r 0.98) and `Fe-KA` (red
  channel, r 0.91); the green channel is unused (|r| < 0.1 with both). It sits on the 1800 x 240 grid of the
  RTX map image. The colour assignment agrees with the RTX display settings (`MapImage1` colour 16711680 =
  0xFF0000, read as BGR = blue; `MapImage8` colour 255 = red).
- `Video Mosaic.png` (948 x 6000) is the RTX `Video Mosaic` image (6000 x 948) **rotated 90 degrees
  clockwise**, r 0.95 to 0.97 per channel against RTX planes 0, 1, 2 respectively. It is not pixel-exact
  (the opposite rotation and the plain transpose give r below 0.7); the cause of the difference (rendering,
  resampling or colour processing) is unknown.

### Limitations / Unresolved

- The PNG files are not identical to any decoded plane, so their processing history is unknown.
- ~~The BCF holds only its own `Video` image (240 x 1800); it exposes no mosaic (section 8.4).~~ Correction (section 9.1): RosettaSciIO **returns** only the video, but the BCF header contains seven images, among them an overview mosaic byte-identical to the RTX mosaics and a `PixelTimes` image.

## 4. BCF / RTX association

### Confirmed

- The two files share an identical filename stem and the calendar date 2026-07-30 (weak evidence).
- The RTX `Mapdaten` image plane 0 (`Video 1`) is **bit-identical** to the BCF `Video` dataset
  (`np.array_equal`; 240 x 1800, `uint16`, 864,000 bytes; section 8.4). This is content-level evidence of a
  common acquisition, not name-level.
- The RTX X/Y calibration of that image equals the BCF axis scales, and its size equals the BCF grid.
- The RTX element planes agree pixel by pixel with counts summed from the BCF spectra (Pearson r on the first 40
  scan lines: Ca 0.9996, Mn 0.9979, Fe 0.9969, Zn 0.9932; over the whole file Ca 0.9996, Fe 0.9995);
  section 8.4.
- The RTX contains **no explicit reference** to the BCF (no file name, path, identifier or GUID element).
- Recorded times: BCF EDX 09:55:52; RTX images 9:38:26, 9:55:51, 9:55:53; RTX header 15:58:47.

### Limitations / Unresolved

- The association is now **established from content** (bit-identical video image, identical calibration,
  pixel-level agreement of element maps with BCF spectra), not from an identifier. Timestamps and filename
  similarity remain separate, weaker evidence.
- The meaning of each timestamp (start, end, save time) and their time zone are unknown. The BCF EDX time lies
  one second from the RTX `Mapdaten` image time and one second from the second mosaic.

## 5. Open questions

1. ~~What does the decoded RTX payload contain?~~ Answered in section 8.
2. ~~Does the RTX duplicate, complement, or contradict the BCF content?~~ Section 8.6: it duplicates the video
   image and the calibration, and complements them with 14 element maps and two mosaics; no contradiction was
   found. Still open: how the RTX element maps were computed.
3. What are the original units of the beam energy? Partial: the RTX annotation reads `HV: 50,0 kV` (section
   8.3), consistent with kV, but the BCF raw metadata unit has still not been read.
4. Do the video image and the PNG files spatially correspond to the EDX grid? Partial: the RTX stores its video
   plane and the element planes on one 1800 x 240 grid, the element planes agree pixelwise with BCF spectra,
   and `CaFe.png` matches the RTX planes (sections 3 and 8.4). `Video Mosaic.png` is on a different grid.
5. Are the spatial and spectral calibrations correct?
6. Why is no elemental-map dataset exposed by the reader — are none stored, or are they not parsed? The RTX
   holds 14 element maps (section 8.3). Answered for this file (section 9.1): the BCF header holds no element maps; it does hold seven images, of which RosettaSciIO returns one.
7. ~~How can the EDS cube be read in bounded memory?~~ Answered for this file in section 7: the unmodified
   compiled decoder can be fed a few scan lines at a time. Still open: files with other pixel encodings
   (`flag` 0/1, extra pulses) and zlib-compressed SFS containers, which this file does not exercise.
8. ~~Is the sequential design also applicable to the RTX payload?~~ Yes, verified: section 8.2 decodes it with
   bounded buffers (54 to 56 MiB peak process working set for a 62.9 MB payload).
9. What are the processing method and units of the RTX element maps (ROI window, net counts, deconvolution)?
10. What do the plane fields `Valid` and `LineCounter` mean, and what distinguishes mosaic 1 from mosaic 2
    beyond time and annotations?
11. Why is `Video Mosaic.png` rotated relative to the RTX mosaic, and are mosaic and map axes parallel? Partly answered (section 9.3): the axes are parallel with the same direction (r 0.9956, inclusive `Map` rectangle); the rotation of the PNG is still unexplained.

## 6. Verified environment and diagnostic run

Recorded 2026-09-29 on Windows 11, 31.5 GB RAM.

| Package | Version |
|---|---|
| Python | 3.12.14 (conda-forge) |
| RosettaSciIO | 0.14.0 |
| Dask | 2026.8.0 |
| NumPy | 2.5.3 |
| Matplotlib | 3.11.2 |
| Pillow | 12.3.0 |
| Pint | 0.26.1 (RosettaSciIO dependency) |
| pytest | 9.1.1 |
| netCDF4 | 1.7.4 (HDF5 1.14.6, netCDF-C 4.9.3), added 2026-09-29 (section 9.9) |
| xarray | 2026.7.0, added 2026-09-29 |
| pandas | 3.0.6 (xarray dependency) |
| cftime | 1.6.6 (netCDF4 dependency) |

- The `microxrf_to_netcdf` environment is Python 3.12 with the packages of `requirements.txt`; the package is
  installed editable (`pip install -e . --no-deps`). A from-scratch `conda env create -f environment.yml` **failed** here with a conda-forge
  SSL certificate verification error, so the from-scratch path is not yet verified
  (`pip install --dry-run -r requirements.txt` resolves and `pip check` reports no broken
  requirements in `microxrf_to_netcdf`).
- `python tools/diagnostics/inspect_bruker.py --data-dir data` exits with status 0. Its output is saved in
  [reference/inspect_bruker_output.txt](reference/inspect_bruker_output.txt). File modification
  times in the inventory reflect when the files were copied and are not stable across copies.
- `pytest` before the RTX work: 27 passed (9 fast, 18 marked `large`); about 1.4 s, peaking at about
  196 MiB working set. After the RTX work (section 8): **67 passed** (32 fast, 35 `large`) in about 8 s
  (`tests/test_rtx_synthetic.py`: 23 fast; `tests/test_rtx_original.py`: 17 `large`).
  Latest full run (2026-09-30, after the unified converter of section 9): **197 passed** (145 fast, 50 `large`, 2 `slow`) in about 148 s. The 145 fast tests also pass in a fresh virtualenv built from `requirements.txt` alone (about 47 s).

## 7. Sequential (streaming) BCF decoding and incremental NetCDF-4 writing

Investigation of 2026-09-29. Random access to the BCF is not needed; the goal is bounded-memory
sequential conversion. The original files and the installed RosettaSciIO 0.14.0 were not modified;
no converter exists yet. Evidence is produced by [probe_bcf_streaming.py](tools/benchmarks/probe_bcf_streaming.py); its
output on this machine (Windows 11, 31.5 GB RAM, `microxrf_to_netcdf` environment) is saved in
[reference/probe_bcf_streaming_output.txt](reference/probe_bcf_streaming_output.txt). The probe is a
diagnostic with a throw-away line framer, **not** the reader or the converter.

### 7.1 How the RosettaSciIO decoder works (Confirmed by reading `rsciio/bruker/_api.py` and `unbcf_fast.pyx`)

- The BCF is an SFS container. `SFS_reader(filename).get_file("EDSDatabase/SpectrumData0")` returns an
  item whose `get_iter_and_properties()` yields the internal file as a **sequential Python iterator of
  byte blocks** (one file open and read per block). That layer is already a streaming interface. Only
  the header item (`EDSDatabase/HeaderData`, 34,006,136 bytes) is read whole; constructing `BCF_reader`
  peaked at about 182 MiB working set.
- In this file the SFS container is **not** zlib-compressed (`SFS_reader.compression == "None"`): the
  spectrum stream is 577,830,250 bytes in 142,183 SFS blocks of at most 4064 bytes. The zlib branch
  (`_iter_read_compr_chunks`) exists in the code but is **not exercised by this file** and was not tested.
- The compiled decoder `unbcf_fast.parse_to_numpy(virtual_file, shape, dtype, downsample=1)` is monolithic:
  it calls `np.zeros(shape)` itself, then `bin_to_numpy` loops over the `height` lines whose count it reads
  from the first 4 bytes of the stream, adding each pixel's spectrum into that array. It has no callback,
  no generator and no line-range argument. `dask.delayed` in `parse_hypermap` merely wraps this one call,
  which is why the Dask array is a single chunk.
- The stream is strictly sequential and self-framing (the decoder never seeks backwards). Layout read by the
  decoder: 8-byte header (`height`, `width` as `uint32`), seek to `0x1A0`, then per line a `uint32` pixel
  count followed by one record per pixel: a 22-byte pixel header (`pixel_x` uint32, two uint16 channel
  fields, 4 skipped bytes, `flag` uint16, `data_size1` uint16, `n_of_pulses` uint16, `data_size2` uint32)
  and a payload. `flag` 0 (16-bit pulses) and 1 (12-bit pulses) consume `data_size2` bytes; other flags use
  the "instructed" packing and consume `data_size2 - 4` bytes plus 4 bytes (or, if `n_of_pulses > 0`,
  `4 + 2 * n_of_pulses` bytes).
- Decoder hazards: it is compiled with `boundscheck(False)` (a stream whose `pixel_x` or line count exceeds
  `shape` writes out of bounds instead of raising); counts are added in place in the chosen `dtype`, so a
  too-small dtype **wraps silently**; the dtype comes from the heuristic `estimate_map_depth` (maximum of
  the sum spectrum divided by the pixel count, times 2); channels at or above `shape[2]` are dropped
  silently.

### 7.2 Confirmed on this file

| Observation | Value |
|---|---|
| Line structure | 240 lines, each announcing 1800 pixels; 432,000 pixels in total; the framer ends **exactly** at the last byte (0 trailing bytes) |
| Pixel encoding present | every pixel has `flag == 2` (instructed packing) and `n_of_pulses == 0`; `pixel_x` max 1799; `data_size2` 216 to 1752 bytes |
| Not exercised | `flag` 0 and 1, and `n_of_pulses > 0`; the framer's rules for them come from reading the code only |
| Largest scan-line record | 2,710,881 bytes (about 2.6 MiB raw) |
| Public-path decode with `downsample=8` | (30, 225, 4096) `uint16`, about 1.9 s, 182 MiB peak working set (includes header parsing) |
| Maximum per-pixel, per-channel count (streamed `uint16` decode, all channels) | 149; total counts 677,896,870 |

### 7.3 Smallest practical unit of sequential decoding (Confirmed)

- **The scan line is the natural unit**: the stream is line-ordered and each line announces its pixel count.
  The unmodified decoder cannot be started mid-line because it reads the line loop from the stream header,
  so the smallest unit reusable without recompiling is **one or more whole lines**.
- Technique verified: a small Python framer walks the record structure without decoding spectra and collects
  the raw bytes of `n` lines. They are prefixed with the original `0x1A0`-byte header whose `height` field
  is patched to `n`, and passed to the **unmodified** compiled `parse_to_numpy` through a duck-typed object
  exposing `get_iter_and_properties()`. The decoder returns an `(n, width, channels)` array for exactly
  those lines. No installed file is altered.
- Memory of one band (`n` lines, 4096 channels, decoded as `uint16` so counts cannot wrap) is
  `n x 1800 x 4096 x 2` bytes = 14.1 MiB per line, plus the raw record of those lines (about 2.6 MiB per
  line) and SFS blocks (4064 bytes each). Measured with `probe_bcf_streaming.py memory`:

| Lines per band `n` | Band size (uint16) | Peak working set of the whole process |
|---|---|---|
| 1 | 14.1 MiB | 85.6 MiB (about 48 MiB before streaming starts) |
| 8 | 112.5 MiB | 319.7 MiB |
| 30 | 421.9 MiB | 1043.6 MiB |

  The peak is about 2.4 x the band size plus a fixed ~50 MiB (band, `np.zeros` and temporaries of
  `band.max()`); this is a measurement, not a guarantee. The full cube as `uint8` is 1,769,472,000 bytes;
  one-line streaming peaks at about 5 % of that.
- Speed: streaming all 240 lines takes about 3.7 s (framing about 1.7 s in Python, the rest decoding),
  against about 1.9 s for the single public call. Decoding is not the bottleneck; writing is (7.5).

### 7.4 Numerical validation of the streamed decode (Confirmed, this file only)

Two comparisons against the unmodified public decoder, both bounded in memory:

1. **All 4096 channels, binned 8 x 8:** streamed one line at a time and binned in `int64`, the result equals
   `BCF_reader.parse_hypermap(downsample=8)` **exactly** (`np.array_equal` is `True`).
2. **Per pixel, first 400 channels, all 432,000 pixels:** streamed bands equal `parse_to_numpy` run on the
   real stream with `shape=(240, 1800, 400)`, exactly.

These checks cannot detect an error that cancels inside an 8 x 8 bin at channels above 399. A stronger
per-pixel check of the full-channel output is listed in [TODO.md](TODO.md) (P3).

### 7.5 Answers to the design questions

1. **Can the decoder yield data progressively?** Not through its public Dask interface or its function
   signature. Yes through the technique in 7.3, at scan-line granularity, with the unmodified decoder.
2. **Smallest practical unit and memory:** one scan line; 14.1 MiB decoded (`uint16`, 4096 channels); 85.6 MiB
   measured peak for the whole process.
3. **Reuse or adaptation.** Options, best first:
   - **A (recommended now): framer plus the unmodified compiled decoder.** New code is only the framing
     walk (about 30 lines, structural, no spectral decoding). No compiler and no fork are needed. Its cost:
     the record-size rules are mirrored from `unbcf_fast.pyx`, so they must be guarded by tests (7.2 lists the
     branches this file does not exercise) and by a RosettaSciIO version pin checked at run time.
   - **B: upstream contribution to RosettaSciIO**: a line-range or generator variant of `parse_to_numpy`
     (yield after each line in `bin_to_numpy`; a small change in the `.pyx`). It removes the framer and the
     duplication but depends on a release cycle. RosettaSciIO is GPL-3.0-or-later.
   - **C: vendored, modified copy of `unbcf_fast.pyx`.** Needs Cython and a C compiler (neither `cl` nor `gcc`
     was found on `PATH`) and duplicates the parser. Not recommended.
   - **Rejected:** Dask rechunking or slicing (one delayed call decodes everything first); `downsample` or
     `cutoff_at_kV` as the conversion route (they discard information); reading the full cube.
4. **Proposed pipeline** (the RTX stage was verified afterwards; see section 8):

```
BCF: SFS block iterator -> line framer -> n raw lines -> unmodified decoder -> uint16 band (n, W, C)
                                                      -> range check -> target dtype (uint8 only if max <= 255)
RTX: Base64 + zlib incremental decode + expat (verified, section 8) -> 3 images, 62.9 MB XML-like payload
                                   |
                scientific data model (P2: dimensions, coordinates, units, provenance)
                                   |
                netCDF4 incremental writer: variable[y0:y1] = band, bands aligned to chunks
                                   |
                new .nc file, reopened lazily with xarray + Dask
```

   Memory bound: `peak ~ baseline + k x n x W x C x 2 bytes` with `k` about 2.4 measured (7.3). `n` is the
   user-facing parameter and must be a multiple of the chunk extent along `y`. Header parsing (about
   180 MiB) is a separate fixed term.
5. **netCDF4 Python library** (netCDF4 1.7.4, xarray 2026.7.0, installed in a scratch `pip --target`
   directory, **not** in `microxrf_to_netcdf`, so `requirements.txt` is unchanged): direct incremental writes work.
   `Variable.__setitem__` with a slice along `y` writes one band at a time; `chunksizes`, `zlib`, `complevel`
   and `shuffle` are per-variable options. For `counts` (`uint8`, 240 x 1800 x 4096), chunks (4, 60, 4096),
   zlib level 4 with shuffle, bands of 4 lines (`probe_bcf_streaming.py netcdf`):

| Quantity | Value |
|---|---|
| Write time | about 22 s for the whole cube |
| Output size | 196.8 MiB (uncompressed `uint8` cube: 1687.5 MiB) |
| Peak working set of the writing process | 260 to 303 MiB over three runs |
| Lazy reopen | `xarray.open_dataset(..., chunks={})` gives a Dask array, `uint8`, chunks (4, 60, 4096) |
| Values | a pixel spectrum, a 3-D window and a whole energy-channel image equal the reference; the whole file binned 8 x 8 equals the public-path reference |

   A scratch run with chunks (4, 20, 512) gave a similar result (188.1 MiB, 21.5 s, 269 MiB peak); the
   default chunk shape must be chosen from access patterns (P3). Whole-file reductions through xarray/Dask
   peaked at about 0.9 to 1.0 GiB working set on reopen (multi-threaded Dask; this is reader memory, not
   converter memory). Observed pitfalls:
   - **Do not set `_FillValue = 0` on count variables.** xarray's CF decoding then turned the variable into
     `float32` and every zero count into `NaN` (observed: values differed and the dtype changed). The probe
     passes `fill_value=False` (no `_FillValue` attribute). The netCDF default fill of an unwritten `ubyte`
     chunk is not 0, so completeness of the write must be verified (P3).
   - netCDF4 1.7.4 with NumPy 2.5.3 emits `DeprecationWarning: Setting the shape on a NumPy array has been
     deprecated in NumPy 2.5` on band assignment. It works today; watch it when pinning versions.
   - The output was written to a scratch location, never inside `data/`.
6. **Disk-backed intermediate representation:** not needed, because direct streaming is technically feasible
   for the BCF (7.3 to 7.5). Do not add one; reconsider only if an RTX payload proves not decodable
   sequentially. This did not happen: `RTData` (Base64 + zlib) was decoded incrementally with
   `zlib.decompressobj` over Base64 chunks and parsed as a stream, peaking at 54 to 56 MiB (section 8.2), so
   no disk-backed intermediate is needed for this file.

### 7.6 Technical limitations and open risks

- All evidence comes from **one** file with one pixel encoding (`flag` 2, no extra pulses) and an
  uncompressed SFS container. The other decoder branches and the zlib SFS branch are unverified. The framer
  must fail loudly on any unseen `flag` or inconsistent length, and needs synthetic test streams for those
  branches.
- Framer correctness rests on mirroring `unbcf_fast.pyx`; a RosettaSciIO release may change it.
- `uint8` suffices for this file only. The converter must decode to a wider type, check the maximum against
  the target dtype and refuse or widen; it must not trust `estimate_map_depth`.
- "Larger than RAM" cannot be demonstrated with this 1.77 GB cube on a 31.5 GB machine. The memory bound was
  measured only up to 30 lines per band. [TODO.md](TODO.md) P3 specifies a capped-memory test.
- `EDSDatabase/HeaderData` (34 MB) is parsed whole by RosettaSciIO; that memory is not chunked.
- Working-set numbers come from Windows `GetProcessMemoryInfo`; they include the interpreter and imported
  libraries, and the method is not portable.

## 8. RTX payload decoding and content

Investigation of 2026-09-29 in the `microxrf_to_netcdf` environment (Windows 11, 31.5 GB RAM), with
[inspect_rtx.py](tools/diagnostics/inspect_rtx.py) (read-only, standalone). Its output is saved in
[reference/inspect_rtx_output.txt](reference/inspect_rtx_output.txt). Tests: `tests/test_rtx_synthetic.py`
(23 fast, synthetic files) and `tests/test_rtx_original.py` (17 `large`). The original files and the installed
RosettaSciIO were not modified; nothing was written into `data/`. Classification follows the legend at the top.

### 8.1 Decoding (Confirmed)

| Step | Observation |
|---|---|
| Outer file | Windows-1252 XML; `RTHeader/RTCompression compressor="zlib" encoder="base64"`; `RTData` text starts at byte 300; after `</RTData>` only `\r\n</TRTProject>\r\n` |
| Base64 | 43,723,624 characters (43,723,950 bytes minus the XML markup), no invalid character, length a multiple of 4 |
| zlib | one complete stream; end-of-stream reached; 0 unused bytes after it |
| Decompressed size | **62,915,687 bytes** (the zlib stream is about 32.8 MB, about 52 % of the payload; derived from the Base64 length) |
| SHA-256 of the decompressed payload | `832f6f4c090e5413033621335127f7095ceb9ff7f9f7c3604bca2ffba00d4b61` |
| First bytes | `3c 43 6f 6d 70 44 61 74 61 3e ...` = the ASCII text `<CompData><ClassInstance Type="TRTProject" Name="Bruker project">` |
| Nature of the payload | **XML-like text, not binary and not a raw array.** It is well formed when parsed as Windows-1252 (the payload has no XML declaration; parsing it as UTF-8 fails at the byte `0xB5` of `µm`, decompressed offset 45,525,498) |
| Well-formedness | 5,305 elements, maximum depth 14, one root `CompData` |

The Windows-1252 reading is confirmed only by the byte `0xB5` in `40000 µm` decoding correctly; the payload
does not declare its own encoding.

### 8.2 Memory and time (Confirmed)

Method: `inspect_rtx.py` reads the file in 1 MiB chunks, decodes Base64 in 4-character groups, feeds
`zlib.decompressobj`, and passes each decompressed chunk to `xml.parsers.expat` (Windows-1252 override).
Neither the encoded nor the decoded payload is held whole and no temporary file is needed. The largest
retained object is one decoded image plane (at most 5,688,000 bytes; its Base64 text of 7,584,000 characters
arrives in expat pieces of at most about 1 MiB).

| Quantity | Value |
|---|---|
| Time | about 1.0 to 1.1 s for the whole file (statistics, SHA-256 and hashing of every plane included) |
| Peak working set of the whole process | 54.0 to 56.2 MiB over runs (Python, NumPy and expat included; `GetProcessMemoryInfo`, Windows only) |
| Payload for comparison | 62.9 MB decompressed, 43.7 MB encoded |
| Test bound | a fresh-process run must peak below 150 MiB (`test_peak_memory_of_decoding_is_far_below_payload_size`) |

The sequential design of section 7 therefore also holds for the RTX: `Base64 -> zlib -> XML` is decodable with
bounded buffers on this file. A file whose single plane is very large would need the plane text handled in
pieces; the current decoder keeps one plane. Nothing else is verified about larger RTX files.

### 8.3 Content (Confirmed)

The payload is a tree of `ClassInstance` elements, each with a `Type` attribute (and often `Name`). Counts by
type: `TRTProject` 1, `TRTImageData` 3, `TRTImageOverlay` 7, `TRTImageConfigurationData` 3,
`TRTRectangleOverlayElement` 6, `TRTLineOverlayElement` 6, `TRTTextOverlayElement` 7, `TRTMapConfigurationData` 1,
`TRTImagePalette` 1, `TRTMapViewSettings` 1. The `TRT` prefix and the names are the file's own; their meaning
beyond the observed content is not assumed. All 86 distinct element names (numbered names collapsed) are
listed in the reference output.

**Three images** (`TRTImageData`). Each stores its pixels as one Base64 text per plane, **not compressed
further**: the decoded length equals `Width x Height x ItemSize` exactly for every plane (`Size` element,
`PlaneCount` and geometry are consistent for all three images).

| # | Name | Payload offset | Width x Height | ItemSize (dtype) | Planes | X/Y calibration | Date, Time |
|---|---|---|---|---|---|---|---|
| 0 | `Video Mosaic` | 86 | 6000 x 948 | 1 (`uint8`) | 3 | 31,308984675955 / 31,3089846759549 | 30.7.2026, 9:38:26 |
| 1 | `Video Mosaic` | 22,761,856 | 6000 x 948 | 1 (`uint8`) | 3 | same | 30.7.2026, 9:55:51 |
| 2 | `Mapdaten` | 45,528,486 | 1800 x 240 | 2 (`uint16`, little-endian assumed) | 15 | 100,014812159301 / 100,0148121593 | 30.7.2026, 9:55:53 |

Numbers use a decimal comma. Plane byte offsets inside the decompressed payload are in the reference output
(for example `Mapdaten` `Data` elements begin at 45,535,007 for plane 0 and 61,745,995 for plane 14).

- **Images 0 and 1 have identical pixels** (SHA-256 of each of the 3 planes equal; plane sums 560,840,716,
  526,036,488 and 525,967,475). They differ in `Time` and in annotations (below). The unit of the mosaic
  calibration is not written next to it; it is micrometres by derivation (8.4).
- **Image 2 (`Mapdaten`) plane names** (`Description`): plane 0 `Video 1`, then the fourteen maps
  `Ca-KA, K-KA, S-KA, Si-K, Ti-KA, Cr-KA, Mn-KA, Fe-KA, Ni-KA, Zn-KA, Sr-KA, Pd-KA, Rh-KA, Al-K`. The names are
  the file's own; the meaning of the suffixes `KA` and `K` is not verified.
- Recorded plane statistics (`uint16`): `Video 1` min 3263, max 65280, sum 11,806,561,069; `Ca-KA` max 1607,
  sum 233,524,498; `Fe-KA` max 2159, sum 44,031,298. Every element plane has 0 as minimum; full table in the
  reference output.
- **Annotations.** Text overlays: image 1 `40000 µm`, `Mosaik`; image 2 `40000 µm`, `Mapdaten`,
  `HV: 50,0 kV`, `Ca`, `Fe`. Rectangle overlays named `Map` mark a region on each mosaic (below). The
  configuration element lists the label fields `MAG,HV,WD,Scale,Name,Comment`, but **no value for magnification
  or working distance is stored** anywhere.
- **Display settings** (image 2): `MapUsed` is 1 only for `MapImage1` (colour 16711680) and `MapImage8` (colour
  255); gamma 1,07407407462597; map filter `None`, width 7; palettes of 256 entries per plane. These are
  presentation parameters for the Ca/Fe composite, not measurements.
- **Opaque fields.** Every plane has `Valid` = 0. `LineCounter` is a comma-separated list with one entry per
  image line: `1` for every line of the video planes, `3` for every line of every element plane. Their meaning
  is unknown; they are recorded, not interpreted.
- **Absent from the payload (negative finding for this file, pinned by a test that scans the element
  vocabulary):** spectra, spectral or energy axes, detector, instrument, software or version fields, sample
  identifiers other than the name `Mapdaten`, stage coordinates, quantification tables, peak lists, file names
  or paths, GUIDs. The only regions are drawn overlay rectangles. The outer `RTHeader` `Creator` and `Comment`
  elements are empty.
- **Beam energy.** The annotation `HV: 50,0 kV` is a text label of image 2; it is consistent with the BCF value
  `50` being kilovolts, but the BCF raw metadata unit remains unverified (TODO.md P1).

### 8.4 Comparison with the BCF and the PNG files (Confirmed unless marked)

| Test | Result |
|---|---|
| RTX `Mapdaten` plane 0 vs BCF `Video` dataset (RosettaSciIO, `select_type="image"`) | **bit-identical** (`np.array_equal`, `uint16`, 240 x 1800) |
| RTX `Mapdaten` X/Y calibration vs BCF height/width axis scale | equal to relative 1e-12 (100.0148121593 µm); width 1800 and height 240 equal the BCF axis sizes |
| RTX element planes vs BCF spectra, first 40 scan lines, streamed one line per band with the unmodified decoder, windows of +-15 channels (+-0.15 keV) around tabulated K-alpha energies converted with the BCF energy axis | Pearson r (the tested 40 lines): Ca 0.9996, Mn 0.9979, Fe 0.9969, Zn 0.9932. Control: `Al-K` plane vs the Fe window r 0.48 (test threshold 0.5, so this control is weak). Whole-file exploratory run, 14 elements: Ca 0.9996, Fe 0.9995, Mn 0.9962, Zn 0.9919, Ni 0.976, Sr 0.972, Ti 0.964, K 0.948, Rh 0.932, S 0.930, Cr 0.919, Pd 0.906, Si 0.865, Al 0.822 |
| Same comparison, totals | RTX sums are not equal to the window sums (Ca 233.5 M vs 252.7 M; Fe 44.0 M vs 46.6 M; Sr 4.34 M vs 4.10 M), so the RTX maps are a **processed product** with an unknown method, not window sums |
| `CaFe.png` vs RTX planes | blue channel with `Ca-KA` r 0.98, red channel with `Fe-KA` r 0.91, green unused |
| `Video Mosaic.png` vs RTX image 1 | equals the mosaic rotated 90° clockwise, r 0.95 to 0.97 (not pixel-exact) |
| Derived: `Map` rectangle of mosaic 1 (Left 192, Top 89, Right 5941, Bottom 855) | width 5750 px x 31,308984675955 = 180,027 µm and height 767 px = 24,013 µm (inclusive rectangle), against the BCF 1800 x 100,0148 = 180,027 µm and 240 x 100,0148 = 24,004 µm: equal within one mosaic pixel (31,3 µm) |
| Derived: scale bar `40000 µm` | its line spans x = 4541 to 5819 (about 1278 px); 1278 x 31,309 = 40,013 µm, so the mosaic calibration is in **micrometres per pixel** |
| Derived: mosaic 0 `Map` rectangle (Left 740, Top 89, Right 5181, Bottom 621) | 4442 x 533 px = about 139 mm x 16.7 mm; it does not match the BCF footprint (a different region, probably an earlier planning step; unverified) |

The exploratory full-file run (14 elements, 8.4 s streaming) was run in a scratch script and not stored; the
test uses 40 lines and 4 elements to stay fast. The "tabulated K-alpha energies" (Ca 3.692, Mn 5.899, Fe 6.404, Zn 8.638 keV) are
from general knowledge, not from a file in this project. The comparison shows that the BCF energy axis puts the
Ca, Mn, Fe and Zn lines within about +-0.15 keV of the tabulated positions; it is **not** a precise validation
of the spectral calibration (TODO.md P1 stays open).

### 8.5 Timestamps (kept separate from identifiers)

| Source | Value |
|---|---|
| RTX outer header | 30.7.2026 15:58:47 |
| RTX mosaic 0 | 30.7.2026 9:38:26 |
| RTX mosaic 1 | 30.7.2026 9:55:51 |
| RTX `Mapdaten` | 30.7.2026 9:55:53 |
| BCF EDX acquisition | 2026-07-30 09:55:52 |

The order (mosaic 0, then mosaic 1, then the BCF acquisition, all within 18 minutes) is a fact; what each time
records and the time zone are unknown.

### 8.6 Relationship: what the RTX has that the BCF lacks (Confirmed content, Interpretation of role)

| Content | In the BCF | In the RTX |
|---|---|---|
| Video image, 240 x 1800 `uint16` | yes (`Video`) | yes, bit-identical (duplicate) |
| Spatial calibration of the map grid | yes | yes, equal |
| Spectrum image, 4096 channels, energy axis, detector, elevation angle, real time, sample name | yes | **no** |
| 14 element maps, 240 x 1800 `uint16`, with element-line names | not returned by RosettaSciIO | **yes** |
| Optical mosaic, 6000 x 948 x 3 `uint8`, 31.309 µm/px, and map footprint rectangle | not returned by RosettaSciIO; **the BCF header contains the same mosaic pixels** (section 9.1), without the footprint rectangle | **yes** |
| Composite display settings (colours, gamma, filter), scale-bar and text annotations, `HV: 50,0 kV` | no | yes |

No contradiction between the two files was found. The BCF `EDSDatabase` was later searched (section 9.1): its header holds no element maps and does hold more images than RosettaSciIO returns.

### 8.7 Implications for the NetCDF-4 data model (proposals only; the schema is decided in P2)

Everything below is a candidate, not a design; nothing is implemented.

- **Element maps (RTX):** 14 `uint16` (y, x) arrays on the same grid as the BCF cube. Candidate: ancillary
  variables on the shared `y`/`x` coordinates, named from the RTX `Description` text, with attributes recording
  the source (RTX), the payload plane index, and that the processing method is **unknown**. They must not be
  presented as counts equal to a window of the cube (8.4).
- **Video image:** duplicate of the BCF `Video`; store once, with provenance stating that both sources were
  verified equal, or store the BCF copy only.
- **Optical mosaic (RTX):** a different grid (6000 x 948, 31.309 µm) that needs its own dimensions and
  coordinates, not the map's. Candidate: a separate group, 3 `uint8` planes, plus the map footprint rectangle
  in mosaic pixels as attributes. Whether the mosaic axes are parallel to the map axes, and the origin, are
  unverified (open question 11), so no georeferencing may be claimed yet.
- **Plane statistics, `LineCounter`, `Valid`, palettes, display settings, overlays:** presentation or opaque
  data. Candidate: preserve verbatim as text attributes or an ancillary text variable, not as scientific
  variables. `HV: 50,0 kV` is annotation text and supports, but does not replace, reading the beam energy from
  the BCF header.
- **Calibration:** both sources agree; keep one set of coordinate values and record both sources.
- **Provenance:** the outer RTX `Date`/`Time`, the three image times and their unknown meaning, the payload
  SHA-256, and the decode parameters belong in the provenance attributes.
- **Bounded memory:** the largest RTX unit is one plane (5.7 MB here), so an RTX stage can run beside the BCF
  band stream at negligible cost. This is verified on this file only.

### 8.8 Limitations

- One RTX file. Other RTX files may hold other `TRT*` types, spectra, ROI tables or compressed planes; the
  module reads unknown elements only as counts and names, and fails loudly on an unsupported `RTCompression`,
  invalid Base64, a truncated or over-long zlib stream and non-well-formed XML (all covered by synthetic tests).
- Little-endian `uint16` for `ItemSize` 2 is an assumption that the equality with the BCF `Video` dataset
  supports for that plane; `ItemSize` values other than 1 and 2 are recorded, not decoded.
- The synthetic fixtures mimic only the element layout observed here.
- Payload offsets are decompressed-stream byte offsets reported by expat; they are valid for this file's
  SHA-256 only.

## 9. Unified BCF + RTX converter: new findings and measured results

Work of 2026-09-29/30 in the `microxrf_to_netcdf` environment (Windows 11, 31.5 GB RAM, Python 3.12.14). The converter is the
package `microxrf_to_netcdf` 0.1.0 (schema in [NETCDF_SCHEMA.md](NETCDF_SCHEMA.md)); the deliverable is
`output/GRF17A_9-29cm_slab3.nc`. The original files, the installed RosettaSciIO, Hyper Fusion and HSI Toolbox were not
modified. Measurements are for this one acquisition and this machine. Saved evidence is in [reference/](reference/):
`conversion_report.json`, `convert_log.txt`, `convert_stdout.json`, `validation_independent.json`,
`benchmark_memory.txt`, `benchmark_chunks.txt`.

### 9.1 The BCF header holds images that RosettaSciIO does not return (Confirmed)

`EDSDatabase/HeaderData` (34,006,136 bytes) is an XML document of the same `TRT*` family as the RTX payload. Decoded
incrementally with the diagnostic decoder of `inspect_rtx.py` from the SFS blocks (0.6 s, 63,153 characters of residual XML
once plane text is removed) it holds **7** `TRTImageData` elements. RosettaSciIO 0.14.0 returns only the first
(`Video`); `_set_images` reads the overview images only when `FileVersion == 2`, and this file has `FileVersion` 1.

| # | Name | Size (W x H) | ItemSize | Planes | Calibration text | Time | Content |
|---|---|---|---|---|---|---|---|
| 0 | (none) | 1800 x 240 | 2 | 1 (`Video`) | 0,0 | 9:55:52 | equals the RTX `Video 1` plane (SHA-256 6841bf7c…) |
| 1 | `Counter` | 1800 x 240 | 2 | **0** | 0,0 | 9:55:52 | no planes |
| 2 | `PixelTimes` | 1800 x 240 | **4** | 1 | 0,0 | 9:55:52 | uint32; minimum 14,500, maximum 31,400, 32 distinct values, sum 12,955,261,600 |
| 3 | `Default` | 6000 x 948 | 1 | 3 | 31,308984675955 | 9:55:51 | **byte-identical to both RTX mosaics** (all three planes) |
| 4 | `Image_0` | 1024 x 768 | 1 | 3 | 14,7528011327442 | 9:55:51 | other overview image |
| 5 | `Image_1` | 1024 x 768 | 1 | 3 | 1,44248358289655 / 1,44001779899416 | 9:55:51 | other overview image |
| 6 | `Image_2` | 752 x 480 | 1 | 3 | 1,08936170212766 / 1,28 | 9:55:51 | other overview image |

- The corresponding statement in section 3 (the BCF "holds only its own Video image") and the "no" entry for the BCF in
  the mosaic row of section 8.6 hold for **what RosettaSciIO returns**, not for what the BCF contains. The BCF
  contains the mosaic (identical pixels to the RTX) and four more images.
- The sum of `PixelTimes` (12,955,261,600) is close to the header real time of 12,960 s expressed in microseconds
  (12,960,000,000; difference 0.04 %). This is a numerical observation. It suggests per-pixel times in microseconds
  but does **not** verify the unit; the file assigns no unit.
- The relationship of `Image_0`, `Image_1` and `Image_2` to the acquisition grid is not established (no footprint and no
  stage coordinates were found). Their calibrations carry no unit.
- No element map or other image exists in the BCF header beyond these seven. The element list that
  RosettaSciIO reads from the header (`Sample.elements`) is empty for this file. This answers, for this file, the "BCF `EDSDatabase`
  searched for element maps" part of TODO P1 follow-up D.
- All seven images, the residual header XML (every element except the image-plane text) and their roles are stored in the
  output file (NETCDF_SCHEMA.md sections 5 to 7). Every header plane was compared again with the diagnostic decoder
  (`14 of 14 planes`).

### 9.2 Header sum spectrum versus the decoded pixels (Confirmed; cause unknown)

The header spectrum (`Channels`, 4,096 values, total 680,903,120) is **not** the sum of the decoded pixel spectra
(677,896,870): the header is larger by 3,006,250 counts (0.442 %), in 3,873 channels, and the pixel sum never exceeds the
header in any channel (largest difference in one channel 19,559; largest relative difference 31 % in a sparse channel).
It is not a wrap-around (a wrap would differ by multiples of 65,536 at the strongest channels). Interpretation, unverified:
the header may include events that were not assigned to a pixel. Consequence: the header sum is used only as an upper
bound (the converter refuses a decode whose channel sums exceed it) and is stored beside the derived pixel sum.

### 9.3 Mosaic registration (Confirmed for RTX mosaic instance 1)

Method: area-average the grey mosaic crop inside the inclusive `Map` rectangle to 240 x 1800 and correlate it with the BCF
`Video` image (Pearson r).

| Comparison | r |
|---|---|
| identity orientation, rectangle (192, 89, 5941, 855) | **0.9956** |
| vertical flip / horizontal flip / rotation by 180° | 0.4165 / 0.3782 / 0.2412 |
| best shift within ±6 mosaic pixels (step 2) | (0, 0), r 0.9956 |
| rectangle read as exclusive (192, 89, 5942, 856) / shifted by one pixel (191, 88, …) | 0.9928 / 0.9869 |

So mosaic and acquisition axes are parallel with the same direction, and the `Map` rectangle is inclusive, for that
instance. RTX image 0 has a `Map` rectangle of another size (4442 x 533 px) and is recorded as not verified; the BCF `Default`
instance has no rectangle. Pixel-corner versus pixel-centre coordinates and the rotation of the supplied `Video Mosaic.png`
(90° clockwise relative to the RTX mosaic, r 0.95 to 0.97 per channel, confirmed again) remain unexplained.
Each mosaic image in the RTX carries **two** `Map` rectangle records: a degenerate one (all four values 0) and the drawn
footprint; the converter ignores degenerate records and requires exactly one real one.

### 9.4 Counts precision: uint8 is safe for this acquisition; the uint16 wrap hazard is real (Confirmed)

- An exact (uint32) decode of the whole cube gives a maximum per-pixel, per-channel count of **149** and 215,280,205
  non-zero values out of 1,769,472,000 (the rest are valid zeros). `uint8` therefore holds every count of this file; the
  output dtype is chosen from that measured maximum, never from RosettaSciIO's heuristic.
- Decoding every band as `uint16` and as `uint32` gives identical arrays for all 60 bands of 4 lines (`python -m microxrf_to_netcdf validate
  --strict`), so no count wrapped.
- The compiled decoder adds in place in the chosen dtype. A synthetic count of 70,000 decoded as `uint16` becomes **4,464**,
  a value no "maximum near the ceiling" test would flag (a first version of the guard that widened bands with maximum at
  least 32,768 failed this test). The converter therefore scans the stream once in `uint32` (exact) and decodes the write
  and read-back passes in `uint16` only if the exact maximum fits, requiring the per-channel sums, total, non-zero count and
  maximum of the two decodes to be identical; otherwise it decodes in `uint32`.

### 9.5 Memory (Measured)

| Quantity | Value |
|---|---|
| Converter peak process working set, complete conversion including scan, write and read-back passes, lines per band 1 / 4 / 8 / 20 | 312.1 / 378.6 / 642.7 / 1324.4 MiB (chunks `(1,60,4096)` for one line, `(4,60,4096)` otherwise; `reference/benchmark_memory.txt`) |
| Final deliverable run (default: 4 lines per band, chunks `(2,30,4096)`) | **453.4 MiB** peak (earlier identical runs 408 to 465 MiB: allocator variation) |
| Stage peaks in profiled runs (MiB, cumulative peak so far) | header parse 190; scan pass 235; RTX + registration 226; write pass 323 to 360; read-back 408 to 448 |
| Preflight estimate for the final run | 495.3 MiB (an **estimate**; it stays above every measured peak of 408 to 465 MiB) |
| Full cube for comparison | 1,769,472,000 bytes as uint8; the peak is about 27 % of it |
| Retention | after the band loops about 180 MiB stays resident (HDF5 chunk cache and allocator); it is included in the estimate |
| Reading the counts back through netCDF4 | a fresh process that only reads all bands peaks at 177 to 277 MiB (1 to 8 lines per call), against a 41 MiB baseline; a library-side cost, not converter-specific, and it is why the read-back stage is the highest peak |
| `import numpy` on this machine | commits **about 685 MiB** (OpenBLAS thread buffers) versus about 17 MiB with `OPENBLAS_NUM_THREADS=1`; the package defaults it to 1 because the converter does no linear algebra |
| Hard cap test (Windows Job Object, committed process memory) with a synthetic 786 MB (48 x 1000 x 16384, uint8) cube streamed on the fly | completes at caps of 500 and 350 MiB (peak working set 337 MiB); fails **cleanly** at 250 MiB with a `MemoryError` wrapped in `ConversionError` and no file left; a 40 MiB control fails |

The claims "bounded memory" and "larger than the memory limit" are therefore demonstrated for a **process-memory cap**,
not for a dataset larger than the physical RAM of this machine (31.5 GB): no such dataset was available.

### 9.6 Chunk shapes and compression (Measured; `benchmark_conversion.py`, `reference/benchmark_chunks.txt`)

Four `counts` chunk shapes converted with the same code (zlib 4, shuffle); read timings on a fresh open of the file
(cold HDF5 chunk cache, warm operating-system cache; indicative). Table and reasoning are in
[NETCDF_SCHEMA.md](NETCDF_SCHEMA.md) section 9. Summary: `2x30x4096` was chosen as default (pixel spectrum 0.61 ms, best
16x16 window 4.32 ms, 212.8 MiB, 7,200 chunks); `1x16x4096` has the fastest spectra (0.27 ms) but 27,120 chunks; `8x64x256`
is ten times faster for one energy-channel image (308 ms against 3,295 ms) and smallest (198.2 MiB) but 7 to 16 times
slower for a pixel spectrum. An earlier version of the read benchmark took the minimum of three repeats on an open file
and measured the HDF5 chunk cache (0.07 ms per spectrum); it was replaced by a fresh open per trial.
Compression: the 1.81 GB logical dataset is 226,985,418 bytes (216.5 MiB), a ratio of 7.97; most counts are zero
(215,280,205 non-zero values out of 1.77 billion).

### 9.7 Final conversion of the real acquisition (Measured)

| Item | Value |
|---|---|
| Command | `python -m microxrf_to_netcdf convert --bcf data/GRF17A_9-29cm_slab3_Elemental_map.bcf --rtx data/GRF17A_9-29cm_slab3_Elemental_map.rtx --out output/GRF17A_9-29cm_slab3.nc --overwrite` |
| Output | `output/GRF17A_9-29cm_slab3.nc`, **226,985,418 bytes**, SHA-256 `e6a723be8a6d4c95d0cd3601a7fab06d391106694a236ba706b7f65c3c90e3ae` (sidecars `.sha256`, `.report.json`) |
| Time | 51.3 s (scan 7.5 s, write pass 28.8 s, the rest RTX/header passes, read-back and finalization) |
| Peak process working set | 453.4 MiB |
| Inputs | BCF 617,320,728 bytes, SHA-256 `5ec18488fcd40db55392e14af570d42568ff4517fc082bec49e519d521253d3a`; RTX 43,723,950 bytes, SHA-256 `c1be8801321d4840b285799f09c957e6cdf407c0503ac6e1c0c280190de0ff82` (decompressed payload `832f6f4c…` as in 8.1); unchanged (size and modification time compared) |
| Compression | zlib level 4 with the shuffle filter on every data variable |
| Actual chunks | `counts` (2, 30, 4096) uint8; `element_maps` (1, 120, 900) uint16; `video` and `pixel_times` (120, 900); `mosaic/pixels` (1, 256, 1024); `overview/Image_0/pixels` and `Image_1` (1, 256, 1024); `Image_2` (1, 256, 752) |
| Estimated logical size (estimate) | 1,808,352,420 bytes; disk check: 382 GB free, worst case 1.84 GB plus margin = 2.38 GB required |
| dtype decision | `counts` as uint8 because the exact maximum is 149; decoding uint32 (scan) then uint16 (write and read-back) |
| Validation recorded in the file | 29 checks, all passed (`/metadata/validation_json`) |

Independent validation of this file (`python -m microxrf_to_netcdf validate --strict --data-dir data`, 21 checks, all passed,
`reference/validation_independent.json`): lazy reopening (xarray and Dask, dtype `uint8`, no `_FillValue`, chunks
`(2, 30, 4096)`); recorded source hashes and sizes equal the input files; **counts binned 8 x 8 over all 4,096 channels
equal the public `BCF_reader.parse_hypermap(downsample=8)`**; **the first 400 channels of every one of the 432,000 pixels
equal the unmodified decoder on the real stream**; the stored pixel sum spectrum equals the sum of the counts read back; every
RTX plane (21) and every BCF header plane (14) equals the diagnostic decoder; the mosaic pixels rotated 90° clockwise
correlate with `Video Mosaic.png` (r 0.970, 0.958, 0.953); the Ca and Fe maps correlate with the blue and red channels of
`CaFe.png` (r 0.982, 0.911); uint16 and uint32 decodes are identical.
Additional large-file tests compare all 240 lines x 4,096 channels with the throw-away framer of
`probe_bcf_streaming.py`, representative pixel spectra, a 3-D window and one energy-channel image after lazy reopening,
the element maps, video and mosaic with the diagnostic RTX decoder, and the calibration with the public reader.

### 9.8 What cannot yet be independently validated

1. Per-pixel counts of **channels 400 to 4,095** against the public decoder: the public path cannot return them without
   materializing (nearly) the whole cube. They are validated (a) in aggregate over all channels by the 8 x 8 binned
   comparison, (b) against a second framer (the prototype) that shares the same compiled decoder, and (c) by the uint16 versus
   uint32 decodes. An error that cancels inside an 8 x 8 block in a high channel and is common to both framers would go
   unnoticed.
2. The **processing and units of the RTX element maps** (unknown; they are not window sums, section 8.4).
3. Spectral calibration and pixel size against an **external reference** (only the internal agreement of BCF and RTX and a
   +-0.15 keV consistency of four K-alpha windows exist).
4. **Units** of beam energy, elevation angle, `PixelTimes` and the overview calibrations; the meaning of `Valid`, both
   `LineCounter` fields and every timestamp; the 0.44 % header excess; corner-versus-centre coordinates; the footprint of
   mosaic image 0; the relationship of the three other overview images to the acquisition grid.
5. Pixel encodings `flag` 0, `flag` 1 and `n_of_pulses > 0` (synthetic streams only), zlib-compressed SFS containers, BCF
   files with several hypermaps or `FileVersion` 2, other RTX layouts, other detectors: none is exercised by the one
   real acquisition. Little-endian 32-bit `PixelTimes` is an assumption.
6. Behaviour on an acquisition larger than the physical RAM (only a process-memory cap was demonstrated).

### 9.9 Environment and tests

- Added to `microxrf_to_netcdf` with pip (verified SSL, no environment deleted or base touched): netCDF4 1.7.4 (HDF5 1.14.6,
  netCDF-C 4.9.3), xarray 2026.7.0, pandas 3.0.6, cftime 1.6.6. `pip check`: no broken requirements. NumPy stays 2.5.3, Dask
  2026.8.0, RosettaSciIO 0.14.0.
- A **fresh virtual environment built from `requirements.txt` alone** installed the same versions and passed `pip check`; the
  fast tests ran in it (TODO P0). The from-scratch `conda env create` still cannot be tested: conda-forge fails on the
  development machine with an SSL certificate verification error (verification was not disabled).
- netCDF4 1.7.4 with NumPy 2.5 emits `DeprecationWarning: Setting the shape on a NumPy array has been deprecated`; it is
  harmless today and is filtered inside the writer.
- Test suite after this work: see section 6.
