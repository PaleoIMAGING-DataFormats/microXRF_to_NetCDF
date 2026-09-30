# microXRF to NetCDF — NetCDF-4 schema 1.0.0

One paired BCF/RTX acquisition produces **one** NetCDF-4 (HDF5) file. This document is the contract of
`microxrf_to_netcdf` version 0.1.0, schema version **1.0.0** (`microxrf_to_netcdf_schema_version` in the file). Different acquisitions
stay separate files. Every fact below was observed on the one available acquisition
(`GRF17A_9-29cm_slab3_Elemental_map.bcf` + `.rtx`, see [FINDINGS.md](FINDINGS.md)); a claim that rests on one file is
marked. A test (`tests/test_header_images.py::test_every_written_variable_is_documented_in_the_schema`) fails if the
converter writes a variable that this document does not name.

## 1. Principles

- The file is self-contained: spectra, element maps, optical images, calibration and metadata of both sources.
- Every piece of content keeps its **source identity** (BCF or RTX, and where in it), **calibration**, **units** and
  **provenance**. A unit that the source does not state is not assigned; the attribute `units_status` says why.
- Raw values and derived values are named differently (`energy_calibration` holds the raw parameters, `energy` is
  the derived coordinate).
- Content that cannot be interpreted is preserved without interpretation (residual XML, verbatim texts).
- The file has several groups and several spatial grids. There is **no limit** on its size. It is written
  sequentially with bounded memory and is meant to be read lazily (spatial and spectral subsets).
- No `_FillValue` is set on any count or image variable. A count of `0` is data (integer zero), never missing.
  Completeness is proven by explicit checks and by `conversion_status`.

## 2. Layout

```
/                                   global identity, versions, source names, sizes, SHA-256, history
├── acquisition/                    the BCF grid (y, x) and the spectral axis (energy)
│   dims  y=240  x=1800  energy=4096  element=14
│   y, x, channel, energy                       coordinates
│   energy_calibration, spatial_calibration     RAW calibration parameters (scalar variables holding attributes)
│   counts            (y, x, energy)  uint8     EDS counts (BCF)
│   video             (y, x)          uint16    video image, stored once (BCF and RTX are bit-identical)
│   pixel_times       (y, x)          uint32    BCF-header image "PixelTimes" (units unknown)
│   element_maps      (element, y, x) uint16    14 RTX element maps
│   element, element_line, element_label, element_map_source_plane, element_map_sha256, element_map_valid_flag
│   sum_spectrum_header (energy) uint64 · sum_spectrum_from_pixels (energy) int64 · line_counter (y) int32
│   instrument                                  instrument/acquisition metadata (scalar variable holding attributes)
├── mosaic/                         the optical mosaic on its own grid
│   dims  mosaic_y=948  mosaic_x=6000  plane=3
│   mosaic_y, mosaic_x                          coordinates
│   pixels            (plane, mosaic_y, mosaic_x) uint8      stored ONCE (pixels_2, ... only if pixels differ)
│   instance_rtx_0/  instance_rtx_1/  instance_bcf_3_Default/    one subgroup per source instance
│       annotations_json, map_footprint         per-instance annotations and footprint
├── overview/                       further images stored in the BCF header (own grids; RosettaSciIO does not return them)
│   Image_0/ Image_1/ Image_2/      dims plane=3 row column; variables row, column, pixels, annotations_json
└── metadata/                       source metadata, provenance, validation, open questions (text variables)
    bcf_original_metadata_json, bcf_header_residual_xml, bcf_image_table_json,
    rtx_outer_prefix_text, rtx_outer_tail_text, rtx_payload_residual_xml, rtx_image_table_json,
    timestamps_json, mapdaten_annotations_json, mapdaten_display_settings_json,
    provenance_json, validation_json, unresolved_properties_json
```

Dimension sizes are those of the reference acquisition; they follow the source. The group and instance names are
generated (`instance_<source>_<index>[_<name>]`, `<name>` from the BCF header image name).

## 3. Global attributes

| Attribute | Content |
|---|---|
| `title`, `structure`, `Conventions` | text; `Conventions` states that only the units and coordinate-variable conventions of CF-1.8 are followed (the file is not CF-compliant as a whole) |
| `microxrf_to_netcdf_schema_version`, `microxrf_to_netcdf_version` | `1.0.0`, `0.1.0` |
| `acquisition_id` | common stem of the source file names |
| `conversion_status`, `conversion_status_note` | `incomplete` while writing; `complete` only after every band and plane was written and verified. A file whose status is not `complete` must not be used |
| `source_bcf_filename`, `source_bcf_size_bytes`, `source_bcf_sha256` | the BCF input |
| `source_rtx_filename`, `source_rtx_size_bytes`, `source_rtx_sha256` | the RTX input |
| `created_utc`, `history` | creation time and one history line |

The SHA-256 of the output file itself cannot be stored inside it; it is written to `<file>.sha256` and to
`<file>.report.json` next to it.

## 4. `/acquisition`

Coordinates and calibration:

| Variable | Dims | dtype | Units | Source and meaning |
|---|---|---|---|---|
| `y` | y | float64 | `um` | `row_index x pixel_size_y_raw`; offset 0 as reported by the source; corner or centre of the pixel is **not verified**; attribute `source_units_label` keeps the reader's original label |
| `x` | x | float64 | `um` | `column_index x pixel_size_x_raw` |
| `channel` | energy | int32 | – | raw zero-based channel index |
| `energy` | energy | float64 | `keV` | **derived**: `calib_abs + calib_lin x channel` |
| `energy_calibration` | – | int32 (value 0) | `keV` | **raw** BCF spectrum-header values as attributes: `calib_abs_raw`, `calib_lin_raw`, `sigma_abs_raw`, `sigma_lin_raw` |
| `spatial_calibration` | – | int32 (value 0) | `um` | **raw** `pixel_size_y_raw`, `pixel_size_x_raw`, the RTX calibration texts, `source_units_label`; states that no external reference exists |

The scalar variables that only carry attributes hold the value 0; their data are not meaningful.

Data:

| Variable | Dims | dtype | Source | Notes |
|---|---|---|---|---|
| `counts` | y, x, energy | uint8 (auto) | BCF `EDSDatabase/SpectrumData0` through the unmodified RosettaSciIO decoder | dtype is the smallest of uint8/uint16/uint32 that holds the exact maximum count (`dtype_selection` attribute); units `counts`; no `_FillValue` |
| `video` | y, x | uint16 | BCF `Video` = RTX `Mapdaten` plane 0 (`also_in` attribute) | stored once because the two are bit-identical (validated by SHA-256 and element-wise); also equal to the BCF header image of the same name |
| `pixel_times` | y, x | uint32 | BCF header image `PixelTimes` (not returned by RosettaSciIO) | no `units` attribute: unknown. `data_min`, `data_max`, `data_sum` are derived while writing. The sum is close to the header real time expressed in microseconds, which suggests microseconds; **not verified** |
| `element_maps` | element, y, x | uint16 | RTX `Mapdaten` planes 1..14 | `processing_method = "unknown"`, `units_status = "unknown"`, no `units` attribute, not raw energy-window sums (FINDINGS 8.4) |
| `element` | element | str | text before `-` of the RTX plane description (`Ca` of `Ca-KA`) | plain text split |
| `element_line` | element | str | text after `-` (`KA`) | the meaning of the suffixes is not verified |
| `element_label` | element | str | RTX plane description, verbatim | |
| `element_map_source_plane` | element | int32 | plane index in the RTX map image | |
| `element_map_sha256` | element | str | SHA-256 of the decoded plane bytes | |
| `element_map_valid_flag` | element | str | RTX plane field `Valid`, verbatim | meaning unknown |
| `sum_spectrum_header` | energy | uint64 | BCF spectrum header `Channels` | not equal to the pixel sums on the reference file (0.442 % more counts; cause unknown) |
| `sum_spectrum_from_pixels` | energy | int64 | **derived** while writing | sum of `counts` over all pixels |
| `line_counter` | y | int32 | BCF header `LineCounter` | meaning not verified |
| `instrument` | – | int32 (value 0) | BCF header | attributes: `sample_name`, `detector_type`, `detector_technology`, `detector_serial`, `bcf_date`, `bcf_time`, `beam_energy_recorded`, `beam_energy_units_status`, `elevation_angle_recorded`, `elevation_angle_units_status`, `real_time_seconds_from_rosettasciio`, `real_time_note`. **Units of beam energy and elevation angle are not stored in the BCF**: no `units` attribute is set; the `..._units_status` attributes give the supporting evidence (RTX label `HV: 50,0 kV`; RosettaSciIO) and say it is not independently verified. The complete header is in `/metadata/bcf_original_metadata_json` |

**Shared grid rule.** `counts`, `video`, `pixel_times` and `element_maps` share `y` and `x` only because the
correspondence was validated before writing: the RTX map image has the BCF height and width; its X/Y calibration
equals the BCF pixel size (relative tolerance 1e-9); one of its planes is bit-identical to the BCF video; the header
image `PixelTimes` comes from the same BCF and has the same size. The validation record is the group attribute
`grid_validation` and `/metadata/validation_json`. Any mismatch stops the conversion (no silent merging).

## 5. `/mosaic`

The optical mosaic keeps its native resolution and its own dimensions.

| Variable | Dims | dtype | Notes |
|---|---|---|---|
| `mosaic_y`, `mosaic_x` | mosaic_y / mosaic_x | float64, `um` | `pixel_index x calibration` (31.309 per pixel on the reference file); the calibration has no unit next to it in the sources: micrometres per pixel by derivation from the `40000 µm` scale-bar annotation (FINDINGS 8.4) |
| `pixels` | plane, mosaic_y, mosaic_x | uint8 | 3 planes; correlate with R, G, B of the supplied `Video Mosaic.png` (r 0.95 to 0.97; that PNG is this array rotated 90° clockwise). Stored once for every instance whose per-plane SHA-256 are equal (`stored_once_for_instances`). If an instance differs, it gets its own `pixels_2`, `pixels_3`, … |

**Instances.** On the reference file three source images have the same pixels: RTX image 0 (time 9:38:26), RTX image 1
(9:55:51) and the BCF-header image `Default` (9:55:51). Equal pixels do **not** make them interchangeable, so each has
a subgroup (`instance_rtx_0`, `instance_rtx_1`, `instance_bcf_3_Default`) with its own source, timestamp, calibration
text, plane hashes, `LineCounter` and `Valid` values, and:

| Variable | Notes |
|---|---|
| `annotations_json` | overlay elements as parsed (type, name, text such as `40000 µm` or `Mosaik`, rectangle, position); the complete source is in the residual XML |
| `map_footprint` | scalar variable whose attributes give the `Map` rectangle (`rect_left`, `rect_top`, `rect_right`, `rect_bottom`, inclusive mosaic pixels), its size in pixels and micrometres against the acquisition extent, `size_matches`, and `registration_status` / `relationship_to_acquisition_grid` / `registration_json` |

**Verified spatial relationship (reference file, instance `rtx_1`).** The BCF video image and the mosaic crop inside
the inclusive `Map` rectangle (192, 89, 5941, 855), area-averaged to 240 × 1800, have Pearson r = **0.9956** in identity
orientation and at most 0.42 for the three other orientations; the best shift within ±6 mosaic pixels is (0, 0). So the
mosaic axes are parallel to the acquisition axes with the same direction, and the rectangle is inclusive. Instance
`rtx_0` has a rectangle of another size (139 mm × 16.7 mm) and is recorded as `not_verified`; the BCF `Default` instance
has no rectangle. Corner-versus-centre conventions are not verified. `registration_status` is `verified` only when the
identity orientation correlates at r ≥ 0.90 and beats every other orientation by 0.30; otherwise it says
`not_verified` with the numbers.

## 6. `/overview`

Images found in the BCF header that RosettaSciIO does not return and that are not the mosaic: on the reference file
`Image_0` (1024 × 768, 14.75 per pixel), `Image_1` (1024 × 768, 1.44 per pixel) and `Image_2` (752 × 480, 1.09 / 1.28 per
pixel), all 8-bit with 3 planes and time 9:55:51. Each is a subgroup with dims `plane`, `row`, `column`, coordinates
`row`, `column` (`index x calibration`, **no unit**: `units_status` says the header carries none), `pixels` (uint8) and
`annotations_json`. Their spatial relationship to the acquisition grid is **not established** (no footprint and no stage
coordinates were found).

## 7. `/metadata`

All are scalar Unicode string variables (`str`, UTF-8 in the file).

| Variable | Content |
|---|---|
| `bcf_original_metadata_json` | every BCF header dictionary RosettaSciIO exposes (Hardware, Detector, Analysis, Spectrum, DSP Configuration, Stage, Microscope, HyperHeader), values verbatim |
| `bcf_header_residual_xml` | the BCF header XML (34 MB on the reference file) with the text of every image-plane `Data` element removed (63 KB on the reference file); everything else is preserved |
| `bcf_image_table_json` | every image found in the BCF header with geometry, timestamp, planes, plane hashes and **where it was stored** (including `Counter`, which has no planes) |
| `rtx_outer_prefix_text`, `rtx_outer_tail_text` | the outer RTX text before `<RTData>` and after `</RTData>`, verbatim |
| `rtx_payload_residual_xml` | the decompressed RTX payload (Windows-1252 decoded to Unicode) with the text of every image-plane `Data` element removed (132 KB on the reference file); unknown fields, `LineCounter`, `Valid`, overlays, palettes and display settings are preserved |
| `rtx_image_table_json` | every RTX image and plane with geometry, calibration text, timestamp, statistics and hashes |
| `timestamps_json` | every timestamp of both sources, verbatim (`iso_naive` added where the text parses); meaning and time zone unknown |
| `mapdaten_annotations_json`, `mapdaten_display_settings_json` | overlays (including `HV: 50,0 kV`) and display settings (colours, gamma, filter) of the RTX map image: presentation, not data |
| `provenance_json` | software versions (microXRF to NetCDF under the key `microxrf_to_netcdf`, RosettaSciIO, netCDF4, HDF5, NumPy, xarray, Dask), all conversion parameters, source names, sizes and hashes, the preflight **estimates** (marked as estimates), the actual chunking and filters as stored, the decoding precision of each pass, scan statistics, measured time and peak memory up to finalization, registration results |
| `validation_json` | every validation check with its status and detail |
| `unresolved_properties_json` | the scientific properties that cannot yet be independently validated |

## 8. Conventions

- **Axis order** `(y, x, energy)`: `y` is the BCF `height` axis (scan line, first line first), `x` the `width` axis.
- **Origin and orientation.** Coordinates start at 0 at the first pixel (source offset 0). `y` and `x` increase with the
  index. The orientation relative to the sample is not defined by the sources except through the verified mosaic
  relationship in section 5.
- **Units.** `um` for spatial coordinates with a calibration in the source; `keV` for the energy coordinate; `counts` for
  the counts. No units are assigned to the element maps, `pixel_times`, the overview images, beam energy or elevation
  angle. The text attribute `units_status` explains each case.
- **Missing data.** None is defined: there is no mask and no `_FillValue`. The RTX plane field `Valid` is copied
  verbatim and is not interpreted.
- **Endianness.** 16- and 32-bit source items are read as little-endian (supported for the 16-bit video by its bit
  identity with the BCF video; assumed for `pixel_times`).
- **Strings** are UTF-8; RTX and BCF header texts (Windows-1252) are decoded to Unicode.

## 9. Chunking and compression

HDF5 chunking and zlib are configurable (`ConversionConfig`, CLI `--chunks`, `--complevel`, `--no-shuffle`). Defaults:
zlib level 4 with the shuffle filter; chunks are clipped to the dimension sizes.

| Variable | Default chunks | Reason |
|---|---|---|
| `counts` | `(2, 30, 4096)` (245 KB uncompressed as uint8) | full spectrum per chunk: a pixel spectrum touches one chunk, a small spatial window touches a few |
| `element_maps` | `(1, 120, 900)` | one map per chunk row: a map is read without the others; 216 KB |
| `video`, `pixel_times` | `(120, 900)` | four chunks per image |
| `mosaic/pixels`, `overview/*/pixels` | `(1, 256, 1024)` | one plane per chunk row; regions of the mosaic read independently |

Measured on the reference file (candidates converted with the same code, read on a warm operating-system cache and a
cold HDF5 chunk cache; machine-specific, indicative; `benchmark_conversion.py`, saved in `reference/`):

| Chunks (y, x, energy) | Size MiB | Pixel spectrum ms | 16×16 window, all energy ms | 8×100 window × 400 channels ms | One energy channel, whole image ms |
|---|---|---|---|---|---|
| 4×60×4096 | 212.4 | 1.92 | 7.77 | 17.65 | 3295 |
| 1×16×4096 | 215.1 | **0.27** | 4.99 | 7.99 | 3495 |
| 8×64×256 | **198.2** | 4.30 | 13.43 | **4.73** | **308** |
| **2×30×4096 (default)** | 212.8 | 0.61 | **4.32** | 7.25 | 3294 |

Trade-offs. Full-energy chunks favour spectra and spatial windows but force a whole-image read for one energy
channel (about 3.3 s here, the whole compressed cube is decompressed); energy-split chunks such as `8 64 256` are ten
times faster for single-channel images and windows over few channels, and slightly smaller, but 7 to 16 times slower for
a pixel spectrum. `1×16×4096` has the fastest spectra but 27,120 chunks instead of 7,200, which grows the chunk index
for larger acquisitions. The default keeps the fast spectrum, the best small-window time and a moderate chunk count.
For work centred on energy-channel images, convert with `--chunks 8 64 256` (`--band-lines` must stay a multiple of
the chunk extent along `y`). Sequential writing: `--band-lines` (default 4) is a multiple of the `y` chunk extent, so
each chunk is written once from one band; the bands are converted and written one at a time.

## 10. Completion, atomic finalization and validation

The converter writes to `.<name>.partial` in the destination directory, verifies it, then renames it to the final name
with one atomic `os.replace`. A failure or interruption removes the partial file; a hard kill leaves only the
dot-prefixed partial file, never a file that looks finished, and the next run refuses to start until it is removed.
Before the rename the converter (a) compares every written band with a fresh decode, (b) compares every element map,
mosaic, overview and header image with the source by SHA-256, (c) checks coordinates against the raw calibration,
(d) reopens the file lazily with xarray and Dask. `python -m microxrf_to_netcdf validate` repeats the source-facing checks with
independent references (the public RosettaSciIO path, the diagnostic RTX decoder, the companion PNG files).

## 11. Combination rules (what goes where)

| Source element | Destination |
|---|---|
| BCF spectrum stream, all pixels and channels | `/acquisition/counts` |
| BCF spectrum header: `CalibAbs`, `CalibLin`, `SigmaAbs`, `SigmaLin` | `energy_calibration` (raw), `energy` (derived) |
| BCF `DX`, `DY` | `spatial_calibration` (raw), `y`, `x` |
| BCF `Channels` (header sum spectrum), `LineCounter` | `sum_spectrum_header`, `line_counter` |
| BCF header image `Video` | `video` (one copy for BCF and RTX) |
| BCF header image `PixelTimes` | `pixel_times` |
| BCF header image `Default` (6000 × 948) | `/mosaic` instance `instance_bcf_<i>_Default`; pixels shared with equal instances |
| BCF header images `Image_0..2` | `/overview/Image_k` |
| BCF header image `Counter` (no planes) | recorded in `bcf_image_table_json` and the residual XML only |
| BCF instrument, detector, stage, DSP, microscope, date, time | `instrument` attributes (selected) and `bcf_original_metadata_json` (all) |
| Every other element of the BCF header | `bcf_header_residual_xml` |
| RTX `Mapdaten` plane 0 (`Video 1`) | validated equal to the BCF video; not stored again |
| RTX `Mapdaten` planes 1..14 | `element_maps` with `element`, `element_line`, `element_label` |
| RTX `Video Mosaic` images 0 and 1 | `/mosaic` instances `instance_rtx_0`, `instance_rtx_1` |
| RTX overlays, `Map` rectangles, texts, display settings | instance `annotations_json`, `map_footprint`, `mapdaten_*_json` |
| RTX `Valid`, `LineCounter`, palettes, timestamps, header fields, every other element | verbatim in `rtx_payload_residual_xml`, `rtx_outer_*`, `timestamps_json`, plane attributes |
| Anything the readers do not decode | refused (unsupported layout) or listed in `unresolved_properties_json`; never dropped silently |

## 12. Scientific properties not independently validated

The file carries them in `/metadata/unresolved_properties_json`; the authoritative list is in FINDINGS.md section 9.
In short: the processing and units of the element maps; the meaning of `Valid`, `LineCounter` and every timestamp;
why the BCF header sum spectrum exceeds the pixel sums by 0.44 %; the units of beam energy, elevation angle, `PixelTimes`
and the overview calibrations; spectral calibration and pixel size against an external reference; pixel corner versus
centre; the mosaic instance `rtx_0` footprint; and pixel encodings, container and file layouts the reference file does
not exercise.

## 13. Opening the file lazily

```python
import xarray as xr

path = "output/GRF17A_9-29cm_slab3.nc"

# Nothing is loaded: every variable is a Dask array backed by the HDF5 chunks.
acq = xr.open_dataset(path, group="acquisition", chunks={})
counts = acq["counts"]                           # (y, x, energy) uint8, chunks (2, 30, 4096)

spectrum = counts.isel(y=100, x=900).values      # one pixel spectrum: reads one chunk
window = counts.isel(y=slice(50, 58), x=slice(1000, 1100), energy=slice(0, 400)).compute()   # a 3-D subset
ca_map = acq["element_maps"].sel(element="Ca").values           # element map (units unknown)
roi = counts.sel(energy=slice(3.5, 3.9)).sum("energy")          # spectral window by keV, still lazy until .compute()

mosaic = xr.open_dataset(path, group="mosaic", chunks={})["pixels"]     # (plane, mosaic_y, mosaic_x) uint8
meta = xr.open_dataset(path, group="metadata")                          # text variables (JSON, XML)
```

Zero counts come back as integer `0` and the dtype stays `uint8`. With the netCDF4 library directly:
`netCDF4.Dataset(path)["acquisition"]["counts"][100, 900, :]`.

## 14. Versioning

`microxrf_to_netcdf_schema_version` follows semantic versioning: a change that renames or removes a variable or changes a
meaning increments the major version; an added optional variable or attribute increments the minor version. Readers
should check the attribute and `conversion_status`.
