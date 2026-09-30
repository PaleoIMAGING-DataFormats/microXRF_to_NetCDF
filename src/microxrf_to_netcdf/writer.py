"""
Authors: Andre L. Belem (https://github.com/andrebelem) and F.R.I.D.A.Y. (https://observatoriooceanografico.org/people/friday-bot)

Incremental NetCDF-4 writer implementing schema 1.0.0 (NETCDF_SCHEMA.md).

The file is created once with every dimension, coordinate, variable and attribute defined; the large variables
are then filled progressively: ``write_band`` stores one BCF band along ``y``, ``rtx_sink`` stores each RTX
plane as it is decoded. Nothing is buffered beyond the band or plane being written. No ``_FillValue`` is set on
any count or image variable (an EDS count of 0 is data, not missing; a ``_FillValue`` of 0 would turn zeros
into NaN in xarray). Completeness is therefore tracked explicitly (``lines_written``, ``planes_written``).

The writer never touches an existing final file; the caller passes a temporary path and finalizes it.
"""

from __future__ import annotations

import hashlib
import json
import warnings
from pathlib import Path
from typing import Any

import netCDF4
import numpy as np

from . import SCHEMA_VERSION, __version__
from .bcf import Band, BCFInfo, original_metadata_json
from .config import ConversionConfig
from .errors import ValidationError
from .model import Layout, iso_from_rtx, split_label
from .preflight import Plan
from .registration import footprint_summary, find_map_rectangle
from .rtx import RTXScan

ELEMENT_MAP_NOTE = ("Element maps are stored exactly as decoded from the RTX. Their processing method (ROI window, "
                    "net counts, deconvolution or other) and their units are UNKNOWN; they are NOT raw energy-window "
                    "sums of the counts variable (FINDINGS.md 8.4). No units attribute is assigned.")
_SHAPE_WARNING = "Setting the shape on a NumPy array has been deprecated"


def _attrs(obj: Any, **values: Any) -> None:
    for key, value in values.items():
        if value is None:
            continue
        if isinstance(value, bool):
            value = int(value)
        obj.setncattr(key, value)


def _text(group: netCDF4.Group, name: str, value: str, **attributes: Any) -> netCDF4.Variable:
    variable = group.createVariable(name, str, ())
    variable[...] = value
    _attrs(variable, **attributes)
    return variable


def _sha256_stream(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def file_sha256(path: Path) -> str:
    return _sha256_stream(Path(path))


class UnifiedWriter:
    """Creates the schema-1.0.0 file at ``path`` and fills it incrementally."""

    def __init__(self, path: Path, info: BCFInfo, scan: RTXScan, layout: Layout, plan: Plan,
                 config: ConversionConfig, identity: dict[str, Any]):
        self.path, self.info, self.scan, self.layout, self.plan, self.config = path, info, scan, layout, plan, config
        self.identity = identity
        self.counts_dtype = np.dtype(plan.counts_dtype)
        self.counts_max = int(np.iinfo(self.counts_dtype).max)
        self.dataset: netCDF4.Dataset | None = None
        self.lines_written: set[int] = set()
        self.planes_written: set[tuple[str, int, int]] = set()
        self._targets: dict[tuple[str, int, int], tuple[str, Any]] = {}
        self._aux_variables: dict[str, Any] = {}
        self._aux_stats: dict[str, dict[str, int]] = {}
        self._footprints: dict[int, tuple[Any, dict[str, int], dict[str, Any]]] = {}

    # ------------------------------------------------------------------------------------------
    # Creation
    # ------------------------------------------------------------------------------------------

    def _compression(self) -> dict[str, Any]:
        on = self.config.complevel > 0
        return {"zlib": on, "complevel": max(self.config.complevel, 1), "shuffle": bool(on and self.config.shuffle),
                "fill_value": False}

    def create(self) -> None:
        info, layout, config = self.info, self.layout, self.config
        dataset = self.dataset = netCDF4.Dataset(self.path, "w", format="NETCDF4")
        try:
            with warnings.catch_warnings():   # netCDF4 1.7.4 with NumPy 2.5: harmless shape DeprecationWarning
                warnings.filterwarnings("ignore", message=_SHAPE_WARNING, category=DeprecationWarning)
                self._create_root(dataset)
                self._create_acquisition(dataset.createGroup("acquisition"))
                self._create_mosaic(dataset.createGroup("mosaic"))
                self._create_overview(dataset.createGroup("overview"))
                self._create_metadata(dataset.createGroup("metadata"))
        except BaseException:
            dataset.close()
            self.dataset = None
            raise

    def _create_root(self, ds: netCDF4.Dataset) -> None:
        identity = self.identity
        _attrs(
            ds,
            title=f"microXRF to NetCDF unified acquisition: {identity['acquisition_id']}",
            Conventions=f"microxrf_to_netcdf schema {SCHEMA_VERSION}; only the units and coordinate-variable "
                        "conventions of CF-1.8 are followed",
            microxrf_to_netcdf_schema_version=SCHEMA_VERSION, microxrf_to_netcdf_version=__version__,
            acquisition_id=identity["acquisition_id"],
            conversion_status="incomplete",
            conversion_status_note="set to 'complete' only after every band and plane was written and verified; "
                                   "a file whose status is not 'complete' must not be used",
            source_bcf_filename=self.info.filename, source_bcf_size_bytes=int(self.info.size_bytes),
            source_bcf_sha256=identity["bcf_sha256"],
            source_rtx_filename=Path(self.scan.path).name, source_rtx_size_bytes=int(self.scan.size_bytes),
            source_rtx_sha256=identity["rtx_sha256"],
            created_utc=identity["created_utc"],
            history=identity["history"],
            structure="groups: /acquisition (BCF grid: counts, video, element maps, pixel times), /mosaic (optical "
                      "mosaic on its own grid, one subgroup per source instance), /overview (other overview images "
                      "stored in the BCF header, each on its own grid), /metadata (source metadata, provenance, "
                      "validation)",
        )

    def _coordinate(self, group, name, values, units, **attributes):
        variable = group.createVariable(name, "f8", (name,), fill_value=False)
        variable[:] = values
        _attrs(variable, units=units, **attributes)
        return variable

    def _create_acquisition(self, g: netCDF4.Group) -> None:
        info, layout, config, plan = self.info, self.layout, self.config, self.plan
        y, x, e = info.height, info.width, info.channels
        maps = layout.element_planes
        _attrs(g, description="Data on the BCF acquisition grid (y, x) with the spectral axis (energy).",
               grid_shared_by="counts, video and element_maps, after validation of the grid correspondence",
               grid_validation=json.dumps(layout.checks))
        for name, size in (("y", y), ("x", x), ("energy", e), ("element", len(maps))):
            g.createDimension(name, size)

        cy = self._coordinate(
            g, "y", np.arange(y, dtype="f8") * info.pixel_size_y, "um", long_name="acquisition row coordinate",
            axis="Y", positive_direction="increasing with row index (rows are scan lines, first line first)",
            coordinate_rule="coordinate = row_index * pixel_size_y_raw; source offset 0; whether the value refers "
                            "to the pixel corner or centre is not verified",
            source_units_label=info.units_label)
        cx = self._coordinate(
            g, "x", np.arange(x, dtype="f8") * info.pixel_size_x, "um", long_name="acquisition column coordinate",
            axis="X", coordinate_rule="coordinate = column_index * pixel_size_x_raw; source offset 0; corner or "
                                      "centre not verified", source_units_label=info.units_label)
        channel = g.createVariable("channel", "i4", ("energy",), fill_value=False)
        channel[:] = np.arange(e, dtype="i4")
        _attrs(channel, long_name="spectral channel index (raw, zero-based)")
        energy = self._coordinate(
            g, "energy", info.calib_abs + info.calib_lin * np.arange(e, dtype="f8"), "keV",
            long_name="channel energy (DERIVED coordinate)", axis="Z",
            derivation="energy_keV = calib_abs + calib_lin * channel, with the raw parameters stored in the "
                       "variable energy_calibration; not validated against an external source (FINDINGS.md)")

        calibration = g.createVariable("energy_calibration", "i4", ())
        calibration.assignValue(0)
        _attrs(calibration, long_name="RAW spectral calibration parameters as read from the BCF spectrum header",
               calib_abs_raw=float(info.calib_abs), calib_lin_raw=float(info.calib_lin),
               sigma_abs_raw=info.sigma_abs, sigma_lin_raw=info.sigma_lin, units="keV",
               units_note="unit of CalibAbs/CalibLin is reported as keV by RosettaSciIO",
               source="BCF EDSDatabase spectrum header (CalibAbs, CalibLin, SigmaAbs, SigmaLin)",
               meaning_of_sigma="not verified")
        spatial = g.createVariable("spatial_calibration", "i4", ())
        spatial.assignValue(0)
        _attrs(spatial, long_name="RAW pixel size parameters as read from the BCF header",
               pixel_size_y_raw=float(info.pixel_size_y), pixel_size_x_raw=float(info.pixel_size_x),
               units="um", source_units_label=info.units_label,
               rtx_x_calibration_raw_text=layout.grid.x_calibration_text,
               rtx_y_calibration_raw_text=layout.grid.y_calibration_text,
               agreement="BCF and RTX values agree to relative 1e-9 (validated)",
               external_validation="none: no independent reference for the pixel size has been used")

        counts = g.createVariable("counts", self.counts_dtype, ("y", "x", "energy"), chunksizes=plan.chunks["counts"], **self._compression())
        counts.set_auto_maskandscale(False)
        _attrs(counts, long_name="EDS X-ray counts per pixel and energy channel", units="counts",
               source="BCF EDSDatabase/SpectrumData0 decoded by the unmodified RosettaSciIO decoder",
               dtype_selection=plan.counts_dtype_reason,
               zero_counts="valid data (integer 0); the variable has no _FillValue",
               axis_order="(y, x, energy) = (BCF height, width, Energy)")
        self.counts = counts

        video = g.createVariable("video", "u2", ("y", "x"), chunksizes=plan.chunks["video"], **self._compression())
        video.set_auto_maskandscale(False)
        _attrs(video, long_name="video (optical) image of the mapped area", source="BCF Video dataset",
               also_in="RTX Mapdaten plane %d (%r), bit-identical to this array (sha256 %s)"
                       % (layout.video_plane.index, layout.video_plane.description, layout.video_plane.sha256),
               stored_once="the two sources are byte-identical, so a single copy is stored",
               units_note="no units are stored in either source")
        video[:] = info.video

        elements = g.createVariable("element", str, ("element",))
        lines = g.createVariable("element_line", str, ("element",))
        labels = g.createVariable("element_label", str, ("element",))
        split = [split_label(p.description or f"plane_{p.index}") for p in maps]
        elements[:] = np.array([s[0] for s in split], dtype=object)
        lines[:] = np.array([s[1] for s in split], dtype=object)
        labels[:] = np.array([p.description or f"plane_{p.index}" for p in maps], dtype=object)
        _attrs(elements, long_name="element symbol (text before the hyphen of the RTX plane description)")
        _attrs(lines, long_name="line identifier (text after the hyphen of the RTX plane description)",
               note="the meaning of the suffixes (for example KA, K) is not verified")
        _attrs(labels, long_name="original RTX plane description, verbatim")
        source_plane = g.createVariable("element_map_source_plane", "i4", ("element",), fill_value=False)
        source_plane[:] = np.array([p.index for p in maps], dtype="i4")
        _attrs(source_plane, long_name="plane index in the RTX map image")
        sha = g.createVariable("element_map_sha256", str, ("element",))
        sha[:] = np.array([p.sha256 for p in maps], dtype=object)
        _attrs(sha, long_name="SHA-256 of the decoded plane bytes (little-endian uint16)")
        valid = g.createVariable("element_map_valid_flag", str, ("element",))
        valid[:] = np.array([p.valid if p.valid is not None else "" for p in maps], dtype=object)
        _attrs(valid, long_name="RTX plane field 'Valid', verbatim", meaning="unknown (opaque source field)")

        element_maps = g.createVariable("element_maps", "u2", ("element", "y", "x"),
                                        chunksizes=plan.chunks["element_maps"], **self._compression())
        element_maps.set_auto_maskandscale(False)
        _attrs(element_maps, long_name="RTX element distribution maps", source="RTX Mapdaten image planes "
               "(little-endian uint16 assumed; supported by the bit-identity of its video plane with the BCF)",
               processing_method="unknown", units_status="unknown; no units attribute is assigned",
               scientific_note=ELEMENT_MAP_NOTE, coordinates="element element_line element_label")
        self.element_maps = element_maps

        sum_header = g.createVariable("sum_spectrum_header", "u8", ("energy",), fill_value=False)
        sum_header[:] = info.sum_spectrum.astype(np.uint64)
        _attrs(sum_header, long_name="sum spectrum as recorded in the BCF header", units="counts",
               source="BCF spectrum header 'Channels'",
               note="NOT equal to the sum of the pixel spectra in the real file (FINDINGS.md section 9): the header "
                    "is larger in every channel where they differ; the cause is not verified")
        self.sum_from_pixels = g.createVariable("sum_spectrum_from_pixels", "i8", ("energy",), fill_value=False)
        _attrs(self.sum_from_pixels, long_name="sum over all pixels of the stored counts (DERIVED, computed while writing)",
               units="counts")
        line_counter = g.createVariable("line_counter", "i4", ("y",), fill_value=False)
        line_counter[:] = np.array(info.line_counter, dtype="i4")
        _attrs(line_counter, long_name="BCF header 'LineCounter' entry per scan line, verbatim", meaning="not verified")

        instrument = g.createVariable("instrument", "i4", ())
        instrument.assignValue(0)
        detector = info.original_metadata.get("Detector", {})
        _attrs(instrument, long_name="instrument and acquisition metadata read from the BCF header",
               sample_name=info.sample_name, detector_type=info.detector_type,
               detector_technology=detector.get("Technology"), detector_serial=str(detector.get("Serial", "")) or None,
               bcf_date=info.date_iso, bcf_time=info.time_iso,
               beam_energy_recorded=info.beam_energy_recorded,
               beam_energy_units_status="NOT stored in the BCF; RosettaSciIO treats the value as kV and the RTX "
                                        "annotation reads 'HV: 50,0 kV'; not independently verified, no units attribute",
               elevation_angle_recorded=info.elevation_angle_recorded,
               elevation_angle_units_status="NOT stored in the BCF; degrees per RosettaSciIO, not independently verified",
               real_time_seconds_from_rosettasciio=info.real_time_seconds,
               real_time_note="computed by RosettaSciIO from LineCounter, LineAverage, PixelAverage and PixelTime "
                              "assuming PixelTime is in microseconds (not verified)",
               complete_header="see /metadata/bcf_original_metadata_json")

        for index in range(len(maps)):
            self._targets[("rtx", layout.grid.index, maps[index].index)] = ("map", index)
        self._targets[("rtx", layout.grid.index, layout.video_plane.index)] = ("video", None)
        if layout.header_video_image is not None:
            hv = layout.header_video_image
            self._targets[("bcf", hv.index, hv.planes[0].index)] = ("header_video", None)
            video.header_image_index = hv.index
        for aux in layout.aux_images:
            self._create_aux(g, aux)

    def _create_aux(self, g: netCDF4.Group, aux) -> None:
        image = aux.image
        dtype = {1: "u1", 2: "u2", 4: "u4"}[image.itemsize]
        variable = g.createVariable(aux.variable, dtype, ("y", "x"), chunksizes=self.plan.chunks["video"],
                                    **self._compression())
        variable.set_auto_maskandscale(False)
        note = ""
        if image.name == "PixelTimes":
            note = ("the recorded sum (see data_sum) is close to the acquisition real time of the header expressed "
                    "in microseconds, which suggests per-pixel times in microseconds; this is NOT verified, so no "
                    "units attribute is assigned")
        _attrs(variable, long_name=f"image stored in the BCF header: {image.name!r} (not returned by RosettaSciIO)",
               source=f"BCF header TRTImageData index {image.index}, name {image.name!r}",
               units_status="unknown; no units attribute is assigned", note=note or None,
               width=image.width, height=image.height, item_size_bytes=image.itemsize,
               dtype_assumption="little-endian unsigned integer (ItemSize) is assumed",
               calibration_raw_text=f"{image.x_calibration_text},{image.y_calibration_text}",
               sha256=image.planes[0].sha256, date=image.date, time=image.time,
               grid_shared_with="counts (the image has the acquisition height and width, from the same BCF header)")
        self._aux_variables[aux.variable] = variable
        self._targets[("bcf", image.index, image.planes[0].index)] = ("aux", aux)

    def _create_mosaic(self, g: netCDF4.Group) -> None:
        layout, plan, info = self.layout, self.plan, self.info
        first = layout.mosaics[0].image if layout.mosaics else None
        _attrs(g, description="Optical mosaic at its native resolution, on its own grid (independent of the "
                              "acquisition grid). One subgroup per source instance (RTX images and BCF header "
                              "images) keeps its own source, timestamp, annotations and footprint.")
        if first is None:
            return
        my, mx = first.height, first.width
        g.createDimension("mosaic_y", my)
        g.createDimension("mosaic_x", mx)
        g.createDimension("plane", 3)
        for name, size, cal, text in (("mosaic_y", my, first.y_calibration, first.y_calibration_text),
                                      ("mosaic_x", mx, first.x_calibration, first.x_calibration_text)):
            self._coordinate(
                g, name, np.arange(size, dtype="f8") * cal, "um", long_name=f"mosaic {name[-1]} coordinate",
                calibration_raw_text=text,
                units_status="no unit is stored next to the calibration; micrometres per pixel by derivation from "
                             "the '40000 um' scale-bar annotation of the RTX (FINDINGS.md 8.4)",
                coordinate_rule="coordinate = pixel_index * calibration; origin at the mosaic's first pixel")
        for entry in layout.mosaics:
            image = entry.image
            if entry.writes_pixels:
                pixels = g.createVariable(entry.pixels_variable, "u1", ("plane", "mosaic_y", "mosaic_x"),
                                          chunksizes=plan.chunks["mosaic"], **self._compression())
                pixels.set_auto_maskandscale(False)
                same = [m.key for m in layout.mosaics if m.pixels_variable == entry.pixels_variable]
                _attrs(pixels, long_name="optical mosaic pixels, 3 planes, native resolution",
                       source="RTX 'Video Mosaic' image planes 0..2 and/or the BCF header 'Overview' image "
                              "(see stored_once_for_instances)",
                       stored_once_for_instances=json.dumps(same),
                       equal_pixels_note="the instances listed are byte-identical (per-plane SHA-256); they are NOT "
                                         "assumed scientifically interchangeable: sources, timestamps and "
                                         "annotations differ and are kept in the instance subgroups",
                       plane_semantics="planes 0, 1, 2 correlate with the R, G, B channels of the supplied "
                                       "'Video Mosaic.png' (r 0.95 to 0.97, FINDINGS.md 3); colour order otherwise "
                                       "not verified",
                       plane_sha256=json.dumps([p.sha256 for p in image.planes]),
                       orientation="rows and columns as stored in the source (row 0 first); 'Video Mosaic.png' is "
                                   "this array rotated 90 degrees clockwise",
                       units_note="no units are stored in the source")
            for plane in image.planes:
                self._targets[(entry.source, image.index, plane.index)] = ("mosaic", entry)
            self._create_instance(g.createGroup(entry.group_name), entry)

    def _instance_attrs(self, g: netCDF4.Group, entry) -> None:
        image = entry.image
        _attrs(g, long_name=f"{entry.source.upper()} image {image.index} ({image.name!r}): per-instance metadata",
               source_file="RTX payload" if entry.source == "rtx" else "BCF header (EDSDatabase/HeaderData)",
               source_image_index=image.index, source_image_name=image.name, source_date=image.date,
               source_time=image.time, timestamp_iso_naive=iso_from_rtx(image.date, image.time) or None,
               timestamp_note="day-first date inferred from 30.7.2026; the time zone and the event the time "
                              "records (start, end, save) are unknown",
               width=image.width, height=image.height, plane_count=image.plane_count,
               x_calibration_raw_text=image.x_calibration_text, y_calibration_raw_text=image.y_calibration_text,
               pixels_variable=f"/mosaic/{entry.pixels_variable}", pixels_duplicate_of=entry.duplicate_of,
               plane_sha256=json.dumps([p.sha256 for p in image.planes]),
               plane_line_counter_values=json.dumps({str(p.index): p.line_counter_values for p in image.planes}),
               plane_valid_flags=json.dumps({str(p.index): p.valid for p in image.planes}))

    def _create_instance(self, g: netCDF4.Group, entry) -> None:
        image, info = entry.image, self.info
        rect = find_map_rectangle(image.overlays)
        self._instance_attrs(g, entry)
        _text(g, "annotations_json", json.dumps(image.overlays, ensure_ascii=False),
              long_name="overlay elements of this image (type, name, text, rectangle, position), as parsed; the "
                        "complete source is in /metadata/rtx_payload_residual_xml or bcf_header_residual_xml")
        footprint = g.createVariable("map_footprint", "i4", ())
        footprint.assignValue(0)
        if rect is None:
            _attrs(footprint, long_name="Map rectangle overlay", present=0,
                   relationship_to_acquisition_grid="not established: no 'Map' rectangle in this instance")
            return
        summary = footprint_summary(rect, (image.x_calibration, image.y_calibration),
                                    (info.height, info.width), (info.pixel_size_y, info.pixel_size_x))
        _attrs(footprint, long_name="'Map' rectangle overlay in mosaic pixel coordinates (inclusive), with its "
                                    "size compared to the acquisition extent", present=1,
               **{k: (float(v) if isinstance(v, float) else (int(v) if isinstance(v, (int, bool)) else v))
                  for k, v in summary.items()})
        self._footprints[entry.key] = (g, rect, summary)

    def _create_overview(self, g: netCDF4.Group) -> None:
        layout = self.layout
        _attrs(g, description="Additional optical/overview images stored in the BCF header that RosettaSciIO does not "
                              "return. Each has its own grid. Their spatial relationship to the acquisition grid is "
                              "not established (no footprint or stage coordinates were found).")
        for overview in layout.overviews:
            image = overview.image
            sub = g.createGroup(overview.group_name)
            sub.createDimension("plane", 3)
            sub.createDimension("row", image.height)
            sub.createDimension("column", image.width)
            for name, size, cal, text in (("row", image.height, image.y_calibration, image.y_calibration_text),
                                          ("column", image.width, image.x_calibration, image.x_calibration_text)):
                variable = sub.createVariable(name, "f8", (name,), fill_value=False)
                variable[:] = np.arange(size, dtype="f8") * cal
                _attrs(variable, long_name=f"overview image {name} coordinate = index * calibration",
                       calibration_raw_text=text,
                       units_status="unknown: the calibration carries no unit in the BCF header; not assumed")
            chunks = tuple(min(a, b) for a, b in zip(self.config.mosaic_chunks, (3, image.height, image.width)))
            pixels = sub.createVariable("pixels", "u1", ("plane", "row", "column"), chunksizes=chunks,
                                        **self._compression())
            pixels.set_auto_maskandscale(False)
            _attrs(pixels, long_name=f"overview image {image.name!r} pixels, 3 planes",
                   source=f"BCF header TRTImageData index {image.index}, name {image.name!r}",
                   plane_sha256=json.dumps([p.sha256 for p in image.planes]),
                   relationship_to_acquisition_grid="not established", units_note="no units are stored")
            _attrs(sub, long_name=f"BCF header image {image.index} ({image.name!r})", source_image_index=image.index,
                   source_image_name=image.name, source_date=image.date, source_time=image.time,
                   timestamp_iso_naive=iso_from_rtx(image.date, image.time) or None, width=image.width,
                   height=image.height, x_calibration_raw_text=image.x_calibration_text,
                   y_calibration_raw_text=image.y_calibration_text)
            _text(sub, "annotations_json", json.dumps(image.overlays, ensure_ascii=False),
                  long_name="overlay elements as parsed; the complete source is in /metadata/bcf_header_residual_xml")
            for plane in image.planes:
                self._targets[("bcf", image.index, plane.index)] = ("overview", (overview, pixels))

    def _create_metadata(self, g: netCDF4.Group) -> None:
        info, scan, layout = self.info, self.scan, self.layout
        _attrs(g, description="Source metadata kept verbatim or parsed, provenance and validation. Fields whose "
                              "meaning is unknown are preserved without interpretation.")
        _text(g, "bcf_original_metadata_json", original_metadata_json(info),
              long_name="every BCF header dictionary exposed by RosettaSciIO (Hardware, Detector, Analysis, Spectrum, "
                        "DSP Configuration, Stage, Microscope, HyperHeader), verbatim values")
        _text(g, "rtx_outer_prefix_text", scan.outer_prefix_text, long_name="RTX text before <RTData>, verbatim",
              rtx_outer_fields=json.dumps(scan.outer_header_fields), rtx_compression=json.dumps(scan.outer_attributes),
              note="the outer Date/Time (30.7.2026 15:58:47) is not one of the image timestamps; meaning unknown")
        _text(g, "rtx_outer_tail_text", scan.outer_tail_text, long_name="RTX text after </RTData>, verbatim")
        _text(g, "rtx_payload_residual_xml", scan.residual_xml,
              long_name="the decompressed RTX payload (Windows-1252 decoded to Unicode) with the text of every "
                        "image-plane Data element removed; all other elements, unknown fields, LineCounter, Valid, "
                        "overlays, palettes and display settings are preserved",
              payload_sha256=scan.payload_sha256, payload_bytes=int(scan.payload_bytes),
              element_count=int(scan.payload_elements),
              class_instance_types=json.dumps(scan.class_instance_types))
        table = [{"rtx_image_index": im.index, "name": im.name, "width": im.width, "height": im.height,
                  "item_size_bytes": im.itemsize, "plane_count": im.plane_count,
                  "x_calibration_raw_text": im.x_calibration_text, "y_calibration_raw_text": im.y_calibration_text,
                  "date": im.date, "time": im.time, "timestamp_iso_naive": iso_from_rtx(im.date, im.time),
                  "role": "acquisition grid image (video + element maps)" if im is layout.grid else "mosaic instance",
                  "planes": [{"index": p.index, "description": p.description, "valid": p.valid, "sha256": p.sha256,
                              "min": p.minimum, "max": p.maximum, "sum": p.total} for p in im.planes]}
                 for im in scan.images]
        _text(g, "rtx_image_table_json", json.dumps(table, ensure_ascii=False),
              long_name="every RTX image and plane with geometry, calibration text, timestamp and statistics")
        stamps = [{"source": "BCF EDX spectrum", "date": info.date_iso, "time": info.time_iso},
                  {"source": "RTX outer header", "date": scan.outer_header_fields.get("Date"),
                   "time": scan.outer_header_fields.get("Time")}]
        stamps += [{"source": f"RTX image {im.index} ({im.name})", "date": im.date, "time": im.time,
                    "iso_naive": iso_from_rtx(im.date, im.time)} for im in scan.images]
        if info.header_scan is not None:
            stamps += [{"source": f"BCF header image {im.index} ({im.name})", "date": im.date, "time": im.time,
                        "iso_naive": iso_from_rtx(im.date, im.time)} for im in info.header_scan.images]
        header = info.header_scan
        if header is not None:
            _text(g, "bcf_header_residual_xml", header.residual_xml,
                  long_name="the BCF header XML (EDSDatabase/HeaderData, Windows-1252 decoded to Unicode) with the "
                            "text of every image-plane Data element removed; every other element is preserved",
                  header_sha256=header.sha256, header_decoded_bytes=int(header.decoded_bytes),
                  element_count=int(header.elements), class_instance_types=json.dumps(header.class_instance_types))
            roles = {}
            for aux in layout.aux_images:
                roles[aux.image.index] = f"acquisition-grid image stored as /acquisition/{aux.variable}"
            for m in layout.mosaics:
                if m.source == "bcf":
                    roles[m.image.index] = f"mosaic instance /mosaic/{m.group_name}"
            for o in layout.overviews:
                roles[o.image.index] = f"overview image /overview/{o.group_name}"
            for im in layout.empty_header_images:
                roles[im.index] = "no planes: recorded here and in the residual XML only"
            if layout.header_video_image is not None:
                roles[layout.header_video_image.index] = "video image (the one RosettaSciIO returns), see /acquisition/video"
            _text(g, "bcf_image_table_json", json.dumps([
                {"header_image_index": im.index, "name": im.name, "width": im.width, "height": im.height,
                 "item_size_bytes": im.itemsize, "plane_count": im.plane_count,
                 "x_calibration_raw_text": im.x_calibration_text, "y_calibration_raw_text": im.y_calibration_text,
                 "date": im.date, "time": im.time, "role": roles.get(im.index, "unclassified"),
                 "planes": [{"index": p.index, "description": p.description, "sha256": p.sha256} for p in im.planes]}
                for im in header.images], ensure_ascii=False),
                long_name="every image found in the BCF header with its geometry, timestamp, planes and where it "
                          "was stored (RosettaSciIO returns only the video)")
        _text(g, "timestamps_json", json.dumps(stamps),
              long_name="all recorded timestamps, verbatim, kept separate; meaning and time zone unknown")
        grid = layout.grid
        _text(g, "mapdaten_annotations_json", json.dumps(grid.overlays, ensure_ascii=False),
              long_name="overlay elements of the RTX map image (including the label 'HV: 50,0 kV')")
        _text(g, "mapdaten_display_settings_json", json.dumps(grid.map_display),
              long_name="RTX map-image display settings (colours, gamma, filter); presentation parameters, not data")
        _text(g, "provenance_json", "", long_name="conversion provenance (filled at the end)")
        _text(g, "validation_json", "", long_name="validation results (filled at the end)")
        _text(g, "unresolved_properties_json", json.dumps(UNRESOLVED), long_name="scientific properties that "
              "cannot yet be independently validated")

    # ------------------------------------------------------------------------------------------
    # Incremental writes
    # ------------------------------------------------------------------------------------------

    def write_band(self, band: Band) -> None:
        data = band.data
        if int(data.max()) > self.counts_max:
            raise OverflowError(f"a decoded count exceeds the output dtype {self.counts_dtype} "
                                f"(max {self.counts_max}); refusing to truncate or wrap")
        first, count = band.first_line, data.shape[0]
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message=_SHAPE_WARNING, category=DeprecationWarning)
            self.counts[first:first + count] = data.astype(self.counts_dtype, copy=False)
        self.lines_written.update(range(first, first + count))

    def write_sum_from_pixels(self, sums: np.ndarray) -> None:
        self.sum_from_pixels[:] = sums.astype(np.int64)

    def rtx_sink(self, image_index: int, plane_index: int, raw: bytes) -> None:
        self._sink("rtx", self.scan.images, image_index, plane_index, raw)

    def bcf_header_sink(self, image_index: int, plane_index: int, raw: bytes) -> None:
        self._sink("bcf", self.info.header_scan.images, image_index, plane_index, raw)

    def _sink(self, source: str, images, image_index: int, plane_index: int, raw: bytes) -> None:
        key = (source, image_index, plane_index)
        if key not in self._targets:
            raise ValidationError(f"{source} plane {key[1:]} was not present in the scan pass")
        image = next(im for im in images if im.index == image_index)
        plane = image.planes[plane_index]
        if hashlib.sha256(raw).hexdigest() != plane.sha256:
            raise ValidationError(f"{source} plane {key[1:]} changed between the scan pass and the write pass")
        kind, target = self._targets[key]
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message=_SHAPE_WARNING, category=DeprecationWarning)
            if kind in ("video", "header_video"):
                array = np.frombuffer(raw, dtype="<u2").reshape(image.height, image.width)
                if not np.array_equal(array, self.info.video):
                    raise ValidationError(f"the {source} video plane differs from the BCF Video image")
            elif kind == "map":
                array = np.frombuffer(raw, dtype="<u2").reshape(image.height, image.width)
                self.element_maps[target, :, :] = array.astype(np.uint16, copy=False)
            elif kind == "aux":
                dtype = np.dtype({1: "u1", 2: "<u2", 4: "<u4"}[image.itemsize])
                array = np.frombuffer(raw, dtype=dtype).reshape(image.height, image.width)
                variable = self._aux_variables[target.variable]
                variable[:, :] = array.astype(dtype.newbyteorder("="), copy=False)
                _attrs(variable, data_min=int(array.min()), data_max=int(array.max()),
                       data_sum=int(array.sum(dtype=np.uint64)), data_stats_note="derived while writing")
            elif kind == "mosaic":
                if target.writes_pixels:
                    array = np.frombuffer(raw, dtype="u1").reshape(image.height, image.width)
                    self.dataset["mosaic"][target.pixels_variable][plane_index, :, :] = array
            elif kind == "overview":
                overview, variable = target
                variable[plane_index, :, :] = np.frombuffer(raw, dtype="u1").reshape(image.height, image.width)
        self.planes_written.add(key)

    def all_planes_written(self) -> bool:
        return self.planes_written == set(self._targets)

    def missing_planes(self) -> list[tuple[str, int, int]]:
        return sorted(set(self._targets) - self.planes_written)

    def all_lines_written(self) -> bool:
        return len(self.lines_written) == self.info.height

    def mosaic_gray_reader(self):
        """Return ``read_gray(top, bottom_excl, left, right_excl)`` over the stored mosaic (None if outside)."""
        group = self.dataset["mosaic"]
        variable = group["pixels"]
        my, mx = variable.shape[1], variable.shape[2]

        def read_gray(top: int, bottom: int, left: int, right: int):
            if top < 0 or left < 0 or bottom > my or right > mx or bottom <= top or right <= left:
                return None
            crop = variable[:, top:bottom, left:right]
            return np.asarray(crop, dtype=np.float32).mean(axis=0)

        return read_gray

    def footprints(self) -> dict[int, tuple[Any, dict[str, int], dict[str, Any]]]:
        return self._footprints

    def record_registration(self, group: netCDF4.Group, result: dict[str, Any]) -> None:
        variable = group["map_footprint"]
        status = result["status"]
        by_orientation = result.get("correlation_by_orientation", {})
        others = max((v for k, v in by_orientation.items() if k != "identity"), default=float("nan"))
        relation = (
            "verified by image correlation: the BCF Video image matches the inclusive Map rectangle of this mosaic, "
            f"area-averaged to the acquisition pixel size, in identity orientation (Pearson r "
            f"{by_orientation.get('identity', float('nan')):.4f}; every other orientation at most {others:.2f}); mosaic "
            "and acquisition axes are therefore parallel with the same direction (no flip, no rotation); whether "
            "coordinates refer to pixel corners or centres is not verified"
            if status == "verified" else f"not verified ({result.get('reason', 'see registration_json')})")
        _attrs(variable, relationship_to_acquisition_grid=relation, registration_status=status,
               registration_json=json.dumps(result))

    def note_footprint_without_registration(self, group: netCDF4.Group, summary: dict[str, Any], why: str) -> None:
        _attrs(group["map_footprint"], registration_status="not_verified",
               relationship_to_acquisition_grid=f"not verified: {why}")

    def actual_layout(self) -> dict[str, Any]:
        """Chunking, filters and shapes as HDF5 actually stored them (read from the open file)."""
        ds = self.dataset
        report: dict[str, Any] = {}
        paths = ["acquisition/counts", "acquisition/video", "acquisition/element_maps", "mosaic/pixels"]
        paths += [f"acquisition/{aux.variable}" for aux in self.layout.aux_images]
        paths += [f"overview/{o.group_name}/pixels" for o in self.layout.overviews]
        for path in paths:
            group, name = ds, path
            *groups, name = path.split("/")
            for part in groups:
                group = group[part]
            if name not in group.variables:
                continue
            variable = group[name]
            chunking = variable.chunking()
            report[path] = {"shape": list(variable.shape), "dtype": str(variable.dtype),
                            "chunks": list(chunking) if isinstance(chunking, (list, tuple)) else chunking,
                            "filters": {k: v for k, v in variable.filters().items() if v}}
        return report

    def close(self) -> None:
        if self.dataset is not None:
            self.dataset.close()
            self.dataset = None


UNRESOLVED = [
    "Processing method and units of the RTX element maps (ROI window, net counts, deconvolution) are unknown.",
    "The meaning of the RTX plane fields 'Valid' and 'LineCounter' and of the BCF 'LineCounter' is unknown.",
    "The BCF header sum spectrum exceeds the sum of the decoded pixel spectra in the real file (about 0.44 % of "
    "total counts, every differing channel); the cause is unknown.",
    "Beam-energy and elevation-angle units are not stored in the BCF; kV and degrees are consistent with an RTX "
    "annotation and RosettaSciIO but not independently verified.",
    "Spectral calibration (CalibAbs, CalibLin) and pixel size have no independent external reference; only K-alpha "
    "windows of Ca, Mn, Fe and Zn fall within +-0.15 keV of tabulated positions (not a precise validation).",
    "Whether the mosaic and acquisition coordinates refer to pixel corners or centres is not verified.",
    "Mosaic instance 0 has a Map rectangle of a different size than the acquisition extent; its meaning is unknown.",
    "The meaning and time zone of every timestamp (BCF, RTX header, three RTX images) are unknown.",
    "BCF pixel encodings flag 0, flag 1 and n_of_pulses > 0, zlib-compressed SFS containers, multi-hypermap BCF and "
    "other RTX layouts are not exercised by the real file (synthetic tests only or unsupported).",
    "The RTX 8-bit mosaic plane colour order and the little-endian reading of 16-bit RTX planes rest on indirect evidence.",
]
