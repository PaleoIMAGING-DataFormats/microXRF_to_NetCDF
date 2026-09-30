"""
Authors: Andre L. Belem (https://github.com/andrebelem) and F.R.I.D.A.Y. (https://observatoriooceanografico.org/people/friday-bot)

Bounded-memory RTX reader. It reuses the verified decoder of ``rtx_payload.py`` (FINDINGS.md section 8):
file -> Base64 (4-character groups) -> ``zlib.decompressobj`` -> expat, holding at most ONE decoded image
plane at a time. Two sequential passes are used, each reading the file once:

1. ``scan_rtx`` decodes everything without keeping planes and returns the geometry, plane statistics and
   SHA-256 values, annotations and display settings, plus a *residual XML* document: the payload with the
   text of every image-plane ``Data`` element removed. The residual keeps every other element verbatim
   (unknown fields, ``LineCounter``, ``Valid``, overlays, palettes, timestamps), so nothing that is not pixel
   data is discarded even where its meaning is unknown.
2. ``stream_rtx_planes`` decodes again and hands each plane to a sink while it is the only plane in memory.

Memory: the read chunk, the zlib chunk, the expat buffer, one plane (5.7 MB in the real file) and the
residual XML (a few hundred kilobytes in the real file).

The BCF header (``EDSDatabase/HeaderData``) is an XML document of the same ``TRT*`` family and holds its own
``TRTImageData`` elements (FINDINGS.md section 9). ``scan_trt_document`` and ``stream_trt_planes`` are the
generic entry points over any iterator of decompressed XML chunks; ``microxrf_to_netcdf.bcf`` uses them for the header.
"""

from __future__ import annotations

import re
import xml.parsers.expat
from collections.abc import Callable, Iterator
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any
from xml.sax.saxutils import escape, quoteattr

from . import rtx_payload as _decoder
from .errors import GridMismatchError, RTXFormatError

PlaneSink = Callable[[int, int, bytes], None]
MAX_RESIDUAL_CHARS = 64 << 20
_PLANE_TAG = re.compile(r"Plane\d+")
ITEMSIZE_DTYPES = {1: "|u1", 2: "<u2", 4: "<u4"}   # little-endian is an assumption (see FINDINGS.md 8.8)


@dataclass
class RTXPlane:
    index: int
    description: str | None
    valid: str | None
    size: int | None
    sha256: str
    decoded_bytes: int
    dtype: str | None
    minimum: int | None
    maximum: int | None
    total: int | None
    line_counter_values: dict[str, int]


@dataclass
class RTXImage:
    index: int
    name: str | None
    width: int
    height: int
    itemsize: int
    plane_count: int
    x_calibration: float               # micrometres per pixel by derivation (FINDINGS 8.4), unit not stored
    y_calibration: float
    x_calibration_text: str
    y_calibration_text: str
    date: str | None
    time: str | None
    planes: list[RTXPlane]
    overlays: list[dict[str, Any]]
    map_display: dict[str, dict[str, str]]
    payload_offset: int
    scalars: dict[str, str] = field(default_factory=dict)


@dataclass
class RTXScan:
    path: str
    size_bytes: int
    outer_prefix_text: str             # everything before <RTData> (the outer header), verbatim
    outer_tail_text: str               # text after </RTData>, verbatim
    outer_attributes: dict[str, str]   # RTCompression attributes
    outer_header_fields: dict[str, str]
    payload_sha256: str
    payload_bytes: int
    payload_root: str | None
    payload_elements: int
    class_instance_types: dict[str, int]
    images: list[RTXImage]
    residual_xml: str
    seconds: float
    peak_working_set_mib: float | None


class _Residual:
    """Re-serializes the payload without the text of image-plane ``Data`` elements (event based)."""

    def __init__(self) -> None:
        self.parser = xml.parsers.expat.ParserCreate(_decoder.PAYLOAD_ENCODING)
        self.parts: list[str] = []
        self.chars = 0
        self.stack: list[str] = []
        self.parser.StartElementHandler = self._start
        self.parser.EndElementHandler = self._end
        self.parser.CharacterDataHandler = self._text

    def _emit(self, text: str) -> None:
        self.chars += len(text)
        if self.chars > MAX_RESIDUAL_CHARS:
            raise RTXFormatError("the non-pixel content of the RTX payload exceeds the residual-XML bound "
                                 f"({MAX_RESIDUAL_CHARS} characters)")
        self.parts.append(text)

    def _start(self, tag: str, attributes: dict[str, str]) -> None:
        rendered = "".join(f" {key}={quoteattr(value)}" for key, value in attributes.items())
        self.stack.append(tag)
        self._emit(f"<{tag}{rendered}>")

    def _end(self, tag: str) -> None:
        self.stack.pop()
        self._emit(f"</{tag}>")

    def _text(self, data: str) -> None:
        if len(self.stack) > 1 and self.stack[-1] == "Data" and _PLANE_TAG.fullmatch(self.stack[-2]):
            return
        self._emit(escape(data))

    def feed(self, chunk: bytes) -> None:
        self.parser.Parse(chunk, False)

    def finish(self) -> str:
        self.parser.Parse(b"", True)
        return "".join(self.parts)


def _outer_fields(prefix: str) -> dict[str, str]:
    fields = {}
    for tag in ("Date", "Time", "Creator", "Comment"):
        match = re.search(rf"<{tag}>(.*?)</{tag}>", prefix, re.S)
        fields[tag] = match.group(1) if match else ""
    return fields


def _to_float(text: str | None, what: str) -> float:
    try:
        return _decoder.parse_decimal(text or "")
    except ValueError as error:
        raise RTXFormatError(f"{what} is not a number: {text!r}") from error


def _build_image(index: int, raw: dict[str, Any], allowed_itemsizes: tuple[int, ...] = (1, 2)) -> RTXImage:
    scalars = raw["scalars"]
    try:
        width, height = int(scalars["Width"]), int(scalars["Height"])
        itemsize, plane_count = int(scalars["ItemSize"]), int(scalars["PlaneCount"])
    except (KeyError, ValueError) as error:
        raise RTXFormatError(f"image {index} lacks a numeric Width, Height, ItemSize or PlaneCount") from error
    if itemsize not in allowed_itemsizes:
        raise RTXFormatError(f"image {index}: ItemSize {itemsize} is not one of {allowed_itemsizes} "
                             "(only those are verified)")
    if not raw.get("geometry_consistent"):
        raise RTXFormatError(f"image {index} ({raw.get('name')!r}): planes do not decode to "
                             "Width x Height x ItemSize bytes as declared")
    planes = []
    for plane_index, plane in sorted(raw["planes"].items()):
        if "decode_error" in plane:
            raise RTXFormatError(f"image {index} plane {plane_index}: {plane['decode_error']}")
        planes.append(RTXPlane(
            index=plane_index, description=plane.get("Description"), valid=plane.get("Valid"),
            size=int(plane["Size"]) if "Size" in plane else None, sha256=plane["sha256"],
            decoded_bytes=plane["decoded_bytes"], dtype=plane.get("dtype") or ITEMSIZE_DTYPES.get(itemsize),
            minimum=plane.get("min"),
            maximum=plane.get("max"), total=plane.get("sum"),
            line_counter_values=plane.get("line_counter_values", {})))
    if [p.index for p in planes] != list(range(plane_count)):
        raise RTXFormatError(f"image {index}: plane indexes {[p.index for p in planes]} do not run 0..{plane_count - 1}")
    return RTXImage(
        index=index, name=raw.get("name"), width=width, height=height, itemsize=itemsize,
        plane_count=plane_count,
        x_calibration=_to_float(scalars.get("XCalibration"), "XCalibration"),
        y_calibration=_to_float(scalars.get("YCalibration"), "YCalibration"),
        x_calibration_text=scalars.get("XCalibration", ""), y_calibration_text=scalars.get("YCalibration", ""),
        date=scalars.get("Date"), time=scalars.get("Time"), planes=planes, overlays=list(raw["overlays"]),
        map_display={key: dict(value) for key, value in raw["map_display"].items()},
        payload_offset=raw["payload_offset"], scalars=dict(scalars))


def _wrap(error: Exception) -> RTXFormatError:
    return RTXFormatError(str(error))


@dataclass
class ImageScan:
    """Images and residual XML of one TRT document (the RTX payload or the BCF header)."""

    images: list[RTXImage]
    residual_xml: str
    root: str | None
    elements: int
    class_instance_types: dict[str, int]
    decoded_bytes: int
    sha256: str
    seconds: float


def scan_trt_document(chunks: Iterator[bytes], allowed_itemsizes: tuple[int, ...] = (1, 2)) -> ImageScan:
    """Decode a whole TRT XML document once (bounded memory, no plane kept) from decompressed chunks."""
    import hashlib
    import time
    started = time.perf_counter()
    residual = _Residual()
    digest = hashlib.sha256()
    size = [0]

    def tee(source: Iterator[bytes]) -> Iterator[bytes]:
        for chunk in source:
            residual.feed(chunk)
            digest.update(chunk)
            size[0] += len(chunk)
            yield chunk

    try:
        payload = _decoder.parse_payload(tee(chunks))
        residual_xml = residual.finish()
    except _decoder.RTXError as error:
        raise _wrap(error) from error
    except xml.parsers.expat.ExpatError as error:
        raise RTXFormatError(f"document is not well-formed XML: {error}") from error
    images = [_build_image(index, raw, allowed_itemsizes) for index, raw in enumerate(payload["images"])]
    return ImageScan(images, residual_xml, payload["root"], payload["elements"], payload["class_instance_types"],
                     size[0], digest.hexdigest(), round(time.perf_counter() - started, 2))


def stream_trt_planes(chunks: Iterator[bytes], sink: PlaneSink) -> None:
    """Decode a TRT document again, calling ``sink(image_index, plane_index, raw_bytes)`` per plane."""
    try:
        _decoder.parse_payload(chunks, sink)
    except _decoder.RTXError as error:
        raise _wrap(error) from error
    except xml.parsers.expat.ExpatError as error:
        raise RTXFormatError(f"document is not well-formed XML: {error}") from error


def scan_rtx(path: str | Path) -> RTXScan:
    """Pass 1: decode the whole RTX once, keeping no plane. Raises RTXFormatError on any structural fault."""
    import time
    from .memory import peak_working_set_mib
    path = Path(path)
    if not path.is_file():
        raise RTXFormatError(f"RTX file not found: {path}")
    started = time.perf_counter()
    stats = _decoder.PayloadStats()
    try:
        outer = _decoder.read_outer_header(path)
        document = scan_trt_document(_decoder.iter_payload(path, stats))
    except _decoder.RTXError as error:
        raise _wrap(error) from error
    if not stats.zlib_finished or stats.zlib_unused_bytes:
        raise RTXFormatError("the zlib stream did not end cleanly")
    with path.open("rb") as source:
        prefix = source.read(outer["prefix_bytes"]).decode(_decoder.PAYLOAD_ENCODING)
    if not document.images:
        raise RTXFormatError("the RTX payload holds no TRTImageData element")
    return RTXScan(
        path=str(path), size_bytes=path.stat().st_size, outer_prefix_text=prefix,
        outer_tail_text=stats.outer_tail, outer_attributes=dict(outer["attributes"]),
        outer_header_fields=_outer_fields(prefix), payload_sha256=stats.sha256, payload_bytes=stats.decoded_bytes,
        payload_root=document.root, payload_elements=document.elements,
        class_instance_types=document.class_instance_types, images=document.images,
        residual_xml=document.residual_xml, seconds=round(time.perf_counter() - started, 2),
        peak_working_set_mib=peak_working_set_mib())


def stream_rtx_planes(path: str | Path, sink: PlaneSink) -> None:
    """Pass 2: decode again and call ``sink(image_index, plane_index, raw_bytes)`` for every plane.

    The sink runs while ``raw_bytes`` is the only plane in memory; it must not keep a reference to it.
    """
    try:
        _decoder.parse_payload(_decoder.iter_payload(Path(path)), sink)
    except _decoder.RTXError as error:
        raise _wrap(error) from error
    except xml.parsers.expat.ExpatError as error:
        raise RTXFormatError(f"decompressed payload is not well-formed XML: {error}") from error


def classify_images(scan: RTXScan, height: int, width: int) -> tuple[RTXImage, list[RTXImage]]:
    """Split the RTX images into the one on the BCF grid and the mosaics.

    The map image is the single image with ``Width x Height`` equal to the BCF grid and 16-bit items. Every
    other image must be an 8-bit, 3-plane image (a mosaic). Anything else is refused rather than dropped.
    """
    on_grid = [image for image in scan.images if (image.width, image.height) == (width, height)]
    if len(on_grid) != 1:
        raise GridMismatchError(
            f"{len(on_grid)} RTX images have the BCF grid {width} x {height} (exactly one is required); "
            f"RTX image sizes: {[(i.name, i.width, i.height) for i in scan.images]}")
    grid = on_grid[0]
    if grid.itemsize != 2:
        raise RTXFormatError(f"the RTX image on the BCF grid has ItemSize {grid.itemsize}, not the verified 2")
    mosaics = [image for image in scan.images if image is not grid]
    for image in mosaics:
        if image.itemsize != 1 or image.plane_count != 3:
            raise RTXFormatError(
                f"RTX image {image.index} ({image.name!r}, {image.width} x {image.height}, ItemSize "
                f"{image.itemsize}, {image.plane_count} planes) is neither the map image nor a "
                "3-plane 8-bit mosaic; refusing to drop it silently")
    return grid, mosaics
