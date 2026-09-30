"""
Authors: Andre L. Belem (https://github.com/andrebelem) and F.R.I.D.A.Y. (https://observatoriooceanografico.org/people/friday-bot)

Bounded-memory decoder of the Bruker RTX container and its TRT payload (read-only).

Verified structure (FINDINGS.md section 2): the RTX is an outer Windows-1252 XML file whose ``RTData``
element holds Base64 text of one zlib stream. The zlib stream decompresses to a second, XML-like
document (``CompData``) with no XML declaration. This module decodes both levels incrementally:

    file (1 MiB reads) -> Base64 (4-character groups) -> zlib.decompressobj -> chunks of bytes
    -> expat (Windows-1252 override) -> per-element summaries; large ``Data`` texts are Base64-decoded
    plane by plane.

Memory bound: the read chunk, the zlib output chunk, the expat buffer and at most ONE decoded image
plane (``Width x Height x ItemSize`` bytes) at a time. The encoded and decoded payloads are never held
whole and no copy of either is written to disk. Nothing is written anywhere and the input is opened
read-only. The full payload is never printed.

The command line report is tools/diagnostics/inspect_rtx.py.
"""

from __future__ import annotations

import base64
import binascii
import collections
import hashlib
import re
import sys
import time
import xml.parsers.expat
import zlib
from collections.abc import Callable, Iterator
from pathlib import Path
from typing import Any

import numpy as np

AUTHORS = "Authors: Andre L. Belem (https://github.com/andrebelem) and F.R.I.D.A.Y. (https://observatoriooceanografico.org/people/friday-bot)"

READ_CHUNK = 1 << 20               # bytes read from the file per step
HEADER_SCAN_BYTES = 1 << 20        # the outer header (about 300 bytes in the real file) must lie within this
RTDATA_OPEN = b"<RTData>"
RTDATA_CLOSE = b"</RTData>"
SUPPORTED_COMPRESSOR = "zlib"
SUPPORTED_ENCODER = "base64"
PAYLOAD_ENCODING = "windows-1252"  # payload has no XML declaration; expat is told explicitly
SIGNATURE_BYTES = 64               # leading decompressed bytes kept as the signature (never more)
MAX_SCALAR_CHARS = 256             # longer element texts are not stored
MAX_LINE_COUNTER_CHARS = 4096      # LineCounter texts (about 2 characters per image line) are stored
ITEM_DTYPES = {1: np.dtype("u1"), 2: np.dtype("<u2")}  # ItemSize -> dtype (little-endian is an assumption)
IMAGE_SCALARS = ("ItemSize", "Width", "Height", "PlaneCount", "MultiImage",
                 "XCalibration", "YCalibration", "Date", "Time")


class RTXError(ValueError):
    """The RTX does not have the structure verified for this format."""


def parse_decimal(text: str) -> float:
    """Parse the payload's decimal-comma numbers (``31,308984675955``)."""
    return float(text.replace(",", "."))


# --------------------------------------------------------------------------------------------------
# Level 1: outer XML -> Base64 -> zlib
# --------------------------------------------------------------------------------------------------


def read_outer_header(path: Path) -> dict[str, Any]:
    """Return the outer-XML text before ``<RTData>`` and the declared compressor and encoder."""
    with path.open("rb") as source:
        head = source.read(HEADER_SCAN_BYTES)
    index = head.find(RTDATA_OPEN)
    if index < 0:
        raise RTXError(f"no <RTData> element in the first {HEADER_SCAN_BYTES} bytes of the file")
    prefix = head[:index].decode(PAYLOAD_ENCODING)
    match = re.search(r"<RTCompression\b([^>]*)/?>", prefix)
    if match is None:
        raise RTXError("no <RTCompression> element before <RTData>")
    attributes = dict(re.findall(r'(\w+)="([^"]*)"', match.group(1)))
    if attributes.get("compressor") != SUPPORTED_COMPRESSOR or attributes.get("encoder") != SUPPORTED_ENCODER:
        raise RTXError(f"unsupported RTCompression {attributes!r}; verified only zlib + base64")
    return {"prefix_bytes": index + len(RTDATA_OPEN), "attributes": attributes}


class PayloadStats:
    """Counters filled while the decoded payload streams by (no payload bytes are retained)."""

    def __init__(self) -> None:
        self.base64_chars = 0
        self.decoded_bytes = 0
        self.zlib_finished = False
        self.zlib_unused_bytes = 0
        self.outer_tail = ""
        self.first_bytes = b""
        self._sha256 = hashlib.sha256()

    @property
    def sha256(self) -> str:
        return self._sha256.hexdigest()


def iter_payload(path: Path, stats: PayloadStats | None = None) -> Iterator[bytes]:
    """Yield the decompressed ``RTData`` payload in chunks; fill ``stats`` as it goes.

    Raises ``RTXError`` on invalid Base64, a truncated or trailing zlib stream, or a missing
    ``</RTData>``. ``stats.outer_tail`` receives the text after ``</RTData>`` (at most 64 characters).
    """
    stats = stats if stats is not None else PayloadStats()
    read_outer_header(path)
    decompressor = zlib.decompressobj()
    pending = b""            # Base64 characters not yet forming a whole 4-character group
    carry = b""              # end of the previous block, searched for a split tag
    started = closed = False
    with path.open("rb") as source:
        while not closed:
            block = source.read(READ_CHUNK)
            if not block:
                break
            block = carry + block
            carry = b""
            if not started:
                index = block.find(RTDATA_OPEN)
                if index < 0:
                    carry = block[-len(RTDATA_OPEN):]
                    continue
                block, started = block[index + len(RTDATA_OPEN):], True
            end = block.find(RTDATA_CLOSE)
            if end >= 0:
                stats.outer_tail = (block[end + len(RTDATA_CLOSE):][:64]).decode(PAYLOAD_ENCODING)
                block, closed = block[:end], True
            elif len(block) >= len(RTDATA_CLOSE):
                carry = block[-len(RTDATA_CLOSE) + 1:]
                block = block[:-len(RTDATA_CLOSE) + 1]
            else:
                carry, block = block, b""
            text = pending + b"".join(block.split())
            usable = len(text) // 4 * 4
            pending = text[usable:]
            stats.base64_chars += usable
            if usable:
                try:
                    encoded = base64.b64decode(text[:usable], validate=True)
                except binascii.Error as error:
                    raise RTXError(f"invalid Base64 near character {stats.base64_chars}: {error}") from error
                try:
                    chunk = decompressor.decompress(encoded)
                except zlib.error as error:
                    raise RTXError(f"zlib decompression failed: {error}") from error
                if decompressor.eof and decompressor.unused_data:
                    raise RTXError("bytes follow the end of the zlib stream")
                if chunk:
                    if len(stats.first_bytes) < SIGNATURE_BYTES:
                        stats.first_bytes += chunk[:SIGNATURE_BYTES - len(stats.first_bytes)]
                    stats.decoded_bytes += len(chunk)
                    stats._sha256.update(chunk)
                    yield chunk
    if not started or not closed:
        raise RTXError("no closed <RTData> element found")
    if pending:
        raise RTXError(f"{len(pending)} Base64 characters left over (not a multiple of 4)")
    try:
        tail = decompressor.flush()
    except zlib.error as error:
        raise RTXError(f"zlib decompression failed: {error}") from error
    if tail:
        stats.decoded_bytes += len(tail)
        stats._sha256.update(tail)
        yield tail
    stats.zlib_finished = decompressor.eof
    stats.zlib_unused_bytes = len(decompressor.unused_data)
    if not decompressor.eof:
        raise RTXError("zlib stream is truncated (no end-of-stream marker)")


# --------------------------------------------------------------------------------------------------
# Level 2: the decompressed CompData document
# --------------------------------------------------------------------------------------------------

PlaneSink = Callable[[int, int, bytes], None]  # (image index, plane index, plane bytes)


class _PlaneDecoder:
    """Incremental Base64 decoder for one ``Data`` element; keeps only the current plane."""

    def __init__(self) -> None:
        self.reset()

    def reset(self) -> None:
        self.chars = 0
        self.parts: list[bytes] = []
        self.carry = ""
        self.error: str | None = None

    def feed(self, text: str) -> None:
        if self.error:
            return
        text = self.carry + "".join(text.split())
        usable = len(text) // 4 * 4
        self.carry = text[usable:]
        self.chars += usable
        if usable:
            try:
                self.parts.append(base64.b64decode(text[:usable], validate=True))
            except binascii.Error as error:
                self.error = str(error)

    def finish(self) -> bytes | None:
        if self.error is None and self.carry:
            self.error = f"{len(self.carry)} Base64 characters left over"
        return None if self.error else b"".join(self.parts)


def _plane_summary(raw: bytes, itemsize: int) -> dict[str, Any]:
    """Statistics of one decoded plane (values only; the plane itself is not kept)."""
    summary: dict[str, Any] = {"sha256": hashlib.sha256(raw).hexdigest()}
    if itemsize in ITEM_DTYPES and len(raw) % itemsize == 0:
        values = np.frombuffer(raw, dtype=ITEM_DTYPES[itemsize])
        summary.update(dtype=ITEM_DTYPES[itemsize].str, min=int(values.min()), max=int(values.max()),
                       sum=int(values.sum(dtype=np.int64)), nonzero=int(np.count_nonzero(values)))
    return summary


def _geometry_consistent(image: dict[str, Any]) -> bool:
    """True if every plane decoded to exactly Width x Height x ItemSize bytes, as its Size element says."""
    scalars = image["scalars"]
    try:
        expected = int(scalars["Width"]) * int(scalars["Height"]) * int(scalars["ItemSize"])
    except (KeyError, ValueError):
        return False
    planes = image["planes"]
    return (len(planes) == int(scalars.get("PlaneCount", -1))
            and all(p.get("decoded_bytes") == expected == int(p.get("Size", -1)) for p in planes.values()))


def _new_image(name: str | None, offset: int) -> dict[str, Any]:
    return {"name": name, "payload_offset": offset, "scalars": {}, "planes": {}, "overlays": [],
            "map_display": collections.defaultdict(dict)}


def parse_payload(chunks: Iterator[bytes], plane_sink: PlaneSink | None = None) -> dict[str, Any]:
    """Stream the decompressed payload through expat and summarize it.

    ``plane_sink``, if given, receives each fully decoded plane as ``(image_index, plane_index, bytes)``
    while it is the only plane in memory; it must not keep a reference unless the caller accepts the cost.
    """
    parser = xml.parsers.expat.ParserCreate(PAYLOAD_ENCODING)
    stack: list[str] = []
    texts: list[list[str] | None] = []          # per open element; None once over MAX_SCALAR_CHARS
    text_length: list[int] = []
    instances: list[tuple[str | None, str | None, int]] = []   # (Type, Name, depth) of open ClassInstances
    overlays: list[dict[str, Any]] = []         # open overlay ClassInstances
    type_counts: collections.Counter[str] = collections.Counter()
    tag_counts: collections.Counter[str] = collections.Counter()
    images: list[dict[str, Any]] = []
    planes = _PlaneDecoder()
    state: dict[str, Any] = {"max_depth": 0, "elements": 0, "root": None, "plane": None}

    def current_image() -> dict[str, Any] | None:
        return images[-1] if images else None

    def start(tag: str, attributes: dict[str, str]) -> None:
        stack.append(tag)
        texts.append([])
        text_length.append(0)
        state["elements"] += 1
        state["max_depth"] = max(state["max_depth"], len(stack))
        tag_counts[re.sub(r"\d+$", "#", tag)] += 1   # Plane0..Plane14 -> Plane#, C0..C255 -> C#
        if state["root"] is None:
            state["root"] = tag
        if tag == "ClassInstance":
            kind, name = attributes.get("Type"), attributes.get("Name")
            instances.append((kind, name, len(stack)))
            type_counts[kind or "(no Type)"] += 1
            if kind == "TRTImageData":
                images.append(_new_image(name, parser.CurrentByteIndex))
            elif kind and "Overlay" in kind and current_image() is not None:
                overlay = {"type": kind, "name": name}
                current_image()["overlays"].append(overlay)
                overlays.append(overlay)
        elif tag == "Data" and re.fullmatch(r"Plane\d+", stack[-2] if len(stack) > 1 else ""):
            planes.reset()
            state["data_offset"] = parser.CurrentByteIndex
        elif tag.startswith("Plane") and tag[5:].isdigit() and current_image() is not None:
            index = int(tag[5:])
            state["plane"] = index
            current_image()["planes"][index] = {"payload_offset": parser.CurrentByteIndex}

    def characters(data: str) -> None:
        if not stack:
            return
        if stack[-1] == "Data" and len(stack) > 1 and re.fullmatch(r"Plane\d+", stack[-2]):
            planes.feed(data)
            return
        if texts[-1] is not None:
            text_length[-1] += len(data)
            limit = MAX_LINE_COUNTER_CHARS if stack[-1] == "LineCounter" else MAX_SCALAR_CHARS
            if text_length[-1] > limit:
                texts[-1] = None
            else:
                texts[-1].append(data)

    def end(tag: str) -> None:
        text_parts = texts.pop()
        text_length.pop()
        text = "".join(text_parts).strip() if text_parts is not None else None
        image = current_image()
        parent = stack[-2] if len(stack) > 1 else ""
        if image is not None and text is not None and text != "":
            record_scalar(image, tag, parent, text)
        if tag == "Data" and re.fullmatch(r"Plane\d+", parent) and image is not None:
            finish_plane(image)
        if tag == "ClassInstance":
            kind = instances.pop()[0]
            if kind and "Overlay" in kind and overlays and overlays[-1].get("type") == kind:
                overlays.pop()
        stack.pop()

    def record_scalar(image: dict[str, Any], tag: str, parent: str, text: str) -> None:
        path = stack
        plane = image["planes"].get(state["plane"]) if state["plane"] is not None else None
        if len(path) >= 2 and re.fullmatch(r"Plane\d+", parent) and plane is not None:
            if tag in ("Description", "Valid", "Size"):
                plane[tag] = text
            elif tag == "LineCounter":
                values = collections.Counter(text.split(","))
                plane["line_counter_values"] = dict(values)
                plane["line_counter_chars"] = len(text)
        elif tag in IMAGE_SCALARS and _is_direct_child(image, path):
            image["scalars"][tag] = text
        elif tag == "Text" and overlays:
            overlays[-1]["text"] = text
        elif tag in ("Left", "Top", "Right", "Bottom") and parent == "Rect" and overlays:
            overlays[-1].setdefault("rect", {})[tag] = int(text)
        elif tag in ("PosX", "PosY") and parent == "Pos" and overlays:
            overlays[-1].setdefault("pos", {})[tag] = int(text)
        elif tag in ("MapUsed", "MapFactor", "MapColor") and re.fullmatch(r"MapImage\d+", parent):
            image["map_display"][parent][tag] = text
        elif tag in ("Options", "Gamma", "ColorMixMethod", "MapFilterType") and parent == "ClassInstance":
            image["map_display"]["_general"][tag] = text

    def _is_direct_child(image: dict[str, Any], path: list[str]) -> bool:
        # The image ClassInstance is the innermost TRTImageData; its scalars sit directly below it.
        for index in range(len(instances) - 1, -1, -1):
            if instances[index][0] == "TRTImageData":
                return len(path) == instances[index][2] + 1
        return False

    def finish_plane(image: dict[str, Any]) -> None:
        index = state["plane"]
        plane = image["planes"][index]
        plane["data_base64_chars"] = planes.chars
        plane["data_payload_offset"] = state["data_offset"]
        raw = planes.finish()
        if raw is None:
            plane["decode_error"] = planes.error
            return
        plane["decoded_bytes"] = len(raw)
        itemsize = int(image["scalars"].get("ItemSize", 0) or 0)
        plane.update(_plane_summary(raw, itemsize))
        if plane_sink is not None:
            plane_sink(len(images) - 1, index, raw)
        planes.reset()

    parser.StartElementHandler = start
    parser.EndElementHandler = end
    parser.CharacterDataHandler = characters
    signature = b""
    try:
        for chunk in chunks:
            if len(signature) < SIGNATURE_BYTES:
                signature += chunk[:SIGNATURE_BYTES - len(signature)]
            parser.Parse(chunk, False)
        parser.Parse(b"", True)
    except xml.parsers.expat.ExpatError as error:
        raise RTXError(f"decompressed payload is not well-formed XML: {error}") from error

    for image in images:
        image["map_display"] = {key: dict(value) for key, value in image["map_display"].items()}
        image["geometry_consistent"] = _geometry_consistent(image)
    return {
        "signature_text": signature.decode("ascii", errors="replace"),
        "root": state["root"],
        "elements": state["elements"],
        "max_depth": state["max_depth"],
        "class_instance_types": dict(type_counts),
        "tag_vocabulary": sorted(tag_counts),
        "images": images,
    }


# --------------------------------------------------------------------------------------------------
# Peak memory (Windows working set, as in tools/benchmarks/probe_bcf_streaming.py)
# --------------------------------------------------------------------------------------------------


def peak_memory_mib() -> float | None:
    """Peak working set of this process in MiB; ``None`` where the Windows API is unavailable."""
    try:
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
        return round(counters.peak / 2**20, 1)
    except (AttributeError, OSError):
        return None


def inspect_rtx_payload(path: Path, plane_sink: PlaneSink | None = None) -> dict[str, Any]:
    """Decode and summarize an RTX. Never raises for a malformed file: sets ``report['error']``."""
    report: dict[str, Any] = {"path": str(path), "file_bytes": path.stat().st_size, "error": None}
    stats = PayloadStats()
    started = time.perf_counter()
    try:
        report["outer"] = read_outer_header(path)
        report["payload"] = parse_payload(iter_payload(path, stats), plane_sink)
    except RTXError as error:
        report["error"] = f"{type(error).__name__}: {error}"
    report["decode"] = {
        "base64_chars": stats.base64_chars,
        "decoded_bytes": stats.decoded_bytes,
        "zlib_finished": stats.zlib_finished,
        "zlib_unused_bytes": stats.zlib_unused_bytes,
        "outer_tail": stats.outer_tail,
        "first_bytes_hex": stats.first_bytes.hex(" "),
        "sha256_decoded": stats.sha256,
    }
    report["seconds"] = round(time.perf_counter() - started, 2)
    report["peak_working_set_mib"] = peak_memory_mib()
    return report
