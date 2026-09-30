"""
Authors: Andre L. Belem (https://github.com/andrebelem) and F.R.I.D.A.Y. (https://observatoriooceanografico.org/people/friday-bot)

Bounded-memory BCF reader (FINDINGS.md section 7, option A): a line framer that walks the record structure
of the SFS spectrum stream without decoding spectra, plus the UNMODIFIED compiled RosettaSciIO decoder
(``rsciio.bruker.unbcf_fast.parse_to_numpy``) applied to a few whole scan lines at a time.

Memory bound (decoding only): ``band_lines x width x channels x itemsize`` for the band, times about 2.4
(measured, FINDINGS.md 7.3), plus one raw record per line (about 2.6 MiB on the real file). The header of the
BCF (34 MB in the real file, about 180 MiB while RosettaSciIO parses it) is read once and released.

The framer mirrors the record layout that ``unbcf_fast.pyx`` follows. It refuses, with BCFStreamError:
an unknown ``flag``, payload sizes that would make the decoder read outside the record, ``pixel_x >= width``
(the decoder has no bounds check), more pixels than the width in a line, oversized records, a short stream
and trailing bytes. Encodings not present in the real file (flag 0, flag 1, n_of_pulses > 0) are covered by
synthetic tests only; they are unverified on real acquisitions.
"""

from __future__ import annotations

import gc
import json
import struct
from collections.abc import Iterable, Iterator
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, NamedTuple, Protocol

import numpy as np

from .errors import BCFStreamError, MicroXRFToNetCDFError
from .rtx import ImageScan, PlaneSink, scan_trt_document, stream_trt_planes

HEADER_BYTES = 0x1A0                      # the decoder seeks here to find the first line
PIXEL_HEADER = struct.Struct("<IHHIHHHI")  # pixel_x, chan1, chan2, 4 skipped bytes, flag, size1, n_pulses, size2
KNOWN_FLAGS = (0, 1, 2)                   # only 2 occurs in the real file; 0 and 1 follow the decoder source
MAX_PIXEL_BYTES = 4 << 20                 # a pixel payload larger than this is treated as corruption
DEFAULT_MAX_LINE_BYTES = 512 << 20        # a scan-line record larger than this is treated as corruption
GUARD_BYTES = 8192                        # zero padding after each band block (see decode_band)
WIDEN_THRESHOLD = 1 << 15                 # a uint16 band reaching this is re-decoded as uint32
SUPPORTED_ROSETTASCIIO = ((0, 14), (0, 14))  # inclusive (major, minor) range verified: 0.14.x
STREAM_PATH = "EDSDatabase/SpectrumData{index}"
HEADER_PATH = "EDSDatabase/HeaderData"
HEADER_ITEMSIZES = (1, 2, 4)              # 16-bit video, 8-bit overview images, 32-bit PixelTimes (little-endian assumed)


class LineRecord(NamedTuple):
    index: int
    data: bytes
    pixels: int


class Band(NamedTuple):
    first_line: int
    data: np.ndarray        # (lines, width, channels), unsigned integer
    widened: bool           # True if the band was re-decoded in uint32 because uint16 counts neared overflow


@dataclass
class BCFInfo:
    """Everything read from the BCF header (never spectra). Filled by ``open_bcf`` or by a synthetic source."""

    filename: str
    size_bytes: int
    height: int
    width: int
    channels: int
    video: np.ndarray                        # (height, width) uint16, the BCF "Video" image
    sum_spectrum: np.ndarray                 # (channels,) uint64 as recorded in the header
    calib_abs: float                         # raw CalibAbs (keV)
    calib_lin: float                         # raw CalibLin (keV per channel)
    sigma_abs: float | None
    sigma_lin: float | None
    pixel_size_y: float                      # raw Microscope.DY
    pixel_size_x: float                      # raw Microscope.DX
    units_label: str                         # the label the reader reports (a micro sign + "m" in the real file)
    sample_name: str
    date_iso: str | None
    time_iso: str | None
    beam_energy_recorded: float | None       # EsmaHeader PrimaryEnergy; units NOT stored in the file
    elevation_angle_recorded: float | None   # EsmaHeader ElevationAngle; units NOT stored in the file
    detector_type: str | None
    real_time_seconds: float | None
    line_counter: list[int]
    original_metadata: dict[str, Any]        # every header dict the reader exposes, JSON-serializable
    video_title: str = "Video"
    software: dict[str, str] = field(default_factory=dict)
    stream_bytes: int | None = None
    sfs_compression: str | None = None
    header_bytes: int | None = None          # size of EDSDatabase/HeaderData (memory estimate only)
    header_blocks: Callable[[], Iterator[bytes]] | None = None   # fresh iterator over the header XML bytes
    header_scan: ImageScan | None = None     # every TRTImageData in the header (RosettaSciIO returns only one)


def attach_header(info: BCFInfo, blocks: Callable[[], Iterator[bytes]]) -> None:
    """Scan the header XML for images (one incremental pass, one plane in memory) and keep the block factory."""
    info.header_blocks = blocks
    info.header_scan = scan_trt_document(blocks(), HEADER_ITEMSIZES)


def stream_header_planes(info: BCFInfo, sink: PlaneSink) -> None:
    """Second pass over the header XML: ``sink(image_index, plane_index, raw_bytes)`` per plane."""
    if info.header_blocks is None:
        return
    stream_trt_planes(info.header_blocks(), sink)


class BCFSource(Protocol):
    """What the converter needs from a BCF: header information and a fresh iterator over stream blocks."""

    info: BCFInfo

    def blocks(self) -> Iterator[bytes]: ...


def check_rosettasciio_version() -> str:
    """Return the installed RosettaSciIO version; raise if it is outside the verified 0.14.x range."""
    from importlib.metadata import version
    installed = version("rosettasciio")
    parts = installed.split(".")
    try:
        major_minor = (int(parts[0]), int(parts[1]))
    except (ValueError, IndexError) as error:
        raise MicroXRFToNetCDFError(f"cannot parse RosettaSciIO version {installed!r}") from error
    low, high = SUPPORTED_ROSETTASCIIO
    if not low <= major_minor <= high:
        raise MicroXRFToNetCDFError(
            f"RosettaSciIO {installed} is outside the verified range 0.14.x: the record-framing rules of "
            f"microxrf_to_netcdf.bcf mirror unbcf_fast.pyx of 0.14.0 and must be re-verified before using it")
    return installed


# ---------------------------------------------------------------------------------------------------
# Framing
# ---------------------------------------------------------------------------------------------------


def _flag1_min_payload(pulses: int) -> int:
    """Smallest payload (bytes) that the 12-bit unpacker of unbcf_fast.pyx reads for ``pulses`` pulses."""
    if pulses == 0:
        return 0
    quotient, remainder = divmod(pulses - 1, 4)
    return 6 * quotient + (1, 3, 5, 5)[remainder] + 1


def frame_lines(blocks: Iterable[bytes], height: int | None = None, width: int | None = None,
                max_line_bytes: int = DEFAULT_MAX_LINE_BYTES) -> Iterator[bytes | LineRecord]:
    """Yield the ``HEADER_BYTES`` header, then one ``LineRecord`` per scan line.

    Memory: one line record plus one SFS block. Raises BCFStreamError on any inconsistency, including
    a stream that ends early or has bytes left after the last line.
    """
    source = iter(blocks)
    buffer = bytearray()
    position = 0

    def need(count: int) -> None:
        while len(buffer) - position < count:
            try:
                buffer.extend(next(source))
            except StopIteration:
                raise BCFStreamError("the spectrum stream ended before the record was complete "
                                     "(truncated file?)") from None

    need(HEADER_BYTES)
    header = bytes(buffer[:HEADER_BYTES])
    position = HEADER_BYTES
    stream_height, stream_width = struct.unpack_from("<II", header)
    if stream_height == 0 or stream_width == 0:
        raise BCFStreamError(f"stream header declares an empty grid ({stream_height} x {stream_width})")
    if height is not None and stream_height != height:
        raise BCFStreamError(f"stream height {stream_height} differs from the header height {height}")
    if width is not None and stream_width != width:
        raise BCFStreamError(f"stream width {stream_width} differs from the header width {width}")
    yield header
    for line in range(stream_height):
        del buffer[:position]
        position = 0
        need(4)
        (pixels,) = struct.unpack_from("<I", buffer, position)
        position += 4
        if pixels > stream_width:
            raise BCFStreamError(f"line {line} announces {pixels} pixels, more than the width {stream_width}")
        for _ in range(pixels):
            need(PIXEL_HEADER.size)
            pixel_x, _c1, _c2, _skip, flag, _size1, pulses, size2 = PIXEL_HEADER.unpack_from(buffer, position)
            position += PIXEL_HEADER.size
            if pixel_x >= stream_width:  # unbcf_fast is compiled with boundscheck(False)
                raise BCFStreamError(f"pixel_x {pixel_x} >= width {stream_width} on line {line}")
            if flag not in KNOWN_FLAGS:
                raise BCFStreamError(f"unknown pixel flag {flag} on line {line}; only {KNOWN_FLAGS} are handled")
            if flag == 0:
                if size2 < 2 * pulses:
                    raise BCFStreamError(f"flag 0 payload {size2} bytes is shorter than {pulses} 16-bit pulses")
                length = size2
            elif flag == 1:
                if size2 < _flag1_min_payload(pulses):
                    raise BCFStreamError(f"flag 1 payload {size2} bytes is too short for {pulses} 12-bit pulses")
                length = size2
            else:
                if size2 < 4:
                    raise BCFStreamError(f"flag {flag} payload size {size2} is smaller than its 4-byte trailer")
                length = size2 - 4 + ((4 + 2 * pulses) if pulses > 0 else 4)
            if length > MAX_PIXEL_BYTES:
                raise BCFStreamError(f"pixel payload of {length} bytes on line {line} exceeds {MAX_PIXEL_BYTES}")
            if position + length > max_line_bytes:
                raise BCFStreamError(f"line {line} record exceeds {max_line_bytes} bytes")
            need(length)
            position += length
        yield LineRecord(line, bytes(buffer[:position]), pixels)
    trailing = len(buffer) - position
    if trailing:
        raise BCFStreamError(f"{trailing} unexpected bytes after the last line")
    try:
        next(source)
    except StopIteration:
        return
    raise BCFStreamError("unexpected blocks after the last line")


class _Block:
    """Duck-typed SFS item: ``parse_to_numpy`` only calls ``get_iter_and_properties()``."""

    def __init__(self, padded: bytes, size: int):
        self.padded, self.size = padded, size

    def get_iter_and_properties(self):
        return iter([self.padded]), self.size, 1


def decode_band(header: bytes, records: list[bytes], width: int, channels: int, dtype: Any) -> np.ndarray:
    """Decode whole scan lines with the unmodified compiled decoder.

    The header is copied with its height field set to the number of lines. ``GUARD_BYTES`` of zeros follow
    the block inside the bytes object but are excluded from the declared size, so that a payload whose
    internal structure is inconsistent (the decoder does not check it) reads zeros instead of foreign memory.
    """
    from rsciio.bruker import unbcf_fast
    patched = bytearray(header)
    struct.pack_into("<I", patched, 0, len(records))
    payload = b"".join([bytes(patched), *records])
    return unbcf_fast.parse_to_numpy(_Block(payload + bytes(GUARD_BYTES), len(payload)),
                                     (len(records), width, channels), np.dtype(dtype))


def iter_bands(source: BCFSource, band_lines: int, dtype: Any = np.uint16,
               max_line_bytes: int = DEFAULT_MAX_LINE_BYTES,
               on_line: Any = None) -> Iterator[Band]:
    """Yield ``Band(first_line, data, widened)`` over the whole stream, ``band_lines`` scan lines at a time.

    ``dtype`` is the decoding type. uint32 is exact for any realistic count. uint16 is exact ONLY while every
    true per-pixel, per-channel count is below 65536: the compiled decoder adds in place in that type, so a larger
    count wraps silently (for example 70000 becomes 4464, a value that looks harmless). A uint16 band whose
    maximum reaches 2**15 is re-decoded as uint32 as a best-effort guard, which cannot catch a wrap to a small
    value. The guarantee therefore does not come from here: the converter scans the stream once in uint32 and
    uses uint16 afterwards only if the exact maximum fits and the per-channel sums of both decodes agree exactly.
    ``on_line`` (optional) receives each ``LineRecord`` (used for statistics).
    """
    if band_lines < 1:
        raise ValueError("band_lines must be at least 1")
    info = source.info
    dtype = np.dtype(dtype)
    if dtype not in (np.dtype(np.uint16), np.dtype(np.uint32)):
        raise ValueError("the decoding dtype must be uint16 or uint32")
    framer = frame_lines(source.blocks(), info.height, info.width, max_line_bytes)
    header = next(framer)
    pending: list[bytes] = []
    first = 0
    for record in framer:
        if on_line is not None:
            on_line(record)
        pending.append(record.data)
        if len(pending) == band_lines or record.index == info.height - 1:
            data = decode_band(header, pending, info.width, info.channels, dtype)
            widened = False
            if dtype == np.dtype(np.uint16) and int(data.max()) >= WIDEN_THRESHOLD:
                data = decode_band(header, pending, info.width, info.channels, np.uint32)
                widened = True
            yield Band(first, data, widened)
            first += len(pending)
            pending = []


# ---------------------------------------------------------------------------------------------------
# Real files
# ---------------------------------------------------------------------------------------------------


def _jsonable(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    return value


def _optional_float(value: Any) -> float | None:
    try:
        return None if value is None else float(value)
    except (TypeError, ValueError):
        return None


class RealBCF:
    """A BCF file opened read-only through RosettaSciIO's SFS container reader."""

    def __init__(self, path: Path, info: BCFInfo, index: int):
        self.path, self.info, self._index = path, info, index

    def blocks(self) -> Iterator[bytes]:
        from rsciio.bruker._api import SFS_reader
        item = SFS_reader(str(self.path)).get_file(STREAM_PATH.format(index=self._index))
        return item.get_iter_and_properties()[0]


def open_bcf(path: str | Path) -> RealBCF:
    """Read the BCF header (not the spectra) and return a source. Releases the parsed header afterwards."""
    from rsciio.bruker._api import BCF_reader, SFS_reader
    path = Path(path)
    version = check_rosettasciio_version()
    if not path.is_file():
        raise MicroXRFToNetCDFError(f"BCF file not found: {path}")
    reader = BCF_reader(str(path))
    try:
        if len(reader.available_indexes) != 1:
            raise MicroXRFToNetCDFError(f"BCF holds spectrum indexes {reader.available_indexes}; only a single "
                                 "hypermap is verified")
        index = reader.def_index
        header = reader.header
        spectrum = header.spectra_data[index]
        images = header.image.images
        if len(images) != 1:
            raise MicroXRFToNetCDFError(f"expected exactly one BCF image plane, found {len(images)}")
        video = np.ascontiguousarray(images[0]["data"])
        if video.dtype != np.uint16:
            raise MicroXRFToNetCDFError(f"BCF Video image dtype {video.dtype} is not the verified uint16")
        item = reader.get_file(STREAM_PATH.format(index=index))
        header_bytes = int(reader.get_file(HEADER_PATH).size)
        original = _jsonable({
            "Hardware": spectrum.hardware_metadata, "Detector": spectrum.detector_metadata,
            "Analysis": spectrum.esma_metadata, "Spectrum": spectrum.spectrum_metadata,
            "DSP Configuration": header.dsp_metadata, "Stage": header.stage_metadata,
            "Microscope": header.sem_metadata,
            "HyperHeader": {"name": header.name, "version": header.version, "date": header.date,
                            "time": header.time, "channel_count": header.channel_count,
                            "detector_count": header.mapping_count, "elements": header.elements,
                            "mode_guess_by_rosettasciio": header.mode, "hv_from_microscope_block": header.hv},
        })
        info = BCFInfo(
            filename=path.name, size_bytes=path.stat().st_size,
            height=int(header.image.height), width=int(header.image.width),
            channels=int(spectrum.data.shape[0]), video=video,
            sum_spectrum=np.asarray(spectrum.data, dtype=np.uint64).copy(),
            calib_abs=float(spectrum.spectrum_metadata["CalibAbs"]),
            calib_lin=float(spectrum.spectrum_metadata["CalibLin"]),
            sigma_abs=_optional_float(spectrum.spectrum_metadata.get("SigmaAbs")),
            sigma_lin=_optional_float(spectrum.spectrum_metadata.get("SigmaLin")),
            pixel_size_y=float(header.y_res), pixel_size_x=float(header.x_res),
            units_label=str(header.units), sample_name=str(header.name),
            date_iso=str(header.date) if header.date else None,
            time_iso=str(header.time) if header.time else None,
            beam_energy_recorded=_optional_float(spectrum.esma_metadata.get("PrimaryEnergy")),
            elevation_angle_recorded=_optional_float(spectrum.esma_metadata.get("ElevationAngle")),
            detector_type=str(spectrum.detector_type) if spectrum.detector_type else None,
            real_time_seconds=float(header.calc_real_time()),
            line_counter=[int(v) for v in header.line_counter],
            original_metadata=original,
            software={"rosettasciio": version},
            stream_bytes=int(item.size), sfs_compression=str(item.sfs.compression),
            header_bytes=header_bytes,
        )
    finally:
        del reader
        gc.collect()
    attach_header(info, lambda: SFS_reader(str(path)).get_file(HEADER_PATH).get_iter_and_properties()[0])
    if info.sfs_compression not in ("None", None):
        raise MicroXRFToNetCDFError(f"SFS compression {info.sfs_compression!r} is not verified (only 'None')")
    if info.channels < 1 or len(info.line_counter) != info.height:
        raise MicroXRFToNetCDFError("inconsistent BCF header (channel count or line counter)")
    if video.shape != (info.height, info.width):
        raise MicroXRFToNetCDFError(f"BCF Video shape {video.shape} differs from the grid {(info.height, info.width)}")
    return RealBCF(path, info, index)


def original_metadata_json(info: BCFInfo) -> str:
    return json.dumps(info.original_metadata, sort_keys=True, ensure_ascii=False)
