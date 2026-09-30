"""
Authors: Andre L. Belem (https://github.com/andrebelem) and F.R.I.D.A.Y. (https://observatoriooceanografico.org/people/friday-bot)

Synthetic BCF record streams and RTX files for the fast tests. Not a test module.

The BCF encoder writes the record layout that ``unbcf_fast.pyx`` reads (FINDINGS.md 7.1) for every pixel
encoding: ``flag`` 0 (16-bit pulses), ``flag`` 1 (12-bit pulses), ``flag`` 2 (instructed packing) with and without
extra pulses. The real file exercises only ``flag`` 2 without extra pulses; the other encodings are therefore
verified here by synthetic streams only, and FINDINGS.md says so.
"""

from __future__ import annotations

import struct
from collections.abc import Callable, Iterator
from pathlib import Path

import numpy as np

from microxrf_to_netcdf.bcf import HEADER_BYTES, BCFInfo, attach_header

from test_rtx_synthetic import b64, write_rtx

PIXEL_HEADER = struct.Struct("<IHHIHHHI")


# ---------------------------------------------------------------------------------------------------
# BCF stream encoder
# ---------------------------------------------------------------------------------------------------


def _pulses_of(spectrum: np.ndarray) -> list[int]:
    return [int(c) for c, n in enumerate(spectrum) for _ in range(int(n))]


def _pack12(channels: list[int]) -> bytes:
    out = bytearray()
    for g in range(0, len(channels), 4):
        c = (channels[g:g + 4] + [0, 0, 0, 0])[:4]
        out += bytes([((c[0] & 0xF) << 4) | (c[1] >> 8), c[0] >> 4, c[2] >> 4, c[1] & 0xFF, c[3] & 0xFF,
                      ((c[2] & 0xF) << 4) | (c[3] >> 8)])
    return bytes(out)


def _instructed(spectrum: np.ndarray) -> bytes:
    """Bunches: zero runs, then nibble (size 1), byte (size 2) or word (size 4, 32-bit gain) runs."""
    out = bytearray()
    i, n = 0, len(spectrum)
    while i < n:
        j = i
        if spectrum[i] == 0:
            while j < n and spectrum[j] == 0 and j - i < 255:
                j += 1
            out += bytes([0, j - i])
        else:
            while j < n and spectrum[j] != 0 and j - i < 255:
                j += 1
            segment = spectrum[i:j].astype(np.int64)
            gain = int(segment.min())
            residual = segment - gain
            top = int(residual.max())
            if top > 65535:  # a range no run can express: one single-channel run per channel (32-bit gain)
                for value in segment:
                    out += bytes([4, 1]) + struct.pack("<I", int(value)) + struct.pack("<H", 0)
            elif top <= 15 and gain <= 255:
                out += bytes([1, len(segment), gain])
                data = bytearray((len(segment) + 1) // 2)
                for k, v in enumerate(residual):
                    data[k // 2] |= int(v) if k % 2 == 0 else int(v) << 4
                out += data
            elif top <= 255 and gain <= 65535:
                out += bytes([2, len(segment)]) + struct.pack("<H", gain) + bytes(int(v) for v in residual)
            else:
                out += bytes([4, len(segment)]) + struct.pack("<I", gain)
                out += b"".join(struct.pack("<H", int(v)) for v in residual)
        i = j
    return bytes(out)


def encode_pixel(x: int, spectrum: np.ndarray, flag: int, extra: int = 0) -> bytes:
    """One pixel record. ``extra`` (flag 2 only) moves that many counts into the trailing pulse list."""
    spectrum = spectrum.copy()
    if flag == 0:
        pulses = _pulses_of(spectrum)
        payload = b"".join(struct.pack("<H", c) for c in pulses)
        return PIXEL_HEADER.pack(x, 0, 0, 0, 0, 0, len(pulses), len(payload)) + payload
    if flag == 1:
        pulses = _pulses_of(spectrum)
        assert all(c < 4096 for c in pulses), "12-bit pulses carry channels below 4096"
        payload = _pack12(pulses)
        return PIXEL_HEADER.pack(x, 0, 0, 0, 1, 0, len(pulses), len(payload)) + payload
    extra_channels: list[int] = []
    for c in np.nonzero(spectrum)[0][:extra]:
        spectrum[c] -= 1
        extra_channels.append(int(c))
    body = _instructed(spectrum)
    tail = (struct.pack("<I", 2 * len(extra_channels)) + b"".join(struct.pack("<H", c) for c in extra_channels)
            if extra_channels else b"\0\0\0\0")
    return PIXEL_HEADER.pack(x, 0, 0, 0, 2, 0, len(extra_channels), len(body) + 4) + body + tail


def encode_line(spectra: np.ndarray, flags: Callable[[int, int], tuple[int, int]], y: int) -> bytes:
    """One scan line of ``spectra`` (width, channels); ``flags(y, x)`` returns (flag, extra_pulses)."""
    out = bytearray(struct.pack("<I", spectra.shape[0]))
    for x in range(spectra.shape[0]):
        flag, extra = flags(y, x)
        out += encode_pixel(x, spectra[x], flag, extra)
    return bytes(out)


def stream_header(height: int, width: int) -> bytes:
    return struct.pack("<II", height, width) + bytes(HEADER_BYTES - 8)


def cycling_flags(y: int, x: int) -> tuple[int, int]:
    """Exercise every encoding: flag 0, flag 1, flag 2, flag 2 with two extra pulses."""
    return [(0, 0), (1, 0), (2, 0), (2, 2)][(x + y) % 4]


def build_stream(counts: np.ndarray, flags: Callable[[int, int], tuple[int, int]] = cycling_flags) -> bytes:
    height, width, _ = counts.shape
    return stream_header(height, width) + b"".join(encode_line(counts[y], flags, y) for y in range(height))


def split_blocks(stream: bytes, size: int) -> list[bytes]:
    return [stream[i:i + size] for i in range(0, len(stream), size)]


def random_counts(height=5, width=7, channels=40, seed=1, high=6, density=0.3, big=None) -> np.ndarray:
    rng = np.random.default_rng(seed)
    counts = rng.integers(0, high, size=(height, width, channels)) * (rng.random((height, width, channels)) < density)
    if big is not None:  # {(y, x, channel): value}
        for (y, x, c), value in big.items():
            counts[y, x, c] = value
    return counts.astype(np.int64)


class FixedBlockFile:
    """Duck-typed SFS item with a FIXED block size, as the public decoder expects."""

    def __init__(self, stream: bytes, block_size: int = 4064):
        self.blocks_, self.size = split_blocks(stream, block_size), block_size

    def get_iter_and_properties(self):
        return iter(self.blocks_), self.size, len(self.blocks_)


# ---------------------------------------------------------------------------------------------------
# Synthetic BCF source
# ---------------------------------------------------------------------------------------------------


class SyntheticBCF:
    """A ``BCFSource`` over an in-memory stream or a block-generator factory."""

    def __init__(self, info: BCFInfo, blocks_factory: Callable[[], Iterator[bytes]]):
        self.info, self._factory = info, blocks_factory

    def blocks(self) -> Iterator[bytes]:
        return self._factory()


def make_info(height: int, width: int, channels: int, video: np.ndarray | None = None,
              sum_spectrum: np.ndarray | None = None, pixel_size: float = 100.0, calib_abs: float = -0.96079607,
              calib_lin: float = 0.010001, name: str = "synthetic") -> BCFInfo:
    if video is None:  # structured noise: unlike a ramp it is not (nearly) symmetric under flips
        video = np.random.default_rng(12345).integers(0, 65536, size=(height, width))
    return BCFInfo(
        filename=f"{name}.bcf", size_bytes=0, height=height, width=width, channels=channels,
        video=np.ascontiguousarray(video, dtype=np.uint16),
        sum_spectrum=sum_spectrum if sum_spectrum is not None else np.full(channels, 2**40, dtype=np.uint64),
        calib_abs=calib_abs, calib_lin=calib_lin, sigma_abs=0.001, sigma_lin=0.0004, pixel_size_y=pixel_size,
        pixel_size_x=pixel_size, units_label="µm", sample_name="Synthetic", date_iso="2026-01-02",
        time_iso="03:04:05", beam_energy_recorded=50.0, elevation_angle_recorded=50.0, detector_type="Synthetic",
        real_time_seconds=1.0, line_counter=[3] * height,
        original_metadata={"Detector": {"Serial": 1, "Technology": "SDD"}, "Spectrum": {"CalibAbs": calib_abs}},
        software={"rosettasciio": "test"}, stream_bytes=None, sfs_compression="None", header_bytes=1000)


def synthetic_source(counts: np.ndarray, block_size: int = 4064, header: bool | bytes = False,
                     **info_kwargs) -> SyntheticBCF:
    """``header=True`` adds a BCF-header-like document with images; bytes are used verbatim; False adds none."""
    height, width, channels = counts.shape
    stream = build_stream(counts)
    info = make_info(height, width, channels,
                     sum_spectrum=info_kwargs.pop("sum_spectrum", counts.sum(axis=(0, 1)).astype(np.uint64) + 5),
                     **info_kwargs)
    if header:
        document = header if isinstance(header, bytes) else header_document(info.video)
        attach_header(info, lambda: iter(split_blocks(document, 4064)))
    return SyntheticBCF(info, lambda: iter(split_blocks(stream, block_size)))


# ---------------------------------------------------------------------------------------------------
# Synthetic RTX
# ---------------------------------------------------------------------------------------------------

MAP_OVERLAY = ('<ClassInstance Type="TRTRectangleOverlayElement" Name="Map"><TRTOverlayElement><Pos><PosX>0</PosX>'
               '<PosY>0</PosY></Pos><Rect><Left>{l}</Left><Top>{t}</Top><Right>{r}</Right><Bottom>{b}</Bottom>'
               '</Rect></TRTOverlayElement></ClassInstance>')


def _plane_xml(index: int, description: str | None, raw: bytes, height: int, valid: str = "0") -> str:
    desc = f"<Description>{description}</Description>" if description else ""
    return (f"<Plane{index}>{desc}<Valid>{valid}</Valid><LineCounter>{','.join(['3'] * height)}</LineCounter>"
            f"<Data>{b64(raw)}</Data><Size>{len(raw)}</Size></Plane{index}>")


def image_xml(name: str, width: int, height: int, itemsize: int, planes: list[tuple[str | None, np.ndarray]],
              calibration: str, date: str, time: str, extra: str = "") -> str:
    dtype = {1: "u1", 2: "<u2", 4: "<u4"}[itemsize]
    body = "".join(_plane_xml(i, d, np.asarray(a, dtype=dtype).tobytes(), height) for i, (d, a) in enumerate(planes))
    return (f'<ClassInstance Type="TRTImageData" Name="{name}"><ItemSize>{itemsize}</ItemSize><Width>{width}</Width>'
            f'<Height>{height}</Height><PlaneCount>{len(planes)}</PlaneCount><MultiImage>0</MultiImage>'
            f'<XCalibration>{calibration}</XCalibration><YCalibration>{calibration}</YCalibration>'
            f'<Date>{date}</Date><Time>{time}</Time>{body}{extra}</ClassInstance>')


def mosaic_canvas(video: np.ndarray, scale: int = 4, margin_top: int = 5, margin_left: int = 10):
    """A mosaic that contains the video, up-sampled ``scale`` times, inside margins; and the Map rectangle."""
    gray = (video.astype(np.uint32) >> 8).astype(np.uint8)
    big = np.kron(gray, np.ones((scale, scale), dtype=np.uint8))
    canvas = np.zeros((big.shape[0] + 2 * margin_top, big.shape[1] + 2 * margin_left), dtype=np.uint8)
    canvas[margin_top:margin_top + big.shape[0], margin_left:margin_left + big.shape[1]] = big
    rect = dict(l=margin_left, t=margin_top, r=margin_left + big.shape[1] - 1, b=margin_top + big.shape[0] - 1)
    return canvas, rect


def header_document(video: np.ndarray, with_mosaic: bool = True, with_overview: bool = True,
                    pixel_times: np.ndarray | None = None, video_override: np.ndarray | None = None,
                    extra_images: str = "", mosaic_scale: int = 4) -> bytes:
    """A BCF-header-like XML with the images the real header holds: Video, Counter (no planes), PixelTimes,
    the overview mosaic 'Default' and further overview images."""
    height, width = video.shape
    shown = video_override if video_override is not None else video
    times = pixel_times if pixel_times is not None else (
        np.arange(height * width, dtype=np.uint32).reshape(height, width) * 7 + 100)
    canvas, _ = mosaic_canvas(video, mosaic_scale)
    mosaic_calibration = f"{100.0 / mosaic_scale:g}".replace(".", ",")
    parts = [image_xml("", width, height, 2, [("Video", shown)], "0,0", "30.7.2026", "9:55:52"),
             image_xml("Counter", width, height, 2, [], "0,0", "30.7.2026", "9:55:52"),
             image_xml("PixelTimes", width, height, 4, [(None, times)], "0,0", "30.7.2026", "9:55:52")]
    if with_mosaic:
        parts.append(image_xml("Default", canvas.shape[1], canvas.shape[0], 1, [(None, canvas)] * 3,
                               mosaic_calibration, "30.7.2026", "9:55:51"))
    if with_overview:
        small = (np.arange(6 * 9, dtype=np.uint8).reshape(6, 9) * 5)
        parts.append(image_xml("Image_0", 9, 6, 1, [(None, small), (None, small + 1), (None, small + 2)], "14,75",
                               "30.7.2026", "9:55:51"))
    text = ('<?xml version="1.0" encoding="WINDOWS-1252" standalone="yes"?><TRTSpectrumDatabase><RTHeader/>'
            '<ClassInstance Type="TRTSpectrumDatabase" Name="synthetic"><Header><Date>30.7.2026</Date>'
            '<Time>9:55:52</Time></Header>' + "".join(parts) + extra_images
            + "</ClassInstance></TRTSpectrumDatabase>")
    return text.encode("windows-1252")


def write_pair_rtx(path: Path, video: np.ndarray, element_names=("Ca-KA", "Fe-KA", "Al-K"), seed: int = 3,
                   calibration: str = "100,0", grid_size: tuple[int, int] | None = None,
                   mosaics: str = "consistent", video_override: np.ndarray | None = None,
                   mosaic_scale: int = 4, element_seed_offset: int = 0) -> dict:
    """Write an RTX matching ``video`` (BCF grid). ``mosaics``: 'consistent' (two identical mosaics, the second with
    a Map rectangle that matches the video), 'flipped' (mosaic mirrored vertically), 'different' (second mosaic has other pixels), 'none'."""
    height, width = video.shape
    if grid_size is not None:
        height, width = grid_size
    rng = np.random.default_rng(seed + element_seed_offset)
    maps = [(n, rng.integers(0, 500, size=(height, width))) for n in element_names]
    shown = video_override if video_override is not None else video
    grid = image_xml("Mapdaten", width, height, 2, [("Video 1", shown)] + maps, calibration, "30.7.2026", "9:55:53")
    canvas, rect = mosaic_canvas(video, mosaic_scale)
    if mosaics == "flipped":  # same rectangle, but the optical mosaic is mirrored vertically
        canvas = np.flipud(canvas)
    mh, mw = canvas.shape
    mosaic_calibration = f"{100.0 / mosaic_scale:g}".replace(".", ",")
    images = []
    if mosaics != "none":
        planes0 = [(None, canvas)] * 3
        other = np.roll(canvas, 3, axis=1)
        planes1 = planes0 if mosaics in ("consistent", "flipped") else [(None, other)] * 3
        images.append(image_xml("Video Mosaic", mw, mh, 1, planes0, mosaic_calibration, "30.7.2026", "9:38:26",
                                extra=MAP_OVERLAY.format(l=0, t=0, r=0, b=0)
                                + MAP_OVERLAY.format(l=1, t=1, r=5, b=5)))
        images.append(image_xml("Video Mosaic", mw, mh, 1, planes1, mosaic_calibration, "30.7.2026", "9:55:51",
                                extra=MAP_OVERLAY.format(l=0, t=0, r=0, b=0) + MAP_OVERLAY.format(**rect)))
    images.append(grid)
    from test_rtx_synthetic import payload_bytes
    write_rtx(path, payload_bytes(*images))
    return {"maps": maps, "canvas": canvas, "rect": rect, "height": height, "width": width}


# ---------------------------------------------------------------------------------------------------
# Large sparse synthetic stream, generated on the fly (bounded-memory tests)
# ---------------------------------------------------------------------------------------------------


def sparse_line(y: int, width: int, channels: int, peaks: int, seed: int) -> tuple[np.ndarray, np.ndarray]:
    """Deterministic ``peaks`` (channel, count) entries per pixel of scan line ``y`` (channels may repeat)."""
    rng = np.random.default_rng(seed * 100003 + y)
    return rng.integers(0, channels, size=(width, peaks)), rng.integers(1, 10, size=(width, peaks))


def sparse_expected_line(y: int, width: int, channels: int, peaks: int, seed: int) -> np.ndarray:
    index, value = sparse_line(y, width, channels, peaks, seed)
    dense = np.zeros((width, channels), dtype=np.int64)
    np.add.at(dense, (np.repeat(np.arange(width), peaks), index.ravel()), value.ravel())
    return dense


def encode_sparse_line(y: int, width: int, channels: int, peaks: int, seed: int) -> bytes:
    """One scan line in ``flag`` 2 encoding; only the non-zero channels are written (zero runs between them)."""
    index, value = sparse_line(y, width, channels, peaks, seed)
    out = bytearray(struct.pack("<I", width))
    for x in range(width):
        merged: dict[int, int] = {}
        for c, v in zip(index[x].tolist(), value[x].tolist()):
            merged[c] = merged.get(c, 0) + v
        body, previous = bytearray(), 0
        for c in sorted(merged):
            gap = c - previous
            while gap > 0:
                step = min(gap, 255)
                body += bytes([0, step])
                gap -= step
            body += bytes([2, 1]) + struct.pack("<H", merged[c]) + b"\0"
            previous = c + 1
        out += PIXEL_HEADER.pack(x, 0, 0, 0, 2, 0, 0, len(body) + 4) + bytes(body) + b"\0\0\0\0"
    return bytes(out)


def sparse_blocks(height: int, width: int, channels: int, peaks: int, seed: int, block: int = 4064) -> Iterator[bytes]:
    """Yield the whole stream in blocks, one scan line in memory at a time."""
    yield stream_header(height, width)
    for y in range(height):
        line = encode_sparse_line(y, width, channels, peaks, seed)
        for start in range(0, len(line), block):
            yield line[start:start + block]
