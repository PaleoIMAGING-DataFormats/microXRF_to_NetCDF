"""
Authors: Andre L. Belem (https://github.com/andrebelem) and F.R.I.D.A.Y. (https://observatoriooceanografico.org/people/friday-bot)

The scientific data model rules that decide how the RTX content and the images inside the BCF header map onto
the BCF acquisition (NETCDF_SCHEMA.md).

The EDS counts and the RTX element maps share the ``y`` and ``x`` dimensions ONLY after the grid
correspondence has been validated here:

* the RTX map image has exactly the BCF height and width;
* its X/Y calibration equals the BCF pixel size (relative difference at most 1e-9);
* one plane of it is bit-identical to the BCF ``Video`` image. A video plane that differs is a
  ValidationError, never a silent merge; a plane identical to the video is stored once (from the BCF) with
  the provenance of both sources.

Everything else on the grid image is an element map. The BCF header also holds images that RosettaSciIO does not
return (FINDINGS.md section 9): they are classified with the same care. A header image on the acquisition grid
with one plane (``PixelTimes``) is stored on the shared grid; 8-bit 3-plane images of the mosaic size join the
mosaic group; other 8-bit 3-plane images become overview images with their own grids; an image without planes is
only recorded. Anything else is refused instead of being dropped.

Mosaic instances (RTX images and BCF header images) whose pixels are byte-identical share one stored pixel array,
but each instance keeps its own source, timestamp, annotations and footprint, because equal pixels do not make two
instances scientifically interchangeable.
"""

from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass, field
from datetime import datetime
from typing import Any

import numpy as np

from .bcf import BCFInfo
from .errors import GridMismatchError, RTXFormatError, ValidationError
from .rtx import RTXImage, RTXPlane, RTXScan, classify_images

CALIBRATION_RTOL = 1e-9


@dataclass
class MosaicInstance:
    source: str                     # "rtx" or "bcf"
    image: RTXImage
    pixels_variable: str            # name of the stored pixel array this instance's pixels equal
    writes_pixels: bool             # False when an earlier instance already stored identical pixels
    duplicate_of: str | None        # "source:index" of the earlier instance holding the same pixels

    @property
    def key(self) -> str:
        return f"{self.source}:{self.image.index}"

    @property
    def group_name(self) -> str:
        suffix = f"_{safe_name(self.image.name)}" if self.source == "bcf" and self.image.name else ""
        return f"instance_{self.source}_{self.image.index}{suffix}"


@dataclass
class AuxImage:
    """A one-plane BCF header image on the acquisition grid (for example ``PixelTimes``)."""

    image: RTXImage
    variable: str


@dataclass
class OverviewImage:
    """An 8-bit 3-plane BCF header image with its own grid (not the mosaic)."""

    image: RTXImage
    group_name: str


@dataclass
class Layout:
    grid: RTXImage
    video_plane: RTXPlane
    element_planes: list[RTXPlane]
    mosaics: list[MosaicInstance]
    aux_images: list[AuxImage] = field(default_factory=list)
    overviews: list[OverviewImage] = field(default_factory=list)
    empty_header_images: list[RTXImage] = field(default_factory=list)
    header_video_image: RTXImage | None = None
    checks: list[dict[str, Any]] = field(default_factory=list)

    @property
    def unique_mosaic_sets(self) -> int:
        return sum(1 for m in self.mosaics if m.writes_pixels)


def safe_name(name: str | None) -> str:
    return re.sub(r"[^A-Za-z0-9_]+", "_", (name or "unnamed").strip()).strip("_") or "unnamed"


def video_sha256(video: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(video, dtype="<u2").tobytes()).hexdigest()


def split_label(label: str) -> tuple[str, str]:
    """``Ca-KA`` -> (``Ca``, ``KA``). A plain text split: the meaning of the suffix is not verified."""
    match = re.fullmatch(r"([A-Z][a-z]?)-(.+)", label.strip())
    return (match.group(1), match.group(2)) if match else (label.strip(), "")


def iso_from_rtx(date: str | None, time: str | None) -> str:
    """``30.7.2026`` + ``9:38:26`` -> ``2026-07-30T09:38:26`` (day-first order inferred from 30.7.2026).

    Returns an empty string if the text does not parse. The time zone is unknown and is not added.
    """
    try:
        return datetime.strptime(f"{date} {time}", "%d.%m.%Y %H:%M:%S").isoformat()
    except (ValueError, TypeError):
        return ""


def _check(name: str, passed: bool, detail: str) -> dict[str, Any]:
    return {"check": name, "status": "passed" if passed else "FAILED", "detail": detail}


def _mosaic_signature(image: RTXImage) -> tuple[str, ...]:
    return tuple(p.sha256 for p in image.planes) + (f"{image.width}x{image.height}",)


def build_layout(info: BCFInfo, scan: RTXScan) -> Layout:
    """Validate the BCF/RTX grid correspondence, classify every image of both sources, decide where planes go."""
    grid, mosaic_images = classify_images(scan, info.height, info.width)
    checks = [_check("rtx_map_image_size_equals_bcf_grid", True,
                     f"{grid.width} x {grid.height} == {info.width} x {info.height}")]
    for axis, rtx_value, bcf_value in (("x", grid.x_calibration, info.pixel_size_x),
                                       ("y", grid.y_calibration, info.pixel_size_y)):
        close = abs(rtx_value - bcf_value) <= CALIBRATION_RTOL * abs(bcf_value)
        checks.append(_check(f"rtx_calibration_equals_bcf_pixel_size_{axis}", close,
                             f"RTX {rtx_value!r} vs BCF {bcf_value!r} (relative tolerance {CALIBRATION_RTOL})"))
        if not close:
            raise GridMismatchError(f"RTX {axis} calibration {rtx_value} differs from the BCF pixel size "
                                    f"{bcf_value}; the RTX maps are not on the BCF grid")
    expected = video_sha256(info.video)
    video_like = [p for p in grid.planes if (p.description or "").lower().startswith("video")]
    identical = [p for p in grid.planes if p.sha256 == expected]
    if not video_like and not identical:
        raise ValidationError("the RTX map image has no plane described as video and none equal to the BCF video")
    if video_like and (len(video_like) != 1 or video_like[0].sha256 != expected):
        raise ValidationError(
            f"the RTX video plane {[p.index for p in video_like]} is not bit-identical to the BCF Video image "
            f"(sha256 {[p.sha256[:16] for p in video_like]} vs {expected[:16]}); refusing to merge the two")
    if len(identical) != 1:
        raise ValidationError(f"{len(identical)} RTX planes equal the BCF Video image; exactly one is required")
    video_plane = identical[0]
    checks.append(_check("rtx_video_plane_bit_identical_to_bcf_video", True,
                         f"plane {video_plane.index} ({video_plane.description!r}) sha256 {expected}"))
    element_planes = [p for p in grid.planes if p is not video_plane]
    if not element_planes:
        raise RTXFormatError("the RTX map image holds no element-map plane besides the video plane")
    for plane in element_planes:
        if plane.dtype != "<u2":
            raise RTXFormatError(f"element plane {plane.index} has dtype {plane.dtype}, not the verified uint16")

    layout = Layout(grid, video_plane, element_planes, [], checks=checks)
    stored: dict[tuple[str, ...], tuple[str, str]] = {}   # signature -> (variable, instance key)

    def add_mosaic(source: str, image: RTXImage) -> None:
        signature = _mosaic_signature(image)
        if signature in stored:
            variable, first = stored[signature]
            layout.mosaics.append(MosaicInstance(source, image, variable, False, first))
        else:
            variable = "pixels" if not stored else f"pixels_{len(stored) + 1}"
            stored[signature] = (variable, f"{source}:{image.index}")
            layout.mosaics.append(MosaicInstance(source, image, variable, True, None))

    for image in mosaic_images:
        add_mosaic("rtx", image)
    shapes = {(m.image.width, m.image.height) for m in layout.mosaics}
    if len(shapes) > 1:
        raise RTXFormatError(f"mosaic instances have different sizes {sorted(shapes)}; not supported")
    _classify_header_images(info, layout, add_mosaic, shapes)
    return layout


def _classify_header_images(info: BCFInfo, layout: Layout, add_mosaic, mosaic_shapes: set) -> None:
    """Classify every image found in the BCF header (those RosettaSciIO does not return included)."""
    scan = info.header_scan
    if scan is None:
        return
    grid_shape = (info.width, info.height)
    expected_video = video_sha256(info.video)
    eight_bit = [im for im in scan.images if im.itemsize == 1 and im.plane_count == 3 and (im.width, im.height) != grid_shape]
    if not mosaic_shapes and eight_bit:
        biggest = max(eight_bit, key=lambda im: im.width * im.height)
        mosaic_shapes = {(biggest.width, biggest.height)}
    seen_names: set[str] = set()
    for image in scan.images:
        shape = (image.width, image.height)
        if image.plane_count == 0:
            layout.empty_header_images.append(image)
        elif shape == grid_shape:
            hashes = {p.sha256 for p in image.planes}
            if expected_video in hashes:
                layout.header_video_image = image
                layout.checks.append(_check("bcf_header_video_image_equals_rosettasciio_video", True,
                                            f"header image {image.index} sha256 {expected_video}"))
                continue
            if any((p.description or "").lower().startswith("video") for p in image.planes):
                raise ValidationError(f"the BCF header image {image.index} holds a plane described as video that "
                                      "differs from the video returned by RosettaSciIO; refusing to merge the two")
            if image.plane_count != 1:
                raise RTXFormatError(f"BCF header image {image.index} ({image.name!r}) on the acquisition grid has "
                                     f"{image.plane_count} planes; only single-plane images are supported")
            name = "pixel_times" if image.name == "PixelTimes" else f"aux_{safe_name(image.name or f'image_{image.index}')}"
            if name in seen_names:
                name = f"{name}_{image.index}"
            seen_names.add(name)
            layout.aux_images.append(AuxImage(image, name))
        elif image.itemsize == 1 and image.plane_count == 3:
            if shape in mosaic_shapes:
                add_mosaic("bcf", image)
            else:
                group = f"{safe_name(image.name or f'image_{image.index}')}"
                if any(o.group_name == group for o in layout.overviews):
                    group = f"{group}_{image.index}"
                layout.overviews.append(OverviewImage(image, group))
        else:
            raise RTXFormatError(f"BCF header image {image.index} ({image.name!r}, {image.width} x {image.height}, "
                                 f"ItemSize {image.itemsize}, {image.plane_count} planes) has an unsupported layout; "
                                 "refusing to drop it silently")
    shapes = {(m.image.width, m.image.height) for m in layout.mosaics}
    if len(shapes) > 1:
        raise RTXFormatError(f"mosaic instances have different sizes {sorted(shapes)}; not supported")
