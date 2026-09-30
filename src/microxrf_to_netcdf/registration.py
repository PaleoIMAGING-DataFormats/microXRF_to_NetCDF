"""
Authors: Andre L. Belem (https://github.com/andrebelem) and F.R.I.D.A.Y. (https://observatoriooceanografico.org/people/friday-bot)

Spatial relationship between an optical mosaic and the acquisition grid (FINDINGS.md 8.4 and section 9).

The RTX draws a rectangle named ``Map`` on each mosaic. Its size can be compared with the acquisition
extent (derived), and its content can be compared with the BCF video image (measured): the mosaic crop
inside the rectangle, resampled by area averaging to the acquisition grid, is correlated with the video for
the four axis orientations (identity, vertical flip, horizontal flip, both) and for small shifts.

The relationship is called ``verified`` only when the identity orientation correlates at r >= 0.9 and beats
every other orientation by at least 0.3. Otherwise the status says exactly what was not established.
"""

from __future__ import annotations

from typing import Any

import numpy as np
from PIL import Image

VERIFIED_R = 0.90
VERIFIED_MARGIN = 0.30
SHIFTS = range(-6, 7, 2)     # mosaic pixels, both axes
SIZE_TOLERANCE_PX = 2        # footprint size versus acquisition extent, in mosaic pixels


def find_map_rectangle(overlays: list[dict[str, Any]]) -> dict[str, int] | None:
    """The single non-degenerate rectangle overlay named ``Map`` (Left, Top, Right, Bottom in mosaic pixels).

    The real RTX holds two ``Map`` rectangle records per mosaic: a degenerate one (all four values 0) and the
    drawn footprint. Degenerate records are ignored (they stay in the preserved annotations); if zero or several
    non-degenerate ones remain the footprint is ambiguous and None is returned.
    """
    found = []
    for overlay in overlays:
        if overlay.get("type") == "TRTRectangleOverlayElement" and overlay.get("name") == "Map" and "rect" in overlay:
            rect = overlay["rect"]
            if all(key in rect for key in ("Left", "Top", "Right", "Bottom")):
                rect = {key: int(rect[key]) for key in ("Left", "Top", "Right", "Bottom")}
                if rect["Right"] > rect["Left"] and rect["Bottom"] > rect["Top"]:
                    found.append(rect)
    return found[0] if len(found) == 1 else None


def footprint_summary(rect: dict[str, int], mosaic_cal: tuple[float, float], grid_size: tuple[int, int],
                      grid_pixel_size: tuple[float, float]) -> dict[str, Any]:
    """Compare the inclusive rectangle with the acquisition extent. ``grid_size`` is (height, width)."""
    width_px = rect["Right"] - rect["Left"] + 1
    height_px = rect["Bottom"] - rect["Top"] + 1
    expected_w = grid_size[1] * grid_pixel_size[1] / mosaic_cal[0]
    expected_h = grid_size[0] * grid_pixel_size[0] / mosaic_cal[1]
    return {
        "rect_left": rect["Left"], "rect_top": rect["Top"], "rect_right": rect["Right"], "rect_bottom": rect["Bottom"],
        "rect_width_px_inclusive": width_px, "rect_height_px_inclusive": height_px,
        "rect_width_um": width_px * mosaic_cal[0], "rect_height_um": height_px * mosaic_cal[1],
        "acquisition_width_um": grid_size[1] * grid_pixel_size[1],
        "acquisition_height_um": grid_size[0] * grid_pixel_size[0],
        "expected_width_px": expected_w, "expected_height_px": expected_h,
        "size_matches": abs(width_px - expected_w) <= SIZE_TOLERANCE_PX and abs(height_px - expected_h) <= SIZE_TOLERANCE_PX,
    }


def _resize_mean(crop: np.ndarray, height: int, width: int) -> np.ndarray:
    return np.asarray(Image.fromarray(np.ascontiguousarray(crop, dtype=np.float32), mode="F")
                      .resize((width, height), Image.BOX))


def _pearson(a: np.ndarray, b: np.ndarray) -> float:
    a, b = a.ravel().astype(np.float64), b.ravel().astype(np.float64)
    if a.std() == 0 or b.std() == 0:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def correlate_footprint(read_gray: Any, rect: dict[str, int], video: np.ndarray,
                        mosaic_shape: tuple[int, int]) -> dict[str, Any]:
    """Correlate the video with the mosaic crop. ``read_gray(top, bottom_excl, left, right_excl)`` returns the
    grey crop (float32, mean of the three planes) or None when the window leaves the mosaic."""
    height, width = video.shape
    orientations = {"identity": lambda a: a, "flip_vertical": np.flipud, "flip_horizontal": np.fliplr,
                    "rotate_180": lambda a: a[::-1, ::-1]}
    base = read_gray(rect["Top"], rect["Bottom"] + 1, rect["Left"], rect["Right"] + 1)
    if base is None:
        return {"status": "not_verified", "reason": "the Map rectangle lies outside the mosaic"}
    resized = _resize_mean(base, height, width)
    correlations = {name: _pearson(fn(resized), video) for name, fn in orientations.items()}
    best_shift, best_r = (0, 0), correlations["identity"]
    for dy in SHIFTS:
        for dx in SHIFTS:
            crop = read_gray(rect["Top"] + dy, rect["Bottom"] + 1 + dy, rect["Left"] + dx, rect["Right"] + 1 + dx)
            if crop is None:
                continue
            value = _pearson(_resize_mean(crop, height, width), video)
            if value > best_r:
                best_shift, best_r = (dx, dy), value
    others = max(v for k, v in correlations.items() if k != "identity")
    verified = (correlations["identity"] >= VERIFIED_R and correlations["identity"] - others >= VERIFIED_MARGIN)
    return {
        "status": "verified" if verified else "not_verified",
        "method": "Pearson correlation of the BCF Video image with the area-averaged grey mosaic crop of the "
                  "inclusive Map rectangle",
        "correlation_by_orientation": correlations,
        "best_shift_mosaic_px_dx_dy": list(best_shift),
        "correlation_at_best_shift": best_r,
        "shift_search_range_px": [min(SHIFTS), max(SHIFTS)], "shift_search_step_px": 2,
    }
