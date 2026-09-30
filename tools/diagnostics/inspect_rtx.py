"""
Authors: Andre L. Belem (https://github.com/andrebelem) and F.R.I.D.A.Y. (https://observatoriooceanografico.org/people/friday-bot)

Read-only command line diagnostic for Bruker RTX files: prints a bounded report of the RTX container and its TRT payload.

Verified structure (FINDINGS.md section 2): the RTX is an outer Windows-1252 XML file whose ``RTData``
element holds Base64 text of one zlib stream. The zlib stream decompresses to a second, XML-like
document (``CompData``) with no XML declaration. The decoder (``microxrf_to_netcdf.rtx_payload``) decodes both levels incrementally:

    file (1 MiB reads) -> Base64 (4-character groups) -> zlib.decompressobj -> chunks of bytes
    -> expat (Windows-1252 override) -> per-element summaries; large ``Data`` texts are Base64-decoded
    plane by plane.

Memory bound: the read chunk, the zlib output chunk, the expat buffer and at most ONE decoded image
plane (``Width x Height x ItemSize`` bytes) at a time. The encoded and decoded payloads are never held
whole and no copy of either is written to disk. Nothing is written anywhere and the input is opened
read-only. The full payload is never printed.

Usage (from the repository root):  python tools/diagnostics/inspect_rtx.py [path-to.rtx]
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from microxrf_to_netcdf.rtx_payload import IMAGE_SCALARS, inspect_rtx_payload  # noqa: E402

DEFAULT_RTX = Path("data/GRF17A_9-29cm_slab3_Elemental_map.rtx")


def print_report(report: dict[str, Any]) -> None:
    """Print a bounded, human-readable report (no payload content beyond a 64-character signature)."""
    print("RTX PAYLOAD REPORT")
    print(f"- file: {report['path']} ({report['file_bytes']} bytes)")
    if report["error"]:
        print(f"- ERROR: {report['error']}")
    outer = report.get("outer")
    if outer:
        print(f"- outer RTCompression: {outer['attributes']}; RTData text starts at byte {outer['prefix_bytes']}")
    decode = report["decode"]
    print(f"- Base64 characters: {decode['base64_chars']}; decompressed bytes: {decode['decoded_bytes']}")
    print(f"- zlib stream finished: {decode['zlib_finished']}; unused bytes after it: {decode['zlib_unused_bytes']}")
    print(f"- text after </RTData>: {decode['outer_tail']!r}")
    print(f"- first decompressed bytes: {decode['first_bytes_hex']}")
    print(f"- sha256 of decompressed payload: {decode['sha256_decoded']}")
    payload = report.get("payload")
    if payload:
        print(f"- payload signature text: {payload['signature_text']!r}")
        print(f"- payload root: {payload['root']}; elements: {payload['elements']}; max depth: {payload['max_depth']}")
        print(f"- ClassInstance types: {payload['class_instance_types']}")
        print(f"- distinct tag names ({len(payload['tag_vocabulary'])}; trailing digits shown as #): {payload['tag_vocabulary']}")
        for number, image in enumerate(payload["images"]):
            scalars = image["scalars"]
            print(f"\nIMAGE {number}: {image['name']!r} at decompressed offset {image['payload_offset']}")
            print("  " + "; ".join(f"{key}={scalars.get(key)}" for key in IMAGE_SCALARS))
            print(f"  planes decoded to Width x Height x ItemSize bytes, PlaneCount matches: {image['geometry_consistent']}")
            for index, plane in sorted(image["planes"].items()):
                stats = " ".join(f"{key}={plane[key]}" for key in ("dtype", "min", "max", "sum", "nonzero") if key in plane)
                print(f"  plane {index:>2} {plane.get('Description', '-')!s:8} Size={plane.get('Size')} "
                      f"decoded={plane.get('decoded_bytes')} data@{plane.get('data_payload_offset')} {stats}")
            counters = {index: plane["line_counter_values"] for index, plane in sorted(image["planes"].items())
                        if "line_counter_values" in plane}
            print(f"  LineCounter distinct values per plane (meaning unknown): {counters}")
            texts = [o["text"] for o in image["overlays"] if "text" in o]
            rects = [(o["name"], o["rect"]) for o in image["overlays"] if o["type"] == "TRTRectangleOverlayElement" and "rect" in o]
            print(f"  overlay texts: {texts}")
            print(f"  rectangle overlays: {rects}")
            used = {key: value for key, value in image["map_display"].items() if value.get("MapUsed") == "1"}
            if image["map_display"]:
                print(f"  map display entries used: {used}")
    print(f"\n- time: {report['seconds']} s; peak working set: {report['peak_working_set_mib']} MiB")


def main(argv: list[str]) -> int:
    path = Path(argv[1]) if len(argv) > 1 else DEFAULT_RTX
    if not path.is_file():
        print(f"RTX file not found: {path}", file=sys.stderr)
        return 2
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    report = inspect_rtx_payload(path)
    print_report(report)
    return 1 if report["error"] else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
