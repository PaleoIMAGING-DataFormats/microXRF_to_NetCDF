"""
Authors: Andre L. Belem (https://github.com/andrebelem) and F.R.I.D.A.Y. (https://observatoriooceanografico.org/people/friday-bot)

Read-only diagnostic report for Bruker BCF and RTX acquisition files.
"""

from __future__ import annotations

import argparse
import importlib
import importlib.metadata
import inspect
import os
import re
import struct
import sys
import xml.sax
from collections.abc import Mapping
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from xml.sax.handler import ContentHandler


AUTHORS = "Authors: Andre L. Belem (https://github.com/andrebelem) and F.R.I.D.A.Y. (https://observatoriooceanografico.org/people/friday-bot)"
REFERENCE_PATTERN = re.compile(r"[^\x00\r\n<>|]{1,512}?\.(?:bcf|rtx|png|jpg|jpeg|tif|tiff|bmp|csv|txt|spx|msa|emsa|xml)\b", re.IGNORECASE)
BASE64_CHUNK_PATTERN = re.compile(r"^[A-Za-z0-9+/\s=]*$")
KEYWORDS = (
    "acquisition", "date", "time", "instrument", "model", "software", "version",
    "sample", "scan", "width", "height", "resolution", "pixel", "coordinate",
    "stage", "position", "element", "spectrum", "map", "image", "file", "path",
)


def local_name(tag: str) -> str:
    """Return an XML tag name without its namespace."""
    return tag.rsplit("}", 1)[-1]


def format_size(size: int) -> str:
    """Format byte counts without losing the exact value."""
    units = ("B", "KiB", "MiB", "GiB", "TiB")
    value = float(size)
    for unit in units:
        if value < 1024 or unit == units[-1]:
            return f"{size} bytes ({value:.2f} {unit})"
        value /= 1024
    return f"{size} bytes"


def file_signature(path: Path, count: int = 32) -> str:
    """Read a short signature only; never load an entire source file."""
    with path.open("rb") as source:
        return source.read(count).hex(" ").upper()


def inventory_files(data_dir: Path) -> list[Path]:
    """Print exact filesystem metadata and bounded file signatures."""
    print("\nFILE INVENTORY")
    files = sorted(path for path in data_dir.iterdir() if path.is_file())
    if not files:
        print("No files found.")
        return files
    for path in files:
        stat = path.stat()
        modified = datetime.fromtimestamp(stat.st_mtime, tz=timezone.utc).isoformat()
        print(f"- Name: {path.name}")
        print(f"  Size: {format_size(stat.st_size)}")
        print(f"  Modified (UTC): {modified}")
        print(f"  Signature (first 32 bytes): {file_signature(path)}")
    return files


def inspect_png(path: Path) -> dict[str, Any] | None:
    """Read PNG IHDR metadata without decoding image pixels."""
    with path.open("rb") as source:
        header = source.read(33)
    if len(header) < 33 or header[:8] != b"\x89PNG\r\n\x1a\n" or header[12:16] != b"IHDR":
        return None
    width, height, bit_depth, color_type, compression, filter_method, interlace = struct.unpack(
        ">IIBBBBB", header[16:29]
    )
    return {
        "width": width,
        "height": height,
        "bit_depth": bit_depth,
        "color_type": color_type,
        "compression": compression,
        "filter": filter_method,
        "interlace": interlace,
    }


def inspect_companion_images(files: list[Path]) -> None:
    """Report concise PNG metadata for the acquisition's companion images."""
    print("\nCOMPANION IMAGE HEADERS")
    found = False
    for path in files:
        if path.suffix.lower() != ".png":
            continue
        found = True
        details = inspect_png(path)
        if details is None:
            print(f"- {path.name}: not a valid PNG IHDR header.")
        else:
            print(
                f"- {path.name}: {details['width']} x {details['height']} pixels, "
                f"bit depth {details['bit_depth']}, color type {details['color_type']}."
            )
    if not found:
        print("No PNG files found.")


def short_text(value: str | None, limit: int = 240) -> str:
    """Normalize text for concise reporting without exposing large payloads."""
    normalized = " ".join((value or "").split())
    if len(normalized) > limit:
        return normalized[:limit] + " [truncated]"
    return normalized


def relevant_metadata(tag: str, attributes: dict[str, str], text: str) -> list[tuple[str, str]]:
    """Extract bounded candidate metadata values from an XML element."""
    candidates: list[tuple[str, str]] = []
    tag_lower = tag.lower()
    if text and any(keyword in tag_lower for keyword in KEYWORDS):
        candidates.append((tag, text))
    for name, value in attributes.items():
        label = f"{tag}@{name}"
        if any(keyword in label.lower() for keyword in KEYWORDS):
            candidates.append((label, short_text(value)))
    return candidates


class RTXHandler(ContentHandler):
    """Collect XML structure while keeping at most 240 text characters per element."""

    def __init__(self, report: dict[str, Any]) -> None:
        super().__init__()
        self.report = report
        self.stack: list[dict[str, Any]] = []

    def startElement(self, name: str, attrs: xml.sax.xmlreader.AttributesImpl) -> None:
        tag = local_name(name)
        if self.report["root"] is None:
            self.report["root"] = tag
        elif len(self.stack) == 1:
            self.report["major_sections"].append(tag)
        if self.stack:
            self.report["hierarchy"][(self.stack[-1]["tag"], tag)] += 1
        self.report["element_counts"][tag] += 1
        self.stack.append(
            {
                "tag": tag,
                "attributes": {local_name(key): value for key, value in attrs.items()},
                "text_length": 0,
                "preview": [],
                "base64_only": True,
            }
        )

    def characters(self, content: str) -> None:
        if not self.stack:
            return
        frame = self.stack[-1]
        frame["text_length"] += len(content)
        remaining = 240 - sum(map(len, frame["preview"]))
        if remaining > 0:
            frame["preview"].append(content[:remaining])
        if not BASE64_CHUNK_PATTERN.fullmatch(content):
            frame["base64_only"] = False

    def endElement(self, name: str) -> None:
        frame = self.stack.pop()
        tag = frame["tag"]
        preview = short_text("".join(frame["preview"]))
        if frame["text_length"] > 240:
            self.report["long_text_elements"][tag] += 1
        if frame["text_length"] >= 256 and frame["base64_only"]:
            self.report["base64_like_elements"][tag] += 1
        self.report["metadata"].extend(relevant_metadata(tag, frame["attributes"], preview))
        if frame["text_length"] <= 240 and (preview or frame["attributes"]):
            self.report["scalar_values"].append((tag, preview, frame["attributes"]))
            for value in (preview, *frame["attributes"].values()):
                for reference in REFERENCE_PATTERN.findall(value):
                    self.report["references"].add(short_text(reference.strip(' \"\'')))


def inspect_rtx(path: Path) -> dict[str, Any]:
    """Incrementally validate and summarize an RTX XML file in binary mode."""
    report: dict[str, Any] = {
        "well_formed": False,
        "error": None,
        "root": None,
        "major_sections": [],
        "element_counts": Counter(),
        "hierarchy": Counter(),
        "metadata": [],
        "scalar_values": [],
        "references": set(),
        "long_text_elements": Counter(),
        "base64_like_elements": Counter(),
        "nul_bytes": 0,
        "bytes_scanned": 0,
    }
    try:
        parser = xml.sax.make_parser()
        parser.setContentHandler(RTXHandler(report))
        with path.open("rb") as source:
            while chunk := source.read(1024 * 1024):
                parser.feed(chunk)
        parser.close()
        report["well_formed"] = True
    except xml.sax.SAXParseException as error:
        report["error"] = f"XML ParseError: line {error.getLineNumber()}, column {error.getColumnNumber()}: {error.getMessage()}"
    except (OSError, UnicodeError) as error:
        report["error"] = f"Read error: {type(error).__name__}: {error}"

    with path.open("rb") as source:
        while chunk := source.read(1024 * 1024):
            report["bytes_scanned"] += len(chunk)
            report["nul_bytes"] += chunk.count(b"\x00")
    return report


def print_rtx_report(path: Path, report: dict[str, Any]) -> None:
    """Print the XML analysis while separating verified facts from hypotheses."""
    print("\nRTX XML REPORT")
    print(f"- File: {path.name}")
    if report["well_formed"]:
        print("- Well-formed XML: verified (full incremental parse completed).")
    else:
        print(f"- Well-formed XML: not verified. {report['error']}")
    print(f"- Root element: {report['root'] or 'not found'}")
    sections = report["major_sections"]
    print(f"- Major root sections: {', '.join(sections) if sections else 'not found'}")
    print(f"- Binary scan: {report['nul_bytes']} NUL byte(s) in {report['bytes_scanned']} byte(s) scanned.")
    if report["nul_bytes"] or report["base64_like_elements"]:
        print("- Embedded binary indication: possible; see NUL/base64-like counts below.")
    else:
        print("- Embedded binary indication: not found by NUL-byte or base64-like text checks.")
    print("- Element frequencies:")
    for tag, count in report["element_counts"].most_common():
        print(f"  {tag}: {count}")
    print("- Nesting hierarchy (parent > child):")
    for (parent, child), count in sorted(report["hierarchy"].items()):
        print(f"  {parent} > {child}: {count}")
    print("- Long text elements (payload not printed):")
    if report["long_text_elements"]:
        for tag, count in report["long_text_elements"].most_common():
            print(f"  {tag}: {count}")
    else:
        print("  not found")
    print("- Base64-like text elements (payload not printed):")
    if report["base64_like_elements"]:
        for tag, count in report["base64_like_elements"].most_common():
            print(f"  {tag}: {count}")
    else:
        print("  not found")
    print("- Short scalar XML values and attributes:")
    if report["scalar_values"]:
        for tag, value, attributes in report["scalar_values"]:
            attributes_text = ", ".join(f"{name}={short_text(value)!r}" for name, value in attributes.items())
            suffix = f"; attributes: {attributes_text}" if attributes_text else ""
            print(f"  {tag}: {value or '[empty]'}{suffix}")
    else:
        print("  not found")
    print("- Candidate metadata (verified XML values; absent requested fields are not found):")
    seen: set[tuple[str, str]] = set()
    for label, value in report["metadata"]:
        item = (label, value)
        if item not in seen:
            seen.add(item)
            print(f"  {label}: {value}")
    if not seen:
        print("  not found")
    scalar_labels = {tag.lower() for tag, _, _ in report["scalar_values"]}
    print("- Requested acquisition fields:")
    requested_fields = {
        "Acquisition date/time": ("date" in scalar_labels or "time" in scalar_labels),
        "Instrument model": any("instrument" in label or "model" in label for label in scalar_labels),
        "Software version": any("software" in label or "version" in label for label in scalar_labels),
        "Sample identity": any("sample" in label or "specimen" in label for label in scalar_labels),
        "Scan dimensions": any(word in label for label in scalar_labels for word in ("width", "height", "dimension", "scan")),
        "Spatial resolution": any(word in label for label in scalar_labels for word in ("resolution", "pixel", "scale")),
        "Coordinate information": any(word in label for label in scalar_labels for word in ("coordinate", "position", "stage", "location")),
    }
    for field, found in requested_fields.items():
        print(f"  {field}: {'verified above' if found else 'not found'}")
    print("- Referenced files or formats:")
    if report["references"]:
        for reference in sorted(report["references"], key=str.lower):
            print(f"  {reference}")
    else:
        print("  not found")


def summarize_signal(signal: dict[str, Any], index: int) -> None:
    """Print only lightweight RosettaSciIO signal metadata."""
    data = signal.get("data")
    shape = getattr(data, "shape", None)
    dtype = getattr(data, "dtype", None)
    data_type = f"{type(data).__module__}.{type(data).__qualname__}"
    print(f"  Dataset {index}: shape={shape}, dtype={dtype}, array_type={data_type}")
    for axis in signal.get("axes", []):
        units = axis.get("units")
        if isinstance(units, str):
            units = units.replace("\u00c2\u00b5m", "um").replace("\u00b5m", "um")
        print(
            "    Axis: "
            f"name={axis.get('name')!r}, size={axis.get('size')}, "
            f"scale={axis.get('scale')}, offset={axis.get('offset')}, units={units!r}"
        )
    metadata = signal.get("metadata", {})
    original_metadata = signal.get("original_metadata", {})
    general = metadata.get("General", {}) if isinstance(metadata, Mapping) else {}
    signal_metadata = metadata.get("Signal", {}) if isinstance(metadata, Mapping) else {}
    title = general.get("title") if isinstance(general, Mapping) else None
    signal_type = signal_metadata.get("signal_type") if isinstance(signal_metadata, Mapping) else None
    if title or signal_type:
        print(f"    Dataset role: title={title!r}, signal_type={signal_type!r}")
    if getattr(data, "ndim", None) == 3 and any(axis.get("name") == "Energy" for axis in signal.get("axes", [])):
        print("    Dataset interpretation: lazy EDS spectrum image; data values were not materialized.")
    elif getattr(data, "ndim", None) == 2:
        print("    Dataset interpretation: two-dimensional image.")
    print(f"    Metadata keys: {', '.join(sorted(metadata)) or 'not found'}")
    print(f"    Original metadata keys: {', '.join(sorted(original_metadata)) or 'not found'}")
    scalar_values = metadata_scalar_values(metadata) + metadata_scalar_values(original_metadata, prefix=("original",))
    print("    Scalar metadata (bounded):")
    if scalar_values:
        for key, value in scalar_values[:50]:
            print(f"      {key}: {value}")
        if len(scalar_values) > 50:
            print(f"      [truncated: {len(scalar_values) - 50} additional scalar values]")
    else:
        print("      not found")


def metadata_scalar_values(value: Any, prefix: tuple[str, ...] = ()) -> list[tuple[str, str]]:
    """Flatten only small scalar BCF metadata values for a concise terminal report."""
    values: list[tuple[str, str]] = []
    if isinstance(value, Mapping):
        for key, child in value.items():
            values.extend(metadata_scalar_values(child, prefix + (str(key),)))
    elif isinstance(value, (str, int, float, bool)) and not isinstance(value, bytes):
        values.append((".".join(prefix), short_text(str(value))))
    return values


def normalize_date(value: str) -> str | None:
    """Normalize the known RTX and BCF date formats for comparison only."""
    for date_format in ("%Y-%m-%d", "%d.%m.%Y"):
        try:
            return datetime.strptime(value, date_format).date().isoformat()
        except ValueError:
            continue
    return None


def inspect_bcf(path: Path) -> dict[str, Any]:
    """Probe RosettaSciIO and request image-only lazy reading when supported."""
    result: dict[str, Any] = {"available": False, "version": None, "api": None, "error": None, "signals": []}
    try:
        result["version"] = importlib.metadata.version("rosettasciio")
    except importlib.metadata.PackageNotFoundError:
        result["error"] = "RosettaSciIO is not installed. Required command: python -m pip install rosettasciio"
        return result
    try:
        bruker = importlib.import_module("rsciio.bruker")
        reader = getattr(bruker, "file_reader")
        signature = inspect.signature(reader)
        result["available"] = True
        result["api"] = str(signature)
        parameters = signature.parameters
        kwargs: dict[str, Any] = {}
        if "lazy" in parameters:
            kwargs["lazy"] = True
        if "select_type" in parameters:
            kwargs["select_type"] = "images"
        result["read_request"] = kwargs
        signals = reader(str(path), **kwargs)
        result["signals"] = signals
    except Exception as error:  # Reader exceptions vary by optional backend and file compatibility.
        result["error"] = f"{type(error).__name__}: {error}"
    return result


def print_bcf_report(path: Path, result: dict[str, Any]) -> None:
    """Report the installed reader API and the outcome of its narrowest read."""
    print("\nBCF ROSETTASCIIO REPORT")
    print(f"- File: {path.name}")
    if not result["available"]:
        print(f"- Reader availability: not available. {result['error']}")
        return
    print(f"- RosettaSciIO version: {result['version']}")
    print(f"- Bruker reader API: rsciio.bruker.file_reader{result['api']}")
    print(f"- Initial read request: {result.get('read_request', {})} (image-only/lazy when supported).")
    if result["error"]:
        print(f"- Compatibility/read result: not verified. {result['error']}")
        return
    print(f"- Compatibility/read result: verified; returned {len(result['signals'])} dataset(s).")
    for index, signal in enumerate(result["signals"], start=1):
        summarize_signal(signal, index)
    elemental_map_titles = []
    for signal in result["signals"]:
        metadata = signal.get("metadata", {})
        general = metadata.get("General", {}) if isinstance(metadata, Mapping) else {}
        title = general.get("title") if isinstance(general, Mapping) else None
        if isinstance(title, str) and any(element in title.lower() for element in ("ca", "fe", "map")):
            elemental_map_titles.append(title)
    print(f"- Explicit elemental-map datasets: {', '.join(elemental_map_titles) if elemental_map_titles else 'not found in the returned BCF datasets'}")


def compare_metadata(rtx: dict[str, Any], bcf: dict[str, Any], files: list[Path]) -> None:
    """Compare only metadata proven to have been obtained by this run."""
    print("\nCROSS-FILE COMPARISON")
    bcf_path = next((path for path in files if path.suffix.lower() == ".bcf"), None)
    rtx_path = next((path for path in files if path.suffix.lower() == ".rtx"), None)
    bcf_name = bcf_path.name if bcf_path else None
    rtx_references = {reference.lower() for reference in rtx["references"]}
    if bcf_name and any(bcf_name.lower() in reference for reference in rtx_references):
        print("- RTX-to-BCF linkage: verified by an explicit RTX reference.")
    else:
        print("- RTX-to-BCF linkage: no explicit RTX reference found.")
    if bcf_path and rtx_path and bcf_path.stem == rtx_path.stem:
        print("- Filename evidence: BCF and RTX have the same filename stem.")
    if bcf.get("signals"):
        bcf_dates = {
            value for signal in bcf["signals"]
            for key, value in metadata_scalar_values(signal.get("metadata", {}))
            if key.endswith("General.date")
        }
        rtx_dates = {value for tag, value, _ in rtx.get("scalar_values", []) if tag == "Date"}
        if bcf_dates and rtx_dates:
            normalized_rtx_dates = {normalized for value in rtx_dates if (normalized := normalize_date(value))}
            normalized_bcf_dates = {normalized for value in bcf_dates if (normalized := normalize_date(value))}
            match = bool(normalized_rtx_dates & normalized_bcf_dates)
            outcome = "same calendar date" if match else "calendar-date match not verified"
            print(f"- Date evidence: RTX={', '.join(sorted(rtx_dates))}; BCF={', '.join(sorted(bcf_dates))}; {outcome} after format normalization.")
        print("- Same-acquisition assessment: probable from matching filename stem, calendar date, and complementary project/data roles; not conclusively verified by an explicit cross-reference.")
    else:
        print("- Metadata comparison: BCF metadata unavailable, so same-acquisition status is a hypothesis only (shared filename stem).")
    cafe = next((path for path in files if path.name.lower() == "cafe.png"), None)
    image_shapes = {tuple(signal.get("data").shape) for signal in bcf.get("signals", []) if getattr(signal.get("data"), "ndim", None) == 2}
    cafe_details = inspect_png(cafe) if cafe else None
    if cafe_details and (cafe_details["height"], cafe_details["width"]) in image_shapes:
        print("- CaFe.png comparison: dimensions match the BCF two-dimensional image grid (1800 x 240); pixel values were not compared.")
    else:
        print("- CaFe.png comparison: no matching BCF image grid verified.")
    print("- Video Mosaic.png comparison: dimensions do not match the BCF image grid; its pixel content was not compared.")
    print("- Duplicate images/metadata: not verified beyond the CaFe.png grid match.")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=Path(__file__).resolve().parents[2] / "data")
    arguments = parser.parse_args()
    data_dir = arguments.data_dir.resolve()
    print(AUTHORS)
    print("Read-only mode: no source data will be modified.")
    print(f"Python: {sys.version.split()[0]} ({sys.executable})")
    if not data_dir.is_dir():
        print(f"ERROR: Data directory not found: {data_dir}")
        return 2

    files = inventory_files(data_dir)
    inspect_companion_images(files)
    rtx_files = [path for path in files if path.suffix.lower() == ".rtx"]
    bcf_files = [path for path in files if path.suffix.lower() == ".bcf"]
    rtx = inspect_rtx(rtx_files[0]) if rtx_files else {"references": set()}
    if rtx_files:
        print_rtx_report(rtx_files[0], rtx)
    else:
        print("\nRTX XML REPORT\n- RTX file: not found")
    bcf = inspect_bcf(bcf_files[0]) if bcf_files else {"available": False, "error": "BCF file not found.", "signals": []}
    if bcf_files:
        print_bcf_report(bcf_files[0], bcf)
    else:
        print("\nBCF ROSETTASCIIO REPORT\n- BCF file: not found")
    compare_metadata(rtx, bcf, files)
    print("\nUNRESOLVED QUESTIONS AND RECOMMENDED NEXT STEPS")
    print("- Import into Hyper Fusion should begin with the BCF reader only after compatibility is verified above.")
    print("- Preserve the RTX as complementary project metadata; map its verified fields to Hyper Fusion import fields.")
    print("- Do not materialize a full hyperspectral cube until available memory and lazy-array behavior are confirmed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())