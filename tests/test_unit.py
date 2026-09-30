"""
Authors: Andre L. Belem (https://github.com/andrebelem) and F.R.I.D.A.Y. (https://observatoriooceanografico.org/people/friday-bot)

Fast unit tests: synthetic inputs only, no original acquisition files.
"""

from __future__ import annotations

import struct
import zlib
from pathlib import Path

import inspect_bruker as ib


def _png_bytes(width: int, height: int) -> bytes:
    def chunk(kind: bytes, payload: bytes) -> bytes:
        body = kind + payload
        return struct.pack(">I", len(payload)) + body + struct.pack(">I", zlib.crc32(body))

    ihdr = struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0)
    return b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", ihdr) + chunk(b"IEND", b"")


def test_authorship_line_is_exact():
    assert ib.AUTHORS == (
        "Authors: Andre L. Belem (https://github.com/andrebelem) and "
        "F.R.I.D.A.Y. (https://observatoriooceanografico.org/people/friday-bot)"
    )


def test_local_name_strips_namespace():
    assert ib.local_name("{urn:x}Tag") == "Tag"
    assert ib.local_name("Tag") == "Tag"


def test_format_size_keeps_exact_bytes():
    assert ib.format_size(617320728).startswith("617320728 bytes (588.72 MiB)")


def test_normalize_date_known_formats():
    assert ib.normalize_date("30.7.2026") == "2026-07-30"
    assert ib.normalize_date("2026-07-30") == "2026-07-30"
    assert ib.normalize_date("July 30") is None


def test_inspect_png_reads_header_only(tmp_path: Path):
    path = tmp_path / "x.png"
    path.write_bytes(_png_bytes(1800, 240))
    details = ib.inspect_png(path)
    assert (details["width"], details["height"]) == (1800, 240)
    assert details["bit_depth"] == 8 and details["color_type"] == 2


def test_inspect_png_rejects_non_png(tmp_path: Path):
    path = tmp_path / "x.png"
    path.write_bytes(b"not a png" * 10)
    assert ib.inspect_png(path) is None


def test_inspect_rtx_on_synthetic_project(tmp_path: Path):
    path = tmp_path / "s.rtx"
    path.write_bytes(
        b'<?xml version="1.0" encoding="WINDOWS-1252"?>'
        b'<TRTProject><RTHeader><ProjectHeader><Date>1.2.2020</Date></ProjectHeader>'
        b'<RTCompression compressor="zlib" encoder="base64"/></RTHeader>'
        b"<RTData>" + b"QUJD" * 100 + b"</RTData></TRTProject>"
    )
    report = ib.inspect_rtx(path)
    assert report["well_formed"] is True
    assert report["root"] == "TRTProject"
    assert report["major_sections"] == ["RTHeader", "RTData"]
    assert report["base64_like_elements"]["RTData"] == 1
    assert report["nul_bytes"] == 0


def test_inspect_rtx_reports_malformed_xml(tmp_path: Path):
    path = tmp_path / "bad.rtx"
    path.write_bytes(b"<a><b></a>")
    report = ib.inspect_rtx(path)
    assert report["well_formed"] is False
    assert "XML ParseError" in report["error"]


def test_metadata_scalar_values_flattens_only_scalars():
    flat = dict(ib.metadata_scalar_values({"a": {"b": 1, "c": [1, 2]}, "d": "x"}))
    assert flat == {"a.b": "1", "d": "x"}
