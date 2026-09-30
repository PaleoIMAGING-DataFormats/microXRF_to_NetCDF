"""
Authors: Andre L. Belem (https://github.com/andrebelem) and F.R.I.D.A.Y. (https://observatoriooceanografico.org/people/friday-bot)

Fast tests of microxrf_to_netcdf.rtx_payload on small synthetic RTX files (no original acquisition files).
The synthetic payload imitates only the structure observed in the real RTX (FINDINGS.md section 8).
"""

from __future__ import annotations

import base64
import zlib
from pathlib import Path

import numpy as np
import pytest

from microxrf_to_netcdf import rtx_payload as rtx

OUTER = ('<?xml version="1.0" encoding="WINDOWS-1252" standalone="yes"?>\r\n<TRTProject>\r\n  <RTHeader>\r\n'
         '  <ProjectHeader><Date>1.2.2026</Date><Time>3:04:05</Time></ProjectHeader>\r\n'
         '    <RTCompression compressor="{compressor}" encoder="{encoder}"/>\r\n  </RTHeader>\r\n'
         '  <RTData>{data}</RTData>\r\n</TRTProject>\r\n')


def b64(raw: bytes) -> str:
    return base64.b64encode(raw).decode("ascii")


def image_xml(name="Img", width=4, height=3, itemsize=2, planes=None, extra="") -> str:
    """One TRTImageData with the element layout seen in the real payload."""
    dtype = "<u2" if itemsize == 2 else "u1"
    planes = planes if planes is not None else [np.arange(width * height, dtype=dtype) * 7]
    body = "".join(
        f"<Plane{i}><Description>P{i}</Description><Valid>0</Valid>"
        f"<LineCounter>{','.join(['3'] * height)}</LineCounter>"
        f"<Data>{b64(np.asarray(p, dtype=dtype).tobytes())}</Data><Size>{width * height * itemsize}</Size></Plane{i}>"
        for i, p in enumerate(planes))
    return (f'<ClassInstance Type="TRTImageData" Name="{name}"><ItemSize>{itemsize}</ItemSize><Width>{width}</Width>'
            f'<Height>{height}</Height><PlaneCount>{len(planes)}</PlaneCount><MultiImage>1</MultiImage>'
            f'<XCalibration>12,5</XCalibration><YCalibration>12,5</YCalibration><Date>1.2.2026</Date><Time>3:04:05</Time>'
            f'{body}{extra}</ClassInstance>')


OVERLAY = ('<ClassInstance Type="TRTTextOverlayElement" Name="T"><TRTOverlayElement><Pos><PosX>1</PosX><PosY>2</PosY></Pos>'
           '<Rect><Left>-1</Left><Top>0</Top><Right>3</Right><Bottom>2</Bottom></Rect></TRTOverlayElement>'
           '<Text>HV: 50,0 kV \xb5m</Text></ClassInstance>')


def payload_bytes(*images: str) -> bytes:
    text = ('<CompData><ClassInstance Type="TRTProject" Name="Bruker project"><ChildClassInstances>'
            + "".join(images) + "</ChildClassInstances></ClassInstance></CompData>")
    return text.encode("windows-1252")


def write_rtx(path: Path, payload: bytes, compressor="zlib", encoder="base64", wrap=0) -> Path:
    data = b64(zlib.compress(payload))
    if wrap:
        data = "\r\n".join(data[i:i + wrap] for i in range(0, len(data), wrap))
    path.write_bytes(OUTER.format(compressor=compressor, encoder=encoder, data=data).encode("windows-1252"))
    return path


@pytest.fixture
def small_rtx(tmp_path):
    return write_rtx(tmp_path / "a.rtx", payload_bytes(image_xml(extra=OVERLAY)))


def test_authorship_line_is_exact():
    assert rtx.AUTHORS.startswith("Authors: Andre L. Belem (https://github.com/andrebelem) and F.R.I.D.A.Y.")
    assert rtx.AUTHORS in Path(rtx.__file__).read_text(encoding="utf-8")


def test_parse_decimal_comma():
    assert rtx.parse_decimal("31,308984675955") == pytest.approx(31.308984675955)
    assert rtx.parse_decimal("2") == 2.0


def test_decodes_synthetic_rtx(small_rtx):
    report = rtx.inspect_rtx_payload(small_rtx)
    assert report["error"] is None
    assert report["outer"]["attributes"] == {"compressor": "zlib", "encoder": "base64"}
    assert report["decode"]["zlib_finished"] is True and report["decode"]["zlib_unused_bytes"] == 0
    assert report["decode"]["outer_tail"].strip() == "</TRTProject>"
    assert report["decode"]["first_bytes_hex"].startswith("3c 43 6f 6d 70 44 61 74 61 3e")  # "<CompData>"
    payload = report["payload"]
    assert payload["root"] == "CompData" and payload["signature_text"].startswith("<CompData><ClassInstance")
    (image,) = payload["images"]
    assert image["name"] == "Img"
    assert image["scalars"]["Width"] == "4" and image["scalars"]["XCalibration"] == "12,5"
    assert image["geometry_consistent"] is True
    plane = image["planes"][0]
    assert plane["Description"] == "P0" and plane["decoded_bytes"] == 24 and plane["dtype"] == "<u2"
    assert (plane["min"], plane["max"], plane["sum"], plane["nonzero"]) == (0, 77, 7 * sum(range(12)), 11)
    assert plane["line_counter_values"] == {"3": 3}


def test_overlay_text_uses_windows_1252_and_rect_is_captured(small_rtx):
    (image,) = rtx.inspect_rtx_payload(small_rtx)["payload"]["images"]
    (overlay,) = image["overlays"]
    assert overlay["text"] == "HV: 50,0 kV \xb5m"
    assert overlay["rect"] == {"Left": -1, "Top": 0, "Right": 3, "Bottom": 2}
    assert overlay["pos"] == {"PosX": 1, "PosY": 2}


def test_plane_sink_receives_exact_values(tmp_path):
    values = np.array([[1, 2, 65535], [0, 300, 4]], dtype="<u2")
    path = write_rtx(tmp_path / "b.rtx", payload_bytes(image_xml(width=3, height=2, planes=[values, values + 1])))
    received = {}
    report = rtx.inspect_rtx_payload(path, lambda i, k, raw: received.__setitem__((i, k), np.frombuffer(raw, "<u2").reshape(2, 3)))
    assert report["error"] is None
    assert np.array_equal(received[(0, 0)], values) and np.array_equal(received[(0, 1)], values + 1)


def test_uint8_itemsize(tmp_path):
    plane = np.arange(6, dtype="u1")
    path = write_rtx(tmp_path / "c.rtx", payload_bytes(image_xml(width=3, height=2, itemsize=1, planes=[plane])))
    (image,) = rtx.inspect_rtx_payload(path)["payload"]["images"]
    assert image["planes"][0]["dtype"] == "|u1" and image["planes"][0]["max"] == 5


def test_result_is_independent_of_read_chunking_and_line_wrapping(tmp_path, monkeypatch):
    payload = payload_bytes(image_xml(width=40, height=30, planes=[np.arange(1200) % 251, np.arange(1200) % 97]),
                            image_xml(name="Second", extra=OVERLAY))
    reference = rtx.inspect_rtx_payload(write_rtx(tmp_path / "ref.rtx", payload))
    # tiny reads split <RTData>, </RTData>, Base64 groups and expat tokens at arbitrary places
    monkeypatch.setattr(rtx, "READ_CHUNK", 7)
    for wrap in (0, 76, 5):
        chunked = rtx.inspect_rtx_payload(write_rtx(tmp_path / f"w{wrap}.rtx", payload, wrap=wrap))
        assert chunked["error"] is None, wrap
        assert chunked["payload"] == reference["payload"]
        assert chunked["decode"]["sha256_decoded"] == reference["decode"]["sha256_decoded"]
        assert chunked["decode"]["decoded_bytes"] == len(payload)


def test_header_larger_than_first_read_is_reported_not_guessed(tmp_path):
    path = tmp_path / "d.rtx"
    path.write_bytes(b"<TRTProject><RTHeader>" + b" " * (rtx.HEADER_SCAN_BYTES + 10) + b"</RTHeader></TRTProject>")
    report = rtx.inspect_rtx_payload(path)
    assert "no <RTData>" in report["error"]


@pytest.mark.parametrize("compressor, encoder", [("lzma", "base64"), ("zlib", "hex"), ("", "")])
def test_unsupported_declared_compression_is_refused(tmp_path, compressor, encoder):
    path = write_rtx(tmp_path / "e.rtx", payload_bytes(image_xml()), compressor=compressor, encoder=encoder)
    report = rtx.inspect_rtx_payload(path)
    assert report["error"] and "unsupported RTCompression" in report["error"]
    assert "payload" not in report


def test_missing_rtdata_and_missing_compression(tmp_path):
    path = tmp_path / "f.rtx"
    path.write_bytes(b"<TRTProject><RTHeader/></TRTProject>")
    assert "no <RTData>" in rtx.inspect_rtx_payload(path)["error"]
    path.write_bytes(b"<TRTProject><RTHeader/><RTData>AAAA</RTData></TRTProject>")
    assert "no <RTCompression>" in rtx.inspect_rtx_payload(path)["error"]


def _replace_data(path: Path, data: str) -> None:
    text = path.read_bytes().decode("windows-1252")
    start, end = text.index("<RTData>") + 8, text.index("</RTData>")
    path.write_bytes((text[:start] + data + text[end:]).encode("windows-1252"))


def test_truncated_zlib_stream_is_an_error(small_rtx):
    stream = zlib.compress(payload_bytes(image_xml()))
    _replace_data(small_rtx, b64(stream[: len(stream) // 2]))
    report = rtx.inspect_rtx_payload(small_rtx)
    assert report["error"] is not None
    assert report["decode"]["zlib_finished"] is False


def test_bytes_after_the_zlib_stream_are_an_error(small_rtx):
    _replace_data(small_rtx, b64(zlib.compress(payload_bytes(image_xml())) + b"JUNK"))
    assert "follow the end of the zlib stream" in rtx.inspect_rtx_payload(small_rtx)["error"]


def test_invalid_base64_is_an_error(small_rtx):
    _replace_data(small_rtx, "!!!!" + b64(zlib.compress(b"<a/>")))
    assert "invalid Base64" in rtx.inspect_rtx_payload(small_rtx)["error"]


def test_base64_length_not_multiple_of_four_is_an_error(small_rtx):
    _replace_data(small_rtx, b64(zlib.compress(b"<a/>"))[:-1].rstrip("="))
    assert rtx.inspect_rtx_payload(small_rtx)["error"] is not None


def test_zlib_data_that_is_not_zlib_is_an_error(small_rtx):
    _replace_data(small_rtx, b64(b"this is not a zlib stream"))
    assert "zlib decompression failed" in rtx.inspect_rtx_payload(small_rtx)["error"]


def test_unclosed_rtdata_is_an_error(tmp_path):
    path = tmp_path / "g.rtx"
    data = b64(zlib.compress(b"<a/>"))
    complete = OUTER.format(compressor="zlib", encoder="base64", data=data)
    path.write_text(complete[: complete.index("</RTData>")], encoding="windows-1252")   # file ends inside RTData
    assert "no closed <RTData>" in rtx.inspect_rtx_payload(path)["error"]
    path.write_text(complete.replace("</RTData>", ""), encoding="windows-1252")           # tail text read as Base64
    assert rtx.inspect_rtx_payload(path)["error"] is not None


def test_malformed_payload_xml_is_reported(tmp_path):
    path = write_rtx(tmp_path / "h.rtx", b"<CompData><ClassInstance></CompData>")
    assert "not well-formed" in rtx.inspect_rtx_payload(path)["error"]


def test_bad_plane_base64_and_size_mismatch_are_recorded(tmp_path):
    xml = image_xml(width=2, height=2, planes=[np.zeros(4, "<u2")]).replace("<Size>8</Size>", "<Size>16</Size>")
    report = rtx.inspect_rtx_payload(write_rtx(tmp_path / "i.rtx", payload_bytes(xml)))
    (image,) = report["payload"]["images"]
    assert image["geometry_consistent"] is False        # Size says 16 bytes, 8 were decoded
    bad = image_xml(width=2, height=2, planes=[np.zeros(4, "<u2")]).replace("<Data>AAAAAAAAAAA=</Data>", "<Data>@@@@</Data>")
    (image,) = rtx.inspect_rtx_payload(write_rtx(tmp_path / "j.rtx", payload_bytes(bad)))["payload"]["images"]
    assert "decode_error" in image["planes"][0] and image["geometry_consistent"] is False


def test_tag_vocabulary_collapses_numbered_tags(tmp_path):
    payload = rtx.inspect_rtx_payload(write_rtx(tmp_path / "k.rtx", payload_bytes(
        image_xml(planes=[np.zeros(12, "<u2")] * 3))))["payload"]
    assert "Plane#" in payload["tag_vocabulary"] and "Plane0" not in payload["tag_vocabulary"]


def test_input_file_is_not_modified(small_rtx):
    before = (small_rtx.stat().st_size, small_rtx.stat().st_mtime_ns, small_rtx.read_bytes())
    rtx.inspect_rtx_payload(small_rtx)
    assert (small_rtx.stat().st_size, small_rtx.stat().st_mtime_ns, small_rtx.read_bytes()) == before


def test_many_planes_are_decoded_one_at_a_time(tmp_path):
    """The sink sees planes sequentially; the module itself keeps no plane bytes."""
    planes = [np.full(12, i, "<u2") for i in range(20)]
    path = write_rtx(tmp_path / "l.rtx", payload_bytes(image_xml(planes=planes)))
    seen = []
    report = rtx.inspect_rtx_payload(path, lambda i, k, raw: seen.append((k, len(raw))))
    assert seen == [(k, 24) for k in range(20)]
    assert all("bytes" not in str(type(v)) for v in report["payload"]["images"][0]["planes"][0].values())
