"""
Authors: Andre L. Belem (https://github.com/andrebelem) and F.R.I.D.A.Y. (https://observatoriooceanografico.org/people/friday-bot)

Fast tests of the bounded-memory BCF line framer and band decoder (microxrf_to_netcdf.bcf) on synthetic record streams.
They cover the pixel encodings the real file does not exercise (flag 0, flag 1, extra pulses), arbitrary block
boundaries, several band sizes, overflow protection and malformed or truncated streams.
"""

from __future__ import annotations

import struct

import numpy as np
import pytest
from rsciio.bruker import unbcf_fast

from microxrf_to_netcdf import bcf
from microxrf_to_netcdf.errors import BCFStreamError

import synthetic_data as sd


def decode_all(source, band_lines, dtype=np.uint16):
    bands = list(bcf.iter_bands(source, band_lines, dtype))
    return np.concatenate([b.data for b in bands]), bands


@pytest.fixture(scope="module")
def counts():
    return sd.random_counts(height=6, width=9, channels=64, seed=11, high=5)


def test_authorship_line_is_exact():
    text = (bcf.Path(bcf.__file__)).read_text(encoding="utf-8")
    assert "Authors: Andre L. Belem (https://github.com/andrebelem) and F.R.I.D.A.Y. (https://observatoriooceanografico.org/people/friday-bot)" in text


def test_encoder_and_decoder_agree_on_every_encoding(counts):
    """The decoded cube equals the counts that were encoded, for flags 0, 1, 2 and flag 2 with extra pulses."""
    source = sd.synthetic_source(counts)
    flags = {sd.cycling_flags(y, x) for y in range(counts.shape[0]) for x in range(counts.shape[1])}
    assert flags == {(0, 0), (1, 0), (2, 0), (2, 2)}
    decoded, _ = decode_all(source, 2)
    assert np.array_equal(decoded, counts)


@pytest.mark.parametrize("band_lines", [1, 2, 4, 6, 100])
def test_band_size_does_not_change_the_result(counts, band_lines):
    decoded, bands = decode_all(sd.synthetic_source(counts), band_lines)
    assert np.array_equal(decoded, counts)
    assert [b.first_line for b in bands] == list(range(0, counts.shape[0], min(band_lines, counts.shape[0])))
    assert all(b.data.shape[0] <= band_lines for b in bands)


@pytest.mark.parametrize("block_size", [1, 7, 100, 4064, 10**9])
def test_arbitrary_block_boundaries(counts, block_size):
    decoded, _ = decode_all(sd.synthetic_source(counts, block_size=block_size), 3)
    assert np.array_equal(decoded, counts)


def test_framer_plus_decoder_equals_the_unmodified_decoder_on_the_whole_stream(counts):
    stream = sd.build_stream(counts)
    whole = unbcf_fast.parse_to_numpy(sd.FixedBlockFile(stream), counts.shape, np.uint16)
    decoded, _ = decode_all(sd.synthetic_source(counts), 2)
    assert np.array_equal(decoded, whole)
    assert np.array_equal(whole, counts)


def test_line_records_partition_the_stream_exactly(counts):
    stream = sd.build_stream(counts)
    framer = bcf.frame_lines(sd.split_blocks(stream, 50), counts.shape[0], counts.shape[1])
    header = next(framer)
    records = list(framer)
    assert len(header) == bcf.HEADER_BYTES and len(records) == counts.shape[0]
    assert len(header) + sum(len(r.data) for r in records) == len(stream)
    assert all(r.pixels == counts.shape[1] for r in records)


def test_uint16_wraparound_is_a_real_hazard_that_uint32_decoding_avoids():
    """A count of 70000 does not fit uint16. The compiled decoder wraps it silently, to a small value that no
    'is the maximum near the ceiling' test would notice; uint32 decoding is exact. The converter therefore scans
    in uint32 and cross-checks the uint16 passes exactly (tests/test_convert_synthetic.py)."""
    counts = sd.random_counts(height=3, width=4, channels=32, seed=5, big={(1, 2, 7): 70000, (2, 0, 3): 40000})
    stream = sd.build_stream(counts, flags=lambda y, x: (2, 0))
    framer = bcf.frame_lines(sd.split_blocks(stream, 64), 3, 4)
    header = next(framer)
    records = [r.data for r in framer]
    wrapped = bcf.decode_band(header, records, 4, 32, np.uint16)
    assert wrapped[1, 2, 7] == 70000 - 65536, "the hazard is real: uint16 wraps silently"
    source = sd.synthetic_source(counts, sum_spectrum=counts.sum(axis=(0, 1)).astype(np.uint64))
    silent, bands = decode_all(source, 1)          # uint16 with the best-effort widening guard
    assert int(silent.max()) < 70000, "the best-effort guard cannot see a wrap to a small value"
    assert bands[2].widened and int(bands[2].data.max()) == 40000    # a count >= 2**15 is caught and widened
    exact, exact_bands = decode_all(source, 1, np.uint32)
    assert exact.dtype == np.uint32 and int(exact.max()) == 70000 and np.array_equal(exact, counts)
    assert not any(b.widened for b in exact_bands)


def test_valid_zero_counts_survive_decoding():
    counts = np.zeros((2, 5, 16), dtype=np.int64)
    counts[0, 1, 3] = 2
    decoded, _ = decode_all(sd.synthetic_source(counts), 2)
    assert int(decoded.sum()) == 2 and np.array_equal(decoded, counts)


def test_invalid_arguments():
    source = sd.synthetic_source(sd.random_counts())
    with pytest.raises(ValueError):
        list(bcf.iter_bands(source, 0))
    with pytest.raises(ValueError):
        list(bcf.iter_bands(source, 1, np.uint8))


# ------------------------------------------------------------------------------------------------
# Malformed and truncated streams
# ------------------------------------------------------------------------------------------------


@pytest.fixture()
def good_stream(counts):
    return sd.build_stream(counts)


def frame(stream, height=6, width=9, block=64):
    return list(bcf.frame_lines(sd.split_blocks(stream, block), height, width))


@pytest.mark.parametrize("cut", [10, bcf.HEADER_BYTES + 2, 1000, -1, -50])
def test_truncated_stream_fails_clearly(good_stream, cut):
    with pytest.raises(BCFStreamError, match="ended|truncated"):
        frame(good_stream[:cut] if cut > 0 else good_stream[:cut])


def test_trailing_bytes_are_refused(good_stream):
    with pytest.raises(BCFStreamError, match="unexpected"):
        frame(good_stream + b"\x00\x01\x02")


def test_extra_block_after_the_last_line_is_refused(good_stream):
    blocks = sd.split_blocks(good_stream, 64) + [b""]
    with pytest.raises(BCFStreamError):
        list(bcf.frame_lines(iter(blocks + [b"x"]), 6, 9))


def test_header_dimension_mismatch_is_refused(good_stream):
    with pytest.raises(BCFStreamError, match="height"):
        frame(good_stream, height=7)
    with pytest.raises(BCFStreamError, match="width"):
        frame(good_stream, width=8)


def test_empty_grid_is_refused():
    with pytest.raises(BCFStreamError, match="empty grid"):
        list(bcf.frame_lines([sd.stream_header(0, 5)], None, None))


def _first_pixel_offset():
    return bcf.HEADER_BYTES + 4


def _patch(stream, offset, fmt, value):
    data = bytearray(stream)
    struct.pack_into(fmt, data, offset, value)
    return bytes(data)


def test_unknown_flag_is_refused(good_stream):
    # pixel header: pixel_x(4) c1(2) c2(2) skip(4) flag(2) ...
    with pytest.raises(BCFStreamError, match="unknown pixel flag"):
        frame(_patch(good_stream, _first_pixel_offset() + 12, "<H", 7))


def test_pixel_x_beyond_width_is_refused(good_stream):
    with pytest.raises(BCFStreamError, match="pixel_x"):
        frame(_patch(good_stream, _first_pixel_offset(), "<I", 9))


def test_more_pixels_than_width_is_refused(good_stream):
    with pytest.raises(BCFStreamError, match="more than the width"):
        frame(_patch(good_stream, bcf.HEADER_BYTES, "<I", 10))


def test_oversized_pixel_payload_is_refused_without_reading_it(good_stream):
    with pytest.raises(BCFStreamError, match="exceeds"):
        frame(_patch(good_stream, _first_pixel_offset() + 18, "<I", 2**31))


def test_flag_zero_payload_shorter_than_its_pulses_is_refused():
    record = struct.pack("<II", 1, 3) + bytes(bcf.HEADER_BYTES - 8)
    record += struct.pack("<I", 1) + sd.PIXEL_HEADER.pack(0, 0, 0, 0, 0, 0, 5, 4) + bytes(4)
    with pytest.raises(BCFStreamError, match="shorter"):
        list(bcf.frame_lines([record], 1, 3))


def test_flag_one_payload_too_short_is_refused():
    record = struct.pack("<II", 1, 3) + bytes(bcf.HEADER_BYTES - 8)
    record += struct.pack("<I", 1) + sd.PIXEL_HEADER.pack(0, 0, 0, 0, 1, 0, 8, 6) + bytes(6)
    with pytest.raises(BCFStreamError, match="too short"):
        list(bcf.frame_lines([record], 1, 3))


def test_instructed_payload_smaller_than_its_trailer_is_refused():
    record = struct.pack("<II", 1, 3) + bytes(bcf.HEADER_BYTES - 8)
    record += struct.pack("<I", 1) + sd.PIXEL_HEADER.pack(0, 0, 0, 0, 2, 0, 0, 2) + bytes(8)
    with pytest.raises(BCFStreamError, match="trailer"):
        list(bcf.frame_lines([record], 1, 3))


def test_a_line_over_the_configured_limit_is_refused(good_stream):
    with pytest.raises(BCFStreamError, match="exceeds"):
        list(bcf.frame_lines(sd.split_blocks(good_stream, 64), 6, 9, max_line_bytes=50))


def test_a_corrupted_length_never_reads_the_rest_of_the_stream_into_memory():
    """A length field of 3 MiB is refused after reading only that much, not the (huge) remainder."""
    consumed = []

    def blocks():
        yield struct.pack("<II", 1, 3) + bytes(bcf.HEADER_BYTES - 8)
        yield struct.pack("<I", 1) + sd.PIXEL_HEADER.pack(0, 0, 0, 0, 2, 0, 0, 5 << 20)
        for _ in range(10**6):
            consumed.append(1)
            yield bytes(4064)

    with pytest.raises(BCFStreamError):
        list(bcf.frame_lines(blocks(), 1, 3))
    assert len(consumed) == 0


def test_rosettasciio_version_gate(monkeypatch):
    assert bcf.check_rosettasciio_version().startswith("0.14")
    monkeypatch.setattr("importlib.metadata.version", lambda name: "0.15.0")
    with pytest.raises(Exception, match="verified range"):
        bcf.check_rosettasciio_version()
