"""
Authors: Andre L. Belem (https://github.com/andrebelem) and F.R.I.D.A.Y. (https://observatoriooceanografico.org/people/friday-bot)

Conversion parameters with documented defaults. The reasons for the defaults are in NETCDF_SCHEMA.md
(chunking) and FINDINGS.md sections 7 and 9 (memory and timings).
"""

from __future__ import annotations

from dataclasses import asdict, dataclass


@dataclass
class ConversionConfig:
    """Every parameter that changes what is written or how much memory is used."""

    # Decoding
    band_lines: int = 4                      # scan lines decoded per band; must be a multiple of counts_chunks[0]
    decode_dtype: str = "uint16"             # "uint16" (default) or "uint32"; wrap-around is guarded (see bcf.py)
    counts_dtype: str = "auto"               # "auto" picks the smallest of uint8/uint16/uint32 holding the maximum
    max_line_bytes: int = 512 << 20          # a scan-line record above this is treated as corruption

    # HDF5 layout (chunks are clipped to the dimension sizes)
    counts_chunks: tuple[int, int, int] = (2, 30, 4096)     # (y, x, energy); chosen from measurements, NETCDF_SCHEMA.md
    map_chunks: tuple[int, int, int] = (1, 120, 900)        # (element, y, x)
    video_chunks: tuple[int, int] = (120, 900)              # (y, x)
    mosaic_chunks: tuple[int, int, int] = (1, 256, 1024)    # (plane, mosaic_y, mosaic_x)
    complevel: int = 4                       # zlib level 1-9; 0 disables compression
    shuffle: bool = True                     # HDF5 byte-shuffle filter (only with compression)

    # Resource preflight
    disk_margin_fraction: float = 0.10       # extra free space required on top of the worst-case size
    disk_margin_bytes: int = 512 << 20       # ... but at least this much
    allow_low_disk: bool = False             # explicit override when free space cannot cover the worst case
    memory_fraction: float = 0.80            # estimated peak must not exceed this fraction of free RAM
    allow_low_memory: bool = False           # explicit override of the memory check

    # Output policy
    overwrite: bool = False                  # replace an existing final file (never an input, never in data/)
    verify_readback: bool = True             # third pass: compare every written band and plane with the source
    check_registration: bool = True          # correlate the BCF video with the mosaic crop of the Map footprint
    registration_max_crop_bytes: int = 256 << 20

    def validate(self) -> None:
        if self.band_lines < 1:
            raise ValueError("band_lines must be at least 1")
        if self.decode_dtype not in ("uint16", "uint32"):
            raise ValueError("decode_dtype must be 'uint16' or 'uint32'")
        if self.counts_dtype not in ("auto", "uint8", "uint16", "uint32"):
            raise ValueError("counts_dtype must be auto, uint8, uint16 or uint32")
        if not 0 <= self.complevel <= 9:
            raise ValueError("complevel must be between 0 and 9")
        for name in ("counts_chunks", "map_chunks", "video_chunks", "mosaic_chunks"):
            if any(int(v) < 1 for v in getattr(self, name)):
                raise ValueError(f"{name} entries must be positive")
        if not 0 < self.memory_fraction <= 1:
            raise ValueError("memory_fraction must be in (0, 1]")

    def as_dict(self) -> dict:
        return asdict(self)
