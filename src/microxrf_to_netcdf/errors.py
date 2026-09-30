"""
Authors: Andre L. Belem (https://github.com/andrebelem) and F.R.I.D.A.Y. (https://observatoriooceanografico.org/people/friday-bot)

Exception hierarchy. Every failure of the converter is one of these, with a message that says what was
observed. None of them is ever swallowed: the converter removes its partial output and re-raises.
"""

from __future__ import annotations


class MicroXRFToNetCDFError(Exception):
    """Base class of all microXRF to NetCDF errors."""


class BCFStreamError(MicroXRFToNetCDFError, ValueError):
    """The BCF spectrum stream is truncated, inconsistent, or uses an unverified encoding."""


class RTXFormatError(MicroXRFToNetCDFError, ValueError):
    """The RTX does not have the structure that was verified (see FINDINGS.md section 8)."""


class GridMismatchError(MicroXRFToNetCDFError, ValueError):
    """The BCF and RTX spatial grids or calibrations do not correspond."""


class ValidationError(MicroXRFToNetCDFError, ValueError):
    """A scientific cross-check failed (for example the RTX video plane differs from the BCF video)."""


class PreflightError(MicroXRFToNetCDFError, RuntimeError):
    """Resources are insufficient or unknown, or the destination is unsafe. Nothing was written."""


class ConversionError(MicroXRFToNetCDFError, RuntimeError):
    """The conversion failed after it started. The partial output was removed."""
