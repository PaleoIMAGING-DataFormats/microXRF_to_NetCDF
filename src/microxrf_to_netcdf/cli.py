"""
Authors: Andre L. Belem (https://github.com/andrebelem) and F.R.I.D.A.Y. (https://observatoriooceanografico.org/people/friday-bot)

Command line: ``python -m microxrf_to_netcdf {plan,convert,validate}``. Run inside the ``microxrf_to_netcdf`` Conda environment.
Exit status: 0 success, 1 conversion or validation failure, 2 usage or preflight error.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from . import __version__
from .config import ConversionConfig
from .errors import MicroXRFToNetCDFError, PreflightError


def _log(message: str) -> None:
    print(message, file=sys.stderr, flush=True)


def _add_source_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--bcf", type=Path, required=True, help="input .bcf (read-only)")
    parser.add_argument("--rtx", type=Path, required=True, help="input .rtx (read-only)")


def _add_config_options(parser: argparse.ArgumentParser) -> None:
    defaults = ConversionConfig()
    parser.add_argument("--band-lines", type=int, default=defaults.band_lines,
                        help="scan lines decoded per band; multiple of the counts chunk along y (default %(default)s)")
    parser.add_argument("--chunks", type=int, nargs=3, metavar=("Y", "X", "ENERGY"), default=defaults.counts_chunks,
                        help="HDF5 chunk shape of the counts (default %(default)s)")
    parser.add_argument("--complevel", type=int, default=defaults.complevel, help="zlib level 0-9 (default %(default)s)")
    parser.add_argument("--no-shuffle", action="store_true", help="disable the HDF5 shuffle filter")
    parser.add_argument("--decode-dtype", choices=("uint16", "uint32"), default=defaults.decode_dtype)
    parser.add_argument("--counts-dtype", choices=("auto", "uint8", "uint16", "uint32"), default=defaults.counts_dtype,
                        help="output dtype of the counts; auto picks the smallest that holds the observed maximum")
    parser.add_argument("--allow-low-disk", action="store_true", help="override an insufficient-disk preflight failure")
    parser.add_argument("--allow-low-memory", action="store_true", help="override the memory preflight failure")


def _config(args: argparse.Namespace) -> ConversionConfig:
    return ConversionConfig(
        band_lines=args.band_lines, counts_chunks=tuple(args.chunks), complevel=args.complevel,
        shuffle=not args.no_shuffle, decode_dtype=args.decode_dtype, counts_dtype=args.counts_dtype,
        allow_low_disk=args.allow_low_disk, allow_low_memory=args.allow_low_memory,
        overwrite=getattr(args, "overwrite", False), verify_readback=not getattr(args, "no_verify", False))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="microxrf_to_netcdf", description=__doc__.strip().splitlines()[0])
    parser.add_argument("--version", action="version", version=f"microxrf_to_netcdf {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    plan = sub.add_parser("plan", help="inspect the inputs and print the resource preflight; writes nothing")
    _add_source_options(plan)
    plan.add_argument("--out", type=Path, required=True, help="intended output .nc (its filesystem is checked)")
    _add_config_options(plan)

    convert = sub.add_parser("convert", help="convert one paired BCF/RTX acquisition into one NetCDF-4 file")
    _add_source_options(convert)
    convert.add_argument("--out", type=Path, required=True, help="output .nc (must not exist unless --overwrite)")
    convert.add_argument("--overwrite", action="store_true")
    convert.add_argument("--no-verify", action="store_true", help="skip the read-back comparison pass")
    convert.add_argument("--quiet", action="store_true")
    _add_config_options(convert)

    validate = sub.add_parser("validate", help="validate a converted file against the sources and reference decoders")
    validate.add_argument("--nc", type=Path, required=True)
    _add_source_options(validate)
    validate.add_argument("--data-dir", type=Path, help="directory with 'Video Mosaic.png' and 'CaFe.png' (optional)")
    validate.add_argument("--strict", action="store_true", help="also decode every band as uint16 and uint32")
    validate.add_argument("--json", type=Path, help="write the check list to this file")

    args = parser.parse_args(argv)
    from .convert import convert as run_convert
    try:
        if args.command in ("plan", "convert"):
            config = _config(args)
            log = _log if not getattr(args, "quiet", False) else (lambda _: None)
            if args.command == "plan":
                print(json.dumps(run_convert(args.bcf, args.rtx, args.out, config, log, dry_run=True), indent=2, default=str))
                return 0
            report = run_convert(args.bcf, args.rtx, args.out, config, log)
            print(json.dumps({k: report[k] for k in ("output", "output_size_bytes", "output_sha256", "conversion_seconds",
                                                     "peak_working_set_mib", "counts_dtype", "compression",
                                                     "actual_layout")}, indent=2))
            return 0
        from .validate import validate_against_sources
        checks = validate_against_sources(args.nc, args.bcf, args.rtx, args.data_dir, args.strict, _log)
        for check in checks:
            print(f"[{check['status']}] {check['check']}: {check['detail']}")
        if args.json:
            args.json.write_text(json.dumps(checks, indent=2, ensure_ascii=False), encoding="utf-8")
        return 0 if all(c["status"] == "passed" for c in checks) else 1
    except PreflightError as error:
        print(f"PREFLIGHT FAILED (nothing was written): {error}", file=sys.stderr)
        return 2
    except MicroXRFToNetCDFError as error:
        print(f"FAILED: {type(error).__name__}: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
