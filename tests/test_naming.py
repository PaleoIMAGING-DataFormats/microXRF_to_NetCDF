"""
Authors: Andre L. Belem (https://github.com/andrebelem) and F.R.I.D.A.Y. (https://observatoriooceanografico.org/people/friday-bot)

Guards the project identity: the package is ``microxrf_to_netcdf`` (src layout), the tools live under ``tools/``,
and no file of the working tree carries an earlier experimental project name. The earlier names are assembled from
fragments here so that this file does not contain them.
"""

import re
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
PACKAGE = ROOT / "src" / "microxrf_to_netcdf"
EARLIER_NAMES = re.compile("|".join([
    "test" + "_micro" + "ct", "micro" + "ct", "bruker" + r"[\s_-]?" + "spy", "bruc" + "ker" + "_spy",
]), re.IGNORECASE)
SKIPPED_PARTS = {"data", "__pycache__", ".pytest_cache", ".git", ".venv"}


def _text_files():
    for path in sorted(ROOT.rglob("*")):
        relative = path.relative_to(ROOT)
        if (not path.is_file() or SKIPPED_PARTS & set(relative.parts) or relative.parts[0] == "output"
                or any(part.endswith(".egg-info") for part in relative.parts)):
            continue
        yield path


def test_layout():
    assert (PACKAGE / "__init__.py").is_file() and (PACKAGE / "__main__.py").is_file()
    assert (ROOT / "pyproject.toml").is_file()
    for script in ("inspect_bruker.py", "inspect_rtx.py"):
        assert (ROOT / "tools" / "diagnostics" / script).is_file()
    for script in ("probe_bcf_streaming.py", "benchmark_conversion.py"):
        assert (ROOT / "tools" / "benchmarks" / script).is_file()
    assert not list(ROOT.glob("*.py")), "scripts belong under tools/ or tests/"


def test_no_text_file_carries_an_earlier_project_name():
    offenders = []
    for path in _text_files():
        try:
            text = path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            continue  # binary artifact
        for number, line in enumerate(text.splitlines(), 1):
            if EARLIER_NAMES.search(line):
                offenders.append(f"{path.relative_to(ROOT).as_posix()}:{number}: {line.strip()[:100]}")
    assert not offenders, "\n".join(offenders)


def test_no_file_or_directory_is_named_after_an_earlier_project_name():
    assert not [p for p in ROOT.rglob("*") if EARLIER_NAMES.search(p.name) and "data" not in p.relative_to(ROOT).parts]


def test_no_compatibility_alias_for_an_earlier_package_name():
    alias = "bruker" + "_spy"
    failed = subprocess.run([sys.executable, "-c", f"import {alias}"], cwd=ROOT, capture_output=True, text=True)
    assert failed.returncode != 0


def test_error_base_class_carries_the_project_name():
    from microxrf_to_netcdf import errors

    assert hasattr(errors, "MicroXRFToNetCDFError")


def test_cli_entry_point():
    result = subprocess.run([sys.executable, "-m", "microxrf_to_netcdf", "--version"], cwd=ROOT, capture_output=True,
                            text=True, timeout=120)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip().startswith("microxrf_to_netcdf ")


def test_documents_and_environment_use_the_project_name():
    assert (ROOT / "environment.yml").read_text(encoding="utf-8").startswith("name: microxrf_to_netcdf")
    for name in ("README.md", "CLAUDE.md", "FINDINGS.md", "TODO.md", "NETCDF_SCHEMA.md"):
        first = (ROOT / name).read_text(encoding="utf-8").splitlines()[0]
        assert "microXRF to NetCDF" in first, name
    assert "python -m microxrf_to_netcdf plan" in (ROOT / "README.md").read_text(encoding="utf-8")


def test_schema_attributes_are_documented_under_the_project_name():
    text = (ROOT / "NETCDF_SCHEMA.md").read_text(encoding="utf-8")
    assert "microxrf_to_netcdf_schema_version" in text
