#!/usr/bin/env python3
"""Validate the committed checkout locally on Python 3.12 and 3.14."""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tarfile
from pathlib import Path
from tempfile import TemporaryDirectory

ROOT = Path(__file__).resolve().parents[1]
VERSIONS = ("3.12", "3.14")
PDF = Path("output/pdf/OPH-FEP-Synergetics_combined.pdf")


def run(*command: str | Path, cwd: Path, env: dict[str, str]) -> None:
    subprocess.run(command, cwd=cwd, env=env, stdin=subprocess.DEVNULL, check=True)


def check_version(version: str, archive: Path, work: Path) -> None:
    source = work / "source"
    source.mkdir(parents=True)
    with tarfile.open(archive) as snapshot:
        snapshot.extractall(source, filter="data")
    if not (source / PDF).is_file():
        raise RuntimeError(f"The committed manuscript PDF is missing: {PDF}")

    environment = work / "venv"
    env = os.environ.copy()
    for name in ("PYTHONHOME", "PYTHONPATH", "VIRTUAL_ENV", "UV_PROJECT", "UV_WORKING_DIR", "UV_WORKING_DIRECTORY"):
        env.pop(name, None)
    env.update(
        UV_PROJECT_ENVIRONMENT=str(environment),
        MPLBACKEND="Agg",
        MPLCONFIGDIR=str(work / "matplotlib"),
        COVERAGE_FILE=str(work / "coverage"),
    )
    binary_dir = environment / ("Scripts" if os.name == "nt" else "bin")
    suffix = ".exe" if os.name == "nt" else ""
    python = binary_dir / f"python{suffix}"
    cli = binary_dir / f"oph-fep-synergetics{suffix}"

    print(f"== Python {version}: locked, noneditable install ==", flush=True)
    run("uv", "sync", "--locked", "--python", version, "--no-editable", "--extra", "dev", "--extra", "figures",
        cwd=source, env=env)
    run(python, "--version", cwd=work, env=env)
    print(f"== Python {version}: installed source register ==", flush=True)
    run(cli, "sources", "--json", cwd=work, env=env)
    print(f"== Python {version}: tests and coverage ==", flush=True)
    run(python, "-m", "pytest", "tests/", "--cov", "--cov-report=term-missing", cwd=source, env=env)
    print(f"== Python {version}: executable probes and generated documents ==", flush=True)
    run(cli, "verify", cwd=work, env=env)
    run(python, "scripts/render_docs.py", "--check", cwd=source, env=env)
    print(f"Python {version} local CI passed", flush=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--python", choices=VERSIONS, help="Run one version. With no selection, run both.")
    args = parser.parse_args()
    for tool in ("git", "uv", "pdftotext"):
        if shutil.which(tool) is None:
            parser.error(f"{tool} is required on PATH")
    git = ("git", "--no-optional-locks")
    try:
        commit = subprocess.check_output((*git, "rev-parse", "--verify", "HEAD"), cwd=ROOT, text=True).strip()
        changes = subprocess.check_output((*git, "status", "--porcelain", "--untracked-files=no"), cwd=ROOT)
        if changes:
            raise RuntimeError("Commit the tracked changes before running the local CI gate.")
        versions = (args.python,) if args.python else VERSIONS
        print(f"Local CI for {commit}", flush=True)
        with TemporaryDirectory(prefix="oph-fep-local-ci-") as temporary:
            work = Path(temporary)
            archive = work / "source.tar"
            with archive.open("wb") as stream:
                subprocess.run((*git, "archive", commit), cwd=ROOT, stdout=stream, stdin=subprocess.DEVNULL, check=True)
            for version in versions:
                check_version(version, archive, work / version)
        print(f"Local CI passed for {commit} on Python {', '.join(versions)}", flush=True)
    except (RuntimeError, OSError, subprocess.CalledProcessError) as error:
        print(f"Local CI failed: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
