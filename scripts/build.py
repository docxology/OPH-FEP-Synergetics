#!/usr/bin/env python3
"""Build the manuscript end to end: figures, tokens, PDF, contract check.

One command that runs the four documented steps in order (README.md, "The
manuscript"):

1. ``scripts/make_figures.py`` -- the figures and their accessibility registry;
2. the token hydration, run by the template checkout's own interpreter, the one
   that passes the template's bytecode audit;
3. the template's ``stage_03_render.py --skip-manuscript-hydration`` (hydration
   already ran);
4. ``scripts/check_render.py``, which reads the render contract back out of the
   artefact.

The template checkout is ``--template-root``, else ``$TEMPLATE_ROOT``, else the
sibling ``template/`` directory beside this project. With ``--deliver`` the
fresh PDF is copied to ``~/Downloads/`` -- the acceptance convention is a
byte-identical copy there -- and both sha256 digests are printed.

Exit codes: 0 every stage passed; 1 unused; 2 a stage failed or no template
checkout was found.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import shutil
import subprocess
import sys
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[1]
PROJECT_NAME = "OPH-FEP-Synergetics"


def _template_root(explicit: str | None) -> Path:
    """Find the template checkout: flag, environment, then the sibling default."""
    candidates = [
        Path(explicit) if explicit else None,
        Path(os.environ["TEMPLATE_ROOT"]) if os.environ.get("TEMPLATE_ROOT") else None,
        PROJECT.parent / "template",
    ]
    for candidate in candidates:
        if candidate is None:
            continue
        python = candidate / ".venv" / "bin" / "python"
        stage = candidate / "scripts" / "pipeline" / "stage_03_render.py"
        if python.is_file() and stage.is_file():
            return candidate
    raise SystemExit(
        "no template checkout found: pass --template-root, set $TEMPLATE_ROOT, "
        "or place a template/ checkout beside this project"
    )


def _run(command: list[str], *, cwd: Path | None = None, env: dict[str, str] | None = None) -> None:
    """Run one stage; a non-zero exit stops the build with exit code 2."""
    completed = subprocess.run(command, cwd=cwd, env=env)
    if completed.returncode != 0:
        raise SystemExit(2)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--template-root", help="the docxology/template checkout")
    parser.add_argument(
        "--deliver",
        action="store_true",
        help="copy the fresh PDF to ~/Downloads and print both sha256 digests",
    )
    args = parser.parse_args(argv)

    template = _template_root(args.template_root)
    template_python = template / ".venv" / "bin" / "python"
    project_python = sys.executable

    print("== 1. figures and accessibility registry ==")
    _run([project_python, str(PROJECT / "scripts" / "make_figures.py")], cwd=PROJECT)

    print("== 2. token hydration with the template's interpreter ==")
    env = dict(os.environ, TEMPLATE_REPO_ROOT=str(template))
    _run(
        [str(template_python), str(PROJECT / "scripts" / "z_generate_manuscript_variables.py")],
        cwd=template,
        env=env,
    )

    print("== 3. the PDF ==")
    _run(
        [
            str(template_python),
            "scripts/pipeline/stage_03_render.py",
            "--project",
            f"working/{PROJECT_NAME}",
            "--skip-manuscript-hydration",
        ],
        cwd=template,
    )

    print("== 4. the render contract ==")
    _run([project_python, str(PROJECT / "scripts" / "check_render.py")], cwd=PROJECT)

    pdf = PROJECT / "output" / "pdf" / f"{PROJECT_NAME}_combined.pdf"
    print(f"rendered: {pdf}")
    if args.deliver:
        target = Path.home() / "Downloads" / pdf.name
        shutil.copyfile(pdf, target)
        rendered = hashlib.sha256(pdf.read_bytes()).hexdigest()
        delivered = hashlib.sha256(target.read_bytes()).hexdigest()
        print(f"sha256 {rendered}  {pdf}")
        print(f"sha256 {delivered}  {target}")
        if rendered != delivered:
            raise SystemExit(2)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
