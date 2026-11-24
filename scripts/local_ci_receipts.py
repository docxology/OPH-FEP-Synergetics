#!/usr/bin/env python3
# OPH-FEP-Synergetics local CI pre-push hook
"""Record local CI passes and guard pushes into refs/heads/main."""

from __future__ import annotations

import json
import os
import re
import subprocess
import sys
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Iterable

PYTHON_VERSION = "3.14"
SCHEMA = 1
HOOK_MARKER = b"# OPH-FEP-Synergetics local CI pre-push hook"
RUN_GATE = f"uv run --no-project --python {PYTHON_VERSION} python scripts/check_local_ci.py"


def _receipt_path(root: Path, commit: str) -> Path:
    path = root.resolve() / "output" / "local-ci" / commit / f"{PYTHON_VERSION}.json"
    if path.resolve() != path:
        raise ValueError(f"Expected an ordinary local CI receipt path: {path}")
    return path


def _atomic_write(path: Path, data: bytes, mode: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile(dir=path.parent, prefix=f".{path.name}-", delete=False) as stream:
        temporary = Path(stream.name)
        try:
            stream.write(data)
            stream.close()
            temporary.chmod(mode)
            os.replace(temporary, path)
        finally:
            temporary.unlink(missing_ok=True)


def is_gate_version(runtime: object) -> bool:
    return isinstance(runtime, str) and re.fullmatch(rf"{re.escape(PYTHON_VERSION)}\.\d+", runtime) is not None


def _is_commit(commit: str) -> bool:
    return re.fullmatch(r"(?:[0-9a-f]{40}|[0-9a-f]{64})", commit) is not None and bool(commit.strip("0"))


def record_pass(root: Path, commit: str, runtime: str) -> Path:
    if not _is_commit(commit) or not is_gate_version(runtime):
        raise ValueError("A local CI pass requires a full commit ID and the requested Python runtime.")
    payload = {"schema": SCHEMA, "commit": commit, "python": runtime, "result": "passed"}
    path = _receipt_path(root, commit)
    _atomic_write(path, (json.dumps(payload, indent=2) + "\n").encode(), 0o600)
    return path


def _passed(root: Path, commit: str) -> bool:
    try:
        record = json.loads(_receipt_path(root, commit).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return False
    if not isinstance(record, dict) or set(record) != {"schema", "commit", "python", "result"}:
        return False
    return (type(record["schema"]) is int and record["schema"] == SCHEMA and record["commit"] == commit
            and record["result"] == "passed" and is_gate_version(record["python"]))


def pre_push(root: Path, updates: Iterable[str]) -> None:
    for line in updates:
        fields = line.split()
        if len(fields) != 4:
            raise ValueError("Malformed Git pre-push ref update.")
        _, commit, destination, _ = fields
        if destination != "refs/heads/main":
            continue
        if not _is_commit(commit):
            raise ValueError("An update to refs/heads/main must publish a tested commit.")
        if not _passed(root, commit):
            raise ValueError(f"{commit} is missing a passing local CI receipt for Python {PYTHON_VERSION}.\n"
                             f"Check out and test the commit before pushing to main.\n  {RUN_GATE}")


def _git_path(root: Path, *args: str) -> Path:
    return Path(subprocess.check_output(("git", "rev-parse", "--path-format=absolute", *args),
                                        cwd=root, text=True).strip())


def install(root: Path) -> Path:
    hook = _git_path(root, "--git-path", "hooks/pre-push")
    default = _git_path(root, "--git-common-dir") / "hooks" / "pre-push"
    if hook != default or hook.parent.is_symlink() or hook.is_symlink():
        raise ValueError(f"Git uses a custom hook location. Review it before installing: {hook}")
    if hook.exists() and (not hook.is_file() or HOOK_MARKER not in hook.read_bytes().splitlines()[:2]):
        raise ValueError(f"An existing pre-push hook is installed: {hook}")
    _atomic_write(hook, Path(__file__).read_bytes(), 0o755)
    return hook


def main() -> int:
    try:
        root = Path.cwd()
        if sys.argv[1:] == ["--install"]:
            print(f"Installed the local CI pre-push hook at {install(root)}")
        else:
            pre_push(root, sys.stdin)
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        print(f"OPH-FEP-Synergetics local CI: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
