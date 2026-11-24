import json
import os
import subprocess

import pytest

from scripts import local_ci_receipts as receipts

MAIN = "refs/heads/main"
PREVIOUS = "a" * 40
CANDIDATE = "b" * 40


def update(commit=CANDIDATE, destination=MAIN):
    return f"refs/heads/topic {commit} {destination} {PREVIOUS}\n"


def test_main_requires_the_python_314_receipt_for_the_pushed_object(tmp_path):
    previous = receipts.record_pass(tmp_path, PREVIOUS, "3.14.10")
    with pytest.raises(ValueError, match=rf"{CANDIDATE}.*Python 3\.14"):
        receipts.pre_push(tmp_path, [update()])

    legacy = previous.parent.parent / CANDIDATE / "3.12.json"
    legacy.parent.mkdir()
    legacy.write_text(json.dumps({"schema": 1, "commit": CANDIDATE, "python": "3.12.16", "result": "passed"}))
    with pytest.raises(ValueError, match=r"Python 3\.14"):
        receipts.pre_push(tmp_path, [update()])
    with pytest.raises(ValueError, match="requested Python runtime"):
        receipts.record_pass(tmp_path, CANDIDATE, "3.15.2")
    other = receipts.record_pass(tmp_path, CANDIDATE, "3.14.10")
    receipts.pre_push(tmp_path, [update()])
    record = json.loads(other.read_text())
    for wrong in ({"commit": PREVIOUS}, {"python": "3.12.16"}, {"result": "failed"}, {"schema": True}):
        other.write_text(json.dumps(record | wrong))
        with pytest.raises(ValueError, match="Python 3.14"):
            receipts.pre_push(tmp_path, [update()])
    other.write_text("{broken")
    with pytest.raises(ValueError, match="Python 3.14"):
        receipts.pre_push(tmp_path, [update()])


def test_hook_filters_remote_destinations_and_rejects_invalid_main_updates(tmp_path):
    receipts.pre_push(tmp_path, [update(destination="refs/heads/review")])
    with pytest.raises(ValueError, match="must publish a tested commit"):
        receipts.pre_push(tmp_path, [update(commit="0" * 40)])
    with pytest.raises(ValueError, match="Malformed Git"):
        receipts.pre_push(tmp_path, ["truncated ref update\n"])


def test_install_keeps_an_unrelated_hook_and_installs_an_executable_copy(tmp_path, monkeypatch):
    monkeypatch.setenv("GIT_CONFIG_GLOBAL", os.devnull)
    monkeypatch.setenv("GIT_CONFIG_NOSYSTEM", "1")
    subprocess.run(["git", "init", "-q", "--template=", str(tmp_path)], check=True)
    hook = tmp_path / ".git" / "hooks" / "pre-push"
    hook.parent.mkdir()
    hook.write_text("# an unrelated hook\n")
    with pytest.raises(ValueError, match="existing pre-push hook"):
        receipts.install(tmp_path)
    assert hook.read_text() == "# an unrelated hook\n"
    hook.unlink()
    assert receipts.install(tmp_path) == hook.resolve()
    assert os.access(hook, os.X_OK)
    assert receipts.install(tmp_path) == hook.resolve()
    subprocess.run(["git", "config", "--local", "core.hooksPath", "custom-hooks"], cwd=tmp_path, check=True)
    with pytest.raises(ValueError, match="custom hook location"):
        receipts.install(tmp_path)
    assert not (tmp_path / "custom-hooks").exists()
