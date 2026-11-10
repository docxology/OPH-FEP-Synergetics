"""Source pins: what each cited program is, and whether the checkout moved."""

from __future__ import annotations

import json
import subprocess

import pytest

from oph_fep_synergetics.analysis import (
    DEFAULT_SOURCES,
    SOURCES_ENV,
    PinCheck,
    SourcePin,
    check_pins,
    head_of,
    load_sources,
)


def test_the_committed_pin_file_is_found():
    assert DEFAULT_SOURCES.is_file()
    assert DEFAULT_SOURCES.name == "sources.json"
    assert SOURCES_ENV == "OPH_FEP_SYNERGETICS_SOURCES"


def test_pins_record_the_three_cited_programs():
    pins = load_sources()
    assert [pin.id for pin in pins] == ["oph", "fep_formal", "quadmath"]
    by_id = {pin.id: pin for pin in pins}
    assert by_id["oph"].spdx == "Apache-2.0 AND CC-BY-NC-SA-4.0"
    assert by_id["fep_formal"].spdx == "CC-BY-4.0"
    assert by_id["quadmath"].spdx == "Apache-2.0"
    for pin in pins:
        assert len(pin.commit) == 40
        assert pin.branch
        assert pin.url.startswith("https://github.com/")
        assert pin.reuse_policy
        assert pin.short_commit == pin.commit[:12]


def test_ophs_licence_is_recorded_as_a_split():
    oph = {pin.id: pin for pin in load_sources()}["oph"]
    scopes = {term.scope: term for term in oph.terms}
    assert len(scopes) == 2
    assert any(term.spdx == "Apache-2.0" for term in oph.terms)
    assert any(term.spdx == "CC-BY-NC-SA-4.0" for term in oph.terms)
    assert all(term.evidence for term in oph.terms)
    assert "code/LICENSE" in {term.evidence for term in oph.terms}
    assert "Do not copy" in oph.reuse_policy or "do not reuse" in oph.reuse_policy


def test_every_cited_identifier_names_its_role():
    for pin in load_sources():
        assert pin.cited
        for cited in pin.cited:
            assert cited.kind
            assert cited.path
            assert cited.symbol
            assert cited.role


def test_identifier_lookup_returns_the_named_entry():
    oph = {pin.id: pin for pin in load_sources()}["oph"]
    entry = oph.identifier("native")
    assert entry.path.endswith("native.py")
    assert entry.lines == "45-51"
    with pytest.raises(KeyError):
        oph.identifier("not-a-symbol")


def test_opts_out_of_a_missing_pin_file():
    with pytest.raises(FileNotFoundError):
        load_sources(DEFAULT_SOURCES.parent / "missing.json")


def test_an_unsupported_schema_is_rejected(tmp_path):
    path = tmp_path / "sources.json"
    path.write_text(json.dumps({"schema": 99, "sources": []}), encoding="utf-8")
    with pytest.raises(ValueError):
        load_sources(path)


def test_without_a_vault_root_every_pin_is_unavailable():
    checks = check_pins(load_sources())
    assert len(checks) == 3
    assert all(check.status == "unavailable" for check in checks)
    assert all(check.observed is None for check in checks)
    assert not any(check.matched for check in checks)


def test_a_missing_checkout_is_absent(tmp_path):
    checks = check_pins(load_sources(), vault_root=tmp_path)
    assert all(check.status == "absent" for check in checks)


def _point_head_at_a_commit(checkout) -> str:
    """Create a repository with one resolvable commit, without an identity.

    The commit bytes are written straight through ``hash-object``, so nothing
    here reads or writes Git's configured author or committer identity.
    """
    subprocess.run(
        ["git", "init", "-q", "-b", "main", str(checkout)],
        check=True,
        stdin=subprocess.DEVNULL,
    )
    subprocess.run(
        ["git", "-C", str(checkout), "symbolic-ref", "HEAD", "refs/heads/main"],
        check=True,
        stdin=subprocess.DEVNULL,
    )
    tree = subprocess.run(
        ["git", "-C", str(checkout), "hash-object", "-t", "tree", "-w", "--stdin"],
        input="",
        capture_output=True,
        text=True,
        check=True,
    ).stdout.strip()
    payload = (
        f"tree {tree}\n"
        "author Fixture Test <fixture@example.invalid> 1000000000 +0000\n"
        "committer Fixture Test <fixture@example.invalid> 1000000000 +0000\n"
        "\n"
        "test fixture\n"
    )
    commit = subprocess.run(
        ["git", "-C", str(checkout), "hash-object", "-t", "commit", "-w", "--stdin"],
        input=payload,
        capture_output=True,
        text=True,
        check=True,
    ).stdout.strip()
    subprocess.run(
        ["git", "-C", str(checkout), "update-ref", "refs/heads/main", commit],
        check=True,
        stdin=subprocess.DEVNULL,
    )
    return commit


def test_a_matching_and_a_moved_checkout_are_distinguished(tmp_path):
    checkout = tmp_path / "quadmath"
    checkout.mkdir(parents=True)
    observed = _point_head_at_a_commit(checkout)
    assert head_of(checkout) == observed

    pin = SourcePin(
        id="quadmath",
        name="QuadMath",
        url="https://github.com/docxology/quadmath",
        checkout="quadmath",
        commit=observed,
        branch="main",
        spdx="Apache-2.0",
        reuse_policy="reuse with notice",
    )
    matched = check_pins((pin,), vault_root=tmp_path)[0]
    assert matched.status == "match"
    assert matched.matched is True

    moved_pin = SourcePin(
        id="quadmath",
        name="QuadMath",
        url="https://github.com/docxology/quadmath",
        checkout="quadmath",
        commit="0" * 40,
        branch="main",
        spdx="Apache-2.0",
        reuse_policy="reuse with notice",
    )
    moved = check_pins((moved_pin,), vault_root=tmp_path)[0]
    assert moved.status == "moved"
    assert moved.observed == observed
    assert moved.matched is False


def test_head_of_a_directory_that_is_not_a_repository_is_none(tmp_path):
    assert head_of(tmp_path) is None
    assert head_of(tmp_path / "does-not-exist") is None


def test_pin_check_matched_is_a_property_of_status():
    assert PinCheck("x", "path", "a", "a", "match").matched is True
    assert PinCheck("x", "path", "a", "b", "moved").matched is False
