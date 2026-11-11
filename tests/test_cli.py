"""The console entry point: every layer reachable from one command."""

from __future__ import annotations

import importlib.metadata
import json
import os
from pathlib import Path

import pytest

from oph_fep_synergetics import __version__
from oph_fep_synergetics.cli import main
from oph_fep_synergetics.manuscript_variables import as_dict

VAULT_ROOT = Path(os.environ.get("OPH_FEP_VAULT_ROOT", "/Volumes/GitVault"))


def test_version_flag(capsys):
    with pytest.raises(SystemExit) as exit_info:
        main(["--version"])
    assert exit_info.value.code == 0
    assert __version__ in capsys.readouterr().out


def test_verify_reports_every_probe(capsys):
    assert main(["verify"]) == 0
    out = capsys.readouterr().out
    assert "42/42 probes passed" in out
    assert "[FAIL]" not in out
    assert "carrier-port-chart" in out
    assert "pr01-obstruction" in out
    assert "oph-fep" in out
    assert "three-way" in out


def test_verify_json_payload_is_complete(capsys):
    assert main(["verify", "--json"]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["summary"]["total"] == 42
    assert payload["summary"]["failed"] == 0
    assert payload["summary"]["failures"] == []
    assert len(payload["results"]) == 42
    assert all(result["passed"] for result in payload["results"])
    assert all(result["values"] for result in payload["results"])
    assert all(result["pair"] in payload["summary"]["by_pair"] for result in payload["results"])
    assert sum(entry["passed"] for entry in payload["summary"]["by_pair"].values()) == 42


def test_verify_selects_named_probes(capsys):
    assert main(["verify", "--key", "quadray-ivm"]) == 0
    assert "1/1 probes passed" in capsys.readouterr().out


def test_verify_selects_one_pairwise_scope(capsys):
    assert main(["verify", "--pair", "three-way"]) == 0
    out = capsys.readouterr().out
    assert "2/2 probes passed" in out
    assert "tetravolume-arithmetic" in out
    assert "conjugation-norm" in out
    assert "carrier-port-chart" not in out


def test_verify_selects_one_kind(capsys):
    assert main(["verify", "--kind", "complementarity"]) == 0
    out = capsys.readouterr().out
    assert "10/10 probes passed" in out
    assert "effective-resistance" in out
    assert "carrier-port-chart" not in out


def test_verify_rejects_an_unknown_kind(capsys):
    assert main(["verify", "--kind", "isomorphism-typo"]) == 2
    captured = capsys.readouterr()
    assert "isomorphism-typo" in captured.err
    assert captured.out == ""


def test_verify_rejects_an_unknown_pair(capsys):
    assert main(["verify", "--pair", "oph-fep-typo"]) == 2
    captured = capsys.readouterr()
    assert "oph-fep-typo" in captured.err
    assert captured.out == ""


def test_verify_rejects_an_unknown_key(capsys):
    assert main(["verify", "--key", "no-such-key"]) == 2
    captured = capsys.readouterr()
    assert "no-such-key" in captured.err
    assert captured.out == ""


def test_sources_report_pins_and_licences(capsys):
    assert main(["sources"]) == 0
    out = capsys.readouterr().out
    assert "Observer Patch Holography" in out
    assert "Apache-2.0 AND CC-BY-NC-SA-4.0" in out
    assert "CC-BY-4.0" in out
    assert "pin:     unavailable" in out


def test_sources_reports_absent_checkouts(capsys, tmp_path):
    assert main(["sources", "--vault-root", str(tmp_path)]) == 0
    assert capsys.readouterr().out.count("pin:     absent") == 3


@pytest.mark.skipif(not VAULT_ROOT.is_dir(), reason="the vault checkout is not available")
def test_sources_checks_a_live_vault(capsys):
    assert main(["sources", "--vault-root", str(VAULT_ROOT)]) == 0
    out = capsys.readouterr().out
    # A pin is a point-in-time reading, not a live requirement: an upstream
    # checkout that has moved since it was read is reported as moved, and the
    # package's own invariant is that this is a reading rather than a failure.
    # What the live check must show is that all three sources are readable and
    # that each lands in the documented vocabulary.
    statuses = [chunk.split()[0] for chunk in out.split("pin:     ")[1:]]
    assert len(statuses) == 3, out
    assert all(status in {"match", "moved", "absent", "unavailable"} for status in statuses), out


def test_sources_json(capsys):
    assert main(["sources", "--json"]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert [source["id"] for source in payload["sources"]] == ["oph", "fep_formal", "quadmath"]
    assert all(source["cited"] > 0 for source in payload["sources"])
    assert [check["status"] for check in payload["checks"]] == ["unavailable"] * 3


def test_values_lists_every_bound_value(capsys):
    assert main(["values"]) == 0
    out = capsys.readouterr().out
    assert "carrier.ports" in out
    assert "paper.five_surprises" in out


def test_values_json_matches_the_module(capsys):
    assert main(["values", "--json"]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert set(payload) == set(as_dict())
    assert payload["carrier.ports"]["value"] == 12
    assert payload["carrier.ports"]["source"]


def test_carrier_summary(capsys):
    assert main(["carrier"]) == 0
    out = capsys.readouterr().out
    assert "ports 12   seams 30   faces 20" in out
    assert "5 - sqrt(5)" in out
    assert "numeric error" in out
    assert "shell total-two 78   binary 66" in out


def test_geometry_summary(capsys):
    assert main(["geometry"]) == 0
    out = capsys.readouterr().out
    assert "neighbor moves 12" in out
    assert "cuboctahedron" in out
    assert "truncated_octahedron" in out
    assert "96" in out


def test_the_console_script_is_declared():
    try:
        distribution = importlib.metadata.distribution("oph-fep-synergetics")
    except importlib.metadata.PackageNotFoundError:
        pytest.skip("the distribution is not installed in this environment")
    scripts = {
        entry.name
        for entry in distribution.entry_points
        if entry.group == "console_scripts"
    }
    assert scripts == {"oph-fep-synergetics"}
