"""The generated documents stay in step with the registry they describe."""

from __future__ import annotations

import importlib.util
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _load_renderer():
    spec = importlib.util.spec_from_file_location(
        "render_docs", ROOT / "scripts" / "render_docs.py"
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_rendered_registry_document_is_current():
    renderer = _load_renderer()
    current = (ROOT / "docs" / "CORRESPONDENCES.md").read_text(encoding="utf-8")
    assert current == renderer.render_registry()


def test_rendered_values_document_is_current():
    renderer = _load_renderer()
    current = (ROOT / "docs" / "MANUSCRIPT_VALUES.md").read_text(encoding="utf-8")
    assert current == renderer.render_values()


def test_rendered_concordance_document_is_current():
    renderer = _load_renderer()
    current = (ROOT / "docs" / "CONCORDANCE.md").read_text(encoding="utf-8")
    assert current == renderer.render_concordance()


def test_the_outline_names_every_registered_claim():
    renderer = _load_renderer()
    outline = (ROOT / "manuscript" / "MANUSCRIPT_STATUS.md").read_text(encoding="utf-8")
    missing = [entry.key for entry in renderer.CORRESPONDENCES if f"`{entry.key}`" not in outline]
    assert missing == []


def test_every_claim_is_named_by_some_section():
    renderer = _load_renderer()
    text = renderer.render_concordance()
    assert "no section names it" not in text
    assert "**no**" not in text


def test_check_mode_reports_stale_documents(tmp_path, monkeypatch, capsys):
    renderer = _load_renderer()
    monkeypatch.setattr(renderer, "REGISTRY_DOC", tmp_path / "missing-registry.md")
    monkeypatch.setattr(renderer, "VALUES_DOC", tmp_path / "missing-values.md")
    monkeypatch.setattr(renderer, "CONCORDANCE_DOC", tmp_path / "missing-concordance.md")
    assert renderer.main(["--check"]) == 1
    captured = capsys.readouterr()
    assert "out of date" in captured.err
    assert not (tmp_path / "missing-registry.md").exists()


def test_write_mode_creates_missing_documents(tmp_path, monkeypatch, capsys):
    renderer = _load_renderer()
    target = tmp_path / "nested" / "registry.md"
    monkeypatch.setattr(renderer, "REGISTRY_DOC", target)
    monkeypatch.setattr(renderer, "VALUES_DOC", tmp_path / "values.md")
    monkeypatch.setattr(renderer, "CONCORDANCE_DOC", tmp_path / "concordance.md")
    assert renderer.main([]) == 0
    assert target.is_file()
    assert target.read_text(encoding="utf-8") == renderer.render_registry()
    assert "wrote" in capsys.readouterr().out


def test_format_value_renders_containers_and_scalars():
    renderer = _load_renderer()
    assert renderer.format_value({"a": 1, "b": (2, 3)}) == "a=1, b=(2, 3)"
    assert renderer.format_value("plain") == "plain"


def test_every_relative_link_in_the_guides_resolves():
    broken = []
    for guide in (ROOT / "README.md", ROOT / "AGENTS.md"):
        text = guide.read_text(encoding="utf-8")
        links = re.findall(r"\]\(([^)]+)\)", text)
        assert links, f"{guide.name} should link somewhere"
        for link in links:
            if link.startswith(("http://", "https://", "#", "mailto:")):
                continue
            target = link.split("#", 1)[0]
            if not (ROOT / target).exists():
                broken.append(f"{guide.name}: {link}")
    assert broken == []


def test_the_guides_name_the_project_commands():
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    for command in (
        "oph-fep-synergetics verify",
        "oph-fep-synergetics sources",
        "oph-fep-synergetics values",
        "scripts/render_docs.py",
    ):
        assert command in readme


def test_the_guides_state_the_collected_test_count(request):
    """A stale test count in the guides is a documentation failure.

    The count is read from the live session, so adding a test without updating
    the guides fails here rather than leaving a number that quietly stops being
    true.
    """
    count = request.session.testscollected
    for guide in (ROOT / "README.md", ROOT / "docs" / "VERIFICATION.md"):
        assert f"{count} tests" in guide.read_text(encoding="utf-8"), guide.name
