"""Index the labelled objects and Pandoc cross-references in this manuscript."""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

_LABEL = r"(?:sec|fig|tbl|eq|def|prop|thm|lem|cor|rmk|exa):[a-z0-9_-]+"
_ATTRIBUTES = re.compile(r"\{([^{}\n]+)\}")
_ANCHOR = re.compile(rf"(?<!\S)#({_LABEL})(?=\s|$)")
_CITATION_GROUP = re.compile(r"\[(-?@[^\[\]]+)\]")
_REFERENCE = re.compile(rf"(?:^|;)\s*-?@({_LABEL})(?=[\s,;]|$)")
_COMMENT = re.compile(r"<!--.*?-->", re.DOTALL)
_FENCE = re.compile(r"^ {0,3}(?P<marker>`{3,}|~{3,})[^\n]*\n.*?^ {0,3}(?P=marker)[ \t]*$",
                    re.MULTILINE | re.DOTALL)
_INLINE_CODE = re.compile(r"(?P<ticks>`+)[^\n]*?(?P=ticks)")
_FIGURE = re.compile(
    r"^!\[(?P<caption>.*?)\]\(\.\./\.\./output/figures/(?P<filename>fig_[a-z0-9_]+\.png)\)"
    r"\{(?P<attrs>[^{}\n]+)\}",
    re.MULTILINE | re.DOTALL,
)


@dataclass(frozen=True)
class Figure:
    """An image declaration with its exact authored label and caption."""

    label: str
    filename: str
    caption: str


@dataclass(frozen=True)
class ManuscriptReferences:
    """Declarations and in-text citations, indexed by the renderer's labels."""

    declared: dict[str, tuple[str, ...]]
    cited_by: dict[str, tuple[str, ...]]
    figures: tuple[Figure, ...]

    def labels(self, kind: str) -> tuple[str, ...]:
        return tuple(label for label in self.declared if label.startswith(f"{kind}:"))

    def problems(self) -> list[str]:
        """Require unique targets and a use in text for each numbered object."""
        repeated = [f"duplicate label {label} in {', '.join(sections)}"
                    for label, sections in self.declared.items() if len(sections) != 1]
        unknown = [f"no target for {label} in {', '.join(sections)}"
                   for label, sections in self.cited_by.items() if label not in self.declared]
        unused = [f"no in-text reference to {label} in {', '.join(sections)}"
                  for label, sections in self.declared.items()
                  if not label.startswith("sec:") and label not in self.cited_by]
        return repeated + unknown + unused


def read_sections(manuscript: Path) -> dict[str, str]:
    """Read numbered sections in render order, including the ontology supplement."""
    return {path.name: path.read_text(encoding="utf-8") for path in sorted(manuscript.glob("[0-9]*.md"))}


def _citations(text: str) -> list[str]:
    groups = _CITATION_GROUP.findall(_INLINE_CODE.sub("", text))
    return [label for group in groups for label in _REFERENCE.findall(group)]


def index_references(sections: Mapping[str, str]) -> ManuscriptReferences:
    """Index actual object attributes and references, excluding code and comments."""
    declared: dict[str, list[str]] = {}
    cited_by: dict[str, list[str]] = {}
    figures: list[Figure] = []
    for name, source in sections.items():
        visible = _FENCE.sub("", _COMMENT.sub("", source))
        references = _citations(_FIGURE.sub("", visible))
        for match in _FIGURE.finditer(visible):
            anchors = _ANCHOR.findall(match["attrs"])
            for label in anchors:
                if label.startswith("fig:"):
                    figures.append(Figure(label, match["filename"], match["caption"]))
            # The figure declaration cannot serve as its own reference.
            references.extend(label for label in _citations(match["caption"]) if label not in anchors)
        labels = _INLINE_CODE.sub("", visible)
        for attrs in _ATTRIBUTES.findall(labels):
            for label in _ANCHOR.findall(attrs):
                declared.setdefault(label, []).append(name)

        for label in dict.fromkeys(references):
            cited_by.setdefault(label, []).append(name)
    return ManuscriptReferences(
        {label: tuple(sections) for label, sections in declared.items()},
        {label: tuple(sections) for label, sections in cited_by.items()},
        tuple(figures),
    )
