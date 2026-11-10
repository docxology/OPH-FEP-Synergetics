#!/usr/bin/env python3
"""Resolve the manuscript tokens and inject them into ``output/manuscript/``.

The rendering stage calls this script before Pandoc runs
(``infrastructure.rendering._manuscript_source.run_manuscript_variable_script``
looks for exactly this path).  All computation lives in
``oph_fep_synergetics.manuscript_variables``; all substitution lives in
``infrastructure.rendering.manuscript_injection``.

The script additionally fails when a token survives hydration.  The rendering
stage only warns about unresolved markdown tokens, so without this check a
mistyped token would ship verbatim into the PDF.

Exit codes:
    0   tokens written, tree injected, no token left unresolved
    1   a token survived hydration, or the infrastructure import failed
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

_TEMPLATE_ROOT = os.environ.get("TEMPLATE_REPO_ROOT")
if _TEMPLATE_ROOT:
    sys.path.insert(0, _TEMPLATE_ROOT)

TOKEN_RE = re.compile(r"\{\{([A-Z][A-Z0-9_]*)\}\}")
CODE_SPAN_RE = re.compile(r"`[^`\n]*`")


def main() -> int:
    """Write the token map, inject the manuscript tree, and verify hydration."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--allow-draft",
        action="store_true",
        help=(
            "accepted for parity with the template exemplars; this project's tokens are "
            "bound to values rather than to analysis outputs, so a draft mode changes nothing"
        ),
    )
    parser.parse_args()

    try:
        from infrastructure.rendering.manuscript_injection import write_resolved_manuscript_tree
    except ImportError as exc:  # pragma: no cover - depends on the caller's environment
        print(
            "cannot import the rendering infrastructure; run this script from the "
            f"docxology/template checkout with TEMPLATE_REPO_ROOT set ({exc})",
            file=sys.stderr,
        )
        return 1

    from oph_fep_synergetics.manuscript_variables import as_tokens

    tokens = as_tokens()
    out_path = PROJECT_ROOT / "output" / "data" / "manuscript_variables.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(tokens, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    resolved_root = write_resolved_manuscript_tree(PROJECT_ROOT, tokens)
    print(out_path)
    print(resolved_root)

    leftovers: dict[str, list[str]] = {}
    for section in sorted(Path(resolved_root).glob("*.md")):
        # A token inside a backtick code span documents the syntax rather than
        # asking for a value, so it is stripped before the scan.  This matches
        # the infrastructure's own config.yaml rule.
        text = CODE_SPAN_RE.sub("", section.read_text(encoding="utf-8"))
        found = sorted(set(TOKEN_RE.findall(text)))
        if found:
            leftovers[section.name] = found
    if leftovers:
        for name, keys in leftovers.items():
            print(f"unresolved token in {name}: {', '.join(keys)}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
