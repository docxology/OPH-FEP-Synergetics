"""Source pins, licence terms and reuse policy for the three mapped programs.

The registry records what each source is, which commit was read, what its
licence permits, and which identifiers this package cites.  Nothing here
relicenses source material: the module reads the pins from ``data/sources.json``
and, when a checkout is reachable, reports whether the recorded commit still
matches that checkout.

A pin that has *moved* is reported as moved rather than as a failure.  A pin is
a point-in-time reading, not a claim about the source's present state.
"""

from __future__ import annotations

import json
import os
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

#: Environment variable that overrides the committed pin file.
SOURCES_ENV = "OPH_FEP_SYNERGETICS_SOURCES"


def _default_sources() -> Path:
    """Locate the committed pin file from either a source or an installed layout.

    Wheels bundle the repository's ``data/sources.json`` inside the package.
    The ancestor search also supports imports directly from a source checkout.
    """
    override = os.environ.get(SOURCES_ENV)
    if override:
        return Path(override)
    for parent in Path(__file__).resolve().parents:
        candidate = parent / "data" / "sources.json"
        if candidate.is_file():
            return candidate
    return Path(__file__).resolve().parents[2] / "data" / "sources.json"


#: The committed pin file, resolved at import time.
DEFAULT_SOURCES = _default_sources()


@dataclass(frozen=True)
class LicenceTerm:
    """One licence term applying to a declared scope of a source."""

    scope: str
    spdx: str
    evidence: str


@dataclass(frozen=True)
class CitedIdentifier:
    """One source identifier this package relies on, with the role it plays."""

    kind: str
    path: str
    symbol: str
    role: str
    lines: str | None = None


@dataclass(frozen=True)
class SourcePin:
    """A pinned source program and the terms under which it is cited."""

    id: str
    name: str
    url: str
    checkout: str
    commit: str
    branch: str
    spdx: str
    reuse_policy: str
    terms: tuple[LicenceTerm, ...] = ()
    cited: tuple[CitedIdentifier, ...] = ()
    raw: dict = field(default_factory=dict, repr=False)

    @property
    def short_commit(self) -> str:
        return self.commit[:12]

    def identifier(self, symbol: str) -> CitedIdentifier:
        """Return the cited identifier with this symbol, or raise ``KeyError``."""
        for entry in self.cited:
            if entry.symbol == symbol:
                return entry
        raise KeyError(f"{self.id} does not cite {symbol!r}")


@dataclass(frozen=True)
class PinCheck:
    """The outcome of comparing a recorded pin with a live checkout."""

    source_id: str
    checkout: str
    expected: str
    observed: str | None
    status: str

    @property
    def matched(self) -> bool:
        return self.status == "match"


def load_sources(path: Path | str = DEFAULT_SOURCES) -> tuple[SourcePin, ...]:
    """Read the committed source pins in declaration order."""
    document = json.loads(Path(path).read_text(encoding="utf-8"))
    if document.get("schema") != 1:
        raise ValueError("unsupported sources schema")
    pins: list[SourcePin] = []
    for entry in document["sources"]:
        licence = entry["license"]
        terms = tuple(
            LicenceTerm(scope=term["scope"], spdx=term["spdx"], evidence=term["evidence"])
            for term in licence.get("split", ())
        )
        cited = tuple(
            CitedIdentifier(
                kind=item["kind"],
                path=item["path"],
                symbol=item["symbol"],
                role=item["role"],
                lines=item.get("lines"),
            )
            for item in entry.get("cited_identifiers", ())
        )
        pins.append(
            SourcePin(
                id=entry["id"],
                name=entry["name"],
                url=entry["url"],
                checkout=entry["checkout"],
                commit=entry["commit"],
                branch=entry["branch"],
                spdx=licence["spdx"],
                reuse_policy=entry["reuse_policy"],
                terms=terms,
                cited=cited,
                raw=entry,
            )
        )
    return tuple(pins)


def head_of(checkout: Path) -> str | None:
    """Return the commit at ``HEAD`` of a checkout, or ``None`` when unavailable."""
    observed, _ = _read_checkout_head(checkout)
    return observed


def _read_checkout_head(checkout: Path) -> tuple[str | None, str]:
    """Distinguish missing Git metadata from a failed observation of it."""
    try:
        (checkout / ".git").lstat()
    except (FileNotFoundError, NotADirectoryError):
        return None, "absent"
    except OSError:
        return None, "unavailable"

    # Git's repository selectors take precedence over -C. Read the indicated
    # checkout even when the caller is running from a different Git operation.
    env = os.environ.copy()
    for name in ("GIT_DIR", "GIT_WORK_TREE", "GIT_COMMON_DIR", "GIT_NAMESPACE"):
        env.pop(name, None)
    env["GIT_NO_LAZY_FETCH"] = "1"
    try:
        completed = subprocess.run(
            ["git", "-C", str(checkout), "--git-dir=.git", "rev-parse", "--verify", "-q", "HEAD^{commit}"],
            capture_output=True,
            text=True,
            timeout=30,
            stdin=subprocess.DEVNULL,
            env=env,
        )
    except (OSError, subprocess.SubprocessError):
        return None, "unavailable"
    observed = completed.stdout.strip()
    if completed.returncode != 0 or not observed:
        return None, "unavailable"
    return observed, "observed"


def check_pins(
    sources: tuple[SourcePin, ...] | None = None,
    *,
    vault_root: Path | str | None = None,
) -> tuple[PinCheck, ...]:
    """Compare each recorded pin with the live checkout beneath ``vault_root``.

    With no ``vault_root`` the check reports ``unavailable`` for every pin
    instead of guessing a location.
    """
    pins = load_sources() if sources is None else sources
    root = Path(vault_root) if vault_root is not None else None
    checks: list[PinCheck] = []
    for pin in pins:
        if root is None:
            checks.append(PinCheck(pin.id, pin.checkout, pin.commit, None, "unavailable"))
            continue
        observed, status = _read_checkout_head(root / pin.checkout)
        if observed is not None:
            status = "match" if observed == pin.commit else "moved"
        checks.append(PinCheck(pin.id, pin.checkout, pin.commit, observed, status))
    return tuple(checks)
