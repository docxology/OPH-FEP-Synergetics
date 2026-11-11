"""Command-line entry point for the correspondence battery."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Sequence

from . import __version__
from .analysis import check_pins, load_sources, summary, verify
from .core import icosahedral_carrier, shell_gap
from .geometry import (
    MAX_SHELL,
    VOLUME_RATIOS,
    cumulative_count,
    neighbor_moves,
    shell_count,
    squared_norm,
)
from .manuscript_variables import MANUSCRIPT_VALUES


def _json_default(value: object) -> object:
    if isinstance(value, (set, frozenset)):
        return sorted(value)
    if hasattr(value, "tolist"):
        return value.tolist()  # type: ignore[attr-defined]
    return str(value)


def _cmd_verify(args: argparse.Namespace) -> int:
    keys = tuple(args.key) if args.key else None
    pairs = tuple(args.pair) if args.pair else None
    kinds = tuple(args.kind) if args.kind else None
    try:
        results = verify(keys=keys, pairs=pairs, kinds=kinds)
    except ValueError as error:
        print(f"error: {error}", file=sys.stderr)
        return 2
    report = summary(results)
    if args.json:
        payload = {
            "summary": report,
            "results": [
                {
                    "key": result.key,
                    "kind": result.kind,
                    "pair": result.pair,
                    "status": result.status,
                    "passed": result.passed,
                    "values": result.values,
                }
                for result in results
            ],
        }
        print(json.dumps(payload, indent=2, default=_json_default))
    else:
        for result in results:
            mark = "PASS" if result.passed else "FAIL"
            print(f"[{mark}] {result.pair:<16} {result.kind:<16} {result.status:<16} {result.key}")
        print(f"\n{report['passed']}/{report['total']} probes passed")
        for pair, bucket in report["by_pair"].items():
            print(f"  {pair:<18} {bucket['passed']} passed, {bucket['failed']} failed")
        for kind, bucket in report["by_kind"].items():
            print(f"  {kind:<18} {bucket['passed']} passed, {bucket['failed']} failed")
        if report["failures"]:
            print("failed: " + ", ".join(report["failures"]))
    return 0 if report["failed"] == 0 else 1


def _cmd_sources(args: argparse.Namespace) -> int:
    pins = load_sources()
    root = Path(args.vault_root) if args.vault_root else None
    checks = check_pins(pins, vault_root=root)
    if args.json:
        print(
            json.dumps(
                {
                    "sources": [
                        {
                            "id": pin.id,
                            "name": pin.name,
                            "url": pin.url,
                            "commit": pin.commit,
                            "spdx": pin.spdx,
                            "reuse_policy": pin.reuse_policy,
                            "terms": [
                                {"scope": term.scope, "spdx": term.spdx, "evidence": term.evidence}
                                for term in pin.terms
                            ],
                            "cited": len(pin.cited),
                        }
                        for pin in pins
                    ],
                    "checks": [
                        {
                            "source_id": check.source_id,
                            "status": check.status,
                            "expected": check.expected,
                            "observed": check.observed,
                        }
                        for check in checks
                    ],
                },
                indent=2,
                default=_json_default,
            )
        )
        return 0
    for pin, check in zip(pins, checks):
        print(f"{pin.id} ({pin.name})")
        print(f"  url:     {pin.url}")
        print(f"  commit:  {pin.commit}")
        print(f"  licence: {pin.spdx}")
        print(f"  reuse:   {pin.reuse_policy}")
        for term in pin.terms:
            print(f"    - {term.spdx:<18} {term.scope}  [{term.evidence}]")
        print(f"  cited:   {len(pin.cited)} identifiers")
        if check.observed and not check.matched:
            print(f"  pin:     {check.status} (live HEAD {check.observed[:12]})")
        else:
            print(f"  pin:     {check.status}")
    return 0


def _cmd_values(args: argparse.Namespace) -> int:
    if args.json:
        print(
            json.dumps(
                {
                    entry.key: {
                        "value": entry.value,
                        "exact": entry.exact,
                        "source": entry.source,
                        "note": entry.note,
                    }
                    for entry in MANUSCRIPT_VALUES
                },
                indent=2,
                default=_json_default,
            )
        )
        return 0
    for entry in MANUSCRIPT_VALUES:
        exact = f"   ({entry.exact})" if entry.exact else ""
        print(f"{entry.key:<38} {entry.value!r}{exact}")
        print(f"{'':<38} source: {entry.source}")
        if entry.note:
            print(f"{'':<38} note:   {entry.note}")
    return 0


def _cmd_carrier(args: argparse.Namespace) -> int:
    carrier = icosahedral_carrier()
    numeric = carrier.spectral_gap()
    exact = 5.0 - math.sqrt(5.0)
    gap = shell_gap()
    exact_gap = exact / 60.0
    print(f"ports {carrier.port_count}   seams {carrier.seam_count}   faces {carrier.face_count}")
    print(f"degrees {sorted({carrier.degree(port) for port in range(carrier.port_count)})}")
    print(
        f"lambda_1 {numeric:.15f}   exact 5 - sqrt(5) = {exact:.15f}   "
        f"numeric error {abs(numeric - exact):.1e}"
    )
    print(
        f"shell gap {gap:.15f}   exact (5 - sqrt(5)) / 60 = {exact_gap:.15f}   "
        f"numeric error {abs(gap - exact_gap):.1e}"
    )
    print(
        "shell total-two "
        f"{len(carrier.shell_configurations(total=2))}   "
        f"binary {len(carrier.binary_shell_configurations())}"
    )
    return 0


def _cmd_geometry(args: argparse.Namespace) -> int:
    moves = neighbor_moves()
    print(f"neighbor moves {len(moves)}   squared norm {float(squared_norm(moves[0])):.1f}")
    print(f"shell counts   {[int(shell_count(k)) for k in range(6)]}")
    print(f"cumulative     {[int(cumulative_count(k)) for k in range(6)]}")
    print(f"max shell      {MAX_SHELL}")
    print("volume ratios (unit IVM tetrahedron = 1)")
    for name, ratio in sorted(VOLUME_RATIOS.items(), key=lambda item: str(item[0])):
        print(f"  {str(name):<24} {ratio}")
    return 0


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="oph-fep-synergetics",
        description=(
            "Executable correspondences, obstructions and complements across "
            "Observer Patch Holography, the Free Energy Principle and Fuller's Synergetics."
        ),
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    verify_parser = subparsers.add_parser(
        "verify", help="run the correspondence registry and report every probe"
    )
    verify_parser.add_argument(
        "--key", action="append", default=[], help="run only this correspondence; repeatable"
    )
    verify_parser.add_argument(
        "--pair",
        action="append",
        default=[],
        help="run only this pairwise scope (oph-fep, fep-synergetics, synergetics-oph, three-way); repeatable",
    )
    verify_parser.add_argument(
        "--kind",
        action="append",
        default=[],
        help="run only this kind of claim (isomorphism, disagreement, complementarity); repeatable",
    )
    verify_parser.add_argument("--json", action="store_true", help="emit machine-readable output")
    verify_parser.set_defaults(handler=_cmd_verify)

    sources_parser = subparsers.add_parser(
        "sources", help="show source pins, licence terms and reuse policy"
    )
    sources_parser.add_argument(
        "--vault-root", default=None, help="check the pins against checkouts under this root"
    )
    sources_parser.add_argument("--json", action="store_true", help="emit machine-readable output")
    sources_parser.set_defaults(handler=_cmd_sources)

    values_parser = subparsers.add_parser(
        "values", help="show the source-bound manuscript values"
    )
    values_parser.add_argument("--json", action="store_true", help="emit machine-readable output")
    values_parser.set_defaults(handler=_cmd_values)

    carrier_parser = subparsers.add_parser("carrier", help="show the carrier invariants")
    carrier_parser.set_defaults(handler=_cmd_carrier)

    geometry_parser = subparsers.add_parser("geometry", help="show the quadray/IVM invariants")
    geometry_parser.set_defaults(handler=_cmd_geometry)

    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Run the command line and return its exit status."""
    parser = _build_parser()
    args = parser.parse_args(argv)
    return int(args.handler(args))


if __name__ == "__main__":  # pragma: no cover - module entry point
    sys.exit(main())
