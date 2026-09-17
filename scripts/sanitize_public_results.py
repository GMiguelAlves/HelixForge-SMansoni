#!/usr/bin/env python3
"""Replace private runtime prefixes before publishing text-based artifacts."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


TEXT_SUFFIXES = {".csv", ".html", ".json", ".md", ".txt", ".tsv", ".yaml", ".yml"}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("--replace", action="append", default=[], metavar="PRIVATE=SYMBOLIC")
    parser.add_argument("--forbid", action="append", default=[])
    parser.add_argument("--report", required=True, type=Path)
    args = parser.parse_args()

    replacements: list[tuple[str, str]] = []
    for item in args.replace:
        if "=" not in item:
            parser.error(f"invalid replacement: {item!r}")
        private, symbolic = item.split("=", 1)
        if not private:
            parser.error("private replacement prefix cannot be empty")
        replacements.append((private, symbolic))

    files_scanned = 0
    files_changed = 0
    replacements_applied = 0
    violations: list[str] = []
    for path in sorted(args.root.rglob("*")):
        if not path.is_file() or path.suffix.lower() not in TEXT_SUFFIXES:
            continue
        files_scanned += 1
        text = path.read_text(encoding="utf-8")
        updated = text
        for private, symbolic in replacements:
            variants = (
                (private, symbolic),
                (private.replace("/", r"\/"), symbolic.replace("/", r"\/")),
            )
            for source, target in variants:
                count = updated.count(source)
                if count:
                    replacements_applied += count
                    updated = updated.replace(source, target)
        if updated != text:
            path.write_text(updated, encoding="utf-8")
            files_changed += 1
        for forbidden in args.forbid:
            if forbidden and forbidden in updated:
                violations.append(str(path.relative_to(args.root)))
                break

    report = {
        "schema_version": "1.0",
        "status": "PASS" if not violations else "FAIL",
        "files_scanned": files_scanned,
        "files_changed": files_changed,
        "replacements_applied": replacements_applied,
        "violations": violations,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if not violations else 1


if __name__ == "__main__":
    raise SystemExit(main())
