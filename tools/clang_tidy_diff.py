#!/usr/bin/env python3
"""Run Pappus's clang-tidy naming policy on only changed source lines."""

from __future__ import annotations

import argparse
import json
import pathlib
import re
import subprocess
import sys
from collections import defaultdict

SOURCE_SUFFIXES = {".c", ".cc", ".cpp", ".cxx"}
HEADER_SUFFIXES = {".h", ".hh", ".hpp", ".hxx"}
HUNK = re.compile(r"^@@ -\d+(?:,\d+)? \+(\d+)(?:,(\d+))? @@")
FORBIDDEN_CONSTANT_PREFIX = re.compile(r"\b(?:inline\s+)?constexpr\b[^;=]*\bk_[A-Za-z0-9_]*\b")


def command(*args: str, cwd: pathlib.Path) -> str:
    return subprocess.check_output(args, cwd=cwd, text=True)


def changed_lines(diff: str, root: pathlib.Path) -> dict[pathlib.Path, list[list[int]]]:
    result: dict[pathlib.Path, list[list[int]]] = defaultdict(list)
    path: pathlib.Path | None = None
    for line in diff.splitlines():
        if line.startswith("+++ "):
            name = line[4:]
            path = None if name == "/dev/null" else root / name.removeprefix("b/")
            continue
        match = HUNK.match(line)
        if match is not None and path is not None:
            first, count = int(match.group(1)), int(match.group(2) or 1)
            if count:
                result[path].append([first, first + count - 1])
    return result


def is_pappus_source(path: pathlib.Path, root: pathlib.Path) -> bool:
    try:
        relative = path.relative_to(root)
    except ValueError:
        return False
    return relative.parts[:2] in (("include", "pappus"), ("test", "source")) and path.suffix in SOURCE_SUFFIXES | HEADER_SUFFIXES


def check_forbidden_constant_prefix(paths: dict[pathlib.Path, list[list[int]]]) -> None:
    violations: list[str] = []
    for path, ranges in paths.items():
        if path.suffix not in SOURCE_SUFFIXES | HEADER_SUFFIXES:
            continue
        lines = path.read_text().splitlines()
        for first, last in ranges:
            for number in range(first, last + 1):
                if FORBIDDEN_CONSTANT_PREFIX.search(lines[number - 1]):
                    violations.append(f"{path}:{number}: forbidden k_ prefix on a constexpr constant")
    if violations:
        print("\n".join(violations), file=sys.stderr)
        raise SystemExit(1)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", required=True, help="Git revision to diff against")
    parser.add_argument("--build-dir", required=True, type=pathlib.Path)
    parser.add_argument("--clang-tidy", default="clang-tidy")
    args = parser.parse_args()

    root = pathlib.Path(__file__).resolve().parents[1]
    diff = command("git", "diff", "--no-ext-diff", "--unified=0", args.base, cwd=root)
    lines = {path: ranges for path, ranges in changed_lines(diff, root).items() if is_pappus_source(path, root)}
    if not lines:
        print("No changed Pappus C++ source or headers to check.")
        return

    check_forbidden_constant_prefix(lines)
    database = json.loads((root / args.build_dir / "compile_commands.json").read_text())
    changed_sources = {path.resolve() for path in lines if path.suffix in SOURCE_SUFFIXES}
    changed_headers = any(path.suffix in HEADER_SUFFIXES for path in lines)
    translation_units = {
        pathlib.Path(entry["file"]).resolve()
        for entry in database
        if pathlib.Path(entry["file"]).resolve() in changed_sources or changed_headers
    }
    if not translation_units:
        raise SystemExit("Changed Pappus code has no matching compile command.")

    line_filter = json.dumps([
        {"name": str(path.relative_to(root)), "lines": ranges}
        for path, ranges in sorted(lines.items())
    ])
    header_filter = ".*"
    for source in sorted(translation_units):
        subprocess.run([
            args.clang_tidy,
            "-p", str(args.build_dir),
            "--config-file", str(root / ".clang-tidy"),
            "--line-filter", line_filter,
            "--header-filter", header_filter,
            str(source),
        ], cwd=root, check=True)


if __name__ == "__main__":
    main()
