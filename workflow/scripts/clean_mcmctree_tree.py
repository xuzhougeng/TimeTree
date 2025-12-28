#!/usr/bin/env python3
"""
Clean a tree file for MCMCtree: remove leading headers (e.g. "5 1") and keep a single Newick line.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path


HEADER_RE = re.compile(r"^\s*(\d+)\s+(\d+)")


def extract_header(lines: list[str]) -> tuple[int | None, int | None]:
    for ln in lines:
        if not ln.strip():
            continue
        m = HEADER_RE.match(ln)
        if m:
            return int(m.group(1)), int(m.group(2))
        break
    return None, None


def extract_newick(text: str) -> str:
    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    newick_parts = []
    started = False
    for ln in lines:
        if not started:
            if "(" not in ln and "[" not in ln:
                continue
            started = True
        newick_parts.append(ln)
    if not newick_parts:
        return text.strip()
    newick = " ".join(newick_parts)
    if not newick.endswith(";"):
        newick += ";"
    return newick


def fix_calibration_syntax(newick: str) -> str:
    """
    Convert IQ-TREE style 'B_41.0_52.0_' to MCMCtree 'B(41.0,52.0)'.
    Handles B/L/U with one or two bounds.
    """
    def repl(match: re.Match[str]) -> str:
        code = match.group(1)
        v1 = match.group(2)
        v2 = match.group(3)
        if v2:
            return f"'{code}({v1},{v2})'"
        else:
            return f"'{code}({v1})'"

    return re.sub(r"'([A-Za-z])_([0-9.eE+-]+)(?:_([0-9.eE+-]+))?_'", repl, newick)


def count_tips(newick: str) -> int:
    """Rough tip count: strip quoted calibrations, then count name tokens."""
    stripped = re.sub(r"'[^']*'", "", newick)
    tokens = re.findall(r"[A-Za-z0-9][A-Za-z0-9_.-]*", stripped)
    return len(tokens)


def main() -> int:
    if len(sys.argv) != 2:
        return 0
    tree_path = Path(sys.argv[1])
    if not tree_path.exists():
        return 0
    text = tree_path.read_text()
    lines = text.splitlines()
    ns, ntree = extract_header(lines)
    newick = extract_newick(text)
    newick = fix_calibration_syntax(newick)
    if ns is None:
        ns = count_tips(newick)
    if ntree is None:
        ntree = 1
    header = f"{ns} {ntree}"
    tree_path.write_text(f"{header}\n{newick}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
