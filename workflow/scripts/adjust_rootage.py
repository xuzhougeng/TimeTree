#!/usr/bin/env python3
"""
Adjust RootAge in MCMCtree ctl to be safely above calibration maxima.
Usage: adjust_rootage.py <ctl_file> <tree_file>
"""

from __future__ import annotations

import re
import sys
from pathlib import Path


CAL_RE = re.compile(r"[BLU]\(\s*([0-9.eE+-]+)(?:\s*,\s*([0-9.eE+-]+))?")


def max_calibration(tree_text: str) -> float | None:
    vals = []
    for m in CAL_RE.finditer(tree_text):
        v1 = float(m.group(1))
        vals.append(v1)
        if m.group(2):
            vals.append(float(m.group(2)))
    return max(vals) if vals else None


def adjust_rootage(ctl_path: Path, tree_path: Path):
    tree_txt = tree_path.read_text()
    max_age = max_calibration(tree_txt)
    if max_age is None:
        return
    target = max_age * 1.5

    lines = ctl_path.read_text().splitlines()
    new_lines = []
    replaced = False
    for ln in lines:
        if ln.strip().lower().startswith("rootage"):
            comment = ""
            if "*" in ln:
                ln, comment = ln.split("*", 1)
                comment = "*" + comment
            new_lines.append(f"RootAge = <{target:.6g} {comment}".rstrip() + "\n")
            replaced = True
        else:
            new_lines.append(ln + "\n")
    if not replaced:
        new_lines.append(f"RootAge = <{target:.6g}\n")
    ctl_path.write_text("".join(new_lines))


def main() -> int:
    if len(sys.argv) != 3:
        return 0
    ctl = Path(sys.argv[1])
    tree = Path(sys.argv[2])
    if ctl.exists() and tree.exists():
        adjust_rootage(ctl, tree)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
