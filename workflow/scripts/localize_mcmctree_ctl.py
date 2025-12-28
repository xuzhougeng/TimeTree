#!/usr/bin/env python3
"""
Rewrite path values in an MCMCtree control file to local basenames.
This makes the ctl runnable from a staging directory.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

KEYS = {"seqfile", "treefile", "mcmcfile", "outfile", "hessianfile", "ckpfile"}


def localize_line(line: str) -> str:
    stripped = line.lstrip()
    lower = stripped.lower()
    for key in KEYS:
        if lower.startswith(key):
            if "=" not in stripped:
                return line
            left, right = stripped.split("=", 1)
            value_part, comment = (right.split("*", 1) + [""])[:2]
            value = value_part.strip()
            if not value:
                return line
            value = value.strip("'\"")
            basename = Path(value).name
            new_right = f" {basename}"
            if comment:
                new_right += f" *{comment}"
            prefix = line[: len(line) - len(stripped)]
            return f"{prefix}{left.strip()} ={new_right}\n"
    return line


def main() -> int:
    if len(sys.argv) != 2:
        return 0
    ctl_path = Path(sys.argv[1])
    if not ctl_path.exists():
        return 0
    text = ctl_path.read_text()
    lines = [localize_line(line) for line in text.splitlines(keepends=True)]
    ctl_path.write_text("".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
