#!/usr/bin/env python3
"""
Extract mcmcfile value from an MCMCtree control file.
Prints the value without a trailing newline, or empty string if not found.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Optional


def parse_mcmcfile(text: str) -> Optional[str]:
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("*"):
            continue
        if line.lower().startswith("mcmcfile"):
            if "=" not in raw:
                return None
            value = raw.split("=", 1)[1]
            value = value.split("*", 1)[0].strip()
            return value or None
    return None


def main() -> int:
    if len(sys.argv) != 2:
        return 0
    ctl_path = Path(sys.argv[1])
    if not ctl_path.exists():
        return 0
    mcmcfile = parse_mcmcfile(ctl_path.read_text())
    if mcmcfile:
        print(mcmcfile, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
