#!/usr/bin/env python3
"""
Ensure an MCMCtree control file has `print = 1`.

IQ-TREE 3 `--dating mcmctree` may generate a ctl without a real `print` line
(e.g. `... FossilErrprint = 1` inside a comment due to a missing newline),
which causes MCMCtree to skip writing `FigTree.tre` and the MCMC log.

This script is idempotent:
- If a non-comment `print = ...` line exists, it is rewritten to `print = 1`.
- Otherwise, `print = 1` is appended.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

PRINT_RE = re.compile(r"^\s*print\s*=", re.IGNORECASE)


def split_comment(line: str) -> tuple[str, str]:
    """
    Split on the first '*' (MCMCtree comment marker).
    Returns (before_comment, comment_including_star_or_empty).
    """
    before, sep, after = line.partition("*")
    if not sep:
        return line, ""
    return before, f"*{after}"


def ensure_print(ctl_path: Path) -> bool:
    lines = ctl_path.read_text().splitlines(keepends=True)
    out_lines: list[str] = []
    found = False

    for raw in lines:
        newline = "\n" if raw.endswith("\n") else ""
        body = raw[:-1] if newline else raw

        before_comment, comment = split_comment(body)
        if PRINT_RE.match(before_comment):
            indent = re.match(r"^(\s*)", before_comment).group(1)  # type: ignore[union-attr]
            out_lines.append(f"{indent}print = 1{(' ' if comment and not comment.startswith(' *') else '')}{comment}{newline}")
            found = True
        else:
            out_lines.append(raw)

    if not found:
        if out_lines and not out_lines[-1].endswith("\n"):
            out_lines[-1] += "\n"
        out_lines.append("print = 1\n")

    ctl_path.write_text("".join(out_lines))
    return True


def main() -> int:
    if len(sys.argv) != 2:
        return 0
    ctl_path = Path(sys.argv[1])
    if not ctl_path.exists():
        return 0
    ensure_print(ctl_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

