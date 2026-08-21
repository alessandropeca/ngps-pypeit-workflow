#!/usr/bin/env python3
"""Apply the pinned, opt-in manual-trace-refit patch to a PypeIt checkout."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


PINNED_COMMIT = "e9ed85c1a237c49626227f4227e323fc390def4b"
PATCH = Path(__file__).resolve().parents[1] / "patches" / "pypeit-manual-trace-refit.patch"


def git(root: Path, *arguments: str, check: bool = True) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ["git", "-C", str(root), *arguments], text=True, capture_output=True, check=check,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pypeit_root", type=Path, help="Path to the pinned PypeIt checkout")
    args = parser.parse_args()
    root = args.pypeit_root.expanduser().resolve()
    if not (root / ".git").exists():
        parser.error(f"Not a Git checkout: {root}")
    if not PATCH.is_file():
        parser.error(f"Missing workflow patch: {PATCH}")
    if git(root, "apply", "--reverse", "--check", str(PATCH), check=False).returncode == 0:
        print("PypeIt manual-trace-refit patch is already applied.")
        return 0
    commit = git(root, "rev-parse", "HEAD").stdout.strip()
    if commit != PINNED_COMMIT:
        parser.error(
            "PypeIt is not at the workflow's pinned commit. "
            f"Expected {PINNED_COMMIT}, found {commit}."
        )
    check = git(root, "apply", "--check", str(PATCH), check=False)
    if check.returncode != 0:
        sys.stderr.write(check.stderr)
        return check.returncode
    git(root, "apply", str(PATCH))
    print("Applied the pinned manual-trace-refit patch to PypeIt.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
