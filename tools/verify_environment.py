#!/usr/bin/env python3
"""Verify that installed upstream checkouts match upstream-lock.yml."""

from __future__ import annotations

import importlib
import subprocess
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
LOCK = yaml.safe_load((ROOT / "upstream-lock.yml").read_text())


def git_output(directory: Path, *args: str) -> str | None:
    result = subprocess.run(
        ["git", "-C", str(directory), *args],
        text=True, capture_output=True, check=False,
    )
    return result.stdout.strip() if result.returncode == 0 else None


def check_module(name: str, locked_commit: str) -> bool:
    module = importlib.import_module(name)
    module_path = Path(module.__file__).resolve()
    checkout = next((p for p in (module_path, *module_path.parents)
                     if (p / ".git").exists()), None)
    print(f"{name}: {module_path}")
    if checkout is None:
        print("  WARNING: installed from a wheel; Git revision cannot be verified.")
        return False
    actual = git_output(checkout, "rev-parse", "HEAD")
    dirty = git_output(checkout, "status", "--porcelain")
    if actual != locked_commit:
        print(f"  ERROR: expected {locked_commit}, found {actual}")
        return False
    if dirty:
        print(f"  ERROR: checkout is dirty: {checkout}")
        return False
    print(f"  OK: {actual}")
    return True


def main() -> int:
    ok = check_module("pypeit", LOCK["upstreams"]["pypeit"]["commit"])
    ok &= check_module("ngps_pipeline", LOCK["upstreams"]["ngps_pipeline"]["commit"])
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())

