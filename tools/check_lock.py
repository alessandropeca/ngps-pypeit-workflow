#!/usr/bin/env python3
"""Validate immutable upstream pins without network access."""

from __future__ import annotations

from pathlib import Path
import re
import sys
import yaml

lock = yaml.safe_load((Path(__file__).parents[1] / "upstream-lock.yml").read_text())
for name, item in lock["upstreams"].items():
    if not re.fullmatch(r"[0-9a-f]{40}", item["commit"]):
        sys.exit(f"{name}: commit must be a full 40-character lowercase SHA")
    if not item["repository"].startswith("https://github.com/"):
        sys.exit(f"{name}: repository must be an explicit GitHub URL")
print("upstream-lock.yml is structurally valid")

