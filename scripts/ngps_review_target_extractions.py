#!/usr/bin/env python3
"""Review every reduced 2D science exposure for one NGPS target.

The review runs one full 2D frame at a time, across the target's U/G/R/I
channels and repeated observations.  It uses the same accept-automatic or
manual-click decision as ngps_interactive_extract.py.
"""

from __future__ import annotations

import argparse
import csv
import os
import subprocess
from pathlib import Path

from ngps_interactive_extract import (
    ask_yes_no,
    create_manual_copy,
    exposure_from_spec2d,
    interactive_select,
    manual_value,
)


def find_spec2d(directory: Path, raw_filename: str) -> Path | None:
    matches = sorted(directory.glob(f"spec2d_{Path(raw_filename).stem}-*.fits"))
    return matches[0] if matches else None


def matching_rows(
    rows: list[dict[str, str]], target: str, channel: str | None, setup: str | None
) -> list[dict[str, str]]:
    matches = []
    for row in rows:
        if "science" not in row["frametype"].lower():
            continue
        if row["target"].casefold() != target.casefold():
            continue
        if channel is not None and row["channel"].lower() != channel:
            continue
        if setup is not None and row["setup"] != setup:
            continue
        matches.append(row)
    return sorted(matches, key=lambda row: (row["channel"], row["setup"], row["mjd"]))


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Review all full 2D NGPS science frames for one target."
    )
    parser.add_argument("date", help="UT date, e.g. 20260623")
    parser.add_argument("--target", required=True, help="Target name from the inventory")
    parser.add_argument("--channel", choices=("u", "g", "r", "i"), help="Optional channel filter")
    parser.add_argument("--setup", help="Optional PypeIt setup filter")
    parser.add_argument("--fwhm", type=float, default=4.0, help="Initial manual FWHM in pixels")
    parser.add_argument("--max-select", type=int, default=3, help="Manual positions allowed: 1 to 3")
    parser.add_argument("--summary", action="store_true", help="List 2D frames without opening them")
    args = parser.parse_args()

    if not 1 <= args.max_select <= 3 or args.fwhm <= 0:
        parser.error("--max-select must be 1 to 3 and --fwhm must be positive")

    root = Path(os.environ.get("NGPS_WORK_ROOT", Path.home() / "ngps_data" / "work")) / args.date
    inventory = root / "science_standard_inventory.csv"
    if not inventory.is_file():
        parser.error(f"Inventory not found: {inventory}; run ngps_inventory_standards.py first")
    with inventory.open() as handle:
        rows = matching_rows(list(csv.DictReader(handle)), args.target, args.channel, args.setup)
    if not rows:
        parser.error("No matching science exposures found")

    review_items: list[tuple[dict[str, str], Path]] = []
    for row in rows:
        science = root / f"manual_setup_{row['channel'].lower()}" / row["setup"] / "Science"
        spec2d = find_spec2d(science, row["filename"])
        if spec2d is None:
            print(f"WARNING: no reduced spec2d for {row['filename']}")
            continue
        review_items.append((row, spec2d))

    print(f"Found {len(review_items)} 2D frame(s) for {args.target}.")
    for index, (row, spec2d) in enumerate(review_items, start=1):
        print(f"{index:2d}. {row['channel']}  {row['setup']}  {row['filename']}  {spec2d.name}")
    if args.summary:
        return 0

    for index, (row, spec2d) in enumerate(review_items, start=1):
        print(
            f"\n{'=' * 78}\n"
            f"{index}/{len(review_items)}  {row['target']} | {row['channel']} | "
            f"{row['setup']} | {row['filename']}\n"
            f"{'=' * 78}"
        )
        result = interactive_select(spec2d, args.fwhm, args.max_select)
        if result.decision == "automatic":
            print("Automatic extraction accepted; no files changed.")
            continue
        if result.decision != "manual":
            print("Review cancelled for this exposure; no files changed.")
            continue

        source_pypeit = next(iter(sorted(spec2d.parent.parent.glob("*.pypeit"))), None)
        if source_pypeit is None:
            print("WARNING: no source .pypeit file found; manual selection was not written.")
            continue
        print("PypeIt manual value:\n  " + manual_value(result.selections))
        if not ask_yes_no("Write this to a copied manual PypeIt setup?"):
            print("No files changed.")
            continue
        try:
            manual_dir, manual_pypeit = create_manual_copy(
                source_pypeit, exposure_from_spec2d(spec2d), result.selections
            )
        except FileExistsError as error:
            print(f"WARNING: {error}")
            continue
        print(f"Manual setup: {manual_dir}\nManual PypeIt file: {manual_pypeit}")
        if ask_yes_no("Run PypeIt on the copied manual setup now?"):
            status = subprocess.run(["run_pypeit", manual_pypeit.name], cwd=manual_dir).returncode
            if status != 0:
                return status

    print("\nTarget extraction review complete.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
