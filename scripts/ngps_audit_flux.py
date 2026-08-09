#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path

from astropy.io import fits


def find_spec1d(science_dir: Path, raw_filename: str):
    stem = Path(raw_filename).stem
    return sorted(science_dir.glob(f"spec1d_{stem}-*.fits"))


def has_flux_columns(path: Path) -> bool:
    try:
        with fits.open(path) as hdul:
            for hdu in hdul[1:]:
                if not hasattr(hdu, "columns") or hdu.columns is None:
                    continue

                names = hdu.columns.names or []

                if "OPT_FLAM" in names or "BOX_FLAM" in names:
                    return True

    except Exception as exc:
        print(f"ERROR reading {path}: {exc}")

    return False


def main() -> int:

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "date",
        help="UT observing date, e.g. 20260623",
    )

    args = parser.parse_args()

    root = Path(os.environ.get(
        "NGPS_WORK_ROOT", Path.home() / "ngps_data" / "work"
    )) / args.date

    inventory = root / "science_standard_inventory.csv"

    with inventory.open() as handle:
        rows = list(csv.DictReader(handle))

    missing_science = []
    missing_standards = []

    total_fluxed = 0
    successfully_fluxed = 0

    print()
    print("=" * 80)
    print("CHECKING REDUCED SPEC1D FILES")
    print("=" * 80)

    for row in rows:

        channel = row["channel"].lower()
        setup = row["setup"]
        frametype = row["frametype"].lower()

        # Ignore Setup A because it was intentionally not reduced.
        if setup.endswith("_A"):
            continue

        science_dir = (
            root
            / f"manual_setup_{channel}"
            / setup
            / "Science"
        )

        matches = find_spec1d(
            science_dir,
            row["filename"],
        )

        if not matches:

            item = (
                row["channel"],
                setup,
                row["frametype"],
                row["target"],
                row["filename"],
            )

            if "standard" in frametype:
                missing_standards.append(item)

            elif "science" in frametype:
                missing_science.append(item)

    print()
    print("MISSING STANDARD SPEC1D FILES")
    print("-" * 80)

    if not missing_standards:
        print("None")
    else:
        for item in missing_standards:
            print(
                f"{item[0]:2s}  "
                f"{item[1]:20s}  "
                f"{item[3]:20s}  "
                f"{item[4]}"
            )

    print()
    print("MISSING SCIENCE SPEC1D FILES")
    print("-" * 80)

    if not missing_science:
        print("None")
    else:
        for item in missing_science:
            print(
                f"{item[0]:2s}  "
                f"{item[1]:20s}  "
                f"{item[3]:20s}  "
                f"{item[4]}"
            )

    print()
    print("=" * 80)
    print("CHECKING FLUXED FILES")
    print("=" * 80)

    for channel in ("r", "g", "i", "u"):

        setup_root = (
            root
            / f"manual_setup_{channel}"
        )

        for fluxed_dir in sorted(
            setup_root.glob("*/Fluxed")
        ):

            files = sorted(
                fluxed_dir.glob("spec1d_*.fits")
            )

            if not files:
                continue

            good = sum(
                has_flux_columns(path)
                for path in files
            )

            total_fluxed += len(files)
            successfully_fluxed += good

            print(
                f"{channel.upper():2s}  "
                f"{fluxed_dir.parent.name:20s}  "
                f"files={len(files):3d}  "
                f"with_FLAM={good:3d}"
            )

    print()
    print("=" * 80)
    print("SUMMARY")
    print("=" * 80)

    print(
        f"Missing standard spec1d: "
        f"{len(missing_standards)}"
    )

    print(
        f"Missing science spec1d:  "
        f"{len(missing_science)}"
    )

    print(
        f"Fluxed-directory files:  "
        f"{total_fluxed}"
    )

    print(
        f"Files containing FLAM:   "
        f"{successfully_fluxed}"
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
