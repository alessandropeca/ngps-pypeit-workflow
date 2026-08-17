#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path

import numpy as np
from astropy.io import fits


def find_spec1d(science_dir: Path, raw_filename: str):
    stem = Path(raw_filename).stem
    return sorted(science_dir.glob(f"spec1d_{stem}-*.fits"))


def flux_quality(path: Path) -> tuple[bool, bool]:
    """Return whether a file has FLAM and at least one usable FLAM uncertainty."""
    has_flam = False
    has_valid_ivar = False
    try:
        with fits.open(path) as hdul:
            for hdu in hdul[1:]:
                if not hasattr(hdu, "columns") or hdu.columns is None:
                    continue

                names = hdu.columns.names or []

                if "OPT_FLAM" in names or "BOX_FLAM" in names:
                    has_flam = True
                for flux_name, ivar_name in (("OPT_FLAM", "OPT_FLAM_IVAR"), ("BOX_FLAM", "BOX_FLAM_IVAR")):
                    if flux_name not in names or ivar_name not in names:
                        continue
                    flux = np.asarray(hdu.data[flux_name], dtype=float)
                    ivar = np.asarray(hdu.data[ivar_name], dtype=float)
                    if np.any(np.isfinite(flux) & np.isfinite(ivar) & (ivar > 0)):
                        has_valid_ivar = True

    except Exception as exc:
        print(f"ERROR reading {path}: {exc}")

    return has_flam, has_valid_ivar


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
    valid_uncertainty = 0

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

            quality = [flux_quality(path) for path in files]
            good = sum(has_flam for has_flam, _ in quality)
            usable = sum(has_ivar for _, has_ivar in quality)

            total_fluxed += len(files)
            successfully_fluxed += good
            valid_uncertainty += usable

            print(
                f"{channel.upper():2s}  "
                f"{fluxed_dir.parent.name:20s}  "
                f"files={len(files):3d}  "
                f"with_FLAM={good:3d}  "
                f"with_valid_FLAM_IVAR={usable:3d}"
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

    print(
        f"Files with valid FLAM IVAR: "
        f"{valid_uncertainty}"
    )

    if valid_uncertainty != total_fluxed:
        print("WARNING: Some fluxed files have no usable FLAM uncertainties.")
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
