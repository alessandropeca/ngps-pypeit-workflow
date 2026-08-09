#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import os
import re
from pathlib import Path


CHANNELS = ("r", "g", "i", "u")


def main() -> int:

    parser = argparse.ArgumentParser(
        description=(
            "Create a chronological inventory of NGPS science "
            "and standard-star exposures for all channels and configurations."
        )
    )

    parser.add_argument(
        "date",
        help="UT observing date in YYYYMMDD format, e.g. 20260623",
    )

    args = parser.parse_args()

    date = args.date

    root = Path(os.environ.get(
        "NGPS_WORK_ROOT", Path.home() / "ngps_data" / "work"
    )) / date

    if not root.exists():
        print(f"ERROR: observing-night directory does not exist:")
        print(root)
        return 1

    rows = []

    for channel in CHANNELS:

        setup_root = root / f"manual_setup_{channel}"

        if not setup_root.exists():
            print(
                f"WARNING: setup directory not found: "
                f"{setup_root}"
            )
            continue

        for pypeit_file in sorted(
            setup_root.glob("*/*.pypeit")
        ):

            text = pypeit_file.read_text()

            # Configuration/binning
            match = re.search(
                r"binning:\s*([^\n]+)",
                text,
            )

            binning = (
                match.group(1).strip()
                if match
                else "unknown"
            )

            setup = pypeit_file.parent.name

            in_data = False

            for line in text.splitlines():

                stripped = line.strip()

                if stripped == "data read":
                    in_data = True
                    continue

                if stripped == "data end":
                    in_data = False
                    continue

                if not in_data:
                    continue

                # Ignore empty and commented-out rows.
                if (
                    not stripped
                    or stripped.startswith("#")
                ):
                    continue

                parts = [
                    value.strip()
                    for value in line.split("|")
                ]

                # Expected columns:
                #
                # filename
                # frametype
                # ra
                # dec
                # target
                # dispname
                # decker
                # binning
                # mjd
                # airmass
                # exptime
                # calib

                if len(parts) < 12:
                    continue

                filename = parts[0]
                frametype = parts[1].lower()

                # Ignore the table header.
                if filename.lower() == "filename":
                    continue

                # Keep only science and standards.
                if (
                    "science" not in frametype
                    and "standard" not in frametype
                ):
                    continue

                try:
                    mjd = float(parts[8])
                except ValueError:
                    continue

                try:
                    airmass = float(parts[9])
                except ValueError:
                    airmass = float("nan")

                try:
                    exptime = float(parts[10])
                except ValueError:
                    exptime = float("nan")

                rows.append(
                    {
                        "channel": channel.upper(),
                        "setup": setup,
                        "binning": binning,
                        "frametype": frametype,
                        "filename": filename,
                        "target": parts[4],
                        "mjd": mjd,
                        "airmass": airmass,
                        "exptime": exptime,
                    }
                )

    # Sort by channel, setup, then time.
    rows.sort(
        key=lambda row: (
            row["channel"],
            row["setup"],
            row["mjd"],
        )
    )

    # ---------------------------------------------------------
    # Print readable summary
    # ---------------------------------------------------------

    for channel in CHANNELS:

        channel_upper = channel.upper()

        print()
        print("=" * 90)
        print(f"CHANNEL {channel_upper}")
        print("=" * 90)

        setups = sorted(
            {
                row["setup"]
                for row in rows
                if row["channel"] == channel_upper
            }
        )

        for setup in setups:

            subset = [
                row
                for row in rows
                if (
                    row["channel"] == channel_upper
                    and row["setup"] == setup
                )
            ]

            if not subset:
                continue

            print()
            print(
                f"{setup}   "
                f"binning={subset[0]['binning']}"
            )

            print("-" * 90)

            for row in subset:

                kind = (
                    "STD"
                    if "standard" in row["frametype"]
                    else "SCI"
                )

                print(
                    f"{kind:3s}  "
                    f"MJD={row['mjd']:.6f}  "
                    f"airmass={row['airmass']:.3f}  "
                    f"{row['target']:20s}  "
                    f"{row['filename']}"
                )

    # ---------------------------------------------------------
    # Write CSV
    # ---------------------------------------------------------

    output = (
        root
        / "science_standard_inventory.csv"
    )

    with output.open(
        "w",
        newline="",
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "channel",
                "setup",
                "binning",
                "frametype",
                "filename",
                "target",
                "mjd",
                "airmass",
                "exptime",
            ],
        )

        writer.writeheader()
        writer.writerows(rows)

    print()
    print("=" * 90)
    print(f"Inventory written to:")
    print(output)

    print(
        f"\nTotal science/standard rows: "
        f"{len(rows)}"
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
