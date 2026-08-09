#!/usr/bin/env python3

from __future__ import annotations

import argparse
import os
import re
import subprocess
from pathlib import Path


CHANNELS = ("r", "g", "i", "u")


def run(cmd: list[str], cwd: Path | None = None) -> bool:
    """Run a command and return True on success."""
    print("\n>>>", " ".join(cmd), flush=True)

    result = subprocess.run(cmd, cwd=cwd)

    if result.returncode != 0:
        print(f"ERROR: command failed with code {result.returncode}")
        return False

    return True


def inspect_pypeit_file(path: Path) -> dict:
    """Inspect active rows in a generated PypeIt file."""

    text = path.read_text()

    match = re.search(r"binning:\s*([^\n]+)", text)
    binning = match.group(1).strip() if match else "unknown"

    frame_types = set()
    in_data_block = False

    for line in text.splitlines():

        stripped = line.strip()

        if stripped == "data read":
            in_data_block = True
            continue

        if stripped == "data end":
            in_data_block = False
            continue

        if not in_data_block:
            continue

        if not stripped or stripped.startswith("#"):
            continue

        # Actual PypeIt table rows contain pipe-separated columns.
        parts = line.split("|")

        if len(parts) < 2:
            continue

        frametype = parts[1].strip().lower()

        if not frametype or frametype == "frametype":
            continue

        for ft in frametype.split(","):
            frame_types.add(ft.strip())

    return {
        "binning": binning,
        "has_science": "science" in frame_types,
        "has_arc": "arc" in frame_types,
        "has_flat": any(
            x in frame_types
            for x in ("pixelflat", "illumflat", "trace")
        ),
        "has_standard": "standard" in frame_types,
        "has_bias": "bias" in frame_types,
    }
  


def main() -> int:

    parser = argparse.ArgumentParser(
        description=(
            "Reduce all valid NGPS PypeIt configurations "
            "for all four channels."
        )
    )

    parser.add_argument(
        "date",
        help="UT date, e.g. 20260623",
    )

    parser.add_argument(
        "--force-setup",
        action="store_true",
        help="Delete and regenerate manual setup directories.",
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Pass --overwrite to run_pypeit.",
    )

    args = parser.parse_args()

    work_root = Path(os.environ.get(
        "NGPS_WORK_ROOT", Path.home() / "ngps_data" / "work"
    )) / args.date
    raw_dir = work_root / "raw"

    if not raw_dir.exists():
        print(f"ERROR: raw directory does not exist:\n{raw_dir}")
        return 1

    print(f"\nNGPS night: {args.date}")
    print(f"Raw data:   {raw_dir}")

    jobs: list[tuple[str, Path, dict]] = []

    # ---------------------------------------------------------
    # STEP 1: Generate all configurations for every channel
    # ---------------------------------------------------------

    for channel in CHANNELS:

        setup_root = work_root / f"manual_setup_{channel}"

        if args.force_setup and setup_root.exists():
            import shutil

            print(f"\nRemoving old setup directory: {setup_root}")
            shutil.rmtree(setup_root)

        if not setup_root.exists():

            ok = run(
                [
                    "pypeit_setup",
                    "-s",
                    f"p200_ngps_{channel}",
                    "-r",
                    str(raw_dir),
                    "-d",
                    str(setup_root),
                    "-c",
                    "all",
                ]
            )

            if not ok:
                print(f"Skipping channel {channel}: setup failed.")
                continue

        # -----------------------------------------------------
        # STEP 2: Discover every generated configuration
        # -----------------------------------------------------

        pypeit_files = sorted(setup_root.glob("*/*.pypeit"))

        print(
            f"\n{'=' * 70}\n"
            f"CHANNEL {channel.upper()}: "
            f"{len(pypeit_files)} configuration(s)\n"
            f"{'=' * 70}"
        )

        for pf in pypeit_files:

            info = inspect_pypeit_file(pf)

            setup_name = pf.parent.name

            print(
                f"\n{setup_name}"
                f"\n  binning:  {info['binning']}"
                f"\n  science:  {info['has_science']}"
                f"\n  arc:      {info['has_arc']}"
                f"\n  flat:     {info['has_flat']}"
                f"\n  bias:     {info['has_bias']}"
                f"\n  standard: {info['has_standard']}"
            )

            # Require science + arc + flat.
            #
            # This automatically skips your incomplete Setup A,
            # which has science but no matching calibration set.
            if (
                info["has_science"]
                and info["has_arc"]
                and info["has_flat"]
            ):
                print("  --> VALID: will reduce")
                jobs.append((channel, pf, info))

            else:
                print("  --> SKIP: incomplete science/calibration setup")

    # ---------------------------------------------------------
    # Summary before running expensive reductions
    # ---------------------------------------------------------

    print("\n")
    print("=" * 70)
    print("CONFIGURATIONS SELECTED FOR REDUCTION")
    print("=" * 70)

    for channel, pf, info in jobs:
        print(
            f"{channel.upper():2s}  "
            f"{pf.parent.name:20s}  "
            f"binning={info['binning']}"
        )

    if not jobs:
        print("\nNo valid configurations found.")
        return 1

    # ---------------------------------------------------------
    # STEP 3: Run PypeIt sequentially
    # ---------------------------------------------------------

    for number, (channel, pf, info) in enumerate(jobs, start=1):

        print("\n")
        print("#" * 70)
        print(
            f"REDUCTION {number}/{len(jobs)}"
            f"  channel={channel.upper()}"
            f"  setup={pf.parent.name}"
            f"  binning={info['binning']}"
        )
        print("#" * 70)

        science_dir = pf.parent / "Science"

        existing = (
            list(science_dir.glob("spec1d_*.fits"))
            if science_dir.exists()
            else []
        )

        if existing and not args.overwrite:

            print(
                f"Found {len(existing)} existing spec1d files."
                " Skipping this setup."
            )
            continue

        cmd = ["run_pypeit", pf.name]

        if args.overwrite:
            cmd.append("--overwrite")

        ok = run(cmd, cwd=pf.parent)

        if not ok:
            print(
                f"\nWARNING: reduction failed for "
                f"{pf.parent.name}."
            )

    # ---------------------------------------------------------
    # Final summary
    # ---------------------------------------------------------

    print("\n")
    print("=" * 70)
    print("REDUCTION SUMMARY")
    print("=" * 70)

    total_spec1d = 0

    for channel, pf, info in jobs:

        science_dir = pf.parent / "Science"

        spec1d = (
            sorted(science_dir.glob("spec1d_*.fits"))
            if science_dir.exists()
            else []
        )

        total_spec1d += len(spec1d)

        print(
            f"{channel.upper():2s}  "
            f"{pf.parent.name:20s}  "
            f"binning={info['binning']:8s}  "
            f"spec1d={len(spec1d)}"
        )

    print(f"\nTotal spec1d files: {total_spec1d}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
