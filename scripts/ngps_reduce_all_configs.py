#!/usr/bin/env python3

from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
from pathlib import Path


CHANNELS = ("r", "g", "i", "u")


def run(cmd: list[str], label: str, log: Path, cwd: Path | None = None) -> bool:
    """Run a PypeIt command quietly, retaining its complete log."""
    log.parent.mkdir(parents=True, exist_ok=True)
    print(f"  {label} ...", end=" ", flush=True)
    with log.open("w") as stream:
        result = subprocess.run(cmd, cwd=cwd, stdout=stream, stderr=subprocess.STDOUT)
    if result.returncode == 0:
        print("done")
        return True

    print(f"FAILED (details: {log})")
    tail = log.read_text(errors="replace").splitlines()[-12:]
    if tail:
        print("  Last messages from PypeIt:")
        for line in tail:
            print(f"    {line}")
    return False


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
        help="Delete and regenerate the PypeIt setup directories.",
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Pass --overwrite to run_pypeit.",
    )

    parser.add_argument(
        "--auto",
        action="store_true",
        help=(
            "Keep PypeIt's automatic extraction, without opening review windows. "
            "A four-channel extraction-review PDF is still saved for every science exposure."
        ),
    )

    parser.add_argument(
        "--no-review",
        action="store_true",
        help="Reduce only; do not create extraction-review PDFs or open review windows.",
    )

    args = parser.parse_args()
    overwrite_products = args.overwrite or args.auto

    work_root = Path(os.environ.get(
        "NGPS_WORK_ROOT", Path.home() / "ngps_data" / "work"
    )) / args.date
    raw_dir = work_root / "raw"

    if not raw_dir.exists():
        print(f"ERROR: raw directory does not exist:\n{raw_dir}")
        return 1

    logs_dir = work_root / "logs"
    if args.no_review:
        mode = "reduction only (no review PDFs)"
    elif args.auto:
        mode = "automatic reduction and review PDFs (existing products overwritten)"
    else:
        mode = "interactive extraction review"
    print(f"\nNGPS reduction — night {args.date}")
    print(f"Raw data: {raw_dir}")
    print(f"Mode: {mode}")
    print("\n1/3 Checking configurations")

    jobs: list[tuple[str, Path, dict]] = []

    # ---------------------------------------------------------
    # STEP 1: Generate all configurations for every channel
    # ---------------------------------------------------------

    for channel in CHANNELS:

        setup_root = work_root / f"manual_setup_{channel}"

        if args.force_setup and setup_root.exists():
            import shutil

            print(f"  {channel.upper()}: rebuilding configurations")
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
                ],
                label=f"{channel.upper()}: creating configurations",
                log=logs_dir / f"pypeit_setup_{channel}.log",
            )

            if not ok:
                print(f"  {channel.upper()}: skipped because setup creation failed.")
                continue

        # -----------------------------------------------------
        # STEP 2: Discover every generated configuration
        # -----------------------------------------------------

        pypeit_files = [
            path for path in sorted(setup_root.glob("*/*.pypeit"))
            if "_manual_" not in path.parent.name and "_auto_" not in path.parent.name
        ]

        valid_for_channel = 0

        for pf in pypeit_files:

            info = inspect_pypeit_file(pf)

            setup_name = pf.parent.name

            # Require science + arc + flat.
            #
            # This automatically skips your incomplete Setup A,
            # which has science but no matching calibration set.
            if (
                info["has_science"]
                and info["has_arc"]
                and info["has_flat"]
            ):
                jobs.append((channel, pf, info))
                valid_for_channel += 1

            else:
                print(f"  {channel.upper()} {setup_name}: skipped (missing arc or flat calibration)")

        print(
            f"  {channel.upper()}: {valid_for_channel} valid configuration(s) "
            f"of {len(pypeit_files)} found"
        )

    # ---------------------------------------------------------
    # Summary before running expensive reductions
    # ---------------------------------------------------------

    if not jobs:
        print("\nNo valid configurations found; nothing was reduced.")
        return 1

    print(f"\n2/3 Reducing {len(jobs)} valid configuration(s)")

    # ---------------------------------------------------------
    # STEP 3: Run PypeIt sequentially
    # ---------------------------------------------------------

    for number, (channel, pf, info) in enumerate(jobs, start=1):

        progress = f"[{number}/{len(jobs)}] {channel.upper()} {pf.parent.name}"

        science_dir = pf.parent / "Science"

        existing = (
            list(science_dir.glob("spec1d_*.fits"))
            if science_dir.exists()
            else []
        )

        if existing and not overwrite_products:

            print(f"  {progress}: kept {len(existing)} existing spec1d file(s)")
            continue

        cmd = ["run_pypeit", pf.name]

        if overwrite_products:
            cmd.append("--overwrite")

        ok = run(
            cmd,
            label=progress,
            log=logs_dir / f"run_pypeit_{channel}_{pf.parent.name}.log",
            cwd=pf.parent,
        )

        if not ok:
            print(f"  WARNING: {progress} was not reduced.")

    # ---------------------------------------------------------
    # Final summary
    # ---------------------------------------------------------

    print("\nReduction summary")

    total_spec1d = 0

    for channel, pf, info in jobs:

        science_dir = pf.parent / "Science"

        spec1d = (
            sorted(science_dir.glob("spec1d_*.fits"))
            if science_dir.exists()
            else []
        )

        total_spec1d += len(spec1d)

        print(f"  {channel.upper()} {pf.parent.name}: {len(spec1d)} spec1d file(s)")

    print(f"Total: {total_spec1d} spec1d file(s)")

    if not args.no_review:
        reviewer = Path(__file__).with_name("ngps_manual_target_extractions.py")
        review_cmd = [sys.executable, str(reviewer), args.date, "--all"]
        if args.auto:
            review_cmd.append("--auto")
        print("\n3/3 " + ("Saving automatic extraction-review PDFs" if args.auto else "Opening extraction-review dashboards"))
        if args.auto:
            print(f"Review PDFs: {work_root / 'ExtractionQA'}")
        else:
            print("Choose Accept automatic, Manual extraction, or Cancel for each exposure.")
        review_status = subprocess.run(review_cmd).returncode
        if review_status != 0:
            print(f"ERROR: extraction review stopped with code {review_status}.")
            return 1

    if args.no_review:
        print("\nFinished. No review PDFs were created (--no-review).")
    elif args.auto:
        print("\nFinished. Inspect the PDFs above before moving to step 4 (flux calibration).")
    elif not args.no_review:
        print("\nFinished. Flux-calibrate after you have accepted the extractions you want to keep.")
    print(f"Detailed PypeIt logs: {logs_dir}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
