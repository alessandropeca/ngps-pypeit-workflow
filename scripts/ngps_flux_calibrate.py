#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
from pathlib import Path


CHANNELS = ("r", "g", "i", "u")


def run_command(
    command: list[str],
    cwd: Path | None = None,
) -> bool:
    """Run an external command and return True on success."""

    print("\n>>>", " ".join(str(x) for x in command), flush=True)

    result = subprocess.run(
        command,
        cwd=cwd,
    )

    if result.returncode != 0:
        print(
            f"ERROR: command failed with exit code "
            f"{result.returncode}"
        )
        return False

    return True


def find_spec1d(
    science_dir: Path,
    raw_filename: str,
) -> Path | None:
    """
    Find the spec1d corresponding to a raw NGPS filename.

    Example:
        ngps_260623_0097.fits

    matches:
        spec1d_ngps_260623_0097-hz44_....fits
    """

    stem = Path(raw_filename).stem

    matches = sorted(
        science_dir.glob(f"spec1d_{stem}-*.fits")
    )

    if not matches:
        return None

    if len(matches) > 1:
        print(
            f"WARNING: multiple spec1d files found for "
            f"{raw_filename}. Using:"
        )
        print(matches[0])

    return matches[0]


def safe_name(text: str) -> str:
    """Convert target names to filesystem-friendly strings."""

    cleaned = []

    for char in text.strip():
        if char.isalnum() or char in ("-", "_"):
            cleaned.append(char)
        else:
            cleaned.append("_")

    return "".join(cleaned).strip("_")


def main() -> int:

    parser = argparse.ArgumentParser(
        description=(
            "Generate sensitivity functions and flux-calibrate "
            "all reduced NGPS science spectra."
        )
    )

    parser.add_argument(
        "date",
        help="UT observing date in YYYYMMDD format",
    )

    parser.add_argument(
        "--run",
        action="store_true",
        help=(
            "Actually run pypeit_sensfunc and "
            "pypeit_flux_calib. Without this option, "
            "only print the planned associations."
        ),
    )

    parser.add_argument(
        "--force-sensfunc",
        action="store_true",
        help="Regenerate sensitivity functions that already exist.",
    )

    parser.add_argument(
        "--overwrite-fluxed",
        action="store_true",
        help=(
            "Replace existing copies in the Fluxed directories."
        ),
    )

    args = parser.parse_args()

    root = Path(os.environ.get(
        "NGPS_WORK_ROOT", Path.home() / "ngps_data" / "work"
    )) / args.date

    inventory_file = (
        root
        / "science_standard_inventory.csv"
    )

    if not inventory_file.exists():
        print("ERROR: inventory file does not exist:")
        print(inventory_file)
        print()
        print(
            "Run ngps_inventory_standards.py first."
        )
        return 1

    # ---------------------------------------------------------
    # Read inventory
    # ---------------------------------------------------------

    with inventory_file.open() as handle:
        rows = list(
            csv.DictReader(handle)
        )

    for row in rows:
        row["mjd"] = float(row["mjd"])

    total_flux_jobs = 0
    total_sensfuncs = 0

    # ---------------------------------------------------------
    # Process each channel/configuration independently
    # ---------------------------------------------------------

    setups = sorted(
        {
            (row["channel"].lower(), row["setup"])
            for row in rows
            if (
                "science" in row["frametype"].lower()
                or "standard" in row["frametype"].lower()
            )
        }
    )

    for channel, setup in setups:

        setup_rows = [
            row
            for row in rows
            if (
                row["channel"].lower() == channel
                and row["setup"] == setup
            )
        ]

        standards = [
            row
            for row in setup_rows
            if "standard" in row["frametype"].lower()
        ]

        science = [
            row
            for row in setup_rows
            if "science" in row["frametype"].lower()
        ]

        # Skip incomplete setups such as Setup A.
        if not standards or not science:
            print()
            print("=" * 78)
            print(
                f"SKIPPING {setup}: "
                f"science={len(science)}, "
                f"standards={len(standards)}"
            )
            print("=" * 78)
            continue

        setup_dir = (
            root
            / f"manual_setup_{channel}"
            / setup
        )

        science_dir = setup_dir / "Science"

        if not science_dir.exists():
            print(
                f"WARNING: Science directory missing:"
            )
            print(science_dir)
            continue

        sens_dir = setup_dir / "Sensfunc"
        fluxed_dir = setup_dir / "Fluxed"
        flux_files_dir = setup_dir / "FluxFiles"

        sens_dir.mkdir(exist_ok=True)
        fluxed_dir.mkdir(exist_ok=True)
        flux_files_dir.mkdir(exist_ok=True)

        print()
        print("#" * 78)
        print(
            f"{channel.upper()}  {setup}"
        )
        print("#" * 78)

        # -----------------------------------------------------
        # Locate standard spec1d files
        # -----------------------------------------------------

        available_standards = []

        for std in standards:

            spec1d = find_spec1d(
                science_dir,
                std["filename"],
            )

            if spec1d is None:
                print(
                    f"WARNING: no spec1d found for standard "
                    f"{std['filename']} "
                    f"({std['target']})"
                )
                continue

            std_name = safe_name(std["target"])

            raw_stem = Path(
                std["filename"]
            ).stem

            sensfile = (
                sens_dir
                / (
                    f"sens_{channel.upper()}_"
                    f"{setup}_"
                    f"{std_name}_"
                    f"{raw_stem}.fits"
                )
            )

            std_info = dict(std)
            std_info["spec1d"] = spec1d
            std_info["sensfile"] = sensfile

            available_standards.append(
                std_info
            )

        if not available_standards:
            print(
                "WARNING: no reduced standard-star "
                "spec1d files found."
            )
            continue

        # -----------------------------------------------------
        # Build one sensitivity function per standard exposure
        # -----------------------------------------------------

        for std in available_standards:

            sensfile = std["sensfile"]

            print()
            print(
                f"STANDARD: {std['target']}"
            )
            print(
                f"  MJD:      {std['mjd']:.6f}"
            )
            print(
                f"  spec1d:   {std['spec1d']}"
            )
            print(
                f"  sensfunc: {sensfile}"
            )

            if (
                sensfile.exists()
                and not args.force_sensfunc
            ):
                print(
                    "  Sensitivity function already exists."
                )

            elif args.run:

                # Do not force UVIS/IR here.
                #
                # Let the NGPS PypeIt spectrograph defaults
                # determine the sensitivity-function algorithm.
                ok = run_command(
                    [
                        "pypeit_sensfunc",
                        str(std["spec1d"]),
                        "-o",
                        str(sensfile),
                    ],
                    cwd=setup_dir,
                )

                if not ok:
                    print(
                        "  WARNING: sensfunc generation failed."
                    )
                    continue

            total_sensfuncs += 1

        # Keep only standards with usable sensitivity functions
        # when actually running.
        if args.run:

            available_standards = [
                std
                for std in available_standards
                if std["sensfile"].exists()
            ]

            if not available_standards:
                print(
                    "No usable sensitivity functions "
                    "for this setup."
                )
                continue

        # -----------------------------------------------------
        # Match each science exposure to nearest standard
        # -----------------------------------------------------

        assignments = []

        for sci in science:

            science_spec1d = find_spec1d(
                science_dir,
                sci["filename"],
            )

            if science_spec1d is None:
                print(
                    f"WARNING: no spec1d found for science "
                    f"{sci['filename']}"
                )
                continue

            nearest = min(
                available_standards,
                key=lambda std: abs(
                    std["mjd"] - sci["mjd"]
                ),
            )

            delta_days = abs(
                nearest["mjd"]
                - sci["mjd"]
            )

            delta_minutes = (
                delta_days
                * 24.0
                * 60.0
            )

            assignments.append(
                {
                    "science": sci,
                    "science_spec1d": science_spec1d,
                    "standard": nearest,
                    "delta_minutes": delta_minutes,
                }
            )

        # -----------------------------------------------------
        # Print association table
        # -----------------------------------------------------

        print()
        print(
            "SCIENCE → STANDARD ASSOCIATIONS"
        )
        print("-" * 78)

        for item in assignments:

            sci = item["science"]
            std = item["standard"]

            print(
                f"{sci['target']:20s}  "
                f"{sci['filename']:25s}  "
                f"→ {std['target']:12s}  "
                f"Δt={item['delta_minutes']:6.1f} min"
            )

        if not assignments:
            continue

        # -----------------------------------------------------
        # Copy science spec1d files into Fluxed/
        # -----------------------------------------------------

        flux_rows = []

        for item in assignments:

            source = item["science_spec1d"]

            destination = (
                fluxed_dir
                / source.name
            )

            if (
                not destination.exists()
                or args.overwrite_fluxed
            ):
                shutil.copy2(
                    source,
                    destination,
                )

            flux_rows.append(
                (
                    destination,
                    item["standard"]["sensfile"],
                )
            )

        # -----------------------------------------------------
        # Write PypeIt .flux file
        # -----------------------------------------------------

        flux_file = (
            flux_files_dir
            / f"{channel}_{setup}.flux"
        )

        with flux_file.open("w") as handle:

            handle.write(
                "# Auto-generated NGPS flux calibration file\n"
            )

            handle.write(
                f"# Date: {args.date}\n"
            )

            handle.write(
                f"# Channel: {channel.upper()}\n"
            )

            handle.write(
                f"# Setup: {setup}\n\n"
            )

            handle.write("flux read\n")
            handle.write(
                "    filename | sensfile\n"
            )

            for science_file, sensfile in flux_rows:

                handle.write(
                    f"    {science_file} | {sensfile}\n"
                )

            handle.write("flux end\n")

        print()
        print(
            f"Flux file written:"
        )
        print(flux_file)

        total_flux_jobs += len(flux_rows)

        # -----------------------------------------------------
        # Apply flux calibration
        # -----------------------------------------------------

        if args.run:

            run_command(
                [
                    "pypeit_flux_calib",
                    str(flux_file),
                ],
                cwd=setup_dir,
            )

    # ---------------------------------------------------------
    # Final summary
    # ---------------------------------------------------------

    print()
    print("=" * 78)

    if args.run:
        print("FLUX CALIBRATION RUN COMPLETE")
    else:
        print("DRY RUN COMPLETE")

    print("=" * 78)

    print(
        f"Potential sensitivity functions: "
        f"{total_sensfuncs}"
    )

    print(
        f"Science spectra assigned: "
        f"{total_flux_jobs}"
    )

    if not args.run:

        print()
        print(
            "No files were flux calibrated."
        )

        print(
            "Review the associations above, then run:"
        )

        print()
        print(
            f"python {Path(__file__)} "
            f"{args.date} --run"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
