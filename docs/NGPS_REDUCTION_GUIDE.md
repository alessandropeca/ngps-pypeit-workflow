NGPS REDUCTION GUIDE — PINNED INSTALLATION THROUGH EXTRACTED 1D SPECTRA
===============================================================================

PROJECT BASELINE
----------------
This is the maintained guide for the ngps-pypeit-workflow project.  It uses
fixed Git commits, never a floating upstream branch.  Before reducing data, set:

    WORKFLOW_ROOT="$HOME/Documents/GitHub/ngps-pypeit-workflow"

The exact versions are recorded in $WORKFLOW_ROOT/upstream-lock.yml.  The NGPS
wrapper is downloaded from Alessandro Peca's maintained fork; the original
Christoffer Fremling repository is retained only as the documented upstream.



This guide reproduces the workflow that successfully reduced the 20260623 NGPS
night through wavelength-calibrated, sky-subtracted, extracted 1D spectra.

WORKING SETUP
-------------
Software directory:
    ~/Software

Repositories:
    ~/Software/PypeIt
    ~/Software/ngps_pipeline

Both are required and are installed separately. PypeIt is the reduction engine
and provides the NGPS instrument support; ngps_pipeline is the operational
wrapper that calls PypeIt. The wrapper does not include PypeIt.

Conda environment:
    ngps

Python:
    3.11

Raw-data work directory:
    ~/ngps_data/work/<DATE>/raw

Example date used here:
    20260623

IMPORTANT RESULT FROM THIS NIGHT
--------------------------------
The standard NGPS wrapper command:

    python -m ngps_pipeline.reduce 20260623 --skip-db-import

was not sufficient for this night because the data contained multiple PypeIt
configurations distinguished by detector binning.

For each channel, PypeIt found four configurations:

    Setup A: binning 1,2
    Setup B: binning 3,2
    Setup C: binning 3,4
    Setup D: binning 4,4

Setup A contained science but no matching calibration set. Setups B, C, and D
contained science plus calibration frames. The successful solution was therefore
to run pypeit_setup with -c all for every channel and reduce all valid setups.

===============================================================================
1. OPTIONAL CLEAN START
===============================================================================

Deactivate any current Conda environment:

    conda deactivate

Remove old dedicated NGPS/PypeIt environments if needed:

    conda env remove -n ngps -y
    conda env remove -n pypeit -y

Remove old clones only if you truly want a fresh installation:

    rm -rf ~/Software/ngps_pipeline
    rm -rf ~/Software/PypeIt

Optional: remove old PypeIt cached/reference data:

    rm -rf ~/.pypeit

Check remaining environments:

    conda env list

Check old PypeIt executables:

    hash -r
    rehash 2>/dev/null

    which -a pypeit_setup
    which -a run_pypeit
    which -a pypeit_show_1dspec

===============================================================================
2. CREATE THE SOFTWARE DIRECTORY
===============================================================================

    mkdir -p ~/Software
    cd ~/Software
    pwd

Expected:

    /Users/xpecax/Software

===============================================================================
3. CLONE THE REPOSITORIES
===============================================================================

Clone the NGPS wrapper:

    git clone https://github.com/alessandropeca/ngps_pipeline.git
    git -C ~/Software/ngps_pipeline checkout 55fa9491eb1683769006118c46b26963bbf33ea2

Clone the NGPS-enabled PypeIt fork and check out its validated commit:

    git clone https://github.com/cfremling/PypeIt.git
    git -C ~/Software/PypeIt checkout e9ed85c1a237c49626227f4227e323fc390def4b

Check:

    ls -ld ~/Software/ngps_pipeline ~/Software/PypeIt

===============================================================================
4. CREATE THE CONDA ENVIRONMENT
===============================================================================

    conda create -n ngps python=3.11 -y
    conda activate ngps

Check:

    which python
    python --version

Expected Python path:

    /opt/anaconda3/envs/ngps/bin/python

===============================================================================
5. UPDATE INSTALLATION TOOLS
===============================================================================

    python -m pip install --upgrade pip setuptools wheel

===============================================================================
6. INSTALL PYPEIT AND THE NGPS WRAPPER
===============================================================================

    python -m pip install -e ~/Software/PypeIt
    python -m pip install -e ~/Software/ngps_pipeline

===============================================================================
7. VERIFY THE INSTALLATION
===============================================================================

    python -c "import sys, pypeit; print('Python:', sys.executable); print('PypeIt:', pypeit.__file__); print('PypeIt version:', pypeit.__version__)"

Expected paths should resemble:

    /opt/anaconda3/envs/ngps/bin/python
    /Users/xpecax/Software/PypeIt/pypeit/__init__.py

Check executables:

    which pypeit_setup
    which run_pypeit
    which pypeit_sensfunc
    which pypeit_flux_calib
    which pypeit_show_1dspec

They should point into:

    /opt/anaconda3/envs/ngps/bin/

===============================================================================
8. VERIFY NGPS SUPPORT
===============================================================================

Do not use pypeit_setup --list with this development version.

Instead:

    python - <<'PY'
from pypeit.spectrographs import p200_ngps
print("NGPS module:", p200_ngps.__file__)
print("NGPS class:", p200_ngps.P200NGPSSpectrograph)
PY

or, in one line:

python -c "from pypeit.spectrographs import p200_ngps; print('NGPS module:', p200_ngps.__file__); print('NGPS class:', p200_ngps.P200NGPSSpectrograph)"

Expected module path:

    /Users/xpecax/Software/PypeIt/pypeit/spectrographs/p200_ngps.py

Check the NGPS wrapper:

    python -c "import ngps_pipeline.reduce; print(ngps_pipeline.reduce.__file__)"

Expected:

    /Users/xpecax/Software/ngps_pipeline/ngps_pipeline/reduce.py

Check options:

    python -m ngps_pipeline.reduce --help

===============================================================================
9. OPTIONAL WRAPPER CONFIG FILE
===============================================================================

The wrapper can work without this file because its default work directory is:

    ~/ngps_data/work

To create it explicitly:

    mkdir -p ~/.config

    cat > ~/.config/ngps_pipeline.toml <<'CFGEOF'
work_dir = "/Users/xpecax/ngps_data/work"
db_path = "/Users/xpecax/ngps_data/ngps_db.sqlite"

[palomar]
longitude_deg = -116.8639
latitude_deg = 33.3563
altitude_m = 1712.0
CFGEOF

===============================================================================
10. PREPARE THE RAW-DATA DIRECTORY
===============================================================================

Set the observing date:

    DATE=20260623

Create the raw directory:

    mkdir -p ~/ngps_data/work/${DATE}/raw

Copy the original NGPS FITS files:

    RAW_DIR="/PATH/TO/RAW/FILES"
    rsync -av "${RAW_DIR}/"*.fits ~/ngps_data/work/${DATE}/raw/

Check the count:

    find ~/ngps_data/work/${DATE}/raw -maxdepth 1 -name "*.fits" | wc -l

For 20260623 there were 151 raw FITS files.

Keep the original NGPS FITS files unchanged. Do not manually split the U, G, R,
and I extensions.

===============================================================================
11. WHY THE STANDARD WRAPPER FAILED FOR THIS NIGHT
===============================================================================

The first attempt was:

    python -m ngps_pipeline.reduce 20260623 --skip-db-import

The wrapper generated setup files but all four run_pypeit jobs failed.
The detailed error was found with:

    tail -n 150 ~/ngps_data/work/20260623/logs/run_r.log

The fatal error was:

    PypeItError: No frames of type=arc provided.

The generated Setup A contained:

    binning: 1,2

and one science frame, but no matching arc, flat, bias, or standard frames.

Running pypeit_setup with -c all showed that the night actually contained four
configurations, and the complete calibration sets were in B, C, and D.

===============================================================================
12. OPTIONAL MANUAL DIAGNOSTIC FOR ONE CHANNEL
===============================================================================

Example for R:

    rm -rf ~/ngps_data/work/20260623/manual_setup_r

    pypeit_setup \
        -s p200_ngps_r \
        -r ~/ngps_data/work/20260623/raw \
        -d ~/ngps_data/work/20260623/manual_setup_r \
        -c all

List generated PypeIt files:

    find ~/ngps_data/work/20260623/manual_setup_r -name "*.pypeit" -print

Inspect their frame assignments:

    for pf in ~/ngps_data/work/20260623/manual_setup_r/*/*.pypeit; do
        echo
        echo "================================================"
        echo "$pf"
        echo "================================================"
        grep -E "Setup |arc|tilt|pixelflat|trace|illumflat|bias|science|standard" "$pf"
    done

For this night:

    A = 1,2   incomplete calibration set
    B = 3,2   valid
    C = 3,4   valid
    D = 4,4   valid

===============================================================================
13. OPTIONAL TEST OF ONE CONFIGURATION
===============================================================================

R/Setup B was tested manually first:

    cd ~/ngps_data/work/20260623/manual_setup_r/p200_ngps_r_B
    run_pypeit p200_ngps_r_B.pypeit

This completed successfully and produced spec1d and spec2d products.

===============================================================================
14. AUTOMATE ALL VALID CONFIGURATIONS
===============================================================================

Create:

    Do not create a new local script. The maintained version is:

    /scripts/ngps_reduce_all_configs.py

The historical inline copy below is retained only as a record of the 20260623
solution; do not edit or run it from ~/Software.

Paste the following script:

-------------------------------------------------------------------------------
BEGIN SCRIPT
-------------------------------------------------------------------------------

#!/usr/bin/env python3

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
from pathlib import Path


CHANNELS = ("r", "g", "i", "u")


def run(cmd: list[str], cwd: Path | None = None) -> bool:
    print("\n>>>", " ".join(cmd), flush=True)
    result = subprocess.run(cmd, cwd=cwd)
    if result.returncode != 0:
        print(f"ERROR: command failed with code {result.returncode}")
        return False
    return True


def inspect_pypeit_file(path: Path) -> dict:
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
        description="Reduce all valid NGPS PypeIt configurations for all four channels."
    )

    parser.add_argument("date", help="UT date, e.g. 20260623")
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

    work_root = Path.home() / "ngps_data" / "work" / args.date
    raw_dir = work_root / "raw"

    if not raw_dir.exists():
        print(f"ERROR: raw directory does not exist:\n{raw_dir}")
        return 1

    print(f"\nNGPS night: {args.date}")
    print(f"Raw data:   {raw_dir}")

    jobs: list[tuple[str, Path, dict]] = []

    for channel in CHANNELS:
        setup_root = work_root / f"manual_setup_{channel}"

        if args.force_setup and setup_root.exists():
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

        pypeit_files = sorted(setup_root.glob("*/*.pypeit"))

        print(
            f"\n{'=' * 70}\n"
            f"CHANNEL {channel.upper()}: {len(pypeit_files)} configuration(s)\n"
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

            if info["has_science"] and info["has_arc"] and info["has_flat"]:
                print("  --> VALID: will reduce")
                jobs.append((channel, pf, info))
            else:
                print("  --> SKIP: incomplete science/calibration setup")

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
                f"Found {len(existing)} existing spec1d files. "
                "Skipping this setup."
            )
            continue

        cmd = ["run_pypeit", pf.name]
        if args.overwrite:
            cmd.append("--overwrite")

        ok = run(cmd, cwd=pf.parent)
        if not ok:
            print(f"\nWARNING: reduction failed for {pf.parent.name}.")

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

-------------------------------------------------------------------------------
END SCRIPT
-------------------------------------------------------------------------------

Save in nano:

    Ctrl+O
    Enter
    Ctrl+X

Run:

    conda activate ngps
    python "/scripts/ngps_reduce_all_configs.py" 20260623

Optional: regenerate setup directories:

    python "/scripts/ngps_reduce_all_configs.py" 20260623 --force-setup

Optional: rerun already reduced configurations:

    python "/scripts/ngps_reduce_all_configs.py" 20260623 --overwrite

===============================================================================
15. SUCCESSFUL RESULT FOR 20260623
===============================================================================

The successful reduction produced:

    R   p200_ngps_r_B   binning=3,2   spec1d=10
    R   p200_ngps_r_C   binning=3,4   spec1d=4
    R   p200_ngps_r_D   binning=4,4   spec1d=15

    G   p200_ngps_g_B   binning=3,2   spec1d=10
    G   p200_ngps_g_C   binning=3,4   spec1d=4
    G   p200_ngps_g_D   binning=4,4   spec1d=28

    I   p200_ngps_i_B   binning=3,2   spec1d=10
    I   p200_ngps_i_C   binning=3,4   spec1d=4
    I   p200_ngps_i_D   binning=4,4   spec1d=28

    U   p200_ngps_u_B   binning=3,2   spec1d=10
    U   p200_ngps_u_C   binning=3,4   spec1d=4
    U   p200_ngps_u_D   binning=4,4   spec1d=28

Total:

    155 spec1d files

At this point you have individual reduced 1D spectra, not yet coadded and not yet
fully flux calibrated.

===============================================================================
16. FIND THE REDUCED 1D SPECTRA
===============================================================================

Example for R:

    find ~/ngps_data/work/20260623/manual_setup_r -path "*/Science/spec1d*.fits"

Across all channels:

    find ~/ngps_data/work/20260623/manual_setup_{r,g,i,u} -path "*/Science/spec1d*.fits"

Find one target:

    find ~/ngps_data/work/20260623/manual_setup_{r,g,i,u} \
        -path "*Science/spec1d*MGC+04-48-002*.fits"

===============================================================================
17. VIEW A 1D SPECTRUM
===============================================================================

Example:

    cd ~/ngps_data/work/20260623/manual_setup_r/p200_ngps_r_B

    pypeit_show_1dspec \
        ./Science/spec1d_ngps_260623_0123-MGC+04-48-002_NGPS_r_20260623T094242.770.fits

If Ginga times out on first launch, start it manually:

    ginga --rcport=11771 --modules=RC,SlitWavelength

Leave Ginga running and rerun pypeit_show_1dspec in another terminal.

===============================================================================
18. WHAT IS COMPLETE AT THIS POINT
===============================================================================

For every valid B/C/D configuration in every U/G/R/I channel:

    [DONE] raw-file organization
    [DONE] frame typing
    [DONE] separation by detector-binning configuration
    [DONE] bias processing
    [DONE] flat processing
    [DONE] trace/slit calibration
    [DONE] wavelength calibration
    [DONE] sky subtraction
    [DONE] reduced 2D spectra
    [DONE] extracted 1D spectra
    [DONE] QA products

===============================================================================
19. WHAT IS STILL MISSING
===============================================================================

1. SENSITIVITY FUNCTION / FLUX CALIBRATION
-----------------------------------------
For each channel and each configuration, identify a suitable reduced standard-star
spec1d file, build a sensitivity function with pypeit_sensfunc, and apply it to
science spec1d files with pypeit_flux_calib.

This still needs a multi-configuration-aware automation script.

2. COADDITION OF REPEATED EXPOSURES
-----------------------------------
Repeated observations of the same target must be grouped by target, channel, and
configuration/binning, then coadded into one spectrum per target per channel.

3. TELLURIC CORRECTION
----------------------
Atmospheric absorption, especially in the redder channels, still needs to be
corrected. The NGPS wrapper has an empirical telluric procedure, but its behavior
must be adapted or checked for multiple configurations.

4. MERGE U + G + R + I
----------------------
After each target has one calibrated/coadded spectrum per channel, merge the four
channels with overlap checks, edge masking, flux-scale consistency checks, and
appropriate inverse-variance weighting.

5. SCIENCE VALIDATION
---------------------
Inspect wavelength accuracy, sky residuals, extraction quality, standard-star
response, flux calibration, channel overlaps, telluric residuals, propagated
uncertainties, and masks before scientific use.

===============================================================================
20. SPECIAL NOTE ABOUT SETUP A
===============================================================================

For 20260623:

    Setup A = binning 1,2

It contained science frames but no calibration frames that PypeIt considered
compatible with that setup, so the automation intentionally skipped it.

This does not prove that no Setup A calibration exposures were ever taken. It only
means that the available files and current PypeIt configuration matching did not
associate a usable calibration set with the 1,2 setup.

===============================================================================
21. CURRENT END POINT
===============================================================================

    151 raw FITS files
        ↓
    12 valid channel/configuration reductions
        ↓
    155 extracted spec1d files
        ↓
    NEXT:
        sensitivity functions
        flux calibration
        coaddition
        telluric correction
        U/G/R/I merging

===============================================================================

===============================================================================
22. FLUX CALIBRATION AND AUDIT — MAINTAINED PROJECT SCRIPTS
===============================================================================

The original guide stopped before flux calibration. The pinned workflow now
tracks the required scripts. After the extracted spectra are present, run:

    python "$WORKFLOW_ROOT/scripts/ngps_inventory_standards.py" 20260623
    python "$WORKFLOW_ROOT/scripts/ngps_flux_calibrate.py" 20260623
    python "$WORKFLOW_ROOT/scripts/ngps_flux_calibrate.py" 20260623 --run
    python "$WORKFLOW_ROOT/scripts/ngps_audit_flux.py" 20260623

The first flux-calibration invocation identifies the proposed associations and
writes the PypeIt configuration files; only `--run` creates sensitivity
functions and applies flux calibration. The scripts preserve unfluxed inputs by
working in Fluxed/ directories. For the known I/C sensitivity-function QA-only
failure, the tracked script retries without QA so that the valid
sensitivity-function FITS product is saved.

===============================================================================
23. INTERACTIVE REVIEW BEFORE COADDING
===============================================================================

Do not confuse repeat exposures with the three image-slicer traces that appear
inside one raw NGPS exposure. The coadd reviewer finds repeat observations of
one target by its inventory name, then groups them by channel and PypeIt setup.
For example, 0121, 0122, and 0123 are three repeat exposures; each panel in the
review window shows the three slicer traces belonging to one of those exposures.
For a good point-source exposure, keep all three slicer traces together: they
are the pieces that NGPS intends to recombine for the full source signal. A
spatially extended or blended source must first be checked in the 2D frame and,
where needed, re-extracted with the interactive extraction tool.

First list target/channel/setup groups without opening a review window:

    python "$WORKFLOW_ROOT/scripts/ngps_interactive_coadd.py" 20260623 \
        --list-groups

Then inspect the groups and proposed inputs for one target:

    python "$WORKFLOW_ROOT/scripts/ngps_interactive_coadd.py" 20260623 \
        --target MGC+04-48-002 \
        --summary

Then run the interactive review:

    python "$WORKFLOW_ROOT/scripts/ngps_interactive_coadd.py" 20260623 \
        --target MGC+04-48-002

If the target has several channel/setup groups, enter the group numbers to
review. The display contains an overlay and one panel plus checkbox per repeat
exposure. Unticking an exposure excludes all three of its slicer traces. To
choose only particular files before the window opens, use one or more
`--exposure` arguments. After acceptance, the script asks separately whether to
write a new coadd setup and whether to run it. It never alters the individual
Fluxed spectra or overwrites a previous coadd directory.

Before this coadd review, check the extraction-review PDF made for every raw
science exposure. To reduce without stopping, while still creating those PDFs,
use:

    python "$WORKFLOW_ROOT/scripts/ngps_reduce_all_configs.py" 20260623 --auto

The PDFs are in `ExtractionQA/<target>/`. Each one presents one exposure as a
four-channel U/G/R/I dashboard: four slicer-aligned 2D diagnostic panels,
coloured spatial profiles, and quick-look 1D spectra. The aligned 2D view is a
review display only, not a science coadd.

To revise a suspicious target, open the same per-exposure dashboards:

    python "$WORKFLOW_ROOT/scripts/ngps_manual_target_extractions.py" 20260623 \
        --target MGC+04-48-002

**Accept automatic** retains PypeIt's result. **Manual extraction** lets you
click a spatial position in any channel panel; that position is mapped to all
three slicers of that channel and linked across U/G/R/I. **Adjust this channel
only** makes the next click move only its selected channel. The
accepted manual review replaces the automatic PDF as the audit
record. Its PypeIt reduction is written only to a copied manual setup, leaving
the automatic detector products unchanged. After a manual setup is reduced,
rerun the inventory and flux-calibration steps before coadding.

Omit `--auto` when running `ngps_reduce_all_configs.py` to open each
per-exposure dashboard as soon as reduction is complete.

The generated setup, selected-file record, and coadd product are kept in:

    $NGPS_WORK_ROOT/20260623/Coadds/<target>_<channel>_<setup>/

Telluric correction and U/G/R/I merging should be reviewed in the same way;
they are not performed by this coadd command. The intended order is: coadd
repeat observations within each channel/setup, correct telluric absorption on
the resulting per-channel spectra, then merge U+G+R+I.

===============================================================================
24. RECORD THE WORKFLOW VERSION
===============================================================================

For each reduction, record the Git commit or release tag of this workflow
repository plus the contents of upstream-lock.yml. Before a new night, run:

    cd "$WORKFLOW_ROOT"
    python tools/verify_environment.py

If this reports a different or dirty upstream checkout, stop and recreate the
pinned environment before continuing. See docs/MAINTENANCE.md before any update.
