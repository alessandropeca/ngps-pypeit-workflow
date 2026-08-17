# NGPS / PypeIt reduction workflow

This repository contains the reproducible NGPS reduction scripts. It installs
the pinned NGPS-enabled [PypeIt fork](https://github.com/cfremling/PypeIt)
separately. 
PypeIt is the reduction engine and
[`ngps_pipeline`](https://github.com/alessandropeca/ngps_pipeline) is its
operational wrapper. The official PypeIt project is
[here](https://github.com/pypeit/PypeIt).

## Install once

Choose folders for your own computer before running these commands. The two
paths below are examples, not requirements; use any locations you prefer and
keep them consistent.

Note: The long strings such as `e9ed85c1a237c49626227f4227e323fc390def4b` are Git
commit IDs: unique labels for the exact tested software version. Do not edit
them during installation.

```bash
export GITHUB_ROOT="$HOME/Documents/GitHub"
export SOFTWARE_ROOT="$HOME/Software"
export WORKFLOW_ROOT="$GITHUB_ROOT/ngps-pypeit-workflow"
mkdir -p "$GITHUB_ROOT" "$SOFTWARE_ROOT"

git clone https://github.com/alessandropeca/ngps-pypeit-workflow.git "$WORKFLOW_ROOT"
cd "$WORKFLOW_ROOT"

conda env create -f environment.yml
conda activate ngps

git clone https://github.com/cfremling/PypeIt.git "$SOFTWARE_ROOT/PypeIt"
git -C "$SOFTWARE_ROOT/PypeIt" checkout e9ed85c1a237c49626227f4227e323fc390def4b
git clone https://github.com/alessandropeca/ngps_pipeline.git "$SOFTWARE_ROOT/ngps_pipeline"
git -C "$SOFTWARE_ROOT/ngps_pipeline" checkout 55fa9491eb1683769006118c46b26963bbf33ea2

python -m pip install -e "$SOFTWARE_ROOT/PypeIt"
python -m pip install -e "$SOFTWARE_ROOT/ngps_pipeline"
python tools/verify_environment.py
```

The final command must report both pinned commits as `OK`.

## Reduction for one night

Set the date and data location:

```bash
conda activate ngps  # if not already active
export WORKFLOW_ROOT="$HOME/Documents/GitHub/ngps-pypeit-workflow"  # your chosen location
cd "$WORKFLOW_ROOT"
export NGPS_WORK_ROOT="$HOME/ngps_data/work"
export DATE=20260623
export NIGHT="$NGPS_WORK_ROOT/$DATE"
```

1. Copy raw FITS files into `$NIGHT/raw/`; do not modify or split them.

   ```bash
   mkdir -p "$NIGHT/raw"
   rsync -av "/PATH/TO/RAW/FILES/"*.fits "$NIGHT/raw/"
   ```

2. Reduce every valid channel/setup and save an automatic extraction-review PDF
   for every science exposure in the entire night. This does not pause for decisions.
   `--auto` keeps existing extracted spectra and refreshes the review PDFs. It
   does not overwrite an extraction. If you omit `--auto`, a dashboard opens
   for each exposure so you can review and change the automatic settings.

   ```bash
   python scripts/ngps_reduce_all_configs.py "$DATE" --auto
   ```

   The PDFs to inspect are in `$NIGHT/ExtractionQA/<target>/` (for this night:
   `~/ngps_data/work/20260623/ExtractionQA/`). On macOS, open that folder with:

   ```bash
   open "$NIGHT/ExtractionQA"
   ```

   These PDFs are derived science-review products. Each one has four
   aligned U/G/R/I 2D panels, coloured spatial profiles, and a quick-look 1D
   panel.
   The aligned 2D panels are for checking the same source in the three slicers, they are not a science coadd.

3. To revise one already-reduced target whose automatic PDF does not look
   right, open the same dashboard interactively. This is a **single-target,
   single-exposure review command**.
   
   The **Manual extraction** button enables a click in a channel
   panel and links the position across U/G/R/I.

   The **Adjust this channel only** button makes the next click move only the panel selected. The saved PDF is
   replaced with your marked version.
   Manual positions inherit PypeIt's measured FWHM separately for every
   channel/slicer. Therefore, you can change pixel position but not extraction width.

   The **Accept Manual** and **Accept Automatic** buttons run an isolated
   one-exposure PypeIt setup. They replace only that exposure's derived
   `spec1d/spec2d` products (images and files) in the baseline setup.

   The **Return to automatic** button clears the manual aperture and restores the original
   profile and quick-look plots.

   The **Cancel** button, or closing the window, leaves the existing PDF and all PypeIt
   products unchanged.

   The **Re-norm U/G/R/I** buttons sets the quick-look y-range from one channel only; it is
   a display aid and does not change the detector counts or extracted spectrum.

   ```bash
   python scripts/ngps_manual_target_extractions.py "$DATE" \
     --target 'MGC+04-48-002' \
     --exposure 0121
   ```

   Alternatively, omit `--auto` during reduction to open the dashboard for each
   source as soon as the reductions finish:

   ```bash
   python scripts/ngps_reduce_all_configs.py "$DATE"
   ```

   If you accept automatic or manual re-extraction after flux calibration,
   repeat step 4 before coadding.

4. Flux-calibrate the 1D products. Only
   the third command changes science products.

   ```bash
   # 1. Read-only inventory: confirm the science and standard-star exposures.
   python scripts/ngps_inventory_standards.py "$DATE"

   # 2. Read-only dry run: print the proposed standard/science associations.
   python scripts/ngps_flux_calibrate.py "$DATE"

   # 3. Run: create sensitivity functions and Fluxed science copies.
   python scripts/ngps_flux_calibrate.py "$DATE" --run

   # 4. Read-only audit: verify that the Fluxed files contain FLAM.
   python scripts/ngps_audit_flux.py "$DATE"
   ```

5. Find repeated observations by target name, then review and coadd them.

   ```bash
   python scripts/ngps_interactive_coadd.py "$DATE" --list-groups
   python scripts/ngps_interactive_coadd.py "$DATE" \
     --target 'MGC+04-48-002'
   ```

   The window has one panel and checkbox per repeat exposure (e.g. 0121,
   0122, 0123). Each panel contains that exposure’s three NGPS slicer traces;
   accepting or rejecting an exposure keeps those three pieces together. The
   script asks before writing and again before coadding.

   To preselect observations:

   ```bash
   python scripts/ngps_interactive_coadd.py "$DATE" \
     --target 'MGC+04-48-002' --channel r --setup p200_ngps_r_B \
     --exposure 0121 --exposure 0123
   ```

## Additional notes

Three slicer traces belong to one raw exposure, while obs numbers such as 0121, 0122, and 0123, might be repeated
exposures of the same source (check that in observation log). 
Coadd repeats within each channel/setup first, then correct telluric
absorption, then merge U+G+R+I. Telluric correction and four-channel merging
are not automated yet, so the current final reproducible product is a reviewed,
flux-calibrated coadd per target/channel/setup.

For the fuller guide, including background and troubleshooting, see
[NGPS_REDUCTION_GUIDE.md](docs/NGPS_REDUCTION_GUIDE.md).
