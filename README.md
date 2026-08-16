# NGPS / PypeIt reduction workflow

This repository contains the reproducible NGPS reduction scripts. It installs
the pinned NGPS-enabled [PypeIt fork](https://github.com/cfremling/PypeIt)
separately; PypeIt is the reduction engine and
[`ngps_pipeline`](https://github.com/alessandropeca/ngps_pipeline) is its
operational wrapper. The official PypeIt project is
[here](https://github.com/pypeit/PypeIt).

Pinned revisions are recorded in [`upstream-lock.yml`](upstream-lock.yml).
Do not update either upstream checkout without following
[maintenance](docs/MAINTENANCE.md).

The long strings such as `e9ed85c1a237c49626227f4227e323fc390def4b` are Git
commit IDs: unique labels for the exact tested software version. Do not edit
them during installation.

## Install once

Choose folders for your own computer before running these commands. The two
paths below are examples, not requirements; use any locations you prefer and
keep them consistent.

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

## Nightly reduction

Set the date and data location:

```bash
conda activate ngps
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

2. Reduce every valid channel/setup.

   ```bash
   python scripts/ngps_reduce_all_configs.py "$DATE"
   ```

3. Optionally review a suspicious 2D extraction. Gold lines are PypeIt’s
   automatic traces. Choose **Accept automatic**, or click up to three manual
   positions and choose **Accept manual positions**. Manual choices write only
   to a copied setup.

   ```bash
   python scripts/ngps_interactive_extract.py \
     "$NIGHT/manual_setup_r/p200_ngps_r_B/Science/spec2d_<EXPOSURE>.fits"
   ```

   If a manual setup is rerun, flux-calibrate its new products before coadding.

4. Inventory, flux-calibrate, and audit the 1D products.

   ```bash
   python scripts/ngps_inventory_standards.py "$DATE"
   python scripts/ngps_flux_calibrate.py "$DATE"
   python scripts/ngps_flux_calibrate.py "$DATE" --run
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

## Product order

Three slicer traces belong to one raw exposure; 0121, 0122, and 0123 are repeat
exposures. Coadd repeats within each channel/setup first, then correct telluric
absorption, then merge U+G+R+I. Telluric correction and four-channel merging
are not automated yet, so the current final reproducible product is a reviewed,
flux-calibrated coadd per target/channel/setup.

For the fuller guide, including background and troubleshooting, see
[NGPS_REDUCTION_GUIDE.md](docs/NGPS_REDUCTION_GUIDE.md).
