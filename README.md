# Reproducible NGPS / PypeIt reductions

This repository contains our NGPS reduction workflow scripts. It does **not**
contain a modified copy of PypeIt: it installs the fixed NGPS-enabled PypeIt
fork [`cfremling/PypeIt`](https://github.com/cfremling/PypeIt) separately.

The main reduction engine is [PypeIt](https://github.com/pypeit/PypeIt), whose documentation is at [pypeit.readthedocs.io](https://pypeit.readthedocs.io/).

## What this workflow produces today

It takes unmodified raw NGPS FITS files through configuration-aware PypeIt
reduction, 2D/1D products, interactive manual-extraction review, flux
calibration, and a reviewed 1D coadd for one target, channel, and setup.

### Two different meanings of “three spectra”

For a source, your two, three, or four spectra are **repeat exposures**: separate
raw observations of the same target, which are the spectra you want to inspect
and potentially coadd.

Within *each one* of those raw NGPS exposures, the instrument's image slicer can
produce three PypeIt 1D **traces**. Those traces are spatial pieces of the same
exposure, not three repeat observations. They can look different, especially for
an extended source, an offset target, or an imperfect extraction.

| What you are deciding | Current workflow behavior |
| --- | --- |
| Which 2–4 repeat exposures of one source to coadd | The coadd window shows one candidate spectrum per exposure; you tick or untick each one. |
| How to extract a source within one exposure | Inspect the 2D image first; use the interactive extraction tool if the automatic trace is unsuitable. |
| Whether to combine the three image-slicer traces within one exposure | **Not automatic.** This is a separate spatial/science decision and must be defined from the 2D data before coadding repeat exposures. |

The current coadd command uses PypeIt's central-slicer trace as the initial
point-source candidate from each repeat exposure. It displays these candidate
exposures and requires a human decision before writing or running a coadd. It
does not silently combine the three slicer traces.

Telluric correction and the final U+G+R+I merge are not yet automated here.
Therefore, the current final reproducible product is one reviewed, fluxed,
coadded spectrum per target/channel/setup—not yet one merged four-channel
science spectrum. Do not describe the merged product as produced by this
repository until those two reviewed stages are implemented and validated.

## One-time installation

Run these steps once on a new computer. They create a dedicated environment and
install **two separate programs** from two clean, pinned upstream checkouts:

1. **PypeIt** — the spectroscopic reduction engine, including NGPS instrument
   support.
2. **`ngps_pipeline`** — the operational NGPS wrapper that uses PypeIt.

`ngps_pipeline` does not contain PypeIt. Both `python -m pip install -e`
commands below are required. Do not install either program from a branch name.

```bash
git clone https://github.com/alessandropeca/ngps-pypeit-workflow.git \
  ~/Documents/GitHub/ngps-pypeit-workflow
cd ~/Documents/GitHub/ngps-pypeit-workflow

conda env create -f environment.yml
conda activate ngps-pinned

mkdir -p ~/Software
git clone https://github.com/cfremling/PypeIt.git ~/Software/PypeIt
git -C ~/Software/PypeIt checkout e9ed85c1a237c49626227f4227e323fc390def4b

git clone https://github.com/alessandropeca/ngps_pipeline.git ~/Software/ngps_pipeline
git -C ~/Software/ngps_pipeline checkout 55fa9491eb1683769006118c46b26963bbf33ea2

python -m pip install -e ~/Software/PypeIt
python -m pip install -e ~/Software/ngps_pipeline

python tools/check_lock.py
python tools/verify_environment.py
```

`verify_environment.py` must report both commits as `OK`. If it reports a
different commit or a dirty checkout, stop and recreate the pinned checkouts;
do not continue with a mixture of versions.

## Reduction checklist: raw FITS to reviewed coadd

Use a new terminal for each night. Replace `20260623` below with the observing
date in `YYYYMMDD` form. The raw files and all products live outside this Git
repository, so they are never committed accidentally.

### 1. Activate the pinned environment and set paths

```bash
conda activate ngps-pinned
cd ~/Documents/GitHub/ngps-pypeit-workflow

export NGPS_WORK_ROOT="$HOME/ngps_data/work"
export DATE=20260623
export NIGHT="$NGPS_WORK_ROOT/$DATE"
```

For an already validated older environment named `ngps`, it is acceptable to
activate `ngps` instead, but first run `python tools/verify_environment.py`.

### 2. Stage raw data without altering it

```bash
mkdir -p "$NIGHT/raw"
rsync -av "/PATH/TO/RAW/FILES/"*.fits "$NIGHT/raw/"
find "$NIGHT/raw" -maxdepth 1 -name '*.fits' | wc -l
```

Keep the original FITS files unchanged. Do not split the NGPS U/G/R/I data
extensions manually.

### 3. Run the configuration-aware reduction

```bash
python scripts/ngps_reduce_all_configs.py "$DATE"
```

The script runs `pypeit_setup -c all` for U, G, R, and I, discovers the
detector-binning configurations, and reduces only those with science plus the
required arc and flat calibration frames. Existing completed setups are skipped.

Use `--overwrite` only when intentionally rerunning an already reduced setup.
Use `--force-setup` only when intentionally discarding and regenerating all
generated setup directories for that night.

The principal PypeIt products are here:

```text
$NIGHT/manual_setup_<channel>/<setup>/Science/spec2d_*.fits
$NIGHT/manual_setup_<channel>/<setup>/Science/spec1d_*.fits
```

### 4. Review or correct an extraction when needed

For a suspicious, extended, blended, or missed source, open the matching 2D
product. This uses Astropy and Matplotlib, not `pypeit_show_1dspec`.

```bash
python scripts/ngps_interactive_extract.py \
  "$NIGHT/manual_setup_r/p200_ngps_r_B/Science/spec2d_<EXPOSURE>.fits"
```

The display shows the sky-subtracted 2D spectrum, PypeIt’s suggested traces,
and any manual extraction positions. You may drag a manual marker, select up to
three positions, then press Enter. The script asks before creating a copied
manual setup and asks again before reducing it. It never edits the original
automatic setup.

### 5. Build the science/standard inventory

```bash
python scripts/ngps_inventory_standards.py "$DATE"
```

This writes the authoritative target names and setup names to:

```text
$NIGHT/science_standard_inventory.csv
```

Use its exact `target` and `setup` values in later coadd commands.

### 6. Plan flux calibration, then run it

```bash
python scripts/ngps_flux_calibrate.py "$DATE"
python scripts/ngps_flux_calibrate.py "$DATE" --run
python scripts/ngps_audit_flux.py "$DATE"
```

The first command is a review/planning pass: it identifies associations and
writes configuration files, but does not apply flux calibration. `--run` builds
the sensitivity functions and writes flux-calibrated copies, preserving the
unfluxed 1D inputs. The audit must show no missing reduced science or standard
files and report `FLAM` for the expected fluxed files.

Flux-calibrated inputs for coaddition are here:

```text
$NIGHT/manual_setup_<channel>/<setup>/Fluxed/spec1d_*.fits
```

### 7. Preview one proposed coadd

First ask the script which input spectra it proposes. This writes nothing and
does not open a window:

```bash
python scripts/ngps_interactive_coadd.py "$DATE" \
  --target 'MGC+04-48-002' \
  --channel r \
  --setup p200_ngps_r_B \
  --summary
```

### 8. Interactively accept or reject the exposures, then coadd

```bash
python scripts/ngps_interactive_coadd.py "$DATE" \
  --target 'MGC+04-48-002' \
  --channel r \
  --setup p200_ngps_r_B
```

The window contains an overlay plus one panel per exposure. Untick any exposure
that should not contribute, then select **Accept selection**. The script asks
before it writes the PypeIt coadd input and asks separately before it launches
the coadd. Select **Cancel** to write nothing.

Accepted inputs, a JSON selection record, the PypeIt `.coadd1d` file, and the
coadded FITS product are kept together in a new directory:

```text
$NIGHT/Coadds/MGC_04-48-002_r_p200_ngps_r_B/
```

Run steps 7–8 separately for every target, channel, and setup that you intend
to combine. Do not combine distinct setups or image-slicer traces merely because
they have the same target name.

### 9. Validate before scientific use

Inspect the 2D data, individual fluxed spectra, selected coadd inputs, masks,
wavelength coverage, sky/telluric residuals, and flux consistency. The
interactive review records what was accepted, but it does not replace scientific
quality control.

## Working with the 2026-06-23 reduction now

That night already has complete extracted and flux-calibrated products. Do
**not** restart raw-data staging, PypeIt reduction, or flux calibration. Start
at step 7 above for one source. For `MGC+04-48-002` in R/Setup B, run the
preview command, then the interactive command, and accept only the exposures
that look scientifically consistent.

## Safe maintenance

An upstream update is a small reproducibility project, not a routine `git
pull`. Follow [docs/MAINTENANCE.md](docs/MAINTENANCE.md): use a separate branch,
pin candidate commits, create a new environment, rerun the smoke test, review
the differences, and tag the accepted workflow revision. Never edit PypeIt or
`ngps_pipeline` inside a working reduction without recording that change.

PypeIt and `ngps_pipeline` remain independent upstream projects. Their source
repositories, licenses, and pinned revisions are recorded in
[`upstream-lock.yml`](upstream-lock.yml).
