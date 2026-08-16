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
| Which 2–4 repeat exposures of one source to coadd | The coadd window groups candidates by exposure: one checkbox for 0121, 0122, or 0123 keeps or rejects all of that exposure's slicer traces together. |
| How to extract a source within one exposure | The 2D review window shows the automatic PypeIt traces. Click **Accept automatic**, or click up to three replacement positions and select **Accept manual positions**. |
| Whether to combine the three image-slicer traces within one exposure | For an integrated point-source spectrum, normally keep all three good slicer traces. Keep them separate only for a deliberately spatially resolved analysis. |

The current coadd command proposes one PypeIt trace from each of the three
slicer slices in every repeat exposure. It shows all of those candidates and
requires a human decision before PypeIt combines the accepted traces. Normally,
accept all three slices for each good point-source exposure. A rejected slice
must have a documented reason, such as a bad extraction or contaminating source.

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

The display shows the sky-subtracted 2D spectrum and PypeIt’s automatic traces
in gold. Click **Accept automatic** if they are correct: no files are changed.
If not, click the trace you want (up to three manual positions), drag a marker
or adjust its FWHM if needed, then click **Accept manual positions**. The script
asks before creating a copied manual setup and asks again before reducing it; it
never edits the original automatic setup.

After running a copied manual setup, repeat steps 5–6 for that copied setup so
its new 1D products are flux-calibrated before they enter coaddition. Do not
coadd an unfluxed manual product.

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

### 7. Find repeat observations of the same target

First list all target/channel/setup groups. This writes nothing and does not
open a window:

```bash
python scripts/ngps_interactive_coadd.py "$DATE" --list-groups
```

Then ask which repeat observations it found for one target. This also writes
nothing and does not open a window:

```bash
python scripts/ngps_interactive_coadd.py "$DATE" \
  --target 'MGC+04-48-002' \
  --summary
```

### 8. Interactively accept or reject the exposures, then coadd

```bash
python scripts/ngps_interactive_coadd.py "$DATE" \
  --target 'MGC+04-48-002'
```

If the target occurs in more than one channel or setup, choose the numbered
groups to review. Each window has one panel and checkbox per **repeat exposure**
(for example 0121, 0122, and 0123). Each panel displays that exposure's three
slicer traces, which stay together: unticking 0122 rejects all three of its
traces. A bad individual trace should instead be corrected in the earlier 2D
extraction-review step. Select **Accept selection** to continue. The script
asks before it writes the PypeIt coadd input and asks separately before it
launches the coadd. Select **Cancel** to write nothing.

To preselect only named raw observations, supply `--exposure` more than once:

```bash
python scripts/ngps_interactive_coadd.py "$DATE" \
  --target 'MGC+04-48-002' \
  --channel r \
  --setup p200_ngps_r_B \
  --exposure 0121 \
  --exposure 0123
```

Accepted inputs, a JSON selection record, the PypeIt `.coadd1d` file, and the
coadded FITS product are kept together in a new directory:

```text
$NIGHT/Coadds/MGC_04-48-002_r_p200_ngps_r_B/
```

The tool always keeps U, G, R, and I separate. Coadd repeat observations within
each channel/setup first; telluric-correct those per-channel coadds next; only
then inspect and merge the U+G+R+I products. Do not combine distinct setups
merely because they have the same target name.

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
