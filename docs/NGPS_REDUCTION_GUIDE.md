## Purpose

This guide reduces one night of Palomar/NGPS data from unmodified raw FITS
files to flux-calibrated, telluric-corrected (R and I), coadded 1D spectra.
It uses the maintained, version-pinned workflow. Do not replace the pinned
software commits with a newer upstream branch during a reduction.

This is a working guide, not a recipe for modifying raw data. Keep the raw
FITS files unchanged. All products are written under the chosen work folder.

## Before you start

You need:

- macOS or Linux, Git, and Conda (Miniconda or Anaconda)
- enough free disk space for raw data, calibrations, and reductions
- the original NGPS FITS files for one observing night
- internet access during the one-time installation

Choose your own folders. The examples below use:

GitHub projects:  $HOME/Documents/GitHub
Software:         $HOME/Software
Work products:    $HOME/ngps_data/work

Do not copy another person's absolute paths. Set the variables below once in
each new terminal session so the commands work on your computer.

# PART A — INSTALL ONCE

1. **Clone the maintained workflow**

---

Open Terminal and run:

```bash
export GITHUB_ROOT="$HOME/Documents/GitHub"
export SOFTWARE_ROOT="$HOME/Software"
export WORKFLOW_ROOT="$GITHUB_ROOT/ngps-pypeit-workflow"
mkdir -p "$GITHUB_ROOT" "$SOFTWARE_ROOT"
git clone https://github.com/alessandropeca/ngps-pypeit-workflow.git "$WORKFLOW_ROOT"
cd "$WORKFLOW_ROOT"
```

2. **Create the Conda environment**

---

```bash
conda env create -f environment.yml
conda activate ngps
```

3. **Install the pinned reduction software**

---

PypeIt performs the reduction. ngps_pipeline is the operational NGPS wrapper.
Both are required. The long strings below are fixed Git commit IDs. Leave them
exactly as written.

```bash
git clone https://github.com/cfremling/PypeIt.git "$SOFTWARE_ROOT/PypeIt"
git -C "$SOFTWARE_ROOT/PypeIt" checkout e9ed85c1a237c49626227f4227e323fc390def4b
git clone https://github.com/alessandropeca/ngps_pipeline.git "$SOFTWARE_ROOT/ngps_pipeline"
git -C "$SOFTWARE_ROOT/ngps_pipeline" checkout 55fa9491eb1683769006118c46b26963bbf33ea2
python -m pip install -e "$SOFTWARE_ROOT/PypeIt"
python -m pip install -e "$SOFTWARE_ROOT/ngps_pipeline"
python tools/apply_pypeit_manual_refit_patch.py "$SOFTWARE_ROOT/PypeIt"
python tools/verify_environment.py
```

The final command must report both pinned commits as OK. If it does not, stop
and ask the supervisor before reducing data.

# PART B — REDUCE ONE REAL NIGHT

The real example below uses the 2026-06-23 NGPS night and the target
MGC+04-48-002. Substitute your own date and target later, but keep the order
of operations unchanged.

1. **Start a new reduction session**

---

```bash
conda activate ngps
export WORKFLOW_ROOT="$HOME/Documents/GitHub/ngps-pypeit-workflow"
export NGPS_WORK_ROOT="$HOME/ngps_data/work"
export DATE=20260623
export NIGHT="$NGPS_WORK_ROOT/$DATE"
cd "$WORKFLOW_ROOT"
```

2. **Copy the raw FITS files**

---

Replace /PATH/TO/RAW/FILES with the folder containing the original FITS files.

```bash
mkdir -p "$NIGHT/raw"
rsync -av "/PATH/TO/RAW/FILES/"*.fits "$NIGHT/raw/"
find "$NIGHT/raw" -maxdepth 1 -name '*.fits' | wc -l
```

Never rename, edit, split, or overwrite files in $NIGHT/raw/.

3. **Reduce and review every science exposure**

---

```bash
python scripts/ngps_reduce_all_configs.py "$DATE"
```

This reduces every valid U/G/R/I channel and instrumental setup, then opens
one extraction-review window at a time. Re-running replaces existing
reduction products and refreshes the review PDFs.

Do not use `--auto` for the student or full science-review workflow. It skips
the windows and saves refreshed automatic-review PDFs only.

The review PDFs are here:

`$NIGHT/ExtractionQA/<target>/`

For the example target, open:

`open "$NIGHT/ExtractionQA/MGC_04-48-002"`

4. **Understand the extraction-review window**

---

Each PDF contains:

- four central-slicer 2D panels: U, G, R, I
- gold curves: PypeIt's automatic traces
- spatial profiles, one colour per channel
- quick-look 1D counts spectra

The three NGPS slicers are extracted separately by PypeIt. The central-slicer
panel is only the clearest place to choose the source position.

If the automatic trace follows the desired source, click **Accept automatic**.
If it does not, click **Manual extraction + refit**, select the desired trace
in the 2D panel, and then click **Accept manual**.

5. **Re-open one exposure for another review**

---

Example:

```bash
python scripts/ngps_manual_target_extractions.py "$DATE" --target 'MGC+04-48-002' --exposure 0121
```

Buttons in the window:

Accept automatic
Rerun and replace only this exposure using PypeIt's automatic choice.

Manual extraction + refit
Click the target in any channel. The same slicer-relative position is
applied to U/G/R/I. PypeIt refits the trace and FWHM independently in
every channel and in all three slicers.

Adjust this channel only
Click a channel to refit only that channel and its three slicers. The
remaining channels retain their previous extracted products.

Return to automatic
Remove manual choices and restore the automatic display.

Accept manual
Rerun and replace only the selected exposure/channel products.

Cancel, or close the window
Make no changes to products or the existing review PDF.

After accepting automatic or manual extraction, wait for the terminal to say
that the re-extraction has finished. If the exposure was already flux
calibrated, repeat Parts 9–12 below before using it in a coadd.

6. **Build and inspect the flux-calibration plan**

---

```bash
python scripts/ngps_inventory_standards.py "$DATE"
python scripts/ngps_flux_calibrate.py "$DATE"
```

These commands do not alter spectra. They create:

`$NIGHT/science_standard_inventory.csv`
`$NIGHT/science_standard_associations.csv`

Inspect science_standard_associations.csv. Each row assigns one standard star
to one consecutive science-exposure group, channel, and setup. If an
association is unsuitable, edit the standard_filename in that row, save it,
and rerun the second command to check the plan.

7. **Run and audit flux calibration**

---

```bash
python scripts/ngps_flux_calibrate.py "$DATE" --run
python scripts/ngps_audit_flux.py "$DATE"
```

Flux-calibrated copies are written in each setup's Fluxed/ folder. The audit
must be read before coadding. A group with no safe standard remains unfluxed
and is excluded from coaddition. Read:

`$NIGHT/sensitivity_review.csv`

8. **Identify repeat observations and review the proposed coadds**

---

```bash
python scripts/ngps_interactive_coadd.py "$DATE" --list-groups
```

This writes `$NIGHT/coadd_review.csv`. It groups repeated observations of the
same target by channel and setup. Review this file and set status to discard
for any exposure flagged in the observing log or visibly unsuitable.

9. **Coadd all safe repeat groups automatically**

---

```bash
python scripts/ngps_interactive_coadd.py "$DATE" --all --auto
```

This saves a review PDF for each coadd and writes the coadded FITS products.
It keeps all three slicer traces from an included exposure together. Outputs:

Review PDFs:  `$NIGHT/CoaddQA/<target>/`
Coadded FITS: `$NIGHT/Coadds/<target_channel_setup>/`
Summary:      `$NIGHT/coadd_run_summary.csv`

For a single target, open the interactive coadd review instead:

```bash
python scripts/ngps_interactive_coadd.py "$DATE" --target 'MGC+04-48-002'
```

Click Accept selection to replace that target/channel/setup coadd. Click
Cancel or close the window to leave it unchanged.

10. **Apply telluric correction to R and I coadds**

---

Install the PypeIt atmospheric model once:

```bash
pypeit_install_telluric TellPCA_3000_26000_R10000.fits
```

First inspect the proposed work:

```bash
python scripts/ngps_telluric_correct.py "$DATE" --all
```

Then run it:

```bash
python scripts/ngps_telluric_correct.py "$DATE" --all --run
```

Telluric correction is applied to R and I only. U and G have no telluric
correction in this workflow. Corrected products and QA PDFs are saved in:

`$NIGHT/Telluric/<target>/
$NIGHT/TelluricQA/<target>/`

Read `$NIGHT/telluric_review.csv`. Failed telluric products are not used for the
final plot.

11. **Make final U/G/R/I plots**

---

For every target:

```bash
python scripts/ngps_plot_final_spectra.py "$DATE" --all --noUGedges
```

For one target, including an interactive plot window:

```bash
python scripts/ngps_plot_final_spectra.py "$DATE" --target 'MGC+04-48-002'
```

Final PDF and PNG files are saved in:

`$NIGHT/FinalQA/<target>/`

The final plot uses telluric-corrected R/I products when available. It keeps
the four channels separate and does not merge them into a single spectrum.
The grey curve is native sampling; the coloured curve is a display-only
inverse-variance rebin of two pixels. Neither changes the FITS data.

12. **Deliverable files**

---

For scientific sharing, provide the final coadded FITS files:

U/G: `$NIGHT/Coadds/<target_channel_setup>/*_coadd.fits`
R/I: `$NIGHT/Telluric/<target>/*_coadd_tellcorr.fits`

Also provide the final QA plot and state clearly:

- observing date
- target name
- total exposure time
- channel/setup
- whether R/I are telluric-corrected
- whether a manual extraction was accepted

# COMMON RULES AND TROUBLESHOOTING

1. Do not use pypeit_show_1dspec in this workflow. Use the saved PDFs and the
final plotting script instead.
2. If an automatic extraction is wrong, revise only the affected exposure.
Do not rerun the full night unless you want to replace all automatic
products.
3. A manual position is not merely shifted. It triggers a new trace and FWHM
fit for the selected source. Linked mode refits all channels. Per-channel
mode refits only the chosen channel.
4. If you change an extraction after flux calibration, you must repeat flux
calibration, coaddition, telluric correction, and final plotting for the
affected target before distributing its final spectrum.
5. If the flux audit reports an unsafe or missing standard, do not force a
coadd. Read sensitivity_review.csv and ask the supervisor.
6. If a target has only one exposure, it is not a repeat-exposure coadd group.
Its flux-calibrated extracted spectrum remains a valid single-exposure
product, but it should be labelled accordingly.
7. Keep the pinned installation. To update software later, follow the
repository's docs/MAINTENANCE.md and validate on test data before using a
new version for science.
