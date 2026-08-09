# NGPS reduction guide — pinned workflow

This guide intentionally refers to this repository, not floating upstream
branches or an editable `~/Software/PypeIt` checkout.

## Start a night

Activate the project environment and choose where observational data live:

```bash
conda activate ngps-pinned
export NGPS_WORK_ROOT="/Users/xpecax/ngps_data/work"
python tools/verify_environment.py
```

The verifier must report the commits in `upstream-lock.yml`. If it does not,
stop: results would not be tied to this documented workflow.

## Reduce and assess

The scripts are project-owned copies of the working 2026-06-23 workflow:

```bash
python scripts/ngps_reduce_all_configs.py 20260623
python scripts/ngps_inventory_standards.py 20260623
python scripts/ngps_flux_calibrate.py 20260623        # prints plan
python scripts/ngps_flux_calibrate.py 20260623 --run  # performs calibration
python scripts/ngps_audit_flux.py 20260623
```

Use `--force-setup` and `--overwrite` only after reviewing existing
products. The default behavior avoids destructive replacement.

## Manual extraction

Inspect a reduced `spec2d` product first. If manual extraction is warranted,
make a separate setup copy; never change the automatic setup. Use the project's
manual-selection procedure to select at most three traces (single, dual, or
rare triplet), show the resulting manual string, and ask before writing the
copied `.pypeit` file. Record the project tag and selected positions with the
resulting product.

## Reproduce a published reduction

Checkout the project tag named in the reduction record, create the matching
environment from `environment.yml`, run the verifier, then use the same data
root and commands above. See [maintenance](MAINTENANCE.md) for upgrades.
