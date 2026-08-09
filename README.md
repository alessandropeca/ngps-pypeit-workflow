# Reproducible NGPS / PypeIt reductions

This repository is the stable, laboratory-owned wrapper for the NGPS reduction
workflow. It does **not** copy or modify PypeIt or `ngps_pipeline`. Instead it
installs the exact upstream commits recorded in `upstream-lock.yml`, keeps our
small orchestration scripts here, and records a testable update procedure.

The validated baseline is the environment used for the 2026-06-23 reduction:

| Component | Upstream | pinned commit |
| --- | --- | --- |
| PypeIt NGPS fork | `cfremling/PypeIt` | `e9ed85c1a237c49626227f4227e323fc390def4b` |
| NGPS operational wrapper | `cfremling/ngps_pipeline` | `55fa9491eb1683769006118c46b26963bbf33ea2` |

No command in this repository follows an upstream branch such as `main` or
`develop`. A later upstream change cannot affect an existing environment until
we deliberately change the lock file, create a review branch, and pass the
smoke tests.

## Install a pinned environment

Create an isolated environment (for example, `conda env create -f
environment.yml`) and activate it. Then verify that it really points at the
pinned checkouts:

```bash
python tools/verify_environment.py
```

The verifier reports installed module paths and, for Git checkouts, rejects an
incorrect commit or a dirty checkout. Keep raw data and reduction products
outside this repository; set `NGPS_WORK_ROOT` to their parent work directory.

## Workflow

1. Generate/reduce configurations with `scripts/ngps_reduce_all_configs.py`.
2. Inventory science and standards with `scripts/ngps_inventory_standards.py`.
3. Flux calibrate using `scripts/ngps_flux_calibrate.py`.
4. Inspect results using `scripts/ngps_audit_flux.py`.
5. Use `scripts/ngps_manual_extract.py` for an optional, isolated manual run.

Every wrapper defaults to dry-run or refuses to overwrite. Manual extraction
copies a setup, updates only the copy, and can select at most three traces.

See [the reduction guide](docs/NGPS_REDUCTION_GUIDE.md) and
[maintenance](docs/MAINTENANCE.md).

## Attribution and licensing

PypeIt and `ngps_pipeline` remain independent upstream projects. Their code is
not vendored here; installation retains their respective licenses. This
repository includes our wrapper code only. Upstream URLs, commits, licenses, and
citation guidance are in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

