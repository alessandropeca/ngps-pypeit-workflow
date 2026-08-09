# Maintenance and upstream updates

An upstream update is a small reproducibility project, not a routine pull.

## Before an update

1. Make sure the baseline is recorded: `python tools/verify_environment.py`.
2. Preserve one representative, non-proprietary smoke-test night (or its
   expected metadata/product manifest) outside Git if it is too large to store.
3. Create a branch, e.g. `upgrade/pypeit-YYYY-MM-DD`. Never update the
   baseline branch in place.

## Update procedure

1. Fetch the two upstream repositories and select exact candidate *commits*;
   do not use branch names.
2. Change only `upstream-lock.yml` and the matching VCS pins in
   `environment.yml`.
3. Build a new named environment; retain the previous environment unchanged.
4. Run `python tools/check_lock.py` and `python tools/verify_environment.py`.
5. On the smoke-test night, run the reducer, inventory, flux-calibration dry
   run, and audit. Compare:
   - generated PypeIt setup count and frame classifications;
   - science/standard inventory associations;
   - number of spec1d/spec2d products;
   - fluxed files with `OPT_FLAM` or `BOX_FLAM`;
   - manual-extraction copy behavior.
6. Review any difference in a pull request. Update this guide's validated
   baseline and add a dated changelog entry only after scientific review.
7. Tag the accepted repository commit (for example,
   `ngps-workflow-2026.08.0`). Put that tag—not an upstream branch—in papers,
   notebooks, and reduction logs.

## Rollback

Use the previous project tag and recreate that exact environment. Because the
old lock contains full object IDs, rollback does not require guessing which
upstream branch used to work.

## Responsibilities

- Upstream fixes belong upstream whenever practical; file an issue/PR with a
  minimal reproducer.
- A local workaround belongs in a clearly documented wrapper, never as an
  unrecorded edit inside a PypeIt checkout.
- Do not commit raw FITS data, spectra, tokens, or personal configuration.

