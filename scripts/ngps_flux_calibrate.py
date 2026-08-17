#!/usr/bin/env python3
"""Plan and run grouped NGPS flux calibration."""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from astropy.io import fits


ASSOCIATION_FIELDS = (
    "channel", "setup", "group_id", "target", "science_filenames",
    "group_mid_mjd", "standard_filename", "standard_target",
    "standard_mjd", "delta_minutes", "assignment_status",
)
REQUIRED_ASSOCIATION_FIELDS = tuple(
    field for field in ASSOCIATION_FIELDS if field != "assignment_status"
)
SENSITIVITY_REVIEW_FIELDS = (
    "channel", "setup", "standard_filename", "standard_target", "status", "detail",
)
MAX_SENSITIVITY_DISAGREEMENT_MAG = 1.0
MAX_STANDARD_REFERENCE_RESIDUAL_DEX = np.log10(2.0)


def run_command(command: list[str], cwd: Path) -> bool:
    """Run one PypeIt command and report whether it succeeded."""
    print("\n>>> " + " ".join(str(item) for item in command), flush=True)
    try:
        result = subprocess.run(command, cwd=cwd)
    except OSError as error:
        print(f"ERROR: could not start command: {error}")
        return False
    if result.returncode == 0:
        return True
    print(f"ERROR: command failed with exit code {result.returncode}")
    return False


def pypeit_command(name: str, *arguments: str) -> list[str]:
    """Run a PypeIt entry point from the active Python environment."""
    executable = Path(sys.executable).with_name(name)
    return [str(executable) if executable.is_file() else name, *arguments]


def validate_sensfunc(path: Path) -> tuple[bool, str]:
    """Reject standard observations that do not reproduce their reference flux."""
    try:
        with fits.open(path, memmap=False) as hdul:
            if "SENS" not in hdul:
                return False, "missing SENS extension"
            data = hdul["SENS"].data
            names = data.dtype.names or ()
            required = {
                "SENS_ZEROPOINT_FIT", "SENS_ZEROPOINT_FIT_GPM",
                "SENS_FLUXED_STD_WAVE", "SENS_FLUXED_STD_FLAM",
                "SENS_FLUXED_STD_FLAM_IVAR", "SENS_FLUXED_STD_MASK",
                "SENS_STD_MODEL_FLAM",
            }
            missing = required - set(names)
            if missing:
                return False, f"missing columns: {','.join(sorted(missing))}"
            fit = np.asarray(data["SENS_ZEROPOINT_FIT"], dtype=float).ravel()
            fit_gpm = np.asarray(data["SENS_ZEROPOINT_FIT_GPM"], dtype=bool).ravel()
            standard_wave = np.asarray(data["SENS_FLUXED_STD_WAVE"], dtype=float).ravel()
            standard_flux = np.asarray(data["SENS_FLUXED_STD_FLAM"], dtype=float).ravel()
            standard_ivar = np.asarray(data["SENS_FLUXED_STD_FLAM_IVAR"], dtype=float).ravel()
            standard_mask = np.asarray(data["SENS_FLUXED_STD_MASK"], dtype=bool).ravel()
            standard_model = np.asarray(data["SENS_STD_MODEL_FLAM"], dtype=float).ravel()
    except Exception as error:
        return False, f"cannot read file: {error}"
    usable = fit_gpm & np.isfinite(fit)
    fit_fraction = usable.sum() / max(fit_gpm.sum(), 1)
    if fit_fraction < 0.95:
        return False, f"only {fit_fraction:.0%} of fitted sensitivity pixels are finite"
    if not np.any(np.isfinite(standard_ivar) & (standard_ivar > 0)):
        return False, "standard-star flux inverse variance is zero everywhere"
    standard_good = (
        standard_mask & np.isfinite(standard_wave) & np.isfinite(standard_flux)
        & np.isfinite(standard_model) & np.isfinite(standard_ivar)
        & (standard_flux > 0) & (standard_model > 0) & (standard_ivar > 0)
    )
    if standard_good.sum() < 20:
        return False, "too few valid reference-spectrum comparison pixels"
    residual = np.abs(np.log10(standard_flux[standard_good] / standard_model[standard_good]))
    residual_90 = float(np.percentile(residual, 90))
    if residual_90 > MAX_STANDARD_REFERENCE_RESIDUAL_DEX:
        return False, (
            "calibrated standard disagrees with its reference spectrum in 10% or more "
            f"of pixels (90th-percentile residual {residual_90:.2f} dex)"
        )
    return True, "valid"


def sensitivity_curve(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Read a sensitivity fit and its fitted-good-pixel mask."""
    with fits.open(path, memmap=False) as hdul:
        data = hdul["SENS"].data
        fit = np.asarray(data["SENS_ZEROPOINT_FIT"], dtype=float).ravel()
        good = np.asarray(data["SENS_ZEROPOINT_FIT_GPM"], dtype=bool).ravel()
    return fit, good & np.isfinite(fit)


def central_sensitivity_offset(first: Path, second: Path) -> float | None:
    """Return the central-wavelength zero-point difference between two standards."""
    fit_a, good_a = sensitivity_curve(first)
    fit_b, good_b = sensitivity_curve(second)
    if len(fit_a) != len(fit_b):
        return None
    interior = np.zeros(len(fit_a), dtype=bool)
    interior[len(fit_a) // 10: 9 * len(fit_a) // 10] = True
    good = interior & good_a & good_b
    if good.sum() < 20:
        return None
    return float(abs(np.median(fit_a[good] - fit_b[good])))


def sensitivity_review(
    plans: list[dict[str, object]], attempted: dict[Path, bool],
) -> tuple[list[dict[str, str]], set[tuple[str, str]]]:
    """Review all standard responses and block setups with strong disagreement."""
    rows: list[dict[str, str]] = []
    blocked: set[tuple[str, str]] = set()
    for plan in plans:
        channel = str(plan["channel"]).upper()
        setup = str(plan["setup"])
        setup_rows: list[dict[str, str]] = []
        valid: list[tuple[dict[str, object], Path]] = []
        for standard in plan["standards"]:
            path = Path(standard["sensfile"])
            status, detail = "valid", "finite fitted sensitivity response"
            if not attempted.get(path, False):
                _, reason = validate_sensfunc(path) if path.exists() else (False, "sensitivity file was not built")
                status, detail = "rejected", reason
            else:
                valid.append((standard, path))
            row = {
                "channel": channel,
                "setup": setup,
                "standard_filename": str(standard["filename"]),
                "standard_target": str(standard["target"]),
                "status": status,
                "detail": detail,
            }
            rows.append(row)
            setup_rows.append(row)

        offsets = [
            offset
            for index, (_, first) in enumerate(valid)
            for _, second in valid[:index]
            if (offset := central_sensitivity_offset(first, second)) is not None
        ]
        if offsets and max(offsets) > MAX_SENSITIVITY_DISAGREEMENT_MAG:
            maximum = max(offsets)
            detail = (
                f"standard responses disagree by up to {maximum:.2f} mag in the central wavelength range"
            )
            for row in setup_rows:
                if row["status"] == "valid":
                    row["status"] = "review required"
                    row["detail"] = detail
            blocked.add((str(plan["channel"]).lower(), setup))
    return rows, blocked


def write_sensitivity_review(path: Path, rows: list[dict[str, str]]) -> None:
    """Write the standard-response quality record for one observing night."""
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SENSITIVITY_REVIEW_FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def find_spec1d(science_dir: Path, raw_filename: str) -> Path | None:
    """Return the extracted 1D product corresponding to one raw file."""
    matches = sorted(science_dir.glob(f"spec1d_{Path(raw_filename).stem}-*.fits"))
    if not matches:
        return None
    if len(matches) > 1:
        print(f"WARNING: multiple spec1d files for {raw_filename}. Using {matches[0].name}")
    return matches[0]


def safe_name(text: str) -> str:
    """Convert a target name to a stable filesystem-safe label."""
    return "".join(char if char.isalnum() or char in "-_" else "_" for char in text.strip()).strip("_")


def group_id(target: str, science: list[dict[str, object]]) -> str:
    """Return a stable identifier for one consecutive science sequence."""
    first = Path(str(science[0]["filename"])).stem
    last = Path(str(science[-1]["filename"])).stem
    return f"{safe_name(target)}__{first}-{last}"


def consecutive_groups(science: list[dict[str, object]]) -> list[list[dict[str, object]]]:
    """Group consecutive science exposures with the same target name."""
    groups: list[list[dict[str, object]]] = []
    for row in sorted(science, key=lambda item: float(item["mjd"])):
        if groups and str(groups[-1][0]["target"]).casefold() == str(row["target"]).casefold():
            groups[-1].append(row)
        else:
            groups.append([row])
    return groups


def nearest_standard(standards: list[dict[str, object]], mjd: float) -> dict[str, object]:
    return min(standards, key=lambda standard: abs(float(standard["mjd"]) - mjd))


def association_key(channel: str, setup: str, identifier: str) -> tuple[str, str, str]:
    return channel.lower(), setup, identifier


def read_associations(path: Path) -> dict[tuple[str, str, str], dict[str, str]]:
    """Read a reviewed association file keyed by channel, setup, and group."""
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows or set(REQUIRED_ASSOCIATION_FIELDS) - set(rows[0]):
        raise RuntimeError(f"{path} is missing required association columns")
    result: dict[tuple[str, str, str], dict[str, str]] = {}
    for row in rows:
        row.setdefault("assignment_status", "")
        key = association_key(row["channel"], row["setup"], row["group_id"])
        if key in result:
            raise RuntimeError(f"duplicate association group: {row['group_id']}")
        result[key] = row
    return result


def write_associations(path: Path, associations: list[dict[str, str]]) -> None:
    """Write a proposed association table for human review."""
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=ASSOCIATION_FIELDS)
        writer.writeheader()
        writer.writerows(associations)


def collect_plans(root: Path, rows: list[dict[str, object]]) -> list[dict[str, object]]:
    """Collect usable science groups and standard-star candidates by setup."""
    grouped_rows: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        frame_type = str(row["frametype"]).lower()
        if "science" in frame_type or "standard" in frame_type:
            grouped_rows[(str(row["channel"]).lower(), str(row["setup"]))].append(row)

    plans: list[dict[str, object]] = []
    for (channel, setup), setup_rows in sorted(grouped_rows.items()):
        if "_manual_" in setup or "_auto_" in setup:
            continue
        standards = [row for row in setup_rows if "standard" in str(row["frametype"]).lower()]
        science = [row for row in setup_rows if "science" in str(row["frametype"]).lower()]
        if not standards or not science:
            print(f"Skipping {channel.upper()} {setup}: science={len(science)}, standards={len(standards)}")
            continue

        setup_dir = root / f"manual_setup_{channel}" / setup
        science_dir = setup_dir / "Science"
        if not science_dir.exists():
            print(f"WARNING: science directory missing: {science_dir}")
            continue

        sens_dir = setup_dir / "Sensfunc"
        candidates: list[dict[str, object]] = []
        for standard in standards:
            spec1d = find_spec1d(science_dir, str(standard["filename"]))
            if spec1d is None:
                print(f"WARNING: no spec1d for standard {standard['filename']} ({standard['target']})")
                continue
            candidate = dict(standard)
            candidate["spec1d"] = spec1d
            candidate["sensfile"] = sens_dir / (
                f"sens_{channel.upper()}_{setup}_{safe_name(str(standard['target']))}_{Path(str(standard['filename'])).stem}.fits"
            )
            candidates.append(candidate)
        if not candidates:
            print(f"WARNING: no usable standards for {channel.upper()} {setup}")
            continue

        usable_science: list[dict[str, object]] = []
        for item in science:
            spec1d = find_spec1d(science_dir, str(item["filename"]))
            if spec1d is None:
                print(f"WARNING: no spec1d for science {item['filename']}")
                continue
            item = dict(item)
            item["spec1d"] = spec1d
            usable_science.append(item)
        if not usable_science:
            continue

        groups = []
        for members in consecutive_groups(usable_science):
            target = str(members[0]["target"])
            midpoint = sum(float(member["mjd"]) for member in members) / len(members)
            default = nearest_standard(candidates, midpoint)
            groups.append({
                "id": group_id(target, members),
                "target": target,
                "members": members,
                "midpoint": midpoint,
                "default": default,
                "standard": default,
                "manual": False,
                "status": "automatic",
            })
        plans.append({
            "channel": channel,
            "setup": setup,
            "setup_dir": setup_dir,
            "science_dir": science_dir,
            "sens_dir": sens_dir,
            "fluxed_dir": setup_dir / "Fluxed",
            "flux_files_dir": setup_dir / "FluxFiles",
            "standards": candidates,
            "groups": groups,
        })
    return plans


def proposal_rows(plans: list[dict[str, object]]) -> list[dict[str, str]]:
    rows = []
    for plan in plans:
        for group in plan["groups"]:
            standard = group["standard"]
            midpoint = float(group["midpoint"])
            rows.append({
                "channel": str(plan["channel"]).upper(),
                "setup": str(plan["setup"]),
                "group_id": str(group["id"]),
                "target": str(group["target"]),
                "science_filenames": " ".join(str(member["filename"]) for member in group["members"]),
                "group_mid_mjd": f"{midpoint:.8f}",
                "standard_filename": str(standard["filename"]),
                "standard_target": str(standard["target"]),
                "standard_mjd": f"{float(standard['mjd']):.8f}",
                "delta_minutes": f"{abs(float(standard['mjd']) - midpoint) * 1440.0:.1f}",
                "assignment_status": str(group["status"]),
            })
    return rows


def apply_reviewed_associations(
    plans: list[dict[str, object]], reviewed: dict[tuple[str, str, str], dict[str, str]],
) -> None:
    """Apply reviewed group choices and reject stale or invalid choices."""
    expected = {
        association_key(str(plan["channel"]), str(plan["setup"]), str(group["id"]))
        for plan in plans for group in plan["groups"]
    }
    missing = expected - set(reviewed)
    extra = set(reviewed) - expected
    if missing or extra:
        raise RuntimeError("association file does not match this inventory. Run the dry run with --reset-associations")

    for plan in plans:
        choices = {str(standard["filename"]): standard for standard in plan["standards"]}
        for group in plan["groups"]:
            key = association_key(str(plan["channel"]), str(plan["setup"]), str(group["id"]))
            choice = reviewed[key]["standard_filename"]
            if choice not in choices:
                raise RuntimeError(
                    f"{group['id']}: {choice} is not a usable standard in "
                    f"{str(plan['channel']).upper()} {plan['setup']}"
                )
            group["standard"] = choices[choice]
            saved_status = reviewed[key].get("assignment_status", "")
            if saved_status == "automatic fallback":
                group["manual"] = False
                group["status"] = saved_status
            elif choice != str(group["default"]["filename"]):
                group["manual"] = True
                group["status"] = "manual"


def print_associations(plans: list[dict[str, object]]) -> int:
    """Print one concise line for every consecutive science group."""
    count = 0
    print("\nScience-group associations")
    for plan in plans:
        for group in plan["groups"]:
            members = group["members"]
            exposures = ",".join(Path(str(member["filename"])).stem[-4:] for member in members)
            standard = group["standard"]
            midpoint = float(group["midpoint"])
            minutes = abs(float(standard["mjd"]) - midpoint) * 1440.0
            source = str(group["status"])
            print(
                f"  {str(plan['channel']).upper()} {plan['setup']} {group['target']} "
                f"[{exposures}] -> {standard['target']} ({minutes:.1f} min, {source})"
            )
            count += len(members)
    return count


def ensure_sensfunc(
    plan: dict[str, object], standard: dict[str, object], force: bool,
    attempted: dict[Path, bool],
) -> bool:
    """Use or build one sensitivity function exactly once per command run."""
    sensfile = Path(standard["sensfile"])
    if sensfile in attempted:
        return attempted[sensfile]
    sensfile.parent.mkdir(parents=True, exist_ok=True)
    if sensfile.exists() and not force:
        valid, reason = validate_sensfunc(sensfile)
        if valid:
            print(f"Sensitivity function kept: {sensfile.name}")
        else:
            print(f"WARNING: rejecting {sensfile.name}: {reason}")
        attempted[sensfile] = valid
        return valid
    print(f"Building sensitivity function: {sensfile.name}")
    success = run_command(
        pypeit_command("pypeit_sensfunc", str(standard["spec1d"]), "-o", str(sensfile)),
        Path(plan["setup_dir"]),
    )
    if not success:
        print(f"WARNING: sensitivity function failed for {standard['target']}")
        attempted[sensfile] = False
        return False
    valid, reason = validate_sensfunc(sensfile)
    if not valid:
        print(f"WARNING: rejecting {sensfile.name}: {reason}")
    attempted[sensfile] = valid
    return valid


def build_candidate_sensfuncs(
    plans: list[dict[str, object]], force: bool,
) -> dict[Path, bool]:
    """Build and validate every available standard in each setup."""
    attempted: dict[Path, bool] = {}
    for plan in plans:
        for standard in plan["standards"]:
            ensure_sensfunc(plan, standard, force, attempted)
    return attempted


def propose_fallbacks(
    plans: list[dict[str, object]], force: bool, attempted: dict[Path, bool],
) -> tuple[int, list[str]]:
    """Replace failed choices with the nearest standard that builds successfully."""
    changed = 0
    unresolved: list[str] = []
    for plan in plans:
        for group in plan["groups"]:
            selected = group["standard"]
            if attempted.get(Path(selected["sensfile"]), Path(selected["sensfile"]).exists()):
                continue
            midpoint = float(group["midpoint"])
            candidates = sorted(
                plan["standards"],
                key=lambda standard: abs(float(standard["mjd"]) - midpoint),
            )
            replacement = None
            for candidate in candidates:
                if str(candidate["filename"]) == str(selected["filename"]):
                    continue
                if ensure_sensfunc(plan, candidate, force, attempted):
                    replacement = candidate
                    break
            if replacement is None:
                unresolved.append(str(group["id"]))
                continue
            group["standard"] = replacement
            group["manual"] = False
            group["status"] = "automatic fallback"
            changed += 1
    return changed, unresolved


def quarantine_fluxed_group(plan: dict[str, object], group: dict[str, object]) -> None:
    """Move stale fluxed copies aside when their standard is not safe to use."""
    quarantine = Path(plan["setup_dir"]) / "Fluxed_invalid_standard"
    for member in group["members"]:
        source = Path(plan["fluxed_dir"]) / Path(member["spec1d"]).name
        if not source.exists():
            continue
        destination = quarantine / source.name
        if destination.exists():
            print(f"WARNING: keeping existing invalid copy: {destination.name}")
            continue
        quarantine.mkdir(parents=True, exist_ok=True)
        shutil.move(source, destination)
        print(f"Moved invalid fluxed copy aside: {destination}")


def flux_calibrate(
    plans: list[dict[str, object]], date: str, skipped_groups: set[str],
) -> int:
    """Copy and calibrate safe groups, keeping unsafe products out of Fluxed/."""
    total = 0
    for plan in plans:
        rows = []
        for group in plan["groups"]:
            if str(group["id"]) in skipped_groups:
                quarantine_fluxed_group(plan, group)
                continue
            standard = group["standard"]
            sensfile = Path(standard["sensfile"])
            if not sensfile.exists():
                print(f"WARNING: skipping {group['id']}. Missing {sensfile.name}")
                continue
            for member in group["members"]:
                source = Path(member["spec1d"])
                destination = Path(plan["fluxed_dir"]) / source.name
                Path(plan["fluxed_dir"]).mkdir(parents=True, exist_ok=True)
                shutil.copy2(source, destination)
                rows.append((destination, sensfile))
        if not rows:
            continue

        flux_files_dir = Path(plan["flux_files_dir"])
        flux_files_dir.mkdir(parents=True, exist_ok=True)
        flux_file = flux_files_dir / f"{plan['channel']}_{plan['setup']}.flux"
        with flux_file.open("w") as handle:
            handle.write("# NGPS flux calibration file\n")
            handle.write(f"# Date: {date}\n")
            handle.write("flux read\n    filename | sensfile\n")
            for science_file, sensfile in rows:
                handle.write(f"    {science_file} | {sensfile}\n")
            handle.write("flux end\n")
        print(f"Flux-calibrating {len(rows)} science spectrum/s: {plan['channel'].upper()} {plan['setup']}")
        if run_command(pypeit_command("pypeit_flux_calib", str(flux_file)), Path(plan["setup_dir"])):
            total += len(rows)
    return total


def main() -> int:
    parser = argparse.ArgumentParser(description="Plan and run grouped NGPS flux calibration.")
    parser.add_argument("date", help="UT observing date in YYYYMMDD format")
    parser.add_argument("--run", action="store_true", help="Build sensitivity functions and flux-calibrate the reviewed plan")
    parser.add_argument("--force-sensfunc", action="store_true", help="Regenerate existing sensitivity functions")
    parser.add_argument("--reset-associations", action="store_true", help="Write a new automatic association proposal during a dry run")
    args = parser.parse_args()
    if args.run and args.reset_associations:
        parser.error("run --reset-associations first, review the file, then run --run")

    root = Path(os.environ.get("NGPS_WORK_ROOT", Path.home() / "ngps_data" / "work")) / args.date
    inventory = root / "science_standard_inventory.csv"
    associations_file = root / "science_standard_associations.csv"
    if not inventory.exists():
        print(f"ERROR: inventory file does not exist: {inventory}")
        print("Run ngps_inventory_standards.py first.")
        return 1
    with inventory.open(newline="") as handle:
        rows: list[dict[str, object]] = list(csv.DictReader(handle))
    for row in rows:
        row["mjd"] = float(str(row["mjd"]))

    plans = collect_plans(root, rows)
    if not plans:
        print("ERROR: no usable science and standard-star groups found.")
        return 1

    if associations_file.exists() and not args.reset_associations:
        try:
            apply_reviewed_associations(plans, read_associations(associations_file))
        except RuntimeError as error:
            print(f"ERROR: {error}")
            return 1
        print(f"Using reviewed associations: {associations_file}")
    else:
        if args.run:
            print(f"ERROR: association file does not exist: {associations_file}")
            print("Run the dry run first, review the file, then run --run.")
            return 1
        write_associations(associations_file, proposal_rows(plans))
        print(f"Created association proposal: {associations_file}")
        print("Each row is one consecutive science group. Edit standard_filename if needed.")

    assigned = print_associations(plans)
    if not args.run:
        print("\nDry run complete. No spectra were flux-calibrated.")
        print(f"Review or edit: {associations_file}")
        print(f"Then run: python scripts/ngps_flux_calibrate.py {args.date} --run")
        return 0

    attempted = build_candidate_sensfuncs(plans, args.force_sensfunc)
    fallbacks, unresolved = propose_fallbacks(plans, args.force_sensfunc, attempted)
    review_rows, blocked_setups = sensitivity_review(plans, attempted)
    sensitivity_review_file = root / "sensitivity_review.csv"
    write_sensitivity_review(sensitivity_review_file, review_rows)
    blocked_groups = [
        str(group["id"])
        for plan in plans
        if (str(plan["channel"]).lower(), str(plan["setup"])) in blocked_setups
        for group in plan["groups"]
    ]
    if fallbacks:
        write_associations(associations_file, proposal_rows(plans))
        print("\nAutomatic fallback associations recorded")
        print(f"Review: {associations_file}")
    skipped_groups = set(unresolved) | set(blocked_groups)
    if skipped_groups:
        print("\nNo safe standard is available. These groups will not be flux-calibrated:")
        for identifier in sorted(skipped_groups):
            print(f"  {identifier}")
        print(f"Review: {sensitivity_review_file}")

    sensfuncs = sum(attempted.values())
    fluxed = flux_calibrate(plans, args.date, skipped_groups)
    skipped_spectra = sum(
        len(group["members"])
        for plan in plans
        for group in plan["groups"]
        if str(group["id"]) in skipped_groups
    )
    expected_fluxed = assigned - skipped_spectra
    if fluxed != expected_fluxed:
        print("\nFlux calibration incomplete")
        print(f"Science spectra assigned: {assigned}")
        print(f"Science spectra flux-calibrated: {fluxed}")
        return 1
    print("\nFlux calibration complete for all groups with a safe standard")
    print(f"Sensitivity functions available: {sensfuncs}")
    print(f"Science spectra assigned: {assigned}")
    print(f"Science spectra flux-calibrated: {fluxed}")
    if skipped_groups:
        print(f"Science spectra skipped: {skipped_spectra}")
        print(f"Review: {sensitivity_review_file}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
