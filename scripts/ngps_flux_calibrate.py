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


ASSOCIATION_FIELDS = (
    "channel", "setup", "group_id", "target", "science_filenames",
    "group_mid_mjd", "standard_filename", "standard_target",
    "standard_mjd", "delta_minutes", "assignment_status",
)
REQUIRED_ASSOCIATION_FIELDS = tuple(
    field for field in ASSOCIATION_FIELDS if field != "assignment_status"
)


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
        print(f"Sensitivity function kept: {sensfile.name}")
        attempted[sensfile] = True
        return True
    print(f"Building sensitivity function: {sensfile.name}")
    success = run_command(
        pypeit_command("pypeit_sensfunc", str(standard["spec1d"]), "-o", str(sensfile)),
        Path(plan["setup_dir"]),
    )
    if not success:
        print(f"WARNING: sensitivity function failed for {standard['target']}")
    attempted[sensfile] = success
    return success


def build_selected_sensfuncs(
    plans: list[dict[str, object]], force: bool,
) -> dict[Path, bool]:
    """Build the standards currently selected in the reviewed association plan."""
    attempted: dict[Path, bool] = {}
    for plan in plans:
        for group in plan["groups"]:
            ensure_sensfunc(plan, group["standard"], force, attempted)
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


def flux_calibrate(plans: list[dict[str, object]], date: str) -> int:
    """Copy selected spectra, write flux files, and run PypeIt calibration."""
    total = 0
    for plan in plans:
        rows = []
        for group in plan["groups"]:
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

    attempted = build_selected_sensfuncs(plans, args.force_sensfunc)
    fallbacks, unresolved = propose_fallbacks(plans, args.force_sensfunc, attempted)
    if fallbacks:
        write_associations(associations_file, proposal_rows(plans))
        print("\nAutomatic fallback proposal written")
        print(f"Review: {associations_file}")
        print("No spectra were flux-calibrated. Run the dry run, then run --run again.")
        return 0
    if unresolved:
        print("\nFlux calibration stopped. No usable fallback standard was found for:")
        for identifier in unresolved:
            print(f"  {identifier}")
        return 1

    sensfuncs = sum(attempted.values())
    fluxed = flux_calibrate(plans, args.date)
    if fluxed != assigned:
        print("\nFlux calibration incomplete")
        print(f"Science spectra assigned: {assigned}")
        print(f"Science spectra flux-calibrated: {fluxed}")
        return 1
    print("\nFlux calibration complete")
    print(f"Sensitivity functions available: {sensfuncs}")
    print(f"Science spectra assigned: {assigned}")
    print(f"Science spectra flux-calibrated: {fluxed}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
