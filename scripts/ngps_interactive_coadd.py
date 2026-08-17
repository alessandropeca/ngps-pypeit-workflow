#!/usr/bin/env python3
"""Interactively review and safely coadd one NGPS target/channel/setup.

NGPS divides each raw exposure into three image-slicer traces.  The reviewer
proposes one source trace from each slice of every repeat exposure and makes the
user accept or reject every proposed trace before PypeIt combines them.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import subprocess
from dataclasses import dataclass, asdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from matplotlib.widgets import Button, CheckButtons


COADD_REVIEW_FIELDS = (
    "target", "channel", "setup", "exposures", "status", "reason", "notes",
)


@dataclass
class Candidate:
    raw_filename: str
    spec1d: str
    obj_id: str
    slit_id: int


@dataclass
class ObservationGroup:
    """Repeat observations of one target in one channel and setup."""

    target: str
    channel: str
    setup: str
    rows: list[dict[str, str]]


def ask_yes_no(question: str, default: bool = False) -> bool:
    suffix = " [Y/n]: " if default else " [y/N]: "
    while True:
        answer = input(question + suffix).strip().lower()
        if not answer:
            return default
        if answer in {"y", "yes"}:
            return True
        if answer in {"n", "no"}:
            return False
        print("Please enter y or n.")


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_")


def find_spec1d(directory: Path, raw_filename: str) -> Path | None:
    matches = sorted(directory.glob(f"spec1d_{Path(raw_filename).stem}-*.fits"))
    return matches[0] if matches else None


def slit_and_spat(name: str) -> tuple[int, int] | None:
    match = re.fullmatch(r"SPAT(\d+)-SLIT(\d+)-DET01", name)
    return (int(match.group(2)), int(match.group(1))) if match else None


def object_metric(hdu) -> float:
    """A stable brightness proxy for choosing a trace within one slit."""
    names = hdu.data.dtype.names or ()
    flux_name = "OPT_FLAM" if "OPT_FLAM" in names else "OPT_COUNTS"
    ivar_name = "OPT_FLAM_IVAR" if "OPT_FLAM_IVAR" in names else "OPT_COUNTS_IVAR"
    flux = np.asarray(hdu.data[flux_name], dtype=float)
    ivar = np.asarray(hdu.data[ivar_name], dtype=float)
    with np.errstate(invalid="ignore"):
        return float(np.nansum(np.abs(flux) * np.sqrt(np.maximum(ivar, 0))))


def central_slit_from_standard(spec1d: Path) -> tuple[int, int]:
    """Use the brightest standard-star trace to identify NGPS's central slit."""
    best = (-np.inf, None, None)
    with fits.open(spec1d, memmap=False) as hdul:
        for hdu in hdul[1:]:
            parsed = slit_and_spat(hdu.name)
            if parsed is None:
                continue
            metric = object_metric(hdu)
            if metric > best[0]:
                slit, spat = parsed
                best = (metric, slit, spat)
    if best[1] is None:
        raise RuntimeError(f"No PypeIt object traces found in standard: {spec1d}")
    return best[1], best[2]


def best_trace_per_slicer(spec1d: Path) -> list[tuple[int, str]]:
    """Return the strongest unique PypeIt trace from each image-slicer slit.

    This creates review *candidates*, not an irreversible science decision. The
    user sees every candidate and can deselect a bad extraction before coadding.
    """
    candidates: dict[int, tuple[int, float, str]] = {}
    with fits.open(spec1d, memmap=False) as hdul:
        for hdu in hdul[1:]:
            parsed = slit_and_spat(hdu.name)
            if parsed is None:
                continue
            slit_id, _ = parsed
            metric = object_metric(hdu)
            # A manually clicked extraction is the explicit user choice.  If
            # it exists, prefer it to an automatic trace on the same slicer.
            manual = int(bool(hdu.header.get("HAND_EXTRACT_FLAG", False)))
            old = candidates.get(slit_id)
            if old is None or (manual, metric) > (old[0], old[1]):
                candidates[slit_id] = (manual, metric, hdu.name)
    return [(slit_id, item[2]) for slit_id, item in sorted(candidates.items())]


def plot_arrays(path: Path, obj_id: str) -> tuple[np.ndarray, np.ndarray]:
    with fits.open(path, memmap=False) as hdul:
        data = hdul[obj_id].data
        wave = np.asarray(data["OPT_WAVE"], dtype=float)
        flux = np.asarray(data["OPT_FLAM"], dtype=float)
        mask = np.asarray(data["OPT_MASK"], dtype=bool)
    good = np.isfinite(wave) & np.isfinite(flux) & mask
    return wave[good], flux[good]


def candidate_label(item: Candidate) -> str:
    exposure = Path(item.raw_filename).stem.rsplit("_", 1)[-1]
    return f"{exposure} | SLIT{item.slit_id:04d}"


def exposure_label(raw_filename: str) -> str:
    return Path(raw_filename).stem.rsplit("_", 1)[-1]


def requested_exposure(row: dict[str, str], requested: list[str]) -> bool:
    """Match an exact filename, its stem, or a four-digit frame number."""
    if not requested:
        return True
    filename = row["filename"]
    stem = Path(filename).stem
    return any(value in {filename, stem} or stem.endswith(f"_{value}") for value in requested)


def find_observation_groups(
    rows: list[dict[str, str]],
    target: str | None = None,
    channel: str | None = None,
    setup: str | None = None,
    exposures: list[str] | None = None,
) -> list[ObservationGroup]:
    """Find compatible repeat-observation groups from the inventory."""
    requested = exposures or []
    groups: dict[tuple[str, str, str], list[dict[str, str]]] = {}
    for row in rows:
        if "science" not in row["frametype"].lower():
            continue
        if target is not None and row["target"].casefold() != target.casefold():
            continue
        if channel is not None and row["channel"].lower() != channel:
            continue
        if setup is not None and row["setup"] != setup:
            continue
        if not requested_exposure(row, requested):
            continue
        key = (row["target"], row["channel"].lower(), row["setup"])
        groups.setdefault(key, []).append(row)
    return [
        ObservationGroup(name, channel_name, setup_name, sorted(items, key=lambda row: row["mjd"]))
        for (name, channel_name, setup_name), items in sorted(groups.items())
    ]


def review_key(group: ObservationGroup) -> tuple[str, str, str]:
    return group.target.casefold(), group.channel.lower(), group.setup


def review_exposures(group: ObservationGroup) -> str:
    return ",".join(exposure_label(row["filename"]) for row in group.rows)


def review_path(root: Path) -> Path:
    return root / "coadd_review.csv"


def group_review_row(root: Path, group: ObservationGroup) -> dict[str, str]:
    """Create the initial review decision for a possible coadd."""
    exposures = review_exposures(group)
    row = {
        "target": group.target,
        "channel": group.channel.upper(),
        "setup": group.setup,
        "exposures": exposures,
        "status": "review",
        "reason": "repeat science exposures",
        "notes": "",
    }
    if len(group.rows) < 2:
        row["status"] = "discard"
        row["reason"] = "only one science exposure"
        return row
    fluxed = root / f"manual_setup_{group.channel}" / group.setup / "Fluxed"
    missing = [
        exposure_label(item["filename"])
        for item in group.rows
        if find_spec1d(fluxed, item["filename"]) is None
    ]
    if missing:
        row["status"] = "discard"
        row["reason"] = f"missing Fluxed spectrum: {','.join(missing)}"
    return row


def update_coadd_review(root: Path, groups: list[ObservationGroup]) -> dict[tuple[str, str, str], dict[str, str]]:
    """Create or update the persistent, user-editable coadd-review table."""
    path = review_path(root)
    existing: dict[tuple[str, str, str], dict[str, str]] = {}
    if path.exists():
        with path.open(newline="") as handle:
            for row in csv.DictReader(handle):
                if set(COADD_REVIEW_FIELDS) - set(row):
                    raise RuntimeError(f"Invalid coadd review file: {path}")
                existing[(row["target"].casefold(), row["channel"].lower(), row["setup"])] = row

    updated: list[dict[str, str]] = []
    for group in groups:
        key = review_key(group)
        generated = group_review_row(root, group)
        saved = existing.get(key)
        if saved is not None and saved["exposures"] == generated["exposures"]:
            generated["status"] = saved["status"]
            generated["reason"] = saved["reason"]
            generated["notes"] = saved["notes"]
        updated.append(generated)

    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COADD_REVIEW_FIELDS)
        writer.writeheader()
        writer.writerows(updated)
    return {
        (row["target"].casefold(), row["channel"].lower(), row["setup"]): row
        for row in updated
    }


def reviewable_groups(
    groups: list[ObservationGroup], review: dict[tuple[str, str, str], dict[str, str]],
) -> list[ObservationGroup]:
    """Return only repeat groups that have not been discarded in the review file."""
    return [
        group for group in groups
        if len(group.rows) >= 2 and review[review_key(group)]["status"].casefold() == "review"
    ]


def print_groups(groups: list[ObservationGroup]) -> None:
    print("\nCompatible repeat-observation groups")
    print("#  target                    channel  setup                 exposures")
    for index, group in enumerate(groups, start=1):
        frames = ", ".join(exposure_label(row["filename"]) for row in group.rows)
        print(
            f"{index:<2} {group.target:<25} {group.channel.upper():<7} "
            f"{group.setup:<21} {frames}"
        )


def choose_groups(groups: list[ObservationGroup]) -> list[ObservationGroup]:
    """Let the user choose channel/setup groups; no files are written here."""
    if len(groups) == 1:
        return groups
    print_groups(groups)
    while True:
        answer = input("Groups to review [all, or comma-separated numbers]: ").strip().lower()
        if not answer or answer == "all":
            return groups
        try:
            selected = sorted({int(value.strip()) for value in answer.split(",")})
        except ValueError:
            selected = []
        if selected and all(1 <= index <= len(groups) for index in selected):
            return [groups[index - 1] for index in selected]
        print("Enter all, or valid group numbers separated by commas.")


def candidates_for_group(
    root: Path,
    rows: list[dict[str, str]],
    group: ObservationGroup,
) -> tuple[list[Candidate], int, int]:
    """Find all three slicer candidates for every repeat observation."""
    setup_rows = [
        row for row in rows
        if row["channel"].lower() == group.channel and row["setup"] == group.setup
    ]
    standard_rows = [row for row in setup_rows if "standard" in row["frametype"].lower()]
    if not standard_rows:
        raise RuntimeError("No standard-star row in this channel/setup")
    setup_dir = root / f"manual_setup_{group.channel}" / group.setup
    standard = find_spec1d(setup_dir / "Science", standard_rows[0]["filename"])
    if standard is None:
        raise RuntimeError("Could not find the reduced standard-star spec1d file")
    central_slit, standard_spat = central_slit_from_standard(standard)
    candidates: list[Candidate] = []
    for row in group.rows:
        spec1d = find_spec1d(setup_dir / "Fluxed", row["filename"])
        if spec1d is None:
            print(f"WARNING: no flux-calibrated spec1d for {row['filename']}")
            continue
        trace_ids = best_trace_per_slicer(spec1d)
        if not trace_ids:
            print(f"WARNING: no PypeIt slicer traces in {spec1d.name}")
            continue
        for slit_id, obj_id in trace_ids:
            candidates.append(Candidate(row["filename"], str(spec1d), obj_id, slit_id))
    return candidates, central_slit, standard_spat


def review(candidates: list[Candidate], target: str, channel: str, setup: str) -> list[Candidate] | None:
    """Review repeat exposures while keeping each exposure's slicer traces together."""
    if not candidates:
        return None
    by_exposure: dict[str, list[Candidate]] = {}
    for item in candidates:
        by_exposure.setdefault(item.raw_filename, []).append(item)
    exposures = list(by_exposure)
    figure, axes = plt.subplots(
        len(exposures) + 1,
        1,
        figsize=(14, 2.75 * (len(exposures) + 1)),
        sharex=True,
    )
    figure.subplots_adjust(left=0.26, bottom=0.06, top=0.93, hspace=0.38)
    axes = np.atleast_1d(axes)
    included = {exposure: True for exposure in exposures}
    plotted = {exposure: [] for exposure in exposures}

    overlay = axes[0]
    overlay.set_title(
        f"Proposed coadd: {target} | {channel.upper()} | {setup}\n"
        "One panel per repeat exposure; its three slicer traces stay together"
    )
    overlay.set_ylabel("Flux")
    line_styles = ("-", "--", ":")
    for index, exposure in enumerate(exposures):
        items = by_exposure[exposure]
        colour = f"C{index % 10}"
        panel = axes[index + 1]
        for trace_index, item in enumerate(items):
            wave, flux = plot_arrays(Path(item.spec1d), item.obj_id)
            label = exposure_label(exposure) if trace_index == 0 else "_nolegend_"
            line = overlay.plot(
                wave, flux, color=colour, ls=line_styles[trace_index % 3], lw=0.8, label=label
            )[0]
            plotted[exposure].append(line)
            line = panel.plot(
                wave,
                flux,
                color=colour,
                ls=line_styles[trace_index % 3],
                lw=0.75,
                label=f"SLIT{item.slit_id:04d}",
            )[0]
            plotted[exposure].append(line)
        panel.axhline(0, color="0.6", lw=0.6)
        panel.set_ylabel("Flux")
        panel.set_title(
            f"Exposure {exposure_label(exposure)}: {len(items)} slicer trace(s)",
            loc="left",
            fontsize=9,
        )
        panel.legend(fontsize=7, loc="upper right")
    overlay.legend(ncol=min(4, len(exposures)), fontsize=8, loc="upper right")
    axes[-1].set_xlabel("Vacuum wavelength (Å)")

    labels = [exposure_label(exposure) for exposure in exposures]
    label_to_exposure = dict(zip(labels, exposures))
    check_axis = figure.add_axes((0.02, 0.27, 0.22, 0.62))
    checks = CheckButtons(check_axis, labels, [included[exposure] for exposure in exposures])
    check_axis.set_title("Include exposure\n(all its slices)", fontsize=9)
    for label in checks.labels:
        label.set_fontsize(7)
    result = {"accepted": False}

    def toggle(label: str) -> None:
        exposure = label_to_exposure[label]
        included[exposure] = not included[exposure]
        for line in plotted[exposure]:
            line.set_alpha(1.0 if included[exposure] else 0.12)
        figure.canvas.draw_idle()

    def accept(event) -> None:
        result["accepted"] = True
        plt.close(figure)

    def cancel(event) -> None:
        plt.close(figure)

    checks.on_clicked(toggle)
    accept_axis = figure.add_axes((0.03, 0.18, 0.18, 0.05))
    cancel_axis = figure.add_axes((0.03, 0.11, 0.18, 0.05))
    Button(accept_axis, "Accept selection").on_clicked(accept)
    Button(cancel_axis, "Cancel").on_clicked(cancel)
    plt.show()
    if not result["accepted"]:
        return None
    return [item for item in candidates if included[item.raw_filename]]


def write_coadd_input(out_dir: Path, target: str, channel: str, setup: str,
                      central_slit: int, candidates: list[Candidate]) -> tuple[Path, Path]:
    """Write immutable review records and PypeIt's input file; never overwrite."""
    stem = f"{safe_name(target)}_{channel}_{safe_name(setup)}"
    out_dir.mkdir(parents=True, exist_ok=False)
    coadd_file = out_dir / f"{stem}.coadd1d"
    output = out_dir / f"{stem}_coadd.fits"
    coadd_file.write_text(
        "[coadd1d]\n"
        f"    coaddfile = '{output}'\n\n"
        "coadd1d read\n"
        "    filename | obj_id\n"
        + "".join(f"    {item.spec1d} | {item.obj_id}\n" for item in candidates)
        + "coadd1d end\n"
    )
    (out_dir / f"{stem}_selection.json").write_text(json.dumps({
        "target": target, "channel": channel, "setup": setup,
        "central_slit_anchor": central_slit,
        "candidates": [asdict(item) for item in candidates],
    }, indent=2) + "\n")
    return coadd_file, output


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Find, review, and optionally coadd repeat flux-calibrated NGPS observations."
    )
    parser.add_argument("date", help="UT date, e.g. 20260623")
    parser.add_argument("--target", help="Target name; matching is case-insensitive")
    parser.add_argument("--channel", choices=("u", "g", "r", "i"), help="Optional channel filter")
    parser.add_argument("--setup", help="Optional PypeIt setup filter, e.g. p200_ngps_r_B")
    parser.add_argument(
        "--exposure",
        action="append",
        default=[],
        help="Limit to a raw filename, stem, or four-digit frame number; repeat as needed",
    )
    parser.add_argument("--list-groups", action="store_true", help="List target/channel/setup groups and exit")
    parser.add_argument("--summary", action="store_true", help="Print candidates without opening the review window")
    args = parser.parse_args()

    root = Path(os.environ.get("NGPS_WORK_ROOT", Path.home() / "ngps_data" / "work")) / args.date
    inventory = root / "science_standard_inventory.csv"
    if not inventory.is_file():
        parser.error(f"Inventory not found: {inventory}; run ngps_inventory_standards.py first")
    with inventory.open() as handle:
        rows = list(csv.DictReader(handle))
    all_groups = find_observation_groups(rows)
    try:
        coadd_review = update_coadd_review(root, all_groups)
    except RuntimeError as error:
        parser.error(str(error))
    if args.list_groups:
        groups = find_observation_groups(rows, channel=args.channel, setup=args.setup)
        groups = reviewable_groups(groups, coadd_review)
        if not groups:
            parser.error("No reviewable repeat-observation groups found")
        print_groups(groups)
        discarded = sum(1 for row in coadd_review.values() if row["status"].casefold() == "discard")
        print(f"\nCoadd review file: {review_path(root)}")
        print(f"Automatically or manually discarded groups: {discarded}")
        return 0

    if args.target is None:
        parser.error("Specify --target, or use --list-groups to discover target names")
    groups = find_observation_groups(
        rows, args.target, args.channel, args.setup, args.exposure
    )
    groups = reviewable_groups(groups, coadd_review)
    if not groups:
        parser.error("No matching reviewable repeat observations found")
    if args.summary:
        print_groups(groups)
        for group in groups:
            candidates, central_slit, standard_spat = candidates_for_group(root, rows, group)
            print(
                f"\n{group.target} | {group.channel.upper()} | {group.setup}"
                f"\nCentral-slit anchor: SLIT{central_slit:04d} "
                f"(standard spatial pixel {standard_spat})"
            )
            for item in candidates:
                print(f"  {candidate_label(item):23s}  {item.obj_id}")
        return 0

    for group in choose_groups(groups):
        candidates, central_slit, standard_spat = candidates_for_group(root, rows, group)
        print(
            f"\nReviewing {group.target} | {group.channel.upper()} | {group.setup}"
            f"\nCentral-slit anchor: SLIT{central_slit:04d} "
            f"(standard spatial pixel {standard_spat})"
        )
        if not candidates:
            continue
        accepted = review(candidates, group.target, group.channel, group.setup)
        if not accepted:
            print("Coadd cancelled; no files were written for this group.")
            continue
        exposure_count = len({item.raw_filename for item in accepted})
        print(
            f"Accepted {exposure_count} of {len(group.rows)} repeat exposure(s), "
            f"containing {len(accepted)} slicer trace(s)."
        )
        if not ask_yes_no("Write this coadd selection and PypeIt input file?"):
            print("No files were written for this group.")
            continue
        out_dir = root / "Coadds" / (
            f"{safe_name(group.target)}_{group.channel}_{safe_name(group.setup)}"
        )
        coadd_file, output = write_coadd_input(
            out_dir, group.target, group.channel, group.setup, central_slit, accepted
        )
        print(f"Coadd input: {coadd_file}\nExpected output: {output}")
        if ask_yes_no("Run PypeIt coaddition now?"):
            status = subprocess.run([
                "pypeit_coadd_1dspec", str(coadd_file),
                "--par_outfile", str(out_dir / "coadd1d.par"),
            ]).returncode
            if status != 0:
                return status
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
