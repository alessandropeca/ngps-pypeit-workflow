#!/usr/bin/env python3
"""Interactively review and safely coadd one NGPS target/channel/setup.

The reviewer intentionally starts with the standard-star-anchored *central*
image-slicer slit.  It never silently combines the three slicer regions: an
integrated-galaxy extraction must be designed from the 2D data first.
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


@dataclass
class Candidate:
    raw_filename: str
    spec1d: str
    obj_id: str


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


def central_trace(spec1d: Path, slit_id: int, reference_spat: int) -> str | None:
    """Pick the source closest to the standard trace on the specified slit."""
    candidates: list[tuple[int, str]] = []
    with fits.open(spec1d, memmap=False) as hdul:
        for hdu in hdul[1:]:
            parsed = slit_and_spat(hdu.name)
            if parsed is not None and parsed[0] == slit_id:
                candidates.append((parsed[1], hdu.name))
    if not candidates:
        return None
    return min(candidates, key=lambda item: abs(item[0] - reference_spat))[1]


def plot_arrays(path: Path, obj_id: str) -> tuple[np.ndarray, np.ndarray]:
    with fits.open(path, memmap=False) as hdul:
        data = hdul[obj_id].data
        wave = np.asarray(data["OPT_WAVE"], dtype=float)
        flux = np.asarray(data["OPT_FLAM"], dtype=float)
        mask = np.asarray(data["OPT_MASK"], dtype=bool)
    good = np.isfinite(wave) & np.isfinite(flux) & mask
    return wave[good], flux[good]


def review(candidates: list[Candidate], target: str, channel: str, setup: str) -> list[Candidate] | None:
    """Display every proposed spectrum and return the accepted subset."""
    if not candidates:
        return None
    figure, axes = plt.subplots(len(candidates) + 1, 1, figsize=(14, 2.5 * (len(candidates) + 1)), sharex=True)
    figure.subplots_adjust(left=0.25, bottom=0.08, top=0.93, hspace=0.35)
    axes = np.atleast_1d(axes)
    included = [True] * len(candidates)
    plotted: list[list] = [[] for _ in candidates]

    overlay = axes[0]
    overlay.set_title(f"Proposed coadd: {target} | {channel.upper()} | {setup}\nCentral-slit traces only")
    overlay.set_ylabel("Flux")
    for index, item in enumerate(candidates):
        wave, flux = plot_arrays(Path(item.spec1d), item.obj_id)
        label = Path(item.raw_filename).stem.rsplit("_", 1)[-1]
        line = overlay.plot(wave, flux, lw=0.8, label=label)[0]
        plotted[index].append(line)
        panel = axes[index + 1]
        line = panel.plot(wave, flux, color=line.get_color(), lw=0.75)[0]
        plotted[index].append(line)
        panel.axhline(0, color="0.6", lw=0.6)
        panel.set_ylabel("Flux")
        panel.set_title(f"{label}: {item.obj_id}", loc="left", fontsize=9)
    overlay.legend(ncol=min(4, len(candidates)), fontsize=8, loc="upper right")
    axes[-1].set_xlabel("Vacuum wavelength (Å)")

    labels = [Path(item.raw_filename).stem.rsplit("_", 1)[-1] for item in candidates]
    check_axis = figure.add_axes((0.02, 0.55, 0.18, 0.28))
    checks = CheckButtons(check_axis, labels, included)
    check_axis.set_title("Include exposure", fontsize=9)
    result = {"accepted": False}

    def toggle(label: str) -> None:
        index = labels.index(label)
        included[index] = not included[index]
        for line in plotted[index]:
            line.set_alpha(1.0 if included[index] else 0.12)
        figure.canvas.draw_idle()

    def accept(event) -> None:
        result["accepted"] = True
        plt.close(figure)

    def cancel(event) -> None:
        plt.close(figure)

    checks.on_clicked(toggle)
    accept_axis = figure.add_axes((0.03, 0.43, 0.16, 0.05))
    cancel_axis = figure.add_axes((0.03, 0.36, 0.16, 0.05))
    Button(accept_axis, "Accept selection").on_clicked(accept)
    Button(cancel_axis, "Cancel").on_clicked(cancel)
    plt.show()
    if not result["accepted"]:
        return None
    return [item for item, keep in zip(candidates, included) if keep]


def write_coadd_input(out_dir: Path, target: str, channel: str, setup: str,
                      slit_id: int, candidates: list[Candidate]) -> tuple[Path, Path]:
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
        "central_slit": slit_id, "candidates": [asdict(item) for item in candidates],
    }, indent=2) + "\n")
    return coadd_file, output


def main() -> int:
    parser = argparse.ArgumentParser(description="Review and optionally coadd one flux-calibrated NGPS target.")
    parser.add_argument("date", help="UT date, e.g. 20260623")
    parser.add_argument("--target", required=True, help="Target name exactly as shown by the inventory")
    parser.add_argument("--channel", required=True, choices=("u", "g", "r", "i"))
    parser.add_argument("--setup", required=True, help="PypeIt setup, e.g. p200_ngps_r_B")
    parser.add_argument("--summary", action="store_true", help="Print candidates without opening the review window")
    args = parser.parse_args()

    root = Path(os.environ.get("NGPS_WORK_ROOT", Path.home() / "ngps_data" / "work")) / args.date
    inventory = root / "science_standard_inventory.csv"
    if not inventory.is_file():
        parser.error(f"Inventory not found: {inventory}; run ngps_inventory_standards.py first")
    with inventory.open() as handle:
        rows = list(csv.DictReader(handle))
    setup_rows = [row for row in rows if row["channel"].lower() == args.channel and row["setup"] == args.setup]
    science_rows = [row for row in setup_rows if "science" in row["frametype"].lower() and row["target"] == args.target]
    standard_rows = [row for row in setup_rows if "standard" in row["frametype"].lower()]
    if not science_rows or not standard_rows:
        parser.error("No matching science rows or no standard-star row in this channel/setup")

    setup_dir = root / f"manual_setup_{args.channel}" / args.setup
    standard = find_spec1d(setup_dir / "Science", standard_rows[0]["filename"])
    if standard is None:
        parser.error("Could not find the reduced standard-star spec1d file")
    slit_id, standard_spat = central_slit_from_standard(standard)
    candidates: list[Candidate] = []
    for row in science_rows:
        spec1d = find_spec1d(setup_dir / "Fluxed", row["filename"])
        if spec1d is None:
            print(f"WARNING: no flux-calibrated spec1d for {row['filename']}")
            continue
        obj_id = central_trace(spec1d, slit_id, standard_spat)
        if obj_id is None:
            print(f"WARNING: no central-slit trace in {spec1d.name}")
            continue
        candidates.append(Candidate(row["filename"], str(spec1d), obj_id))
    print(f"Central slit: SLIT{slit_id:04d} (standard spatial pixel {standard_spat})")
    for item in candidates:
        print(f"  {item.raw_filename}  {item.obj_id}")
    if not candidates or args.summary:
        return 0

    accepted = review(candidates, args.target, args.channel, args.setup)
    if not accepted:
        print("Coadd cancelled; no files were written.")
        return 0
    print(f"Accepted {len(accepted)} of {len(candidates)} exposure(s).")
    if not ask_yes_no("Write this coadd selection and PypeIt input file?"):
        print("No files were written.")
        return 0
    out_dir = root / "Coadds" / f"{safe_name(args.target)}_{args.channel}_{safe_name(args.setup)}"
    coadd_file, output = write_coadd_input(out_dir, args.target, args.channel, args.setup, slit_id, accepted)
    print(f"Coadd input: {coadd_file}\nExpected output: {output}")
    if ask_yes_no("Run PypeIt coaddition now?"):
        return subprocess.run(["pypeit_coadd_1dspec", str(coadd_file), "--par_outfile", str(out_dir / "coadd1d.par")]).returncode
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
