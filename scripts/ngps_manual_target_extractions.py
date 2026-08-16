#!/usr/bin/env python3
"""Review NGPS extractions as four-channel, per-exposure dashboards.

The dashboard is a quality-assurance display.  It aligns the three image-slicer
slits around PypeIt's trace and combines them only for viewing.  It never
replaces the detector images or performs a science coadd.
"""

from __future__ import annotations

import argparse
import os
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from matplotlib.gridspec import GridSpec
from matplotlib.widgets import Button

from ngps_interactive_extract import (
    Selection,
    create_manual_copy,
    exposure_from_spec2d,
    manual_value,
    read_traces,
)


CHANNELS = ("u", "g", "r", "i")
COLOURS = {"u": "tab:purple", "g": "tab:green", "r": "tab:red", "i": "tab:orange"}
NAME_RE = re.compile(r"spec2d_ngps_\d+_(?P<exposure>\d{4})-(?P<target>.+?)_NGPS_(?P<channel>[ugri])_")


@dataclass
class Frame:
    channel: str
    target: str
    exposure: str
    spec2d: Path


def safe_name(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.+-]+", "_", text).strip("_")


def parse_frame(path: Path) -> Frame | None:
    match = NAME_RE.search(path.name)
    if match is None:
        return None
    return Frame(match["channel"], match["target"], match["exposure"], path)


def discover_frames(root: Path) -> list[Frame]:
    result = []
    for path in sorted(root.glob("manual_setup_*/*/Science/spec2d_*.fits")):
        # The automatic dashboard reviews the original PypeIt products.  Manual
        # copies are kept separate so a revision never silently becomes baseline.
        if "_manual_" in path.parent.parent.name:
            continue
        frame = parse_frame(path)
        if frame is not None:
            result.append(frame)
    return result


def group_frames(frames: list[Frame]) -> dict[tuple[str, str], dict[str, Frame]]:
    groups: dict[tuple[str, str], dict[str, Frame]] = {}
    for frame in frames:
        groups.setdefault((frame.target.casefold(), frame.exposure), {})[frame.channel] = frame
    return groups


def frame_arrays(frame: Frame) -> tuple[np.ndarray, np.ndarray, list[tuple[int, np.ndarray, np.ndarray]], list[tuple[str, np.ndarray, np.ndarray]]]:
    with fits.open(frame.spec2d, memmap=False) as hdul:
        image = np.asarray(hdul["DET01-SCIIMG"].data, dtype=float) - np.asarray(hdul["DET01-SKYMODEL"].data, dtype=float)
        ivar = np.asarray(hdul["DET01-IVARMODEL"].data, dtype=float) if "DET01-IVARMODEL" in hdul else np.ones_like(image)
        slits = [(int(row["spat_id"]), np.asarray(row["left_init"]), np.asarray(row["right_init"])) for row in hdul["DET01-SLITS"].data]
    spec1d = frame.spec2d.parent / frame.spec2d.name.replace("spec2d_", "spec1d_", 1)
    return image, ivar, slits, read_traces(spec1d if spec1d.is_file() else None)


def slit_id(trace_name: str) -> int | None:
    match = re.search(r"SLIT(\d+)", trace_name)
    return int(match.group(1)) if match else None


def centre_for_slit(slit: tuple[int, np.ndarray, np.ndarray], traces: list[tuple[str, np.ndarray, np.ndarray]], rows: np.ndarray) -> np.ndarray:
    identifier, left, right = slit
    choices = [(spatial, spectral) for name, spatial, spectral in traces if slit_id(name) == identifier]
    if choices:
        spatial, spectral = choices[0]
        order = np.argsort(spectral)
        return np.interp(rows, spectral[order], spatial[order])
    return (left + right) / 2


def aligned_image(frame: Frame, half_width: int = 38) -> tuple[np.ndarray, np.ndarray, list[Selection]]:
    """Return a slicer-aligned diagnostic image and manual positions at offset zero."""
    image, ivar, slits, traces = frame_arrays(frame)
    rows = np.arange(image.shape[0])
    offsets = np.arange(-half_width, half_width + 1)
    numerator = np.zeros((len(rows), len(offsets)))
    weights = np.zeros_like(numerator)
    automatic: list[Selection] = []
    for slit in slits:
        centre = centre_for_slit(slit, traces, rows)
        middle = len(rows) // 2
        automatic.append(Selection(float(centre[middle]), float(rows[middle]), 4.0))
        x = np.rint(centre[:, None] + offsets[None, :]).astype(int)
        valid = (x >= 0) & (x < image.shape[1])
        values = np.full_like(numerator, np.nan, dtype=float)
        value_weights = np.zeros_like(numerator)
        yy = np.broadcast_to(rows[:, None], x.shape)
        values[valid] = image[yy[valid], x[valid]]
        local_ivar = ivar[yy[valid], x[valid]]
        value_weights[valid] = np.where(np.isfinite(local_ivar) & (local_ivar > 0), local_ivar, 0.0)
        valid_values = np.isfinite(values) & (value_weights > 0)
        numerator[valid_values] += values[valid_values] * value_weights[valid_values]
        weights[valid_values] += value_weights[valid_values]
    combined = np.divide(numerator, weights, out=np.full_like(numerator, np.nan), where=weights > 0)
    return combined, offsets, automatic


def quicklook_spectrum(frame: Frame) -> tuple[np.ndarray, np.ndarray] | None:
    spec1d = frame.spec2d.parent / frame.spec2d.name.replace("spec2d_", "spec1d_", 1)
    if not spec1d.is_file():
        return None
    with fits.open(spec1d, memmap=False) as hdul:
        candidates = [hdu for hdu in hdul[1:] if hdu.name.startswith("SPAT") and getattr(hdu, "data", None) is not None]
        if not candidates:
            return None
        # The middle slit is the least vignetted quick-look choice.  This is not a final coadd.
        hdu = candidates[len(candidates) // 2]
        names = hdu.data.dtype.names or ()
        wave_key = "OPT_WAVE" if "OPT_WAVE" in names else "BOX_WAVE"
        flux_key = "OPT_FLAM" if "OPT_FLAM" in names else ("OPT_COUNTS" if "OPT_COUNTS" in names else "BOX_COUNTS")
        if wave_key not in names or flux_key not in names:
            return None
        return np.asarray(hdu.data[wave_key]), np.asarray(hdu.data[flux_key])


def audit_path(root: Path, target: str, exposure: str) -> Path:
    output = root / "ExtractionQA" / safe_name(target)
    output.mkdir(parents=True, exist_ok=True)
    return output / f"ngps_extraction_review_{exposure}.pdf"


def selections_for_offsets(frame: Frame, offsets: list[float]) -> list[Selection]:
    _, _, base = aligned_image(frame)
    return [Selection(item.spatial + offset, item.spectral, item.fwhm) for offset in offsets for item in base]


def review_group(root: Path, target: str, exposure: str, frames: dict[str, Frame], interactive: bool, maximum: int = 3) -> tuple[str, list[float]]:
    """Save a dashboard.  In interactive mode return the chosen extraction decision."""
    figure = plt.figure(figsize=(17, 11))
    grid = GridSpec(3, 3, figure=figure, width_ratios=(1, 1, 1.05), height_ratios=(1, 1, .75))
    axes = {"u": figure.add_subplot(grid[0, 0]), "g": figure.add_subplot(grid[0, 1]), "r": figure.add_subplot(grid[1, 0]), "i": figure.add_subplot(grid[1, 1])}
    profile_axis = figure.add_subplot(grid[:2, 2])
    spectrum_axis = figure.add_subplot(grid[2, :2])
    control_axis = figure.add_subplot(grid[2, 2])
    control_axis.axis("off")
    figure.suptitle(f"{target}  |  exposure {exposure}  |  NGPS extraction review", fontsize=15)
    selected: list[float] = []
    state = {"decision": "automatic" if not interactive else "cancel", "manual": False}

    for channel in CHANNELS:
        axis = axes[channel]
        frame = frames.get(channel)
        if frame is None:
            axis.text(.5, .5, "No reduced frame", ha="center", va="center", transform=axis.transAxes)
            axis.set_title(channel.upper())
            continue
        image, offsets, _ = aligned_image(frame)
        finite = image[np.isfinite(image)]
        limits = np.percentile(finite, (5, 99)) if finite.size else (-1, 1)
        axis.imshow(image, origin="lower", aspect="auto", cmap="gray", vmin=limits[0], vmax=limits[1], extent=(offsets[0], offsets[-1], 0, image.shape[0] - 1))
        axis.axvline(0, color="gold", lw=1.2, label="PypeIt automatic centre")
        axis.set_title(f"{channel.upper()}: aligned three-slicer diagnostic")
        axis.set_xlabel("Offset from automatic trace (pixels)")
        axis.set_ylabel("Spectral pixel")
        profile = np.nanmedian(image, axis=0)
        scale = np.nanmax(np.abs(profile))
        if np.isfinite(scale) and scale > 0:
            profile_axis.plot(profile / scale, offsets, color=COLOURS[channel], label=channel.upper())
        spectrum = quicklook_spectrum(frame)
        if spectrum is not None:
            wave, flux = spectrum
            finite_spec = np.isfinite(wave) & np.isfinite(flux)
            spectrum_axis.plot(wave[finite_spec], flux[finite_spec], color=COLOURS[channel], lw=.65, label=channel.upper())

    profile_axis.axvline(0, color="0.7", lw=.8)
    profile_axis.set_title("Spatial profiles (one colour per channel)")
    profile_axis.set_xlabel("Normalised sky-subtracted profile")
    profile_axis.set_ylabel("Offset from automatic trace (pixels)")
    profile_axis.legend(loc="best")
    spectrum_axis.set_title("Quick-look central-slicer 1D spectra (not a coadd)")
    spectrum_axis.set_xlabel("Vacuum wavelength (Å)")
    spectrum_axis.set_ylabel("FLAM or counts")
    spectrum_axis.legend(loc="best", ncol=4)

    selection_artists: list[object] = []
    def redraw_selections() -> None:
        for artist in selection_artists:
            artist.remove()
        selection_artists.clear()
        for index, offset in enumerate(selected):
            for axis in axes.values():
                selection_artists.append(axis.axvspan(offset - 2, offset + 2, color="tab:red", alpha=.24))
                selection_artists.append(axis.axvline(offset, color="tab:red", lw=.9))
            selection_artists.append(profile_axis.axhline(offset, color="tab:red", lw=.9, alpha=.8))
        figure.canvas.draw_idle()

    def click(event) -> None:
        if not state["manual"] or event.inaxes not in axes.values() or event.xdata is None:
            return
        if len(selected) >= maximum:
            print(f"Maximum of {maximum} source component(s) reached.")
            return
        selected.append(float(event.xdata))
        redraw_selections()

    def accept_auto(event) -> None:
        state["decision"] = "automatic"
        plt.close(figure)

    def begin_manual(event) -> None:
        state["manual"] = True
        control_message.set_text("Click a red extraction band in any channel panel.\nUse Add another component for a dual/triplet; then Accept manual.")
        figure.canvas.draw_idle()

    def add_component(event) -> None:
        state["manual"] = True
        control_message.set_text(f"Click the next component (up to {maximum}).")
        figure.canvas.draw_idle()

    def accept_manual(event) -> None:
        if not selected:
            control_message.set_text("Choose at least one position first.")
            figure.canvas.draw_idle()
            return
        state["decision"] = "manual"
        plt.close(figure)

    def cancel(event) -> None:
        state["decision"] = "cancel"
        plt.close(figure)

    control_message = control_axis.text(.5, .94, "Gold line: PypeIt automatic centre.\nRed bands: your manual choices.", ha="center", va="top", wrap=True, transform=control_axis.transAxes)
    button_widgets: list[Button] = []
    if interactive:
        button_specs = [
            ("Accept automatic", accept_auto), ("Manual extraction", begin_manual),
            ("Add another component", add_component), ("Accept manual", accept_manual), ("Cancel", cancel),
        ]
        for index, (label, callback) in enumerate(button_specs):
            button_axis = figure.add_axes((.72, .12 + (.055 * (4 - index)), .21, .04))
            button = Button(button_axis, label)
            button.on_clicked(callback)
            button_widgets.append(button)
    else:
        control_axis.text(.5, .45, "Automatic run\n\nThis PDF records PypeIt's automatic extraction.\nUse ngps_manual_target_extractions.py\nto revise it.", ha="center", va="center", wrap=True, transform=control_axis.transAxes)
    figure.canvas.mpl_connect("button_press_event", click)
    figure.subplots_adjust(left=.06, right=.98, bottom=.08, top=.88, wspace=.17, hspace=.32)
    output = audit_path(root, target, exposure)
    figure.savefig(output)
    print(f"Saved review PDF: {output}")
    if interactive:
        plt.show()
        # Save again so the accepted manual bands replace the automatic PDF.
        figure.savefig(output)
    plt.close(figure)
    return state["decision"], selected


def write_manual(root: Path, frames: dict[str, Frame], offsets: list[float], run_manual: bool) -> int:
    for channel, frame in sorted(frames.items()):
        source = next(iter(sorted(frame.spec2d.parent.parent.glob("*.pypeit"))), None)
        if source is None:
            print(f"WARNING: no PypeIt file for {frame.spec2d}")
            continue
        selections = selections_for_offsets(frame, offsets)
        print(f"{channel.upper()} PypeIt manual value:\n  {manual_value(selections)}")
        try:
            manual_dir, manual_pypeit = create_manual_copy(source, exposure_from_spec2d(frame.spec2d), selections)
        except FileExistsError as error:
            print(f"WARNING: {error}")
            continue
        print(f"Manual setup: {manual_dir}")
        if run_manual:
            status = subprocess.run(["run_pypeit", manual_pypeit.name], cwd=manual_dir).returncode
            if status != 0:
                return status
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description="Review NGPS 2D extractions in four-channel per-exposure dashboards.")
    parser.add_argument("date", help="UT date, e.g. 20260623")
    parser.add_argument("--target", help="Target name from the science inventory")
    parser.add_argument("--channel", choices=CHANNELS, help="Optional channel filter")
    parser.add_argument("--auto", action="store_true", help="Save PDFs only; do not open review windows")
    parser.add_argument("--all", action="store_true", help="Review every reduced science exposure (used by the reduction driver)")
    parser.add_argument("--run-manual", action="store_true", help="Run each copied manual setup after acceptance")
    parser.add_argument("--max-components", type=int, default=3, help="Maximum source components: 1 to 3")
    args = parser.parse_args()
    if not args.target and not args.all:
        parser.error("provide --target, or use --all")
    if not 1 <= args.max_components <= 3:
        parser.error("--max-components must be 1 to 3")
    root = Path(os.environ.get("NGPS_WORK_ROOT", Path.home() / "ngps_data" / "work")) / args.date
    frames = discover_frames(root)
    if args.target:
        frames = [frame for frame in frames if frame.target.casefold() == args.target.casefold()]
    if args.channel:
        frames = [frame for frame in frames if frame.channel == args.channel]
    groups = group_frames(frames)
    if not groups:
        parser.error("No reduced science spec2d files matched")
    for group_key in sorted(groups):
        _, exposure = group_key
        group = groups[group_key]
        target = next(iter(group.values())).target
        print(f"\n{'=' * 76}\n{target} | exposure {exposure} | channels: {', '.join(channel.upper() for channel in sorted(group))}\n{'=' * 76}")
        decision, offsets = review_group(root, target, exposure, group, not args.auto, args.max_components)
        if decision == "manual":
            # The dashboard PDF above has already replaced the automatic PDF.  Detector products
            # remain protected in copied manual setups.
            status = write_manual(root, group, offsets, args.run_manual)
            if status != 0:
                return status
        elif decision == "automatic":
            print("Automatic extraction retained; only the review PDF was written.")
        else:
            print("Review cancelled; the saved PDF remains as an audit record.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
