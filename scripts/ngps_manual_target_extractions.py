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
        left, right = slit[1], slit[2]
        # Exclude pixels beyond the curved slicer edges from the aligned view.
        valid = ((x >= 0) & (x < image.shape[1])
                 & (x >= left[:, None]) & (x <= right[:, None]))
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
    # Extraction review is intentionally before flux calibration.  Always use
    # the original Science product so the display matches the extraction that
    # the user is accepting or revising.
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


def manual_quicklook_spectrum(
    frame: Frame, offset: float, half_width: float = 2.0
) -> tuple[np.ndarray, np.ndarray] | None:
    """Make a live, display-only 1D aperture sum at an aligned offset."""
    image, offsets, _ = aligned_image(frame)
    aperture = np.abs(offsets - offset) <= half_width
    if not np.any(aperture):
        return None
    flux = np.nansum(image[:, aperture], axis=1)
    reference = quicklook_spectrum(frame)
    if reference is None:
        return None
    wave, _ = reference
    if len(wave) != len(flux):
        wave = np.interp(
            np.arange(len(flux)),
            np.linspace(0, len(flux) - 1, len(wave)),
            wave,
        )
    return wave, flux


def normalise_for_display(flux: np.ndarray) -> np.ndarray:
    """Scale one channel robustly for an extraction-quality quick-look."""
    finite = np.isfinite(flux)
    if not np.any(finite):
        return flux
    # The 90th percentile keeps a broad throughput rise in one channel from
    # compressing the other three, while remaining insensitive to rare spikes.
    scale = np.nanpercentile(np.abs(flux[finite]), 90)
    return flux / scale if np.isfinite(scale) and scale > 0 else flux


def audit_path(root: Path, target: str, exposure: str) -> Path:
    output = root / "ExtractionQA" / safe_name(target)
    output.mkdir(parents=True, exist_ok=True)
    return output / f"ngps_extraction_review_{exposure}.pdf"


def selections_for_offsets(frame: Frame, offsets: list[float]) -> list[Selection]:
    _, _, base = aligned_image(frame)
    return [Selection(item.spatial + offset, item.spectral, item.fwhm) for offset in offsets for item in base]


def review_group(root: Path, target: str, exposure: str, frames: dict[str, Frame], interactive: bool) -> tuple[str, list[float]]:
    """Save a dashboard.  In interactive mode return the chosen extraction decision."""
    figure = plt.figure(figsize=(22, 10.5))
    grid = GridSpec(2, 5, figure=figure, width_ratios=(1, 1, 1, 1, 1.08), height_ratios=(1, 1))
    axes = {channel: figure.add_subplot(grid[0, index]) for index, channel in enumerate(CHANNELS)}
    profile_axis = figure.add_subplot(grid[0, 4])
    spectrum_axis = figure.add_subplot(grid[1, :4])
    control_axis = figure.add_subplot(grid[1, 4])
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
        axis.imshow(image, origin="lower", aspect="auto", cmap="viridis", vmin=limits[0], vmax=limits[1], extent=(offsets[0], offsets[-1], 0, image.shape[0] - 1))
        axis.set_box_aspect(1)
        axis.axvline(0, color="gold", lw=1.2, label="PypeIt automatic centre")
        axis.set_title(f"{channel.upper()}: aligned three-slicer diagnostic")
        axis.set_xlabel("Offset from automatic trace (pixels)")
        axis.set_ylabel("Spectral pixel")
        profile = np.nanmedian(image, axis=0)
        scale = np.nanmax(np.abs(profile))
        if np.isfinite(scale) and scale > 0:
            profile_axis.plot(offsets, profile / scale, color=COLOURS[channel], label=channel.upper())
    profile_axis.axvline(0, color="0.7", lw=.8)
    profile_axis.set_title("Spatial profiles (one colour per channel)")
    profile_axis.set_xlabel("Offset from automatic trace (pixels)")
    profile_axis.set_ylabel("Normalised sky-subtracted profile")
    profile_axis.legend(loc="best")
    spectrum_axis.set_xlabel("Wavelength (Å)")

    spectrum_lines: list[object] = []

    def redraw_spectra(manual_offset: float | None = None) -> None:
        for line in spectrum_lines:
            line.remove()
        spectrum_lines.clear()
        display_values: list[np.ndarray] = []
        for channel in CHANNELS:
            frame = frames.get(channel)
            if frame is None:
                continue
            spectrum = (manual_quicklook_spectrum(frame, manual_offset)
                        if manual_offset is not None else quicklook_spectrum(frame))
            if spectrum is None:
                continue
            wave, flux = spectrum
            display_flux = normalise_for_display(flux)
            finite_spec = np.isfinite(wave) & np.isfinite(display_flux)
            if np.any(finite_spec):
                display_values.append(display_flux[finite_spec])
                spectrum_lines.append(spectrum_axis.plot(
                    wave[finite_spec], display_flux[finite_spec], color=COLOURS[channel],
                    lw=.65, label=channel.upper())[0])
        spectrum_axis.set_title(
            "Manual-aperture quick-look spectra"
            if manual_offset is not None else
            "Quick-look central-slicer 1D spectra, individually scaled"
        )
        if manual_offset is not None:
            spectrum_axis.set_ylabel("Normalised detector counts")
        else:
            spectrum_axis.set_ylabel("Normalised detector counts")
        if display_values:
            values = np.concatenate(display_values)
            lower, upper = np.nanpercentile(values, (1, 99))
            span = upper - lower
            if np.isfinite(span) and span > 0:
                spectrum_axis.set_ylim(lower - .08 * span, upper + .08 * span)
        if spectrum_lines:
            spectrum_axis.legend(loc="best", ncol=4)

    redraw_spectra()

    selection_artists: list[object] = []
    def redraw_selections() -> None:
        for artist in selection_artists:
            artist.remove()
        selection_artists.clear()
        for offset in selected:
            for axis in axes.values():
                selection_artists.append(axis.axvspan(offset - 2, offset + 2, color="tab:red", alpha=.24))
                selection_artists.append(axis.axvline(offset, color="tab:red", lw=.9))
            selection_artists.append(profile_axis.axvspan(offset - 2, offset + 2, color="tab:red", alpha=.18))
            selection_artists.append(profile_axis.axvline(offset, color="tab:red", lw=.9, alpha=.8))
        profile_axis.set_title(
            "Spatial profiles with manual aperture" if selected
            else "Spatial profiles (one colour per channel)"
        )
        redraw_spectra(selected[0] if selected else None)
        figure.canvas.draw_idle()

    def click(event) -> None:
        if not state["manual"] or event.inaxes not in axes.values() or event.xdata is None:
            return
        # One component for now: clicking again moves the same manual aperture.
        selected[:] = [float(event.xdata)]
        redraw_selections()

    def accept_auto(event) -> None:
        state["decision"] = "automatic"
        plt.close(figure)

    def begin_manual(event) -> None:
        state["manual"] = True

    def return_to_automatic(event) -> None:
        state["manual"] = False
        selected.clear()
        redraw_selections()

    def accept_manual(event) -> None:
        if not selected:
            print("Choose a position in a channel panel before accepting manual extraction.")
            return
        state["decision"] = "manual"
        plt.close(figure)

    def cancel(event) -> None:
        state["decision"] = "cancel"
        plt.close(figure)

    button_widgets: list[Button] = []
    if interactive:
        button_specs = [
            ("Accept automatic", accept_auto, "#D7F2DF", "#BCE8CA"),
            ("Manual extraction", begin_manual, "#D7E9FF", "#BCD8F5"),
            ("Return to automatic", return_to_automatic, "#E8E1FF", "#D3C9F2"),
            ("Accept manual", accept_manual, "#FFE6B3", "#F5D296"),
            ("Cancel", cancel, "#FFD9D9", "#F2BFBF"),
        ]
        for index, (label, callback, colour, hover_colour) in enumerate(button_specs):
            button_axis = figure.add_axes((.835, .19 + (.075 * (4 - index)), .135, .05))
            button = Button(button_axis, label, color=colour, hovercolor=hover_colour)
            button.on_clicked(callback)
            button_widgets.append(button)
    figure.canvas.mpl_connect("button_press_event", click)
    figure.subplots_adjust(left=.045, right=.985, bottom=.09, top=.86, wspace=.22, hspace=.36)
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
    args = parser.parse_args()
    if not args.target and not args.all:
        parser.error("provide --target, or use --all")
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
        decision, offsets = review_group(root, target, exposure, group, not args.auto)
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
