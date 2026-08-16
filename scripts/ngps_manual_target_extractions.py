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
    create_target_copy,
    exposure_from_spec2d,
    manual_value,
    read_traces,
    replace_target_products,
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


def automatic_fwhm_by_slit(frame: Frame) -> dict[int, float]:
    """Read PypeIt's measured FWHM for each automatically found slicer trace."""
    spec1d = frame.spec2d.parent / frame.spec2d.name.replace("spec2d_", "spec1d_", 1)
    if not spec1d.is_file():
        return {}
    result: dict[int, float] = {}
    with fits.open(spec1d, memmap=False) as hdul:
        for hdu in hdul[1:]:
            identifier = slit_id(hdu.name)
            value = hdu.header.get("FWHM")
            if identifier is None or value is None:
                continue
            try:
                fwhm = float(value)
            except (TypeError, ValueError):
                continue
            if np.isfinite(fwhm) and fwhm > 0:
                result[identifier] = fwhm
    return result


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
    fwhms = automatic_fwhm_by_slit(frame)
    for slit in slits:
        centre = centre_for_slit(slit, traces, rows)
        middle = len(rows) // 2
        automatic.append(Selection(
            float(centre[middle]), float(rows[middle]), fwhms.get(slit[0], 4.0)
        ))
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
    frame: Frame, offset: float, fwhm: float
) -> tuple[np.ndarray, np.ndarray] | None:
    """Make a live, display-only 1D aperture sum at an aligned offset."""
    image, offsets, _ = aligned_image(frame)
    aperture = np.abs(offsets - offset) <= fwhm / 2
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


def display_fwhm(frame: Frame) -> float:
    """One representative width for the aligned dashboard display only."""
    _, _, selections = aligned_image(frame)
    widths = [selection.fwhm for selection in selections if selection.fwhm > 0]
    return float(np.median(widths)) if widths else 4.0


def display_quality_mask(flux: np.ndarray) -> np.ndarray:
    """Mask implausible endpoint outliers for every channel's display spectrum."""
    mask = np.isfinite(flux)
    edge = max(5, int(.01 * len(flux)))
    core = flux[edge:-edge][np.isfinite(flux[edge:-edge])]
    if core.size < 10:
        return mask
    centre = np.nanmedian(core)
    scatter = 1.4826 * np.nanmedian(np.abs(core - centre))
    amplitude = np.nanpercentile(np.abs(core - centre), 95)
    limit = max(20 * scatter, 20 * amplitude, 1.0)
    endpoints = np.zeros(len(flux), dtype=bool)
    endpoints[:edge] = True
    endpoints[-edge:] = True
    mask[endpoints & (np.abs(flux - centre) > limit)] = False
    return mask


def audit_path(root: Path, target: str, exposure: str) -> Path:
    output = root / "ExtractionQA" / safe_name(target)
    output.mkdir(parents=True, exist_ok=True)
    return output / f"ngps_extraction_review_{exposure}.pdf"


def selections_for_offsets(frame: Frame, offsets: list[float]) -> list[Selection]:
    _, _, base = aligned_image(frame)
    return [Selection(item.spatial + offset, item.spectral, item.fwhm) for offset in offsets for item in base]


def review_group(root: Path, target: str, exposure: str, frames: dict[str, Frame], interactive: bool) -> tuple[str, dict[str, float]]:
    """Save a dashboard.  In interactive mode return the chosen extraction decision."""
    # Fits a typical laptop display at Matplotlib's default 100 dpi while
    # preserving a readable four-channel review layout.
    figure = plt.figure(figsize=(14, 7.6))
    grid = GridSpec(2, 5, figure=figure, width_ratios=(1, 1, 1, 1, 1.08), height_ratios=(1, 1))
    axes = {channel: figure.add_subplot(grid[0, index]) for index, channel in enumerate(CHANNELS)}
    profile_axis = figure.add_subplot(grid[0, 4])
    spectrum_axis = figure.add_subplot(grid[1, :4])
    control_axis = figure.add_subplot(grid[1, 4])
    control_axis.axis("off")
    figure.suptitle(f"{target}  |  exposure {exposure}  |  NGPS extraction review", fontsize=15)
    selected: dict[str, float] = {}
    channel_fwhm: dict[str, float] = {}
    state = {
        "decision": "automatic" if not interactive else "cancel",
        "manual": False,
        "channel_only": False,
        "focus_channel": None,
    }

    def add_final_mode_label(label: str) -> None:
        """Add a mode badge only to the final saved audit record."""
        manual = label == "MANUAL MODE"
        figure.text(
            .9025, .445, label, ha="center", va="center", fontsize=10,
            fontweight="bold", color="#8A4B00" if manual else "#1D5F34",
            bbox={
                "boxstyle": "round,pad=.45",
                "facecolor": "#FFE6B3" if manual else "#D7F2DF",
                "edgecolor": "#E8B65D" if manual else "#92C9A5",
            },
        )

    for channel in CHANNELS:
        axis = axes[channel]
        frame = frames.get(channel)
        if frame is None:
            axis.text(.5, .5, "No reduced frame", ha="center", va="center", transform=axis.transAxes)
            axis.set_title(channel.upper())
            continue
        image, offsets, _ = aligned_image(frame)
        channel_fwhm[channel] = display_fwhm(frame)
        finite = image[np.isfinite(image)]
        limits = np.percentile(finite, (5, 99)) if finite.size else (-1, 1)
        axis.imshow(image, origin="lower", aspect="auto", cmap="viridis", vmin=limits[0], vmax=limits[1], extent=(offsets[0], offsets[-1], 0, image.shape[0] - 1))
        axis.set_box_aspect(1)
        axis.axvline(0, color="gold", lw=1.2, label="PypeIt automatic centre")
        axis.set_title(f"{channel.upper()}: aligned three-slicer", fontsize=9)
        axis.set_xlabel("Offset from automatic trace (pixels)", fontsize=8)
        if channel == "u":
            axis.set_ylabel("Spectral pixel", fontsize=8)
        axis.tick_params(labelsize=8)
        profile = np.nanmedian(image, axis=0)
        scale = np.nanmax(np.abs(profile))
        if np.isfinite(scale) and scale > 0:
            profile_axis.plot(offsets, profile / scale, color=COLOURS[channel], label=channel.upper())
    profile_axis.axvline(0, color="0.7", lw=.8)
    profile_axis.set_title("Spatial profiles\none colour per channel", fontsize=9)
    profile_axis.set_xlabel("Offset from automatic trace (pixels)", fontsize=8)
    profile_axis.set_ylabel("Normalised profile", fontsize=8)
    profile_axis.tick_params(labelsize=8)
    profile_axis.legend(loc="best")
    spectrum_axis.set_xlabel("Wavelength (Å)")

    spectrum_lines: list[object] = []

    def redraw_spectra(manual_offsets: dict[str, float] | None = None) -> None:
        for line in spectrum_lines:
            line.remove()
        spectrum_lines.clear()
        displayed: dict[str, np.ndarray] = {}
        for channel in CHANNELS:
            frame = frames.get(channel)
            if frame is None:
                continue
            spectrum = (manual_quicklook_spectrum(
                            frame, manual_offsets[channel], channel_fwhm[channel]
                        )
                        if manual_offsets is not None and channel in manual_offsets
                        else quicklook_spectrum(frame))
            if spectrum is None:
                continue
            wave, flux = spectrum
            finite_spec = np.isfinite(wave) & display_quality_mask(flux)
            if np.any(finite_spec):
                displayed[channel] = flux[finite_spec]
                spectrum_lines.append(spectrum_axis.plot(
                    wave[finite_spec], flux[finite_spec], color=COLOURS[channel],
                    lw=.65, label=channel.upper())[0])
        spectrum_axis.set_title(
            "Manual-aperture quick-look spectra"
            if manual_offsets else
            "Quick-look central-slicer 1D spectra"
        )
        spectrum_axis.set_ylabel("Detector counts")
        if spectrum_lines:
            spectrum_axis.relim()
            spectrum_axis.autoscale_view()
            focus = state["focus_channel"]
            if focus in displayed:
                values = displayed[focus]
                lower, upper = np.nanpercentile(values, (1, 99))
                span = upper - lower
                if np.isfinite(span) and span > 0:
                    spectrum_axis.set_ylim(lower - .08 * span, upper + .08 * span)
                    spectrum_axis.set_title(
                        f"Quick-look central-slicer 1D spectra - y-scale: {focus.upper()}"
                    )
            spectrum_axis.legend(loc="best", ncol=4)

    redraw_spectra()

    selection_artists: list[object] = []
    def redraw_selections() -> None:
        for artist in selection_artists:
            artist.remove()
        selection_artists.clear()
        for channel, offset in selected.items():
            axis = axes[channel]
            width = channel_fwhm[channel]
            selection_artists.append(axis.axvspan(offset - width / 2, offset + width / 2, color="tab:red", alpha=.24))
            selection_artists.append(axis.axvline(offset, color="tab:red", lw=.9))
            selection_artists.append(profile_axis.axvspan(offset - width / 2, offset + width / 2, color=COLOURS[channel], alpha=.14))
            selection_artists.append(profile_axis.axvline(offset, color=COLOURS[channel], lw=.9, alpha=.9))
        profile_axis.set_title(
            "Spatial profiles\nmanual apertures" if selected
            else "Spatial profiles\none colour per channel"
        )
        redraw_spectra(selected if selected else None)
        figure.canvas.draw_idle()

    def click(event) -> None:
        if not state["manual"] or event.inaxes not in axes.values() or event.xdata is None:
            return
        channel = next(name for name, axis in axes.items() if axis is event.inaxes)
        if state["channel_only"]:
            selected[channel] = float(event.xdata)
        else:
            # Default: one linked sky position across all four channels.
            selected.update({name: float(event.xdata) for name in frames})
        redraw_selections()

    def accept_auto(event) -> None:
        state["decision"] = "automatic"
        plt.close(figure)

    def begin_manual(event) -> None:
        state["manual"] = True
        state["channel_only"] = False

    def adjust_this_channel(event) -> None:
        state["manual"] = True
        state["channel_only"] = True
        print("Click a channel panel to move only that channel's manual aperture.")

    def return_to_automatic(event) -> None:
        state["manual"] = False
        state["channel_only"] = False
        state["focus_channel"] = None
        selected.clear()
        redraw_selections()

    def renormalise(channel: str) -> None:
        state["focus_channel"] = channel
        redraw_spectra(selected if selected else None)
        figure.canvas.draw_idle()

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
            ("Adjust this channel only", adjust_this_channel, "#D7F4F2", "#B8E8E4"),
            ("Return to automatic", return_to_automatic, "#E8E1FF", "#D3C9F2"),
            ("Accept manual", accept_manual, "#FFE6B3", "#F5D296"),
            ("Cancel", cancel, "#FFD9D9", "#F2BFBF"),
        ]
        for index, (label, callback, colour, hover_colour) in enumerate(button_specs):
            button_axis = figure.add_axes((.835, .16 + (.052 * (5 - index)), .135, .041))
            button = Button(button_axis, label, color=colour, hovercolor=hover_colour)
            button.on_clicked(callback)
            button_widgets.append(button)
        for index, channel in enumerate(CHANNELS):
            column, row = index % 2, index // 2
            button_axis = figure.add_axes((.835 + (.07 * column), .095 + (.03 * (1 - row)), .065, .024))
            button = Button(button_axis, f"Re-norm {channel.upper()}", color="#F1F3F5", hovercolor="#DDE4EA")
            button.on_clicked(lambda event, selected_channel=channel: renormalise(selected_channel))
            button_widgets.append(button)
    figure.canvas.mpl_connect("button_press_event", click)
    figure.subplots_adjust(left=.045, right=.985, bottom=.09, top=.86, wspace=.22, hspace=.36)
    if interactive:
        plt.show()
        # Cancel (including closing the window) is deliberately a true no-op:
        # preserve the previous audit PDF and every PypeIt product exactly.
        if state["decision"] == "cancel":
            plt.close(figure)
            return state["decision"], selected
        # Button axes are useful only in the live review.  Remove them before
        # saving the accepted result as an uncluttered scientific record.
        for button in button_widgets:
            button.ax.remove()
        if state["decision"] in {"automatic", "manual"}:
            add_final_mode_label(
                "MANUAL MODE" if state["decision"] == "manual" else "AUTO MODE"
            )
    else:
        # --auto writes an explicitly labelled automatic audit record.
        add_final_mode_label("AUTO MODE")
    # There is exactly one audit image per target/exposure.  In interactive
    # mode it is written only after the window closes, so the accepted manual
    # bands and live quick-look spectrum replace the earlier automatic review.
    figure.canvas.draw()
    output = audit_path(root, target, exposure)
    figure.savefig(output)
    print(f"Saved review PDF: {output}")
    plt.close(figure)
    return state["decision"], selected


def rerun_selected_exposure(
    frames: dict[str, Frame], offsets: dict[str, float] | None = None,
) -> int:
    """Rerun and install only the accepted exposure's channel products."""
    manual = offsets is not None
    for channel, frame in sorted(frames.items()):
        if manual and channel not in offsets:
            print(f"{channel.upper()}: automatic extraction retained.")
            continue
        source = next(iter(sorted(frame.spec2d.parent.parent.glob("*.pypeit"))), None)
        if source is None:
            print(f"WARNING: no PypeIt file for {frame.spec2d}")
            continue
        selections = selections_for_offsets(frame, [offsets[channel]]) if manual else None
        if selections is not None:
            print(f"{channel.upper()} PypeIt manual value:\n  {manual_value(selections)}")
        try:
            run_dir, target_pypeit = create_target_copy(
                source, exposure_from_spec2d(frame.spec2d), selections,
            )
        except (OSError, RuntimeError) as error:
            print(f"ERROR: could not create one-exposure setup: {error}")
            return 1
        mode = "manual" if manual else "automatic"
        print(f"{channel.upper()} one-exposure {mode} setup: {run_dir}")
        status = subprocess.run(["run_pypeit", target_pypeit.name], cwd=run_dir).returncode
        if status != 0:
            return status
        try:
            count = replace_target_products(run_dir, source.parent, frame.exposure)
        except (OSError, RuntimeError) as error:
            print(f"ERROR: reduction succeeded but products were not installed: {error}")
            return 1
        print(f"{channel.upper()}: replaced {count} derived product(s) for exposure {frame.exposure}.")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description="Review NGPS 2D extractions in four-channel per-exposure dashboards.")
    parser.add_argument("date", help="UT date, e.g. 20260623")
    parser.add_argument("--target", help="Target name from the science inventory")
    parser.add_argument("--channel", choices=CHANNELS, help="Optional channel filter")
    parser.add_argument("--exposure", help="Review one four-digit exposure, e.g. 0121")
    parser.add_argument("--auto", action="store_true", help="Save PDFs only; do not open review windows")
    parser.add_argument("--all", action="store_true", help="Internal: review every reduced exposure (used by ngps_reduce_all_configs.py)")
    args = parser.parse_args()
    if not args.target and not args.all:
        parser.error("provide --target, or use --all")
    if args.target and not args.all and not args.exposure:
        parser.error("--target review requires --exposure, e.g. --exposure 0121")
    if args.exposure and not re.fullmatch(r"\d{4}", args.exposure):
        parser.error("--exposure must be four digits, e.g. 0121")
    root = Path(os.environ.get("NGPS_WORK_ROOT", Path.home() / "ngps_data" / "work")) / args.date
    frames = discover_frames(root)
    if args.target:
        frames = [frame for frame in frames if frame.target.casefold() == args.target.casefold()]
    if args.channel:
        frames = [frame for frame in frames if frame.channel == args.channel]
    if args.exposure:
        frames = [frame for frame in frames if frame.exposure == args.exposure]
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
            # The dashboard PDF above has already replaced the automatic PDF.
            # Manual detector products are rebuilt for this exposure only.
            status = rerun_selected_exposure(group, offsets)
            if status != 0:
                return status
        elif decision == "automatic":
            status = rerun_selected_exposure(group)
            if status != 0:
                return status
        else:
            print("Review cancelled; existing PDFs and PypeIt products were not changed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
