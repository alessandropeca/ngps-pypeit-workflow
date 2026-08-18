#!/usr/bin/env python3
"""Inspect an NGPS 2D frame and create a safe manual-extraction copy.

This viewer intentionally avoids PypeIt's Ginga display utilities.  It reads
the reduced spec2d FITS file with Astropy and uses Matplotlib for selection.
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
from dataclasses import asdict, dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from matplotlib.widgets import Button, Slider


@dataclass
class Selection:
    spatial: float
    spectral: float
    fwhm: float


@dataclass
class ReviewResult:
    decision: str
    selections: list[Selection]


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


def manual_value(selections: list[Selection]) -> str:
    return ";".join(
        f"1:{item.spatial:.1f}:{item.spectral:.1f}:{item.fwhm:.1f}"
        for item in selections
    )


def spec1d_for(spec2d: Path) -> Path | None:
    matches = sorted(spec2d.parent.glob(spec2d.name.replace("spec2d_", "spec1d_", 1)))
    return matches[0] if matches else None


def read_traces(spec1d: Path | None) -> list[tuple[str, np.ndarray, np.ndarray]]:
    if spec1d is None:
        return []
    traces = []
    with fits.open(spec1d, memmap=False) as hdul:
        for hdu in hdul[1:]:
            if not hdu.name.startswith("SPAT") or getattr(hdu, "data", None) is None:
                continue
            names = hdu.data.dtype.names or ()
            if "TRACE_SPAT" in names and "trace_spec" in names:
                traces.append((hdu.name, np.asarray(hdu.data["TRACE_SPAT"]),
                               np.asarray(hdu.data["trace_spec"])))
    return traces


def interactive_select(spec2d: Path, initial_fwhm: float, maximum: int) -> ReviewResult:
    """Review automatic traces or define up to three replacement positions."""
    with fits.open(spec2d, memmap=False) as hdul:
        image = (np.asarray(hdul["DET01-SCIIMG"].data, dtype=float)
                 - np.asarray(hdul["DET01-SKYMODEL"].data, dtype=float))
        slit_data = hdul["DET01-SLITS"].data
        slits = [(int(row["spat_id"]), np.asarray(row["left_init"]),
                  np.asarray(row["right_init"])) for row in slit_data]

    # Rectify the display only.  NGPS slicers are tilted on the detector, so a
    # source at a fixed slicer position can look slanted.  We remove the median
    # slit-centre tilt for viewing, then convert each click back to raw detector
    # coordinates before writing any PypeIt manual extraction value.
    rows = np.arange(image.shape[0])
    slit_centres = np.array([(left + right) / 2 for _, left, right in slits])
    row_shift = np.nanmedian(slit_centres, axis=0)
    row_shift -= row_shift[len(row_shift) // 2]
    raw_columns = np.arange(image.shape[1])
    display_image = np.full_like(image, np.nan, dtype=float)
    for row, shift in enumerate(row_shift):
        source = np.rint(raw_columns + shift).astype(int)
        valid = (source >= 0) & (source < image.shape[1])
        display_image[row, valid] = image[row, source[valid]]

    finite = display_image[np.isfinite(display_image)]
    if finite.size == 0:
        raise RuntimeError("The sky-subtracted image has no finite pixels.")
    vmin, vmax = np.percentile(finite, (5, 99))
    figure, axis = plt.subplots(figsize=(13, 8))
    figure.subplots_adjust(bottom=0.22)
    axis.imshow(display_image, origin="lower", aspect="auto", cmap="gray",
                vmin=vmin, vmax=vmax, interpolation="nearest")
    axis.set_xlabel("Rectified spatial detector pixel")
    axis.set_ylabel("Spectral detector pixel")
    axis.set_title(
        "Rectified sky-subtracted 2D spectrum\n"
        "A fixed slicer position is vertical. "
        "Gold: PypeIt automatic traces. Accept them, or left-click a replacement position (max 3). "
        "Drag coloured marker: move it. Right-click: remove. "
        "FWHM slider applies to active marker."
    )
    for slit_id, left, right in slits:
        axis.plot(left - row_shift, rows, color="deepskyblue", lw=0.6, alpha=0.65)
        axis.plot(right - row_shift, rows, color="deepskyblue", lw=0.6, alpha=0.65)
        axis.text(np.nanmedian((left + right) / 2 - row_shift), image.shape[0] - 25,
                  f"SLIT{slit_id:04d}", color="deepskyblue", fontsize=8,
                  ha="center", va="top")
    for name, spatial, spectral in read_traces(spec1d_for(spec2d)):
        display_spatial = spatial - np.interp(spectral, rows, row_shift)
        axis.plot(display_spatial, spectral, color="gold", lw=0.8, alpha=0.9)
        middle = len(spatial) // 2
        axis.text(display_spatial[middle], spectral[middle], name, color="gold",
                  fontsize=7, rotation=90, va="bottom")

    selections: list[Selection] = []
    artists: list[tuple[object, object, object]] = []
    active: int | None = None
    dragging = False
    decision = {"value": "cancel"}
    slider_axis = figure.add_axes((0.2, 0.12, 0.6, 0.03))
    slider = Slider(slider_axis, "Active FWHM (pixels)", 1.0, 20.0,
                    valinit=initial_fwhm, valstep=0.5)

    def redraw() -> None:
        nonlocal artists
        for marker, left, right in artists:
            marker.remove()
            left.remove()
            right.remove()
        artists = []
        for index, item in enumerate(selections):
            colour = "tab:orange" if index == active else "tab:red"
            display_spatial = item.spatial - np.interp(item.spectral, rows, row_shift)
            marker = axis.plot(display_spatial, item.spectral, "o",
                               color=colour, markersize=7)[0]
            left = axis.axvline(display_spatial - item.fwhm / 2, color=colour, lw=0.9)
            right = axis.axvline(display_spatial + item.fwhm / 2, color=colour, lw=0.9)
            artists.append((marker, left, right))
        figure.canvas.draw_idle()

    def nearest(event) -> int | None:
        if event.xdata is None or event.ydata is None:
            return None
        for index, item in enumerate(selections):
            display_spatial = item.spatial - np.interp(item.spectral, rows, row_shift)
            dx, dy = (event.xdata - display_spatial) / 20, (event.ydata - item.spectral) / 60
            if dx * dx + dy * dy < 1:
                return index
        return None

    def on_press(event) -> None:
        nonlocal active, dragging
        if event.inaxes != axis:
            return
        if event.button == 3:
            index = nearest(event)
            if index is not None:
                selections.pop(index)
                active = min(index, len(selections) - 1) if selections else None
                redraw()
            return
        if event.button != 1 or event.xdata is None or event.ydata is None:
            return
        index = nearest(event)
        if index is not None:
            active, dragging = index, True
            slider.set_val(selections[index].fwhm)
            redraw()
        elif len(selections) < maximum:
            row = int(np.clip(round(event.ydata), 0, len(rows) - 1))
            selections.append(Selection(
                float(event.xdata + row_shift[row]), float(event.ydata), float(slider.val)
            ))
            active = len(selections) - 1
            redraw()
        else:
            print(f"Maximum of {maximum} manual selections reached.")

    def on_motion(event) -> None:
        if dragging and active is not None and event.inaxes == axis and event.xdata is not None and event.ydata is not None:
            row = int(np.clip(round(event.ydata), 0, len(rows) - 1))
            selections[active].spatial = float(event.xdata + row_shift[row])
            selections[active].spectral = float(event.ydata)
            redraw()

    def on_release(event) -> None:
        nonlocal dragging
        dragging = False

    def on_slider(value: float) -> None:
        if active is not None:
            selections[active].fwhm = float(value)
            redraw()

    def on_key(event) -> None:
        if event.key == "a":
            decision["value"] = "automatic"
            plt.close(figure)
        elif event.key in {"enter", "return"}:
            if selections:
                decision["value"] = "manual"
                plt.close(figure)
            else:
                print("Click Accept automatic, or add at least one manual position first.")
        elif event.key == "q":
            decision["value"] = "cancel"
            plt.close(figure)

    def accept_automatic(event) -> None:
        decision["value"] = "automatic"
        plt.close(figure)

    def accept_manual(event) -> None:
        if not selections:
            print("Add at least one manual position before accepting manual extraction.")
            return
        decision["value"] = "manual"
        plt.close(figure)

    def cancel(event) -> None:
        decision["value"] = "cancel"
        plt.close(figure)

    figure.canvas.mpl_connect("button_press_event", on_press)
    figure.canvas.mpl_connect("motion_notify_event", on_motion)
    figure.canvas.mpl_connect("button_release_event", on_release)
    figure.canvas.mpl_connect("key_press_event", on_key)
    slider.on_changed(on_slider)
    auto_axis = figure.add_axes((0.08, 0.035, 0.22, 0.05))
    manual_axis = figure.add_axes((0.39, 0.035, 0.24, 0.05))
    cancel_axis = figure.add_axes((0.72, 0.035, 0.16, 0.05))
    # Keep widget instances alive for the entire interactive session.  Without
    # these references, Matplotlib can garbage-collect the buttons and their
    # click callbacks stop responding.
    buttons = [
        Button(auto_axis, "Accept automatic (a)"),
        Button(manual_axis, "Accept manual positions"),
        Button(cancel_axis, "Cancel (q)"),
    ]
    buttons[0].on_clicked(accept_automatic)
    buttons[1].on_clicked(accept_manual)
    buttons[2].on_clicked(cancel)
    plt.show()
    return ReviewResult(decision["value"], selections)


def exposure_from_spec2d(spec2d: Path) -> str:
    match = re.search(r"spec2d_ngps_\d+_(\d{4})-", spec2d.name)
    if match is None:
        raise RuntimeError(f"Cannot identify exposure from {spec2d.name}")
    return match.group(1)


def write_manual_column(source: Path, destination: Path, exposure: str, value: str) -> None:
    output, in_data, manual_index, matched = [], False, None, 0
    for line in source.read_text().splitlines():
        stripped = line.strip()
        if stripped == "data read":
            in_data = True
            output.append(line)
            continue
        if stripped == "data end":
            in_data = False
            output.append(line)
            continue
        if not in_data or not stripped or stripped.startswith("#"):
            output.append(line)
            continue
        parts = [part.strip() for part in line.split("|")]
        if parts[0].lower() == "filename":
            lowered = [part.lower() for part in parts]
            if "manual" in lowered:
                manual_index = lowered.index("manual")
            else:
                parts.append("manual")
                manual_index = len(parts) - 1
            output.append(" | ".join(parts))
            continue
        if manual_index is None:
            output.append(line)
            continue
        while len(parts) <= manual_index:
            parts.append("")
        if parts[0].endswith(".fits") and f"_{exposure}.fits" in parts[0]:
            matched += 1
            parts[manual_index] = value
        output.append(" | ".join(parts))
    if manual_index is None or matched != 1:
        raise RuntimeError(f"Expected one PypeIt table row for exposure {exposure}; found {matched}.")
    destination.write_text("\n".join(output) + "\n")


def write_target_pypeit(
    source: Path, destination: Path, exposure: str, manual: str | None = None,
) -> None:
    """Write one science exposure, retaining all calibration and standard rows."""
    output, in_data, matched = [], False, 0
    filename_index = frametype_index = manual_index = None
    for line in source.read_text().splitlines():
        stripped = line.strip()
        if stripped == "data read":
            in_data = True
            output.append(line)
            continue
        if stripped == "data end":
            in_data = False
            output.append(line)
            continue
        if not in_data or not stripped or stripped.startswith("#"):
            output.append(line)
            continue
        parts = [part.strip() for part in line.split("|")]
        if parts[0].lower() == "filename":
            lowered = [part.lower() for part in parts]
            if "frametype" not in lowered:
                raise RuntimeError("PypeIt data table has no frametype column.")
            filename_index = lowered.index("filename")
            frametype_index = lowered.index("frametype")
            manual_index = lowered.index("manual") if "manual" in lowered else None
            if manual is not None and manual_index is None:
                parts.append("manual")
                manual_index = len(parts) - 1
            output.append(" | ".join(parts))
            continue
        if filename_index is None or frametype_index is None:
            output.append(line)
            continue
        while manual_index is not None and len(parts) <= manual_index:
            parts.append("")
        filename = parts[filename_index]
        frametype = parts[frametype_index].lower()
        selected = filename.endswith(".fits") and f"_{exposure}.fits" in filename
        if "science" in frametype and not selected:
            continue
        if selected:
            matched += 1
            if manual_index is not None:
                parts[manual_index] = manual or ""
        output.append(" | ".join(parts))
    if matched != 1:
        raise RuntimeError(
            f"Expected one science row for exposure {exposure}; found {matched}."
        )
    destination.write_text("\n".join(output) + "\n")


def create_target_copy(
    source_pypeit: Path, exposure: str, selections: list[Selection] | None = None,
) -> tuple[Path, Path]:
    """Create an isolated PypeIt setup for one reviewed science exposure."""
    setup_dir = source_pypeit.parent
    mode = "manual" if selections is not None else "automatic"
    run_root = setup_dir / ".ngps_target_runs"
    run_dir = run_root / f"{setup_dir.name}_{mode}_{exposure}"
    if run_dir.exists():
        shutil.rmtree(run_dir)
    shutil.copytree(
        setup_dir, run_dir,
        ignore=shutil.ignore_patterns(
            "Science", "QA", "ExtractionQA", "Fluxed", "FluxFiles", "Sensfunc",
            ".ngps_target_runs", "*.log", "*.pdf", "*.png",
        ),
    )
    copied = run_dir / source_pypeit.name
    destination = run_dir / f"{source_pypeit.stem}_{mode}_{exposure}.pypeit"
    write_target_pypeit(
        copied, destination, exposure,
        manual_value(selections) if selections is not None else None,
    )
    copied.unlink()
    if selections is not None:
        (run_dir / f"manual_selection_{exposure}.json").write_text(
            json.dumps([asdict(item) for item in selections], indent=2) + "\n"
        )
    return run_dir, destination


def replace_target_products(run_dir: Path, setup_dir: Path, exposure: str) -> int:
    """Replace only one exposure's detector-level products in its base setup."""
    source_science = run_dir / "Science"
    destination_science = setup_dir / "Science"
    patterns = (
        f"spec1d_*_{exposure}-*.fits",
        f"spec2d_*_{exposure}-*.fits",
        f"spec1d_*_{exposure}-*.txt",
    )
    spec1d = sorted(source_science.glob(patterns[0]))
    spec2d = sorted(source_science.glob(patterns[1]))
    if not spec1d or not spec2d:
        raise RuntimeError(
            f"One-exposure reduction did not produce spec1d and spec2d files for {exposure}."
        )
    destination_science.mkdir(parents=True, exist_ok=True)
    copied = 0
    for pattern in patterns:
        for previous in destination_science.glob(pattern):
            previous.unlink()
        for product in source_science.glob(pattern):
            shutil.copy2(product, destination_science / product.name)
            copied += 1
    return copied


def create_manual_copy(
    source_pypeit: Path, exposure: str, selections: list[Selection],
    replace_existing: bool = False,
) -> tuple[Path, Path]:
    setup_dir = source_pypeit.parent
    manual_dir = setup_dir.with_name(f"{setup_dir.name}_manual_{exposure}")
    if manual_dir.exists():
        if not replace_existing:
            raise FileExistsError(f"Manual setup already exists: {manual_dir}")
        # This directory contains only derived manual products for this exact
        # setup/exposure.  Rebuild it from the untouched automatic setup.
        shutil.rmtree(manual_dir)
    shutil.copytree(setup_dir, manual_dir, ignore=shutil.ignore_patterns(
        "Science", "QA", "ExtractionQA", "Fluxed", "FluxFiles", "Sensfunc", "*.log", "*.pdf", "*.png"))
    copied = manual_dir / source_pypeit.name
    destination = manual_dir / f"{source_pypeit.stem}_manual_{exposure}.pypeit"
    write_manual_column(copied, destination, exposure, manual_value(selections))
    copied.unlink()
    (manual_dir / f"manual_selection_{exposure}.json").write_text(
        json.dumps([asdict(item) for item in selections], indent=2) + "\n")
    return manual_dir, destination


def main() -> int:
    parser = argparse.ArgumentParser(description="Interactively position safe NGPS manual extractions.")
    parser.add_argument("spec2d", type=Path, help="Reduced PypeIt spec2d FITS file")
    parser.add_argument("--pypeit", type=Path, help="Source PypeIt file; auto-found when unambiguous")
    parser.add_argument("--fwhm", type=float, default=4.0, help="Initial FWHM in detector pixels")
    parser.add_argument("--max-select", type=int, default=3, help="Maximum selections: 1 to 3")
    parser.add_argument("--write-manual", action="store_true", help="Write copied manual setup without confirmation")
    parser.add_argument("--run-manual", action="store_true", help="Run copied setup; implies --write-manual")
    args = parser.parse_args()
    spec2d = args.spec2d.expanduser().resolve()
    if not spec2d.is_file():
        parser.error(f"File not found: {spec2d}")
    if not 1 <= args.max_select <= 3 or args.fwhm <= 0:
        parser.error("--max-select must be 1 to 3 and --fwhm must be positive")

    result = interactive_select(spec2d, args.fwhm, args.max_select)
    if result.decision == "automatic":
        print("Automatic PypeIt extraction accepted. Nothing was written.")
        return 0
    if result.decision != "manual":
        print("Extraction review cancelled. Nothing was written.")
        return 0
    selections = result.selections
    print("\nPypeIt manual value:\n  " + manual_value(selections))
    write = args.write_manual or args.run_manual or ask_yes_no("Write this to a copied manual PypeIt setup?")
    if not write:
        print("No files were written.")
        return 0
    if args.pypeit is not None:
        source = args.pypeit.expanduser().resolve()
    else:
        candidates = sorted(spec2d.parent.parent.glob("*.pypeit"))
        if len(candidates) != 1:
            raise RuntimeError("Use --pypeit because the setup does not have one unique .pypeit file.")
        source = candidates[0]
    manual_dir, manual_pypeit = create_manual_copy(source, exposure_from_spec2d(spec2d), selections)
    print(f"Manual setup: {manual_dir}\nManual PypeIt file: {manual_pypeit}")
    if args.run_manual or ask_yes_no("Run PypeIt on the copied manual setup now?"):
        return subprocess.run(["run_pypeit", manual_pypeit.name], cwd=manual_dir).returncode
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
