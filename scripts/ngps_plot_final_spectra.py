#!/usr/bin/env python3
"""Plot the four separate final NGPS channel coadds for one target."""

from __future__ import annotations

import argparse
import csv
import os
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits


CHANNELS = ("u", "g", "r", "i")
COLOURS = {"u": "#9467bd", "g": "#2ca02c", "r": "#d62728", "i": "#ff7f0e"}


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_")


def completed_coadds(root: Path, target: str) -> dict[str, dict[str, Path]]:
    """Find completed coadds keyed by setup letter and channel."""
    coadds = root / "Coadds"
    name = safe_name(target)
    found: dict[str, dict[str, Path]] = {}
    for channel in CHANNELS:
        prefix = f"{name}_{channel}_p200_ngps_{channel}_"
        for directory in sorted(coadds.glob(f"{prefix}*")):
            if not directory.is_dir():
                continue
            configuration = directory.name.removeprefix(prefix).upper()
            filename = directory / f"{directory.name}_coadd.fits"
            if filename.is_file():
                found.setdefault(configuration, {})[channel] = filename
    return found


def complete_configurations(found: dict[str, dict[str, Path]]) -> dict[str, dict[str, Path]]:
    """Return only configurations with all four final channel coadds."""
    return {key: value for key, value in found.items() if set(value) == set(CHANNELS)}


def choose_configuration(
    found: dict[str, dict[str, Path]], requested: str | None,
) -> tuple[str, dict[str, Path]]:
    """Choose one complete U/G/R/I configuration, without guessing among several."""
    complete = complete_configurations(found)
    if requested is not None:
        configuration = requested.upper()
        if configuration not in complete:
            present = ", ".join(sorted(found.get(configuration, {}))) or "none"
            raise ValueError(
                f"Configuration {configuration} does not have all U/G/R/I coadds. "
                f"Available channels: {present}"
            )
        return configuration, complete[configuration]
    if len(complete) == 1:
        return next(iter(complete.items()))
    if not complete:
        details = "; ".join(
            f"{key}: {','.join(sorted(value)) or 'none'}" for key, value in sorted(found.items())
        ) or "none"
        raise ValueError(f"No complete U/G/R/I coadd set found for this target. Found: {details}")
    choices = ", ".join(sorted(complete))
    raise ValueError(f"More than one complete configuration is available: {choices}. Use --configuration.")


def completed_targets(root: Path) -> list[str]:
    """Read target names from the completed coadd-review table."""
    review = root / "coadd_review.csv"
    if not review.is_file():
        raise ValueError(f"Coadd review file not found: {review}")
    with review.open(newline="") as handle:
        return sorted({
            row["target"] for row in csv.DictReader(handle)
            if row.get("status", "").casefold() == "coadded"
        })


def read_spectrum(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Read valid flux-calibrated OneSpec samples from PypeIt's final coadd."""
    with fits.open(path, memmap=False) as hdul:
        data = hdul["SPECTRUM"].data
        fluxed = bool(hdul["SPECTRUM"].header.get("FLUXED", False))
        if not fluxed:
            raise ValueError(f"Coadd is not flux calibrated: {path}")
        wave = np.asarray(data["wave"], dtype=float)
        flux = np.asarray(data["flux"], dtype=float)
        mask = np.asarray(data["mask"], dtype=bool)
    good = mask & np.isfinite(wave) & np.isfinite(flux)
    if not np.any(good):
        raise ValueError(f"No valid samples in: {path}")
    return wave[good], flux[good]


def display_limits(fluxes: list[np.ndarray]) -> tuple[float, float]:
    """Choose robust linear limits so isolated endpoint spikes do not hide the spectrum."""
    values = np.concatenate(fluxes)
    low, high = np.nanpercentile(values, (1.0, 99.0))
    if not np.isfinite(low) or not np.isfinite(high) or low == high:
        low, high = np.nanmin(values), np.nanmax(values)
    padding = 0.08 * (high - low) if high > low else max(abs(high) * 0.1, 1.0)
    return low - padding, high + padding


def save_plot(
    root: Path, target: str, configuration: str, paths: dict[str, Path], show: bool,
) -> tuple[Path, Path]:
    """Save one complete U/G/R/I QA plot and optionally open its interactive window."""
    spectra = {channel: read_spectrum(paths[channel]) for channel in CHANNELS}

    figure, axes = plt.subplots(len(CHANNELS), 1, figsize=(10.0, 7.0))
    for axis, channel in zip(axes, CHANNELS):
        wave, flux = spectra[channel]
        axis.plot(wave, flux, lw=0.8, color=COLOURS[channel])
        axis.axhline(0, color="0.65", lw=0.7)
        axis.set_xlim(np.min(wave), np.max(wave))
        axis.set_ylim(*display_limits([flux]))
        axis.set_ylabel(r"Flux ($10^{-17}$ erg s$^{-1}$ cm$^{-2}$ Å$^{-1}$)")
        axis.set_title(f"{channel.upper()} coadd", loc="left", color=COLOURS[channel])
        axis.set_xlabel("Wavelength (Å)")
    figure.suptitle(
        f"{target} | final separate U/G/R/I coadds | configuration {configuration}",
        fontsize=14,
    )
    figure.tight_layout(rect=(0, 0, 1, 0.97))

    output_dir = root / "FinalQA" / safe_name(target)
    output_dir.mkdir(parents=True, exist_ok=True)
    stem = f"{safe_name(target)}_UGRI_{configuration}_coadds"
    pdf = output_dir / f"{stem}.pdf"
    png = output_dir / f"{stem}.png"
    figure.savefig(pdf)
    figure.savefig(png, dpi=180)
    if show:
        print("Opening final-spectrum plot. Close the window when you are finished zooming or panning.")
        plt.show()
    plt.close(figure)

    return pdf, png


def print_saved_plot(target: str, configuration: str, paths: dict[str, Path], pdf: Path, png: Path) -> None:
    """Print one concise plot result and its source coadds."""
    print(f"\nFinal U/G/R/I coadd plot saved: {target} | configuration {configuration}")
    for channel in CHANNELS:
        print(f"{channel.upper()}: {paths[channel]}")
    print(f"PDF: {pdf}")
    print(f"PNG: {png}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Save U/G/R/I flux-versus-wavelength plots from final NGPS coadds."
    )
    parser.add_argument("date", help="UT date, e.g. 20260623")
    selection = parser.add_mutually_exclusive_group(required=True)
    selection.add_argument("--target", help="Save and open one target plot, e.g. MGC+04-48-002")
    selection.add_argument("--all", action="store_true", help="Save plots for every complete target/configuration")
    parser.add_argument(
        "--configuration", help="NGPS setup letter for one target when multiple complete sets exist, e.g. B",
    )
    args = parser.parse_args()
    if args.all and args.configuration:
        parser.error("--configuration applies only with --target")

    root = Path(os.environ.get("NGPS_WORK_ROOT", Path.home() / "ngps_data" / "work")) / args.date
    if args.target:
        try:
            configuration, paths = choose_configuration(
                completed_coadds(root, args.target), args.configuration
            )
            pdf, png = save_plot(root, args.target, configuration, paths, show=True)
        except ValueError as error:
            parser.error(str(error))
        print_saved_plot(args.target, configuration, paths, pdf, png)
        return 0

    try:
        targets = completed_targets(root)
    except ValueError as error:
        parser.error(str(error))
    completed = 0
    failed = 0
    for target in targets:
        for configuration, paths in sorted(complete_configurations(completed_coadds(root, target)).items()):
            try:
                pdf, png = save_plot(root, target, configuration, paths, show=False)
            except ValueError as error:
                print(f"Failed to plot {target} | configuration {configuration}: {error}")
                failed += 1
                continue
            print_saved_plot(target, configuration, paths, pdf, png)
            completed += 1
    print(f"\nFinal U/G/R/I plots saved: {completed}")
    if failed:
        print(f"Final U/G/R/I plots failed: {failed}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
