#!/usr/bin/env python3
"""Plot the four separate final NGPS channel coadds for one target."""

from __future__ import annotations

import argparse
import csv
import hashlib
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


def source_sha256(path: Path) -> str:
    """Return a stable checksum used to reject stale telluric products."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


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


def completed_tellurics(root: Path, target: str) -> dict[str, dict[str, Path]]:
    """Find validated telluric products whose source coadds have not changed."""
    review = root / "telluric_review.csv"
    if not review.is_file():
        return {}
    found: dict[str, dict[str, Path]] = {}
    with review.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        if row.get("target", "").casefold() != target.casefold():
            continue
        if row.get("status", "").casefold() != "completed":
            continue
        channel = row.get("channel", "").casefold()
        source = Path(row.get("source_coadd", ""))
        corrected = Path(row.get("corrected_spectrum", ""))
        if channel not in CHANNELS or not source.is_file() or not corrected.is_file():
            continue
        if not row.get("source_sha256") or source_sha256(source) != row["source_sha256"]:
            continue
        setup = row.get("setup", "")
        if not setup:
            continue
        found.setdefault(setup.rsplit("_", 1)[-1].upper(), {})[channel] = corrected
    return found


def final_products(root: Path, target: str) -> dict[str, dict[str, Path]]:
    """Use validated telluric products for R/I while retaining the source coadds elsewhere."""
    found = {configuration: dict(paths) for configuration, paths in completed_coadds(root, target).items()}
    for configuration, corrected in completed_tellurics(root, target).items():
        if configuration in found:
            found[configuration].update(corrected)
    return found


def choose_configuration(
    found: dict[str, dict[str, Path]], requested: str | None,
) -> tuple[str, dict[str, Path]]:
    """Choose one configuration that has at least one final channel coadd."""
    if requested is not None:
        configuration = requested.upper()
        if configuration not in found:
            present = ", ".join(sorted(found.get(configuration, {}))) or "none"
            raise ValueError(
                f"Configuration {configuration} has no completed coadd. "
                f"Available channels: {present}"
            )
        return configuration, found[configuration]
    if len(found) == 1:
        return next(iter(found.items()))
    if not found:
        details = ", ".join(
            f"{key}: {','.join(sorted(value)) or 'none'}" for key, value in sorted(found.items())
        ) or "none"
        raise ValueError(f"No final channel coadds found for this target. Found: {details}")
    choices = ", ".join(sorted(found))
    raise ValueError(f"More than one configuration is available: {choices}. Use --configuration.")


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
    """Return full linear limits for the samples shown in the plot."""
    values = np.concatenate(fluxes)
    low, high = np.nanmin(values), np.nanmax(values)
    padding = 0.08 * (high - low) if high > low else max(abs(high) * 0.1, 1.0)
    return low - padding, high + padding


def display_spectrum(
    wave: np.ndarray, flux: np.ndarray, channel: str, no_ug_edges: bool,
    manual_ranges: dict[str, tuple[float, float]],
) -> tuple[np.ndarray, np.ndarray, int]:
    """Apply optional wavelength display ranges without changing FITS data."""
    if channel in manual_ranges:
        minimum, maximum = manual_ranges[channel]
        keep = (wave >= minimum) & (wave <= maximum)
    elif no_ug_edges and channel == "u":
        minimum, maximum = 3100.0, 4350.0
        keep = (wave >= 3100.0) & (wave <= 4350.0)
    elif no_ug_edges and channel == "g":
        minimum, maximum = 4280.0, float(np.max(wave))
        keep = wave >= 4280.0
    else:
        return wave, flux, 0
    if not np.any(keep):
        raise ValueError(
            f"{channel.upper()} display range contains no samples: "
            f"{minimum:.0f}-{maximum:.0f} A"
        )
    return wave[keep], flux[keep], int((~keep).sum())


def save_plot(
    root: Path, target: str, configuration: str, paths: dict[str, Path], show: bool,
    no_ug_edges: bool, manual_ranges: dict[str, tuple[float, float]],
) -> tuple[Path, Path]:
    """Save available channel coadds, using an empty panel for an unavailable channel."""
    spectra: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for channel, path in paths.items():
        wave, flux = read_spectrum(path)
        wave, flux, removed = display_spectrum(
            wave, flux, channel, no_ug_edges, manual_ranges,
        )
        if removed:
            print(f"{channel.upper()} {target} {configuration}: omitted {removed} sample(s) outside the selected display range")
        spectra[channel] = wave, flux

    figure, axes = plt.subplots(len(CHANNELS), 1, figsize=(10.0, 7.0))
    figure.subplots_adjust(left=0.12, right=0.98, bottom=0.10, top=0.90, hspace=0.24)
    for axis, channel in zip(axes, CHANNELS):
        label = f"{channel.upper()} {'telluric-corrected' if path_is_telluric(paths.get(channel)) else 'coadd'}"
        label_colour = COLOURS[channel]
        axis.text(
            0.015, 0.88, label, transform=axis.transAxes, ha="left", va="top",
            color=label_colour, fontsize=10, fontweight="bold",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78, "pad": 1.8},
        )
        if channel not in spectra:
            axis.set_facecolor("white")
            axis.set_xlim(0, 1)
            axis.set_ylim(0, 1)
            axis.set_xticks([])
            axis.set_yticks([])
            axis.text(
                0.5, 0.5, "No flux-calibrated coadd", transform=axis.transAxes,
                ha="center", va="center", color="0.45",
            )
            continue
        wave, flux = spectra[channel]
        axis.plot(wave, flux, lw=0.8, color=COLOURS[channel])
        axis.axhline(0, color="0.65", lw=0.7)
        axis.set_xlim(np.min(wave), np.max(wave))
        axis.set_ylim(*display_limits([flux]))
    axes[-1].set_xlabel("Wavelength (Å)")
    figure.text(
        0.04, 0.5, r"Flux ($10^{-17}$ erg s$^{-1}$ cm$^{-2}$ Å$^{-1}$)",
        ha="center", va="center", rotation="vertical",
    )
    figure.suptitle(
        f"{target} | final available U/G/R/I spectra | configuration {configuration}",
        fontsize=14,
    )

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


def path_is_telluric(path: Path | None) -> bool:
    return path is not None and path.name.endswith("_tellcorr.fits")


def print_saved_plot(target: str, configuration: str, paths: dict[str, Path], pdf: Path, png: Path) -> None:
    """Print one concise plot result and its source coadds."""
    print(f"\nFinal U/G/R/I coadd plot saved: {target} | configuration {configuration}")
    for channel in CHANNELS:
        print(f"{channel.upper()}: {paths.get(channel, 'unavailable')}")
    print(f"PDF: {pdf}")
    print(f"PNG: {png}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Save U/G/R/I flux-versus-wavelength plots from final NGPS coadds."
    )
    parser.add_argument("date", help="UT date, e.g. 20260623")
    selection = parser.add_mutually_exclusive_group(required=True)
    selection.add_argument("--target", help="Save and open one target plot, e.g. MGC+04-48-002")
    selection.add_argument("--all", action="store_true", help="Save plots for every target/configuration with a completed coadd")
    parser.add_argument(
        "--configuration", help="NGPS setup letter for one target when multiple configurations exist, e.g. B",
    )
    parser.add_argument(
        "--noUGedges", "--no-ug-edges", dest="no_ug_edges", action="store_true",
        help="Plot U from 3100 to 4350 A and G from 4280 A onward without changing FITS data.",
    )
    parser.add_argument(
        "--manual", action="store_true",
        help="Use one or more manual wavelength display ranges supplied with --U, --G, --R, or --I.",
    )
    for channel in CHANNELS:
        parser.add_argument(
            f"--{channel.upper()}", dest=f"{channel}_range", nargs=2, type=float,
            metavar=("MIN", "MAX"), help=f"Manual {channel.upper()} display range in Angstrom.",
        )
    args = parser.parse_args()
    if args.all and args.configuration:
        parser.error("--configuration applies only with --target")
    manual_ranges = {
        channel: tuple(getattr(args, f"{channel}_range"))
        for channel in CHANNELS
        if getattr(args, f"{channel}_range") is not None
    }
    if args.manual and not manual_ranges:
        parser.error("--manual requires at least one of --U, --G, --R, or --I")
    if manual_ranges and not args.manual:
        parser.error("Use --manual when supplying --U, --G, --R, or --I")
    if args.manual and args.no_ug_edges:
        parser.error("Choose either --manual or --noUGedges")
    for channel, (minimum, maximum) in manual_ranges.items():
        if minimum >= maximum:
            parser.error(f"{channel.upper()} range must have MIN < MAX")

    root = Path(os.environ.get("NGPS_WORK_ROOT", Path.home() / "ngps_data" / "work")) / args.date
    if args.target:
        try:
            configuration, paths = choose_configuration(
                final_products(root, args.target), args.configuration
            )
            pdf, png = save_plot(
                root, args.target, configuration, paths, show=True,
                no_ug_edges=args.no_ug_edges, manual_ranges=manual_ranges,
            )
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
        for configuration, paths in sorted(final_products(root, target).items()):
            try:
                pdf, png = save_plot(
                    root, target, configuration, paths, show=False,
                    no_ug_edges=args.no_ug_edges, manual_ranges=manual_ranges,
                )
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
