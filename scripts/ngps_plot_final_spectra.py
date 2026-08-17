#!/usr/bin/env python3
"""Plot the four separate final NGPS channel coadds for one target."""

from __future__ import annotations

import argparse
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


def choose_configuration(
    found: dict[str, dict[str, Path]], requested: str | None,
) -> tuple[str, dict[str, Path]]:
    """Choose one complete U/G/R/I configuration, without guessing among several."""
    complete = {key: value for key, value in found.items() if set(value) == set(CHANNELS)}
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


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Save one U/G/R/I flux-versus-wavelength plot from final NGPS coadds."
    )
    parser.add_argument("date", help="UT date, e.g. 20260623")
    parser.add_argument("--target", required=True, help="Target name, e.g. MGC+04-48-002")
    parser.add_argument(
        "--configuration", help="NGPS setup letter when multiple complete channel sets exist, e.g. B",
    )
    args = parser.parse_args()

    root = Path(os.environ.get("NGPS_WORK_ROOT", Path.home() / "ngps_data" / "work")) / args.date
    try:
        configuration, paths = choose_configuration(
            completed_coadds(root, args.target), args.configuration
        )
        spectra = {channel: read_spectrum(paths[channel]) for channel in CHANNELS}
    except ValueError as error:
        parser.error(str(error))

    figure, axes = plt.subplots(len(CHANNELS), 1, figsize=(12, 10.5))
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
        f"{args.target} | final separate U/G/R/I coadds | configuration {configuration}",
        fontsize=14,
    )
    figure.tight_layout(rect=(0, 0, 1, 0.97))

    output_dir = root / "FinalQA" / safe_name(args.target)
    output_dir.mkdir(parents=True, exist_ok=True)
    stem = f"{safe_name(args.target)}_UGRI_{configuration}_coadds"
    pdf = output_dir / f"{stem}.pdf"
    png = output_dir / f"{stem}.png"
    figure.savefig(pdf)
    figure.savefig(png, dpi=180)
    plt.close(figure)

    print("Final U/G/R/I coadd plot saved")
    for channel in CHANNELS:
        print(f"{channel.upper()}: {paths[channel]}")
    print(f"PDF: {pdf}")
    print(f"PNG: {png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
