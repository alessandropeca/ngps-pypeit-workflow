#!/usr/bin/env python3
"""Telluric-correct completed NGPS R/I coadds with the pinned PypeIt."""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import re
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits


DEFAULT_CHANNELS = ("r", "i")
TELLURIC_GRID = "TellPCA_3000_26000_R10000.fits"
REVIEW_FIELDS = (
    "target", "channel", "setup", "source_coadd", "source_sha256",
    "corrected_spectrum", "telluric_model", "status", "reason",
    "min_transmission", "low_transmission_fraction", "qa_pdf", "qa_png",
)


@dataclass(frozen=True)
class Product:
    """One completed coadd eligible for telluric correction."""

    target: str
    channel: str
    setup: str
    source: Path


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_")


def configuration(setup: str) -> str:
    """Return the NGPS setup letter from a PypeIt setup name."""
    return setup.rsplit("_", 1)[-1].upper()


def source_sha256(path: Path) -> str:
    """Return a stable checksum so stale corrections cannot be selected later."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def coadd_path(root: Path, target: str, channel: str, setup: str) -> Path:
    stem = f"{safe_name(target)}_{channel}_{safe_name(setup)}"
    return root / "Coadds" / stem / f"{stem}_coadd.fits"


def output_paths(root: Path, product: Product) -> dict[str, Path]:
    """Name immutable-source telluric products separately from Coadds/."""
    directory = root / "Telluric" / safe_name(product.target)
    stem = product.source.stem
    qa_directory = root / "TelluricQA" / safe_name(product.target)
    return {
        "directory": directory,
        "corrected": directory / f"{stem}_tellcorr.fits",
        "model": directory / f"{stem}_tellmodel.fits",
        "parameters": directory / f"{stem}_telluric.par",
        "log": directory / f"{stem}_telluric.log",
        "qa_pdf": qa_directory / f"{stem}_telluric_qa.pdf",
        "qa_png": qa_directory / f"{stem}_telluric_qa.png",
    }


def read_coadd_review(root: Path) -> list[dict[str, str]]:
    path = root / "coadd_review.csv"
    if not path.is_file():
        raise ValueError(f"Coadd review file not found: {path}")
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    needed = {"target", "channel", "setup", "status"}
    if any(needed - set(row) for row in rows):
        raise ValueError(f"Invalid coadd review file: {path}")
    return rows


def selected_products(
    root: Path, target: str | None, channels: tuple[str, ...], setup_filter: str | None,
) -> tuple[list[Product], list[str]]:
    """Return completed, fluxed coadds and concise reasons for exclusions."""
    products: list[Product] = []
    skipped: list[str] = []
    for row in read_coadd_review(root):
        channel = row["channel"].casefold()
        if row["status"].casefold() != "coadded":
            continue
        if target is not None and row["target"].casefold() != target.casefold():
            continue
        if channel not in channels:
            continue
        if setup_filter is not None and configuration(row["setup"]) != setup_filter.upper():
            continue
        source = coadd_path(root, row["target"], channel, row["setup"])
        if not source.is_file():
            skipped.append(
                f"{row['target']} | {channel.upper()} | {row['setup']}: missing final coadd"
            )
            continue
        try:
            with fits.open(source, memmap=False) as hdul:
                fluxed = bool(hdul["SPECTRUM"].header.get("FLUXED", False))
        except (OSError, KeyError) as error:
            skipped.append(f"{row['target']} | {channel.upper()} | {row['setup']}: {error}")
            continue
        if not fluxed:
            skipped.append(
                f"{row['target']} | {channel.upper()} | {row['setup']}: coadd is not flux calibrated"
            )
            continue
        products.append(Product(row["target"], channel, row["setup"], source))
    return sorted(products, key=lambda item: (item.target.casefold(), item.setup, item.channel)), skipped


def read_spectrum(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Read one PypeIt coadd or telluric-corrected OneSpec product."""
    with fits.open(path, memmap=False) as hdul:
        data = hdul["SPECTRUM"].data
        names = set(data.dtype.names or ())
        required = {"wave", "flux", "ivar", "mask"}
        if required - names:
            raise ValueError(f"Missing spectrum columns in {path}")
        wave = np.asarray(data["wave"], dtype=float)
        flux = np.asarray(data["flux"], dtype=float)
        ivar = np.asarray(data["ivar"], dtype=float)
        mask = np.asarray(data["mask"], dtype=bool)
    return wave, flux, ivar, mask


def validate_output(corrected: Path, model: Path) -> tuple[float, float]:
    """Reject incomplete or non-physical PypeIt telluric products."""
    wave, flux, ivar, mask = read_spectrum(corrected)
    with fits.open(corrected, memmap=False) as hdul:
        data = hdul["SPECTRUM"].data
        if "telluric" not in (data.dtype.names or ()):
            raise ValueError("corrected spectrum has no telluric transmission column")
        transmission = np.asarray(data["telluric"], dtype=float)
    valid = mask & np.isfinite(wave) & np.isfinite(flux) & np.isfinite(ivar) & np.isfinite(transmission)
    if int(valid.sum()) < 10:
        raise ValueError("corrected spectrum has fewer than ten valid samples")
    transmission = transmission[valid]
    if np.any(transmission < 0) or np.any(transmission > 1.05):
        raise ValueError("telluric transmission is outside the physical range")
    with fits.open(model, memmap=False) as hdul:
        table = hdul["MODEL"].data
        names = set(table.dtype.names or ())
        if "SUCCESS" not in names or not bool(np.all(table["SUCCESS"])):
            raise ValueError("PypeIt did not mark the telluric fit successful")
    return float(np.min(transmission)), float(np.mean(transmission < 0.2))


def plot_limits(values: list[np.ndarray]) -> tuple[float, float]:
    finite = np.concatenate([item[np.isfinite(item)] for item in values])
    low, high = np.nanpercentile(finite, (0.5, 99.5))
    padding = 0.08 * (high - low) if high > low else max(abs(high) * 0.1, 1.0)
    return float(low - padding), float(high + padding)


def save_qa(source: Path, corrected: Path, output: dict[str, Path], product: Product) -> None:
    """Save a compact before/after and transmission QA plot."""
    source_wave, source_flux, _, source_mask = read_spectrum(source)
    wave, flux, _, mask = read_spectrum(corrected)
    with fits.open(corrected, memmap=False) as hdul:
        transmission = np.asarray(hdul["SPECTRUM"].data["telluric"], dtype=float)
    source_good = source_mask & np.isfinite(source_wave) & np.isfinite(source_flux)
    corrected_good = mask & np.isfinite(wave) & np.isfinite(flux)
    transmission_good = mask & np.isfinite(wave) & np.isfinite(transmission)
    figure, (spectrum_axis, transmission_axis) = plt.subplots(2, 1, figsize=(10.0, 5.6), sharex=True)
    figure.subplots_adjust(left=0.12, right=0.98, bottom=0.12, top=0.88, hspace=0.12)
    spectrum_axis.plot(source_wave[source_good], source_flux[source_good], color="0.45", lw=0.7, label="fluxed coadd")
    spectrum_axis.plot(wave[corrected_good], flux[corrected_good], color="#d62728", lw=0.75, label="telluric-corrected")
    spectrum_axis.axhline(0, color="0.70", lw=0.6)
    spectrum_axis.set_ylim(*plot_limits([source_flux[source_good], flux[corrected_good]]))
    spectrum_axis.set_ylabel(r"Flux ($10^{-17}$ erg s$^{-1}$ cm$^{-2}$ Å$^{-1}$)")
    spectrum_axis.legend(loc="upper right", fontsize=8)
    transmission_axis.plot(wave[transmission_good], transmission[transmission_good], color="#1f77b4", lw=0.8)
    transmission_axis.axhline(0.2, color="#d62728", ls="--", lw=0.8, label="low-transmission threshold")
    transmission_axis.set_ylim(-0.03, 1.05)
    transmission_axis.set_ylabel("Transmission")
    transmission_axis.set_xlabel("Wavelength (Å)")
    transmission_axis.legend(loc="lower left", fontsize=8)
    figure.suptitle(f"{product.target} | {product.channel.upper()} | {product.setup} | telluric QA", fontsize=13)
    output["qa_pdf"].parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output["qa_pdf"])
    figure.savefig(output["qa_png"], dpi=180)
    plt.close(figure)


def tellfit_command(source: Path, parameters: Path) -> list[str]:
    runner = Path(sys.executable).with_name("pypeit_tellfit")
    executable = str(runner) if runner.is_file() else "pypeit_tellfit"
    return [executable, str(source), "--objmodel", "poly", "--par_outfile", str(parameters), "-v", "1"]


def require_telluric_grid() -> None:
    """Require an explicit one-time install instead of an implicit network download."""
    try:
        from pypeit import dataPaths
        path = Path(dataPaths.telgrid.get_file_path(TELLURIC_GRID))
    except ImportError as error:
        raise ValueError("PypeIt is not installed in the active environment") from error
    if not path.is_file():
        raise ValueError(
            f"Telluric PCA model is missing. Install it once with: "
            f"pypeit_install_telluric {TELLURIC_GRID}"
        )


def run_product(root: Path, product: Product) -> dict[str, str]:
    """Run PypeIt in staging and atomically publish only validated results."""
    output = output_paths(root, product)
    output["directory"].mkdir(parents=True, exist_ok=True)
    checksum = source_sha256(product.source)
    with tempfile.TemporaryDirectory(prefix=".telluric_", dir=output["directory"]) as temporary:
        staging = Path(temporary)
        parameters = staging / output["parameters"].name
        command = tellfit_command(product.source, parameters)
        result = subprocess.run(command, cwd=staging, capture_output=True, text=True)
        log = result.stdout + result.stderr
        (staging / output["log"].name).write_text(log)
        expected_corrected = staging / output["corrected"].name
        expected_model = staging / output["model"].name
        if result.returncode != 0:
            return review_row(product, output, checksum, "failed", f"PypeIt returned status {result.returncode}")
        if not expected_corrected.is_file() or not expected_model.is_file():
            return review_row(product, output, checksum, "failed", "PypeIt did not write both telluric FITS products")
        try:
            minimum, low_fraction = validate_output(expected_corrected, expected_model)
        except (OSError, KeyError, ValueError) as error:
            return review_row(product, output, checksum, "failed", str(error))
        for key, staged in (
            ("corrected", expected_corrected), ("model", expected_model),
            ("parameters", parameters), ("log", staging / output["log"].name),
        ):
            os.replace(staged, output[key])
    save_qa(product.source, output["corrected"], output, product)
    return review_row(product, output, checksum, "completed", "PypeIt telluric fit completed", minimum, low_fraction)


def review_row(
    product: Product, output: dict[str, Path], checksum: str, status: str, reason: str,
    minimum: float | None = None, low_fraction: float | None = None,
) -> dict[str, str]:
    return {
        "target": product.target,
        "channel": product.channel.upper(),
        "setup": product.setup,
        "source_coadd": str(product.source),
        "source_sha256": checksum,
        "corrected_spectrum": str(output["corrected"]),
        "telluric_model": str(output["model"]),
        "status": status,
        "reason": reason,
        "min_transmission": "" if minimum is None else f"{minimum:.5f}",
        "low_transmission_fraction": "" if low_fraction is None else f"{low_fraction:.5f}",
        "qa_pdf": str(output["qa_pdf"]),
        "qa_png": str(output["qa_png"]),
    }


def load_review(path: Path) -> dict[tuple[str, str, str], dict[str, str]]:
    if not path.is_file():
        return {}
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if any(set(REVIEW_FIELDS) - set(row) for row in rows):
        raise ValueError(f"Invalid telluric review file: {path}")
    return {(row["target"].casefold(), row["channel"].casefold(), row["setup"]): row for row in rows}


def write_review(path: Path, rows: dict[tuple[str, str, str], dict[str, str]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=REVIEW_FIELDS)
        writer.writeheader()
        writer.writerows(sorted(rows.values(), key=lambda row: (row["target"].casefold(), row["setup"], row["channel"])))


def print_plan(products: list[Product], skipped: list[str], run: bool) -> None:
    print("\nTelluric correction plan")
    print("target                    channel  setup                 source coadd")
    for product in products:
        print(f"{product.target:<25} {product.channel.upper():<7} {product.setup:<21} {product.source}")
    if skipped:
        print("\nSkipped")
        for item in skipped:
            print(item)
    if not run:
        print("\nNo files were changed. Add --run to create separate Telluric products.")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Telluric-correct completed NGPS R/I coadds without changing the source coadds."
    )
    parser.add_argument("date", help="UT date, e.g. 20260623")
    selection = parser.add_mutually_exclusive_group(required=True)
    selection.add_argument("--target", help="Process one target")
    selection.add_argument("--all", action="store_true", help="Process every completed coadd")
    parser.add_argument("--channel", action="append", choices=("u", "g", "r", "i"), help="Channel to process. Repeat as needed. Default: R and I")
    parser.add_argument("--configuration", help="Setup letter for one target, e.g. B")
    parser.add_argument("--run", action="store_true", help="Run PypeIt and replace prior telluric products only after validation")
    args = parser.parse_args()
    if args.all and args.configuration:
        parser.error("--configuration applies only with --target")
    channels = tuple(dict.fromkeys(args.channel or DEFAULT_CHANNELS))
    root = Path(os.environ.get("NGPS_WORK_ROOT", Path.home() / "ngps_data" / "work")) / args.date
    try:
        products, skipped = selected_products(root, args.target, channels, args.configuration)
    except ValueError as error:
        parser.error(str(error))
    if not products:
        parser.error("No completed, flux-calibrated coadds match this selection")
    print_plan(products, skipped, args.run)
    if not args.run:
        return 0
    try:
        require_telluric_grid()
    except ValueError as error:
        parser.error(str(error))
    review_file = root / "telluric_review.csv"
    try:
        review = load_review(review_file)
    except ValueError as error:
        parser.error(str(error))
    completed = 0
    failed = 0
    for product in products:
        print(f"\nTelluric correction: {product.target} | {product.channel.upper()} | {product.setup}")
        row = run_product(root, product)
        review[(product.target.casefold(), product.channel, product.setup)] = row
        print(f"{row['status']}: {row['reason']}")
        if row["status"] == "completed":
            completed += 1
        else:
            failed += 1
    write_review(review_file, review)
    print(f"\nTelluric products completed: {completed}")
    print(f"Telluric products failed: {failed}")
    print(f"Review: {review_file}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
