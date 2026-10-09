"""Fixed slicer-relative boxcar extraction, independent of source tracing.

The input is PypeIt's already reduced, sky-subtracted detector image. No
object finder, source trace, profile fit, or optimal extraction is used here.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from io import BytesIO
import json
from pathlib import Path
import shutil
import tempfile

import numpy as np
from astropy.io import fits
from matplotlib.widgets import Button, CheckButtons, TextBox
from pypeit.spec2dobj import Spec2DObj
from pypeit.specobj import SpecObj
from pypeit.specobjs import SpecObjs
from pypeit.spectrographs.util import load_spectrograph
from pypeit.core import parse


@dataclass
class FixedAperture:
    offset: float
    half_width: float
    channels: tuple[str, ...]
    reference_channel: str | None = None
    pixel_scales: dict[str, float] = field(default_factory=dict)
    qa_pdf: bytes | None = field(default=None, repr=False)

    def validate(self) -> None:
        if not np.isfinite(self.offset):
            raise ValueError("The aperture centre must be finite.")
        if not np.isfinite(self.half_width) or self.half_width <= 0:
            raise ValueError("Half-width must be a positive, finite number of pixels.")
        if not self.channels or any(c not in ("u", "g", "r", "i") for c in self.channels):
            raise ValueError("Select at least one available channel.")
        if self.pixel_scales:
            required = set(self.channels) | {self.reference_channel}
            if any(c not in self.pixel_scales or not np.isfinite(self.pixel_scales[c])
                   or self.pixel_scales[c] <= 0 for c in required):
                raise ValueError("Valid spatial pixel scales are required for the selected channels.")

    def parameters(self, channel: str) -> tuple[float, float]:
        ratio = (self.pixel_scales[self.reference_channel] / self.pixel_scales[channel]
                 if self.pixel_scales else 1.)
        return self.offset * ratio, self.half_width * ratio


@dataclass
class ApertureSpectrum:
    wave: np.ndarray
    counts: np.ndarray
    ivar: np.ndarray
    mask: np.ndarray
    npix: np.ndarray
    fraction: np.ndarray
    sky: np.ndarray
    centre: np.ndarray


def sum_fixed_band(image, ivar, waveimg, good_pixels, sky, left, right,
                   offset: float, half_width: float) -> ApertureSpectrum:
    """Sum fractional pixel overlaps and propagate independent-pixel variance.

    Limits are constant in rectified, geometric-slicer coordinates. Detector
    curvature is described only by the slit edges, never by a source trace.
    A row is flagged bad if any of its aperture is outside the slit/detector
    or masked. Partial sums are retained for diagnostics but get zero IVAR.
    """
    FixedAperture(offset, half_width, ("r",)).validate()
    image, ivar, waveimg, sky = [np.asarray(a, dtype=float) for a in (image, ivar, waveimg, sky)]
    good_pixels = np.asarray(good_pixels, dtype=bool)
    if image.ndim != 2 or any(a.shape != image.shape for a in (ivar, waveimg, sky, good_pixels)):
        raise ValueError("Image, wavelength, variance and mask shapes must agree.")
    left, right = np.asarray(left, dtype=float), np.asarray(right, dtype=float)
    if left.shape != (image.shape[0],) or right.shape != left.shape:
        raise ValueError("Slit edges must cover every spectral row.")
    if not np.all(np.isfinite(left) & np.isfinite(right) & (right > left)):
        raise ValueError("The slicer edges are invalid.")
    centre = (left + right) / 2 + offset
    low, high = centre - half_width, centre + half_width
    x = np.arange(image.shape[1], dtype=float)[None, :]
    # Pixel centres are integers. Boundaries are at x +/- 0.5.
    weights = np.clip(np.minimum(x + .5, high[:, None]) - np.maximum(x - .5, low[:, None]), 0, 1)
    inside_slit = (x >= left[:, None]) & (x <= right[:, None])
    valid = (good_pixels & inside_slit & np.isfinite(image) & np.isfinite(sky)
             & np.isfinite(ivar) & (ivar > 0) & np.isfinite(waveimg) & (waveimg > 0))
    used = weights * valid
    npix = used.sum(axis=1)
    fraction = npix / (2 * half_width)
    variance_pixels = np.divide(1., ivar, out=np.zeros_like(ivar), where=valid)
    variance = (used ** 2 * variance_pixels).sum(axis=1)
    counts = (used * np.where(valid, image, 0)).sum(axis=1)
    sky_counts = (used * np.where(valid, sky, 0)).sum(axis=1)
    # Use this aperture's wavelength map, not another object's 1D spectrum.
    wavelength = np.divide(
        (used * np.where(valid, waveimg, 0)).sum(axis=1), npix,
        out=np.zeros_like(npix), where=npix > 0,
    )
    mask = (np.isclose(fraction, 1., rtol=0, atol=1e-7)
            & np.isfinite(variance) & (variance > 0)
            & np.isfinite(wavelength) & (wavelength > 0))
    out_ivar = np.divide(1., variance, out=np.zeros_like(variance), where=mask)
    return ApertureSpectrum(wavelength, counts, out_ivar, mask, npix, fraction, sky_counts, centre)


def default_half_width(spec2d: Path) -> float:
    """Use the central slicer's automatic boxcar radius, not its FWHM."""
    spec1d = spec2d.with_name(spec2d.name.replace("spec2d_", "spec1d_", 1))
    if spec1d.is_file():
        with fits.open(spec1d, memmap=False) as hdul:
            objects = sorted((h for h in hdul[1:] if h.name.startswith("SPAT")),
                             key=lambda h: h.header.get("SLITID", 0))
            if objects:
                value = objects[len(objects) // 2].header.get("BOX_R_PIX")
                if value is not None and np.isfinite(float(value)) and float(value) > 0:
                    return float(value)
    return 4.0


def spatial_pixel_scale(spec2d: Path) -> float:
    reduced = Spec2DObj.from_file(str(spec2d), "DET01")
    _, spatial_binning = parse.parse_binning(reduced.detector.binning)
    value = float(reduced.detector.platescale * spatial_binning)
    if not np.isfinite(value) or value <= 0:
        raise ValueError("Missing spatial pixel scale. Cannot link channel positions safely.")
    return value


def extract_fixed_slicers(spec2d: Path, offset: float, half_width: float) -> SpecObjs:
    """Create genuine BOX spectra for all slicers without copying any trace."""
    reduced = Spec2DObj.from_file(str(spec2d), "DET01")
    spec1d = spec2d.with_name(spec2d.name.replace("spec2d_", "spec1d_", 1))
    header = fits.getheader(spec1d if spec1d.is_file() else spec2d).copy()
    spectrograph = load_spectrograph(header["PYP_SPEC"])
    objects = SpecObjs(header=header)
    image = reduced.sciimg - reduced.skymodel
    good_pixels = reduced.bpmmask.mask == 0
    slits = reduced.slits
    # Match the geometric centre used in the dashboard's rectified image.
    order = np.argsort(np.nanmedian((slits.left_init + slits.right_init) / 2, axis=0))
    bad_slits = slits.bitmask.flagged(slits.mask, and_not=slits.bitmask.exclude_for_reducing)
    for index in order:
        left, right = slits.left_init[:, index], slits.right_init[:, index]
        spectrum = sum_fixed_band(image, reduced.ivarmodel, reduced.waveimg,
                                  good_pixels, reduced.skymodel, left, right, offset, half_width)
        if bad_slits[index]:
            spectrum.mask[:] = False
            spectrum.ivar[:] = 0
        obj = SpecObj(spectrograph.pypeline, "DET01", OBJTYPE="science",
                      SLITID=int(slits.spat_id[index]), PYP_SPEC=header["PYP_SPEC"])
        obj.DETECTOR = reduced.detector
        obj.TRACE_SPAT = spectrum.centre  # Geometric band centre, not a fitted source trace.
        obj.trace_spec = np.arange(len(spectrum.wave))
        middle = len(spectrum.wave) // 2
        obj.SPAT_PIXPOS = float(spectrum.centre[middle])
        obj.SPAT_PIXPOS_ID = int(np.rint(obj.SPAT_PIXPOS))
        obj.SPAT_FRACPOS = float((spectrum.centre[middle] - left[middle]) / (right[middle] - left[middle]))
        obj.OBJID = 1
        obj.hand_extract_flag = True
        obj.BOX_R_PIX = float(half_width)
        obj.BOX_R_ASEC = float(half_width * obj.platescale)
        obj.BOX_WAVE, obj.BOX_COUNTS = spectrum.wave, spectrum.counts
        obj.BOX_COUNTS_IVAR, obj.BOX_MASK = spectrum.ivar, spectrum.mask
        obj.BOX_COUNTS_SIG = np.sqrt(np.divide(1., spectrum.ivar, out=np.zeros_like(spectrum.ivar), where=spectrum.mask))
        obj.BOX_COUNTS_SKY = spectrum.sky
        obj.BOX_NPIX, obj.BOX_FRAC_USE = spectrum.npix, spectrum.fraction
        obj.S2N = float(np.median(spectrum.counts[spectrum.mask] * np.sqrt(spectrum.ivar[spectrum.mask]))) if np.any(spectrum.mask) else 0.
        obj.RA, obj.DEC = header.get("RA"), header.get("DEC")
        obj.VEL_TYPE, obj.VEL_CORR = reduced.vel_type, reduced.vel_corr
        # waveimg already includes PypeIt's global flexure and velocity correction.
        # Do not apply the old bright source's local spectral-flexure correction.
        obj.set_name()
        objects.add_sobj(obj)
    if not objects.nobj or not any(np.any(obj.BOX_MASK) for obj in objects):
        raise ValueError("The aperture contains no fully valid spectral rows in any slicer.")
    objects.header["NGPSMODE"] = "FIXED"
    objects.header["FIXOFF"] = float(offset)
    objects.header["FIXHALF"] = float(half_width)
    objects.header["FIXCOORD"] = "SLICER_RELATIVE"
    objects.header["HISTORY"] = "NGPS fixed aperture: geometric slicer centre, no source trace or profile fit."
    return objects


def install_fixed_exposure(root: Path, target: str, exposure: str,
                           paths: dict[str, tuple[Path, Path]], aperture: FixedAperture,
                           qa_path: Path) -> None:
    """Stage and validate all selected channels, then replace with rollback.

    ``paths`` maps channel to (input spec2d, baseline setup directory).
    Existing spec2d data are unchanged. Only accepted channels' spec1d/text
    products and stale Fluxed copies are replaced/invalidated.
    """
    aperture.validate()
    if any(c not in paths for c in aperture.channels):
        raise ValueError("An unavailable channel was selected.")
    with tempfile.TemporaryDirectory(prefix=".ngps_fixed_", dir=root) as temporary:
        staging = Path(temporary)
        replacements: dict[Path, Path] = {}
        removals: set[Path] = set()
        provenance = {"mode": "fixed", "target": target, "exposure": exposure,
                      "offset_pixels": aperture.offset, "half_width_pixels": aperture.half_width,
                      "channels": list(aperture.channels), "coordinate": "geometric slicer-relative",
                      "reference_channel": aperture.reference_channel,
                      "pixel_scales_arcsec": aperture.pixel_scales,
                      "channel_parameters_pixels": {c: dict(zip(("offset", "half_width"), aperture.parameters(c))) for c in aperture.channels},
                      "mask_policy": "any missing aperture pixel masks the 1D row",
                      "variance": "sum of squared fractional weights times pixel variance"}
        for channel in aperture.channels:
            spec2d, setup = paths[channel]
            offset, half_width = aperture.parameters(channel)
            products = extract_fixed_slicers(spec2d, offset, half_width)
            destination = setup / "Science" / spec2d.name.replace("spec2d_", "spec1d_", 1)
            channel_stage = staging / channel
            channel_stage.mkdir()
            staged = channel_stage / destination.name
            products.header["FIXCHANS"] = ",".join(aperture.channels)
            products.write_to_fits(products.header, str(staged))
            # Test the actual downstream reader before touching existing spectra.
            checked = SpecObjs.from_fitsfile(str(staged))
            if len(checked) != len(products) or any(o.OPT_COUNTS is not None for o in checked):
                raise ValueError("Fixed-aperture FITS validation failed.")
            replacements[destination] = staged
            summary = staged.with_suffix(".txt")
            summary.write_text(
                f"Fixed aperture: centre offset {offset:g} px, half-width {half_width:g} px\n"
                "slit | object | valid_rows | total_rows\n" + "".join(
                    f"{o.SLITID} | {o.NAME} | {np.count_nonzero(o.BOX_MASK)} | {len(o.BOX_MASK)}\n" for o in checked)
            )
            replacements[destination.with_suffix(".txt")] = summary
            removals.update((setup / "Fluxed").glob(f"spec1d_*_{exposure}-*.fits"))
        manifest = qa_path.with_name(f"ngps_fixed_aperture_{exposure}.json")
        staged_manifest = staging / "selection.json"
        staged_manifest.write_text(json.dumps(provenance, indent=2) + "\n")
        replacements[manifest] = staged_manifest
        if aperture.qa_pdf is not None:
            staged_pdf = staging / "review.pdf"
            staged_pdf.write_bytes(aperture.qa_pdf)
            replacements[qa_path] = staged_pdf
        backup = staging / "backup"
        backup.mkdir()
        old: dict[Path, Path | None] = {}
        for index, destination in enumerate(set(replacements) | removals):
            saved = backup / str(index) if destination.is_file() else None
            if saved is not None:
                shutil.copy2(destination, saved)
            old[destination] = saved
        try:
            for destination, source in replacements.items():
                destination.parent.mkdir(parents=True, exist_ok=True)
                source.replace(destination)
            for destination in removals:
                destination.unlink()
        except BaseException:
            for destination, saved in old.items():
                if saved is None:
                    destination.unlink(missing_ok=True)
                else:
                    saved.replace(destination)
            raise


def pdf_bytes(figure) -> bytes:
    stream = BytesIO()
    figure.savefig(stream, format="pdf")
    return stream.getvalue()


class FixedApertureDialog:
    """A popup panel on the same canvas, avoiding a second native GUI loop."""

    def __init__(self, figure, offset, half_width, channels, on_extract, on_close):
        self.figure, self.on_extract, self.on_close = figure, on_extract, on_close
        self.offset, self.channels = offset, tuple(channels)
        self.axes = []
        def panel(rect):
            axis = figure.add_axes(rect, zorder=100)
            self.axes.append(axis)
            return axis
        background = panel((.26, .22, .48, .53))
        background.set_facecolor("#F7F9FC")
        background.set_xticks([])
        background.set_yticks([])
        background.text(.5, .92, "Fixed aperture", ha="center", fontsize=14, fontweight="bold")
        background.text(.10, .80, f"Center: {offset:.2f} px", fontsize=11)
        background.text(.10, .53, "Apply to:", fontsize=10)
        background.text(.10, .10, "Extract previews the band. Accept fixed aperture saves it.", fontsize=9)
        self.error = background.text(.10, .02, "", fontsize=9, color="#A32020")
        self.width = TextBox(panel((.42, .56, .15, .045)), "Half-width (px): ", initial=f"{half_width:g}")
        self.checks = CheckButtons(panel((.42, .32, .15, .17)), [c.upper() for c in channels], [True] * len(channels))
        self.extract = Button(panel((.60, .45, .10, .05)), "Extract", color="#D7F2DF")
        self.cancel = Button(panel((.60, .37, .10, .05)), "Cancel", color="#FFD9D9")
        self.extract.on_clicked(self.submit)
        self.cancel.on_clicked(lambda event: self.close())
        figure.canvas.draw_idle()

    def submit(self, event=None):
        try:
            selection = FixedAperture(float(self.offset), float(self.width.text),
                                     tuple(c for c, enabled in zip(self.channels, self.checks.get_status()) if enabled))
            selection.validate()
            self.on_extract(selection)
        except (ValueError, OSError, RuntimeError) as error:
            self.error.set_text(str(error))
            self.figure.canvas.draw_idle()
            return
        self.close()

    def close(self):
        for axis in self.axes:
            axis.remove()
        self.axes.clear()
        self.on_close()
        self.figure.canvas.draw_idle()
