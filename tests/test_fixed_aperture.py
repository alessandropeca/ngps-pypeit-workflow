"""Numerical, FITS, UI and downstream regressions for fixed apertures."""

import importlib
from pathlib import Path
import sys
from unittest.mock import patch

import numpy as np
import pytest
from astropy.io import fits
from astropy.table import Table
from matplotlib.backend_bases import MouseEvent
from pypeit.images.detector_container import DetectorContainer
from pypeit.images.imagebitmask import ImageBitMaskArray
from pypeit.slittrace import SlitTraceSet
from pypeit.spec2dobj import Spec2DObj
from pypeit.specobjs import SpecObjs

sys.path.insert(0, str(Path(__file__).parents[1] / "scripts"))
fixed = importlib.import_module("ngps_fixed_aperture")
review = importlib.import_module("ngps_manual_target_extractions")
coadd = importlib.import_module("ngps_interactive_coadd")


def arrays():
    image = np.full((8, 30), 10.)
    ivar = np.full_like(image, .25)
    wave = 5000 + np.arange(8)[:, None] * 2 + np.arange(30)[None, :] * .01
    return image, ivar, wave, np.ones_like(image, dtype=bool), np.full_like(image, 2.), np.full(8, 2.), np.full(8, 26.)


def test_fixed_band_fractional_counts_variance_and_own_wavelength():
    result = fixed.sum_fixed_band(*arrays(), offset=-3.2, half_width=2.)
    # Centre=10.8. Edges=8.8,12.8. Weights are .7,1,1,1,.3.
    np.testing.assert_allclose(result.counts, 40.)
    np.testing.assert_allclose(result.ivar, 1 / (4 * (.7**2 + 3 + .3**2)))
    np.testing.assert_allclose(result.wave, 5000 + np.arange(8) * 2 + 10.8 * .01)
    np.testing.assert_allclose(result.centre, 10.8)
    assert result.mask.all()


def test_fixed_band_ignores_bright_trace_and_follows_only_geometric_edges():
    values = list(arrays())
    curve = np.linspace(-.3, .3, 8)
    values[5] += curve
    values[6] += curve
    # A moving bright object, outside our chosen aperture, must have no effect.
    values[0][np.arange(8), np.arange(8) + 19] = 1e8
    result = fixed.sum_fixed_band(*values, offset=-3., half_width=2.)
    np.testing.assert_allclose(result.centre - (values[5] + values[6]) / 2, -3.)
    np.testing.assert_allclose(result.counts, 40.)


def test_fixed_band_propagates_bad_pixels_nans_zero_ivar_and_edges():
    values = list(arrays())
    values[3][0, 11] = False
    values[0][1, 11] = np.nan
    values[1][2, 11] = 0
    values[2][3, 11] = np.nan
    result = fixed.sum_fixed_band(*values, offset=-3., half_width=2.)
    assert not result.mask[:4].any()
    assert result.mask[4:].all()
    assert (result.ivar[:4] == 0).all()
    assert np.isfinite(result.counts).all()
    assert (result.fraction[:4] < 1).all()
    truncated = fixed.sum_fixed_band(*arrays(), offset=-12., half_width=4.)
    assert not truncated.mask.any()
    assert (truncated.ivar == 0).all()


@pytest.mark.parametrize("width", [0., -1., np.nan, np.inf])
def test_fixed_band_rejects_invalid_widths(width):
    with pytest.raises(ValueError):
        fixed.sum_fixed_band(*arrays(), offset=0., half_width=width)


def test_fixed_link_uses_spatial_pixel_scale():
    aperture = fixed.FixedAperture(-7.2, 4., ("u", "r"), "u", {"u": .2, "r": .4})
    aperture.validate()
    assert aperture.parameters("u") == (-7.2, 4.)
    assert aperture.parameters("r") == (-3.6, 2.)


@pytest.fixture
def synthetic_night(tmp_path):
    root = tmp_path / "20260623"
    frames = {}
    for channel in "ugri":
        setup = root / f"manual_setup_{channel}" / f"p200_ngps_{channel}_B"
        science = setup / "Science"
        science.mkdir(parents=True)
        path = science / f"spec2d_ngps_260623_0121-target_NGPS_{channel}_time.fits"
        shape = (24, 90)
        left = np.tile([2., 32., 62.], (shape[0], 1))
        slits = SlitTraceSet(left, left + 22., "MultiSlit", nspat=shape[1], PYP_SPEC=f"p200_ngps_{channel}")
        flex = Table({"spat_id": slits.spat_id, "sci_spec_flexure": np.zeros(3)})
        detector = DetectorContainer(dataext=0, specaxis=0, specflip=False, spatflip=False,
                                     platescale=.2, saturation=1e5, mincounts=-1e9,
                                     nonlinear=.9, numamplifiers=1, gain=np.array([1.]),
                                     ronoise=np.array([3.]), det=1, binning="1,1")
        wave = np.broadcast_to(5000 + np.arange(shape[0])[:, None] * 2., shape).copy()
        reduced = Spec2DObj(sciimg=np.full(shape, 12.), ivarraw=np.full(shape, .25),
                            skymodel=np.full(shape, 2.), bkg_redux_skymodel=None,
                            objmodel=np.full(shape, 10.), ivarmodel=np.full(shape, .25),
                            scaleimg=np.ones(shape), waveimg=wave, bpmmask=ImageBitMaskArray(shape),
                            detector=detector, sci_spat_flexure=0., sci_spec_flexure=flex,
                            vel_type="heliocentric", vel_corr=1., slits=slits, wavesol=None,
                            tilts=np.ones(shape), maskdef_designtab=None)
        reduced.to_file(str(path))
        with fits.open(path, mode="update") as hdul:
            hdul[0].header["PYP_SPEC"] = f"p200_ngps_{channel}"
            hdul[0].header["PYPELINE"] = "MultiSlit"
            hdul[0].header["TARGET"] = "target"
            hdul[0].header["EXPTIME"] = 600.
            hdul[0].header["AIRMASS"] = 1.2
            hdul[0].header["DISPNAME"] = "VPH"
        frames[channel] = review.Frame(channel, "target", "0121", path)
    return root, frames


def test_fixed_slicers_round_trip_as_box_only(synthetic_night):
    root, frames = synthetic_night
    frame = frames["u"]
    objects = fixed.extract_fixed_slicers(frame.spec2d, -3., 2.)
    output = root / "test_spec1d.fits"
    objects.write_to_fits(objects.header, str(output))
    checked = SpecObjs.from_fitsfile(str(output))
    assert len(checked) == 3
    for obj in checked:
        assert coadd.slit_and_spat(obj.NAME) is not None
        assert obj.OPT_COUNTS is None
        np.testing.assert_allclose(obj.BOX_COUNTS, 40.)
        np.testing.assert_allclose(obj.BOX_WAVE, 5000 + np.arange(24) * 2)
        assert obj.BOX_MASK.all()
    assert checked.header["NGPSMODE"] == "FIXED"


def test_fixed_install_changes_only_checked_channel_and_preserves_2d(synthetic_night):
    root, frames = synthetic_night
    before = {c: f.spec2d.read_bytes() for c, f in frames.items()}
    for frame in frames.values():
        product = frame.spec2d.with_name(frame.spec2d.name.replace("spec2d_", "spec1d_"))
        product.write_bytes(b"original 1d")
        fluxed = frame.spec2d.parent.parent / "Fluxed" / product.name
        fluxed.parent.mkdir()
        fluxed.write_bytes(b"old fluxed")
    # Remove U's dummy input. It is optional metadata, not an extraction reference.
    frames["u"].spec2d.with_name(frames["u"].spec2d.name.replace("spec2d_", "spec1d_")).unlink()
    qa = review.audit_path(root, "target", "0121")
    aperture = fixed.FixedAperture(-3., 2., ("u",), qa_pdf=b"test QA")
    paths = {c: (f.spec2d, review.base_setup_dir(f)) for c, f in frames.items()}
    fixed.install_fixed_exposure(root, "target", "0121", paths, aperture, qa)
    assert qa.read_bytes() == b"test QA"
    assert qa.with_name("ngps_fixed_aperture_0121.json").is_file()
    for channel, frame in frames.items():
        assert frame.spec2d.read_bytes() == before[channel]
        product = frame.spec2d.with_name(frame.spec2d.name.replace("spec2d_", "spec1d_"))
        fluxed = frame.spec2d.parent.parent / "Fluxed" / product.name
        if channel == "u":
            assert product.read_bytes().startswith(b"SIMPLE")
            assert not fluxed.exists()
        else:
            assert product.read_bytes() == b"original 1d"
            assert fluxed.read_bytes() == b"old fluxed"


def test_fixed_stage_failure_keeps_all_old_products(synthetic_night):
    root, frames = synthetic_night
    qa = review.audit_path(root, "target", "0121")
    qa.write_bytes(b"previous PDF")
    paths = {c: (f.spec2d, review.base_setup_dir(f)) for c, f in frames.items()}
    with patch.object(fixed, "extract_fixed_slicers", side_effect=ValueError("bad aperture")):
        with pytest.raises(ValueError):
            fixed.install_fixed_exposure(root, "target", "0121", paths, fixed.FixedAperture(0., 4., ("u",)), qa)
    assert qa.read_bytes() == b"previous PDF"
    assert not qa.with_name("ngps_fixed_aperture_0121.json").exists()


def test_fixed_commit_failure_rolls_back_spectra_fluxed_and_qa(synthetic_night):
    root, frames = synthetic_night
    old = {}
    for channel in ("u", "g"):
        frame = frames[channel]
        product = frame.spec2d.with_name(frame.spec2d.name.replace("spec2d_", "spec1d_"))
        objects = fixed.extract_fixed_slicers(frame.spec2d, 0., 3.)
        objects.write_to_fits(objects.header, str(product))
        old[product] = product.read_bytes()
        fluxed = product.parent.parent / "Fluxed" / product.name
        fluxed.parent.mkdir()
        fluxed.write_bytes(b"old calibrated spectrum")
        old[fluxed] = fluxed.read_bytes()
    qa = review.audit_path(root, "target", "0121")
    qa.write_bytes(b"old QA")
    old[qa] = qa.read_bytes()
    destination_g = frames["g"].spec2d.with_name(frames["g"].spec2d.name.replace("spec2d_", "spec1d_"))
    replace = Path.replace
    failures = []
    def fail_once(source, destination, *args, **kwargs):
        if Path(destination) == destination_g and not failures:
            failures.append(True)
            raise OSError("simulated interrupted install")
        return replace(source, destination, *args, **kwargs)
    paths = {c: (f.spec2d, review.base_setup_dir(f)) for c, f in frames.items()}
    with patch.object(Path, "replace", autospec=True, side_effect=fail_once):
        with pytest.raises(OSError):
            fixed.install_fixed_exposure(root, "target", "0121", paths,
                                         fixed.FixedAperture(-3., 2., ("u", "g"), qa_pdf=b"new QA"), qa)
    assert failures
    for product, content in old.items():
        assert product.read_bytes() == content


@pytest.mark.parametrize("channels,available", [(("u",), "ugri"), (("u", "g", "r", "i"), "ugri"), (("r", "i"), "ri")])
def test_fixed_ui_single_click_channel_selection_no_refits(synthetic_night, channels, available):
    root, frames = synthetic_night
    frames = {c: frame for c, frame in frames.items() if c in available}
    callbacks, dialogs = {}, []
    real_button = review.Button
    class FakeButton:
        def __init__(self, axis, label, **kwargs):
            self.ax, self.label = axis, label
            self.widget = real_button(axis, label, **kwargs)
        def on_clicked(self, callback):
            callbacks[self.label] = callback
            self.widget.on_clicked(callback)
    original_dialog = fixed.FixedApertureDialog
    def popup(*args, **kwargs):
        dialog = original_dialog(*args, **kwargs)
        dialogs.append(dialog)
        return dialog
    def interact():
        callbacks["Fixed aperture"](None)
        figure = review.plt.gcf()
        axis = next(a for a in figure.axes if a.get_title().startswith(f"{available[0].upper()}:"))
        x, y = axis.transData.transform((-3., 12.))
        event = MouseEvent("button_press_event", figure.canvas, x, y, button=1)
        figure.canvas.callbacks.process("button_press_event", event)
        assert len(dialogs) == 1
        dialog = dialogs[0]
        figure.savefig(root / "fixed-popup.png")
        dialog.width.set_val("2")
        for index, channel in enumerate(dialog.channels):
            if channel not in channels:
                dialog.checks.set_active(index)
        dialog.submit()
        assert not dialog.axes
        figure.savefig(root / "fixed-preview.png")
        callbacks["Accept fixed aperture"](None)
    with patch.object(review, "Button", FakeButton), \
         patch.object(review, "FixedApertureDialog", side_effect=popup), \
         patch.object(review, "refit_central_trace", side_effect=AssertionError("No fixed trace fitting")), \
         patch.object(review, "refit_central_fwhm", side_effect=AssertionError("No fixed width fitting")), \
         patch.object(review.plt, "show", side_effect=interact):
        decision, selection = review.review_group(root, "target", "0121", frames, True)
    assert decision == "fixed"
    assert selection.channels == channels
    assert selection.offset == pytest.approx(-3.)
    assert selection.half_width == 2.
    assert selection.qa_pdf.startswith(b"%PDF")
    assert not review.audit_path(root, "target", "0121").exists()


def test_popup_invalid_width_empty_channels_and_cancel_preserve_selection():
    figure = review.plt.figure()
    accepted, closed = [], []
    dialog = fixed.FixedApertureDialog(figure, -7.2, 4., ("u", "r"), accepted.append, lambda: closed.append(True))
    dialog.width.set_val("nan")
    dialog.submit()
    assert dialog.axes and not accepted
    dialog.width.set_val("4")
    dialog.checks.set_active(0)
    dialog.checks.set_active(1)
    dialog.submit()
    assert dialog.axes and not accepted
    dialog.close()
    assert closed == [True] and not accepted
    review.plt.close(figure)


def test_coadd_uses_genuine_box_fields(synthetic_night):
    root, frames = synthetic_night
    objects = fixed.extract_fixed_slicers(frames["u"].spec2d, -3., 2.)
    # Supply a simple already-fluxed spectrum, without a sensitivity-function test.
    for obj in objects:
        obj.BOX_FLAM = obj.BOX_COUNTS * .1
        obj.BOX_FLAM_IVAR = obj.BOX_COUNTS_IVAR * 100
        obj.BOX_FLAM_SIG = obj.BOX_COUNTS_SIG * .1
    path = root / "fluxed.fits"
    objects.write_to_fits(objects.header, str(path))
    candidates = [coadd.Candidate("ngps_260623_0121.fits", str(path), obj.NAME, obj.SLITID) for obj in objects]
    assert coadd.coadd_extraction(candidates) == "BOX"
    assert len(coadd.best_trace_per_slicer(path)) == 3
    wave, flux = coadd.plot_arrays(path, objects[0].NAME)
    np.testing.assert_allclose(flux, 4.)
    input_file, _ = coadd.write_coadd_input(root / "coadd", "target", "u", "B", 43, candidates)
    assert "ex_value = BOX" in input_file.read_text()
    # PypeIt's downstream coadd reader must accept the BOX-only data.
    checked = SpecObjs.from_fitsfile(str(path))
    checked[[0]].unpack_object(ret_flam=True, extract_type="BOX")


def test_pypeit_flux_calibrates_fixed_box_spectra(synthetic_night):
    _, frames = synthetic_night
    objects = fixed.extract_fixed_slicers(frames["r"].spec2d, 0., 3.)
    wave = np.linspace(4900., 5200., 100)
    for obj in objects:
        obj.apply_flux_calib(wave, np.full_like(wave, 25.), exptime=600.)
        assert obj.OPT_FLAM is None
        assert np.all(np.isfinite(obj.BOX_FLAM[obj.BOX_MASK]))
        assert (obj.BOX_FLAM_IVAR[obj.BOX_MASK] > 0).all()


def test_mixed_coadd_uses_box_for_automatic_and_fixed_inputs(synthetic_night):
    root, frames = synthetic_night
    paths = []
    for automatic in (False, True):
        objects = fixed.extract_fixed_slicers(frames["r"].spec2d, 0., 2.)
        for obj in objects:
            obj.BOX_FLAM = obj.BOX_COUNTS * .1
            obj.BOX_FLAM_IVAR = obj.BOX_COUNTS_IVAR * 100
            if automatic:
                obj.OPT_WAVE = obj.BOX_WAVE.copy()
                obj.OPT_FLAM = obj.BOX_FLAM * 2
                obj.OPT_FLAM_IVAR = obj.BOX_FLAM_IVAR.copy()
                obj.OPT_MASK = obj.BOX_MASK.copy()
        path = root / f"mixed_{automatic}.fits"
        objects.write_to_fits(objects.header, str(path))
        paths.append(coadd.Candidate("raw.fits", str(path), objects[0].NAME, objects[0].SLITID))
    mode = coadd.coadd_extraction(paths)
    assert mode == "BOX"
    for item in paths:
        _, flux = coadd.plot_arrays(Path(item.spec1d), item.obj_id, mode)
        np.testing.assert_allclose(flux, 4.)
