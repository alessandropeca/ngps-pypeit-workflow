import csv
import importlib.util
import os
from pathlib import Path
import py_compile
import subprocess
import sys
import tempfile
from unittest.mock import patch

ROOT = Path(__file__).parents[1]


def test_lock_is_valid():
    result = subprocess.run(
        [sys.executable, "tools/check_lock.py"], cwd=ROOT, capture_output=True, text=True
    )
    assert result.returncode == 0, result.stderr


def test_wrappers_compile():
    # Keep test artifacts out of the tracked scripts directory.
    with tempfile.TemporaryDirectory() as temporary:
        output_directory = Path(temporary)
        for script in (ROOT / "scripts").glob("*.py"):
            py_compile.compile(
                str(script),
                cfile=str(output_directory / f"{script.stem}.pyc"),
                doraise=True,
            )


def test_interactive_extraction_tool_is_present():
    assert (ROOT / "scripts" / "ngps_interactive_extract.py").is_file()
    assert (ROOT / "scripts" / "ngps_review_target_extractions.py").is_file()
    assert (ROOT / "scripts" / "ngps_manual_target_extractions.py").is_file()


def test_manual_trace_refit_patch_and_target_config():
    source = ROOT / "scripts" / "ngps_interactive_extract.py"
    spec = importlib.util.spec_from_file_location("ngps_interactive_extract", source)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)

    patch = ROOT / "patches" / "pypeit-manual-trace-refit.patch"
    assert patch.is_file()
    assert "manual_refit_trace" in patch.read_text()

    with tempfile.TemporaryDirectory() as temporary:
        destination = Path(temporary) / "target.pypeit"
        destination.write_text("[rdx]\nsetup read\nsetup end\n")
        module.enable_manual_trace_refit(destination)
        text = destination.read_text()

    assert "manual_refit_trace = True" in text
    assert "trace_maxshift = 3.0" in text


def test_manual_selection_can_be_linked_or_channel_only():
    source = ROOT / "scripts" / "ngps_manual_target_extractions.py"
    spec = importlib.util.spec_from_file_location("ngps_manual_target_extractions", source)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.path.insert(0, str(ROOT / "scripts"))
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    frames = {channel: object() for channel in ("u", "g", "r", "i")}
    selected = {"r": -1.0}

    module.apply_manual_selection(selected, frames, "g", 2.5, True)
    assert selected == {"r": -1.0, "g": 2.5}

    module.apply_manual_selection(selected, frames, "i", -3.0, False)
    assert selected == {"u": -3.0, "g": -3.0, "r": -3.0, "i": -3.0}


def test_discover_frames_falls_back_to_a_target_run(tmp_path):
    source = ROOT / "scripts" / "ngps_manual_target_extractions.py"
    spec = importlib.util.spec_from_file_location("ngps_manual_target_extractions_discovery", source)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.path.insert(0, str(ROOT / "scripts"))
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    science = (tmp_path / "manual_setup_r" / "p200_ngps_r_C" / ".ngps_target_runs"
               / "p200_ngps_r_C_manual_0106" / "Science")
    science.mkdir(parents=True)
    (science / "spec2d_ngps_260727_0106-target_NGPS_r_20260727T045817.156.fits").touch()
    frames = module.discover_frames(tmp_path)
    assert len(frames) == 1
    assert frames[0].channel == "r"


def test_target_run_frame_resolves_to_its_baseline_setup():
    source = ROOT / "scripts" / "ngps_manual_target_extractions.py"
    spec = importlib.util.spec_from_file_location("ngps_manual_target_extractions", source)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.path.insert(0, str(ROOT / "scripts"))
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    setup = Path("/tmp/night/manual_setup_r/p200_ngps_r_C")
    frame = module.Frame(
        "r", "target", "0106",
        setup / ".ngps_target_runs" / "p200_ngps_r_C_manual_0106" / "Science" / "spec2d_test.fits",
    )
    assert module.base_setup_dir(frame) == setup


def test_interactive_coadd_tool_is_present():
    assert (ROOT / "scripts" / "ngps_interactive_coadd.py").is_file()


def test_telluric_plan_uses_only_completed_fluxed_ri_coadds():
    source = ROOT / "scripts" / "ngps_telluric_correct.py"
    spec = importlib.util.spec_from_file_location("ngps_telluric_correct", source)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)

    class FakeSpectrum:
        header = {"FLUXED": True}

    class FakeHdul:
        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return False

        def __getitem__(self, _name):
            return FakeSpectrum()

    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        rows = []
        for channel in ("u", "r", "i"):
            setup = f"p200_ngps_{channel}_B"
            rows.append({
                "target": "target_a", "channel": channel.upper(), "setup": setup,
                "status": "coadded",
            })
            path = module.coadd_path(root, "target_a", channel, setup)
            path.parent.mkdir(parents=True)
            path.touch()
        with (root / "coadd_review.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=("target", "channel", "setup", "status"))
            writer.writeheader()
            writer.writerows(rows)

        with patch.object(module.fits, "open", return_value=FakeHdul()):
            products, skipped = module.selected_products(root, None, ("r", "i"), None)

    assert skipped == []
    assert [product.channel for product in products] == ["i", "r"]


def test_flux_plan_groups_consecutive_exposures_and_allows_a_manual_standard():
    with tempfile.TemporaryDirectory() as temporary:
        work_root = Path(temporary)
        night = work_root / "20260101"
        setup = night / "manual_setup_r" / "p200_ngps_r_B"
        science = setup / "Science"
        science.mkdir(parents=True)
        fields = [
            "channel", "setup", "binning", "frametype", "filename", "target",
            "mjd", "airmass", "exptime",
        ]
        rows = [
            ["R", "p200_ngps_r_B", "3,2", "standard", "ngps_0000.fits", "std_a", "1.0", "1.0", "10"],
            ["R", "p200_ngps_r_B", "3,2", "science", "ngps_0001.fits", "target_a", "1.10", "1.0", "600"],
            ["R", "p200_ngps_r_B", "3,2", "science", "ngps_0002.fits", "target_a", "1.20", "1.0", "600"],
            ["R", "p200_ngps_r_B", "3,2", "standard", "ngps_0003.fits", "std_b", "1.25", "1.0", "10"],
        ]
        with (night / "science_standard_inventory.csv").open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(fields)
            writer.writerows(rows)
        for filename in ("ngps_0000", "ngps_0001", "ngps_0002", "ngps_0003"):
            (science / f"spec1d_{filename}-test.fits").touch()

        environment = dict(os.environ, NGPS_WORK_ROOT=str(work_root))
        command = [sys.executable, "scripts/ngps_flux_calibrate.py", "20260101"]
        result = subprocess.run(command, cwd=ROOT, env=environment, capture_output=True, text=True)
        assert result.returncode == 0, result.stdout + result.stderr

        associations = night / "science_standard_associations.csv"
        with associations.open(newline="") as handle:
            plan = list(csv.DictReader(handle))
        assert len(plan) == 1
        assert plan[0]["science_filenames"] == "ngps_0001.fits ngps_0002.fits"
        assert plan[0]["standard_filename"] == "ngps_0003.fits"
        assert plan[0]["assignment_status"] == "automatic"

        plan[0]["standard_filename"] = "ngps_0000.fits"
        with associations.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=plan[0].keys())
            writer.writeheader()
            writer.writerows(plan)
        result = subprocess.run(command, cwd=ROOT, env=environment, capture_output=True, text=True)
        assert result.returncode == 0, result.stdout + result.stderr
        assert "target_a [0001,0002] -> std_a" in result.stdout
        assert "manual" in result.stdout


def test_flux_plan_proposes_a_reviewable_fallback_after_a_standard_failure():
    source = ROOT / "scripts" / "ngps_flux_calibrate.py"
    spec = importlib.util.spec_from_file_location("ngps_flux_calibrate", source)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    failed = {"filename": "hz44.fits", "target": "hz44", "mjd": 1.0, "sensfile": Path("failed.fits")}
    fallback = {"filename": "wolf.fits", "target": "wolf", "mjd": 1.2, "sensfile": Path("fallback.fits")}
    group = {
        "id": "target__0001-0003", "midpoint": 1.1, "standard": failed,
        "status": "automatic", "manual": False,
    }
    plan = {"standards": [failed, fallback], "groups": [group]}

    def fake_ensure(_plan, standard, _force, attempted):
        success = standard["filename"] == "wolf.fits"
        attempted[Path(standard["sensfile"])] = success
        return success

    with patch.object(module, "ensure_sensfunc", side_effect=fake_ensure):
        changed, unresolved = module.propose_fallbacks([plan], False, {Path("failed.fits"): False})

    assert changed == 1
    assert unresolved == []
    assert group["standard"] is fallback
    assert group["status"] == "automatic fallback"


def test_coadd_review_discards_singletons_and_keeps_fluxed_repeats():
    source = ROOT / "scripts" / "ngps_interactive_coadd.py"
    spec = importlib.util.spec_from_file_location("ngps_interactive_coadd", source)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)

    singleton = module.ObservationGroup(
        "target_a", "r", "p200_ngps_r_A",
        [{"filename": "ngps_0001.fits", "mjd": "1.0"}],
    )
    repeated = module.ObservationGroup(
        "target_b", "r", "p200_ngps_r_B",
        [
            {"filename": "ngps_0002.fits", "mjd": "2.0"},
            {"filename": "ngps_0003.fits", "mjd": "3.0"},
        ],
    )
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        fluxed = root / "manual_setup_r" / "p200_ngps_r_B" / "Fluxed"
        fluxed.mkdir(parents=True)
        (fluxed / "spec1d_ngps_0002-test.fits").touch()
        (fluxed / "spec1d_ngps_0003-test.fits").touch()
        review = module.update_coadd_review(root, [singleton, repeated])

    assert review[module.review_key(singleton)]["status"] == "discard"
    assert review[module.review_key(singleton)]["reason"] == "only one science exposure"
    assert module.reviewable_groups([singleton, repeated], review) == [repeated]
    review[module.review_key(repeated)]["status"] = "coadded"
    assert module.reviewable_groups([singleton, repeated], review) == []
    assert module.reviewable_groups([singleton, repeated], review, include_coadded=True) == [repeated]
