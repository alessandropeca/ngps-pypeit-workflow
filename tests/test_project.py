import csv
import os
from pathlib import Path
import py_compile
import subprocess
import sys
import tempfile

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


def test_interactive_coadd_tool_is_present():
    assert (ROOT / "scripts" / "ngps_interactive_coadd.py").is_file()


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

        plan[0]["standard_filename"] = "ngps_0000.fits"
        with associations.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=plan[0].keys())
            writer.writeheader()
            writer.writerows(plan)
        result = subprocess.run(command, cwd=ROOT, env=environment, capture_output=True, text=True)
        assert result.returncode == 0, result.stdout + result.stderr
        assert "target_a [0001,0002] -> std_a" in result.stdout
        assert "manual" in result.stdout
