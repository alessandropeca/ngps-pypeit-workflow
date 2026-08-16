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


def test_interactive_coadd_tool_is_present():
    assert (ROOT / "scripts" / "ngps_interactive_coadd.py").is_file()
