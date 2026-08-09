from pathlib import Path
import py_compile
import subprocess
import sys

ROOT = Path(__file__).parents[1]


def test_lock_is_valid():
    result = subprocess.run(
        [sys.executable, "tools/check_lock.py"], cwd=ROOT, capture_output=True, text=True
    )
    assert result.returncode == 0, result.stderr


def test_wrappers_compile():
    for script in (ROOT / "scripts").glob("*.py"):
        py_compile.compile(str(script), doraise=True)

