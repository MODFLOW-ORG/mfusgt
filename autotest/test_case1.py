"""
Tests the binary concentration output of GSI Documented Example 1, Case 1 (1D).

Each available executable in bin/ runs the case once with a single concentration
file and once with the BCT MULTIFILE option.

Cases:
  - test_case1 : both runs terminate normally, FloPy reads every record of both
                 concentration files, and the two files are byte-identical.
"""

import shutil
import subprocess
import sys
from pathlib import Path

import flopy
import pytest

ROOT = Path(__file__).resolve().parents[1]
DATA = Path(__file__).resolve().parent / "data" / "case1"
EXT = ".exe" if sys.platform == "win32" else ""
TARGETS = ["mfusgt", "mfusgtdbl"]
NTIMES = 1000

# BCT line 1 is read as 19 fixed fields of 10 characters; options follow them
OPTION_COLUMN = 190


def run_case(exe, ws, multifile):
    shutil.copytree(DATA, ws)
    if multifile:
        fpth = ws / "case1.bct"
        lines = fpth.read_text().splitlines()
        lines[0] = lines[0].rstrip().ljust(OPTION_COLUMN) + " MULTIFILE mc"
        fpth.write_text("\n".join(lines) + "\n")
    proc = subprocess.run(
        [str(exe), "case1.nam"], cwd=ws, capture_output=True, text=True
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "Normal termination" in proc.stdout


def read_conc(fpth):
    cobj = flopy.utils.HeadFile(fpth, text="CONC")
    times = cobj.get_times()
    assert len(times) == NTIMES
    return cobj.get_data(totim=times[-1])


@pytest.mark.parametrize("target", TARGETS)
def test_case1(tmp_path, target):
    exe = ROOT / "bin" / f"{target}{EXT}"
    if not exe.is_file():
        pytest.skip(f"{exe.name} is not built")

    run_case(exe, tmp_path / "single", multifile=False)
    run_case(exe, tmp_path / "multi", multifile=True)

    single = tmp_path / "single" / "case1.con"
    multi = tmp_path / "multi" / "MC CONC.CON"
    conc = read_conc(single)
    assert conc.shape == (2, 2, 101)
    read_conc(multi)
    assert single.read_bytes() == multi.read_bytes()
