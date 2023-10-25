import sys
from unittest.mock import patch

import pytest

from acpype.cli import _chk_py_ver, init_main
from acpype.utils import _getoutput


def test_no_ac(janitor, capsys):
    binaries = {"ac_bin": "no_ac", "obabel_bin": "obabel"}
    msg = f"ERROR: no '{binaries['ac_bin']}' executable... aborting!"
    inp = "AAA.mol2"
    with pytest.raises(SystemExit) as e_info:
        init_main(argv=["-di", inp, "-c", "gas"], binaries=binaries)
    captured = capsys.readouterr()
    assert msg in captured.out
    assert e_info.typename == "SystemExit"
    assert e_info.value.code == 19


def test_only_ac(janitor, capsys):
    binaries = {"ac_bin": "antechamber", "obabel_bin": "no_obabel"}
    msg1 = f"WARNING: no '{binaries['obabel_bin']}' executable, no PDB file can be used as input!"
    msg2 = "Total time of execution:"
    msg3 = "WARNING: No Openbabel python module, no chiral groups"
    inp = "AAA.mol2"
    temp_base = "vir_temp"
    init_main(argv=["-di", inp, "-c", "gas", "-b", temp_base], binaries=binaries)
    captured = capsys.readouterr()
    assert msg1 in captured.out
    assert msg2 in captured.out
    assert msg3 in captured.out
    _getoutput(f"rm -vfr {temp_base}* .*{temp_base}*")


@pytest.mark.parametrize(
    ("inp", "msg"),
    [
        ("AAA.pdb", "ERROR: no 'no_obabel' executable; you need it if input is PDB or SMILES"),
        ("cccc", "WARNING: your input may be a SMILES but"),
    ],
)
def test_no_obabel(janitor, capsys, inp, msg):
    binaries = {"ac_bin": "antechamber", "obabel_bin": "no_obabel"}
    with pytest.raises(SystemExit) as e_info:
        init_main(argv=["-di", inp, "-c", "gas"], binaries=binaries)
    captured = capsys.readouterr()
    assert msg in captured.out
    assert e_info.typename == "SystemExit"
    assert e_info.value.code == 19


def test_amb2gmx_no_bins(janitor, capsys):
    binaries = {"ac_bin": "no_ac", "obabel_bin": "no_obabel"}
    argv = ["-x", "Base.inpcrd", "-p", "Base.prmtop"]
    temp_base = "vir_temp"
    init_main(argv=argv + ["-b", temp_base], binaries=binaries)
    captured = capsys.readouterr()
    assert "Total time of execution:" in captured.out
    _getoutput(f"rm -vfr {temp_base}* .*{temp_base}*")


def test_chk_py_ver_python_3_9_or_higher():
    # This should not raise an exception for Python 3.9 or higher
    _chk_py_ver()


def test_chk_py_ver_python_3_8():
    # Mock sys.version_info to mimic Python 3.8
    with patch.object(sys, "version_info", (3, 8)):
        with pytest.raises(Exception, match="Sorry, you need python 3.9 or higher"):
            # This should raise an exception for Python 3.8
            # Ensure that the exception message matches the expected one
            _chk_py_ver()
