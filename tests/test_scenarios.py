import os
import pytest
from acpype import cli
from acpype.utils import _getoutput


def test_no_ac(capsys):
    binaries = {"ac_bin": "no_ac", "obabel_bin": "obabel"}
    msg = f"ERROR: no '{binaries['ac_bin']}' executable... aborting!"
    inp = "AAA.mol2"
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    with pytest.raises(SystemExit) as e_info:
        cli.init_main(argv=["-di", inp, "-c", "gas"], binaries=binaries)
    captured = capsys.readouterr()
    assert msg in captured.out
    assert e_info.typename == "SystemExit"
    assert e_info.value.code == 19


def test_only_ac(capsys):
    binaries = {"ac_bin": "antechamber", "obabel_bin": "no_obabel"}
    msg1 = f"WARNING: no '{binaries['obabel_bin']}' executable, no PDB file can be used as input!"
    msg2 = "Total time of execution:"
    msg3 = "WARNING: No Openbabel python module, no chiral groups"
    inp = "AAA.mol2"
    temp_base = "vir_temp"
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    cli.init_main(argv=["-di", inp, "-c", "gas", "-b", temp_base], binaries=binaries)
    captured = capsys.readouterr()
    assert msg1 in captured.out
    assert msg2 in captured.out
    assert msg3 in captured.out
    _getoutput(f"rm -fr {temp_base}*")


@pytest.mark.parametrize(
    ("inp", "msg"),
    [
        ("AAA.pdb", "ERROR: no 'no_obabel' executable; you need it if input is PDB or SMILES"),
        ("cccc", "WARNING: your input may be a SMILES but"),
    ],
)
def test_no_obabel(capsys, inp, msg):
    binaries = {"ac_bin": "antechamber", "obabel_bin": "no_obabel"}
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    with pytest.raises(SystemExit) as e_info:
        cli.init_main(argv=["-di", inp, "-c", "gas"], binaries=binaries)
    captured = capsys.readouterr()
    assert msg in captured.out
    assert e_info.typename == "SystemExit"
    assert e_info.value.code == 19


def test_amb2gmx_no_bins(capsys):
    binaries = {"ac_bin": "no_ac", "obabel_bin": "no_obabel"}
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    argv = ["-x", "Base.inpcrd", "-p", "Base.prmtop"]
    temp_base = "vir_temp"
    cli.init_main(argv=argv + ["-b", temp_base], binaries=binaries)
    captured = capsys.readouterr()
    assert "Total time of execution:" in captured.out
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    _getoutput(f"rm -fr {temp_base}*")
