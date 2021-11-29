import os
import ujson as json
from acpype_lib.acs_api import acpype_api

file_types = (
    ("file_name"),
    ("em_mdp"),
    ("AC_frcmod"),
    ("AC_inpcrd"),
    ("AC_lib"),
    ("AC_prmtop"),
    ("mol2"),
    ("CHARMM_inp"),
    ("CHARMM_prm"),
    ("CHARMM_rtf"),
    ("CNS_inp"),
    ("CNS_par"),
    ("CNS_top"),
    ("GMX_OPLS_itp"),
    ("GMX_OPLS_top"),
    ("GMX_gro"),
    ("GMX_itp"),
    ("GMX_top"),
    ("NEW_pdb"),
    ("md_mdp"),
)


def test_json_simple():
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    jj = json.loads(acpype_api(inputFile="benzene.pdb", debug=True))
    assert min([len(jj.get(x)) for x in file_types]) >= 7
    assert jj.get("file_name") == "benzene"


def test_json_failed():
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    jj = json.loads(acpype_api(inputFile="_fake_", debug=True))
    assert "ERROR: [Errno 2] No such file or directory" in jj.get("file_name")
    assert "tests/_fake_" in jj.get("file_name")


def get_json():
    json_output = acpype_api(inputFile="AAA.mol2", chargeType="gas", atomType="gaff2", debug=True, basename="AAA")
    return json.loads(json_output)


def test_json():
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    jj = get_json()
    for ft in file_types:
        assert len(jj.get(ft)) > 2
