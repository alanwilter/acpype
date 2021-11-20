#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Requirements:
* pytest
* pytest-html

To run:
pytest -s --html=report.html

"""
import os
import pytest
import ujson as json
from acpype_lib.acs_api import acpype_api

os.chdir(os.path.dirname(os.path.abspath(__file__)))
json_output = acpype_api(
    inputFile="AAA.mol2", chargeType="gas", atomType="gaff2", debug=True, basename="AAA", chiral=False
)
jj = json.loads(json_output)


@pytest.mark.parametrize(
    ("atype"),
    [
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
    ],
)
def test_json(atype):
    assert len(jj.get(atype)) > 2
