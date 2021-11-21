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
import ujson as json
from acpype_lib.acs_api import acpype_api


def get_json():
    json_output = acpype_api(
        inputFile="AAA.mol2", chargeType="gas", atomType="gaff2", debug=True, basename="AAA", chiral=False
    )
    return json.loads(json_output)


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


def test_json():
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    jj = get_json()
    for ft in file_types:
        assert len(jj.get(ft)) > 2
