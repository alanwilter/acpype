from __future__ import absolute_import, print_function
from pytraj.core import *
from pytraj.math import *
from pytraj.externals import *
from pytraj.trajectory.c_traj import *
from pytraj.analysis.hbond_analysis import *
from pytraj.datafiles.load_samples import load_sample_data
from pytraj.trajectory.trajectory import Trajectory

from pytraj import *
from pytraj.datasets import *
from pytraj.all_actions import *

from pytraj.core import c_dict
from pytraj.utils.misc import get_atts


def run_tests():
    print("try to load sample data")
    load_sample_data()
    load_sample_data('tz2')

    print("try to make all action objects")
    failed_list = [
        'createreservoir',
    ]
    DatasetList()
    print("try to make all analysis objects")
    from pytraj import analdict
    failed_list = []

    for key in analdict.keys():
        if key not in failed_list:
            analdict[key]

    print("try to make all dataset stuff")
    DatasetDouble()
    DatasetFloat()
    DatasetInteger()
    DatasetString()
    DatasetMatrixDouble()
    DatasetGridFloat()
    DatasetMatrixFloat()
    DatasetVector()
    DatasetMatrix3x3()
    DatasetCoords()
    DatasetCoordsRef()
    DatasetCoordsCRD()

    print("try to make structure-related objects")
    Topology()
    Molecule()
    Residue()
    Atom()
    Box()
    Frame()

    print("try to create Trajectory-like objects")
    TrajectoryIterator()
    Trajectory()

    print("other stuff. throw all tests don't belong anywhere else here")
    keys = get_atts(c_dict)
    cdict = c_dict.__dict__

    for key in keys:
        if isinstance(cdict[key], dict):
            assert cdict[key].keys() is not None

    # other objects
    CpptrajState()

    print("OK")


if __name__ == '__main__':
    run_tests()
