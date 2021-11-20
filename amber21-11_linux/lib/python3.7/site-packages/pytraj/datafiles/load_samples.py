from __future__ import absolute_import
import os
from ..trajectory.trajectory_iterator import TrajectoryIterator

__all__ = [
    'load_sample_data',
    'load_rna',
    'load_tz2_ortho',
    'load_tz2',
    'load_ala3',
    'load_dpdp',
    'load_trpcage',
    'load_remd_ala2',
]


def load_sample_data(data_name=None):
    """
    Return TrajectoryIterator instance for Ala3 or tz2 data

    Parameters
    ----------
    data_name : str, {'ala3', 'tz2', 'rna', 'trpcage'}, default 'ala3'

    Notes
    -----
    tz2 dataset : $AMBERHOME/AmberTools/test/cpptraj/
        explicit water, ortho box
    """
    data_dict = {
        'ala3': ["ala3/Ala3.crd", "ala3/Ala3.top"],
        'tz2': ["tz2/tz2.ortho.nc", "tz2/tz2.ortho.parm7"],
        'tz2_dry': ["tz2/tz2.nc", "tz2/tz2.parm7"],
        'rna': ["rna.pdb", "rna.pdb"],
        'trpcage': ["trpcage/trpcage.pdb.gz", "trpcage/trpcage.pdb.gz"],
        'dpdp': ["dpdp/DPDP.nc", "dpdp/DPDP.parm7"],
        'remd_ala2': [[
            "remd_ala2/rem.nc.000", "remd_ala2/rem.nc.001",
            "remd_ala2/rem.nc.002", "remd_ala2/rem.nc.003"
        ], "remd_ala2/ala2.parm7"],
    }

    mydir = os.path.dirname(os.path.abspath(__file__))
    if data_name is None:
        data_name = 'ala3'
    if data_name == 'remd_ala2':
        crd = [os.path.join(mydir, fn) for fn in data_dict[data_name][0]]
    else:
        crd = os.path.join(mydir, data_dict[data_name][0])
    top = os.path.join(mydir, data_dict[data_name][1])
    return TrajectoryIterator(crd, top)


def load_rna():
    '''return pytraj.TrajectoryIterator for an RNA trajectory with 3 frames
    '''
    return load_sample_data('rna')


def load_dpdp():
    '''return pytraj.TrajectoryIterator for an RNA trajectory with 3 frames
    '''
    return load_sample_data('dpdp')


def load_trpcage():
    '''return pytraj.TrajectoryIterator for an trp-cage, 38 frames, NMR (pdb: 1l2y)

    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_trpcage()
    >>> traj.n_frames
    38
    >>> traj.n_atoms
    304

    '''
    return load_sample_data('trpcage')


def load_remd_ala2():
    return load_sample_data('remd_ala2')


def load_tz2_ortho():
    return load_sample_data('tz2')


def load_tz2():
    return load_sample_data('tz2_dry')


def load_ala3():
    return load_sample_data('ala3')
