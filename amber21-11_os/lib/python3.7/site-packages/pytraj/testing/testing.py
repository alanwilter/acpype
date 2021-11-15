from __future__ import absolute_import, print_function
import os

from pytraj.datafiles.load_samples import load_sample_data
from pytraj.utils import eq
from pytraj.utils import Timer
from pytraj.utils import tempfolder
from pytraj.utils import duplicate_traj
from numpy.testing import assert_almost_equal
from functools import partial

aa_eq = partial(assert_almost_equal, decimal=4)

__all__ = [
    'load_sample_data', 'eq', 'aa_eq', 'duplicate_traj', 'Timer', 'tempfolder',
    'amberhome', 'cpptraj_test_dir', 'get_fn', 'get_remd_fn',
    'assert_equal_topology'
]

# find cpptraj test dir
# find in CPPTRAJHOME first
# if not exist CPPTRAJHOME, find in AMBERHOME
# if not exist CPPTRAJHOME and AMBERHOME, find in ../cpptraj/ folder
# (supposed we run pytraj's test in pytraj/tests folder)

# if you are adding more test files to cpptraj master branch on github,
# you should export CPPTRAJHOME

cpptrajhome = os.environ.get('CPPTRAJHOME', '')
amberhome = os.environ.get('AMBERHOME', '')

if cpptrajhome:
    cpptraj_test_dir = os.path.join(cpptrajhome, 'test')
else:
    if amberhome:
        cpptraj_test_dir = os.path.join(amberhome, 'AmberTools', 'src',
                                        'cpptraj', 'test')
    else:
        cpptrajhome = ''
        amberhome = ''
        cpptraj_test_dir = ''

DEFAULT_PATH = os.path.dirname(__file__) + "/../../cpptraj/test"

if os.path.exists(DEFAULT_PATH):
    cpptraj_test_dir = os.path.abspath(DEFAULT_PATH)

header_doc = '''
import pytraj as pt
import pytraj.all_actions as pyca
traj = pt.load_sample_data("tz2")
'''


def assert_equal_topology(top, new_top, traj=None):
    import pytraj as pt
    assert new_top.n_atoms == top.n_atoms, 'same n_atoms'
    assert new_top.n_residues == top.n_residues, 'same n_residues'
    assert new_top.n_mols == top.n_mols, 'same n_mols'
    # there are inverted bond indices [5292 5291] vs [5291 5292]
    # so use distance to assert
    if traj is not None:
        aa_eq(
            pt.distance(traj, new_top.bond_indices),
            pt.distance(traj, top.bond_indices))
        # same for dihedral_indices
        aa_eq(
            pt.dihedral(traj, new_top.dihedral_indices),
            pt.dihedral(traj, top.dihedral_indices))
    aa_eq(new_top.dihedral_indices, top.dihedral_indices)
    aa_eq(new_top.mass, top.mass)
    aa_eq(new_top.charge, top.charge)
    aa_eq(new_top.box.values, top.box.values)

    assert [res.name for res in top.residues] == [
        res.name for res in new_top.residues
    ], 'equal resnames'
    assert [atom.name for atom in top.atoms] == [
        atom.name for atom in new_top.atoms
    ], 'equal resnames'

    for res, res_new in zip(top.residues, new_top.residues):
        assert res.first_atom_index == res_new.first_atom_index, 'first atom'
        assert res.last_atom_index == res_new.last_atom_index, 'last atom'


def get_fn(txt):
    '''get absolute path for trajectory and topology samples. Legit text = 'ala3', 'tz2',
    'rna'. Mostly for testing purpose.

    Examples
    --------
    >>> # get trajectory file name
    >>> fname = get_fn('tz2')[0]
    >>> fname.split('/')[-1]
    'tz2.ortho.nc'
    '''
    from pytraj import load_sample_data
    traj = load_sample_data(txt)
    return traj.filename, traj.top.filename


def get_remd_fn(txt):
    """

    >>> fnlist, tn = get_remd_fn('remd_ala2')
    """
    from pytraj import load_sample_data
    traj = load_sample_data(txt)
    return traj.filelist, traj.top.filename


if __name__ == "__main__":
    print(amberhome)
    print(cpptraj_test_dir)
