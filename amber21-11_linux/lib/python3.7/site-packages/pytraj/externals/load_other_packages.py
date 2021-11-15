from __future__ import absolute_import
from ..utils.get_common_objects import _load_Topology
from ..utils.context import tempfolder


def load_parmed(parm, traj=True, **kwd):
    """return pytraj's Topology or Trajectory objects

    Parameters
    ----------
    parm : ParmEd's Structure object
    traj: bool, default True
        if True, return pytraj.Trajectory
        if False, return Topology

    >>> import parmed as pmd
    >>> import pytraj as pt
    >>> p = pmd.download_PDB("1l2y")
    >>> traj = pt.load_parmed(p)
    """
    with tempfolder():
        fname = 'tmp.parm7'
        parm.save(fname)
        top = _load_Topology(fname)
    if traj:
        from pytraj import Trajectory
        coords = parm.get_coordinates()
        traj = Trajectory(xyz=coords, top=top)
        traj.unitcells = parm.get_box()
        return traj
    else:
        return top

