from .shared_trajectory import SharedTrajectory
from .shared_methods import _xyz, _box
from .frame import Frame


class StrippedTrajectoryIterator(SharedTrajectory):
    """Out-of-core, indexable StrippedTrajectory
    This class is used for NGLView, nothing else.

    versionadded: 1.0.7

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> straj = traj.autoimage().superpose('@CA').strip(":WAT")
    >>> view = straj.view()
    """

    def __init__(self, origtraj, mask):
        top = origtraj.top.copy()
        top.strip(mask)
        self.top = top
        self._traj = origtraj
        self._atm = self._traj.top(mask)
        self._atm.invert_mask()
        self._mask = mask

    def __iter__(self):
        for frame in self._traj:
            yield Frame(frame, self._atm)

    def __getitem__(self, index):
        traj = self._traj[index]
        if isinstance(traj, Frame):
            return Frame(traj, self._atm)
        else:
            return traj.strip(self._mask)

    @property
    def xyz(self):
        return _xyz(self)

    @property
    def unitcells(self):
        return _box(self)

    @property
    def n_frames(self):
        return self._traj.n_frames

    @property
    def n_atoms(self):
        return self.top.n_atoms
