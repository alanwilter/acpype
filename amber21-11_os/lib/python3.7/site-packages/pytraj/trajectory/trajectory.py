from __future__ import absolute_import

import numpy as np
from ..core.box import Box
from .frame import Frame
from ..utils.check_and_assert import is_int, is_frame_iter
from ..utils.convert import array_to_cpptraj_atommask
from ..core.c_core import AtomMask
from .shared_trajectory import SharedTrajectory

# use absolute import here
from pytraj.utils.get_common_objects import get_topology

from ..topology.topology import Topology
from .shared_methods import iterframe_master
from ..utils.cyutils import _fast_iterptr, _fast_iterptr_withbox
from .frameiter import FrameIterator

__all__ = ['Trajectory']


class Trajectory(SharedTrajectory):
    """Simple in-memory Trajectory. It has only information about 3D coordinates
    and unitcells (no time, no velocity, no force, ...)

    Parameters
    ----------
    filename: str, trajectory filename
    top : Topology
    xyz: 3D-array
        if filename is not given, pytraj will construct Trajectory from given
        Topology and given xyz array.

    Examples
    --------
    >>> import pytraj as pt
    >>> from pytraj.testing import get_fn
    >>> fn, tn = get_fn('ala3')

    >>> # load from filename and topology name
    >>> traj = pt.Trajectory(fn, tn)
    >>> # load from a list of filenames
    >>> traj = pt.Trajectory([fn, fn], tn)

    >>> # load from 3D array-like
    >>> xyz = traj.xyz
    >>> traj_1 = pt.Trajectory(xyz=xyz, top=traj.top)

    >>> # get new Trajectory with only CA atoms
    >>> traj['@CA'].xyz[:, :, 0]
    array([[  3.970048 ,   7.6400076,  10.1610562]])

    >>> # iterate
    >>> for frame in traj: pass
    """

    def __init__(self,
                 filename=None,
                 top=None,
                 xyz=None,
                 velocity=None,
                 force=None,
                 time=None):
        self._top = get_topology(filename, top)
        if velocity is not None:
            velocity = np.asarray(velocity, dtype='f8')
        if force is not None:
            force = np.asarray(force, dtype='f8')
        if time is not None:
            time = np.asarray(time, dtype='f8')
        self.velocities = velocity
        self.forces = force
        self.time = time

        if self._top is None:
            self._top = Topology()

        self._xyz = None
        self._boxes = None

        # use those to keep lifetime of Frame
        self._life_holder = None
        self._frame_holder = None

        if filename is None or filename == "":
            if xyz is not None:
                if self.top.is_empty():
                    raise ValueError("must have a non-empty Topology")
                else:
                    assert self.top.n_atoms == xyz.shape[
                        1], "must have the same n_atoms"
                self._xyz = np.asarray(xyz)
            else:
                self._xyz = None
        elif hasattr(filename, 'xyz'):
            # make sure to use `float64`
            self._xyz = filename.xyz.astype(np.float64)
        elif isinstance(filename, (str, list, tuple)):
            self.load(filename)
        else:
            raise ValueError('filename must be None, a Trajectory or a string')

        if hasattr(self._xyz, 'shape'):
            assert self.top.n_atoms == self._xyz.shape[
                1], "must have the same n_atoms"

        if hasattr(filename, 'unitcells'):
            self._boxes = filename.unitcells

    @property
    def top(self):
        '''Topology (having info about atom, residue, ...)

        See also
        --------
        pytraj.Trajectory.topology
        '''
        return self._top

    @top.setter
    def top(self, value):
        '''Topology
        '''
        # always make a copy
        self._top = value.copy()

    @property
    def topology(self):
        '''longer name for ``top``

        >>> import pytraj as pt
        >>> traj = pt.datafiles.load_ala3()[:]
        >>> traj.topology.n_residues
        3
        >>> ca_indices = pt.select('@CA', traj.topology)
        '''
        return self.top

    @top.setter
    def topology(self, value):
        '''Topology
        '''
        self.top = value

    @property
    def xyz(self):
        '''Trajectory's coordinates, shape=(n_frames, n_frames, 3), dtype='f8'

        Examples
        --------
        >>> import pytraj as pt
        >>> traj0 = pt.datafiles.load_ala3()
        >>> traj1 = pt.Trajectory(xyz=np.empty((traj0.n_frames, traj0.n_atoms, 3), dtype='f8'), top=traj0.top)
        >>> traj1.xyz = traj0.xyz.copy()
        >>> # autoconvert from fortran order to c order
        >>> xyz = np.asfortranarray(traj0.xyz)
        >>> traj1.xyz = xyz
        '''
        return self._xyz

    @xyz.setter
    def xyz(self, values):
        '''assign new coordinates for Trajectory
        '''
        if self.shape[1] and self.n_atoms != values.shape[1]:
            raise ValueError("must have the same number of atoms")
        # make sure dtype='f8'
        values = np.asarray(values, dtype='f8')
        if not values.flags['C_CONTIGUOUS']:
            # autoconvert
            values = np.ascontiguousarray(values)
        self._xyz = values

    def copy(self):
        '''return a deep copy of trajectory

        Examples
        --------
        >>> import pytraj as pt
        >>> t0 = pt.datafiles.load_rna()[:]
        >>> isinstance(t0.copy(), pt.Trajectory)
        True
        '''
        traj = self.__class__()
        traj.top = self.top.copy()
        traj.xyz = self._xyz.copy()
        return traj

    @property
    def shape(self):
        '''(n_frames, n_atoms, 3)
        '''
        try:
            return self._xyz.shape
        except AttributeError:
            return (None, None, 3)

    @property
    def n_atoms(self):
        '''n_atoms
        '''
        return self.top.n_atoms

    @property
    def n_frames(self):
        '''n_frames
        '''
        try:
            n_frames = self._xyz.shape[0]
        except (AttributeError, IndexError):
            n_frames = 0
        return n_frames

    def __iter__(self):
        """return a Frame view of coordinates

        Notes
        -----
        update frame view will update Trajectory.xyz too
        if want to use listcomp, need to make copy for every frame `[frame.copy() for frame in traj]`

        Examples
        --------
        >>> import pytraj as pt
        >>> traj = pt.load_sample_data('tz2')[:]
        >>> for frame in traj: pass
        """
        indices = range(self.n_frames)
        return self._iterframe_indices(indices)

    def _iterframe_indices(self, indices):
        """return a Frame view of coordinates

        Notes
        -----
        update frame view will update Trajectory.xyz too

        Examples
        --------
        >>> import pytraj as pt
        >>> traj = pt.load_sample_data('tz2')[:]
        >>> for frame in traj._iterframe_indices([3, 5 ,7]): pass
        """

        if self._boxes is None:
            iter_obj = _fast_iterptr(self.xyz, self.n_atoms, indices, self.top)
        else:
            iter_obj = _fast_iterptr_withbox(self.xyz, self._boxes,
                                             self.n_atoms, indices, self.top)

        for (index, frame) in iter_obj:
            if self.forces is not None:
                frame.force = self.forces[index]
            if self.velocities is not None:
                frame.velocity = self.velocities[index]
            if self.time is not None:
                frame.time = self.time[index]
            yield frame

    def _handle_setting_box_force_velocity(self, frame, index):
        if self._boxes is not None:
            frame.box = Box(self._boxes[index])
        if self.velocities is not None:
            frame.velocity = self.velocities[index]
        if self.forces is not None:
            frame.force = self.forces[index]
        if self.time is not None:
            frame.time = self.time[index]

    def __getitem__(self, index):
        """return a view or copy of coordinates (follow numpy's rule)

        Examples
        --------
        >>> import pytraj as pt
        >>> import pytraj as pt
        >>> traj = pt.load_sample_data('tz2')[:]
        >>> assert isinstance(traj, pt.Trajectory)
        >>> # create mutable trajectory from TrajectoryIterator
        >>> t0 = traj[:]

        >>> # get a Frame view

        >>> # get a Trajetory view
        >>> t0 = traj[0:8:2]

        >>> # get a copy of Trajetory
        >>> t0 = traj[[0, 4, 6]]

        >>> # get a copy, keep only CA atoms
        >>> t0 = traj['@CA']

        >>> # get a copy, keep only CA atoms for 3 frames
        >>> t0 = traj[:3, '@CA']

        >>> # get a new stripped Frame
        >>> t0 = traj[0, '@CA']

        >>> atm = traj.top('@CA')
        >>> t1 = traj[atm]
        """
        if self.n_frames == 0:
            raise IndexError("Your Trajectory is empty, how can I index it?")

        if is_int(index):
            # traj[0]
            # return a single Frame as a view
            arr0 = self._xyz[index]
            frame = Frame(self.n_atoms, arr0, _as_ptr=True)
            self._handle_setting_box_force_velocity(frame, index)
            self._life_holder = frame
        else:
            # return a new Trajectory
            traj = self.__class__()
            atm = None
            arr0 = None

            if isinstance(index, (str, AtomMask)):
                # return a copy
                # traj['@CA']
                atm = self.top(index) if isinstance(index,
                                                    str) else index
                traj.top = self.top._modify_state_by_mask(atm)
                arr0 = self._xyz[:, atm.indices]
                if self.forces is not None:
                    forces = self.forces[:, atm.indices]
                    traj.forces = forces.copy()
                if self.velocities is not None:
                    velocities = self.velocities[:, atm.indices]
                    traj.velocities = velocities.copy()
                # make copy to create contigous memory block
                if self._boxes is not None:
                    # always make a copy in this case
                    traj._boxes = self._boxes.copy()
                if self.time is not None:
                    # always make a copy in this case
                    traj.time = self.time.copy()
                traj._xyz = arr0.copy()
            elif not isinstance(index, tuple):
                # might return a view or a copy
                # based on numpy array rule
                # traj.xyz[index]
                traj.top = self.top
                traj._xyz = self._xyz[index]
                if self._boxes is not None:
                    traj._boxes = self._boxes[index]
                if self.forces is not None:
                    forces = self.forces[index]
                    traj.forces = forces
                if self.velocities is not None:
                    velocities = self.velocities[index]
                    traj.velocities = velocities
                if self.time is not None:
                    traj.time = self.time[index]
            else:
                # is a tuple
                if len(index) == 1:
                    traj = self[index[0]]
                elif len(index) == 2 and is_int(index[0]) and isinstance(
                        index[1], str):
                    # traj[0, '@CA']: return a stripped Frame
                    frame = self[index[0]].copy()
                    # make AtomMask object
                    atm = self.top(index[1])
                    # create new Frame with AtomMask
                    self._life_holder = Frame(frame, atm)
                    return self._life_holder
                else:
                    self._life_holder = self[index[0]]
                    if isinstance(self._life_holder, Frame):
                        self._frame_holder = self._life_holder
                    traj = self._life_holder[index[1:]]
            self._life_holder = traj
        return self._life_holder

    def __setitem__(self, index, other):
        if self.n_frames == 0:
            raise IndexError("Your Trajectory is empty, how can I index it?")

        if other is None:
            raise ValueError("why bothering assign None?")
        if is_int(index):
            if hasattr(other, 'xyz') or isinstance(other, Frame):
                # traj[1] = frame
                self._xyz[index] = other.xyz
            else:
                # traj[1] = xyz
                # check shape?
                self._xyz[index] = other
        elif index == '*':
            # why need this?
            # traj.xyz = xyz
            # update all atoms, use fast version
            self._xyz[:] = other  # xyz
        elif isinstance(index, str):
            # update xyz for mask
            # traj['@CA'] = xyz
            atm = self.top(index)
            if isinstance(other, Trajectory):
                indices = atm.indices

                for i in range(self.n_frames):
                    for j, k in enumerate(indices):
                        self.xyz[i, k] = other.xyz[i, j]
            else:
                view3d = other
                int_view = atm.indices.astype('i4')
                # loop all frames
                for i in range(view3d.shape[0]):
                    self._xyz[:, int_view] = view3d[:]
        else:
            # really need this?
            # example: self[0, 0, 0] = 100.
            self._xyz[index] = other

    def append_xyz(self, xyz):
        '''append 3D numpy array

        Examples
        --------
        >>> import pytraj as pt
        >>> traj = pt.load_sample_data('tz2')
        >>> t0 = pt.Trajectory(top=traj.top)
        >>> t0.append_xyz(traj.xyz)
        >>> t0.n_frames
        10
        >>> t0.append_xyz(traj.xyz)
        >>> t0.n_frames
        20

        Notes
        -----
        This method is not well optimized for speed.
        '''
        # make sure 3D
        if xyz.ndim != 3:
            raise ValueError("ndim must be 3")

        if self.shape == (None, None, 3):
            self._xyz = xyz
        else:
            self._xyz = np.vstack((self._xyz, xyz))

    def _append_unitcells(self, box):
        '''append unitcells

        Examples
        --------
        >>> import pytraj as pt
        >>> traj = pt.datafiles.load_tz2_ortho()
        >>> traj2 =  pt.Trajectory(top=traj.top)
        >>> traj2._append_unitcells(traj.unitcells)
        >>> traj2.unitcells.shape
        (10, 6)

        >>> traj3 =  pt.Trajectory(top=traj.top)
        >>> clen, cangle = traj.unitcells[:, :3], traj.unitcells[:, 3:]
        >>> traj3._append_unitcells((clen, cangle))
        >>> traj3.unitcells.shape
        (10, 6)

        >>> traj4 = pt.Trajectory(top=traj.top)
        >>> traj4._append_unitcells(traj.unitcells)
        >>> traj4._append_unitcells(traj.unitcells)
        '''
        if isinstance(box, tuple):
            clen, cangle = box
            data = np.hstack((clen, cangle))
            if self._boxes is None:
                self._boxes = np.asarray([data], dtype='f8')
            else:
                self._boxes = np.vstack((self._boxes, data))

        else:
            if self._boxes is None:
                self._boxes = np.asarray(box, dtype='f8')
            else:
                self._boxes = np.vstack((self._boxes, box))

        if self._boxes.ndim == 3:
            shape = self._boxes.shape
            n_frames = int(shape[0] * shape[1] * shape[2] / 6)
            self._boxes = self._boxes.reshape((n_frames, 6))

    def append(self, other):
        """other: xyz, Frame, Trajectory, ...

        Examples
        --------
        >>> import pytraj as pt
        >>> import numpy as np
        >>> traj = pt.load_sample_data('tz2')[:]
        >>> t0 = pt.Trajectory(top=traj.top)
        >>> t0.n_frames
        0
        >>> f0 = traj[0]
        >>> t0.append(f0)
        >>> t0.n_frames
        1
        >>> t0.append(np.array([traj[3].xyz,]))
        >>> t0.n_frames
        2
        >>> t0.append(traj)
        >>> t0.n_frames
        12
        >>> t0.append(traj())
        >>> t0.n_frames
        22

        >>> t1 = pt.Trajectory(top=traj.top)
        >>> t1.append(traj)

        >>> t2 = pt.Trajectory(top=traj.top)
        >>> t2.append(traj.xyz)

        Notes
        -----
            - Can not append TrajectoryIterator object since we use Trajectory in TrajectoryIterator class
            - This method is not well optimized for speed.
        """
        if isinstance(other, Frame):
            arr0 = other.xyz.reshape((1, other.n_atoms, 3))
            barr = other.box.values.reshape((1, 6))
            if self._xyz is None:
                self._xyz = arr0.copy()
                self._boxes = barr
            else:
                self._xyz = np.vstack((self._xyz, arr0))
                self._boxes = np.vstack((self._boxes, barr))
        elif isinstance(other, np.ndarray) and other.ndim == 3:
            if self._xyz is None:
                self._xyz = other
                self._boxes = np.empty((other.shape[0], 6))
            else:
                self._xyz = np.vstack((self._xyz, other))

                fake_box_arr = np.empty((other.shape[0], 6))
                self._boxes = np.vstack((self._boxes, fake_box_arr))
        elif hasattr(other, 'n_frames') and hasattr(other, 'xyz'):
            # assume Trajectory-like object
            if self._xyz is None:
                self._xyz = other.xyz[:]
                self._boxes = other.unitcells
            else:
                self._xyz = np.vstack((self._xyz, other.xyz))
                self._boxes = np.vstack((self._boxes, other.unitcells))
        elif is_frame_iter(other):
            for frame in other:
                self.append(frame)
        else:
            # try to iterate to get Frame
            for frame in iterframe_master(other):
                self.append(frame)

    def __call__(self, *args, **kwd):
        '''shortcut of ``iterframe`` method

        Examples
        --------
        >>> import pytraj as pt
        >>> traj = pt.load_sample_data()[:]
        >>> for f in traj(0, 8, 2): pass
        >>> # same as
        >>> for f in traj.iterframe(0, 8, 2): pass

        '''
        return self.iterframe(*args, **kwd)

    def load(self, filename=''):
        '''This is for internal use. User should always use
        ``pytraj.load`` (or ``iterload``) method

        Examples
        --------
        >>> import pytraj as pt
        >>> from pytraj.testing import get_fn
        >>> fname, tname = get_fn('tz2')
        >>> traj = pt.Trajectory()
        >>> traj.top = pt.load_topology(tname)
        >>> traj.load(fname)
        >>> traj.n_atoms
        5293

        >>> # from list/tuple
        >>> traj.load([fname, fname])
        >>> traj.load((fname, fname))

        Notes
        -----
        It's better to use ``pytraj.load`` method

        >>> traj = pt.load(fname, tname)
        >>> traj.n_atoms
        5293
        '''
        if self.top is None or self.top.is_empty():
            raise RuntimeError('Must have a valid Topology')

        # always use self.top
        from pytraj import TrajectoryIterator
        ts = TrajectoryIterator()
        ts.top = self.top.copy()
        ts._load(filename)
        ts = ts[:]
        self._xyz = ts.xyz
        self.unitcells = ts.unitcells
        self.forces = ts.forces
        self.velocities = ts.velocities

    def autoimage(self, command=''):
        '''perform autoimage

        Returns
        -------
        self

        Examples
        --------
        >>> import pytraj as pt; from pytraj.testing import get_fn
        >>> t0 = pt.load(*get_fn('tz2'))
        >>> t0.top.has_box()
        True
        >>> t0 = t0.autoimage()
        '''
        from pytraj.analysis.c_action import c_action

        act = c_action.Action_AutoImage()
        act(command, self, top=self.top)
        return self

    def rotate(self, command=''):
        '''do rotation

        Returns
        -------
        self

        Examples
        --------
        >>> import pytraj as pt
        >>> traj = pt.load_sample_data('ala3')[:]
        >>> traj = traj.rotate('@CA x 20')
        '''
        from pytraj.analysis.c_action import c_action

        act = c_action.Action_Rotate()
        act(command, self, top=self.top)
        return self

    def translate(self, command=''):
        '''do translation

        Returns
        -------
        self

        Examples
        --------
        >>> import pytraj as pt
        >>> traj = pt.load_sample_data('ala3')[:]
        >>> traj = traj.translate('@CA x 1.2')
        '''
        from pytraj.analysis.c_action import c_action

        act = c_action.Action_Translate()
        act(command, self, top=self.top)
        return self

    def scale(self, command=''):
        '''do scaling

        Returns
        -------
        self

        Examples
        --------
        >>> import pytraj as pt
        >>> traj = pt.load_sample_data('ala3')[:]
        >>> traj = traj.scale('@CA x 1.2')
        '''
        from pytraj.analysis.c_action import c_action

        act = c_action.Action_Scale()
        act(command, self, top=self.top)
        return self

    def center(self, *args, **kwargs):
        '''do centering

        Returns
        -------
        self

        Examples
        --------
        >>> import pytraj as pt
        >>> traj = pt.load_sample_data('ala3')[:]
        >>> traj = traj.center('@CA origin')

        See also
        --------
        pytraj.all_actions.center
        '''
        from pytraj.all_actions import center
        return center(self, *args, **kwargs)

    def align_principal_axis(self, command=''):
        """align principal axis

        Examples
        --------
        >>> import pytraj as pt
        >>> traj = pt.load_sample_data('ala3')[:]
        >>> traj = traj.align_principal_axis()
        """
        from pytraj.analysis.c_action import c_action
        act = c_action.Action_Principal()

        command += " dorotation"
        act(command, self, top=self.top)
        return self

    def transform(self, commands, frame_indices=None):
        '''apply a series of cpptraj commands to trajectory

        Returns
        -------
        self

        Examples
        --------
        >>> import pytraj as pt
        >>> traj = pt.datafiles.load_tz2_ortho()[:]
        >>> traj.xyz[0, 0]
        array([ 15.55458927,  28.54844856,  17.18908691])
        >>> traj = traj.transform(['autoimage', 'center @CA origin', 'translate x 1.2'])
        >>> traj.xyz[0, 0]
        array([-1.19438073,  8.75046229, -1.82742397])

        >>> # which is similiar to below:
        >>> traj2 = pt.datafiles.load_tz2_ortho()[:]
        >>> traj2.xyz[0, 0] # before transforming
        array([ 15.55458927,  28.54844856,  17.18908691])
        >>> traj = traj2.autoimage().center('@CA origin').translate('x 1.2')

        >>> traj2.xyz[0, 0] # after transforming
        array([-1.19438073,  8.75046229, -1.82742397])
        '''
        from pytraj.analysis.c_action.actionlist import pipe
        fi = pipe(self, commands, frame_indices=frame_indices)

        for _ in fi:
            pass
        return self

    @property
    def unitcells(self):
        '''return 2D ndarray, shape=(n_frames, 6)

        Examples
        --------
        >>> import pytraj as pt
        >>> traj = pt.load_sample_data('tz2')[:]
        >>> traj.unitcells[0]
        array([ 35.26277966,  41.84554768,  36.16862953,  90.        ,
                90.        ,  90.        ])
        '''
        return self._boxes

    @unitcells.setter
    def unitcells(self, values):
        self._boxes = values

    def rmsfit(self, *args, **kwd):
        """do the fitting to reference Frame by rotation and translation

        Parameters
        ----------
        ref : {Frame, int}, default=None (first Frame)
            Reference
        mask : str or AtomMask object, default='*' (fit all atoms)

        Examples
        --------
        >>> traj.rmsfit(0) # fit to 1st frame # doctest: +SKIP
        >>> traj.rmsfit(-1, '@CA') # fit to last frame using @CA atoms # doctest: +SKIP

        Returns
        -------
        self

        Notes
        -----
        this is alias of superpose

        """
        return self.superpose(*args, **kwd)

    def superpose(self,
                  mask="*",
                  ref=None,
                  ref_mask='',
                  frame_indices=None,
                  mass=False):
        """do the fitting to reference Frame by rotation and translation

        Parameters
        ----------
        mask : str or AtomMask object, default='*' (fit all atoms)
        ref : {Frame object, int, str}, default=None
            Reference
        ref_mask : str, default ''
            if not given, use `mask`
            if given, use it
        mass : bool, default False
            if True, use mass-weighted
        frame_indices : array-like, default None, optional
            if not None, only do fitting for specific frames

        Returns
        -------
        self

        Examples
        --------
        >>> import pytraj as pt
        >>> from pytraj.testing import get_fn
        >>> traj = pt.load(*get_fn('tz2'))
        >>> traj = traj.superpose() # fit to 1st frame
        >>> traj = traj.superpose(ref=0) # fit to 1st frame, explitly specify
        >>> traj = traj.superpose(ref=-1, mask='@CA') # fit to last frame using @CA atoms
        """
        from pytraj.all_actions import superpose

        return superpose(
            self,
            mask=mask,
            ref=ref,
            ref_mask=ref_mask,
            mass=mass,
            frame_indices=frame_indices)

    def _allocate(self, n_frames, n_atoms):
        '''allocate (n_frames, n_atoms, 3) coordinates
        '''
        self._xyz = np.zeros((n_frames, n_atoms, 3), dtype='f8')

    def strip(self, mask):
        '''strip atoms with given mask

        Examples
        --------
        >>> import pytraj as pt
        >>> traj = pt.load_sample_data()[:]
        >>> traj.n_atoms
        34
        >>> t0 = traj.strip('!@CA') # keep only CA atoms
        >>> isinstance(t0, pt.Trajectory)
        True
        >>> t0.n_atoms
        3
        '''
        # AtomMask
        atm = self.top(mask)
        atm.invert_mask()
        self.top.strip(mask)

        self._xyz = np.ascontiguousarray(self._xyz[:, atm.indices])
        return self

    def iterframe(self,
                  start=0,
                  stop=None,
                  step=1,
                  mask=None,
                  autoimage=False,
                  frame_indices=None,
                  rmsfit=None):
        '''iterate trajectory with given frame_indices or given (start, stop, step)

        Parameters
        ----------
        start : int, default 0
        stop : {None, int}, default None
            if None, iterate to final frame
        step : int, default 1
        mask : {None, str}, default None
            if None, use all atoms. If not None, use given mask
        autoimage : bool, default False
            if True, perform autoimage for each frame
        rmsfit : {None, int, tuple}, default None
            if not None, perform superpose each Frame to to reference.
        frame_indices : {None, array-like}
            if not None, iterate trajectory for given indices. If frame_indices is given,
            (start, stop, step) will be ignored.


        Examples
        --------
        >>> import pytraj as pt
        >>> from pytraj.testing import get_fn
        >>> traj = pt.load(*get_fn('tz2'))
        >>> for frame in traj.iterframe(0, 8, 2): pass
        >>> for frame in traj.iterframe(0, 8, 2, autoimage=True): pass
        >>> # use negative index
        >>> traj.n_frames
        10
        >>> fi = traj.iterframe(0, -1, 2, autoimage=True)
        >>> fi.n_frames
        5
        >>> # mask is atom indices
        >>> fi = traj.iterframe(0, -1, 2, mask=range(100), autoimage=True)
        >>> fi.n_atoms
        100
        '''

        if mask is None:
            _top = self.top
        else:
            if isinstance(mask, str):
                mask = mask
                _top = self.top._get_new_from_mask(mask)
            else:
                mask = array_to_cpptraj_atommask(mask)
                _top = self.top._get_new_from_mask(mask)

        if rmsfit is not None:
            if isinstance(rmsfit, tuple):
                assert len(rmsfit) == 2, (
                    "rmsfit must be a tuple of one (frame,) "
                    "or two elements (frame, mask)")
            elif isinstance(rmsfit, (int, Frame)):
                rmsfit = (rmsfit, '*')
            else:
                raise ValueError("rmsfit must be a tuple or an integer")

            if is_int(rmsfit[0]):
                index = rmsfit[0]
                rmsfit = ([self[index], rmsfit[1]])

        # check how many frames will be calculated
        if frame_indices is None:
            start, stop, step = slice(start, stop, step).indices(self.n_frames)
            # make sure `range` return iterator
            indices = range(start, stop, step)
            n_frames = len(indices)
        else:
            # frame_indices is not None
            start, stop, step = None, None, None
            try:
                n_frames = len(frame_indices)
            except TypeError:
                # itertools.chain
                n_frames = None
            indices = frame_indices

        frame_iter_super = self._iterframe_indices(indices)

        return FrameIterator(
            frame_iter_super,
            original_top=self.top,
            new_top=_top,
            start=start,
            stop=stop,
            step=step,
            mask=mask,
            autoimage=autoimage,
            rmsfit=rmsfit,
            n_frames=n_frames,
            frame_indices=frame_indices,
            copy=False)

    @property
    def _estimated_GB(self):
        """esimated GB of data will be loaded to memory

        Examples
        --------
        >>> import pytraj as pt
        >>> traj = pt.load_sample_data('tz2')[:]
        >>> traj._estimated_GB
        0.0011830776929855347
        """
        return self.n_frames * self.n_atoms * 3 * 8 / (1024**3)

    @classmethod
    def from_iterable(cls, iterables, top=None):
        '''create a new Trajectory from iterable (produce Frame)

        Examples
        --------
        >>> import pytraj as pt
        >>> traj = pt.load_sample_data('tz2')
        >>> t0 = pt.Trajectory.from_iterable(traj(3, 8, 2))
        >>> from pytraj import pipe
        >>> fi = pipe(traj, ['autoimage', 'rms'])
        >>> t0 = pt.Trajectory.from_iterable(fi, top=traj.top)
        >>> t0.n_frames
        10
        >>> pt.radgyr(t0)
        array([ 18.90953437,  18.93564662,  18.85415458,  18.90994856,
                18.85884218,  18.88551081,  18.9364612 ,  18.89353463,
                18.91772124,  18.87070283])

        '''
        if top is None or top.is_empty():
            if hasattr(iterables, 'top'):
                top = iterables.top
            else:
                raise ValueError("must provide non-empty Topology")

        fa = Trajectory()
        fa.top = top

        if hasattr(iterables, 'n_frames'):
            _n_frames = iterables.n_frames
        else:
            try:
                _n_frames = len(iterables)
            except TypeError:
                _n_frames = None

        # faster
        if _n_frames is None:
            xyz = np.array([data.xyz.copy() for data in iterables])
            fa.xyz = xyz
        else:
            fa._allocate(_n_frames, fa.top.n_atoms)
            fa._boxes = np.empty((_n_frames, 6), dtype='f8')
            for index, frame in enumerate(iterables):
                fa._xyz[index] = frame.xyz
                fa._boxes[index] = frame.box.data
        return fa

    def __len__(self):
        return self.n_frames

    def __del__(self):
        self._xyz = None
        self._boxes = None

    def __add__(self, other):
        '''merge two trajectories together. Order matter.

        Notes
        -----
        this is convenient method and it is not really optimized for memory and speed.

        Examples
        --------
        >>> import pytraj as pt
        >>> traj1 = pt.datafiles.load_ala3()[:1]
        >>> traj2 = pt.datafiles.load_tz2_ortho()[:1]
        >>> traj3 = traj1 + traj2
        >>> traj1.n_atoms
        34
        >>> traj2.n_atoms
        5293
        >>> traj3.n_atoms
        5327
        >>> traj1.xyz[0, 0] == traj3.xyz[0, 0]
        array([ True,  True,  True], dtype=bool)
        >>> traj2.xyz[0, -1] == traj3.xyz[0, -1]
        array([ True,  True,  True], dtype=bool)

        See also
        --------
        pytraj.tools.merge_trajs
        '''
        if self.n_frames != other.n_frames:
            raise ValueError('two trajs must have the same n_frames')

        traj = self.__class__()
        traj._allocate(self.n_frames, self.n_atoms + other.n_atoms)
        traj.top = self.top + other.top

        for f1, f2, frame in zip(self, other, traj):
            frame.xyz[:] = np.vstack((f1.xyz, f2.xyz))
        return traj

    @property
    def crdinfo(self):
        return {'has_velocity': self.velocities is not None,
         'has_temperature': False,
         'has_time': self.time is not None,
         'has_force': self.forces is not None,
         'has_box': self.unitcells is not None,
         'has_replcica_dims': False,
         'n_frames': self.n_frames,
         'n_atoms': self.top.n_atoms,
         'box_type': self.unitcells is not None and self.top.box.type or None}
