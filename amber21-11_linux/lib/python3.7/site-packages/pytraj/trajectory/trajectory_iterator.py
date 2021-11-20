"""out-of-core TrajectoryIterator
"""
from __future__ import absolute_import
import os
import re
from glob import glob
import numpy as np
from .c_traj.c_trajectory import TrajectoryCpptraj
from .shared_trajectory import SharedTrajectory
from ..topology.topology import Topology
from .frame import Frame
from ..utils import is_int
from ..utils.cyutils import get_positive_idx
from .frameiter import FrameIterator
from ..utils.get_common_objects import _load_Topology
from ..utils import split_range
from ..utils.convert import array_to_cpptraj_atommask
from pytraj.utils.get_common_objects import get_reference

__all__ = [
    'TrajectoryIterator',
]


# tryint, alphanum_key, sort_filename_by_number are adapted from
# http://nedbatchelder.com/blog/200712/human_sorting.html
def tryint(s):
    try:
        return int(s)
    except:
        return s


def alphanum_key(s):
    """ Turn a string into a list of string and number chunks.
        "z23a" -> ["z", 23, "a"]
    """
    return [tryint(c) for c in re.split('([0-9]+)', s)]


def sort_filename_by_number(filelist):
    """ Sort the given list in the way that humans expect.
    """
    return sorted(filelist, key=alphanum_key)


def _make_frame_slices(n_files, original_frame_slice):
    '''
    >>> _make_frame_slices(2, (0, -1))
    [(0, -1), (0, -1)]

    >>> _make_frame_slices(3, [(0, -1), (0, -2)],)
    [(0, -1), (0, -2), (0, -1, 1)]

    >> # raise
    >>> _make_frame_slices(3, None)
    Traceback (most recent call last):
        ...
    ValueError: must be a tuple of integer values or a list of tuple of integer values

    >>> _make_frame_slices(2, [(0, -1), (0, -2), (0, -1, 1)])
    Traceback (most recent call last):
        ...
    ValueError: len of frame_slice tuple-list must be smaller or equal number of files
    '''
    if isinstance(original_frame_slice, tuple):
        return [original_frame_slice for i in range(n_files)]
    elif isinstance(original_frame_slice, list):
        fs_len = len(original_frame_slice)
        if fs_len < n_files:
            new_list = original_frame_slice[:] + [
                (0, -1, 1) for _ in range(fs_len, n_files)
            ]
        elif fs_len == n_files:
            new_list = original_frame_slice
        else:
            raise ValueError(
                "len of frame_slice tuple-list must be smaller or equal number of files"
            )
        return new_list
    else:
        raise ValueError(
            "must be a tuple of integer values or a list of tuple of integer values"
        )


class TrajectoryIterator(TrajectoryCpptraj, SharedTrajectory):
    '''out-of-core trajectory holder.

    Examples
    --------
    >>> import pytraj as pt
    >>> from pytraj.testing import get_fn
    >>> traj_name, top_name = get_fn('tz2')
    >>> traj = pt.TrajectoryIterator(traj_name, top_name)

    >>> # user should always use :method:`pytraj.iterload` to load TrajectoryIterator
    >>> traj = pt.iterload(['remd.x.000', 'remd.x.001'], 'input.parm7') # doctest: +SKIP

    Notes
    -----
    It's a bit tricky to pickle this class. As default, new TrajectoryIterator will
    use original trj filename and top filename.
    '''

    def __init__(self, filename=None, top=None, *args, **kwd):
        self._force_load = False
        # use self._chunk to store `chunk` in iterchunk
        # to deallocate memory
        self._chunk = None
        # only allow to load <= 1 GB
        self._size_limit_in_GB = 1
        super(TrajectoryIterator, self).__init__()

        if not top:
            self.top = Topology()
        elif isinstance(top, str):
            self.top = _load_Topology(top)
        elif isinstance(top, Topology):
            self.top = top.copy()
        else:
            raise ValueError("Topology must be None/string/Topology")

        self._frame_slice_list = []

        if filename:
            if self.top.is_empty():
                raise ValueError(
                    'First argument is always a trajectory filename'
                    ' or a list of filenames'
                    'must have a non-empty Topology')
            self._load(filename, self.top, *args, **kwd)

        self.__dict__.update({
            '_top_filename': self.top.filename,
            'filelist': self.filelist,
            '_frame_slice_list': self._frame_slice_list,
        })

    def __setstate__(self, state):
        self.__dict__ = state
        self.top = _load_Topology(state['_top_filename'])
        self._load(state['filelist'], frame_slice=state['_frame_slice_list'])

        _transform_commands = state.get('_transform_commands')
        if _transform_commands:
            self._transform_commands = _transform_commands
            self._reset_transformation()

    def __getstate__(self):
        '''

        Examples
        --------
        >>> import pytraj as pt
        >>> traj = pt.datafiles.load_ala3()
        >>> # pickle by reloading Topology from filename
        >>> pt.to_pickle(traj, 'output/test.pk')

        >>> # pickle by rebuilding Topology
        >>> pt.to_pickle(traj, 'output/test.pk')
        '''
        # slow
        # Topology is pickable
        sdict = self.__dict__
        sdict['_transform_commands'] = self._transform_commands
        return sdict

    def __iter__(self):
        '''do not make a frame copy here
        '''
        for frame in super(TrajectoryIterator, self).__iter__():
            yield frame

    def copy(self):
        '''return a deep copy. Use this method with care since the copied traj just reuse
        the filenames
        '''
        other = self.__class__()
        other.top = self.top.copy()

        for fname, frame_slice in zip(self.filelist, self._frame_slice_list):
            other._load(fname, frame_slice=frame_slice)
        return other

    def _load(self,
              filename=None,
              top=None,
              frame_slice=(0, -1, 1),
              stride=None):
        """load trajectory/trajectories from filename/filenames
        with a single frame_slice or a list of frame_slice

        Notes
        -----
        if stride is not None, frame_slice will be ignored
        """
        if not top:
            top_ = self.top
        else:
            top_ = top

        frame_slice_ = frame_slice if stride is None else (0, -1, stride)

        if isinstance(filename, str) and os.path.exists(filename):
            super(TrajectoryIterator, self)._load(filename, top_, frame_slice_)
            self._frame_slice_list.append(frame_slice_)
        elif isinstance(filename,
                        str) and not os.path.exists(filename):

            flist = sort_filename_by_number(glob(filename))
            if not flist:
                if "\\" in filename:
                    # filename with white space
                    flist = [filename]
                if not flist:
                    raise ValueError(
                        "Must provie a filename or list of filenames or file pattern"
                    )
                frame_slice_ = [
                    (0, -1, stride),
                ] * len(flist) if stride is not None else frame_slice_
            self._load(flist, top=top, frame_slice=frame_slice_)
        elif isinstance(filename, (list, tuple)):
            flist = filename

            if stride is None:
                full_frame_slice = _make_frame_slices(len(flist), frame_slice)
            else:
                full_frame_slice = [
                    (0, -1, stride),
                ] * len(flist)

            for fname, fslice in zip(flist, full_frame_slice):
                self._frame_slice_list.append(frame_slice)
                super(TrajectoryIterator, self)._load(
                    fname, top_, frame_slice=fslice)
        else:
            raise ValueError("filename must a string or a list of string")

    @property
    def topology(self):
        """alias of traj.top

        Examples
        --------
        >>> import pytraj as pt
        >>> from pytraj.testing import get_fn
        >>> fname, tname = get_fn('ala3')
        >>> traj = pt.iterload(fname, tname)
        >>> traj.topology
        <Topology: 34 atoms, 3 residues, 1 mols, non-PBC>
        >>> new_traj = pt.TrajectoryIterator()
        >>> new_traj.topology = traj.topology
        >>> new_traj._load(traj.filename)
        """
        return self.top

    @topology.setter
    def topology(self, newtop):
        self.top = newtop

    @property
    def _estimated_GB(self):
        """esimated GB of data will be loaded to memory
        """
        return self.n_frames * self.n_atoms * 3 * 8 / (1024**3)

    @property
    def xyz(self):
        '''return 3D array of coordinates'''
        size_in_GB = self._estimated_GB
        # check if larger than size_limit_in_GB
        if size_in_GB > self._size_limit_in_GB and not self._force_load:
            raise MemoryError(
                "you are loading %s GB, larger than size_limit %s GB. "
                "Please increase traj._size_limit_in_GB or set traj._force_load to True"
                % (size_in_GB, self._size_limit_in_GB))
        return super(TrajectoryIterator, self).xyz

    def iterframe(self,
                  start=0,
                  stop=None,
                  step=1,
                  mask=None,
                  autoimage=False,
                  rmsfit=None,
                  copy=False,
                  frame_indices=None):
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
        >>> traj = pt.load_sample_data('tz2')
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
            top_ = self.top
        else:
            if isinstance(mask, str):
                mask = mask
                top_ = self.top._get_new_from_mask(mask)
            else:
                mask = array_to_cpptraj_atommask(mask)
                top_ = self.top._get_new_from_mask(mask)

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
            # only check if does not have frame_indices
            if stop is None or stop >= self.n_frames:
                stop = self.n_frames
            elif stop < 0:
                stop = get_positive_idx(stop, self.n_frames)
            n_frames = len(range(start, stop, step))
            frame_iter_super = super(TrajectoryIterator, self).iterframe(
                start, stop, step)
        else:
            stop = None
            start = None
            step = None
            try:
                n_frames = len(frame_indices)
            except TypeError:
                # itertools.chain
                n_frames = None
            frame_iter_super = super(TrajectoryIterator,
                                     self)._iterframe_indices(frame_indices)

        return FrameIterator(
            frame_iter_super,
            original_top=self.top,
            new_top=top_,
            start=start,
            stop=stop,
            step=step,
            mask=mask,
            autoimage=autoimage,
            rmsfit=rmsfit,
            n_frames=n_frames,
            copy=copy,
            frame_indices=frame_indices)

    def iterchunk(self,
                  chunksize=2,
                  start=0,
                  stop=-1,
                  autoimage=False,
                  rmsfit=None):
        """iterate trajectory by chunk

        Parameters
        ----------
        chunk : int, default=2
            size of each chunk. Notes: final chunk's size might be changed
        start : int, default=0 (first frame)
        stop : int, default=-1 (last frame)
        autoimage : bool, default=False
        rmsfit : None | tuple/list of (reference frame, mask)

        Examples
        --------
        >>> import pytraj as pt
        >>> traj = pt.load_sample_data('tz2')
        >>> ref = traj[3]
        >>> for chunk in traj.iterchunk(3, autoimage=True, rmsfit=(ref, '@CA')): pass

        Notes
        -----
        if using 'autoimage` with reference frame for rms-fit, make sure to `autoimage`
        ref first
        """
        if rmsfit is not None:
            ref, mask_for_rmsfit = rmsfit
            need_align = True
        else:
            need_align = False
            ref, mask_for_rmsfit = None, None

        for chunk in super(TrajectoryIterator, self).iterchunk(
                chunksize, start, stop):
            # always perform autoimage before doing fitting
            # chunk is `Trajectory` object, having very fast `autoimage` and
            # `rmsfit` methods
            if autoimage:
                chunk.autoimage()
            if need_align:
                chunk.superpose(ref=ref, mask=mask_for_rmsfit)
            # free memory
            # if not, memory will be quicly accumulated.
            if self._chunk:
                self._chunk.__del__()

            self._chunk = chunk
            yield self._chunk

    @property
    def filename(self):
        '''return 1st filename in filelist. For testing only
        '''
        return self.filelist[0]

    @property
    def shape(self):
        '''(n_frames, n_atoms, 3)

        >>> import pytraj as pt
        >>> traj = pt.datafiles.load_tz2_ortho()
        >>> traj.shape
        (10, 5293, 3)
        '''
        return (self.n_frames, self.n_atoms, 3)

    def superpose(self, mask='*', ref=None, ref_mask='', mass=False):
        """register to superpose to reference frame when iterating. 
        To turn off superposing, set traj._being_transformed = False

        Notes
        -----
        This method is different from ``superpose`` in pytraj.Trajectory.
        It does not change the coordinates of TrajectoryCpptraj/TrajectoryIterator itself but 
        changing the coordinates of copied Frame.

        This method is mainly for NGLView in Jupyter notebook, to view out-of-core data.
        It's good to do translation and rotation on the fly.


        Examples
        --------
        >>> import pytraj as pt
        >>> traj = pt.datafiles.load_tz2()
        >>> isinstance(traj, pt.TrajectoryIterator)
        True
        >>> traj[0].xyz[0]
        array([-1.88900006,  9.1590004 ,  7.56899977])

        >>> # turn on superpose
        >>> _ = traj.superpose(ref=-1, mask='@CA')
        >>> traj[0].xyz[0]
        array([ 6.97324167,  8.82901548,  1.31844696])

        >>> # turn off superpose
        >>> traj._being_transformed = False
        >>> traj[0].xyz[0]
        array([-1.88900006,  9.1590004 ,  7.56899977])

        Examples for NGLView::

            import pytraj as pt, nglview as nv
            traj = pt.datafiles.load_tz2()
            traj.superpose(ref=0, mask='@CA')
            view = nv.show_pytraj(traj)
            view
        """
        ref = get_reference(self, ref)
        super(TrajectoryIterator, self).superpose(
            mask=mask, ref=ref, ref_mask=ref_mask, mass=mass)
        return self

    def _split_iterators(self,
                         n_chunks=1,
                         start=0,
                         stop=-1,
                         step=1,
                         mask=None,
                         autoimage=False,
                         rmsfit=None,
                         rank=0,
                         func=None):
        """simple splitting `self` to n_chunks FrameIterator objects

        Parameters
        ----------
        n_chunks : int, default 1
            number of chunks
        start, stop, step: int, default (0, -1, 1)
        mask : {None, str}, default None
            if given, iterate only atom in `mask`
        autoimage : bool, default False
            if True, do autoimage
        rmsfit : {None, Frame, int}, default None
            if given, perform rmsfit to reference. If autoimage=True and rmsfit is given, always perform autoimage first
        rank : rank
        func : {None, function}, default None
            if given, coordinates will be modified by this func.

        Examples
        --------
        >>> import pytraj as pt
        >>> traj = pt.load_sample_data('tz2')
        >>> list(traj._split_iterators(n_chunks=4, mask='@CA'))
        [<Frame with 12 atoms>, <Frame with 12 atoms>]
        >>> isinstance(traj._split_iterators(n_chunks=4, mask='@CA', rank=-1), list)
        True

        >>> # reset stop value to max n_framaes if this number looks 'weird'
        >>> fi = traj._split_iterators(n_chunks=4, mask='@CA', stop=-100)
        >>> fi = traj._split_iterators(n_chunks=4, mask='@CA', stop=traj.n_frames+100)
        """

        assert 0 <= start <= self.n_frames, "0 <= start <= self.n_frames"

        # if stop <= 0 or stop > self.n_frames:
        if not (0 < stop <= self.n_frames):
            stop = self.n_frames

        if rank >= 0:
            _start, _stop = split_range(
                n_chunks=n_chunks, start=start, stop=stop)[rank]
            return self.iterframe(
                start=_start,
                stop=_stop,
                step=step,
                mask=mask,
                autoimage=autoimage,
                rmsfit=rmsfit)
        else:
            list_of_iterators = []
            for (_start, _stop) in split_range(
                    n_chunks=n_chunks, start=start, stop=stop):
                list_of_iterators.append(
                    self.iterframe(
                        start=_start,
                        stop=_stop,
                        step=step,
                        mask=mask,
                        autoimage=autoimage,
                        rmsfit=rmsfit))
            return list_of_iterators

    @property
    def temperatures(self):
        """return 1D array of temperatures
        """
        return np.array([frame.temperature for frame in self])

    def at(self, index):
        """same as traj[index]
        """
        return self[index]

    def strip(self, mask):
        from pytraj.trajectory.stripped_trajectory import StrippedTrajectoryIterator
        return StrippedTrajectoryIterator(self, mask)
