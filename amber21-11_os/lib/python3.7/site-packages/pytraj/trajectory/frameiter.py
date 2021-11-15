from __future__ import absolute_import
from pytraj.trajectory.frame import Frame
from .shared_methods import iterframe_master

__all__ = ['iterframe', 'iterchunk', 'FrameIterator']


def iterframe(traj, *args, **kwd):
    """create frame iterator with given indices, mask or some iter_options

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> for frame in pt.iterframe(traj, 0, 8, 2): pass
    >>> for frame in pt.iterframe(traj, 4, mask='@CA'): print(frame)
    <Frame with 12 atoms>
    <Frame with 12 atoms>
    <Frame with 12 atoms>
    <Frame with 12 atoms>
    <Frame with 12 atoms>
    <Frame with 12 atoms>

    # create frame iterator for given indices
    >>> for frame in pt.iterframe(traj, frame_indices=[0, 7, 3]): print(frame)
    <Frame with 5293 atoms>
    <Frame with 5293 atoms>
    <Frame with 5293 atoms>


    >>> fi = pt.iterframe(traj)
    >>> # iterframe its self
    >>> fi = pt.iterframe(fi)

    See also
    --------
    pytraj.TrajectoryIterator.iterframe
    """
    if hasattr(traj, 'iterframe'):
        return traj.iterframe(*args, **kwd)
    else:
        return iterframe_master(traj)


def iterchunk(traj, *args, **kwd):
    """iterate ``traj`` by chunk

    Parameters
    ----------
    traj : TrajectoryIterator
    chunksize : int
        the number of frames in each chunk
    start : int, default 0
        start frame to iterate
    start : int, default -1 (last frame)
        stop frame
    autoimage : bool, default False
        if True, do autoimage for chunk

    Return
    ------
    pytraj.Trajectory, n_frames=chunksize
        The final chunk might not have the n_frames=chunksize

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> for frame in pt.iterchunk(traj, 4): pass
    >>> for frame in pt.iterchunk(traj, chunksize=4, start=2): pass
    >>> for frame in pt.iterchunk(traj, chunksize=4, start=2, stop=9): pass
    >>> for frame in pt.iterchunk(traj, chunksize=4, autoimage=True): pass

    See also
    --------
    pytraj.TrajectoryIterator.iterchunk
    """
    return traj.iterchunk(*args, **kwd)


class FrameIterator(object):
    """
    create this class to hold all iterating information. This class is for internal use.

    Parameters
    ----------
    top : new Topology
    original_top : original Topology
    start, stop, step : int
    mask : str or None, default None (all atoms)
        only take atom with given mask
    frame_indices : iterable, default: None
        if frame_indices is not None: ignore (start, stop, step)
    autoimage : bool, default: False
        if autoimage, perform autoimage
    rmsfit : int or a tuple, default False
        if rmsfit, perform rms fit to reference. If ``rmsfit`` is an integer, perform
        rms fit to indicated frame for all atoms. If ``rmsfit`` is a tuple, perform rmsfit
        to given frame with given mask. if both ``autoimage`` and ``rmsfit`` are specified,
        do ``autoimage`` first.
    n_frames : total number of frame. read-only
    copy : bool, defaul: True
        if True, always make a copy of Frame when iterating.

    Notes
    -----
    if 'autoimage' and 'rmsfit', reference frame is also autoimaged

    Examples
    --------
    >>> # short cut:
    >>> # create FrameIterator with start=0, stop=8, step=2
    >>> import pytraj as pt
    >>> traj = pt.load_sample_data('tz2')
    >>> fi = traj(0, 8, 2)
    >>> # perform radgyr calculation with FrameIterator
    >>> pt.radgyr(traj(0, 8, 2))
    array([ 18.91114428,  18.84969884,  18.8568644 ,  18.9430491 ])

    >>> # create FrameIterator with start, stop, step = 0, 8, 2
    >>> # autoimage=False, rmsfit=False
    >>> fi = traj.iterframe(0, 8, 2)

    >>> # create FrameIterator with start, stop, step = 2, 8, 1
    >>> # autoimage=False, rmsfit=False
    >>> fi = traj.iterframe(2, 8)

    >>> # create FrameIterator with start, stop, step = 2, 8, 1
    >>> # autoimage=False, rmsfit=False, mask='@CA'
    >>> fi = traj.iterframe(2, 8, mask='@CA')

    >>> # create FrameIterator with start, stop, step = 2, 8, 1
    >>> # autoimage=True, rmsfit=False, mask='@CA'
    >>> for frame in traj.iterframe(2, 8, autoimage=True, mask='@CA'): print(frame)
    <Frame with 12 atoms>
    <Frame with 12 atoms>
    <Frame with 12 atoms>
    <Frame with 12 atoms>
    <Frame with 12 atoms>
    <Frame with 12 atoms>

    >>> # with rmsfit
    >>> for frame in traj.iterframe(2, 8, autoimage=True, rmsfit=0): pass
    >>> for frame in traj.iterframe(2, 8, autoimage=True, rmsfit=(1, '!@H=')): pass

    >>> # rmsfit
    >>> fi = traj.iterframe(2, 8, rmsfit=(0, '@CA'))
    >>> fi.n_frames
    6

    >>> fi = traj.iterframe(2, 8, mask='@1,2,3,4,5')
    >>> fi.n_atoms
    5

    >>> # make copy of each Frame
    >>> fi = traj.iterframe(2, 8, mask='@1,2,3,4,5', copy=True)

    >>> fi = traj.iterframe(2, 8, rmsfit=3)
    >>> fi = traj.iterframe(2, 8, mask='@1,2,3,4,5')

    >>> # explit use copy=True to give different Frame with list
    # >>> fi = traj.iterframe(2, 8)
    # >>> fi.copy = True
    # >>> pt.radgyr(list(fi), top=traj.top)
    # array([ 18.84969884,  18.90449256,  18.8568644 ,  18.88917208,
    #         18.9430491 ,  18.88878079])
    """

    def __init__(self,
                 fi_generator,
                 original_top=None,
                 new_top=None,
                 start=0,
                 stop=-1,
                 step=1,
                 mask="",
                 autoimage=False,
                 rmsfit=None,
                 n_frames=None,
                 copy=True,
                 frame_indices=None):
        self.top = new_top
        self.original_top = original_top
        self.frame_iter = fi_generator
        self.start = start
        self.stop = stop
        self.step = step
        self.mask = mask
        self.autoimage = autoimage
        self.rmsfit = rmsfit
        # use `copy_frame` for TrajectoryIterator
        self._n_frames = n_frames
        self.copy = copy
        self.frame_indices = frame_indices

    @property
    def n_frames(self):
        return self._n_frames

    @property
    def n_atoms(self):
        return self.top.n_atoms

    @property
    def __name__(self):
        '''for inspecting
        '''
        return "FrameIterator"

    def __str__(self):
        root_msg = '<FrameIterator with '
        root_msg2 = 'start=%s, stop=%s, step=%s, n_frames=%s, \n' % (
            self.start, self.stop, self.step, self.n_frames)
        root_msg3 = 'frame_indices=%s, \n' % self.frame_indices

        more_msg = 'mask=%s, autoimage=%s, rmsfit=%s, copy=%s> \n' % (
            self.mask, self.autoimage, self.rmsfit, self.copy)
        return "".join((root_msg, root_msg2, root_msg3, more_msg))

    def __repr__(self):
        return self.__str__()

    def save(self, filename='', overwrite=False, options='', *args, **kwd):
        '''save to different file format.

        Notes
        -----
        FrameIterator will be exhausted since this is an iterator.

        Examples
        --------
        >>> import pytraj as pt
        >>> traj = pt.load_sample_data('tz2')
        >>> fi = traj(2, 8, 2, mask='@CA')
        >>> fi.save('output/test.nc', overwrite=True)
        >>> # short version
        >>> traj(2, 8, 2, mask='@CA').save('output/test.nc', overwrite=True)
        '''
        from pytraj.io import write_traj
        write_traj(
            filename=filename,
            traj=self,
            frame_indices=None,
            overwrite=overwrite,
            options=options,
            *args,
            **kwd)

    def __iter__(self):
        from pytraj.analysis.c_action import c_action
        # do not import c_action in the top to avoid circular importing
        if self.autoimage:
            image_act = c_action.Action_AutoImage()
            image_act.read_input("", top=self.original_top)
            image_act.setup(self.original_top)
        if self.rmsfit is not None:
            ref, mask_for_rmsfit = self.rmsfit
            need_align = True
            if self.autoimage:
                # need to do autoimage for ref too
                # make a copy to avoid changing ref
                ref = ref.copy()
                image_act.compute(ref)
            rmsd_act = c_action.Action_Rmsd()
            rmsd_act.read_input(mask_for_rmsfit, top=self.original_top)
            rmsd_act.setup(self.original_top)
            # creat first frame to trick cpptraj to align to this.
            rmsd_act.compute(ref)
        else:
            need_align = False
            ref, mask_for_rmsfit = None, None

        if self.mask is not None:
            mask = self.mask
            atm = self.original_top(mask)

        for frame0 in self.frame_iter:
            if self.copy:
                # use copy for TrajectoryIterator
                # so [f for f in traj()] will return a list of different
                # frames
                frame = frame0.copy()
            else:
                frame = frame0
            if self.autoimage:
                # from pytraj.c_action.c_action import Action_AutoImage
                # Action_AutoImage()("", frame, self.top)
                image_act.compute(frame)
            if need_align:
                # trick cpptraj to fit to 1st frame (=ref)
                rmsd_act.compute(frame)
            if self.mask is not None:
                frame2 = Frame(frame, atm)
                yield frame2
            else:
                yield frame
