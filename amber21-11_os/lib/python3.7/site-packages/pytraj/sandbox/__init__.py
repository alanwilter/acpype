'''sandbox for Julialang and other stuff.
Codes in this module might not be tested and they might be not correct, at all.
'''
import numpy as np
from functools import wraps
from pytraj.analysis.c_action import c_action


def take(traj, indices):
    return traj[indices]


def itake(traj, indices):
    return traj.iterframe(frame_indices=indices)


def get_top(traj):
    return traj.top


def set_top(traj, top):
    traj.top = top


def translate(traj, mask, anchor, to_point=None):
    '''translate a group of atoms in mask by moving anchor mask to given point

    Returns
    -------
    updated traj
    '''
    if to_point is None:
        to_point = [0.0, 0.0, 0.0]
    indices = traj.top.select(mask)

    for frame in traj:
        diff = np.asarray(to_point) - np.asarray(anchor)
        frame.xyz[indices] += diff
    return traj


def write_traj(filename, traj=None, mode='', frame_indices=None):
    '''
    Parameters
    ----------

    Notes
    -----
    cpptraj will detect file format based on extension for writting.


    Examples
    --------
    '''

    command = ' '.join((filename, mode))
    fi = traj if frame_indices is None else traj.iterframe(
        frame_indices=frame_indices)

    act = c_action.Action_Outtraj()
    _top = traj.top
    act(command, fi, top=_top)


def write_mask(traj, command):
    """very simple wrapping.
    """
    act = c_action.Action_Mask()
    act(command, traj, top=traj.top)


class TrajectoryWriter:
    # give wrong n_atoms for last frame. Why?
    '''
    Examples
    --------

    >>> with TrajectoryWriter('test.pdb', mode='model', top=top) trajout:
    >>>     for frame in traj:
    >>>         trajout.write(frame)

    '''

    def __init__(self, filename, mode='', top=None):
        from pytraj.analysis.c_action.c_action import Action_Outtraj
        self._outtraj = Action_Outtraj()
        command = ' '.join((filename, mode))
        self._outtraj.read_input(command, top=top)
        self._outtraj.process(top)

    def write(self, frame):
        self._outtraj.do_action(frame)

    def __enter__(self):
        return self

    def __exit__(self, *args):
        pass


# for testing
from pytraj.utils.get_common_objects import (get_data_from_dtype, get_topology,
                                             get_reference, get_fiterator)
from pytraj.analysis.c_action import c_action
from pytraj.datasets import CpptrajDatasetList
from pytraj.utils.convert import array_to_cpptraj_atommask


def _dispatch_traj_ref_top_frame_indices(f):
    @wraps(f)
    def inner(*args, **kwd):
        args = list(args)
        traj = kwd.get('traj', args[0])
        frame_indices = kwd.get('frame_indices')
        ref = kwd.get('ref')
        top = kwd.get('top')

        if 'mask' in kwd.keys():
            mask = kwd.get('mask')
        else:
            mask = args[1]

        # overwrite
        kwd['top'] = get_topology(traj, top)
        if ref is not None:
            kwd['ref'] = get_reference(traj, ref)
        if 'traj' in kwd.keys():
            kwd['traj'] = get_fiterator(traj, frame_indices)
        else:
            args[0] = get_fiterator(traj, frame_indices)
        if not isinstance(mask, str):
            mask = array_to_cpptraj_atommask(mask)
        if 'mask' in kwd.keys():
            kwd['mask'] = mask
        else:
            args[1] = mask
        return f(*args, **kwd)

    return inner


@_dispatch_traj_ref_top_frame_indices
def _toy_radgyr(traj,
                mask="",
                top=None,
                dtype='ndarray',
                nomax=True,
                frame_indices=None):
    '''
    Examples
    --------
    >>> import pytraj as pt
    >>> pt.mindist(traj, '@CA @H')
    '''
    act = c_action.Action_Radgyr()
    dslist = CpptrajDatasetList()

    _nomax = 'nomax' if nomax else ''
    mask = ' '.join((mask, _nomax))

    act(mask, traj, top=top, dslist=dslist)
    return get_data_from_dtype(dslist, dtype=dtype)


def calc_linear_interaction_energy(traj=None,
                                   mask="",
                                   top=None,
                                   dtype='dataset',
                                   frame_indices=None,
                                   *args,
                                   **kwd):
    command = mask
    act = c_action.Action_LIE()

    dslist = CpptrajDatasetList()
    act(command, traj, top=top, dslist=dslist, *args, **kwd)
    return get_data_from_dtype(dslist, dtype)


# alias
calc_LIE = calc_linear_interaction_energy


def _dbscan(traj=None,
            mask='*',
            minpoints=None,
            epsilon=0.,
            sievetoframe=False,
            random_sieveseed=1,
            kdist=None,
            kfile=None,
            sieve=1,
            metric='rms',
            top=None,
            options=''):
    '''perform clustering and return cluster index for each frame

    Parameters
    ----------
    traj : Trajectory-like or iterable that produces Frame
    mask : str, default: * (all atoms)
    n_clusters: int, default: 10
    random_point : bool, default: True
    maxit : int, default: 100
        max iterations
    metric : str, {'rms', 'dme'}
        distance metric
    top : Topology, optional, default: None
        only need to provide this Topology if ``traj`` does not have one
    options : option to save data to files.


    Returns
    -------
    1D numpy array of frame indices
    '''

    # don't need to get_topology
    _top = get_topology(traj, top)
    _clusters = 'dbscan minpoints ' + str(minpoints)
    _mask = mask
    _epsilon = 'epsilon ' + str(epsilon)
    _sievetoframe = 'sievetoframe' if sievetoframe else ''
    _kdist = 'kdist' + str(kdist) if kdist is not None else ''
    _kfile = kfile if kfile is not None else ''
    _sieve = 'sieve ' + str(sieve)
    _metric = metric
    _random_sieveseed = 'random ' + str(random_sieveseed)
    _output = options
    command = ' '.join((_clusters, _epsilon, _sievetoframe, _kdist, _sieve,
                        _kfile, _metric, _random_sieveseed, _mask, _output))
    return _cluster(traj, command, top=_top, dtype='ndarray')


#@_super_dispatch()


def calc_density(traj=None,
                 command="",
                 top=None,
                 dtype='ndarray',
                 frame_indices=None):
    # NOTE: trick cpptraj to write to file first and the reload
    '''

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> pt.density(traj, 'charge')
    '''

    with tempfolder():

        def _calc_density(traj, command):
            # TODO: update this method if cpptraj save data to
            # CpptrajDatasetList
            dflist = DataFileList()

            tmp_filename = "tmp_pytraj_out.txt"
            command = "out " + tmp_filename + " " + command
            act = c_action.Action_Density()
            act(command, traj, top=top, dflist=dflist)
            act.post_process()
            dflist.write_all_datafiles()
            absolute_path_tmp = os.path.abspath(tmp_filename)
            return absolute_path_tmp

        dslist = CpptrajDatasetList()
        fname = _calc_density(traj, command)
        dslist.read_data(fname)
        return get_data_from_dtype(dslist, dtype)


def _rotdif(arr,
            nvecs=1000,
            rvecin=None,
            rseed=80531,
            order=2,
            ncorr=-1,
            tol=1E-6,
            d0=0.03,
            nmesh=2,
            dt=0.002,
            ti=0.0,
            tf=-1,
            itmax=-1.,
            dtype='ndarray'):
    '''

    Parameters
    ----------
    arr : array-like, shape (n_frames, 3, 3) or (n_frames, 9)
    '''
    _nvecs = 'nvecs ' + str(nvecs)
    _rvecin = 'rvecin ' + rvecin if rvecin is not None else ''
    _rseed = 'rseed ' + str(rseed)
    _order = 'order ' + str(order)
    _ncorr = 'ncorr ' + str(ncorr) if ncorr > 0 else ''
    _tol = 'tol ' + str(tol)
    _d0 = 'd0 ' + str(d0)
    _nmesh = 'nmesh ' + str(nmesh)
    _dt = 'dt ' + str(dt)
    _ti = 'ti ' + str(ti)
    _tf = 'tf ' + str(tf)
    _itmax = 'itmax ' + str(itmax)
    _rmatrix = 'rmatrix mymat'

    act = CpptrajAnalyses.Analysis_Rotdif()
    dslist = CpptrajDatasetList()
    dslist.add_set('mat3x3', 'mymat')

    arr = np.asarray(arr, dtype='f8')
    msg = 'array must have shape=(n_frames, 9) or (n_frames, 3, 3)'
    shape = arr.shape
    if arr.ndim == 2:
        assert shape[1] == 9, msg
    elif arr.ndim == 3:
        assert shape[1:3] == (3, 3), msg
        # need to reshape to (n_frames, 9)
        arr = arr.reshape(shape[0], shape[1] * shape[2])
    else:
        raise ValueError(msg)


def lifetime(data, cut=0.5, rawcurve=False, more_options='', dtype='ndarray'):
    """lifetime (adapted lightly from cpptraj doc)

    Parameters
    ----------
    data : 1D-array or 2D array-like
    cut : cutoff to use when determining if data is 'present', default 0.5
    more_options : str, more cpptraj's options. Check cpptraj's manual.
    """
    data = np.asarray(data)
    if data.ndim == 1:
        data_ = [
            data,
        ]
    else:
        data_ = data

    _outname = 'name lifetime_'
    _cut = 'cut ' + str(cut)
    _rawcurve = 'rawcurve' if rawcurve else ''
    # do not sorting dataset's names. We can accessing by indexing them.
    _nosort = 'nosort'

    namelist = []
    cdslist = CpptrajDatasetList()
    for idx, arr in enumerate(data_):
        # create datasetname so we can reference them
        name = 'data_' + str(idx)
        if 'int' in arr.dtype.name:
            cdslist.add_set("integer", name)
        else:
            cdslist.add_set("double", name)
        cdslist[-1].data = np.asarray(arr)
        namelist.append(name)

    act = CpptrajAnalyses.Analysis_Lifetime()
    _cm = ' '.join(namelist)
    command = " ".join((_cm, _outname, _cut, _rawcurve, _nosort, more_options))
    act(command, dslist=cdslist)

    for name in namelist:
        cdslist.remove_set(cdslist[name])
    return get_data_from_dtype(cdslist, dtype=dtype)


# parallel


def to_parmed(traj, all_coords=False):
    import parmed as pmd

    parm = pmd.Structure()
    for atom in traj.top.simplify().atoms:
        p_atom = pmd.Atom(
            name=atom.name,
            type=atom.type,
            atomic_number=atom.atomic_number,
            charge=atom.charge,
            mass=atom.mass)
        parm.add_atom(p_atom, resname=atom.resname, resnum=atom.resid)
        # chain=str(atom.molnum))
    if all_coords:
        parm.coordinates = traj.xyz
    else:
        parm.coordinates = traj.xyz[0]
    parm.box = traj.top.box
    return parm
