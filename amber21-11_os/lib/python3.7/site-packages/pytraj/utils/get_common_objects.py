from __future__ import absolute_import
# TODO: rename this file
from functools import wraps

# do not import anything else here.
from pytraj.utils.convert import array_to_cpptraj_atommask
from pytraj.trajectory.shared_methods import iterframe_master


def _load_Topology(filename):
    from pytraj import Topology, ParmFile
    top = Topology()
    parm = ParmFile()
    parm.read(filename, top)
    return top


def get_topology(traj, top):
    '''

    >>> import pytraj as pt
    >>> from pytraj import Topology
    >>> traj = pt.datafiles.load_rna()
    >>> top_filename = traj.top.filename
    >>> isinstance(get_topology(traj, None), Topology)
    True
    >>> isinstance(get_topology(traj, top_filename), Topology)
    True
    >>> top = traj.top
    >>> isinstance(get_topology(traj, top), Topology)
    True
    >>> # find Topology in a list
    >>> isinstance(get_topology([traj, traj], None), Topology)
    True
    >>> get_topology(None, None) is None
    True
    '''
    if isinstance(top, str):
        # if provide a filename, load to Topology
        top_ = _load_Topology(top)
    elif top is None:
        # if user does not provide Topology, try to find it in traj
        if hasattr(traj, 'top'):
            top_ = traj.top
        else:
            # list, tuple of traj objects
            try:
                for tmp in traj:
                    if hasattr(tmp, 'top'):
                        top_ = tmp.top
                        break
            except TypeError:
                top_ = None
    else:
        top_ = top
    return top_


def get_data_from_dtype(d0, dtype='dataset'):
    from pytraj.datasets.datasetlist import DatasetList as DSL

    if (dtype is None or dtype == 'dataset') and hasattr(d0, 'set_own_memory'):
        d0.set_own_memory(False)

    dtype = dtype.lower()
    if dtype == 'dataset':
        return DSL(d0)
    elif dtype == 'ndarray':
        return d0.to_ndarray()
    elif dtype == 'dict':
        return d0.to_dict()
    elif dtype == 'dataframe':
        if hasattr(d0, 'key'):
            d0.key = d0.key.replace(':', '_')
            d0.key = d0.key.replace('-', '_')
        else:
            for _d in d0:
                _d.key = _d.key.replace(':', '_')
                _d.key = _d.key.replace('-', '_')
        return d0.to_dataframe()
    elif dtype == 'cpptraj_dataset':
        return d0
    else:
        raise NotImplementedError()


def get_list_of_commands(mask_or_commands):
    '''

    Examples
    --------
    >>> get_list_of_commands('@CA')
    ['@CA']
    >>> get_list_of_commands(('@CA'))
    ['@CA']
    >>> get_list_of_commands(100)
    Traceback (most recent call last):
        ...
    ValueError: must be string or list/tuple of strings
    '''
    if isinstance(mask_or_commands, str):
        return [
            mask_or_commands,
        ]
    elif isinstance(mask_or_commands, (list, tuple)):
        return list(mask_or_commands)
    else:
        raise ValueError("must be string or list/tuple of strings")


def get_matrix_from_dataset(dset, mat_type='full'):
    '''return full or half matrix

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> dset = pt.matrix.dist(traj, '@CA', dtype='cpptraj_dataset')[0]
    >>> get_matrix_from_dataset(dset, mat_type='full').shape
    (12, 12)
    >>> get_matrix_from_dataset(dset, mat_type='half').shape
    (78,)
    '''
    # dset in DatasetMatrixDouble object
    if mat_type == 'full':
        return dset.values
    elif mat_type == 'half':
        return dset._to_cpptraj_sparse_matrix()
    else:
        raise ValueError()


def get_reference(traj, ref):
    '''try best to get reference

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> frame = get_reference(traj, 3)
    >>> isinstance(frame, pt.Frame)
    True
    >>> frame = get_reference(traj, None)
    >>> isinstance(frame, pt.Frame)
    True
    >>> ref = traj[5]
    >>> frame = get_reference(traj, ref)
    '''
    if isinstance(ref, int):
        try:
            ref_ = traj[ref]
            ref_.top = traj.top
        except TypeError:
            raise TypeError("%s does not support indexing" % str(traj))
    elif ref is None:
        try:
            ref_ = traj[0]
            ref_.top = traj.top
        except TypeError:
            raise TypeError(
                "If reference is an integer, %s must support indexing" %
                traj.__str__())
    elif 'Trajectory' in ref.__class__.__name__:
        assert ref.n_frames == 1, "only support 1-frame Trajectory as reference"
        ref_ = ref[0]
        ref_.top = ref.top
    else:
        ref_ = ref

    return ref_


def get_fiterator(traj, frame_indices=None):
    if frame_indices is None:
        return traj
    else:
        return traj.iterframe(frame_indices=frame_indices)


def get_resrange(resrange):
    '''return resrange as a string

    Examples
    --------
    >>> get_resrange('1-3')
    'resrange 1-3'
    >>> get_resrange(0)
    'resrange 1'
    >>> get_resrange(range(3))
    'resrange 1,2,3'
    >>> get_resrange([2, 5, 7])
    'resrange 3,6,8'
    >>> get_resrange(None)
    ''
    '''
    from pytraj.utils import convert, is_int

    if resrange is not None:
        if is_int(resrange):
            resrange = [
                resrange,
            ]
        if isinstance(resrange, str):
            _resrange = "resrange " + resrange
        else:
            _resrange = convert.array_to_cpptraj_range(resrange)
            _resrange = "resrange " + str(_resrange)
    else:
        _resrange = ""
    return _resrange


class super_dispatch(object):
    # TODO: more descriptive method name?
    '''apply a series of functions to ``f``'s args and kwargs

    - get Topology from a given traj (Trajectory, frame iterator, ...) and top
        get_topology(traj, top)
    - create frame iterator from traj and frame_indices
        get_fiterator(traj, frame_indices)
    - create Amber mask from atom index array
        array_to_cpptraj_atommask(mask)
    - convert int ref to Frame ref
    '''

    def __init__(self, refindex=None):
        self.refindex = refindex

    def __call__(self, f):
        import inspect

        try:
            args_spec = inspect.getfullargspec(f)
        except AttributeError:
            # py2
            args_spec = inspect.getargspec(f)
        n_default = len(args_spec.defaults) if args_spec.defaults else 0
        try:
            kwargs_spec = dict(
                (k, v)
                for (
                    k,
                    v) in zip(args_spec.args[-n_default:], args_spec.defaults))
        except TypeError:
            kwargs_spec = {}

        has_ref_arg = 'ref' in args_spec.args
        has_mask_arg = 'mask' in args_spec.args

        @wraps(f)
        def inner(*args, **kwargs):
            args = list(args)
            # traj is always 1st argument
            if 'traj' in kwargs:
                traj = kwargs.get('traj')
                has_traj_arg = True
            else:
                traj = args[0]
                has_traj_arg = False

            mask = kwargs.get('mask', kwargs_spec.get('mask'))
            ref = kwargs.get('ref', kwargs_spec.get('ref'))
            frame_indices = kwargs.get('frame_indices')
            top = kwargs.get('top')

            if has_mask_arg and isinstance(mask, str):
                if mask == '':
                    if has_traj_arg:
                        try:
                            mask = args[0]
                        except IndexError:
                            mask = ''
                    else:
                        try:
                            mask = args[1]
                        except IndexError:
                            mask = ''

            if has_ref_arg:
                if ref is None:
                    try:
                        ref = args[
                            self.
                            refindex] if self.refindex is not None else args[2]
                    except IndexError:
                        ref = 0
                ref = get_reference(traj, ref)

            # update traj to args or kwargs
            if 'traj' in kwargs:
                kwargs['traj'] = get_fiterator(traj, frame_indices)
            else:
                args[0] = get_fiterator(traj, frame_indices)

            # update topology to kwargs
            kwargs['top'] = get_topology(traj, top)

            # update reference to args or kwargs
            if has_ref_arg:
                kwargs['ref'] = get_reference(traj, ref)

            # update mask to args or kwargs
            if has_mask_arg and not isinstance(mask, str):
                mask = array_to_cpptraj_atommask(mask)
            if 'mask' in kwargs:
                kwargs['mask'] = mask
            else:
                if has_mask_arg:
                    if 'traj' not in kwargs:
                        try:
                            args[1] = mask
                        except IndexError:
                            args.append(mask)
            return f(*args, **kwargs)

        inner._is_super_dispatched = True
        return inner


def get_iterator_from_dslist(traj,
                             mask,
                             frame_indices,
                             top,
                             crdname='dataset_coords'):
    from pytraj import Trajectory, TrajectoryIterator
    from pytraj.datasets import CpptrajDatasetList

    dslist = CpptrajDatasetList()
    dslist.add("coords", crdname)
    # need to set "rmsout" to trick cpptraj not giving error
    # need " " (space) before crdset too

    if isinstance(traj, (Trajectory, TrajectoryIterator)):
        # we do atom stripping here before copying to DatasetCoordsCRD to save memory if
        # loading from TrajectoryIterator
        fi = traj.iterframe(mask=mask, frame_indices=frame_indices)
        command = ''
        # use Topology from fi (could be stripped to save memory)
        dslist[0].top = fi.top
        top_ = fi.top
    else:
        # ignore frame_indices
        fi = iterframe_master(traj)
        command = mask
        top_ = get_topology(traj, top)
        dslist[0].top = top_
    for frame in fi:
        dslist[0].append(frame)
    return dslist, top_, command
