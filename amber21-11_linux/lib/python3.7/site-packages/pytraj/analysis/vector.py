from __future__ import print_function, absolute_import
import numpy as np
from ..utils.decorators import register_pmap
from .c_action import do_action, c_action
from .c_action.actionlist import ActionList
from ..utils.get_common_objects import (
    get_topology,
    get_data_from_dtype,
    get_list_of_commands,
    get_fiterator,
    super_dispatch, )
from ..datasets.c_datasetlist import DatasetList as CpptrajDatasetList

SUPPORTED_TYPES = [
    x
    for x in
    'minimage dipole center corrplane box boxcenter ucellx ucelly ucellz principal'.
    split()
]

__all__ = ['vector', 'vector_mask'] + SUPPORTED_TYPES[:]


def _2darray_to_atommask_groups(seq):
    '''
    >>> list(_2darray_to_atommask_groups([[0, 3], [4, 7]]))
    ['@1 @4', '@5 @8']
    '''
    for arr in seq:
        # example: arr = [0, 3]; turns ot '@1 @4'
        yield '@' + str(arr[0] + 1) + ' @' + str(arr[1] + 1)


@register_pmap
def vector(traj=None,
           command="",
           frame_indices=None,
           dtype='ndarray',
           top=None):
    """perform vector calculation. See example below. Same as 'vector' command in cpptraj.

    Parameters
    ----------
    traj : Trajectory-like or iterable that produces :class:`pytraj.Frame`
    command : str or a list of strings, cpptraj command
    frame_indices : array-like, optional, default None
        only perform calculation for given frame indices
    dtype : output's dtype, default 'ndarray'
    top : Topology, optional, default None

    Returns
    -------
    out : numpy ndarray, shape (n_frames, 3) if command is a string
          numpy ndarray, shape (n_vectors, n_frames, 3) if command is a list of strings

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> data = pt.vector.vector(traj, "@CA @CB")
    >>> data = pt.vector.vector(traj, [("@CA @CB"),])
    >>> data = pt.vector.vector(traj, "principal z")
    >>> data = pt.vector.vector(traj, "principal x")
    >>> data = pt.vector.vector(traj, "ucellx")
    >>> data = pt.vector.vector(traj, "boxcenter")
    >>> data = pt.vector.vector(traj, "box")

    Notes
    -----
    It's faster to calculate with a list of commands.
    For example, if you need to perform 3 calculations for 'ucellx', 'boxcenter', 'box'
    like below:

    >>> data = pt.vector.vector(traj, "ucellx")
    >>> data = pt.vector.vector(traj, "boxcenter")
    >>> data = pt.vector.vector(traj, "box")

    You should use a list of commands for faster calculation.
    >>> comlist = ['ucellx', 'boxcenter', 'box']
    >>> data = pt.vector.vector(traj, comlist)
    """
    c_dslist = CpptrajDatasetList()
    top_ = get_topology(traj, top)
    list_of_commands = get_list_of_commands(command)
    fi = get_fiterator(traj, frame_indices)
    actlist = ActionList()

    for command in list_of_commands:
        act = c_action.Action_Vector()
        actlist.add(act, command, top_, dslist=c_dslist)
    actlist.compute(fi)

    return get_data_from_dtype(c_dslist, dtype=dtype)


@super_dispatch()
def multivector(traj,
                resrange,
                names,
                top=None,
                dtype='dataset',
                frame_indices=None):
    '''

    Parameters
    ----------
    traj : Trajectory-like
    resrange : str, residue range
    names : {str, tuple of str}
    top : Topology, optional
    dtype : str, default 'dataset'
    frame_indices : {None, 1D array-like}, optional, default None

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> vecs = pt.multivector(traj, resrange='1-5', names=('C', 'N'))
    >>> vecs = pt.multivector(traj, resrange='1-5', names='C N')
    '''
    _resrange = 'resrange ' + resrange
    if 'name1' in names or 'name2' in names:
        # cpptraj style
        _names = names
    else:
        if isinstance(names, str):
            name1, name2 = names.split()
        else:
            # try to unpack
            name1, name2 = names
        _names = ' '.join(('name1', name1, 'name2', name2))
    command = ' '.join((_resrange, _names))

    c_dslist, _ = do_action(traj, command, c_action.Action_MultiVector)
    return get_data_from_dtype(c_dslist, dtype)


@register_pmap
def vector_mask(traj=None,
                mask="",
                frame_indices=None,
                dtype='ndarray',
                top=None):
    """compute vector between two maskes

    Parameters
    ----------
    traj : Trajectory-like or iterable that produces Frame
    mask: str or array of string or array of intergers, shape (n_vectors, 2)
        vector maskes
    frame_indices : array-like or iterable that produces integer number
        frame indices
    dtype : str, default 'ndarray'
        output dtype
    top : Topology, optional, default None

    Returns
    -------
    if mask is a string, return 2D ndarray, shape (n_frames, 3)
    if mask is a list of strings or a 2D ndarray, return 3D ndarray, shape (n_vectors, n_frames, 3)

    Examples
    --------
    >>> # calcualte N-H vector
    >>> import pytraj as pt
    >>> import numpy as np
    >>> traj = pt.load_sample_data('tz2')
    >>> from pytraj import vector as va
    >>> n_indices = pt.select_atoms('@N', traj.top)
    >>> h_indices = n_indices + 1

    >>> # create n-h pair for vector calculation
    >>> n_h_pairs = np.array(list(zip(n_indices, h_indices)))
    >>> data_vec = va.vector_mask(traj, n_h_pairs, dtype='ndarray')

    >>> # compute vectors for specific frame indices (0, 4)
    >>> data_vec = va.vector_mask(traj, n_h_pairs, frame_indices=[0, 4], dtype='ndarray')
    """
    from pytraj.utils.get_common_objects import get_topology, get_data_from_dtype, get_fiterator
    from pytraj.utils.get_common_objects import get_list_of_commands
    from pytraj.datasets.c_datasetlist import DatasetList as CpptrajDatasetList
    from pytraj.analysis.c_action.c_action import Action_Vector
    from pytraj.analysis.c_action.actionlist import ActionList

    fi = get_fiterator(traj, frame_indices)
    _top = get_topology(fi, top)
    dslist = CpptrajDatasetList()
    template_command = ' mask '

    cm_arr = np.asarray(mask)
    if cm_arr.dtype.kind != 'i':
        list_of_commands = get_list_of_commands(mask)
    else:
        if cm_arr.ndim != 2:
            raise ValueError(
                'if mask is a numpy.ndarray, it must have ndim = 2')
        list_of_commands = _2darray_to_atommask_groups(cm_arr)

    actlist = ActionList()

    for command in list_of_commands:
        act = Action_Vector()
        _command = command + template_command
        actlist.add(act, _command, _top, dslist=dslist)
    actlist.compute(fi)
    return get_data_from_dtype(dslist, dtype=dtype)


_template = '''
@register_pmap
def %s(traj=None, command="", frame_indices=None, dtype='ndarray', top=None):
    """
    Parameters
    ----------
    traj : Trajectory-like
    command : str or a list-like of strings
    frame_indices : array-like, default None
        if specified, only perform calculation with given frames
    top : {str, Topology}, optional, default None
    """
    from pytraj.utils.get_common_objects import get_topology, get_data_from_dtype, get_fiterator
    from pytraj.utils.get_common_objects import get_list_of_commands
    from pytraj.datasets.c_datasetlist import DatasetList as CpptrajDatasetList
    from pytraj.analysis.c_action.c_action import Action_Vector
    from pytraj.analysis.c_action.actionlist import ActionList

    fi = get_fiterator(traj, frame_indices)
    _top = get_topology(fi, top)
    dslist = CpptrajDatasetList()
    template_command = ' %s '

    list_of_commands = get_list_of_commands(command)
    actlist = ActionList()

    for command in list_of_commands:
        act = Action_Vector()
        _command = command + template_command
        actlist.add(act, _command, _top, dslist=dslist)
    actlist.compute(fi)
    return get_data_from_dtype(dslist, dtype=dtype)
'''

for _key in SUPPORTED_TYPES:
    _my_func_str = _template % (_key, _key)
    exec(_my_func_str)

del _key
