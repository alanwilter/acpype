"""seperate module, only use stdlib
If want to use external package, import it inside the function

This module stores all useful functions that does not fit to anywhere else.
"""

import os
from itertools import islice
from collections import OrderedDict, defaultdict
import numpy as np
from functools import reduce


class WrapBareIterator(object):
    def __init__(self, obj, top):
        # obj : any iterable object
        # top : Topology
        self.obj = obj
        self.top = top

    def __iter__(self):
        return iter(self.obj)


def estimate_size(n_frames, n_atoms, dtype='f8'):
    '''return MB

    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> estimate_size(traj.n_frames, traj.n_atoms, 'f8')
    1.2114715576171875
    >>> estimate_size(traj.n_frames, traj.n_atoms, 'f4')
    0.6057357788085938
    '''
    if dtype == 'f8':
        n_bytes = 8
    if dtype == 'f4':
        n_bytes = 4
    return n_frames * n_atoms * 3 * n_bytes / (1024**2)


def groupby(key, seq):
    # lightly adapted from `toolz` package.
    # see license in $PYTRAJHOME/licenses/externals/toolz.txt
    '''
    Examples
    --------
    >>> names = ['Alice', 'Bob', 'Charlie', 'Dan', 'Edith', 'Frank']
    >>> _ = groupby(len, names)
    '''
    d = defaultdict(lambda: seq.__class__().append)
    for item in seq:
        d[key(item)](item)
    rv = {}
    for k, v in iteritems(d):
        rv[k] = v.__self__
    return rv


def _array_to_cpptraj_range(seq):
    # use "i+1" since cpptraj use 1-based index for mask
    '''
    Examples
    --------
    >>> _array_to_cpptraj_range([2, 4])
    '3,5'
    '''
    return ",".join((str(i + 1) for i in seq))


def iteritems(d, **kw):
    """Return an iterator over the (key, value) pairs of a dictionary.

    Examples
    --------
    >>> for k, v in iteritems({'x': 3, 'y': 4}): print(k, v) : # doctest: +SKIP
    x 3
    y 4
    """
    return iter(d.items())(**kw)


# this module gathers commonly used functions
# from toolz, stackoverflow, ... and from myself
# should make this independent from pytraj


def split(data, n_chunks):
    """split `self.data` to n_chunks

    Examples
    --------
    >>> for data in split(range(30), 3): print(data)
    [0 1 2 3 4 5 6 7 8 9]
    [10 11 12 13 14 15 16 17 18 19]
    [20 21 22 23 24 25 26 27 28 29]

    Notes
    -----
    same as numpy.array_split
    """
    return np.array_split(data, n_chunks)


def block_average(self, n_chunk):
    '''average by chunk

    Examples
    --------
    >>> block_average(range(30), 3)
    array([  4.5,  14.5,  24.5])
    '''

    return np.array(list(map(np.mean, split(self, n_chunk))))


def moving_average(data, n):
    """moving average

    Examples
    --------
    >>> moving_average([1, 2, 3, 4, 6], 2)
    array([ 0.5,  1.5,  2.5,  3.5,  5. ])

    Notes
    -----
    from `stackoverflow <http://stackoverflow.com/questions/11352047/finding-moving-average-from-data-points-in-python>`_
    """
    window = np.ones(int(n)) / float(n)
    return np.convolve(data, window, 'same')


def _compose2(f, g):
    # copied from pandas
    # see license in pytraj/license/
    """Compose 2 callables"""
    return lambda *args, **kwargs: f(g(*args, **kwargs))


def compose(*funcs):
    """
    Notes: copied from pandas (added pytraj's example)
    see license in pytraj/license/

    Compose 2 or more callables

    Examples
    --------
    >>> import pytraj as pt
    >>> from pytraj.testing import get_fn
    >>> func = compose(pt.calc_radgyr, pt.iterload)
    >>> fname, tname = get_fn('tz2')
    >>> func(fname, tname)
    array([ 18.91114428,  18.93654996,  18.84969884,  18.90449256,
            18.8568644 ,  18.88917208,  18.9430491 ,  18.88878079,
            18.91669565,  18.87069722])
    """
    assert len(funcs) > 1, 'At least 2 callables must be passed to compose'
    return reduce(_compose2, funcs)


def grep_key(self, key):
    """grep key

    Examples
    --------
    >>> import pytraj as pt
    >>> traj  = pt.load_sample_data('tz2')
    >>> dslist = pt.calc_multidihedral(traj, dtype='dataset')
    >>> pt.tools.grep_key(dslist, 'psi').values[0]
    array([ 176.6155643 ,  166.82129574,  168.79510009,  167.42561927,
            151.18334989,  134.17610997,  160.99207908,  165.1126967 ,
            147.94332109,  145.42901383])
    """
    new_self = self.__class__()
    for d in self:
        if key in d.key:
            new_self.append(d)
    return new_self


def flatten(x):
    """Returns a single, flat list which contains all elements retrieved
    from the sequence and all recursively contained sub-sequences
    (iterables).

    Notes
    -----
    from: http://kogs-www.informatik.uni-hamburg.de/~meine/python_tricks

    Examples
    --------
    >>> [1, 2, [3,4], (5,6)]
    [1, 2, [3, 4], (5, 6)]
    >>> flatten([[[1,2,3], (42,None)], [4,5], [6], 7, (8,9,10)])
    [1, 2, 3, 42, None, 4, 5, 6, 7, 8, 9, 10]"""

    result = []
    for el in x:
        # if isinstance(el, (list, tuple)):
        if hasattr(el, "__iter__") and not isinstance(el, str):
            result.extend(flatten(el))
        else:
            result.append(el)
    return result


def n_grams(a, n):
    """n_grams

    Parameters
    ----------
    a : sequence
    n : number of elements
    asarray : bool, default False
        if False: return an iterator
        if True: return a numpy array

    Examples
    --------
    >>> list(n_grams([2, 3, 4 ,5], 2))
    [(2, 3), (3, 4), (4, 5)]

    Notes
    -----
    adapted from: http://sahandsaba.com/thirty-python-language-features-and-tricks-you-may-not-know.html
    """

    z = (islice(a, i, None) for i in range(n))
    return zip(*z)


def dict_to_ndarray(dict_of_array):
    """convert OrderedDict to numpy array

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.load_sample_data('tz2')
    >>> dslist = pt.multidihedral(traj, dihedral_types='phi psi', resrange='2', dtype='dict')
    >>> list(dslist.keys())
    ['phi:2', 'psi:2']
    >>> dict_to_ndarray(dslist)
    array([[-128.72617304, -109.44321317, -130.93278259, ..., -146.70146067,
            -121.58263643, -112.74485175],
           [ 150.11249102,  142.52303293,  131.11609265, ...,  123.44883266,
             141.18992429,  120.03168126]])
    """
    if not isinstance(dict_of_array, OrderedDict):
        raise NotImplementedError("support only OrderedDict")
    return np.array([v for _, v in dict_of_array.items()])


def concat_dict(iterables):
    """concat dict

    Examples
    --------
    >>> dict_0 = {'x' : [1, 2, 3,]}
    >>> dict_1 = {'x' : [4, 5]}
    >>> concat_dict((dict_0, dict_1))
    OrderedDict([('x', array([1, 2, 3, 4, 5]))])
    """
    new_dict = OrderedDict()
    for i, d in enumerate(iterables):
        if i == 0:
            # make a copy of first dict
            new_dict.update(d)
        else:
            for k, v in new_dict.items():
                new_dict[k] = np.concatenate((new_dict[k], d[k]))
    return new_dict


def merge_coordinates(iterables):
    """merge_coordinates from frames

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.load_sample_data('tz2')
    >>> merge_coordinates(traj(0, 3)) # doctest: +SKIP
    array([[ 15.55458927,  28.54844856,  17.18908691],
           [ 16.20579147,  29.07935524,  17.74959946],
           [ 14.95065975,  29.27651787,  16.83513069],
           ...,
           [ 34.09399796,   7.88915873,  15.6500845 ],
           [ 34.4160347 ,   8.53098011,  15.01716137],
           [ 34.29132462,   8.27471733,  16.50368881]])
    """
    return np.vstack([f.xyz.copy() for f in iterables])


def merge_frames(iterables):
    """merge from frames to a single Frame. Order matters.
    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.load_sample_data('tz2')
    >>> traj[0]
    <Frame with 5293 atoms>
    >>> merge_frames(traj(0, 3))
    <Frame with 15879 atoms>
    """
    from pytraj import Frame
    xyz = np.vstack([f.xyz.copy() for f in iterables])
    frame = Frame()
    frame.append_xyz(xyz)
    return frame


def merge_frame_from_trajs(trajlist):
    """
    Examples
    --------
    >>> import numpy as np
    >>> import pytraj as pt
    >>> traj0 = pt.load_sample_data('tz2')[:3]
    >>> traj1 = pt.load_sample_data('tz2')[3:6]
    >>> traj2 = pt.load_sample_data('tz2')[6:9]
    >>> print(traj0.n_atoms, traj1.n_atoms, traj2.n_atoms)
    5293 5293 5293
    >>> for frame in pt.tools.merge_frame_from_trajs((traj0, traj1, traj2)): print(frame)
    <Frame with 15879 atoms>
    <Frame with 15879 atoms>
    <Frame with 15879 atoms>
   """
    for iterables in zip(*trajlist):
        yield merge_frames(iterables)


def rmsd_1darray(a1, a2):
    '''rmsd of a1 and a2

    Examples
    --------
    >>> a0 = [1, 3, 4]
    >>> a1 = [1.4, 3.5, 4.2]
    >>> rmsd_1darray(a0, a1)
    0.3872983346207417

    >>> rmsd_1darray(a0, [3, 4, 5, 7, 8])
    Traceback (most recent call last):
        ...
    ValueError: must have the same shape
    '''
    import numpy as np
    from math import sqrt
    arr1 = np.asarray(a1)
    arr2 = np.asarray(a2)

    if len(arr1.shape) > 1 or len(arr2.shape) > 1:
        raise ValueError("1D array only")

    if arr1.shape != arr2.shape:
        raise ValueError("must have the same shape")

    tmp = sum((arr1 - arr2)**2)
    return sqrt(tmp / arr1.shape[0])


def rmsd(a1, a2, flatten=True):
    """rmsd for two array with the same shape

    Parameters
    ----------
    a1, a2: np.ndarray
    flatten : bool, default True
        if True: always flatten two input arrays

    Examples
    --------
    >>> import pytraj as pt
    >>> t0 = pt.load_sample_data('ala3')
    >>> t1 = t0[:]
    >>> t1.xyz += 1.
    >>> rmsd(t0.xyz, t1.xyz)
    1.0

    Notes
    -----
    This method is different from ``pytraj.rmsd``
    """
    import numpy as np
    a1 = np.asarray(a1)
    a2 = np.asarray(a2)
    if a1.shape != a2.shape and not flatten:
        raise ValueError("must have the same shape")
    return rmsd_1darray(a1.flatten(), a2.flatten())


def mean_and_error(a1, a2):
    """calculate mean and error from two 1D array-like

    Examples
    --------
    >>> import pytraj as pt
    >>> a0 = [2, 4, 6]
    >>> a1 = [3, 5, 7]
    >>> mean_and_error(a0, a1)
    (4.5, 0.5)
    """
    import numpy as np
    mean = np.mean

    a1 = np.asarray(a1)
    a2 = np.asarray(a2)
    assert len(a1.shape) == len(a2.shape) == 1, "1D array"
    return (mean(a1 + a2) / 2, mean(np.abs(a1 - a2)) / 2)


def split_traj_by_residues(traj, start=0, stop=-1, step=1):
    '''return a generator

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_rna()
    >>> g = pt.tools.split_traj_by_residues(traj)
    >>> t0 = next(g)
    >>> print(t0.top.n_residues)
    1
    '''
    from pytraj.utils.cyutils import get_positive_idx

    _stop = get_positive_idx(stop, traj.top.n_residues)

    for i in range(start, _stop, step):
        j = ':' + str(i + 1)
        # example: traj[':3']
        yield traj[j]


def read_gaussian_output(filename=None, top=None):
    """return a `pytraj.trajectory.Trajectory` object

    Parameters
    ----------
    fname : str, filename
    top : {str, Topology}, optional, default None
        pytraj.Topology or a filename or None
        if None, use `antechamber` to generate mol2 file, need set $AMBERHOME env

    Requires
    --------
    cclib (``pip install cclib``)

    >>> import pytraj as pt
    >>> pt.tools.read_gaussian_output("gau.out", "mytest.pdb") # doctest: +SKIP
    """
    import cclib
    from pytraj.trajectory.trajectory import Trajectory
    from pytraj.utils.context import tempfolder
    from pytraj.utils.get_common_objects import get_topology

    _top = get_topology(None, top)
    gau = cclib.parser.Gaussian(filename)
    go = gau.parse()

    if _top is None:  # pragma: no cover
        try:
            amberhome = os.environ['AMBERHOME']
        except KeyError:
            raise KeyError("must set AMBERHOME")

        fpath = os.path.abspath(filename)

        with tempfolder():
            at = amberhome + "/bin/antechamber"
            out = "-i %s -fi gout -o tmp.mol2 -fo mol2 -at amber" % fpath
            cm = " ".join((at, out))
            os.system(cm)

            return Trajectory(xyz=go.atomcoords, top="tmp.mol2")
    else:
        return Trajectory(xyz=go.atomcoords, top=_top)


def merge_trajs(traj1, traj2, start_new_mol=True, n_frames=None):
    """

    Examples
    --------
    >>> import pytraj as pt
    >>> import numpy as np
    >>> traj1 = pt.load_sample_data('ala3')[:1]
    >>> traj2 = pt.load_sample_data('tz2')[:1]
    >>> traj3 = merge_trajs(traj1, traj2)
    >>> # from frame_iter for saving memory
    >>> traj3 = merge_trajs((traj1(0, 10, 2), traj1.top), (traj2(100, 110, 2), traj2.top), n_frames=6)

    >>> # raise error if not having the same n_frames
    >>> traj4 = pt.load_sample_data('tz2')[:]
    >>> traj4.n_frames
    10
    >>> traj1.n_frames
    1
    >>> merge_trajs(traj1, traj4)
    Traceback (most recent call last):
        ...
    ValueError: must have the same n_frames


    Notes
    -----
    Code might be changed
    """
    from pytraj import Trajectory
    import numpy as np

    if isinstance(traj1, (list, tuple)):
        n_frames_1 = n_frames
        top1 = traj1[1]
        _traj1 = traj1[0]
    else:
        n_frames_1 = traj1.n_frames
        top1 = traj1.top
        _traj1 = traj1

    if isinstance(traj2, (list, tuple)):
        n_frames_2 = n_frames
        top2 = traj2[1]  # example: (traj(0, 5), traj.top)
        _traj2 = traj2[0]
    else:
        n_frames_2 = traj2.n_frames
        top2 = traj2.top
        _traj2 = traj2

    if n_frames_1 != n_frames_2:
        raise ValueError("must have the same n_frames")

    traj = Trajectory()
    traj._allocate(n_frames_1, top1.n_atoms + top2.n_atoms)

    # merge Topology
    top = top1.copy()
    #if start_new_mol:
    #    top.start_new_mol()
    top.join(top2)
    traj.top = top

    # update coords
    for f1, f2, frame in zip(_traj1, _traj2, traj):
        frame.xyz = np.vstack((f1.xyz, f2.xyz))

    return traj


def as_2darray(traj_or_xyz):
    '''reshape traj.xyz to 2d array, shape=(n_frames, n_atoms * 3)

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.load_sample_data('tz2')
    >>> traj.xyz.shape
    (10, 5293, 3)
    >>> as_2darray(traj).shape
    (10, 15879)
    >>> as_2darray(traj.xyz).shape
    (10, 15879)

    Notes
    -----
    if ``traj`` is mutable, this method return a view of its coordinates.
    '''
    import numpy as np

    if hasattr(traj_or_xyz, 'xyz'):
        traj = traj_or_xyz
        # Trajectory-like
        return traj.xyz.reshape(traj.n_frames, traj.n_atoms * 3)
    else:
        # array-like, assume 3D
        xyz = np.asarray(traj_or_xyz)
        assert xyz.ndim == 3, 'xyz must has ndim=3'
        shape = xyz.shape
        return xyz.reshape(shape[0], shape[1] * shape[2])


def as_3darray(xyz):
    '''reshape xyz to 3d array, shape=(n_frames, n_atoms, 3)

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.load_sample_data('tz2')
    >>> traj.xyz.shape
    (10, 5293, 3)
    >>> xyz_2d = as_2darray(traj)
    >>> xyz_2d.shape
    (10, 15879)
    >>> as_3darray(xyz_2d).shape
    (10, 5293, 3)
    >>> as_3darray(traj.xyz)
    Traceback (most recent call last):
        ...
    ValueError: ndim must be 2
    '''
    shape = xyz.shape
    if len(shape) != 2:
        raise ValueError('ndim must be 2')
    new_shape = (shape[0], int(shape[1] / 3), 3)
    return xyz.reshape(new_shape)


def split_and_write_traj(self,
                         n_chunks=None,
                         root_name="trajx",
                         ext='nc',
                         *args,
                         **kwd):
    '''
    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.load_sample_data('tz2')
    >>> split_and_write_traj(traj, n_chunks=3, root_name='output/trajx', overwrite=True)
    '''

    chunksize = self.n_frames // n_chunks
    for idx, traj in enumerate(self.iterchunk(chunksize=chunksize)):
        fname = ".".join((root_name, str(idx), ext))
        traj.save(fname, *args, **kwd)


def read_to_array(fname):
    '''read text from file to numpy array'''
    import numpy as np
    with open(fname, 'r') as fh:
        arr0 = np.array([[x for x in line.split()] for line in fh.readlines()])
        return np.array(flatten(arr0), dtype='f8')


def make_fake_topology(n_atoms):
    '''make fake Topology, just for writing xyz array to supported formats
    (netcdf, dcd, trr, ...)

    >>> import pytraj as pt
    >>> top = pt.tools.make_fake_topology(100)
    >>> top.n_atoms
    100
    >>> isinstance(top, pt.Topology)
    True
    >>> import numpy as np
    >>> xyz = np.random.rand(10*100*3).reshape(10, 100, 3)
    >>> traj0 = pt.Trajectory(xyz=xyz, top=top)
    >>> pt.write_traj('output/test.nc', traj0, overwrite=True)
    >>> traj = pt.iterload('output/test.nc', top=top)
    >>> traj.n_atoms
    100
    '''
    from pytraj import Atom, Residue, Topology

    top = Topology()

    for _ in range(n_atoms):
        atom = Atom(name='X', type='X', charge=0., mass=0., resid=0)
        residue = Residue('Y', 0)
        top.add_atom(atom, residue)
    return top


def dir_(obj):
    '''return a list of obj's attributes with no private method
    '''
    return [x for x in dir(obj) if not x.startswith('_')]
