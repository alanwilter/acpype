import time
from .check_and_assert import file_exist, is_generator
from .check_and_assert import eq
from .check_and_assert import _import, is_int
from .check_and_assert import has_, is_array
from .check_and_assert import ensure_not_none_or_string
from .Timer import Timer
from .context import tempfolder
from . import convert
from .tools import WrapBareIterator


def fn(name):
    # return absolute dir of pytraj/tests/data/name
    import pytraj
    base = pytraj.__path__[0]
    return base + '/../tests/data/' + name


class Timer:
    def __init__(self):
        self._t0 = self.value = None

    def __enter__(self):
        self._t0 = time.time()
        return self

    def __exit__(self, *args, **kwargs):
        self.value = time.time() - self._t0



def duplicate_traj(orig_traj, n_times):
    '''
    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.load_sample_data('ala3')
    >>> traj.n_frames
    1
    >>> duplicate_traj(traj, 3).n_frames
    3
    >>> t0 = pt.Trajectory(xyz=traj.xyz, top=traj.top)
    >>> duplicate_traj(t0, 3).n_frames
    4
    '''
    traj = orig_traj.copy()
    for _ in range(n_times - 1):
        if 'Iter' in orig_traj.__class__.__name__:
            # TrajectoryIterator
            traj._load(orig_traj.filelist)
        else:
            # Trajectory
            traj.append(traj.copy())
    return traj


def join_mask(m, res=None):
    """
    Examples
    --------
    >>> join_mask(('CA', 'CB'), res='1')
    ':1@CA :1@CB'
    >>> join_mask('CA CB', res='1')
    ':1@CA :1@CB'
    >>> join_mask('CA CB', res=0)
    ':1@CA :1@CB'
    """

    if is_int(res):
        res = str(res + 1)
    else:
        res = res

    if isinstance(m, str):
        # 'CA CB' to ['CA', 'CB']
        m = m.split()
    elif not isinstance(m, (list, tuple)):
        raise ValueError("must be a list/tuple")

    return " ".join(':' + res + '@' + s for s in m)


def split_range(n_chunks, start, stop):
    '''split a given range to n_chunks

    Examples
    --------
    >>> split_range(3, 0, 10)
    [(0, 3), (3, 6), (6, 10)]
    '''
    list_of_tuple = []
    chunksize = (stop - start) // n_chunks
    for i in range(n_chunks):
        if i < n_chunks - 1:
            _stop = start + (i + 1) * chunksize
        else:
            _stop = stop
        list_of_tuple.append((start + i * chunksize, _stop))
    return list_of_tuple
