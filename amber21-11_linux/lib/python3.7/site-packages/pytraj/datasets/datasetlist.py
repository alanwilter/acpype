from __future__ import absolute_import
import numpy as np
from collections import OrderedDict
from pytraj.datasets.c_datasetlist import DatasetList as CpptrajDatasetList
from pytraj.utils import is_int, is_array, is_generator
from pytraj.datafiles import DataFile
from pytraj.core.c_core import ArgList
from .array import DataArray

__all__ = ['load_datafile', 'stack', 'DatasetList']


def _groupby(self, key):
    # adapted from `toolz` package.
    # see license in $PYTRAJHOME/licenses/externals/toolz.txt
    import collections
    d = collections.defaultdict(lambda: self.__class__().append)
    for item in self:
        d[key(item)](item)
    rv = {}
    for k, v in d.items():
        rv[k] = v.__self__
    return rv


def load_datafile(filename):
    """load cpptraj's output
    >>> d = load_datafile('data/tc5b.native_contacts.dat')
    """
    ds = DatasetList()
    ds.read_data(filename)
    return ds


def stack(args):
    """return a new DatasetList by joining (vstack)

    Parameters
    ----------
    args : list/tuple of DatasetList

    Notes
    -----
        similiar to numpy.vstack

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.load_sample_data('tz2')
    >>> traj1 = traj[:3]
    >>> traj2 = traj[5:]

    >>> d1 = pt.dssp(traj1, dtype='dataset')
    >>> len(d1[0])
    3
    >>> d2 = pt.dssp(traj2, dtype='dataset')
    >>> len(d2[0])
    5

    >>> d3 = stack((d1, d2))
    >>> len(d3[0])
    8

    >>> d4 = d1.copy()
    >>> d4[0].key = 'dasffafa'
    >>> d5 = stack((d1, d4))
    Traceback (most recent call last):
        ...
    KeyError: "Don't support stack different key"
    """
    is_subcriptable = not (isinstance(args, map) or is_generator(args))

    if is_subcriptable:
        dslist0 = args[0].copy()
    else:
        dslist0 = next(args)

    dslist_iter = args[1:] if is_subcriptable else args

    for dslist in dslist_iter:
        for d0, d in zip(dslist0, dslist):
            if d0.dtype != d.dtype:
                raise TypeError("Dont support stack different dtype together")
            if d0.key != d.key:
                raise KeyError("Don't support stack different key")
            d0.append(d.copy())
    return dslist0


concat_datasetlist = stack


class DatasetList(list):
    '''similiar to python's list but the data is labeled.
    Think as a OrderedDict-like and list-like object. This class is suitable for small
    datasets. For high performance, user should use pandas' DataFrame.

    Examples
    --------
    >>> import pytraj as pt
    >>> from pytraj import DatasetList
    >>> traj = pt.load_sample_data('tz2')[:2]
    >>> dslist = pt.multidihedral(traj, dtype='dataset')
    >>> dslist['phi:2'].values
    array([-128.72617304, -109.44321317])

    >>> # make a copy
    >>> d2 = dslist.copy()

    >>> d3 = DatasetList({'x': [0, 3, 7], 'y': [5, 6]})
    '''

    def __init__(self, dslist=None, copy=False):
        if dslist:
            if isinstance(dslist, dict):
                # {'x': [1, 3, 5], 'y': [4, 7, 8]}
                for key, values in dslist.items():
                    self.append(DataArray({key: values}), copy=copy)
            else:
                for d0 in dslist:
                    # always make a copy
                    # from DataArray(d0)
                    # so we set copy=False here
                    # to avoid copying twice
                    self.append(DataArray(d0), copy=copy)

    def copy(self):
        dslist = self.__class__()
        for d0 in self:
            dslist.append(d0, copy=True)
        return dslist

    def dtypes(self):
        '''
        >>> import pytraj as pt
        >>> dslist = pt.multidihedral(pt.load_sample_data('ala3'), dtype='dataset')
        >>> dslist.dtypes()
        [dtype('float64'), dtype('float64'), dtype('float64'), dtype('float64'), dtype('float64'), dtype('float64')]
        '''
        return [d.dtype for d in self]

    def __str__(self):
        safe_msg = "<pytraj.DatasetList with %s datasets>\n" % self.size
        if self.size == 0:
            return safe_msg
        msg = "\n\n".join("\n".join((d.key, d.values.__str__())) for d in self)
        str_first_3 = "\n\n".join("\n".join((d.key, d.values.__str__()))
                                  for d in self[:3])
        str_last_2 = "\n\n".join("\n".join((d.key, d.values.__str__()))
                                 for d in self[-2:])

        if self.size <= 5:
            return safe_msg + msg
        else:
            return safe_msg + str_first_3 + "\n...\n\n" + str_last_2

    def __repr__(self):
        return self.__str__()

    @property
    def size(self):
        return len(self)

    def __getslice__(self, i, j):
        return self.__getitem__(slice(i, j))

    def __getitem__(self, idx):
        """return a DataSet instance. Memory view is applied (which mean this new insance
        is just alias of self[idx]). Should we use a copy instead?

        Examples
        --------
        >>> import pytraj as pt
        >>> traj = pt.datafiles.load_tz2_ortho()
        >>> dslist = pt.multidihedral(traj)
        >>> d0 = dslist[0]
        >>> d1 = dslist['phi:3']
        >>> d2 = dslist[:6:2]
        >>> d3 = dslist[[0, 3, 8]]
        >>> d4 = dslist.__getslice__(0, 3)

        >>> d5 = d3.__class__()
        >>> d5[0]
        Traceback (most recent call last):
            ...
        ValueError: size = 0: can not index

        >> # dummy
        >>> d6 = dslist[pt.Frame]
        Traceback (most recent call last):
            ...
        ValueError: index must be int, string, slice or array-like
        """
        if self.size == 0:
            raise ValueError("size = 0: can not index")

        if is_int(idx):
            return super(DatasetList, self).__getitem__(idx)
        elif isinstance(idx, str):
            for d0 in self:
                if d0.key.upper() == idx.upper():
                    d0._base = self
                    return d0
        elif isinstance(idx, slice):
            # return new view of `self`
            start, stop, step = idx.indices(self.size)
            new_dslist = self.__class__()
            for _idx in range(start, stop, step):
                new_dslist.append(self[_idx], copy=False)
            return new_dslist
        elif is_array(idx) or isinstance(
                idx, list) and not isinstance(idx[0], bool):
            new_dslist = self.__class__()
            for _idx in idx:
                new_dslist.append(self[_idx], copy=False)
            return new_dslist
        else:
            raise ValueError('index must be int, string, slice or array-like')

    def keys(self):
        """return a list

        >>> import pytraj as pt
        >>> dslist = pt.multidihedral(pt.load_sample_data('tz2'), dtype='dataset')
        >>> keys = dslist.keys()
        """
        tmp_list = []
        for d0 in self:
            tmp_list.append(d0.key)
        return tmp_list

    def grep(self, key, mode='key', copy=False):
        """"return a new DatasetList object as a view of `self`

        Parameters
        ----------
        key : str or list
            keyword for searching
        mode: str, default='key'
            mode = 'key' | 'name' | 'dtype' | 'aspect'

        Examples
        --------
        >>> import pytraj as pt
        >>> traj = pt.load_sample_data('tz2')
        >>> dslist = pt.multidihedral(traj, dtype='dataset')
        >>> sub_dslist = dslist.grep('phi')
        >>> sub_dslist = dslist.grep(['phi', 'psi'])

        >>> dslist.grep(3)
        Traceback (most recent call last):
            ...
        ValueError: support string or list/tuple of strings
        """
        import re

        # use __class__ so we can `filter` return the same class
        # we subclass this Cython class to python level
        dtmp = self.__class__()

        # dont free mem here
        for d0 in self:
            att = getattr(d0, mode)
            if isinstance(key, str):
                if re.search(key, att):
                    dtmp.append(d0, copy=copy)
            elif isinstance(key, (list, tuple)):
                for _key in key:
                    if re.search(_key, att):
                        dtmp.append(d0, copy=copy)
            else:
                raise ValueError("support string or list/tuple of strings")
        return dtmp

    def to_dict(self):
        """return a dict object with key=key, value=list"""
        return OrderedDict((d0.key, d0.to_ndarray()) for d0 in self)

    @property
    def values(self):
        """return read-only ndarray"""
        return self.to_ndarray()

    def to_ndarray(self):
        """
        >>> DatasetList({'x': [0, 2]}).values
        array([0, 2])
        """
        if self.size == 1:
            return self[0].values
        else:
            # more than one set
            return np.array([x.values for x in self])

    def to_dataframe(self):
        """return pandas' DataFrame

        Requires
        --------
        pandas
        """
        import pandas
        my_dict = OrderedDict((d0.key, d0.to_ndarray(copy=True))
                              for d0 in self)
        return pandas.DataFrame(my_dict)

    @classmethod
    def read_data(cls, filename, arg=""):
        '''
        >>> from pytraj.datasets.datasetlist import DatasetList
        >>> d = DatasetList.read_data('data/tc5b.native_contacts.dat')
        '''
        df = DataFile()
        dslist = CpptrajDatasetList()
        df.read_data(filename, ArgList(arg), dslist)
        return DatasetList(dslist)

    def append(self, dset, copy=True):
        """append anythin having `key` attribute
        >>> d = DatasetList()
        >>> d.append({'x' : [0, 3]})

        >>> import pytraj as pt
        >>> dslist = pt.multidihedral(pt.load_sample_data('tz2'), dtype='dataset')
        >>> dslist[0].key
        'psi:1'
        >>> d.append(dslist[0])
        >>> d.append(dslist[1], copy=False)
        >>> d.append(dslist[1])
        Traceback (most recent call last):
            ...
        KeyError: 'must have different key'

        >>> d.append(100)
        Traceback (most recent call last):
            ...
        ValueError: must have key or be a dict
        """
        if copy:
            try:
                d0 = dset.copy()
            except AttributeError:
                d0 = dset
        else:
            d0 = dset

        def check_key(self, key):
            for k0 in self.keys():
                if k0 == key:
                    raise KeyError("must have different key")

        if hasattr(dset, 'key'):
            check_key(self, dset.key)
            super(DatasetList, self).append(d0)
        elif isinstance(dset, dict):
            for key, values in dset.items():
                check_key(self, key)
                super(DatasetList, self).append(DataArray({key: values}))
        else:
            raise ValueError('must have key or be a dict')

    def groupby(self, func_or_key):
        '''
        Examples
        --------
        >>> import pytraj as pt
        >>> traj = pt.load_sample_data('tz2')
        >>> dslist = pt.multidihedral(traj, dtype='dataset')
        >>> x = dslist.groupby(lambda x: 'psi' in x.key)
        '''
        return _groupby(self, func_or_key)
