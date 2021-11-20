"""perform nucleic acid analysis
"""

import numpy as np
from ..utils.get_common_objects import get_topology, get_resrange
from ..utils.get_common_objects import get_reference, get_fiterator

__all__ = ['nastruct', 'nupars']


def _group(self, key):
    # adapted from `toolz` package.
    # see license in $PYTRAJHOME/licenses/externals/toolz.txt
    import collections
    # d = collections.defaultdict(lambda: self.__class__().append)
    d = collections.defaultdict(lambda: [].append)
    for item in self:
        d[key(item)](item)
    rv = {}
    for k, v in d.items():
        rv[k] = v.__self__
    return rv


def nastruct(traj=None,
             ref=0,
             resrange=None,
             resmap=None,
             hbcut=3.5,
             frame_indices=None,
             pucker_method='altona',
             dtype='nupars',
             baseref=None,
             groove_3dna=True,
             top=None):
    """compute nucleic acid parameters. (adapted from cpptraj doc)

    Parameters
    ----------
    traj : Trajectory-like
    ref : {Frame, int}, default 0 (first frame)
    resrange : None, str or array-like of integers
    resmap : residue map, example: 'AF2:A'
    hbcut : float, default=3.5 Angstrong
        Distance cutoff for determining basepair hbond
    pucker_method : str, {'altona', 'cremer'}, default 'altona'
        'altona' : Use method of Altona & Sundaralingam to calculate sugar pucker
        'cremer' : Use method of Cremer and Pople to calculate sugar pucker'
    frame_indices : array-like, default None (all frames)
    baseref : {None, str}, default None
        if given filename, use it for base reference (versionadded: 1.0.6)
    groove_3dna : bool, default True
        if True, major and minor groove will match 3DNA's output.
    dtype : str, {'nupars', 'cpptraj_dataset'}, default 'nupars'

    Returns
    -------
    out : nupars object. One can assess different values (major groove width, xdips values
    ...) by accessing its attribute. See example below.

    Examples
    --------
    >>> import pytraj as pt
    >>> import numpy as np
    >>> traj = pt.datafiles.load_rna()
    >>> data = pt.nastruct(traj, groove_3dna=False)
    >>> data.keys()[:5] # doctest: +SKIP
    ['buckle', 'minor', 'major', 'xdisp', 'stagger']
    >>> # get minor groove width values for each pairs for each snapshot
    >>> # data.minor is a tuple, first value is a list of basepairs, seconda value is
    >>> # numpy array, shape=(n_frames, n_pairs)

    >>> data.minor # doctest: +SKIP
    (['1G16C', '2G15C', '3G14C', '4C13G', '5G12C', '6C11G', '7C10G', '8C9G'],
     array([[ 13.32927036,  13.403409  ,  13.57159901, ...,  13.26655865,
             13.43054485,  13.4557209 ],
           [ 13.32002068,  13.45918751,  13.63253593, ...,  13.27066231,
             13.42743683,  13.53450871],
           [ 13.34087658,  13.53778553,  13.57062435, ...,  13.29017353,
             13.38542843,  13.46101475]]))

    >>> data.twist # doctest: +SKIP
    (['1G16C-2G15C', '2G15C-3G14C', '3G14C-4C13G', '4C13G-5G12C', '5G12C-6C11G', '6C11G-7C10G', '7C10G-8C9G'],
    array([[ 34.77773666,  33.98158646,  30.18647003, ...,  35.14608765,
             33.9628334 ,  33.13056946],
           [ 33.39176178,  32.68476105,  28.36385536, ...,  36.59774399,
             30.20827484,  26.48732948],
           [ 36.20665359,  32.58955002,  27.47707367, ...,  33.42843246,
             30.90047073,  33.73724365]]))

    >>> # change dtype
    >>> data = pt.nastruct(traj, dtype='cpptraj_dataset')
    """
    from pytraj.datasets.c_datasetlist import DatasetList as CpptrajDatasetList
    from .c_action import c_action
    from pytraj.datasets.array import DataArray

    _resrange = get_resrange(resrange)

    fi = get_fiterator(traj, frame_indices)
    ref = get_reference(traj, ref)
    _top = get_topology(traj, top)
    _resmap = "resmap " + resmap if resmap is not None else ""
    _hbcut = "hbcut " + str(hbcut) if hbcut is not None else ""
    _pucker_method = pucker_method
    _groove_3dna = 'groovecalc 3dna' if groove_3dna else ''
    baseref_ = 'baseref {}'.format(baseref) if baseref is not None else ''

    command = " ".join((_resrange, _resmap, _hbcut, _pucker_method,
                        _groove_3dna, baseref_))

    dslist = CpptrajDatasetList()

    # need to construct 3 steps so we can pickle this method for parallel
    # not sure why?
    act = c_action.Action_NAstruct(command, top=_top, dslist=dslist)
    act.compute(ref)
    act.compute(fi)
    act.post_process()

    if dtype == 'cpptraj_dataset':
        return dslist
    elif dtype == 'nupars':
        dslist_py = []
        for d in dslist:
            dslist_py.append(DataArray(d))
            dslist_py[-1].values = dslist_py[-1].values[1:]
        return nupars(_group(dslist_py, lambda x: x.aspect))
    else:
        raise ValueError("only support dtype = {'nupars', 'cpptraj_dataset'}")


class nupars(object):
    '''class holding data for nucleic acid.

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_rna()
    >>> nu = pt.nastruct(traj, groove_3dna=False)
    >>> nu.major # doctest: +SKIP
    (['1G16C', '2G15C', '3G14C', '4C13G', '5G12C', '6C11G', '7C10G', '8C9G'],
     array([[  0.        ,  18.60012245,  18.7782402 , ...,  18.45940208,
              18.78943062,   0.        ],
            [  0.        ,  18.82868767,  19.09348869, ...,  18.86692429,
              18.89782524,   0.        ],
            [  0.        ,  18.16820908,  18.77448463, ...,  19.1559639 ,
              18.38820267,   0.        ]]))
    >>> nu.keys()[:10] # doctest: +SKIP
    ['shift',
     'prop',
     'buckle',
     'rise',
     'hrise',
     'major',
     'twist',
     'minor',
     'zp',
     'hb']
    >>> import numpy as np
    >>> nu._summary(np.mean, keys=['major'])
    {'major': array([ 14.36732316,  14.3188405 ,  14.28000116])}
    '''

    def __init__(self, adict):
        self._dict = adict

    def __str__(self):
        return '<nupars, keys = %s>' % str(self.keys())

    def __repr__(self):
        return str(self)

    def __getitem__(self, key):
        '''self['minor'], ...
        '''
        return self.__getattr__(key)

    def __getattr__(self, aspect):
        '''self.minor, ...
        '''
        # data is a list of DataArray
        data = self._dict[aspect]
        arr = np.empty((len(data), len(data[0])), dtype='f8')

        keylist = []
        for idx, arr0 in enumerate(data):
            keylist.append(arr0.key)
            arr[idx] = arr0.values
        return keylist, arr.T

    def keys(self):
        return list(self._dict)

    def __dir__(self):
        '''for autocompletion in ipython
        '''
        return self.keys() + ['_summary', '_explain']

    def _summary(self, ops, keys=None, indices=None):
        '''
        Parameters
        op : numpy method
        keys: optional
        indices : optional

        Examples
        --------
        self._summary(np.mean, indices=range(2, 8))
        '''
        _keys = keys if keys is not None else self.keys()
        if isinstance(_keys, str):
            _keys = [
                _keys,
            ]

        sumlist = []
        ops = [
            ops,
        ] if not isinstance(ops, (list, tuple)) else ops

        for op in ops:
            sumdict = {}
            for k in _keys:
                values = self[k][1]
                if indices is None:
                    sumdict[k] = op(values, axis=1)
                else:
                    sumdict[k] = op(values[indices], axis=1)
            sumlist.append(sumdict)
        if len(ops) == 1:
            return sumlist[0]
        else:
            return sumlist

    @classmethod
    def _explain(cls):
        '''copied from cpptraj doc
        '''
        return '''
        [shear] Base pair shear.
        [stretch] Base pair stretch.
        [stagger] Base pair stagger.
        [buckle] Base pair buckle.
        [prop] Base pair propeller.
        [open] Base pair opening.
        [hb] Number of hydrogen bonds between bases in base pair.
        [pucker] Base sugar pucker.
        [major] Rough estimate of major groove width, calculated between P atoms of each
        base.
        [minor] Rough estimate of minor groove width, calculated between O4 atoms of
        each base.
        [shift] Base pair step shift.
        [slide] Base pair step slide.
        [rise] Base pair step rise.
        [title] Base pair step tilt.
        [roll] Base pair step roll.
        [twist] Base pair step twist.
        [xdisp] Helical X displacement.
        [ydisp] Helical Y displacement.
        [hrise] Helical rise.
        [incl] Helical inclination.
        [tip] Helical tip.
        [htwist] Helical twist.
        '''

    def __setstate__(self, state):
        self._dict = state

    def __getstate__(self):
        return self._dict
