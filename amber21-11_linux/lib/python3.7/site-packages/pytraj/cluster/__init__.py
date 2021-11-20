from __future__ import absolute_import
from collections import Counter
import numpy as np
from pytraj.utils.get_common_objects import get_topology, get_data_from_dtype
from pytraj.utils.decorators import register_pmap, register_openmp
from pytraj.utils.context import capture_stdout
from pytraj.utils.get_common_objects import super_dispatch, get_iterator_from_dslist
from pytraj.analysis.c_analysis import c_analysis
from pytraj.datasets.c_datasetlist import DatasetList as CpptrajDatasetList
from pytraj.datafiles.datafiles import DataFileList

__all__ = [
    'cluster',
    'kmeans',
    'dbscan',
    'hieragglo'
    'cluster_dataset',
]


class ClusteringDataset(object):
    '''
    Notes
    -----
    Unstable API

    Parameters
    ----------
    cpp_out : Tuple[CpptrajDatasetList, str]

    Attributes
    ----------
    cluster_index : np.ndarray[int]
        Cluster index of each frame
    n_frames : int
        Total frame
    population : Counter([int, int])
        Number of frames for each cluster
    fraction : np.ndarray[float]
        Fraction of each cluster
    centroids : np.ndarray[int]
        Representative frame index for each cluster
    '''

    def __init__(self, cpp_out):
        self._cpp_out = cpp_out

    def summary(self):
        return self._cpp_out[1]

    @property
    def cluster_index(self):
        return self._cpp_out[0]

    @property
    def n_frames(self):
        return sum(val for _, val in self.population.items())

    @property
    def population(self):
        return Counter(self.cluster_index)

    @property
    def fraction(self):
        return np.array(
            sorted(
                [
                    float(val) / self.n_frames
                    for _, val in self.population.items()
                ],
                reverse=True))

    @property
    def centroids(self):
        words = '#Representative frames:'
        for line in self.summary().split('\n'):
            if line.startswith(words):
                line = line.strip(words)
                return np.array(line.split(), dtype='i4')


def _cluster(traj,
             algorithm,
             mask="",
             frame_indices=None,
             dtype='dataset',
             top=None,
             options=''):
    """clustering. Limited support.

    Parameters
    ----------
    traj : Trajectory-like or any iterable that produces Frame
    mask : str
        atom mask
    dtype : str
        return data type
    top : Topology, optional
    options: str
        more cpptraj option

    Notes
    -----
    Call `pytraj._verbose()` to see more output. Turn it off by `pytraj._verbose(False)`


    cpptraj manual::

        Algorithms:
          [hieragglo [epsilon <e>] [clusters <n>] [linkage|averagelinkage|complete]
            [epsilonplot <file>]]
          [dbscan minpoints <n> epsilon <e> [sievetoframe] [kdist <k> [kfile <prefix>]]]
          [dpeaks epsilon <e> [noise] [dvdfile <density_vs_dist_file>]
            [choosepoints {manual | auto}]
            [distancecut <distcut>] [densitycut <densitycut>]
            [runavg <runavg_file>] [deltafile <file>] [gauss]]
          [kmeans clusters <n> [randompoint [kseed <seed>]] [maxit <iterations>]
          [{readtxt|readinfo} infofile <file>]
        Distance metric options: {rms | srmsd | dme | data}
          { [[rms | srmsd] [<mask>] [mass] [nofit]] | [dme [<mask>]] |
             [data <dset0>[,<dset1>,...]] }
          [sieve <#> [random [sieveseed <#>]]] [loadpairdist] [savepairdist] [pairdist <name>]
          [pairwisecache {mem | none}]
        Output options:
          [out <cnumvtime>] [gracecolor] [summary <summaryfile>] [info <infofile>]
          [summarysplit <splitfile>] [splitframe <comma-separated frame list>]
          [clustersvtime <filename> cvtwindow <window size>]
          [cpopvtime <file> [normpop | normframe]] [lifetime]
          [sil <silhouette file prefix>]
        Coordinate output options:
          [ clusterout <trajfileprefix> [clusterfmt <trajformat>] ]
          [ singlerepout <trajfilename> [singlerepfmt <trajformat>] ]
          [ repout <repprefix> [repfmt <repfmt>] [repframe] ]
          [ avgout <avgprefix> [avgfmt <avgfmt>] ]
        Experimental options:
          [[drawgraph | drawgraph3d] [draw_tol <tolerance>] [draw_maxit <iterations]]
        Cluster structures based on coordinates (RMSD/DME) or given data set(s).
        <crd set> can be created with the 'createcrd' command.
    """
    # Note: do not use super_dispatch here. We use get_iterator_from_dslist

    ana = c_analysis.Analysis_Clustering()
    # need to creat `dslist` here so that every time `do_clustering` is called,
    # we will get a fresh one (or will get segfault)
    crdname = 'DEFAULT_NAME'
    dslist, _top, mask2 = get_iterator_from_dslist(
        traj, mask, frame_indices, top, crdname=crdname)
    dflist = DataFileList()

    if 'summary' not in options.split():
        options += ' summary'

    # do not output cluster info to STDOUT
    command = ' '.join((algorithm, mask2, "crdset {0}".format(crdname),
                       options))

    with capture_stdout() as (out, _):
        ana(command, dslist, dflist=dflist)

    # remove frames in dslist to save memory
    dslist.remove_set(dslist[crdname])
    dflist.write_all_datafiles()
    return ClusteringDataset((get_data_from_dtype(dslist[:1], dtype='ndarray'),
                              out.read()))


def hieragglo(traj=None, mask="", options='', dtype='dataset'):
    return _cluster(
        traj=traj,
        algorithm='hieragglo',
        mask=mask,
        dtype=dtype,
        top=traj.top,
        options=options)


def dbscan(traj=None, mask="", options='', dtype='dataset'):
    return _cluster(
        traj=traj,
        algorithm='dbscan',
        mask=mask,
        dtype=dtype,
        top=traj.top,
        options=options)


def dpeaks(traj=None, mask="", options='', dtype='dataset'):
    return _cluster(
        traj=traj,
        algorithm='dpeaks',
        mask=mask,
        dtype=dtype,
        top=traj.top,
        options=options)


def readinfo(traj=None, mask="", options='', dtype='dataset'):
    return _cluster(
        traj=traj,
        algorithm='readinfo',
        mask=mask,
        dtype=dtype,
        top=traj.top,
        options=options)


hieragglo.__doc__ = _cluster.__doc__ + """

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2()
    >>> data = pt.cluster.hieragglo(traj, mask='@CA', options='')
"""

dbscan.__doc__ = _cluster.__doc__ + """

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2()
    >>> data = pt.cluster.dbscan(traj, mask='@CA', options='epsilon 1.7 minpoints 5')
"""

dpeaks.__doc__ = _cluster.__doc__
readinfo.__doc__ = _cluster.__doc__


@super_dispatch()
@register_openmp
def kmeans(traj=None,
           mask='*',
           n_clusters=10,
           random_point=True,
           kseed=1,
           maxit=100,
           metric='rms',
           top=None,
           frame_indices=None,
           options='',
           dtype='ndarray'):
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
    frame_indices : {None, 1D array-like}, optional
        if not None, only perform clustering for given indices. Notes that this is
        different from ``sieve`` keywords.
    options : str, optional
        extra cpptraj options controlling output, sieve, ...

    Sieve options::

        [sieve <#> [random [sieveseed <#>]]]

    Output options::

        [out <cnumvtime>] [gracecolor] [summary <summaryfile>] [info <infofile>]
        [summarysplit <splitfile>] [splitframe <comma-separated frame list>]
        [clustersvtime <filename> cvtwindow <window size>]
        [cpopvtime <file> [normpop | normframe]] [lifetime]
        [sil <silhouette file prefix>]

    Coordinate output options::

        [ clusterout <trajfileprefix> [clusterfmt <trajformat>] ]
        [ singlerepout <trajfilename> [singlerepfmt <trajformat>] ]
        [ repout <repprefix> [repfmt <repfmt>] [repframe] ]
        [ avgout <avgprefix> [avgfmt <avgfmt>] ]

    Notes
    -----
    - if the distance matrix is large (get memory Error), should add sieve number to
    ``options`` (check example)
    - install ``libcpptraj`` with ``-openmp`` flag to speed up this calculation.

    Returns
    -------
    1D numpy array of frame indices

    Examples
    --------
    >>> import pytraj as pt
    >>> from pytraj.cluster import kmeans
    >>> traj = pt.datafiles.load_tz2()
    >>> # use default options
    >>> cluster_data = kmeans(traj)
    >>> cluster_data.cluster_index
    array([8, 8, 6, ..., 0, 0, 0], dtype=int32)
    >>> cluster_data.centroids
    array([95, 34, 42, 40, 71, 10, 12, 74,  1, 64], dtype=int32)
    >>> # update n_clusters
    >>> data = kmeans(traj, n_clusters=5)
    >>> # update n_clusters with CA atoms
    >>> data = kmeans(traj, n_clusters=5, mask='@CA')
    >>> # specify distance metric
    >>> data = kmeans(traj, n_clusters=5, mask='@CA', kseed=100, metric='dme')
    >>> # add sieve number for less memory
    >>> data = kmeans(traj, n_clusters=5, mask='@CA', kseed=100, metric='rms', options='sieve 5')
    >>> # add sieve number for less memory, and specify random seed for sieve
    >>> data = kmeans(traj, n_clusters=5, mask='@CA', kseed=100, metric='rms', options='sieve 5 sieveseed 1')
    '''
    # don't need to get_topology
    _clusters = 'kmeans clusters ' + str(n_clusters)
    _random_point = 'randompoint' if random_point else ''
    _kseed = 'kseed ' + str(kseed)
    _maxit = str(maxit)
    _metric = metric
    # turn of cpptraj's cluster info
    _output = options
    command = ' '.join((_clusters, _random_point, _kseed, _maxit, _metric,
                        _output))
    return _cluster(
        traj,
        mask,
        frame_indices=frame_indices,
        top=top,
        dtype=dtype,
        options=command)


def cluster_dataset(array_like, options=''):
    '''cluster dataset

    Parameters
    ----------
    array_like : array_like
    options : str, cpptraj options

    Returns
    -------
    cluster index for each data point

    Examples
    --------
    >>> import pytraj as pt
    >>> import numpy as np
    >>> array_like = np.random.randint(0, 10, 1000)
    >>> data = pt.cluster.cluster_dataset(array_like, 'clusters 10 epsilon 3.0')
    '''
    import numpy as np
    c_dslist = CpptrajDatasetList()
    c_dslist.add('double', '__array_like')
    c_dslist[0].resize(len(array_like))
    c_dslist[0].values[:] = array_like
    act = c_analysis.Analysis_Clustering()
    command = 'data __array_like ' + options
    act(command, dslist=c_dslist)

    return np.array(c_dslist[-2])
