from __future__ import absolute_import
import numpy as np
from pytraj import DatasetList, tools
from ..utils.get_common_objects import get_data_from_dtype, super_dispatch, get_topology
from ..utils.decorators import register_openmp
from ..datasets.c_datasetlist import DatasetList as CpptrajDatasetList
from .c_action.c_action import Action_DSSP

__all__ = ['dssp', 'dssp_allatoms', 'dssp_allresidues']


@register_openmp
@super_dispatch()
def dssp(traj=None,
         mask="",
         frame_indices=None,
         dtype='ndarray',
         simplified=False,
         top=None):
    """return dssp profile for frame/traj

    Parameters
    ----------
    traj : Trajectory-like
    mask: str
        atom mask
    frame_indices : {None, array-like}, default None, optional
        specify frame numbers for calculation.
        if None, do all frames
    dtype : str, default 'ndarray'
        return data type, for regular user, just use default one (ndarray).
        use dtype='dataset' if wanting to get secondary structure in integer format
    simplified : bool, default False
        if True, use simplified codes, only has 'H', 'E' and 'C'
        if False, use all DSSP codes

    Returns
    -------
    out_0: ndarray, shape=(n_residues,)
        residue names
    out_1: ndarray, shape=(n_frames, n_residues)
        DSSP for each residue
    out_2 : pytraj.DatasetList
        average value for each secondary structure type

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.load_pdb_rcsb('1l2y')
    >>> residues, ss, _ = pt.dssp(traj, ":2-10")
    >>> residues # doctest: +SKIP
    array(['LEU:2', 'TYR:3', 'ILE:4', 'GLN:5', 'TRP:6', 'LEU:7', 'LYS:8',
           'ASP:9', 'GLY:10'],
          dtype='<U6')
    >>> ss # doctest: +SKIP
    array([['0', 'H', 'H', ..., 'H', 'T', '0'],
           ['0', 'H', 'H', ..., 'H', 'T', '0'],
           ['0', 'H', 'H', ..., 'H', 'T', '0'],
           ...,
           ['0', 'H', 'H', ..., 'H', 'T', '0'],
           ['0', 'H', 'H', ..., 'H', 'H', '0'],
           ['0', 'H', 'H', ..., 'H', 'T', '0']],
          dtype='<U1')

    >>> residues, ss, _ = pt.dssp(traj, mask=range(100))

    >>> traj = pt.fetch_pdb('1l2y')
    >>> residues, ss, _ = pt.dssp(traj, simplified=True)
    >>> ss[0].tolist() # first frame
    ['C', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'C', 'C', 'H', 'H', 'H', 'H', 'C', 'C', 'C', 'C', 'C', 'C']


    Notes
    -----
    ========= ======= ========= =======================
    Character Integer DSSP_Char Secondary structure type
    ========= ======= ========= =======================
    0         0       '0'       None
    b         1       'E'       Parallel Beta-sheet
    B         2       'B'       Anti-parallel Beta-sheet
    G         3       'G'       3-10 helix
    H         4       'H'       Alpha helix
    I         5       'I'       Pi (3-14) helix
    T         6       'T'       Turn
    S         7       'S'       Bend
    ========= ======= ========= =======================

    Simplified codes::

        - 'H': include 'H', 'G', 'I' (helix)
        - 'E': include 'E', 'B' (strand)
        - 'C': include 'T', 'S' or '0' (coil)

    Simplified codes will be mostly used for visualization in other packages.
    """

    command = mask

    dslist = CpptrajDatasetList()

    Action_DSSP()(command, traj, top=top, dslist=dslist)

    # replace key to something nicer
    for key, dset in dslist.iteritems():
        if 'DSSP' in key:
            key = key.replace("DSSP_00000[", "")
            key = key.replace("]", "_avg")
            dset.key = key.lower()
    dtype = dtype.lower()

    if dtype == 'ndarray':
        # get all dataset from DatSetList if dtype == integer
        arr0 = dslist.grep("integer", mode='dtype').values
        keys = dslist.grep("integer", mode='dtype').keys()
        avg_dict = DatasetList(dslist.grep('_avg'))
        ss_array = np.asarray([
            _to_string_secondary_structure(arr, simplified=simplified)
            for arr in arr0
        ]).T
        return np.asarray(keys), ss_array, avg_dict
    else:
        return get_data_from_dtype(dslist, dtype=dtype)


# _s0 = ['None', 'Para', 'Anti', '3-10', 'Alpha', 'Pi', 'Turn', 'Bend']
_s1 = ["0", "b", "B", "G", "H", "I", "T", "S"]


def _to_string_secondary_structure(arr0, simplified=False):
    """
    arr0 : ndarray
    """
    if not simplified:
        ssdict = {
            0: '0',
            1: 'b',
            2: 'B',
            3: 'G',
            4: 'H',
            5: 'I',
            6: 'T',
            7: 'S'
        }
    else:
        ssdict = {
            0: 'C',
            1: 'E',
            2: 'E',
            3: 'H',
            4: 'H',
            5: 'H',
            6: 'C',
            7: 'C'
        }

    return np.vectorize(lambda key: ssdict[key])(arr0)


def get_ss_per_frame(arr, top, res_indices, simplified=False, all_atoms=False):
    if simplified:
        symbol = 'C'

    for idx, res in enumerate(top.residues):
        if idx in res_indices:
            ss = arr[res_indices.index(idx)]
            if all_atoms:
                yield [
                    ss
                    for _ in range(res.first_atom_index, res.last_atom_index)
                ]
            else:
                # only residues
                yield [
                    ss,
                ]
        else:
            if all_atoms:
                yield [
                    symbol
                    for _ in range(res.first_atom_index, res.last_atom_index)
                ]
            else:
                yield [
                    symbol,
                ]


def dssp_allatoms(traj, *args, **kwd):
    '''calculate dssp for all atoms

    Returns
    -------
    ndarray, shape=(n_frames, n_atoms)

    Notes
    -----
    this method is not well optimized for speed.

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.fetch_pdb('1l2y')
    >>> x = pt.dssp_allatoms(traj, simplified=True)
    >>> x[0, :3].tolist()
    ['C', 'C', 'C']

    See also
    --------
    dssp
    '''
    res_labels, data = dssp(traj, *args, **kwd)[:2]
    top = get_topology(traj, kwd.get('top'))
    res_indices = [int(x.split(':')[-1]) - 1 for x in res_labels]

    new_data = np.empty((traj.n_frames, traj.n_atoms), dtype='U2')
    simplified = kwd.get('simplified', False)
    for fid, arr in enumerate(data):
        new_data[fid][:] = tools.flatten(
            get_ss_per_frame(
                arr, top, res_indices, simplified, all_atoms=True))
    return new_data


def dssp_allresidues(traj, *args, **kwd):
    '''calculate dssp for all residues. Mostly used for visualization.

    Returns
    -------
    ndarray, shape=(n_frames, n_residues)

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_dpdp()
    >>> x = pt.dssp_allresidues(traj, simplified=True)
    >>> x[0].tolist()
    ['C', 'E', 'E', 'E', 'E', 'C', 'C', 'C', 'C', 'E', 'E', 'E', 'E', 'E', 'C', 'C', 'E', 'E', 'E', 'E', 'C', 'C']
    >>> len(x[0]) == traj.top.n_residues
    True

    >>> # load trajectory having waters
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> x = pt.dssp_allresidues(traj, simplified=True)
    >>> len(x[0]) == traj.top.n_residues
    True
    >>> len(x[0])
    1704
    >>> # only calculate protein residues, use `pytraj.dssp`
    >>> y = pt.dssp(traj, simplified=True)
    >>> len(y[0])
    13

    Notes
    -----
    this method is not well optimized for speed.

    See also
    --------
    dssp
    '''
    res_labels, data = dssp(traj, *args, **kwd)[:2]
    top = get_topology(traj, kwd.get('top', None))

    # do not need to compute again if there is no solvent or weird residues
    if len(res_labels) == top.n_residues:
        return data

    res_indices = [int(x.split(':')[-1]) - 1 for x in res_labels]

    new_data = np.empty((traj.n_frames, traj.top.n_residues), dtype='U2')

    simplified = kwd.get('simplified', False)
    for fid, arr in enumerate(data):
        new_data[fid][:] = tools.flatten(
            get_ss_per_frame(
                arr, top, res_indices, simplified, all_atoms=False))
    return new_data
