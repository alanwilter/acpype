import numpy as np

from ..utils.get_common_objects import (
    get_topology, get_data_from_dtype, get_matrix_from_dataset, get_reference,
    get_fiterator, super_dispatch, get_iterator_from_dslist)
from ..utils.convert import array_to_cpptraj_atommask
from ..utils.decorators import register_pmap, register_openmp
from .c_action import c_action
from .c_action import do_action
from .c_analysis import c_analysis
from .c_action.actionlist import ActionList
from ..datasets.datasetlist import DatasetList
from ..datasets.c_datasetlist import DatasetList as CpptrajDatasetList

__all__ = [
    'rotation_matrix',
    'pairwise_rmsd',
    'rmsd_perres',
    'rmsd_nofit',
    'rmsd',
    'symmrmsd',
]


@register_pmap
@super_dispatch()
def rotation_matrix(traj=None,
                    mask="",
                    ref=0,
                    mass=False,
                    frame_indices=None,
                    top=None,
                    with_rmsd=False):
    ''' Compute rotation matrix with/without rmsd

    Parameters
    ----------
    traj : Trajectory-like
    ref : {int, Frame}, default 0 (first Frame)
        reference
    mask : str, default all atoms
    mass : bool, default False
        if True, rmsfit with mass weighted
    frame_indices : {None, array-like}
        if not None, compute for given indices
    top : Topology, optional
    with_rmsd : bool, default False
        - if False, return only rotation matrix.
        - if True, return rotation matrix and rmsd values

    Returns
    -------
    out : if with_rmsd=False, return numpy array, shape (n_frames, 3, 3)
          if with_rmsd=True, return a tuple (mat, rmsd)

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2()
    >>> mat = pt.calc_rotation_matrix(traj, mask='@CA')
    >>> mat.shape
    (101, 3, 3)
    '''
    c_dslist = CpptrajDatasetList()
    mass_ = 'mass' if mass else ''

    command = ' '.join(('tmp', mask, 'savematrices', mass_))

    act = c_action.Action_Rmsd()
    act(command, [ref, traj], top=top, dslist=c_dslist)
    mat = c_dslist[-1].values
    # exclude data for reference
    if with_rmsd:
        return mat[1:], np.array(c_dslist[0].values[1:])
    else:
        return mat[1:]


@register_openmp
def pairwise_rmsd(traj=None,
                  mask="",
                  metric='rms',
                  top=None,
                  dtype='ndarray',
                  mat_type='full',
                  frame_indices=None):
    """ Calculate pairwise rmsd with different metrics.

    Parameters
    ----------

    traj : Trajectory-like or iterable object
    mask : mask
        if mask is "", use all atoms
    metric : {'rms', 'dme', 'srmsd', 'nofit'}
        if 'rms', perform rms fit
        if 'dme', use distance RMSD
        if 'srmsd', use symmetry-corrected RMSD
        if 'nofit', perform rmsd without fitting
    top : Topology, optional, default=None
    dtype: ndarray
        return type
    mat_type : str, {'full', 'half'}
        if 'full': return 2D array, shape=(n_frames, n_frames)
        if 'half': return 1D array, shape=(n_frames*(n_frames-1)/2, )

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> arr = pt.pairwise_rmsd(traj(0, 1000, mask='@CA'))
    >>> # calculate pairwise rmsd for all frames using CA atoms, use `dme` (distance RMSD)
    >>> # convert to numpy array
    >>> arr_np = pt.pairwise_rmsd(traj, "@CA", metric="dme", dtype='ndarray')
    >>> # calculate pairwise rmsd for all frames using CA atoms, nofit for RMSD
    >>> # convert to numpy array
    >>> arr_np = pt.pairwise_rmsd(traj, "@CA", metric="nofit", dtype='ndarray')
    >>> # calculate pairwise rmsd for all frames using CA atoms
    >>> # use symmetry-corrected RMSD, convert to numpy array
    >>> arr_np = pt.pairwise_rmsd(traj, "@CA", metric="srmsd", dtype='ndarray')
    >>> # use different dtype
    >>> arr_np = pt.pairwise_rmsd(traj, "@CA", metric="srmsd", dtype='dataset')

    Notes
    -----
    Install ``libcpptraj`` with ``openmp`` to get benefit from parallel
    """
    # we copy Frame coordinates to DatasetCoordsCRD first

    if not isinstance(mask, str):
        mask = array_to_cpptraj_atommask(mask)

    act = c_analysis.Analysis_Rms2d()

    crdname = 'default_coords'
    c_dslist, top_, command = get_iterator_from_dslist(
        traj, mask, frame_indices, top, crdname=crdname)

    command = ' '.join((command, metric,
                        "crdset {} rmsout mycrazyoutput".format(crdname)))

    act(command, dslist=c_dslist)
    # remove dataset coords to free memory
    c_dslist.remove_set(c_dslist[0])

    if dtype == 'ndarray':
        return get_matrix_from_dataset(c_dslist[0], mat_type)
    else:
        return get_data_from_dtype(c_dslist, dtype)


calc_pairwise_rmsd = pairwise_rmsd


@register_pmap
def rmsd_perres(traj=None,
                mask="",
                ref=0,
                mass=False,
                resrange=None,
                perres_mask=None,
                perres_center=False,
                perres_invert=False,
                frame_indices=None,
                top=None,
                dtype='dataset',
                **kwd):
    """superpose ``traj`` to ``ref`` with `mask`, then calculate nofit rms for residues
    in `resrange` with given `perresmask`

    Returns
    -------
    out : pytraj.DatasetList, shape=(1+n_residues, n_frames)
        out[0]: regular rmsd
        out[1:]: perres rmsd for all given residues
        `out.values` will return corresponding numpy array
    """
    range_ = 'range %s ' % resrange
    perresmask_ = 'perresmask ' + perres_mask if perres_mask is not None else ''
    perrestcenter_ = 'perrescenter' if perres_center else ''
    perrestinvert_ = 'perresinvert' if perres_invert else ''

    cm = " ".join((mask, 'perres', range_, perresmask_, perrestcenter_,
                   perrestinvert_))
    return rmsd(
        traj=traj,
        mask=cm,
        ref=ref,
        nofit=False,
        mass=mass,
        frame_indices=frame_indices,
        top=top,
        dtype=dtype,
        **kwd)


@register_pmap
def rmsd_nofit(traj=None,
               mask="",
               ref=0,
               mass=False,
               frame_indices=None,
               top=None,
               dtype='ndarray',
               **kwd):
    '''compute rmsd without fitting (translating and rotating)

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
    ref : Frame or int
    mass : bool, default False
        if True, use mass-weighted
    frame_indices : 1D array-like, default None
        if given, only perform calculation for those frames

    Notes
    -----
    This method is equal to pytraj.rmsd(traj, mask, ref, nofit=True, ...)

    When comparing the same structures (e.g. small ligands), the atoms need to be in the exact same order in the trajectory and reference frames.
    The function `atom_map(traj, ref, rmsfit=False)` can be used to attempt to reorder the atoms in the correct way before a RMSD calculation.
    '''
    return rmsd(
        traj=traj,
        mask=mask,
        ref=ref,
        mass=mass,
        nofit=True,
        frame_indices=frame_indices,
        top=top,
        dtype=dtype,
        **kwd)


calc_rmsd_nofit = rmsd_nofit


@register_pmap
def rmsd(traj=None,
         mask="",
         ref=0,
         ref_mask='',
         nofit=False,
         mass=False,
         update_coordinate=True,
         frame_indices=None,
         top=None,
         dtype='ndarray'):
    """compute rmsd

    Parameters
    ----------
    traj : Trajectory-like
    mask : str or 1D array-like of string or 1D or 2D array-like
        Atom mask/indices
    ref : {Frame, int}, default=0 (first frame)
        Reference frame or index.
    ref_mask: str, optional
        if given, use it instead of `mask`
    nofit : bool, default False
        if False, perform fitting (rotation and translation).
        if ``traj`` is mutable, its coordinates will be updated
        if True, not fitting.
    mass : bool, default False
        if True, compute mass-weighted rmsd
    update_coordinate : bool, default True
        if True, coordinates will be updated. But this only apply to mutable Trajectory
        if False (same as `nomod` in cpptraj), no modification
    frame_indices : int 1D array-like, default None
        if not None, only compute rmsd for given frame indices
    top : {Topology, str}, default None, optional
    dtype : return data type, default='ndarray'

    Notes
    -----
    - if traj and ref has diffrent n_atoms, make sure to update ref.top
    - you can use `pytraj.rmsd` to superpose structure (use update_coordinate=True)

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_trpcage()
    >>> # all atoms, do fitting, using ref=traj[-3]
    >>> data = pt.rmsd(traj, ref=-3)
    >>> # rmsd for 3 maskes, do fitting, using ref=traj[0] (defaul)
    >>> data = pt.rmsd(traj, mask=['@CA', '@C', ':3-18@CA'], dtype='dataset')
    >>> # rmsd to first frame, use mass ':3-13' but do not perorm fitting
    >>> data= pt.rmsd(traj, ref=traj[0], mask=':3-13', nofit=True)
    >>> # use atom indices for mask
    >>> data= pt.rmsd(traj, ref=traj[0], mask=range(40), nofit=True)
    >>> # compute rmsd (and align) with reference having different atoms
    >>> trpcage_traj = pt.datafiles.load_trpcage()[:]
    >>> tz2_traj = pt.datafiles.load_tz2()[:1]
    >>> data = pt.rmsd(trpcage_traj, mask='@1-10', ref=tz2_traj, ref_mask='@11-20')
    >>> data
    array([ 2.16203842,  2.28859396,  2.15817654, ...,  2.20767189,
            2.30087764,  1.92654945])
    Notes
    -----
    if ``traj`` is mutable and update_coordinate=True, its coordinates will be updated.

    When comparing the same structures (e.g. small ligands), the atoms need to be in the exact same order in the trajectory and reference frames.
    The function `atom_map(traj, ref, rmsfit=False)` can be used to attempt to reorder the atoms in the correct way before a RMSD calculation.
    """

    nofit_ = 'nofit' if nofit else ''
    mass_ = 'mass' if mass else ''
    nomod_ = 'nomod' if not update_coordinate else ''
    options = ' '.join((nofit_, mass_, nomod_))

    if ref_mask:
        if not mask:
            raise ValueError('mask must be provided if ref_mask is given')
        if not isinstance(ref_mask, str):
            ref_mask = array_to_cpptraj_atommask(ref_mask)

    if isinstance(mask, str):
        command = [
            mask,
        ]
    else:
        try:
            cmd = np.asarray(mask)
        except ValueError:
            raise ValueError("don't mix different types")
        dname = cmd.dtype.name
        if 'str' in dname:
            command = cmd
        elif 'int' in dname:
            if cmd.ndim == 1:
                command = [
                    array_to_cpptraj_atommask(mask),
                ]
            else:
                # assume ndim==2
                command = [array_to_cpptraj_atommask(x) for x in mask]
        elif 'object' in dname:
            # different array lens or mix type
            # dangerous: assume array of two array
            command = [array_to_cpptraj_atommask(x) for x in mask]
        else:
            raise ValueError("not supported")

    top_ = get_topology(traj, top)

    ref = get_reference(traj, ref)
    fi = get_fiterator(traj, frame_indices)

    alist = ActionList()
    c_dslist = CpptrajDatasetList()

    ref_top = ref.top if ref.top else top_
    c_dslist.add('reference', name='myref')
    c_dslist[-1].top = ref_top
    c_dslist[-1].add_frame(ref)

    ref_mask = ' '.join((ref_mask, 'ref myref'))

    for cm in command:
        cm_ = ' '.join((cm, ref_mask, options))
        if 'savematrices' in cm_ and dtype not in [
                'dataset', 'cpptraj_dataset'
        ]:
            raise ValueError('if savematrices, dtype must be "dataset"')
        alist.add(c_action.Action_Rmsd(), cm_, top=top_, dslist=c_dslist)

    alist.compute(fi)

    # pop Reference Dataset
    c_dslist._pop(0)

    dnew = DatasetList(c_dslist)
    return get_data_from_dtype(dnew, dtype=dtype)


@super_dispatch()
def symmrmsd(traj,
             mask='',
             ref=0,
             ref_mask=None,
             fit=True,
             remap=False,
             mass=False,
             top=None,
             dtype='ndarray',
             frame_indices=None):
    """Compute symmetry-corrected RMSD

    Parameters
    ----------
    traj : Trajectory-like
    mask : str, default '' (all atoms)
    ref : {int, Frame}, default 0 (first frame)
    ref_mask : {str, None}, default None
        if None, use traj's mask
        if given, use it
    fit : Bool, default True
        if True, do fitting
        if False, nofit
    mass : Bool, default False
        if True, compute mass-weighted rmsd
        if False, no mas-weighted
    remap : Bool, default False
        if True, frames will be modifed for symmetry as well
    dtype : str, default 'ndarray'
        return data type
    frame_indices : {None, array-like}, default None
       if given, only compute RMSD for those

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.load("TYR.nc", "TYR.parm7") # doctest: +SKIP
    >>> data = pt.symmrmsd(traj, ref=0) # doctest: +SKIP
    Notes
    -----
    versionadded: 1.0.6

    When comparing the same structures (e.g. small ligands), the atoms need to be in the exact same order in the trajectory and reference frames.
    The function `atom_map(traj, ref, rmsfit=False)` can be used to attempt to reorder the atoms in the correct way before a RMSD calculation.
    """

    mask_ = mask
    refmask_ = ref_mask if ref_mask is not None else ''
    reftop = ref.top if hasattr(ref, 'top') else top
    nofit_ = 'nofit' if not fit else ''
    mass_ = 'mass' if mass else ''
    remap_ = 'remap' if remap else ''

    refname = 'myref'
    ref_command_ = 'ref {}'.format(refname)

    command = ' '.join((mask_, refmask_, nofit_, mass_, remap_, ref_command_))

    if reftop is None:
        reftop = traj.top

    c_dslist = CpptrajDatasetList()
    c_dslist.add('reference', name=refname)
    c_dslist[0].top = reftop
    c_dslist[0].add_frame(ref)

    act = c_action.Action_SymmetricRmsd()
    act.read_input(command, top=top, dslist=c_dslist)
    act.setup(top)

    for frame in traj:
        new_frame = act.compute(frame, get_new_frame=remap)
        if remap:
            frame.xyz[:] = new_frame.xyz[:]
    act.post_process()

    # remove ref
    c_dslist._pop(0)

    return get_data_from_dtype(c_dslist, dtype=dtype)


@register_pmap
@super_dispatch()
def distance_rmsd(traj=None,
                  mask='',
                  ref=0,
                  top=None,
                  dtype='ndarray',
                  frame_indices=None):
    '''compute distance rmsd between traj and reference

    Parameters
    ----------
    traj : Trajectory-like or iterator that produces Frame
    ref : {int, Frame}, default 0 (1st Frame)
    mask : str
    top : Topology or str, optional, default None
    dtype : return dtype, default 'ndarray'

    Returns
    -------
    1D ndarray if dtype is 'ndarray' (default)

    Examples
    --------
    >>> import pytraj as pt
    >>> # compute distance_rmsd to last frame
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> data = pt.distance_rmsd(traj, ref=-1)
    >>> # compute distance_rmsd to first frame with mask = '@CA'
    >>> data = pt.distance_rmsd(traj, ref=0, mask='@CA')
    '''
    command = mask
    c_dslist, _ = do_action(
        [ref, traj], command, c_action.Action_DistRmsd, top=top)
    # exclude ref value
    for d in c_dslist:
        d.data = d.data[1:]
    return get_data_from_dtype(c_dslist, dtype=dtype)
