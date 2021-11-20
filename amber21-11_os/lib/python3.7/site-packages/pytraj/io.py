import os
from pathlib import Path
import tempfile
import numpy as np

from .core.c_core import _load_batch
from .serialize.serialize import to_pickle, read_pickle
from .datafiles.load_samples import load_sample_data
from .core.c_options import set_error_silent
from .topology.topology import Topology, ParmFile
from .trajectory.shared_methods import iterframe_master
from .trajectory.trajectory import Trajectory
from .trajectory.trajectory_iterator import TrajectoryIterator
from .trajectory.frameiter import iterframe
from .trajectory.c_traj.c_trajout import TrajectoryWriter
from .utils.decorators import ensure_exist
from .utils.context import tempfolder

try:
    from urllib.request import urlopen
except ImportError:
    from urllib import urlopen

__all__ = [
    'load',
    'iterload',
    'load_remd',
    'iterload_remd',
    'load_pdb_rcsb',
    'load_sample_data',
    'load_parmed',
    'load_leap',
    'load_topology',
    'write_parm',
    'save',
    'write_traj',
    'read_pickle',
    'to_pickle',
    'select_atoms',
]


def load(filename, top=None, frame_indices=None, mask=None, stride=None):
    """try loading and returning appropriate values. See example below.

    Parameters
    ----------
    filename : str, Trajectory filename
    top : Topology filename or a Topology
    frame_indices : {None, array_like}, default None
        only load frames with given number given in frame_indices
    stride : {None, int}, default None
        if given, frame will be skip every `stride`.
        Note: if both frame_indices and stride are given, `frame_indices` will be ignored.
    mask : {str, None}, default None
        if None: load coordinates for all atoms
        if string, load coordinates for given atom mask

    Returns
    -------
    pytraj.Trajectory

    Notes
    -----
    - For further slicing options, see pytraj.TrajectoryIterator (created by ``pytraj.iterload``)

    - Also see `pytraj.iterload` for loading a series of trajectories that don't fit to
      memory

    See also
    --------
    iterload

    Examples
    --------
    >>> import pytraj as pt
    >>> # load netcdf file with given amber parm file
    >>> traj = pt.load('traj.nc', '2koc.parm7') # doctest: +SKIP

    >>> # load mol2 file
    >>> traj = pt.load('traj.mol2') # # doctest: +SKIP

    >>> # load pdb file
    >>> traj = pt.load('traj.pdb') # doctest: +SKIP

    >>> # load given frame numbers
    >>> from pytraj.testing import get_fn
    >>> fn, tn = get_fn('tz2_dry')
    >>> traj = pt.load(fn, top=tn, frame_indices=[0, 3, 5, 12, 20])
    >>> traj = pt.load(fn, top=tn, frame_indices=[0, 3, 5, 12, 20], mask='!@H=')
    >>> # load last frame
    >>> traj = pt.load(fn, top=tn, frame_indices=[-1])

    >>> # load with frame slice
    >>> traj = pt.load(fn, tn, frame_indices=slice(0, 10, 2))
    >>> # which is equal to:
    >>> traj = pt.load(fn, tn, frame_indices=range(0, 10, 2))

    >>> # load with stride
    >>> traj = pt.load(fn, tn)
    >>> traj.n_frames
    101
    >>> traj = pt.load(fn, tn, stride=5)
    >>> traj.n_frames
    21

    >>> # load with stride for more than one filename
    >>> traj = pt.load([fn, fn], tn, stride=5)
    >>> traj.n_frames
    42
    >>> traj.n_atoms
    223

    >>> # load with stride for more than one filename, and with mask
    >>> traj = pt.load([fn, fn], tn, stride=5, mask='@CA')
    >>> traj.n_frames
    42
    >>> traj.n_atoms
    12
    """
    # load to TrajectoryIterator object first
    # do not use frame_indices_ here so we can optimize the slicing speed
    traj = load_traj(filename, top, stride=stride)

    # do the slicing and other things if needed.
    if stride is not None:
        if mask is None:
            return traj[:]
        else:
            return traj[mask]
    else:
        frame_indices_ = frame_indices
        if isinstance(frame_indices_, tuple):
            frame_indices_ = list(frame_indices_)
        if frame_indices_ is None and mask is None:
            # load all
            return traj[:]
        elif frame_indices_ is None and mask is not None:
            # load all frames with given mask
            # eg: traj['@CA']
            return traj[mask]
        elif frame_indices_ is not None and mask is None:
            # eg. traj[[0, 3, 7]]
            return traj[frame_indices_]
        else:
            # eg. traj[[0, 3, 7], '@CA']
            return traj[frame_indices_, mask]


def iterload(*args, **kwd):
    """return TrajectoryIterator object

    Parameters
    ----------
    filename: {str, list-like of filenames, pattern}
        input trajectory filename(s). You can use a single filename, a list of filenames
        or a pattern.
    top : {str, Topology}
        input topology. If str, pytraj will load from disk to Topology first
    frame_slice: tuple or list of tuple
        specify start, stop, step for each trajectory you want to read.

        cpptraj input::

            trajin traj0.nc 1 10
            trajin traj1.nc

        In pytraj, corresponding frame_slice=[(0, 10), (0, -1)]
    stride : {None, int}, default None
        if not None, trajectories will be strided.
        Note: if both stride and frame_slice are not None, frame_slice will be ignored

    Returns
    -------
    pytraj.TrajectoryIterator

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.iterload('traj.nc', '2koc.parm7') # doctest: +SKIP

    >>> # load from a list of filenames
    >>> traj = pt.iterload(['traj0.nc', 'traj1.nc'], '2koc.parm7') # doctest: +SKIP

    >>> # load all files with a given pattern (sorted)
    >>> traj = pt.iterload('./traj*.nc', '2koc.parm7') # doctest: +SKIP

    >>> # load from a list of files with given frame step
    >>> # for each file, only register to load from frame 0 to 9 (skip 10), skip every 2 frames
    >>> traj = pt.iterload(['traj0.nc', 'traj1.nc'], '2koc.parm7', frame_slice=[(0, 10, 2),]*2) # doctest: +SKIP

    >>> # load from frame 0 to 9 for `traj0.nc`
    >>> # load all frames from `traj1.nc`
    >>> traj = pt.iterload(['traj0.nc', 'traj1.nc'], '2koc.parm7', frame_slice=[(0, 10), (0, -1)]) # doctest: +SKIP

    >>> # use stride, skip every 2 frames
    >>> from pytraj.testing import get_remd_fn
    >>> filenames, topology_filename = get_remd_fn('remd_ala2')
    >>> [fn.split('/')[-1] for fn in filenames]
    ['rem.nc.000', 'rem.nc.001', 'rem.nc.002', 'rem.nc.003']
    >>> traj = pt.iterload(filenames, topology_filename, stride=2)

    Notes
    -----
    Unlike `pytraj.load`, you can not arbitarily set `frame_indices`. If you want to do
    so, first load trajectories to TrajectoryIterator object, then do fancy slicing

    >>> import pytraj as pt
    >>> # register to load traj.nc from 0-th to 99-th frame
    >>> traj = pt.iterload('traj.nc', 'prmtop', frame_slice=(0, 100)]) # doctest: +SKIP
    >>> # do fancy indexing to load specific frames to memory
    >>> traj[[0, 8, 3, 50, 7]] # doctest: +SKIP

    >>> # load to disk with given mask
    >>> traj[[0, 8, 3, 50, 7], '!@H='] # doctest: +SKIP
    """
    if kwd and 'frame_indices' in kwd.keys():
        raise ValueError(
            "do not support indices for TrajectoryIterator loading")
    return load_traj(*args, **kwd)


def load_traj(filename=None, top=None, *args, **kwd):
    """load trajectory from filename

    Parameters
    ----------
    filename : str
    top : {str, Topology}
    frame_indices : {None, list, array ...}
    *args, **kwd: additional arguments, depending on `engine`

    Returns
    -------
    TrajectoryIterator : if frame_indices is None
    Trajectory : if there is indices
    """
    if isinstance(top, str):
        top = load_topology(top)
    if top is None or top.is_empty():
        top = load_topology(filename)
    ts = TrajectoryIterator(top=top)

    if 'stride' in kwd:
        ts._load(filename, stride=kwd['stride'])
    elif 'frame_slice' in kwd:
        ts._load(filename, frame_slice=kwd['frame_slice'])
    else:
        ts._load(filename)

    return ts


def _load_from_frame_iter(iterable, top=None):
    '''
    '''
    return Trajectory.from_iterable(iterable, top)


def iterload_remd(filename, top=None, T="300.0"):
    """Load temperature remd trajectory for single temperature.
    e.g: Suppose you have replica trajectoris remd.x.00{1-4}.
    You want to load and extract only frames at 300 K, use this method

    Parameters
    ----------
    filename : str
    top : {str, Topology}
    T : {float, str}, default=300.0

    Returns
    -------
    pytraj.traj.TrajectoryCpptraj

    Notes
    -----

    """
    from pytraj.core.c_core import CpptrajState, Command

    state = CpptrajState()

    # add keyword 'remdtraj' to trick cpptraj
    trajin = ' '.join(('trajin', filename, 'remdtraj remdtrajtemp', str(T)))
    if isinstance(top, str):
        top = load_topology(top)
    else:
        top = top
    state.data.add('topology', 'remdtop')
    # set topology
    state.data['remdtop']._top = top

    # load trajin
    with Command() as cm:
        cm.dispatch(state, trajin)
        cm.dispatch(state, 'loadtraj name remdtraj')

    # state.data.remove_set(state.data['remdtop'])
    traj = state.data[-1]

    # assign state to traj to avoid memory free
    traj._base = state
    return traj


def load_remd(filename, top=None, T="300.0"):
    return iterload_remd(filename, top, T)[:]


def _files_exist(filename, n_frames, options):
    """ Will be used for `write_traj` with options that will write
    multilple frames; example: x.pdb -> x.pdb.{1, 2, 3, ...}
    
    Note: does not support format 'x.{1, 2, ...}.pdb'

    Parameters
    ----------
    filename : str
    n_frames : int
    options : str
        cpptraj options
    """
    option_set = set([s for s in options.split() if s])
    exists = []
    if 'multi' in option_set:
        # 'multi' is only available with pdb format.
        for index in range(n_frames):
            if 'keepext' in option_set:
                # e.g: x.1.pdb
                words = filename.split('.')
                if len(words) == 1:
                    fn = filename + '.' + str(index+1)
                else:
                    ext = words[-1]
                    fn = "{}{}.{}".format(filename.strip(ext), index+1, ext)
            else:
                # e.g: x.pdb.1
                fn = filename + '.' + str(index+1)

            if os.path.exists(fn):
                exists.append(fn)
    else:
        if os.path.exists(filename):
            exists.append(filename)
    return exists


def write_traj(filename,
               traj,
               format='infer',
               frame_indices=None,
               overwrite=False,
               force=False,
               velocity=False,
               time=False,
               options=""):
    """

    Parameters
    ----------
    filename : str
    traj : Trajectory
    format : str, default 'infer'
        if 'inter', detect format based on extension. If can not detect, use amber mdcrd format.
    frame_indices: array-like or iterator that produces integer, default: None
        If not None, only write output for given frame indices
    overwrite: bool, default: False
        Note: does not respect options='keepext'
    velocity : bool, default False
        if True, write velocity. Make sure your trajectory or Frame does have velocity
    force : bool, default False
        if True, write force. Make sure your trajectory or Frame does have force
    time: bool, default False
        if True, write time.
    options : str, additional cpptraj keywords

    Notes
    -----
    ===================  =========
    Format               Extension
    ===================  =========
    Amber Trajectory     .crd
    Amber NetCDF         .nc
    Amber Restart        .rst7
    Amber NetCDF         .ncrst
    Charmm DCD           .dcd
    PDB                  .pdb
    Mol2                 .mol2
    Scripps              .binpos
    Gromacs              .trr
    SQM Input            .sqm
    ===================  =========

    'options' for writing to pdb format (cptraj manual)::

        dumpq:       Write atom charge/GB radius in occupancy/B-factor columns (PQR format)."
        parse:       Write atom charge/PARSE radius in occupancy/B-factor columns (PQR format)."
        vdw:         Write atom charge/VDW radius in occupancy/B-factor columns (PQR format)."
        pdbres:      Use PDB V3 residue names."
        pdbatom:     Use PDB V3 atom names."
        pdbv3:       Use PDB V3 residue/atom names."
        teradvance:  Increment record (atom) # for TER records (default no)."
        terbyres:    Print TER cards based on residue sequence instead of molecules."
        model:       Write to single file separated by MODEL records."
        multi:       Write each frame to separate files."
        chainid <c>: Write character 'c' in chain ID column."
        sg <group>:  Space group for CRYST1 record, only if box coordinates written."
        include_ep:  Include extra points."
        conect:      Write CONECT records using bond information.");

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> pt.write_traj("output/t.nc", traj, overwrite=True) # write to amber netcdf file

    >>> # write to multi pdb files (t.pdb.1, t.pdb.2, ...)
    >>> pt.write_traj("output/t.pdb", traj, overwrite=True, options='multi')

    >>> # write all frames to single pdb file and each frame is seperated by "MODEL" word
    >>> pt.write_traj("output/t.pdb", traj, overwrite=True, options='model')

    >>> # write to DCD file
    >>> pt.write_traj("output/test.dcd", traj, overwrite=True)

    'options' for writing to amber netcdf format (cptraj manual)::

        remdtraj: Write temperature to trajectory (makes REMD trajectory)."
        velocity: Write velocities to trajectory."
        force: Write forces to trajectory.");

    'options' for writing to amber netcdf restart format(cptraj manual)::

        novelocity: Do not write velocities to restart file."
        notime:     Do not write time to restart file."
        remdtraj:   Write temperature to restart file."
        time0:      Time for first frame (default 1.0)."
        dt:         Time step for subsequent frames, t=(time0+frame)*dt; (default 1.0)");
        keepext     Keep filename extension; write '<name>.<num>.<ext>' instead (example: myfile.1.rst7)

    'options' for writing to mol2 format (cptraj manual)::

        single   : Write to a single file."
        multi    : Write each frame to a separate file."
        sybyltype: Convert Amber atom types (if present) to SYBYL types.");

    'options'  for other formats::

        please check http://ambermd.org/doc12/Amber15.pdf
    """
    existing_files = _files_exist(filename, traj.n_frames, options)
    if existing_files and not overwrite:
        for fn in existing_files:
            print("{} exists. Use overwrite=True or remove the file".format(fn))
        raise IOError()

    if hasattr(traj, '_crdinfo'):
        crdinfo = traj._crdinfo
    else:
        crdinfo = dict()

    crdinfo['has_force'] = force
    crdinfo['has_velocity'] = velocity
    crdinfo['has_time'] = time

    with TrajectoryWriter(
            filename=filename,
            top=traj.top,
            format=format,
            crdinfo=crdinfo,
            options=options) as writer:
        for frame in iterframe(traj, frame_indices=frame_indices):
            writer.write(frame)


def write_parm(filename=None, top=None, format='amberparm', overwrite=False):
    if os.path.exists(filename) and not overwrite:
        raise RuntimeError(
            '{0} exists, must set overwrite=True'.format(filename))
    parm = ParmFile()
    parm.write(filename=filename, top=top, format=format)


def load_topology(filename, option=''):
    """load Topology from a filename or from url or from ParmEd object. Adapted from cpptraj doc.

    Parameters
    ----------
    filename : str, Amber prmtop, pdb, mol2, psf, cif, gromacs topology, sdf, tinker formats
    option : cpptraj options.
        if filename is a pdb file, option = {'pqr', 'noconnect'}.

        - pqr     : Read atomic charge/radius from occupancy/B-factor columns.
        - noconect: Do not read CONECT records if present.

    Notes
    -----
    if cpptraj/pytraj does not support specific file format, you still can convert to PDB
    file. cpptraj will do the bond guess based on distance.

    Examples
    --------
    >>> import pytraj as pt
    >>> # from a filename
    >>> pt.load_topology("data/tz2.ortho.parm7")
    <Topology: 5293 atoms, 1704 residues, 1692 mols, PBC with box type = ortho>

    >>> # read with option
    >>> pt.load_topology('1KX5.pdb', 'bondsearch 0.2') # doctest: +SKIP
    """
    top = Topology()

    # always read box info from pdb
    option = ' '.join(('readbox', option))

    if isinstance(filename, str):
        parm = ParmFile()
        set_error_silent(True)
        parm.read(filename=filename, top=top, option=option)
        set_error_silent(False)
    else:
        raise ValueError('filename must be a string')

    if top.n_atoms == 0:
        raise RuntimeError(
            'n_atoms = 0: make sure to load correct Topology filename '
            'or load supported topology (pdb, amber parm, psf, ...)')
    return top


def load_parmed(parm, traj=True, **kwd):
    """return pytraj's Topology or Trajectory objects

    Parameters
    ----------
    parm : ParmEd's Structure object
    traj: bool, default True
        if True, return pytraj.Trajectory
        if False, return Topology

    >>> import parmed as pmd
    >>> import pytraj as pt
    >>> p = pmd.download_PDB("1l2y")
    >>> traj = pt.load_parmed(p)
    """
    with tempfolder():
        fname = 'tmp.parm7'
        parm.save(fname)
        top = load_topology(fname)
    if traj:
        coords = parm.get_coordinates()
        traj = Trajectory(xyz=coords, top=top)
        traj.unitcells = parm.get_box()
        return traj
    else:
        return top


def loadpdb_rcsb(pdbid):
    """load pdb file from rcsb website

    Parameters
    ----------
    pdbid : str

    Examples
    --------
        io.loadpdb_rcsb("2KOC") # popular RNA hairpin
    """
    pdbid = pdbid.upper()
    url = 'http://files.rcsb.org/download/{pdbid}.pdb'.format(pdbid=pdbid)
    return _make_traj_from_remote_file(url)


def load_url(url):
    """
    
    versionadded: 1.0.7
    """
    return _make_traj_from_remote_file(url)


def _make_traj_from_remote_file(remote_file):
    _, fname = tempfile.mkstemp()
    content = urlopen(remote_file).read()
    Path(fname).write_text(content.decode())
    return load(fname)


def download_pdb(pdbid, location="./"):
    """download pdb to local disk
    """
    fname = location + pdbid + ".pdb"
    url = 'http://www.rcsb.org/pdb/files/%s.pdb' % pdbid
    content = urlopen(url).read().decode()
    Path(fname).write_text(content)


# create alias
download_PDB = download_pdb
load_pdb_rcsb = loadpdb_rcsb


@ensure_exist
def load_single_frame(filename=None, top=None, index=0):
    """load a single Frame"""
    return iterload(filename, top)[index]


load_frame = load_single_frame

# creat alias
save_traj = write_traj


def save(filename, obj, *args, **kwd):
    '''an universal method

    Parameters
    ----------
    filename : output filename
    obj : Topology or Trajetory-like
        if Topology, write a new Topology to disk
        if Trajetory-like, write a trajectory to disk
    *args, **kwd: additional options

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_ala3()
    >>> isinstance(traj, pt.TrajectoryIterator)
    True
    >>> top = traj.top
    >>> isinstance(top, pt.Topology)
    True
    >>> # save Topology to a new Topology
    >>> pt.save('output/prmtop', top, overwrite=True)
    >>> isinstance(pt.load_topology('output/prmtop'), pt.Topology)
    True
    >>> # save TrajectoryIterator to disk
    >>> pt.save('output/traj.nc', traj, overwrite=True)
    >>> isinstance(pt.iterload('output/traj.nc', traj.top), pt.TrajectoryIterator)
    True

    See also
    --------
    write_traj
    write_parm
    '''
    if isinstance(obj, Topology):
        write_parm(filename, obj, *args, **kwd)
    else:
        write_traj(filename, obj, *args, **kwd)


def get_coordinates(iterable,
                    autoimage=None,
                    rmsfit=None,
                    mask=None,
                    frame_indices=None):
    '''return 3D-ndarray coordinates of `iterable`, shape=(n_frames, n_atoms, 3). This method is more memory
    efficient if use need to perform autoimage and rms fit to reference before loading all coordinates
    from disk.

    This method is good (fast, memory efficient) if you just want to get raw numpy array
    to feed to external package, such as sciki-learn, ...

    Parameters
    ----------
    iterable : could be anything that produces Frame when iterating
               (Trajectory, TrajectoryIterator, FrameIterator, Python's generator, ...)

    Notes
    -----
    - if using both ``autoimage`` and ``rmsfit``, autoimage will be always processed before doing rmsfit.
    - You will get faster speed if ``iterable`` has attribute ``n_frames``

    Examples
    --------
    >>> import pytraj as pt
    >>> from pytraj.testing import get_fn
    >>> fn, tn = get_fn('tz2')
    >>> # load to out-of-core pytraj.TrajectoryIterator
    >>> traj = pt.iterload(fn, tn)
    >>> traj.n_frames, traj.n_atoms
    (10, 5293)

    >>> # simple
    >>> xyz = pt.get_coordinates(traj)  # same as traj.xyz
    >>> xyz.shape
    (10, 5293, 3)

    >>> # load coordinates to memory and doing autoimage
    >>> xyz = pt.get_coordinates(traj, autoimage=True)
    >>> xyz.shape
    (10, 5293, 3)

    >>> # load coordinates to memory for given mask
    >>> xyz = pt.get_coordinates(traj, mask='@CA')
    >>> xyz.shape
    (10, 12, 3)

    >>> # load coordinates of specific frame to memory and doing autoimage
    >>> xyz = pt.get_coordinates(traj, autoimage=True, frame_indices=[3, 6, 2, 5])
    >>> xyz.shape
    (4, 5293, 3)

    >>> # create frame iterator with some given cpptraj's commands
    >>> fi = pt.pipe(traj, ['autoimage', 'rms', 'center :1-6 origin'])
    >>> xyz = pt.get_coordinates(fi)
    >>> xyz.shape
    (10, 5293, 3)

    >>> # make your own out-of-core method
    >>> def my_method(traj):
    ...     for frame in traj:
    ...         frame.xyz += 2.
    ...         yield frame
    >>> fi = my_method(traj)
    >>> fi.__class__.__name__
    'generator'
    >>> xyz = pt.get_coordinates(fi)
    >>> xyz.shape
    (10, 5293, 3)
    '''
    has_any_iter_options = any(
        x is not None for x in (autoimage, rmsfit, mask, frame_indices))
    # try to iterate to get coordinates
    if isinstance(iterable, (Trajectory, TrajectoryIterator)):
        fi = iterable.iterframe(
            autoimage=autoimage,
            rmsfit=rmsfit,
            mask=mask,
            frame_indices=frame_indices)
    else:
        if has_any_iter_options:
            raise ValueError(
                'only support autoimage, rmsfit or mask for Trajectory and TrajectoryIterator'
            )
        fi = iterframe_master(iterable)
    if hasattr(fi, 'n_frames') and hasattr(fi, 'n_atoms'):
        # faster
        n_frames = fi.n_frames
        shape = (n_frames, fi.n_atoms, 3)
        arr = np.empty(shape, dtype='f8')
        for idx, frame in enumerate(fi):
            # real calculation
            arr[idx] = frame.xyz
        return arr
    else:
        # slower
        return np.array(
            [frame.xyz.copy() for frame in iterframe_master(iterable)],
            dtype='f8')


def load_batch(traj, txt):
    '''perform calculation for traj with cpptraj's batch style. This is for internal use.

    Parameters
    ----------
    traj : pytraj.TrajectoryIterator
    txt : text or a list of test
        cpptraj's commands

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> text = """
    ... autoimage
    ... radgyr @CA nomax
    ... molsurf !@H=
    ... """
    >>> state = pt.load_batch(traj, text)
    >>> state = state.run()
    >>> state.data
    <pytraj.datasets.CpptrajDatasetList - 3 datasets>

    >>> # raise if not TrajectoryIterator
    >>> traj2 = pt.Trajectory(xyz=traj.xyz, top=traj.top)
    >>> not isinstance(traj2, pt.TrajectoryIterator)
    True
    >>> pt.load_batch(traj2, text)
    Traceback (most recent call last):
        ...
    ValueError: only support TrajectoryIterator
    '''
    if not isinstance(traj, TrajectoryIterator):
        raise ValueError('only support TrajectoryIterator')
    return _load_batch(txt, traj=traj)


def read_data(filename, options=''):
    """same as readdata in cpptraj

    Returns
    -------
    out : CpptrajDatasetList
    """
    from pytraj.datasets import CpptrajDatasetList

    cdslist = CpptrajDatasetList()
    cdslist.read_data(filename, options)
    return cdslist


def write_data(filename, data):
    """same as writedata in cpptraj. (this is work in progress, only support 1D-array for now)

    Parameters
    ----------
    filename : str
        output filename. File format is deteced by its name.
    data : ndarray
    dtype : str, {'double', 'float', 'grid_float', 'grid_double'}, default 'double'

    Examples
    --------
    >>> # xmgrace format
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2()
    >>> rmsd_data = pt.rmsd(traj, ref=0)
    >>> pt.io.write_data('my.agr', rmsd_data)

    >>> # gnu format
    >>> pt.io.write_data('my.gnu', rmsd_data)
    """
    data = np.asarray(data, dtype='f8')

    dtype_map = {1: 'double', 3: 'grid_float'}

    from pytraj.core.c_core import CpptrajState, Command

    dtype = dtype_map[data.ndim]

    if data.ndim == 3:
        data = data.astype('f4')

    state = CpptrajState()
    dataset_name = 'my_data'
    d0 = state.data.add(dtype=dtype, name=dataset_name)
    d0.data = np.asarray(data)

    cm = 'writedata {filename} {dataset_name}'.format(
        filename=filename, dataset_name=dataset_name)

    with Command() as command:
        command.dispatch(state, cm)


def _format_convert(input_filename, output_filename):
    """

    Examples
    --------
    import pytraj as pt
    pt.io._format_convert('test.cpp4', 'test.dx')
    """
    from pytraj.core.c_core import CpptrajState, Command

    state = CpptrajState()

    with Command() as command:
        command.dispatch(state,
                         'readdata {} name mydata'.format(input_filename))
        command.dispatch(state, 'writedata {} mydata '.format(output_filename))


def _get_amberhome():
    amberhome = os.getenv('AMBERHOME', '')
    if not amberhome:
        raise EnvironmentError("must set AMBERHOME")
    return amberhome


def load_leap(command, verbose=False):
    """create pytraj.Trajectory from tleap's command.

    Notes
    -----
    If you load pdb file, make sure to use absolute dir.
    This method is not extensively tested. Use it with your own risk.
    """
    import subprocess

    command = command.strip()

    amberhome = _get_amberhome()
    tleap = amberhome + '/bin/tleap'

    lines = command.split('\n')

    if 'quit' not in lines[-1]:
        command = command + '\n' + 'quit'

    for line in lines:
        if line.lower().strip().startswith('saveamberparm'):
            parm, crd = line.split()[-2:]

    with tempfolder():
        leapin = 'tmp.leap'

        with open(leapin, 'w') as fh:
            fh.write(command)

        with open(os.devnull, 'wb') as devnull:
            if not verbose:
                subprocess.check_call(
                    [tleap, ' -f {}'.format(leapin)],
                    stdout=devnull,
                    stderr=subprocess.STDOUT)
            else:
                subprocess.check_call([tleap, ' -f {}'.format(leapin)])

        return load(crd, parm)


def select_atoms(mask, topology):
    '''return atom indices

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> atom_indices = pt.select_atoms('@CA', traj.top)
    >>> atom_indices
    array([  4,  15,  39, ..., 159, 173, 197])
    >>> pt.select_atoms(traj.top, '@CA')
    array([  4,  15,  39, ..., 159, 173, 197])
    '''
    if isinstance(mask, Topology) and isinstance(topology, str):
        mask, topology = topology, mask
    return topology.select(mask)
