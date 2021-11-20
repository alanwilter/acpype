from __future__ import absolute_import
from ..trajectory.trajectory import Trajectory
from ..trajectory.frame import Frame
from ..analysis.c_action import c_action
from ..analysis.c_action.actionlist import ActionList
from ..datasets.c_datasetlist import DatasetList as CpptrajDatasetList


def make_structure(traj, command="", ref=None):
    """limited support for make_structure
    
    Parameters
    ----------
    traj : Trajectory
    command : cpptraj command or a list of cpptraj commands
    ref : None or a Frame or a Trajectory with one Frame
        if not None, use `ref` as reference. If ref.top is None, use `traj.top`, else use ref.top

        If `ref` is given, the `command` should be command = 'ref:<range>:<ref range>'.
        Please note that this syntax is different from cpptraj, you do not need to provide
        reference name.
        If ref is a Trajectory having more than 1 frame, only first frame is considered.
      
    Returns
    -------
    traj : itself
  
    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2()[:]

    >>> # Make alpha helix for residue 1 to 12)
    >>> traj = pt.make_structure(traj, "alpha:1-12")

    >>> # Make hairpin for residue 1 to 5 and make alpha helix for residue 6 to 12
    >>> traj = pt.make_structure(traj, ["hairpin:1-5", "alpha:6-12"])

    >>> # From cpptraj example:
    >>> # Make residues 1 and 12 'extended', residues 6 and 7 a type I' turn, and two
    >>> # custom assignments, one (custom1) for residues 2-5, the other (custom2) for residues 8-11:
    >>> traj = pt.make_structure(traj, ["extended:1,12",
    ...                                 "custom1:2-5:-80.:130.:-130.:140.",
    ...                                 "typeI':6-7",
    ...                                 "custom2:8-11:-140.:170.:-100.:140."])
    >>>

    >>> # Make new structure from reference
    >>> def make_new_by_using_ref(): # doctest: +SKIP
    ...     tz2_parm7 = 'tz2.parm7'
    ...     ref_rst7 = 'tz2.rst7'
    ...     trajin_rst7 = pp2.rst7.save'
    ...     traj = pt.load(trajin_rst7, top=tz2_parm7)
    ...     ref = pt.load(ref_rst7, top=tz2_parm7)
    ...     # make dihedral angle for targeted residues 1 to 13 by using
    ...     # reference residues 1 to 13
    ...     pt.make_structure(traj, "ref:1-13:1-13", ref=ref)


    Notes
    -----
    Apply dihedrals to specified residues using arguments found in <List of Args>,
    where an argument is 1 or more of the following arg types:

    '<sstype>:<res range>'

        Apply SS type (phi/psi) to residue range.

        <sstype> standard = alpha, left, pp2, hairpin, extended

        <sstype> turn = typeI, typeII, typeVIII, typeI', typeII, typeVIa1, typeVIa2, typeVIb
            Turns are applied to 2 residues at a time, so resrange must be divisible by 4.

    '<custom ss>:<res range>:<phi>:<psi>'

        Apply custom <phi>/<psi> to residue range.

    '<custom turn>:<res range>:<phi1>:<psi1>:<phi2>:<psi2>'

        Apply custom turn <phi>/<psi> pair to residue range.

    '<custom dih>:<res range>:<dih type>:<angle>'

        Apply <angle> to dihedrals in range.

        <dih type> = alpha beta gamma delta epsilon zeta nu1 nu2 h1p c2p chin phi psi chip omega

    '<custom dih>:<res range>:<at0>:<at1>:<at2>:<at3>:<angle>[:<offset>]'

          Apply <angle> to dihedral defined by atoms <at1>, <at2>, <at3>, and <at4>.

          Offset -2=<a0><a1> in previous res, -1=<a0> in previous res,
                  0=All <aX> in single res,
                  1=<a3> in next res, 2=<a2><a3> in next res.
    """
    assert isinstance(command, list) or isinstance(
        command, str), 'command must be a string or a list of string'
    assert isinstance(traj, Trajectory), 'traj must be a Trajectory object'

    cmlist = [
        command,
    ] if isinstance(command, str) else command

    if ref is not None:
        ref_name = 'myref'
        for index, cm in enumerate(cmlist):
            if cm.startswith('ref'):
                tokens = cm.split(':')
                assert len(
                    tokens
                ) == 3, 'command must follow format ref:<range>:<ref range>'
                # insert ref_name
                cm_ = ':'.join(tokens[:2] + [ref_name] + tokens[-1:])
                cmlist[index] = cm_
        c_dslist = CpptrajDatasetList()
        frame_dset = c_dslist.add('reference', name=ref_name)
        frame_dset.top = ref.top if ref.top is not None else traj.top
        ref_ = ref if isinstance(ref, Frame) else ref[0]
        frame_dset.append(ref_)
    else:
        c_dslist = None
    actlist = ActionList()

    for cm in cmlist:
        if ref is not None:
            actlist.add(
                c_action.Action_MakeStructure(),
                cm,
                top=traj.top,
                dslist=c_dslist)
        else:
            actlist.add(c_action.Action_MakeStructure(), cm, top=traj.top)

    for frame in traj:
        actlist.compute(frame)
    return traj
