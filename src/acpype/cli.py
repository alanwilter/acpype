#!/usr/bin/env python3
"""Command line interface for ACPYPE, built on typer."""

import os
import sys
import time
from enum import StrEnum
from shutil import rmtree
from typing import Annotated

import typer

from acpype.logger import copy_log, tmpLogFile
from acpype.logger import set_logging_conf as logger
from acpype.params import MAXTIME, binaries, epilog
from acpype.topol import AbstractTopol, ACTopol, MolTopol, header
from acpype.utils import elapsedTime, while_replace


class ChargeMethod(StrEnum):
    """Charge methods antechamber can apply."""

    gas = "gas"
    bcc = "bcc"
    abcg2 = "abcg2"
    user = "user"


class AtomType(StrEnum):
    """Supported atom type sets."""

    gaff = "gaff"
    amber = "amber"
    gaff2 = "gaff2"
    amber2 = "amber2"


class QProg(StrEnum):
    """Quantum programs usable for am1-bcc charges.

    Only ``sqm`` is usable out of the box. ``mopac`` and ``divcon`` are kept for
    compatibility but need a full external AmberTools: antechamber drives them
    through its own ``mopac.sh``/``divcon`` wrapper scripts, which are absent from
    the AmberTools bundled here and from conda-forge's ``ambertools`` package.
    """

    mopac = "mopac"
    sqm = "sqm"
    divcon = "divcon"


class OutTopol(StrEnum):
    """Topology flavours ACPYPE can emit."""

    all = "all"
    gmx = "gmx"
    cns = "cns"
    charmm = "charmm"


app = typer.Typer(
    add_completion=False,
    context_settings={"help_option_names": ["-h", "--help"]},
)


def _chk_py_ver():
    # `requires-python` already blocks installs below the floor, but acpype is often
    # run straight from a source checkout, where no such check applies.
    if sys.version_info < (3, 12):  # noqa: UP036
        msg = "Sorry, you need python 3.12 or higher"
        logger().error(msg)
        raise Exception(msg)


def _handle_exception(level):
    _exceptionType, exceptionValue, _exceptionTraceback = sys.exc_info()
    logger(level).exception(f"ACPYPE FAILED: {exceptionValue}")
    return True


@app.command(epilog=epilog.strip())
def main(
    ctx: typer.Context,
    input: Annotated[
        str | None,
        typer.Option("-i", "--input", help="input file ('.pdb', '.mdl', '.mol2') or SMILES string"),
    ] = None,
    basename: Annotated[
        str | None, typer.Option("-b", "--basename", help="a basename for the project (folder and output files)")
    ] = None,
    inpcrd: Annotated[
        str | None, typer.Option("-x", "--inpcrd", help="amber inpcrd file name (always used with -p)")
    ] = None,
    prmtop: Annotated[
        str | None, typer.Option("-p", "--prmtop", help="amber prmtop file name (always used with -x)")
    ] = None,
    charge_method: Annotated[
        ChargeMethod,
        typer.Option(
            "-c",
            "--charge_method",
            help="charge method: abcg2 is recommended for GAFF2, user reads charges from the mol2",
        ),
    ] = ChargeMethod.bcc,
    net_charge: Annotated[
        int | None, typer.Option("-n", "--net_charge", help="net molecular charge; guessed when not given")
    ] = None,
    multiplicity: Annotated[int, typer.Option("-m", "--multiplicity", help="multiplicity (2S+1)")] = 1,
    predindex: Annotated[
        int,
        typer.Option(
            "-r",
            "--predindex",
            min=0,
            max=5,
            help=(
                "antechamber atom/bond type prediction index (its -j flag): 0 none, "
                "1 atom type, 2 full bond types, 3 part bond types, 4 atom and full "
                "bond type, 5 atom and part bond type. ACPYPE needs atom types, so "
                "only 1, 4 and 5 can produce a topology"
            ),
        ),
    ] = 4,
    atom_type: Annotated[
        AtomType,
        typer.Option("-a", "--atom_type", help="atom type; amber is AMBER14SB, amber2 is AMBER14SB + GAFF2"),
    ] = AtomType.gaff2,
    qprog: Annotated[
        QProg,
        typer.Option(
            "-q",
            "--qprog",
            help="am1-bcc engine; only sqm is bundled, mopac and divcon need an external AmberTools",
        ),
    ] = QProg.sqm,
    keyword: Annotated[str | None, typer.Option("-k", "--keyword", help="mopac or sqm keyword, inside quotes")] = None,
    force: Annotated[bool, typer.Option("-f", "--force", help="force topologies recalculation anew")] = False,
    debug: Annotated[
        bool,
        typer.Option("-d", "--debug", help="keep any temporary file created (not allowed with -w)"),
    ] = False,
    verboseless: Annotated[bool, typer.Option("-w", "--verboseless", help="print nothing (not allowed with -d)")] = (
        False
    ),
    outtop: Annotated[OutTopol, typer.Option("-o", "--outtop", help="output topologies")] = OutTopol.all,
    gmx4: Annotated[bool, typer.Option("-z", "--gmx4", help="write RB dihedrals for old GMX 4.0")] = False,
    cnstop: Annotated[
        bool, typer.Option("-t", "--cnstop", help="write CNS topology with allhdg-like parameters (experimental)")
    ] = False,
    max_time: Annotated[
        int, typer.Option("-s", "--max_time", help="max time (in sec) tolerance for sqm/mopac")
    ] = MAXTIME,
    ipython: Annotated[bool, typer.Option("-y", "--ipython", help="start iPython interpreter")] = False,
    merge: Annotated[
        bool,
        typer.Option("-g", "--merge", help="merge lower and uppercase atomtypes in GMX top if parameters match"),
    ] = False,
    split_molecules: Annotated[
        bool,
        typer.Option(
            "-S",
            "--split_molecules",
            help="in 'amb2gmx' mode, write one [ moleculetype ] per molecule AMBER found, e.g. protein and ligand apart",
        ),
    ] = False,
    direct: Annotated[
        bool,
        typer.Option("-u", "--direct", help="in 'amb2gmx' mode, do a direct conversion for any solvent (EXPERIMENTAL)"),
    ] = False,
    is_sorted: Annotated[bool, typer.Option("-l", "--sorted", help="sort atoms for GMX ordering")] = False,
    chiral: Annotated[
        bool, typer.Option("-j", "--chiral", help="create improper dihedral parameters for chiral atoms in CNS")
    ] = False,
    version: Annotated[bool, typer.Option("-v", "--version", help="show the ACPYPE version and exit")] = False,
) -> None:
    """Generate topologies for chemical compounds, using Antechamber.

    Topologies can be written for GROMACS, CNS/XPLOR, CHARMM and AMBER. Given an AMBER
    prmtop and inpcrd pair instead ('-p' and '-x'), ACPYPE runs in 'amb2gmx' mode.
    """
    state = ctx.obj if isinstance(ctx.obj, dict) else {}
    ac_binaries = state.get("binaries", binaries)

    at0 = time.time()
    amb2gmxF = False

    if version:
        print(header)
        sys.exit(0)

    # argparse enforced this through a mutually exclusive group; typer has no
    # equivalent, so it is checked by hand.
    if debug and verboseless:
        raise typer.BadParameter("'-d' and '-w' are mutually exclusive")

    # `verboseless` inverts: the historic flag switched verbosity off.
    verbose = not verboseless

    level = 20
    if verboseless:
        level = 100
    if debug:
        level = 10

    logger(level).info(header)

    if not input:
        amb2gmxF = True
        if not inpcrd or not prmtop:
            raise typer.BadParameter("missing input files")
    elif inpcrd or prmtop:
        raise typer.BadParameter("either '-i' or ('-p', '-x'), but not both")

    logger(level).debug(f"CLI: {' '.join(state.get('argv', sys.argv[1:]))}")
    texta = f"Python Version {sys.version}"
    logger(level).debug(while_replace(texta))

    if direct and not amb2gmxF:
        raise typer.BadParameter("option -u is only meaningful in 'amb2gmx' mode (args '-p' and '-x')")

    acpypeFailed = False
    if amb2gmxF:
        logger(level).info("Converting Amber input files to Gromacs ...")
        try:
            molecule: AbstractTopol = MolTopol(
                acFileXyz=inpcrd or "",
                acFileTop=prmtop or "",
                amb2gmx=True,
                debug=debug,
                basename=basename,
                verbose=verbose,
                gmx4=gmx4,
                merge=merge,
                direct=direct,
                split_molecules=split_molecules,
                is_sorted=is_sorted,
                chiral=chiral,
            )
        except Exception:
            acpypeFailed = _handle_exception(level)
        if not acpypeFailed:
            try:
                molecule.writeGromacsTopolFiles()
                molecule.printDebug("prmtop and inpcrd files parsed")
            except Exception:
                acpypeFailed = _handle_exception(level)

    else:
        try:
            molecule = ACTopol(
                input,
                binaries=ac_binaries,
                chargeType=charge_method.value,
                chargeVal=net_charge,
                debug=debug,
                multiplicity=multiplicity,
                predIndex=predindex,
                atomType=atom_type.value,
                force=force,
                outTopol=outtop.value,
                allhdg=cnstop,
                basename=basename,
                timeTol=max_time,
                qprog=qprog.value,
                ekFlag=f'''"{keyword}"''',
                verbose=verbose,
                gmx4=gmx4,
                merge=merge,
                direct=direct,
                split_molecules=split_molecules,
                is_sorted=is_sorted,
                chiral=chiral,
                amb2gmx=False,
            )
        except Exception:
            acpypeFailed = _handle_exception(level)
        if not acpypeFailed:
            try:
                molecule.createACTopol()
            except Exception:
                acpypeFailed = _handle_exception(level)
        if not acpypeFailed:
            try:
                molecule.createMolTopol()
            except Exception:
                acpypeFailed = _handle_exception(level)

    execTime = round(time.time() - at0)
    amsg = "less than a second" if execTime == 0 else elapsedTime(execTime)
    logger(level).info(f"Total time of execution: {amsg}")

    if ipython:
        try:
            import IPython

            IPython.embed(colors="neutral")
        except ModuleNotFoundError:
            logger(level).exception("No 'ipython' installed")

    if not debug:
        try:
            rmtree(molecule.tmpDir)
        except Exception:
            logger(level).debug("No tmp folder left to be removed")
    else:
        try:
            if molecule.tmpDir:
                logger(level).debug(f"Keeping folder '{molecule.tmpDir}' for possible helping debugging")
        except Exception:
            logger(level).debug("No tmp folder left to be removed")

    try:
        copy_log(molecule)
    except UnboundLocalError:
        print(f"Log tmp location: {tmpLogFile}")

    if acpypeFailed:
        sys.exit(19)

    os.chdir(molecule.rootDir)

    if not amb2gmxF and molecule.obabelExe and molecule.checkSmiles():
        afile = "smiles_molecule.mol2"
        if os.path.exists(afile):
            os.remove(afile)

    # Tells init_main this was a complete run, so the exit status typer raises on
    # success can be swallowed while `--version` keeps exiting explicitly.
    state["completed"] = True


def init_main(binaries: dict[str, str] = binaries, argv: list[str] | None = None):
    """Orchestrate the command line usage for ACPYPE with its all input arguments.

    Args:
        binaries (dict[str, str], optional): Mostly used for debug and testing. Defaults to ``acpype.params.binaries``.
        argv (list[str] | None, optional): Mostly used for debug and testing. Defaults to None.

    Returns:
        SystemExit(status): 0 or 19 (failed)
    """
    _chk_py_ver()
    if argv is None:
        argv = sys.argv[1:]

    state: dict = {"binaries": binaries, "argv": list(argv), "completed": False}
    command = typer.main.get_command(app)
    try:
        command(args=list(argv), prog_name="acpype", obj=state, standalone_mode=True)
    except SystemExit as exc:
        # typer exits 0 after any successful invocation, but callers of init_main
        # expect a plain return from a normal run; `--version` still exits.
        if exc.code == 0 and state["completed"]:
            return
        raise


if __name__ == "__main__":
    init_main()  # necessary for to call in anaconda package;
