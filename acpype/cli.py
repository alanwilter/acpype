#!/usr/bin/env python3

import os
import sys
import time
from shutil import rmtree
from typing import Dict, List, Optional

from acpype.logger import copy_log
from acpype.logger import set_logging_conf as logger
from acpype.logger import tmpLogFile
from acpype.params import binaries
from acpype.parser_args import get_option_parser
from acpype.topol import AbstractTopol, ACTopol, MolTopol, header
from acpype.utils import elapsedTime, set_for_pip, while_replace


def chk_py_ver():
    if sys.version_info < (3, 6):
        msg = "Sorry, you need python 3.6 or higher"
        logger().error(msg)
        raise Exception(msg)


def handle_exception(level):
    _exceptionType, exceptionValue, _exceptionTraceback = sys.exc_info()
    logger(level).exception(f"ACPYPE FAILED: {exceptionValue}")
    return True


def init_main(binaries: Dict[str, str] = binaries, argv: Optional[List[str]] = None):

    """
    Orchestrate the command line usage for ACPYPE with its all input arguments.

    Args:
        binaries (Dict[str, str], optional): Mostly used for debug and testing. Defaults to ``acpype.params.binaries``.
        argv (Optional[List[str]], optional): Mostly used for debug and testing. Defaults to None.

    Returns:
        SystemExit(status): 0 or 19 (failed)
    """
    chk_py_ver()
    set_for_pip(binaries)
    if argv is None:
        argv = sys.argv[1:]

    parser = get_option_parser()
    args = parser.parse_args(argv)

    at0 = time.time()

    amb2gmxF = False

    if args.version:
        print(header)
        sys.exit(0)

    level = 20
    if not args.verboseless:
        level = 100
    if args.debug:
        level = 10

    logger(level).info(header)

    if not args.input:
        amb2gmxF = True
        if not args.inpcrd or not args.prmtop:
            parser.error("missing input files")
    elif args.inpcrd or args.prmtop:
        parser.error("either '-i' or ('-p', '-x'), but not both")

    logger(level).debug(f"CLI: {' '.join(argv)}")
    texta = f"Python Version {sys.version}"
    logger(level).debug(while_replace(texta))

    if args.direct and not amb2gmxF:
        parser.error("option -u is only meaningful in 'amb2gmx' mode (args '-p' and '-x')")

    acpypeFailed = False
    if amb2gmxF:
        logger(level).info("Converting Amber input files to Gromacs ...")
        try:
            molecule: AbstractTopol = MolTopol(
                acFileXyz=args.inpcrd,
                acFileTop=args.prmtop,
                amb2gmx=True,
                debug=args.debug,
                basename=args.basename,
                verbose=args.verboseless,
                gmx4=args.gmx4,
                merge=args.merge,
                direct=args.direct,
                is_sorted=args.sorted,
                chiral=args.chiral,
            )
        except Exception:
            acpypeFailed = handle_exception(level)
        if not acpypeFailed:
            try:
                molecule.writeGromacsTopolFiles()
                molecule.printDebug("prmtop and inpcrd files parsed")
            except Exception:
                acpypeFailed = handle_exception(level)

    else:
        try:
            molecule = ACTopol(
                args.input,
                binaries=binaries,
                chargeType=args.charge_method,
                chargeVal=args.net_charge,
                debug=args.debug,
                multiplicity=args.multiplicity,
                atomType=args.atom_type,
                force=args.force,
                outTopol=args.outtop,
                allhdg=args.cnstop,
                basename=args.basename,
                timeTol=args.max_time,
                qprog=args.qprog,
                ekFlag=f'''"{args.keyword}"''',
                verbose=args.verboseless,
                gmx4=args.gmx4,
                merge=args.merge,
                direct=args.direct,
                is_sorted=args.sorted,
                chiral=args.chiral,
                amb2gmx=False,
            )
        except Exception:
            acpypeFailed = handle_exception(level)
        if not acpypeFailed:
            try:
                molecule.createACTopol()
            except Exception:
                acpypeFailed = handle_exception(level)
        if not acpypeFailed:
            try:
                molecule.createMolTopol()
            except Exception:
                acpypeFailed = handle_exception(level)

    execTime = int(round(time.time() - at0))
    if execTime == 0:
        amsg = "less than a second"
    else:
        amsg = elapsedTime(execTime)
    logger(level).info(f"Total time of execution: {amsg}")

    if args.ipython:
        try:
            import IPython

            IPython.embed(colors="neutral")
        except ModuleNotFoundError:
            logger(level).exception("No 'ipython' installed")

    if not args.debug:
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

    if not amb2gmxF and molecule.obabelExe:
        if molecule.checkSmiles():
            afile = "smiles_molecule.mol2"
            if os.path.exists(afile):
                os.remove(afile)


if __name__ == "__main__":
    init_main()  # necessary for to call in anaconda package;
