#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
    Requirements: Python 3.6 or higher
                  Antechamber (from AmberTools preferably)
                  OpenBabel (optional, but strongly recommended)

    This code is released under GNU General Public License V3.

          <<<  NO WARRANTY AT ALL!!!  >>>

    It was inspired by:

    - amb2gmx.pl (Eric Sorin, David Mobley and John Chodera)
      and depends on Antechamber and Openbabel

    - YASARA Autosmiles:
      http://www.yasara.org/autosmiles.htm (Elmar Krieger)

    - topolbuild (Bruce Ray)

    - xplo2d (G.J. Kleywegt)

    For Non-uniform 1-4 scale factor conversion (e.g. if using GLYCAM06), please cite:

    BERNARDI, A., FALLER, R., REITH, D., and KIRSCHNER, K. N. ACPYPE update for
    nonuniform 1-4 scale factors: Conversion of the GLYCAM06 force field from AMBER
    to GROMACS. SoftwareX 10 (2019), 100241. doi: 10.1016/j.softx.2019.100241

    For Antechamber, please cite:
    1. WANG, J., WANG, W., KOLLMAN, P. A., and CASE, D. A. Automatic atom type and
       bond type perception in molecular mechanical calculations. Journal of Molecular
       Graphics and Modelling 25, 2 (2006), 247-260. doi: 10.1016/j.jmgm.2005.12.005
    2. WANG, J., WOLF, R. M., CALDWELL, J. W., KOLLMAN, P. A., and CASE, D. A.
       Development and testing of a General Amber Force Field. Journal of Computational
       Chemistry 25, 9 (2004), 1157-1174. doi: 10.1002/jcc.20035

    If you use this code, I am glad if you cite:

    SOUSA DA SILVA, A. W. & VRANKEN, W. F.
    ACPYPE - AnteChamber PYthon Parser interfacE.
    BMC Research Notes 5 (2012), 367 doi: 10.1186/1756-0500-5-367
    http://www.biomedcentral.com/1756-0500/5/367

    BATISTA, P. R.; WILTER, A.; DURHAM, E. H. A. B. & PASCUTTI, P. G. Molecular
    Dynamics Simulations Applied to the Study of Subtypes of HIV-1 Protease.
    Cell Biochemistry and Biophysics 44 (2006), 395-404. doi: 10.1385/CBB:44:3:395

    Alan Silva, D.Sc.
    alanwilter _at_ gmail _dot_ com
"""

import os
import sys
import time
import traceback
from shutil import rmtree

from acpype.logger import copy_log
from acpype.logger import set_logging_conf as logger
from acpype.params import binaries
from acpype.parser_args import get_option_parser
from acpype.topol import ACTopol, MolTopol, header
from acpype.utils import elapsedTime, set_for_pip, while_replace


def chk_py_ver():
    if sys.version_info < (3, 6):
        msg = "Sorry, you need python 3.6 or higher"
        logger().exception(msg)
        raise Exception(msg)


def init_main(binaries=binaries, argv=None):

    """
    Main function, to satisfy Conda
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
    texta = "Python Version %s" % sys.version
    logger(level).debug(while_replace(texta))

    if args.direct and not amb2gmxF:
        parser.error("option -u is only meaningful in 'amb2gmx' mode (args '-p' and '-x')")

    try:
        if amb2gmxF:
            logger(level).info("Converting Amber input files to Gromacs ...")
            system = MolTopol(
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

            system.printDebug("prmtop and inpcrd files parsed")
            system.writeGromacsTopolFiles()
            copy_log(system)
        else:
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
                ekFlag='''"%s"''' % args.keyword,
                verbose=args.verboseless,
                gmx4=args.gmx4,
                merge=args.merge,
                direct=args.direct,
                is_sorted=args.sorted,
                chiral=args.chiral,
                amb2gmx=False,
            )

            molecule.createACTopol()
            molecule.createMolTopol()
            copy_log(molecule)
        acpypeFailed = False
    except Exception:
        _exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
        logger(level).error("ACPYPE FAILED: %s" % exceptionValue)
        if args.debug:
            traceback.print_tb(exceptionTraceback, file=sys.stdout)
        acpypeFailed = True

    execTime = int(round(time.time() - at0))
    if execTime == 0:
        amsg = "less than a second"
    else:
        amsg = elapsedTime(execTime)
    logger(level).info("Total time of execution: %s" % amsg)

    if args.ipython:
        import IPython  # pylint: disable=import-outside-toplevel

        IPython.embed(colors="neutral")

    try:
        rmtree(molecule.tmpDir)
    except Exception:
        pass

    if acpypeFailed:
        sys.exit(19)
    try:
        os.chdir(molecule.rootDir)
    except Exception:
        pass

    if not amb2gmxF and molecule.obabelExe:
        if molecule.checkSmiles():
            afile = "smiles_molecule.mol2"
            if os.path.exists(afile):
                os.remove(afile)


if __name__ == "__main__":
    init_main()  # necessary for to call in anaconda package;
