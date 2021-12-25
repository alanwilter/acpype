import argparse

from acpype.params import MAXTIME, epilog, outTopols, usage


def get_option_parser():
    # not used yet: -e -l -r
    parser = argparse.ArgumentParser(usage=usage + epilog)
    group = parser.add_mutually_exclusive_group()
    parser.add_argument(
        "-i",
        "--input",
        action="store",
        dest="input",
        help="input file type like '.pdb', '.mdl', '.mol2' or SMILES string (mandatory if -p and -x not set)",
    )
    parser.add_argument(
        "-b",
        "--basename",
        action="store",
        dest="basename",
        help="a basename for the project (folder and output files)",
    )
    parser.add_argument(
        "-x",
        "--inpcrd",
        action="store",
        dest="inpcrd",
        help="amber inpcrd file name (always used with -p)",
    )
    parser.add_argument(
        "-p",
        "--prmtop",
        action="store",
        dest="prmtop",
        help="amber prmtop file name (always used with -x)",
    )
    parser.add_argument(
        "-c",
        "--charge_method",
        choices=["gas", "bcc", "user"],
        action="store",
        default="bcc",
        dest="charge_method",
        help="charge method: gas, bcc (default), user (user's charges in mol2 file)",
    )
    parser.add_argument(
        "-n",
        "--net_charge",
        action="store",
        type=int,
        default=None,
        dest="net_charge",
        help="net molecular charge (int), for gas default is 0",
    )
    parser.add_argument(
        "-m",
        "--multiplicity",
        action="store",
        type=int,
        default=1,
        dest="multiplicity",
        help="multiplicity (2S+1), default is 1",
    )
    parser.add_argument(
        "-a",
        "--atom_type",
        choices=["gaff", "amber", "gaff2", "amber2"],
        action="store",
        default="gaff2",
        dest="atom_type",
        help="atom type, can be gaff, gaff2 (default), amber (AMBER14SB) or amber2 (AMBER14SB + GAFF2)",
    )
    parser.add_argument(
        "-q",
        "--qprog",
        choices=["mopac", "sqm", "divcon"],
        action="store",
        default="sqm",
        dest="qprog",
        help="am1-bcc flag, sqm (default), divcon, mopac",
    )
    parser.add_argument(
        "-k",
        "--keyword",
        action="store",
        dest="keyword",
        help="mopac or sqm keyword, inside quotes",
    )
    parser.add_argument(
        "-f",
        "--force",
        action="store_true",
        dest="force",
        help="force topologies recalculation anew",
    )
    group.add_argument(
        "-d",
        "--debug",
        action="store_true",
        dest="debug",
        help="for debugging purposes, keep any temporary file created (not allowed with arg -w)",
    )
    group.add_argument(
        "-w",
        "--verboseless",
        action="store_false",
        default=True,
        dest="verboseless",
        help="print nothing (not allowed with arg -d)",
    )
    parser.add_argument(
        "-o",
        "--outtop",
        choices=["all"] + outTopols,
        action="store",
        default="all",
        dest="outtop",
        help="output topologies: all (default), gmx, cns or charmm",
    )
    parser.add_argument(
        "-z",
        "--gmx4",
        action="store_true",
        dest="gmx4",
        help="write RB dihedrals old GMX 4.0",
    )
    parser.add_argument(
        "-t",
        "--cnstop",
        action="store_true",
        dest="cnstop",
        help="write CNS topology with allhdg-like parameters (experimental)",
    )
    parser.add_argument(
        "-s",
        "--max_time",
        type=int,
        action="store",
        default=MAXTIME,
        dest="max_time",
        help="max time (in sec) tolerance for sqm/mopac, default is %i hours" % (MAXTIME // 3600),
    )
    parser.add_argument(
        "-y",
        "--ipython",
        action="store_true",
        dest="ipython",
        help="start iPython interpreter",
    )
    parser.add_argument(
        "-g",
        "--merge",
        action="store_true",
        dest="merge",
        help="Merge lower and uppercase atomtypes in GMX top file if identical parameters",
    )
    parser.add_argument(
        "-u",
        "--direct",
        action="store_true",
        dest="direct",
        help="for 'amb2gmx' mode, does a direct conversion, for any solvent (EXPERIMENTAL)",
        # NOTE: when solvent is present, gmx mdrun is not working, lack solvent topology
    )
    parser.add_argument(
        "-l",
        "--sorted",
        action="store_true",
        dest="sorted",
        help="sort atoms for GMX ordering",
    )
    parser.add_argument(
        "-j",
        "--chiral",
        action="store_true",
        dest="chiral",
        help="create improper dihedral parameters for chiral atoms in CNS",
    )
    parser.add_argument(
        "-v",
        "--version",
        action="store_true",
        dest="version",
        help="Show the Acpype version and exit",
    )
    return parser
