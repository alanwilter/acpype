import abc
import array
import collections
import math
import os
import pickle
import re
import signal
import subprocess as sub
from datetime import datetime
from shutil import copy2, rmtree, which

from acpype import __version__ as version
from acpype.logger import set_logging_conf as logger
from acpype.mol import Angle, Atom, AtomType, Bond, Dihedral
from acpype.params import (
    FF14SB_RETYPES,
    MAXTIME,
    PROTEIN_FF,
    PROTEIN_FF_LIBS,
    TLEAP_TEMPLATE,
    binaries,
    cal,
    dictAtomTypeAmb2OplsGmxCode,
    dictAtomTypeGaff2OplsGmxCode,
    diffTol,
    ionOrSolResNameList,
    maxDist,
    maxDist2,
    minDist,
    minDist2,
    oplsCode2AtomTypeDict,
    outTopols,
    qConv,
    qDict,
    radPi,
    specialGaffAtoms,
)
from acpype.utils import (
    _getoutput,
    bundled_amber_dir,
    distanceAA,
    elapsedTime,
    find_bin,
    imprDihAngle,
    job_pids_family,
    mergeFrcmodGaps,
    parmMerge,
    readAmberLibTypes,
    readParmAtomTypes,
    retypeMol2Atoms,
    set_for_pip,
    unknownMol2Types,
    while_replace,
)

year = datetime.today().year
tag = version

lineHeader = f"""
| ACPYPE: AnteChamber PYthon Parser interfacE v. {tag} (c) {year} AWSdS |
"""
frameLine = (len(lineHeader) - 2) * "="
header = f"{frameLine}{lineHeader}{frameLine}"

#    TODO:
#        Howto Charmm and Amber with NAMD
#        Howto build topology for a modified amino acid
#        CYANA topology files

head = "%s created by acpype (v: " + tag + ") on %s\n"

date = datetime.now().ctime()

pid: int


class Topology_14:
    """Amber topology abstraction for non-uniform 1-4 scale factors."""

    def __init__(self) -> None:
        self.pointers = array.array("d")
        self.charge = array.array("d")
        self.atom_type_index = array.array("d")
        self.nonbonded_parm_index = array.array("d")
        self.scee_scale_factor = array.array("d")
        self.scnb_scale_factor = array.array("d")
        self.dihedral_force_constants = array.array("d")
        self.dihedral_periodicity = array.array("d")
        self.dihedral_phase = array.array("d")
        self.dihedral_yes_H = array.array("d")
        self.dihedral_no_H = array.array("d")
        self.lennard_jones_acoef = array.array("d")
        self.lennard_jones_bcoef = array.array("d")

    def read_amber_topology(self, buff):
        """Read AMBER topology file."""
        flag_strings = [
            "%FLAG POINTERS",
            "%FLAG CHARGE",
            "%FLAG ATOM_TYPE_INDEX",
            "%FLAG NONBONDED_PARM_INDEX",
            "%FLAG SCEE_SCALE_FACTOR",
            "%FLAG SCNB_SCALE_FACTOR",
            "%FLAG DIHEDRAL_FORCE_CONSTANT",
            "%FLAG DIHEDRAL_PERIODICITY",
            "%FLAG DIHEDRAL_PHASE",
            "%FLAG DIHEDRALS_INC_HYDROGEN",
            "%FLAG DIHEDRALS_WITHOUT_HYDROGEN",
            "%FLAG LENNARD_JONES_ACOEF",
            "%FLAG LENNARD_JONES_BCOEF",
        ]
        attributes = [
            "pointers",
            "charge",
            "atom_type_index",
            "nonbonded_parm_index",
            "scee_scale_factor",
            "scnb_scale_factor",
            "dihedral_force_constants",
            "dihedral_periodicity",
            "dihedral_phase",
            "dihedral_yes_H",
            "dihedral_no_H",
            "lennard_jones_acoef",
            "lennard_jones_bcoef",
        ]
        for i, _item in enumerate(attributes):
            try:
                setattr(self, attributes[i], self.p7_array_read(buff, flag_strings[i]))
            except Exception:
                logger().exception(f"Skipping non-existent attributes {attributes[i]} {flag_strings[i]}")

    @staticmethod
    def skipline(buff, index) -> int:
        """Skip line."""
        while buff[index] != "\n":
            index += 1
        index += 1
        return index

    def p7_array_read(self, buff, flag_string):
        """Convert AMBER topology data to python array."""
        myarray = array.array("d")
        i = buff.index(flag_string)
        i = self.skipline(buff, i)
        i = self.skipline(buff, i)
        while 1:
            while buff[i] == " " or buff[i] == "\t" or buff[i] == "\n":
                i += 1
            j = i
            if buff[i] == "%":
                break
            while buff[i] != " " and buff[i] != "\t" and buff[i] != "\n":
                i += 1
            myarray.append(float(buff[j:i]))
        return myarray

    def print_gmx_pairs(self):
        """Generate non-bonded pairs list."""
        pair_list = []
        pair_buff = "[ pairs_nb ]\n;       ai         aj  funct         qi         qj           sigma         epsilon\n"
        pair_list.append(pair_buff)
        dihedrals = self.dihedral_yes_H + self.dihedral_no_H
        dih_number = len(dihedrals)
        j = 0
        while j < dih_number:
            if dihedrals[j + 2] > 0:
                parm_idx = int(dihedrals[j + 4]) - 1
                scee_scale_factor = self.scee_scale_factor[parm_idx]
                if scee_scale_factor == 0:
                    scee_scale_factor = 1.2
                ai = int(abs(dihedrals[j]) / 3)
                al = int(abs(dihedrals[j + 3]) / 3)
                qi = self.charge[ai] / qConv
                ql = self.charge[al] / qConv / scee_scale_factor
                ntypes = int(self.pointers[1])
                ai_index = int(self.atom_type_index[ai])
                al_index = int(self.atom_type_index[al])
                nb_parm_index = int(self.nonbonded_parm_index[ntypes * (ai_index - 1) + al_index - 1]) - 1
                scnb_scale_factor = self.scnb_scale_factor[parm_idx]
                if scnb_scale_factor == 0:
                    scnb_scale_factor = 2
                lj_acoeff = self.lennard_jones_acoef[nb_parm_index] / scnb_scale_factor
                lj_bcoeff = self.lennard_jones_bcoef[nb_parm_index] / scnb_scale_factor
                if lj_bcoeff != 0:
                    sigma6 = lj_acoeff / lj_bcoeff
                else:
                    sigma6 = 1  # arbitrary and doesn't matter
                epsilon = lj_bcoeff / 4 / sigma6 * 4.184
                sigma = sigma6 ** (1 / 6) / 10
                pair_buff = (
                    f"{ai + 1:>10.0f} {al + 1:>10.0f} {1:>6.0f} "
                    + f"{qi:>10.6f} {ql:>10.6f} "
                    + f"{sigma:>15.5e} {epsilon:>15.5e}\n"
                )
                pair_list.append(pair_buff)
            j += 5
        return "".join(pair_list)

    def hasNondefault14(self):
        """Check non-uniform 1-4 scale factor."""
        if any(val not in (0, 1.2) for val in self.scee_scale_factor):
            return True
        return any(val not in (0, 2) for val in self.scnb_scale_factor)

    def patch_gmx_topol14(self, gmx_init_top):
        """Patch GMX topology file for non-uniform 1-4 scale factor."""
        pair_buff = self.print_gmx_pairs()
        jdefault = gmx_init_top.index("\n[ atomtypes ]")
        ipair = gmx_init_top.index("[ pairs ]")
        jpair = gmx_init_top.index("\n[ angles ]")
        init_buff = (
            "\n\n[ defaults ]\n"
            + "; nbfunc        comb-rule       gen-pairs       \n"
            + "1               2               no              \n"
        )
        return (
            gmx_init_top.splitlines()[0]
            + init_buff
            + gmx_init_top[jdefault:ipair]
            + pair_buff
            + gmx_init_top[jpair : len(gmx_init_top)]
        )


class AbstractTopol(abc.ABC):
    """Abstract super class to build topologies"""

    # Assigned by the concrete subclasses rather than here, but used by methods on
    # this class, so they are declared for the benefit of type checkers only.
    level: int
    binaries: dict[str, str]
    _parent: "ACTopol"

    @abc.abstractmethod
    def __init__(self):
        self.debug: bool = False
        self.verbose: bool = False
        self.chargeVal: str | None = None
        self.tmpDir: str = ""
        self.absInputFile: str = ""
        self.chargeType: str = ""
        self.obabelExe: str = ""
        self.baseName: str = ""
        self.acExe: str = ""
        self.force: bool = False
        self.acBaseName: str = ""
        self.atomType: str = ""
        self.acMol2FileName: str = ""
        self.multiplicity: str = ""
        self.predIndex: int = 4
        self.qFlag: int = 0
        self.ekFlag: str = ""
        self.timeTol: int = 0
        self.acXyzFileName: str = ""
        self.acTopFileName: str = ""
        self.acParDict: dict[str, str] = {}
        self.tleapExe: str = ""
        self.parmchkExe: str = ""
        self.acFrcmodFileName: str = ""
        self.gaffDatfile: str = ""
        self.homeDir: str = ""
        self.rootDir: str = ""
        self.extOld: str = ""
        self.direct: bool = False
        self.splitMolecules: bool = False
        self.proteinFF: str = "ff14SB"
        self.gmxSplitNames: list = []
        self.merge: bool = False
        self.gmx4: bool = False
        self.sorted: bool = False
        self.chiral: bool = False
        self.outTopols: list = []
        self.ext: str = ""
        self.xyzFileData: list = []
        self.charmmBase: str = ""
        self.allhdg: bool = False
        self.topo14Data: Topology_14 = Topology_14()
        self.atomPairs: set | list = []
        self.properDihedralsGmx45: list = []
        self.properDihedralsAlphaGamma: list = []
        self.properDihedralsCoefRB: list = []
        self.resName: str = ""
        self.acLog: str = ""
        self.tleapLog: str = ""
        self.parmchkLog: str = ""
        self.inputFile: str = ""
        self.obabelLog: str = ""
        self.absHomeDir: str = ""
        self.molTopol: MolTopol | None = None
        self.topFileData: list = []
        self.residueLabel: list = []
        self._atomTypeNameList: list = []
        self.atomTypeSystem: str = ""
        self.totalCharge: int = 0
        self.atoms: list = []
        self.atomTypes: list = []
        self.pbc: list = []
        self.bonds: list = []
        self.angles: list = []
        self.properDihedrals: list = []
        self.improperDihedrals: list = []
        self.condensedProperDihedrals: list = []
        self.chiralGroups: list = []
        self.excludedAtoms: list = []
        self.atomsGromacs: list = []
        self.atomTypesGromacs: list = []
        self.CnsTopFileName: str = ""
        self.CnsInpFileName: str = ""
        self.CnsParFileName: str = ""
        self.CnsPdbFileName: str = ""
        self.is_smiles: bool = False
        self.smiles: str = ""
        self.amb2gmx: bool = False

    set_for_pip(binaries)  # check and call environments

    def printDebug(self, text=""):
        """Debug log level."""
        logger(self.level).debug(f"{while_replace(text)}")

    def printWarn(self, text=""):
        """Warn log level."""
        logger(self.level).warning(f"{while_replace(text)}")

    def printError(self, text=""):
        """Error log level."""
        logger(self.level).error(f"{while_replace(text)}")

    def printMess(self, text=""):
        """Info log level."""
        logger(self.level).info(f"==> {while_replace(text)}")

    def printDebugQuoted(self, text=""):
        """Print quoted messages."""
        logger(self.level).debug(10 * "+" + "start_quote" + 59 * "+")
        logger(self.level).debug(while_replace(text))
        logger(self.level).debug(10 * "+" + "end_quote" + 61 * "+")

    def printErrorQuoted(self, text=""):
        """Print quoted messages."""
        logger(self.level).error(10 * "+" + "start_quote" + 59 * "+")
        logger(self.level).error(while_replace(text))
        logger(self.level).error(10 * "+" + "end_quote" + 61 * "+")

    def search(self, name: str = "", alist: bool = False):
        """Returns a list with all atomName matching 'name' or just the first case."""
        ll = [x for x in self.atoms if x.atomName == name.upper()]
        if ll and not alist:
            ll = ll[0]
        return ll

    def checkSmiles(self):
        """Check if input arg is a SMILES string.

        Returns:
            bool: True/False
        """
        if find_bin(self.binaries["obabel_bin"]):
            # openbabel-wheel >= 3.1.1.19 is a hard dependency, so only the 3.x API is used.
            from openbabel import openbabel as ob
            from openbabel import pybel

            ob.cvar.obErrorLog.StopLogging()
        else:
            logger(self.level).warning("WARNING: your input may be a SMILES but")
            logger(self.level).warning("         without OpenBabel, this functionality won't work")
            return False

        # Check if input is a smiles string
        try:
            ob.obErrorLog.SetOutputLevel(0)
            pybel.readstring("smi", self.smiles)
            return True
        except Exception:
            ob.obErrorLog.SetOutputLevel(0)

            return False

    def guessCharge(self):
        """Guess the charge of a system based on antechamber.

        Returns None in case of error.
        """
        done = False
        error = False
        charge = self.chargeVal
        localDir = os.path.abspath(".")
        if not os.path.exists(self.tmpDir):
            os.mkdir(self.tmpDir)
        if not os.path.exists(os.path.join(self.tmpDir, self.inputFile)):
            copy2(self.absInputFile, self.tmpDir)
        os.chdir(self.tmpDir)

        if self.chargeType == "user":
            if self.ext == ".mol2":
                self.printMess("Reading user's charges from mol2 file...")
                charge = self.readMol2TotalCharge(self.inputFile)
                done = True
            else:
                self.chargeType = "bcc"
                self.printWarn("cannot read charges from a PDB file")
                self.printWarn("using now 'bcc' method for charge")
        if self.chargeVal is None and not done:
            self.printWarn("no charge value given, trying to guess one...")
            mol2FileForGuessCharge = self.inputFile
            if self.ext == ".pdb":
                cmd = f"'{self.obabelExe}' -ipdb '{self.inputFile}' -omol2 -O '{self.baseName}.mol2'"
                self.printDebug(f"guessCharge: {cmd}")
                out = _getoutput(cmd)
                self.printDebug(out)
                mol2FileForGuessCharge = os.path.abspath(f"{self.baseName}.mol2")
                in_mol = "mol2"
            else:
                in_mol = self.ext[1:]
                if in_mol == "mol":
                    in_mol = "mdl"

            cmd = f"'{self.acExe}' -dr no -i '{mol2FileForGuessCharge}' -fi {in_mol} -o tmp -fo mol2 -c gas -pf n -j {self.predIndex}"

            logger(self.level).debug(while_replace(cmd))

            log = _getoutput(cmd).strip()

            if os.path.exists("tmp"):
                charge = self.readMol2TotalCharge("tmp")
            else:
                try:
                    charge = float(
                        log.strip()
                        .split("equal to the total charge (")[-1]
                        .split(") based on Gasteiger atom type, exit")[0]
                    )
                except Exception:
                    error = True

            if not charge:
                error = True
                charge = 0
            if error:
                self.printError("guessCharge failed")
                os.chdir(localDir)
                rmtree(self.tmpDir)
                self.printErrorQuoted(log)
                self.printMess("Trying with net charge = 0 ...")
        if charge is None:
            msg = "guessCharge failed: no net charge could be determined"
            self.printError(msg)
            raise Exception(msg)
        charge = float(charge)
        charge2 = round(charge)
        drift = abs(charge2 - charge)
        self.printDebug(f"Net charge drift '{drift:3.6f}'")
        if drift > diffTol:
            self.printWarn(f"Net charge drift '{drift:3.5f}' bigger than tolerance '{diffTol:3.5f}'")
            if not self.force and self.chargeType != "user":
                msg = "Error with calculated charge"
                logger(self.level).error(msg)
                rmtree(self.tmpDir)
                raise Exception(msg)
        self.chargeVal = str(charge2)
        self.printMess(f"... charge set to {charge2}")
        os.chdir(localDir)

    def setResNameCheckCoords(self):
        """Set a 3 letter residue name and check coords for issues
        like duplication, atoms too close or too sparse.
        """
        exit_ = False
        localDir = os.path.abspath(".")
        if not os.path.exists(self.tmpDir):
            os.mkdir(self.tmpDir)

        copy2(self.absInputFile, self.tmpDir)
        os.chdir(self.tmpDir)

        exten = self.ext[1:]
        if self.ext == ".pdb":
            tmpFile = open(self.inputFile)
        else:
            if exten == "mol":
                exten = "mdl"
            cmd = f"'{self.acExe}' -dr no -i '{self.inputFile}' -fi {exten} -o tmp -fo ac -pf y -j {self.predIndex}"
            self.printDebug(cmd)
            out = _getoutput(cmd)
            if not out.isspace():
                self.printDebug(out)
            try:
                tmpFile = open("tmp")
            except Exception:
                rmtree(self.tmpDir)
                raise

        tmpData = tmpFile.readlines()
        residues = set()
        coords = {}
        for line in tmpData:
            if "ATOM  " in line or "HETATM" in line:
                residues.add(line[17:20])
                at = line[0:17]
                cs = line[30:54]
                if cs in coords:
                    coords[cs].append(at)
                else:
                    coords[cs] = [at]

        if len(residues) > 1:
            self.printError(f"more than one residue detected '{residues!s}'")
            self.printError(f"verify your input file '{self.inputFile}'. Aborting ...")
            msg = "Only ONE Residue is allowed for ACPYPE to work"
            logger(self.level).error(msg)
            raise Exception(msg)

        dups = ""
        shortd = ""
        longd = ""
        longSet = set()
        id_ = 0
        items = list(coords.items())
        ll = len(items)
        for item in items:
            id_ += 1
            if len(item[1]) > 1:  # if True means atoms with same coordinates
                for i in item[1]:
                    dups += f"{i} {item[0]}\n"
            for id2 in range(id_, ll):
                item2 = items[id2]
                c1 = list(map(float, [item[0][i : i + 8] for i in range(0, 24, 8)]))
                c2 = list(map(float, [item2[0][i : i + 8] for i in range(0, 24, 8)]))
                dist2 = distanceAA(c1, c2)
                if dist2 < minDist2:
                    dist = math.sqrt(dist2)
                    shortd += f"{dist:8.5f}       {item[1]} {item2[1]}\n"
                if dist2 < maxDist2:  # and not longOK:
                    longSet.add(str(item[1]))
                    longSet.add(str(item2[1]))
            if str(item[1]) not in longSet and ll > 1:
                longd += f"{item[1]}\n"

        if dups:
            self.printError(f"Atoms with same coordinates in '{self.inputFile}'!")
            self.printErrorQuoted(dups[:-1])
            exit_ = True

        if shortd:
            self.printError(f"Atoms TOO close (< {minDist} Ang.)")
            self.printErrorQuoted(f"Dist (Ang.)    Atoms\n{shortd[:-1]}")
            exit_ = True

        if longd:
            self.printError(f"Atoms TOO scattered (> {maxDist} Ang.)")
            self.printErrorQuoted(longd[:-1])
            exit_ = True

        if exit_:
            if self.force:
                self.printWarn("You chose to proceed anyway with '-f' option. GOOD LUCK!")
            else:
                self.printError("Use '-f' option if you want to proceed anyway. Aborting ...")
                if not self.debug:
                    rmtree(self.tmpDir)
                msg = "Coordinates issues with your system"
                logger(self.level).error(msg)
                rmtree(self.tmpDir)
                raise Exception(msg)

        # escape resname list index out of range: no RES name in pdb for example
        resname = next(iter(residues)).strip()
        if not resname:
            resname = "LIG"
            self.printWarn("No residue name identified, using default resname: 'LIG'")
        newresname = resname

        # To avoid resname likes: 001 (all numbers), 1e2 (sci number), ADD : reserved terms for leap
        leapWords = [
            "_cmd_options_",
            "_types_",
            "add",
            "addAtomTypes",
            "addIons",
            "addIons2",
            "addPath",
            "addPdbAtomMap",
            "addPdbResMap",
            "alias",
            "alignAxes",
            "bond",
            "bondByDistance",
            "center",
            "charge",
            "check",
            "clearPdbAtomMap",
            "clearPdbResMap",
            "clearVariables",
            "combine",
            "copy",
            "createAtom",
            "createParmset",
            "createResidue",
            "createUnit",
            "crossLink",
            "debugOff",
            "debugOn",
            "debugStatus",
            "deleteBond",
            "deleteOffLibEntry",
            "deleteRestraint",
            "desc",
            "deSelect",
            "displayPdbAtomMap",
            "displayPdbResMap",
            "edit",
            "flip",
            "groupSelectedAtoms",
            "help",
            "impose",
            "list",
            "listOff",
            "loadAmberParams",
            "loadAmberPrep",
            "loadMol2",
            "loadOff",
            "loadPdb",
            "loadPdbUsingSeq",
            "logFile",
            "matchVariables",
            "measureGeom",
            "quit",
            "relax",
            "remove",
            "restrainAngle",
            "restrainBond",
            "restrainTorsion",
            "saveAmberParm",
            "saveAmberParmPert",
            "saveAmberParmPol",
            "saveAmberParmPolPert",
            "saveAmberPrep",
            "saveMol2",
            "saveOff",
            "saveOffParm",
            "savePdb",
            "scaleCharges",
            "select",
            "sequence",
            "set",
            "setBox",
            "solvateBox",
            "solvateCap",
            "solvateDontClip",
            "solvateOct",
            "solvateShell",
            "source",
            "transform",
            "translate",
            "verbosity",
            "zMatrix",
        ]
        isLeapWord = False
        for word in leapWords:
            if resname.upper().startswith(word.upper()):
                self.printDebug(f"Residue name is a reserved word: '{word.upper()}'")
                isLeapWord = True
        try:
            float(resname)
            self.printDebug(f"Residue name is a 'number': '{resname}'")
            isNumber = True
        except ValueError:
            isNumber = False

        if resname[0].isdigit() or isNumber or isLeapWord:
            newresname = "R" + resname
        if not resname.isalnum():
            newresname = "MOL"
        if newresname != resname:
            self.printWarn(
                f"In {self.acBaseName}.lib, residue name will be '{newresname}' instead of '{resname}' elsewhere"
            )

        self.resName = newresname

        os.chdir(localDir)
        self.printDebug("setResNameCheckCoords done")

    def readMol2TotalCharge(self, mol2File):
        """Reads the charges in given mol2 file and returns the total."""
        charge = 0.0
        ll = []
        cmd = f"'{self.acExe}' -dr no -i '{mol2File}' -fi mol2 -o tmp -fo mol2 -c wc -cf tmp.crg -pf n -j {self.predIndex}"

        self.printDebug(cmd)

        log = _getoutput(cmd)

        if os.path.exists("tmp.crg"):
            tmpFile = open("tmp.crg")
            tmpData = tmpFile.readlines()
            for line in tmpData:
                ll += line.split()
            charge = sum(map(float, ll))
        if not log.isspace():
            self.printDebugQuoted(log)

        self.printDebug("readMol2TotalCharge: " + str(charge))

        return charge

    def execAntechamber(self, chargeType=None, atomType=None) -> bool:
        r"""To call Antechamber and execute it.

        Args:
            chargeType ([str], optional): bcc, gas or user. Defaults to None/bcc.
            atomType ([str], optional): gaff, amber, gaff2, amber2. Defaults to None/gaff2.

        Returns:
            bool: True if failed.

        ::

            Welcome to antechamber 21.0: molecular input file processor.

            Usage: antechamber \
                    -i     input file name
                    -fi    input file format
                    -o     output file name
                    -fo    output file format
                    -c     charge method
                    -cf    charge file name
                    -nc    net molecular charge (int)
                    -a     additional file name
                    -fa    additional file format
                    -ao    additional file operation
                            crd   : only read in coordinate
                            crg   : only read in charge
                            radius: only read in radius
                            name  : only read in atom name
                            type  : only read in atom type
                            bond  : only read in bond type
                    -m     multiplicity (2S+1), default is 1
                    -rn    residue name, overrides input file, default is MOL
                    -rf    residue topology file name in prep input file,
                            default is molecule.res
                    -ch    check file name for gaussian, default is 'molecule'
                    -ek    mopac or sqm keyword, inside quotes; overwrites previous ones
                    -gk    gaussian job keyword, inside quotes, is ignored when both -gopt and -gsp are used
                    -gopt  gaussian job keyword for optimization, inside quotes
                    -gsp   gaussian job keyword for single point calculation, inside quotes
                    -gm    gaussian memory keyword, inside quotes, such as "%mem=1000MB"
                    -gn    gaussian number of processors keyword, inside quotes, such as "%nproc=8"
                    -gdsk  gaussian maximum disk usage keyword, inside quotes, such as "%maxdisk=50GB"
                    -gv    add keyword to generate gesp file (for Gaussian 09 only)
                            1    : yes
                            0    : no, the default
                    -ge    gaussian esp file generated by iop(6/50=1), default is g09.gesp
                    -tor   torsional angle list, inside a pair of quotes, such as "1-2-3-4:0,5-6-7-8"
                            ':1' or ':0' indicates the torsional angle is frozen or not
                    -df    am1-bcc precharge flag, 2 - use sqm(default); 0 - use mopac
                    -at    atom type
                            gaff : the default
                            gaff2: for gaff2 (beta-version)
                            amber: for PARM94/99/99SB
                            bcc  : bcc
                            sybyl: sybyl
                    -du    fix duplicate atom names: yes(y)[default] or no(n)
                    -bk    component/block Id, for ccif
                    -an    adjust atom names: yes(y) or no(n)
                            the default is 'y' for 'mol2' and 'ac' and 'n' for the other formats
                    -j     atom type and bond type prediction index, default is 4
                            0    : no assignment
                            1    : atom type
                            2    : full  bond types
                            3    : part  bond types
                            4    : atom and full bond type
                            5    : atom and part bond type
                    -s     status information: 0(brief), 1(default) or 2(verbose)
                    -eq    equalizing atomic charge, default is 1 for '-c resp' and '-c bcc' and 0
                            for the other charge methods
                            0    : no use
                            1    : by atomic paths
                            2    : by atomic paths and structural information, i.e. E/Z configurations
                    -pf    remove intermediate files: yes(y) or no(n)[default]
                    -pl    maximum path length to determin equivalence of atomic charges for resp and bcc,
                            the smaller the value, the faster the algorithm, default is -1 (use full length),
                            set this parameter to 10 to 30 if your molecule is big (# atoms >= 100)
                    -seq   atomic sequence order changable: yes(y)[default] or no(n)
                    -dr    acdoctor mode: yes(y)[default] or no(n)
                    -i -o -fi and -fo must appear; others are optional
                    Use 'antechamber -L' to list the supported file formats and charge methods

                                List of the File Formats

                    file format type  abbre. index | file format type abbre. index
                    --------------------------------------------------------------
                    Antechamber        ac       1  | Sybyl Mol2         mol2    2
                    PDB                pdb      3  | Modified PDB       mpdb    4
                    AMBER PREP (int)   prepi    5  | AMBER PREP (car)   prepc   6
                    Gaussian Z-Matrix  gzmat    7  | Gaussian Cartesian gcrt    8
                    Mopac Internal     mopint   9  | Mopac Cartesian    mopcrt 10
                    Gaussian Output    gout    11  | Mopac Output       mopout 12
                    Alchemy            alc     13  | CSD                csd    14
                    MDL                mdl     15  | Hyper              hin    16
                    AMBER Restart      rst     17  | Jaguar Cartesian   jcrt   18
                    Jaguar Z-Matrix    jzmat   19  | Jaguar Output      jout   20
                    Divcon Input       divcrt  21  | Divcon Output      divout 22
                    SQM Input          sqmcrt  23  | SQM Output         sqmout 24
                    Charmm             charmm  25  | Gaussian ESP       gesp   26
                    Component cif      ccif    27  | GAMESS dat         gamess 28
                    Orca input         orcinp  29  | Orca output        orcout 30
                    --------------------------------------------------------------
                    AMBER restart file can only be read in as additional file.

                                List of the Charge Methods

                    charge method     abbre. index | charge method    abbre. index
                    --------------------------------------------------------------
                    RESP               resp     1  |  AM1-BCC           bcc     2
                    CM1                cm1      3  |  CM2               cm2     4
                    ESP (Kollman)      esp      5  |  Mulliken          mul     6
                    Gasteiger          gas      7  |  Read in charge    rc      8
                    Write out charge   wc       9  |  Delete Charge     dc     10
                    --------------------------------------------------------------
        """
        global pid

        self.printMess("Executing Antechamber...")

        self.makeDir()

        ct = chargeType or self.chargeType
        at = atomType or self.atomType
        if "amber2" in at:
            at = "amber"

        if ct == "user":
            ct = ""
        else:
            ct = f"-c {ct}"

        exten = self.ext[1:]
        if exten == "mol":
            exten = "mdl"

        cmd = f"'{self.acExe}' -dr no -i '{self.inputFile}' -fi {exten} -o '{self.acMol2FileName}' -fo mol2 {ct} -nc {self.chargeVal} -m {self.multiplicity} -j {self.predIndex} -s 2 -df {self.qFlag} -at {at} -pf n {self.ekFlag}"

        self.printDebug(cmd)

        if os.path.exists(self.acMol2FileName) and not self.force:
            self.printMess("AC output file present... doing nothing")
        else:
            try:
                os.remove(self.acMol2FileName)
            except Exception:
                self.printDebug("No file left to be removed")
            signal.signal(signal.SIGALRM, self.signal_handler)
            signal.alarm(self.timeTol)
            p = sub.Popen(cmd, shell=True, stderr=sub.STDOUT, stdout=sub.PIPE)
            pid = p.pid

            out = str(p.communicate()[0].decode())  # p.stdout.read()
            self.acLog = out

        if os.path.exists(self.acMol2FileName):
            self.printMess("* Antechamber OK *")
        else:
            self.printErrorQuoted(self.acLog)
            return True
        return False

    def signal_handler(self, _signum, _frame):  # , pid = 0):
        """Signal handler."""
        global pid
        pids = job_pids_family(pid)
        self.printDebug(f"PID: {pid}, PIDS: {pids}")
        self.printMess(f"Timed out! Process {pids} killed, max exec time ({self.timeTol}s) exceeded")
        # os.system('kill -15 %s' % pids)
        for i in pids.split():
            os.kill(int(i), 15)
        msg = "Semi-QM taking too long to finish... aborting!"
        logger(self.level).error(msg)
        raise Exception(msg)

    def delOutputFiles(self):
        """Delete temporary output files."""
        delFiles = [
            "mopac.in",
            "tleap.in",
            "fixbo.log",
            "addhs.log",
            "ac_tmp_ot.mol2",
            "frcmod.ac_tmp",
            "fragment.mol2",
            self.tmpDir,
        ]  # , 'divcon.pdb', 'mopac.pdb', 'mopac.out'] #'leap.log'
        self.printMess("Removing temporary files...")
        for file_ in delFiles:
            file_ = os.path.join(self.absHomeDir, file_)
            if os.path.exists(file_):
                if os.path.isdir(file_):
                    rmtree(file_)
                else:
                    os.remove(file_)

    def checkXyzAndTopFiles(self):
        """Check XYZ and TOP files."""
        fileXyz = self.acXyzFileName
        fileTop = self.acTopFileName
        return os.path.exists(fileXyz) and os.path.exists(fileTop)

    def execTleap(self):
        """Execute tleap."""
        fail = False

        self.makeDir()

        if self.ext == ".pdb":
            self.printMess("... converting pdb input file to mol2 input file")
            if self.convertPdbToMol2():
                self.printError("convertPdbToMol2 failed")

        if self.execAntechamber():
            self.printError("Antechamber failed")
            fail = True
        else:
            self.retypeForProteinFF()
            if self.checkAtomTypes():
                fail = True
        if not fail and self.execParmchk():
            self.printError("Parmchk failed")
            fail = True

        if fail:
            return True

        tleapScpt = TLEAP_TEMPLATE % self.acParDict

        fp = open("tleap.in", "w")
        fp.write(tleapScpt)
        fp.close()

        cmd = f"'{self.tleapExe}' -f tleap.in"

        if self.checkXyzAndTopFiles() and not self.force:
            self.printMess("Topologies files already present... doing nothing")
        else:
            try:
                os.remove(self.acTopFileName)
                os.remove(self.acXyzFileName)
            except Exception:
                self.printDebug("No crd or prm files left to be removed")
            self.printMess("Executing Tleap...")
            self.printDebug(cmd)
            self.tleapLog = _getoutput(cmd)
            self.checkLeapLog(self.tleapLog)

        if self.checkXyzAndTopFiles():
            self.printMess("* Tleap OK *")
        else:
            missing = sorted(set(re.findall(r"Could not find (\w+) parameter for atom types:\s*(.+)", self.tleapLog)))
            if missing:
                shown = "; ".join(f"{kind} {types.strip()}" for kind, types in missing[:6])
                if len(missing) > 6:
                    shown += f"; ... ({len(missing)} in total)"
                self.printError(f"tleap has no parameters for: {shown}")
            self.printErrorQuoted(self.tleapLog)
            return True
        return False

    def checkLeapLog(self, log):
        """Check Leap log."""
        log = log.splitlines(True)
        check = ""
        block = False
        for line in log:
            # print "*"+line+"*"
            if "Checking '" in line:
                # check += line
                block = True
            if "Checking Unit." in line:
                block = False
            if block:
                check += line
        self.printDebugQuoted(check[:-1])

    def locateDat(self, aFile):
        """Locate a file pertinent to $AMBERHOME/dat/leap/parm/."""
        amberhome = os.environ.get("AMBERHOME")
        if amberhome:
            aFileF = os.path.join(amberhome, "dat/leap/parm", aFile)
            if os.path.exists(aFileF):
                return aFileF
        aFileF = os.path.join(os.path.dirname(self.acExe), "../dat/leap/parm", aFile)
        if os.path.exists(aFileF):
            return aFileF
        return None

    def locateLib(self, aFile):
        """Locate a file pertinent to $AMBERHOME/dat/leap/lib/."""
        amberhome = os.environ.get("AMBERHOME")
        if amberhome:
            aFileF = os.path.join(amberhome, "dat/leap/lib", aFile)
            if os.path.exists(aFileF):
                return aFileF
        aFileF = os.path.join(os.path.dirname(self.acExe), "../dat/leap/lib", aFile)
        if os.path.exists(aFileF):
            return aFileF
        return None

    def retypeForProteinFF(self):
        """Apply the protein force field's residue-specific atom types after antechamber.

        ff14SB types the backbone C-alpha as CX and several side-chain carbons as CO,
        2C, 3C or C8, and keys its torsions on them. antechamber's AMBER typer predates
        those types and emits CT throughout, so every ff14SB backbone torsion would
        otherwise miss, fall back to a zero wildcard, and be back-filled by parmchk2
        from GAFF. Atoms whose residue and atom names match an amino-acid template take
        the template type; anything else keeps what antechamber assigned.
        """
        libs = PROTEIN_FF_LIBS[self.proteinFF]
        if not libs or "amber" not in self.atomType:
            return
        paths = [path for path in (self.locateLib(name) for name in libs) if path]
        if len(paths) != len(libs):
            self.printWarn(f"{self.proteinFF} residue templates not found; keeping antechamber's atom types")
            return
        templates = readAmberLibTypes(paths)
        with open(self.acMol2FileName) as fh:
            lines = fh.readlines()
        lines, changes = retypeMol2Atoms(lines, templates, FF14SB_RETYPES)
        if not changes:
            self.printDebug("no residue template atom types to apply")
            return
        with open(self.acMol2FileName, "w") as fh:
            fh.writelines(lines)
        counts = collections.Counter(f"{old}->{new}" for _, _, old, new in changes)
        summary = ", ".join(f"{change} x{n}" for change, n in sorted(counts.items()))
        self.printMess(f"Applied {self.proteinFF} residue atom types to {len(changes)} atoms: {summary}")

    def checkAtomTypes(self):
        """Fail early, and clearly, when antechamber typed atoms outside the force field.

        antechamber writes ``DU`` for an atom it cannot type, and with ``-at amber`` it
        also emits types such as ``NO`` for a nitro nitrogen that neither the AMBER
        protein set nor GAFF define. Left alone, parmchk2 aborts with an empty frcmod,
        which used to pass as "Parmchk OK", and tleap then fails once per missing bond
        and angle. This names the atoms instead.

        Returns:
            bool: True when the topology cannot be built with the requested atom types.
        """
        if "amber" in self.atomType:
            ff = PROTEIN_FF[self.proteinFF]
            names = (ff["parm"], ff["frcmod"], self.gaffDatfile)
            hint = "This chemistry is outside the AMBER protein atom types; try '-a gaff2'."
        else:
            names = (self.gaffDatfile,)
            hint = f"antechamber could not assign {self.atomType} types to it."
        files = [path for path in (self.locateDat(name) for name in names) if path]
        if len(files) != len(names):
            return False  # cannot check without the parameter files; let tleap decide
        with open(self.acMol2FileName) as fh:
            offending = unknownMol2Types(fh.readlines(), readParmAtomTypes(files))
        if not offending:
            return False
        shown = ", ".join(f"{name}({atomType})" for name, atomType in offending[:8])
        if len(offending) > 8:
            shown += f", ... ({len(offending)} in total)"
        self.printError(f"No parameters for atom(s) {shown}. {hint}")
        return True

    def execParmchk(self):
        """Execute parmchk."""
        self.makeDir()
        cmdBase = f"'{self.parmchkExe}' -i '{self.acMol2FileName}' -f mol2"

        if "amber" in self.atomType:
            ff = PROTEIN_FF[self.proteinFF]
            gaffFile = self.locateDat(self.gaffDatfile)
            parmfile = self.locateDat(ff["parm"])
            frcmodSB = self.locateDat(ff["frcmod"])
            # Two passes, because parmchk2 does not only fill gaps. Given a parameter
            # file with GAFF merged in it also replaces AMBER terms it already has with
            # GAFF analogues its correspondence table scores as equivalent -- the amide
            # X-C-N-X torsion by c3-c-n-c3, the carboxylate H-CT-C-O2 by hc-c3-c-o --
            # and tleap loads the frcmod last, so those analogues win. The first pass,
            # against AMBER alone, decides what is genuinely missing; the second, against
            # the merge, supplies GAFF values for exactly those entries and nothing else.
            amberOnly = parmMerge(parmfile, frcmodSB, frcmod=True)
            withGaff = parmMerge(parmMerge(parmfile, gaffFile), frcmodSB, frcmod=True)
            gapsFrcmod = self.acFrcmodFileName + ".amber"
            gaffFrcmod = self.acFrcmodFileName + ".gaff"
            self.parmchkLog = ""
            for out, parm in ((gapsFrcmod, amberOnly), (gaffFrcmod, withGaff)):
                cmd = f"{cmdBase} -o '{out}' -p '{parm}'"
                self.printDebug(cmd)
                self.parmchkLog += _getoutput(cmd)
            if os.path.exists(gapsFrcmod) and os.path.exists(gaffFrcmod):
                with open(gapsFrcmod) as gaps, open(gaffFrcmod) as filled:
                    merged = mergeFrcmodGaps(gaps.readlines(), filled.readlines())
                with open(self.acFrcmodFileName, "w") as fh:
                    fh.writelines(merged)
                if not self.debug:
                    os.remove(gapsFrcmod)
                    os.remove(gaffFrcmod)
        else:
            cmd = f"{cmdBase} -o '{self.acFrcmodFileName}'"
            if "gaff2" in self.atomType:
                cmd += " -s 2"
            self.printDebug(cmd)
            self.parmchkLog = _getoutput(cmd)

        unknown = sorted(set(re.findall(r"Atom type of (\S+) does not exist in PARMCHK.DAT", self.parmchkLog)))
        if unknown:
            self.printError(f"parmchk2 does not know atom type(s) {', '.join(unknown)}; its frcmod would be incomplete")
            return True

        if os.path.exists(self.acFrcmodFileName):
            check = self.checkFrcmod()
            if check:
                self.printWarn("Couldn't determine all parameters:")
                self.printMess(f"From file '{self.acFrcmodFileName + check}'\n")
            else:
                self.printMess("* Parmchk OK *")
        else:
            self.printErrorQuoted(self.parmchkLog)
            return True
        return False

    def checkFrcmod(self):
        """Check FRCMOD file."""
        check = ""
        frcmodContent = open(self.acFrcmodFileName).readlines()
        for line in frcmodContent:
            if "ATTN, need revision" in line:
                check += line
        return check

    def convertPdbToMol2(self):
        """Convert PDB to MOL2 by using obabel."""
        if self.ext == ".pdb":
            if self.execObabel():
                self.printError(f"convert pdb to mol2 via {binaries['obabel_bin']} failed")
                return True
        return False

    def convertSmilesToMol2(self):
        """Convert Smiles to MOL2 by using obabel."""
        # if not self.obabelExe:
        #     msg = "SMILES needs OpenBabel python module"
        #     logger(self.level).error(msg)
        #     raise Exception(msg)

        from openbabel import pybel

        try:
            mymol = pybel.readstring("smi", str(self.smiles))
            mymol.addh()
            mymol.make3D()
            mymol.write(self.ext.replace(".", ""), self.absInputFile, overwrite=True)
            return True
        except Exception:
            return False

    def execObabel(self):
        """Execute obabel."""
        self.makeDir()

        cmd = f"'{self.obabelExe}' -ipdb '{self.inputFile}' -omol2 -O '{self.baseName}.mol2'"
        self.printDebug(cmd)
        self.obabelLog = _getoutput(cmd)
        self.ext = ".mol2"
        self.inputFile = self.baseName + self.ext
        self.acParDict["ext"] = "mol2"
        if os.path.exists(self.inputFile):
            self.printMess("* Babel OK *")
        else:
            self.printErrorQuoted(self.obabelLog)
            return True
        return False

    def makeDir(self):
        """Make Dir."""
        os.chdir(self.rootDir)
        self.absHomeDir = os.path.abspath(self.homeDir)
        if not os.path.exists(self.homeDir):
            os.mkdir(self.homeDir)
        os.chdir(self.homeDir)
        if self.absInputFile:
            copy2(self.absInputFile, ".")
        return True

    def createACTopol(self):
        """Build the AMBER prmtop/inpcrd pair.

        Returns:
            bool: True when any step failed, so callers can stop before converting.
        """
        failed = self.execTleap()
        if failed and self.tleapLog:
            # Only when tleap itself ran; earlier steps report their own failure.
            self.printError("Tleap failed")
        if not self.debug:
            self.delOutputFiles()
        return failed

    def createMolTopol(self):
        """Create MolTopol obj."""
        self.topFileData = open(self.acTopFileName).readlines()
        self.molTopol = MolTopol(
            self,  # acTopolObj
            verbose=self.verbose,
            debug=self.debug,
            gmx4=self.gmx4,
            merge=self.merge,
            direct=self.direct,
            is_sorted=self.sorted,
            chiral=self.chiral,
        )
        if self.outTopols:
            if "cns" in self.outTopols:
                self.molTopol.writeCnsTopolFiles()
            if "gmx" in self.outTopols:
                self.molTopol.writeGromacsTopolFiles()
            if "charmm" in self.outTopols:
                self.writeCharmmTopolFiles()
        try:  # scape the pickle save error
            self.pickleSave()
        except Exception:
            self.printError("pickleSave failed")

        if not self.debug:
            self.delOutputFiles()  # required to use on Jupyter Notebook
        os.chdir(self.rootDir)

    def pickleSave(self):
        """Create portable serialized representations of System's Python objects.

        Example:
            to restore:

            .. code-block:: python

                from acpype import *

                # import cPickle as pickle
                import pickle

                mol = pickle.load(open("DDD.pkl", "rb"))
        """
        pklFile = self.baseName + ".pkl"
        dumpFlag = False
        if not os.path.exists(pklFile):
            mess = "Writing pickle file %s" % pklFile
            dumpFlag = True
        elif self.force:
            mess = "Overwriting pickle file %s" % pklFile
            dumpFlag = True
        else:
            mess = "Pickle file %s already present... doing nothing" % pklFile
        self.printMess(mess)
        if dumpFlag:
            with open(pklFile, "wb") as f:  # for python 3.3 or higher
                pickle.dump(self, f)

    def getFlagData(self, flag):
        """For a given acFileTop flag, return a list of the data related."""

        def proc_line(line):
            # data need format
            data = line.rstrip()
            sdata = [data[i : i + f].strip() for i in range(0, len(data), f)]
            if "+" and "." in data and flag != "RESIDUE_LABEL":  # it's a float
                ndata = list(map(float, sdata))
            elif flag != "RESIDUE_LABEL":
                try:  # try if it's integer
                    ndata = list(map(int, sdata))
                except Exception:
                    ndata = sdata
            else:
                ndata = sdata

            return ndata

        block = False
        tFlag = "%FLAG " + flag
        ndata = []

        if not self.topFileData:
            msg = "PRMTOP file empty?"
            logger(self.level).error(msg)
            raise Exception(msg)

        for rawLine in self.topFileData:
            line = rawLine.replace("\r", "").replace("\n", "")
            if tFlag in line:
                block = True
                continue
            if block and "%FLAG " in line:
                break
            if block:
                if "%FORMAT" in line:
                    line = line.strip().strip("%FORMAT()").split(".")[0]
                    for c in line:
                        if c.isalpha():
                            f = int(line.split(c)[1])
                            break
                    continue
                ndata += proc_line(line)

        if flag == "AMBER_ATOM_TYPE":
            nn = []
            ll = set()
            prefixed = False
            for ii in ndata:
                prefixed = True
                if ii[0].isdigit():
                    ll.add(ii)
                    ii = "A" + ii
                nn.append(ii)
            if prefixed and ll:
                self.printDebug("GMX does not like atomtype starting with Digit")
                self.printDebug("prefixing AtomType %s with 'A'." % list(ll))
            ndata = nn
        return ndata  # a list

    def getResidueLabel(self):
        """Get a 3 capital letters code from acFileTop.

        Returns a list.
        """
        residueLabel = self.getFlagData("RESIDUE_LABEL")
        residueLabel = list(map(str, residueLabel))
        if residueLabel[0] != residueLabel[0].upper():
            self.printWarn(f"residue label '{residueLabel[0]}' in '{self.inputFile}' is not all UPPERCASE")
            self.printWarn("this may raise problem with some applications like CNS")
        self.residueLabel = residueLabel

    def getCoords(self):
        """For a given acFileXyz file, return a list of coords as::

        [[x1,y1,z1],[x2,y2,z2], etc.]
        """
        if not self.xyzFileData:
            msg = "INPCRD file empty?"
            logger(self.level).error(msg)
            raise Exception(msg)
        data = ""
        for rawLine in self.xyzFileData[2:]:
            line = rawLine.replace("\r", "").replace("\n", "")
            data += line
        ll = len(data)
        ndata = list(map(float, [data[i : i + 12] for i in range(0, ll, 12)]))

        gdata = []
        for i in range(0, len(ndata), 3):
            gdata.append([ndata[i], ndata[i + 1], ndata[i + 2]])

        self.printDebug("getCoords done")

        return gdata

    def getAtoms(self):
        """Set a list with all atoms objects built from dat in acFileTop.

        Set also atomType system is gaff or amber, list of atomTypes, resid
        and system's total charge.
        """
        atomNameList = self.getFlagData("ATOM_NAME")
        atomTypeNameList = self.getFlagData("AMBER_ATOM_TYPE")
        self._atomTypeNameList = atomTypeNameList
        massList = self.getFlagData("MASS")
        chargeList = self.getFlagData("CHARGE")
        resIds = [*self.getFlagData("RESIDUE_POINTER"), 0]
        coords = self.getCoords()
        ACOEFs, BCOEFs = self.getABCOEFs()

        atoms = []
        atomTypes = []
        tmpList = []  # a list with unique atom types
        totalCharge = 0.0
        countRes = 0
        id_ = 0
        FirstNonSoluteId = None
        for atomName in atomNameList:
            if atomName != atomName.upper():
                self.printDebug("atom name '%s' HAS to be all UPPERCASE... Applying this here." % atomName)
                atomName = atomName.upper()
            atomTypeName = atomTypeNameList[id_]
            if id_ + 1 == resIds[countRes]:
                resid = countRes
                countRes += 1
            resName = self.residueLabel[resid]
            if resName in ionOrSolResNameList and not FirstNonSoluteId:
                FirstNonSoluteId = id_
            mass = massList[id_]
            charge = chargeList[id_]
            chargeConverted = charge / qConv
            totalCharge += charge
            coord = coords[id_]
            ACOEF = ACOEFs[id_]
            BCOEF = BCOEFs[id_]
            atomType = AtomType(atomTypeName, mass, ACOEF, BCOEF)
            if atomTypeName not in tmpList:
                tmpList.append(atomTypeName)
                atomTypes.append(atomType)
            atom = Atom(atomName, atomType, id_ + 1, resid, mass, chargeConverted, coord)
            atoms.append(atom)
            id_ += 1

        soluteRanges = []
        if self.splitMolecules:
            soluteRanges = self.amberMoleculeRanges(FirstNonSoluteId or len(chargeList))
        if soluteRanges:
            self.balanceChargesPerMolecule(chargeList, atoms, soluteRanges)
            balanceChargeList = chargeList
        else:
            balanceChargeList, balanceValue, balanceIds = self.balanceCharges(chargeList, FirstNonSoluteId)

            for id_ in balanceIds:
                atoms[id_].charge = balanceValue / qConv

        if atomTypeName[0].islower():
            self.atomTypeSystem = "gaff"
        else:
            self.atomTypeSystem = "amber"
        self.printDebug("Balanced TotalCharge %13.10f" % float(sum(balanceChargeList) / qConv))

        self.totalCharge = round(totalCharge / qConv)

        self.atoms = atoms
        self.atomTypes = atomTypes

        self.pbc = []
        if len(coords) == len(atoms) + 2 or len(coords) == len(atoms) * 2 + 2:
            self.pbc = [coords[-2], coords[-1]]
        self.printDebug("PBC = %s" % self.pbc)
        self.printDebug("getAtoms done")

    def amberMoleculeRanges(self, nSoluteAtoms):
        """Split the solute into the molecules AMBER already identified.

        ``ATOMS_PER_MOLECULE`` is written by tleap from bonded connectivity, and its
        atoms are contiguous by construction, which is exactly what GROMACS needs: the
        coordinate file order has to match the order of ``[ molecules ]``. Working from
        that flag means no connectivity analysis and no reordering here.

        The flag only exists when the prmtop carries a box, so a vacuum system yields
        nothing and the caller keeps its single moleculetype.

        Args:
            nSoluteAtoms: how many leading atoms belong to the solute, i.e. everything
                before the first water or ion handled by its own template.

        Returns:
            list: ``(lo, hi)`` pairs of 0-based, half-open atom ranges covering exactly
            the solute, or ``[]`` when the split cannot be made safely.
        """
        perMolecule = self.getFlagData("ATOMS_PER_MOLECULE")
        if not perMolecule:
            return []
        ranges = []
        lo = 0
        for count in perMolecule:
            if lo >= nSoluteAtoms:
                break
            ranges.append((lo, lo + count))
            lo += count
        # The molecule boundaries must land exactly on the solute boundary. If they do
        # not, the solvent split disagrees with AMBER's and splitting would silently
        # move atoms between moleculetypes.
        if lo != nSoluteAtoms:
            self.printDebug("ATOMS_PER_MOLECULE does not align with the solute; not splitting")
            return []
        return ranges

    def balanceChargesPerMolecule(self, chargeList, atoms, ranges):
        """Round each molecule's charge to an integer on its own.

        AMBER stores charges to a limited precision, so a molecule taken out of a
        system is typically off by ~1e-3: ILDN's two chains carry +2.001 each and its
        DMP ligand -0.002, summing to an exact +4. Balancing globally, as
        :meth:`balanceCharges` does, hides that inside one moleculetype. Once the
        molecules are written separately each has to be integral by itself, or grompp
        reports a non-integer charge for every one of them.
        """
        for lo, hi in ranges:
            subList = chargeList[lo:hi]
            balanced, value, ids = self.balanceCharges(subList)
            chargeList[lo:hi] = balanced
            for id_ in ids:
                atoms[lo + id_].charge = value / qConv

    def getBonds(self):
        """Get Bonds."""
        uniqKbList = self.getFlagData("BOND_FORCE_CONSTANT")
        uniqReqList = self.getFlagData("BOND_EQUIL_VALUE")
        bondCodeHList = self.getFlagData("BONDS_INC_HYDROGEN")
        bondCodeNonHList = self.getFlagData("BONDS_WITHOUT_HYDROGEN")
        bondCodeList = bondCodeHList + bondCodeNonHList
        bonds = []
        for i in range(0, len(bondCodeList), 3):
            idAtom1 = bondCodeList[i] // 3  # remember python starts with id 0
            idAtom2 = bondCodeList[i + 1] // 3
            bondTypeId = bondCodeList[i + 2] - 1
            atom1 = self.atoms[idAtom1]
            atom2 = self.atoms[idAtom2]
            kb = uniqKbList[bondTypeId]
            req = uniqReqList[bondTypeId]
            atoms = [atom1, atom2]
            bond = Bond(atoms, kb, req)
            bonds.append(bond)
        self.bonds = bonds
        self.printDebug("getBonds done")

    def getAngles(self):
        """Get Angles."""
        uniqKtList = self.getFlagData("ANGLE_FORCE_CONSTANT")
        uniqTeqList = self.getFlagData("ANGLE_EQUIL_VALUE")
        # for list below, true atom number = index/3 + 1
        angleCodeHList = self.getFlagData("ANGLES_INC_HYDROGEN")
        angleCodeNonHList = self.getFlagData("ANGLES_WITHOUT_HYDROGEN")
        angleCodeList = angleCodeHList + angleCodeNonHList
        angles = []
        for i in range(0, len(angleCodeList), 4):
            idAtom1 = angleCodeList[i] // 3  # remember python starts with id 0
            idAtom2 = angleCodeList[i + 1] // 3
            idAtom3 = angleCodeList[i + 2] // 3
            angleTypeId = angleCodeList[i + 3] - 1
            atom1 = self.atoms[idAtom1]
            atom2 = self.atoms[idAtom2]
            atom3 = self.atoms[idAtom3]
            kt = uniqKtList[angleTypeId]
            teq = uniqTeqList[angleTypeId]  # angle given in rad in prmtop
            atoms = [atom1, atom2, atom3]
            angle = Angle(atoms, kt, teq)
            angles.append(angle)
        self.angles = angles
        self.printDebug("getAngles done")

    def getDihedrals(self):
        """Get dihedrals (proper and imp), condensed list of prop dih and
        atomPairs.
        """
        uniqKpList = self.getFlagData("DIHEDRAL_FORCE_CONSTANT")
        uniqPeriodList = self.getFlagData("DIHEDRAL_PERIODICITY")
        uniqPhaseList = self.getFlagData("DIHEDRAL_PHASE")
        # for list below, true atom number = abs(index)/3 + 1
        dihCodeHList = self.getFlagData("DIHEDRALS_INC_HYDROGEN")
        dihCodeNonHList = self.getFlagData("DIHEDRALS_WITHOUT_HYDROGEN")
        dihCodeList = dihCodeHList + dihCodeNonHList
        properDih = []
        improperDih = []
        condProperDih = []  # list of dihedrals condensed by the same quartet
        # atomPairs = []
        atomPairs = set()
        for i in range(0, len(dihCodeList), 5):
            idAtom1 = dihCodeList[i] // 3  # remember python starts with id 0
            idAtom2 = dihCodeList[i + 1] // 3
            # 3 and 4 indexes can be negative: if id3 < 0, end group interations
            # in amber are to be ignored; if id4 < 0, dihedral is improper
            idAtom3raw = dihCodeList[i + 2] // 3  # can be negative -> exclude from 1-4vdw
            idAtom4raw = dihCodeList[i + 3] // 3  # can be negative -> Improper
            idAtom3 = abs(idAtom3raw)
            idAtom4 = abs(idAtom4raw)
            dihTypeId = dihCodeList[i + 4] - 1
            atom1 = self.atoms[idAtom1]
            atom2 = self.atoms[idAtom2]
            atom3 = self.atoms[idAtom3]
            atom4 = self.atoms[idAtom4]
            kPhi = uniqKpList[dihTypeId]  # already divided by IDIVF
            period = int(uniqPeriodList[dihTypeId])  # integer
            phase = uniqPhaseList[dihTypeId]  # angle given in rad in prmtop
            if phase == kPhi == 0:
                period = 0  # period is set to 0
            atoms = [atom1, atom2, atom3, atom4]
            dihedral = Dihedral(atoms, kPhi, period, phase)
            if idAtom4raw > 0:
                try:
                    atomsPrev = properDih[-1].atoms
                except Exception:
                    atomsPrev = []
                properDih.append(dihedral)
                if idAtom3raw < 0 and atomsPrev == atoms:
                    condProperDih[-1].append(dihedral)
                else:
                    condProperDih.append([dihedral])
                pair = (atom1, atom4)
                # if atomPairs.count(pair) == 0 and idAtom3raw > 0:
                if idAtom3raw > 0:
                    atomPairs.add(pair)
            else:
                improperDih.append(dihedral)
        if self.sorted:
            atomPairs = sorted(atomPairs, key=lambda x: (x[0].id, x[1].id))
        self.properDihedrals = properDih
        self.improperDihedrals = improperDih
        self.condensedProperDihedrals = condProperDih  # [[],[],...]
        self.atomPairs = atomPairs  # set((atom1, atom2), ...)
        self.printDebug("getDihedrals done")

    def getChirals(self):
        """Get chiral atoms (for CNS only!).

        Plus its 4 neighbours and improper dihedral angles
        to store non-planar improper dihedrals.
        """
        if not self._parent.obabelExe:
            self.printWarn("No Openbabel python module, no chiral groups")
            self.chiralGroups = []
            return

        from openbabel import openbabel as ob
        from openbabel import pybel

        self.printMess("Using OpenBabel v." + ob.OBReleaseVersion() + "\n")

        " obchiral script - replace the obchiral executable"
        out = []
        _filename, file_extension = os.path.splitext(self.inputFile)
        mol = pybel.readfile(file_extension.replace(".", ""), self.inputFile)
        for ml in mol:
            for at in ml:
                if ob.OBStereoFacade(ml.OBMol).HasTetrahedralStereo(at.idx):
                    out.append(at.idx)
        " end of obchiral script "

        chiralGroups = []
        for id_ in out:
            atChi = self.atoms[id_]
            quad = []
            for bb in self.bonds:
                bAts = bb.atoms[:]
                if atChi in bAts:
                    bAts.remove(atChi)
                    quad.append(bAts[0])
            if len(quad) != 4:
                if self.chiral:
                    self.printWarn(f"Atom {atChi} has less than 4 connections to 4 different atoms. It's NOT Chiral!")
                continue
            v1, v2, v3, v4 = (x.coords for x in quad)
            chiralGroups.append((atChi, quad, imprDihAngle(v1, v2, v3, v4)))
        self.chiralGroups = chiralGroups

    def sortAtomsForGromacs(self):
        """Re-sort atoms for gromacs, which expects all hydrogens to immediately
        follow the heavy atom they are bonded to and belong to the same charge
        group.

        Currently, atom mass < 1.2 is taken to denote a proton.  This behaviour
        may be changed by modifying the 'is_hydrogen' function within.

        JDC 2011-02-03
        """
        # Build dictionary of bonded atoms.
        bonded_atoms = dict()
        for atom in self.atoms:
            bonded_atoms[atom] = list()
        for bond in self.bonds:
            [atom1, atom2] = bond.atoms
            bonded_atoms[atom1].append(atom2)
            bonded_atoms[atom2].append(atom1)

        # Define hydrogen and heavy atom classes.
        def is_hydrogen(atom):
            """Check for H."""
            return atom.mass < 1.2

        def is_heavy(atom):
            """Check for non H."""
            return not is_hydrogen(atom)

        # Build list of sorted atoms, assigning charge groups by heavy atom.
        sorted_atoms = list()
        cgnr = 1  # charge group number: each heavy atoms is assigned its own charge group
        # First pass: add heavy atoms, followed by the hydrogens bonded to them.
        for atom in self.atoms:
            if is_heavy(atom):
                # Append heavy atom.
                atom.cgnr = cgnr
                sorted_atoms.append(atom)
                # Append all hydrogens.
                for bonded_atom in bonded_atoms[atom]:
                    if is_hydrogen(bonded_atom) and bonded_atom not in sorted_atoms:
                        # Append bonded hydrogen.
                        bonded_atom.cgnr = cgnr
                        sorted_atoms.append(bonded_atom)
                cgnr += 1

        # Second pass: Add any remaining atoms.
        if len(sorted_atoms) < len(self.atoms):
            for atom in self.atoms:
                if atom not in sorted_atoms:
                    atom.cgnr = cgnr
                    sorted_atoms.append(atom)
                    cgnr += 1

        # Replace current list of atoms with sorted list.
        self.atoms = sorted_atoms

        # Renumber atoms in sorted list, starting from 1.
        for index, atom in enumerate(self.atoms):
            atom.id = index + 1

    def balanceCharges(self, chargeList, FirstNonSoluteId=None):
        """Spread charges fractions among atoms to balance system's total charge.

        Note that python is very annoying about floating points.
        Even after balance, there will always be some residue of order :math:`e^{-12}`
        to :math:`e^{-16}`, which is believed to vanished once one writes a topology
        file, say, for CNS or GMX, where floats are represented with 4 or 5
        maximum decimals.
        """
        limIds = []
        total = sum(chargeList)
        totalConverted = total / qConv
        self.printDebug("charge to be balanced: total %13.10f" % (totalConverted))
        maxVal = max(chargeList[:FirstNonSoluteId])
        minVal = min(chargeList[:FirstNonSoluteId])
        if abs(maxVal) >= abs(minVal):
            lim = maxVal
        else:
            lim = minVal
        nLims = chargeList.count(lim)
        diff = totalConverted - round(totalConverted)
        fix = lim - diff * qConv / nLims
        id_ = 0
        for c in chargeList:
            if c == lim:
                limIds.append(id_)
                chargeList[id_] = fix
            id_ += 1
        self.printDebug("balanceCharges done")
        return chargeList, fix, limIds

    def getABCOEFs(self):
        """Get non-bonded coefficients."""
        uniqAtomTypeIdList = self.getFlagData("ATOM_TYPE_INDEX")
        nonBonIdList = self.getFlagData("NONBONDED_PARM_INDEX")
        rawACOEFs = self.getFlagData("LENNARD_JONES_ACOEF")
        rawBCOEFs = self.getFlagData("LENNARD_JONES_BCOEF")
        ACOEFs = []
        BCOEFs = []
        ntypes = max(uniqAtomTypeIdList)
        for id_ in range(len(self._atomTypeNameList)):
            atomTypeId = uniqAtomTypeIdList[id_]
            index = ntypes * (atomTypeId - 1) + atomTypeId
            nonBondId = nonBonIdList[index - 1]
            ACOEFs.append(rawACOEFs[nonBondId - 1])
            BCOEFs.append(rawBCOEFs[nonBondId - 1])
        self.printDebug("getABCOEFs done")
        return ACOEFs, BCOEFs

    def setProperDihedralsCoef(self):
        """Create proper dihedrals list with Ryckaert-Bellemans coefficients.

        It takes self.condensedProperDihedrals and returns
        self.properDihedralsCoefRB, a reduced list of quartet atoms + RB.
        Coefficients ready for GMX (multiplied by 4.184)::

            self.properDihedralsCoefRB = [[atom1, ..., atom4], C[0:5]]

        For proper dihedrals: a quartet of atoms may appear with more than
        one set of parameters and to convert to GMX they are treated as RBs.

        The resulting coefficients calculated here may look slighted different
        from the ones calculated by amb2gmx.pl because python is taken full float
        number from prmtop and not rounded numbers from rdparm.out as amb2gmx.pl does.
        """
        properDihedralsCoefRB = []
        properDihedralsAlphaGamma = []
        properDihedralsGmx45 = []
        for item in self.condensedProperDihedrals:
            V = 6 * [0.0]
            C = 6 * [0.0]
            for dih in item:
                period = dih.period  # Pn
                kPhi = dih.kPhi  # in rad
                phaseRaw = dih.phase * radPi  # in degree
                phase = int(phaseRaw)  # in degree
                if period > 4 and self.gmx4:
                    rmtree(self.absHomeDir)
                    msg = "Likely trying to convert ILDN to RB, DO NOT use option '-z'"
                    logger(self.level).error(msg)
                    raise Exception(msg)
                if phase in [0, 180]:
                    properDihedralsGmx45.append([item[0].atoms, phaseRaw, kPhi, period])
                    if self.gmx4:
                        if kPhi != 0:
                            V[period] = 2 * kPhi * cal
                        if period == 1:
                            C[0] += 0.5 * V[period]
                            if phase == 0:
                                C[1] -= 0.5 * V[period]
                            else:
                                C[1] += 0.5 * V[period]
                        elif period == 2:
                            if phase == 180:
                                C[0] += V[period]
                                C[2] -= V[period]
                            else:
                                C[2] += V[period]
                        elif period == 3:
                            C[0] += 0.5 * V[period]
                            if phase == 0:
                                C[1] += 1.5 * V[period]
                                C[3] -= 2 * V[period]
                            else:
                                C[1] -= 1.5 * V[period]
                                C[3] += 2 * V[period]
                        elif period == 4:
                            if phase == 180:
                                C[2] += 4 * V[period]
                                C[4] -= 4 * V[period]
                            else:
                                C[0] += V[period]
                                C[2] -= 4 * V[period]
                                C[4] += 4 * V[period]
                else:
                    properDihedralsAlphaGamma.append([item[0].atoms, phaseRaw, kPhi, period])
                    # print phaseRaw, kPhi, period
            if phase in [0, 180]:
                properDihedralsCoefRB.append([item[0].atoms, C])

        self.printDebug("setProperDihedralsCoef done")

        self.properDihedralsCoefRB = properDihedralsCoefRB
        self.properDihedralsAlphaGamma = properDihedralsAlphaGamma
        self.properDihedralsGmx45 = properDihedralsGmx45

    def writeCharmmTopolFiles(self):
        """Write CHARMM topology files."""
        self.printMess("Writing CHARMM files\n")

        at = self.atomType
        self.getResidueLabel()
        res = self.resName

        cmd = f"'{self.acExe}' -dr no -i '{self.acMol2FileName}' -fi mol2 -o '{self.charmmBase}' \
            -fo charmm -s 2 -at {at} -pf n -rn {res} -j {self.predIndex}"

        self.printDebug(cmd)

        log = _getoutput(cmd)
        self.printDebugQuoted(log)

    def writePdb(self, afile):
        """Write a new PDB file with the atom names defined by Antechamber.

        The format generated here use is slightly different from:

            old: http://www.wwpdb.org/documentation/file-format-content/format23/sect9.html
            latest: http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html

        respected to atom name. Using GAFF2 atom types::

            CU/Cu Copper, CL/cl Chlorine, BR/br Bromine

        Args:
            afile ([str]): file path name
        """
        # TODO: assuming only one residue ('1')
        pdbFile = open(afile, "w")
        fbase = os.path.basename(afile)
        pdbFile.write("REMARK " + head % (fbase, date))
        id_ = 1
        for atom in self.atoms:
            # id_ = self.atoms.index(atom) + 1
            aName = atom.atomName
            if atom.atomType.atomTypeName.upper() in specialGaffAtoms:
                s = atom.atomType.atomTypeName.upper()
            else:
                s = atom.atomType.atomTypeName[0].upper()

            rName = self.residueLabel[0]
            x = atom.coords[0]
            y = atom.coords[1]
            z = atom.coords[2]
            line = "%-6s%5d %4s %3s Z%4d%s%8.3f%8.3f%8.3f%6.2f%6.2f%s%2s\n" % (
                "ATOM",
                id_,
                aName,
                rName,
                1,
                4 * " ",
                x,
                y,
                z,
                1.0,
                0.0,
                10 * " ",
                s,
            )
            pdbFile.write(line)
            id_ += 1
        pdbFile.write("END\n")

    def writeGromacsTopolFiles(self):
        """Write GMX topology Files.

        ::

            # from ~/Programmes/amber10/dat/leap/parm/gaff.dat
            #atom type        atomic mass        atomic polarizability        comments
            ca                12.01                 0.360                    Sp2 C in pure aromatic systems
            ha                1.008                 0.135                    H bonded to aromatic carbon

            #bonded atoms        harmonic force kcal/mol/A^2       eq. dist. Ang.  comments
            ca-ha                  344.3*                           1.087**         SOURCE3  1496    0.0024    0.0045
            * for gmx: 344.3 * 4.184 * 100 * 2 = 288110 kJ/mol/nm^2 (why factor 2?)
            ** convert Ang to nm ( div by 10) for gmx: 1.087 A = 0.1087 nm
            # CA HA      1    0.10800   307105.6 ; ged from 340. bsd on C6H6 nmodes; PHE,TRP,TYR (from ffamber99bon.itp)
            # CA-HA  367.0    1.080       changed from 340. bsd on C6H6 nmodes; PHE,TRP,TYR (from parm99.dat)

            # angle        HF kcal/mol/rad^2    eq angle degrees     comments
            ca-ca-ha        48.5*             120.01                SOURCE3 2980   0.1509   0.2511
            * to convert to gmx: 48.5 * 4.184 * 2 = 405.848 kJ/mol/rad^2 (why factor 2?)
            # CA  CA  HA           1   120.000    418.400 ; new99 (from ffamber99bon.itp)
            # CA-CA-HA    50.0      120.00 (from parm99.dat)

            # dihedral    idivf        barrier hight/2 kcal/mol  phase degrees       periodicity     comments
            X -ca-ca-X    4           14.500*                    180.000             2.000           intrpol.bsd.on C6H6
            *convert 2 gmx: 14.5/4 * 4.184 * 2 (?) (yes in amb2gmx, not in topolbuild, why?) = 30.334 or 15.167 kJ/mol
            # X -CA-CA-X    4   14.50        180.0     2.         intrpol.bsd.on C6H6 (from parm99.dat)
            # X CA CA  X    3   30.334       0.000   -30.33400     0.000     0.000     0.000   ; intrpol.bsd.on C6H6
            ;propers treated as RBs in GMX to use combine multiple AMBER torsions per quartet (from ffamber99bon.itp)

            # impr. dihedral        barrier hight/2      phase degrees       periodicity     comments
            X -X -ca-ha             1.1*                  180.                      2.       bsd.on C6H6 nmodes
            * to convert to gmx: 1.1 * 4.184 = 4.6024 kJ/mol/rad^2
            # X -X -CA-HA         1.1          180.          2.           bsd.on C6H6 nmodes (from parm99.dat)
            # X   X   CA  HA       1      180.00     4.60240     2      ; bsd.on C6H6 nmodes
            ;impropers treated as propers in GROMACS to use correct AMBER analytical function (from ffamber99bon.itp)

            # 6-12 parms     sigma = 2 * r * 2^(-1/6)    epsilon
            # atomtype        radius Ang.                    pot. well depth kcal/mol      comments
            ha                  1.4590*                      0.0150**                         Spellmeyer
            ca                  1.9080                    0.0860                            OPLS
            * to convert to gmx:
                sigma = 1.4590 * 2^(-1/6) * 2 = 2 * 1.29982 Ang. = 2 * 0.129982 nm  = 1.4590 * 2^(5/6)/10 =  0.259964 nm
            ** to convert to gmx: 0.0150 * 4.184 = 0.06276 kJ/mol
            # amber99_3    CA     0.0000  0.0000  A   3.39967e-01  3.59824e-01 (from ffamber99nb.itp)
            # amber99_22   HA     0.0000  0.0000  A   2.59964e-01  6.27600e-02 (from ffamber99nb.itp)
            # C*          1.9080  0.0860             Spellmeyer
            # HA          1.4590  0.0150             Spellmeyer (from parm99.dat)
            # to convert r and epsilon to ACOEF and BCOEF
            # ACOEF = sqrt(e1*e2) * (r1 + r2)^12 ; BCOEF = 2 * sqrt(e1*e2) * (r1 + r2)^6 = 2 * ACOEF/(r1+r2)^6
            # to convert ACOEF and BCOEF to r and epsilon
            # r = 0.5 * (2*ACOEF/BCOEF)^(1/6); ep = BCOEF^2/(4*ACOEF)
            # to convert ACOEF and BCOEF to sigma and epsilon (GMX)
            # sigma = (ACOEF/BCOEF)^(1/6) * 0.1 ; epsilon = 4.184 * BCOEF^2/(4*ACOEF)
            #   ca   ca       819971.66        531.10
            #   ca   ha        76245.15        104.66
            #   ha   ha         5716.30         18.52

        For proper dihedrals: a quartet of atoms may appear with more than
        one set of parameters and to convert to GMX they are treated as RBs;
        use the algorithm:

        .. code-block:: c++

            for(my $j=$i;$j<=$lines;$j++){
                my $period = $pn{$j};
                if($pk{$j}>0) {
                $V[$period] = 2*$pk{$j}*$cal;
                }
                # assign V values to C values as predefined #
                if($period==1){
                $C[0]+=0.5*$V[$period];
                if($phase{$j}==0){
                    $C[1]-=0.5*$V[$period];
                }else{
                    $C[1]+=0.5*$V[$period];
                }
                }elsif($period==2){
                if(($phase{$j}==180)||($phase{$j}==3.14)){
                    $C[0]+=$V[$period];
                    $C[2]-=$V[$period];
                }else{
                    $C[2]+=$V[$period];
                }
                }elsif($period==3){
                $C[0]+=0.5*$V[$period];
                if($phase{$j}==0){
                    $C[1]+=1.5*$V[$period];
                    $C[3]-=2*$V[$period];
                }else{
                    $C[1]-=1.5*$V[$period];
                    $C[3]+=2*$V[$period];
                }
                }elsif($period==4){
                if(($phase{$j}==180)||($phase{$j}==3.14)){
                    $C[2]+=4*$V[$period];
                    $C[4]-=4*$V[$period];
                }else{
                    $C[0]+=$V[$period];
                    $C[2]-=4*$V[$period];
                    $C[4]+=4*$V[$period];
                }
                }
            }
        """
        if self.amb2gmx:
            os.chdir(self.absHomeDir)

        self.printMess("Writing GROMACS files\n")

        self.setAtomType4Gromacs()

        self.writeGroFile()

        self.writePosreFile()

        self.writeGromacsTop()

        self.writeMdpFiles()

        if self.amb2gmx:
            os.chdir(self.rootDir)

    def setAtomType4Gromacs(self):
        """Set atom types for Gromacs.

        Atom types names in Gromacs TOP file are not case sensitive;
        this routine will append a '_' to lower case atom type.

        Example:
            >>> CA and ca -> CA and ca_
        """
        if self.merge:
            self.printMess("Merging identical lower and uppercase atomtypes in GMX top file.\n")
            atNames = [at.atomTypeName for at in self.atomTypes]
            delAtomTypes = []
            modAtomTypes = []
            atomTypesGromacs = []
            dictAtomTypes = {}
            for at in self.atomTypes:
                atName = at.atomTypeName
                dictAtomTypes[atName] = at
                if atName.islower() and atName.upper() in atNames:
                    atUpper = self.atomTypes[atNames.index(atName.upper())]
                    if at.ACOEF == atUpper.ACOEF and at.BCOEF == atUpper.BCOEF and at.mass == atUpper.mass:
                        delAtomTypes.append(atName)
                    else:
                        newAtName = atName + "_"
                        modAtomTypes.append(atName)
                        atomType = AtomType(newAtName, at.mass, at.ACOEF, at.BCOEF)
                        atomTypesGromacs.append(atomType)
                        dictAtomTypes[newAtName] = atomType
                else:
                    atomTypesGromacs.append(at)

            atomsGromacs = []
            for a in self.atoms:
                atName = a.atomType.atomTypeName
                if atName in delAtomTypes:
                    atom = Atom(
                        a.atomName,
                        dictAtomTypes[atName.upper()],
                        a.id,
                        a.resid,
                        a.mass,
                        a.charge,
                        a.coords,
                    )
                    atom.cgnr = a.cgnr
                    atomsGromacs.append(atom)
                elif atName in modAtomTypes:
                    atom = Atom(
                        a.atomName,
                        dictAtomTypes[atName + "_"],
                        a.id,
                        a.resid,
                        a.mass,
                        a.charge,
                        a.coords,
                    )
                    atom.cgnr = a.cgnr
                    atomsGromacs.append(atom)
                else:
                    atomsGromacs.append(a)

            self.atomTypesGromacs = atomTypesGromacs
            self.atomsGromacs = atomsGromacs
            return

        self.printMess("Disambiguating lower and uppercase atomtypes in GMX top file, even if identical.\n")
        self.atomTypesGromacs = self.atomTypes
        self.atomsGromacs = self.atoms

    def shiftLineIds(self, line, nIds, offset):
        """Renumber the leading atom-id columns of a topology line, keeping its width."""
        out = ""
        rest = line
        for _ in range(nIds):
            match = re.match(r"(\s*)(\d+)", rest)
            if not match:
                break
            out += "%*d" % (len(match.group(0)), int(match.group(2)) - offset)
            rest = rest[match.end() :]
        return out + rest

    def gmxMoleculeName(self, lo, hi, used):
        """Name one solute molecule, preferring its residue name when it has just one."""
        resids = {self.atoms[id_].resid for id_ in range(lo, hi)}
        if len(resids) == 1:
            name = self.residueLabel[resids.pop()]
        else:
            name = "%s_%i" % (self.baseName, len(used) + 1)
        # GROMACS matches moleculetypes by name, so two blocks may never share one.
        # Identical molecules are written out separately rather than collapsed into an
        # nmols count: mistaking two different molecules for one corrupts the system,
        # while repeating an identical one only costs a few lines.
        stem, count = name, 2
        while name in used:
            name = "%s_%i" % (stem, count)
            count += 1
        used.append(name)
        return name

    def gmxSoluteBlocks(self, headMoleculetype, sections):
        """Render the solute as one ``[ moleculetype ]`` per AMBER molecule.

        Args:
            headMoleculetype: the ``[ moleculetype ]`` template, taking a name.
            sections: ``(header, lines, nIds)`` triples in the order they were built,
                where ``nIds`` is how many leading columns of each line are atom ids.

        Returns:
            list: the lines to write, or ``None`` when the solute cannot be split and
            the caller should fall back to a single block.
        """
        atomLines = sections[0][1]
        ranges = self.amberMoleculeRanges(len(atomLines))
        if len(ranges) < 2:
            return None

        text = []
        names = []
        for lo, hi in ranges:
            text.append(headMoleculetype % self.gmxMoleculeName(lo, hi, names))
            for header, lines, nIds in sections:
                # A bonded term never spans two molecules -- that is what makes them
                # separate molecules -- so the first id decides where the line belongs.
                picked = [line for line in lines if lo < int(line.split()[0]) <= hi]
                if not picked:
                    continue
                text.append(header)
                if nIds == 1:
                    qtot = 0.0
                    for line in picked:
                        qtot += float(line.split()[6])
                        shifted = self.shiftLineIds(line, nIds, lo)
                        text.append(re.sub(r"(; qtot )\S+", "\\g<1>%1.3f" % qtot, shifted))
                else:
                    text += [self.shiftLineIds(line, nIds, lo) for line in picked]
        self.printMess("Split the solute into %i moleculetypes: %s\n" % (len(names), ", ".join(names)))
        self.gmxSplitNames = names
        return text

    def writeGromacsTop(self):
        """Write GMX topology file."""
        if self.atomTypeSystem == "amber":
            d2opls = dictAtomTypeAmb2OplsGmxCode
        else:
            d2opls = dictAtomTypeGaff2OplsGmxCode

        topText = []
        itpText = []
        oitpText = []
        otopText = []
        top = self.baseName + "_GMX.top"
        itp = self.baseName + "_GMX.itp"
        posre = "posre_" + self.baseName + ".itp"
        otop = self.baseName + "_GMX_OPLS.top"
        oitp = self.baseName + "_GMX_OPLS.itp"

        headDefault = """
[ defaults ]
; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ
1               2               yes             0.5     0.8333333333
"""
        headItp = """
; Include %s topology
#include "%s"
"""
        headLigPosre = """
; Ligand position restraints
#ifdef POSRES_LIG
#include "%s"
#endif
"""
        headOpls = """
; Include forcefield parameters
#include "ffoplsaa.itp"
"""
        headSystem = """
[ system ]
 %s
"""
        headMols = """
[ molecules ]
; Compound        nmols
"""
        headAtomtypes = """
[ atomtypes ]
;name   bond_type     mass     charge   ptype   sigma         epsilon       Amb
"""
        headAtomtypesOpls = """
; For OPLS atomtypes manual fine tuning
; AC_at:OPLS_at:OPLScode: Possible_Alternatives (see ffoplsaa.atp and ffoplsaanb.itp)
"""
        headMoleculetype = """
[ moleculetype ]
;name            nrexcl
 %-16s 3
"""
        headAtoms = """
[ atoms ]
;   nr  type  resi  res  atom  cgnr     charge      mass       ; qtot   bond_type
"""
        headBonds = """
[ bonds ]
;   ai     aj funct   r             k
"""
        headPairs = """
[ pairs ]
;   ai     aj    funct
"""
        headAngles = """
[ angles ]
;   ai     aj     ak    funct   theta         cth
"""
        headProDih = """
[ dihedrals ] ; propers
; treated as RBs in GROMACS to use combine multiple AMBER torsions per quartet
;    i      j      k      l   func    C0         C1         C2         C3         C4         C5
"""

        headProDihAlphaGamma = """; treated as usual propers in GROMACS since Phase angle diff from 0 or 180 degrees
;    i      j      k      l   func   phase     kd      pn
"""

        headProDihGmx45 = """
[ dihedrals ] ; propers
; for gromacs 4.5 or higher, using funct 9
;    i      j      k      l   func   phase     kd      pn
"""

        headImpDih = """
[ dihedrals ] ; impropers
; treated as propers in GROMACS to use correct AMBER analytical function
;    i      j      k      l   func   phase     kd      pn
"""

        # NOTE: headTopWaterTip3p and headTopWaterSpce actually do NOTHING
        # ==============================================================================================================
        #         _headTopWaterTip3p = """
        # [ bondtypes ]
        #   ; i    j      func       b0          kb
        #   OW    HW         1    0.09572   462750.4 ; TIP3P water
        #   HW    HW         1    0.15139   462750.4 ; TIP3P water
        #
        # [ angletypes ]
        #   ;  i    j    k  func       th0       cth
        #   HW  OW  HW           1   104.520    836.800 ; TIP3P water
        #   HW  HW  OW           1   127.740      0.000 ; (found in crystallographic water with 3 bonds)
        # """
        #
        #         _headTopWaterSpce = """
        # [ bondtypes ]
        #   ; i    j      func       b0          kb
        #   OW    HW         1    0.1       462750.4 ; SPCE water
        #   HW    HW         1    0.1633    462750.4 ; SPCE water
        #
        # [ angletypes ]
        #   ;  i    j    k  func       th0       cth
        #   HW  OW  HW           1   109.47      836.800 ; SPCE water
        #   HW  HW  OW           1   125.265     0.000 ; SPCE water
        # """
        # ==============================================================================================================

        headNa = """
[ moleculetype ]
  ; molname       nrexcl
  NA+             1

[ atoms ]
  ; id_    at type res nr  residue name     at name  cg nr  charge   mass
    1       %s      1          NA+         NA+       1      1     22.9898
"""
        headCl = """
[ moleculetype ]
  ; molname       nrexcl
  CL-             1

[ atoms ]
  ; id_    at type res nr  residue name     at name  cg nr  charge   mass
    1       %s      1         CL-           CL-      1     -1     35.45300
"""
        headK = """
[ moleculetype ]
  ; molname       nrexcl
  K+             1

[ atoms ]
  ; id_    at type res nr  residue name     at name  cg nr  charge   mass
    1       %s       1          K+         K+       1      1     39.100
"""
        headWaterTip3p = """
[ moleculetype ]
; molname       nrexcl ; TIP3P model
  WAT             2

[ atoms ]
;   nr   type  resnr residue  atom   cgnr     charge       mass
     1     OW      1     WAT     O      1     -0.834   16.00000
     2     HW      1     WAT    H1      1      0.417    1.00800
     3     HW      1     WAT    H2      1      0.417    1.00800

#ifdef FLEXIBLE
[ bonds ]
; i j   funct   length  force.c.
1   2   1   0.09572   462750.4 0.09572   462750.4
1   3   1   0.09572   462750.4 0.09572   462750.4

[ angles ]
; i j   k   funct   angle   force.c.
2   1   3   1   104.520    836.800  104.520    836.800
#else
[ settles ]
; i j   funct   length
1   1   0.09572 0.15139

[ exclusions ]
1   2   3
2   1   3
3   1   2
#endif
"""

        headWaterSpce = """
[ moleculetype ]
; molname       nrexcl ; SPCE model
  WAT             2

[ atoms ]
;   nr   type  resnr residue  atom   cgnr     charge       mass
     1     OW      1     WAT     O      1     -0.8476  15.99940
     2     HW      1     WAT    H1      1      0.4238   1.00800
     3     HW      1     WAT    H2      1      0.4238   1.00800

#ifdef FLEXIBLE
[ bonds ]
; i j   funct   length  force.c.
1   2   1   0.1 462750.4  0.1     462750.4
1   3   1   0.1 462750.4  0.1     462750.4

[ angles ]
; i j   k   funct   angle   force.c.
2   1   3   1   109.47  836.800 109.47  836.800
#else
[ settles ]
; OW    funct   doh dhh
1   1   0.1 0.16330

[ exclusions ]
1   2   3
2   1   3
3   1   2
#endif
"""
        if self.direct and self.amb2gmx:
            self.printMess("Converting directly from AMBER to GROMACS (EXPERIMENTAL).\n")

        # Dict of ions dealt by acpype emulating amb2gmx
        ionsDict = {"Na+": headNa, "Cl-": headCl, "K+": headK}
        ionsSorted = []
        # NOTE: headWaterTip3p and headWaterSpce actually do the real thing
        #      so, skipping headTopWaterTip3p and headWaterTip3p
        # headTopWater = headTopWaterTip3p
        headWater = headWaterTip3p

        nWat = 0
        topText.append("; " + head % (top, date))
        otopText.append("; " + head % (otop, date))
        topText.append(headDefault)

        nSolute = 0
        # Section lines are collected rather than written straight out, so that a split
        # solute can partition them per molecule without duplicating any formatting.
        soluteSections = []
        if not self.amb2gmx:
            topText.append(headItp % (itp, itp))
            topText.append(headLigPosre % posre)
            otopText.append(headOpls)
            otopText.append(headItp % (itp, itp))
            otopText.append(headLigPosre % posre)
            itpText.append("; " + head % (itp, date))
            oitpText.append("; " + head % (oitp, date))

        self.printDebug("atomTypes %i" % len(self.atomTypesGromacs))
        temp = []
        otemp = []
        for aType in self.atomTypesGromacs:
            aTypeName = aType.atomTypeName
            oaCode = d2opls.get(aTypeName, ["x", "0"])[:-1]
            aTypeNameOpls = oplsCode2AtomTypeDict.get(oaCode[0], "x")
            A = aType.ACOEF
            B = aType.BCOEF
            # one cannot infer sigma or epsilon for B = 0, assuming 0 for them
            if B == 0.0:
                sigma, epsilon, r0, epAmber = 0, 0, 0, 0
            else:
                r0 = 0.5 * math.pow((2 * A / B), (1.0 / 6))
                epAmber = 0.25 * B * B / A
                sigma = 0.1 * math.pow((A / B), (1.0 / 6))
                epsilon = cal * epAmber
            if aTypeName == "OW":
                if A == 629362.166 and B == 625.267765:
                    # headTopWater = headTopWaterSpce
                    headWater = headWaterSpce
            # OW 629362.166 625.267765 spce
            # OW 581935.564 594.825035 tip3p
            #       print aTypeName, A, B
            line = (
                " %-8s %-11s %3.5f  %3.5f   A   %13.5e %13.5e"
                % (
                    aTypeName,
                    aTypeName,
                    aType.mass,
                    0.0,
                    sigma,
                    epsilon,
                )
                + f" ; {r0:4.2f}  {epAmber:1.4f}\n"
            )
            oline = f"; {aTypeName}:{aTypeNameOpls}:opls_{oaCode[0]}: {oaCode[1:]!r}\n"
            # tmpFile.write(line)
            temp.append(line)
            otemp.append(oline)
        if self.amb2gmx:
            topText.append(headAtomtypes)
            topText += temp
            nWat = self.residueLabel.count("WAT")
            for ion in ionsDict:
                nIon = self.residueLabel.count(ion)
                if nIon > 0:
                    idIon = self.residueLabel.index(ion)
                    ionType = self.search(name=ion).atomType.atomTypeName
                    ionsSorted.append((idIon, nIon, ion, ionType))
            ionsSorted.sort()
        else:
            itpText.append(headAtomtypes)
            itpText += temp
            oitpText.append(headAtomtypesOpls)
            oitpText += otemp
        self.printDebug("GMX atomtypes done")

        if len(self.atoms) > 3 * nWat + sum(x[1] for x in ionsSorted):
            nSolute = 1

        if nWat:
            # topText.append(headTopWater)
            self.printDebug("type of water '%s'" % headWater[43:48].strip())

        if nSolute:
            if not self.amb2gmx:
                itpText.append(headMoleculetype % self.baseName)
                oitpText.append(headMoleculetype % self.baseName)

        self.printDebug("atoms %i" % len(self.atoms))
        qtot = 0.0
        count = 1
        temp = []
        otemp = []
        id2oplsATDict = {}
        for atom in self.atomsGromacs:
            resid = atom.resid
            resname = self.residueLabel[resid]
            if not self.direct:
                if resname in [*list(ionsDict), "WAT"]:
                    break
            aName = atom.atomName
            aType = atom.atomType.atomTypeName
            oItem = d2opls.get(aType, ["x", 0])
            oplsAtName = oplsCode2AtomTypeDict.get(oItem[0], "x")
            id_ = atom.id
            id2oplsATDict[id_] = oplsAtName
            oaCode = "opls_" + str(oItem[0])
            cgnr = id_
            if self.sorted:
                cgnr = atom.cgnr  # JDC
            charge = atom.charge
            mass = atom.mass
            omass = float(oItem[-1])
            qtot += charge
            resnr = resid + 1
            line = "%6d %4s %5d %5s %5s %4d %12.6f %12.5f ; qtot %1.3f\n" % (
                id_,
                aType,
                resnr,
                resname,
                aName,
                cgnr,
                charge,
                mass,
                qtot,
            )  # JDC
            oline = "%6d %4s %5d %5s %5s %4d %12.6f %12.5f ; qtot % 3.3f  %-4s\n" % (
                id_,
                oaCode,
                resnr,
                resname,
                aName,
                cgnr,
                charge,
                omass,
                qtot,
                oplsAtName,
            )  # JDC
            count += 1
            temp.append(line)
            otemp.append(oline)
        if temp:
            if self.amb2gmx:
                soluteSections.append((headAtoms, temp, 1))
            else:
                itpText.append(headAtoms)
                itpText += temp
                oitpText.append(headAtoms)
                oitpText += otemp
        self.printDebug("GMX atoms done")

        # remove bond of water
        self.printDebug("bonds %i" % len(self.bonds))
        temp = []
        otemp = []
        for bond in self.bonds:
            res1 = self.residueLabel[bond.atoms[0].resid]
            res2 = self.residueLabel[bond.atoms[0].resid]
            if "WAT" in [res1, res2]:
                continue
            a1Name = bond.atoms[0].atomName
            a2Name = bond.atoms[1].atomName
            id1 = bond.atoms[0].id
            id2 = bond.atoms[1].id
            oat1 = id2oplsATDict.get(id1)
            oat2 = id2oplsATDict.get(id2)
            line = "%6i %6i %3i %13.4e %13.4e ; %6s - %-6s\n" % (
                id1,
                id2,
                1,
                bond.rEq * 0.1,
                bond.kBond * 200 * cal,
                a1Name,
                a2Name,
            )
            oline = "%6i %6i %3i ; %13.4e %13.4e ; %6s - %-6s %6s - %-6s\n" % (
                id1,
                id2,
                1,
                bond.rEq * 0.1,
                bond.kBond * 200 * cal,
                a1Name,
                a2Name,
                oat1,
                oat2,
            )
            temp.append(line)
            otemp.append(oline)
        temp.sort()
        otemp.sort()
        if temp:
            if self.amb2gmx:
                soluteSections.append((headBonds, temp, 2))
            else:
                itpText.append(headBonds)
                itpText += temp
                oitpText.append(headBonds)
                oitpText += otemp
        self.printDebug("GMX bonds done")

        self.printDebug("atomPairs %i" % len(self.atomPairs))
        temp = []
        for pair in self.atomPairs:
            # if not printed:
            #    tmpFile.write(headPairs)
            #    printed = True
            a1Name = pair[0].atomName
            a2Name = pair[1].atomName
            id1 = pair[0].id
            id2 = pair[1].id
            # id1 = self.atoms.index(pair[0]) + 1
            # id2 = self.atoms.index(pair[1]) + 1
            line = "%6i %6i %6i ; %6s - %-6s\n" % (id1, id2, 1, a1Name, a2Name)
            temp.append(line)
        temp.sort()
        if temp:
            if self.amb2gmx:
                soluteSections.append((headPairs, temp, 2))
            else:
                itpText.append(headPairs)
                itpText += temp
                oitpText.append(headPairs)
                oitpText += temp
        self.printDebug("GMX pairs done")

        self.printDebug("angles %i" % len(self.angles))
        temp = []
        otemp = []
        for angle in self.angles:
            a1 = angle.atoms[0].atomName
            a2 = angle.atoms[1].atomName
            a3 = angle.atoms[2].atomName
            id1 = angle.atoms[0].id
            id2 = angle.atoms[1].id
            id3 = angle.atoms[2].id
            oat1 = id2oplsATDict.get(id1)
            oat2 = id2oplsATDict.get(id2)
            oat3 = id2oplsATDict.get(id3)
            line = "%6i %6i %6i %6i %13.4e %13.4e ; %6s - %-6s - %-6s\n" % (
                id1,
                id2,
                id3,
                1,
                angle.thetaEq * radPi,
                2 * cal * angle.kTheta,
                a1,
                a2,
                a3,
            )
            oline = "%6i %6i %6i %6i ; %13.4e %13.4e ; %6s - %-4s - %-6s %4s - %+4s - %-4s\n" % (
                id1,
                id2,
                id3,
                1,
                angle.thetaEq * radPi,
                2 * cal * angle.kTheta,
                a1,
                a2,
                a3,
                oat1,
                oat2,
                oat3,
            )
            temp.append(line)
            otemp.append(oline)
        temp.sort()
        otemp.sort()
        if temp:
            if self.amb2gmx:
                soluteSections.append((headAngles, temp, 3))
            else:
                itpText.append(headAngles)
                itpText += temp
                oitpText.append(headAngles)
                oitpText += otemp
        self.printDebug("GMX angles done")

        self.setProperDihedralsCoef()
        self.printDebug("properDihedralsCoefRB %i" % len(self.properDihedralsCoefRB))
        self.printDebug("properDihedralsAlphaGamma %i" % len(self.properDihedralsAlphaGamma))
        self.printDebug("properDihedralsGmx45 %i" % len(self.properDihedralsGmx45))
        temp = []
        otemp = []
        if self.gmx4:
            self.printMess("Writing RB dihedrals for old GMX 4.\n")
            for dih in self.properDihedralsCoefRB:
                a1 = dih[0][0].atomName
                a2 = dih[0][1].atomName
                a3 = dih[0][2].atomName
                a4 = dih[0][3].atomName
                id1 = dih[0][0].id
                id2 = dih[0][1].id
                id3 = dih[0][2].id
                id4 = dih[0][3].id
                oat1 = id2oplsATDict.get(id1)
                oat2 = id2oplsATDict.get(id2)
                oat3 = id2oplsATDict.get(id3)
                oat4 = id2oplsATDict.get(id4)
                c0, c1, c2, c3, c4, c5 = dih[1]
                line = "%6i %6i %6i %6i %6i %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f" % (
                    id1,
                    id2,
                    id3,
                    id4,
                    3,
                    c0,
                    c1,
                    c2,
                    c3,
                    c4,
                    c5,
                ) + " ; %6s-%6s-%6s-%6s\n" % (a1, a2, a3, a4)
                oline = "%6i %6i %6i %6i %6i ; %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f" % (
                    id1,
                    id2,
                    id3,
                    id4,
                    3,
                    c0,
                    c1,
                    c2,
                    c3,
                    c4,
                    c5,
                ) + " ; %6s-%6s-%6s-%6s    %4s-%4s-%4s-%4s\n" % (a1, a2, a3, a4, oat1, oat2, oat3, oat4)
                temp.append(line)
                otemp.append(oline)
            temp.sort()
            otemp.sort()
            if temp:
                if self.amb2gmx:
                    soluteSections.append((headProDih, temp, 4))
                else:
                    itpText.append(headProDih)
                    itpText += temp
                    oitpText.append(headProDih)
                    oitpText += otemp
            self.printDebug("GMX proper dihedrals done")
        else:
            self.printMess("Writing GMX dihedrals for GMX 4.5 and higher.\n")
            funct = 9  # 9
            for dih in self.properDihedralsGmx45:
                a1 = dih[0][0].atomName
                a2 = dih[0][1].atomName
                a3 = dih[0][2].atomName
                a4 = dih[0][3].atomName
                id1 = dih[0][0].id
                id2 = dih[0][1].id
                id3 = dih[0][2].id
                id4 = dih[0][3].id
                ph = dih[1]  # phase already in degree
                kd = dih[2] * cal  # kPhi PK
                pn = dih[3]  # .period
                line = "%6i %6i %6i %6i %6i %8.2f %9.5f %3i ; %6s-%6s-%6s-%6s\n" % (
                    id1,
                    id2,
                    id3,
                    id4,
                    funct,
                    ph,
                    kd,
                    pn,
                    a1,
                    a2,
                    a3,
                    a4,
                )
                oline = "%6i %6i %6i %6i %6i ; %8.2f %9.5f %3i ; %6s-%6s-%6s-%6s\n" % (
                    id1,
                    id2,
                    id3,
                    id4,
                    funct,
                    ph,
                    kd,
                    pn,
                    a1,
                    a2,
                    a3,
                    a4,
                )
                temp.append(line)
                otemp.append(oline)
            temp.sort()
            otemp.sort()
            if temp:
                if self.amb2gmx:
                    soluteSections.append((headProDihGmx45, temp, 4))
                else:
                    itpText.append(headProDihGmx45)
                    itpText += temp
                    oitpText.append(headProDihGmx45)
                    oitpText += otemp

        # for properDihedralsAlphaGamma
        if not self.gmx4:
            funct = 4  # 4
        else:
            funct = 1
        temp = []
        otemp = []
        for dih in self.properDihedralsAlphaGamma:
            a1 = dih[0][0].atomName
            a2 = dih[0][1].atomName
            a3 = dih[0][2].atomName
            a4 = dih[0][3].atomName
            id1 = dih[0][0].id
            id2 = dih[0][1].id
            id3 = dih[0][2].id
            id4 = dih[0][3].id
            ph = dih[1]  # phase already in degree
            kd = dih[2] * cal  # kPhi PK
            pn = dih[3]  # .period
            line = "%6i %6i %6i %6i %6i %8.2f %9.5f %3i ; %6s-%6s-%6s-%6s\n" % (
                id1,
                id2,
                id3,
                id4,
                funct,
                ph,
                kd,
                pn,
                a1,
                a2,
                a3,
                a4,
            )
            oline = "%6i %6i %6i %6i %6i ; %8.2f %9.5f %3i ; %6s-%6s-%6s-%6s\n" % (
                id1,
                id2,
                id3,
                id4,
                funct,
                ph,
                kd,
                pn,
                a1,
                a2,
                a3,
                a4,
            )
            temp.append(line)
            otemp.append(oline)
        temp.sort()
        otemp.sort()
        if temp:
            if self.amb2gmx:
                soluteSections.append((headProDihAlphaGamma, temp, 4))
            else:
                itpText.append(headProDihAlphaGamma)
                itpText += temp
                oitpText.append(headProDihAlphaGamma)
                oitpText += otemp
        self.printDebug("GMX special proper dihedrals done")

        self.printDebug("improperDihedrals %i" % len(self.improperDihedrals))
        temp = []
        otemp = []
        for dih in self.improperDihedrals:
            a1 = dih.atoms[0].atomName
            a2 = dih.atoms[1].atomName
            a3 = dih.atoms[2].atomName
            a4 = dih.atoms[3].atomName
            id1 = dih.atoms[0].id
            id2 = dih.atoms[1].id
            id3 = dih.atoms[2].id
            id4 = dih.atoms[3].id
            kd = dih.kPhi * cal
            pn = dih.period
            ph = dih.phase * radPi
            line = "%6i %6i %6i %6i %6i %8.2f %9.5f %3i ; %6s-%6s-%6s-%6s\n" % (
                id1,
                id2,
                id3,
                id4,
                funct,
                ph,
                kd,
                pn,
                a1,
                a2,
                a3,
                a4,
            )
            oline = "%6i %6i %6i %6i %6i ; %8.2f %9.5f %3i ; %6s-%6s-%6s-%6s\n" % (
                id1,
                id2,
                id3,
                id4,
                funct,
                ph,
                kd,
                pn,
                a1,
                a2,
                a3,
                a4,
            )
            temp.append(line)
            otemp.append(oline)
        temp.sort()
        otemp.sort()
        if temp:
            if self.amb2gmx:
                soluteSections.append((headImpDih, temp, 4))
            else:
                itpText.append(headImpDih)
                itpText += temp
                oitpText.append(headImpDih)
                oitpText += otemp
        self.printDebug("GMX improper dihedrals done")

        if self.amb2gmx and nSolute:
            blocks = self.gmxSoluteBlocks(headMoleculetype, soluteSections) if self.splitMolecules else None
            if blocks is None:
                blocks = [headMoleculetype % self.baseName]
                for header, lines, _ in soluteSections:
                    blocks.append(header)
                    blocks += lines
            topText += blocks

        if not self.direct:
            for ion in ionsSorted:
                topText.append(ionsDict[ion[2]] % ion[3])

            if nWat:
                topText.append(headWater)

        topText.append(headSystem % (self.baseName))
        topText.append(headMols)
        otopText.append(headSystem % (self.baseName))
        otopText.append(headMols)

        if nSolute > 0:
            for name in self.gmxSplitNames or [self.baseName]:
                topText.append(" %-16s %-6i\n" % (name, nSolute))
            otopText.append(" %-16s %-6i\n" % (self.baseName, nSolute))

        if not self.direct:
            for ion in ionsSorted:
                topText.append(" %-16s %-6i\n" % (ion[2].upper(), ion[1]))

            if nWat:
                topText.append(" %-16s %-6i\n" % ("WAT", nWat))

        if self.topo14Data.hasNondefault14():
            citation = (
                "     BERNARDI, A., FALLER, R., REITH, D., and KIRSCHNER, K. N. ACPYPE update for\n"
                + "     nonuniform 1-4 scale factors: Conversion of the GLYCAM06 force field from AMBER\n"
                + '     to GROMACS. SoftwareX 10 (2019), 100241. doi: 10.1016/j.softx.2019.100241"\n'
            )

            msg = "Non-default 1-4 scale parameters detected.  Converting individually. Please cite:\n\n" + citation

            self.printMess(msg)
            topText = self.topo14Data.patch_gmx_topol14("".join(topText))

        gmxDir = os.path.abspath(".")
        topFileName = os.path.join(gmxDir, top)
        topFile = open(topFileName, "w")
        topFile.writelines(topText)
        self.topText = topText

        if not self.amb2gmx:
            itpFileName = os.path.join(gmxDir, itp)
            itpFile = open(itpFileName, "w")
            itpFile.writelines(itpText)
            oitpFileName = os.path.join(gmxDir, oitp)
            oitpFile = open(oitpFileName, "w")
            oitpFile.writelines(oitpText)
            otopFileName = os.path.join(gmxDir, otop)
            otopFile = open(otopFileName, "w")
            otopFile.writelines(otopText)

    def writeGroFile(self):
        """Write GRO files."""
        # print "Writing GROMACS GRO file\n"
        self.printDebug("writing GRO file")
        gro = self.baseName + "_GMX.gro"
        gmxDir = os.path.abspath(".")
        groFileName = os.path.join(gmxDir, gro)
        groFile = open(groFileName, "w")
        groFile.write(head % (gro, date))
        groFile.write(" %i\n" % len(self.atoms))
        count = 1
        for atom in self.atoms:
            coords = [c * 0.1 for c in atom.coords]
            resid = atom.resid
            line = "%5d%5s%5s%5d%8.3f%8.3f%8.3f\n" % (
                resid + 1,
                self.residueLabel[resid],
                atom.atomName,
                count,
                coords[0],
                coords[1],
                coords[2],
            )
            count += 1
            if count == 100000:
                count = 0
            groFile.write(line)
        if self.pbc:
            boxX = self.pbc[0][0] * 0.1
            boxY = self.pbc[0][1] * 0.1
            boxZ = self.pbc[0][2] * 0.1
            vX = self.pbc[1][0]
            # vY = self.pbc[1][1]
            # vZ = self.pbc[1][2]
            if vX == 90.0:
                self.printDebug("PBC triclinic")
                text = f"{boxX:11.5f} {boxY:11.5f} {boxZ:11.5f}\n"
            elif round(vX, 2) == 109.47:
                self.printDebug("PBC octahedron")
                f1 = 0.471405  # 1/3 * sqrt(2)
                f2 = 0.333333 * boxX
                v22 = boxY * 2 * f1
                v33 = boxZ * f1 * 1.73205  # f1 * sqrt(3)
                v21 = v31 = v32 = 0.0
                v12 = f2
                v13 = -f2
                v23 = f1 * boxX
                text = f"{boxX:11.5f} {v22:11.5f} {v33:11.5f} {v21:11.5f} {v31:11.5f} {v12:11.5f} {v32:11.5f} {v13:11.5f} {v23:11.5f}\n"
        else:
            self.printDebug("Box size estimated")
            X = [a.coords[0] * 0.1 for a in self.atoms]
            Y = [a.coords[1] * 0.1 for a in self.atoms]
            Z = [a.coords[2] * 0.1 for a in self.atoms]
            boxX = max(X) - min(X)  # + 2.0 # 2.0 is double of rlist
            boxY = max(Y) - min(Y)  # + 2.0
            boxZ = max(Z) - min(Z)  # + 2.0
            text = f"{boxX * 20.0:11.5f} {boxY * 20.0:11.5f} {boxZ * 20.0:11.5f}\n"
        groFile.write(text)

    def writePosreFile(self, fc=1000):
        """Write file with positional restraints for heavy atoms.

        http://www.mdtutorials.com/gmx/complex/06_equil.html
        """
        self.printDebug("writing POSRE file")
        posre = "posre_" + self.baseName + ".itp"
        gmxDir = os.path.abspath(".")
        posreFileName = os.path.join(gmxDir, posre)
        posreFile = open(posreFileName, "w")
        posreFile.write("; " + head % (posre, date))
        posreFile.write("\n[ position_restraints ]\n; atom  type      fx      fy      fz\n")
        for atom in self.atoms:
            if not atom.atomType.atomTypeName.upper().startswith("H"):
                posreFile.write(f"{atom.id:>6d}     1 {fc:>5d} {fc:>5d} {fc:>5d}\n")

    def writeMdpFiles(self):
        """Write MDP for test with GROMACS."""
        emMdp = f"""; to test
; echo 0 | gmx editconf -f {self.baseName}_GMX.gro -bt octahedron -d 10 -c -princ
; gmx grompp -f em.mdp -c out.gro -p {self.baseName}_GMX.top -o em.tpr -v
; gmx mdrun -ntmpi 1 -v -deffnm em

; Parameters describing what to do, when to stop and what to save
integrator      = steep     ; Algorithm (steep = steepest descent minimization)
nsteps          = 500       ; Maximum number of (minimization) steps to perform
nstxout         = 10

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist         = 1             ; Frequency to update the neighbour list and long range forces
cutoff-scheme   = Verlet
rlist           = 1.2           ; Cut-off for making neighbour list (short range forces)
coulombtype     = PME           ; Treatment of long range electrostatic interactions
rcoulomb        = 1.2           ; long range electrostatic cut-off
vdw-type        = cutoff
vdw-modifier    = force-switch
rvdw-switch     = 1.0
rvdw            = 1.2           ; long range Van der Waals cut-off
pbc             = xyz           ; Periodic Boundary Conditions
DispCorr        = no
; vmd em.gro em.trr
"""
        mdMdp = f"""; to test
; gmx grompp -f md.mdp -c em.gro -p {self.baseName}_GMX.top -o md.tpr
; gmx mdrun -ntmpi 1 -v -deffnm md
; define                   = -DPOSRES_LIG
integrator               = md
nsteps                   = 10000
nstxout                  = 10
cutoff-scheme            = verlet
coulombtype              = PME
constraints              = h-bonds
vdwtype                  = cutoff
vdw-modifier             = force-switch
rlist                    = 1.0
rvdw                     = 1.0
rvdw-switch              = 0.9
rcoulomb                 = 1.1
DispCorr                 = EnerPres
lincs-iter               = 2
fourierspacing           = 0.25
gen-vel                  = yes
; vmd md.gro md.trr
"""
        rungmx = f"""
echo 0 | gmx editconf -f {self.baseName}_GMX.gro -bt octahedron -d 10 -c -princ
gmx grompp -f em.mdp -c out.gro -p {self.baseName}_GMX.top -o em.tpr -v
gmx mdrun -ntmpi 1 -v -deffnm em
gmx grompp -f md.mdp -c em.gro -p {self.baseName}_GMX.top -o md.tpr -r em.gro
gmx mdrun -ntmpi 1 -v -deffnm md
"""
        emMdpFile = open("em.mdp", "w")
        mdMdpFile = open("md.mdp", "w")
        runGmxFile = open("rungmx.sh", "w")
        emMdpFile.write(emMdp)
        mdMdpFile.write(mdMdp)
        runGmxFile.write(rungmx)
        os.chmod("rungmx.sh", 0o744)

    def writeCnsTopolFiles(self):
        """Write CNS topology files."""
        if self.amb2gmx:
            os.chdir(self.absHomeDir)

        autoAngleFlag = True
        autoDihFlag = True
        cnsDir = os.path.abspath(".")

        pdb = self.baseName + "_NEW.pdb"
        par = self.baseName + "_CNS.par"
        top = self.baseName + "_CNS.top"
        inp = self.baseName + "_CNS.inp"

        pdbFileName = os.path.join(cnsDir, pdb)
        parFileName = os.path.join(cnsDir, par)
        topFileName = os.path.join(cnsDir, top)
        inpFileName = os.path.join(cnsDir, inp)

        self.CnsTopFileName = topFileName
        self.CnsInpFileName = inpFileName
        self.CnsParFileName = parFileName
        self.CnsPdbFileName = pdbFileName

        parFile = open(parFileName, "w")
        topFile = open(topFileName, "w")
        inpFile = open(inpFileName, "w")

        self.printMess("Writing NEW PDB file\n")
        self.writePdb(pdbFileName)

        self.printMess("Writing CNS/XPLOR files\n")

        # print "Writing CNS PAR file\n"
        parFile.write("Remarks " + head % (par, date))
        parFile.write("\nset echo=false end\n")

        parFile.write("\n{ Bonds: atomType1 atomType2 kb r0 }\n")
        lineSet = []
        for bond in self.bonds:
            a1Type = bond.atoms[0].atomType.atomTypeName + "_"
            a2Type = bond.atoms[1].atomType.atomTypeName + "_"
            kb = 1000.0
            if not self.allhdg:
                kb = bond.kBond
            r0 = bond.rEq
            line = "BOND %5s %5s %8.1f %8.4f\n" % (a1Type, a2Type, kb, r0)
            lineRev = "BOND %5s %5s %8.1f %8.4f\n" % (a2Type, a1Type, kb, r0)
            if line not in lineSet:
                if lineRev not in lineSet:
                    lineSet.append(line)
        for item in lineSet:
            parFile.write(item)

        parFile.write("\n{ Angles: aType1 aType2 aType3 kt t0 }\n")
        lineSet = []
        for angle in self.angles:
            a1 = angle.atoms[0].atomType.atomTypeName + "_"
            a2 = angle.atoms[1].atomType.atomTypeName + "_"
            a3 = angle.atoms[2].atomType.atomTypeName + "_"
            kt = 500.0
            if not self.allhdg:
                kt = angle.kTheta
            t0 = angle.thetaEq * radPi
            line = "ANGLe %5s %5s %5s %8.1f %8.2f\n" % (a1, a2, a3, kt, t0)
            lineRev = "ANGLe %5s %5s %5s %8.1f %8.2f\n" % (a3, a2, a1, kt, t0)
            if line not in lineSet:
                if lineRev not in lineSet:
                    lineSet.append(line)
        for item in lineSet:
            parFile.write(item)

        parFile.write(
            "\n{ Proper Dihedrals: aType1 aType2 aType3 aType4 kt per\
iod phase }\n"
        )
        lineSet = set()
        for item in self.condensedProperDihedrals:
            seq = ""
            id_ = 0
            for dih in item:
                # id_ = item.index(dih)
                ll = len(item)
                a1 = dih.atoms[0].atomType.atomTypeName + "_"
                a2 = dih.atoms[1].atomType.atomTypeName + "_"
                a3 = dih.atoms[2].atomType.atomTypeName + "_"
                a4 = dih.atoms[3].atomType.atomTypeName + "_"
                kp = 750.0
                if not self.allhdg:
                    kp = dih.kPhi
                p = dih.period
                ph = dih.phase * radPi
                if ll > 1:
                    if id_ == 0:
                        line = (
                            "DIHEdral %5s %5s %5s %5s  MULT %1i %7.3f %4i %8\
.2f\n"
                            % (a1, a2, a3, a4, ll, kp, p, ph)
                        )
                    else:
                        line = "%s %7.3f %4i %8.2f\n" % (40 * " ", kp, p, ph)
                else:
                    line = "DIHEdral %5s %5s %5s %5s %15.3f %4i %8.2f\n" % (
                        a1,
                        a2,
                        a3,
                        a4,
                        kp,
                        p,
                        ph,
                    )
                seq += line
                id_ += 1
            lineSet.add(seq)
        for item in lineSet:
            parFile.write(item)

        parFile.write(
            "\n{ Improper Dihedrals: aType1 aType2 aType3 aType4 kt p\
eriod phase }\n"
        )
        lineSet = set()
        for idh in self.improperDihedrals:
            a1 = idh.atoms[0].atomType.atomTypeName + "_"
            a2 = idh.atoms[1].atomType.atomTypeName + "_"
            a3 = idh.atoms[2].atomType.atomTypeName + "_"
            a4 = idh.atoms[3].atomType.atomTypeName + "_"
            kp = 750.0
            if not self.allhdg:
                kp = idh.kPhi
            p = idh.period
            ph = idh.phase * radPi
            line = "IMPRoper %5s %5s %5s %5s %13.1f %4i %8.2f\n" % (
                a1,
                a2,
                a3,
                a4,
                kp,
                p,
                ph,
            )
            lineSet.add(line)
        if self.chiral:
            for idhc in self.chiralGroups:
                _atc, neig, angle = idhc
                a1 = neig[0].atomType.atomTypeName + "_"
                a2 = neig[1].atomType.atomTypeName + "_"
                a3 = neig[2].atomType.atomTypeName + "_"
                a4 = neig[3].atomType.atomTypeName + "_"
                kp = 11000.0
                p = 0
                ph = angle
                line = "IMPRoper %5s %5s %5s %5s %13.1f %4i %8.2f\n" % (
                    a1,
                    a2,
                    a3,
                    a4,
                    kp,
                    p,
                    ph,
                )
                lineSet.add(line)

        for item in lineSet:
            parFile.write(item)

        parFile.write("\n{ Nonbonded: Type Emin sigma; (1-4): Emin/2 sigma }\n")
        for at in self.atomTypes:
            A = at.ACOEF
            B = at.BCOEF
            atName = at.atomTypeName + "_"
            if B == 0.0:
                sigma = epAmber = ep2 = sig2 = 0.0
            else:
                epAmber = 0.25 * B * B / A
                ep2 = epAmber / 2.0
                sigma = math.pow((A / B), (1.0 / 6))
                sig2 = sigma
            line = "NONBonded %5s %11.6f %11.6f %11.6f %11.6f\n" % (
                atName,
                epAmber,
                sigma,
                ep2,
                sig2,
            )
            parFile.write(line)
        parFile.write("\nset echo=true end\n")

        # print "Writing CNS TOP file\n"
        topFile.write("Remarks " + head % (top, date))
        topFile.write("\nset echo=false end\n")
        topFile.write(f"\nautogenerate angles={autoAngleFlag} dihedrals={autoDihFlag} end\n")
        topFile.write("\n{ atomType  mass }\n")

        for at in self.atomTypes:
            atType = at.atomTypeName + "_"
            mass = at.mass
            line = "MASS %-5s %8.3f\n" % (atType, mass)
            topFile.write(line)

        topFile.write("\nRESIdue %s\n" % self.residueLabel[0])
        topFile.write("\nGROUP\n")

        topFile.write("\n{ atomName  atomType  Charge }\n")

        for at in self.atoms:
            atName = at.atomName
            atType = at.atomType.atomTypeName + "_"
            charge = at.charge
            line = "ATOM %-5s TYPE= %-5s CHARGE= %8.4f END\n" % (atName, atType, charge)
            topFile.write(line)

        topFile.write("\n{ Bonds: atomName1  atomName2 }\n")
        for bond in self.bonds:
            a1Name = bond.atoms[0].atomName
            a2Name = bond.atoms[1].atomName
            line = "BOND %-5s %-5s\n" % (a1Name, a2Name)
            topFile.write(line)

        if not autoAngleFlag or 1:  # noqa: SIM222 # generating angles anyway
            topFile.write("\n{ Angles: atomName1 atomName2 atomName3}\n")
            for angle in self.angles:
                a1Name = angle.atoms[0].atomName
                a2Name = angle.atoms[1].atomName
                a3Name = angle.atoms[2].atomName
                line = "ANGLe %-5s %-5s %-5s\n" % (a1Name, a2Name, a3Name)
                topFile.write(line)

        if not autoDihFlag or 1:  # noqa: SIM222 # generating dihedrals anyway
            topFile.write("\n{ Proper Dihedrals: name1 name2 name3 name4 }\n")
            for item in self.condensedProperDihedrals:
                for dih in item:
                    a1Name = dih.atoms[0].atomName
                    a2Name = dih.atoms[1].atomName
                    a3Name = dih.atoms[2].atomName
                    a4Name = dih.atoms[3].atomName
                    line = "DIHEdral %-5s %-5s %-5s %-5s\n" % (
                        a1Name,
                        a2Name,
                        a3Name,
                        a4Name,
                    )
                    break
                topFile.write(line)

        topFile.write("\n{ Improper Dihedrals: aName1 aName2 aName3 aName4 }\n")
        for dih in self.improperDihedrals:
            a1Name = dih.atoms[0].atomName
            a2Name = dih.atoms[1].atomName
            a3Name = dih.atoms[2].atomName
            a4Name = dih.atoms[3].atomName
            line = "IMPRoper %-5s %-5s %-5s %-5s\n" % (a1Name, a2Name, a3Name, a4Name)
            topFile.write(line)

        if self.chiral:
            for idhc in self.chiralGroups:
                _atc, neig, angle = idhc
                a1Name = neig[0].atomName
                a2Name = neig[1].atomName
                a3Name = neig[2].atomName
                a4Name = neig[3].atomName
                line = "IMPRoper %-5s %-5s %-5s %-5s\n" % (
                    a1Name,
                    a2Name,
                    a3Name,
                    a4Name,
                )
                topFile.write(line)

        topFile.write("\nEND {RESIdue %s}\n" % self.residueLabel[0])

        topFile.write("\nset echo=true end\n")

        inpFile.write("Remarks " + head % (inp, date))
        inpData = """
topology
  @%(CNS_top)s
end

parameters
  @%(CNS_par)s
  nbonds
      atom cdie shift eps=1.0  e14fac=0.4   tolerance=0.5
      cutnb=9.0 ctonnb=7.5 ctofnb=8.0
      nbxmod=5 vswitch wmin 1.0
  end
  remark dielectric constant eps set to 1.0
end

flags exclude elec ? end

segment name="    "
  chain
   coordinates @%(NEW_pdb)s
  end
end
coordinates @%(NEW_pdb)s
coord copy end

! Remarks If you want to shake up the coordinates a bit ...
 vector do (x=x+6*(rand()-0.5)) (all)
 vector do (y=y+6*(rand()-0.5)) (all)
 vector do (z=z+6*(rand()-0.5)) (all)
 write coordinates output=%(CNS_ran)s end

! Remarks RMS diff after randomisation and before minimisation
coord rms sele=(known and not hydrogen) end

print threshold=0.02 bonds
print threshold=3.0 angles
print threshold=3.0 dihedrals
print threshold=3.0 impropers

! Remarks Do Powell energy minimisation
minimise powell
  nstep=250 drop=40.0
end

write coordinates output=%(CNS_min)s end
write structure   output=%(CNS_psf)s end

! constraints interaction (not hydro) (not hydro) end

print threshold=0.02 bonds
print threshold=3.0 angles
print threshold=3.0 dihedrals
print threshold=3.0 impropers

flags exclude * include vdw end energy end
distance from=(not hydro) to=(not hydro) cutoff=2.6 end

! Remarks RMS fit after minimisation
coord fit sele=(known and not hydrogen) end

stop
"""
        dictInp = {}
        dictInp["CNS_top"] = top
        dictInp["CNS_par"] = par
        dictInp["NEW_pdb"] = pdb
        dictInp["CNS_min"] = self.baseName + "_NEW_min.pdb"
        dictInp["CNS_psf"] = self.baseName + "_CNS.psf"
        dictInp["CNS_ran"] = self.baseName + "_rand.pdb"
        line = inpData % dictInp
        inpFile.write(line)
        if not self.amb2gmx:
            self.printDebug("chiralGroups %i" % len(self.chiralGroups))
        else:
            os.chdir(self.rootDir)


class ACTopol(AbstractTopol):
    """Class to build the AC topologies (Antechamber AmberTools)."""

    def __init__(
        self,
        inputFile,
        binaries=binaries,
        chargeType="bcc",
        chargeVal=None,
        multiplicity="1",
        predIndex=4,
        atomType="gaff2",
        force=False,
        basename=None,
        debug=False,
        outTopol="all",
        allhdg=False,
        timeTol=MAXTIME,
        qprog="sqm",
        ekFlag=None,
        verbose=True,
        gmx4=False,
        merge=False,
        direct=False,
        split_molecules=False,
        protein_ff="ff14SB",
        is_sorted=False,
        chiral=False,
        amb2gmx=False,
        level=20,
    ):
        super().__init__()
        self.binaries = binaries
        self.amb2gmx = amb2gmx
        self.debug = debug
        self.verbose = verbose
        self.gmx4 = gmx4
        self.merge = merge
        self.direct = direct
        self.splitMolecules = split_molecules
        self.sorted = is_sorted
        self.chiral = chiral

        if not self.verbose:
            level = 100
        if self.debug:
            level = 10
        self.level = level or 20

        self.acExe = find_bin(binaries["ac_bin"])
        if not os.path.exists(self.acExe):
            self.printError(f"no '{binaries['ac_bin']}' executable... aborting! ")
            hint1 = "HINT1: is 'AMBERHOME' environment variable set?"
            hint2 = (
                f"HINT2: is '{binaries['ac_bin']}' in your $PATH?"
                + f"\n    What 'which {binaries['ac_bin']}' in your terminal says?"
                + "\n    'alias' doesn't work for ACPYPE."
            )
            self.printMess(hint1)
            self.printMess(hint2)
            if bundled_amber_dir() is None:
                self.printMess(
                    "HINT3: this ACPYPE installation carries no AmberTools of its own."
                    "\n    Bundled binaries ship only in the Linux x86_64 and macOS Apple Silicon"
                    "\n    wheels. On any other platform, or when installing from the source"
                    "\n    distribution, install AmberTools yourself, for example:"
                    "\n        conda install -c conda-forge ambertools"
                    "\n    or use the conda/docker route described in the README."
                )
            msg = "Missing ANTECHAMBER"
            logger(self.level).error(msg)
            raise Exception(msg)
        self.inputFile = os.path.basename(inputFile)
        self.rootDir = os.path.abspath(".")
        self.absInputFile = os.path.abspath(inputFile)

        if not os.path.exists(self.absInputFile) and not re.search(r"\.mol2$|\.mdl$|\.pdb$", self.inputFile):
            self.smiles = inputFile
            if self.checkSmiles():
                self.is_smiles = True
                if not basename:
                    self.inputFile = "smiles_molecule.mol2"
                else:
                    self.inputFile = f"{basename}.mol2"
                self.absInputFile = os.path.abspath(self.inputFile)
            else:
                self.is_smiles = False
                self.smiles = ""
        elif not os.path.exists(self.absInputFile):
            msg = f"Input file {inputFile} DOES NOT EXIST"
            logger(self.level).error(msg)
            raise Exception(msg)
        baseOriginal, ext = os.path.splitext(self.inputFile)
        base = basename or baseOriginal
        self.baseOriginal = baseOriginal
        self.ext = ext
        self.baseName = base  # name of the input file without ext.
        self.obabelExe = find_bin(binaries["obabel_bin"])
        if not os.path.exists(self.obabelExe):
            if self.ext != ".mol2" and self.ext != ".mdl":
                self.printError(f"no '{binaries['obabel_bin']}' executable; you need it if input is PDB or SMILES")
                self.printError("otherwise use only MOL2 or MDL file as input ... aborting!")
                msg = "Missing OBABEL"
                logger(self.level).error(msg)
                raise Exception(msg)
            else:
                self.printWarn(f"no '{binaries['obabel_bin']}' executable, no PDB file can be used as input!")
        if self.is_smiles:
            self.convertSmilesToMol2()
        self.timeTol = timeTol
        self.printDebug("Max execution time tolerance is %s" % elapsedTime(self.timeTol))
        # ekFlag e.g. (default used by sqm):
        # acpype -i cccc -k "qm_theory='AM1', grms_tol=0.0005, scfconv=1.d-10, ndiis_attempts=700, qmcharge=0"
        if ekFlag == '"None"' or ekFlag is None:
            self.ekFlag = ""
        else:
            self.ekFlag = "-ek %s" % ekFlag
        self.extOld = ext
        self.homeDir = self.baseName + ".acpype"
        self.chargeType = chargeType
        self.chargeVal = chargeVal
        self.multiplicity = multiplicity
        self.predIndex = predIndex
        self.atomType = atomType
        self.gaffDatfile = "gaff.dat"
        if protein_ff not in PROTEIN_FF:
            raise ValueError(f"protein_ff must be one of {sorted(PROTEIN_FF)}, not {protein_ff!r}")
        self.proteinFF = protein_ff
        leapGaffFile = "leaprc.gaff"
        if "2" in self.atomType:
            leapGaffFile = "leaprc.gaff2"
            self.gaffDatfile = "gaff2.dat"
        self.force = force
        self.allhdg = allhdg
        self.tleapExe = which("tleap") or ""
        self.parmchkExe = which("parmchk2") or ""
        acBase = base + "_AC"
        self.acBaseName = acBase
        self.acXyzFileName = acBase + ".inpcrd"
        self.acTopFileName = acBase + ".prmtop"
        self.acFrcmodFileName = acBase + ".frcmod"
        self.tmpDir = os.path.join(self.rootDir, ".acpype_tmp_%s" % os.path.basename(base))
        self.setResNameCheckCoords()
        self.guessCharge()
        acMol2FileName = f"{base}_{chargeType}_{atomType}.mol2"
        self.acMol2FileName = acMol2FileName
        self.charmmBase = "%s_CHARMM" % base
        self.qFlag = qDict[qprog]
        self.outTopols = [outTopol]
        if outTopol == "all":
            self.outTopols = outTopols
        self.acParDict = {
            "base": base,
            "ext": ext[1:],
            "acBase": acBase,
            "acMol2FileName": acMol2FileName,
            "res": self.resName,
            "leapAmberFile": PROTEIN_FF[protein_ff]["leaprc"],
            "baseOrg": self.baseOriginal,
            "leapGaffFile": leapGaffFile,
        }


class MolTopol(AbstractTopol):
    """Class to write topologies and parameters files for several applications.

    https://ambermd.org/FileFormats.php

    Parser, take information in AC xyz and top files and convert to objects.

    Args:
        acFileXyz
        acFileTop

    Returns:
        molTopol obj or None
    """

    def __init__(
        self,
        acTopolObj=None,
        acFileXyz: str = "",
        acFileTop: str = "",
        debug=False,
        basename=None,
        verbose=True,
        gmx4=False,
        merge=False,
        direct=False,
        split_molecules=False,
        is_sorted=False,
        chiral=False,
        amb2gmx=False,
        level=20,
    ):
        super().__init__()
        self.amb2gmx = amb2gmx
        self.chiral = chiral
        self.allhdg = False
        self.debug = debug
        self.level = level
        self.gmx4 = gmx4
        self.merge = merge
        self.direct = direct
        self.splitMolecules = split_molecules
        self.sorted = is_sorted
        self.verbose = verbose
        self.inputFile = acFileTop

        if not self.verbose:
            level = 100
        if self.debug:
            level = 10
        self.level = level

        if acTopolObj:
            if not acFileXyz:
                acFileXyz = acTopolObj.acXyzFileName
            if not acFileTop:
                acFileTop = acTopolObj.acTopFileName
            self._parent = acTopolObj
            self.allhdg = self._parent.allhdg
            self.debug = self._parent.debug
            self.inputFile = self._parent.inputFile
        elif not self.amb2gmx:
            self.amb2gmx = True
        if not os.path.exists(acFileXyz) or not os.path.exists(acFileTop):
            self.printError(f"Files '{acFileXyz}' and/or '{acFileTop}' don't exist")
            self.printError("molTopol object won't be created")

        self.xyzFileData = open(acFileXyz).readlines()
        self.topFileData = [x for x in open(acFileTop).readlines() if not x.startswith("%COMMENT")]
        self.topo14Data = Topology_14()
        self.topo14Data.read_amber_topology("".join(self.topFileData))
        self.printDebug("prmtop and inpcrd files loaded")

        self.getResidueLabel()
        if len(self.residueLabel) > 1:
            self.baseName = basename or os.path.splitext(os.path.basename(acFileTop))[0]  # 'solute'
        else:
            self.baseName = basename or self.residueLabel[0]  # 3 caps letters
        if acTopolObj:
            self.baseName = basename or acTopolObj.baseName
        self.printDebug("basename defined = '%s'" % self.baseName)

        self.getAtoms()

        self.getBonds()

        self.getAngles()

        self.getDihedrals()

        if self.amb2gmx:
            self.rootDir = os.path.abspath(".")
            self.homeDir = f"{self.baseName}.amb2gmx"
            self.makeDir()
        else:
            self.getChirals()

        # Sort atoms for gromacs output. # JDC
        if self.sorted:
            self.printMess("Sorting atoms for gromacs ordering.\n")
            self.sortAtomsForGromacs()
