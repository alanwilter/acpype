import sys, time, json
from acpype_lib.acpype import ACTopol, MAXTIME, while_replace, header, elapsedTime, traceback
import io

em_mdp = io.StringIO()
AC_frcmod = io.StringIO()
AC_inpcrd = io.StringIO()
AC_lib = io.StringIO()
AC_prmtop = io.StringIO()
bcc_gaff_mol2 = io.StringIO()
CHARMM_inp = io.StringIO()
CHARMM_prm = io.StringIO()
CHARMM_rtf = io.StringIO()
CNS_inp = io.StringIO()
CNS_par = io.StringIO()
CNS_top = io.StringIO()
GMX_OPLS_itp = io.StringIO()
GMX_OPLS_top = io.StringIO()
GMX_gro = io.StringIO()
GMX_itp = io.StringIO()
GMX_top = io.StringIO()
NEW_pdb = io.StringIO()
md_mdp = io.StringIO()
sqm_in = io.StringIO()
sqm_out = io.StringIO()

filesInMemory = [(em_mdp, 'em.mdp'),(AC_frcmod,'_AC.frcmod'),(AC_inpcrd,'_AC.inpcrd'),(AC_lib,'_AC.lib'), (AC_prmtop,'_AC.prmtop'),(bcc_gaff_mol2,'_bcc_gaff.mol2'), (CHARMM_inp,'_CHARMM.inp'),(CHARMM_prm,'_CHARMM.prm'),(CHARMM_rtf,'_CHARMM.rtf'), (CNS_inp, '_CNS.inp'),
(CNS_par,'_CNS.par'),(CNS_top,'_CNS.top'),(GMX_OPLS_itp,'_GMX_OPLS.itp'),(GMX_OPLS_top,'_GMX_OPLS.top'),(GMX_gro,'_GMX.gro'),(GMX_itp,'_GMX.itp'),(GMX_top,'_GMX.top'),(NEW_pdb,'_NEW.pdb'), (md_mdp,'md.mdp'), (sqm_in,'sqm.in'),(sqm_out,'sqm.out')]



def clearFileInMemory():
    for files in filesInMemory:
        files[0].seek(0)
        files[0].truncate(0)

def readFiles(basename):
    for files in filesInMemory:
        if files[1] == 'em.mdp' or files[1] == 'md.mdp' or files[1] =='sqm.in' or files[1] == 'sqm.out':
            filename = files[1]
        else:
            filename = basename+files[1]
        readfile = tuple(open(filename, 'r'))
        for line in readfile:
            files[0].write(line)

def acpype_api(
    inputFile,
    chargeType="bcc",
    chargeVal=None,
    multiplicity="1",
    atomType="gaff",
    force=False,
    basename=None,
    debug=False,
    outTopol="all",
    engine="tleap",
    allhdg=False,
    timeTol=MAXTIME,
    qprog="sqm",
    ekFlag=None,
    verbose=True,
    gmx4=False,
    disam=False,
    direct=False,
    is_sorted=False,
    chiral=False,
    is_smiles = False):

    at0 = time.time()
    print(header)

    if debug:
        texta = "Python Version %s" % sys.version
        print("DEBUG: %s" % while_replace(texta))
    try:
        molecule = ACTopol(inputFile=inputFile, chargeType=chargeType, chargeVal=chargeVal, debug=debug, multiplicity=multiplicity, atomType=atomType, force=force,
            outTopol=outTopol, engine=engine, allhdg=allhdg, basename=basename, timeTol=timeTol, qprog=qprog, ekFlag=ekFlag, verbose=verbose, gmx4=gmx4, disam=disam,
            direct=direct, is_sorted=is_sorted, chiral=chiral)

        molecule.createACTopol()
        molecule.createMolTopol()
        if not basename:
            file_name = 'smiles_molecule'
        else:
            file_name = basename

        "Output in JSON format"
        readFiles(file_name)
        output={
            "file_name":file_name,
            "em_mdp" : em_mdp.getvalue(),
            "AC_frcmod" : AC_frcmod.getvalue(),
            "AC_inpcrd" : AC_inpcrd.getvalue(),
            "AC_lib" : AC_lib.getvalue(),
            "AC_prmtop" : AC_prmtop.getvalue(),
            "bcc_gaff_mol2" : bcc_gaff_mol2.getvalue(),
            "CHARMM_inp" : CHARMM_inp.getvalue(),
            "CHARMM_prm" : CHARMM_prm.getvalue(),
            "CHARMM_rtf" : CHARMM_rtf.getvalue(),
            "CNS_inp" : CNS_inp.getvalue(),
            "CNS_par" : CNS_par.getvalue(),
            "CNS_top" : CNS_top.getvalue(),
            "GMX_OPLS_itp" : GMX_OPLS_itp.getvalue(),
            "GMX_OPLS_top" : GMX_OPLS_top.getvalue(),
            "GMX_gro" : GMX_gro.getvalue(),
            "GMX_itp" : GMX_itp.getvalue(),
            "GMX_top" : GMX_top.getvalue(),
            "NEW_pdb" : NEW_pdb.getvalue(),
            "md_mdp" : md_mdp.getvalue(),
            "sqm_in" : sqm_in.getvalue(),
            "sqm_out" : sqm_out.getvalue()}

    except Exception:
        _exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
        print("ACPYPE FAILED: %s" % exceptionValue)
        if debug:
            traceback.print_tb(exceptionTraceback, file=sys.stdout)
            output = {'file_name': "error"}

    execTime = int(round(time.time() - at0))
    if execTime == 0:
        amsg = "less than a second"
    else:
        amsg = elapsedTime(execTime)
    print("Total time of execution: %s" % amsg)
    clearFileInMemory()
    return json.dumps(output)
