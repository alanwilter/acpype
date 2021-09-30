import sys, time, json
from acpype_lib.acpype import ACTopol, MAXTIME, while_replace, mdatFileInMemory, sleapFileInMemory, pdbFileInMemory, topFileInMemory, itpFileInMemory, oitpFileInMemory,otopFileInMemory,groFileInMemory, emMdpFileInMemory, mdMdpFileInMemory,parFileInMemory,inpFileInMemory, header, elapsedTime

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
    InMemory = True,
    is_smiles = False):

    at0 = time.time()
    print(header)

    if debug:
        texta = "Python Version %s" % sys.version
        print("DEBUG: %s" % while_replace(texta))
    try:
        molecule = ACTopol(inputFile=inputFile, chargeType=chargeType, chargeVal=chargeVal, debug=debug, multiplicity=multiplicity, atomType=atomType, force=force,
            outTopol=outTopol, engine=engine, allhdg=allhdg, basename=basename, timeTol=timeTol, qprog=qprog, ekFlag=ekFlag, verbose=verbose, gmx4=gmx4, disam=disam,
            direct=direct, is_sorted=is_sorted, chiral=chiral, InMemory=InMemory)

        molecule.createACTopol()
        molecule.createMolTopol()
        if not basename:
            file_name = 'smiles_molecule.mol2'
        else:
            file_name = basename

        "Output in JSON format"
        output={
            "file_name":file_name,
            "pdbFile":pdbFileInMemory.getvalue(),
            "topFile":topFileInMemory.getvalue(),
            "itpFile":itpFileInMemory.getvalue(),
            "oitpFile":oitpFileInMemory.getvalue(),
            "otopFile":otopFileInMemory.getvalue(),
            "groFile":groFileInMemory.getvalue(),
            "emMdpFile":emMdpFileInMemory.getvalue(),
            "mdMdpFile":mdMdpFileInMemory.getvalue(),
            "parFile":parFileInMemory.getvalue(),
            "inpFile":inpFileInMemory.getvalue()}

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

    return json.dumps(output)
