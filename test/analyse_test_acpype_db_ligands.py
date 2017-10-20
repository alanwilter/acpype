#!/usr/bin/env python
from commands import getoutput

# Gmm with 200 atoms, biggest OK

# usage: ./analyse_test_acpype_db_ligands.py (opt: n=10 or "['001','Rca']")

# semi-QM = mopac (AT 1.2) or sqm (AT 1.3)

import os, sys, fnmatch

global id
id = 1

atomicDetailed = True

checkChargeComparison = True

results = {}

ccpCodes = os.listdir('other')
ccpCodes.sort()

groupResults = [
               ["dirs totally empty or at least missing one *.mol2 input files for either PDB or IDEAL)", "Dirs missing *.mol2 input files",
                []],
               ["still running presumably", "Mols still running presumably",
                ['0 E, 0 W, ET , WT ']],
               ["mols clean, no erros or warnings","Mols clean",
                ['0 E, 1 W, ET , WT _0','0 E, 2 W, ET , WT _0_1',
                 '0 E, 2 W, ET , WT _0_2','0 E, 3 W, ET , WT _0_1_2',
                 '0 E, 2 W, ET , WT _0_3','0 E, 3 W, ET , WT _0_1_3',
                 '0 E, 3 W, ET , WT _0_2_3','0 E, 4 W, ET , WT _0_1_2_3',
                 '0 E, 2 W, ET , WT _0_8']],
               ["only guessCharge failed, but running semi-QM with charge = 0 finished fine","Mols only guessCharge failed",
                ['1 E, 1 W, ET _0, WT _0', '1 E, 2 W, ET _0, WT _0_1',
                 '1 E, 2 W, ET _0, WT _0_2', '1 E, 3 W, ET _0, WT _0_1_2']],
               ["atoms in close contact", "Mols have atoms in close contact",
                ['0 E, 2 W, ET , WT _0_7', '0 E, 3 W, ET , WT _0_1_7',
                 '0 E, 3 W, ET , WT _0_2_7', '0 E, 3 W, ET , WT _0_3_7',
                 '0 E, 4 W, ET , WT _0_1_2_7', '0 E, 4 W, ET , WT _0_1_3_7',
                 '0 E, 4 W, ET , WT _0_2_3_7', '0 E, 5 W, ET , WT _0_1_2_3_7']],
               ["irregular bonds", "Mols have irregular bonds",
                ['0 E, 2 W, ET , WT _0_5', '0 E, 3 W, ET , WT _0_1_5',
                 '0 E, 3 W, ET , WT _0_2_5', '0 E, 4 W, ET , WT _0_1_2_5']],
               ["irregular bonds and atoms in close contact", "Mols have irregular bonds and atoms in close contact",
                ['0 E, 3 W, ET , WT _0_5_7', '0 E, 4 W, ET , WT _0_1_5_7',
                 '0 E, 4 W, ET , WT _0_2_5_7', '0 E, 5 W, ET , WT _0_1_2_5_7',
                 '0 E, 6 W, ET , WT _0_1_2_3_5_7', '0 E, 4 W, ET , WT _0_3_5_7']],
               ["couldn't determine all parameters", "Mols have missing parameters",
                ['0 E, 3 W, ET , WT _0_2_4', '0 E, 3 W, ET , WT _0_1_4', '0 E, 2 W, ET , WT _0_4']],
               ["missing parameters and atoms in close contact", "Mols have missing parameters and atoms in close contact",
                ['0 E, 4 W, ET , WT _0_1_4_7']],
               ["missing parameters, irregular bonds and atoms in close contact", "Mols have missing parameters, irregular bonds and atoms in close contact",
                ['0 E, 5 W, ET , WT _0_2_4_5_7']],
               ["missing parameters, irregular bonds, maybe wrong atomtype and atoms in close contact", "Mols have missing parameters, irregular bonds, maybe wrong atomtype and atoms in close contact",
                ['0 E, 5 W, ET , WT _0_4_5_6_7']],
               ["no 'tmp', acpype did nothing at all", "Mols have no 'tmp'",
                ['1 E, 0 W, ET _7, WT ']],
               ["atoms with same coordinates", "Mols have duplicated coordinates",
                ['1 E, 0 W, ET _1, WT ']],
               ["maybe wrong atomtype", "Mols with maybe wrong atomtype",
                ['0 E, 2 W, ET , WT _0_6','0 E, 3 W, ET , WT _0_1_6','0 E, 4 W, ET , WT _0_1_3_6','0 E, 3 W, ET , WT _0_3_6']],
               ["maybe wrong atomtype and atoms in close contact", "Mols with maybe wrong atomtype and atoms in close contact",
                ['0 E, 4 W, ET , WT _0_1_6_7', '0 E, 3 W, ET , WT _0_6_7', '0 E, 4 W, ET , WT _0_3_6_7']],
               ["irregular bonds, maybe wrong atomtype and atoms in close contact", "Mols with irregular bonds, maybe wrong atomtype and atoms in close contact",
                ['0 E, 5 W, ET , WT _0_1_5_6_7']],
               ["guessCharge failed and atoms in close contact", "Mols have guessCharge failed and atoms in close contact",
                ['1 E, 3 W, ET _0, WT _0_1_7', '1 E, 3 W, ET _0, WT _0_2_7',
                 '1 E, 4 W, ET _0, WT _0_1_2_7']],
               ["guessCharge failed and missing parameters", "Mols have guessCharge failed and missing parameters",
                ['1 E, 2 W, ET _0, WT _0_4', '1 E, 3 W, ET _0, WT _0_1_4',
                 '1 E, 3 W, ET _0, WT _0_2_4', '1 E, 4 W, ET _0, WT _0_1_2_4']],
               ["guessCharge failed and maybe wrong atomtype", "Mols have guessCharge failed and maybe wrong atomtype",
                ['1 E, 2 W, ET _0, WT _0_6', '1 E, 3 W, ET _0, WT _0_1_6',
                 '1 E, 3 W, ET _0, WT _0_2_6', '1 E, 4 W, ET _0, WT _0_1_2_6']],
               ["guessCharge failed and irregular bonds", "Mols have guessCharge failed and irregular bonds",
                ['1 E, 3 W, ET _0, WT _0_2_5']],
               ["guessCharge failed, missing parameters and maybe wrong atomtype", "Mols have guessCharge failed, missing parameters and maybe wrong atomtype",
                ['1 E, 3 W, ET _0, WT _0_4_6', '1 E, 4 W, ET _0, WT _0_1_4_6',
                 '1 E, 4 W, ET _0, WT _0_2_4_6']],
               ["guessCharge failed, irregular bonds and maybe wrong atomtype", "Mols have guessCharge failed, irregular bonds and maybe wrong atomtype",
                ['1 E, 4 W, ET _0, WT _0_2_5_6', '0 E, 3 W, ET , WT _0_5_6']],
               ["guessCharge failed, missing parameters and atoms in close contact", "Mols have guessCharge failed, missing parameters and atoms in close contact",
                ['1 E, 3 W, ET _0, WT _0_4_7', '1 E, 4 W, ET _0, WT _0_2_4_7']],
               ["guessCharge failed, maybe wrong atomtype and atoms in close contact", "Mols have guessCharge failed, maybe wrong atomtype and atoms in close contact",
                ['1 E, 3 W, ET _0, WT _0_6_7', '1 E, 4 W, ET _0, WT _0_1_6_7',
                 '1 E, 4 W, ET _0, WT _0_2_6_7', '1 E, 5 W, ET _0, WT _0_1_2_6_7']],
               ["guessCharge failed, irregular bonds and atoms in close contact", "Mols have guessCharge failed, irregular bonds and atoms in close contact",
                ['1 E, 4 W, ET _0, WT _0_2_5_7']],
               ["guessCharge failed, irregular bonds, maybe wrong atomtype and atoms in close contact", "Mols have guessCharge failed, irregular bonds, maybe wrong atomtype and atoms in close contact",
                ['1 E, 5 W, ET _0, WT _0_1_5_6_7', '1 E, 5 W, ET _0, WT _0_2_5_6_7']],
               ["atoms too close", "Mols have atoms too close",
                ['1 E, 0 W, ET _2, WT ']],
               ["atoms too alone", "Mols have atoms too alone",
                ['1 E, 0 W, ET _3, WT ']],
               ["tleap failed", "Mols have tleap failed",
                ['3 E, 1 W, ET _4_5_6, WT _0', '3 E, 2 W, ET _4_5_6, WT _0_1']],
               ["tleap failed, maybe wrong atomtype", "Mols have tleap failed and maybe wrong atomtype",
                ['3 E, 3 W, ET _4_5_6, WT _0_1_6', '3 E, 2 W, ET _4_5_6, WT _0_6']],
               ["semi-QM timeout", "Mols have semi-QM timeout",
                ['1 E, 2 W, ET _9, WT _0_1', '1 E, 1 W, ET _9, WT _0']],
               ["semi-QM timeout and maybe wrong atomtype", "Mols have semi-QM timeout and maybe wrong atomtype",
                ['1 E, 3 W, ET _9, WT _0_1_6', '1 E, 2 W, ET _9, WT _0_6']],
               ["guessCharge and tleap failed", "Mols have guessCharge and tleap failed",
                ['4 E, 1 W, ET _0_4_5_6, WT _0', '4 E, 2 W, ET _0_4_5_6, WT _0_1']],
               ["guessCharge and tleap failed, maybe wrong atomtype", "Mols have guessCharge and tleap failed, maybe wrong atomtype",
                ['4 E, 2 W, ET _0_4_5_6, WT _0_6', '4 E, 3 W, ET _0_4_5_6, WT _0_1_6']],
               ["guessCharge failed and semi-QM timeout", "Mols have guessCharge failed and semi-QM timeout",
                ['2 E, 1 W, ET _0_9, WT _0']],
               ["atoms with same coordinates and maybe wrong atomtype", "Mols have atoms with same coordinates and maybe wrong atomtype",
                ['1 E, 1 W, ET _1, WT _6']],
               ["atoms too close and maybe wrong atomtype", "Mols have atoms too close and maybe wrong atomtype",
                ['1 E, 1 W, ET _2, WT _6']],
               ["atoms too alone and maybe wrong atomtype", "Mols have atoms too alone and maybe wrong atomtype",
                ['1 E, 1 W, ET _3, WT _6']],
               ["atoms with same coordinates and too close", "Mols have atoms with same coordinates and too close",
                ['2 E, 0 W, ET _1_2, WT ']],
               ["atoms with same coordinates and too alone", "Mols have atoms with same coordinates and too alone",
                ['2 E, 0 W, ET _1_3, WT ']],
               ["atoms with same coordinates, too alone and maybe wrong atomtype", "Mols have atoms with same coordinates, too alone and maybe wrong atomtype",
                ['2 E, 1 W, ET _1_3, WT _6']],
               ["atoms with same coordinates, too close and too alone", "Mols have atoms with same coordinates, too close and too alone",
                ['3 E, 0 W, ET _1_2_3, WT ']],
               ]

error_warn_messages = \
'''
    warnType 0 = 'no charge value given, trying to guess one...'
    warnType 1 = "In ..., residue name will be 'R instead of ... elsewhere"
    warnType 2 = 'residue label ... is not all UPPERCASE'

    # mild warning (need to be sure the charge guessed is correct)
    warnType 3 = "The unperturbed charge of the unit ... is not zero"

    # serious warnings (topology not reliable):
    warnType 4 = "Couldn't determine all parameters"
    warnType 5 = 'There is a bond of ... angstroms between'
    warnType 6 = 'atom type may be wrong'
    warnType 7 = 'Close contact of ... angstroms between ...'

    #warnType 8 = 'no 'babel' executable, no PDB file as input can be used!'
    warnType 8 = "residue name will be 'MOL' instead of"
    warnType 9 = 'UNKNOWN WARN'

    errorType 0 = 'guessCharge failed' # (not so serious if only err    or because semi-QM worked with charge Zero)
    errorType 1 = 'Atoms with same coordinates in'
    errorType 2 = 'Atoms TOO close'
    errorType 3 = 'Atoms TOO alone'
    errorType 4 = 'Antechamber failed'
    errorType 5 = 'Parmchk failed'
    errorType 6 = 'Tleap failed'
    errorType 7 = "No such file or directory: 'tmp'" # can be bondtyes wrong or wrong frozen atom type
    errorType 8 = 'syntax error'
    errorType 9 = 'Semi-QM taking too long to finish'
    errorType 10 = 'UNKNOWN ERROR'
'''

totalPdbOkCount = 0
totalIdealOkCount = 0
totalIdealMol2Count = 0
totalPdbMol2Count = 0

emptyDir = []
failedPdb = []
failedIdeal = []
failedBoth = []
missPdbMol2 = []
missIdealMol2 = []
multIdealMol2 = []
multPdbMol2 = []

ET0 = []
ET1 = []
ET2 = []
ET3 = []
ET4 = []
ET5 = []
ET6 = []
ET7 = set()
ET8 = set()
ET9 = []

WT3 = set()
WT4 = []
WT5 = set()
WT6 = set()
WT7 = set()

mapResults = {}

execTime = {}

compareCharges = set()
compareChargesOK = set()
# ==> Trying with net charge = 0 ...
# ==> ... charge set to 0
# Overall charge: 0

SCFfailedList = []

def analyseFile(mol, structure, file):
    """
         mol = string molDir, e.g. 'Gnp'
         structure = string 'pdb' or 'ideal'
         file = string e.g. '10a/10a.none_neutral.pdb.out'
         returns e.g. 'pdb: 0 E, 3 W, ET , WT _0_1_7'
    """
    _flagChargeType = None # if charge was guesed correctly ('Y') or tried with 0 ('Z' for Zero)
    guessChargeValue = None
    warnTypes = ''
    errorTypes = ''
    tmpFile = open(file, 'r')
    content = tmpFile .readlines()
    for line in content:
        if 'WARNING: ' in line:
            if "no charge value given" in line:
                warnTypes += '_0'
            elif "residue name will be 'R" in line:
                warnTypes += '_1'
            elif "not all UPPERCASE" in line:
                warnTypes += '_2'
            elif "applications like CNS" in line:
                pass
            elif "The unperturbed charge of the" in line:
                warnTypes += '_3'
                charge = round(float(line[44:54]))
                WT3.add('%s: %s' % (mol, charge))
            elif "Couldn't determine all parameters" in line:
                warnTypes += '_4'
                WT4.append('%s_%s'% (mol, structure))
            elif "There is a bond of" in line:
                warnTypes += '_5'
                _dist = round(float(line[27:37]))
                WT5.add('%s_%s' % (mol, structure))
            elif ' atom type of ' in line:
                warnTypes += '_6'
                WT6.add('%s_%s'% (mol, structure))
            #elif "no 'babel' executable, no PDB" in line:
            elif "residue name will be 'MOL' instead of" in line:
                warnTypes += '_8'
            else:
                print "UNKNOWN WARN:", file, line
                warnTypes += '_9'
        if 'Warning: Close contact of ' in line:
            warnTypes += '_7'
            WT7.add('%s_%s'% (mol, structure))

        if 'ERROR: ' in line:
            if 'guessCharge failed' in line:
                errorTypes += '_0'
                ET0.append('%s_%s'% (mol, structure))
            elif 'Atoms with same coordinates in' in line:
                errorTypes += '_1'
                ET1.append('%s_%s'% (mol, structure))
            elif 'Atoms TOO close' in line:
                errorTypes += '_2'
                ET2.append('%s_%s'% (mol, structure))
            elif 'Atoms TOO alone' in line:
                errorTypes += '_3'
                ET3.append('%s_%s'% (mol, structure))
            elif 'Antechamber failed' in line:
                errorTypes += '_4'
                ET4.append('%s_%s'% (mol, structure))
            elif 'Parmchk failed' in line:
                errorTypes += '_5'
                ET5.append('%s_%s'% (mol, structure))
            elif 'Tleap failed' in line:
                errorTypes += '_6'
                ET6.append('%s_%s'% (mol, structure))
            elif 'syntax error' in line:
                errorTypes += '_8'
                ET8.add('%s_%s'% (mol, structure))
            elif "Use '-f' option if you want to proceed anyway. Aborting" in line:
                pass
            else:
                print "UNKNOWN ERROR:", file, line
                errorTypes += '_10'
        if "No such file or directory: 'tmp'" in line:
            errorTypes += '_7'
            ET7.add('%s_%s'% (mol, structure))
        if "Semi-QM taking too long to finish" in line:
            errorTypes += '_9'
            ET9.append('%s_%s'% (mol, structure))
        if "Total time of execution:" in line:
            if execTime.has_key(mol):
                #execTime[mol].append({structure:line[:-1].split(':')[1]})
                execTime[mol][structure] = line[:-1].split(':')[1]
            else:
                #execTime[mol] = [{structure:line[:-1].split(':')[1]}]
                execTime[mol] = {structure:line[:-1].split(':')[1]}
    # to compare charges from acpype out with input mol2
        if "==> ... charge set to" in line:
            guessChargeValue = int(line.split('to')[1])
            _flagChargeType = 'Y'
            #print "*%s*" % guessChargeValue
        elif "Trying with net charge =" in line:
            guessChargeValue = 0
            _flagChargeType = 'Z'

    if checkChargeComparison:
        mol2FileName = file.replace('.out','.mol2')
        mol2File = open(mol2FileName,'r').readlines()
        for line in mol2File:
            if "# Overall charge:" in line:
                mol2Charge = int(line.split(':')[1])
        if guessChargeValue:
            if mol2Charge != guessChargeValue:
                compareCharges.add("%s_%i_%i" % (mol,mol2Charge, guessChargeValue))
            else:
                compareChargesOK.add("%s_%i" % (mol, mol2Charge))
            # this excpetion should happen only if acpype early aborted for both outs
            #if mol not in `compareChargesOK.union(compareCharges)`:
            if mol not in `compareChargesOK` and mol not in `compareCharges`:
                print "!!!! Unable to compare. Failed to guess charge for", mol, structure

    out = parseSummurisedLine(warnTypes, errorTypes)
    if mapResults.has_key(out):
        if mapResults[out].has_key(mol):
            mapResults[out][mol].append(structure)
        else:
            mapResults[out][mol] = [structure]
    else:
        mapResults[out] = {mol:[structure]}
    return out

def parseSummurisedLine(warnTypes, errorTypes):
    wt = list(set(warnTypes.split('_')))
    if wt != ['']:
        wt = wt[1:]
        wt.sort(cmp=lambda x,y: int(x)-int(y))
    warnTypes = ''
    for i in wt:
        if i: warnTypes += '_'+i
    countWarn = warnTypes.count('_')
    et = list(set(errorTypes.split('_')))
    et.sort()
    errorTypes = ''
    for j in et:
        if j: errorTypes += '_'+j
    countError = errorTypes.count('_')
    return "%i E, %i W, ET %s, WT %s" % (countError, countWarn, errorTypes, warnTypes)

def parseChargeList(WT, dList):
    listOk = []
    for wt in WT:
        if wt.split(':')[0] in dList:
            listOk.append(wt)
    listOk.sort()
    return listOk

def myComp(vx,vy):
    ix = int(vx.split()[0])
    iy = int(vy.split()[0])
    tx = vx[-1:]
    ty = vy[-1:]
    if tx.isdigit(): x = int(tx)
    else: x = 0
    if ty.isdigit(): y = int(ty)
    else: y = 0

    if ix>iy:
            return 1
    elif ix==iy:
        if x>y:
            return 1
        elif x==y:
            return 0
        else:
            return -1
    else:
        if x>y:
            return 1
        elif x==y:
            return 0
        else:
            return -1

def sortList(lista,typeMess):
    for mol, l in mapResults[typeMess].items():
        if len(l) == 2:
            lista[0].append(mol)
        elif len(l) == 1:
            if l[0] == 'pdb': lista[1].append(mol)
            elif l[0] == 'ideal': lista[2].append(mol)
        else:
            print "problem with", typeMess, mol, l
    return lista

def printResults(lista, subHead, header=None):
    '''
    print results as
    --------------------------------------------------------------------------------

    *** For results [1], [2], [3], totally empty dirs (NO mol2 input files for either PDB or IDEAL):

    [1] Dirs missing pdb.mol2 input files for both PDB and Ideal:
    0     []

    [2] Dirs missing pdb.mol2 input files with PDB ONLY, besides [1]:
    0     []

    [3] Dirs missing pdb.mol2 input files with IDEAL ONLY, besides [1]:
    1     ['006']
    PDB Total: 0 of 20 (0.00%)
    IDEAL Total: 1 of 20 (5.00%)
    Total: 1 of 20 (5.00%)
    --------------------------------------------------------------------------------
    '''
    global id
    dList, pList, iList = lista
    if not dList and not pList and not iList:
        return 0
    dList.sort()
    pList.sort()
    iList.sort()
    id1 = id + 1
    id2 = id + 2
    print 80*'-'
    if header:
        print "\n*** For results [%i], [%i], [%i], %s:" % (id,id1,id2,header)
    print '\n[%i] %s for both PDB and Ideal:\n%i\t %s' % (id, subHead, len(dList), str(dList))
    print '\n[%i] %s with PDB ONLY, besides [%i]:\n%i\t %s' % (id1, subHead, id, len(pList), str(pList))
    print '\n[%i] %s with IDEAL ONLY, besides [%i]:\n%i\t %s' % (id2, subHead, id, len(iList), str(iList))
    pTotal = len(dList)+len(pList)
    iTotal = len(dList)+len(iList)
    total = pTotal+iTotal
    per = total / (allTotal * 0.01)
    perPdb = pTotal / (allTotal * 0.01)
    perIdeal = iTotal / (allTotal * 0.01)
    print "PDB Total: %i of %i (%3.2f%%)" % (pTotal, allTotal, perPdb)
    print "IDEAL Total: %i of %i (%3.2f%%)" % (iTotal, allTotal, perIdeal)
    print "Total: %i of %i (%3.2f%%)" % (total, allTotal, per)
    print 80*'-'
    id += 3
    return total

def elapsedTime(seconds, suffixes=['y','w','d','h','m','s'], add_s=False, separator=' '):
    """
    Takes an amount of seconds and turns it into a human-readable amount of time.
    """
    # the formatted time string to be returned
    if seconds == 0:
        return '0s'
    time = []

    # the pieces of time to iterate over (days, hours, minutes, etc)
    # - the first piece in each tuple is the suffix (d, h, w)
    # - the second piece is the length in seconds (a day is 60s * 60m * 24h)
    parts = [(suffixes[0], 60 * 60 * 24 * 7 * 52),
          (suffixes[1], 60 * 60 * 24 * 7),
          (suffixes[2], 60 * 60 * 24),
          (suffixes[3], 60 * 60),
          (suffixes[4], 60),
          (suffixes[5], 1)]

    # for each time piece, grab the value and remaining seconds, and add it to
    # the time string
    for suffix, length in parts:
        value = seconds / length
        if value > 0:
            seconds = seconds % length
            time.append('%s%s' % (str(value),
                           (suffix, (suffix, suffix + 's')[value > 1])[add_s]))
        if seconds < 1:
            break

    return separator.join(time)

def convertStringTime2Seconds(stringTime):
    if 'less' in stringTime:
        return 0
    vecTime = stringTime.split()
    totalSec = 0
    for n in vecTime:
        if 'h' in n:
            totalSec += eval(n[:-1])*3600
        elif 'm' in n:
            totalSec += eval(n[:-1])*60
        elif 's' in n:
            totalSec += eval(n[:-1])
    return totalSec

def locate(pattern, root=os.curdir):
    '''Locate all files matching supplied filename pattern in and below
    supplied root directory.'''
    for path, _dirs, files in os.walk(os.path.abspath(root)):
        for filename in fnmatch.filter(files, pattern):
            yield os.path.join(path, filename)

os.chdir('other')

# Only run this on 'other'!
if len(sys.argv) > 1:
    args = sys.argv[1:]
    if args[0].startswith('n='):
        num = int(args[0][2:])
        ccpCodes = ccpCodes[:num]
    else:
        args = list(set(eval(args[0])))
        args.sort()
        ccpCodes = args

for molDir in ccpCodes:
#    files = os.listdir(molDir)
    files = []
#    for dirpath, dirnames, filenames in os.walk(molDir):
#        files += filenames
    files = list(locate("*",molDir))
    results[molDir] = []
    pdb = False
    ideal = False
    pdbMol2 = False
    idealMol2 = False
    countPdbMol2 = 1
    countIdealMol2 = 1
    out1 = ''
    out2 = ''
    pass1 = False
    pass2 = False
    if not files:
        emptyDir.append(molDir)
    for file in files:
        if 'pdb_NEW.pdb' in file:
            pdb = True
            results[molDir].append('pdb_OK')
            totalPdbOkCount += 1
        elif 'ideal_NEW.pdb' in file:
            ideal = True
            results[molDir].append('ideal_OK')
            totalIdealOkCount += 1
        elif 'ideal.out' in file:
            out1 = analyseFile(molDir, 'ideal', os.path.join(molDir, file))
            results[molDir].append('ideal: '+out1)
        elif 'pdb.out' in file:
            out2 = analyseFile(molDir, 'pdb', os.path.join(molDir, file))
            results[molDir].append('pdb: '+out2)
        elif '.ideal.mol2' in file:
            if idealMol2:
                countIdealMol2 +=1
            else:
                idealMol2 = True
                totalIdealMol2Count += 1
        elif '.pdb.mol2' in file:
            if pdbMol2:
                countPdbMol2 +=1
            else:
                pdbMol2 = True
                totalPdbMol2Count += 1
        elif ".ideal.acpype/sqm.out" in file:
            content = open(file, 'r').read()
            if "No convergence in SCF after" in content:
                SCFfailedList.append("%s_%s" % (molDir, 'ideal'))
        elif ".pdb.acpype/sqm.out" in file:
            content = open(file, 'r').read()
            if "No convergence in SCF after" in content:
                SCFfailedList.append("%s_%s" % (molDir, 'pdb'))
    if countIdealMol2 > 2: multIdealMol2.append(molDir)
    if countPdbMol2 > 2: multPdbMol2.append(molDir)
    if not pdbMol2: missPdbMol2.append(molDir)
    if not idealMol2: missIdealMol2.append(molDir)
    if not pdb: failedPdb.append(molDir)
    if not ideal: failedIdeal.append(molDir)
    if not pdb and not ideal: failedBoth.append(molDir)

a, b, c = set(emptyDir), set(missPdbMol2), set(missIdealMol2)
pdbList = list(b.difference(a))
idealList = list(c.difference(a))
groupResults[0].append([emptyDir, pdbList, idealList])

c = 1
while c < len(groupResults):
    groupResults[c].append([[],[],[]])
    c += 1

keys = mapResults.keys()
keys.sort()
keys.sort(cmp=myComp)

# Print messages that was not classified yet
groupMess = []
for m in groupResults:
    groupMess += m[2]
for typeMess in keys:
    if typeMess not in groupMess:
        size = len(typeMess)
        txt = str(mapResults[typeMess])
        Nboth = txt.count("['pdb', 'ideal']") + txt.count("['ideal', 'pdb']")
        Npdb = txt.count("['pdb']")
        Nideal = txt.count("['ideal']")
        print '*%s*%s%i %i %i' % (typeMess,(40-size)*' ', 2*Nboth, Npdb, Nideal)

# Append to groupResults[index] a list of Mols belonging to this group
for typeMess in keys:
    index = 0
    for group in groupResults:
        msg = group[2]
        if typeMess in msg:
            groupResults[index][3] = sortList(groupResults[index][3], typeMess)
        index += 1

allTotal = len(ccpCodes) *2
jobsOK = [] # list of Mols that have at least a PDB or IDEAL clean: [[both],[pdb,ideal]]
totalTxt = ''
for group in groupResults:
    header, subHead, dummy, lista = group
    if "Mols clean" == subHead:
        jobsOK = lista
    subTot = printResults(lista, subHead, header)
    if not subTot: continue
    totalTxt += "+ %i " % subTot
totalTxt = totalTxt[2:]
sumVal = eval(totalTxt)
print "%s= %s" % (totalTxt, sumVal)

#print results

#print 'execTime', execTime

#print 'jobsOK', jobsOK

if atomicDetailed:
    print "\n>>> Detailed report per atom <<<\n"

    if ET1:
        print "=>Mols have duplicated coordinates"
        for molLabel in ET1:
            mol,structure = molLabel.split('_')
            cmd = "grep -e '^ATOM' %s/*%s.out" % (mol, structure)
            out = getoutput(cmd)
            print "# %s #" % molLabel
            print out,"\n"

    if WT7:
        print "=>Mols have atoms in close contact"
        for molLabel in WT7:
            mol,structure = molLabel.split('_')
            cmd = "grep -e '^Warning: Close contact of' %s/*%s.out" % (mol, structure)
            out = getoutput(cmd)
            print "# %s #" % molLabel
            print out,"\n"

    if WT5:
        print "=>Mols have irregular bonds"
        for molLabel in WT5:
            mol,structure = molLabel.split('_')
            cmd = "grep -A 1 -e '^WARNING: There is a bond of' %s/*%s.out" % (mol, structure)
            out = getoutput(cmd)
            print "# %s #" % molLabel
            print out,"\n"

    if ET2:
        print "=>Mols have atoms too close"
        for molLabel in ET2:
            mol,structure = molLabel.split('_')
            cmd = "grep -e '^ .*ATOM' %s/*%s.out" % (mol, structure)
            out = getoutput(cmd)
            print "# %s #" % molLabel
            print out,"\n"

    if ET3:
        print "=>Mols have atoms too alone"
        for molLabel in ET3:
            mol,structure = molLabel.split('_')
            cmd = '''grep -e "^\['ATOM" %s/*%s.out''' % (mol, structure)
            out = getoutput(cmd)
            print "# %s #" % molLabel
            print out,"\n"

# mols with MOL2 charge different from ACPYPE guessed charge
if compareCharges:
    print "\n>>> Mols with MOL2 charge different from ACPYPE guessed charge <<<\n"
    lCC = list(compareCharges)
    lCC.sort()
    print len(lCC), lCC

maxExecTime = 0
firstPass = True
maxMolTime = None
minMolTime = None
totalCleanExecTime = 0
nJobs = 0
listMolTime = []
for mol in jobsOK[0]:
    listMolTime.append([mol+'_pdb',execTime[mol]['pdb']])
    listMolTime.append([mol+'_ideal',execTime[mol]['ideal']])
for mol in jobsOK[1]:
    listMolTime.append([mol+'_pdb',execTime[mol]['pdb']])
for mol in jobsOK[2]:
    listMolTime.append([mol+'_ideal',execTime[mol]['ideal']])

#print 'listMolTime', listMolTime
#print jobsOK

if SCFfailedList:
    SCFfailedList.sort()
    print "\n>>>Mol Jobs whose sqm.out has 'No convergence in SCF'<<<\n"
    print len(SCFfailedList), SCFfailedList
    molsOKinSCFfailedList = []
    for item in SCFfailedList:
        if item in [x[0] for x in listMolTime]:
            molsOKinSCFfailedList.append(item)
    if molsOKinSCFfailedList:
        molsOKinSCFfailedList.sort()
        print "\n>>>Mol Jobs whose sqm.out has 'No convergence in SCF' but finished OK<<<\n"
        print len(molsOKinSCFfailedList), molsOKinSCFfailedList

for item in listMolTime:
    mol, jtime = item
    tSec = convertStringTime2Seconds(jtime)
#    print 'tSec', tSec
    if firstPass:
        minExecTime = tSec
        firstPass = False
#        print 'firstPass'
    if tSec >= maxExecTime:
        maxExecTime = tSec
        maxMolTime = mol
    if tSec <= minExecTime:
        minExecTime = tSec
        minMolTime = mol
    totalCleanExecTime += tSec
    nJobs += 1

print "\n>>> Time Job Execution Summary <<<\n"
if listMolTime:
    #print maxExecTime, minExecTime
    print "Number of clean jobs:", nJobs
    print "Longest job: Mol='%s', time= %s" % (maxMolTime, elapsedTime(maxExecTime))
    print "Fatest job: Mol='%s', time= %s" % (minMolTime, elapsedTime(minExecTime))
    print "Average time of execution per clean job: %s" % elapsedTime(totalCleanExecTime/nJobs)
else:
    print "NO time stats available for clean jobs"

# Global average exec time
totalGlobalExecTime = 0
nGJobs = 0
for item in execTime.items():
    mol, tDict = item
    if tDict.has_key('pdb'):
        totalGlobalExecTime += convertStringTime2Seconds(tDict['pdb'])
        nGJobs += 1
    if tDict.has_key('ideal'):
        totalGlobalExecTime += convertStringTime2Seconds(tDict['ideal'])
        nGJobs += 1
print "\nTotal number of jobs:", nGJobs
if nGJobs:
    print "Global average time of execution per job: %s" % elapsedTime(totalGlobalExecTime/nGJobs)

# mols with charge not 0
#print WT3
