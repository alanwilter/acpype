"This module is for GAMESS"

import linecache
import numpy
from parmed.periodic_table import AtomicNum
from pymsmt.mol.constants import B_TO_A

#------------------------------------------------------------------------------
#--------------------------Write GAMESS input file-----------------------------
#------------------------------------------------------------------------------

def write_gms_optf(goptf2, smchg, SpinNum, gatms, signum=3):

    ##GAMESS OPT file
    optf2 = open(goptf2, 'w')
    print(" $SYSTEM MEMDDI=400 MWORDS=200 $END", file=optf2)
    print(" $CONTRL DFTTYP=B3LYP RUNTYP=OPTIMIZE ICHARG=%d MULT=%d $END" %(smchg, SpinNum), file=optf2)
    print(" $STATPT NSTEP=1000 $END", file=optf2)
    print(" $BASIS GBASIS=N31 NGAUSS=6 NDFUNC=1 $END", file=optf2)
    print(" $DATA", file=optf2)
    print("Cluster/6-31G", file=optf2)
    print("C1", file=optf2)
    optf2.close()

    for gatmi in gatms:
        write_gmsatm(gatmi, goptf2, signum)

    ##Print the last line in GAMESS input file
    ##Geometry Optimization file
    optf2 = open(goptf2, 'a')
    print(" $END", file=optf2)
    optf2.close()

def write_gms_fcf(gfcf2, smchg, SpinNum):

    ##GAMESS FC file
    fcf2 = open(gfcf2, 'w')
    print(" $SYSTEM MEMDDI=400 MWORDS=200 $END", file=fcf2)
    print(" $CONTRL DFTTYP=B3LYP RUNTYP=HESSIAN ICHARG=%d MULT=%d $END" %(smchg, SpinNum), file=fcf2)
    print(" $BASIS GBASIS=N31 NGAUSS=6 NDFUNC=1 $END", file=fcf2)
    print(" $DATA", file=fcf2)
    print("Cluster/6-31G", file=fcf2)
    print("C1", file=fcf2)
    print(" ", file=fcf2)
    print(" $END", file=fcf2)
    fcf2.close()

def write_gms_mkf(gmsf, lgchg, SpinNum, gatms, signum=3):

    #For GAMESS MK Charge file
    w_gmsf = open(gmsf, 'w')
    print(" $SYSTEM MEMDDI=400 MWORDS=200 $END", file=w_gmsf)
    print(" $CONTRL DFTTYP=B3LYP ICHARG=%d MULT=%d $END" %(lgchg, SpinNum), file=w_gmsf)
    print(" $ELPOT IEPOT=1 WHERE=PDC $END", file=w_gmsf)
    print(" $PDC PTSEL=CONNOLLY CONSTR=NONE $END", file=w_gmsf)
    print(" $BASIS GBASIS=N31 NGAUSS=6 NDFUNC=1 $END", file=w_gmsf)
    print(" $DATA", file=w_gmsf)
    print("Cluster/6-31G(d)", file=w_gmsf)
    print("C1", file=w_gmsf)
    w_gmsf.close()

    #For GAMESS file
    for gatmi in gatms:
        write_gmsatm(gatmi, gmsf, signum)

    #Print the end character for GAMESS input file
    w_gmsf = open(gmsf, 'a')
    print(' $END', file=w_gmsf)
    w_gmsf.close()

def write_gmsatm(gmsatm, fname, signum=3):

    wf = open(fname, 'a')
    element = gmsatm.element
    nuchg = AtomicNum[element]
    nuchg = round(nuchg, 1)
    if signum == 3:
        print("%-6s  %6.1f  %8.3f %8.3f %8.3f" %(gmsatm.element, \
                 nuchg, gmsatm.crdx, gmsatm.crdy, gmsatm.crdz), file=wf)
    elif signum == 4:
        print("%-6s  %6.1f  %9.4f %9.4f %9.4f" %(gmsatm.element, \
                 nuchg, gmsatm.crdx, gmsatm.crdy, gmsatm.crdz), file=wf)
    wf.close()

#------------------------------------------------------------------------------
#----------------------Read info from GAMESS output file-----------------------
#------------------------------------------------------------------------------

def get_crds_from_gms(logfile):

    unit = 'angs' #Coordinates will use angs unit in default

    bln = None
    ln = 1
    fp = open(logfile, 'r')
    for line in fp:
        if ' ATOM      ATOMIC                      COORDINATES (BOHR)' in line:
            bln = ln + 2
            unit = 'bohr'    #Means using Bohr unit
        elif ' ATOM      ATOMIC                      COORDINATES (ANGS' in line:
            bln = ln + 2
            unit = 'angs'
        elif '          INTERNUCLEAR DISTANCES' in line:
            eln = ln - 2
        ln = ln + 1
    fp.close()

    if bln is None:
        raise pymsmtError('There is no atomic coordinates found in the '
                          'GAMESS-US output file. Please check whether '
                          'the GAMESS-US jobs are finished normally, and '
                          'whether you are using the correct output file.')

    crdl = []
    for i in range(bln, eln+1):
        line = linecache.getline(logfile, i)
        line = line.strip('\n')
        line = line.split()
        if unit == 'bohr':
            crdl.append(float(line[2]))
            crdl.append(float(line[3]))
            crdl.append(float(line[4]))
        elif unit == 'angs':
            crdl.append(float(line[2])/B_TO_A)
            crdl.append(float(line[3])/B_TO_A)
            crdl.append(float(line[4])/B_TO_A)
    linecache.clearcache()

    return crdl

def get_matrix_from_gms(logfile, msize):

    ln = 1
    hasfc = 0
    fp = open(logfile, 'r')
    for line in fp:
        if 'CARTESIAN FORCE CONSTANT MATRIX' in line:
            hasfc = hasfc + 1
            bln = ln + 6
        ln = ln + 1
    fp.close()

    if hasfc > 0:
        pass
    else:
        raise pymsmtError('There is no \'CARTESIAN FORCE CONSTANT MATRIX\' '
                          'found in the GAMESS-US output file. Please check '
                          'whether the GAMESS-US jobs are finished normally, '
                          'and whether you are using the correct output file.')

    fcmatrix = numpy.array([[float(0) for x in range(msize)] for x in range(msize)])

    cycles = msize//6

    for i in range(0, cycles): #To see how many cycles need
        for j in range(0, msize-i*6): #How many lines in the cycle
            line = linecache.getline(logfile, bln)
            line = line.strip('\n')
            for k in range(0, 6): #There are 6 values in one line
                fcmatrix[j+i*6][k+i*6] = float(line[20+9*k:29+9*k])
            bln = bln + 1
        bln = bln + 4

    #If there is one more section remaining
    if (msize%6 == 3):
        for j in range(0, 3):
            line = linecache.getline(logfile, bln)
            line = line.strip('\n')
            for k in range(0, 3): #There are 6 values in one line
                fcmatrix[j + cycles * 6][k + cycles * 6] = float(line[20+9*k:29+9*k])
            bln = bln + 1

    linecache.clearcache()

    #To complete the matrix
    for j in range(0, msize):
        for k in range(0, j+1):
            fcmatrix[k][j] = fcmatrix[j][k]

    return fcmatrix

def get_esp_from_gms(logfile, espfile):

    #ESP file uses Bohr unit

    unit = 'bohr' #In default ESP coordinates will use bohr unit

    #---------For Coordinates--------
    crdl = []

    hasesp = 0
    bln1 = 0
    bln2 = None
    ln = 1
    fp = open(logfile, 'r')
    for line in fp:
        if 'COORDINATES OF ALL ATOMS ARE (ANGS)' in line:
            bln1 = ln + 3
            unit = 'angs'
        elif 'COORDINATES OF ALL ATOMS ARE (BOHR)' in line:
            bln1 = ln + 3
            unit = 'bohr'
        elif ' ATOM      ATOMIC                      COORDINATES (BOHR)' in line:
            bln2 = ln + 2
            unit = 'bohr'
        elif ' ATOM      ATOMIC                      COORDINATES (ANGS' in line:
            bln2 = ln + 2
            unit = 'angs'
        elif 'TOTAL NUMBER OF ATOMS' in line:
            atnumber = int(line.strip('\n').split()[-1])
        elif 'ELECTROSTATIC POTENTIAL' in line:
            espbln0 = ln
            hasesp = hasesp + 1
        elif 'NUMBER OF POINTS SELECTED FOR FITTING' in line:
            espnumber = int(line.strip('\n').split()[-1])
            espbln = ln + 1
            espeln = ln + espnumber

        ln = ln + 1
    fp.close()

    if (bln1 == 0) and (bln2 is None):
        raise pymsmtError('There is no atomic coordinates found in the '
                          'GAMESS-US output file. Please check whether '
                          'the GAMESS-US jobs are finished normally, and '
                          'whether you are using the correct output file.')

    if bln1 == 0:
        bln = bln2
    else:
        bln = bln1

    eln = bln + atnumber - 1

    for i in range(bln, eln+1):
        line = linecache.getline(logfile, i)
        line = line.strip('\n')
        line = line.split()
        if unit == 'angs':
            crd = (float(line[2])/B_TO_A, float(line[3])/B_TO_A, float(line[4])/B_TO_A)
        elif unit == 'bohr':
            crd = (float(line[2]), float(line[3]), float(line[4]))
        crdl.append(crd)

    linecache.clearcache()

    #---------For ESP points---------
    #ESP files use Bohr and Atomic Unit Charge

    if hasesp > 0:
        pass
    else:
        raise pymsmtError('There is no \'ELECTROSTATIC POTENTIAL\' '
                          'found in the GAMESS-US output file. Please check '
                          'whether the GAMESS-US jobs are finished normally, '
                          'and whether you are using the correct output file.')

    for i in range(espbln0, espbln):
        line = linecache.getline(logfile, i)
        if 'VDW RADIUS' in line:
            print("CAUTION: " + line.strip('\n') + " IN THE GAMESS CALCULATION.")

    espdict = {}
    for i in range(espbln, espeln+1):
        line = linecache.getline(logfile, i)
        line = line.strip('\n')
        line = line.split()
        val = (float(line[1]), float(line[2]), float(line[3]), float(line[6]))
        espdict[int(line[0])] = val

    linecache.clearcache()

    w_espf = open(espfile, 'w')
    print("%5d%5d%5d" %(len(crdl), len(espdict), 0), file=w_espf)
    for i in crdl:
        print("%16s %15.7E %15.7E %15.7E" %(' ', i[0], i[1], i[2]), file=w_espf)
    for i in range(1, len(espdict)+1):
        val = espdict[i]
        print("%16.7E %15.7E %15.7E %15.7E" %(val[3], val[0], val[1], val[2]), file=w_espf)
    w_espf.close()
