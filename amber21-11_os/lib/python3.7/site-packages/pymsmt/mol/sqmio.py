"This module for SQM"

import linecache
from pymsmt.mol.mol import gauatm
from parmed.periodic_table import AtomicNum

#------------------------------------------------------------------------------
#------------------------------Write SQM input file----------------------------
#------------------------------------------------------------------------------

def write_sqm_optf(siopf, smchg, gatms):

    sqm_scf = open(siopf, 'w')
    print("Run semi-empirical minimization", file=sqm_scf)
    print(" &qmmm", file=sqm_scf)
    print(" qm_theory='PM6', grms_tol=0.0002,", file=sqm_scf)
    print(" tight_p_conv=1, scfconv=1.d-10, qmcharge=%d," %smchg, file=sqm_scf)
    print(" /", file=sqm_scf)
    for gatmi in gatms:
        nuchg = int(AtomicNum[gatmi.element])
        print("%-2s %5s %10.4f %10.4f %10.4f" \
        %(nuchg, gatmi.element, gatmi.crdx, gatmi.crdy, gatmi.crdz), file=sqm_scf)
    sqm_scf.close()

#------------------------------------------------------------------------------
#--------------------------Read info from SQM output file----------------------
#------------------------------------------------------------------------------

def get_crdinfo_from_sqm(outfile):

    gauatms = []

    ln = 1
    fp = open(outfile, 'r')
    for line in fp:
        if " Final Structure" in line:
            bln = ln + 4
        elif "--------- Calculation Completed ----------" in line:
            eln = ln - 2
        ln = ln + 1
    fp.close()

    for i in range(bln, eln+1):
        line = linecache.getline(outfile, i)
        line = line.strip('\n')
        line = line.split()
        if line[0] == 'QMMM:':
            atm = gauatm(line[3], float(line[4]), float(line[5]), float(line[6]))
            gauatms.append(atm)

    linecache.clearcache()
    return gauatms
