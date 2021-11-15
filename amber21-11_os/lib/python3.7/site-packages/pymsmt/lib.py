"""
This module is used for getting the parameter information from mol2,
dat, and frcmod files.
"""

from pymsmt.mol.mol2io import get_atominfo
from pymsmt.exp import *
import os
import numpy
from scipy.optimize import curve_fit
import linecache
import re

#-----------------------------------------------------------------------------
#Get the charge library file and the atom type definition for the amino acid
#atoms
#-----------------------------------------------------------------------------

#-------------Regular representation of different parameters------------------
# Float: (\s+\-?\d+\.\d*)
# Int: (\s+\d+)
# Note: (\s+\w.*|.*)

# Example: "C  12.01         0.616  !            sp2 C carbonyl group"
_massre = re.compile(r'^(\w.)(\s+\-?\d+\.\d*)(\s+\w.*|.*)$')

# Example: "OW-HW  553.0    0.9572    ! TIP3P water"
_bondre = re.compile(r'^(\w.\-\w.)(\s+\-?\d+\.\d*)(\s+\-?\d+\.\d*)(\s+\w.*|.*)$')

# Example: "HW-OW-HW    100.      104.52    TIP3P water"
_angre = re.compile(r'^(\w.\-\w.\-\w.)(\s+\-?\d+\.\d*)(\s+\-?\d+\.\d*)(\s+\w.*|.*)$')

# Example: "X -C -C -X    4   14.50        180.0             2.         Junmei et al, 1999"
#          "br-c3-c3-br   1    0.500         0.000          -3          m9 GA AUE=0.9626 RMSE=1.1958 TorType=2"
_dihre = re.compile(r'^(\w.\-\w.\-\w.-\w.)(\s+\d+)(\s+\-?\d+\.\d*)(\s+\-?\d+\.\d*)(\s+\-?\d\.?\d*)(\s+\w.*|.*)$')

# Example: "X -X -C -O          10.5         180.          2.           JCC,7,(1986),230"
_impre = re.compile(r'^(\w.\-\w.\-\w.-\w.)(\s+\-?\d+\.\d*)(\s+\-?\d+\.\d*)(\s+\-?\d\.\d*)(\s+\w.*|.*)$')

# Example: "  H           0.6000  0.0157            !Ferguson base pair geom."
# _nbre = re.compile(r'^(\s*\w.{1,3})(\s+\-?\d+\.\d*)(\s+\-?\d+\.\d*)(\s+\w.*|.*)$')

# Example: "  OW  O3           1.775931    0.162750    1.860500    0.210000"
_ljedre = re.compile(r'(\s*\w.)(\s+\w.)(\s+\-?\d+\.\d*)(\s+\-?\d+\.\d*)(\s+\-?\d+\.\d*)(\s+\-?\d+\.\d*)(\s+\w.*|.*)')

#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Define the force field dictionary
#-----------------------------------------------------------------------------

#Test amberhome or msmthome
amberhome = os.getenv('AMBERHOME')
if amberhome is None:
    raise pymsmtError('Could not perform modeling without setting $AMBERHOME '
                      'in the computer setting.')
else:
    cmdadd = amberhome + '/dat/leap/cmd/'
    libadd = amberhome + '/dat/pymsmt/'
    parmadd = amberhome + '/dat/leap/parm/'

class force_field:
    def __init__(self, sleaprcf, mol2f, datf, frcmodfs=[]):
        self.sleaprcf = sleaprcf
        self.lleaprcf = cmdadd + sleaprcf
        self.mol2f = libadd + mol2f
        self.datf = parmadd + datf
        if frcmodfs != []:
            self.frcmodfs = [parmadd + i for i in frcmodfs]
        else:
            self.frcmodfs = frcmodfs

#Define the avaliable force fields

# Old force fields
ff94 = force_field('oldff/leaprc.ff94', 'parm94.mol2', 'parm94.dat')

ff99 = force_field('oldff/leaprc.ff99', 'parm94.mol2', 'parm99.dat')

ff99SB = force_field('oldff/leaprc.ff99', 'parm94.mol2', 'parm99.dat',
    ['frcmod.ff99SB'])

ff03 = force_field('oldff/leaprc.ff03', 'parm03.mol2', 'parm99.dat',
    ['frcmod.ff03'])

ff03_r1 = force_field('leaprc.protein.ff03.r1', 'parm03_r1.mol2',
    'parm99.dat', ['frcmod.ff03'])

ff10 = force_field('oldff/leaprc.ff10', 'parm10.mol2', 'parm10.dat')

ff14ipq = force_field('oldff/leaprc.ff14ipq', 'parm14ipq.mol2',
    'parm14ipq.dat', ['frcmod.tip4pew'])

ff14SB = force_field('oldff/leaprc.ff14SB', 'parm12.mol2', 'parm10.dat',
    ['frcmod.ff14SB'])

# New force fields
ff14SB_redq = force_field('leaprc.ff14SB.redq', 'parm12_redq.mol2', 'parm10.dat',
    ['frcmod.ff14SB'])

ff14SBonlysc = force_field('leaprc.protein.ff14SBonlysc', 'parm12.mol2', 'parm10.dat',
    ['frcmod.ff14SB', 'frcmod.ff99SB14'])

ff19SB = force_field('leaprc.protein.ff19SB', 'parm19.mol2', 'parm19.dat',
    ['frcmod.ff19SB'])

ff15ipq = force_field('leaprc.protein.ff15ipq', 'parm15ipq_10.0.mol2',
    'parm15ipq_10.3.dat')

ff15ipq_vac = force_field('leaprc.protein.ff15ipq-vac', 'parm15ipq-vac_10.0.mol2',
    'parm15ipq_10.3.dat')

fb15 = force_field('leaprc.protein.fb15', 'parm_fb15.mol2', 'parm99.dat',
    ['frcmod.fb15', 'frcmod.tip3pfb'])

FF_DICT = {'ff94': ff94, 'ff99': ff99, 'ff99SB': ff99SB, 'ff03': ff03,
    'ff03.r1': ff03_r1, 'ff10': ff10, 'ff14ipq': ff14ipq,'ff14SB': ff14SB,
    'ff14SB.redq': ff14SB_redq, 'ff14SBonlysc': ff14SBonlysc, 'ff19SB': ff19SB,
    'ff15ipq': ff15ipq, 'ff15ipq-vac': ff15ipq_vac, 'fb15': fb15}

#-----------------------------------------------------------------------------
# About the force field lib parameters
#-----------------------------------------------------------------------------

def get_lib_dict(ff_choice):

    if ff_choice in list(FF_DICT.keys()):
        mol, atids, resids = get_atominfo(FF_DICT[ff_choice].mol2f)
    else:
        mol, atids, resids = get_atominfo(ff_choice)

    libdict = {} #resname + atname : atom type, atom charge
    chargedict = {} #resname : charge

    for i in resids:
        charge = 0.0
        for j in mol.residues[i].resconter:
            #key is residue name and atom name
            key =  mol.residues[i].resname + '-' +  mol.atoms[j].atname
            if len(mol.atoms[j].atomtype) == 1:
                mol.atoms[j].atomtype = mol.atoms[j].atomtype + ' '
            #value is atom type and atom charge
            val =  (mol.atoms[j].atomtype, mol.atoms[j].charge)
            libdict[key] = val
            #cummulative charge of the residue
            charge = charge + mol.atoms[j].charge
        chargedict[mol.residues[i].resname] = charge

        #Alias HN as H
        #atnames = [mol.atoms[k].atname for k in mol.residues[i].resconter]
        #if set(['H', 'N', 'C', 'O']) < set(atnames):
        #  libdict[mol.residues[i].resname + '-HN'] = \
        #  libdict[mol.residues[i].resname + '-H']
    return libdict, chargedict

#-----------------------------------------------------------------------------
# About the force field params parameters
#-----------------------------------------------------------------------------

class Parms:
    def __init__(self, mass, bond, ang, dih, imp, nb, ljed):
        self.mass = mass
        self.bond = bond
        self.ang = ang
        self.dih = dih
        self.imp = imp
        self.nb = nb
        self.ljed = ljed

    def combine(self, Parms2):

        #Mass
        self.mass.update(Parms2.mass)

        #Bond
        for i in list(Parms2.bond.keys()):
            if (i in list(self.bond.keys())) or (i[::-1] in list(self.bond.keys())):
                self.bond[i] = Parms2.bond[i]
            else:
                self.bond[i] = Parms2.bond[i]

        #Angle
        for i in list(Parms2.ang.keys()):
            if (i in list(self.ang.keys())) or (i[::-1] in list(self.ang.keys())):
                self.ang[i] = Parms2.ang[i]
            else:
                self.ang[i] = Parms2.ang[i]

        #Dih
        for i in list(Parms2.dih.keys()):
            if (i in list(self.dih.keys())) or (i[::-1] in list(self.dih.keys())):
                self.dih[i] = Parms2.dih[i]
            else:
                self.dih[i] = Parms2.dih[i]

        #Imp
        for i in list(Parms2.imp.keys()):
            if (i in list(self.imp.keys())) or ((i[0], i[3], i[2], i[1]) in list(self.imp.keys())) \
              or ((i[1], i[0], i[2], i[3]) in list(self.imp.keys())) \
              or ((i[1], i[3], i[2], i[0]) in list(self.imp.keys())) \
              or ((i[3], i[0], i[2], i[1]) in list(self.imp.keys())) \
              or ((i[3], i[1], i[2], i[0]) in list(self.imp.keys())):
                self.imp[i] = Parms2.imp[i]
            else:
                self.imp[i] = Parms2.imp[i]

        #Nb
        self.nb.update(Parms2.nb)

        #LJEd
        for i in list(Parms2.ljed.keys()):
            if (i in list(self.ljed.keys())) or (i[::-1] in list(self.ljed.keys())):
                self.ljed[i] = Parms2.ljed[i]
            else:
                self.ljed[i] = Parms2.ljed[i]

def readmass(massparms, line):
    attyp = line[0:2]
    mass = line[2:]
    massparms[attyp] = mass
    return massparms

def readbond(bondparms, line):
    at1 = line[0:2]
    at2 = line[3:5]
    bondparm = line[5:]
    bondparms[(at1, at2)] = bondparm
    return bondparms

def readang(angparms, line):
    at1 = line[0:2]
    at2 = line[3:5]
    at3 = line[6:8]
    angparm = line[8:]
    angparms[(at1, at2, at3)] = angparm
    return angparms

def readdih(dihparms, line):
    at1 = line[0:2]
    at2 = line[3:5]
    at3 = line[6:8]
    at4 = line[9:11]
    dihtyp = (at1, at2, at3, at4)

    terms = re.search(_dihre, line)
    sdihtyp, n, vn, pha, pero, annot = \
        [t(s) for t,s in zip((str, int, float, float, float, str), terms.groups())]

    nvnp = ' ' + str(n).rjust(4) + ' ' + str(round(vn, 5)).rjust(10) + ' ' + \
           str(round(pha, 5)).rjust(10) + ' '
    pero = int(round(pero, 0))

    dihparm = [nvnp, pero, annot]

    if dihtyp in list(dihparms.keys()):
        has_pero = dihparms[dihtyp][1::3]
        if dihparm[1] not in has_pero:
            dihparm = dihparms[dihtyp] + dihparm
    elif dihtyp[::-1] in list(dihparms.keys()):
        has_pero = dihparms[dihtyp[::-1]][1::3]
        if dihparm[1] not in has_pero:
            dihparm = dihparms[dihtyp[::-1]] + dihparm
    dihparms[dihtyp] = dihparm

    return dihparms

def readimp(impparms, line):
    at1 = line[0:2]
    at2 = line[3:5]
    at3 = line[6:8]
    at4 = line[9:11]
    impparm = line[11:]
    if (at1, at2, at3, at4) in list(impparms.keys()):
        print(line[0:11])

    impparms[(at1, at2, at3, at4)] = impparm
    return impparms

def readeqnb(eqdict, line):
    eqatms = line.split()
    for i in eqatms:
        if len(i) == 1:
            i = i + ' '
    eqdict[eqatms[0]] = eqatms[1:]
    return eqdict

def readnb(nbparms, line):
    at1 = line.split()[0]
    if len(at1) == 1:
        at1 = at1 + ' '
    nbparm = line[2:]
    nbparms[at1] = nbparm
    return nbparms

def readljed(ljedparms, line):
    terms = re.search(_ljedre, line)
    at1, at2, r1, e1, r2, e2, annot = \
        [t(s) for t,s in zip((str, str, str, str, str, str, str), terms.groups())]

    at1 = at1.strip()
    at2 = at2.strip()
    ljedparm = r1 + e1 + r2 + e2 + annot
    if len(at1) == 1:
        at1 = at1 + ' '
    if len(at2) == 1:
        at2 = at2 + ' '
    ljedparms[(at1, at2)] = ljedparm
    return ljedparms

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------

def read_dat_file(datf):

    #Read the parameter into dicts
    massparms = {}
    bondparms = {}
    angparms = {}
    dihparms = {}
    impparms = {}
    nbparms = {}
    ljedparms = {}
    hasljed = False

    count = len(linecache.getlines(datf))
    for i in range(2, count+1):

        rline = linecache.getline(datf, i)
        line = rline.strip()

        if (line[0:4] == 'MOD4') and (line.split()[1] == 'RE'):
            nbbln = i
        if (line[0:4] == 'LJED'):
            hasljed = True
            ljedbln = i
        if line[0:3] == 'END':
            nbeln = i
            break

    massln = []
    bondln = []
    angln = []
    dihln = []
    impln = []
    
    # Read from mass to improper parameters
    for i in range(2, nbbln):

        rline = linecache.getline(datf, i)
        line = rline.strip()

        mass_match = _massre.match(line)
        bond_match = _bondre.match(line)
        ang_match = _angre.match(line)
        dih_match = _dihre.match(line)
        imp_match = _impre.match(line)

        if line:       
            if mass_match:
                massparms = readmass(massparms, line)
                massln.append(i)
            elif bond_match:
                bondparms = readbond(bondparms, line)
                bondln.append(i)
            elif ang_match:
                angparms = readang(angparms, line)
                angln.append(i)
            elif dih_match:
                dihparms = readdih(dihparms, line)
                dihln.append(i)
            elif imp_match:
                impparms = readimp(impparms, line)
                impln.append(i)

    massln0, massln1 = (min(massln), max(massln))
    bondln0, bondln1 = (min(bondln), max(bondln))
    angln0, angln1 = (min(angln), max(angln))
    dihln0, dihln1 = (min(dihln), max(dihln))
    impln0, impln1 = (min(impln), max(impln))

    if (bondln0 - massln1 == 3) and (angln0 - bondln1 == 2) \
        and (dihln0 - angln1 == 2) and (impln0 - dihln1 == 2) and (nbbln - impln1 > 2):
        pass
    else:
        raise pymsmtError('Error of reading the .dat file! Please check it whether ' 
            'it has different parameter types mixed in one section!')

    # Read the NB
    if hasljed is True:
         for i in range(nbbln+1, ljedbln):
            rline = linecache.getline(datf, i)
            line = rline.strip()
            if line:
                nbparms = readnb(nbparms, line)
         for i in range(ljedbln+1, nbeln):
            rline = linecache.getline(datf, i)
            line = rline.strip()
            if line:
                ljedparms = readljed(ljedparms, line)
    else:
        for i in range(nbbln+1, nbeln):
            rline = linecache.getline(datf, i)
            line = rline.strip()
            if line:
                nbparms = readnb(nbparms, line)

    # Deal with the equil atoms
    eqdict = {}
    for i in range(nbbln-3, nbbln):
        rline = linecache.getline(datf, i)
        line = rline.strip()
        if line and rline[0] != ' ':
            eqdict = readeqnb(eqdict, line)

    for i in list(eqdict.keys()):
        for j in eqdict[i]:
            if len(i) == 1:
                nbparms[j] = nbparms[i + ' ']
            else:
                nbparms[j] = nbparms[i]

    # Merge all the parameters into one dict
    parmdict = Parms(massparms, bondparms, angparms, dihparms, impparms, nbparms, ljedparms)
    linecache.clearcache()

    return parmdict

def read_frcmod_file(frcmodf):

    #Get range of each parameter part in the frcmodf
    rfrcmodf = open(frcmodf, 'r')
    cardlist = ['MASS', 'BOND', 'ANGL', 'DIHE', 'IMPR', 'NONB', 'LJED']
    lnlist1 = []
    lnlist2 = []
    ln = 1
    for line in rfrcmodf:
        for card in cardlist:
            if line[0:len(card)] == card:
                lnlist1.append(card)
                lnlist2.append(ln)
        ln = ln + 1
    tln = ln - 1 # Terminal line number
    rfrcmodf.close()

    lndict = {}
    for i in range(0, len(lnlist1)-1):
        lndict[lnlist1[i]] = (lnlist2[i]+1, lnlist2[i+1])

    lndict[lnlist1[-1]] = (lnlist2[-1]+1, tln)

    #Read the parameter into dicts
    massparms = {}
    bondparms = {}
    angparms = {}
    dihparms = {}
    impparms = {}
    nbparms = {}
    ljedparms = {}

    for i in list(lndict.keys()):
        if i == "MASS":
            for j in range(lndict[i][0],lndict[i][1]):
                rline = linecache.getline(frcmodf, j)
                line = rline.strip()
                if line:
                    massparms = readmass(massparms, line)
        elif i == "BOND":
            for j in range(lndict[i][0], lndict[i][1]):
                rline = linecache.getline(frcmodf, j)
                line = rline.strip()
                if line:
                    bondparms = readbond(bondparms, line)
        elif i == "ANGL":
            for j in range(lndict[i][0], lndict[i][1]):
                rline = linecache.getline(frcmodf, j)
                line = rline.strip()
                if line:
                    angparms = readang(angparms, line)
        elif i == "DIHE":
            for j in range(lndict[i][0], lndict[i][1]):
                rline = linecache.getline(frcmodf, j)
                line = rline.strip()
                if line:
                    dihparms = readdih(dihparms, line)
        elif i == "IMPR":
            for j in range(lndict[i][0], lndict[i][1]):
                rline = linecache.getline(frcmodf, j)
                line = rline.strip()
                if line:
                    impparms = readimp(impparms, line)
        elif i == "NONB":
            for j in range(lndict[i][0], lndict[i][1]):
                rline = linecache.getline(frcmodf, j)
                line = rline.strip()
                if line:
                    nbparms = readnb(nbparms, line)
        elif i == "LJED":
            for j in range(lndict[i][0], lndict[i][1]):
                rline = linecache.getline(frcmodf, j)
                line = rline.strip()
                if line:
                    ljedparms = readljed(ljedparms, line)

    parmdict = Parms(massparms, bondparms, angparms, dihparms, impparms, nbparms, ljedparms)
    linecache.clearcache()

    return parmdict

def get_parm_dict(ff_choice, gaff, frcmodfs):

    #1. Read the parm*.dat file
    parmdict = read_dat_file(FF_DICT[ff_choice].datf)

    #2. Read the frcmod file for each force field
    for i in FF_DICT[ff_choice].frcmodfs:
        parmdict1 = read_frcmod_file(i)
        parmdict.combine(parmdict1)

    #3. GAFF
    if gaff == 1:
        parmf2 = parmadd + 'gaff.dat'
        parmdict2 =  read_dat_file(parmf2)
        parmdict.combine(parmdict2)
    elif gaff == 2:
        parmf2 = parmadd + 'gaff2.dat'
        parmdict2 =  read_dat_file(parmf2)
        parmdict.combine(parmdict2)

    #4. Additional frcmod file
    for i in frcmodfs:
        parmdict3 = read_frcmod_file(i)
        parmdict.combine(parmdict3)

    return parmdict

#------------------------------------------------------------------------------
# About Empirical Fitting
#------------------------------------------------------------------------------

def expf(x, a, b, c):
    return a * numpy.exp(-b * x) + c

def getfc(fname, dis):

    lengthl = []
    fcl = []
    fcf = open(fname, 'r')
    for line in fcf:
        length, fc = line.split()[:2]
        length = float(length)
        lengthl.append(length)
        fc = float(fc)
        fcl.append(fc)
    fcf.close()

    initguess = [5E6, 5.5, 0.0]
    lengthl = numpy.array(lengthl)
    fcl = numpy.array(fcl)

    optparas, convar = curve_fit(expf, lengthl, fcl, p0=initguess, maxfev=10000)
    a, b, c = optparas

    val = expf(dis, a, b, c)
    val = round(val, 1)
    return val
