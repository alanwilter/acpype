"""
This module is written for reading the atom and bond information from mol2
file.
"""

from pymsmt.mol.mol import Atom, Residue, Molecule
from pymsmt.mol.element import ionnamel, METAL_PDB
from pymsmt.exp import *
import sys
import linecache

def get_atominfo(fname):

    #Detect the line numbers of each part information
    fp = open(fname, 'r')
    lnum = 1
    for line in fp:
        if ("@<TRIPOS>ATOM" in line):
            atbgin = lnum + 1
        elif ("@<TRIPOS>BOND" in line):
            atend = lnum
        lnum = lnum + 1
    fp.close()

    Atoms = {}
    Residues = {}

    atids = []
    resids = []
    resnamedict = {}
    conterdict = {}

    for i in range(atbgin, atend):
        atid, atname, crdx, crdy, crdz, atomtype, resid, resname, charge = \
        linecache.getline(fname, i).split()[:9]

        #for atom part
        gtype = "ATOM"
        atid = int(atid)
        atids.append(atid)
        crd = (float(crdx),float(crdy),float(crdz))
        charge = float(charge)
        resid = int(resid)

        if (resname, atname) in list(METAL_PDB.keys()):
            element = METAL_PDB[(resname, atname)][0]
        elif atname[0:2].upper() in ['CL', 'BR']:
            element = atname[0].upper() + atname[1].lower()
        else:
            element = atname[0]

        if atid not in list(Atoms.keys()):
            Atoms[atid] = Atom(gtype, atid, atname, element, atomtype, crd, charge, resid, resname)
        else:
            raise pymsmtError('There are more than one atom with atom id '
                              '%d in the mol2 file : %s .' %(atid, fname))

        #for the residue part
        if resid not in resids:
            resids.append(resid)
        if resid not in list(resnamedict.keys()):
            resnamedict[resid] = resname

    #clean the memory
    linecache.clearcache()

    resids.sort()

    for i in resids:
        preconter = []
        for j in atids:
            if (Atoms[j].resid == i) and (j not in preconter):
                preconter.append(j)
        preconter.sort()
        conterdict[i] = preconter

    for i in resids:
        resname = resnamedict[i]
        resconter = conterdict[i]
        Residues[i] = Residue(i, resname, resconter)

    del resnamedict
    del conterdict

    mol = Molecule(Atoms, Residues)

    return mol, atids, resids

def get_bondinfo(fname):
    fp = open(fname, 'r')
    lnum = 1
    for line in fp:
        if ("@<TRIPOS>BOND" in line):
            bdbgin = lnum + 1
        elif ("@<TRIPOS>SUBSTRUCTURE" in line):
            bdend = lnum
        lnum = lnum + 1
    fp.close()

    blist = []

    for i in range(bdbgin, bdend):
        a, b, c, d = linecache.getline(fname, i).split()[:4]

        if (int(b) > int(c)):
            blist.append((int(c), int(b), d))
        else:
            blist.append((int(b), int(c), d))

    linecache.clearcache()
    return blist


