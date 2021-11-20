"""
This module is written for reading the atom and bond information from car
file.
"""

from pymsmt.mol.mol import Atom, Residue, Molecule
from pymsmt.mol.element import ionnamel, METAL_PDB
from pymsmt.mol.getlist import get_blist
from pymsmt.exp import *
import sys
import linecache

def read_carf(fname):

    end = 0
    pbc = 0
    pbc_size = ' '
    #Detect the line numbers of each part information
    fp = open(fname, 'r')
    lnum = 1
    for line in fp:
        if ("PBC=YES" in line) or ("PBC=ON" in line):
            pbc = 1
        if pbc == 1:
            if ("PBC " in line):
                atbgin = lnum + 1
                pbc_size = line.lstrip("PBC")
                pbc_size = pbc_size.strip('\n')
                pbc_size = pbc_size.replace('(', '')
                pbc_size = pbc_size.replace(')', '')
        else:
            if ("!DATE" in line):
                atbgin = lnum + 1
        if ("end" in line) and end == 0:
            end = 1
            atend = lnum
        lnum = lnum + 1
    fp.close()

    Atoms = {}
    Residues = {}

    atids = []
    resids = []
    resnamedict = {}
    conterdict = {}

    #print(atbgin, atend)

    atid = 0
    for i in range(atbgin, atend):
        atname, crdx, crdy, crdz, resname, resid, atomtype, element, charge = \
        linecache.getline(fname, i).split()[:9]

        #for atom part
        gtype = "ATOM"
        atid = atid + 1
        atids.append(atid)
        crd = (float(crdx),float(crdy),float(crdz))
        charge = float(charge)
        resid = int(resid)

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

    return mol, atids, resids, pbc_size

def print_mol2f(mol, atids, outfile, resname2, attyp_dict):

    #iddict1: atom id
    #mol: atname, crd
    #sddict: atomtype
    #stdict: atom charge
    #blist_each: bond information

    print('***Generating the ' + outfile + ' file...')

    blist = get_blist(mol, atids)

    mol2f = open(outfile, 'w')

    ##1. molecule information
    print("@<TRIPOS>MOLECULE", file=mol2f)
    print(resname2, file=mol2f)
    print('%5d%6d%6d%6d%6d' %(len(atids), len(blist),
                                       1, 0, 0), file=mol2f) #atom number and bond number
    print('SMALL', file=mol2f)
    print('Car File Charge', file=mol2f)
    print(' ', file=mol2f)
    print(' ', file=mol2f)

    ##2. atom information
    print('@<TRIPOS>ATOM', file=mol2f)
    for atm in atids:
        #new atom id
        atid = mol.atoms[atm].atid
        atname = mol.atoms[atm].atname
        crd = mol.atoms[atm].crd
        atomtype = mol.atoms[atm].atomtype
        if atomtype in list(attyp_dict.keys()):
            atomtype = attyp_dict[atomtype]
        chg = mol.atoms[atm].charge
        resid = mol.atoms[atm].resid
        print('%7d %-4s    %10.4f%10.4f%10.4f %-4s %6d %-4s %12.6f'\
                        %(atid, atname, crd[0], crd[1], crd[2], atomtype,
                         resid, resname2, chg), file=mol2f)

    ##3. bond information
    print('@<TRIPOS>BOND', file=mol2f)
    for bonds in range(0, len(blist)):
        print('%6d%5d%5d%2d' %(bonds+1, blist[bonds][0],
                        blist[bonds][1], 1), file=mol2f) #all as single bond

    ##4. substructure information
    print('@<TRIPOS>SUBSTRUCTURE', file=mol2f)
    print('     1', resname2, '        1 TEMP' + \
                    '              0 ****  ****    0 ROOT', file=mol2f)
    mol2f.close()

def print_pdbf(mol, atids, fname, pbc_size):

    print('***Generating the ' + fname + ' file...')

    wf = open(fname, 'w')
    print('REMARK, BUILD BY MCPB.PY', file=wf)

    if pbc_size != ' ':
        print('CRYST1' + pbc_size, file=wf)

    for i in atids:
        atm = mol.atoms[i]
        gtype = atm.gtype
        atid = atm.atid
        if len(atm.atname) == 3:
            atname = atm.atname
        else:
            atname = atm.atname.center(4)
        crd = atm.crd
        resid = atm.resid
        resname = atm.resname
        print("%-6s%5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f" \
                 %(gtype, atid, atname, resname, 'A', resid, crd[0], crd[1], crd[2], 1.00, 0.00), file=wf)
    print('END', file=wf)
    wf.close()


