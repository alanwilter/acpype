"""
This is the code for reading and writting pdb files.
"""

from pymsmt.mol.mol import Atom, Residue, Molecule, get_reslist
from pymsmt.mol.getlist import get_pure_type, get_pure_num
from pymsmt.mol.element import ionnamel, CoRadiiDict, METAL_PDB
from pymsmt.exp import *

def get_atominfo_fpdb(fname):
    Atoms = {}
    Residues = {}

    atids = []
    resids = []
    resnamedict = {}
    conterdict = {}

    fp = open(fname, 'r')

    for line in fp:
        if (line[0:4] == "ATOM") or (line[0:6] == "HETATM"):
            gtype = line[0:6].strip(" ")
            atid = int(line[6:11])
            atids.append(atid)
            atname = line[12:16].strip(" ")
            allocind = line[16:17]
            resname = line[17:20].strip(" ")
            chainid = line[21:22]
            resid = int(line[22:26])
            codeinsert = line[26:27]
            crdx = float(line[30:38])
            crdy = float(line[38:46])
            crdz = float(line[46:54])
            crd = (crdx, crdy, crdz)
            occp = line[54:60]
            tempfac = line[60:66]
            atomtype = line[76:78].strip(" ")
            charge = line[78:80]

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
                                  '%d in the PDB file : %s .' %(atid, fname))

            if resid not in resids:
                resids.append(resid)
            if resid not in list(resnamedict.keys()):
                resnamedict[resid] = resname

    fp.close()

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

def writepdb(mol, atids, fname):

    wf = open(fname, 'w')
    print('REMARK, BUILD BY MCPB.PY', file=wf)

    resids = []
    for i in atids:
        if mol.atoms[i].resid not in resids:
            resids.append(mol.atoms[i].resid)

    reslist = get_reslist(mol, resids)

    terlist = reslist.nterm + reslist.nonstd + reslist.water

    terlist = list(set(terlist))

    for i in resids:
        if (i in terlist) and (i != min(resids)):
            print('TER', file=wf)
        for j in mol.residues[i].resconter:
            atm = mol.atoms[j]
            gtype = atm.gtype
            atid = atm.atid
            if len(atm.atname) == 3:
                atname = atm.atname
            else:
                atname = atm.atname.center(4)
            crd = atm.crd
            resid = atm.resid
            resname = mol.residues[resid].resname
            print("%-6s%5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f" \
                     %(gtype, atid, atname, resname, 'A', resid, crd[0], crd[1], crd[2], 1.00, 0.00), file=wf)
    print('END', file=wf)
    wf.close()

def writepdbatm(pdbatm, fname):
    wf = open(fname, 'a')
    if len(pdbatm.atname) == 3:
        print("%-6s%5d %4s %3s %1s%4d   %8.3f%8.3f%8.3f%6.2f%6.2f" %(pdbatm.tiker, pdbatm.atid,\
          pdbatm.atname, pdbatm.resname, pdbatm.chainid, pdbatm.resid, pdbatm.crdx, pdbatm.crdy, \
          pdbatm.crdz, pdbatm.occp, pdbatm.tempfac), file=wf)
    else:
        print("%-6s%5d %4s %3s %1s%4d   %8.3f%8.3f%8.3f%6.2f%6.2f" %(pdbatm.tiker, pdbatm.atid,\
         pdbatm.atname.center(4), pdbatm.resname, pdbatm.chainid, pdbatm.resid, pdbatm.crdx, pdbatm.crdy, \
         pdbatm.crdz, pdbatm.occp, pdbatm.tempfac), file=wf)
