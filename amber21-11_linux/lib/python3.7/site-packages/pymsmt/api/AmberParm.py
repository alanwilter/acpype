
from parmed.amber import AmberParm
from pymsmt.mol.element import AtnumRev
from pymsmt.mol.rstio import read_rstf
from pymsmt.mol.mol import *

def read_amber_prm(pfile, cfile):

    prmtop = AmberParm(pfile)
    crds = read_rstf(cfile)

    atids = list(range(1, len(prmtop.parm_data['ATOM_NAME'])+1))
    resids = list(range(1, len(prmtop.parm_data['RESIDUE_LABEL'])+1))

    if len(prmtop.parm_data['ATOM_NAME']) != len(crds):
        raise ReadError('The toplogy and coordinates file are not \
                          consistent in the atom numbers.')

    residues = {}
    atoms = {}
    for i in range(0, len(prmtop.parm_data['RESIDUE_LABEL'])):

        resid = i + 1
        resname = prmtop.parm_data['RESIDUE_LABEL'][i]

        if i < len(prmtop.parm_data['RESIDUE_LABEL'])-1:
            resconter = list(range(prmtop.parm_data['RESIDUE_POINTER'][i], \
                         prmtop.parm_data['RESIDUE_POINTER'][i+1]))
        else:
            resconter = list(range(prmtop.parm_data['RESIDUE_POINTER'][i], \
                         len(prmtop.parm_data['ATOM_NAME'])+1))

        residues[resid] = Residue(resid, resname, resconter)

        for j in resconter:
            gtype = 'ATOM'
            atid = j
            atname = prmtop.parm_data['ATOM_NAME'][j-1]
            element = prmtop.parm_data['ATOMIC_NUMBER'][j-1]
            atomtype = prmtop.parm_data['AMBER_ATOM_TYPE'][j-1]
            crd = crds[j-1]
            charge = prmtop.parm_data['CHARGE'][j-1]

            try:
                element = AtnumRev[element]
            except:
                print(resname, atname, atomtype)
                element = 'X'

            atoms[j] = Atom(gtype, atid, atname, element, atomtype, crd, charge, \
                           resid, resname)

    mol = Molecule(atoms, residues)

    return prmtop, mol, atids, resids
