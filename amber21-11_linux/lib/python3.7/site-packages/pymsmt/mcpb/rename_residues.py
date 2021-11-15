"""
This module is written for detecting the disulfide bond and renaming the
residues.
"""

from pymsmt.mol.cal import calc_bond
from pymsmt.mol.mol import get_reslist

def get_diS_bond(mol, atids):

    #Residue IDs
    resids = []
    for i in atids:
        if mol.atoms[i].resid not in resids:
            resids.append(mol.atoms[i].resid)

    disul = []
    for i in resids:
        if mol.residues[i].resname in ['CYX', 'CYS']:
            #for atom in CYX or CYS residue
            for j in mol.residues[i].resconter:
                if mol.atoms[j].atname == 'SG':
                    sgcrd = mol.atoms[j].crd

                    #for every atom
                    for k in atids:
                        if mol.atoms[k].resname in ['CYX', 'CYS']:
                            if mol.atoms[k].atname == 'SG' and k != j:
                                atkcrd = mol.atoms[k].crd
                                dis = calc_bond(sgcrd, atkcrd)
                                if dis <= 2.50:
                                    if (j < k) and ((j, k) not in disul):
                                        disul.append((j, k))
                                    elif ((k, j) not in disul):
                                        disul.append((k, j))

    return disul

def rename_res(mol, atids):

    #Residue IDs
    resids = []
    for i in atids:
        if mol.atoms[i].resid not in resids:
            resids.append(mol.atoms[i].resid)

    #Correct the names of the HIS, ASP, GLU, LYS, CYS
    for i in resids:
        #HIS
        if mol.residues[i].resname == 'HIS':
            hasatoms = []
            for j in mol.residues[i].resconter:
                atname = mol.atoms[j].atname
                hasatoms.append(atname)
            if ('HD1' in hasatoms) and ('HE2' in hasatoms):
                mol.residues[i].resname = 'HIP'
            elif ('HD1' in hasatoms):
                mol.residues[i].resname = 'HID'
            elif ('HE2' in hasatoms):
                mol.residues[i].resname = 'HIE'
        #ASP
        elif mol.residues[i].resname == 'ASP':
            hasatoms = []
            for j in mol.residues[i].resconter:
                atname = mol.atoms[j].atname
                hasatoms.append(atname)
            if ('HD1' in hasatoms) or ('HD2' in hasatoms):
                mol.residues[i].resname = 'ASH'
        #GLU
        elif mol.residues[i].resname == 'GLU':
            hasatoms = []
            for j in mol.residues[i].resconter:
                atname = mol.atoms[j].atname
                hasatoms.append(atname)
            if ('HE1' in hasatoms) or ('HE2' in hasatoms):
                mol.residues[i].resname = 'GLH'
        #LYS
        elif mol.residues[i].resname == 'LYS':
            hasatoms = []
            for j in mol.residues[i].resconter:
                atname = mol.atoms[j].atname
                hasatoms.append(atname)
            if ('HZ1' not in hasatoms):
                mol.residues[i].resname = 'LYN'
        #CYS
        elif mol.residues[i].resname == 'CYS':
            hasatoms = []
            for j in mol.residues[i].resconter:
                atname = mol.atoms[j].atname
                hasatoms.append(atname)
            if ('HG' not in hasatoms): ##There are two different situations
                #for atom in CYS residue
                for j in mol.residues[i].resconter:
                    if mol.atoms[j].atname == 'SG':
                        sgcrd = mol.atoms[j].crd
                        #for every atom
                        for k in atids:
                            if mol.atoms[k].atname == 'SG' and k != j:
                                atkcrd = mol.atoms[k].crd
                                dis = calc_bond(sgcrd, atkcrd)
                                if dis <= 2.50:
                                    mol.residues[i].resname = 'CYX'
                                else:
                                    mol.residues[i].resname = 'CYM'

    reslist = get_reslist(mol, resids)

    #rename the HN atom to H atom in amino acid residues
    for i in resids:
        if i in reslist.std:
            for j in mol.residues[i].resconter:
                if mol.atoms[j].atname == 'HN':
                    mol.atoms[j].atname = 'H'

    return mol
