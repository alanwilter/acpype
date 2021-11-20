"""
This module was used to generate the pdb, fingerprint files of the small,
standard and large models and Gaussian input file of small (for bond and
angle parameter fitting) and large models(for RESP charge fitting).
"""

from pymsmt.mol.pdbio import get_atominfo_fpdb, writepdbatm
from pymsmt.mol.cal import calc_bond
from pymsmt.mol.mol import pdbatm, gauatm, get_reslist
from pymsmt.mol.getlist import get_blist
from pymsmt.mol.element import (Atnum, CoRadiiDict,
                              get_ionljparadict, AtnumRev, bdld)
from pymsmt.mol.gauio import (write_gauatm, write_gauatm_opth, write_gau_optf,
                             write_gau_fcf, write_gau_mkf)
from pymsmt.mol.gmsio import (write_gmsatm, write_gms_optf, write_gms_fcf,
                             write_gms_mkf)
from pymsmt.mol.sqmio import get_crdinfo_from_sqm, write_sqm_optf
from pymsmt.lib import get_lib_dict
import os

H_NAMES = ['HH31', 'HH32', 'HH33']  #hydrogen names for ACE and NME methyl group
SH_NAMES = ['H1', 'H2', 'H3'] #The names of the three Hs in the methyl group
GH_NAMES = ['HA2', 'HA3'] #GLY hydrogen names
SH_NAMES2 = ['H1', 'H2']
BIND_ATOMS = ['N', 'O', 'S', 'F', 'Cl', 'Br', 'I']

def del_files(fnamel):
    for fname in fnamel:
        if os.path.exists(fname):
            os.system("rm %s" %fname)

def count_lines(fname):
    ln = 0
    fp = open(fname, 'r')
    for line in fp:
        ln = ln + 1
    return ln
    fp.close()

#-------------------Get metal center residue names-----------------------------
def get_ms_resnames(pdbfile, ionids, cutoff, addres, addbpairs):

    global BIND_ATOMS

    mol, atids, resids = get_atominfo_fpdb(pdbfile)
    ionids = ionids #metal ion atom id
    metresids = [] #metal ion residue id

    #Get the metal ion id
    for i in ionids:
        resid = mol.atoms[i].resid
        metresids.append(resid)

    msresids = [] #metal site residues
    msresids = msresids + metresids + addres

    #Get the atoms which is in the cutoff of metal ion
    for met in ionids:
        for i in atids:
            if (i != met):
                dis = calc_bond(mol.atoms[met].crd, mol.atoms[i].crd)
                if (dis <= cutoff) and mol.atoms[i].element in BIND_ATOMS:
                    if (mol.atoms[i].resid not in msresids):
                        msresids.append(mol.atoms[i].resid)

    #Add the residue numbers of additional bonded atoms
    for i in addbpairs:
        if (i[0] in ionids) and (mol.atoms[i[1]].resid not in msresids):
            msresids.append(mol.atoms[i[1]].resid)
        elif (i[1] in ionids) and (mol.atoms[i[0]].resid not in msresids):
            msresids.append(mol.atoms[i[0]].resid)

    msresids.sort()

    mcresnames = [] #New names of the metal site residues
    tmpl = []
    for i in msresids:
        if len(mol.residues[i].resname) == 3:
            nresname = mol.residues[i].resname[0:3:2]
        elif len(mol.residues[i].resname) in [1, 2]:
            nresname = mol.residues[i].resname

        counter = 1
        for j in tmpl:
            if j == nresname:
                counter = counter + 1

        tmpl.append(nresname)
        nresname = nresname + str(counter)
        mcresnames.append(nresname)

    #Residue names of metal site ligating groups
    mcresnames0 = [mol.residues[i].resname for i in \
                   list(set(msresids)-set(metresids))]
    return mcresnames0, mcresnames

#--------------------Get metal site bonded atom ids--------------------------
def get_ms_ids(mol, atids, ionids, cutoff, addbpairs):

    global BIND_ATOMS

    bdatmids = []
    bdatnams = []

    #Get the atoms which is in the cutoff of metal ion
    for met in ionids:
        for i in atids:
            if (i != met):
                dis = calc_bond(mol.atoms[met].crd, mol.atoms[i].crd)
                if (dis <= cutoff) and mol.atoms[i].element in BIND_ATOMS:
                    if i not in bdatmids:
                        bdatmids.append(i)
                        bdatnams.append(mol.atoms[i].atname)

    #Add the residue numbers of additional bonded atoms
    for i in addbpairs:
        if (i[0] in ionids):
            bdatmids.append(i[1])
            bdatnams.append(mol.atoms[i[1]].atname)
        elif (i[1] in ionids):
            bdatmids.append(i[0])
            bdatnams.append(mol.atoms[i[0]].atname)

    return bdatmids, bdatnams

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#---------------------Write ACE residue into the PDB file---------------------
def write_ace(mol, i, gatms, pdbf, fpf=None):

    """
    ACE group
    a. C, O atoms are kept.
    b. CA --> CH3
    c. HA/HA2, CB/HA3, N --> Three Hs bonds to CH3 with adapting the \
       bond length.
    d. other atoms are deleted.

    Speical case:
    If there is a PRO was treated as ACE, there will be no influence.
    """

    global H_NAMES, SH_NAMES, GH_NAMES

    print("Creating the residue " + str(i) + '-' + \
          mol.residues[i].resname +  " into ACE...")

    #get the coordinates of the CA atom
    for j in mol.residues[i].resconter:
        atname = mol.atoms[j].atname
        if atname == 'CA':
            cacrd = mol.atoms[j].crd

    #rename the atom names to get the large model
    atnames = []
    hdict = {}
    for j in mol.residues[i].resconter:
        atname = mol.atoms[j].atname
        if (atname in ['HA', 'CB', 'N', 'HA2', 'HA3']):
            atnames.append(atname)

    for j in range(0, len(atnames)):
        hdict[atnames[j]] = H_NAMES[j]

    for j in mol.residues[i].resconter:
        atname = mol.atoms[j].atname
        #C, O atoms will be still while CA is CH3, HA, CB and N
        #are three HH3s

        #Only for the backbone things
        if (atname in ['C', 'O', 'CA', 'HA', 'CB', 'N', 'HA2', 'HA3']):
            if (atname == 'C') or (atname == 'O'):
                crdx = mol.atoms[j].crd[0]
                crdy = mol.atoms[j].crd[1]
                crdz = mol.atoms[j].crd[2]
                element = atname
            elif (atname == 'CA'):
                atname = 'CH3'
                crdx = mol.atoms[j].crd[0]
                crdy = mol.atoms[j].crd[1]
                crdz = mol.atoms[j].crd[2]
                element = 'C'
            else: #left were HA, CB, C, HA2, HA3
                atname = hdict[atname]
                bvec = calc_bond(cacrd, mol.atoms[j].crd)
                crdx = cacrd[0] + bdld['CH'] * (mol.atoms[j].crd[0] - cacrd[0])/bvec
                crdy = cacrd[1] + bdld['CH'] * (mol.atoms[j].crd[1] - cacrd[1])/bvec
                crdz = cacrd[2] + bdld['CH'] * (mol.atoms[j].crd[2] - cacrd[2])/bvec
                element = 'H'

            crdx = round(crdx, 3)
            crdy = round(crdy, 3)
            crdz = round(crdz, 3)

            #gausssian file
            gatms.append(gauatm(element, crdx, crdy, crdz))

            #assign other parameters to it
            tiker = mol.atoms[j].gtype
            atid = mol.atoms[j].atid
            chainid = 'A'
            resid = mol.atoms[j].resid
            resname = 'ACE'
            occp = 1.00
            tempfac = 0.00

            #pdb file
            atmi = pdbatm(tiker, atid, atname, resname, chainid, resid,
                          crdx, crdy, crdz, occp, tempfac)
            writepdbatm(atmi, pdbf)

            #fingerprint file
            if fpf is not None:
                fpff = open(fpf, 'a')
                print(str(resid) + '-' + 'ACE-' + atname, file=fpff)
                fpff.close()

#---------------------Write CH3NH2 residue into the PDB file-------------------
def write_ant(mol, i, gatms, pdbf, fpf=None):

    """
    ANT group
    a. N, H1/HN1, H2/HN2, H3/HN3 atoms are kept.
    b. CA --> CH3
    c. C, HA/HA2, CB/HA3 --> Three Hs bonds to CH3 with adapting the bond length.
    d. other atoms are deleted

    If there is a PRO was treated as ACT, there will be no influence.
    """

    global H_NAMES, SH_NAMES, GH_NAMES

    #If the resname is not PRO, keep the NH2 or NH3 group, delete the sidechain
    #and change CA to CH3, HA, CB and C to HH31, HH32, HH33
    if mol.residues[i].resname != 'PRO':

        print("Creating the residue " + str(i) + '-' + \
              mol.residues[i].resname + " into ANT...")

        #get the coordinates of the CA atom
        for j in mol.residues[i].resconter:
            atname = mol.atoms[j].atname
            if atname == 'CA':
                cacrd = mol.atoms[j].crd

        atnames = []
        hdict = {}

        #rename the atom names to get the large model
        for j in mol.residues[i].resconter:
            atname = mol.atoms[j].atname
            if (atname in ['C', 'HA', 'HA2', 'CB', 'HA3']):
                atnames.append(atname)

        for j in range(0, len(atnames)):
            hdict[atnames[j]] = H_NAMES[j]

        for j in mol.residues[i].resconter:
            atname = mol.atoms[j].atname
            if (atname in ['N', 'H1', 'H2', 'H3', 'HN1', 'HN2', 'HN3', 'CA', 'HA', 'HA2', 'CB', 'HA3']):
                if (atname == 'N'):
                    crdx = mol.atoms[j].crd[0]
                    crdy = mol.atoms[j].crd[1]
                    crdz = mol.atoms[j].crd[2]
                    element = atname
                elif (atname in ['HN1', 'H1']):
                    atname = 'H1'
                    crdx = mol.atoms[j].crd[0]
                    crdy = mol.atoms[j].crd[1]
                    crdz = mol.atoms[j].crd[2]
                    element = 'H'
                elif (atname in ['HN2', 'H2']):
                    atname = 'H2'
                    crdx = mol.atoms[j].crd[0]
                    crdy = mol.atoms[j].crd[1]
                    crdz = mol.atoms[j].crd[2]
                    element = 'H'
                elif (atname in ['HN3', 'H3']):
                    atname = 'H3'
                    crdx = mol.atoms[j].crd[0]
                    crdy = mol.atoms[j].crd[1]
                    crdz = mol.atoms[j].crd[2]
                    element = 'H'
                elif (atname == 'CA'):
                    atname = 'CH3'
                    crdx = mol.atoms[j].crd[0]
                    crdy = mol.atoms[j].crd[1]
                    crdz = mol.atoms[j].crd[2]
                    element = 'C'
                else: #left were C, HA/HA2, CB/HA3
                    atname = hdict[atname]
                    bvec = calc_bond(cacrd, mol.atoms[j].crd)
                    crdx = cacrd[0] + bdld['CH'] * (mol.atoms[j].crd[0] - cacrd[0])/bvec
                    crdy = cacrd[1] + bdld['CH'] * (mol.atoms[j].crd[1] - cacrd[1])/bvec
                    crdz = cacrd[2] + bdld['CH'] * (mol.atoms[j].crd[2] - cacrd[2])/bvec
                    element = 'H'

                crdx = round(crdx, 3)
                crdy = round(crdy, 3)
                crdz = round(crdz, 3)

                #gaussian file
                gatms.append(gauatm(element, crdx, crdy, crdz))

                #assign other parameters to it
                tiker = mol.atoms[j].gtype
                atid = mol.atoms[j].atid
                chainid = 'A'
                resid = mol.atoms[j].resid
                resname = 'ANT'
                occp = 1.00
                tempfac = 0.00

                #pdb file
                atmi = pdbatm(tiker, atid, atname, resname, chainid, resid,
                              crdx, crdy, crdz, occp, tempfac)
                writepdbatm(atmi, pdbf)

                #fingerprint file
                if fpf is not None:
                    fpff = open(fpf, 'a')
                    print(str(resid) + '-' + 'ANT-' + atname, file=fpff)
                    fpff.close()

    #If the resname is PRO, change HA, C to HA2, HA3 and keep the ring
    elif mol.residues[i].resname == 'PRO':

        print("Creating the residue " + str(i) + '-' + \
              mol.residues[i].resname + " into ANT...")

        #get the coordinates of the CA atom
        for j in mol.residues[i].resconter:
            atname = mol.atoms[j].atname
            if atname == 'CA':
                cacrd = mol.atoms[j].crd

        atnames = []
        hdict = {}

        #rename the atom names to get the large model
        for j in mol.residues[i].resconter:
            atname = mol.atoms[j].atname
            if (atname in ['HA', 'C']):
                atnames.append(atname)

        for j in range(0, len(atnames)):
            hdict[atnames[j]] = GH_NAMES[j]

        for j in mol.residues[i].resconter:
            atname = mol.atoms[j].atname
            if (atname != 'O'):
                if (atname in ['HN1', 'H1']):
                    atname = 'H1'
                    crdx = mol.atoms[j].crd[0]
                    crdy = mol.atoms[j].crd[1]
                    crdz = mol.atoms[j].crd[2]
                    element = 'H'
                elif (atname in ['HN2', 'H2']):
                    atname = 'H2'
                    crdx = mol.atoms[j].crd[0]
                    crdy = mol.atoms[j].crd[1]
                    crdz = mol.atoms[j].crd[2]
                    element = 'H'
                elif (atname in ['HN3', 'H3']):
                    atname = 'H3'
                    crdx = mol.atoms[j].crd[0]
                    crdy = mol.atoms[j].crd[1]
                    crdz = mol.atoms[j].crd[2]
                    element = 'H'
                elif (atname in ['C', 'HA']):
                    #left were C, HA
                    atname = hdict[atname]
                    bvec = calc_bond(cacrd, mol.atoms[j].crd)
                    crdx = cacrd[0] + bdld['CH'] * (mol.atoms[j].crd[0] - cacrd[0])/bvec
                    crdy = cacrd[1] + bdld['CH'] * (mol.atoms[j].crd[1] - cacrd[1])/bvec
                    crdz = cacrd[2] + bdld['CH'] * (mol.atoms[j].crd[2] - cacrd[2])/bvec
                    element = 'H'
                else:
                    crdx = mol.atoms[j].crd[0]
                    crdy = mol.atoms[j].crd[1]
                    crdz = mol.atoms[j].crd[2]
                    element = mol.atoms[j].element

                crdx = round(crdx, 3)
                crdy = round(crdy, 3)
                crdz = round(crdz, 3)

                #gaussian file
                gatms.append(gauatm(element, crdx, crdy, crdz))

                #assign other parameters to it
                tiker = mol.atoms[j].gtype
                atid = mol.atoms[j].atid
                chainid = 'A'
                resid = mol.atoms[j].resid
                resname = 'PNT'
                occp = 1.00
                tempfac = 0.00

                #pdb file
                atmi = pdbatm(tiker, atid, atname, resname, chainid, resid,
                              crdx, crdy, crdz, occp, tempfac)
                writepdbatm(atmi, pdbf)

                #fingerprint file
                if fpf is not None:
                    fpff = open(fpf, 'a')
                    print(str(resid) + '-' + 'PNT-' + atname, file=fpff)
                    fpff.close()

#---------------------Write CH3CO2- residue into the PDB file---------------------
def write_act(mol, i, gatms, pdbf, fpf=None):

    """
    ACT group
    a. C, O, OXT atoms are kept.
    b. CA --> CH3
    c. N, HA/HA2, CB/HA3 --> Three Hs bonds to CH3 with adapting the bond length.
    d. other atoms are deleted

    If there is a PRO was treated as ACT, there will be no influence.
    """

    global H_NAMES, SH_NAMES, GH_NAMES

    print("Creating the residue " + str(i) + '-' + \
          mol.residues[i].resname + " into ACT...")

    #get the coordinates of the CA atom
    for j in mol.residues[i].resconter:
        atname = mol.atoms[j].atname
        if atname == 'CA':
            cacrd = mol.atoms[j].crd

    atnames = []
    hdict = {}

    #rename the atom names to get the large model
    for j in mol.residues[i].resconter:
        atname = mol.atoms[j].atname
        if (atname in ['N', 'HA', 'HA2', 'CB', 'HA3']):
            atnames.append(atname)

    for j in range(0, len(atnames)):
        hdict[atnames[j]] = H_NAMES[j]

    for j in mol.residues[i].resconter:
        atname = mol.atoms[j].atname
        if (atname in ['N', 'HA', 'HA2', 'CB', 'HA3', 'C', 'O', 'OXT', 'CA']):
            if (atname == 'C') or (atname == 'O'):
                crdx = mol.atoms[j].crd[0]
                crdy = mol.atoms[j].crd[1]
                crdz = mol.atoms[j].crd[2]
                element = atname
            elif (atname == 'OXT'):
                crdx = mol.atoms[j].crd[0]
                crdy = mol.atoms[j].crd[1]
                crdz = mol.atoms[j].crd[2]
                element = 'O'
            elif (atname == 'CA'):
                atname = 'CH3'
                crdx = mol.atoms[j].crd[0]
                crdy = mol.atoms[j].crd[1]
                crdz = mol.atoms[j].crd[2]
                element = 'C'
            else: #left were HA, CB, C, HA2, HA3
                atname = hdict[atname]
                bvec = calc_bond(cacrd, mol.atoms[j].crd)
                crdx = cacrd[0] + bdld['CH'] * (mol.atoms[j].crd[0] - cacrd[0])/bvec
                crdy = cacrd[1] + bdld['CH'] * (mol.atoms[j].crd[1] - cacrd[1])/bvec
                crdz = cacrd[2] + bdld['CH'] * (mol.atoms[j].crd[2] - cacrd[2])/bvec
                element = 'H'

            crdx = round(crdx, 3)
            crdy = round(crdy, 3)
            crdz = round(crdz, 3)

            #gaussian file
            gatms.append(gauatm(element, crdx, crdy, crdz))

            #assign other parameters to it
            tiker = mol.atoms[j].gtype
            atid = mol.atoms[j].atid
            chainid = 'A'
            resid = mol.atoms[j].resid
            resname = 'ACT'
            occp = 1.00
            tempfac = 0.00

            #pdb file
            atmi = pdbatm(tiker, atid, atname, resname, chainid, resid,
                          crdx, crdy, crdz, occp, tempfac)
            writepdbatm(atmi, pdbf)

            #fingerprint file
            if fpf is not None:
                fpff = open(fpf, 'a')
                print(str(resid) + '-' + 'ACT-' + atname, file=fpff)
                fpff.close()

#---------------------Write NME residue into the PDB file---------------------
def write_nme(mol, i, gatms, pdbf, fpf=None):

    """
    NME group
    a. N, H atoms are kept.
    b. CA --> CH3
    c. HA/HA2, CB/HA3, C --> Three Hs bonds to CH3 with adapting the \
       bond length.
    d. other atoms are deleted.

    If there is a PRO was treated as NME, the atom CD need to change to \
       atom H.
    """

    global H_NAMES, SH_NAMES, GH_NAMES

    print("Creating the residue " + str(i) + '-' + \
          mol.residues[i].resname + " into NME...")

    #If the resname is not PRO
    if mol.residues[i].resname != 'PRO':

        #get the coordinates of the CA atom
        for j in mol.residues[i].resconter:
            atname = mol.atoms[j].atname
            if atname == 'CA':
                cacrd = mol.atoms[j].crd

        atnames = []
        hdict = {}

        #rename the atom names to get the large model
        for j in mol.residues[i].resconter:
            atname = mol.atoms[j].atname
            if (atname in ['HA', 'CB', 'C', 'HA2', 'HA3']):
                atnames.append(atname)

        for j in range(0, len(atnames)):
            hdict[atnames[j]] = H_NAMES[j]

        for j in mol.residues[i].resconter:
            atname = mol.atoms[j].atname
            if (atname in ['N', 'H', 'HN', 'CA', 'HA', 'CB', 'C', 'HA2', 'HA3']):
                #N, H atoms will be still
                #while CA is CH3, HA, CB and C are three HH3s
                if (atname == 'N') or (atname == 'H'):
                    crdx = mol.atoms[j].crd[0]
                    crdy = mol.atoms[j].crd[1]
                    crdz = mol.atoms[j].crd[2]
                    element = atname
                elif (atname == 'HN'):
                    atname = 'H'
                    crdx = mol.atoms[j].crd[0]
                    crdy = mol.atoms[j].crd[1]
                    crdz = mol.atoms[j].crd[2]
                    element = 'H'
                elif (atname == 'CA'):
                    atname = 'CH3'
                    crdx = mol.atoms[j].crd[0]
                    crdy = mol.atoms[j].crd[1]
                    crdz = mol.atoms[j].crd[2]
                    element = 'C'
                else: #left were HA, CB, C, HA2, HA3
                    atname = hdict[atname]
                    bvec = calc_bond(cacrd, mol.atoms[j].crd)
                    crdx = cacrd[0] + bdld['CH'] * (mol.atoms[j].crd[0] - cacrd[0])/bvec
                    crdy = cacrd[1] + bdld['CH'] * (mol.atoms[j].crd[1] - cacrd[1])/bvec
                    crdz = cacrd[2] + bdld['CH'] * (mol.atoms[j].crd[2] - cacrd[2])/bvec
                    element = 'H'

                crdx = round(crdx, 3)
                crdy = round(crdy, 3)
                crdz = round(crdz, 3)

                #gaussian file
                gatms.append(gauatm(element, crdx, crdy, crdz))

                #assign other parameters to it
                tiker = mol.atoms[j].gtype
                atid = mol.atoms[j].atid
                chainid = 'A'
                resid = mol.atoms[j].resid
                resname = 'NME'
                occp = 1.00
                tempfac = 0.00

                #pdb file
                atmi = pdbatm(tiker, atid, atname, resname, chainid, resid,
                              crdx, crdy, crdz, occp, tempfac)
                writepdbatm(atmi, pdbf)

                #fingerprint file
                if fpf is not None:
                    fpff = open(fpf, 'a')
                    print(str(resid) + '-' + 'NME-' + atname, file=fpff)
                    fpff.close()

    #If the resname is PRO
    else:

        #get the coordinates of the CA atom
        for j in mol.residues[i].resconter:
            atname = mol.atoms[j].atname
            if atname == 'CA':
                cacrd = mol.atoms[j].crd
            if atname == 'N':
                ncrd = mol.atoms[j].crd

        atnames = []
        hdict = {}

        #rename the atom names to get the large model
        for j in mol.residues[i].resconter:
            atname = mol.atoms[j].atname
            if (atname in ['HA', 'CB', 'C']):
                atnames.append(atname)

        for j in range(0, len(atnames)):
            hdict[atnames[j]] = H_NAMES[j]

        for j in mol.residues[i].resconter:
            atname = mol.atoms[j].atname
            if (atname in ['N', 'CA', 'HA', 'CB', 'C', 'CD']):
                #N atom will be still, CD will change to H
                #while CA is CH3, HA, CB and C are three HH3s
                if (atname == 'N'):
                    crdx = mol.atoms[j].crd[0]
                    crdy = mol.atoms[j].crd[1]
                    crdz = mol.atoms[j].crd[2]
                    element = atname
                elif (atname == 'CA'):
                    atname = 'CH3'
                    crdx = mol.atoms[j].crd[0]
                    crdy = mol.atoms[j].crd[1]
                    crdz = mol.atoms[j].crd[2]
                    element = 'C'
                elif (atname == 'CD'):
                    atname = 'H'
                    bvec = calc_bond(ncrd, mol.atoms[j].crd)
                    crdx = ncrd[0] + bdld['CH'] * (mol.atoms[j].crd[0] - ncrd[0])/bvec
                    crdy = ncrd[1] + bdld['CH'] * (mol.atoms[j].crd[1] - ncrd[1])/bvec
                    crdz = ncrd[2] + bdld['CH'] * (mol.atoms[j].crd[2] - ncrd[2])/bvec
                    element = 'H'
                else: #left were HA, CB, C
                    atname = hdict[atname]
                    bvec = calc_bond(cacrd, mol.atoms[j].crd)
                    crdx = cacrd[0] + bdld['CH'] * (mol.atoms[j].crd[0] - cacrd[0])/bvec
                    crdy = cacrd[1] + bdld['CH'] * (mol.atoms[j].crd[1] - cacrd[1])/bvec
                    crdz = cacrd[2] + bdld['CH'] * (mol.atoms[j].crd[2] - cacrd[2])/bvec
                    element = 'H'

                crdx = round(crdx, 3)
                crdy = round(crdy, 3)
                crdz = round(crdz, 3)

                #gaussian file
                gatms.append(gauatm(element, crdx, crdy, crdz))

                #assign other parameters to it
                tiker = mol.atoms[j].gtype
                atid = mol.atoms[j].atid
                chainid = 'A'
                resid = mol.atoms[j].resid
                resname = 'NME'
                occp = 1.00
                tempfac = 0.00

                #pdb file
                atmi = pdbatm(tiker, atid, atname, resname, chainid, resid,
                              crdx, crdy, crdz, occp, tempfac)
                writepdbatm(atmi, pdbf)

                #fingerprint file
                if fpf is not None:
                    fpff = open(fpf, 'a')
                    print(str(resid) + '-' + 'NME-' + atname, file=fpff)
                    fpff.close()

#---------------------Write GLY residue into the PDB file---------------------
def write_gly(mol, i, gatms, pdbf, fpf=None):

    """
    GLY group
    a. N, H, C, O, CA are kept.
    b. HA, CB --> Two Hs bond to CA in GLY.
    c. other atoms are deleted.
    If there is a PRO was treated as GLY, the atom CD need to change to \
    atom H as well.
    """

    global H_NAMES, SH_NAMES, GH_NAMES

    #get the coordinates of the CA atom
    print("Creating the residue " + str(i) + '-' + \
          mol.residues[i].resname + " into GLY...")

    #If the resname is not PRO
    if mol.residues[i].resname != 'PRO':

        for j in mol.residues[i].resconter:
            atname = mol.atoms[j].atname
            if atname == 'CA':
                cacrd = mol.atoms[j].crd

        #rename the atom names to get the large model
        atnames = []
        hdict = {}

        for j in mol.residues[i].resconter:
            atname = mol.atoms[j].atname
            if (atname in ['HA', 'CB']):
                atnames.append(atname)

        for j in range(0, len(atnames)):
            hdict[atnames[j]] = GH_NAMES[j]

        for j in mol.residues[i].resconter:
            atname = mol.atoms[j].atname
            if (atname in ['N', 'NH', 'H', 'HN', 'CA', 'HA', 'CB', 'C', 'O','HA2', 'HA3', 'H1', 'H2', 'H3', 'HN1', 'HN2', 'HN3']):
                if atname in ['N', 'H', 'CA', 'C', 'O', 'HA2', 'HA3']:
                    crdx = mol.atoms[j].crd[0]
                    crdy = mol.atoms[j].crd[1]
                    crdz = mol.atoms[j].crd[2]
                    element = mol.atoms[j].element
                elif (atname == 'NH'):
                    atname = 'N'
                    crdx = mol.atoms[j].crd[0]
                    crdy = mol.atoms[j].crd[1]
                    crdz = mol.atoms[j].crd[2]
                    element = 'N'
                elif (atname == 'HN'):
                    atname = 'H'
                    crdx = mol.atoms[j].crd[0]
                    crdy = mol.atoms[j].crd[1]
                    crdz = mol.atoms[j].crd[2]
                    element = 'H'
                elif (atname == 'HN1'):
                    atname = 'H1'
                    crdx = mol.atoms[j].crd[0]
                    crdy = mol.atoms[j].crd[1]
                    crdz = mol.atoms[j].crd[2]
                    element = 'H'
                elif (atname == 'HN2'):
                    atname = 'H2'
                    crdx = mol.atoms[j].crd[0]
                    crdy = mol.atoms[j].crd[1]
                    crdz = mol.atoms[j].crd[2]
                    element = 'H'
                elif (atname == 'HN3'):
                    atname = 'H3'
                    crdx = mol.atoms[j].crd[0]
                    crdy = mol.atoms[j].crd[1]
                    crdz = mol.atoms[j].crd[2]
                    element = 'H'
                elif (atname in ['H1', 'H2', 'H3']):
                    crdx = mol.atoms[j].crd[0]
                    crdy = mol.atoms[j].crd[1]
                    crdz = mol.atoms[j].crd[2]
                    element = 'H'
                else: #the left were HA, CB
                    atname = hdict[atname]
                    bvec = calc_bond(cacrd, mol.atoms[j].crd)
                    crdx = cacrd[0] + bdld['CH'] * (mol.atoms[j].crd[0] - cacrd[0])/bvec
                    crdy = cacrd[1] + bdld['CH'] * (mol.atoms[j].crd[1] - cacrd[1])/bvec
                    crdz = cacrd[2] + bdld['CH'] * (mol.atoms[j].crd[2] - cacrd[2])/bvec
                    element = 'H'

                crdx = round(crdx, 3)
                crdy = round(crdy, 3)
                crdz = round(crdz, 3)

                #gaussian file
                gatms.append(gauatm(element, crdx, crdy, crdz))

                #assign other parameters to it
                tiker = mol.atoms[j].gtype
                atid = mol.atoms[j].atid
                chainid = 'A'
                resid = mol.atoms[j].resid
                resname = 'GLY'
                occp = 1.00
                tempfac = 0.00

                #pdb file
                atmi = pdbatm(tiker, atid, atname, resname, chainid, resid,
                              crdx, crdy, crdz, occp, tempfac)
                writepdbatm(atmi, pdbf)

                #fingerprint file
                if fpf is not None:
                    fpff = open(fpf, 'a')
                    print(str(resid) + '-' + 'GLY-' + atname, file=fpff)
                    fpff.close()

    #If the resname is PRO
    else:

        for j in mol.residues[i].resconter:
            atname = mol.atoms[j].atname
            if atname == 'CA':
                cacrd = mol.atoms[j].crd
            if atname == 'N':
                ncrd = mol.atoms[j].crd

        #rename the atom names to get the large model
        atnames = []
        hdict = {}
        for j in mol.residues[i].resconter:
            atname = mol.atoms[j].atname
            if (atname in ['HA', 'CB']):
                atnames.append(atname)

        for j in range(0, len(atnames)):
            hdict[atnames[j]] = GH_NAMES[j]

        for j in mol.residues[i].resconter:
            atname = mol.atoms[j].atname
            if (atname in ['N', 'HN', 'CD', 'CA', 'HA', 'CB', 'C', 'O', 'H2', 'H3', 'HN2', 'HN3']):
                if atname in ['N', 'CA', 'C', 'O']:
                    crdx = mol.atoms[j].crd[0]
                    crdy = mol.atoms[j].crd[1]
                    crdz = mol.atoms[j].crd[2]
                    element = mol.atoms[j].element
                elif atname in ['HA', 'CB']:
                    atname = hdict[atname]
                    bvec = calc_bond(cacrd, mol.atoms[j].crd)
                    crdx = cacrd[0] + bdld['CH'] * (mol.atoms[j].crd[0] - cacrd[0])/bvec
                    crdy = cacrd[1] + bdld['CH'] * (mol.atoms[j].crd[1] - cacrd[1])/bvec
                    crdz = cacrd[2] + bdld['CH'] * (mol.atoms[j].crd[2] - cacrd[2])/bvec
                    element = 'H'
                elif (atname == 'NH'):
                    atname = 'N'
                    crdx = mol.atoms[j].crd[0]
                    crdy = mol.atoms[j].crd[1]
                    crdz = mol.atoms[j].crd[2]
                    element = 'N'
                elif (atname == 'HN2'):
                    atname = 'H2'
                    crdx = mol.atoms[j].crd[0]
                    crdy = mol.atoms[j].crd[1]
                    crdz = mol.atoms[j].crd[2]
                    element = 'H'
                elif (atname == 'HN3'):
                    atname = 'H3'
                    crdx = mol.atoms[j].crd[0]
                    crdy = mol.atoms[j].crd[1]
                    crdz = mol.atoms[j].crd[2]
                    element = 'H'
                elif atname in ['H2', 'H3']:
                    crdx = mol.atoms[j].crd[0]
                    crdy = mol.atoms[j].crd[1]
                    crdz = mol.atoms[j].crd[2]
                    element = 'H'
                else: #the only left is atom CD
                    atname = 'H'
                    bvec = calc_bond(ncrd, mol.atoms[j].crd)
                    crdx = ncrd[0] + bdld['NH'] * (mol.atoms[j].crd[0] - ncrd[0])/bvec
                    crdy = ncrd[1] + bdld['NH'] * (mol.atoms[j].crd[1] - ncrd[1])/bvec
                    crdz = ncrd[2] + bdld['NH'] * (mol.atoms[j].crd[2] - ncrd[2])/bvec
                    element = 'H'

                crdx = round(crdx, 3)
                crdy = round(crdy, 3)
                crdz = round(crdz, 3)

                #gaussian file
                gatms.append(gauatm(element, crdx, crdy, crdz))

                #assign other parameters to it
                tiker = mol.atoms[j].gtype
                atid = mol.atoms[j].atid
                chainid = 'A'
                resid = mol.atoms[j].resid
                resname = 'GLY'
                occp = 1.00
                tempfac = 0.00

                #pdb file
                atmi = pdbatm(tiker, atid, atname, resname, chainid, resid,
                              crdx, crdy, crdz, occp, tempfac)
                writepdbatm(atmi, pdbf)

                #fingerprint file
                if fpf is not None:
                    fpff = open(fpf, 'a')
                    print(str(resid) + '-' + 'GLY-' + atname, file=fpff)
                    fpff.close()

#---------------------Write normal residue into the PDB file---------------------
def write_normal(mol, reslist, i, gatms, pdbf, fpf=None):

    print("It contains the residue " + str(i) + '-' + \
          mol.residues[i].resname + " as normal.")

    for j in mol.residues[i].resconter:
        tiker = mol.atoms[j].gtype
        atid = mol.atoms[j].atid

        atname = mol.atoms[j].atname
        if (i in reslist.std) and (atname == 'HN'):
            atname = 'H'

        element = mol.atoms[j].element
        chainid = 'A'
        resid = mol.atoms[j].resid
        resname = mol.residues[resid].resname

        crdx = mol.atoms[j].crd[0]
        crdy = mol.atoms[j].crd[1]
        crdz = mol.atoms[j].crd[2]
        occp = 1.00
        tempfac = 0.00

        crdx = round(crdx, 3)
        crdy = round(crdy, 3)
        crdz = round(crdz, 3)

        #Gaussian file
        gatms.append(gauatm(element, crdx, crdy, crdz))

        #PDB file
        atmi = pdbatm(tiker, atid, atname, resname, chainid, resid,
                      crdx, crdy, crdz, occp, tempfac)
        writepdbatm(atmi, pdbf)

        #Fingerprint file
        if fpf is not None:
            fpff = open(fpf, 'a')
            print(str(resid) + '-' + resname + '-' + atname, file=fpff)
            fpff.close()

#-----------------------Write Sidechain residues-------------------------------
def write_sc(mol, i, gatms, smpdbf):

    global H_NAMES, SH_NAMES, GH_NAMES

    print("It contains the residue " + str(i) + '-' + \
          mol.residues[i].resname + " as sidechain coordinated.")

    #get the coordinates of the Ca atom
    for j in mol.residues[i].resconter:
        atname = mol.atoms[j].atname
        if atname == 'CA':
            cacrd = mol.atoms[j].crd

    #N, C, HA are three Hs in the small model
    atnames = []
    hdict = {}

    for j in mol.residues[i].resconter:
        atname = mol.atoms[j].atname
        if (atname in ['N', 'C', 'HA']):
            atnames.append(atname)

    for j in range(0, len(atnames)):
        hdict[atnames[j]] = SH_NAMES[j]

    for j in mol.residues[i].resconter:
        atname = mol.atoms[j].atname
        if (atname not in ['H', 'HN', 'O', 'OXT', 'H1', 'H2', 'H3', 'HN1', 'HN2', 'HN3']): #HN is a alias of H
            #These two backbone atoms will be deleted in the sidechain modeling
            resname = mol.residues[i].resname
            if (atname == 'CA'):
                atname = 'CH3'
                crdx = mol.atoms[j].crd[0]
                crdy = mol.atoms[j].crd[1]
                crdz = mol.atoms[j].crd[2]
                element = 'C'
            elif (atname in ['N', 'C', 'HA']):
                atname = hdict[atname]
                bvec = calc_bond(cacrd, mol.atoms[j].crd)
                crdx = cacrd[0] + bdld['CH'] * (mol.atoms[j].crd[0] - cacrd[0])/bvec
                crdy = cacrd[1] + bdld['CH'] * (mol.atoms[j].crd[1] - cacrd[1])/bvec
                crdz = cacrd[2] + bdld['CH'] * (mol.atoms[j].crd[2] - cacrd[2])/bvec
                element = 'H'
            else: #Which is in the sidechain
                crdx = mol.atoms[j].crd[0]
                crdy = mol.atoms[j].crd[1]
                crdz = mol.atoms[j].crd[2]
                element = mol.atoms[j].element

            crdx = round(crdx, 3)
            crdy = round(crdy, 3)
            crdz = round(crdz, 3)

            #Gaussian file
            gatms.append(gauatm(element, crdx, crdy, crdz))

            tiker = mol.atoms[j].gtype
            atid = mol.atoms[j].atid
            chainid = 'A'
            resid = mol.atoms[j].resid
            resname = mol.residues[resid].resname
            occp = 1.00
            tempfac = 0.00

            #PDB file
            atmi = pdbatm(tiker, atid, atname, resname, chainid, resid,
                          crdx, crdy, crdz, occp, tempfac)
            writepdbatm(atmi, smpdbf)

#-----------------------Write Sidechain residues2-------------------------------

#By keeping the N, H group inside the residue, this residue is coordinate
#to the metal ion and it is beside the residue which use backbone Oxygen
#to coordinate with the metal ion

def write_sc_knh(mol, i, gatms, smpdbf):

    global H_NAMES, SH_NAMES, GH_NAMES, SH_NAMES2

    print("It contains the residue " + str(i) + '-' + \
          mol.residues[i].resname + " as keeping sidechain and NH group.")

    #get the coordinates of the Ca atom
    for j in mol.residues[i].resconter:
        atname = mol.atoms[j].atname
        if atname == 'CA':
            cacrd = mol.atoms[j].crd

    #C, HA are two Hs in the small model
    atnames = []
    hdict = {}

    for j in mol.residues[i].resconter:
        atname = mol.atoms[j].atname
        if (atname in ['C', 'HA']):
            atnames.append(atname)

    for j in range(0, len(atnames)):
        hdict[atnames[j]] = SH_NAMES2[j]

    for j in mol.residues[i].resconter:
        atname = mol.atoms[j].atname
        #Get rid of O
        #Change CA to CH3, C, HA to Hs, and keep remaining
        if (atname not in ['O', 'OXT']):
            resname = mol.residues[i].resname
            if (atname == 'CA'):
                atname = 'CH3'
                crdx = mol.atoms[j].crd[0]
                crdy = mol.atoms[j].crd[1]
                crdz = mol.atoms[j].crd[2]
                element = 'C'
            elif (atname in ['C', 'HA']):
                atname = hdict[atname]
                bvec = calc_bond(cacrd, mol.atoms[j].crd)
                crdx = cacrd[0] + bdld['CH'] * (mol.atoms[j].crd[0] - cacrd[0])/bvec
                crdy = cacrd[1] + bdld['CH'] * (mol.atoms[j].crd[1] - cacrd[1])/bvec
                crdz = cacrd[2] + bdld['CH'] * (mol.atoms[j].crd[2] - cacrd[2])/bvec
                element = 'H'
            elif atname == 'HN':
                atname = 'H'
                crdx = mol.atoms[j].crd[0]
                crdy = mol.atoms[j].crd[1]
                crdz = mol.atoms[j].crd[2]
                element = 'H'
            else:
                crdx = mol.atoms[j].crd[0]
                crdy = mol.atoms[j].crd[1]
                crdz = mol.atoms[j].crd[2]
                element = mol.atoms[j].element

            crdx = round(crdx, 3)
            crdy = round(crdy, 3)
            crdz = round(crdz, 3)

            #Gaussian file
            gatms.append(gauatm(element, crdx, crdy, crdz))

            tiker = mol.atoms[j].gtype
            atid = mol.atoms[j].atid
            chainid = 'A'
            resid = mol.atoms[j].resid
            resname = mol.residues[resid].resname
            occp = 1.00
            tempfac = 0.00

            #PDB file
            atmi = pdbatm(tiker, atid, atname, resname, chainid, resid,
                          crdx, crdy, crdz, occp, tempfac)
            writepdbatm(atmi, smpdbf)

#-----------------------Write Sidechain residues3-------------------------------

#By keeping the C, O group inside the residue, this residue is coordinate
#to the metal ion and it is beside the residue which use backbone Nitrogen
#to coordinate with the metal ion

def write_sc_kco(mol, i, gatms, smpdbf):

    global H_NAMES, SH_NAMES, GH_NAMES, SH_NAMES2

    print("It contains the residue " + str(i) + '-' + \
          mol.residues[i].resname + " as keeping sidechain and CO group.")

    #get the coordinates of the Ca atom
    for j in mol.residues[i].resconter:
        atname = mol.atoms[j].atname
        if atname == 'CA':
            cacrd = mol.atoms[j].crd

    #C, HA are two Hs in the small model
    atnames = []
    hdict = {}

    for j in mol.residues[i].resconter:
        atname = mol.atoms[j].atname
        if (atname in ['N', 'HA']):
            atnames.append(atname)

    for j in range(0, len(atnames)):
        hdict[atnames[j]] = SH_NAMES2[j]

    for j in mol.residues[i].resconter:
        atname = mol.atoms[j].atname
        if (atname not in ['H', 'HN', 'H1', 'H2', 'H3', 'HN1', 'HN2', 'HN3']):
        #Get rid of the N and H atoms
        #Change CA to CH3, change N, HA to Hs and keep CB and sidechain
            resname = mol.residues[i].resname
            if (atname == 'CA'):
                atname = 'CH3'
                crdx = mol.atoms[j].crd[0]
                crdy = mol.atoms[j].crd[1]
                crdz = mol.atoms[j].crd[2]
                element = 'C'
            elif (atname in ['N', 'HA']):
                atname = hdict[atname]
                bvec = calc_bond(cacrd, mol.atoms[j].crd)
                crdx = cacrd[0] + bdld['CH'] * (mol.atoms[j].crd[0] - cacrd[0])/bvec
                crdy = cacrd[1] + bdld['CH'] * (mol.atoms[j].crd[1] - cacrd[1])/bvec
                crdz = cacrd[2] + bdld['CH'] * (mol.atoms[j].crd[2] - cacrd[2])/bvec
                element = 'H'
            else:
                crdx = mol.atoms[j].crd[0]
                crdy = mol.atoms[j].crd[1]
                crdz = mol.atoms[j].crd[2]
                element = mol.atoms[j].element

            crdx = round(crdx, 3)
            crdy = round(crdy, 3)
            crdz = round(crdz, 3)

            #Gaussian file
            gatms.append(gauatm(element, crdx, crdy, crdz))

            tiker = mol.atoms[j].gtype
            atid = mol.atoms[j].atid
            chainid = 'A'
            resid = mol.atoms[j].resid
            resname = mol.residues[resid].resname
            occp = 1.00
            tempfac = 0.00

            #PDB file
            atmi = pdbatm(tiker, atid, atname, resname, chainid, resid,
                          crdx, crdy, crdz, occp, tempfac)
            writepdbatm(atmi, smpdbf)

def write_only_c2h(mol, i, gatms, smpdbf):

    print("Only keep the backbone C atom of residue " + str(i) + '-' + \
          mol.residues[i].resname + " as a capped H.")

    for j in mol.residues[i+1].resconter:
        atname = mol.atoms[j].atname
        if atname == 'N':
            ncrd = mol.atoms[j].crd

    for j in mol.residues[i].resconter:
        atname = mol.atoms[j].atname
        if (atname == 'C'):
            atname = 'H*'
            element = 'H'
            bvec = calc_bond(ncrd, mol.atoms[j].crd)
            crdx = ncrd[0] + bdld['NH'] * (mol.atoms[j].crd[0] - ncrd[0])/bvec
            crdy = ncrd[1] + bdld['NH'] * (mol.atoms[j].crd[1] - ncrd[1])/bvec
            crdz = ncrd[2] + bdld['NH'] * (mol.atoms[j].crd[2] - ncrd[2])/bvec

            crdx = round(crdx, 3)
            crdy = round(crdy, 3)
            crdz = round(crdz, 3)

            #Gaussian file
            gatms.append(gauatm(element, crdx, crdy, crdz))

            tiker = mol.atoms[j].gtype
            atid = mol.atoms[j].atid
            chainid = 'A'
            resid = mol.atoms[j].resid
            resname = mol.residues[resid].resname
            occp = 1.00
            tempfac = 0.00

            #PDB file
            atmi = pdbatm(tiker, atid, atname, resname, chainid, resid,
                          crdx, crdy, crdz, occp, tempfac)
            writepdbatm(atmi, smpdbf)

def write_only_n2h(mol, i, gatms, smpdbf):

    print("Only keep the backbone N atom of residue " + str(i) + '-' + \
          mol.residues[i].resname + " as a capped H.")

    for j in mol.residues[i-1].resconter:
        atname = mol.atoms[j].atname
        if atname == 'C':
            ccrd = mol.atoms[j].crd

    for j in mol.residues[i].resconter:
        atname = mol.atoms[j].atname
        if (atname == 'N'):
            atname = 'H#'
            element = 'H'
            bvec = calc_bond(ccrd, mol.atoms[j].crd)
            crdx = ccrd[0] + bdld['CH'] * (mol.atoms[j].crd[0] - ccrd[0])/bvec
            crdy = ccrd[1] + bdld['CH'] * (mol.atoms[j].crd[1] - ccrd[1])/bvec
            crdz = ccrd[2] + bdld['CH'] * (mol.atoms[j].crd[2] - ccrd[2])/bvec

            crdx = round(crdx, 3)
            crdy = round(crdy, 3)
            crdz = round(crdz, 3)

            #Gaussian file
            gatms.append(gauatm(element, crdx, crdy, crdz))

            tiker = mol.atoms[j].gtype
            atid = mol.atoms[j].atid
            chainid = 'A'
            resid = mol.atoms[j].resid
            resname = mol.residues[resid].resname
            occp = 1.00
            tempfac = 0.00

            #PDB file
            atmi = pdbatm(tiker, atid, atname, resname, chainid, resid,
                          crdx, crdy, crdz, occp, tempfac)
            writepdbatm(atmi, smpdbf)

#------------------------------Sidechain------------------------------------
def build_small_model(mol, reslist, ionids, cutoff, smresids, smresace,
    smresnme, smresgly, smresant, smresact, smresknh, smreskco, smchg,
    smspin, addred, outf, sqmopt):

    """
    For building the small model
    1) Backbone O, H/HN was deleted.
    2) Backbone N, C, HA --> three Hs bonds to CA atom, with distance equals
       CH bond length.
    3) CA atom --> CH3 atom
    4) Other atoms is kept.

    Speical case: PRO does not have H atom, with CD atom instead, so there
    will be no influence since H was not considered in this modeling.
    """

    #Sidechain model file
    smpdbf = outf + '_small.pdb'

    #Gaussian
    goptf = outf + '_small_opt.com'
    gfcf = outf + '_small_fc.com'

    #GAMESS
    goptf2 = outf + '_small_opt.inp'
    gfcf2 = outf + '_small_fc.inp'

    #SQM
    siopf = outf + '_small_sqm.in'
    soopf = outf + '_small_sqm.out'

    #Delete the possible existing file
    del_files([smpdbf, gfcf, goptf])

    #-------------------------------------------------------------------------
    ###############################Sidechain model############################
    #-------------------------------------------------------------------------

    print("***Creating the small model...")

    gatms = [] #gaussian atom list

    for i in smresids:
        #1) For residue switching to ACE
        if i in smresace:
            write_ace(mol, i, gatms, smpdbf)
        #2) For residue switching to NME
        elif i in smresnme:
            write_nme(mol, i, gatms, smpdbf)
        #3) For residue switching to GLY
        elif i in smresgly:
            write_gly(mol, i, gatms, smpdbf)
        #4) For residue switching to CH3NH3+
        elif i in smresant:
            write_ant(mol, i, gatms, smpdbf)
        #5) For residue switching to CH3CO2-
        elif i in smresact:
            write_act(mol, i, gatms, smpdbf)
        #6) For residue which keep N and H in the model
        elif i in smresknh:
            write_sc_knh(mol, i, gatms, smpdbf)
        #7) For residue which keep C and O in the model
        elif i in smreskco:
            write_sc_kco(mol, i, gatms, smpdbf)
        #8) For normal amino acid residues, keep the small
        elif i in reslist.std:
            write_sc(mol, i, gatms, smpdbf)
        #9) For speical residue
        else:
            write_normal(mol, reslist, i, gatms, smpdbf)

    ln = count_lines(smpdbf)
    print("Totally there are " + str(ln) + " atoms in the small model.")

    if addred == 0:
        naddred_list = []
    else:
        # Get the metal site bond and angle list, and add them to the Gaussian
        # input file as additional reduncdant
        sm_mol, sm_atids, sm_resids = get_atominfo_fpdb(smpdbf)
        natids_dict = {}
        for i in range(0, len(sm_atids)):
            natids_dict[sm_atids[i]] = i + 1

        sm_blist = get_blist(sm_mol, sm_atids)
        bdatmids = [] #Just to prevent duplicated things happen
        addred_list = []
        naddred_list = []
        #---Get the atoms which is in the cutoff of metal ion
        for met in ionids:
            for i in sm_atids:
                if (i != met):
                    dis = calc_bond(mol.atoms[met].crd, mol.atoms[i].crd)
                    if (dis <= cutoff) and mol.atoms[i].element in BIND_ATOMS:
                        if i not in bdatmids:
                            bdatmids.append(i)
                            if i > met:
                                addred_list.append((met, i))
                            else:
                                addred_list.append((i, met))
 
                            for bond in sm_blist:
                                if (bond[0] == met) and (bond[1]!=i):
                                    if bond[1] > i:
                                        addred_list.append((i, met, bond[1]))
                                    else:
                                        addred_list.append((bond[1], met, i))
                                elif (bond[1] == met) and (bond[0]!=i):
                                    if bond[0] > i:
                                        addred_list.append((i, met, bond[0]))
                                    else:
                                        addred_list.append((bond[0], met, i))
                                elif (bond[0] == i) and (bond[1]!=met):
                                    if bond[1] > met:
                                        addred_list.append((met, i, bond[1]))
                                    else:
                                        addred_list.append((bond[1], i, met))
                                elif (bond[1] == i) and (bond[0]!=met):
                                    if bond[0] > met:
                                        addred_list.append((met, i, bond[0]))
                                    else:
                                        addred_list.append((bond[0], i, met))
 
        addred_list = sorted(list(set(addred_list)))
 
        for i in addred_list:
            if len(i) == 2:
                naddred_list.append((natids_dict[i[0]], natids_dict[i[1]]))
 
        for i in addred_list:
            if len(i) == 3:
                naddred_list.append((natids_dict[i[0]], natids_dict[i[1]], natids_dict[i[2]]))
        #print(naddred_list)

    #Calculate the spin number and print it into gaussian file
    gaelemts = 0
    for gatm in gatms:
        AtNum = Atnum[gatm.element]
        gaelemts = gaelemts + AtNum

    ElecNum = gaelemts - smchg
    ElecNum = int(round(ElecNum, 0))
    print("Totally there are " + str(ElecNum) + " electrons in the small model.")

    if smspin == -99:
        if ElecNum%2 == 0:
            SpinNum = 1
        else:
            SpinNum = 2
    else:
        SpinNum = smspin

    #Gaussian
    if addred == 0:
        write_gau_optf(outf, goptf, smchg, SpinNum, gatms, [])
        write_gau_fcf(outf, gfcf, [])
    if addred == 1:
        write_gau_optf(outf, goptf, smchg, SpinNum, gatms, naddred_list)
        write_gau_fcf(outf, gfcf, [])
    elif addred == 2:
        write_gau_optf(outf, goptf, smchg, SpinNum, gatms, [])
        write_gau_fcf(outf, gfcf, naddred_list)

    #GAMESS
    write_gms_optf(goptf2, smchg, SpinNum, gatms)
    write_gms_fcf(gfcf2, smchg, SpinNum)

    #Perform the SQM calcualtion under PM6 first
    if (sqmopt == 1) or (sqmopt == 3):
        #Delete the possible existing file
        del_files([siopf, soopf])
        write_sqm_optf(siopf, smchg, gatms)
        if SpinNum == 1:
            print("Performing SQM optimization of small model, please wait...")
            #Run SQM to optimize the coordinates
            os.system("sqm -i %s -o %s" %(siopf, soopf))
            gatms2 = get_crdinfo_from_sqm(soopf)
            write_gau_optf(outf, goptf, smchg, SpinNum, gatms2, naddred_list, 4)
            write_gms_optf(goptf2, smchg, SpinNum, gatms2, 4)
        else:
            print("Could not perform SQM optimization for the small model " + \
                  "with spin number not equal to 1.")

#------------------------------------Standard model---------------------------
def build_standard_model(mol, reslist, cutoff, msresids, outf, ionids,
                         bdedatms, addbpairs, libdict, autoattyp):

    #Standard model file
    stf = outf + '_standard.pdb'
    stpf = outf + '_standard.fingerprint'
    del_files([stf, stpf])

    #-------------------------------------------------------------------------
    ###############################Standard model#############################
    #-------------------------------------------------------------------------

    print("***Creating the standard model...")

    for i in msresids:
        print("It contains the residue " + str(i) + '-' + \
              mol.residues[i].resname + " as normal.")

        for j in mol.residues[i].resconter:

            atname = mol.atoms[j].atname
            if i in reslist.std and atname == 'HN':
                atname = 'H'

            crdx = mol.atoms[j].crd[0]
            crdy = mol.atoms[j].crd[1]
            crdz = mol.atoms[j].crd[2]
            tiker = mol.atoms[j].gtype
            atid = mol.atoms[j].atid
            chainid = 'A'
            resid = mol.atoms[j].resid
            resname = mol.residues[resid].resname
            occp = 1.00
            tempfac = 0.00

            crdx = round(crdx, 3)
            crdy = round(crdy, 3)
            crdz = round(crdz, 3)

            if i in reslist.nterm:
                attype = libdict['N'+resname + '-' + atname][0]
            elif i in reslist.cterm:
                attype = libdict['C'+resname + '-' + atname][0]
            else:
                attype = libdict[resname + '-' + atname][0]

            atmi = pdbatm(tiker, atid, atname, resname, chainid, resid,
                          crdx, crdy, crdz, occp, tempfac)

            writepdbatm(atmi, stf)
            stpff = open(stpf, 'a')

            #assign new atom types to atoms inside the metal site
            attype2 = attype

            if autoattyp == 1:
                if j in ionids:
                    ionresn = mol.atoms[j].resname
                    ionatn = mol.atoms[j].atname
                    ionelmt = mol.atoms[j].element
                    ionchg = libdict[ionresn + '-' + ionatn][1]
                    ionchg = int(ionchg)
                    if ionchg > 0:
                        attype2 = ionelmt + str(ionchg) + '+'
                    elif ionchg < 0:
                        attype2 = ionelmt + str(ionchg) + '-'
                    else:
                        raise pymsmtError('Ion Should not have 0 charge on it.')
            elif autoattyp == 2:
                #for metal ion
                for k in range(0, len(ionids)):
                    if (j == ionids[k]) and (k < 9):
                        k2 = k + 1
                        attype2 = 'M' + str(k2)
                    elif (j == ionids[k]) and (k >= 9):
                        k3 = k - 9 + 1
                        attype2 = 'X' + str(k3)
                #for bonded atoms
                for k in range(0, len(bdedatms)):
                    if (j == bdedatms[k]) and (k < 9):
                        k4 = k + 1
                        attype2 = 'Y' + str(k4)
                    elif (j == bdedatms[k]) and (k >= 9) and (k < 18):
                        k5 = k - 9 + 1
                        attype2 = 'Z' + str(k5)
                    elif (j == bdedatms[k]) and (k >= 18) and (k < 27):
                        k6 = k - 18 + 1
                        attype2 = 'A' + str(k6)
                    elif (j == bdedatms[k]) and (k >= 27):
                        k7 = k - 27 + 1
                        attype2 = 'B' + str(k7)

            print(str(resid) + '-' + resname + '-' + atname,
                  str(atid), attype, '->', attype2, file=stpff)
            stpff.close()

    #Print the link information into standard fingerprint file
    stpff = open(stpf, 'a')
    for met in ionids:
        for i in bdedatms:
            dis = calc_bond(mol.atoms[met].crd, mol.atoms[i].crd)
            if (dis <= cutoff):
                print("LINK", str(met)+'-'+mol.atoms[met].atname,
                      str(i)+'-'+mol.atoms[i].atname, file=stpff)
    #Add the residue numbers of additional bonded atoms
    for i in addbpairs:
        if (i[0] in ionids):
            print("LINK", str(i[0])+'-'+mol.atoms[i[0]].atname,
                      str(i[1])+'-'+mol.atoms[i[1]].atname, file=stpff)
        elif (i[1] in ionids):
            print("LINK", str(i[1])+'-'+mol.atoms[i[1]].atname,
                  str(i[0])+'-'+mol.atoms[i[0]].atname, file=stpff)
    stpff.close()

    ln = count_lines(stf)
    print("Totally there are " + str(ln) + " atoms in the standard model.")

#---------------------------------Large model---------------------------------
def build_large_model(mol, reslist, lmsresids, lmsresace, lmsresnme,
                      lmsresgly, ionids, chargedict, lgchg, lgspin, outf,
                      watermodel, largeopt, sqmopt):

    #Large model file
    lgpdbf = outf + '_large.pdb'
    lfpf = outf + '_large.fingerprint'
    gmkf = outf + '_large_mk.com'
    gmsf = outf + '_large_mk.inp'
    simkf = outf + '_large_sqm.in'
    somkf = outf + '_large_sqm.out'
    del_files([lgpdbf, lfpf, gmkf])

    #-------------------------------------------------------------------------
    ###############################Large model################################
    #-------------------------------------------------------------------------

    print("***Creating the large model...")
    gatms = []
    for i in lmsresids:
        #1) for atoms in ACE ---------------------------------------------------
        if i in lmsresace:
            write_ace(mol, i, gatms, lgpdbf, lfpf)
        #2) for atoms in NME ---------------------------------------------------
        elif i in lmsresnme:
            write_nme(mol, i, gatms, lgpdbf, lfpf)
        #3) for atoms in GLY ---------------------------------------------------
        elif i in lmsresgly:
            write_gly(mol, i, gatms, lgpdbf, lfpf)
        #4) for atoms in other residues ----------------------------------------
        else:
            write_normal(mol, reslist, i, gatms, lgpdbf, lfpf)

    ln = count_lines(lgpdbf)
    print("Totally there are " + str(ln) + " atoms in the large model.")

    #Calculate the spin number and print it into gaussian file
    gaelemts = 0
    for gatm in gatms:
        AtNum = Atnum[gatm.element]
        gaelemts = gaelemts + AtNum

    ElecNum = gaelemts - lgchg
    ElecNum = int(round(ElecNum, 0))
    print("Totally there are " + str(ElecNum) + " electrons in the large model.")

    if lgspin == -99:
        if ElecNum%2 == 0:
            SpinNum = 1
        else:
            SpinNum = 2
    else:
        SpinNum = lgspin

    #For Gaussian file
    IonLJParaDict = get_ionljparadict(watermodel)
    ionnames = [mol.atoms[i].atname for i in ionids]
    ionnames = list(set(ionnames))
    write_gau_mkf(outf, gmkf, lgchg, SpinNum, gatms, ionnames,
                  chargedict, IonLJParaDict, largeopt)

    #For GAMESS file
    write_gms_mkf(gmsf, lgchg, SpinNum, gatms)

    #-------------------------------------------------------------------------
    # Doing SQM Optimization
    #-------------------------------------------------------------------------
    if (sqmopt == 2) or (sqmopt == 3):
        del_files([simkf, somkf])
        write_sqm_optf(simkf, lgchg, gatms)
        if SpinNum == 1:
            print("Performing SQM optimization of large model, please wait...")
            os.system("sqm -i %s -o %s" %(simkf, somkf))
            gatms2 = get_crdinfo_from_sqm(somkf)
            write_gau_mkf(outf, gmkf, lgchg, SpinNum, gatms, ionnames,
                          chargedict, IonLJParaDict, largeopt, 4)
            write_gms_mkf(gmsf, lgchg, SpinNum, gatms2, 4)
        else:
            print("Could not perform SQM optimization for the large model " + \
                  "with spin number not equal to 1.")

def gene_model_files(pdbfile, ionids, addres, addbpairs, outf, ffchoice, naamol2f, cutoff, \
        watermodel, autoattyp, largeopt, sqmopt, smchg, smspin, addred, lgchg, lgspin):

    mol, atids, resids = get_atominfo_fpdb(pdbfile)

    reslist = get_reslist(mol, resids)

    libdict, chargedict = get_lib_dict(ffchoice)

    for mol2f in naamol2f:
        libdict1, chargedict1 = get_lib_dict(mol2f)
        libdict.update(libdict1)
        chargedict.update(chargedict1)

    #-------------------------------------------------------------------------
    # Get the residues in the metal site
    #-------------------------------------------------------------------------

    print("******************************************************************")
    print("*                                                                *")
    print("*=======================Metal Site Information===================*")
    print("*                                                                *")
    print("******************************************************************")

    #1. Metal ions information
    metresids = [] #metal ion residue id
    for i in ionids:
        resid = mol.atoms[i].resid
        metresids.append(resid)
        print("***Selected Metal ion " + mol.atoms[i].atname + " is atom " + \
              str(i) + " in residue " + str(mol.atoms[i].resid) + '-' + \
              mol.residues[resid].resname)
    ionids = list(set(ionids))
    ionids.sort()

    #2. Metal site residues information
    msresids = [] #metal site residues
    bdedatms, bdedatnams = get_ms_ids(mol, atids, ionids, cutoff, addbpairs)

    #3. Get the metal site containing residues
    for i in bdedatms:
        print(str(mol.atoms[i].resid) + '-' + \
              mol.atoms[i].resname + \
              '@' + mol.atoms[i].atname + ' is in ' + str(cutoff) + \
              ' Angstrom of or set bonded (in the input file) to (one of)'
              ' these metal ions')
                    #+ \
                    #str(mol.atoms[i].resid) + '-' + \
                    #mol.atoms[i].resname + '@' + \
                    #mol.atoms[i].atname
        if mol.atoms[i].resid not in msresids:
            msresids.append(mol.atoms[i].resid)

    msresids = msresids + metresids + addres
    msresids = list(set(msresids))
    msresids.sort()

    print("***The following residues are in the Metal Site:")
    totchg = 0.0
    for i in msresids:
        print("Residue " + str(i) + '-' + mol.residues[i].resname)
        if i in reslist.nterm:
            totchg = totchg + chargedict['N' + mol.residues[i].resname]
        elif i in reslist.cterm:
            totchg = totchg + chargedict['C' + mol.residues[i].resname]
        else:
            totchg = totchg + chargedict[mol.residues[i].resname]
    totchg = int(round(totchg, 0))

    if smchg == -99:
        smchg = totchg

    if lgchg == -99:
        lgchg = totchg

    #-------------------------------------------------------------------------
    # Get the residues for building the small model and print
    #-------------------------------------------------------------------------

    smresids = msresids
    smresace = []
    smresnme = []
    smresgly = []
    smresact = []
    smresant = []
    smresknh = [] #Residues to keep N and H, which is connect to the residue
                  #which has backbone oxygen bond to the ion and also bond to
                  #to the ion but with sidechain
    smreskco = [] #Residues to keep C and O, which is connect to the residue
                  #which has backbone nitrogen bond to the ion and also bond to
                  #to the ion but with sidechain

    bdedresids = []
    bdedresdict = {}

    for i in range(0, len(bdedatms)):
        atm = bdedatms[i]
        atname = bdedatnams[i]
        resid = mol.atoms[atm].resid
        if resid not in bdedresids:
            bdedresids.append(resid)

    for resid in bdedresids:
        resatns = []
        for i in bdedatms:
            if mol.atoms[i].resid == resid:
                resatns.append(mol.atoms[i].atname)
        bdedresdict[resid] = resatns

    for resid in bdedresids: #Here the bdedresids is equal to msresids
        resatns = bdedresdict[resid]

        #1. If residue is a n terminal residue
        if resid in reslist.nterm:
            if (set(['N3', 'O']) == set(resatns)) or (set(['N', 'O']) == set(resatns)):
                smresgly.append(resid)
                if resid+1 in smresids:
                    smresknh.append(resid+1)
                else:
                    smresnme.append(resid+1)
            elif ('N3' in resatns) or ('N' in resatns):
                if len(resatns) == 1: #If only nitrogen bond to
                    smresant.append(resid)
                else:
                    smresknh.append(resid)
            elif 'O' in resatns:
                if len(resatns) == 1: #If only oxygen bond to
                    smresace.append(resid)
                else:
                    smreskco.append(resid)
                if resid+1 in smresids:
                    smresknh.append(resid+1)
                else:
                    smresnme.append(resid+1)
        #2. If residue is a C terminal residue
        elif (resid in reslist.cterm):
            if (set(['N', 'O']) == set(resatns)) or (set(['N', 'OXT']) == set(resatns)):
                smresgly.append(resid)
                if resid-1 in smresids:
                    smreskco.append(resid-1)
                else:
                    smresace.append(resid-1)
            elif set(['O', 'OXT']) == set(resatns):
                smresact.append(resid)
            elif ('O' in resatns) or ('OXT' in resatns):
                if len(resatns) == 1: #If only oxygen bond to
                    smresact.append(resid)
                else:
                    smreskco.append(resid)
            elif ('N' in resatns):
                if len(resatns) == 1: #If only nitrogen bond to
                    smresnme.append(resid)
                else:
                    smresknh.append(resid)
                if (resid-1) in smresids:
                    smreskco.append(resid-1)
                else:
                    smresace.append(resid-1)

        #3. If residue is a standard residue but with backbone oxygen and/or nitrogen
        elif (resid in reslist.std):
            if set(['N', 'O']) < set(resatns):
                smresgly.append(resid)
                if (resid-1 in smresids) and (resid+1 in smresids):
                    smreskco.append(resid-1)
                    smresknh.append(resid+1)
                elif (resid-1 in smresids) and (resid+1 not in smresids):
                    smreskco.append(resid-1)
                    smresnme.append(resid+1)
                elif (resid+1 in smresids) and (resid-1 not in smresids):
                    smresknh.append(resid+1)
                    smresace.append(resid-1)
                else:
                    smresace.append(resid-1)
                    smresnme.append(resid+1)
            elif ('O' in resatns):
                if len(resatns) == 1: #If only O bond to, no sidechain atom bond to
                    smresace.append(resid)
                else:
                    smreskco.append(resid)
                if resid+1 in smresids:
                    smresknh.append(resid+1)
                else:
                    smresnme.append(resid+1)
            elif ('N' in resatns):
                if len(resatns) == 1: #If only N bond to, no sidechain atom bond to
                    smresnme.append(resid)
                else:
                    smresknh.append(resid)
                if (resid-1 in smresids):
                    smreskco.append(resid-1)
                else:
                    smresace.append(resid-1)

    smresids = smresids + smresace + smresnme + smresgly + smresknh + \
               smreskco + smresant + smresact
    smresids = list(set(smresids))
    smresids.sort()

    smresf = outf + '_small.res'

    w_smresf = open(smresf, 'w')
    for i in smresace:
        print('ACE-', i, file=w_smresf)
    for i in smresnme:
        print('NME-', i, file=w_smresf)
    for i in smresgly:
        print('GLY-', i, file=w_smresf)
    for i in smresknh:
        print('KNH-', i, file=w_smresf)
    for i in smreskco:
        print('KCO-', i, file=w_smresf)
    for i in smresant:
        print('ANT-', i, file=w_smresf)
    for i in smresact:
        print('ACT-', i, file=w_smresf)
    w_smresf.close()

    print("***The small model contains the following residues: ")
    print(smresids)

    #-------------------------------------------------------------------------
    # Get the residues for building the large model and print
    #-------------------------------------------------------------------------

    lmsresids = [] #large model metal site residues
    lmsresace = [] #large model metal site ACE
    lmsresnme = [] #large model metal site NME
    lmsresgly = [] #large model metal site GLY

    resamino = reslist.nterm + reslist.cterm + reslist.std
    resamino.sort()

    #First several residues
    for i in range(0, len(msresids)-1):
        resi = msresids[i]
        resj = msresids[i+1]
        if resi in resamino:
            if (set(range(resi, resj+1)) < set(resamino)) and (resj - resi <= 5):
            #If two residues within 5 residues apart, treat the residue as GLY
                for j in resids:
                    if (j > resi) and (j < resj) and (j not in lmsresgly):
                        lmsresgly.append(j) #GLY

    for i in range(0, len(msresids)):
        resid = msresids[i]
        resname = mol.residues[resid].resname
        if resid in reslist.nterm:
            if (resid+1 not in lmsresgly) and (resid+1 not in msresids):
                lmsresnme.append(resid+1) #NME
        elif resid in reslist.cterm:
            if (resid-1 not in lmsresgly) and (resid-1 not in msresids):
                lmsresace.append(resid-1) #ACE
        elif resid in reslist.std:
            if (resid-1 not in lmsresgly) and (resid-1 not in msresids):
                lmsresace.append(resid-1) #ACE
            if (resid+1 not in lmsresgly) and (resid+1 not in msresids):
                lmsresnme.append(resid+1) #NME

    lmsresids = msresids + lmsresace + lmsresnme + lmsresgly #Combine the residues
    lmsresids = list(set(lmsresids)) #Delete repeat elements
    lmsresids.sort() #Sort the list

    print("***The large model contains the following residues: ")
    print(lmsresids)

    #-------------------------------------------------------------------------
    # Generate model files
    #-------------------------------------------------------------------------

    print("******************************************************************")
    print("*                                                                *")
    print("*=======================Building models==========================*")
    print("*                                                                *")
    print("******************************************************************")

    build_small_model(mol, reslist, ionids, cutoff, smresids, smresace, smresnme,
                    smresgly, smresant, smresact, smresknh, smreskco, smchg,
                    smspin, addred, outf, sqmopt)

    build_standard_model(mol, reslist, cutoff, msresids, outf, ionids,
                         bdedatms, addbpairs, libdict, autoattyp)

    build_large_model(mol, reslist, lmsresids, lmsresace, lmsresnme, lmsresgly,
              ionids, chargedict, lgchg, lgspin, outf, watermodel, largeopt, sqmopt)

    #Using the automatically detect bond method for the backup
    #else:
    #  for met in ionids:
    #    Radiusmet = CoRadiiDict[mol.atoms[met].element]
    #    for i in atids:
    #      Radiusi = CoRadiiDict[mol.atoms[i].element]
    #      if (i != met):
    #       dis = calc_bond(mol.atoms[met].crd, mol.atoms[i].crd)
    #        cutoff = Radiusmet + Radiusi + 0.40
    #        if (dis <= cutoff) and (dis > 0.1):
    #          print(str(mol.atoms[i].resid) + '-' + \
    #                mol.atoms[i].resname + \
    #                '@' + mol.atoms[i].atname + ' is in ' + str(cutoff) + \
    #                ' Angstrom of ' + str(mol.atoms[met].resid) + '-' + \
    #                mol.atoms[met].resname + '@' + \
    #                mol.atoms[met].atname)
    #          if (mol.atoms[i].resid not in msresids):
    #            msresids.append(mol.atoms[i].resid)
