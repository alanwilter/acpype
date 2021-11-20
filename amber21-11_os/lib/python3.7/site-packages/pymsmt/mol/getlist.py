"""
This module use the bond linkage information to generate all the linkage
(angle, dihedral and nonbonded) information.
"""

from pymsmt.mol.mol import Linklist
from pymsmt.mol.cal import calc_bond
from pymsmt.mol.element import CoRadiiDict

#sort the barray diction to a list
#barrayc = sorted(barray.iteritems(), key=lambda d:d[0])

"""
Use to generate all the linkage information:
Usuage:
  get_all_list(Molecule, Bondlist)
"""

def get_all_list(mol, blist, atids, cutoff):

    ###1. Bond list
    blist = sorted(blist)

    ###2. Angle list
    alist = []

    for i in range(0, len(blist)):
        ati1 = blist[i][0]
        ati2 = blist[i][1]
        for j in range(i+1, len(blist)):
            atj1 = blist[j][0]
            atj2 = blist[j][1]
            if (ati1 == atj1):
                at1 = atj2
                at2 = ati1
                at3 = ati2
                angats = (at1, at2, at3)
                alist.append(angats)
            elif (ati1 == atj2):
                at1 = atj1
                at2 = ati1
                at3 = ati2
                angats = (at1, at2, at3)
                alist.append(angats)
            elif (ati2 == atj1):
                at1 = ati1
                at2 = ati2
                at3 = atj2
                angats = (at1, at2, at3)
                alist.append(angats)
            elif (ati2 == atj2):
                at1 = ati1
                at2 = ati2
                at3 = atj1
                angats = (at1, at2, at3)
                alist.append(angats)

    alist = get_pure_type(alist)

    ###3. Dihedral list
    dlist = []

    for i in range(0, len(alist)):
        ati1 = alist[i][0]
        ati2 = alist[i][1]
        ati3 = alist[i][2]
        for j in range(i + 1, len(alist)):
            atj1 = alist[j][0]
            atj2 = alist[j][1]
            atj3 = alist[j][2]
            if (ati2 == atj1) & (ati3 == atj2) & (ati1 != atj3):
                at1 = ati1
                at2 = ati2
                at3 = ati3
                at4 = atj3
                dihats = (at1, at2, at3, at4)
                dlist.append(dihats)
            elif (ati2 == atj3) & (ati3 == atj2) & (ati1 != atj1):
                at1 = ati1
                at2 = ati2
                at3 = ati3
                at4 = atj1
                dihats = (at1, at2, at3, at4)
                dlist.append(dihats)
            elif (ati1 == atj2) & (ati2 == atj3) & (atj1 != ati3):
                at1 = atj1
                at2 = ati1
                at3 = ati2
                at4 = ati3
                dihats = (at1, at2, at3, at4)
                dlist.append(dihats)
            elif (ati1 == atj2) & (ati2 == atj1) & (atj3 != ati3):
                at1 = atj3
                at2 = ati1
                at3 = ati2
                at4 = ati3
                dihats = (at1, at2, at3, at4)
                dlist.append(dihats)

    dlist = get_pure_type(dlist)


    ###4. Improper torsion list
    ilist = [] #Second is the centeral atom
    for i in range(0, len(alist)):
        ati1 = alist[i][0]
        ati2 = alist[i][1]
        ati3 = alist[i][2]
        for j in range(0, len(blist)):
            atj1 = blist[j][0]
            atj2 = blist[j][1]
            if (ati2 == atj1) and (atj2 != ati1) and (atj2 != ati3):
                at1 = ati1
                at2 = ati3
                at3 = ati2
                at4 = atj2
                impats = (at1, at2, at3, at4)
                ilist.append(impats)
            elif (ati2 == atj2) and (atj1 != ati1) and (atj1 != ati3):
                at1 = ati1
                at2 = ati3
                at3 = ati2
                at4 = atj1
                impats = (at1, at2, at3, at4)
                ilist.append(impats)

    ilist2 = [] #Third is the centeral atoms
    for imp in ilist:
        imp1 = (imp[0], imp[1], imp[2], imp[3])
        imp2 = (imp[0], imp[3], imp[2], imp[1])
        imp3 = (imp[1], imp[0], imp[2], imp[3])
        imp4 = (imp[1], imp[3], imp[2], imp[0])
        imp5 = (imp[3], imp[0], imp[2], imp[1])
        imp6 = (imp[3], imp[1], imp[2], imp[0])
        if (imp1 not in ilist2) and (imp2 not in ilist2) and (imp3 not in ilist2) \
          and (imp4 not in ilist2) and (imp5 not in ilist2) and (imp6 not in ilist2):
            ilist2.append(imp)

    ###5. nonbonded array
    ##get bonded atom list
    bondedatomlist = []

    #bond
    for i in range(0, len(blist)):
        atm1 = blist[i][0]
        atm2 = blist[i][1]
        if (atm1 < atm2):
            bondedatomlist.append((atm1, atm2))
        else:
            bondedatomlist.append((atm2, atm1))

    #angle
    for i in range(0, len(alist)):
        atm1 = alist[i][0]
        atm2 = alist[i][-1]
        if (atm1 < atm2):
            bondedatomlist.append((atm1, atm2))
        else:
            bondedatomlist.append((atm2, atm1))

    #dihedral
    for i in range(0, len(dlist)):
        atm1 = dlist[i][0]
        atm2 = dlist[i][-1]
        if (atm1 < atm2):
            bondedatomlist.append((atm1, atm2))
        else:
            bondedatomlist.append((atm2, atm1))

    bondedatomlist = set(bondedatomlist)

    ##get total atom list
    totlist = []
    for i in range(0, len(atids)):
        for j in range(i+1, len(atids)):
            atm1 = atids[i]
            atm2 = atids[j]
            if (atm1 < atm2):
                totlist.append((atm1, atm2))
            else:
                totlist.append((atm2, atm1))

    totlist = set(totlist)

    ##Get total nb list
    nblist = totlist - bondedatomlist
    nblist = sorted(list(nblist))

    fnblist = []

    for i in range(0, len(nblist)):
        atm1 = nblist[i][0]
        atm2 = nblist[i][1]

        crd1 = mol.atoms[atm1].crd
        crd2 = mol.atoms[atm2].crd

        dis = calc_bond(crd1, crd2)

        if (dis <= cutoff):
            fnblist.append(nblist[i])

    del nblist
    del totlist
    del bondedatomlist

    all_list = Linklist(blist, alist, dlist, ilist2, fnblist)

    return all_list

def get_alist(mol, blist):
    blist = sorted(blist)
    alist = []

    for i in range(0, len(blist)):
        ati1 = blist[i][0]
        ati2 = blist[i][1]
        for j in range(i+1, len(blist)):
            atj1 = blist[j][0]
            atj2 = blist[j][1]
            if (ati1 == atj1):
                at1 = atj2
                at2 = ati1
                at3 = ati2
                angats = (at1, at2, at3)
                alist.append(angats)
            elif (ati1 == atj2):
                at1 = atj1
                at2 = ati1
                at3 = ati2
                angats = (at1, at2, at3)
                alist.append(angats)
            elif (ati2 == atj1):
                at1 = ati1
                at2 = ati2
                at3 = atj2
                angats = (at1, at2, at3)
                alist.append(angats)
            elif (ati2 == atj2):
                at1 = ati1
                at2 = ati2
                at3 = atj1
                angats = (at1, at2, at3)
                alist.append(angats)

    alist = get_pure_type(alist)
    return alist

def get_blist(mol, atids):
    blist = []
    for i in range(0, len(atids)):
        crdi = mol.atoms[atids[i]].crd
        ati = mol.atoms[atids[i]].element
        if (len(ati) == 2):
            ati = ati[0] + ati[1].lower()
        radiusi = CoRadiiDict[ati]
        for j in range(i+1, len(atids)):
            crdj = mol.atoms[atids[j]].crd
            atj = mol.atoms[atids[j]].element
            if (len(atj) == 2):
                atj = atj[0] + atj[1].lower()
            radiusj = CoRadiiDict[atj]
            radiusij = radiusi + radiusj + 0.40
            dis = calc_bond(crdi, crdj)
            if (dis > 0.1) and (dis <= radiusij):
                blist.append((atids[i], atids[j], 1))
    return blist

def get_mc_blist(mol, atids, ionids, fpf):
    "Get Metal Site Bond List"

    blist = []

    rlnk = open(fpf, 'r')
    for line in rlnk:
        if line[0:4] == "LINK":
            line = line.strip('\n')
            line = line.split()
            ati = line[1].split('-')
            atidi = int(ati[0])
            atj = line[2].split('-')
            atidj = int(atj[0])
            if atidi < atidj:
                blist.append((atidi, atidj, 1))
            else:
                blist.append((atidj, atidi, 1))
    rlnk.close()

    for i in range(0, len(atids)):
        crdi = mol.atoms[atids[i]].crd
        ati = mol.atoms[atids[i]].element
        if (len(ati) == 2):
            ati = ati[0] + ati[1].lower()
        radiusi = CoRadiiDict[ati]
        for j in range(i+1, len(atids)):
            crdj = mol.atoms[atids[j]].crd
            atj = mol.atoms[atids[j]].element
            if (len(atj) == 2):
                atj = atj[0] + atj[1].lower()
            radiusj = CoRadiiDict[atj]
            radiusij = radiusi + radiusj + 0.40
            dis = calc_bond(crdi, crdj)

            if list(set([atids[i], atids[j]]) & set(ionids)) == []:
            #If there is no metal ion in the bond
                if (dis > 0.1) and (dis <= radiusij):
                    blist.append((atids[i], atids[j], 1))

    blist = sorted(blist)

    return blist

def get_pure_type(onelist):
    newlist = []
    for i in onelist:
        if (not i in newlist) & (not i[::-1] in newlist):
            newlist.append(i)
    return newlist

def get_pure_num(onelist):
    newlist = []
    for i in onelist:
        if not i in newlist:
            newlist.append(i)
    return newlist

