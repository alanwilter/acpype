
from pymsmt.mol.element import resnamel

"""
This is the code for defining the molecule class
  Molecule class that contains biomolecular data:
  * 0. atids : atom ID list, e.g. [1, 2, 3, ...], int
  * 1. atnames : sequential list of atoms in an array, e.g. [H1, C1, O1,...], char
  * 2. elements : sequential list of elements, e.g. [H, C, O, ...]], char
  * 3. atomtypes : sequential list of atom types, e.g. [H, CT, O2,...], char
  # 4 - 6 is combined together
  * 4. resnames : residue names list, e.g. [CYS, HIS, ...], char
  * 5. resids : residue ID of each residue, e.g. [1, 2, 3, ...], int
  * 6. resconters : the atom ID in each residue, e.g.[(1, 2, 3),(4, 5, 6), ...], int
  # the data information about each molecule
  * 7. crds : coordinates, e.g. [(x1, y1, z1), (x2, y2, z2),...], float
  * 8. charges : sequential list of atom charge in an array, e.g. [-0.382, 0.532, ...], float
  """

class Atom:
    def __init__(self, gtype, atid, atname, element, atomtype, crd, charge, resid, resname):
        self.gtype = gtype
        self.atid = atid
        self.atname = atname
        self.element = element
        self.atomtype = atomtype
        self.crd = crd
        self.charge = charge
        self.resid = resid
        self.resname = resname

class Residue:
    def __init__(self, resid, resname, resconter):
        self.resid = resid
        self.resname = resname
        self.resconter = resconter

class Molecule:

    def __init__(self, atoms, residues):
        self.atoms = atoms
        self.residues = residues

    def renum(self):
        #atom id and resid dict
        atomdict = {}
        resdict = {}

        #dict to store the data
        Atoms = {}
        Residues = {}

        natids = list(self.atoms.keys())
        natids.sort()
        nresids = list(self.residues.keys())
        nresids.sort()

        #renumber the atoms
        for i in range(1, len(self.atoms)+1):
            Atoms[i] = self.atoms[natids[i-1]]

        #renumber the residues
        for i in range(1, len(self.residues)+1):
            Residues[i] = self.residues[nresids[i-1]]

        #return the molecule
        mol1 = Molecule(Atoms, Residues)
        return mol1

    def delwaterion(self):
        #for each residue
        for i in range(1, len(self.residues)+1):
            if (self.residues[i].resname in ['WAT', 'HOH']) or \
              (self.residues[i].resname[-1] in ['+', '-']):
                #for each atom in the residue
                for j in self.residues[i].resconter:
                    del self.atoms[j]
                del self.residues[i]
        #atids = self.atoms.keys().sort()
        #resids = self.residues.keys().sort()

    def delwater(self):
        #for each residue
        for i in range(1, len(self.residues)+1):
            if self.residues[i].resname in ['WAT', 'HOH']:
                #for each atom in the residue
                for j in self.residues[i].resconter:
                    del self.atoms[j]
                del self.residues[i]
        #atids = self.atoms.keys().sort()
        #resids = self.residues.keys().sort()
        #return self

    def delion(self):
        #for each residue
        for i in range(1, len(self.residues)+1):
            if self.residues[i].resname[-1] in ['+', '-']:
                #for each atom in the residue
                for j in self.residues[i].resconter:
                    del self.atoms[j]
                del self.residues[i]
        #atids = self.atoms.keys().sort()
        #resids = self.residues.keys().sort()
        #return self, atids, resids

    def keepaas(self):
        for i in range(1, len(self.residues)+1):
            if (self.residues[i].resname not in resnamel):
                #for each atom in the residue
                for j in self.residues[i].resconter:
                    del self.atoms[j]
                del self.residues[i]

class Linklist:
    def __init__(self, bondlist, anglist, dihlist, implist, nblist):
        self.bondlist = bondlist
        self.anglist = anglist
        self.dihlist = dihlist
        self.implist = implist
        self.nblist = nblist

class pdbatm:
    def __init__(self, tiker, atid, atname, resname, chainid, resid, crdx, crdy, crdz, occp, tempfac):
        self.tiker = tiker
        self.atid = atid
        self.atname = atname
        self.resname = resname
        self.chainid = chainid
        self.resid = resid
        self.crdx = crdx
        self.crdy = crdy
        self.crdz = crdz
        self.occp = occp
        self.tempfac = tempfac

class gauatm:
    def __init__(self, element, crdx, crdy, crdz):
        self.element = element
        self.crdx = crdx
        self.crdy = crdy
        self.crdz = crdz

class XYZatom:
    """Class for the xyz file containing atoms ::: XYZatom(element, coordinates)"""

    def __init__(self, element, crd):
        self.element = element
        self.crd = crd

class residuelist:
    def __init__(self, cterm, nterm, std, nonstd, water):
        self.cterm = cterm
        self.nterm = nterm
        self.std = std
        self.nonstd = nonstd
        self.water = water

def get_reslist(mol, resids):
    cterm = []
    nterm = []
    std = []
    nonstd = []
    water = []

    for i in resids:
        resnamei = mol.residues[i].resname
        atnames = []

        for j in mol.residues[i].resconter:
            atnamej = mol.atoms[j].atname
            atnames.append(atnamej)

        if (set(['CA', 'N', 'C', 'O', 'OXT', 'H2', 'H3']) < set(atnames)) or \
           (set(['CA', 'N', 'C', 'O', 'OXT', 'HN2', 'HN3']) < set(atnames)):
        # If there is a isolated residue
            nonstd.append(i)
        elif (set(['CA', 'N', 'C', 'O', 'H2', 'H3']) < set(atnames)) or \
             (set(['CA', 'N', 'C', 'O', 'HN2', 'HN3']) < set(atnames)):
        # If it is a N-terminal residue
            nterm.append(i)
        elif set(['CA', 'N', 'C', 'O', 'OXT']) < set(atnames):
        # If it is a C-terminal residue
            cterm.append(i)
        elif set(['CA', 'N', 'C', 'O']) < set(atnames):
            if i == min(resids):
                nterm.append(i)
            elif i == max(resids):
                cterm.append(i)
            else:
                std.append(i)
        else:
            nonstd.append(i)

        if resnamei in ['WAT', 'HOH']:
            water.append(i)

    reslist = residuelist(cterm, nterm, std, nonstd, water)
    return reslist
