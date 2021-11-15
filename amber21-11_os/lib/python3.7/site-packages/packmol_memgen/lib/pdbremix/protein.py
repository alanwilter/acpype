# encoding: utf-8

__doc__ = """

Common protein manipulations.

This is a collection of some common PDB and Soup manipulations
of protein strucutres.
"""


import string

from . import data
from . import pdbatoms
from . import util
from . import pdbtext
from . import v3
from . import spacehash


def strip_solvent_pdb(pdb):
  new_pdb = util.fname_variant(pdb)
  txt = pdbtext.strip_solvent(open(pdb).read())
  open(new_pdb, 'w').write(txt)
  return new_pdb


def find_ca_of_resname(atoms, resname):
  for atom in atoms:
    if pdbatoms.split_tag(resname) == \
       (atom.chain_id, atom.res_num, atom.res_insert):
      if "CA" == atom.type:
        return atom
  raise IndexError("Can't find atom %s" % resname)


def get_pdb_transform(pdb, center_res, top_res):
  """
  Returns a transformation matrix that centers pdb to 
  center_res on the z-axis and moves top_res above center_res
  on the y-axis
  """
  soup = pdbatoms.Soup(pdb)
  atoms = soup.atoms()
  soup_center = pdbatoms.get_center(atoms)
  translation = v3.translation(-soup_center)
  soup.transform(translation)
  result = translation

  view = v3.vector(0, 0, 1)

  if center_res is not None:
    center_atom = find_ca_of_resname(soup.atoms(), center_res)
    axis = v3.cross(view, center_atom.pos)
    angle = v3.vec_dihedral(view, axis, center_atom.pos)
    rotation1 = v3.rotation(axis, angle)
    soup.transform(rotation1)
    result = v3.combine(rotation1, result)

  if top_res is not None:
    top_atom = find_ca_of_resname(soup.atoms(), top_res)
    top_dir = v3.vector(0, 1, 0)
    axis = v3.vector(view)
    angle = v3.vec_dihedral(top_dir, axis, top_atom.pos)
    rotation2 = v3.rotation(axis, angle)
    result = v3.combine(rotation2, result)

  return result


def transform_pdbs_to_residues_of_first_pdb(pdbs, center_res, top_res):
  transform = get_pdb_transform(pdbs[0], center_res, top_res)
  new_pdbs = []
  for pdb in pdbs:
    new_pdb = util.fname_variant(pdb)
    soup = pdbatoms.Soup(pdb)
    soup.transform(transform)
    soup.write_pdb(new_pdb)
    new_pdbs.append(new_pdb)
  return new_pdbs


def transformed_soup_from_pdb(
    pdb, center_res=None, top_res=None, 
    width=None, height=None, frame_residues=None):
  soup = pdbatoms.Soup(pdb)
  if center_res and top_res:
    transform = get_pdb_transform(pdb, center_res, top_res)
    soup.transform(transform)
  if frame_residues:
    resnames = [pymol_id_from_res_tag(r) for r in frame_residues]
    soup.frame_pymol_script = "zoom (%s)\n" % ' or '.join(resnames)
  if width: soup.width = width
  if height: soup.height = height
  return soup


def is_peptide_connected(res1, res2):
  if res1.has_atom('CA') and \
     res1.has_atom('C') and \
     res2.has_atom('N') and \
     res2.has_atom('CA'):
     d = v3.distance(res1.atom('C').pos, res2.atom('N').pos)
     if d < 2.0:
       return True
  return False


def find_chains(soup):
  residues = soup.residues()
  n = len(residues)
  if n == 0:
    return
  i_chain = 0
  chain_id = string.ascii_uppercase[i_chain]
  for i in range(0, n):
    res = residues[i]
    if i == 0:
      is_connected_to_prev = False
    else:
      prev_res = residues[i-1]
      is_connected_to_prev = is_peptide_connected(prev_res, res)
    if i < n-1:
      next_res = residues[i+1]
      is_connected_to_next = is_peptide_connected(res, next_res)
    else:
      is_connected_to_next = False
    if is_connected_to_prev or is_connected_to_next:
      res.set_chain_id(chain_id)
      if not is_connected_to_next:
        i_chain += 1
        chain_id = string.ascii_uppercase[i_chain]


def is_connected(i, j, soup, cutoff=3.5):
  if i == j:
    return False
  min_dist = 1000.0
  for atom_i in soup.residue(i).atoms():
    for atom_j in soup.residue(j).atoms():
      dist = v3.distance(atom_i.pos, atom_j.pos)
      if dist < min_dist:
        min_dist = dist
  return min_dist < cutoff


backbone = ['CA', 'HA', 'N', 'H', 'O', 'C']
def is_sidechain_connected(i, j, soup, cutoff=3.5):
  if abs(i-j) <= 2:
    return False
  min_dist = 1000.0
  sidechain_atoms_i = [a for a in soup.residue(i).atoms() 
                       if a.type not in backbone]
  for atom_i in sidechain_atoms_i:
    for atom_j in soup.residue(j).atoms():
      dist = v3.distance(atom_i.pos, atom_j.pos)
      if dist < min_dist:
        min_dist = dist
  return min_dist < cutoff


def find_bb_hbonds(residues):

    cutoff_d_of_n_o = 3.5

    vertices = []
    atoms = []
    for i_residue, residue in enumerate(residues):
        residue.i = i_residue
        for atom in residue.atoms():
            atom.residue = residue
        if residue.has_atom('O'):
            atom = residue.atom('O')
            atoms.append(atom)
            vertices.append(atom.pos)
        if residue.has_atom('N'):
            atom = residue.atom('N')
            atoms.append(atom)
            vertices.append(atom.pos)
    
    for i, j in spacehash.SpaceHash(vertices).close_pairs():
        if abs(i - j) < 3:
            continue

        if atoms[i].element == 'O' and atoms[j].element == 'N':
            o = atoms[i]
            n = atoms[j]
        elif atoms[i].element == 'N' and atoms[j].element == 'O':
            n = atoms[i]
            o = atoms[j]
        else:
            continue

        if v3.distance(o.pos, n.pos) < cutoff_d_of_n_o:
            o.residue.co_partners.append(n.residue.i)
            n.residue.nh_partners.append(o.residue.i)


def unique_append(a_list, item):
    if item not in a_list:
        a_list.append(item)
        a_list.sort()


def find_ss_by_bb_hbonds(soup):
    """
    Analyzes a protein soup and adds the fields to each residue
    to indicate secondary-structure:

      -  res.co_partners = []
      -  res.nh_partners = []
      -  res.beta_contacts = []
      -  res.alpha_contacts = []
      -  res.ss = "C"
      -  res.ss_contacts = []

    The key field is ss_contacts which lists the indices of the
    residues that are naturally in contact due to secondary-structure
    geometry.
    """

    residues = soup.residues()
    n_res = len(residues)

    for res in residues:
        res.co_partners = []
        res.nh_partners = []
        res.beta_contacts = []
        res.alpha_contacts = []
        res.ss_contacts = []
        res.ss = "C"

    find_bb_hbonds(residues)

    def is_conh(i_res, j_res):
        if not (0 <= i_res < n_res):
            return False
        if not (0 <= j_res < n_res):
            return False
        return j_res in residues[i_res].co_partners

    def make_alpha_contacts(i_res, j_res):
        unique_append(residues[i_res].alpha_contacts, j_res)
        unique_append(residues[j_res].alpha_contacts, i_res)

    for i_res in range(n_res):
        if is_conh(i_res - 1, i_res + 3) and is_conh(i_res, i_res + 4):
            # alpha-helix
            for j_res in range(i_res, i_res + 4):
                residues[j_res].ss = 'H'
            make_alpha_contacts(i_res + 3, i_res)
            make_alpha_contacts(i_res + 4, i_res)

    def make_beta_contact(i_res, j_res):
        unique_append(residues[i_res].beta_contacts, j_res)
        unique_append(residues[i_res].beta_contacts, j_res-1)
        unique_append(residues[i_res].beta_contacts, j_res+1)

    def make_beta_contacts(i_res, j_res):
        make_beta_contact(i_res, j_res)
        make_beta_contact(j_res, i_res)

    for i_res in range(n_res):
        for j_res in range(n_res):
            is_beta = False
            if abs(i_res - j_res) <= 2:
                # can't have beta contacts 2 or less apart
                is_beta = False
            elif is_conh(i_res, j_res) and is_conh(j_res, i_res):
                # anti-parallel beta-sheet h-bonded pair
                is_beta = True
            elif is_conh(i_res - 1, j_res + 1) and is_conh(i_res + 1, j_res - 1):
                # anti-parallel beta-sheet non-h-bonded pair
                is_beta = True
            elif is_conh(i_res, j_res - 1) and is_conh(j_res - 1, i_res):
                # parallel beta sheet pairs
                is_beta = True
            if is_beta:
                make_beta_contacts(i_res, j_res)
                residues[i_res].ss = "E"
                residues[j_res].ss = "E"
    
    def check_piece_is_e(i_res, j_res):
        for k_res in range(i_res, j_res+1):
            if k_res < 0 or k_res >= n_res:
                return False
            if residues[k_res].ss != "E":
                return False
        return True

    # add cross-strand i->j-2 contacts as beta
    # TODO: distinguish hb and non-hb pairs
    for i_res in range(n_res):
        if residues[i_res].ss == "E":
            for j_res in residues[i_res].beta_contacts:
                if check_piece_is_e(j_res - 2, j_res):
                    unique_append(residues[i_res].beta_contacts, j_res - 2)
                    unique_append(residues[j_res - 2].beta_contacts, i_res)
                if check_piece_is_e(j_res, j_res + 2):
                    unique_append(residues[i_res].beta_contacts, j_res + 2)
                    unique_append(residues[j_res + 2].beta_contacts, i_res)

    for i_res in range(n_res):
      residue = residues[i_res]
      residue.ss_contacts.extend(residue.alpha_contacts)
      residue.ss_contacts.extend(residue.beta_contacts)




