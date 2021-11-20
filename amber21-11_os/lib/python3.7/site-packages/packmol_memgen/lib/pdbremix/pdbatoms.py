# encoding: utf-8

__doc__ = """

Provides the Soup object to manipulate protein structures.

The Soup object contains a list of Atom objects, which are also
grouped into a list of Residues. The Residues provide a
convenient way to search and access Atoms.

Specifically, chain_ids are not used to organize the data
structures. In the author's experience, for the amount of work to
maintain chain structure, not much utility is gained. As well,
chain_id has only loose semantics that are not strictly
hierchical to residues. By sticking to a Soup as a group of 
residues, the resultant data structure is much cleaner. Chain
analysis can easily be done on a case-by-case basis.
"""


from . import v3
import copy
import string

from . import data


class Atom:
  """
  This is the basic object to hold Atom information.

  The attributes are basically those of a PDB atom field.
  However, pos and vel are proper vectors that can be manipulated
  with the v3 3d-vector geometry library.  
  """

  def __init__(self, pos=None, atom_type="", res_num=None):
    """
    Normally initialized as an empty container, and filled
    up progressively as fields are read by parsers.
    """
    self.is_hetatm = False
    self.pos = v3.vector() if pos is None else pos
    self.vel = v3.vector()
    self.mass = 0.0
    self.charge = 0.0
    self.type = ""
    self.element = ""
    self.chain_id = " "
    self.res_type = ""
    self.res_num = ""
    self.res_insert = ""
    self.bfactor = 0.0
    self.occupancy = 0.0
    self.num = 0
    self.alt_conform = " "

  def copy(self):
    return copy.deepcopy(self)

  def type_str(self):
    """
    Format atom_type to write to a PDB file's atom line.
    """
    atom_type = self.type.strip()
    if len(atom_type) == 1:
      atom_type = " %s  " % atom_type
    elif len(atom_type) == 2:
      if atom_type[0].isdigit():
        atom_type = "%s  " % atom_type
      else:
        atom_type = " %s " % atom_type
    elif len(atom_type) == 3:
      if atom_type[0].isdigit():
        atom_type = "%s " % atom_type
      else:
        atom_type = " %s" % atom_type
    return atom_type    

  def pdb_str(self):
    """
    Returns a string for output to an PDB file.
    """
    if self.is_hetatm:
      field = "HETATM"
    else:
      field = "ATOM  "
    x, y, z = self.pos
    s = "%6s%5s %4s %-4s%1s%4s%1s   %8.3f%8.3f%8.3f%6.2f%6.2f" % \
        (field, 
         str(self.num)[-5:], 
         self.type_str(),
         self.res_type, 
         self.chain_id,
         str(self.res_num)[-4:], 
         self.res_insert,
         x, y, z,
         self.occupancy, 
         self.bfactor)
    return s
               
  def res_tag(self):
    tag = ""
    if self.chain_id != " " and self.chain_id != "":
      tag += self.chain_id + ":"
    tag += str(self.res_num)
    if self.res_insert:
      tag += self.res_insert
    return tag  

  def __str__(self):
    x, y, z = self.pos
    return "%s:%s-%s" %  \
        (self.res_tag(), self.res_type, self.type)

  def transform(self, matrix):
    """
    Transforms the pos vector by a v3.transform matrix.
    """
    new_pos = v3.transform(matrix, self.pos)
    v3.set_vector(self.pos, new_pos)


def AtomFromPdbLine(line):
  """
  Returns an Atom object from an atom line in a pdb file.
  """
  atom = Atom()
  if line.startswith('HETATM'):
    atom.is_hetatm = True
  else:
    atom.is_hetatm = False
  atom.num = int(line[6:11])
  atom.type = line[12:16].strip(" ")
  atom.alt_conform = line[16]
  atom.res_type = line[17:21].strip()
  atom.element = data.guess_element(atom.res_type, atom.type)
  atom.chain_id = line[21]
  atom.res_num = int(line[22:26])
  atom.res_insert = line[26]
  if atom.res_insert == " ":
    atom.res_insert = ""
  x = float(line[30:38])
  y = float(line[38:46])
  z = float(line[46:54])
  v3.set_vector(atom.pos, x, y, z)
  try:
    atom.occupancy = float(line[54:60])
  except:
    atom.occupancy = 100.0
  try:
    atom.bfactor = float(line[60:66])
  except:
    atom.bfactor = 0.0
  return atom
  

# The following functions is for handling lists of atoms
  
def cmp_atom(a1):
  """
  Sorting operator for atoms
  """
  return a1.num


def add_radii(atoms):
  """
  Lookup and assign atom.radius for atoms.
  """
  for atom in atoms:
    if atom.element in data.radii:
      atom.radius = data.radii[atom.element]
    else:
      atom.radius = data.radii['.']


def get_center(atoms):
  """
  Returns the geometric center position vector of atoms.
  """
  center = v3.vector()
  for atom in atoms:
    center += atom.pos
  result = v3.scale(center, 1.0/float(len(atoms)))
  return result


def get_width(atoms, center=None):
  """
  Returns twice the longest distance from the center.
  """
  max_diff = 0
  if center is None:
    center = get_center(atoms)
  for atom in atoms:
    diff = v3.distance(atom.pos, center)
    if diff > max_diff:
      max_diff = diff
  return 2*max_diff


def read_pdb(fname):
  """
  Reads a list of Atoms from a PDB file.
  """
  atoms = []
  for line in open(fname, 'r'):
    if line.startswith(("ENDMDL", "END")):
      break
    if line.startswith(("ATOM", "HETATM")):
      atoms.append(AtomFromPdbLine(line))
  return atoms


def write_pdb(atoms, pdb):
  """
  Writes a list of atoms to a PDB file.
  """
  with open(pdb, 'w') as f:
    for atom in sorted(atoms, key=cmp_atom):
      f.write(atom.pdb_str() + '\n')


# Introducing the Residue structure for organizing atoms

def split_tag(tag):
  """
  Returns (chain_id, res_num, insert). Empty chain_id=" " and
  empty insert="".
  """
  words = tag.split(":")
  if len(words) > 2:
    raise Exception("Too many : in res tag %s" % tag)
  res_num = words[-1]
  insert = ""
  while not res_num[-1].isdigit():
    insert += res_num[-1] 
    res_num = res_num[:-1]
  res_num = int(res_num)
  if len(words) == 1:
    chain_id = " "
  else:
    chain_id = words[0]
    if len(chain_id) > 1:
      raise Exception("chain_id in res tag %s too long" % tag)
  return (chain_id, res_num, insert)


class Residue:
  """
  Class to collect atoms in a residue together. Allows group
  searching where each atom in a residue must have a unique
  atom_type.
  """

  def __init__(self, in_type, in_chain_id, in_num, in_insert=''):
    self.type = in_type
    self.chain_id = in_chain_id
    self.num = in_num
    self.insert = in_insert
    self._atom_dict = {}
 
  def tag(self):
    """
    Returns a name e.g. "A:12" that combines the chain_id and
    residue number. This is a unique tag that can be used to
    identify a residue in a Soup through get_i_residue().
    """
    tag = ""
    if self.chain_id != " " and self.chain_id != "":
      tag += self.chain_id + ":"
    tag += str(self.num)
    if self.insert:
      tag += self.insert
    return tag  

  def __str__(self):
    atom_name_list = [a.type for a in self.atoms()]
    atom_name = " ".join(atom_name_list)
    return "%s-%s { %s }" % (self.type, self.num, atom_name)

  def copy(self):
    return copy.deepcopy(self)
  
  def n_atom(self):
    return len(self._atom_dict)
    
  def atom(self, atom_type):
    return self._atom_dict[atom_type]
    
  def has_atom(self, atom_type):
    return atom_type in list(self._atom_dict.keys())
    
  def change_atom_type(self, atom_type1, atom_type2):
    if not self.has_atom(atom_type1):
      return
    atom = self._atom_dict[atom_type1]
    atom.type = atom_type2
    del self._atom_dict[atom_type1]
    self._atom_dict[atom_type2] = atom

  def atoms(self):
    return list(self._atom_dict.values())
  
  def atom_name(self, atom_type):
    return self.type + self.num + ":" + atom_type

  def insert_atom(self, atom):
    self._atom_dict[atom.type] = atom
    atom.chain_id = self.chain_id
    atom.res_num = self.num
    atom.res_type = self.type
  
  def erase_atom(self, atom_type):
    del self._atom_dict[atom_type]
    
  def set_num(self, i, insert=""):
    self.num = i
    self.insert = insert
    for atom in self.atoms():
      atom.res_num = self.num
      atom.res_insert = insert
     
  def inc_num(self):
    self.set_num(self.num+1, self.insert)

  def dec_num(self):
    self.set_num(self.num-1, self.insert)
    
  def dec_insert(self):
    l = self.insert;
    if l == "A" or l == "a":
      self.insert = ''
    else:
      i = string.ascii_letters.find(l)
      self.insert = string.ascii_letters[i-1]

  def transform(self, matrix):
     for atom in self.atoms():
       atom.transform(matrix)

  def set_chain_id(self, chain_id):
    self.chain_id = chain_id
    for a in self.atoms():
      a.chain_id = chain_id

  def set_type(self, res_type):
    self.type = res_type
    for a in self.atoms():
      a.res_type = res_type

  def load_bfactor(self, bfactor):
    for atom in self.atoms():
      atom.bfactor = bfactor


class Soup():
  """
  The major class that holds a list of atoms and references them
  to a list of residues.

  The methods residues() and atoms() provide access to the
  data structures. 

  Inserting of residues should be done here, as Soup will 
  administer both atom and residue lists.
  """

  def __init__(self, fname=""):
    self._residues = []
    self._atoms = []
    if fname:
      self.read_pdb(fname)

  def clear(self):
    del self._residues[:]
    for atom in self._atoms:
      del atom
    del self._atoms[:]
    
  def copy(self):
    return copy.deepcopy(self)

  def n_atom(self):
    return len(self._atoms)

  def atoms(self):
    return self._atoms

  def atom(self, i):
    return _atoms[i]

  def insert_atom(self, i, atom):
    self._atoms.append(atom)
    self.residue(i).insert_atom(atom)
    
  def erase_atom(self, i, atom_type):
    atom = self.residue(i).atom(atom_type)
    self.residue(i).erase_atom(atom_type)
    for _atom in self._atoms:
      if _atom == atom:
        self._atoms.remove(atom)
        del atom
        break
    
  def transform(self, matrix):
    for atom in self._atoms:
      atom.transform(matrix)

  def residues(self):
    return self._residues

  def residue(self, i):
    return self._residues[i]
    
  def get_i_residue(self, tag):
    """
    Returns the index of residue with tag, or -1 on failure.
    """
    for i, residue in enumerate(self.residues()):
      if split_tag(tag) == (residue.chain_id, residue.num, residue.insert):
        return i
    raise -1
  
  def residue_by_tag(self, tag):
    i = self.get_i_residue(tag)
    if i >= 0:
      return self.residue(i)
    else:
      raise None

  def n_residue(self):
    return len(self._residues)
    
  def insert_residue(self, i, res):
    is_insertion = False
    if i < self.n_residue()-1:
      save_res_num = self.residue(i).num
      if self.residue(i+1).num == save_res_num:
        is_insertion = True

    if self.n_residue() == 0:
      res.set_num(res.num, res.insert)
    elif i < self.n_residue():
      res.set_num(self.residue(i).num, self.residue(i).insert)
    else:
      res.set_num(self.residue(i-1).num, "")
      res.inc_num()

    self._residues.insert(i, res)
    for atom in res.atoms():
      self.insert_atom(i, atom)

    for j in range(i+1, self.n_residue()):
      self.residue(j).inc_num()

    if is_insertion:
      while self.residue(i+1).insert:
        for j in range(i+1, self.n_residue()):
          if self.residue(j).res_num == save_res_num:
            self.residue(k).dec_insert()
    
  def append_residue(self, res):
    self._residues.append(res)
    for atom in res.atoms():
      self.insert_atom(self.n_residue()-1, atom)

  def erase_residue(self, i):  
    save_res_num = self.residue(i).num

    for atom in self.residue(i).atoms():
      self._atoms.remove(atom)
      del atom
    self._residues.pop(i)  
    
    if i < self.n_residue():
      if self.residue(i).num == save_res_num:
        # erasing residue in an insertion
        for j in range(i, self.n_residue()):
          if self.residue(j).num == erase_res_num_int:
            self.residue(j).dec_insert()
      else:
        for j in range(i, self.n_residue()):
          self.residue(j).dec_num()
    
  def extract_soup(self, i, j):
    extract = Soup()
    for res in self.residues()[i:j]:
      extract.append_residue(res.copy())
    return extract
 
  def insert_soup(self, i, insert):
    for res in reversed(insert.residues()):
      self.insert_residue(i, res.copy())
    
  def chain_ids(self):
    chain_id = [r.chain_id for r in self.residues()]
    return list(set(chain_id))

  def extract_chain(self, chain_id):
    extract = Soup()
    for res in self.residues():
      if res.chain_id == chain_id:
        extract.append_residue(res.copy())
    return extract
 
  def load_residue_bfactors(self, res_bfactors):
    for r, b in zip(self.residues(), res_bfactors):
      r.load_bfactor(b)

  def __str__(self):
    res_name_list = [str(res) for res in self._residues]
    return "\n".join(res_name_list)
 
  def read_pdb(self, fname):
    self.clear()
    res_num = -1
    res_insert = " "
    for line in open(fname, 'r').readlines():
      if line.startswith("ATOM") or line.startswith("HETATM"):
        atom = AtomFromPdbLine(line);
        if (res_num != atom.res_num) or \
           (res_insert != atom.res_insert):
          residue = Residue(
              atom.res_type, atom.chain_id,
              atom.res_num, atom.res_insert)
          self.append_residue(residue)
          res_num = atom.res_num
          res_insert = atom.res_insert
        self.insert_atom(-1, atom)
      if line.startswith(("END", "ENDMDL")):
        return

  def write_pdb(self, pdb):
    f = open(pdb, 'w')
    n_atom = 0
    for res in self.residues():
      res_atoms = res.atoms()
      res_atoms.sort(key=cmp_atom)
      for atom in res_atoms:
        f.write(atom.pdb_str() + '\n')
    f.close()
