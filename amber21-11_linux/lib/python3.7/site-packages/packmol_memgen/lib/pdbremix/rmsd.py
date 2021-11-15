# encoding: utf-8

__doc__ = """ 

RMSD calculators and optimal rotation of sets of coordinates.

Two algorithms are provided: 

1. The standard SVD decomposition using the numpy library.
2. The qcp quaternion-based method that doesn't need numpy.

The qcp method is particularly useful because it allows the
library to function without numpy such as when run with pypy.
"""

from . import v3
from . import pdbatoms
import math
from .lib import pyqcprot
try:
  import numpy
  is_numpy = True
except:
  is_numpy = False

standard_residues = {"CYS","MET","HIS","HSD","HIE","SER","GLN","ASP","GLU","TYR","THR","ALA","LEU","ILE","PHE","TRP","ARG","ASN","LYS","VAL","PRO","GLY"}

def numpy_svd_rmsd_rot(in_crds1, in_crds2):
  """
  Returns rmsd and optional rotation between 2 sets of [nx3] arrays.
  
  This requires numpy for svd decomposition.
  The transform direction: transform(m, ref_crd) => target_crd.
  """

  crds1 = numpy.array(in_crds1)
  crds2 = numpy.array(in_crds2)
  assert(crds1.shape[1] == 3)
  assert(crds1.shape == crds2.shape)
  

  n_vec = numpy.shape(crds1)[0]
  correlation_matrix = numpy.dot(numpy.transpose(crds1), crds2)
  v, s, w = numpy.linalg.svd(correlation_matrix)
  is_reflection = (numpy.linalg.det(v) * numpy.linalg.det(w)) < 0.0

  if is_reflection:
    s[-1] = - s[-1]
  E0 = sum(sum(crds1 * crds1)) + \
       sum(sum(crds2 * crds2))
  rmsd_sq = (E0 - 2.0*sum(s)) / float(n_vec)
  rmsd_sq = max([rmsd_sq, 0.0])
  rmsd = numpy.sqrt(rmsd_sq)

  if is_reflection:
    v[-1,:] = -v[-1,:]
  rot33 = numpy.dot(v, w)
  m = v3.identity()
  m[:3,:3] = rot33.transpose()

  return rmsd, m


def pyqcprot_rmsd_rot(crds1, crds2):
  """
  Returns rmsd and optional rotation between 2 sets of [nx3] arrays.
  
  This requires Joshua Adelman's pyqcrot library for quaternion-based
  calculation of Theobauld. The transform direction: 
      transform(m, ref_crd) => target_crd.
  """
  rms, rot9 = lib.pyqcprot.calc_rms_rot(crds1, crds2)
  matrix = v3.identity()
  for i in range(3):
    for j in range(3):
       v3.matrix_elem(matrix, i, j, rot9[i*3+j])
  return rms, matrix


def calc_rmsd_rot(crds1, crds2):
  """
  Returns the rmsd and optimal rotation and chooses the method
  depending on whether numpy is installed.
  """
  if is_numpy:
    return numpy_svd_rmsd_rot(crds1, crds2)
  else:
    return pyqcprot_rmsd_rot(crds1, crds2)


def sum_rmsd(crds1, crds2):
  """
  Returns the direct rmsd between two sets of vectors *without*
  doing any optimal rotation. If calculated between optimal sets,
  should give the proper RMSD.
  """
  sum_squared = 0.0
  for crd1, crd2 in zip(crds1, crds2):
    sum_squared += v3.distance(crd1, crd2)**2
  return math.sqrt(sum_squared/float(len(crds1)))
  

def get_superposable_atoms(soup, segments, atom_types, standard=False):
  """
  Returns a list of atom indices to a soup, built from segments.

  Args:
    segments (list): list of pairs of residue names in the soup,
                     such as ['A:1','A:3'], interpreted as the 
                     two ends of a fragment in soup that we want
                     the atom index of
    atom_types (list): list of atom_types in the residues that
                     we want to generate the indices from.  
  """
  result = []
  allowed_i = []
  residues = soup.residues()
  if segments:
    for res_tag_i, res_tag_j in segments:
      i = soup.get_i_residue(str(res_tag_i))
      j = soup.get_i_residue(str(res_tag_j))
      if i > j:
        i, j = j, i
      allowed_i.extend(list(range(i,j+1)))
  else:
    allowed_i = list(range(len(residues)))
  for i, residue in enumerate(residues):
    if i in allowed_i:
      if standard:
        result.extend([a for a in residue.atoms()
                       if a.type in atom_types and residue.type in standard_residues])
      else:
        result.extend([a for a in residue.atoms()
                       if a.type in atom_types])
  return result


def rmsd_of_soups(
    soup1, soup2, segments1=[], segments2=[], 
    atom_types=['CA'], transform_pdb1=None, standard=False):
  """
  Returns the RMSD between two PDB structures and optionally
  writes the best transformed structure of pdb1 in transform_pdb.

  By default, it chooses the CA atoms in the soup.

  Args:
    segments1 (list): list of pairs of residue names in pdb1,
                     such as ['A:1','A:3'], interpreted as the 
                     two ends of a fragment in soup that we want
                     the atom index of
    segments2 (list): same as above but for pdb2
    atom_types (list): list of atom_types in the residues that
                       we want to generate the indices from.  
  """
  atoms1 = get_superposable_atoms(soup1, segments1, atom_types, standard)
  atoms2 = get_superposable_atoms(soup2, segments2, atom_types, standard)

  crds1 = [a.pos for a in atoms1]
  crds2 = [a.pos for a in atoms2]
 
  center1 = v3.get_center(crds1)
  center2 = v3.get_center(crds2)

  soup1.transform(v3.translation(-center1))
  soup2.transform(v3.translation(-center2))

  rmsd, transform_1_to_2 = calc_rmsd_rot(crds1, crds2)

  if not transform_pdb1:
    return rmsd

  soup1.transform(transform_1_to_2)

  soup1.transform(v3.translation(center2))
  soup2.transform(v3.translation(center2))

  soup1.write_pdb(transform_pdb1)
  return sum_rmsd(crds1, crds2)


def rmsd_of_pdbs(
    pdb1, pdb2, segments1=[], segments2=[], 
    atom_types=['CA'], transform_pdb1=None, standard=False):
  """
  Returns the RMSD between two PDB structures and optionally
  writes the best transformed structure of pdb1 in transform_pdb.

  Args:
    segments1 (list): list of pairs of residue names in pdb1,
                     such as ['A:1','A:3'], interpreted as the 
                     two ends of a fragment in soup that we want
                     the atom index of
    segments2 (list): same as above but for pdb2
    atom_types (list): list of atom_types in the residues that
                       we want to generate the indices from.  
  """
  soup1 = pdbatoms.Soup(pdb1)
  soup2 = pdbatoms.Soup(pdb2)
  return rmsd_of_soups(
    soup1, soup2, segments1, segments2, 
    atom_types, transform_pdb1, standard)
  
