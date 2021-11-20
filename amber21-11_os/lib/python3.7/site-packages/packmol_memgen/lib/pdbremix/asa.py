# encoding: utf-8

__doc__ = """ 

Calculate the accessible-surface area of atoms.

Uses the simple Shrake-Rupley algorithm, that generates a
relatively uniform density of dots over every atoms and
eliminates those within the sphere of another atom. The remaining
dots is used to calculate the area.

Reference: A. Shrake & J. A. Rupley. "Environment and Exposure to
Solvent of Protein Atoms. Lysozyme and Insulin." J Mol Biol. 79
(1973) 351- 371. """

import math

from . import v3
from . import pdbatoms
from collections import defaultdict


def generate_sphere_points(n):
  """
  Returns list of coordinates on a sphere using the Golden-
  Section Spiral algorithm.
  """
  points = []
  inc = math.pi * (3 - math.sqrt(5))
  offset = 2 / float(n)
  for k in range(int(n)):
    y = k * offset - 1 + (offset / 2)
    r = math.sqrt(1 - y*y)
    phi = k * inc
    points.append(v3.vector(math.cos(phi)*r, y, math.sin(phi)*r))
  return points


def find_neighbor_indices(atoms, probe, k):
  """
  Returns list of indices of atoms within probe distance to atom k. 
  """
  neighbor_indices = []
  atom_k = atoms[k]
  radius = atom_k.radius + probe + probe
  indices = list(range(k))
  indices.extend(list(range(k+1, len(atoms))))
  for i in indices:
    atom_i = atoms[i]
    dist = v3.distance(atom_k.pos, atom_i.pos)
    if dist < radius + atom_i.radius:
      neighbor_indices.append(i)
  return neighbor_indices


def calculate_asa(atoms, probe, n_sphere_point=960):
  """
  Returns the accessible-surface areas of the atoms, by rolling a
  ball with probe radius over the atoms with their radius
  defined.
  """
  sphere_points = generate_sphere_points(n_sphere_point)
 
  const = 4.0 * math.pi / len(sphere_points)
  areas = []
  for i, atom_i in enumerate(atoms):
    
    neighbor_indices = find_neighbor_indices(atoms, probe, i)
    n_neighbor = len(neighbor_indices)
    j_closest_neighbor = 0
    radius = probe + atom_i.radius
    
    n_accessible_point = 0
    for point in sphere_points:
      is_accessible = True
      test_point = v3.scale(point, radius) + atom_i.pos
      cycled_indices = list(range(j_closest_neighbor, n_neighbor))
      cycled_indices.extend(list(range(j_closest_neighbor)))
      
      for j in cycled_indices:
        atom_j = atoms[neighbor_indices[j]]
        r = atom_j.radius + probe
        diff = v3.distance(atom_j.pos, test_point)
        if diff*diff < r*r:
          j_closest_neighbor = j
          is_accessible = False
          break
      if is_accessible:
        n_accessible_point += 1
    
    area = const*n_accessible_point*radius*radius 
    areas.append(area)

  return areas


def make_boxes(a, d_max):
    '''
    Returns dictionary which keys are indecies of boxes (regions)
    with d_max length side and values
    are indicies of atoms belonging to these boxes
    '''
    b = defaultdict(list) # space divided into boxes
    for i in range(len(a)):
        atom = a[i]
        box_coor = tuple(int(math.floor(x / d_max)) for x in atom.pos)
        b[box_coor].append(i)
    return b


def add_bond(a, a1, a2, conn, d_max):
    '''
    If distance between atoms a1 and a2 is less than d_max (neighboring atoms),
    add atoms a1 and a2 in adjacency list conn to each other
    '''
    atom1 = a[a1]
    atom2 = a[a2]
    if v3.mag2(atom1.pos - atom2.pos) <= d_max * d_max:  # connected
        conn[a1].append(a2)
        conn[a2].append(a1)


def neighbor_atoms(b, box):
    '''
    Returns list of atoms from half of neighbouring boxes of the box
    another half is accounted when symmetric (opposite) boxes considered
    '''
    na = [] # list for neighboring atoms
    x, y, z = box # coordinates of the box
    # top layer consisting of 9 boxes
    if (x + 1, y + 1, z +1) in b: na.extend(b[(x + 1, y + 1, z +1)])
    if (x, y + 1, z +1) in b: na.extend(b[(x, y + 1, z +1)])
    if (x + 1, y, z +1) in b: na.extend(b[(x + 1, y, z +1)])
    if (x, y, z +1) in b: na.extend(b[(x, y, z +1)])
    if (x - 1, y + 1, z +1) in b: na.extend(b[(x - 1, y + 1, z +1)])
    if (x + 1, y - 1, z +1) in b: na.extend(b[(x + 1, y - 1, z +1)])
    if (x, y - 1, z +1) in b: na.extend(b[(x, y - 1, z +1)])
    if (x - 1, y, z +1) in b: na.extend(b[(x - 1, y, z +1)])
    if (x - 1, y - 1, z +1) in b: na.extend(b[(x - 1, y - 1, z +1)])
    # half of the middle layer excluding the box itself (4 boxes)
    if (x + 1, y + 1, z) in b: na.extend(b[(x + 1, y + 1, z)])
    if (x, y + 1, z) in b: na.extend(b[(x, y + 1, z)])
    if (x + 1, y, z) in b: na.extend(b[(x + 1, y, z)])
    if (x + 1, y - 1, z) in b: na.extend(b[(x + 1, y - 1, z)])
    return na


def adjacency_list(a, d_max):
    '''
    Returns adjacency list from coordinate file
    in O(len(a)) time
    '''
    b = make_boxes(a, d_max) # put atoms into the boxes with dmax length side
    # now go on boxes and check connections inside 3x3 superboxes
    conn = [[] for i in range(len(a))] # list of bond lengths each atom implicated
    for box in b:
        lb = len(b[box])
        for i in range(lb):
            a1 = b[box][i]
            # check possible connections inside the box
            for j in range(i+1, lb):
                a2 = b[box][j]
                add_bond(a, a1, a2, conn, d_max)
            # check connections with atoms from neighbouring boxes
            na = neighbor_atoms(b, box) # list of such atoms
            for a2 in na:
                add_bond(a, a1, a2, conn, d_max)
    return conn


def find_neighbor_indices_modified(atoms, indices, probe, k):
  """
  Returns list of indices of atoms within probe distance to atom k. 
  """
  neighbor_indices = []
  atom_k = atoms[k]
  radius = atom_k.radius + probe + probe
  for i in indices:
    if i == k: continue
    atom_i = atoms[i]
    dist2 = v3.mag2(atom_k.pos - atom_i.pos) # ToAn
    if dist2 < (radius + atom_i.radius) ** 2: # ToAn
      neighbor_indices.append(i)
  return neighbor_indices


def calculate_asa_optimized(atoms, probe, n_sphere_point=960):
  """
  Returns the accessible-surface areas of the atoms, by rolling a
  ball with probe radius over the atoms with their radius
  defined.
  """
  sphere_points = generate_sphere_points(n_sphere_point)
 
  const = 4.0 * math.pi / len(sphere_points)
  areas = []
  neighbor_list = adjacency_list(atoms, 2 * (probe + max(atoms, key=lambda p: p.radius).radius))
  for i, atom_i in enumerate(atoms):
    
    neighbor_indices = [neig for neig in neighbor_list[i]]
    neighbor_indices = find_neighbor_indices_modified(atoms, neighbor_indices, probe, i) # even further narrow diapazon
    n_neighbor = len(neighbor_indices)
    j_closest_neighbor = 0
    radius = probe + atom_i.radius
    
    n_accessible_point = 0
    for point in sphere_points:
      is_accessible = True
      test_point = v3.scale(point, radius) + atom_i.pos
      cycled_indices = list(range(j_closest_neighbor, n_neighbor))
      cycled_indices.extend(list(range(j_closest_neighbor)))
      
      for j in cycled_indices:
        atom_j = atoms[neighbor_indices[j]]
        r = atom_j.radius + probe
        diff2 = v3.mag2(atom_j.pos - test_point)
        if diff2 < r*r:
          j_closest_neighbor = j
          is_accessible = False
          break
      if is_accessible:
        n_accessible_point += 1
    
    area = const*n_accessible_point*radius*radius 
    areas.append(area)

  return areas
