# encoding: utf-8

__doc__ = """ 
Volume calculator for a list of atoms.
"""

import math
import array

from . import pdbatoms
from . import v3


class Grid:
  """
  Class to implement grid-point elimination.

  A large grid of equidistant points is built. Functions are 
  provided to eliminate any points found within the radius of
  a probe atom. The number of eliminated points is used to
  estimate volume.

  Attributes:
    width (float)
    center (vector)
    spacing (float)

  Methods:
    __init__
    reset
    exclude_sphere
    n_excluded
    write_pdb
  """

  def __init__(self, grid_spacing, width, center):
    self.width = float(width)
    half_width = self.width / 2.0
    self.center = v3.vector(center)
    self.spacing = float(grid_spacing)
    self.inv_spacing = 1.0 / self.spacing

    self.n = 1
    cover = 0
    self.low = v3.vector()
    while cover < half_width:
      self.n += 1
      half_n_point = int(self.n / 2)
      self.low[0] = self.center[0] - half_n_point*self.spacing
      self.low[1] = self.center[1] - half_n_point*self.spacing
      self.low[2] = self.center[2] - half_n_point*self.spacing
      width_1 = abs(self.center[0] - self.low[0])
      high_x = self.low[0] + self.n*self.spacing
      width_2 = abs(high_x - self.center[0])
      cover = min(width_1, width_2)
      
    self.actual_width = self.n*self.spacing
    self.total_volume = self.actual_width**3
    self.total_point = self.n**3
    
    self.array = array.array('b')
    for i in range(self.total_point):
      self.array.append(0)
    self.n_sq = self.n**2
      
    self.x = [self.low[0] + i*self.spacing for i in range(self.n)]
    self.y = [self.low[1] + j*self.spacing for j in range(self.n)]
    self.z = [self.low[2] + k*self.spacing for k in range(self.n)]
      
  def reset(self):
    for i in range(self.total_point):
      self.array[i] = 0
    
  def indices(self, pos):
    return ((pos[0]-self.low[0])*self.inv_spacing,
            (pos[1]-self.low[1])*self.inv_spacing,
            (pos[2]-self.low[2])*self.inv_spacing)

  def pos(self, i, j, k):
    return v3.vector(self.x[i], self.y[j], self.z[k])

  def is_grid_point_near_sphere(self, i, j, k, pos, r_sq):
    d_x = self.x[i] - pos[0]
    d_y = self.y[j] - pos[1]
    d_z = self.z[k] - pos[2]
    return (d_x*d_x + d_y*d_y + d_z*d_z) < r_sq
    
  def int_range(self, low_f, high_f):
    low = max(0, int(math.floor(low_f-1)))
    high = min(self.n, int(math.ceil(high_f) + 2))
    return list(range(low, high))

  def exclude_sphere(self, pos, r):
    low = v3.vector(pos[0] - r, pos[1] - r, pos[2] - r)
    low_i, low_j, low_k = self.indices(low)
    high = v3.vector(pos[0] + r, pos[1] + r, pos[2] + r)
    high_i, high_j, high_k = self.indices(high)
    r_sq = r*r
    for i in self.int_range(low_i, high_i):
      for j in self.int_range(low_j, high_j):
        for k in self.int_range(low_k, high_k):
          l = i*self.n_sq + j*self.n + k 
          if self.array[l] == 0:
            if self.is_grid_point_near_sphere(i, j, k, pos, r_sq):
              self.array[l] = 1

  def n_excluded(self):
    return sum(self.array)
    
  def write_pdb(
      self, pdb, res_type="HOH", atom_type="O", element="O"):
    i_res = 1
    with open(pdb, 'w') as f:
      for i in range(self.n):
        for j in range(self.n):
          for k in range(self.n):
            l = i*self.n_sq + j*self.n + k 
            if self.array[l]:
              atom = pdbatoms.Atom()
              atom.pos = self.pos(i,j,k)
              atom.type = atom_type
              atom.is_hetatm = True
              atom.element = element
              atom.res_type = res_type
              atom.res_num = i_res
              atom.num = i_res
              f.write(atom.pdb_str() + '\n')


def volume(atoms, grid_spacing, pdb="", verbose=True):
  """
  Returns the volume of a list of atoms, and writes a PDB file of
  fictious atoms to fill the volume.
  """
  center = pdbatoms.get_center(atoms)
  width = pdbatoms.get_width(atoms, center) + 4.0
  grid = Grid(grid_spacing, width, center)
  if verbose:
    print("%d atoms, grid %d x %d x %d points, width %.2f angstroms" % \
          (len(atoms), grid.n, grid.n, grid.n, grid.actual_width))
  for atom in atoms:
    grid.exclude_sphere(atom.pos, atom.radius)
  d_volume = float(grid_spacing)**3
  volume = grid.n_excluded()*d_volume
  if verbose:
    print("Volume %.1f angstroms^3 (%d x %.3f angstroms^3)" \
            % (volume, grid.n_excluded(), d_volume))
  if pdb:
    if verbose:
      print("Warning: there will probably be more residues/atoms " \
            "than PDB can uniquely number")
      print(pdb)
    grid.write_pdb(pdb)
  return volume
    
