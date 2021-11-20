# encoding: utf-8

__doc__ = """

Interface to read molecular-dynamics trajectories. 

Simulation trajectories are read into a Trajectory object
that provides a Soup object, which is updated everytime a
frame is requested. All operations on Soup are available,
including the writing of PDB files, and extracting fragments
for further interrogation.

Trajectory object:
  Attributes:
    basename
    n_frame
    i_frame
    soup
  Methods:
    load_frame

Currently supported trajectories:
- AMBER: .top .trj .vel.trj
- GROMACS: .top .gro .trr
- NAMD: .psf .dcd .vel.dcd
"""


import os
import copy

from . import util
from . import v3
from . import pdbatoms
from . import data
from . import rmsd
from . import force

from . import namd
from . import amber
from . import gromacs



def open_trajectory(basename):
  """
  Returns Trajectory object by guessing the type from basename.
  """
  if os.path.isfile(basename + '.trr'):
    print("Detected GROMACS trajectory")
    return gromacs.Trajectory(basename)
  if os.path.isfile(basename + '.psf'):
    print("Detected NAMD trajectory")
    return namd.Trajectory(basename)
  if os.path.isfile(basename + '.top'):
    print("Detected AMBER trajectory")
    return amber.Trajectory(basename)
  raise ValueError("No trajectories found at all!")


def make_equilibrium_pdb(trajectory, pdb, start=None, end=None):
  """
  Calculates the equilibrium structure between start and end
  frames. Default: last half of trajectory.
  """
  if start is None:
    start = trajectory.n_frame // 2
  if end is None:
    end = trajectory.n_frame
  
  if start < 0 or start >= trajectory.n_frame:
    raise IndexError("Start frame out of range")
  if end < 0 or end > trajectory.n_frame+1:
    raise IndexError("End frame out of range")
  
  # sum_soup stores cumulative values
  sum_soup = trajectory.soup.copy()
  for a in sum_soup.atoms():
    v3.set_vector(a.pos, 0, 0, 0)
  
  n_frame = 0
  for i in range(start, end):
    trajectory.load_frame(i)
    n_frame += 1
    for sum_atom, atom in zip(
        sum_soup.atoms(), trajectory.soup.atoms()):
      sum_atom.pos += atom.pos
  
  for a in sum_soup.atoms():
    x = a.pos.x / float(n_frame)
    y = a.pos.y / float(n_frame)
    z = a.pos.z / float(n_frame)
    v3.set_vector(a.pos, x, y, z)
  
  trajectory.soup.write_pdb(pdb)


def make_pdb_from_trajectory(
    basename, i_frame, out_pdb, res_bfactors=None):
  """
  Make a PDB from the ith frame of a trajectory described by
  basename. Optional residue bfactors are loaded into residues.
  """
  trajectory = open_trajectory(basename)
  trajectory.load_frame(i_frame)
  if res_bfactors is not None:
    trajectory.soup.load_residue_bfactors(res_bfactors)
  trajectory.soup.write_pdb(out_pdb)
  del trajectory



class TrajectoryAnalyzer(object):
  """
  Abstract strategy object to analyze the frames of a trajectory
  and write the results to a text file.
  """
  var_name = 'override'

  def __init__(self, trj, n_frame_per_ps, ref_pdb):
    self.trj = trj
    self.soup = self.trj.soup
    self.n_frame_per_ps = n_frame_per_ps
    
    # for reference frame
    if ref_pdb:
      self.ref_soup = pdbatoms.Soup(ref_pdb)
    else:
      self.ref_soup = self.soup.copy()

    # output files
    fname = trj.basename + '.' + self.var_name + '.per_frame'
    self.file_per_frame = open(fname, 'w')
    fname = trj.basename + '.' + self.var_name + '.per_ps'
    self.file_per_ps = open(fname, 'w')

    # cumul results for ps save
    self.cumul_results = None
    self.n_cumul_frame = 0

  def calculate_results(self):
    """
    To be overriden: return list of results to print out
    """
    return []
  
  def process_frame_on_ps(self):
    for i in range(len(self.cumul_results)):
      self.cumul_results[i] /= float(self.n_cumul_frame)

    # save results for ps calculation
    s = ' '.join(map(str, self.cumul_results))
    t = round(self.trj.i_frame / float(self.n_frame_per_ps))
    self.file_per_ps.write('%.f %s\n' % (t, s))

    # clear for next ps
    if self.trj.i_frame > 0:
      for i in range(len(self.cumul_results)):
        self.cumul_results[i] = 0.0
    self.n_cumul_frame = 0

  def process_frame(self):
    results = self.calculate_results()

    # write results to file_per_frame
    s = " ".join(map(str, results))
    frame = self.trj.i_frame
    self.file_per_frame.write("%d %s\n" % (frame, s))

    # initialize on first frame
    if self.cumul_results is None:
      self.cumul_results = [0.0 for i in range(len(results))]

    # add results to cumul_results
    for i in range(len(results)):
      self.cumul_results[i] += results[i]
    self.n_cumul_frame += 1

  def close(self):
    self.file_per_frame.close()
    self.file_per_ps.close()



class CaRmsdAnalyzer(TrajectoryAnalyzer):
  """
  TrajectoryAnalyzer to calculate C-alpha RMSD.
  """
  var_name = 'rmsd'
  def calculate_results(self):
    return [rmsd.rmsd_of_soups(self.ref_soup, self.soup)]


class KineticEnergyAnalyzer(TrajectoryAnalyzer):
  """
  TrajectoryAnalyzer to calculate kinetic energy of residues.
  """
  var_name = 'kin'
  def calculate_results(self):
    results = []
    for residue in self.soup.residues():
      if residue.type not in data.solvent_res_types:
        atoms = residue.atoms()
        energy = force.kinetic_energy(atoms)/float(len(atoms))
        results.append(energy)
    return results


class ResidueDeviationAnalyzer(TrajectoryAnalyzer):
  """
  TrajectoryAnalyzer to calculate kinetic energy of residues.
  """
  var_name = 'dev'
  def calculate_results(self):
    ref_crds = []
    crds = []
    for ref_residue, residue in zip(self.ref_soup.residues(), self.soup.residues()):
      if residue.type not in data.solvent_res_types:
        if residue.has_atom('CA'):
          ref_crds.append(ref_residue.atom('CA').pos.copy())
          crds.append(residue.atom('CA').pos.copy())
    center = v3.get_center(crds)
    crds = [c - center for c in crds]
    ref_center = v3.get_center(ref_crds)
    ref_crds = [c - ref_center for c in ref_crds]
    rmsd_val, transform_ref_to_this = rmsd.calc_rmsd_rot(ref_crds, crds)
    ref_crds = [v3.transform(transform_ref_to_this, c) for c in ref_crds]
    return [v3.distance(crd, ref_crd) for crd, ref_crd in zip(crds, ref_crds)]


def guess_n_frame_per_ps(basename):
  """
  Returns the n_frame_per_ps of a trajectory by reading any
  .config files that would have been generated using simualte.py.
  """
  config = basename + ".config"
  try:
    params = util.read_dict(config)
    # assuming 1fs time step
    n_step_per_ps = 1000 
    if 'n_step_per_snapshot' in params:
      n_step_per_snapshot = params['n_step_per_snapshot']
    n_frame_per_ps = n_step_per_ps / n_step_per_snapshot
  except:
    n_frame_per_ps = 50
  return n_frame_per_ps
  
  

def analyze_trajectory(
    basename, n_frame_per_ps=None, 
    analyzer_classes=None, ref_pdb=None):
  """"
  Calculates trajectory variables for the residues in the 
  system using the TrajectoryAnalyzer straegy objects in 
  analyzers.
  """
  trj = open_trajectory(basename)

  # Instantiate analyzers
  if n_frame_per_ps is None:
    n_frame_per_ps = guess_n_frame_per_ps(basename)
  if analyzer_classes is None:
    analyzer_classes = \
        [KineticEnergyAnalyzer, CaRmsdAnalyzer, ResidueDeviationAnalyzer]
  analyzers = [a(trj, n_frame_per_ps, ref_pdb) \
               for a in analyzer_classes]

  # Go through frame by frame
  for i_frame in range(trj.n_frame):
    print("Processing frame %d/%d" % (i_frame, trj.n_frame))
    trj.load_frame(i_frame)
    for analyzer in analyzers:
      analyzer.process_frame()
      # A whole ps has been save to file_per_ps
      if (i_frame+1) % n_frame_per_ps == 0 or i_frame == 0:
        analyzer.process_frame_on_ps()

  # At last frame
  trj.soup.write_pdb(basename + '.last.pdb')
  
  for analyzer in analyzers:
    analyzer.close()



