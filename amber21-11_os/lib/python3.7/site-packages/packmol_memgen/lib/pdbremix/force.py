# encoding: utf-8

__doc__ = """

Generate forces for steered molecular-dynamics using PUFF.

This module provides functions to generate forces by inducing
velocity changes in Soup objects. This can be saved to restart
files of MD simulations for AMBER, NAMD and GROMACS.

There are four functions that are to be used:

1. make_atd_fn(i_residue, heating_temperature, backbone_atoms)

2. make_puff_fn(
    domain1, domain2, target_val, dt=0.1, temperature=None, 
    is_backbone_only=False, is_first_domain_only=False, 
    force_fname='md.puff.out')

3. make_puff_acc_fn(
    domain1, domain2, target_val, dt=0.1, temperature=None, 
    is_backbone_only=False, is_first_domain_only=False, 
    force_fname='md.puff.out')

4. make_rip_fn(i_res, heating_temperature)

These are all function factories to generate a function in the
form:

   def pulse_fn(soup): 
     # change soup velocities
     
The calculations of velocities and energies used here are 
mainly:

- velocity: Ångstrom/ps
- work/energy: Da*angs/ps^2
- mass: Da

However, various conversions are needed to connect with various
thermostatistical identites as well as for output. 
"""


import os
import random
import math
import copy

from . import util
from . import pdbatoms
from . import v3
from . import data


##########################################################
# Unit conversions

# force = Da⋅Å/ps/ps 
#       = 1.66E-13 kg⋅m/s/s 
#       = 1.66E-13 N 
#       = 1.66E-1 pN

# Work/energy conversions
#      = Da⋅Å⋅Å/ps/ps 
#      = Da⋅Å/ps/ps ⋅ Å 
#      = 1.66E-1 pNÅ
work_DaAngSqPerPsSq_to_pNAng = 1.66E-1
work_Nm_to_kcal = 0.000238846
avogadro = 6.02214179E23
work_pNAng_to_kcalPerMol = 1E-12*1E-10*work_Nm_to_kcal*avogadro

# molecular-dynamics integration time-step
timestep_in_ps = 0.001 

# Velocity conversions
vel_mPerS_to_AngsPerPs = 0.01 
velsq_mPerS_to_AngsPerPs = 1.0E-4 

# Boltzmann constant
#    = 1.3806488E-23 J/K 
#    = 1.3806488E-23 kg ⋅ m^2/s^2/K 
#    = 1.3806488E-23 1/1.66E-27Da ⋅ m^2/s^2/K
#    = 8314.47148 Da⋅m^2/s^2/K
boltzmann_in_DaMSqPerSSqPerK = 8314.47148  


##########################################################
# Heating functions


def average_vel(atoms):
  """
  Returns the mass-averaged velocity of atoms.
  """
  momentum = v3.vector()
  mass = 0.0
  for a in atoms:
    momentum += v3.scale(a.vel, a.mass)
    mass += a.mass
  return v3.scale(momentum, 1.0/mass)


def add_vel_to_atoms(atoms, vel_diff):
  """
  Adds vel_diff to the vel vector of atoms.
  """
  for a in atoms:
    a.vel += vel_diff


def maxwell_velocity(temperature, mass):
  """
  Returns a velocity (in angs/ps) sampled from a Maxwell velocity
  distribution determined by mass and temperature. 
  """
  velsq_ave = boltzmann_in_DaMSqPerSSqPerK * temperature / mass 
  return random.gauss(0, math.sqrt(velsq_ave)) * \
         vel_mPerS_to_AngsPerPs


def mean_energy(temperature, n_degree_of_freedom):
  """
  Returns the average energy (Da*angs/ps^2) of n degree of 
  freedom at temperature. if n_degree_of_freedom = 3,
  this is the average energy of a point particle.
  """
  return 0.5 * n_degree_of_freedom * temperature * \
         boltzmann_in_DaMSqPerSSqPerK * \
         velsq_mPerS_to_AngsPerPs 


def random_energy(temperature, n_degree_of_freedom):
  """
  Returns an energy (Da*angs/ps^2) sampled from a Maxwellian
  distribution of energies at temperature.
  """
  average = mean_energy(temperature, n_degree_of_freedom);
  std_dev = math.sqrt(average)
  return random.gauss(average, std_dev)

  
def kinetic_energy(atoms):
  """
  Returns the kinetic energy (Da*angs/ps^2) of the atoms.
  """
  en = 0.0
  for a in atoms:
    vel = v3.mag(a.vel)
    en += 0.5 * a.mass * vel * vel
  return en


def gas_randomize(atoms, temperature):
  """
  Randomly assigns a velocity to atoms based on a Maxwellian
  distribution at temperature.
  """
  for atom in atoms:
    v3.set_vector(
        atom.vel,
        maxwell_velocity(temperature, atom.mass),
        maxwell_velocity(temperature, atom.mass),
        maxwell_velocity(temperature, atom.mass))


def make_atd_fn(i_residue, heating_temperature, backbone_atoms):
  """
  Returns pulse_fn that locally heats a sidechain.
  """
  # This method is called the Anisotropic Thermal Diffusion and
  # was originally desiged by Ota & Agard "Intramolecular
  # Signaling Pathways Revealed by Modeling Anisotropic Thermal
  # Diffusion" JMB (2005) 12:345.

  def gas_heat_sidechain(
      soup, i_residue, heating_temperature, backbone_atoms):
    atoms = [a for a in soup.residue(i_residue).atoms() 
             if a.type not in backbone_atoms]
    gas_randomize(atoms, heating_temperature)

  return lambda soup: gas_heat_sidechain(
      soup, i_residue, heating_temperature, backbone_atoms)


def anderson_velocity_scale(atoms, temperature, n_degree_of_freedom):
  """
  Scales the velocity of atoms such that average energy
  is consistent with the temperature.
  """
  # This is the classic Anderson approach to temperature
  # regulation. Whilst deterministic, can be easily trapped in
  # local minima.
  target_energy = mean_energy(temperature, n_degree_of_freedom)
  kin = kinetic_energy(atoms)
  if v3.is_similar_mag(kin, 0):
    gas_randomize(atoms, temperature)
  else:
    scaling_factor = math.sqrt(target_energy / kin)
    for atom in atoms:
      v3.set_vector(atom.vel, v3.scale(atom.vel, scaling_factor))



##########################################################
# Pushing functions


def get_atoms_of_residues(soup, res_indices, atom_types=None):
  """
  Return atoms of soup that belong to residues indicated by
  res_indices and in atom_types.
  """
  atoms = []
  for i in res_indices:
    res_atoms = soup.residue(i).atoms()
    if atom_types:
      res_atoms = [a for a in res_atoms if a.type in atom_types]
    atoms.extend(res_atoms)
  return atoms
  

class PushApartByVel():
  """
  Strategy to push apart two domains in a Soup.

  This class is designed to be initialized by make_puff_fn(),
  which returns the apply method of this object to the pulse()
  function of simulate.py 
  """

  def __init__(
      self, domain1, domain2, target_val, dt=0.1,
      temperature=None, is_backbone_only=False, 
      is_first_domain_only=True,
      force_fname='md.puff.out'):
    self.domain1 = domain1
    self.domain2 = domain2
    self.dt = dt
    self.target_val = target_val
    self.temperature = temperature
    self.force_fname = os.path.abspath(force_fname)
    self.is_backbone_only = is_backbone_only  
    self.is_first_domain_only = is_first_domain_only  

  def setup_domains(self):
    # scale the temperature (necessary at high pulling speeds)
    if self.temperature:
      atoms = self.soup.atoms()
      anderson_velocity_scale(atoms, self.temperature, 3*len(atoms))

    # select the atoms from the domains definition
    selection = ['CA'] if self.is_backbone_only else None
    self.atoms1 = get_atoms_of_residues(self.soup, self.domain1, selection)
    self.atoms2 = get_atoms_of_residues(self.soup, self.domain2, selection)

    # get direction vectors based on domains
    self.disp2to1 = pdbatoms.get_center(self.atoms1) \
                  - pdbatoms.get_center(self.atoms2)
    self.axis2to1 = v3.norm(self.disp2to1)

    # calculate relative velocities between domains
    self.vel2to1 = average_vel(self.atoms1) - average_vel(self.atoms2)
    self.axis_vel2to1 = v3.parallel(self.vel2to1, self.axis2to1)
    self.vel = v3.dot(self.axis_vel2to1, self.axis2to1)

    if self.is_first_domain_only:
      self.move_atoms = self.atoms1
    else:
      self.move_atoms = self.atoms1 + self.atoms2

  def change_vels(self):
    # calculate the vel diff vector 
    self.old_kinetic_energy = kinetic_energy(self.move_atoms)

    vel_target = self.target_val*self.axis2to1
    if self.is_first_domain_only:
      change_sets = [(self.atoms1, vel_target)]
    else:
      change_sets = [
        (self.atoms1,  0.5*vel_target),
        (self.atoms2, -0.5*vel_target)]

    # now change velocities of movable atoms
    for move_atoms, vel_target in change_sets:
      for a in move_atoms:
        vel_axis = v3.parallel(a.vel, self.axis2to1)
        vel_diff = vel_target - vel_axis
        a.vel += vel_diff

    self.kinetic_energy = kinetic_energy(self.move_atoms)

  def calculate_output(self):
    work_applied = self.kinetic_energy - self.old_kinetic_energy
    work_applied *= work_DaAngSqPerPsSq_to_pNAng
    work_applied *= work_pNAng_to_kcalPerMol
    self.output_dict = {
      'separation': v3.mag(self.disp2to1),
      'mass': sum(a.mass for a in self.move_atoms),
      'target_vel': self.target_val,
      'vel': self.vel,
      'work_applied': work_applied
    }

  def append_output(self):
    with open(self.force_fname, 'a') as f:
      out_s = str(self.output_dict)
      if not out_s.endswith('\n'):
        out_s += '\n'
      f.write(out_s)

  def apply(self, soup):
    """
    Apply velocity changes to the soup that will induce the
    pulling. This is the key method that will be exported to the
    pulse() function of simulate.py.
    """
    self.soup = soup
    self.setup_domains()
    self.change_vels()
    self.calculate_output()
    self.append_output()


def make_puff_fn(
    domain1, domain2, target_val, dt=0.1, temperature=None, 
    is_backbone_only=False, is_first_domain_only=False, 
    force_fname='md.puff.out'):
  """
  Returns a pulse_fn implemented by PushApartByVel. 
  """
  strategy = PushApartByVel(
      domain1, domain2, target_val, dt, temperature, 
      is_backbone_only, is_first_domain_only, 
      force_fname)
  pulse_fn = lambda soup: strategy.apply(soup)
  return pulse_fn


def read_puff_out(md_dir):
  """
  Yields a dictionary representing the properties of 
  PushApartByVel at each frame of a pulsed simulation from the
  specified md.puff.out file.
  """
  # get time in ps, typical MD step is 0.001 ps = 1 fs
  config = os.path.join(md_dir, 'md.puff.config')
  parms = util.read_dict(config)
  dt = 0.001*parms['n_step_per_pulse']
  time = 0.0
  for line in open(os.path.join(md_dir, 'md.puff.out')):
    entry = eval(line)
    entry['time'] = time
    yield entry
    time += dt


class PushApartByAcc(PushApartByVel):
  """
  Strategy to push apart two domains with constant acceleration.

  This is designed to be instantiated by make_puff_acc_fn().
  """
  def change_vels(self):
    # calculate the vel diff vector to apply
    diff_vel = self.target_val * self.dt
    diff_axis_vel2to1 = v3.scale(self.axis2to1, diff_vel)
    self.vel_diff = v3.dot(diff_axis_vel2to1, self.axis2to1)
    if self.is_first_domain_only:
      add_vel_to_atoms(self.atoms1, diff_axis_vel2to1)
    else:
      # apply half of vel_diff to each domain
      diff_axis_vel2to1 = v3.scale(diff_axis_vel2to1, 0.5)
      add_vel_to_atoms(self.atoms1, diff_axis_vel2to1)
      add_vel_to_atoms(self.atoms2, -diff_axis_vel2to1)

  def calculate_output(self):
    PushApartByVel.calculate_output(self)
    self.output_dict['target_acc'] = self.target_val
    del self.output_dict['target_vel']


def make_puff_acc_fn(
    domain1, domain2, target_val, dt=0.1, temperature=None, 
    is_backbone_only=False, is_first_domain_only=False, 
    force_fname='md.puff.out'):
  """
  Returns pulse_fn that implements PushApartByAcc.
  """
  strategy = PushApartByAcc(
      domain1, domain2, target_val, dt, temperature, 
      is_backbone_only, is_first_domain_only, force_fname)
  pulse_fn = lambda soup: strategy.apply(soup)
  return pulse_fn


##########################################################
# 3. Rotation functions

# Rotational Units
# --------------------------------------------------------
# rotational-velocity:
#    1 º/ps = E+12 º/s 
# rotational-acceleration:
#    1 º/ps/ps = E+24º/s/s
# moment-of-inertia = 1 Da⋅Å⋅Å
#                   = 1.66E-27 kg⋅E-10m⋅E-10m 
#                   = 1.66E-47 kg⋅m^2
# torque = moment-of-inertia * rotational-acceleration
#        = Da⋅Å⋅Å⋅º/ps/ps  
#        = 1.66E-47 kg⋅m^2 ⋅ E+24⋅º/s/s
#        = 1.66E-23 º⋅m⋅kg⋅m/s/s 
#        = 1.66E-23 º⋅m⋅N
# force = torque/radius
#       = Da⋅Å⋅Å⋅º/ps/ps / Å 
#       = 1.66E-23 m⋅N/E-10m
#       = 1.66E-13 N
#       = 1.66 E-1 pN


def moment_of_inertia(atom, axis, anchor):
  """
  Returns the moment (DaAng^^2) of atom around axis at anchor.
  """
  r = atom.pos - anchor
  r_perp = v3.perpendicular(r, axis)
  r_len = v3.mag(r_perp)
  return atom.mass * r_len * r_len


def total_moment_of_inertia(atoms, axis, anchor):
  """
  Returns the total moment (DaAng^^2) of atomss around axis at
  anchor.
  """
  moments = [moment_of_inertia(atom, axis, anchor)
             for atom in atoms]
  return sum(moments)
    

def rotational_velocity(atom, axis, anchor):
  """
  Returns the rotational velocity (rad/ps) of the atom connected
  to anchor around axis.
  """
  r = atom.pos - anchor
  r_perp = v3.perpendicular(r, axis)
  vel_perp = v3.perpendicular(atom.vel, axis)
  vel_tang = v3.perpendicular(vel_perp, r_perp)
  pos_ref = v3.cross(axis, r_perp)
  if v3.dot(vel_tang, pos_ref) < 0.0:
    sign = -1.0
  else:
    sign = 1.0
  if v3.is_similar_mag(v3.mag(r_perp), 0):
    result = 0.0
  else:
    result = sign * v3.mag(vel_tang) / v3.mag(r_perp)
  return result
  

def weighted_rotational_velocity(atoms, axis, anchor):
  """
  Returns the average rotational velocity (rad/ps) of a bunch of
  a atoms weighted by the rotational moments.
  """
  moments = [ \
      moment_of_inertia(atom, axis, anchor) for atom in atoms]
  total_moment = sum(moments)
  weights = [moment / total_moment for moment in moments]
  rot_vels = [ \
      rotational_velocity(atom, axis, anchor) for atom in atoms]
  weighted_rot_vels = [ \
      rot_vel*weight for rot_vel, weight in zip(rot_vels, weights)]
  return sum(weighted_rot_vels)


def add_rotational_velocity(atoms, rot_vel, axis, anchor):
  """
  Adds the rot_vel to the vel vector of atoms with respect
  to the rotation around axis and attached to anchor.
  """
  for atom in atoms:
    r_perp = v3.perpendicular(atom.pos - anchor, axis)
    v_tang_dir = v3.cross(axis, r_perp)
    v_tang_dir_len = v3.mag(v_tang_dir)
    if v3.is_similar_mag(v_tang_dir_len, 0):
      v_tang = v3.vector()
    else:
      v_new_len = rot_vel * v3.mag(r_perp)
      v_tang = v3.scale(v_tang_dir, v_new_len/v_tang_dir_len)
    atom.vel += v_tang
  

def get_n_chi(residue):
  """
  Returns the number of chi angles of residue.
  """
  return len(data.get_res_chi_topology(residue.type))


def calculate_chi(residue, i_chi):
  """
  Returns the angle for the i_chi dihedral angle of residue.
  """
  res_chi_topology = data.get_res_chi_topology(residue.type)
  if i_chi < len(res_chi_topology):
    atom_types = res_chi_topology[i_chi]
    crds = [residue.atom(t).pos for t in atom_types]
    return v3.normalize_angle(v3.dihedral(*crds))
  raise ValueError("No Chi%d angle for residue %d" % (i_chi, i))


def get_axis_anchor(residue, i_chi):
  """
  Returns the axis of rotation and an anchor point of i_chi
  dihedral of residue.
  """
  res_chi_topology = data.get_res_chi_topology(residue.type)
  p = [residue.atom(a).pos for a in res_chi_topology[i_chi]]
  axis = p[2] - p[1]
  anchor = p[2]
  return axis, anchor
  
    
def atoms_affected_by_chi(residue, i_chi):
  """
  Returns the atoms in residue that will be rotated if the i_chi
  dihedral is rotated.
  """
  def sidechain_nesting(atom_type):
    label = atom_type
    while label[-1].isdigit():
      label = label[:-1]
    while label[0].isdigit():
      label = label[1:]
    if len(label) < 2:
      nesting = -1
    else:
      nesting = "ABGDEZH".find(label[1]) - 2
      if label[0] == "H":
        nesting += 1
    return nesting
  return [a for a in residue.atoms()
          if sidechain_nesting(a.type) >= i_chi]


def get_rot_vel_chi(residue, i_chi):
  """
  Returns the weighted rotational velocity of the atoms 
  that are rotated by the i_chi dihedral.
  """
  axis, anchor = get_axis_anchor(residue, i_chi)    
  atoms = atoms_affected_by_chi(residue, i_chi)
  return weighted_rotational_velocity(atoms, axis, anchor)


def get_random_chi_rot_vel(residue, i_chi, temperature):
  """
  Returns a random energy from a Maxwellian energy distribution
  consistent with the number of degrees of freedom of the
  atoms involved in the i_chi dihedral.
  """
  axis, anchor = get_axis_anchor(residue, i_chi)
  atoms = atoms_affected_by_chi(residue, i_chi)
  moment = total_moment_of_inertia(atoms, axis, anchor)
  energy = random_energy(temperature, 3*len(atoms))
  return math.sqrt(2 * energy / moment)


def add_rot_vel_to_chi(residue, i_chi, target_rot_vel):
  """
  Adds target_rot_vel to the atoms affected by i_chi around the
  chi axis.
  """
  axis, anchor = get_axis_anchor(residue, i_chi)
  atoms = atoms_affected_by_chi(residue, i_chi)
  add_rotational_velocity(atoms, target_rot_vel, axis, anchor)


class Rip:
  """
  Strategy to directly rotate a residue via chi angles.

  This is designed to be instantiated by make_rip_fn().
  """
  def __init__(self, i_res, heating_temperature):
    self.i_res = i_res
    self.heating_temperature = heating_temperature
    self.mean_chis = None
    self.max_delta_chi = v3.radians(60)

  def apply(self, soup):
    residue = soup.residue(self.i_res)
    atoms = residue.atoms()
    n_chi = get_n_chi(residue)

    if self.mean_chis is None:
      self.mean_chis = [calculate_chi(residue, i) for i in range(n_chi)]

    rot_vels = [get_rot_vel_chi(residue, i) for i in range(n_chi)]
    for atom in atoms:
      v3.set_vector(atom.vel, 0.0, 0.0, 0.0)

    for i_chi in reversed(list(range(n_chi))):
      chi = calculate_chi(residue, i_chi)
      delta_chi = v3.normalize_angle(chi - self.mean_chis[i_chi])
      target_rot_vel = get_random_chi_rot_vel(
          residue, i_chi, self.heating_temperature)
      if abs(delta_chi) > self.max_delta_chi:
        if delta_chi > self.max_delta_chi:
          target_rot_vel = -target_rot_vel
      else:
        if rot_vels[i_chi] < 0.0:
          target_rot_vel *= -target_rot_vel
      add_rot_vel_to_chi(residue, i_chi, target_rot_vel)

    anderson_velocity_scale(atoms, self.heating_temperature, 3*len(atoms))


def make_rip_fn(i_res, heating_temperature):
  """
  Returns pulse_fn that implements Rip.
  """
  rip = Rip(i_res, heating_temperature)
  return lambda soup: rip.apply(soup)

