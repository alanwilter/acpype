# encoding: utf-8

__doc__ = """

Common Interface for molecular dynamics packages.

It provides a common API to 3 modules that wrap standard
relatively-free molecular-dynamics packages:

- amber.py
- gromacs.py
- namd.py

The routines are divided into roughly three sections:

1. Read and write restart files
2. Generate restart files from PDB
3. Run simulations from restart files

"""

import os
import shutil
import copy 

from . import util

from . import amber
from . import namd
from . import gromacs



# ##########################################################

# 0. Housekeeping function to figure out simulation package

def get_md_module(force_field):
  """
  Returns the specific interface module that is referenced by
  force_field.
  """
  if force_field.startswith('GROMACS'):
    return gromacs
  elif force_field.startswith('AMBER'):
    return amber
  elif force_field.startswith('NAMD'):
    return namd
  else:
    raise ValueError("unrecognized force-field" + force_field)



# ##########################################################

# 1. Reading and writing restart files

# In PDBREMIX, all simulations require restart files to run.

# 1. topology file (top)
# 2. coordinates file (crds)
# 3. velocities file (vels)

# When read into a Soup object, it is assumed that the units
# in the atom.pos and atom.vel vectors are:

# - positions: angstroms
# - velocities: angstroms/picosecond


def expand_restart_files(force_field, basename):
  """
  Returns expanded restart files (top, crds, vels) with basename
  for the package implied by force_field. No file checking.
  """
  md_module = get_md_module(force_field)
  return md_module.expand_restart_files(basename)


def get_restart_files(basename):
  """
  Returns restart files (top, crds, vels) only if they exist,
  otherwise raises Exception. Will deduce package by 
  file extensions attached to the basename.
  """
  for module in [amber, gromacs, namd]:
    try:
      return module.get_restart_files(basename)
    except util.FileException:
      pass
  raise util.FileException("Couldn't find restart files for " + basename)
  

def soup_from_restart_files(basename, skip_solvent=True):
  """
  Reads a Soup from the restart files.
  """
  for module in [amber, gromacs, namd]:
    try:
      top, crds, vels = module.get_restart_files(basename)
      return module.soup_from_restart_files(
          top, crds, vels, skip_solvent)
    except util.FileException:
      pass
  raise util.FileException("Couldn't find restart files for " + basename)


def write_soup_to_crds_and_vels(force_field, soup, basename):
  """
  From soup, writes out the coordinate/velocities for a given
  packaged that is deduced from the force_field.
  """
  md_module = get_md_module(force_field)
  return md_module.write_soup_to_crds_and_vels(soup, basename)


def convert_restart_to_pdb(basename, pdb):
  """
  Converts restart files with basename into PDB file.
  """
  # Will now try to guess restart files by trying each module in
  # turn. Since the functions throws a FileException if any
  # expected files are missing, catching these will determine
  # failure
  for module in [amber, gromacs, namd]:
    try:
      return module.convert_restart_to_pdb(basename, pdb)
    except util.FileException:
      pass
  raise util.FileException("Couldn't find restart files for " + basename)


# # 2. Generate restart files from PDB

# The restart files used for PDBREMIX assumes a consistent file
# naming.

# To generate a topology file from the PDB file:
# - handles multiple protein chains
# - hydrogens are removed and then regenerated
# - disulfide bonds are identified 
# - charged residue protonation states are auto-detected
# - explicit water in cubic box with 10.0 angstrom buffer
# - counterions to neutralize the system


def pdb_to_top_and_crds(
    force_field, raw_pdb, basename, solvent_buffer=10.0):
  """
  Creates and returns the absolute pathnames to topology and
  coordinate files required to start an MD simulation using the
  package implied in the force_field.
  """
  md_module = get_md_module(force_field)
  top, crd = md_module.pdb_to_top_and_crds(
      force_field, raw_pdb, basename, solvent_buffer)
  return os.path.abspath(top), os.path.abspath(crd)
    

# ##########################################################

# # 3. Run simulations from restart files

# Simulation approach for implicit solvent:
# - optional positional constraints: 100 kcal/mol/angs**2 
# - Langevin thermostat for constant temperature

# Simulation approach for explict water: 
# - optional positional restraints: 100 kcal/mol/angs**2
# - periodic box with PME electrostatics
# - Langevin thermostat for constant temperature
# - Nose-Hoover barometer with flexible box size

# Each package maintains its own files, but all required
# will share a common basename with standard extensions


def fetch_simulation_parameters(
    force_field, top, crds, restraint_pdb, 
    simulation_type, basename, restraint_force=None):
  """
  Returns a dictionary that contains all the high level
  parameters needed to run an MD simulation using the package
  implied in force_field.

  Options for simulation_type:
  1. 'minimization'
  2. 'constant_energy'
  3. 'langevin_thermometer'
  """
  # use a bit of magic to get the dictionary. Pesumably the
  # dictionary with the correct name is present in each
  # simulation module
  md_module = get_md_module(force_field)
  parms_dict_name = '%s_parms' % simulation_type
  parms = getattr(md_module, parms_dict_name).copy()
  # Several fields in parms is common across all packages, these
  # are taken from the parameters of this function:
  parms.update({
    'force_field': force_field,
    'topology': top,
    'input_crds': crds,
    'output_basename': basename,
  })
  if restraint_pdb:
    parms['restraint_pdb'] = restraint_pdb
    if restraint_force is not None:
      parms['restraint_force'] = restraint_force
  return parms


def run_simulation_with_parameters(parms):
  """
  Carries out simulations based on parms
  """
  # For housekeeping, the parms dictionary is written to a
  # .config file. As an Exception is thrown if simulation failed,
  # the existence of an equivalent .config file is an indicator
  # that the simulation has already successfully run.
  config = parms['output_basename'] + ".config"
  if util.is_same_dict_in_file(parms, config):
    print("Skipping: simulation already run.")
    return
  md_module = get_md_module(parms['force_field'])
  md_module.run(parms)
  # No exceptions were thrown - write .config file.
  util.write_dict(config, parms)


def minimize(
    force_field, in_basename, basename, 
    restraint_pdb="", restraint_force=None, n_step=200):
  """
  Runs an energy minimization on the restart files top & crd.

  This is a crucial step as minimization is more robust than
  dynamics calculations. An initial minimization will find a good
  local energy minimum conformation for a  dynamics simulation to
  start. This will avoid spurious initial energy fluctuations for
  future dynamics.
  """
  md_module = get_md_module(force_field)
  top, crds, vels = md_module.get_restart_files(in_basename)
  parms = fetch_simulation_parameters(
      force_field, top, crds, restraint_pdb, 
      'minimization', basename, restraint_force)
  parms['n_step_minimization'] = n_step
  run_simulation_with_parameters(parms)


def langevin_thermometer(
    force_field, in_basename, n_step, temp, basename, 
    n_step_per_snapshot=50, restraint_pdb="", restraint_force=None,
    random_seed=2343):
  """
  Runs a constant temperature simulation using a Langevin
  thermometer.

  There are two basic thermometers in most packages. Anderson
  and Langevin. Anderson simply rescales the energy to satisfy
  the average kinetic energy equation wheras Langevin adds a
  little random force. Langevin thus avoids getting trapped in
  unintended energy minima for the cost of a bit of stochasity.
  """
  md_module = get_md_module(force_field)
  top, crds, vels = md_module.get_restart_files(in_basename)
  parms = fetch_simulation_parameters(
      force_field, top, crds, restraint_pdb, 
      'langevin_thermometer', basename, restraint_force)
  parms['input_vels'] = vels
  parms['random_seed'] = random_seed
  parms['n_step_dynamics'] = n_step
  parms['n_step_per_snapshot'] = n_step_per_snapshot
  parms['temperature_thermometer'] = "%.1f" % temp
  parms['temperature_initial_velocities'] = "%.1f" % temp
  run_simulation_with_parameters(parms)


def constant_energy(
    force_field, in_basename, n_step, basename, 
    n_step_per_snapshot=50, restraint_pdb="", restraint_force=None):
  """
  Runs a constant energy simulation.

  Constant energy simulation are useful if you want to capture
  energy transfers exactly. Normally such simulations are 
  preceeded by a period of thermal regulation using a Langevin
  thermometer.
  """
  md_module = get_md_module(force_field)
  top, crds, vels = md_module.get_restart_files(in_basename)
  parms = fetch_simulation_parameters(
      force_field, top, crds, restraint_pdb, 
      'constant_energy', basename, restraint_force)
  parms['input_vels'] = vels
  parms['n_step_dynamics'] = n_step
  parms['n_step_per_snapshot'] = n_step_per_snapshot
  assert 'temperature_thermometer' not in parms
  assert 'temp_init' not in parms
  run_simulation_with_parameters(parms)


def merge_trajectories(force_field, basename, src_basenames):
  md_module = get_md_module(force_field)
  md_module.merge_trajectories(basename, src_basenames)


def merge_simulations(force_field, basename, sim_dirs):
  """
  Splices together a bunch of simulations, all with the same
  basename, into one large simulation in the current directory.
  """
  if not sim_dirs:
    return
  src_basenames = [os.path.join(s, basename) for s in sim_dirs]
  merge_trajectories(force_field, basename, src_basenames)
    

def pulse(
    force_field, in_basename, basename, n_step, pulse_fn, 
    n_step_per_pulse=100, restraint_pdb="", restraint_force=None):
  """
  Runs a pulse simulation that uses the restart-file modification
  strategy to manage a steered-molecular dynamics simulation.

  The pulsed approacha pplies external forces in pulses, which is
  practically carried out be running short constant-energy
  simulations and directly modifying the restart velocities
  between each simulation.

  Pulse simulations has certain advantages: for instance, the
  system can respond to the forces between pulses, and the
  incredibly flexibility in applying forces. The disadvantage is
  the costly setup which is hopefully, mitigated by this library.

  Reference: Bosco K. Ho and David A. Agard (2010) "An improved
  strategy for generating forces in steered molecular dynamics:
  the mechanical unfolding of titin, e2lip3 and ubiquitin" PLoS
  ONE 5(9):e13068.
  """

  # Grab the simulation prameters for a constant energy
  # simulation. Constant energy is preferred as we want to ensure
  # energy changes come only from our velocity modification.
  top, crds, vels = get_restart_files(in_basename)
  # use dummy top and crds, which will be overriden 
  overall_config_parms = fetch_simulation_parameters(
      force_field, top, crds, restraint_pdb, 
      'constant_energy', basename, restraint_force)
  overall_config_parms.update({
    'input_md_name': in_basename,
    'input_vels': vels,
    'n_step_dynamics': n_step,
    'n_step_per_snapshot': n_step_per_pulse // 2,
    'n_step_per_pulse': n_step_per_pulse
  }) 

  # Check if the simulation has already run as the config
  # file is not written until the very end
  config = basename + ".config"
  if util.is_same_dict_in_file(overall_config_parms, config):
    print("Skipping: pulsing simulation already run.")
    return

  # The overall_config_parms will be written out at the end.
  # We make a copy for internal use
  pulse_parms = copy.deepcopy(overall_config_parms)

  # Calculate steps for each pulse, esp for last step
  n_pulse = pulse_parms['n_step_dynamics'] / n_step_per_pulse
  n_step_list = [n_step_per_pulse for i in range(n_pulse)]
  n_excess_step = pulse_parms['n_step_dynamics'] % n_step_per_pulse
  if n_excess_step > 0:
    n_pulse += 1
    n_step_list.append(n_excess_step)
  
  # Prepare restart files for first step
  pulse_parms['topology'] = os.path.abspath(pulse_parms['topology'])
  in_basename = pulse_parms['input_md_name']
  pulse_parms['input_md_name'] = os.path.abspath(in_basename)

  # Now loop through pulses
  timer = util.Timer()
  save_dir = os.getcwd()
  pulses = ["pulse%d" % i for i in range(n_pulse)]
  for pulse, n_step in zip(pulses, n_step_list):
    print("Pulse: %s/%d" % (pulse, n_pulse))

    os.chdir(save_dir)
    util.goto_dir(pulse)

    pulse_parms['n_step_dynamics'] = n_step

    soup = soup_from_restart_files(pulse_parms['input_md_name'])

    # Apply forces by modifying the velocities directly 
    pulse_fn(soup)

    crds, vels = write_soup_to_crds_and_vels(
        force_field, soup, basename + '.pulse.in')
    pulse_parms['input_crds'] = crds
    pulse_parms['input_vels'] = vels

    run_simulation_with_parameters(pulse_parms)

    # Setup new restart files based on just-finished pulse
    pulse_parms['input_md_name'] = os.path.abspath(basename)

  os.chdir(save_dir)

  merge_simulations(force_field, basename, pulses)

  # cleanup pulses after merging
  util.clean_fname(*pulses)

  # everything worked, no exceptions thrown
  open(basename+'.time', 'w').write(timer.str()+'\n')
  util.write_dict(config, overall_config_parms)


