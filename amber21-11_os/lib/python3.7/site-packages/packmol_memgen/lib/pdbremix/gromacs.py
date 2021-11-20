# encoding: utf-8

__doc__ = """

Interface to the GROMACS molecular-dynamics package.

The library is split into three sections:

1. Read and write restart files
2. Generate restart files from PDB
3. Run simulations from restart files
4. Read trajectories with some post-processing

Copyright (C) 2009, 2014, Bosco K. Ho
"""

import os
import shutil
import copy
import glob
import xdrlib

from . import util
from . import pdbatoms
from . import v3
from . import data
from . import pdbtext
from . import protein


# ##########################################################

# 1. Reading and writing restart files

# In PDBREMIX, restart files for GROMACS are assumed to
# have the naming scheme:

# 1. topology file: sim.top
# 2. coordinate/velocity file: sim.gro

# Parsers have been written to read .top and .gro files into
# Python structures, and to write these back into .top and .gro
# files, and to convert them into .pdb files

# The units used in GROMACS are:
# - positions: nanometers
# - velocities: nanometers/picosecond
# - force constant: kJ/mol/nm^2

# which will be converted into
# - angstroms
# - picoseconds
# - kilocalories


def read_top(top):
  """
  Returns a list of (mass, charge, chain_id) for the atoms 
  in the topology file.
  """
  util.check_output(top)
  lines = open(top).readlines()

  atoms = []
  is_chain_topologies = False
  chain=" "
  top_dir = os.path.dirname(top)
  for l in lines:
    if not is_chain_topologies:
      if 'chain topologies' in l:
        is_chain_topologies = True
      continue
    if l.startswith("#include"):
      itp = l.split()[1][1:-1]
      itp = os.path.join(top_dir, itp)
      if os.path.isfile(itp):
        full_chain_name = os.path.splitext(itp)[0]
        chain = full_chain_name.split('_')[-1]
        these_atoms = read_top(itp, chain)
        atoms.extend(these_atoms)
    if l.startswith(";"):
      break

  is_atoms = False
  qtot = None
  for l in lines:
    if not is_atoms:
      if '[ atoms ]' in l:
        is_atoms = True
      continue
    if l.startswith('['):
      break
    if l.startswith(";"):
      continue    
    if not l.strip():
      continue
    words = l.split()
    n = int(words[0])
    res_num = int(words[2])
    res_type = words[3]
    q = float(words[6])
    mass = float(words[7])
    atoms.append((mass, q, chain))

  return atoms
  

def AtomFromGroLine(line):
  """
  Returns an Atom object from a .gro atom line.
  """
  atom = pdbatoms.Atom()
  atom.res_num = int(line[0:5])
  atom.res_type = line[5:8].strip()
  atom.type = line[10:15].strip(" ")
  atom.element = data.guess_element(
      atom.res_type, line[12:15])
  atom.num = int(line[15:20])
  # 10 x multiplier converts from nm to angstroms
  x = 10.0*float(line[20:28])
  y = 10.0*float(line[28:36])
  z = 10.0*float(line[36:44])
  v3.set_vector(atom.pos, x, y, z)
  if len(line) > 62:
    # 10 x multiplier converts from nm to angstroms
    x = 10.0*float(line[44:52])
    y = 10.0*float(line[52:60])
    z = 10.0*float(line[60:68])
    v3.set_vector(atom.vel, x, y, z)
  return atom


def convert_to_pdb_atom_names(soup):
  """
  For the soup structure, converts residue 
  atom names peculiar to .gro into standard PDB names
  for interacting with other systems.
  """
  for res in soup.residues():
    if res.type == "ILE":
      if res.has_atom("CD"):
        res.change_atom_type("CD", "CD1")
    if res.has_atom("OC2"):
      res.change_atom_type("OC2", "OXT")
    if res.has_atom("OC1"):
      res.change_atom_type("OC1", "O")
    if res.type == "SOL":
      res.set_type("HOH")
      res.change_atom_type("HW1", "1H")
      res.change_atom_type("HW2", "2H")
      res.change_atom_type("OW", "O")
    if res.type in data.solvent_res_types:
      for a in res.atoms():
        a.is_hetatm = True
    for atom in res.atoms():
      if atom.type[-1].isdigit() and atom.type[0] == "H":
        new_atom_type = atom.type[-1] + atom.type[:-1]
        res.change_atom_type(atom.type, new_atom_type)


def convert_to_gromacs_atom_names(soup):
  """
  For writing into .gro files, converts PDB atom names
  into GROMACS specific atom names.
  """
  for res in soup.residues():
    if res.type == "ILE":
      if res.has_atom("CD1"):
        res.change_atom_type("CD1", "CD")
    if res.has_atom("OXT"):
      res.change_atom_type("OXT", "OC2")
      if res.has_atom("O"):
        res.change_atom_type("O", "OC1")
    if res.type == "HOH":
      res.set_type("SOL")
      res.change_atom_type("1H", "HW1")
      res.change_atom_type("2H", "HW2")
      res.change_atom_type("O", "OW")
    if res.type in data.solvent_res_types:
      for a in res.atoms():
        a.is_hetatm = False
    for atom in res.atoms():
      if atom.type[0].isdigit() and atom.type[1] == "H":
        new_atom_type = atom.type[1:] + atom.type[0]
        res.change_atom_type(atom.type, new_atom_type)


def soup_from_top_gro(top, gro, skip_solvent=False):
  """
  Returns a Soup built from GROMACS restart files.
  If skip_solvent=True, will skip all solvent molecules.
  """
  util.check_output(top)
  util.check_output(gro)

  soup = pdbatoms.Soup()
  soup.remaining_text = ""
  soup.n_remaining_text = 0

  atoms = []

  # Read from .gro because .top does not contain water
  # residue information, which is "inferred"
  lines = open(gro, 'r').readlines()
  for i_line, line in enumerate(lines[2:-1]):
    atom = AtomFromGroLine(line)
    if skip_solvent and atom.res_type == "SOL":
      soup.remaining_text = "".join(lines[i_line+2:-1])
      soup.n_remaining_text = len(lines[i_line+2:-1])
      break
    atoms.append(atom)
  soup.box = [float(w) for w in lines[-1].split()]

  for atom, (mass, q, chain_id) in zip(atoms, read_top(top)):
    atom.mass = mass
    atom.charge = q

  curr_res_num = -1
  for a in atoms:
    if curr_res_num != a.res_num:
      res = pdbatoms.Residue(
          a.res_type, a.chain_id, a.res_num)
      soup.append_residue(res.copy())
      curr_res_num = a.res_num
    soup.insert_atom(-1, a)

  convert_to_pdb_atom_names(soup)
  protein.find_chains(soup)

  return soup


def write_soup_to_gro(in_soup, gro):
  soup = in_soup.copy()
  convert_to_gromacs_atom_names(soup)
  f = open(gro, 'w')
  f.write("Generated by gromacs.py\n")
  atoms = soup.atoms()
  n_atom = len(atoms) + soup.n_remaining_text
  f.write(str(n_atom) + '\n')
  for a in atoms: 
    # GRO doesn't care about numbering so wrap when full
    res_num = a.res_num % 100000
    atom_num = a.num % 100000
    # 0.1 x multiplier converts from angs back to nm
    s = "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n" % \
        (res_num, a.res_type, a.type, atom_num, 
         a.pos[0]*0.1, a.pos[1]*0.1, a.pos[2]*0.1,
         a.vel[0]*0.1, a.vel[1]*0.1, a.vel[2]*0.1)
    f.write(s)
  if soup.remaining_text:
    f.write(soup.remaining_text)
  f.write("%10.5f%10.5f%10.5f\n" % (soup.box[0], soup.box[1], soup.box[2]))
  f.close()
  

# The following functions wrap the above functions into a
# standard API that does not explicitly reference GROMACS

def expand_restart_files(basename):
  """Returns expanded restart files based on basename"""
  top = os.path.abspath(basename + '.top')
  crds = os.path.abspath(basename + '.gro')
  vels = '' # dummy file for cross-package interface
  return top, crds, vels


def get_restart_files(basename):
  """Returns restart files only if they exist"""
  top, crds, vels = expand_restart_files(basename)
  util.check_files(top, crds)
  return top, crds, vels


def soup_from_restart_files(top, crds, vels, skip_solvent=False):
  """Reads pdbatoms.Soup object from restart files."""
  return soup_from_top_gro(top, crds, skip_solvent)


def write_soup_to_crds_and_vels(soup, basename):
  """From soup, writes out the coordinate/velocities, used for pulsing"""
  write_soup_to_gro(soup, basename + '.gro')
  return basename + '.gro', ''


def convert_restart_to_pdb(basename, pdb):
  """Converts restart files with basename into PDB file"""
  top, crds, vels = get_restart_files(basename)
  soup = soup_from_restart_files(top, crds, vels)
  soup.write_pdb(pdb)


# ##########################################################

# # 2. Generate restart files from PDB

# The restart files used for PDBREMIX assumes a consistent file naming. 
# For a given basename `sim`, the files are:
# 1. topology file: sim.top
# 2. coordinate/velocity file: sim.gro
# 3. restraint file: sim.pos_re.itp

# To generate a topology file from the PDB file:
# - handles multiple protein chains
# - hydrogens are removed and then regenerated by GROMACS
# - disulfide bonds are auto-detected
# - charged residue protonation states are auto-detected
# - explicit water in cubic box with 10.0 angstrom buffer
# - counterions to neutralize the system
# - GROMACS4.5: AMBER99 force-field
# - GROMACS4.0: GROMOS96 force-field

# Binaries used to generate restart files:
# 1. pdb2gmx - create topologies for atoms in PDB
# 2. editconf - define cubic box
# 3. genbox - add explicit water molecules
# 4. grompp - build .tpr file for energy terms
# 5. genion - add counterions based on .tpr file


def delete_backup_files(tag):
  util.clean_fname(*util.re_glob('*', '^#' + tag))
    

force_field_mdp = """
; Bond parameters
constraints     = hbonds        ; bonds from heavy atom to H, constrained
continuation    = yes           ; first dynamics run
; constraints     = all-bonds     ; all bonds (even heavy atom-H bonds) constrained
; constraint-algorithm = lincs    ; holonomic constraints 
; lincs_iter      = 1             ; accuracy of LINCS
; lincs_order     = 4             ; also related to accuracy

; Neighborsearching
ns_type         = grid          ; search neighboring grid cels
nstlist         = 5             ; 10 fs
rlist           = 1.0           ; short-range neighborlist cutoff (in nm)

; Periodic boundary conditions
pbc             = xyz           ; 3-dimensional periodic boundary conditions (xyz|no)

; Electrostatics
coulombtype     = PME           ; Particle Mesh Ewald for long-range electrostatics
pme_order       = 4             ; cubic interpolation
fourierspacing  = 0.16          ; grid spacing for FFT
rcoulomb        = 1.0           ; short-range electrostatic cutoff (in nm)
rvdw            = 1.0           ; short-range van der Waals cutoff (in nm)

; Add also refcoord-scaling to work with position restraints and pressure coupling
refcoord-scaling = all

; Dispersion correction
DispCorr        = EnerPres      ; account for cut-off vdW scheme
"""


ions_mdp = """
; ions.mdp - used as input into grompp to generate ions.tpr
; Parameters describing what to do, when to stop and what to save
integrator  = steep   ; Algorithm (steep = steepest descent minimization)
emtol       = 1000.0  ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep      = 0.01    ; Energy step size
nsteps      = 50000   ; Maximum number of (minimization) steps to perform
"""


def neutralize_system_with_salt(
    in_top, in_gro, basename, force_field):
  """
  Takes a .top file and adds counterions to neutralize the overall
  charge of system, and saves to `basename.gro`.
  """
  # Calculate overall charege in the .top file
  qtot = sum([q for mass, q, chain in read_top(in_top)])
  counter_ion_charge = -int(round(qtot))
  if counter_ion_charge == 0:
    shutil.copy(in_gro, basename + '.gro')
    return

  # Create a .tpr paramater file for genion to find low-energy sites
  in_mdp = basename + '.salt.grompp.mdp'
  open(in_mdp, 'w').write(ions_mdp + force_field_mdp)
  top = basename + '.top'
  if in_top != top:
    shutil.copy(in_top, top)
  tpr = basename + '.salt.tpr'
  out_mdp = basename + '.mdp'
  data.binary(
      'grompp',
      '-f %s -po %s -c %s -p %s -o %s' \
            % (in_mdp, out_mdp, in_gro, top, tpr),
      basename + '.salt.grompp')
  util.check_files(tpr)

  # Use genion to generate a gro of system with counterions 
  gro = basename + '.gro'
  # Genion requires user input "SOL" to choose solvent for replacement
  input_fname = basename + '.salt.genion.in'
  open(input_fname, 'w').write('SOL')
  # Different versions of Gromacs use different counterions
  charge_str = ""
  if 'GROMACS4.5' in force_field:
    charge_str = " -pname NA -nname CL "
  elif 'GROMACS4.0' in force_field:
    charge_str = " -pname NA+ -nname CL- "
  else:
    raise ValueError("Cannot recognize force_field " + force_field)
  if counter_ion_charge > 0:
    charge_str += " -np %d " % counter_ion_charge
  else:
    charge_str += " -nn %d " % abs(counter_ion_charge)
  log = basename + '.salt.genion.log'
  data.binary(
      'genion', 
      '-g %s -s %s -o %s -p %s -neutral %s' % \
          (log, tpr, gro, top, charge_str),
      basename + '.salt.genion',
      input_fname)
  util.check_files(gro)


def pdb_to_top_and_crds(force_field, pdb, basename, solvent_buffer=10):
  """
  Converts a PDB file into GROMACS topology and coordinate files,
  and fully converted PDB file. These constitute the restart files
  of a GROMACS simulation.
  """
  util.check_files(pdb)
  full_pdb = os.path.abspath(pdb)
  save_dir = os.getcwd()

  # All intermediate files placed into a subdirectory
  util.goto_dir(basename + '.solvate')

  # Remove all but protein heavy atoms in a single clean conformation
  pdb = basename + '.clean.pdb'
  pdbtext.clean_pdb(full_pdb, pdb)

  # Generate protein topology in pdb2gmx_gro using pdb2gmx
  pdb2gmx_gro = basename + '.pdb2gmx.gro'
  top = basename + '.top'
  itp = basename + '_posre.itp'
  # Choose force field based on GROMACS version
  if 'GROMACS4.5' in force_field:
    ff = 'amber99' 
  elif 'GROMACS4.0' in force_field:
    ff = 'G43a1' 
  else:
    raise ValueError("Couldn't work out pdb2gmx for " + force_field)
  args = '-ignh -ff %s -water spc -missing -f %s -o %s -p %s -i %s -chainsep id_or_ter -merge all' \
          % (ff, pdb, pdb2gmx_gro, top, itp)
  data.binary('pdb2gmx', args, basename+'.pdb2gmx')
  util.check_files(pdb2gmx_gro)

  # Now add a box with editconf
  box_gro = basename + '.box.gro'
  solvent_buffer_in_nm = solvent_buffer/10.0 
  data.binary(
      'editconf', 
      '-f %s -o %s -c -d %f -bt cubic' \
          % (pdb2gmx_gro, box_gro, solvent_buffer_in_nm),
      basename+'.box')
  util.check_files(box_gro)

  # Given box dimensions, can now populate with explict waters
  solvated_gro = basename + '.solvated.gro'
  data.binary(
      'genbox',
      '-cp %s -cs spc216.gro -o %s -p %s' \
          % (box_gro, solvated_gro, top),
       '%s.solvated' % basename)
  util.check_files(solvated_gro)

  # Neutralize with counterions using genion to place ions 
  # based on energy parameters processed by grompp 
  gro = basename + '.gro'
  neutralize_system_with_salt(top, solvated_gro, basename, force_field)
  util.check_files(gro)

  # Make a reference PDB file from restart files for viewing and restraints
  convert_restart_to_pdb(basename, basename+'.pdb')

  # Copy finished restart files back into original directory
  fnames = util.re_glob(
      '*', os.path.basename(basename) + r'[^\.]*\.(pdb|itp|gro|mdp|top)$')
  for fname in fnames:
    shutil.copy(fname, save_dir)

  # Cleanup
  delete_backup_files(basename)
  os.chdir(save_dir)

  return top, gro


# ##########################################################

# # 3. Run simulations from restart files

# Simulation approach: 
# - cubic periodic box 
# - optional positional restraints: 100 kcal/mol/angs**2  
# - PME electrostatics on the periodic box
# - Langevin thermostat for constant temperature
# - Nose-Hoover barometer with flexible periodic box size
# - constraints on hydrogen atoms bonded to heavy atoms

# Binaries used:
# 1. grompp - process topology files for .tpr file
# 2. mdrun - run MD using .tpr on the .gro file
#    mpi version - mpiexec -np 8 /...dir.../gromacs-4.0.7/bin/mdrun

# Files for trajectories:
# 1. coordinate trajectory: md.trr
# 2. restart coordinate/velocity: md.gro


minimization_parms = { 
  'topology' : 'in.top', 
  'input_crds' : 'in.gro', 
  'output_basename' : 'min', 
  'force_field': 'GROMACS',
  'restraint_pdb': '',
  'restraint_force': 100.0,
  'n_step_minimization' : 100, 
} 

constant_energy_parms = { 
  'topology' : 'in.top', 
  'input_crds' : 'in.gro', 
  'output_basename' : 'md', 
  'force_field': 'GROMACS',
  'solvent_state': 2,
  'surface_area': 1,
  'restraint_pdb': '',
  'restraint_force': 100.0,
  'n_step_per_snapshot' : 50, 
  'n_step_dynamics' : 1000, 
} 

langevin_thermometer_parms = { 
  'topology' : 'in.top', 
  'input_crds' : 'in.gro', 
  'output_basename' : 'md', 
  'force_field': 'GROMACS',
  'restraint_pdb': '',
  'restraint_force': 100.0,
  'random_seed' : 2342, 
  'temperature_thermometer' : 300.0, 
  'temperature_initial_velocities': 0.0, # ignored if it is 0.0
  'n_step_per_snapshot' : 50, 
  'n_step_dynamics' : 1000, 
} 

minimization_mdp = """
; template .mdp file used as input into grompp to generate energy minimization for mdrun

; Parameters describing what to do, when to stop and what to save
integrator  = steep    ; Algorithm (steep = steepest descent minimization)
emtol       = 1000.0   ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep      = 0.01     ; Energy step size
nsteps      = %(n_step_minimization)s    ; Maximum number of (minimization) steps to perform
"""

dynamics_mdp = """
title           = Template for constant temperature/pressure

; Run parameters
integrator      = md            ; leap-frog integrator
nsteps          = %(n_step_dynamics)s  ; time = n_step_dynamics*dt
dt              = 0.001         ; time-step in fs

; Output control
nstxout         = %(n_step_per_snapshot)s  ; save coordinates every 0.05 ps
nstvout         = %(n_step_per_snapshot)s  ; save velocities every 0.05 ps
nstenergy       = %(n_step_per_snapshot)s  ; save energies every 0.05 ps
nstlog          = %(n_step_per_snapshot)s  ; update log file every 0.05 ps

; Pressure coupling is on
pcoupl          = Parrinello-Rahman   ; Pressure coupling on in NPT
pcoupltype      = isotropic     ; uniform scaling of box vectors
tau_p           = 2.0           ; time constant, in ps
ref_p           = 1.0           ; reference pressure, in bar
compressibility = 4.5e-5        ; isothermal compressibility of water, bar^-1
"""

temp_mdp = """
; Temperature coupling is on
tcoupl          = V-rescale     ; modified Berendsen thermostat
tc-grps         = Protein Non-Protein   ; two coupling groups - more accurate
tau_t           = 0.1   0.1     ; time constant, in ps
ref_t           = %(temperature_thermometer)s  %(temperature_thermometer)s   ; reference temperature, one for each group, in K
ld-seed         = %(random_seed)s
"""

vel_mdp = """
; Velocity generation
gen_vel          = yes       ; assign velocities from Maxwell distribution
gen_temp         = %(temperature_initial_velocities)s     ; temperature for Maxwell distribution
gen_seed         = -1        ; generate a random seed
"""

restraint_mdp = """
; Restaints turned on
define            = -DPOSRES  ; position restrain the protein
"""

def make_mdp(parms):
  """
  Using GROMACS, we will always solvate, use periodic 
  boundaries and Particle Ewald Mesh.
  """
  if 'n_step_minimization' in parms:
    mdp = minimization_mdp 
  else:
    mdp = dynamics_mdp 
    if 'temperature_thermometer' in parms:
      mdp += temp_mdp 
    else:
      mdp += "; Temperature coupling is off\n"
      mdp += "tcoupl          = no\n"
    if 'temperature_initial_velocities' in parms and parms['temperature_initial_velocities'] > 0.0:
      mdp += vel_mdp 
    else:        
      mdp += "; Velocity generation\n"
      mdp += "gen_vel         = no            ; Velocity generation is off \n"
  if 'restraint_pdb' in parms and parms['restraint_pdb']:
    mdp += restraint_mdp
  mdp += force_field_mdp
  return mdp % parms
  
  
def replace_include_file(chain, r1, r2):
  lines = open(chain, 'r').readlines()
  new_lines = []
  for l in lines:
    if l.startswith('#include'):
      new_lines.append(l.replace(r1, r2))
    else:
      new_lines.append(l)
  open(chain, 'w').write(''.join(new_lines))


restraint_header = """
; In this topology include file, you will find position restraint
; entries for all the heavy atoms in your original pdb file.
; This means that all the protons which were added by pdb2gmx are
; not restrained.

[ position_restraints ]
; atom  type      fx      fy      fz
"""

def make_restraint_itp(restraint_pdb, force):
  txt = restraint_header
  atoms = pdbatoms.Soup(restraint_pdb).atoms()
  for i, atom in enumerate(atoms):
    if atom.bfactor > 0.0:
      txt += "%6s     1 %5.f %5.f %5.f\n" % (i+1, force, force, force)
  return txt


def run(in_parms):
  """
  Run a GROMACS simulations using the PDBREMIX parms dictionary.
  """
  parms = copy.deepcopy(in_parms)
  basename = parms['output_basename']

  # Copies across topology and related *.itp files, with appropriate
  # filename renaming in #includes
  top = basename + '.top'
  in_top = parms['topology']
  shutil.copy(in_top, top)
  in_name = os.path.basename(in_top).replace('.top', '')
  in_dir = os.path.dirname(in_top)
  file_tag = "%s/%s_*itp" % (in_dir, in_name)
  new_files = [top]
  for f in glob.glob(file_tag):
    new_f = os.path.basename(f)
    new_f = new_f.replace(in_name, basename)
    shutil.copy(f, new_f)
    new_files.append(new_f)
  for f in new_files:
    replace_include_file(f, in_name + "_", basename + "_")

  # Copy over input coordinates/velocities
  in_gro = basename + '.in.gro'
  shutil.copy(parms['input_crds'], in_gro)

  # Generates a postiional-restraint topology file
  if parms['restraint_pdb']:
    # 1kcal*mol*A**-2 = 4.184 kJ*mol*(0.1 nm)**-2 
    kcalmolang2_to_kJmolnm2 = 400.184
    open(basename + '_posre.itp', 'w').write(
        make_restraint_itp(
            parms['restraint_pdb'], 
            parms['restraint_force'] * kcalmolang2_to_kJmolnm2))

  # Generate .mdp file based on parms
  in_mdp = basename + '.grompp.mdp'
  open(in_mdp, 'w').write(make_mdp(parms))

  # Now run .grompp to generate this .tpr file 
  tpr = basename + '.tpr'
  # .mdp to save complete set of parameters
  mdp = basename + '.mdrun.mdp'
  data.binary(
      'grompp',
      '-f %s -po %s -c %s -p %s -o %s' \
          % (in_mdp, mdp, in_gro, top, tpr),
      basename + '.grompp')
  util.check_files(tpr)

  # Run simulation with the .tpr file
  data.binary(
      'mdrun',
      '-v -deffnm %s' % (basename),
      basename + '.mdrun')
  top, crds, vels = get_restart_files(basename)
  util.check_output(top)
  util.check_output(crds)
  
  # Cleanup
  delete_backup_files(basename)


# ##########################################################

# # 4. Read trajectories with some post-processing

# The units used in these files are:

# - positions: nanometers
# - velocities: nanometers/picosecond
# - force constant: kJ/mol/nm^2


n_dim = 3

class TrrReader:
  """
  Class to read the coordinates of a GROMACS .trr file.

  Attributes:
    trr (str) - name of trajectory file
    file (file) - file object to trajectory
    n_atom (int) - number of atoms simulated
    n_frame (int) - number of frames in trajectory
    i_frame (int) - index of current frame
    frame (array) - container of coordinates of current frame

  Methods:
    __init__ - initializes trajectory and loads 1st frame
    load_frame(i) - loads the i'th frame
    __getitem__ - returns the frame
    save_to_crd - save current frame to a .crd file
    __repr__ - string representation
  """

  def __init__(self, trr):
    self.trr = trr
    self.file = open(self.trr, 'r')
    self.read_header()
    self.calc_precision()
    self.calc_frame_info()
    self.i_frame = None
    self.load_frame(0)

  def read_header(self):
    self.u = xdrlib.Unpacker(self.file.read(200))
    self.magic = self.u.unpack_int()
    self.version = self.u.unpack_string()
    self.size_ir = self.u.unpack_int()
    self.size_e = self.u.unpack_int()
    self.size_box = self.u.unpack_int()
    self.size_vir = self.u.unpack_int()
    self.size_pres = self.u.unpack_int()
    self.size_top = self.u.unpack_int()
    self.size_sym = self.u.unpack_int()
    self.size_x = self.u.unpack_int()
    self.size_v = self.u.unpack_int()
    self.size_f = self.u.unpack_int()
    self.n_atom = self.u.unpack_int()
    self.step = self.u.unpack_int()
    self.nre = self.u.unpack_int()
    self.t = self.u.unpack_float()
    self.lamb = self.u.unpack_float()
    self.pos_after_header = self.u.get_position()

  def calc_precision(self):
    "Returns 4 for single precision, and 8 for double precision"
    if self.size_box:
      self.precision = self.size_box/n_dim/n_dim
    elif self.size_x:
      self.precision = self.size_x/n_dim
    elif self.size_v:
      self.precision = self.size_v/n_dim
    elif self.size_f:
      self.precision = self.size_f/n_dim
    else:
      raise ValueError("Cannot determine precision")
    if self.precision not in [4, 8]:
      raise ValueError("Precision not single or double!")
    
  def calc_frame_info(self):
    """Given header info and precision, can calculate frame & n_frame"""
    n_vec = 0
    if self.size_box: n_vec += n_dim
    if self.size_vir: n_vec += n_dim
    if self.size_pres: n_vec += n_dim
    if self.size_x: n_vec += self.n_atom
    if self.size_v: n_vec += self.n_atom
    if self.size_f: n_vec += self.n_atom
    self.size_frame = n_vec*n_dim*self.precision + self.pos_after_header

    # Calculates n_frame from end of file
    self.file.seek(0, 2)
    size_file = self.file.tell()
    self.n_frame = size_file / self.size_frame

  def __repr__(self):
    return "< Gromacs TRR Coord file %s with %d frames of %d atoms >" % \
             (self.trj, self.n_frame, self.n_atom)

  def next_3_reals(self):
    if self.precision == 4:
      return [self.u.unpack_float() for i in range(3)]
    if self.precision == 8:
      return [self.u.unpack_double() for i in range(3)]

  def load_frame(self, i_frame):
    if i_frame < - 1*self.n_frame or i_frame >= self.n_frame:
      raise IndexError
    if i_frame < 0:
      i_frame = self.n_frame + i_frame
    if i_frame == self.i_frame:
      return
    self.file.seek(self.pos_after_header + i_frame*self.size_frame)
    self.u = xdrlib.Unpacker(self.file.read(self.size_frame))
    box, positions, velocities, forces = None, None, None, None
    if self.size_box:
      box = [self.next_3_reals() for i in range(n_dim)]
    if self.size_vir:
      dummy = [self.next_3_reals() for i in range(n_dim)]
    if self.size_pres:
      dummy = [self.next_3_reals() for i in range(n_dim)]
    if self.size_x:
      positions = [self.next_3_reals() for i in range(self.n_atom)]
    if self.size_v:
      velocities = [self.next_3_reals() for i in range(self.n_atom)]
    if self.size_f:
      forces = [self.next_3_reals() for i in range(self.n_atom)]
    self.frame = box, positions, velocities, forces
    self.i_frame = i_frame

  def __getitem__(self, i_frame):
    self.load_frame(i_frame)
    return self.frame


class SoupTrajectory():
  """
  Class to interact with an GROMACS trajctory using soup.
  
  Attributes:
    soup (Soup) - Soup object holding current coordinates/velocities
    trr (str) - coordinate/velocity trajectory file
    trr_reader (TrrReader) - the reader of the frames
    n_frame (int) - number of frames in trajectory
    i_frame (int) - index of current frame

  Methods:
    __init__ - load coordinate and velocity trajectories and build soup
    load_frame - loads new frame into soup
  """

  def __init__(self, soup, trr):
    self.soup = soup
    self.trr = trr
    self.trr_reader = TrrReader(self.trr)
    self.n_frame = self.trr_reader.n_frame
    self.atoms = self.soup.atoms()
    self.load_frame(0)

  def __repr__(self):
    return "< Gromacs Trajectory object with %d frames of %d atoms >" % \
             (self.trj, self.n_frame, self.n_atom)

  def load_frame(self, i_frame):
    box, positions, velocities, forces = self.trr_reader[i_frame]
    for i, atom in enumerate(self.atoms):
      v3.set_vector(
          atom.pos,
          positions[i][0]*10,
          positions[i][1]*10,
          positions[i][2]*10)
      v3.set_vector(
          atom.vel,
          velocities[i][0]*10,
          velocities[i][1]*10,
          velocities[i][2]*10)
    self.i_frame = self.trr_reader.i_frame


class Trajectory(object):
  """
  Class to interact with an GROMACS trajctory using soup.
  
  Attributes:
    basename (str) - basename used to guess all required files
    top (str) - topology file of trajectory
    gro (str) - coordinate/velocity restart file
    trr (str) - coordinate/velocity trajectory file
    trr_reader (TrrReader) - the reader of the frames
    n_frame (int) - number of frames in trajectory
    i_frame (int) - index of current frame
    soup (Soup) - Soup object holding current coordinates/velocities

  Methods:
    __init__ - load coordinate and velocity trajectories and build soup
    load_frame - loads new frame into soup
  """

  def __init__(self, basename):
    self.basename = basename
    self.top = basename + '.top'
    self.gro = basename + '.gro'
    self.soup = soup_from_top_gro(self.top, self.gro)
    self.trr = basename + '.trr'
    self.soup_trj = SoupTrajectory(self.soup, self.trr)
    self.n_frame = self.soup_trj.n_frame
    self.i_frame = None
    self.load_frame(0)

  def load_frame(self, i):
    self.soup_trj.load_frame(i)
    self.i_frame = self.soup_trj.i_frame


def merge_trajectories(basename, traj_basenames):
  """
  Given a bunch of directories with consecutive trajectories, all
  with the same basename, the function will splice them into a
  single  trajctory with basename in the current directory.
  """
  save_dir = os.getcwd()

  trr_fname = basename + '.trr'
  trr_list = [b + '.trr' for b in traj_basenames]
  util.check_files(*trr_list)

  f = open(trr_fname, 'w')
  for trr in trr_list:
    trr = TrrReader(trr)
    for i_frame in range(trr.n_frame-1):
      trr.file.seek(trr.size_frame*i_frame)
      txt = trr.file.read(trr.size_frame)
      f.write(txt)
  f.close()

  # Copy parameters of last pulse into current directory
  traj_basename = traj_basenames[-1]
  for ext in ['.top', '.itp', '.tpr', '.mdrun.mdp', 
              '.grompp.mdp', '.gro']:
    for f in glob.glob('%s*%s' % (traj_basename, ext)):
      g = f.replace(traj_basename, basename)
      shutil.copy(f, g)
      if g.endswith('.top'):
        replace_include_file(g, traj_basename + "_", basename + "_")


  os.chdir(save_dir)
    




