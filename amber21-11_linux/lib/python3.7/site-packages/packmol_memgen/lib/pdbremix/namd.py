# encoding: utf-8

__doc__ = """

Interface to the NAMD molecular-dynamics package.

The library is split into three sections:

1. Read and write restart files
2. Generate restart files from PDB
3. Run simulations from restart files
4. Read trajectories with some post-processing

Copyright (C) 2009, 2014, Bosco K. Ho
"""

import os
import sys
import copy
import struct
import shutil

from . import util
from . import data
from . import pdbatoms
from . import v3
from . import pdbtext


# ##########################################################

# 1. Reading and writing restart files

# In PDBREMIX, restart files for NAMD are assumed to
# have the naming scheme, as chosen by namd2 restart
# conventions:

# 1. topology file: sim.psf
# 2. coordinate file: sim.coor
# 2. velocity file: sim.vel

# Parsers have been written to read .psf. The .coor and .vel
# files are simply PDB files. These are the same file formats
# used by CHARMM.

# Standard units used in NAMD are:
# - positions: angstroms
# - velocities: angstroms/picosecond
# - mass: Da
# - charge: electron-charge


def soup_from_psf(psf):
  """
  Returns a Soup from a .psf file
  """
  soup = pdbatoms.Soup()
  curr_res_num = None
  is_header = True
  for line in open(psf):
    if is_header:
      if "NATOM" in line:
        n_atom = int(line.split()[0])
        is_header = False
      continue
    words = line.split()
    atom_num = int(words[0])
    chain_id = words[1]
    res_num = int(words[2])
    res_type = words[3]
    atom_type = words[4]
    charge = float(words[6])
    mass = float(words[7])
    if chain_id.startswith('WT') or chain_id.startswith('ION'):
      is_hetatm = True
      chain_id = " "
    else:
      is_hetatm = False
      chain_id = chain_id[0]
    if curr_res_num != res_num:
      res = pdbatoms.Residue(res_type, chain_id, res_num)
      soup.append_residue(res)
      curr_res_num = res_num
    atom = pdbatoms.Atom()
    atom.vel = v3.vector()
    atom.chain_id = chain_id
    atom.is_hetatm = is_hetatm
    atom.num = atom_num
    atom.res_num = res_num
    atom.res_type = res_type
    atom.type = atom_type
    atom.mass = mass
    atom.charge = charge
    atom.element = data.guess_element(res_type, atom_type)
    soup.insert_atom(-1, atom)
    if len(soup.atoms()) == n_atom:
      break
  convert_to_pdb_atom_names(soup)
  return soup


def convert_to_pdb_atom_names(soup):
  """
  For the soup structure, converts residue 
  atom names peculiar to .coord into standard PDB names
  for interacting with other systems.
  """
  for res in soup.residues():
    if res.type == "ILE":
      if res.has_atom('CD'):
        res.change_atom_type('CD', 'CD1')
    if res.has_atom('OT2'):
      res.change_atom_type('OT2', 'OXT')
      if res.has_atom('OT1'):
        res.change_atom_type('OT1', 'O')
    if res.type == "TIP3":
      res.set_type("HOH")
      res.set_chain_id(" ")
    if res.type == "HSE":
      res.set_type("HIS")
    for atom in res.atoms():
      if atom.type[-1].isdigit() and atom.type[0] == "H":
        new_atom_type = atom.type[-1] + atom.type[:-1]
        res.change_atom_type(atom.type, new_atom_type)


def convert_to_namd_atom_names(soup):
  """
  For writing into .coor files, converts PDB atom names
  into NAMD specific atom names.
  """
  for res in soup.residues():
    if res.type == "ILE" and res.has_atom('CD1'):
      res.change_atom_type('CD1', 'CD')
    if res.has_atom('OXT'):
      res.change_atom_type('OXT', 'OT2')
      if res.has_atom('O'):
        res.change_atom_type('O', 'OT1')
    for atom in res.atoms():
      if atom.type[0].isdigit() and atom.type[1] == "H":
        new_atom_type = atom.type[1:] + atom.type[0]
        res.change_atom_type(atom.type, new_atom_type)
    if res.type == "HOH":
      res.set_type("TIP3")
      res.set_chain_id("W")
    if res.type == "HIS":
      res.set_type("HSE")


# The following functions wrap the above functions into a
# standard API that does not explicitly reference GROMACS


def expand_restart_files(basename):
  """
  Returns expanded restart files based on basename
  """
  psf = os.path.abspath(basename + '.psf')
  coor = os.path.abspath(basename + '.coor')
  vel = os.path.abspath(basename + '.vel')
  return psf, coor, vel


def get_restart_files(basename):
  """
  Returns restart files only if they exist
  """
  psf, coor, vel = expand_restart_files(basename)
  util.check_files(psf, coor)
  if not os.path.isfile(vel):
    vel = ''
  return psf, coor, vel


def soup_from_restart_files(psf, in_coor, in_vel='', skip_solvent=False):
  """
  Reads a Soup from restart files.
  """
  soup = soup_from_psf(psf)
  coord_soup = pdbatoms.Soup(in_coor)
  for atom, coord_atom in zip(soup.atoms(), coord_soup.atoms()):
    p = coord_atom.pos
    v3.set_vector(atom.pos, p[0], p[1], p[2])
  if in_vel:
    vel_soup = pdbatoms.Soup(in_vel)
    for atom, vel_atom in zip(soup.atoms(), vel_soup.atoms()):
      v = vel_atom.pos
      v3.set_vector(atom.vel, v[0], v[1], v[2])
  return soup


def write_soup_to_crds_and_vels(in_soup, basename):
  """
  From soup, writes out the coordinate/velocities, used for pulsing
  """
  soup = in_soup.copy()
  convert_to_namd_atom_names(soup)
  coor = basename + '.coor'
  soup.write_pdb(coor)
  for atom in soup.atoms():
    v3.set_vector(atom.pos, atom.vel[0], atom.vel[1], atom.vel[2])
  vel = basename + '.vel'
  soup.write_pdb(vel)
  return coor, vel


def convert_restart_to_pdb(basename, pdb):
  """
  Converts restart files with basename into PDB file.
  """
  psf, coor, vel = get_restart_files(basename)
  soup = soup_from_restart_files(psf, coor, vel)
  soup.write_pdb(pdb)
    


# ##########################################################

# # 2. Generate restart files from PDB

# The restart files used for PDBREMIX assumes a consistent file 
# naming. For a given basename `sim`, the files are:

# 1. topology file: sim.psf
# 2. coordinate file: sim.coor
# 3. optional velocity: sim.vel
# 4. optional box: sim.xsc

# To generate a topology file from the PDB file:

# - handles multiple protein chains
# - hydrogens are removed and then regenerated by GROMACS
# - disulfide bonds are identified and written into scripts
# - charged residue protonation states are auto-detected
# - explicit water in cubic box with 10.0 angstrom buffer
# - counterions to neutralize the system
# - use the CHARM22 topology/parameters included in this library

# Binaries used to generate restart files:

# 1. psfgen - create topologies for atoms in PDB
# 2. vmd - to add waters and counterions
 

module_load_script = """
package require psfgen 
topology %(topology)s
alias residue HIS HSE 
alias atom ILE CD1 CD 
"""

fixed_waters_script = """
pdbalias residue HOH TIP3
pdbalias residue WAT TIP3
segment wat { 
 auto none
 pdb %(water_pdb)s
} 
pdbalias atom HOH O OH2 
pdbalias atom WAT O OH2 
coordpdb %(water_pdb)s wat
"""

protein_chain_script = """
segment %(chain_id)s { 
 pdb %(chain_pdb)s 
} 
coordpdb %(chain_pdb)s %(chain_id)s
"""

write_script = """
guesscoord 
writepdb %(out_pdb)s
writepsf %(out_psf)s
"""


def make_chain_loading_script(pdb, basename):
  script = ""
  soup = pdbatoms.Soup(pdb)
  for chain_id in soup.chain_ids():
    if chain_id == ' ':
      chain_id = 'A'
    chain_pdb = '%s.chain.%s.pdb' % (basename, chain_id) 
    chain = soup.extract_chain(chain_id).write_pdb(chain_pdb)
    script += protein_chain_script % { 
      'chain_id': chain_id, 'chain_pdb': chain_pdb }
  return script


def make_disulfide_script(pdb):
  """
  Returns the psfgen script for disulfide bonds.

  This function opens in_pdb in a soup object, and searches for
  CYS residues where the SG-SG distance < 3 angs. These residues
  are then renamed to CYX and written to out_pdb. The disulfide bonds
  are then returned in a .tleap script fragment.
  """
  soup = pdbatoms.Soup(pdb)
  n = len(soup.residues())

  # First generate the residue names recognized by psfgen
  res_names = []
  chain_id = None
  i_res = None
  for i in range(n):
    res = soup.residue(i)
    if res.chain_id != chain_id:
      chain_id = res.chain_id
      i_res = 1
    res_names.append("%s:%s" % (chain_id, i_res))
    i_res += 1

  # Then search through for all CYS-CYS pairs and identify disulfide bonds
  script = ""
  for i in range(n):
    for j in range(i+1, n):
      if soup.residue(i).type in 'CYS' and soup.residue(j).type in 'CYS':
        sg1 = soup.residue(i).atom('SG')
        sg2 = soup.residue(j).atom('SG')
        if v3.distance(sg1.pos, sg2.pos) < 3.0:
          script += "patch DISU %s %s\n" % (res_names[i], res_names[j])
  if script:
     script = "# disulfide bonds\n" + script + "\n"

  return script


solvate_vmd_script = """
# Solvate system

# Set minimum padding
set pad %(solvent_buffer)s

# Run solvate with automatic padding option
package require solvate
resetpsf
solvate %(in_psf)s %(in_pdb)s -o %(name)s.vmd -rotate -rotinc 5 -t $pad

# Find the periodic box size.
mol load psf %(name)s.vmd.psf
set mol1 [molinfo top get id]
animate read pdb %(name)s.vmd.pdb $mol1

set sel1 [atomselect $mol1 all]
set cent1 [measure center $sel1]
set size1 [measure minmax $sel1]
mol delete $mol1

set sizemin [lindex $size1 0]
set sizemax [lindex $size1 1]

set xsize [expr [lindex $sizemax 0]-[lindex $sizemin 0]]
set ysize [expr [lindex $sizemax 1]-[lindex $sizemin 1]]
set zsize [expr [lindex $sizemax 2]-[lindex $sizemin 2]]

# Find the padding in each direction to make the box cubic
set lmax $xsize
if {$ysize > $lmax} {
	set lmax $ysize
}
if {$zsize > $lmax} {
	set lmax $zsize
}

# I like my boxe size to be a nice even number
set maxsize [expr int(ceil($lmax))]
if {$maxsize%%2 != 0} { set maxsize [expr $maxsize +1] }

# Calculate additional padding
set xpad [expr {$pad+0.5*($maxsize-$xsize)}]
set ypad [expr {$pad+0.5*($maxsize-$ysize)}]
set zpad [expr {$pad+0.5*($maxsize-$zsize)}]

puts ":Padding: $xpad $ypad $zpad"
puts ":Box size: $lmax"

# Adjust padding for nonzero center of mass. These are used to manually set the padding in each direction (both + and -)
set xplus [expr $xpad - [lindex $cent1 0]]
set xmin [expr $xpad + [lindex $cent1 0]]
set yplus [expr $ypad - [lindex $cent1 1]]
set ymin [expr $ypad + [lindex $cent1 1]]
set zplus [expr $zpad - [lindex $cent1 2]]
set zmin [expr $zpad + [lindex $cent1 2]]

# Rerun solvate on the original structure using calculated padding to make the box cubic
resetpsf
solvate %(name)s.psf %(name)s.pdb -o %(name)s.vmd -rotate -rotinc 5 -x $xmin +x $xplus -y $ymin +y $yplus -z $zmin +z $zplus

# Check that it worked
mol delete all
mol load psf %(name)s.vmd.psf
set mol1 [molinfo top get id]
animate read pdb %(name)s.vmd.pdb $mol1

set sel1 [atomselect $mol1 all]
set cent1 [measure center $sel1]
set size1 [measure minmax $sel1]

set sizemin [lindex $size1 0]
set sizemax [lindex $size1 1]

set xsize [expr [lindex $sizemax 0]-[lindex $sizemin 0]]
set ysize [expr [lindex $sizemax 1]-[lindex $sizemin 1]]
set zsize [expr [lindex $sizemax 2]-[lindex $sizemin 2]]

puts ":Final size: $xsize, $ysize, $zsize"
puts ":Length: $sizemax"
puts ":Center: $cent1"
puts ":Size: $size1"

# Ionize
package require autoionize
# Find original charge
set all [atomselect $mol1 all]
set q [measure sumweights $all weight charge]

# Determine the number and type of ions to use
set natom [expr round($q)]
if {$natom < 0} {
	set natom [expr abs($natom)]
	autoionize -psf %(name)s.vmd.psf -pdb %(name)s.vmd.pdb -nna $natom -ncl 0 -o %(name)s
} elseif {$natom > 0} {
	autoionize -psf %(name)s.vmd.psf -pdb %(name)s.vmd.pdb -nna 0 -ncl $natom -o %(name)s
} elseif {$natom == 0} {
	exec cp %(name)s.vmd.psf %(name)s.psf
	exec cp %(name)s.vmd.pdb %(name)s.pdb
}

# Check that it worked
mol delete all
mol load psf %(name)s.psf
set mol1 [molinfo top get id]
animate read pdb %(name)s.pdb $mol1
set all [atomselect $mol1 all]
puts ":Old charge: $q, New charge: [measure sumweights $all weight charge]"
mol delete all
"""

def solvate_psf(in_psf, in_pdb, basename, solvent_buffer=10.0):
  """
  Uses VMD to add explicit waters to a .psf topology file
  """
  parms = {
    'in_psf': in_psf,
    'in_pdb': in_pdb,
    'name': basename,
    'solvent_buffer': solvent_buffer,
  }
  tcl = basename + '.vmd.tcl'
  open(tcl, 'w').write(solvate_vmd_script % parms)
  data.binary('vmd', '-dispdev text -eofexit', basename+'.vmd', tcl)
  util.check_output(basename+'.vmd.pdb')
  util.check_output(basename+'.pdb')
 

def pdb_to_top_and_crds(force_field, pdb, basename, solvent_buffer=10.0): 
  """
  Creates CHARMM .coor and .psf file for NAMD simulation.
  """
  solv_dir = basename + '.solvate'
  save_dir = os.getcwd()

  pdb = os.path.abspath(pdb)

  util.goto_dir(solv_dir)

  # Remove all but protein heavy atoms in a single clean conformation
  stripped_pdb = basename + '.clean.pdb' 
  pdbtext.clean_pdb(pdb, stripped_pdb)

  # Make input script for psfgen
  psfgen_psf = basename+'.psfgen.psf'
  psfgen_pdb = basename+'.psfgen.pdb'
  script = module_load_script 
  script += make_chain_loading_script(stripped_pdb, basename)
  script += make_disulfide_script(stripped_pdb)
  script += write_script 
  script = script % {
    # load the included CHARMM2 atom topologies
    'topology': os.path.join(data.data_dir, 'charmm22.topology'),
    'out_pdb': psfgen_pdb,
    'out_psf': psfgen_psf
  }

  psfgen_in = basename+".psfgen.in"
  open(psfgen_in, "w").write(script)

  data.binary('psfgen', psfgen_in, basename+'.psfgen')
  util.check_output(psfgen_psf)
  util.check_output(psfgen_pdb)

  solvate_psf(psfgen_psf, psfgen_pdb, basename, solvent_buffer)

  psf = basename+'.psf'
  coor = basename+'.coor'
  pdb = basename+'.pdb'
  os.rename(pdb, coor)
  convert_restart_to_pdb(basename, pdb)
  shutil.copy(psf, save_dir)
  shutil.copy(coor, save_dir)
  shutil.copy(pdb, save_dir)
  
  os.chdir(save_dir)

  return psf, coor


# ##########################################################

# # 3. Run simulations from restart files

# Simulation approach for implicit solvent:
# - optional positional constraints: 100 kcal/mol/angs**2 
# - Langevin thermostat for constant temperature

# Simulation approach for explict water: 
# - periodic box with PME electrostatics
# - optional positional restraints: 100 kcal/mol/angs**2
# - Langevin thermostat for constant temperature
# - Nose-Hoover barometer with flexible periodic box size

# Binaries used:
# 1. namd2

# Files for trajectories:
# 1. coordinate trajectory: md.dcd
# 2. velocitiy trajectory: md.vel.dcd
# 3. restart coordinate: md.coor
# 4. restart velocity: md.vel
# 5. restart box: md.xsc

minimization_parms = { 
  'topology' : 'in.psf', 
  'input_crds' : 'in.pdb', 
  'output_basename' : 'min', 
  'force_field': 'NAMD', 
  'restraint_pdb': '',
  'restraint_force': 100.0,
  'n_step_minimization' : 100, 
} 

constant_energy_parms = { 
  'topology' : 'in.top', 
  'input_crds': 'in.coor',
  'input_vels': 'in.vel',
  'output_basename' : 'md', 
  'random_seed' : 2342, 
  'force_field': 'NAMD', 
  'n_step_per_snapshot' : 5, 
  'n_step_dynamics' : 1000, 
  'restraint_pdb': '',
  'restraint_force': 100.0,
} 

langevin_thermometer_parms = { 
  'topology' : 'in.top', 
  'input_crds': 'in.coor',
  'input_vels': '',
  'temperature_thermometer' : 300.0, 
  'force_field': 'NAMD', 
  'output_basename' : 'md', 
  'random_seed' : 2342, 
  'n_step_per_snapshot' : 5, 
  'n_step_dynamics' : 1000, 
  'restraint_pdb': '',
  'restraint_force': 100.0,
} 


io_script = """
# load topology and initial coordinates
structure %(topology)s
coordinates %(input_crds)s

# set output basenames
outputName %(output_basename)s
binaryOutput no
"""

simulation_parameters_script = """
# forcefield parameters
%(psf_type)s
parameters %(parameter)s

# switching parameters
exclude	scaled1-4
1-4scaling 1.0
dielectric 1.0
switching on
switchDist 8.0
cutoff 12.0
pairListDist 15.0
molly off

# Bond constraints using SHAKE
rigidBonds water

# Periodic Boundary Conditions and PME Electrostatics
wrapWater           on
wrapAll             on
wrapNearest         on
PME                 yes
PMEGridSpacing      1.0
"""

restraint_script="""
# restraints, called constraints here
constraints on
consExp 2
consKFile %(restraint_pdb)s
consKCol B
consRef %(restraint_pdb)s
constraintScaling %(restraint_force)s
"""

molecular_dynamics_script = """
# DCD output
dcdFile %(output_basename)s.dcd
dcdFreq %(n_step_per_snapshot)s
velDcdFile %(output_basename)s.vel.dcd
velDcdFreq %(n_step_per_snapshot)s

# Molecular Dynamics Timesteps
numSteps %(n_step_dynamics)s
timeStep 1  # in fs
firstTimeStep 0
stepsPerCycle 20
nonBondedFreq 1
fullElectFrequency 2
"""

extended_periodic_box_script = """    
# Periodic Box from file
extendedSystem %(xsc)s
"""

new_periodic_box_script = """
# Periodic Box Definition
cellBasisVector1    %(len_x)s     0.0           0.0
cellBasisVector2    0.0           %(len_y)s     0.0
cellBasisVector3    0.0           0.0           %(len_z)s
cellOrigin          %(x_origin)s  %(y_origin)s  %(z_origin)s
"""

def calculate_periodic_box_script(parms):
  """
  Returns namd2 input fragment to parameterize the periodic
  box for the protein. The requires loading the protein
  and directly calculating a good bounding box.
  """
  script = new_periodic_box_script
  p = pdbatoms.Soup(parms['input_crds'])
  atoms = p.atoms()
  parms = {}
  for i_axis, axis in enumerate(['x', 'y', 'z']):
    vals = [a.pos[i_axis] for a in atoms]
    axis_min, axis_max = min(vals), max(vals)
    parms["len_"+axis] = axis_max - axis_min + 0.5
    parms[axis+"_origin"] = sum(vals)/float(len(vals))
  return script % parms


import_velocities_script = """
# Import velocities
velocities %(input_vels)s
"""

generate_velocities_script = """
# Generate velocities from temperature
temperature %(temperature_initial_velocities)s
seed %(random_seed)s
"""

temperature_coupling_script = """
# Temperature Coupling with Langevin Thermometer
langevin on
langevinHydrogen off
langevinTemp %(temperature_thermometer)s
langevinDamping 5.0
"""

constant_pressure_script = """
# Constant Pressure Control (variable volume)
useGroupPressure      yes   # yes needed for rigidBonds
useFlexibleCell       no
useConstantArea       no

# Nose-Hoover Langevin barometer
langevinPiston        on
langevinPistonTarget  1.01325   #  in bar -> 1 atm
langevinPistonPeriod  200.0
langevinPistonDecay   100.0
langevinPistonTemp    %(temperature_thermometer)s
"""
  
minimization_script = """
# Minimization Parameters
minimization on
numsteps %(n_step_minimization)s
"""

def make_namd_input_file(parms):
  """
  Make namd2 input script to run simulations.
  """
  script = io_script % parms
  script += simulation_parameters_script % parms
  if parms['restraint_pdb']:
    script += restraint_script % parms
  if parms['xsc']:
    script += extended_periodic_box_script % parms
  else:
    script += calculate_periodic_box_script(parms)
  if 'n_step_dynamics' in parms:
    script += molecular_dynamics_script % parms
    if parms['input_vels']:
      script += import_velocities_script % parms
    elif 'temperature_initial_velocities' in parms and parms['temperature_initial_velocities'] > 0.0:
      script += generate_velocities_script % parms
    else:
      raise IndexError("No initial velocity information for dynamics run")
    if 'temperature_thermometer' in parms:
      script += temperature_coupling_script % parms
      script += constant_pressure_script % parms
  else:
    script += minimization_script % parms
  return script


def run(in_parms):
  """
  Runs a NAMD simulation using the PDBREMIX in_parms dictionary.
  """
  parms = copy.deepcopy(in_parms)
  name = parms['output_basename']

  # load the included CHARMM2 energy parameters
  parms['parameter'] = os.path.join(data.data_dir, 'charmm22.parameter')
  parms['psf_type'] =  'paraTypeCharmm on'

  # copy over input xsc and topology files (same basename)
  xsc = parms['topology'].replace('.psf', '.xsc')
  if os.path.isfile(xsc):
    shutil.copy(xsc, name + '.in.xsc')
    parms['xsc'] = name + '.in.xsc'
  else:
    parms['xsc'] = ''
  shutil.copy(parms['topology'], name + '.psf')
  parms['topology'] = name + '.psf'

  # copy over coordinates
  shutil.copy(parms['input_crds'], name + '.in.coor')
  parms['input_crds'] = name + '.in.coor'

  # copy over velocities
  if 'input_vels' in parms and parms['input_vels']:
    shutil.copy(parms['input_vels'], name + '.in.vel')
    parms['input_vels'] = name + '.in.vel'
  else:
    parms['input_vels'] = ''

  # copy over restraint coordinates
  if 'restraint_pdb' in parms and parms['restraint_pdb']:
    shutil.copy(parms['restraint_pdb'], name + '.restraint.coor')
    parms['restraint_pdb'] = name + '.restraint.coor'
  else:
    parms['restraint_pdb'] = ''
    
  namd_in = name + ".namd2.in"
  open(namd_in, "w").write(make_namd_input_file(parms))
  
  data.binary('namd2', namd_in, name + '.namd2')

  top, crds, vels = get_restart_files(name)
  util.check_output(top)
  util.check_output(crds)



# ##########################################################

# # 4. Read trajectories with some post-processing

# The units used in these files are:
# - positions: angstroms
# - velocities: angs/ps


def check_dcd_byte_order(dcd):
  """
  Uses flipdcd to DCD trajectories have matching OS endianess.
  """
  if sys.byteorder in 'big':
    option = '-B'
  elif sys.byteorder in 'little':
    option = '-L'
  else:
    raise Exception("Couldn't figure out system byte order %s" % sys.byteorder)
  data.binary('flipdcd', '%s %s' % (option, dcd), dcd+'.flipdcd')


class DcdReader:
  """
  Class to read CHARMM .dcd files.

  Attributes: 
    dcd (str) - name of .dcd file
    n_frame (int) - number of frames in trajectory
    remarks (list) - title of file 
    pos_after_header (int) - position of frames in file
    size_frame (int) - the size of the frame in bytes
    n_atom (int) - number of atoms simulated
    n_frame (int) - number of frames in trajectory
    i_frame (int) - index of current frame
    frame (int) - container of coordinates of current frame

  Methods:
    __init__ - open dcd, read header and load first frame
    load_frame(i) - loads the coordinates of the i'th frame
    __getitme__ - return a frame as a list (3xn) of floats.
  """

  def __init__(self, dcd):
    self.dcd = dcd

    check_dcd_byte_order(self.dcd)    
    
    self.file = open(dcd, 'rb')

    # Read heading block
    if self._read_fmt_val('i') != 84:
      raise Exception("DCD: first int of heading is not 84")
    if self._read_fmt_vals('4c') != ('C', 'O', 'R', 'D') :
     raise Exception("DCD: Missing 'CORD' tag")
    self.n_frame = self._read_fmt_val('i') 
    self.i_start = self._read_fmt_val('i')
    self.n_step_save = self._read_fmt_val('i')
    self._read_fmt_val('5i')
    self.n_fixed_atom = self._read_fmt_val('i')
    if self.n_fixed_atom > 0:
      raise Exception("DCD: this reader doesn't handle fixed atoms")
    self.timeStep = self._read_fmt_val('d')
    self._read_fmt_vals('9i')
    if self._read_fmt_val('i') != 84 :
      raise Exception("DCD: couldn't find ending 84 of the heading")

    # Read title block
    size = self._read_fmt_val('i')
    if (size - 4) % 80 != 0 :
      raise Exception("DCD: title block size is wrong")
    self.remarks = []
    n_line = self._read_fmt_val('i')
    for i in range(0, n_line):
      s = "".join(self._read_fmt_vals('80c'))
      self.remarks.append(s.strip())
    if self._read_fmt_val('i') != size:
      raise Exception("DCD: title block end != block start")

    # Read n_atom block
    if self._read_fmt_val('i') != 4 :
      raise Exception("DCD: n_atom field start != 4")
    self.n_atom = self._read_fmt_val('i')
    if self._read_fmt_val('i') != 4 :
      raise Exception("DCD: n_atom field end != 4")

    # Store end of header position
    self.pos_after_header = self.file.tell()

    # Calculate size_frame from the size of the rest of file
    n_free_atom = self.n_atom - self.n_fixed_atom
    self.size_frame = struct.calcsize('%df6i' % (3 * n_free_atom))
    self.file.seek(0, 2)
    last_pos = self.file.tell()
    size_rest_of_file = last_pos - self.pos_after_header
    implied_size_frame = size_rest_of_file / self.n_frame
    self.extra_block_size = implied_size_frame - self.size_frame 
    if self.extra_block_size > 0:
      self.size_frame += self.extra_block_size

    self.i_frame = None
    self.load_frame(0)

  def _read_fmt_vals(self, fmt):
    return struct.unpack(fmt, self.file.read(struct.calcsize(fmt)))

  def _read_fmt_val(self, fmt):
    return self._read_fmt_vals(fmt)[0]
    
  def load_frame(self, i):
    """
    Reads self.frame from loaded trajectory file. The frame
    is (x_vals, y_vals, z_vals) for each atom.
    """
    if i < - 1*self.n_frame or i >= self.n_frame:
      raise IndexError
    if i < 0:
      i = self.n_frame + i
    if i == self.i_frame:
      return

    frame_coords_pos = self.pos_after_header + i*self.size_frame
    frame_coords_pos += self.extra_block_size
    self.file.seek(frame_coords_pos)

    coords_size_fmt = "%df" % self.n_atom
    size = struct.calcsize(coords_size_fmt)
    if size != self._read_fmt_val('i'):
      raise Exception("DCD: x_vals block start != size")
    x_vals = self._read_fmt_vals(coords_size_fmt)
    if size != self._read_fmt_val('i'):
      raise Exception("DCD: x_vals block end != size")
    if size != self._read_fmt_val('i'):
      raise Exception("DCD: y_vals block start != size")
    y_vals = self._read_fmt_vals(coords_size_fmt)
    if size != self._read_fmt_val('i'):
      raise Exception("DCD: y_vals block end != size")
    if size != self._read_fmt_val('i'):
      raise Exception("DCD: z_vals block start != size")
    z_vals = self._read_fmt_vals(coords_size_fmt)
    if size != self._read_fmt_val('i'):
      raise Exception("DCD: z_vals block end != size")
    self.frame = (x_vals, y_vals, z_vals)
    self.i_frame = i

  def __getitem__(self, i):
    self.load_frame(i)
    return self.frame

  def __repr__(self):
    return "< DCD %s with %d frames of %d atoms >" % \
             (self.dcd, self.n_frame, self.n_atom)


class SoupTrajectory:
  """
  Class to interact with an CHARMM/NAMD DCD trajctory using soup.
  
  Attributes:
    basename (str) - basename used to guess all required files
    psf (str) - topology file of trajectory
    dcd (str) - coordinate trajectory file
    vel_dcd (str) - velocity trajectory file
    coor_dcd_reader (DcdReader) - the reader of the coordinates
    vel_dcd_reader (DcdReader) - the reader of the velocitiies
    n_frame (int) - number of frames in trajectory
    i_frame (int) - index of current frame
    soup (Soup) - Soup object holding current coordinates/velocities

  Methods:
    __init__ - load coordinate and velocity trajectories and build soup
    load_frame - loads new frame into soup
  """  

  def __init__(self, soup, dcd, vel_dcd=''):
    self.soup = soup
    self.dcd = dcd
    self.coor_dcd_reader = DcdReader(self.dcd)
    self.vel_dcd = vel_dcd
    if vel_dcd:
      self.vel_dcd_reader = DcdReader(self.vel_dcd)
    else:
      self.vel_dcd_reader = None
    self.n_frame = self.coor_dcd_reader.n_frame
    self.i_frame = 0
    self.load_frame(0)
    
  def load_frame(self, i_frame):
    x, y, z = self.coor_dcd_reader[i_frame]
    atoms = self.soup.atoms()
    for i in range(len(atoms)):
      v3.set_vector(atoms[i].pos, x[i], y[i], z[i])

    if self.vel_dcd_reader is not None:
      x, y, z = self.vel_dcd_reader[i_frame]
      for i in range(len(atoms)):
          v3.set_vector(atoms[i].vel, x[i], y[i], z[i])

    self.i_frame = self.coor_dcd_reader.i_frame


class Trajectory:
  """
  Class to interact with an CHARMM/NAMD DCD trajctory using soup.
  
  Attributes:
    basename (str) - basename used to guess all required files
    psf (str) - topology file of trajectory
    dcd (str) - coordinate trajectory file
    vel_dcd (str) - velocity trajectory file
    coor_dcd_reader (DcdReader) - the reader of the coordinates
    vel_dcd_reader (DcdReader) - the reader of the velocitiies
    n_frame (int) - number of frames in trajectory
    i_frame (int) - index of current frame
    soup (Soup) - Soup object holding current coordinates/velocities

  Methods:
    __init__ - load coordinate and velocity trajectories and build soup
    load_frame - loads new frame into soup
  """  

  def __init__(self, basename):
    self.basename = basename
    self.psf = basename + '.psf'
    self.soup = soup_from_psf(self.psf)
    self.dcd = basename + '.dcd'    
    self.vel_dcd = basename + '.vel.dcd'    
    if not os.path.isfile(self.vel_dcd):
      self.vel_dcd = ''
    self.soup_trj = SoupTrajectory(self.soup, self.dcd, self.vel_dcd)
    self.n_frame = self.soup_trj.n_frame
    self.i_frame = 0
    self.load_frame(0)
    
  def load_frame(self, i_frame):
    self.soup_trj.load_frame(i_frame)
    self.i_frame = self.soup_trj.i_frame


def merge_dcds(psf, dcds, out_dcd):
  """
  Given a list of traj filenames (trajs), merges them into one complete
  trajectory (out_traj) using top to work out the number of atoms, and
  hence the size of the frame of the trajectory.
  """
  dcd_reader = DcdReader(dcds[0])
  pos_after_header = dcd_reader.pos_after_header  
  size_frame = dcd_reader.size_frame
  del dcd_reader

  shutil.copy(dcds[0], out_dcd)

  merge_dcd_file = open(out_dcd, "ab+")
  for dcd in dcds[1:]:
    dcd_file = open(dcd, "rb")
    dcd_file.seek(-1, 2)
    eof = dcd_file.tell()

    dcd_file.seek(pos_after_header)
    while dcd_file.tell() < eof:
      merge_dcd_file.write(dcd_file.read(size_frame)) 

    dcd_file.close()
  merge_dcd_file.close()


def merge_trajectories(basename, traj_basenames):
  """
  Splices together a bunch of simulations, all with the same
  basename, into one large simulation in the current directory.
  """
  for ext in ['.psf', '.coor', '.vel', '.xsc']:
    f = traj_basenames[-1] + ext
    g = basename + ext
    shutil.copy(f, g)
  trajs = [b + '.dcd' for b in traj_basenames]
  merge_dcds(basename + '.psf', trajs, basename + '.dcd')
  vels = [b + '.vel.dcd' for b in traj_basenames]
  merge_dcds(basename + '.psf', vels, basename + '.vel.dcd')


# def merge_simulations(basename, pulses):
#   """
#   Splices together a bunch of simulations, all with the same
#   basename, into one large simulation in the current directory.
#   """
#   for ext in ['.psf', '.coor', '.vel', '.xsc']:
#     fname = '%s%s' % (basename, ext)
#     shutil.copy('%s/%s' % (pulses[-1], fname), fname)
#   trajs = [os.path.join(pulse, basename + '.dcd') for pulse in pulses]
#   merge_trajectories(basename + '.psf', trajs, basename + '.dcd')
#   vels = [os.path.join(pulse, basename + '.vel.dcd') for pulse in pulses]
#   merge_trajectories(basename + '.psf', vels, basename + '.vel.dcd')
  



