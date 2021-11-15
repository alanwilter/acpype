# encoding: utf-8

__doc__ = """ 

Stores data for PDBREMIX.

- directory of explici data directory
- mappings of residue names to chars
- names and locations of binaries used
- backbone atom types
- solvent residue types
- element radii
- chi dihedral angle topologies
"""


import os
from . import util


module_dir = os.path.dirname(__file__)
data_dir = os.path.join(module_dir, 'data')

binaries = {
  "pymol": "",
  "pymol_batch": "",
  "chimera": "",

  "theseus": "",
  "mafft": "",

  "sander": "",
  "tleap": "",

  "mdrun": "",
  "pdb2gmx": "",
  "trjconv": "",
  "grompp": "",
  "editconf": "",
  "genion": "",
  "genbox": "",

  "vmd": "",
  "psfgen": "",
  "namd2": "",
  "flipdcd": "",

  "mod9v8": ""
}
home_dir = os.path.expanduser('~')
binaries_fname = os.path.join(home_dir, '.pdbremix.config')
backup_binaries_fname = os.path.join(os.environ['AMBERHOME'], 'dat') if 'AMBERHOME' in os.environ else None
if not os.path.isfile(binaries_fname):
  try:
    util.write_dict(binaries_fname, binaries)
  except IOError:
    # May not have permission here, especially in a docker container. If so,
    # fall back on AMBERHOME
    if backup_binaries_fname is None:
      raise
    binaries_fname = backup_binaries_fname
else:
  binaries = util.read_dict(binaries_fname)


def binary(bin, arg_str='', out_name=None, in_fname=None):
  """
  Runs an external binary, handles arguments, writes out
  equivalent .sh file, log file, and can pipe in in_fname.
  """
  if bin in binaries and binaries[bin]:
    bin = binaries[bin]
  else:
    util.check_program(bin)
  if arg_str:
    util.run_with_output_file(
        '%s %s' % (bin, arg_str), out_name, in_fname)
  return '"%s"' % bin


def invert_dict(d):
  """
  Returns dictionaries with swapped key-values.
  """
  return dict((v, k) for k,v in list(d.items()))


res_name_to_char = {
    "ALA":"A", "CYS":"C", "ASP":"D",
    "GLU":"E", "PHE":"F", "GLY":"G",
    "HIS":"H", "ILE":"I", "LYS":"K",
    "LEU":"L", "MET":"M", "ASN":"N",
    "PRO":"P", "GLN":"Q", "ARG":"R",
    "SER":"S", "THR":"T", "VAL":"V",
    "TRP":"W", "TYR":"Y", "ACE":">",
    "NME":"<",
}

res_char_to_name = invert_dict(res_name_to_char)

res_name_to_char.update({"HID":"H", "HIE":"H", "HIP":"H", "HSE":"H", "LYP":"K", "CYM":"C", "CYX":"C", "CYN":"C"})

# recognized atom types for protein backbone 
backbone_atoms = [
    "OXT", # C-terminal carboxyl group
    "H1", "H2", "H3",  # N-terminal charged group
    "C", "O", # peptide-bond carbonyl group
    "H", "HN", "N",  # peptide-bond amide group
    "CA", "HA" # main-chain C-alpha and alkyl group
    ]
    
solvent_res_types = [
    'HOH', 'WAT', 'TIP', 'SOL',
    'CLA', 'SOD', 'NA', 'CL', 
    'NA+', 'CL-', 'Na', 'Cl',
    'Na+', 'Cl-']

radii = { 
 'H':  1.20,
 'N':  1.55,
 'NA': 2.27,
 'CU': 1.40,
 'CL': 1.75,
 'C':  1.70,
 'O':  1.52,
 'I':  1.98,
 'P':  1.80,
 'B':  1.85,
 'BR': 1.85,
 'S':  1.80,
 'SE': 1.90,
 'F':  1.47,
 'FE': 1.80,
 'K':  2.75,
 'MN': 1.73,
 'MG': 1.73,
 'ZN': 1.39,
 'HG': 1.80,
 'XE': 1.80,
 'AU': 1.80,
 'LI': 1.80,
 '.':  1.80
}

masses = {'H'  : 1.008,   'HE' : 4.003,   'LI' : 6.941,   'BE' : 9.012,
          'B'  : 10.811,  'C'  : 12.011,  'N'  : 14.007,  'O'  : 15.999,
          'F'  : 18.998,  'NE' : 20.180,  'NA' : 22.990,  'MG' : 24.305,
          'AL' : 26.982,  'SI' : 28.086,  'P'  : 30.974,  'S'  : 32.066,
          'CL' : 35.453,  'AR' : 39.948,  'K'  : 39.098,  'CA' : 40.078,
          'SC' : 44.956,  'TI' : 47.867,  'V'  : 50.942,  'CR' : 51.996,
          'MN' : 54.938,  'FE' : 55.845,  'CO' : 58.933,  'NI' : 58.693,
          'CU' : 63.546,  'ZN' : 65.38,   'GA' : 69.723,  'GE' : 72.631,
          'AS' : 74.922,  'SE' : 78.971,  'BR' : 79.904,  'KR' : 84.798,
          'RB' : 84.468,  'SR' : 87.62,   'Y'  : 88.906,  'ZR' : 91.224,
          'NB' : 92.906,  'MO' : 95.95,   'TC' : 98.907,  'RU' : 101.07,
          'RH' : 102.906, 'PD' : 106.42,  'AG' : 107.868, 'CD' : 112.414,
          'IN' : 114.818, 'SN' : 118.711, 'SB' : 121.760, 'TE' : 126.7,
          'I'  : 126.904, 'XE' : 131.294, 'CS' : 132.905, 'BA' : 137.328,
          'LA' : 138.905, 'CE' : 140.116, 'PR' : 140.908, 'ND' : 144.243,
          'PM' : 144.913, 'SM' : 150.36,  'EU' : 151.964, 'GD' : 157.25,
          'TB' : 158.925, 'DY' : 162.500, 'HO' : 164.930, 'ER' : 167.259,
          'TM' : 168.934, 'YB' : 173.055, 'LU' : 174.967, 'HF' : 178.49,
          'TA' : 180.948, 'W'  : 183.84,  'RE' : 186.207, 'OS' : 190.23,
          'IR' : 192.217, 'PT' : 195.085, 'AU' : 196.967, 'HG' : 200.592,
          'TL' : 204.383, 'PB' : 207.2,   'BI' : 208.980, 'PO' : 208.982,
          'AT' : 209.987, 'RN' : 222.081, 'FR' : 223.020, 'RA' : 226.025,
          'AC' : 227.028, 'TH' : 232.038, 'PA' : 231.036, 'U'  : 238.029,
          'NP' : 237,     'PU' : 244,     'AM' : 243,     'CM' : 247, 'BK' : 247,
          'CT' : 251,     'ES' : 252,     'FM' : 257,     'MD' : 258, 'NO' : 259,
          'LR' : 262,     'RF' : 261,     'DB' : 262,     'SG' : 266, 'BH' : 264,
          'HS' : 269,     'MT' : 268,     'DS' : 271,     'RG' : 272, 'CN' : 285,
          'NH' : 284,     'FL' : 289,     'MC' : 288,     'LV' : 292, 'TS' : 294,
          'OG' : 294}

two_char_elements = [e for e in list(radii.keys()) if len(e) == 2]


def strip_numbers(s):
  result = ""
  for c in s:
    if not c.isdigit() and c != " ":
      result += c
  return result


def guess_element(res_type, atom_type):
  """
  Returns the element type using a dirty heuristic guess.
  """
  atom_type = strip_numbers(atom_type)
  if len(atom_type) == 0:
    raise Exception("PDB line with no atom name found")
  if res_type in res_name_to_char:
    return atom_type[0]
  if len(atom_type) == 2 and atom_type in two_char_elements:
    return atom_type
  return atom_type[0]  


chi_topology = {
  'ARG': [ ['N', 'CA', 'CB', 'CG'],
           ['CA', 'CB', 'CG', 'CD'],
           ['CB', 'CG', 'CD', 'NE'],
           ['CG', 'CD', 'NE', 'CZ']],
  'ASN': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'OD1']],
  'ASP': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'OD1']],
  'CYS': [['N', 'CA', 'CB', 'SG']],
  'GLN': [ ['N', 'CA', 'CB', 'CG'],
           ['CA', 'CB', 'CG', 'CD'],
           ['CB', 'CG', 'CD', 'OE1']],
  'GLU': [ ['N', 'CA', 'CB', 'CG'],
           ['CA', 'CB', 'CG', 'CD'],
           ['CB', 'CG', 'CD', 'OE1']],
  'HIS': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'ND1']],
  'ILE': [['N', 'CA', 'CB', 'CG1'], ['CA', 'CB', 'CG1', 'CD1']],
  'LEU': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'CD1']],
  'LYN': [ ['N', 'CA', 'CB', 'CG'],
           ['CA', 'CB', 'CG', 'CD'],
           ['CB', 'CG', 'CD', 'CE'],
           ['CG', 'CD', 'CE', 'NZ']],
  'LYP': [ ['N', 'CA', 'CB', 'CG'],
           ['CA', 'CB', 'CG', 'CD'],
           ['CB', 'CG', 'CD', 'CE'],
           ['CG', 'CD', 'CE', 'NZ']],
  'LYS': [ ['N', 'CA', 'CB', 'CG'],
           ['CA', 'CB', 'CG', 'CD'],
           ['CB', 'CG', 'CD', 'CE'],
           ['CG', 'CD', 'CE', 'NZ']],
  'MET': [ ['N', 'CA', 'CB', 'CG'],
           ['CA', 'CB', 'CG', 'SD'],
           ['CB', 'CG', 'SD', 'CE']],
  'PHD': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'OD1']],
  'PHE': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'CD1']],
  'PRO': [ ['N', 'CA', 'CB', 'CG'],
           ['CA', 'CB', 'CG', 'CD'],
           ['CB', 'CG', 'CD', 'N'],
           ['CG', 'CD', 'N', 'CA']],
  'SER': [['N', 'CA', 'CB', 'OG']],
  'THR': [['N', 'CA', 'CB', 'OG1']],
  'TRP': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'CD1']],
  'TYR': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'CD1']],
  'VAL': [['N', 'CA', 'CB', 'CG1']]}


def get_res_chi_topology(res_type):
  """
  Returns the chi topology for a given residue, which is a list of
  atoms that are affected if one rotates the chi0, chi1... 
  dihedral angle.
  """
  # Some common residue renamings in AMBER and GROMACS
  if res_type in ["HID", "HIE", "HIP", "HSE"] or "HIS" in res_type:
    res_type = "HIS"
  if res_type in ["LYP"]:
    res_type = "LYS"
  if res_type in ["CYM", "CYX", "CYN"]:
    res_type = "CYS"
  if res_type not in chi_topology:
    return []
  else:
    return chi_topology[res_type]



