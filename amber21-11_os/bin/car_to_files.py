#!/Users/runner/miniforge3/conda-bld/ambertools_1635195607525/_h_env_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehol/bin/python
# Filename: car_to_files.py

from pymsmt.mol.itfc_dict_table import atom_type_dict3 as attyp_dict
from pymsmt.mol.cario import read_carf, print_mol2f, print_pdbf
from optparse import OptionParser

parser = OptionParser("Usage: car_to_files.py -i input_file -m mol2_file -p pdb_file -r residue_name \n")

parser.add_option("-i", dest="input_file", type='string',
                  help="Input file name")
parser.add_option("-m", dest="mol2_file", type='string',
                  help="Output mol2 file name")
parser.add_option("-p", dest="pdb_file", type='string',
                  help="Output PDB file name")
parser.add_option("-r", dest="resname", type='string',
                  help="Residue name")

(options, args) = parser.parse_args()

mol, atids, resids, pbc_size = read_carf(options.input_file)

print_mol2f(mol, atids, options.mol2_file, options.resname, attyp_dict)

for i in atids:
    mol.atoms[i].resname = options.resname

print_pdbf(mol, atids, options.pdb_file, pbc_size)

