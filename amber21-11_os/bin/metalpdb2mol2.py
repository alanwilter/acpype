#!/Users/runner/miniforge3/conda-bld/ambertools_1635195607525/_h_env_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehol/bin/python
# Filename: metalpdb2mol2.py
"""
A script to convert metal ion PDB files to mol2 files, specifically written
for MCPB.py users.
The PDB file for conversion should have a single metal ion and only.
"""

from optparse import OptionParser
from pymsmt.mol.pdbio import get_atominfo_fpdb

parser = OptionParser("Usage: metalpdb2mol2.py -i pdb_file -o mol2_file -c charge")
parser.add_option("-i", dest="input_file", type='string',
                  help="Input PDB file")
parser.add_option("-o", dest="output_file", type='string',
                  help="Output mol2 file")
parser.add_option("-c", dest="charge", type='float',
                  help="Charge of the metal ion")
(options, args) = parser.parse_args()

def write_mol2file():

    mol2f = open(output_file, 'w')

    print('***Generating the ' + output_file + ' file...')

    ##1. molecule information
    print("@<TRIPOS>MOLECULE", file=mol2f)
    print(resname, file=mol2f)
    print('%5d%6d%6d%6d%6d' %(1, 0, 1, 0, 0), file=mol2f) #atom number and bond number
    print('SMALL', file=mol2f)
    print('User Assigned Charge', file=mol2f)
    print(' ', file=mol2f)
    print(' ', file=mol2f)

    ##2. atom information
    print('@<TRIPOS>ATOM', file=mol2f)
    for atm in atids:
        #new atom id
        atid = mol.atoms[atm].atid
        resid = mol.atoms[atm].resid
        atname = mol.atoms[atm].atname
        crd = mol.atoms[atm].crd
        atomtype = mol.atoms[atm].atname
        print('%7d %-4s    %10.4f%10.4f%10.4f %-4s %6d %-4s %12.6f'\
                        %(atm, atname, crd[0], crd[1], crd[2], atomtype,
                         resid, resname, chg), file=mol2f)

    ##3. bond information
    print('@<TRIPOS>BOND', file=mol2f)

    ##4. substructure information
    print('@<TRIPOS>SUBSTRUCTURE', file=mol2f)
    print('     1', resname, '         1 TEMP' + \
                    '              0 ****  ****    0 ROOT', file=mol2f)
    mol2f.close()

mol, atids, resids = get_atominfo_fpdb(options.input_file)
output_file = options.output_file
chg = options.charge
resname = mol.atoms[atids[0]].resname

write_mol2file()

quit()

