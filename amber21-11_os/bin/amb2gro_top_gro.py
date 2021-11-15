#!/Users/runner/miniforge3/conda-bld/ambertools_1635195607525/_h_env_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehol/bin/python
# Filename: amb2gro_top_gro.py
#
# This is a program written by Pengfei Li to convert the AMBER prmtop and
# inpcrd files to the GROMACS top and gro files
#
import parmed as pmd
from optparse import OptionParser
import os

parser = OptionParser("Usage: amb2gro_top_gro.py -p prmtop -c inpcrd -t top\n"
                      "                          -g gro -b pdb")
parser.add_option("-p", dest="prmtop", type='string',
                  help="Prmtop file")
parser.add_option("-c", dest="inpcrd", type='string',
                  help="Inpcrd file")
parser.add_option("-t", dest="top", type='string',
                  help="GROMACS top file")
parser.add_option("-g", dest="gro", type='string',
                  help="GROMACS gro file")
parser.add_option("-b", dest="pdb", type='string',
                  help="A PDB file to generate")
(options, args) = parser.parse_args()

# Load the AMBER prmtop and inpcrd files
amber = pmd.load_file(options.prmtop, xyz=options.inpcrd)

# Save GROMACS top and gro files, along with a PDB file
amber.save(options.top, overwrite=True, format='gromacs')
amber.save(options.gro, overwrite=True, format='gro')
amber.save(options.pdb, overwrite=True, format='pdb')

quit()

