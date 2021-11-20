#!/Users/runner/miniforge3/conda-bld/ambertools_1635195607525/_h_env_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehol/bin/python
# Filename: espgen.py
"""
This is the espgen.py program written by Pengfei Li in Merz Research
Group at Michigan State University.
It is a program for generating the esp points from the Gaussian/GAMESS-US
output file.
"""
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
from pymsmt.mol.gauio import get_esp_from_gau
from pymsmt.mol.gmsio import get_esp_from_gms
from optparse import OptionParser
from pymsmt.title import print_title

parser = OptionParser("Usage: espgen.py -i input_file -o output_file "
                              "[-v software]")

parser.set_defaults(softversion='gau')

parser.add_option("-i", dest="inputfile", type='string',
                  help="Input file name")
parser.add_option("-o", dest="outputfile", type='string',
                  help="Output file name")
parser.add_option("-v", dest="softversion", type='string',
                  help="Software version [Default is gau (means Gaussian), \n"
                       "           other option is gms (means GAMESS-US)]")
(options, args) = parser.parse_args()

# Print the title of the program
version = '1.0'
print_title('espgen.py', version)

if options.softversion == 'gau':
    get_esp_from_gau(options.inputfile, options.outputfile)
elif options.softversion == 'gms':
    get_esp_from_gms(options.inputfile, options.outputfile)

quit()
