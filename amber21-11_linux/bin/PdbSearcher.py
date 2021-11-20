#!/home/conda/feedstock_root/build_artifacts/ambertools_1635195478443/_h_env_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_/bin/python
# Filename: PdbSearcher.py
"""
This is the PdbSeacher.py program written by Pengfei Li in Merz Research Group
at Michigan State University.
It is a re-written python version of Pdbseacher in MTK++.
It is designed to find the metal center in the PDB files and collecting the
information and generate the metal center complex pdb files for each center
(with metal ion and ligating residues).
"""

from pymsmt.mol.pdbio import get_atominfo_fpdb, writepdbatm
from pymsmt.mol.element import METAL_PDB, CoRadiiDict, resdict
from pymsmt.mol.mol import pdbatm
from pymsmt.mol.cal import calc_bond, det_geo
from pymsmt.title import print_title
from pymsmt.exp import *
from optparse import OptionParser
import warnings
import os

#==============================================================================
# Setting the options
#==============================================================================

parser = OptionParser("Usage: PdbSearcher.py -i/--ion ionname "
                      "-l/--list input_file \n"
                      "                      -e/--env environment_file "
                      "-s/--sum summary_file \n"
                      "                      [-c/--cut cutoff]")
parser.add_option("-i", "--ion", type='string', dest="ionname",
                  help="Element symbol of ion, e.g. Zn")
parser.add_option("-l", "--list", type='string', dest="inputf",
                  help="List file name, list file contains one PDB file name "
                       "per line")
parser.add_option("-e", "--env", type='string', dest='envrmtf',
                  help="Environment file name. An environment file is used to "
                       "store the metal center environment information such "
                       "as ligating atoms, distance, geometry etc. For each "
                       "bond, there is a record.")
parser.add_option("-s", "--sum", type='string', dest='sumf',
                  help="Summary file name. A summary file is used to store "
                       "the metal center summary information such as metal "
                       "center geometry, ligating residues etc. For each "
                       "metal center there is a record.")
parser.add_option("-c", "--cut", type='float', dest='cutoff',
                  help="Optional. The cut off value used to detect the bond "
                       "between metal ion and ligating atoms. The unit is "
                       "Angstroms. If there is no value specified, the "
                       "default algorithm will be used. The default algorithm "
                       "recognizes the bond when its distance is no less than "
                       "0.1 (smaller than 0.1 usually indicates a low quality "
                       "structure) and no bigger than the covalent radius sum "
                       "of the two atoms with a tolerance of 0.4.")
(options, args) = parser.parse_args()

#==============================================================================
# Print the title
#==============================================================================
version = '1.0'
print_title('PdbSearcher.py', version)

#==============================================================================
# Read in th pdb file names from input file
#==============================================================================
#pdb file name list
pdbfnl = []

#read the input file list
fp = open(options.inputf, 'r')
for line in fp:
    line = line.strip('\n').strip()
    pdbfnl.append(line)
fp.close()

#==============================================================================
# Print the title of each file
#==============================================================================

ionname = options.ionname

#Transfer the metal ion name
if len(ionname) == 2:
    ionname = ionname[0] + ionname[1:].lower()

print("The ionname you chosen is : " + ionname)

if options.cutoff != None:
    print("The cutoff is: " + str(options.cutoff) + ' Angstrom.')
else:
    print("Using the default method to determine the bond exists.")

#summary file
sf = open(options.sumf, 'w')
print('PDB_ID,', 'EXP_TECH,', 'RESOLUTION,', 'ATOM_NUMBER,', \
             'ION_NUMBER,', 'RES_ID,', 'RES_NAME,', 'ATOM_ID,', 'ATOM_NAME,', \
             'COORD_SPHERE,', 'GEOMETRY,','GEO_RMS', file=sf)

#environment file
ef = open(options.envrmtf, 'w')
print('PDB,', 'ION_RESID,', 'ION_RESNAME,', 'ION_ATOM_ID,', \
             'ION_ATOM_NAME,', 'RESID,', 'RESNAME,', 'ATOM_ID,', 'ATOM_NAME,',\
             'DISTANCE,', 'GEOMETRY,', 'GEO_RMS,', 'COORDINATE_SPHERE,', \
             'EXP_TECH,', 'RESOLUTION', file=ef)

#==============================================================================
# Do analysis for each pdb file
#==============================================================================

for fname in pdbfnl:
    print("***Performing the " + fname + " file")

    #get the metal list
    mol, atids, resids = get_atominfo_fpdb(fname)

    #Get the resolution and method
    fp1 = open(fname, 'r')
    for line in fp1:
        if 'RESOLUTION.' in line:
            line = line.split()
            try:
                reso = float(line[-1])
            except:
                try:
                    reso = float(line[-2])
                except:
                    reso = 'UNKNOWN'
        elif 'EXPERIMENT TYPE' in line:
            line = line.split()
            exptyp = line[-1]
            if line[-1] == 'DIFFRACTION' and line[-2] == 'X-RAY':
                exptyp = 'X-RAY'
    fp1.close()

    #Get the metal ion which is the ion user want to process
    metallist = []
    for i in atids:
        resname = mol.residues[mol.atoms[i].resid].resname
        atname = mol.atoms[i].atname
        if (resname, atname) in list(METAL_PDB.keys()):
            if METAL_PDB[(resname, atname)][0] == ionname:
                metallist.append(i)

    #for each metal ion in the metal list, print the metal center
    for i in metallist:

        mccrds = [] #The crds of metal site
        crdi = mol.atoms[i].crd
        elmti = mol.atoms[i].element
        residi = mol.atoms[i].resid
        atnamei = mol.atoms[i].atname
        resnamei = mol.residues[residi].resname
        radiusi = CoRadiiDict[elmti]
        mcresids = [] #MetalCenter residue IDs

        #Get the residues which is the metal site
        for j in atids:
            if j != i:
                atnamej = mol.atoms[j].atname
                crdj = mol.atoms[j].crd
                residj = mol.atoms[j].resid
                resnamej = mol.residues[residj].resname
                elmtj = mol.atoms[j].element
                radiusj = CoRadiiDict[elmtj]
                radiusij = radiusi + radiusj + 0.40
                disij = calc_bond(crdi, crdj)

                addon = 0
                if options.cutoff == None:
                    if (disij >= 0.1) and (disij <= radiusij) \
                       and (elmtj != 'H'):
                           addon = 1
                else:
                    if (disij >= 0.1) and (disij <= options.cutoff) \
                       and (elmtj != 'H'):
                           addon = 1

                if addon == 1:
                    #Warning of ligating atoms
                    if elmtj not in ['O', 'N', 'S', 'F', 'Cl', 'Br', 'I']:
                        if options.cutoff == None:
                            warnings.warn('Element %s was found ligating to %s '
                                      'with distance %5.3f, may need to '
                                      'specify the cut off value.'
                                      %(elmtj, elmti, disij), pymsmtWarning)
                        else:
                            warnings.warn('Element %s was found ligating to %s '
                                      'with distance %5.3f, the cut off value '
                                      '%5.3f may need to change.'
                              %(elmtj, elmti, disij, options.cutoff), pymsmtWarning)
                    mccrds.append(crdi)
                    mccrds.append(crdj)
                    if (residj not in mcresids):
                        mcresids.append(residj)

        #Getting the ligating reidue letters
        reslets = ''
        for j in mcresids:
            resname = mol.residues[j].resname
            if resname in list(resdict.keys()):
                reslet = resdict[resname]
            else:
                reslet = 'X'
            reslets = reslets + reslet
        nospace = ''
        reslets = nospace.join(sorted(reslets))
        print('   Find metal center', reslets)

        #Get the geometry and geometry rms
        try:
            geo, georms = det_geo(mccrds)
        except:
            if options.cutoff == None:
                warnings.warn('No atoms were found coordinated to the metal! '
                          'Suggest to specify explicit cut off value.'
                          , pymsmtWarning)
            else:
                warnings.warn('No atoms were found coordinated to the metal! '
                          'The cut off value %5.3f may need to change.'
                          %options.cutoff, pymsmtWarning)

        #add the metal ions into the mcresids
        if mol.atoms[i].resid not in mcresids:
            mcresids.append(mol.atoms[i].resid)

        mcpdbfn = fname.strip('.pdb') + '_res_' + str(i) + '_MetalCenter.pdb'
        if os.path.isfile(mcpdbfn):
            print("Overwritting the metal center pdb file: " + mcpdbfn)
            os.system("rm %s" %mcpdbfn)

        #print the residue which is in the cut off into the pdb file
        for j in mcresids:
            for k in mol.residues[j].resconter:
                tiker = mol.atoms[k].gtype
                atid = mol.atoms[k].atid
                atname = mol.atoms[k].atname
                resname = mol.atoms[k].resname
                chainid = 'A'
                resid = mol.atoms[k].resid
                crdx = round(mol.atoms[k].crd[0], 3)
                crdy = round(mol.atoms[k].crd[1], 3)
                crdz = round(mol.atoms[k].crd[2], 3)
                occp = 1.00
                tempfac = 0.00
                atmj = pdbatm(tiker, atid, atname, resname, chainid, resid,
                              crdx, crdy, crdz, occp, tempfac)
                writepdbatm(atmj, mcpdbfn)

        #Print the environment
        for j in atids:
            atnamej = mol.atoms[j].atname
            crdj = mol.atoms[j].crd
            residj = mol.atoms[j].resid
            resnamej = mol.residues[residj].resname
            elmtj = mol.atoms[j].element
            radiusj = CoRadiiDict[elmtj]
            radiusij = radiusi + radiusj + 0.40
            disij = calc_bond(crdi, crdj)

            if options.cutoff == None:
                if (disij >= 0.1) and (disij <= radiusij) \
                   and (elmtj != 'H'):
                    #print the environment file, for each bond in the metal site
                    print(fname.strip('.pdb'),',', residi,',', resnamei, \
                             ',', i, ',', atnamei,',', residj,',', resnamej,\
                             ',', j,',', atnamej, ',', round(disij, 3),',', geo,\
                             ',', round(georms, 3),',', reslets,',', exptyp,\
                             ',', reso, file=ef)
            else:
                if (disij >= 0.1) and (disij <= options.cutoff) \
                   and (elmtj != 'H'):
                    #print the environment file, for each bond in the metal site
                    print(fname.strip('.pdb'),',', residi,',', resnamei, \
                             ',', i, ',', atnamei,',', residj,',', resnamej,\
                             ',', j,',', atnamej, ',', round(disij, 3),',', geo,\
                             ',', round(georms, 3),',', reslets,',', exptyp,\
                             ',', reso, file=ef)

        #print the summary file, for each metal site
        print(fname.strip('.pdb'),',', exptyp,',', reso, \
                 ',', len(atids), ',', len(metallist),',', residi,\
                 ',', resnamei,',', i,',', atnamei,',', reslets,',', geo,\
                 ',', round(georms, 3), file=sf)

sf.close()
ef.close()
