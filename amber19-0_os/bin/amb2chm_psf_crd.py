#!/opt/anaconda1anaconda2anaconda3/bin/python
# Filename: amb2chm_psf_crd.py
#
# This is a program written by Pengfei Li to convert
# the AMBER prmtop&inpcrd files to the CHARMM psf&crd files
#
import parmed as pmd
from optparse import OptionParser
import os
from pymsmt.mol.amb2chm_typ import (AA_RES_LIST, NA_RES_LIST, NAAA_RES_ATOM_DICT,
                                    ATOM_TYPE_DICT, NA_RES_DICT)

parser = OptionParser("Usage: amb2chm_psf_crd.py -p prmtop -c inpcrd -f psf\n"
                      "                          -d crd -b pdb [--dict dict_file]")
parser.set_defaults(dicf='')
parser.add_option("-p", dest="prmtop", type='string',
                  help="Prmtop file")
parser.add_option("-c", dest="inpcrd", type='string',
                  help="Inpcrd file")
parser.add_option("-f", dest="psf", type='string',
                  help="PSF file")
parser.add_option("-d", dest="crd", type='string',
                  help="CRD file")
parser.add_option("-b", dest="pdb", type='string',
                  help="A PDB file to generate")
parser.add_option("--dict", dest="dicf", type='string',
                  help="Dictionary file name")
(options, args) = parser.parse_args()

amber = pmd.load_file(options.prmtop, options.inpcrd)

# Define the added residue-atom dictionary for metal site residues

ADD_RES_ATOM_DICT = {}

if options.dicf != '':
    with open(options.dicf, 'r') as f:
        for rline in f:
            line = rline.strip('\n')
            line = line.split()
            ADD_RES_ATOM_DICT[(line[0], line[1])] = (line[2], line[3])

"""
ADD_RES_ATOM_DICT = {
            ('HD1', 'H') : ('HD1', 'HN'),   #HD1
            ('HD1', 'HB2') : ('HD1', 'HB1'),
            ('HD1', 'HB3') : ('HD1', 'HB2'),
            ('HD2', 'H') : ('HD2', 'HN'),   #HD2
            ('HD2', 'HB2') : ('HD2', 'HB1'),
            ('HD2', 'HB3') : ('HD2', 'HB2'),
            ('HD3', 'H') : ('HD3', 'HN'),   #HD3
            ('HD3', 'HB2') : ('HD3', 'HB1'),
            ('HD3', 'HB3') : ('HD3', 'HB2'),
       	    ('AN1', 'H')  : ('AN1', 'HN'),    #AN1
	        ('AN1', 'HB2') : ('AN1', 'HB1'),
            ('AN1', 'HB3') : ('AN1', 'HB2'),
	        ('IE1', 'H') : ('IE1', 'HN'),     #IE1
    	    ('IE1', 'CD1') : ('IE1', 'CD'),
    	    ('IE1', 'HD11') : ('IE1', 'HD1'),
    	    ('IE1', 'HD12') : ('IE1', 'HD2'),
    	    ('IE1', 'HD13') : ('IE1', 'HD3'),
    	    ('IE1', 'HG12') : ('IE1', 'HG11'),
            ('IE1', 'HG13') : ('IE1', 'HG12'),
            ('IE1', 'O') : ('IE1', 'OT1'),  #About the terminal group
            ('IE1', 'OXT') : ('IE1', 'OT2'), #About the terminal group
            }
"""

NAAA_RES_ATOM_DICT.update(ADD_RES_ATOM_DICT)

# For the normal residues and atom types
for i in xrange(0, len(amber.residues)):
    resname = amber.residues[i].name

    resid_atid_list = []
    for j in xrange(0, len(amber.residues[i].atoms)):
        atname = amber.residues[i].atoms[j].name
        attype = amber.residues[i].atoms[j].type
        if (resname, atname) in NAAA_RES_ATOM_DICT:
            resid_atid_list.append((i, j))

    if len(resid_atid_list) != 0:
        maxnum = max(resid_atid_list)
    else:
        maxnum = (-1, -1)

    for j in xrange(0, len(amber.residues[i].atoms)):
        atname = amber.residues[i].atoms[j].name
        attype = amber.residues[i].atoms[j].type

        # Change the atom names
        if (resname, atname) in NAAA_RES_ATOM_DICT:
            amber.residues[i].atoms[j].name = NAAA_RES_ATOM_DICT[(resname, atname)][1]

        # Change the residue name
        if ((resname, atname) in NAAA_RES_ATOM_DICT) and ((i, j) == maxnum):
            amber.residues[i].name = NAAA_RES_ATOM_DICT[(resname, atname)][0]

        # Change the atom type
        if attype in ATOM_TYPE_DICT:
            amber.residues[i].atoms[j].type = ATOM_TYPE_DICT[attype]

# For specific residue names (N-terminal and C-terminal amino acids
# & nucleic acids)

NA_atom_set = ['C1\'', 'C2\'', 'C3\'', 'C4\'', 'C5\'']

for i in xrange(0, len(amber.residues)):
    resname = amber.residues[i].name

    # For the N-temrinal and C-terminal amino acids
    if resname in AA_RES_LIST:

        atnames = []
        # Collect the atom names in the residues
        for j in xrange(0, len(amber.residues[i].atoms)):
            atname = amber.residues[i].atoms[j].name
            atnames.append(atname)

        # If it is a N-terminal residue
        if 'H1' in atnames and 'H2' in atnames and 'H3' in atnames:
            amber.residues[i].name = 'N' + resname
            for j in xrange(0, len(amber.residues[i].atoms)):
                atname = amber.residues[i].atoms[j].name
                if atname == 'H1':
                    amber.residues[i].atoms[j].name = 'HT1'
                elif atname == 'H2':
                    amber.residues[i].atoms[j].name = 'HT2'
                elif atname == 'H3':
                    amber.residues[i].atoms[j].name = 'HT3'

        # If it is a C-terminal residue
        elif 'O' in atnames and 'OXT' in atnames:
            amber.residues[i].name = 'C' + resname
            for j in xrange(0, len(amber.residues[i].atoms)):
               atname = amber.residues[i].atoms[j].name
               if atname == 'O':
                   amber.residues[i].atoms[j].name = 'OT1'
               elif atname == 'OXT':
                   amber.residues[i].atoms[j].name = 'OT2'

    # For nucleic acid residues
    elif resname in NA_RES_LIST:

        atnames = []
        # Collect the atom names in the residues
        for j in xrange(0, len(amber.residues[i].atoms)):
            atname = amber.residues[i].atoms[j].name
            atnames.append(atname)

        # If it is a nucleic acid, replace the residue names
        if  set(atnames).intersection(set(NA_atom_set)) == set(NA_atom_set):
            amber.residues[i].name = NA_RES_DICT[resname]

# Save a CHARMM PSF and crd file
amber.save(options.psf, overwrite=True, format='psf')
amber.save(options.crd, overwrite=True, format='charmmcrd')
amber.save(options.pdb, overwrite=True, format='pdb')

quit()

