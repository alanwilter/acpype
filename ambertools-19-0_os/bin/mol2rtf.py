#!/opt/anaconda1anaconda2anaconda3/bin/python
# Filename: mol2rtf.py
#
# This is a program written by Pengfei Li to convert
# the mol2 file to CHARMM rtf file
#
from __future__ import print_function
from optparse import OptionParser
from pymsmt.mol.mol2io import get_atominfo, get_bondinfo
from pymsmt.mol.amb2chm_typ import ATOM_TYPE_DICT

parser = OptionParser("Usage: mol2rtf.py -i mol2_file -o rtf_file -r residue_name\n"
                      "                  -n new_resname [--ref reference_rtf_file]")
parser.set_defaults(ref_rtf='')
parser.add_option("-i", dest="input_file", type='string',
                  help="Input mol2 file")
parser.add_option("-o", dest="output_file", type='string',
                  help="Output RTF file")
parser.add_option("-r", dest="resname", type='string',
                  help="Original residue name")
parser.add_option("-n", dest="new_resname", type='string',
                  help="New residue name")
parser.add_option("--ref", dest="ref_rtf", type='string',
                  help="Reference RTF file")

(options, args) = parser.parse_args()

residues = ["ALA", "ARG", "ASN", "ASP", "CYS", "CYM", "CYX", "GLN", "GLU", "GLH",
            "GLY", "HID", "HIE", "HIP", "ILE", "LEU", "LYS", "MET", "PHE", "PRO",
            "SER", "THR", "TRP", "TYR", "VAL", "ACE", "NME"]

resname = options.resname
new_resname = options.new_resname
tot_chg = 0.0

mol, atids, resids = get_atominfo(options.input_file)
blist = get_bondinfo(options.input_file)

# To read the RTF file for normal amino acids
typ_dict = {}
for i in atids:
    atname = mol.atoms[i].atname
    if atname == 'H':
        atname = 'HN'
    elif atname == 'HB2':
        atname = 'HB1'
    elif atname == 'HB3':
        atname = 'HB2'

    charge = mol.atoms[i].charge
    attyp = mol.atoms[i].atomtype
    typ_dict[atname] = (attyp, charge)
    tot_chg = tot_chg + charge

def read_rtf(fname):
    global resname
    rtf_data = []
    print_yn = 0
    raw_rdf = open(fname)
    for i in raw_rdf:
        if 'RESI ' + resname in i:
            print_yn = 1
        if print_yn == 1:
            rtf_data.append(i)
            #print(i, file=options.output_file)
        if 'PATCH' in i:
            print_yn = 0
    return rtf_data

ref_data = ''
if options.ref_rtf != '':
    ref_data = read_rtf(options.ref_rtf)

print("Old residue name:", resname, "New residue name:", new_resname)

w_output = open(options.output_file, 'w')
if resname in residues:
    # For normal amino acids
    for i in ref_data:
        if 'RESI' in i:
            #RESI HD1          0.00
            #i = i.replace(resname, new_resname)
            #print(i, file=w_output)
            print("%4s %-4s             %-9.6f" %("RESI", new_resname, tot_chg), file=w_output)
        elif 'ATOM' == i[0:4]:
            j = i.split()
            atmtyp = typ_dict[j[1]][0]
            atmchg = typ_dict[j[1]][1]
            if '!' in i:
                com_ind = i.index("!")
                print("%4s %-4s %-2s      %-9.6f  %-s" %('ATOM', j[1], atmtyp, round(atmchg,6), i[com_ind:]), end='', file=w_output)
            else:
                print("%4s %-4s %-2s      %-9.6f  " %('ATOM', j[1], atmtyp, round(atmchg,6)), file=w_output)
        else:
            print(i, end='', file=w_output)
else:
    # For other residues
    #print('RESI', new_resname, file=w_output)
    print("%4s %-4s             %-9.6f" %("RESI", new_resname, tot_chg), file=w_output)
    print('GROUP', file=w_output)
    for i in atids:
        atname = mol.atoms[i].atname
        attype = mol.atoms[i].atomtype
        if attype in ATOM_TYPE_DICT.keys():
            attype = ATOM_TYPE_DICT[attype]
        charge = mol.atoms[i].charge
        print('%4s %-4s %-2s      %-9.6f' %('ATOM', atname, attype, round(charge, 6)), file=w_output)
    for i in blist:
        at1 = i[0]
        at2 = i[1]
        atname1 = mol.atoms[at1].atname
        atname2 = mol.atoms[at2].atname
        print('%4s %-4s %-4s' %('BOND', atname1, atname2), file=w_output)
    print(' ', file=w_output)
w_output.close()

