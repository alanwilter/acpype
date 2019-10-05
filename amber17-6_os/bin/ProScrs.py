#!/opt/anaconda1anaconda2anaconda3/bin/python
# Filename: ProScrs.py

"""
This is the ProScrs.py program written by Pengfei Li
at University of Illinois at Urbana-Champaign
"""
from pymsmt.mol.element import Atnum
from pymsmt.mol.readpdb import get_atominfo_fpdb
from pymsmt.mcpb.gene_model_files import (write_ace, write_ant, write_act,
     write_nme, write_gly, write_normal, write_sc, write_sc_knh,
     write_sc_kco, write_only_c2h, write_only_n2h, count_lines)
from pymsmt.mol.gauio import (write_gau_cls_optf, write_gau_cls_fcf)
from pymsmt.mol.gmsio import (write_gms_optf, write_gms_fcf)
from pymsmt.mol.mol import get_reslist
from optparse import OptionParser
from title import print_title

parser = OptionParser("Usage: ProScrs.py -i input_file -p PDB_file\n"
                      "                  [-s file_name_prefix] [--symcrd  symcrd]\n"
                      "                  [-c charge] [--fix fix] [--crd0 crd0]")

parser.set_defaults(pre='MOL', chg='0', opth='0', fix='0', crd0='0', symcrd='0')

parser.add_option("-i", dest="inputfile", type='string',
                  help="Input file name")
parser.add_option("-p", dest="pdbfile", type='string',
                  help="PDB file name")
parser.add_option("-s", dest="pre", type='string',
                  help="File name prefix (default: MOL)")
parser.add_option("-c", dest="chg", type='int',
                  help="Charge (default: 0)")
parser.add_option("--symcrd", dest="symcrd", type='int',
                  help="Use symbolic Cartesian coordinates (default: 0)")
parser.add_option("--fix", dest="fix", type='int',
                  help="Fix heavy atoms or not (default: 0): " +
                       "0 means no, " +
                       "1 means only backbone N, CA, C, and O atoms, " +
                       "2 means backbone N, CA, C, O, and sidechain beta atoms, " +
                       "3 means all heavy atoms.")
parser.add_option("--crd0", dest="crd0", type='int',
                  help="Reassign the coordinates with first atom as 0 (default: 0)")

(options, args) = parser.parse_args()

#------------------------------------------------------------------------------
#ANT = CH3NH2
#ACT = CH3CO2-
#SC_KNH = sidechain keep N and H
#SC_KCO = sidechain keep C and O

# Print the title of the program
version = '1.0'
group = 'shs'
print_title('ProScrs.py', version, group)

mol, atids, resids = get_atominfo_fpdb(options.pdbfile)

reslist = get_reslist(mol, resids)

# Define original variables
resace = []
resant = []
resact = []
resnme = []
resgly = []
resnormal = []
res_sc = []
res_sc_knh = []
res_sc_kco = []
res_c2h = []
res_n2h = []

gatms = []

#------------------------------------------------------------------------------
# Read the input file
inputf = open(options.inputfile, 'r')
for line in inputf:
    line = line.split()
    if '\n' in line:
        line.remove('\n')
    if ',' in line:
        line.remove(',')
    if '' in line:
        line.remove('')
    if ' ' in line:
        line.remove(' ')
    if ':' in line:
        line.remove(':')
    # Blank line
    if (len(line) == 0):
        continue
    # Comment
    elif (line[0][0] == '#'):
        continue
    #ACE
    elif line[0].lower() == 'ace':
        if len(line) >= 2:
            try:
                resace = line[1:]
                resace = [int(i) for i in resace]
            except:
                raise pymsmtError('ace need to be integer number(s).')
        else:
            raise pymsmtError('ace need to be provided.')
    #ANT: CH3NH2
    elif line[0].lower() == 'ant':
        if len(line) >= 2:
            try:
                resant = line[1:]
                resant = [int(i) for i in resant]
            except:
                raise pymsmtError('ant need to be integer number(s).')
        else:
            raise pymsmtError('ant need to be provided.')
    #ACT: CH3CO2-
    elif line[0].lower() == 'act':
        if len(line) >= 2:
            try:
                resact = line[1:]
                resact = [int(i) for i in resact]
            except:
                raise pymsmtError('act need to be integer number(s).')
        else:
            raise pymsmtError('act need to be provided.')
    #NME
    elif line[0].lower() == 'nme':
        if len(line) >= 2:
            try:
                resnme = line[1:]
                resnme = [int(i) for i in resnme]
            except:
                raise pymsmtError('nme need to be integer number(s).')
        else:
            raise pymsmtError('nme need to be provided.')
    #GLY
    elif line[0].lower() == 'gly':
        if len(line) >= 2:
            try:
                resgly = line[1:]
                resgly = [int(i) for i in resgly]
            except:
                raise pymsmtError('gly need to be integer number(s).')
        else:
            raise pymsmtError('gly need to be provided.')
    #Keep
    elif line[0].lower() == 'keep':
        if len(line) >= 2:
            try:
                reskeep = line[1:]
                reskeep = [int(i) for i in reskeep]
            except:
                raise pymsmtError('keep need to be integer number(s).')
        else:
            raise pymsmtError('keep need to be provided.')
    #Sidechain 
    elif line[0].lower() == 'sc':
        if len(line) >= 2:
            try:
                res_sc = line[1:]
                res_sc = [int(i) for i in res_sc]
            except:
                raise pymsmtError('sidechain need to be integer number(s).')
        else:
            raise pymsmtError('sidechain need to be provided.')
    #Sidechain keep NH group
    elif line[0].lower() == 'sc_knh':
        if len(line) >= 2:
            try:
                res_sc_knh = line[1:]
                res_sc_knh = [int(i) for i in res_sc_knh]
            except:
                raise pymsmtError('sc_knh need to be integer number(s).')
        else:
            raise pymsmtError('sc_knh need to be provided.')
    #Sidechain keep CO group
    elif line[0].lower() == 'sc_kco':
        if len(line) >= 2:
            try:
                res_sc_kco = line[1:]
                res_sc_kco = [int(i) for i in res_sc_kco]
            except:
                raise pymsmtError('sc_kco need to be integer number(s).')
        else:
            raise pymsmtError('sc_kco need to be provided.')
    #Only keep backbone C and change it to H
    elif line[0].lower() == 'c2h':
        if len(line) >= 2:
            try:
                res_c2h = line[1:]
                res_c2h = [int(i) for i in res_c2h]
            except:
                raise pymsmtError('c2h need to be integer number(s).')
        else:
            raise pymsmtError('c2h need to be provided.')
    #Only keep backbone N and change it to H
    elif line[0].lower() == 'n2h':
        if len(line) >= 2:
            try:
                res_n2h = line[1:]
                res_n2h = [int(i) for i in res_n2h]
            except:
                raise pymsmtError('n2h need to be integer number(s).')
        else:
            raise pymsmtError('n2h need to be provided.')
inputf.close()

#------------------------------------------------------------------------------
# Define the names
outf = options.pre
chg = options.chg

# PDB file
pdbf = outf + '_cluster.pdb'

# Gaussian file names
goptf = outf + '_cluster_opt.com'
gfcf = outf + '_cluster_fc.com'

# GAMESS-US file names
goptf2 = outf + '_cluster_opt.inp'
gfcf2 = outf + '_cluster_fc.inp'

#------------------------------------------------------------------------------

#Write the PDB file
wpdbf = open(pdbf, 'w')

for i in resids:
    if i in resace:
        write_ace(mol, i, gatms, pdbf)
    if i in resant:
        write_ant(mol, i, gatms, pdbf)
    if i in resact:
        write_act(mol, i, gatms, pdbf)
    if i in resnme:
        write_nme(mol, i, gatms, pdbf)
    if i in resgly:
        write_gly(mol, i, gatms, pdbf)
    if i in reskeep:
        write_normal(mol, reslist, i, gatms, pdbf)
    if i in res_sc:
        write_sc(mol, i, gatms, pdbf)
    if i in res_sc_knh:
        write_sc_knh(mol, i, gatms, pdbf)
    if i in res_sc_kco:
        write_sc_kco(mol, i, gatms, pdbf)
    if i in res_n2h:
        write_only_n2h(mol, i, gatms, pdbf)
    if i in res_c2h:
        write_only_c2h(mol, i, gatms, pdbf)
wpdbf.close()

ln = count_lines(pdbf)
print("Totally there are " + str(ln) + " atoms in the model.")

# Get the fix list during the optimization
fixlist = []

mol, atids, resids = get_atominfo_fpdb(pdbf)

natid = 0
for i in atids:
    if options.fix == 1:
        if mol.atoms[i].atname in ['N', 'CA', 'C', 'O']:
            fixlist.append(natid)
    elif options.fix == 2:
        if mol.atoms[i].atname in ['N', 'CA', 'C', 'O', 'CB', 'SB']:
            fixlist.append(natid)
    elif options.fix == 3:
        if mol.atoms[i].element != 'H':
            fixlist.append(natid)
    elif options.fix not in [0, 1, 2, 3]:
        raise pymsmtError('Only supports fix equals 0, 1, 2, or 3.')
    natid = natid + 1

# Whether to locate all the coordinates with first atom as {0 0 0}
if options.crd0 == 1:
    crd0x = gatms[0].crdx
    crd0y = gatms[0].crdy
    crd0z = gatms[0].crdz
    for gatmi in gatms:
        gatmi.crdx = gatmi.crdx - crd0x
        gatmi.crdy = gatmi.crdy - crd0y
        gatmi.crdz = gatmi.crdz - crd0z

# Calculate the spin number and print it into gaussian file
gaelemts = 0
for gatm in gatms:
    AtNum = Atnum[gatm.element]
    gaelemts = gaelemts + AtNum

SpinNum = gaelemts - chg
SpinNum = int(round(SpinNum, 0))
print("Totally there are " + str(SpinNum) + " electrons in the model.")

if SpinNum%2 == 0:
    SpinNum = 1
else:
    SpinNum = 2

# Write Gaussian files and GAMESS-US files
# Gaussian
write_gau_cls_optf(outf, goptf, chg, SpinNum, gatms, fixlist, options.symcrd)
write_gau_cls_fcf(outf, gfcf)

# GAMESS-US
write_gms_optf(goptf2, chg, SpinNum, gatms)
write_gms_fcf(gfcf2, chg, SpinNum)

quit()

