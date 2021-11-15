#!/Users/runner/miniforge3/conda-bld/ambertools_1635195607525/_h_env_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehol/bin/python
# Filename: amb2chm_par.py
#
# This is a program written by Pengfei Li to convert
# the AMBER dat/frcmod file to CHARMM par file
#

import parmed as pmd
from optparse import OptionParser
from pymsmt.mol.amb2chm_typ import ATOM_TYPE_DICT

parser = OptionParser("Usage: amb2chm_par.py -i input_file [-f input_file_option]\n"
                      "                      -o output_file [--nat use_new_attype]")
parser.set_defaults(fopt=2, newtype=1)
parser.add_option("-i", dest="inputf", type='string',
                  help="The input file")
parser.add_option("-f", dest="fopt", type='int',
                  help="The input file is a parameter file (1) or just contains "
                       "file names (2) [default: 2]")
parser.add_option("-o", dest="outputf", type='string',
                  help="The output file")
parser.add_option("--nat", dest="newtype", type='int',
                  help="Whether to perform atom type transfer [0 means no, 1 "
                       "means yes, default: 1]")
(options, args) = parser.parse_args()

def get_new_attype(i):
    global ATOM_TYPE_DICT
    try:
        ni = ATOM_TYPE_DICT[i]
    except:
        ni = i
    return ni

def params_update(cparams, aparams, newattyp, imp_rearg):

    for bond in list(aparams.bond_types.keys()):
        i = bond[0]
        j = bond[1]
        if newattyp == 0:
            ni = i
            nj = j
        elif newattyp == 1:
            ni = get_new_attype(i)
            nj = get_new_attype(j)
        if min((ni, nj), (nj, ni)) not in list(cparams.bond_types.keys()):
            cparams.bond_types[min((ni, nj), (nj, ni))] = aparams.bond_types[bond]

    for angle in list(aparams.angle_types.keys()):
        i = angle[0]
        j = angle[1]
        k = angle[2]
        if newattyp == 0:
            ni = i
            nj = j
            nk = k
        elif newattyp == 1:
            ni = get_new_attype(i)
            nj = get_new_attype(j)
            nk = get_new_attype(k)
        if min((ni, nj, nk), (nk, nj, ni)) not in list(cparams.angle_types.keys()):
            cparams.angle_types[min((ni, nj, nk), (nk, nj, ni))] = aparams.angle_types[angle]

    for dih in list(aparams.dihedral_types.keys()):
        i = dih[0]
        j = dih[1]
        k = dih[2]
        l = dih[3]
        if newattyp == 0:
            ni = i
            nj = j
            nk = k
            nl = l
        elif newattyp == 1:
            ni = get_new_attype(i)
            nj = get_new_attype(j)
            nk = get_new_attype(k)
            nl = get_new_attype(l)
        if min((ni, nj, nk, nl), (nl, nk, nj, ni)) not in list(cparams.dihedral_types.keys()):
             cparams.dihedral_types[min((ni, nj, nk, nl), (nl, nk, nj, ni))] = aparams.dihedral_types[dih]

    #for dih in cparams.dihedral_types.keys():
        #print(dih.scee)
        #print(cparams.dihedral_types[dih].scee)
        #cparams.dihedral_types[dih].scee = 1.2
        #cparams.dihedral_types[dih].scnb = 2.0

    for imp in list(aparams.improper_periodic_types.keys()):
        i = imp[0]
        j = imp[1]
        k = imp[2]
        l = imp[3]

        if newattyp == 0:
            ni = i
            nj = j
            nk = k
            nl = l
        elif newattyp == 1:
            ni = get_new_attype(i)
            nj = get_new_attype(j)
            nk = get_new_attype(k)
            nl = get_new_attype(l)

        # Re-arrange the improper torsions by alphabetical order
        if imp_rearg == 1:
            # Assume at3 is the central atom, which must be the case for the
            # AMBER force field
            assert nk != 'X', 'Error about the sequence of the improper torsions, ' + \
                              'X can not be the central atom!'
            attyp = [ni, nj, nk, nl]
            attyp1 = [ni, nj, nl]
            xlist = [attyp[i] for i in range(4) if (attyp[i] == 'X')]
            rattyp1 = [attyp1[i] for i in range(3) if (attyp1[i] != 'X')]

            if len(xlist) == 2:
                  ni = 'X'
                  nj = 'X'
                  nl = rattyp1[0]
            elif len(xlist) == 1:
                  ni = 'X'
                  nj, nl = sorted(rattyp1)
            else:
                  ni, nj, nl = sorted(rattyp1)

        # These are several impropers in the current AMBER FFs which are not
        # follow the normal improper naming scheme - they have the ni, nj,
        # and nk are not in the alphabetical order. These cases are dealt
        # as the following
        imp_nonreg = {  ('2C', 'CA', 'CA', 'CA'): ('CA', 'CA', 'CA', '2C'), #CA-CA-CA-2C in frcmod.ff14SB
                        ('C5', 'CB', 'NG', 'CT'): ('CB', 'C5', 'NG', 'CT'), #CB-C5-N*-CT in parm10.dat
                        ('CC', 'CR', 'NA', 'P') : ('CR', 'CC', 'NA', 'P'), #CR-CC-NA-P in parm10.dat
                        ('C2', 'CB', 'NG', 'CT'): ('CB', 'C2', 'NG', 'CT'), #CB-C2-N*-CT in frcmod.DNA.OL15
                        ('C1', 'CB', 'NG', 'CT'): ('CB', 'C1', 'NG', 'CT'), #CB-C1-N*-CT in frcmod.parmbsc1
                        ('GAC', 'GAC2', 'GAC2', 'GAC3'): ('GAC2', 'GAC', 'GAC2', 'GAC3'),  #c2-c -c2-c3 in GAFF
                        ('GAC2', 'GAC3', 'GAC2', 'GAN2') : ('GAC3', 'GAC2', 'GAC2', 'GAN2'), #c3-c2-c2-n2 in GAFF
                        ('GAC2', 'GAC3', 'GAC2', 'GANA') : ('GAC3', 'GAC2', 'GAC2', 'GANA'), #c3-c2-c2-na in GAFF
                        ('GAC2', 'GACA', 'GACA', 'GACA') : ('GACA', 'GACA', 'GACA', 'GAC2'), #ca-ca-ca-c2 in GAFF
                        ('GAC3', 'GACA', 'GACA', 'GACA') : ('GACA', 'GACA', 'GACA', 'GAC3'), #ca-ca-ca-c3 in GAFF
                        ('GAC3', 'GACA', 'GANA', 'GACA') : ('GACA', 'GACA', 'GANA', 'GAC3'), #ca-ca-na-c3 in GAFF
                        ('GAC2', 'GAN2', 'GACA', 'GAN2') : ('GAN2', 'GAC2', 'GACA', 'GAN2'), #n2-c2-ca-n2 in GAFF
                        ('GACA', 'GAN2', 'GACA', 'GAN2') : ('GAN2', 'GACA', 'GACA', 'GAN2'), #n2-ca-ca-n2 in GAFF
                        ('GAN2', 'GAN2', 'GACA', 'GANA') : ('GANA', 'GAN2', 'GACA', 'GAN2'), #na-n2-ca-n2 in GAFF
                        }

        if (ni, nj, nk, nl) in imp_nonreg:
            (ni, nj, nk, nl) = imp_nonreg[(ni, nj, nk, nl)]

        if (ni, nj, nk, nl) not in list(cparams.improper_periodic_types.keys()):
            cparams.improper_periodic_types[(ni, nj, nk, nl)] = aparams.improper_periodic_types[imp]

    for atom in list(aparams.atom_types.keys()):
        if newattyp == 0:
            natom = atom
        elif newattyp == 1:
            natom = get_new_attype(atom)
        if natom not in list(cparams.atom_types.keys()):
            cparams.atom_types[natom] = aparams.atom_types[atom]
            cparams.atom_types[natom].name = natom

    for nbfix in list(aparams.nbfix_types.keys()):
        i = nbfix[0]
        j = nbfix[1]
        if newattyp == 0:
            ni = i
            nj = i
        elif newattyp == 1:
            ni = get_new_attype(i)
            nj = get_new_attype(j)
        if min((ni, nj), (nj, ni)) not in list(cparams.nbfix_types.keys()):
            cparams.nbfix_types[min((ni, nj), (nj, ni))] = aparams.nbfix_types[nbfix]

    return cparams

# If transfer CHARMM to AMBER, then the dict should be reversed (Not supported now)
#if options.fi.lower() == 'charmm' and options.fo.lower() == 'amber':
#    ATOM_TYPE_DICT2 = {}
#    for i in ATOM_TYPE_DICT.keys():
#        j = ATOM_TYPE_DICT[i]
#        ATOM_TYPE_DICT2[j] = i
#    ATOM_TYPE_DICT = ATOM_TYPE_DICT2
#    del ATOM_TYPE_DICT2

if options.fopt not in [1, 2]:
    raise ValueError('-f should be followed by 1 or 2!')
if options.newtype not in [0, 1]:
    raise ValueError('--nat should be followed by 0 or 1!')

# Setting the option to re-arrange the improper torsion order
imp_rearg = 1

# Generate an empty CHARMM parameter set
params = pmd.charmm.CharmmParameterSet()

# Update the parameter files
if options.fopt == 1:
    params1 = pmd.amber.AmberParameterSet(options.inputf)
    params = params_update(params, params1, options.newtype, imp_rearg)
elif options.fopt == 2:
    frc_files = []
    with open(options.inputf, 'r') as f:
        for rline in f:
            line = rline.strip('\n')
            line = line.strip()
            frc_files.append(line)
        for i in frc_files:
            params1 = pmd.amber.AmberParameterSet(i)
            params = params_update(params, params1, options.newtype, imp_rearg)

# Write the output parameter file
params.write(par=options.outputf)

quit()

