#!/home/conda/feedstock_root/build_artifacts/ambertools_1635195478443/_h_env_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_/bin/python
# Filename: MCPB.py
"""
This is the MCPB.py program written by Pengfei Li in Merz Research Group
at Michigan State University.
It is written to assist the metal center modeling in mixed systems
(especially the protein system).
It is a re-written python program of MCPB in MTK++.
It optimize the workflow of the MCPB and has better supports of different
ions and force feilds.
It supports modeling of more than 80 ions from +1 to +8 oxidation states.
It supports a series of AMBER force fields (ff94, ff99, ff99SB, ff03, ff03.r1,
ff10, ff14ipq, ff14SB, ff14SB.redq, ff14SBonlysc, ff19SB, ff15ipq, ff15ipq-vac,
fb15, GAFF, and GAFF2)

To cite the program please use:
** P. Li, K. M. Merz, JCIM, 2016, 56, 599-604.

The parameterization scheme is adapted from:
** M. B. Peters, Y. Yang, B. Wang, L. Fusti-Molnar, M. N. Weaver, K. M. Merz,
   JCTC, 2010, 6, 2935-2947.

Helpful tips can be found in the following book chapter:
** Pengfei Li and Kenneth M. Merz, Jr., Structural Genomics: General
   Applications, 2021, 257-275.

The Seminario method is from:
** J. M. Seminario, IJQC, 1996, 30, 1271-1277.
"""
#==============================================================================
# Load the MCPB module
#==============================================================================

from pymsmt.mcpb.gene_model_files import get_ms_resnames, gene_model_files
from pymsmt.mcpb.resp_fitting import resp_fitting
from pymsmt.mcpb.gene_pre_frcmod_file import gene_pre_frcmod_file
from pymsmt.mcpb.gene_final_frcmod_file import (gene_by_empirical_way,
          gene_by_QM_fitting_sem, gene_by_QM_fitting_zmatrix, gene_blank)
from pymsmt.mcpb.amber_modeling import gene_leaprc
from pymsmt.mol.element import resnamel
from pymsmt.lib import FF_DICT
from pymsmt.title import print_title
from pymsmt.exp import *
import warnings
import os
from optparse import OptionParser

parser = OptionParser("Usage: MCPB.py -i input_file -s/--step step_number \n"
                      "               [--logf Gaussian/GAMESS-US output logfile] \n"
                      "               [--fchk Gaussian fchk file]")
parser.add_option("-i", dest="inputfile", type='string',
                  help="Input file name")
parser.add_option("-s", "--step", dest="step", type='string',
                  help="Step number")
parser.add_option("--logf", dest="logfile", type='string',
                  help="Gaussian/GAMESS-US output logfile")
parser.add_option("--fchk", dest="fchkfile", type='string',
                  help="Gaussian fchk file")
(options, args) = parser.parse_args()

#==============================================================================
# Get the input variables
#==============================================================================
# Print the title of the program
version = '7.0'
print_title('MCPB.py', version)
options.step = options.step.lower()

# Default values
addred = 0
addbpairs = []
addres = []
anglefc_avg = 0
bondfc_avg = 0
chgfix_resids = []
cutoff = 2.8
ff_choice = 'ff19SB'
frcmodfs = []
g0x = 'gau'
gaff = 1
gname = 'MOL'
ioninfo = []
largeopt = 0
lgchg = -99
lgspin = -99
naamol2fs = []
paraset = '12_6'
scalef = 1.000
smchg = -99
smspin = -99
sqmopt = 0
watermodel = 'opc'
xstru = 0

if options.step not in ['1', '1n', '1m', '1a', '2', '2e', '2s', '2z', '2b',
                        '3', '3a', '3b', '3c', '3d', '4', '4b', '4n1',
                       '4n2']:
    raise pymsmtError('Invalid step number chosen. please choose among the '
                      'following values: 1, 1n, 1m, 1a, 2, 2e, 2s, 2z, 3, '
                      '3a, 3b, 3c, 3d, 4, 4b, 4n1, 4n2')

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
    #---------------------------Required Variables-----------------------------
    #ionids
    elif line[0].lower() == 'ion_ids':
        if len(line) >= 2:
            try:
                ionids = line[1:]
                ionids = [int(i) for i in ionids]
            except:
                raise pymsmtError('ion_ids need to be integer number(s).')
        else:
            raise pymsmtError('ion_ids need to be provided.')
    #ioninfo
    elif line[0].lower() == 'ion_info':
        if (len(line)-1)%4 == 0:
            try:
                line[4::4] = [int(i) for i in line[4::4]]
                ioninfo = line[1:]
            except:
                raise pymsmtError('The charge of the ion in the ion_info '
                                  'should be integer number.')
        else:
            raise pymsmtError('Wrong ion_info format should be: residue name, '
                              'atom name, element, charge.')
    #ionmol2fs
    elif line[0].lower() == 'ion_mol2files':
        if len(line) >= 2:
            ionmol2fs = line[1:]
            for i in ionmol2fs:
                if os.path.exists(i):
                    continue
                else:
                    raise pymsmtError('File %s does not exists.' %i)
        else:
            raise pymsmtError('Need to provide the mol2 file(s) for '
                              'ion_mol2files.')
    #orpdbf
    elif line[0].lower() == 'original_pdb':
        if len(line) == 2:
            orpdbf = line[1]
            if os.path.exists(orpdbf):
                continue
            else:
                raise pymsmtError('File %s does not exists.' %orpdbf)
        else:
            raise pymsmtError('%d pdb files provided for original_pdb, only '
                              'need one.' %(len(line)-1))
    #---------------------------Optional Variables-----------------------------
    # Following options are in the alphabetical order
    #addbpairs
    elif line[0].lower() == 'add_bonded_pairs':
        if len(line) >= 2:
            addbpairs0 = line[1:]
            for i in addbpairs0:
                pairs = i.split('-')
                try:
                    pairs = [int(i) for i in pairs]
                except:
                    raise pymsmtError('Should be integer numbers in the '
                                      'two ends of dash symbol.')
                if len(pairs) != 2:
                    raise pymsmtError('Should be only two numbers in the pairs!')
                pairs = tuple(pairs)
                addbpairs.append(pairs)
        else:
            raise pymsmtError('add_bonded_pairs need to be provided.')
    #addred
    elif line[0].lower() == 'add_redcrd':
        if len(line) == 2:
            try:
                addred = int(line[1])
                if addred not in [0, 1, 2]:
                    raise pymsmtError('add_redcrd varible needs to be 0, 1, or 2. '
                                      '0 means not using, 1 means only add '
                                      'redundant coordinates for optimization of '
                                      'small model. 2 means only do that for force '
                                      'constant calculation of small model.')
            except:
                raise pymsmtError('add_redcrd value is not integer value.')
        elif len(line) == 1:
            warnings.warn('No add_redcrd parameter provided. '
                                'Default value %d is used.' %addred, pymsmtWarning)
        else:
            raise pymsmtError('More than one add_redcrd parameter provided, '
                              'need one.')
    #addres
    elif line[0].lower() == 'additional_resids':
        if len(line) >= 2:
            try:
                addres = line[1:]
                addres = [int(i) for i in addres]
            except:
                raise pymsmtError('additional_resids need to be integer number(s).')
        else:
            raise pymsmtError('additional_resids need to be provided.')
    #anglefc_avg
    elif line[0].lower() == 'anglefc_avg':
        if len(line) == 2:
            try:
                anglefc_avg = int(line[1])
                if anglefc_avg not in [0, 1]:
                    raise pymsmtError('anglefc_avg varible needs to be 0 or 1, '
                                      '0 means not using, 1 means using.')
            except:
                raise pymsmtError('anglefc_avg value is not integer value.')
        elif len(line) == 1:
            warnings.warn('No anglefc_avg parameter provided. '
                          'Default value %d is used.'
                          %anglefc_avg, pymsmtWarning)
        else:
            raise pymsmtError('More than one anglefc_avg parameter provided, '
                              'need one.')
    #bondfc_avg
    elif line[0].lower() == 'bondfc_avg':
        if len(line) == 2:
            try:
                bondfc_avg = int(line[1])
                if bondfc_avg not in [0, 1]:
                    raise pymsmtError('bondfc_avg varible needs to be 0 or 1, '
                                      '0 means not using, 1 means using.')
            except:
                raise pymsmtError('bondfc_avg value is not integer value.')
        elif len(line) == 1:
            warnings.warn('No bondfc_avg parameter provided. '
                          'Default value %d is used.'
                          %bondfc_avg, pymsmtWarning)
        else:
            raise pymsmtError('More than one bondfc_avg parameter provided, '
                              'need one.')
    #chgfix_resids
    elif line[0].lower() == 'chgfix_resids':
        if len(line) >= 2:
            try:
                chgfix_resids = line[1:]
                chgfix_resids = [int(i) for i in chgfix_resids]
            except:
                raise pymsmtError('chgfix_resids need to be integer number(s).')
        else:
            warnings.warn('No chgfix_resids parameter provided. '
                          'Default value %d will be used.'
                          %chgfix_resids, pymsmtWarning)
    #cutoff
    elif line[0].lower() == 'cut_off':
        if len(line) == 2:
            try:
                cutoff = float(line[1])
            except:
                raise pymsmtError('Please provide an float number for the '
                                  'cut_off parameter.')
        elif len(line) == 1:
            warnings.warn('No cut_off parameter provided, Default value '
                          '%5.1f is used.' %cutoff, pymsmtWarning)
        else:
            raise pymsmtError('More than one cut_off values are provided, '
                              'need one.')
    #ff
    elif line[0].lower() == 'force_field':
        if len(line) == 2:
            ff_choice = line[1]
            if ff_choice not in list(FF_DICT.keys()):
                raise pymsmtError('Not support %s force field in current '
                                  'version. Only supports %r.'
                                  %(ff_choice, list(FF_DICT.keys())))
        elif len(line) == 1:
            warnings.warn('No force_field parameter provided. '
                          'Default value %s is used.'
                          %ff_choice, pymsmtWarning)
        else:
            raise pymsmtError('More than one force_field parameter provided, '
                            'need one in %r.' %list(FF_DICT.keys()))
    #frcmodfs
    elif line[0].lower() == 'frcmod_files':
        if len(line) >= 2:
            frcmodfs = line[1:]
            for i in frcmodfs:
                if os.path.exists(i):
                    continue
                else:
                    raise pymsmtError('File %s does not exists.' %i)
        else:
            warnings.warn('No frcmod files is provided for '
                          'frcmod_files.', pymsmtWarning)
    #gaff
    elif line[0].lower() == 'gaff':
        if len(line) == 2:
            try:
                gaff = int(line[1])
                if gaff not in [0, 1, 2]:
                    raise pymsmtError('gaff varible needs to be 0, 1, or 2. '
                                      '0 means not using, 1 means using gaff, '
                                      '2 means using gaff2.')
            except:
                raise pymsmtError('gaff value is not integer value.')
        elif len(line) == 1:
            warnings.warn('No gaff parameter provided. '
                          'Default value %d is used.'
                          %gaff, pymsmtWarning)
        else:
            raise pymsmtError('More than one gaff parameter provided, '
                              'need one.')
    #group: fname
    elif line[0].lower() == 'group_name' :
        if len(line) == 2:
            gname = line[1]
        elif len(line) == 1:
            warnings.warn('None group_name provided in the input file, '
                          'using %s as default.' %gname, pymsmtWarning)
        else:
            raise pymsmtError('More than one group_name provided, need one.')
    #paraset
    elif line[0].lower() == 'ion_paraset':
        if len(line) == 2:
            paraset = line[1].lower()
            if paraset not in ['hfe', 'cm', 'iod', '12_6', '12_6_4']:
                raise pymsmtError('Do not have %s ion parameter set. Only '
                                  'have HFE, CM, IOD, 12_6, and 12_6_4 '
                                  'parameter sets, where the 12_6 set '
                                  'equivalent to the CM set.'
                                  %paraset.upper())
        elif len(line) == 1:
            warnings.warn('No ion_paraset parameter provided. '
                          'Default value %s is used.'
                          %paraset.upper(), pymsmtWarning)
        else:
            raise pymsmtError('More than one ion_paraset parameter provided, '
                              'only need one: HFE, CM, IOD, 12_6, or 12_6_4, '
                              'where the 12_6 set is equivalent to the CM set.')
    #large_opt
    elif line[0].lower() == 'large_opt':
        if len(line) == 2:
            try:
                largeopt = int(line[1])
                if largeopt not in [0, 1, 2]:
                    raise pymsmtError('large_opt varible needs to be 0, 1 or '
                                      '2. 0 means not using optimization in '
                                      'Guassian input file for large model. 1 '
                                      'means only optimize hydrogen positions '
                                      'in Gaussian input file for large model. '
                                      '2 means optimize the whole structure '
                                      'in Gaussian input file for large model.')
            except:
                raise pymsmtError('large_opt value is not integer value.')
        elif len(line) == 1:
            warnings.warn('No large_opt parameter provided. '
                                'Default value %d is used.' %largeopt, pymsmtWarning)
        else:
            raise pymsmtError('More than one large_opt parameter provided, '
                              'need one.')
    #lgchg
    elif line[0].lower() == 'lgmodel_chg':
        if len(line) == 2:
            try:
                lgchg = int(line[1])
            except:
                raise pymsmtError('Please provide an int number for the '
                                  'lgmodel_chg parameter.')
        elif len(line) == 1:
            warnings.warn('No lgmodel_chg parameter provided, program '
                          'will assign a charge automatically.', pymsmtWarning)
        else:
            raise pymsmtError('More than one lgmodel_chg values are provided, '
                              'need one.')
    #lgspin
    elif line[0].lower() == 'lgmodel_spin':
        if len(line) == 2:
            try:
                lgspin = int(line[1])
                if lgspin <= 0:
                    raise pymsmtError('Please provide an positive int number for '
                                      'the lgmodel_spin parameter.')
            except:
                raise pymsmtError('Please provide an int number for the '
                                  'lgmodel_spin parameter.')
        elif len(line) == 1:
            warnings.warn('No lgmodel_spin parameter provided, program '
                          'will assign a spin automatically.', pymsmtWarning)
        else:
            raise pymsmtError('More than one lgmodel_spin values are provided, '
                              'need one.')
    #naamol2fs
    elif line[0].lower() == 'naa_mol2files':
        if len(line) >= 2:
            naamol2fs = line[1:]
            for i in naamol2fs:
                if os.path.exists(i):
                    continue
                else:
                    raise pymsmtError('File %s does not exists.' %i)
        else:
            warnings.warn('No mol2 file is provided for '
                          'naa_mol2files.', pymsmtWarning)
    #scale_factor
    elif line[0].lower() == 'scale_factor':
        if len(line) == 2:
            try:
                scalef = float(line[1])
            except:
                raise pymsmtError('Please provide an float number for the '
                                  'scale_factor parameter.')
        elif len(line) == 1:
            warnings.warn('No scale_factor parameter provided, Default value '
                          '%5.1f is used.' %scalef, pymsmtWarning)
        else:
            raise pymsmtError('More than one scale_factor values are provided, '
                              'need one.')
    #smchg
    elif line[0].lower() == 'smmodel_chg':
        if len(line) == 2:
            try:
                smchg = int(line[1])
            except:
                raise pymsmtError('Please provide an int number for the '
                                  'smmodel_chg parameter.')
        elif len(line) == 1:
            warnings.warn('No smmodel_chg parameter provided, program '
                          'will assign a charge automatically.', pymsmtWarning)
        else:
            raise pymsmtError('More than one smmodel_chg values are provided, '
                              'need one.')
    #smspin
    elif line[0].lower() == 'smmodel_spin':
        if len(line) == 2:
            try:
                smspin = int(line[1])
                if smspin <= 0:
                    raise pymsmtError('Please provide an positive int number '
                                      'for the smmodel_spin parameter.')
            except:
                raise pymsmtError('Please provide an int number for the '
                                  'smmodel_spin parameter.')
        elif len(line) == 1:
            warnings.warn('No smmodel_spin parameter provided, program '
                          'will assign a spin automatically.', pymsmtWarning)
        else:
            raise pymsmtError('More than one smmodel_spin values are provided, '
                              'need one.')
    #software_version
    elif line[0].lower() == 'software_version':
        if len(line) == 2:
            g0x = line[1].lower()
            if g0x not in ['g03', 'g09', 'g16', 'gau', 'gms']:
                raise pymsmtError('Please use either g03, g09, g16, gau, or '
                                  'gms, other software/versions are not '
                                  'guaranteed to be supported.')
        elif len(line) == 1:
            warnings.warn('No g0x parameter provided. Default value '
                          '%s is used.' %g0x, pymsmtWarning)
        else:
            raise pymsmtError('More than one software_version parameters '
                              'provided, need one.')
    #sqmopt
    elif line[0].lower() == 'sqm_opt':
        if len(line) == 2:
            try:
                sqmopt = int(line[1])
                if sqmopt not in [0, 1, 2, 3]:
                    raise pymsmtError('sqm_opt varible needs to be 0, 1, 2 or '
                                      '3, 0 means not using, 1 means only do '
                                      'optimization for small model. 2 '
                                      'means only do optimization for large '
                                      'model. 3 means do optimization for '
                                      'both models.')
            except:
                raise pymsmtError('sqm_opt value is not integer value.')
        elif len(line) == 1:
            warnings.warn('No sqm_opt parameter provided. '
                                'Default value %d is used.' %sqmopt, pymsmtWarning)
        else:
            raise pymsmtError('More than one sqm_opt parameter provided, '
                              'need one.')
    #watermodel
    elif line[0].lower() == 'water_model':
        if len(line) == 2:
            watermodel = line[1].lower()
            if watermodel not in ['tip3p', 'spce', 'tip4pew', 'opc3', 'opc', 'fb3', 'fb4']:
                raise pymsmtError('Not support %s water model. Only support '
                                  'TIP3P, SPCE, TIP4PEW, OPC3, OPC, FB3, and FB4.' %watermodel.upper())
        elif len(line) == 1:
            warnings.warn('No water_model parameter provided. Default '
                          'value %s is used.'
                          %watermodel.upper(), pymsmtWarning)
        else:
            raise pymsmtError('More than one water_model parameter provided, '
                              'only need one: TIP3P, SPCE, TIP4PEW, OPC3, OPC, FB3, or FB4.')
    #xstru
    elif line[0].lower() == 'xstru':
        if len(line) == 2:
            try:
                xstru = int(line[1])
                if xstru not in [0, 1]:
                    raise pymsmtError('xstru varible needs to be 0 or 1, '
                                      '0 means not using, 1 means using.')
            except:
                raise pymsmtError('xstru value is not integer value.')
        elif len(line) == 1:
            warnings.warn('No xstru parameter provided. '
                          'Default value %d is used.'
                          %xstru, pymsmtWarning)
        else:
            raise pymsmtError('More than one xstru parameter provided, '
                              'need one.')
inputf.close()

#Print the input variables
print("The input file you are using is : %s" %options.inputfile)
print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
print("The following is the input variable you have:")

#-----------------------------Required Variables-------------------------------
try:
    print('The variable ion_ids is : ', ionids)
except:
    raise pymsmtError('ion_ids needs to be provided.')

if options.step in ['4n2']:
    if ioninfo == []:
        raise pymsmtError('The variable ion_info need to be provided in step '
                          '%s.' %options.step)
else:
    print('The variable ion_info is : ', ioninfo)

try:
    print('The variable ion_mol2files is : ', ionmol2fs)
except:
    raise pymsmtError('ion_mol2files needs to be provided.')

try:
    print('The variable original_pdb is : ', orpdbf)
except:
    raise pymsmtError('original_pdb needs to be provided.')

#-----------------------------Optional Variables-------------------------------
print('The variable add_bonded_pairs is : ', addbpairs)
print('The variable add_redcrd is : ', addred)
print('The variable additional_resids is : ', addres)
print('The variable anglefc_avg is : ', anglefc_avg)
print('The variable bondfc_avg is : ', bondfc_avg)
print('The variable chgfix_resids is : ', chgfix_resids)
print('The variable cut_off is : ', cutoff)
print('The variable force_field is : ', ff_choice)
print('The variable frcmod_files is : ', frcmodfs)
print('The variable gaff is : ', gaff)
print('The variable group_name is : ', gname)
print('The variable ion_paraset is : ', paraset.upper(),
      "(Only for the ions using the nonbonded model).")
print('The variable large_opt is : ', largeopt)
print('The variable lgmodel_chg is : ', lgchg)
print('The variable lgmodel_spin is : ', lgspin)
print('             -99 means program will assign a charge automatically.')
print('The variable naa_mol2files is : ', naamol2fs)
print('The variable scale_factor is : ', scalef)
print('             ATTENTION: This is the scale factor of frequency. The ')
print('             force constants will be scaled by multiplying the square')
print('             of scale_factor.')
print('The variable smmodel_chg is : ', smchg)
print('The variable smmodel_spin is : ', smspin)
print('             -99 means program will assign a charge automatically.')
print('The variable software_version is : ', g0x)
print('The variable sqm_opt is : ', sqmopt)
print('The variable water_model is : ', watermodel.upper())
print('The variable xstru is : ', xstru)
#==============================================================================
# Related define
#==============================================================================
#Get the renamed residue name
mcresname0, mcresname = get_ms_resnames(orpdbf, ionids, cutoff, addres,
                                        addbpairs)
for i in mcresname0:
    if (i not in resnamel) and (i+'.mol2' not in naamol2fs):
        raise pymsmtError('%s is required in naa_mol2files but not '
                          'provided.' %i)

#The mol2 file used in pre-generated model
premol2fs = ionmol2fs + naamol2fs

##pdb files
smpdbf = gname + '_small.pdb'
stpdbf = gname + '_standard.pdb'
lgpdbf = gname + '_large.pdb'
fipdbf = gname + '_mcpbpy.pdb'

##finger print files
stfpf = gname + '_standard.fingerprint'
lgfpf = gname + '_large.fingerprint'

##residue information file
smresf = gname + '_small.res'

##frcmod files
prefcdf = gname + '_mcpbpy_pre.frcmod'
finfcdf = gname + '_mcpbpy.frcmod'

##log file
if options.logfile is not None:
    fclogf = options.logfile
    mklogf = options.logfile
else:
    fclogf = gname + '_small_fc.log'
    mklogf = gname + '_large_mk.log'

##checkpoint file
if options.fchkfile is not None:
    fcfchkf = options.fchkfile
else:
    fcfchkf = gname + '_small_opt.fchk'

##tleap input file
ileapf = gname + '_tleap.in'

#==============================================================================
# Step 1 General_modeling
#==============================================================================
#1. Generate the modeling files:
#Pdb files for small, standard and large model
#Gaussian input files for small, large model
#Fingerprint files for standard and large model
#Three options:
#1n) Don't rename any of the atom types in the fingerprint file of standard
#    model
#1m) Just rename the metal ion to the AMBER ion atom type style
#1a) Default. Automatically rename the atom type of the atoms in the metal
#    complex.
if (options.step == '1n'):
    gene_model_files(orpdbf, ionids, addres, addbpairs, gname, ff_choice,
        premol2fs, cutoff, watermodel, 0, largeopt, sqmopt, smchg, smspin,
        addred, lgchg, lgspin)
elif (options.step == '1m'):
    gene_model_files(orpdbf, ionids, addres, addbpairs, gname, ff_choice,
        premol2fs, cutoff, watermodel, 1, largeopt, sqmopt, smchg, smspin,
        addred, lgchg, lgspin)
elif (options.step in ['1', '1a']): #Default
    gene_model_files(orpdbf, ionids, addres, addbpairs, gname, ff_choice,
        premol2fs, cutoff, watermodel, 2, largeopt, sqmopt, smchg, smspin,
        addred, lgchg, lgspin)
#==============================================================================
# Step 2 Frcmod file generation
#==============================================================================
#Mass, dihedral, improper, VDW and metal ion non-related bond and metal
#ion non-related angle parameters are generated first. While the metal ion
#related bond and angle parameters are generated later while they could
#generated by using different methods.
#2e) Empirical method developed by Pengfei Li and co-workers in Merz group
#2s) Default. Seminario method developed by Seminario in 1990s
#2z) Z-matrix method
elif (options.step in ['2', '2s', '2e', '2z', '2b']):
    gene_pre_frcmod_file(ionids, premol2fs, stpdbf, stfpf, smresf, prefcdf,
            ff_choice, gaff, frcmodfs, watermodel)
    if options.step == '2e':
        gene_by_empirical_way(smpdbf, ionids, stfpf, prefcdf, finfcdf)
    elif (options.step in ['2', '2s']): #Default
        gene_by_QM_fitting_sem(smpdbf, ionids, xstru, stfpf, prefcdf, finfcdf,
            fcfchkf, fclogf, g0x, scalef, bondfc_avg, anglefc_avg)
    elif (options.step == '2z'):
        gene_by_QM_fitting_zmatrix(smpdbf, ionids, xstru, stfpf, prefcdf,
            finfcdf, fclogf, scalef)
    elif (options.step == '2b'):
        gene_blank(smpdbf, ionids, xstru, stfpf, prefcdf, finfcdf)
#==============================================================================
# Step 3 Doing the RESP charge fitting and generate the mol2 files
#==============================================================================
#3. Generate mol2 files with the charge parameters after resp charge fitting
#3a) All all the charges of the ligating residues could change
#3b) Default. Restrains the charges of backbone heavy atoms according to the
#    force field chosen
#3c) Restrains the charges of backbone atoms according to the force field
#    chosen
#3d) Restrains the charges of backbone atoms and CB atom in the sidechain
#    according to force field chosen
elif (options.step == '3a'):
    resp_fitting(stpdbf, lgpdbf, stfpf, lgfpf, mklogf, ionids, ff_choice,
        premol2fs, mcresname, 0, chgfix_resids, g0x, lgchg)
elif (options.step in ['3', '3b']): #Default
    resp_fitting(stpdbf, lgpdbf, stfpf, lgfpf, mklogf, ionids, ff_choice,
        premol2fs, mcresname, 1, chgfix_resids, g0x, lgchg)
elif (options.step == '3c'):
    resp_fitting(stpdbf, lgpdbf, stfpf, lgfpf, mklogf, ionids, ff_choice,
        premol2fs, mcresname, 2, chgfix_resids, g0x, lgchg)
elif (options.step == '3d'):
    resp_fitting(stpdbf, lgpdbf, stfpf, lgfpf, mklogf, ionids, ff_choice,
        premol2fs, mcresname, 3, chgfix_resids, g0x, lgchg)
#==============================================================================
# Step 4 Prepare the modeling file for leap
#==============================================================================
#4. Prepare the final modeling file
#4b) Default. Bonded model
#4n1) Nonbonded model with refitting the charge in the protein complex
#4n2) Normal Nonbonded model (12-6 nonbonded model) without re-fitting charges
elif (options.step in ['4', '4b']): #bonded model, Default
    gene_leaprc(gname, orpdbf, fipdbf, stpdbf, stfpf, ionids, ionmol2fs,
        ioninfo, mcresname, naamol2fs, ff_choice, gaff, frcmodfs, finfcdf,
        ileapf, 1, watermodel, paraset)
elif (options.step == '4n1'): #nonbonded model with refitting the charge
    gene_leaprc(gname, orpdbf, fipdbf, stpdbf, stfpf, ionids, ionmol2fs,
        ioninfo, mcresname, naamol2fs, ff_choice, gaff, frcmodfs, finfcdf,
        ileapf, 2, watermodel, paraset)
elif (options.step == '4n2'): #normal nonbonded model
    gene_leaprc(gname, orpdbf, fipdbf, stpdbf, stfpf, ionids, ionmol2fs,
        ioninfo, mcresname, naamol2fs, ff_choice, gaff, frcmodfs, finfcdf,
        ileapf, 3, watermodel, paraset)

#Print the reference
print("="*66)
print("To cite MCPB.py please use the following reference:")
print("  Pengfei Li and Kenneth M. Merz, Jr., \"MCPB.py: A Python Based")
print("  Metal Center Parameter Builder\", J. Chem. Inf. Model., 2016,")
print("  56, 599-604.")
print("Helpful tips can be found in the following book chapter:")
print("  Pengfei Li and Kenneth M. Merz, Jr., \"Parameterization of a Dioxygen")
print("  Binding Metal Site Using the MCPB.py Program\", Methods Mol. Biol.,")
print("  2021, 2199, 257-275.")

quit()
