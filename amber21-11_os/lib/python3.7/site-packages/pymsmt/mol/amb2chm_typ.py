# Filename: amb2chm_typ

###############################################################################
#                            FOR NUCLEIC ACIDS
###############################################################################

# All the nucleic acids will have following names (ADE, THY, CYT, GUA, URA) in
# the CHARMM PSF/PDB/CRD file, while CHARMM needs special handles for the
# nucleic acids in the simulation control file
RNA_RES = ['A', 'A3', 'A5', 'AN', 'C', 'C3', 'C5', 'CN',
           'G', 'G3', 'G5', 'GN', 'U', 'U3', 'U5', 'UN',
           'RA', 'RA3', 'RA5', 'RAN', 'RC', 'RC3', 'RC5', 'RCN',
           'RG', 'RG3', 'RG5', 'RGN', 'RU', 'RU3', 'RU5', 'RUN']

DNA_RES = ['DA', 'DA3', 'DA5', 'DAN', 'DC', 'DC3', 'DC5', 'DCN',
           'DG', 'DG3', 'DG5', 'DGN', 'DT', 'DT3', 'DT5', 'DTN']

NA_RES_LIST = RNA_RES + DNA_RES

NA_RES_DICT = {}
# For RNA
for i in RNA_RES:
    if 'A' in i:
        NA_RES_DICT[i] = 'ADE'
    elif 'C' in i:
        NA_RES_DICT[i] = 'CYT'
    elif 'G' in i:
        NA_RES_DICT[i] = 'GUA'
    elif 'U' in i:
        NA_RES_DICT[i] = 'URA'
# For DNA
for i in DNA_RES:
    if 'A' in i:
        NA_RES_DICT[i] = 'ADE'
    elif 'C' in i:
        NA_RES_DICT[i] = 'CYT'
    elif 'G' in i:
        NA_RES_DICT[i] = 'GUA'
    elif 'T' in i:
        NA_RES_DICT[i] = 'THY'

NA_RES_ATOM_DICT = {}
# For RNA, which does have Oxygen at C2' position
for i in RNA_RES:
    new_res = NA_RES_DICT[i]
    NA_RES_ATOM_DICT[(i, 'HO\'2')] = (new_res, 'H2\'')  #@ff99SB
    NA_RES_ATOM_DICT[(i, 'H2\'1')] = (new_res, 'H2\'\'') #@ff99SB
    NA_RES_ATOM_DICT[(i, 'HO2\'')] = (new_res, 'H2\'') #@ff14SB
    NA_RES_ATOM_DICT[(i, 'H2\'')] = (new_res, 'H2\'\'') #@ff14SB

    NA_RES_ATOM_DICT[(i, 'OP1')] = (new_res, 'O1P') #@ff99SB
    NA_RES_ATOM_DICT[(i, 'OP2')] = (new_res, 'O2P') #@ff99SB
    NA_RES_ATOM_DICT[(i, 'H5\'1')] = (new_res, 'H5\'') #@ff99SB
    NA_RES_ATOM_DICT[(i, 'H5\'2')] = (new_res, 'H5\'\'') #@ff99SB
    NA_RES_ATOM_DICT[(i, 'HO3\'')] = (new_res, 'H3T') #@ff14SB
    NA_RES_ATOM_DICT[(i, 'HO5\'')] = (new_res, 'H5T') #@ff14SB

# For DNA, which does not have Oxygen at C2' position
for i in DNA_RES:
    new_res = NA_RES_DICT[i]
    NA_RES_ATOM_DICT[(i, 'H2\'1')] = (new_res, 'H2\''), #@ff99SB
    NA_RES_ATOM_DICT[(i, 'H2\'2')] = (new_res, 'H2\'\''), #@ff99SB

    NA_RES_ATOM_DICT[(i, 'OP1')] = (new_res, 'O1P') #@ff99SB
    NA_RES_ATOM_DICT[(i, 'OP2')] = (new_res, 'O2P') #@ff99SB
    NA_RES_ATOM_DICT[(i, 'H5\'1')] = (new_res, 'H5\'') #@ff99SB
    NA_RES_ATOM_DICT[(i, 'H5\'2')] = (new_res, 'H5\'\'') #@ff99SB
    NA_RES_ATOM_DICT[(i, 'HO3\'')] = (new_res, 'H3T') #@ff14SB
    NA_RES_ATOM_DICT[(i, 'HO5\'')] = (new_res, 'H5T') #@ff14SB

###############################################################################
#                            FOR AMINO ACIDS
###############################################################################

# In total there are 28 normal AAs, 24 N-temrinal AAs (containing ACE),
# and 26 C-terminal AAs (containg NHE and NME) in the AMBER ff14SB
# for protein

AA_RES_LIST = {'ALA', 'ARG', 'ASH', 'ASN', 'ASP', #There is no NASH, CASH
               'CYM', 'CYS', 'CYX',  #There is no NCYM, CCYM
               'GLH', 'GLN', 'GLU', 'GLY', #There is no NGLH, CGLH
               'HID', 'HIE', 'HIP', 'HYP', #There is no NHYP
               'ILE', 'LEU', 'LYN', 'LYS', #There is no NLYN, CLYN
               'MET', 'PHE', 'PRO', 'SER',
               'THR', 'TRP', 'TYR', 'VAL',
               'ACE', 'NHE', 'NME'}

# ASH is protonated ASP
# GLH is protonated GLU
# CYM is unprotonated CYS
# CYX is S-S bonded CYS
# HID, HIE, and HIP are HISs with different protonated states
# HYP is Hydro-PRO
# LYN is unprotonated LYS
# ACE is CH3-CO
# NHE is NH2
# NME is CH3-NH

AA_RES_ATOM_DICT = {
        # ALA
	    ('ALA', 'H')   : ('ALA', 'HN'), 
        # ARG
	    ('ARG', 'H')   : ('ARG', 'HN'), 
	    ('ARG', 'HB2') : ('ARG', 'HB1'), 
        ('ARG', 'HB3') : ('ARG', 'HB2'),
	    ('ARG', 'HG2') : ('ARG', 'HG1'),
        ('ARG', 'HG3') : ('ARG', 'HG2'),
        ('ARG', 'HD2') : ('ARG', 'HD1'),
        ('ARG', 'HD3') : ('ARG', 'HD2'),
        # ASH
        # There is no ASH in CHARMM parm14sb_all.rtf
        # ASN
     	('ASN', 'H')   : ('ASN', 'HN'),  
	    ('ASN', 'HB2') : ('ASN', 'HB1'),
        ('ASN', 'HB3') : ('ASN', 'HB2'),
        # ASP
	    ('ASP', 'H')   : ('ASP', 'HN'),  
	    ('ASP', 'HB2') : ('ASP', 'HB1'),
        ('ASP', 'HB3') : ('ASP', 'HB2'),
        # CYM
	    ('CYM', 'H')   : ('CYM', 'HN'),  
	    ('CYM', 'HB2') : ('CYM', 'HB1'),
	    ('CYM', 'HB3') : ('CYM', 'HB2'),
        # CYS
	    ('CYS', 'H')   : ('CYS', 'HN'),  
        ('CYS', 'HB2') : ('CYS', 'HB1'),
        ('CYS', 'HB3') : ('CYS', 'HB2'),
	    ('CYS', 'HG')  : ('CYS', 'HG1'),
        # CYX
	    ('CYX', 'H')   : ('CYX', 'HN'),  
        ('CYX', 'HB2') : ('CYX', 'HB1'),
        ('CYX', 'HB3') : ('CYX', 'HB2'),
        # GLH
	    ('GLH', 'H')   : ('GLH', 'HN'),  
	    ('GLH', 'HB2') : ('GLH', 'HB1'),
	    ('GLH', 'HB3') : ('GLH', 'HB2'),
    	('GLH', 'HG2') : ('GLH', 'HG1'),
    	('GLH', 'HG3') : ('GLH', 'HG2'),
        # GLN
	    ('GLN', 'H')   : ('GLN', 'HN'),  
	    ('GLN', 'HB2') : ('GLN', 'HB1'),
	    ('GLN', 'HB3') : ('GLN', 'HB2'),
    	('GLN', 'HG2') : ('GLN', 'HG1'),
    	('GLN', 'HG3') : ('GLN', 'HG2'),
        # GLU
    	('GLU', 'H')   : ('GLU', 'HN'),  
        ('GLU', 'HB2') : ('GLU', 'HB1'),
        ('GLU', 'HB3') : ('GLU', 'HB2'),
        ('GLU', 'HG2') : ('GLU', 'HG1'),
        ('GLU', 'HG3') : ('GLU', 'HG2'),
        # GLY
    	('GLY', 'H')   : ('GLY', 'HN'),  
    	('GLY', 'HA2') : ('GLY', 'HA1'),
        ('GLY', 'HA3') : ('GLY', 'HA2'),
        # HID
    	('HID', 'H')   : ('HID', 'HN'),  
    	('HID', 'HB2') : ('HID', 'HB1'),
        ('HID', 'HB3') : ('HID', 'HB2'),
        # HIE
    	('HIE', 'H')   : ('HIE', 'HN'),  
    	('HIE', 'HB2') : ('HIE', 'HB1'),
        ('HIE', 'HB3') : ('HIE', 'HB2'),
        # HIP
    	('HIP', 'H')   : ('HIP', 'HN'),  
    	('HIP', 'HB2') : ('HIP', 'HB1'),
        ('HIP', 'HB3') : ('HIP', 'HB2'),
        # HYP
        # There is no HYP in CHARMM parm14sb_all.rtf
        # ILE
    	('ILE', 'H')   : ('ILE', 'HN'), 
    	('ILE', 'CD1') : ('ILE', 'CD'),
    	('ILE', 'HD11'): ('ILE', 'HD1'),
    	('ILE', 'HD12'): ('ILE', 'HD2'),
    	('ILE', 'HD13'): ('ILE', 'HD3'),
    	('ILE', 'HG12'): ('ILE', 'HG11'),
        ('ILE', 'HG13'): ('ILE', 'HG12'),
        # LEU
        ('LEU', 'H')   : ('LEU', 'HN'), 
        ('LEU', 'HB2') : ('LEU', 'HB1'),
        ('LEU', 'HB3') : ('LEU', 'HB2'),
        # LYN
        # There is no LYN in CHARMM parm14sb_all.rtf
        # LYS
    	('LYS', 'H')   : ('LYS', 'HN'),
    	('LYS', 'HB2') : ('LYS', 'HB1'),
        ('LYS', 'HB3') : ('LYS', 'HB2'),
    	('LYS', 'HG2') : ('LYS', 'HG1'),
        ('LYS', 'HG3') : ('LYS', 'HG2'),
    	('LYS', 'HD2') : ('LYS', 'HD1'),
        ('LYS', 'HD3') : ('LYS', 'HD2'),
    	('LYS', 'HE2') : ('LYS', 'HE1'),
        ('LYS', 'HE3') : ('LYS', 'HE2'),
        # MET
    	('MET', 'H')   : ('MET', 'HN'),
        ('MET', 'HB2') : ('MET', 'HB1'),
        ('MET', 'HB3') : ('MET', 'HB2'),
    	('MET', 'HG2') : ('MET', 'HG1'),
        ('MET', 'HG3') : ('MET', 'HG2'),
        # PHE
    	('PHE', 'H')   : ('PHE', 'HN'),
        ('PHE', 'HB2') : ('PHE', 'HB1'),
        ('PHE', 'HB3') : ('PHE', 'HB2'),
        # PRO
    	('PRO', 'HD2') : ('PRO', 'HD1'),
        ('PRO', 'HD3') : ('PRO', 'HD2'),
    	('PRO', 'HG2') : ('PRO', 'HG1'),
        ('PRO', 'HG3') : ('PRO', 'HG2'),
    	('PRO', 'HB2') : ('PRO', 'HB1'),
        ('PRO', 'HB3') : ('PRO', 'HB2'),
        # SER
    	('SER', 'H')   : ('SER', 'HN'),  
    	('SER', 'HB2') : ('SER', 'HB1'),
        ('SER', 'HB3') : ('SER', 'HB2'),
    	('SER', 'HG')  : ('SER', 'HG1'),
        # THR
    	('THR', 'H')   : ('THR', 'HN'),  
        # TRP
    	('TRP', 'H')   : ('TRP', 'HN'),  
    	('TRP', 'HB2') : ('TRP', 'HB1'), 
        ('TRP', 'HB3') : ('TRP', 'HB2'),
        # TYR
    	('TYR', 'H')   : ('TYR', 'HN'),  
    	('TYR', 'HB2') : ('TYR', 'HB1'),
        ('TYR', 'HB3') : ('TYR', 'HB2'),
        # VAL
    	('VAL', 'H')   : ('VAL', 'HN'),  
        # NHE
        # There is no NHE in CHARMM parm14sb_all.rtf
        # NME
        ('NME', 'H')   : ('NME', 'HN'),
        ('NME', 'CH3') : ('NME', 'CAT'),
        ('NME', 'HH31'): ('NME', 'HT1'),
        ('NME', 'HH32'): ('NME', 'HT2'),
        ('NME', 'HH33'): ('NME', 'HT3'),
        # ACE
        ('ACE', 'CH3') : ('ACE', 'CAY'),
        ('ACE', 'HH31'): ('ACE', 'HY1'),
        ('ACE', 'HH32'): ('ACE', 'HY2'),
        ('ACE', 'HH33'): ('ACE', 'HY3'),
        # WAT, Na+, Cl-
        ('WAT', 'O')   : ('TIP3', 'OH2'),
        ('WAT', 'H1')  : ('TIP3', 'H1'),
        ('WAT', 'H2')  : ('TIP3', 'H2'),  
    	('Na+', 'Na+') : ('SOD', 'SOD'),
        ('K+', 'K+')   : ('POT', 'POT'),
    	('Cl-', 'Cl-') : ('CLA', 'CLA'),
}

# Combine the dicts for nucleic acids and amino acids
NAAA_RES_ATOM_DICT = {}
NAAA_RES_ATOM_DICT.update(NA_RES_ATOM_DICT)
NAAA_RES_ATOM_DICT.update(AA_RES_ATOM_DICT)

"""
# Terminal AA atom dict
TER_AA_RES_ATOM_DICT = {}
for aa in AA_RES_LIST:
    if aa not in ['ASH', 'CYM', 'GLH', 'HYP', 'LYN', 'ACE', 'NHE', 'NME']:
        TER_AA_RES_ATOM_DICT[(aa,'H1')] = ('N' + aa, 'HT1')
        TER_AA_RES_ATOM_DICT[(aa,'H2')] = ('N' + aa, 'HT2')
        TER_AA_RES_ATOM_DICT[(aa,'H3')] = ('N' + aa, 'HT3')
    if aa not in ['ASH', 'CYM', 'GLH', 'LYN', 'ACE', 'NHE', 'NME']:
        TER_AA_RES_ATOM_DICT[(aa,'O')] = ('C' + aa, 'OT1')
        TER_AA_RES_ATOM_DICT[(aa,'OXT')] = ('C' + aa, 'OT2')
"""

###############################################################################
#                            FOR ATOM TYPES IN FF14SB
###############################################################################

"""
ATOM_TYPE_DICT = {'c' : 'C',
                  'c3': 'CT',
                  'o' : 'O2',
                  'c2': 'CM',
                  'hc': 'HC',
                  'ha': 'HA'}
"""

# About AMBER atom types
ATOM_TYPE_DICT = {'Na+': 'SOD',
                  'K+' : 'POT',
                  'Cl-': 'CLA',
                  'C*' : 'CG',
                  'N*' : 'NG',
                  'Cl' : 'CL',
                  'Br' : 'BR',
                  'Zn' : 'ZN'
                 }

###############################################################################
#                               FOR GAFF ATOM TYPES
###############################################################################

# About GAFF atom types
GAFF_ATOM_TYPE_LIST = [ 'c', 'c1', 'c2', 'c3', 'ca', 'cp', 'cq', 'cc', 'cd',
                       'ce', 'cf', 'cg', 'ch', 'cx', 'cy', 'cu', 'cv', 'cz',
                       'h1', 'h2', 'h3', 'h4', 'h5', 'ha', 'hc', 'hn', 'ho',
                       'hp', 'hs', 'hw', 'hx',
                        'f', 'cl', 'br',  'i',
                        'n', 'n1', 'n2', 'n3', 'n4', 'na', 'nb', 'nc', 'nd',
                       'ne', 'nf', 'nh', 'no', 'ni', 'nj', 'nk', 'nl', 'nm',
                       'nn', 'np', 'nq',
                        'o', 'oh', 'os', 'op', 'oq', 'ow',
                       'p2', 'p3', 'p4', 'p5', 'pb', 'pc', 'pd', 'pe', 'pf',
                       'px', 'py',
                        's', 's2', 's4', 's6', 'sh', 'ss', 'sp', 'sq', 'sx', 'sy']

GAFF_ATOM_TYPE_DICT = {}
for i in GAFF_ATOM_TYPE_LIST:
    GAFF_ATOM_TYPE_DICT[i] = 'GA' + i.upper()

ATOM_TYPE_DICT.update(GAFF_ATOM_TYPE_DICT)

# About wild card for the atom type
ATOM_TYPE_DICT['X'] = 'X'

