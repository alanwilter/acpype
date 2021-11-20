"""
module for define the normally used data about elements, amino acids et al.
"""

from parmed.periodic_table import AtomicNum

#-----------------------------------------------------------------------------
#The basic information of ion
#1.Ion names which has the VDW parameters
#2.Mass, making it compatible with AMBER FFs
#3.PDB names, the usually resname and atname for the metal ions
#-----------------------------------------------------------------------------

OrganicMass = { 'H' : 1.008,
                'C' : 12.01,
                'N' : 14.01,
                'O' : 16.00,
                'S' : 32.06,
                'P' : 30.97,
              }

METAL_PDB = {
                ('LI', 'LI')  :  ('Li', 3, 6.94),
                ('BE', 'BE')  :  ('Be', 4, 9.01),
                ('F' , 'F' )  :  ('F' , 9, 19.00),
                ('NA', 'NA')  :  ('Na', 11, 22.99),
                ('MG', 'MG')  :  ('Mg', 12, 24.305),
                ('AL', 'AL')  :  ('Al', 13, 26.98),
                ('SI', 'SI')  :  ('Si', 14, 28.09), #Add on
                ('CL', 'CL')  :  ('Cl', 17, 35.45),
                ('K', 'K')    :  ('K',  19, 39.10),
                ('CA', 'CA')  :  ('Ca', 20, 40.08),
                ('SC', 'SC')  :  ('Sc', 21, 44.96), #Add
                ('TI', 'TI')  :  ('Ti', 22, 47.87), #Add
                ('V', 'V')    :  ('V',  23, 50.94),
                ('CR', 'CR')  :  ('Cr', 24, 52.00),
                ('MN', 'MN')  :  ('Mn', 25, 54.94),
                ('FE2', 'FE2'):  ('Fe', 26, 55.85), #Fe2+
                ('FE', 'FE')  :  ('Fe', 26, 55.85), #Fe3+
                ('CO', 'CO')  :  ('Co', 27, 58.93),
                ('NI', 'NI')  :  ('Ni', 28, 58.69),
                ('CU1', 'CU') :  ('Cu', 29, 63.55), #Cu+
                ('CU', 'CU')  :  ('Cu', 29, 63.55), #Cu2+
                ('ZN', 'ZN')  :  ('Zn', 30, 65.4),
                ('GA', 'GA')  :  ('Ga', 31, 69.72),
                ('GE', 'GE')  :  ('Ge', 32, 72.63), #Add on
                ('ARS', 'AS') :  ('As', 33, 74.92),
                ('SE', 'SE')  :  ('Se', 34, 78.97),
                ('BR', 'BR')  :  ('Br', 35, 79.90),
                ('RB', 'RB')  :  ('Rb', 37, 85.47),
                ('SR', 'SR')  :  ('Sr', 38, 87.62),
                ('Y', 'Y')    :  ('Y',  39, 88.91),
                ('ZR', 'ZR')  :  ('Zr', 40, 91.22), #Add on
                ('NB', 'NB')  :  ('Nb', 41, 92.91), #Add on
                ('MO', 'MO')  :  ('Mo', 42, 95.96),
                ('TC', 'TC')  :  ('Tc', 43, 98.0), #Add on
                ('RU', 'RU')  :  ('Ru', 44, 101.07),
                ('RH1', 'RH') :  ('Rh', 45, 102.91), #Rh+
                ('RH', 'RH')  :  ('Rh', 45, 102.91),
                ('RH3', 'RH') :  ('Rh', 45, 102.91), #Rh3+
                ('PD', 'PD')  :  ('Pd', 46, 106.42),
                ('AG', 'AG')  :  ('Ag', 47, 107.87),
                ('CD', 'CD')  :  ('Cd', 48, 112.41),
                ('IN', 'IN')  :  ('In', 49, 114.82),
                ('SN', 'SN')  :  ('Sn', 50, 118.7), #Add on
                ('SB', 'SB')  :  ('Sb', 51, 121.8), #Add on
                ('TE', 'TE')  :  ('Te', 52, 127.6), #Add on
                ('IOD', 'I')  :  ('I' , 53, 126.9),
                ('CS', 'CS')  :  ('Cs', 55, 132.91),
                ('BA', 'BA')  :  ('Ba', 56, 137.33),
                ('LA', 'LA')  :  ('La', 57, 138.91),
                ('CE', 'CE')  :  ('Ce', 58, 140.12),
                ('PR', 'PR')  :  ('Pr', 59, 140.91),
                ('ND', 'ND')  :  ('Nd', 60, 144.2), #Add on
                ('PM', 'PM')  :  ('Pm', 61, 145.0), #Add on
                ('SM', 'SM')  :  ('Sm', 62, 150.36),
                ('EU', 'EU')  :  ('Eu', 63, 151.96), #Eu2+
                ('EU3', 'EU') :  ('Eu', 63, 151.96), #Eu3+
                ('GD3', 'GD') :  ('Gd', 64, 157.25), #Gd3+
                ('GD', 'GD')  :  ('Gd', 64, 157.25),
                ('TB', 'TB')  :  ('Tb', 65, 158.93),
                ('DY', 'DY')  :  ('Dy', 66, 162.5), #Add on
                ('HO3', 'HO') :  ('Ho', 67, 164.93), #Ho3+
                ('HO', 'HO')  :  ('Ho', 67, 164.93),
                ('ER', 'ER')  :  ('Er', 68, 167.3), #Add on
                ('TM', 'TM')  :  ('Tm', 69, 168.9),
                ('YB2', 'YB2'):  ('Yb', 70, 173.05), #Yb2+
                ('YB', 'YB')  :  ('Yb', 70, 173.05), #Yb3+
                ('LU', 'LU')  :  ('Lu', 71, 174.97),
                ('HF', 'HF')  :  ('Hf', 72, 178.5), #Add on
                ('TA', 'TA')  :  ('Ta', 73, 180.9), #Add on
                ('W', 'W')    :  ('W',  74, 183.84),
                ('RE', 'RE')  :  ('Re', 75, 186.21),
                ('OS4', 'OS') :  ('Os', 76, 190.23), #Os4+
                ('OS', 'OS')  :  ('Os', 76, 190.23),
                ('IR3', 'IR') :  ('Ir', 77, 192.22), #Ir3+
                ('IR', 'IR')  :  ('Ir', 77, 192.22),
                ('PT', 'PT')  :  ('Pt', 78, 195.08),
                ('AU', 'AU')  :  ('Au', 79, 197.0),
                ('HG', 'HG')  :  ('Hg', 80, 200.59),
                ('TL', 'TL')  :  ('Tl', 81, 204.38),
                ('PB', 'PB')  :  ('Pb', 82, 207.2),
                ('BI', 'BI')  :  ('Bi', 83, 209.0), #Add on
                ('PO', 'PO')  :  ('Po', 84, 209.0), #Add on
                ('AT', 'AT')  :  ('At', 85, 210.0), #Add on
                ('FR', 'FR')  :  ('Fr', 87, 223.0), #Add on
                ('RA', 'RA')  :  ('Ra', 88, 226.0), #Add on
                ('AC', 'AC')  :  ('Ac', 89, 227.0), #Add on
                ('TH', 'TH')  :  ('Th', 90, 232.04),
                ('PA', 'PA')  :  ('Pa', 91, 231.0), #Add on
                ('U', 'U')    :  ('U',  92, 238.0), #Add on
                ('NP', 'NP')  :  ('Np', 93, 237.0), #Add on
                ('PU', 'PU')  :  ('Pu', 94, 244.0), #Add on
                ('AM', 'AM')  :  ('Am', 95, 243.0),
                ('CM', 'CM')  :  ('Cm', 96, 247.0), #Add on
                ('BK', 'BK')  :  ('Bk', 97, 247.0), #Add on
                ('CF', 'CF')  :  ('Cf', 98, 251.0), #Add on
                ('ES', 'ES')  :  ('Es', 99, 252.0), #Add on
                ('FM', 'FM')  :  ('Fm', 100, 257.0), #Add on
                ('MD', 'MD')  :  ('Md', 101, 258.0), #Add on
                ('NO', 'NO')  :  ('No', 102, 259.0), #Add on
                ('LR', 'LR')  :  ('Lr', 103, 262.0), #Add on
            }

Mass = OrganicMass
for i in METAL_PDB:
    if METAL_PDB[i][0] not in Mass:
        Mass[METAL_PDB[i][0]] = METAL_PDB[i][2]

ionnamel = []
for i in METAL_PDB:
    if METAL_PDB[i][0] not in ionnamel:
        ionnamel.append(METAL_PDB[i][0])
ionnamel2 = [i[0] + i[1].upper() for i in ionnamel if len(i) > 1]
ionnamel = ionnamel + ionnamel2

#-----------------------------------------------------------------------------

#Residue names and their letters

#-----------------------------------------------------------------------------

resdict = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'ASX': 'B',
           'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLX': 'Z', 'GLY': 'G',
           'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M',
           'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'SEC': 'U',
           'TRP': 'W', 'TYR': 'Y', 'XAA': 'X', 'VAL': 'V'}

#-----------------------------------------------------------------------------

#Residue names avaiable in the AMBER FFs

#-----------------------------------------------------------------------------

resnamel = ['ALA', 'ARG', 'ASH', 'ASN', 'ASP', 'CYM', 'CYS', 'CYX',
            'GLH', 'GLN', 'GLU', 'GLY', 'HIS', 'HID', 'HIE', 'HIP',
            'ILE', 'LEU', 'LYN', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
            'THR', 'TRP', 'TYR', 'VAL']

resnamecl = ['C' + i for i in resnamel]
resnamenl = ['N' + i for i in resnamel]
resnamel = resnamel + resnamecl + resnamenl

reschg = {'ALA': 0, 'ARG': 1, 'ASH': 0, 'ASN': 0, 'ASP': -1, 'CYM': -1,
          'CYS': 0, 'CYX': 0, 'GLH': 0, 'GLN': 0, 'GLU': -1, 'GLY':  0,
          'HIS': 0, 'HID': 0, 'HIE': 0, 'HIP': 1, 'ILE':  0, 'LEU':  0,
          'LYN': 0, 'LYS': 1, 'MET': 0, 'PHE': 0, 'PRO': 0,  'SER':  0,
          'THR': 0, 'TRP': 0, 'TYR': 0, 'VAL': 0}

reschgnt = dict([('N' + k, v+1) for k, v in list(reschg.items())])
reschgct = dict([('C' + k, v-1) for k, v in list(reschg.items())])

ResChgDict = reschg
ResChgDict.update(reschgnt)
ResChgDict.update(reschgct)

"""
Metal1pdb = {'F-':  ('F', 'F'), 'Br-': ('BR', 'BR'), 'Cl-': ('CL', 'CL'),
             'I-':  ('IOD', 'I'), 'Li+': ('LI', 'LI'), 'Na+': ('NA', 'NA'),
             'Rb+': ('RB', 'RB'), 'Tl+': ('TL', 'TL'), 'Cs+': ('CS', 'CS'),
             'K+':  ('K', 'K'), 'Cu+': ('CU1', 'CU'), 'Ag+': ('AG', 'AG'),
             'Au+': ('AU', 'AU')}

Metal2pdb = {'Cu2+': ('CU', 'CU'), 'Ni2+': ('NI', 'NI'), 'Pt2+': ('PT', 'PT'),
             'Zn2+': ('ZN', 'ZN'), 'Co2+': ('CO', 'CO'), 'Pd2+': ('PD', 'PD'),
             'Fe2+': ('FE2', 'FE2'), 'Mg2+': ('MG', 'MG'), 'Mn2+': ('MN', 'MN'),
             'Hg2+': ('HG', 'HG'), 'Cd2+': ('CD', 'CD'), 'Yb2+': ('YB2', 'YB2'),
             'Ca2+': ('CA', 'CA'), 'Pb2+': ('PD', 'PD'), 'Eu2+': ('EU', 'EU'),
             'Sr2+': ('SR', 'SR'), 'Ba2+': ('BA', 'BA')}

Metal3pdb = {'Al3+': ('AL', 'AL'), 'Fe3+': ('FE', 'FE'), 'Cr3+': ('CR', 'CR'),
             'In3+': ('IN', 'IN'), 'Y3+': ('Y', 'Y'), 'La3+': ('LA', 'LA'),
             'Ce3+': ('CE', 'CE'), 'Pr3+': ('PR', 'PR'), 'Sm3+': ('SM', 'SM'),
             'Eu3+': ('EU3', 'EU'), 'Gd3+': ('GD3', 'GD'), 'Tb3+': ('TB', 'TB'),
             'Lu3+': ('LU', 'LU'), 'V3+': ('V', 'V'), 'As3+': ('ARS', 'AS'),
             'Ru3+': ('RU', 'RU')}

Metal4pdb = {'Ir4+': ('IR', 'IR'), 'Mo4+': ('MO', 'MO')}
"""

#-----------------------------------------------------------------------------
#Covalent radii are used for determine the bonds between two atoms.
#The tolerance is 0.4 Angstrom. (<= Coradius1 + Coradius2 + 0.4)
#From Elaine C. Meng and Richard A. Lewis, Journal of Computational Chemistry,
#1991, 12(7), 891-898
#-----------------------------------------------------------------------------

CoRadiiDict = { 'H': 0.23, 'He': 1.50, 'Li': 0.68, 'Be': 0.35,  'B': 0.83,
                'C': 0.68,  'N': 0.68,  'O': 0.68,  'F': 0.64, 'Ne': 1.50,
               'Na': 0.97, 'Mg': 1.10, 'Al': 1.35, 'Si': 1.20,  'P': 1.05,
                'S': 1.02, 'Cl': 0.99, 'Ar': 1.51,  'K': 1.33, 'Ca': 0.99,
               'Sc': 1.44, 'Ti': 1.47,  'V': 1.33, 'Cr': 1.35, 'Mn': 1.35,
               'Fe': 1.34, 'Co': 1.33, 'Ni': 1.50, 'Cu': 1.52, 'Zn': 1.45,
               'Ga': 1.22, 'Ge': 1.17, 'As': 1.21, 'Se': 1.22, 'Br': 1.21,
               'Kr': 1.50, 'Rb': 1.47, 'Sr': 1.12,  'Y': 1.78, 'Zr': 1.56,
               'Nb': 1.48, 'Mo': 1.47, 'Tc': 1.35, 'Ru': 1.40, 'Rh': 1.45,
               'Pd': 1.50, 'Ag': 1.59, 'Cd': 1.69, 'In': 1.63, 'Sn': 1.46,
               'Sb': 1.46, 'Te': 1.47,  'I': 1.40, 'Cs': 1.67, 'Ba': 1.34,
               'La': 1.87, 'Ce': 1.83, 'Pr': 1.82, 'Nd': 1.81, 'Pm': 1.80,
               'Sm': 1.80, 'Eu': 1.99, 'Gd': 1.79, 'Tb': 1.76, 'Dy': 1.75,
               'Ho': 1.74, 'Er': 1.73, 'Tm': 1.72, 'Yb': 1.94, 'Lu': 1.72,
               'Hf': 1.57, 'Ta': 1.43,  'W': 1.37, 'Re': 1.35, 'Os': 1.37,
               'Ir': 1.32, 'Pt': 1.50, 'Au': 1.50, 'Hg': 1.70, 'Tl': 1.55,
               'Pb': 1.54, 'Bi': 1.54, 'Po': 1.68, 'Ra': 1.90, 'Ac': 1.88,
               'Th': 1.79, 'Pa': 1.61,  'U': 1.58, 'Np': 1.55, 'Pu': 1.53,
               'Am': 1.51} #91 elements

#-----------------------------------------------------------------------------

#VDW radii
#From http://periodictable.com/Properties/A/VanDerWaalsRadius.an.html

#-----------------------------------------------------------------------------

#Neutral
VdwRadiiDict = { 'H': 1.20, 'He': 1.40, 'Li': 1.82,  'C': 1.70,  'N': 1.55,
                 'O': 1.52,  'F': 1.47, 'Ne': 1.54, 'Na': 2.27, 'Mg': 1.73,
                'Si': 2.10,  'P': 1.80,  'S': 1.80, 'Cl': 1.75, 'Ar': 1.88,
                 'K': 2.75, 'Ni': 1.63, 'Cu': 1.40, 'Zn': 1.39, 'Ga': 1.87,
                'As': 1.85, 'Se': 1.90, 'Br': 1.85, 'Kr': 2.02, 'Pd': 1.63,
                'Ag': 1.72, 'Cd': 1.58, 'In': 1.93, 'Sn': 2.17, 'Te': 2.06,
                 'I': 1.98, 'Xe': 2.16, 'Pt': 1.75, 'Au': 1.66, 'Hg': 1.55,
                'Tl': 1.96, 'Pb': 2.02,  'U': 1.86} #38 elements

#LJ parameters for the ions in bonded model
def get_ionljparadict(watermodel):

    #Monovalent ions from IOD parameter set
    if watermodel in ['tip3p', 'spce', 'tip4pew']:
        monoljpara =  {
                    'Li1' : (1.315, 0.00594975, 'IOD set for Li+ ion from Li et al. JCTC, 2015, 11, 1645'),
                    'Na1' : (1.465, 0.02909167, 'IOD set for Na+ ion from Li et al. JCTC, 2015, 11, 1645'),
                    'K1'  : (1.745, 0.17018074, 'IOD set for K+ ion from Li et al. JCTC, 2015, 11, 1645'),
                    'Rb1' : (1.820, 0.22962229, 'IOD set for Rb+ ion from Li et al. JCTC, 2015, 11, 1645'),
                    'Cs1' : (2.000, 0.38943250, 'IOD set for Cs+ ion from Li et al. JCTC, 2015, 11, 1645'),
                    'Tl1' : (1.870, 0.27244486, 'IOD set for Tl+ ion from Li et al. JCTC, 2015, 11, 1645'),
                    'Cu1' : (1.214, 0.00139196, 'IOD set for Cu+ ion from Li et al. JCTC, 2015, 11, 1645'),
                    'Ag1' : (1.500, 0.03899838, 'IOD set for Ag+ ion from Li et al. JCTC, 2015, 11, 1645'),
                    'H1'  : (0.925, 0.00000147, 'IOD set for H+ ion from Li et al. JCTC, 2015, 11, 1645'),
                    'F-1' : (1.739, 0.16573832, 'IOD set for F- ion from Li et al. JCTC, 2015, 11, 1645'),
                    'Cl-1': (2.162, 0.53154665, 'IOD set for Cl- ion from Li et al. JCTC, 2015, 11, 1645'),
                    'Br-1': (2.331, 0.65952968, 'IOD set for Br- ion from Li et al. JCTC, 2015, 11, 1645'),
                    'I-1' : (2.590, 0.80293907, 'IOD set for I- ion from Li et al. JCTC, 2015, 11, 1645'),
                       }
    elif watermodel == 'opc3':
        monoljpara =  {
                    'Li1' : (1.321, 0.00641580, 'IOD set for Li+ ion for the OPC3 water model from Sengupta et al., JCIM, 2021, 61, 869'),
                    'Na1' : (1.470, 0.03038310, 'IOD set for Na+ ion for the OPC3 water model from Sengupta et al., JCIM, 2021, 61, 869'),
                    'K1 ' : (1.743, 0.16869420, 'IOD set for K+ ion for the OPC3 water model from Sengupta et al., JCIM, 2021, 61, 869'),
                    'Rb1' : (1.810, 0.22132374, 'IOD set for Rb+ ion for the OPC3 water model from Sengupta et al., JCIM, 2021, 61, 869'),
                    'Cs1' : (1.990, 0.38035199, 'IOD set for Cs+ ion for the OPC3 water model from Sengupta et al., JCIM, 2021, 61, 869'),
                    'Tl1' : (1.866, 0.26894857, 'IOD set for Tl+ ion for the OPC3 water model from Sengupta et al., JCIM, 2021, 61, 869'),
                    'Cu1' : (1.213, 0.00136949, 'IOD set for Cu+ ion for the OPC3 water model from Sengupta et al., JCIM, 2021, 61, 869'),
                    'Ag1' : (1.502, 0.03962711, 'IOD set for Ag+ ion for the OPC3 water model from Sengupta et al., JCIM, 2021, 61, 869'),
                    'F-1' : (1.737, 0.16426906, 'IOD set for F- ion for the OPC3 water model from Sengupta et al., JCIM, 2021, 61, 869'),
                    'Cl-1': (2.160, 0.52988504, 'IOD set for Cl- ion for the OPC3 water model from Sengupta et al., JCIM, 2021, 61, 869'),
                    'Br-1': (2.315, 0.64855145, 'IOD set for Br- ion for the OPC3 water model from Sengupta et al., JCIM, 2021, 61, 869'),
                    'I-1' : (2.590, 0.80293907, 'IOD set for I- ion for the OPC3 water model from Sengupta et al., JCIM, 2021, 61, 869'),
                       }
    elif watermodel == 'opc':
        monoljpara =  {
                    'Li1' : (1.305, 0.00523385, 'IOD set for Li+ ion for the OPC water model from Sengupta et al., JCIM, 2021, 61, 869'),
                    'Na1' : (1.440, 0.02322071, 'IOD set for Na+ ion for the OPC water model from Sengupta et al., JCIM, 2021, 61, 869'),
                    'K1 ' : (1.738, 0.16500296, 'IOD set for K+ ion for the OPC water model from Sengupta et al., JCIM, 2021, 61, 869'),
                    'Rb1' : (1.802, 0.21475916, 'IOD set for Rb+ ion for the OPC water model from Sengupta et al., JCIM, 2021, 61, 869'),
                    'Cs1' : (1.990, 0.38035199, 'IOD set for Cs+ ion for the OPC water model from Sengupta et al., JCIM, 2021, 61, 869'),
                    'Tl1' : (1.845, 0.25078000, 'IOD set for Tl+ ion for the OPC water model from Sengupta et al., JCIM, 2021, 61, 869'),
                    'Cu1' : (1.201, 0.00112300, 'IOD set for Cu+ ion for the OPC water model from Sengupta et al., JCIM, 2021, 61, 869'),
                    'Ag1' : (1.489, 0.03566355, 'IOD set for Ag+ ion for the OPC water model from Sengupta et al., JCIM, 2021, 61, 869'),
                    'F-1' : (1.720, 0.15202035, 'IOD set for F- ion for the OPC water model from Sengupta et al., JCIM, 2021, 61, 869'),
                    'Cl-1': (2.150, 0.52153239, 'IOD set for Cl- ion for the OPC water model from Sengupta et al., JCIM, 2021, 61, 869'),
                    'Br-1': (2.312, 0.64646527, 'IOD set for Br- ion for the OPC water model from Sengupta et al., JCIM, 2021, 61, 869'),
                    'I-1' : (2.573, 0.79541413, 'IOD set for I- ion for the OPC water model from Sengupta et al., JCIM, 2021, 61, 869'),
                       }
    elif watermodel == 'fb3':
        monoljpara =  {
                    'Li1' : (1.320, 0.00633615, 'IOD set for Li+ ion for the TIP3P-FB water model from Sengupta et al., JCIM, 2021, 61, 869'),
                    'Na1' : (1.460, 0.02784010, 'IOD set for Na+ ion for the TIP3P-FB water model from Sengupta et al., JCIM, 2021, 61, 869'),
                    'K1 ' : (1.741, 0.16721338, 'IOD set for K+ ion for the TIP3P-FB water model from Sengupta et al., JCIM, 2021, 61, 869'),
                    'Rb1' : (1.810, 0.22132374, 'IOD set for Rb+ ion for the TIP3P-FB water model from Sengupta et al., JCIM, 2021, 61, 869'),
                    'Cs1' : (1.988, 0.37853483, 'IOD set for Cs+ ion for the TIP3P-FB water model from Sengupta et al., JCIM, 2021, 61, 869'),
                    'Tl1' : (1.865, 0.26807617, 'IOD set for Tl+ ion for the TIP3P-FB water model from Sengupta et al., JCIM, 2021, 61, 869'),
                    'Cu1' : (1.214, 0.00139196, 'IOD set for Cu+ ion for the TIP3P-FB water model from Sengupta et al., JCIM, 2021, 61, 869'),
                    'Ag1' : (1.505, 0.04058327, 'IOD set for Ag+ ion for the TIP3P-FB water model from Sengupta et al., JCIM, 2021, 61, 869'),
                    'F-1' : (1.743, 0.16869420, 'IOD set for F- ion for the TIP3P-FB water model from Sengupta et al., JCIM, 2021, 61, 869'),
                    'Cl-1': (2.166, 0.53486081, 'IOD set for Cl- ion for the TIP3P-FB water model from Sengupta et al., JCIM, 2021, 61, 869'),
                    'Br-1': (2.330, 0.65885086, 'IOD set for Br- ion for the TIP3P-FB water model from Sengupta et al., JCIM, 2021, 61, 869'),
                    'I-1' : (2.597, 0.80596674, 'IOD set for I- ion for the TIP3P-FB water model from Sengupta et al., JCIM, 2021, 61, 869'),
                       }
    elif watermodel == 'fb4':
        monoljpara =  {
                    'Li1' : (1.306, 0.00530214, 'IOD set for Li+ ion for the TIP4P-FB water model from Sengupta et al., JCIM, 2021, 61, 869'),
                    'Na1' : (1.450, 0.02545423, 'IOD set for Na+ ion for the TIP4P-FB water model from Sengupta et al., JCIM, 2021, 61, 869'),
                    'K1 ' : (1.737, 0.16426906, 'IOD set for K+ ion for the TIP4P-FB water model from Sengupta et al., JCIM, 2021, 61, 869'),
                    'Rb1' : (1.810, 0.22132374, 'IOD set for Rb+ ion for the TIP4P-FB water model from Sengupta et al., JCIM, 2021, 61, 869'),
                    'Cs1' : (2.000, 0.38943250, 'IOD set for Cs+ ion for the TIP4P-FB water model from Sengupta et al., JCIM, 2021, 61, 869'),
                    'Tl1' : (1.860, 0.26372453, 'IOD set for Tl+ ion for the TIP4P-FB water model from Sengupta et al., JCIM, 2021, 61, 869'),
                    'Cu1' : (1.208, 0.00126172, 'IOD set for Cu+ ion for the TIP4P-FB water model from Sengupta et al., JCIM, 2021, 61, 869'),
                    'Ag1' : (1.481, 0.03336723, 'IOD set for Ag+ ion for the TIP4P-FB water model from Sengupta et al., JCIM, 2021, 61, 869'),
                    'F-1' : (1.740, 0.16647513, 'IOD set for F- ion for the TIP4P-FB water model from Sengupta et al., JCIM, 2021, 61, 869'),
                    'Cl-1': (2.166, 0.53486081, 'IOD set for Cl- ion for the TIP4P-FB water model from Sengupta et al., JCIM, 2021, 61, 869'),
                    'Br-1': (2.340, 0.66559495, 'IOD set for Br- ion for the TIP4P-FB water model from Sengupta et al., JCIM, 2021, 61, 869'),
                    'I-1' : (2.588, 0.80206648, 'IOD set for I- ion for the TIP4P-FB water model from Sengupta et al., JCIM, 2021, 61, 869'),
                       }

    IonLJParaDict = monoljpara

    #Divalent ions from IOD parameter set
    if watermodel in ['tip3p', 'spce', 'tip4pew']:
        diljpara = {
                  'Be2': (1.168, 0.00063064, 'IOD set for Be2+ ion from Li et al. JCTC, 2013, 9, 2733'),
                  'Cu2': (1.409, 0.01721000, 'IOD set for Cu2+ ion from Li et al. JCTC, 2013, 9, 2733'),
                  'Ni2': (1.373, 0.01179373, 'IOD set for Ni2+ ion from Li et al. JCTC, 2013, 9, 2733'),
                  'Zn2': (1.395, 0.01491700, 'IOD set for Zn2+ ion from Li et al. JCTC, 2013, 9, 2733'),
                  'Co2': (1.404, 0.01636246, 'IOD set for Co2+ ion from Li et al. JCTC, 2013, 9, 2733'),
                  'Cr2': (1.388, 0.01386171, 'IOD set for Cr2+ ion from Li et al. JCTC, 2013, 9, 2733'),
                  'Fe2': (1.409, 0.01721000, 'IOD set for Fe2+ ion from Li et al. JCTC, 2013, 9, 2733'),
                  'Mg2': (1.395, 0.01491700, 'IOD set for Mg2+ ion from Li et al. JCTC, 2013, 9, 2733'),
                  'V2':  (1.476, 0.03198620, 'IOD set for V2+ ion from Li et al. JCTC, 2013, 9, 2733'),
                  'Mn2': (1.467, 0.02960343, 'IOD set for Mn2+ ion from Li et al. JCTC, 2013, 9, 2733'),
                  'Hg2': (1.575, 0.06751391, 'IOD set for Hg2+ ion from Li et al. JCTC, 2013, 9, 2733'),
                  'Cd2': (1.506, 0.04090549, 'IOD set for Cd2+ ion from Li et al. JCTC, 2013, 9, 2733'),
                  'Ca2': (1.608, 0.08337961, 'IOD set for Ca2+ ion from Li et al. JCTC, 2013, 9, 2733'),
                  'Sn2': (1.738, 0.16500296, 'IOD set for Sn2+ ion from Li et al. JCTC, 2013, 9, 2733'),
                  'Sr2': (1.753, 0.17618319, 'IOD set for Sr2+ ion from Li et al. JCTC, 2013, 9, 2733'),
                  'Ba2': (1.913, 0.31060194, 'IOD set for Ba2+ ion from Li et al. JCTC, 2013, 9, 2733'),
                    }
    elif watermodel == 'opc3':
         diljpara = {
                  'Be2': (1.162, 0.00056491, 'IOD set for Be2+ ion for the OPC3 water model from Li et al. JCTC, 2020, 16, 4429'),
                  'Cu2': (1.413, 0.01791152, 'IOD set for Cu2+ ion for the OPC3 water model from Li et al. JCTC, 2020, 16, 4429'),
                  'Ni2': (1.373, 0.01179373, 'IOD set for Ni2+ ion for the OPC3 water model from Li et al. JCTC, 2020, 16, 4429'),
                  'Zn2': (1.400, 0.01570749, 'IOD set for Zn2+ ion for the OPC3 water model from Li et al. JCTC, 2020, 16, 4429'),
                  'Co2': (1.406, 0.01669760, 'IOD set for Co2+ ion for the OPC3 water model from Li et al. JCTC, 2020, 16, 4429'),
                  'Cr2': (1.391, 0.01430674, 'IOD set for Cr2+ ion for the OPC3 water model from Li et al. JCTC, 2020, 16, 4429'),
                  'Fe2': (1.413, 0.01430674, 'IOD set for Fe2+ ion for the OPC3 water model from Li et al. JCTC, 2020, 16, 4429'),
                  'Mg2': (1.400, 0.01570749, 'IOD set for Mg2+ ion for the OPC3 water model from Li et al. JCTC, 2020, 16, 4429'),
                  'V2':  (1.480, 0.03308772, 'IOD set for V2+ ion for the OPC3 water model from Li et al. JCTC, 2020, 16, 4429'),
                  'Mn2': (1.462, 0.02833599, 'IOD set for Mn2+ ion for the OPC3 water model from Li et al. JCTC, 2020, 16, 4429'),
                  'Hg2': (1.584, 0.07163727, 'IOD set for Hg2+ ion for the OPC3 water model from Li et al. JCTC, 2020, 16, 4429'),
                  'Cd2': (1.526, 0.04772212, 'IOD set for Cd2+ ion for the OPC3 water model from Li et al. JCTC, 2020, 16, 4429'),
                  'Ca2': (1.617, 0.08806221, 'IOD set for Ca2+ ion for the OPC3 water model from Li et al. JCTC, 2020, 16, 4429'),
                  'Sn2': (1.746, 0.17092614, 'IOD set for Sn2+ ion for the OPC3 water model from Li et al. JCTC, 2020, 16, 4429'),
                  'Sr2': (1.762, 0.18304100, 'IOD set for Sr2+ ion for the OPC3 water model from Li et al. JCTC, 2020, 16, 4429'),
                  'Ba2': (1.918, 0.31509345, 'IOD set for Ba2+ ion for the OPC3 water model from Li et al. JCTC, 2020, 16, 4429'),
                    }
    elif watermodel == 'opc':
         diljpara = {
                  'Be2': (1.136, 0.00034392, 'IOD set for Be2+ ion for the OPC water model from Li et al. JCTC, 2020, 16, 4429'),
                  'Cu2': (1.391, 0.01430674, 'IOD set for Cu2+ ion for the OPC water model from Li et al. JCTC, 2020, 16, 4429'),
                  'Ni2': (1.345, 0.00858042, 'IOD set for Ni2+ ion for the OPC water model from Li et al. JCTC, 2020, 16, 4429'),
                  'Zn2': (1.373, 0.01179373, 'IOD set for Zn2+ ion for the OPC water model from Li et al. JCTC, 2020, 16, 4429'),
                  'Co2': (1.382, 0.01300356, 'IOD set for Co2+ ion for the OPC water model from Li et al. JCTC, 2020, 16, 4429'),
                  'Cr2': (1.364, 0.01067299, 'IOD set for Cr2+ ion for the OPC water model from Li et al. JCTC, 2020, 16, 4429'),
                  'Fe2': (1.391, 0.01430674, 'IOD set for Fe2+ ion for the OPC water model from Li et al. JCTC, 2020, 16, 4429'),
                  'Mg2': (1.373, 0.01179373, 'IOD set for Mg2+ ion for the OPC water model from Li et al. JCTC, 2020, 16, 4429'),
                  'V2':  (1.456, 0.02686716, 'IOD set for V2+ ion for the OPC water model from Li et al. JCTC, 2020, 16, 4429'),
                  'Mn2': (1.444, 0.02409615, 'IOD set for Mn2+ ion for the OPC water model from Li et al. JCTC, 2020, 16, 4429'),
                  'Hg2': (1.565, 0.06311131, 'IOD set for Hg2+ ion for the OPC water model from Li et al. JCTC, 2020, 16, 4429'),
                  'Cd2': (1.506, 0.04090549, 'IOD set for Cd2+ ion for the OPC water model from Li et al. JCTC, 2020, 16, 4429'),
                  'Ca2': (1.590, 0.07447106, 'IOD set for Ca2+ ion for the OPC water model from Li et al. JCTC, 2020, 16, 4429'),
                  'Sn2': (1.715, 0.14850170, 'IOD set for Sn2+ ion for the OPC water model from Li et al. JCTC, 2020, 16, 4429'),
                  'Sr2': (1.731, 0.15989650, 'IOD set for Sr2+ ion for the OPC water model from Li et al. JCTC, 2020, 16, 4429'),
                  'Ba2': (1.883, 0.28387745, 'IOD set for Ba2+ ion for the OPC water model from Li et al. JCTC, 2020, 16, 4429'),
                    }
    elif watermodel == 'fb3':
         diljpara = {
                  'Be2': (1.163, 0.00057544, 'IOD set for Be2+ ion for the TIP3P-FB water model from Li et al. JCTC, 2020, 16, 4429'),
                  'Cu2': (1.413, 0.01791152, 'IOD set for Cu2+ ion for the TIP3P-FB water model from Li et al. JCTC, 2020, 16, 4429'),
                  'Ni2': (1.373, 0.01179373, 'IOD set for Ni2+ ion for the TIP3P-FB water model from Li et al. JCTC, 2020, 16, 4429'),
                  'Zn2': (1.400, 0.01570749, 'IOD set for Zn2+ ion for the TIP3P-FB water model from Li et al. JCTC, 2020, 16, 4429'),
                  'Co2': (1.406, 0.01669760, 'IOD set for Co2+ ion for the TIP3P-FB water model from Li et al. JCTC, 2020, 16, 4429'),
                  'Cr2': (1.391, 0.01430674, 'IOD set for Cr2+ ion for the TIP3P-FB water model from Li et al. JCTC, 2020, 16, 4429'),
                  'Fe2': (1.413, 0.01791152, 'IOD set for Fe2+ ion for the TIP3P-FB water model from Li et al. JCTC, 2020, 16, 4429'),
                  'Mg2': (1.400, 0.01570749, 'IOD set for Mg2+ ion for the TIP3P-FB water model from Li et al. JCTC, 2020, 16, 4429'),
                  'V2':  (1.480, 0.03308772, 'IOD set for V2+ ion for the TIP3P-FB water model from Li et al. JCTC, 2020, 16, 4429'),
                  'Mn2': (1.462, 0.02833599, 'IOD set for Mn2+ ion for the TIP3P-FB water model from Li et al. JCTC, 2020, 16, 4429'),
                  'Hg2': (1.584, 0.07163727, 'IOD set for Hg2+ ion for the TIP3P-FB water model from Li et al. JCTC, 2020, 16, 4429'),
                  'Cd2': (1.520, 0.04560206, 'IOD set for Cd2+ ion for the TIP3P-FB water model from Li et al. JCTC, 2020, 16, 4429'),
                  'Ca2': (1.617, 0.08806221, 'IOD set for Ca2+ ion for the TIP3P-FB water model from Li et al. JCTC, 2020, 16, 4429'),
                  'Sn2': (1.746, 0.17092614, 'IOD set for Sn2+ ion for the TIP3P-FB water model from Li et al. JCTC, 2020, 16, 4429'),
                  'Sr2': (1.762, 0.18304100, 'IOD set for Sr2+ ion for the TIP3P-FB water model from Li et al. JCTC, 2020, 16, 4429'),
                  'Ba2': (1.918, 0.31509345, 'IOD set for Ba2+ ion for the TIP3P-FB water model from Li et al. JCTC, 2020, 16, 4429'),
                    }
    elif watermodel == 'fb4':
         diljpara = {
                  'Be2': (1.150, 0.00045105, 'IOD set for Be2+ ion for the TIP4P-FB water model from Li et al. JCTC, 2020, 16, 4429'),
                  'Cu2': (1.400, 0.01570749, 'IOD set for Cu2+ ion for the TIP4P-FB water model from Li et al. JCTC, 2020, 16, 4429'),
                  'Ni2': (1.358, 0.00997323, 'IOD set for Ni2+ ion for the TIP4P-FB water model from Li et al. JCTC, 2020, 16, 4429'),
                  'Zn2': (1.383, 0.01314367, 'IOD set for Zn2+ ion for the TIP4P-FB water model from Li et al. JCTC, 2020, 16, 4429'),
                  'Co2': (1.392, 0.01445748, 'IOD set for Co2+ ion for the TIP4P-FB water model from Li et al. JCTC, 2020, 16, 4429'),
                  'Cr2': (1.375, 0.01205473, 'IOD set for Cr2+ ion for the TIP4P-FB water model from Li et al. JCTC, 2020, 16, 4429'),
                  'Fe2': (1.400, 0.01570749, 'IOD set for Fe2+ ion for the TIP4P-FB water model from Li et al. JCTC, 2020, 16, 4429'),
                  'Mg2': (1.383, 0.01314367, 'IOD set for Mg2+ ion for the TIP4P-FB water model from Li et al. JCTC, 2020, 16, 4429'),
                  'V2':  (1.465, 0.02909167, 'IOD set for V2+ ion for the TIP4P-FB water model from Li et al. JCTC, 2020, 16, 4429'),
                  'Mn2': (1.459, 0.02759452, 'IOD set for Mn2+ ion for the TIP4P-FB water model from Li et al. JCTC, 2020, 16, 4429'),
                  'Hg2': (1.572, 0.06617338, 'IOD set for Hg2+ ion for the TIP4P-FB water model from Li et al. JCTC, 2020, 16, 4429'),
                  'Cd2': (1.511, 0.04254294, 'IOD set for Cd2+ ion for the TIP4P-FB water model from Li et al. JCTC, 2020, 16, 4429'),
                  'Ca2': (1.600, 0.07934493, 'IOD set for Ca2+ ion for the TIP4P-FB water model from Li et al. JCTC, 2020, 16, 4429'),
                  'Sn2': (1.731, 0.15989650, 'IOD set for Sn2+ ion for the TIP4P-FB water model from Li et al. JCTC, 2020, 16, 4429'),
                  'Sr2': (1.746, 0.17092614, 'IOD set for Sr2+ ion for the TIP4P-FB water model from Li et al. JCTC, 2020, 16, 4429'),
                  'Ba2': (1.900, 0.29896986, 'IOD set for Ba2+ ion for the TIP4P-FB water model from Li et al. JCTC, 2020, 16, 4429'),
                    }

    IonLJParaDict.update(diljpara)

    #Divalent ions from CM parameter set
    if watermodel == 'tip3p':
        diljpara2 = {
                   'Pt2': (1.266, 0.00307642, 'CM set for Pt2+ ion in TIP3P water from Li et al. JCTC, 2013, 9, 2733'),
                   'Pd2': (1.303, 0.00509941, 'CM set for Pd2+ ion in TIP3P water from Li et al. JCTC, 2013, 9, 2733'),
                   'Ag2': (1.336, 0.00770969, 'CM set for Ag2+ ion in TIP3P water from Li et al. JCTC, 2013, 9, 2733'),
                   'Yb2': (1.642, 0.10185975, 'CM set for Yb2+ ion in TIP3P water from Li et al. JCTC, 2013, 9, 2733'),
                   'Pb2': (1.745, 0.17018074, 'CM set for Pb2+ ion in TIP3P water from Li et al. JCTC, 2013, 9, 2733'),
                   'Eu2': (1.802, 0.21475916, 'CM set for Eu2+ ion in TIP3P water from Li et al. JCTC, 2013, 9, 2733'),
                   'Sm2': (1.819, 0.22878796, 'CM set for Sm2+ ion in TIP3P water from Li et al. JCTC, 2013, 9, 2733'),
                   'Ra2': (2.019, 0.40664608, 'CM set for Ra2+ ion in TIP3P water from Li et al. JCTC, 2013, 9, 2733'),
                     }
    elif watermodel == 'spce':
        diljpara2 = {
                   'Pt2': (1.272, 0.00334975, 'CM set for Pt2+ ion in SPC/E water from Li et al. JCTC, 2013, 9, 2733'),
                   'Pd2': (1.305, 0.00523385, 'CM set for Pd2+ ion in SPC/E water from Li et al. JCTC, 2013, 9, 2733'),
                   'Ag2': (1.337, 0.00780282, 'CM set for Ag2+ ion in SPC/E water from Li et al. JCTC, 2013, 9, 2733'),
                   'Yb2': (1.634, 0.09731901, 'CM set for Yb2+ ion in SPC/E water from Li et al. JCTC, 2013, 9, 2733'),
                   'Pb2': (1.731, 0.15989650, 'CM set for Pb2+ ion in SPC/E water from Li et al. JCTC, 2013, 9, 2733'),
                   'Eu2': (1.786, 0.20184160, 'CM set for Eu2+ ion in SPC/E water from Li et al. JCTC, 2013, 9, 2733'),
                   'Sm2': (1.800, 0.21312875, 'CM set for Sm2+ ion in SPC/E water from Li et al. JCTC, 2013, 9, 2733'),
                   'Ra2': (1.980, 0.37126402, 'CM set for Ra2+ ion in SPC/E water from Li et al. JCTC, 2013, 9, 2733'),
                     }
    elif watermodel == 'tip4pew':
        diljpara2 = {
                   'Pt2': (1.251, 0.00247282, 'CM set for Pt2+ ion in TIP4P/EW water from Li et al. JCTC, 2013, 9, 2733'),
                   'Pd2': (1.288, 0.00417787, 'CM set for Pd2+ ion in TIP4P/EW water from Li et al. JCTC, 2013, 9, 2733'),
                   'Ag2': (1.323, 0.00657749, 'CM set for Ag2+ ion in TIP4P/EW water from Li et al. JCTC, 2013, 9, 2733'),
                   'Yb2': (1.654, 0.10888937, 'CM set for Yb2+ ion in TIP4P/EW water from Li et al. JCTC, 2013, 9, 2733'),
                   'Pb2': (1.758, 0.17997960, 'CM set for Pb2+ ion in TIP4P/EW water from Li et al. JCTC, 2013, 9, 2733'),
                   'Eu2': (1.823, 0.23213110, 'CM set for Eu2+ ion in TIP4P/EW water from Li et al. JCTC, 2013, 9, 2733'),
                   'Sm2': (1.838, 0.24480038, 'CM set for Sm2+ ion in TIP4P/EW water from Li et al. JCTC, 2013, 9, 2733'),
                   'Ra2': (2.050, 0.43454345, 'CM set for Ra2+ ion in TIP4P/EW water from Li et al. JCTC, 2013, 9, 2733'),
                     }
    elif watermodel == 'opc3':
        diljpara2 = {
                   'Pt2': (1.272, 0.00334975, 'CM set for Pt2+ ion in OPC3 water from Li et al. JCTC, 2020, 16, 4429'),
                   'Pd2': (1.309, 0.00551135, 'CM set for Pd2+ ion in OPC3 water from Li et al. JCTC, 2020, 16, 4429'),
                   'Ag2': (1.339, 0.00799176, 'CM set for Ag2+ ion in OPC3 water from Li et al. JCTC, 2020, 16, 4429'),
                   'Yb2': (1.622, 0.09072908, 'CM set for Yb2+ ion in OPC3 water from Li et al. JCTC, 2020, 16, 4429'),
                   'Pb2': (1.723, 0.15415012, 'CM set for Pb2+ ion in OPC3 water from Li et al. JCTC, 2020, 16, 4429'),
                   'Eu2': (1.772, 0.19078645, 'CM set for Eu2+ ion in OPC3 water from Li et al. JCTC, 2020, 16, 4429'),
                   'Sm2': (1.784, 0.20024770, 'CM set for Sm2+ ion in OPC3 water from Li et al. JCTC, 2020, 16, 4429'),
                   'Ra2': (1.966, 0.35853865, 'CM set for Ra2+ ion in OPC3 water from Li et al. JCTC, 2020, 16, 4429'),
                     }
    elif watermodel == 'opc':
        diljpara2 = {
                   'Pt2': (1.219, 0.00150903, 'CM set for Pt2+ ion in OPC water from Li et al. JCTC, 2020, 16, 4429'),
                   'Pd2': (1.269, 0.00321068, 'CM set for Pd2+ ion in OPC water from Li et al. JCTC, 2020, 16, 4429'),
                   'Ag2': (1.305, 0.00523385, 'CM set for Ag2+ ion in OPC water from Li et al. JCTC, 2020, 16, 4429'),
                   'Yb2': (1.602, 0.08034231, 'CM set for Yb2+ ion in OPC water from Li et al. JCTC, 2020, 16, 4429'),
                   'Pb2': (1.707, 0.14295367, 'CM set for Pb2+ ion in OPC water from Li et al. JCTC, 2020, 16, 4429'),
                   'Eu2': (1.753, 0.17618319, 'CM set for Eu2+ ion in OPC water from Li et al. JCTC, 2020, 16, 4429'),
                   'Sm2': (1.766, 0.18612361, 'CM set for Sm2+ ion in OPC water from Li et al. JCTC, 2020, 16, 4429'),
                   'Ra2': (1.960, 0.35308749, 'CM set for Ra2+ ion in OPC water from Li et al. JCTC, 2020, 16, 4429'),
                     }
    elif watermodel == 'fb3':
        diljpara2 = {
                   'Pt2': (1.274, 0.00344520, 'CM set for Pt2+ ion in TIP3P-FB water from Li et al. JCTC, 2020, 16, 4429'),
                   'Pd2': (1.308, 0.00544088, 'CM set for Pd2+ ion in TIP3P-FB water from Li et al. JCTC, 2020, 16, 4429'),
                   'Ag2': (1.339, 0.00799176, 'CM set for Ag2+ ion in TIP3P-FB water from Li et al. JCTC, 2020, 16, 4429'),
                   'Yb2': (1.629, 0.09454081, 'CM set for Yb2+ ion in TIP3P-FB water from Li et al. JCTC, 2020, 16, 4429'),
                   'Pb2': (1.730, 0.15917293, 'CM set for Pb2+ ion in TIP3P-FB water from Li et al. JCTC, 2020, 16, 4429'),
                   'Eu2': (1.782, 0.19865859, 'CM set for Eu2+ ion in TIP3P-FB water from Li et al. JCTC, 2020, 16, 4429'),
                   'Sm2': (1.795, 0.20907204, 'CM set for Sm2+ ion in TIP3P-FB water from Li et al. JCTC, 2020, 16, 4429'),
                   'Ra2': (1.983, 0.37399087, 'CM set for Ra2+ ion in TIP3P-FB water from Li et al. JCTC, 2020, 16, 4429'),
                     }
    elif watermodel == 'fb4':
        diljpara2 = {
                   'Pt2': (1.229, 0.00176831, 'CM set for Pt2+ ion in TIP4P-FB water from Li et al. JCTC, 2020, 16, 4429'),
                   'Pd2': (1.282, 0.00384964, 'CM set for Pd2+ ion in TIP4P-FB water from Li et al. JCTC, 2020, 16, 4429'),
                   'Ag2': (1.316, 0.00602547, 'CM set for Ag2+ ion in TIP4P-FB water from Li et al. JCTC, 2020, 16, 4429'),
                   'Yb2': (1.621, 0.09019198, 'CM set for Yb2+ ion in TIP4P-FB water from Li et al. JCTC, 2020, 16, 4429'),
                   'Pb2': (1.728, 0.15773029, 'CM set for Pb2+ ion in TIP4P-FB water from Li et al. JCTC, 2020, 16, 4429'),
                   'Eu2': (1.782, 0.19865859, 'CM set for Eu2+ ion in TIP4P-FB water from Li et al. JCTC, 2020, 16, 4429'),
                   'Sm2': (1.796, 0.20988115, 'CM set for Sm2+ ion in TIP4P-FB water from Li et al. JCTC, 2020, 16, 4429'),
                   'Ra2': (1.992, 0.38216886, 'CM set for Ra2+ ion in TIP4P-FB water from Li et al. JCTC, 2020, 16, 4429'),
                     }

    IonLJParaDict.update(diljpara2)

    #Trivalent and tetravalent ions from IOD parameter set
    if watermodel == 'tip3p':
        higljpara = {
                   'Al3': (1.297, 0.00471279, 'IOD set for Al3+ ion in TIP3P water from Li et al. JPCB, 2015, 119, 883'),
                   'Fe3': (1.386, 0.01357097, 'IOD set for Fe3+ ion in TIP3P water from Li et al. JPCB, 2015, 119, 883'),
                   'Cr3': (1.344, 0.00848000, 'IOD set for Cr3+ ion in TIP3P water from Li et al. JPCB, 2015, 119, 883'),
                   'In3': (1.461, 0.02808726, 'IOD set for In3+ ion in TIP3P water from Li et al. JPCB, 2015, 119, 883'),
                   'Tl3': (1.513, 0.04321029, 'IOD set for Tl3+ ion in TIP3P water from Li et al. JPCB, 2015, 119, 883'),
                   'Y3':  (1.602, 0.08034231, 'IOD set for Y3+ ion in TIP3P water from Li et al. JPCB, 2015, 119, 883'),
                   'La3': (1.718, 0.15060822, 'IOD set for La3+ ion in TIP3P water from Li et al. JPCB, 2015, 119, 883'),
                   'Ce3': (1.741, 0.16721338, 'IOD set for Ce3+ ion in TIP3P water from Li et al. JPCB, 2015, 119, 883'),
                   'Pr3': (1.733, 0.16134811, 'IOD set for Pr3+ ion in TIP3P water from Li et al. JPCB, 2015, 119, 883'),
                   'Nd3': (1.681, 0.12564307, 'IOD set for Nd3+ ion in TIP3P water from Li et al. JPCB, 2015, 119, 883'),
                   'Sm3': (1.659, 0.11189491, 'IOD set for Sm3+ ion in TIP3P water from Li et al. JPCB, 2015, 119, 883'),
                   'Eu3': (1.666, 0.11617738, 'IOD set for Eu3+ ion in TIP3P water from Li et al. JPCB, 2015, 119, 883'),
                   'Gd3': (1.623, 0.09126804, 'IOD set for Gd3+ ion in TIP3P water from Li et al. JPCB, 2015, 119, 883'),
                   'Tb3': (1.630, 0.09509276, 'IOD set for Tb3+ ion in TIP3P water from Li et al. JPCB, 2015, 119, 883'),
                   'Dy3': (1.609, 0.08389240, 'IOD set for Dy3+ ion in TIP3P water from Li et al. JPCB, 2015, 119, 883'),
                   'Er3': (1.602, 0.08034231, 'IOD set for Er3+ ion in TIP3P water from Li et al. JPCB, 2015, 119, 883'),
                   'Tm3': (1.602, 0.08034231, 'IOD set for Tm3+ ion in TIP3P water from Li et al. JPCB, 2015, 119, 883'),
                   'Lu3': (1.588, 0.07351892, 'IOD set for Lu3+ ion in TIP3P water from Li et al. JPCB, 2015, 119, 883'),
                   'Hf4': (1.499, 0.03868661, 'IOD set for Hf4+ ion in TIP3P water from Li et al. JPCB, 2015, 119, 883'),
                   'Zr4': (1.519, 0.04525501, 'IOD set for Zr4+ ion in TIP3P water from Li et al. JPCB, 2015, 119, 883'),
                   'Ce4': (1.684, 0.12758274, 'IOD set for Ce4+ ion in TIP3P water from Li et al. JPCB, 2015, 119, 883'),
                   'U4':  (1.684, 0.12758274, 'IOD set for U4+ ion in TIP3P water from Li et al. JPCB, 2015, 119, 883'),
                   'Pu4': (1.662, 0.11371963, 'IOD set for Pu4+ ion in TIP3P water from Li et al. JPCB, 2015, 119, 883'),
                   'Th4': (1.708, 0.14364160, 'IOD set for Th4+ ion TIP3P water from Li et al. JPCB, 2015, 119, 883'),
                     }
    elif watermodel == 'spce':
        higljpara = {
                   'Al3': (1.296, 0.00465074, 'IOD set for Al3+ ion in SPC/E water from Li et al. JPCB, 2015, 119, 883'),
                   'Fe3': (1.386, 0.01357097, 'IOD set for Fe3+ ion in SPC/E water from Li et al. JPCB, 2015, 119, 883'),
                   'Cr3': (1.343, 0.00838052, 'IOD set for Cr3+ ion in SPC/E water from Li et al. JPCB, 2015, 119, 883'),
                   'In3': (1.461, 0.02808726, 'IOD set for In3+ ion in SPC/E water from Li et al. JPCB, 2015, 119, 883'),
                   'Tl3': (1.513, 0.04321029, 'IOD set for Tl3+ ion in SPC/E water from Li et al. JPCB, 2015, 119, 883'),
                   'Y3':  (1.602, 0.08034231, 'IOD set for Y3+ ion in SPC/E water from Li et al. JPCB, 2015, 119, 883'),
                   'La3': (1.718, 0.15060822, 'IOD set for La3+ ion in SPC/E water from Li et al. JPCB, 2015, 119, 883'),
                   'Ce3': (1.741, 0.16721338, 'IOD set for Ce3+ ion in SPC/E water from Li et al. JPCB, 2015, 119, 883'),
                   'Pr3': (1.734, 0.16207614, 'IOD set for Pr3+ ion in SPC/E water from Li et al. JPCB, 2015, 119, 883'),
                   'Nd3': (1.681, 0.12564307, 'IOD set for Nd3+ ion in SPC/E water from Li et al. JPCB, 2015, 119, 883'),
                   'Sm3': (1.659, 0.11189491, 'IOD set for Sm3+ ion in SPC/E water from Li et al. JPCB, 2015, 119, 883'),
                   'Eu3': (1.666, 0.11617738, 'IOD set for Eu3+ ion in SPC/E water from Li et al. JPCB, 2015, 119, 883'),
                   'Gd3': (1.623, 0.09126804, 'IOD set for Gd3+ ion in SPC/E water from Li et al. JPCB, 2015, 119, 883'),
                   'Tb3': (1.630, 0.09509276, 'IOD set for Tb3+ ion in SPC/E water from Li et al. JPCB, 2015, 119, 883'),
                   'Dy3': (1.609, 0.08389240, 'IOD set for Dy3+ ion in SPC/E water from Li et al. JPCB, 2015, 119, 883'),
                   'Er3': (1.602, 0.08034231, 'IOD set for Er3+ ion in SPC/E water from Li et al. JPCB, 2015, 119, 883'),
                   'Tm3': (1.602, 0.08034231, 'IOD set for Tm3+ ion in SPC/E water from Li et al. JPCB, 2015, 119, 883'),
                   'Lu3': (1.588, 0.07351892, 'IOD set for Lu3+ ion in SPC/E water from Li et al. JPCB, 2015, 119, 883'),
                   'Hf4': (1.501, 0.03931188, 'IOD set for Hf4+ ion in SPC/E water from Li et al. JPCB, 2015, 119, 883'),
                   'Zr4': (1.521, 0.04595090, 'IOD set for Zr4+ ion in SPC/E water from Li et al. JPCB, 2015, 119, 883'),
                   'Ce4': (1.689, 0.13084945, 'IOD set for Ce4+ ion in SPC/E water from Li et al. JPCB, 2015, 119, 883'),
                   'U4':  (1.689, 0.13084945, 'IOD set for U4+ ion in SPC/E water from Li et al. JPCB, 2015, 119, 883'),
                   'Pu4': (1.666, 0.11617738, 'IOD set for Pu4+ ion in SPC/E water from Li et al. JPCB, 2015, 119, 883'),
                   'Th4': (1.713, 0.14710519, 'IOD set for Th4+ ion in SPC/E water from Li et al. JPCB, 2015, 119, 883'),
                     }
    elif watermodel == 'tip4pew':
        higljpara = {
                   'Al3': (1.285, 0.00401101, 'IOD set for Al3+ ion TIP4P/EW water from Li et al. JPCB, 2015, 119, 883'),
                   'Fe3': (1.375, 0.01205473, 'IOD set for Fe3+ ion TIP4P/EW water from Li et al. JPCB, 2015, 119, 883'),
                   'Cr3': (1.333, 0.00743559, 'IOD set for Cr3+ ion TIP4P/EW water from Li et al. JPCB, 2015, 119, 883'),
                   'In3': (1.450, 0.02545423, 'IOD set for In3+ ion TIP4P/EW water from Li et al. JPCB, 2015, 119, 883'),
                   'Tl3': (1.502, 0.03962711, 'IOD set for Tl3+ ion TIP4P/EW water from Li et al. JPCB, 2015, 119, 883'),
                   'Y3':  (1.590, 0.07447106, 'IOD set for Y3+ ion TIP4P/EW water from Li et al. JPCB, 2015, 119, 883'),
                   'La3': (1.707, 0.14295367, 'IOD set for La3+ ion TIP4P/EW water from Li et al. JPCB, 2015, 119, 883'),
                   'Ce3': (1.729, 0.15845086, 'IOD set for Ce3+ ion TIP4P/EW water from Li et al. JPCB, 2015, 119, 883'),
                   'Pr3': (1.722, 0.15343866, 'IOD set for Pr3+ ion TIP4P/EW water from Li et al. JPCB, 2015, 119, 883'),
                   'Nd3': (1.669, 0.11803919, 'IOD set for Nd3+ ion TIP4P/EW water from Li et al. JPCB, 2015, 119, 883'),
                   'Sm3': (1.647, 0.10475707, 'IOD set for Sm3+ ion TIP4P/EW water from Li et al. JPCB, 2015, 119, 883'),
                   'Eu3': (1.655, 0.10948690, 'IOD set for Eu3+ ion TIP4P/EW water from Li et al. JPCB, 2015, 119, 883'),
                   'Gd3': (1.612, 0.08544204, 'IOD set for Gd3+ ion TIP4P/EW water from Li et al. JPCB, 2015, 119, 883'),
                   'Tb3': (1.619, 0.08912336, 'IOD set for Tb3+ ion TIP4P/EW water from Li et al. JPCB, 2015, 119, 883'),
                   'Dy3': (1.597, 0.07786298, 'IOD set for Dy3+ ion TIP4P/EW water from Li et al. JPCB, 2015, 119, 883'),
                   'Er3': (1.590, 0.07447106, 'IOD set for Er3+ ion TIP4P/EW water from Li et al. JPCB, 2015, 119, 883'),
                   'Tm3': (1.590, 0.07447106, 'IOD set for Tm3+ ion TIP4P/EW water from Li et al. JPCB, 2015, 119, 883'),
                   'Lu3': (1.577, 0.06841702, 'IOD set for Lu3+ ion TIP4P/EW water from Li et al. JPCB, 2015, 119, 883'),
                   'Hf4': (1.483, 0.03393126, 'IOD set for Hf4+ ion TIP4P/EW water from Li et al. JPCB, 2015, 119, 883'),
                   'Zr4': (1.503, 0.03994409, 'IOD set for Zr4+ ion TIP4P/EW water from Li et al. JPCB, 2015, 119, 883'),
                   'Ce4': (1.667, 0.11679623, 'IOD set for Ce4+ ion TIP4P/EW water from Li et al. JPCB, 2015, 119, 883'),
                   'U4':  (1.667, 0.11679623, 'IOD set for U4+ ion TIP4P/EW water from Li et al. JPCB, 2015, 119, 883'),
                   'Pu4': (1.645, 0.10359269, 'IOD set for Pu4+ ion TIP4P/EW water from Li et al. JPCB, 2015, 119, 883'),
                   'Th4': (1.690, 0.13150785, 'IOD set for Th4+ ion TIP4P/EW water from Li et al. JPCB, 2015, 119, 883'),
                     }
    elif watermodel == 'opc3':
        higljpara = {
                   'Al3': (1.287, 0.00412163, 'IOD set for Al3+ ion OPC3 water from Li et al. JCTC, 2021, 17, 2342'),
                   'Fe3': (1.419, 0.01900380, 'IOD set for Fe3+ ion OPC3 water from Li et al. JCTC, 2021, 17, 2342'),
                   'Cr3': (1.364, 0.01067299, 'IOD set for Cr3+ ion OPC3 water from Li et al. JCTC, 2021, 17, 2342'),
                   'In3': (1.453, 0.02615377, 'IOD set for In3+ ion OPC3 water from Li et al. JCTC, 2021, 17, 2342'),
                   'Tl3': (1.496, 0.03776169, 'IOD set for Tl3+ ion OPC3 water from Li et al. JCTC, 2021, 17, 2342'),
                   'Y3':  (1.600, 0.07934493, 'IOD set for Y3+ ion OPC3 water from Li et al. JCTC, 2021, 17, 2342'),
                   'La3': (1.733, 0.16134811, 'IOD set for La3+ ion OPC3 water from Li et al. JCTC, 2021, 17, 2342'),
                   'Ce3': (1.758, 0.17997960, 'IOD set for Ce3+ ion OPC3 water from Li et al. JCTC, 2021, 17, 2342'),
                   'Pr3': (1.750, 0.17392181, 'IOD set for Pr3+ ion OPC3 water from Li et al. JCTC, 2021, 17, 2342'),
                   'Nd3': (1.692, 0.13282966, 'IOD set for Nd3+ ion OPC3 water from Li et al. JCTC, 2021, 17, 2342'),
                   'Sm3': (1.667, 0.11679623, 'IOD set for Sm3+ ion OPC3 water from Li et al. JCTC, 2021, 17, 2342'),
                   'Eu3': (1.675, 0.12180998, 'IOD set for Eu3+ ion OPC3 water from Li et al. JCTC, 2021, 17, 2342'),
                   'Gd3': (1.625, 0.09235154, 'IOD set for Gd3+ ion OPC3 water from Li et al. JCTC, 2021, 17, 2342'),
                   'Tb3': (1.633, 0.09675968, 'IOD set for Tb3+ ion OPC3 water from Li et al. JCTC, 2021, 17, 2342'),
                   'Dy3': (1.608, 0.08337961, 'IOD set for Dy3+ ion OPC3 water from Li et al. JCTC, 2021, 17, 2342'),
                   'Er3': (1.600, 0.07934493, 'IOD set for Er3+ ion OPC3 water from Li et al. JCTC, 2021, 17, 2342'),
                   'Tm3': (1.600, 0.07934493, 'IOD set for Tm3+ ion OPC3 water from Li et al. JCTC, 2021, 17, 2342'),
                   'Lu3': (1.583, 0.07117158, 'IOD set for Lu3+ ion OPC3 water from Li et al. JCTC, 2021, 17, 2342'),
                   'Hf4': (1.489, 0.03566355, 'IOD set for Hf4+ ion OPC3 water from Li et al. JCTC, 2021, 17, 2342'),
                   'Zr4': (1.508, 0.04155519, 'IOD set for Zr4+ ion OPC3 water from Li et al. JCTC, 2021, 17, 2342'),
                   'Ce4': (1.708, 0.14364160, 'IOD set for Ce4+ ion OPC3 water from Li et al. JCTC, 2021, 17, 2342'),
                   'U4':  (1.708, 0.14364160, 'IOD set for U4+ ion OPC3 water from Li et al. JCTC, 2021, 17, 2342'),
                   'Pu4': (1.682, 0.12628793, 'IOD set for Pu4+ ion OPC3 water from Li et al. JCTC, 2021, 17, 2342'),
                   'Th4': (1.721, 0.15272873, 'IOD set for Th4+ ion OPC3 water from Li et al. JCTC, 2021, 17, 2342'),
                     }
    elif watermodel == 'opc':
        higljpara = {
                   'Al3': (1.250, 0.00243637, 'IOD set for Al3+ ion OPC water from Li et al. JCTC, 2021, 17, 2342'),
                   'Fe3': (1.400, 0.01570749, 'IOD set for Fe3+ ion OPC water from Li et al. JCTC, 2021, 17, 2342'),
                   'Cr3': (1.336, 0.00770969, 'IOD set for Cr3+ ion OPC water from Li et al. JCTC, 2021, 17, 2342'),
                   'In3': (1.440, 0.02322071, 'IOD set for In3+ ion OPC water from Li et al. JCTC, 2021, 17, 2342'),
                   'Tl3': (1.483, 0.03393126, 'IOD set for Tl3+ ion OPC water from Li et al. JCTC, 2021, 17, 2342'),
                   'Y3':  (1.575, 0.06751391, 'IOD set for Y3+ ion OPC water from Li et al. JCTC, 2021, 17, 2342'),
                   'La3': (1.700, 0.13818331, 'IOD set for La3+ ion OPC water from Li et al. JCTC, 2021, 17, 2342'),
                   'Ce3': (1.725, 0.15557763, 'IOD set for Ce3+ ion OPC water from Li et al. JCTC, 2021, 17, 2342'),
                   'Pr3': (1.717, 0.14990448, 'IOD set for Pr3+ ion OPC water from Li et al. JCTC, 2021, 17, 2342'),
                   'Nd3': (1.662, 0.11371963, 'IOD set for Nd3+ ion OPC water from Li et al. JCTC, 2021, 17, 2342'),
                   'Sm3': (1.638, 0.09957472, 'IOD set for Sm3+ ion OPC water from Li et al. JCTC, 2021, 17, 2342'),
                   'Eu3': (1.646, 0.10417397, 'IOD set for Eu3+ ion OPC water from Li et al. JCTC, 2021, 17, 2342'),
                   'Gd3': (1.600, 0.07934493, 'IOD set for Gd3+ ion OPC water from Li et al. JCTC, 2021, 17, 2342'),
                   'Tb3': (1.608, 0.08337961, 'IOD set for Tb3+ ion OPC water from Li et al. JCTC, 2021, 17, 2342'),
                   'Dy3': (1.583, 0.07117158, 'IOD set for Dy3+ ion OPC water from Li et al. JCTC, 2021, 17, 2342'),
                   'Er3': (1.575, 0.06751391, 'IOD set for Er3+ ion OPC water from Li et al. JCTC, 2021, 17, 2342'),
                   'Tm3': (1.575, 0.06751391, 'IOD set for Tm3+ ion OPC water from Li et al. JCTC, 2021, 17, 2342'),
                   'Lu3': (1.558, 0.06014121, 'IOD set for Lu3+ ion OPC water from Li et al. JCTC, 2021, 17, 2342'),
                   'Hf4': (1.461, 0.02808726, 'IOD set for Hf4+ ion OPC water from Li et al. JCTC, 2021, 17, 2342'),
                   'Zr4': (1.488, 0.03537062, 'IOD set for Zr4+ ion OPC water from Li et al. JCTC, 2021, 17, 2342'),
                   'Ce4': (1.683, 0.12693448, 'IOD set for Ce4+ ion OPC water from Li et al. JCTC, 2021, 17, 2342'),
                   'U4':  (1.683, 0.12693448, 'IOD set for U4+ ion OPC water from Li et al. JCTC, 2021, 17, 2342'),
                   'Pu4': (1.658, 0.11129023, 'IOD set for Pu4+ ion OPC water from Li et al. JCTC, 2021, 17, 2342'),
                   'Th4': (1.704, 0.14089951, 'IOD set for Th4+ ion OPC water from Li et al. JCTC, 2021, 17, 2342'),
                     }
    elif watermodel == 'fb3':
        higljpara = {
                   'Al3': (1.287, 0.00412163, 'IOD set for Al3+ ion TIP3P-FB water from Li et al. JCTC, 2021, 17, 2342'),
                   'Fe3': (1.419, 0.01900380, 'IOD set for Fe3+ ion TIP3P-FB water from Li et al. JCTC, 2021, 17, 2342'),
                   'Cr3': (1.364, 0.01067299, 'IOD set for Cr3+ ion TIP3P-FB water from Li et al. JCTC, 2021, 17, 2342'),
                   'In3': (1.445, 0.02431873, 'IOD set for In3+ ion TIP3P-FB water from Li et al. JCTC, 2021, 17, 2342'),
                   'Tl3': (1.496, 0.03776169, 'IOD set for Tl3+ ion TIP3P-FB water from Li et al. JCTC, 2021, 17, 2342'),
                   'Y3':  (1.600, 0.07934493, 'IOD set for Y3+ ion TIP3P-FB water from Li et al. JCTC, 2021, 17, 2342'),
                   'La3': (1.731, 0.15989650, 'IOD set for La3+ ion TIP3P-FB water from Li et al. JCTC, 2021, 17, 2342'),
                   'Ce3': (1.754, 0.17693975, 'IOD set for Ce3+ ion TIP3P-FB water from Li et al. JCTC, 2021, 17, 2342'),
                   'Pr3': (1.746, 0.17092614, 'IOD set for Pr3+ ion TIP3P-FB water from Li et al. JCTC, 2021, 17, 2342'),
                   'Nd3': (1.692, 0.13282966, 'IOD set for Nd3+ ion TIP3P-FB water from Li et al. JCTC, 2021, 17, 2342'),
                   'Sm3': (1.667, 0.11679623, 'IOD set for Sm3+ ion TIP3P-FB water from Li et al. JCTC, 2021, 17, 2342'),
                   'Eu3': (1.675, 0.12180998, 'IOD set for Eu3+ ion TIP3P-FB water from Li et al. JCTC, 2021, 17, 2342'),
                   'Gd3': (1.625, 0.09235154, 'IOD set for Gd3+ ion TIP3P-FB water from Li et al. JCTC, 2021, 17, 2342'),
                   'Tb3': (1.633, 0.09675968, 'IOD set for Tb3+ ion TIP3P-FB water from Li et al. JCTC, 2021, 17, 2342'),
                   'Dy3': (1.608, 0.08337961, 'IOD set for Dy3+ ion TIP3P-FB water from Li et al. JCTC, 2021, 17, 2342'),
                   'Er3': (1.600, 0.07934493, 'IOD set for Er3+ ion TIP3P-FB water from Li et al. JCTC, 2021, 17, 2342'),
                   'Tm3': (1.600, 0.07934493, 'IOD set for Tm3+ ion TIP3P-FB water from Li et al. JCTC, 2021, 17, 2342'),
                   'Lu3': (1.583, 0.07117158, 'IOD set for Lu3+ ion TIP3P-FB water from Li et al. JCTC, 2021, 17, 2342'),
                   'Hf4': (1.489, 0.03566355, 'IOD set for Hf4+ ion TIP3P-FB water from Li et al. JCTC, 2021, 17, 2342'),
                   'Zr4': (1.508, 0.04155519, 'IOD set for Zr4+ ion TIP3P-FB water from Li et al. JCTC, 2021, 17, 2342'),
                   'Ce4': (1.705, 0.14158262, 'IOD set for Ce4+ ion TIP3P-FB water from Li et al. JCTC, 2021, 17, 2342'),
                   'U4':  (1.705, 0.14158262, 'IOD set for U4+ ion TIP3P-FB water from Li et al. JCTC, 2021, 17, 2342'),
                   'Pu4': (1.682, 0.12628793, 'IOD set for Pu4+ ion TIP3P-FB water from Li et al. JCTC, 2021, 17, 2342'),
                   'Th4': (1.718, 0.15060822, 'IOD set for Th4+ ion TIP3P-FB water from Li et al. JCTC, 2021, 17, 2342'),
                     }
    elif watermodel == 'fb4':
        higljpara = {
                   'Al3': (1.267, 0.00312065, 'IOD set for Al3+ ion TIP4P-FB water from Li et al. JCTC, 2021, 17, 2342'),
                   'Fe3': (1.404, 0.01636246, 'IOD set for Fe3+ ion TIP4P-FB water from Li et al. JCTC, 2021, 17, 2342'),
                   'Cr3': (1.345, 0.00858042, 'IOD set for Cr3+ ion TIP4P-FB water from Li et al. JCTC, 2021, 17, 2342'),
                   'In3': (1.443, 0.02387506, 'IOD set for In3+ ion TIP4P-FB water from Li et al. JCTC, 2021, 17, 2342'),
                   'Tl3': (1.488, 0.03537062, 'IOD set for Tl3+ ion TIP4P-FB water from Li et al. JCTC, 2021, 17, 2342'),
                   'Y3':  (1.583, 0.07117158, 'IOD set for Y3+ ion TIP4P-FB water from Li et al. JCTC, 2021, 17, 2342'),
                   'La3': (1.715, 0.14850170, 'IOD set for La3+ ion TIP4P-FB water from Li et al. JCTC, 2021, 17, 2342'),
                   'Ce3': (1.738, 0.16500296, 'IOD set for Ce3+ ion TIP4P-FB water from Li et al. JCTC, 2021, 17, 2342'),
                   'Pr3': (1.731, 0.15989650, 'IOD set for Pr3+ ion TIP4P-FB water from Li et al. JCTC, 2021, 17, 2342'),
                   'Nd3': (1.675, 0.12180998, 'IOD set for Nd3+ ion TIP4P-FB water from Li et al. JCTC, 2021, 17, 2342'),
                   'Sm3': (1.650, 0.10651723, 'IOD set for Sm3+ ion TIP4P-FB water from Li et al. JCTC, 2021, 17, 2342'),
                   'Eu3': (1.658, 0.11129023, 'IOD set for Eu3+ ion TIP4P-FB water from Li et al. JCTC, 2021, 17, 2342'),
                   'Gd3': (1.608, 0.08337961, 'IOD set for Gd3+ ion TIP4P-FB water from Li et al. JCTC, 2021, 17, 2342'),
                   'Tb3': (1.617, 0.08806221, 'IOD set for Tb3+ ion TIP4P-FB water from Li et al. JCTC, 2021, 17, 2342'),
                   'Dy3': (1.592, 0.07543075, 'IOD set for Dy3+ ion TIP4P-FB water from Li et al. JCTC, 2021, 17, 2342'),
                   'Er3': (1.583, 0.07117158, 'IOD set for Er3+ ion TIP4P-FB water from Li et al. JCTC, 2021, 17, 2342'),
                   'Tm3': (1.583, 0.07117158, 'IOD set for Tm3+ ion TIP4P-FB water from Li et al. JCTC, 2021, 17, 2342'),
                   'Lu3': (1.567, 0.06397679, 'IOD set for Lu3+ ion TIP4P-FB water from Li et al. JCTC, 2021, 17, 2342'),
                   'Hf4': (1.467, 0.02960343, 'IOD set for Hf4+ ion TIP4P-FB water from Li et al. JCTC, 2021, 17, 2342'),
                   'Zr4': (1.495, 0.03745682, 'IOD set for Zr4+ ion TIP4P-FB water from Li et al. JCTC, 2021, 17, 2342'),
                   'Ce4': (1.692, 0.13282966, 'IOD set for Ce4+ ion TIP4P-FB water from Li et al. JCTC, 2021, 17, 2342'),
                   'U4':  (1.692, 0.13282966, 'IOD set for U4+ ion TIP4P-FB water from Li et al. JCTC, 2021, 17, 2342'),
                   'Pu4': (1.667, 0.11679623, 'IOD set for Pu4+ ion TIP4P-FB water from Li et al. JCTC, 2021, 17, 2342'),
                   'Th4': (1.709, 0.14433113, 'IOD set for Th4+ ion TIP4P-FB water from Li et al. JCTC, 2021, 17, 2342'),
                     }

    IonLJParaDict.update(higljpara)

    otherljpara = {
                'Si0': (2.1475, 0.402, 'Adopted from atom type Si3 from UFF (Rappe et al. JACS, 114, 10024)'),
                'Sc3': (1.6475, 0.019, 'Adopted from atom type Sc3+3 from UFF (Rappe et al. JACS, 114, 10024)'),
                'Ti4': (1.5875, 0.017, 'Adopted from atom type Ti3+4/Ti6+4 from UFF (Rappe et al. JACS, 114, 10024)'),
                'Ga3': (2.1915, 0.415, 'Adopted from atom type Ga3+3 from UFF (Rappe et al. JACS, 114, 10024)'),
                'Ge0': (2.1400, 0.379, 'Adopted from atom type Ge3 from UFF (Rappe et al. JACS, 114, 10024)'),
                'As3': (2.1150, 0.309, 'Adopted from atom type As3+3 from UFF (Rappe et al. JACS, 114, 10024)'),
                'Se2': (2.1025, 0.291, 'Adopted from atom type Se3+2 from UFF (Rappe et al. JACS, 114, 10024)'),
                'Nb5': (1.5825, 0.059, 'Adopted from atom type Nb3+5 from UFF (Rappe et al. JACS, 114, 10024)'),
                'Mo6': (1.5060, 0.056, 'Adopted from atom type Mo3+6/Mo6+6 from UFF (Rappe et al. JACS, 114, 10024)'),
                'Tc5': (1.4990, 0.048, 'Adopted from atom type Tc6+5 from UFF (Rappe et al. JACS, 114, 10024)'),
                'Ru2': (1.4815, 0.056, 'Adopted from atom type Ru6+2 from UFF (Rappe et al. JACS, 114, 10024)'),
                'Rh3': (1.4645, 0.053, 'Adopted from atom type Rh6+3 from UFF (Rappe et al. JACS, 114, 10024)'),
                'Sb3': (2.2100, 0.449, 'Adopted from atom type Sb3+3 from UFF (Rappe et al. JACS, 114, 10024)'),
                'Te2': (2.2350, 0.398, 'Adopted from atom type Te3+2 from UFF (Rappe et al. JACS, 114, 10024)'),
                'Pm3': (1.7735, 0.009, 'Adopted from atom type Pm6+3 from UFF (Rappe et al. JACS, 114, 10024)'),
                'Ho3': (1.7045, 0.007, 'Adopted from atom type Ho6+3 from UFF (Rappe et al. JACS, 114, 10024)'),
                'Ta5': (1.5850, 0.081, 'Adopted from atom type Ta3+5 from UFF (Rappe et al. JACS, 114, 10024)'),
                'W4' : (1.5345, 0.067, 'Adopted from atom type W3+4 from UFF (Rappe et al. JACS, 114, 10024)'),
                'W6' : (1.5345, 0.067, 'Adopted from atom type W3+6/W6+6 from UFF (Rappe et al. JACS, 114, 10024)'),
                'Re5': (1.4770, 0.066, 'Adopted from atom type Re6+5 from UFF (Rappe et al. JACS, 114, 10024)'),
                'Re7': (1.4770, 0.066, 'Adopted from atom type Re3+7 from UFF (Rappe et al. JACS, 114, 10024)'),
                'Os6': (1.5600, 0.037, 'Adopted from atom type Os6+6 from UFF (Rappe et al. JACS, 114, 10024)'),
                'Ir3': (1.4200, 0.073, 'Adopted from atom type Ir6+3 from UFF (Rappe et al. JACS, 114, 10024)'),
                'Au3': (1.6465, 0.039, 'Adopted from atom type Au4+3 from UFF (Rappe et al. JACS, 114, 10024)'),
                'Bi3': (2.1850, 0.518, 'Adopted from atom type Bi3+3 from UFF (Rappe et al. JACS, 114, 10024)'),
                'Po2': (2.3545, 0.325, 'Adopted from atom type Po3+2 from UFF (Rappe et al. JACS, 114, 10024)'),
                'At0': (2.3750, 0.284, 'Adopted from atom type At from UFF (Rappe et al. JACS, 114, 10024)'),
                'Fr0': (2.4500, 0.050, 'Adopted from atom type Fr from UFF (Rappe et al. JACS, 114, 10024)'),
                'Ac3': (1.7390, 0.330, 'Adopted from atom type Ac6+3 from UFF (Rappe et al. JACS, 114, 10024)'),
                'Pa4': (1.7120, 0.022, 'Adopted from atom type Pa6+4 from UFF (Rappe et al. JACS, 114, 10024)'),
                'Np4': (1.7120, 0.019, 'Adopted from atom type Np6+4 from UFF (Rappe et al. JACS, 114, 10024)'),
                'Am4': (1.6905, 0.014, 'Adopted from atom type Am6+4 from UFF (Rappe et al. JACS, 114, 10024)'),
                'Cm3': (1.6630, 0.013, 'Adopted from atom type Cm6+3 from UFF (Rappe et al. JACS, 114, 10024)'),
                'Bk3': (1.6695, 0.013, 'Adopted from atom type Bk6+3 from UFF (Rappe et al. JACS, 114, 10024)'),
                'Cf3': (1.6565, 0.013, 'Adopted from atom type Cf6+3 from UFF (Rappe et al. JACS, 114, 10024)'),
                'Es3': (1.6495, 0.012, 'Adopted from atom type Es6+3 from UFF (Rappe et al. JACS, 114, 10024)'),
                'Fm3': (1.6430, 0.012, 'Adopted from atom type Fm6+3 from UFF (Rappe et al. JACS, 114, 10024)'),
                'Md3': (1.6370, 0.011, 'Adopted from atom type Md6+3 from UFF (Rappe et al. JACS, 114, 10024)'),
                'No3': (1.6240, 0.011, 'Adopted from atom type No6+3 from UFF (Rappe et al. JACS, 114, 10024)'),
                'Lr3': (1.6180, 0.011, 'Adopted from atom type Lw6+3 (Lw is the same as Lr) from UFF (Rappe et al. JACS, 114, 10024)'),
                  }

    IonLJParaDict.update(otherljpara)

    return IonLJParaDict

#LJ parameter for nonbonded model
IonLJparal = [ 'Li1', 'Na1',  'K1', 'Rb1', 'Cs1', 'Tl1', 'Cu1', 'Ag1',  'H1', 'F-1', 
              'Cl-1','Br-1', 'I-1', 'Be2', 'Cu2', 'Ni2', 'Zn2', 'Co2', 'Cr2', 'Fe2',  
               'Mg2',  'V2', 'Mn2', 'Hg2', 'Cd2', 'Ca2', 'Sn2', 'Sr2', 'Ba2', 'Pt2',  
               'Pd2', 'Ag2', 'Yb2', 'Pb2', 'Eu2', 'Sm2', 'Ra2', 'Al3', 'Fe3', 'Cr3',  
               'In3', 'Tl3',  'Y3', 'La3', 'Ce3', 'Pr3', 'Nd3', 'Sm3', 'Eu3', 'Gd3',  
               'Tb3', 'Dy3', 'Er3', 'Tm3', 'Lu3', 'Hf4', 'Zr4', 'Ce4',  'U4', 'Pu4',  
               'Th4']

IonCMparal = ['Be2', 'Cu2', 'Ni2', 'Zn2', 'Co2', 'Cr2', 'Fe2', 'Mg2',  'V2', 'Mn2',
              'Hg2', 'Cd2', 'Ca2', 'Sn2', 'Sr2', 'Ba2', 'Pt2', 'Pd2', 'Ag2', 'Yb2',
              'Pb2', 'Eu2', 'Sm2', 'Ra2']

IonHFEparal = IonLJparal
IonHFEparal.remove('H1')

IonIODparal = list(set(IonLJparal) - \
              set(['Pt2', 'Pd2', 'Ag2', 'Yb2', 'Pb2', 'Eu2', 'Sm2', 'Ra2']))

#-----------------------------------------------------------------------------

#Atom numbers

#-----------------------------------------------------------------------------
Atnum = AtomicNum

AtnumRev = dict([ (v, k) for k, v in list(Atnum.items())])

bdld = {'CH': 1.090, 'NH': 1.010}
