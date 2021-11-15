atom_type_dict = {
 'SC4': 'sA', 'OC23': 'oA', 'OC24': 'oB', 'HOY': 'hA',  'NA+': 'NAp', #Silica
  'K+': 'Kp',  'NA+': 'NAp', 'SY1': 'sB', 'SY2': 'sC', 'AYT1': 'aA',
'AYT2': 'aB',  'AY1': 'aC',  'AY2': 'aD', 'OY1': 'oC',  'OY2': 'oD',
 'OY3': 'oE',  'OY4': 'oF',  'OY5': 'oG', 'OY6': 'oH',  'OY7': 'oI',
 'OY8': 'oJ',  'OY9': 'oK',  'HOY': 'hB', 'HOK': 'hC', #Clay minerals
'CA++': 'CApp','SC1': 'sD',  'OC1': 'oL', 'OC2': 'oM', 'CA+A': 'CApA',
 'AC1': 'aE',  'OC3': 'oN',  'OC4': 'oO', 'OC5': 'oP',  'HOC': 'hD', #Cement minerals
'CA+H': 'CApH','PAP': 'pA', 'OAP1': 'oQ','OAP2': 'oR',  'HOP': 'hE', #Apatite
  'AL': 'Al0',  'NI': 'Ni0',  'CU': 'Cu0', 'PD': 'Pd0', 'AG' : 'Ag0', #Metal
  'PT': 'Pt0',  'AU': 'Au0',  'PB': 'Pb0'
}

for i in list(atom_type_dict.keys()):
    atom_type_dict[i.lower()] = atom_type_dict[i]

#for i in sorted(list(atom_type_dict.values())):
#    print i

atom_type_dict2 = {
 'SC4': '01', 'OC23': '02', 'OC24': '03', 'HOY': '04',  'NA+': 'NAp', #Silica
  'K+': 'Kp',  'NA+': 'NAp', 'SY1': '05', 'SY2': '06', 'AYT1': '07',
'AYT2': '08',  'AY1': '09',  'AY2': '10', 'OY1': '11',  'OY2': '12',
 'OY3': '13',  'OY4': '14',  'OY5': '15', 'OY6': '16',  'OY7': '17',
 'OY8': '18',  'OY9': '19',  'HOY': '20', 'HOK': '21', #Clay minerals
'CA++': 'CApp','SC1': '22',  'OC1': '23', 'OC2': '24', 'CA+A': 'CApA',
 'AC1': '25',  'OC3': '26',  'OC4': '27', 'OC5': '28',  'HOC': '29', #Cement minerals
'CA+H': 'CApH','PAP': '30', 'OAP1': '31','OAP2': '32',  'HOP': '33', #Apatite
  'AL': 'Al0',  'NI': 'Ni0',  'CU': 'Cu0', 'PD': 'Pd0', 'AG' : 'Ag0', #Metal
  'PT': 'Pt0',  'AU': 'Au0',  'PB': 'Pb0'
}

for i in list(atom_type_dict2.keys()):
    atom_type_dict2[i.lower()] = atom_type_dict2[i]

#for i in sorted(list(atom_type_dict2.values())):
#    print i

atom_type_dict3 = {
 'HOC': '01',  'HOK': '02',  'HOP': '03', 'HOY': '04',
'OAP1': '11', 'OAP2': '12',  'OC1': '13', 'OC2': '14',  'OC3': '15',
 'OC4': '16',  'OC5': '17', 'OC23': '18','OC24': '19',  'OY1': '20',
 'OY2': '21',  'OY3': '22',  'OY4': '23', 'OY5': '24',  'OY6': '25',
 'OY7': '26',  'OY8': '27',  'OY9': '28', 
 'AC1': '41',  'AY1': '42',  'AY2': '43','AYT1': '44', 'AYT2': '45',
 'SC1': '51',  'SC4': '52',  'SY1': '53', 'SY2': '54',
 'PAP': '61',
 'NA+': 'NAp',  'K+': 'Kp', 'CA++': 'CApp', 'CA+A': 'CApA', 'CA+H': 'CApH',
  'AL': 'Al0',  'NI': 'Ni0',  'CU': 'Cu0', 'PD': 'Pd0', 'AG' : 'Ag0',
  'PT': 'Pt0',  'AU': 'Au0',  'PB': 'Pb0'
}

for i in list(atom_type_dict3.keys()):
    atom_type_dict3[i.lower()] = atom_type_dict3[i]

