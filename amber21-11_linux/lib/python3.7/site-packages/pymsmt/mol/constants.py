"This module is for the constants"

B_TO_A = 0.529177249 #Bohr to Angstrom
H_TO_KCAL_MOL = 627.5095 #Hatree to kcal/mol
#HB2_TO_KCAL_MOL_A2 = 2240.87 #2240.87 is Hatree/(Bohr^2) to kcal/(mol*angstrom^2)
HB2_TO_KCAL_MOL_A2 = H_TO_KCAL_MOL/(B_TO_A**2)
