clearVariables
logFile ff91.log
addPath ../prep
addPath ../parm
addAtomTypes {
	{ "H"   "H" "sp3" }
	{ "HO"  "H" "sp3" }
	{ "HS"  "H" "sp3" }
	{ "H1"  "H" "sp3" }
	{ "H2"  "H" "sp3" }
	{ "H3"  "H" "sp3" }
	{ "H4"  "H" "sp3" }
	{ "H5"  "H" "sp3" }
	{ "HW"  "H" "sp3" }
	{ "HC"  "H" "sp3" }
	{ "HA"  "H" "sp3" }
	{ "HP"  "H" "sp3" }
	{ "OH"  "O" "sp3" }
	{ "OS"  "O" "sp3" }
	{ "O"   "O" "sp2" }
	{ "O2"  "O" "sp2" }
	{ "OW"  "O" "sp3" }
	{ "CT"  "C" "sp3" }
	{ "CH"  "C" "sp3" }
	{ "C2"  "C" "sp3" }
	{ "C3"  "C" "sp3" }
	{ "C"   "C" "sp2" }
	{ "C*"  "C" "sp2" }
	{ "CA"  "C" "sp2" }
	{ "CB"  "C" "sp2" }
	{ "CC"  "C" "sp2" }
	{ "CN"  "C" "sp2" }
	{ "CM"  "C" "sp2" }
	{ "CK"  "C" "sp2" }
	{ "CQ"  "C" "sp2" }
	{ "CD"  "C" "sp2" }
	{ "CE"  "C" "sp2" }
	{ "CF"  "C" "sp2" }
	{ "CG"  "C" "sp2" }
	{ "CP"  "C" "sp2" }
	{ "CI"  "C" "sp2" }
	{ "CJ"  "C" "sp2" }
	{ "CW"  "C" "sp2" }
	{ "CV"  "C" "sp2" }
	{ "CR"  "C" "sp2" }
	{ "CA"  "C" "sp2" }
	{ "CY"  "C" "sp2" }
	{ "C0"  "Ca" "sp3" }
	{ "MG"  "Mg" "sp3" }
	{ "N"   "N" "sp2" }
	{ "NA"  "N" "sp2" }
	{ "N2"  "N" "sp2" }
	{ "N*"  "N" "sp2" }
	{ "NP"  "N" "sp2" }
	{ "NQ"  "N" "sp2" }
	{ "NB"  "N" "sp2" }
	{ "NC"  "N" "sp2" }
	{ "NT"  "N" "sp3" }
	{ "N3"  "N" "sp3" }
	{ "S"   "S" "sp3" }
	{ "SH"  "S" "sp3" }
	{ "P"   "P" "sp3" }
	{ "LP"  ""  "sp3" }
	{ "F"   "F" "sp3" }
	{ "CL"  "Cl" "sp3" }
	{ "BR"  "Br" "sp3" }
	{ "I"   "I"  "sp3" }
	{ "FE"  "Fe" "sp3" }
# glycam
	{ "OG"  "O" "sp3" }
	{ "OL"  "O" "sp3" }
	{ "AC"  "C" "sp3" }
	{ "EC"  "C" "sp3" }
}
#
#	leap .cmd script for building the residue
#	libraries for the 1991 version of the Weiner
#	et al. force field.
#
#	United atom residues are omitted (see ff91_uni.cmd)
#	since discrepancies between leap and prep/link/edit/parm
#	were not resolved (improper dihedrals).
#
#
################################################
################################################
################################################
################################################
######
######    AMBER 4 prep.in files
######
################################################
################################################
################################################
################################################
#---------------------------------------------------
#
#
#       ALL ATOM FORCE FIELD
#
#
#
# Extract the amino acids from all.in
#
loadAmberPrep  all.in 
a = { ALA GLY SER THR LEU ILE VAL ASN GLN ARG 
	HID HIE HIP TRP PHE TYR GLU ASP LYS PRO CYS CYX MET }
set a   restype   protein

set CYX.1   disulphide  CYX.1.SG

saveOff a   ./all_amino91.lib 

set NME   restype   protein
set NME   tail  null
set NME   head  NME.1.N
set NME.1   connect0  NME.1.N
saveOff NME   ./all_aminoct91.lib 

set NHE   restype   protein
set NHE   tail  null
set NHE   head  NHE.1.N
set NHE.1   connect0  NHE.1.N
saveOff NHE   ./all_aminoct91.lib 

set ACE   restype   protein
set ACE   head  null
set ACE   tail  ACE.1.C
set ACE.1   connect1  ACE.1.C
saveOff ACE   ./all_aminont91.lib 

#
# Save the DNA and RNA residues
#

a = { RADE RCYT RGUA ROHE RPOM RURA }
set a   restype   nucleic
saveOff a   ./all_nucleic91.lib 
a = { DADE DCYT DGUA DOHE DPOM DTHY }
set a   restype   nucleic
saveOff a   ./all_nucleic91.lib 

set HB   head  null
set HB.1   connect0  null
set HE   tail  null
set HE.1   connect1  null

a = { HB HE }
set a   restype   nucleic
saveOff a   ./all_nucleic91.lib 

#
# Extract the N terminus residues
#
clearVariables

loadAmberPrep  allnt.in    N
a = { NALA NGLY NSER NTHR NLEU NILE NVAL NASN NGLN NARG 
NHID NHIE NHIP NTRP NPHE
NTYR NGLU NASP NLYS NPRO NCYS NCYX NMET }
set a   head  null

set NALA.1   nend  null
set NGLY.1   nend  null
set NSER.1   nend  null
set NTHR.1   nend  null
set NLEU.1   nend  null
set NILE.1   nend  null
set NVAL.1   nend  null
set NASN.1   nend  null
set NGLN.1   nend  null
set NARG.1   nend  null
set NHID.1   nend  null
set NHIE.1   nend  null
set NHIP.1   nend  null
set NTRP.1   nend  null
set NPHE.1   nend  null
set NTYR.1   nend  null
set NGLU.1   nend  null
set NASP.1   nend  null
set NLYS.1   nend  null
set NPRO.1   nend  null
set NCYS.1   nend  null
set NCYX.1   nend  null
set NMET.1   nend  null

set a   restype   protein

set NCYX.1   disulphide  NCYX.1.SG

saveOff a   ./all_aminont91.lib 

#
# Extract the C terminus residues
#
loadAmberPrep  allct.in    C
a = { CALA CGLY CSER CTHR CLEU CILE CVAL CASN CGLN CARG 
CHID CHIE CHIP CTRP CPHE
CTYR CGLU CASP CLYS CPRO CCYS CCYX CMET }
set a   tail  null

set CALA.1   cend  null
set CGLY.1   cend  null
set CSER.1   cend  null
set CTHR.1   cend  null
set CLEU.1   cend  null
set CILE.1   cend  null
set CVAL.1   cend  null
set CASN.1   cend  null
set CGLN.1   cend  null
set CARG.1   cend  null
set CHID.1   cend  null
set CHIE.1   cend  null
set CHIP.1   cend  null
set CTRP.1   cend  null
set CPHE.1   cend  null
set CTYR.1   cend  null
set CGLU.1   cend  null
set CASP.1   cend  null
set CLYS.1   cend  null
set CPRO.1   cend  null
set CCYS.1   cend  null
set CCYX.1   cend  null
set CMET.1   cend  null

set a   restype   protein

set CCYX.1   disulphide  CCYX.1.SG

saveOff a   ./all_aminoct91.lib 
quit
