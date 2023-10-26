clearVariables
logFile ions91.log
#
#	Monovalent, monoatomic ions using '91 atom types.
#	The Li+..Cs+ series uses types derived from the
#	work of Aqvist (see force field documentation),
#	while IB is a 'fat ion w/ waters' for vacuum use.
#

i = createAtom   Li+  QL  1.0
set i    element Li
set i    position { 0 0 0 }
r = createResidue Li+
add r i
Li+ = createUnit Li+
add Li+ r
saveOff Li+ ./ions91.lib

i = createAtom   Na+  QN  1.0
set i    element Na
set i    position { 0 0 0 }
r = createResidue Na+
add r i
Na+ = createUnit Na+
add Na+ r
saveOff Na+ ./ions91.lib

i = createAtom   K+   QK  1.0
set i    element K
set i    position { 0 0 0 }
r = createResidue K+
add r i
K+ = createUnit K+
add K+ r
saveOff K+ ./ions91.lib

i = createAtom   Rb+  QR  1.0
set i    element Rb
set i    position { 0 0 0 }
r = createResidue Rb+
add r i
Rb+ = createUnit Rb+
add Rb+ r
saveOff Rb+ ./ions91.lib

i = createAtom   Cs+  QC  1.0
set i    element Cs
set i    position { 0 0 0 }
r = createResidue Cs+
add r i
Cs+ = createUnit Cs+
add Cs+ r
saveOff Cs+ ./ions91.lib

i = createAtom   Cl-  IM  -1.0
set i    element Cl
set i    position { 0 0 0 }
r = createResidue Cl-
add r i
Cl- = createUnit Cl-
add Cl- r
saveOff Cl- ./ions91.lib

i = createAtom   IB   IB  1.0
set i    element Cs
set i    position { 0 0 0 }
r = createResidue IB
add r i
IB = createUnit IB
add IB r
saveOff IB ./ions91.lib


quit
