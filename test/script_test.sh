#!/bin/sh

rm -fr temp_test

mkdir -p temp_test

cd temp_test

pdb2gmx="/sw/bin/pdb2gmx"
editconf="/sw/bin/editconf"
genbox="/sw/bin/genbox"
genion="/sw/bin/genion"
grompp="/sw/bin/grompp"
mdrun="/sw/bin/mdrun"
mrun="" #"${mrun}"

echo "\n#=#=# Test DMP from 1BVG for GROMACS #=#=#\n"

wget -c "http://www.pdbe.org/download/1BVG" -O 1BVG.pdb
grep 'ATOM  ' 1BVG.pdb>| Protein.pdb
grep 'HETATM' 1BVG.pdb>| Ligand.pdb

echo "\n#=#=# CHECK 1BVG.pdb" >| diff_out.log
diff 1BVG.pdb ../Data/1BVG.pdb >> diff_out.log

# Edit Protein.pdb according to ffAMBER (http://ffamber.cnsm.csulb.edu/#usage)
sed s/PRO\ A\ \ \ 1/NPROA\ \ \ 1/g Protein.pdb | sed s/PRO\ B\ \ \ 1/NPROB\ \ \ 1/g \
| sed s/PHE\ A\ \ 99/CPHEA\ \ 99/g | sed s/PHE\ B\ \ 99/CPHEB\ \ 99/g \
| sed s/O\ \ \ CPHE/OC1\ CPHE/g | sed s/OXT\ CPHE/OC2\ CPHE/g \
| sed s/HIS\ /HID\ /g | sed s/LYS\ /LYP\ /g | sed s/CYS\ /CYN\ /g >| ProteinAmber.pdb

\cp Protein.pdb ProteinAmber.pdb

# Process with pdb2gmx and define water
${pdb2gmx} -ff amber99sb -f ProteinAmber.pdb -o Protein2.pdb -p Protein.top -water spce -ignh

# Generate Ligand topology file with acpype (GAFF)
acpype -i Ligand.pdb

# Merge Protein2.pdb + updated Ligand_NEW.pdb -> Complex.pdb
grep -h ATOM Protein2.pdb Ligand.acpype/Ligand_NEW.pdb >| Complex.pdb

# Edit Protein.top -> Complex.top
\cp Ligand.acpype/Ligand_GMX.itp Ligand.itp
\cp Protein.top Complex.top
#  '#include "Ligand.itp"' has to be inserted right after ffamber**.itp line and before Protein_*.itp line in Complex.top.
cat Complex.top | sed '/forcefield\.itp\"/a\
#include "Ligand.itp"
' >| Complex2.top
echo "Ligand   1" >> Complex2.top
\mv Complex2.top Complex.top

echo "\n#=#=# CHECK parm99gaffff99SBparmbsc0File" >> diff_out.log
# Generate Ligand topology file with acpype (AMBERbsc0)
acpype -i Ligand.pdb -a amber -b Ligand_Amber -c gas
diff -w /tmp/parm99gaffff99SB.dat ../../ffamber_additions/parm99SBgaff.dat >> diff_out.log

echo "\n#=#=# CHECK Ligand.itp" >> diff_out.log
diff Ligand.itp ../Data/Ligand.itp >> diff_out.log

# Setup the box and add water
${editconf} -bt triclinic -f Complex.pdb -o Complex.pdb -d 1.0
${genbox} -cp Complex.pdb -cs ffamber_tip3p.gro -o Complex_b4ion.pdb -p Complex.top

# Create em.mdp file
cat << EOF >| EM.mdp
define                   = -DFLEXIBLE
integrator               = cg ; steep
nsteps                   = 200
constraints              = none
emtol                    = 1000.0
nstcgsteep               = 10 ; do a steep every 10 steps of cg
emstep                   = 0.01 ; used with steep
nstcomm                  = 1
coulombtype              = PME
ns_type                  = grid
rlist                    = 1.0
rcoulomb                 = 1.0
rvdw                     = 1.4
Tcoupl                   = no
Pcoupl                   = no
gen_vel                  = no
nstxout                  = 0 ; write coords every # step
optimize_fft             = yes
EOF

# Create md.mdp file
cat << EOF >| MD.mdp
integrator               = md
nsteps                   = 1000
dt                       = 0.002
constraints              = all-bonds
nstcomm                  = 1
ns_type                  = grid
rlist                    = 1.2
rcoulomb                 = 1.1
rvdw                     = 1.0
vdwtype                  = shift
rvdw-switch              = 0.9
coulombtype              = PME-Switch
Tcoupl                   = v-rescale
tau_t                    = 0.1 0.1
tc-grps                  = protein non-protein
ref_t                    = 300 300
Pcoupl                   = parrinello-rahman
Pcoupltype               = isotropic
tau_p                    = 0.5
compressibility          = 4.5e-5
ref_p                    = 1.0
gen_vel                  = yes
nstxout                  = 2 ; write coords every # step
lincs-iter               = 2
DispCorr                 = EnerPres
optimize_fft             = yes
EOF

# Setup ions
${grompp} -f EM.mdp -c Complex_b4ion.pdb -p Complex.top -o Complex_b4ion.tpr
\cp Complex.top Complex_ion.top
echo 15| ${genion} -s Complex_b4ion.tpr -o Complex_b4em.pdb -neutral -conc 0.15 -p Complex_ion.top -norandom
\mv Complex_ion.top Complex.top

# Run minimisaton
#${grompp} -f EM.mdp -c Complex_b4em.pdb -p Complex.top -o em.tpr
#${mdrun} -v -deffnm em

# Run a short simulation
#${grompp} -f MD.mdp -c em.gro -p Complex.top -o md.tpr
#${mdrun} -v -deffnm md

# or with openmpi, for a dual core
${grompp} -f EM.mdp -c Complex_b4em.pdb -p Complex.top -o em.tpr
${mrun} ${mdrun} -v -deffnm em
${grompp} -f MD.mdp -c em.gro -p Complex.top -o md.tpr
${mrun} ${mdrun} -v -deffnm md
# vmd md.gro md.trr

echo "\n#=#=# Test 'amb2gmx' Function on 1BVG #=#=#\n"

# Edit 1BVG.pdb and create ComplexAmber.pdb
sed s/HIS\ /HID\ /g 1BVG.pdb | sed s/H1\ \ PRO/H3\ \ PRO/g >| ComplexAmber.pdb

# Create a input file with commands for tleap and run it
cat << EOF >| leap.in
verbosity 1
source leaprc.ff99SB
source leaprc.gaff
loadoff Ligand.acpype/Ligand_AC.lib
loadamberparams Ligand.acpype/Ligand_AC.frcmod
complex = loadpdb ComplexAmber.pdb
solvatebox complex TIP3PBOX 10.0
addions complex Na+ 23
addions complex Cl- 27
saveamberparm complex ComplexAmber.prmtop ComplexAmber.inpcrd
savepdb complex ComplexNAMD.pdb
quit
EOF
tleap -f leap.in >| leap.out

# convert AMBER to GROMACS
acpype -p ComplexAmber.prmtop -x ComplexAmber.inpcrd

echo "\n#=#=# CHECK ComplexAmber_GMX.top & ComplexAmber_GMX.gro" >> diff_out.log
diff ComplexAmber_GMX.top ../Data/ComplexAmber_GMX.top >> diff_out.log
diff ComplexAmber_GMX.gro ../Data/ComplexAmber_GMX.gro >> diff_out.log

# Run EM and MD
${grompp} -f EM.mdp -c ComplexAmber_GMX.gro -p ComplexAmber_GMX.top -o em.tpr
${mrun} ${mdrun} -v -deffnm em
${grompp} -f MD.mdp -c em.gro -p ComplexAmber_GMX.top -o md.tpr
${mrun} ${mdrun} -v -deffnm md
# vmd md.gro md.trr

echo "############################"
echo "####### DIFF SUMMARY #######"
echo "############################"
cat diff_out.log
