#summary An example illustrating how to use ACPYPE with GROMACS
#labels Phase-Deploy
= Tutorial Using ACPYPE for GROMACS =

== Introduction ==

This tutorial is to show how to prepare a system to run on GROMACS, starting
with a PDB file for a complex protein/ligand.

It is a mere proof of concept. One should be aware of how ACPYPE is trying to do
it. If you have suggestions about how to improve this tutorial, please send a
comment to alanwilter@gmail.com.

*NB:* Besides *acpype*, *antechamber* and *babel*, you will need GROMACS, which
comes with AMBER force fields now.

== Getting GROMACS ==

Install [http://www.gromacs.org/ GROMACS].
Something like:

  * `conda install -c bioconda gromacs` # if you use conda, or

  * `sudo apt-get install gromacs` # if you use Ubuntu Linux, or

  * `brew install gromacs` # if you use Mac

should do the trick.

== Running an Example ==

This is for protein 1BVG.pdb (get it at [http://www.pdb.org PDB]), a homodimer
(HIV protease) with a ligand called DMP. We will use force field Amber99SB.

Luckily, this pdb file has all hydrogens for the ligand, which is necessary for
*antechamber*. One can use either, e.g., `babel -h _mol_w/o_H_.pdb _mol_with_H.pdb`
or [http://www.yasara.org YASARA View] to automatically add missing hydrogens to
your compound. The former just puts 'H' for atom names while the latter puts
more meaningful atom name, e.g., 'HCA' for a H bonded to a CA and not a simply
'H' as *babel* does.

In a script-like way:
{{{
# Assuming Complex.pdb (= 1BVG.pdb), split it in Protein.pdb and Ligand.pdb
wget http://www.ebi.ac.uk/pdbe/entry-files/download/pdb1bvg.ent -O 1BVG.pdb

grep 'ATOM  ' 1BVG.pdb>| Protein.pdb
grep 'HETATM' 1BVG.pdb>| Ligand.pdb

cp Protein.pdb ProteinAmber.pdb

# Process with pdb2gmx and define water
gmx pdb2gmx -ff amber99sb -f ProteinAmber.pdb -o Protein2.pdb -p Protein.top -water spce -ignh

# Generate Ligand topology file with acpype (GAFF)
acpype -di Ligand.pdb -c gas

# Merge Protein2.pdb + updated Ligand_NEW.pdb -> Complex.pdb
grep -h ATOM Protein2.pdb Ligand.acpype/Ligand_NEW.pdb >| Complex.pdb

# Edit Protein.top -> Complex.top
cp Ligand.acpype/Ligand_GMX.itp Ligand.itp
cp Protein.top Complex.top

# `#include "Ligand.itp"` has to be inserted right after
# `#include "amber99sb.ff/forcefield.itp"`
# line and before `Protein_*.itp` line in _Complex.top_.

cat Complex.top | sed '/forcefield\.itp\"/a\
#include "Ligand.itp"
' >| Complex2.top

echo "Ligand   1" >> Complex2.top
mv Complex2.top Complex.top

# Setup the box and add water
gmx editconf -bt triclinic -f Complex.pdb -o Complex.pdb -d 1.0
gmx solvate -cp Complex.pdb -cs spc216.gro -o Complex_b4ion.pdb -p Complex.top

# Create ions.mdp file
cat << EOF >| ions.mdp
; ions.mdp - used as input into grompp to generate ions.tpr
; Parameters describing what to do, when to stop and what to save
integrator  = steep         ; Algorithm (steep = steepest descent minimization)
emtol       = 1000.0        ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep      = 0.01          ; Minimization step size
nsteps      = 200           ; Maximum number of (minimization) steps to perform

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist         = 1         ; Frequency to update the neighbor list and long range forces
cutoff-scheme	= Verlet    ; Buffered neighbor searching
ns_type         = grid      ; Method to determine neighbor list (simple, grid)
coulombtype     = cutoff    ; Treatment of long range electrostatic interactions
rcoulomb        = 1.0       ; Short-range electrostatic cut-off
rvdw            = 1.0       ; Short-range Van der Waals cut-off
pbc             = xyz       ; Periodic Boundary Conditions in all 3 dimensions
EOF

# Create em.mdp file
cat << EOF >| em.mdp
; em.mdp - used as input into grompp to generate em.tpr
; Parameters describing what to do, when to stop and what to save
integrator  = steep         ; Algorithm (steep = steepest descent minimization)
emtol       = 1000.0        ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep      = 0.01          ; Minimization step size
nsteps      = 200           ; Maximum number of (minimization) steps to perform (should be 50000)

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist         = 1         ; Frequency to update the neighbor list and long range forces
cutoff-scheme   = Verlet    ; Buffered neighbor searching
ns_type         = grid      ; Method to determine neighbor list (simple, grid)
coulombtype     = PME       ; Treatment of long range electrostatic interactions
rcoulomb        = 1.0       ; Short-range electrostatic cut-off
rvdw            = 1.0       ; Short-range Van der Waals cut-off
pbc             = xyz       ; Periodic Boundary Conditions in all 3 dimensions
EOF

# Create md.mdp file
cat << EOF >| md.mdp
;define                  = -DPOSRES  ; position restrain the protein
; Run parameters
integrator              = md        ; leap-frog integrator
nsteps                  = 1000      ; 2 * 1000 = 2 ps (should be 50000: 100 ps)
dt                      = 0.002     ; 2 fs
; Output control
nstxout                 = 2       ; save coordinates every 1.0 ps
nstvout                 = 2       ; save velocities every 1.0 ps
nstenergy               = 2       ; save energies every 1.0 ps
nstlog                  = 2       ; update log file every 1.0 ps
; Bond parameters
continuation            = no        ; first dynamics run
constraint_algorithm    = lincs     ; holonomic constraints
constraints             = h-bonds   ; bonds involving H are constrained
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy
; Nonbonded settings
cutoff-scheme           = Verlet    ; Buffered neighbor searching
ns_type                 = grid      ; search neighboring grid cells
nstlist                 = 10        ; 20 fs, largely irrelevant with Verlet
rcoulomb                = 1.0       ; short-range electrostatic cutoff (in nm)
rvdw                    = 1.0       ; short-range van der Waals cutoff (in nm)
DispCorr                = EnerPres  ; account for cut-off vdW scheme
; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
pme_order               = 4         ; cubic interpolation
fourierspacing          = 0.16      ; grid spacing for FFT
; Temperature coupling is on
tcoupl                  = V-rescale             ; modified Berendsen thermostat
tc-grps                 = Protein Non-Protein   ; two coupling groups - more accurate
tau_t                   = 0.1     0.1           ; time constant, in ps
ref_t                   = 300     300           ; reference temperature, one for each group, in K
; Pressure coupling is off
pcoupl                  = no        ; no pressure coupling in NVT
; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC
; Velocity generation
gen_vel                 = yes       ; assign velocities from Maxwell distribution
gen_temp                = 300       ; temperature for Maxwell distribution
gen_seed                = -1        ; generate a random seed
EOF

# Setup ions
gmx grompp -f ions.mdp -c Complex_b4ion.pdb -p Complex.top -o Complex_b4ion.tpr
cp Complex.top Complex_ion.top

gmx genion -s ions.tpr -o 1AKI_solv_ions.gro -p topol.top -pname NA -nname CL -neutral # 13
echo 15| gmx genion -s Complex_b4ion.tpr -o Complex_b4em.pdb -neutral -conc 0.15 -p Complex_ion.top

mv Complex_ion.top Complex.top

# Run minimisaton
gmx grompp -f em.mdp -c Complex_b4em.pdb -p Complex.top -o em.tpr
gmx mdrun -v -deffnm em

# Run a short simulation
gmx grompp -f md.mdp -c em.gro -p Complex.top -o md.tpr # -r em.gro
gmx mdrun -v -deffnm md

# Visualise with VMD
vmd md.gro md.trr
}}}

Voila!
