# Tutorial Using ACPYPE for GROMACS

#### An example illustrating how to use ACPYPE with GROMACS


## Introduction

This tutorial is to show how to prepare a system to run on GROMACS, starting
with a PDB file for a complex protein/ligand.

It is a mere proof of concept. One should be aware of how ACPYPE works.
If you have suggestions about how to improve this tutorial, please send a
comment to alanwilter@gmail.com.

**NB:** Besides **acpype**, **antechamber** and **babel**, you will need GROMACS, which
comes with AMBER force fields now.


```bash
cd
git clone https://github.com/alanwilter/acpype.git
alias acpype='~/acpype/acpype/acpype.py'
tar xfz amber18.tar.gz
source amber18/amber.sh
```

## Getting GROMACS

Install [GROMACS](http://www.gromacs.org/).
Something like:

  * `conda install -c bioconda gromacs # if you use conda, or`

  * `sudo apt-get install gromacs # if you use Ubuntu Linux, or`

  * `brew install gromacs # if you use Mac`

should do the trick.

## Running an Example

This is for protein 1BVG.pdb (get it at [PDB](http://www.pdb.org)), a homodimer
(HIV protease) with a ligand called DMP. We will use force field Amber99SB.

Luckily, this pdb file has all hydrogens for the ligand, which is necessary for
**antechamber**. One can use either, e.g., `babel -h _mol_w/o_H_.pdb _mol_with_H.pdb`
or [YASARA View](http://www.yasara.org) to automatically add missing hydrogens to
your compound. The former just puts 'H' for atom names while the latter puts
more meaningful atom name, e.g., 'HCA' for a H bonded to a CA and not a simply
'H' as **babel** does.

In a script-like way:
```bash
mkdir acpype_tutorial
cd acpype_tutorial

# Assuming Complex.pdb (= 1BVG.pdb), split it in Protein.pdb and Ligand.pdb
wget http://www.ebi.ac.uk/pdbe/entry-files/download/pdb1bvg.ent -O 1BVG.pdb

grep 'ATOM  ' 1BVG.pdb>| Protein.pdb
grep 'HETATM' 1BVG.pdb>| Ligand.pdb

cp Protein.pdb ProteinAmber.pdb

# Process with pdb2gmx and define water
gmx pdb2gmx -ff amber99sb -f ProteinAmber.pdb -o Protein2.pdb -p Protein.top -water spce -ignh

antechamber -i Ligand.pdb -o Ligand.mol2 -fi pdb -fo mol2 -c gas

# Generate Ligand topology file with acpype (GAFF)
acpype -di Ligand.mol2 -c gas

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
cutoff-scheme   = Verlet    ; Buffered neighbor searching
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
nstxout-compressed      = 2        ; save compressed coordinates every 20 fs
nstxout                 = 0        ; save coordinates
nstvout                 = 0        ; save velocities
nstenergy               = 10       ; save energies every 20 fs
nstlog                  = 10       ; update log file every 1.0 ps
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

echo 15| gmx genion -s Complex_b4ion.tpr -o Complex_b4em.pdb -neutral -conc 0.15 -p Complex_ion.top

mv Complex_ion.top Complex.top

# Run minimisaton
gmx grompp -f em.mdp -c Complex_b4em.pdb -p Complex.top -o em.tpr
gmx_d mdrun -v -deffnm em

# Run a short simulation
gmx grompp -f md.mdp -c em.gro -p Complex.top -o md.tpr # -r em.gro
gmx_d mdrun -v -deffnm md

# Create vmd.tcl file
cat << EOF >| vmd.tcl
display projection Orthographic
display rendermode GLSL
mol modselect 0 0 protein
mol modstyle 0 0 NewCartoon 0.300000 10.000000 4.100000 0
mol modcolor 0 0 Structure
mol addrep 0
mol modselect 1 0 noh and resname DMP
mol modstyle 1 0 Licorice 0.300000 12.000000 12.000000
mol modcolor 1 0 Type
color Type C white
mol addrep 0
mol modselect 2 0 noh same resid as within 8 of resname DMP
mol modstyle 2 0 Licorice 0.200000 12.000000 12.000000
mol addrep 0
mol modselect 3 0 noh same resid as within 8 of resname DMP
mol representation Licorice 0.200000 12.000000 12.000000
mol modstyle 3 0 HBonds 3.000000 20.000000 1.000000
mol modcolor 3 0 ColorID 0
mol modstyle 3 0 HBonds 3.000000 20.000000 6.000000
mol modcolor 3 0 ColorID 4
mol smoothrep 0 0 5
mol smoothrep 0 1 5
mol smoothrep 0 2 5
mol smoothrep 0 3 5
EOF

# Visualise with VMD
vmd md.gro md.xtc -e vmd.tcl
```

Voila!
