#!/Users/runner/miniforge3/conda-bld/ambertools_1635195607525/_h_env_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehol/bin/python
# Filename: OptC4.py
"""
This is the OptC4.py program, written by Pengfei Li in Merz Research Group,
Michigan State University. It needs the OpenMM and ParmEd installed and also
needs the SciPy module.
It was designed to optimize the C4 terms for the metal complex in the protein
system.

The 12-6-4 LJ-type nonbonded model was proposed in:
** P. Li, K. M. Merz, JCTC, 2014, 10, 289-297
"""


# OpenMM Imports
import simtk.openmm as mm  #about force field
import simtk.openmm.app as app #about algorithm

# ParmEd imports
from parmed import unit as u
from parmed.amber import AmberParm
from parmed.amber.mask import AmberMask
from parmed.openmm.reporters import (StateDataReporter, NetCDFReporter,
                                     RestartReporter)

# pyMSMT Imports
from pymsmt.mol.getlist import get_blist, get_all_list
from pymsmt.mol.cal import calc_bond, calc_angle, calc_dih
from pymsmt.mol.rstio import read_rstf
from pymsmt.mol.element import Atnum, CoRadiiDict
from pymsmt.api.AmberParm import read_amber_prm
from pymsmt.title import print_title

# Other Imports
from optparse import OptionParser
import os
import sys
import numpy

#-----------------------------------------------------------------------------#
# Functions
#-----------------------------------------------------------------------------#

def get_typ_dict(typinds, typs):
    #Key is the ATOM_TYPE_INDEX, value is the AMBER_ATOM_TYPE
    typdict = {}
    for i in range(0, len(typinds)):
        if typinds[i] not in list(typdict.keys()):
            typdict[typinds[i]] = [typs[i]]
        elif typs[i] in typdict[typinds[i]]:
            continue
        else:
            typdict[typinds[i]].append(typs[i])
    return typdict

def get_rmsd(initparas):

    global idxs, mcresids2, atompairs

    #Modify the C4 terms in the prmtop file
    for i in range(0, len(idxs)):
        prmtop.parm_data['LENNARD_JONES_CCOEF'][idxs[i]] = initparas[i]

    #Overwrite the prmtop file
    prmtop.write_parm('OptC4.top')

    #Perform the OpenMM optimization
    #Use AmberParm function to transfer the topology and
    #coordinate file to the object OpenMM can use
    Ambermol = AmberParm('OptC4.top', options.cfile)

    # Create the OpenMM system
    print('Creating OpenMM System')
    if options.simupha == 'gas':
        system = Ambermol.createSystem(nonbondedMethod=app.NoCutoff)
    elif options.simupha == 'liquid':
        system = Ambermol.createSystem(nonbondedMethod=app.PME,
                                       nonbondedCutoff=8.0*u.angstroms,
                                       constraints=app.HBonds,)

    #Add restraints
    force = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    force.addGlobalParameter("k", 200.0)
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")
    #for i in range(0, len(Ambermol.atoms)):
    for i, atom_crd in enumerate(Ambermol.positions):
        #if (Ambermol.atoms[i].residue.number+1 not in mcresids2) and \
        if (i+1 not in mcresids2) and \
          (Ambermol.atoms[i].residue.name not in ['WAT', 'HOH']) and \
          (Ambermol.atoms[i].name in ['CA', 'C', 'N']):
            force.addParticle(i, atom_crd.value_in_unit(u.nanometers))
    system.addForce(force)

    # Create the integrator to do Langevin dynamics
    # Temperature of heat bath, Friction coefficient, Time step
    integrator = mm.LangevinIntegrator(300*u.kelvin, 1.0/u.picoseconds,
                                       1.0*u.femtoseconds,)

    # Define the platform to use; CUDA, OpenCL, CPU, or Reference. Or do not
    # specify the platform to use the default (fastest) platform
    # Create the Simulation object
    if options.platf == 'ref':
        platform = mm.Platform.getPlatformByName('Reference')
        sim = app.Simulation(Ambermol.topology, system, integrator, platform)
    elif options.platf == 'cpu':
        platform = mm.Platform.getPlatformByName('CPU')
        sim = app.Simulation(Ambermol.topology, system, integrator, platform)
    elif options.platf == 'cuda':
        platform = mm.Platform.getPlatformByName('CUDA')
        prop = dict(CudaPrecision=options.presn)
        sim = app.Simulation(Ambermol.topology, system, integrator, platform,
                             prop)
    elif options.platf == 'opencl':
        platform = mm.Platform.getPlatformByName('OpenCL')
        prop = dict(OpenCLPrecision=options.presn)
        sim = app.Simulation(Ambermol.topology, system, integrator, platform,
                             prop)

    # Set the particle positions
    sim.context.setPositions(Ambermol.positions)

    # Output the rst file
    restrt = RestartReporter(options.rfile, 100, write_velocities=False)
    sim.reporters.append(restrt)

    # Minimize the energy
    print('Minimizing energy ' + str(options.maxsteps) + ' steps.')
    sim.minimizeEnergy(maxIterations=options.maxsteps)

    # Overwrite the final file
    state = sim.context.getState(getPositions=True, enforcePeriodicBox=True)
    restrt.report(sim, state)

    val_aft_min = []
    crds_aft_min = read_rstf(options.rfile)
    for i in atompairs:
        if len(i) == 2:
            crd1 = crds_aft_min[i[0]-1]
            crd2 = crds_aft_min[i[1]-1]
            bond = calc_bond(crd1, crd2)
            val_aft_min.append(('bond', bond))
        elif len(i) == 3:
            crd1 = crds_aft_min[i[0]-1]
            crd2 = crds_aft_min[i[1]-1]
            crd3 = crds_aft_min[i[2]-1]
            angle = calc_angle(crd1, crd2, crd3)
            val_aft_min.append(('angle', angle))
        elif len(i) == 4:
            crd1 = crds_aft_min[i[0]-1]
            crd2 = crds_aft_min[i[1]-1]
            crd3 = crds_aft_min[i[2]-1]
            crd4 = crds_aft_min[i[3]-1]
            dih = calc_dih(crd1, crd2, crd3, crd4)
            val_aft_min.append(('dih', dih))

    valdiffs = []
    for i in range(0, len(atompairs)):
        if val_bf_min[i][0] == 'bond':
            valdiff = abs(val_aft_min[i][1] - val_bf_min[i][1]) * 1.0 / 100.0
        elif val_bf_min[i][0] == 'angle':
            valdiff = abs(val_aft_min[i][1] - val_bf_min[i][1]) * 1.0 / 2.0
        elif val_bf_min[i][0] == 'dih':
            valdiff = abs(val_aft_min[i][1] - val_bf_min[i][1])
            if (360.0 - valdiff < valdiff):
                valdiff = 360.0 - valdiff
        valdiffs.append(valdiff)

    fnldiff = numpy.sum(valdiffs)
    print(fnldiff)

    return fnldiff

    #print('Calculate the RMSD')
    # Perform the RMSD calcualtion, using ptraj in AmberTools
    #os.system("cpptraj -p OptC4.top -i OptC4_ptraj.in > OptC4_ptraj.out")

    #ptrajof = open('OptC4_rmsd.txt', 'r')
    #ln = 1
    #for line in ptrajof:
    #    if ln == 3:
    #        rmsd = float(line[12:21])
    #    ln += 1
    #ptrajof.close()

    #print('RMSD is: ', rmsd)
    #return rmsd

#-----------------------------------------------------------------------------#
# Main Program
#-----------------------------------------------------------------------------#

parser = OptionParser("Usage: OptC4.py -m amber_mask -p topology_file "
                      " -c coordinate_file -r restart_file \n"
                      "                [--maxsteps maxsteps] "
                      "[--phase simulation_phase] \n"
                      "                [--size optimization_step_size] "
                      "[--method optimization_method] \n"
                      "                [--platform device_platform] "
                      "[--model metal_complex_model]")

parser.set_defaults(simupha='gas', maxsteps=1000, stepsize=10.0, minm='bfgs',
                    platf='cpu', presn='single', model=1)

parser.add_option("-m", dest="ion_mask", type='string', help="Amber mask of "
                  "the center metal ion")
parser.add_option("-p", dest="pfile", type='string', help="Topology file")
parser.add_option("-c", dest="cfile", type='string', help="Coordinate file")
parser.add_option("-r", dest="rfile", type='string', help="Restart file")
parser.add_option("--maxsteps", dest="maxsteps", type='int', \
                  help="Maximum minimization steps performed by OpenMM "
                       "in each parameter optimization cycle. "
                       "[Default: 1000]")
parser.add_option("--phase", dest="simupha", type='string', \
                  help="Simulation phase, either gas or liquid. "
                       "[Default: gas]")
parser.add_option("--size", dest="stepsize", type='float', \
                  help="Step size chosen by the user for the C4 value "
                       "during parameter searching. [Default: 10.0]")
parser.add_option("--method", dest="minm", type='string', \
                  help="Optimization method of the C4 terms. The options "
                       "are: powell, cg, bfgs, or slsqp. [Default: bfgs] "
                       "Please check the website: "
                       "http://docs.scipy.org/doc/scipy/reference/optimize.htm"
                       " for more information if interested.")
parser.add_option("--platform", dest="platf", type='string', \
                  help="Platform used. The options are: reference, cpu, cuda "
                       "or opencl. [Default: cpu] Here we use the OpenMM "
                       "software to perform the structure minimization. "
                       "Please check OpenMM user guide for more information "
                       "if interested.")
parser.add_option("--presn", dest="presn", type='string', \
                  help="Precision used. The options are: single, mixed, or "
                       "double. This option is only valid when using the CUDA "
                       "or OpenCL platform. [Default: single]")
parser.add_option("--model", dest="model", type='int', \
                  help="The metal ion complex model chosen to calculate the "
                  "sum of unsigned average errors of bond lengths, angles, "
                  "and dihedrals (the units of them are angstrom, degree and "
                  "degree respectively while the weights of them are 1/100, "
                  "1/2 and 1 respectively). This sum is the criterion for the "
                  "optimization (with a smaller value, better the "
                  "parameters). The options are: 1 or 2. "
                  "1 means a small model (only contains the metal ion and "
                  "binding heavy atoms) while 2 means a big (contains the "
                  "metal ion and heavy atoms in the ligating residues). "
                  "[Default: 1]")
(options, args) = parser.parse_args()

# Print the title of the program
version = '1.1'
print_title('OptC4.py', version)

#lowercase the input variables
options.simupha=options.simupha.lower()
options.minm=options.minm.lower()
options.platf=options.platf.lower()

#Get the metal center information from prmtop and coordinate files
prmtop, mol, atids, resids = read_amber_prm(options.pfile, options.cfile)
mask = AmberMask(prmtop, options.ion_mask)

blist = get_blist(mol, atids)
all_list = get_all_list(mol, blist, atids, 8.0)
alist = all_list.anglist
dlist = all_list.dihlist

#Get the metal ion ids
metids = []       #Atom IDs
mettyps = []      #Amber Atom Type
mettypinds = []   #Atom Type Index
for i in mask.Selected():
    j = i + 1
    metids.append(j)
    atyp = prmtop.parm_data['AMBER_ATOM_TYPE'][i]
    mettyps.append(atyp)
    mettypind = prmtop.parm_data['ATOM_TYPE_INDEX'][i]
    mettypinds.append(mettypind)

smcids = [] #Metal site ligating atom IDs
mcresids = []  #Metal Site Residue IDs
atompairs = [] #Distance pair
for i in metids:
    crdi = mol.atoms[i].crd
    atmi = mol.atoms[i].element
    radiusi = CoRadiiDict[atmi]
    for j in atids:
        if j != i:
            crdj = mol.atoms[j].crd
            atmj = mol.atoms[j].element
            dis = calc_bond(crdi, crdj)
            radiusj = CoRadiiDict[atmj]
            radiusij = radiusi + radiusj
            if (dis <= radiusij + 0.4) and (dis >= 0.1) and (atmj != 'H'):
                smcids.append(j)
                atompairs.append((i, j))
                if mol.atoms[j].resid not in mcresids:
                    mcresids.append(mol.atoms[j].resid)

for i in alist:
    if len(i) != 3:
        raise ValueError('More than 3 atoms in one angle! ' + i)
    j = [k-1 for k in i]
    atnums = [prmtop.parm_data['ATOMIC_NUMBER'][l] for l in j]
    if (list(set(metids) & set(i)) != []) and (1 not in atnums):
        atompairs.append(i)

for i in dlist:
    if len(i) != 4:
        raise ValueError('More than 4 atoms in one dihedral! ' + i)
    j = [k-1 for k in i]
    atnums = [prmtop.parm_data['ATOMIC_NUMBER'][l] for l in j]
    if (list(set(metids) & set(i)) != []) and (1 not in atnums):
        atompairs.append(i)

#Add metal ion IDs to the Metal Site Residue IDs
mcresids2 = mcresids
for i in metids:
    mcresids2.append(mol.atoms[i].resid)
mcresids2 = list(set(mcresids2))
mcresids2.sort()

print('Residues in the metal site: ', mcresids2)

mcids = []  #Metal site atom IDs
if options.model == 1: #Small model
    mcids = smcids
elif options.model == 2: #Big model
    for i in mcresids:
        for j in mol.residues[i].resconter:
            if mol.atoms[j].element != 'H':
                mcids.append(j)

#Calculate the distances between metal ion and ligating atoms
val_bf_min = []
crds_bf_min = read_rstf(options.cfile)
for i in atompairs:
    if len(i) == 2:
        crd1 = crds_bf_min[i[0]-1]
        crd2 = crds_bf_min[i[1]-1]
        bond = calc_bond(crd1, crd2)
        val_bf_min.append(('bond', bond))
    elif len(i) == 3:
        crd1 = crds_bf_min[i[0]-1]
        crd2 = crds_bf_min[i[1]-1]
        crd3 = crds_bf_min[i[2]-1]
        angle = calc_angle(crd1, crd2, crd3)
        val_bf_min.append(('angle', angle))
    elif len(i) == 4:
        crd1 = crds_bf_min[i[0]-1]
        crd2 = crds_bf_min[i[1]-1]
        crd3 = crds_bf_min[i[2]-1]
        crd4 = crds_bf_min[i[3]-1]
        dih = calc_dih(crd1, crd2, crd3, crd4)
        val_bf_min.append(('dih', dih))

#print("Bond, angle and dihedral before minimization...")
#print(val_bf_min)

#-----------------------------------------------------------------------------#
#Get the Amber mask of the metal center complex and print it into ptraj.in file
#-----------------------------------------------------------------------------#

#Print the parmed input file, add new LJ types to the bonded atoms
maskns = []
for i in smcids:
    maskn = str(mol.atoms[i].resid) + '@' + mol.atoms[i].atname
    maskns.append(maskn)

w_parmedf = open('OptC4_parmed.in', 'w')
print("add12_6_4 " + options.ion_mask, file=w_parmedf)
print("outparm %s.c4" %options.pfile, file=w_parmedf)
print("quit", file=w_parmedf)
w_parmedf.close()

os.system("parmed -O -i OptC4_parmed.in -p %s -c %s" %(options.pfile, options.cfile))

#Get the new molecule
prmtop, mol, atids, resids = read_amber_prm(options.pfile + '.c4', options.cfile)
c4terms = prmtop.parm_data['LENNARD_JONES_CCOEF']

mctyps = [] #Metal Site Atom Type
mctypinds = [] #Metal Site Atom Type Index
num = 0
for j in mcids:
    k = j - 1
    #assign new Amber atom types
    #prmtop.parm_data['AMBER_ATOM_TYPE'][k] = 'X' + str(num)
    atyp = prmtop.parm_data['AMBER_ATOM_TYPE'][k]
    mctyps.append(atyp)
    mctypind = prmtop.parm_data['ATOM_TYPE_INDEX'][k]
    mctypinds.append(mctypind)
    num = num + 1

#Get the atom type dictionary for people to see
mettypdict = get_typ_dict(mettypinds, mettyps)
mctypdict = get_typ_dict(mctypinds, mctyps)

print('The following is the dictionary of the atom types: ')
print(mettypdict)
print(mctypdict)

#Delete the repeat ATOM_TYPE_INDEX
mettypinds = list(set(mettypinds))
mettypinds.sort()
mctypinds = list(set(mctypinds))
mctypinds.sort()

#Detect the C4 terms which relates to the metal center complex
ntyps = prmtop.pointers['NTYPES']

idxs = [] #Index of C4 terms which needs to modify
iddict = {} #Dictionary of idx corresponding to ATOM_TYPE_INDEX pair
for i in mettypinds:
    j = i - 1
    for k in mctypinds:
        l = k - 1
        if j < l:
            idx = prmtop.parm_data['NONBONDED_PARM_INDEX'][ntyps*j+l] - 1
        else:
            idx = prmtop.parm_data['NONBONDED_PARM_INDEX'][ntyps*l+j] - 1
        idxs.append(idx)
        iddict[idx] = (i, k)

idxs.sort()
initparas = [c4terms[i] for i in idxs]

print('Initial C4 parameters are : ', initparas)

#Doing optimization of the parameters, initial was the normal C4 term

if options.minm == 'powell':
    from scipy.optimize import fmin_powell as fmin
    xopt = fmin(get_rmsd, initparas)
elif options.minm == 'cg':
    from scipy.optimize import fmin_cg as fmin
    xopt = fmin(get_rmsd, initparas, epsilon=options.stepsize)
elif options.minm == 'bfgs':
    from scipy.optimize import fmin_bfgs as fmin
    xopt = fmin(get_rmsd, initparas, epsilon=options.stepsize)
elif options.minm == 'slsqp':
    from scipy.optimize import fmin_slsqp as fmin
    xopt = fmin(get_rmsd, initparas, epsilon=options.stepsize)

print("Final parameters...")
print(xopt)

