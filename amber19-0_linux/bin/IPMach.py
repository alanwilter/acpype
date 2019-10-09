#!/opt/anaconda1anaconda2anaconda3/bin/python
from __future__ import print_function
from pymsmt.ipmach.calhfe import OneStep_sTI, TwoStep_sTI, OneStep_pTI, TwoStep_pTI
from pymsmt.ipmach.caliod import MD_simulation
from pymsmt.ipmach.tifiles import gene_topcrd
from pymsmt.title import print_title
from pymsmt.exp import *
from scipy.optimize import fmin #as fmin
import os
from optparse import OptionParser
import numpy
import math

#------------------------------------------------------------------------------
# Subprogram functions
#------------------------------------------------------------------------------
class ION:
    def __init__(self, resname, atname, attype, element, charge, rmin, ep):
        self.resname = resname
        self.atname = atname
        self.attype = attype
        self.element = element
        self.charge = charge
        self.rmin = rmin
        self.ep = ep

def get_ep(rmin):
    rmin = round(rmin, 3)
    ep = 10**(-57.36* math.exp(-2.471 * rmin))
    ep = round(ep, 8)
    return rmin, ep

def get_HFE():
    if prog == 'sander':
        if TISteps == 1:
            dG = OneStep_sTI(ti_windows, ti_window_steps, ti_sample_steps,
                             minexe, exe, md_prmtop, md_inpcrd, md0_prmtop,
                             md0_inpcrd, ti_min_steps, ti_nvt_steps,
                             ti_npt_steps, rev, ifc4)
        elif TISteps == 2:
            dG = TwoStep_sTI(ti_vdw_windows, ti_chg_windows, vdw_window_steps,
                             chg_window_steps, vdw_sample_steps,
                             chg_sample_steps, minexe, exe, md0_prmtop,
                             md0_inpcrd, mdv_prmtop, mdv_inpcrd, md_prmtop,
                             ti_min_steps, ti_nvt_steps, ti_npt_steps, rev,
                             ifc4)
    elif prog == 'pmemd':
        if TISteps == 1:
            dG = OneStep_pTI(ti_windows, ti_window_steps,
                             ti_sample_steps, exe, pti_prmtop, pti_inpcrd,
                             ti_min_steps, ti_nvt_steps, ti_npt_steps, rev)
        elif TISteps == 2:
            dG = TwoStep_pTI(ti_vdw_windows, ti_chg_windows, vdw_window_steps,
                             chg_window_steps, vdw_sample_steps,
                             chg_sample_steps, exe, ptiv_prmtop, pti_prmtop,
                             pti_inpcrd, ti_min_steps, ti_nvt_steps,
                             ti_npt_steps, rev)
    return dG

def get_IOD():
    iod, cn = MD_simulation(exe, md_prmtop, md_inpcrd, md_min_steps,
                            md_nvt_steps, md_npt_steps, md_md_steps, ifc4)
    return iod, cn

def get_hfe_iod_for_1264(params, sp=0):

    rmin, c4v = params[0], params[1]
    ion1.rmin, ion1.ep = get_ep(rmin)
    gene_topcrd(ion0, ion1, watermodel, distance, 1, c4v)

    dG = get_HFE()
    iod, cn = get_IOD()

    print("This result of this cycle:")
    print("    rmin=%5.3f (Angstrom), ep=%10.8f (Kcal/mol), c4=%5.0f (Kcal/mol*A^4)" %(ion1.rmin, ion1.ep, c4v))
    print("    HFE=%6.1f (Kcal/mol), IOD=%5.2f (Angstrom), CN=%3.1f" %(dG, iod, cn))
    print("    Exp:HFE=%6.1f (Kcal/mol), IOD=%5.2f (Angstrom)" %(HFE_VAL, IOD_VAL))

    if sp == 1:
        quit()

    HFEerr = abs(dG - HFE_VAL)
    IODerr = abs(iod - IOD_VAL)

    TOTerr = HFEerr + IODerr * 100.0

    if HFEerr <= hfetol and IODerr <= iodtol:
        print("Find the parameters!")
        print("    rmin=%5.3f (Angstrom), ep=%10.8f (Kcal/mol)" %(ion1.rmin, ion1.ep))
        print("    HFE=%6.1f (Kcal/mol), IOD=%5.2f (Angstrom), CN=%3.1f" %(dG, iod, cn))
        print("    Exp:HFE=%6.1f (Kcal/mol), IOD=%5.2f (Angstrom)" %(HFE_VAL, IOD_VAL))
        quit()

    return TOTerr

def get_hfe_for_126(rmin, sp=0):

    ion1.rmin, ion1.ep = get_ep(rmin)
    gene_topcrd(ion0, ion1, watermodel, distance)
    dG = get_HFE()

    print("This result of this cycle:")
    print("    Rmin/2=%5.3f (Angstrom), ep=%10.8f (Kcal/mol)" %(ion1.rmin, ion1.ep))
    print("    HFE=%6.1f (Kcal/mol)" %dG)

    if sp == 1:
        quit()

    HFEerr = abs(dG - HFE_VAL)

    if HFEerr <= hfetol:
        print("######################Find the parameters!####################")
        iod, cn = get_IOD()
        print("    Rmin/2=%5.3f (Angstrom), ep=%10.8f (Kcal/mol)" %(ion1.rmin, ion1.ep))
        print("    HFE=%6.1f (Kcal/mol), IOD=%5.2f (Angstrom), CN=%3.1f" %(dG, iod, cn))
        print("    Exp:HFE=%6.1f (Kcal/mol), IOD=%5.2f (Angstrom)" %(HFE_VAL, IOD_VAL))
        quit()

    return HFEerr

def get_iod_for_126(rmin, sp=0):

    ion1.rmin, ion1.ep = get_ep(rmin)
    gene_topcrd(ion0, ion1, watermodel, distance)
    iod, cn = get_IOD()

    print("This result of this cycle:")
    print("    Rmin/2=%5.3f (Angstrom), ep=%10.8f (Kcal/mol)" %(ion1.rmin, ion1.ep))
    print("    IOD=%5.2f (Angstrom), CN=%3.1f" %(iod, cn))

    if sp == 1:
        quit()

    IODerr = abs(iod - IOD_VAL)

    if IODerr <= iodtol:
        print("######################Find the parameters!####################")
        dG = get_HFE()
        print("    Rmin/2=%5.3f (Angstrom), ep=%10.8f (Kcal/mol)" %(ion1.rmin, ion1.ep))
        print("    HFE=%6.1f (Kcal/mol), IOD=%5.2f (Angstrom), CN=%3.1f" %(dG, iod, cn))
        print("    Exp:HFE=%6.1f (Kcal/mol), IOD=%5.2f (Angstrom)" %(HFE_VAL, IOD_VAL))
        quit()

    return IODerr

def get_sp_for_126(rmin):

    ion1.rmin, ion1.ep = get_ep(rmin)
    gene_topcrd(ion0, ion1, watermodel, distance)
    dG = get_HFE()
    iod, cn = get_IOD()

    print("This result of this cycle:")
    print("    Rmin/2=%5.3f (Angstrom), ep=%10.8f (Kcal/mol)" %(ion1.rmin, ion1.ep))
    print("    HFE=%6.1f (Kcal/mol), IOD=%5.2f (Angstrom), CN=%3.1f" %(dG, iod, cn))
    print("    Exp:HFE=%6.1f (Kcal/mol), IOD=%5.2f (Angstrom)" %(HFE_VAL, IOD_VAL))

#----------------------------------------------------------------------------#
#                            Main Program                                    #
#----------------------------------------------------------------------------#
parser = OptionParser("Usage: IPMach.py -i inputfile")
parser.add_option("-i", dest="inputf", type='string',
                  help="Input file name")
(options, args) = parser.parse_args()

# Print the title of the program
version = '1.0'
print_title('IPMach.py', version)

#---------------------------Default values------------------------------------

# About the program and steps of TI running
prog = 'sander'
cpus = 2
gpus = 0
mode = 'normal'
TIsteps = 2
rev = 1
ti_windows = 7
ti_vdw_windows = 3
ti_chg_windows = 7

# About the running type
opt = 0
para_set = 'HFE'
cal_type = 'OPT'
max_iternum = 100
HFE_VAL = -100.0
IOD_VAL = 2.0
hfetol = 1.0  #Hydration free energy accuray tolerance, unit kcal/mol
iodtol = 0.01 #Ion-oxygen distance accuray tolerance, unit Angstrom

# About the intial parameters
watermodel = 'TIP3P'
distance = 13.0
rmin = 1.500
c4v = 0.0
ifc4 = 0

#----------------------------Read the input files------------------------------
rinput = open(options.inputf, 'r')
for line in rinput:
    line = line.split()
    # About the residue and atname information
    if line[0].lower() == "resname":
        resname = line[1]
    elif line[0].lower() == "atname":
        atname = line[1]
    elif line[0].lower() == "element":
        element = line[1]
    elif line[0].lower() == "attype":
        attype = line[1]
    elif line[0].lower() == "charge":
        charge = round(float(line[1]), 1)
    # About the program and steps of TI running
    elif line[0].lower() == "program":
        prog = line[1].lower()
    elif line[0].lower() == "cpus":
        cpus = int(line[1])
    elif line[0].lower() == "gpus":
        gpus = int(line[1])
    elif line[0].lower() == "mode":
        mode = line[1].lower()
    elif line[0].lower() == "tisteps":
        TISteps = int(line[1])
    elif line[0].lower() == "rev":
        rev = int(line[1])
    elif line[0].lower() == "ti_windows":
        ti_windows = int(line[1])
    elif line[0].lower() == "vdw_windows":
        ti_vdw_windows = int(line[1])
    elif line[0].lower() == "chg_windows":
        ti_chg_windows = int(line[1])
    # About the running type
    elif line[0].lower() == "opt":
        opt = int(line[1])
    elif line[0].lower() == "set":
        para_set = line[1].upper()
    elif line[0].lower() == "cal_type":
        cal_type = line[1].upper()
    elif line[0].lower() == "maxiter":
        max_iternum = int(line[1])
    elif line[0].lower() == "hfe":
        HFE_VAL = round(float(line[1]), 1)
    elif line[0].lower() == "iod":
        IOD_VAL = round(float(line[1]), 2)
    elif line[0].lower() == "hfetol":
        hfetol = float(line[1])
    elif line[0].lower() == "iodtol":
        iodtol = float(line[1])
    # About the intial parameters
    elif line[0].lower() == "watermodel":
        watermodel = line[1].upper()
    elif line[0].lower() == "rmin":
        rmin = round(float(line[1]), 3)
    elif line[0].lower() == "c4":
        c4v = float(line[1])
        c4v = round(c4v, 1)
        if c4v != 0.0:
            ifc4 = 1
            if para_set != '1264':
                para_set = '1264'
                print('The nonzero c4v variable is found, so the para_set '
                      'varaible was set to 1264 automatically.')
    elif line[0].lower() == "distance":
        distance = float(line[1])
        distance = round(distance, 1)
rinput.close()

#-------------Define the residue information and modeling files----------------
ion0 = ION('M0', 'M0', 'M0', element, 0.0, 0.0, 0.0)
rmin, ep = get_ep(rmin)
ion1 = ION(resname, atname, attype, element, charge, rmin, ep)

# TI prmtop and inpcrd files for PMEMD, should have two residues merged
pti_prmtop = element + "_wat_pti.prmtop"
pti_inpcrd = element + "_wat_pti.inpcrd"
ptiv_prmtop = element + "_wat_pvdw.prmtop"

# TI prmtop and inpcrd files for sander
md0_prmtop = element + "_wat_s0.prmtop"
md0_inpcrd = element + "_wat_s0.inpcrd"

mdv_prmtop = element + "_wat_sv.prmtop"
mdv_inpcrd = element + "_wat_sv.inpcrd"

# MD prmtop and inpcrd for sander/pmemd (also use as TI prmtop file for sander)
md_prmtop = element + "_wat_svc.prmtop"
md_inpcrd = element + "_wat_svc.inpcrd"

#-------------------Setting for the program and steps use----------------------
if gpus == 1:
    exe = 'pmemd.cuda'
elif gpus > 1:
    raise ValueError('Only one gpu can be used in this case.')
else:
    if cpus == 1:
        exe = prog
        if exe == 'sander':
            raise ValueError('Could not perform sander TI calculation with only 1 cpu.')
    elif cpus > 1:
        exe = 'mpirun -np %d %s.MPI' %(cpus, prog)
        if prog == 'sander':
            minexe = 'mpirun -np 2 %s.MPI' %(prog)

if (para_set == '1264') and ('pmemd' in exe):
    raise ValueError('Could not perform 12-6-4 TI calculation with pmemd!')

if (para_set == '1264') and (c4v == 0.0):
    raise ValueError('It is not worth to do a 12-6-4 calculation without a C4 term!')

if mode == "test":
    ti_min_steps = 500
    ti_nvt_steps = 500
    ti_npt_steps = 500
    ti_window_steps = 2000
    ti_sample_steps = ti_window_steps * 3/4
    vdw_window_steps = 2000
    vdw_sample_steps = vdw_window_steps * 3/4
    chg_window_steps = 2000
    chg_sample_steps = chg_window_steps * 3/4
    md_min_steps = 500
    md_nvt_steps = 500
    md_npt_steps = 500
    md_md_steps = 2000
elif mode == "normal":
    ti_min_steps = 2000
    ti_nvt_steps = 500000
    ti_npt_steps = 500000
    ti_window_steps = 200000
    ti_sample_steps = ti_window_steps * 3/4
    vdw_window_steps = 200000
    vdw_sample_steps = vdw_window_steps * 3/4
    chg_window_steps = 200000
    chg_sample_steps = chg_window_steps * 3/4
    md_min_steps = 2000
    md_nvt_steps = 500000
    md_npt_steps = 500000
    md_md_steps = 2000000
elif mode == "scan":
    ti_min_steps = 1000
    ti_nvt_steps = 200000
    ti_npt_steps = 200000
    ti_window_steps = 80000
    ti_sample_steps = ti_window_steps * 3/4
    vdw_window_steps = 80000
    vdw_sample_steps = vdw_window_steps * 3/4
    chg_window_steps = 80000
    chg_sample_steps = chg_window_steps * 3/4
    md_min_steps = 1000
    md_nvt_steps = 200000
    md_npt_steps = 200000
    md_md_steps = 800000

#--------------------------------Print the variables---------------------------
print("The input file you are using is : %s" %options.inputf)
print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")

# Necessary inputs
try:
    print('The variable resname is : ', resname)
except:
    raise pymsmtError('resname needs to be provided.')

try:
    print('The variable atname is : ', atname)
except:
    raise pymsmtError('atname needs to be provided.')

try:
    print('The variable element is : ', element)
except:
    raise pymsmtError('element needs to be provided.')

try:
    print('The variable attype is : ', attype)
except:
    raise pymsmtError('attype needs to be provided.')

try:
    print('The variable charge is : ', charge)
except:
    raise pymsmtError('charge needs to be provided.')

# About the program and steps of TI running
print('The variable program is : ', prog)
print('The variable cpus is : ', cpus)
print('The variable gpus is : ', gpus)
print('The variable mode is : ', mode)
print('The variable tisteps is : ', TISteps)
print('The variable rev is : ', rev)
print('The variable ti_windows is : ', ti_windows)
print('The variable vdw_windows is : ', ti_vdw_windows)
print('The variable chg_windows is : ', ti_chg_windows)

# About the running type
print('The variable set is : ', para_set)
print('The variable cal_type is : ', cal_type)
print('The variable maxiter is : ', max_iternum)
print('The variable hfe is : ', HFE_VAL)
print('The variable iod is : ', IOD_VAL)
print('The variable hfetol is : ', hfetol)
print('The variable iodtol is : ', iodtol)

# About the intial parameters
print('The variable watermodel is : ', watermodel)
print('The variable rmin is : ', rmin)
print('The variable c4 is : ', c4v)
print('The variable distance is : ', distance)
print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")

#--------------------------------Doing the optimization------------------------
initparams = [rmin, c4v]

if para_set == "1264":
    if cal_type == 'OPT':
        rminc4opt = fmin(get_hfe_iod_for_1264, initparams)
    elif cal_type == 'SP':
        get_hfe_iod_for_1264(initparams, 1)
elif para_set == "HFE":
    if cal_type == 'OPT':
        rminopt = fmin(get_hfe_for_126, rmin)
    elif cal_type == 'SP':
        get_hfe_for_126(rmin, 1)
elif para_set == "IOD":
    if cal_type == 'OPT':
        rminopt = fmin(get_iod_for_126, rmin)
    elif cal_type == 'SP':
        get_iod_for_126(rmin, 1)


