#------------------------------------------------------------------------------
# This module is for IOD and CN calculation
#------------------------------------------------------------------------------

import numpy
import os

def cal_iod(rs, grs, maxind):
    rlbind = maxind - 10
    rubind = maxind + 11
    xl = rs[rlbind:rubind]
    yl = grs[rlbind:rubind]
    a, b, c = numpy.polyfit(xl, yl, 2)
    iod = -b/(2*a)
    return iod

def print_md_inputf(file_name, nxt, window_steps, ifc4):

    md_mdf = open(file_name, 'w')

    if nxt in ['min']:
        print("min for IOD calculation", file=md_mdf)
    elif nxt in ['nvt']:
        print("nvt for IOD calculation", file=md_mdf)
    elif nxt in ['npt']:
        print("npt for IOD calculation", file=md_mdf)
    elif nxt in ['md']:
        print("md for IOD calculation", file=md_mdf)

    print(" &cntrl", file=md_mdf)

    if nxt in ['min']:
        print("  imin=1,", file=md_mdf)
        print("  maxcyc=%d," %window_steps, file=md_mdf)
        ncyc = window_steps//2
        print("  ncyc=%d," %ncyc, file=md_mdf)
    elif nxt in ['nvt', 'npt', 'md']:
        print("  imin=0,", file=md_mdf)
        print("  irest=0,", file=md_mdf)
        print("  ntx=1,", file=md_mdf)
        print("  ig=-1,", file=md_mdf)
        print("  nstlim=%d," %window_steps, file=md_mdf)
        print("  dt=0.001,", file=md_mdf)

    print("  cut=10.0,", file=md_mdf)
    print("  ntc=2,", file=md_mdf)
    print("  ntf=2,", file=md_mdf)

    if nxt == 'min':
        print("  ntb=1,", file=md_mdf)
    elif nxt == 'nvt':
        print("  ntb=1,", file=md_mdf)
        print("  ntt=3,", file=md_mdf)
        print("  gamma_ln=5.0,", file=md_mdf)
        print("  tempi=0.0,", file=md_mdf)
        print("  temp0=300.0,", file=md_mdf)
    elif nxt in ['npt', 'md']:
        print("  ntp=1,", file=md_mdf)
        print("  pres0=1.01325,", file=md_mdf)
        print("  ntt=3,", file=md_mdf)
        print("  gamma_ln=5.0,", file=md_mdf)
        print("  tempi=300.0,", file=md_mdf)
        print("  temp0=300.0,", file=md_mdf)

    print("  ntpr=500,", file=md_mdf)
    print("  ntwr=500,", file=md_mdf)

    if nxt == 'md':
        ntwx = window_steps//2
        print("  ntwx=%d," %ntwx, file=md_mdf)
    else:
        print("  ntwx=2000,", file=md_mdf)

    if nxt != 'min':
        print("  ntwv=-1,", file=md_mdf)

    print("  ioutfm=1,", file=md_mdf)
    print("  iwrap=1,", file=md_mdf)

    if ifc4 == 1:
        print(" lj1264=1,", file=md_mdf)

    print("/", file=md_mdf)
    md_mdf.close()

def MD_simulation(exe, md_prmtop, md_inpcrd, md_min_steps, md_nvt_steps, md_npt_steps, md_md_steps, ifc4):

    prog = exe.split()[-1]

    #Normal MD simulation-IOD, CN
    print("****Perform MD simulation using %s program..." %prog)

    #MIN
    print("Perform md_min %d steps..." %md_min_steps)
    print_md_inputf('md_min.in', 'min', md_min_steps, ifc4)
    os.system("%s -O -i md_min.in -o md_min.out -p %s -c %s -r md_min.rst -x md_min.netcdf" %(exe, md_prmtop, md_inpcrd))

    #NVT
    print("Perform md_nvt %d steps..." %md_nvt_steps)
    print_md_inputf('md_nvt.in', 'nvt', md_nvt_steps, ifc4)
    os.system("%s -O -i md_nvt.in -o md_nvt.out -p %s -c md_min.rst -r md_nvt.rst -x md_nvt.netcdf" %(exe, md_prmtop))

    #NPT
    print("Perform md_npt %d steps..." %md_npt_steps)
    print_md_inputf('md_npt.in', 'npt', md_npt_steps, ifc4)
    os.system("%s -O -i md_npt.in -o md_npt.out -p %s -c md_nvt.rst -r md_npt.rst -x md_npt.netcdf" %(exe, md_prmtop))

    #MD
    print("Perform md_md %d steps..." %md_md_steps)
    print_md_inputf('md_md.in', 'md', md_md_steps, ifc4)
    os.system("%s -O -i md_md.in -o md_md.out -p %s -c md_npt.rst -r md_md.rst -x md_md.netcdf" %(exe, md_prmtop))

    #Cpptraj input file
    cpptrajf = open('cpptraj.in', 'w')
    print("trajin md_md.netcdf 1 100000 1", file=cpptrajf)
    print("radial M_O 0.01 5.0 :1 :WAT@O volume intrdf M_O", file=cpptrajf)
    cpptrajf.close()

    os.system("cpptraj -p %s -i cpptraj.in > cpptraj.out" %md_prmtop)
    #Get the IOD
    rdff = open('M_O', 'r')
    dr = 0.005
    rs = []
    grs = []
    ccns = []
    for line in rdff:
        if line[0] != "#":
            r, gr, ccn = line.split()
            r = float(r)
            gr = float(gr)
            ccn = float(ccn)
            rs.append(r)
            grs.append(gr)
            ccns.append(ccn)
    rdff.close()

    #Find the first peak
    maxind = grs.index(max(grs))
    iod = cal_iod(rs, grs, maxind)
    drs = [abs(i-iod) for i in rs]

    #Find the first peak again
    maxind2 = drs.index(min(drs))
    iod = cal_iod(rs, grs, maxind2)
    iod = round(iod, 2)

    #Find the second peak
    grs2 = grs[maxind+100:]
    maxind3 = grs2.index(max(grs2))
    grs3 = grs[maxind:maxind+maxind3]

    #Find the first minimum
    minind = grs3.index(min(grs3))
    cn = ccns[maxind+minind]
    cn = round(cn, 1)

    return iod, cn
