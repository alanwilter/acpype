# This module is for calculating the free energy using TI with sander/pmemd


import os

def get_lamadas(windows):
    if windows == 12:
        lamadas = [0.00922, 0.04794, 0.11505, 0.20634, 0.31608, 0.43738,
                   0.56262, 0.68392, 0.79366, 0.88495, 0.95206, 0.99078]
    elif windows == 9:
        lamadas = [0.01592, 0.08198, 0.19331, 0.33787, 0.50000, 0.66213,
                   0.80669, 0.91802, 0.98408]
    elif windows == 7:
        lamadas = [0.02544, 0.12923, 0.29707, 0.50000, 0.70292, 0.87076,
                   0.97455]
    elif windows == 5:
        lamadas = [0.04691, 0.23076, 0.50000, 0.76923, 0.95308]
    elif windows == 3:
        lamadas = [0.11270, 0.50000, 0.88729]
    elif windows == 2:
        lamadas = [0.21132, 0.78867]
    elif windows == 1:
        lamadas = [0.50000]
    return lamadas

def get_weights(windows):
    if windows == 12:
        weights = [0.02359, 0.05347, 0.08004, 0.10158, 0.11675, 0.12457, 0.12457,
                   0.11675, 0.10158, 0.08004, 0.05347, 0.02359]
    elif windows == 9:
        weights = [0.04064, 0.09032, 0.13031, 0.15617, 0.16512, 0.15617, 0.13031,
                   0.09032, 0.04064]
    elif windows == 7:
        weights = [0.06474, 0.13985, 0.19091, 0.20897, 0.19091, 0.13985, 0.06474]
    elif windows == 5:
        weights = [0.11846, 0.23931, 0.28444, 0.23931, 0.11846]
    elif windows == 3:
        weights = [0.27777, 0.44444, 0.27777]
    elif windows == 2:
        weights = [0.5, 0.5]
    elif windows == 1:
        weights = [1.0]
    return weights

def print_message(type, steps, note, window=0):
    if type.lower() in ['min', 'heating', 'equ']:
        print("Perform %s %d steps (%s)..." %(type.upper(), steps, note))
    elif type.lower() in ['tot', 'vdw', 'chg']:
        print("Perform window %d of %s TI %d steps (%s)" %(window, type.upper(), steps, note))

def get_results(ti_windows, output_name, ti_sample_steps):
    weights = get_weights(ti_windows)
    dG = 0.0
    for i in range(1, ti_windows+1):
        dvdl = os.popen("grep -A 10 -i '%s, AVERAGES OVER     500 STEPS' %s | grep -i '%s  =' "
                        "| tail -n %d | "
                        "awk 'BEGIN {sum=0} {sum+=$3} END {print sum/NR}'"
                        %("DV\/DL", str(i)+output_name, "DV\/DL", ti_sample_steps/500)).read()
        dvdl = float(dvdl) * weights[i-1]
        dG += dvdl
    dG = round(dG, 2)
    return dG

#------------------------------------------------------------------------------
#                               Input file template
#------------------------------------------------------------------------------

def print_md_inputf(file_name, nxt, window_steps, lamada):

    ti_mdf = open(file_name, 'w')

    if nxt in ['min']:
        print("min for deltaG calculation", file=ti_mdf)
    elif nxt in ['nvt']:
        print("nvt for deltaG calculation", file=ti_mdf)
    elif nxt in ['npt']:
        print("npt for deltaG calculation", file=ti_mdf)
    elif nxt in ['ti']:
        print("ti for deltaG calculation", file=ti_mdf)

    print(" &cntrl", file=ti_mdf)

    if nxt in ['min']:
        print("  imin=1,", file=ti_mdf)
        print("  irest=0,", file=ti_mdf)
        print("  ntx=1,", file=ti_mdf)
        print("  ntmin=2,", file=ti_mdf)
        print("  maxcyc=%d," %window_steps, file=ti_mdf)
    elif nxt in ['nvt', 'npt', 'ti']:
        print("  imin=0,", file=ti_mdf)
        print("  irest=0,", file=ti_mdf)
        print("  ntx=1,", file=ti_mdf)
        print("  ig=-1,", file=ti_mdf)
        print("  nstlim=%d," %window_steps, file=ti_mdf)
        print("  dt=0.001,", file=ti_mdf)

    print("  cut=10.0,", file=ti_mdf)

    if nxt == 'min':
        print("  ntb=1,", file=ti_mdf)
    elif nxt == ['nvt', 'ti']:
        print("  ntb=1,", file=ti_mdf)
        print("  tempi=0.0,", file=ti_mdf)
        print("  temp0=300.0,", file=ti_mdf)
        print("  ntt=3,", file=ti_mdf)
        print("  gamma_ln=5.0,", file=ti_mdf)
    elif nxt in ['npt']:
        print("  ntp=1,", file=ti_mdf)
        print("  pres0=1.01325,", file=ti_mdf)
        print("  tempi=300.0,", file=ti_mdf)
        print("  temp0=300.0,", file=ti_mdf)
        print("  ntt=3,", file=ti_mdf)
        print("  gamma_ln=5.0,", file=ti_mdf)

    print("  ntc=1,", file=ti_mdf)
    print("  ntf=1,", file=ti_mdf)

    print("  ntpr=500,", file=ti_mdf)
    print("  ntwr=500,", file=ti_mdf)
    print("  ntwx=10000,", file=ti_mdf)
    print("  ntave=500,", file=ti_mdf)
    print("  ioutfm=1,", file=ti_mdf)
    if nxt in ['nvt', 'npt', 'ti']:
        print("  ntwv=-1,", file=ti_mdf)
    print("  iwrap=1,", file=ti_mdf)

    print("  icfe=1,", file=ti_mdf)
    print("  clambda=%6.5f," %lamada, file=ti_mdf)
    ti_mdf.close()

#------------------------------------------------------------------------------
#                               1 step TI using pmemd
#------------------------------------------------------------------------------

def print_pmemd_1step_mdf(file_name, nxt, window_steps, lamada, fwdorbwd):

    print_md_inputf(file_name, nxt, window_steps, lamada)

    ti_mdf = open(file_name, 'a')
    if fwdorbwd == 1:
        print("  ifsc=1,", file=ti_mdf)
        print("  timask1 = \":1\",", file=ti_mdf)
        print("  scmask1 = \":1\",", file=ti_mdf)
        print("  timask2 = \":2\",", file=ti_mdf)
        print("  scmask2 = \":2\",", file=ti_mdf)
    elif fwdorbwd == -1:
        print("  ifsc=1,", file=ti_mdf)
        print("  timask1 = \":2\",", file=ti_mdf)
        print("  scmask1 = \":2\",", file=ti_mdf)
        print("  timask2 = \":1\",", file=ti_mdf)
        print("  scmask2 = \":1\",", file=ti_mdf)
    print("/", file=ti_mdf)
    ti_mdf.close()

def OneStep_pTI(ti_windows, ti_window_steps, ti_sample_steps, exe,
                ti_prmtop, ti_inpcrd, ti_min_steps, ti_nvt_steps,
                ti_npt_steps, rev):

    prog = exe.split()[-1]

    print("****Perform TI calculation using one step method and %s program" %prog)

    #Get lamadas
    lamadas = get_lamadas(ti_windows)

    print("TOT lamadas are: ", lamadas)

    # 1. Total appearing--------------------------------------------------------
    for i in range(1, ti_windows+1):
        lamada = lamadas[i-1]
        print_pmemd_1step_mdf(str(i)+'_fwd.in', 'ti', ti_window_steps, lamada, 1)

        if i == 1:
            # 1.1 Minimization
            print_message('min', ti_min_steps, 'preparation')
            print_pmemd_1step_mdf('ti_min.in', 'min', ti_min_steps, lamada, 1)
            os.system("%s -O -i ti_min.in -o ti_min.out -p %s -c %s -r ti_min.rst -x ti_min.netcdf" %(exe, ti_prmtop, ti_inpcrd))

            # 1.2 NVT heating
            print_message('heating', ti_nvt_steps, 'preparation')
            print_pmemd_1step_mdf('ti_nvt.in', 'nvt', ti_nvt_steps, lamada, 1)
            os.system("%s -O -i ti_nvt.in -o ti_nvt.out -p %s -c ti_min.rst -r ti_nvt.rst -x ti_nvt.netcdf" %(exe, ti_prmtop))

            # 1.3 NPT equil
            print_message('equ', ti_npt_steps, 'preparation')
            print_pmemd_1step_mdf('ti_npt.in', 'npt', ti_npt_steps, lamada, 1)
            os.system("%s -O -i ti_npt.in -o ti_npt.out -p %s -c ti_nvt.rst -r ti_npt.rst -x ti_npt.netcdf" %(exe, ti_prmtop))

            # 1.4 TI
            print_message('tot', ti_window_steps, 'forward', i)
            os.system("%s -O -i %d_fwd.in -o %d_fwd.out -p %s -c ti_npt.rst -r %d_fwd.rst -x %d_fwd.netcdf" %(exe, i, i, ti_prmtop, i, i))
        else:
            print_message('tot', ti_window_steps, 'forward', i)
            j = i - 1
            os.system("%s -O -i %d_fwd.in -o %d_fwd.out -p %s -c %d_fwd.rst -r %d_fwd.rst -x %d_fwd.netcdf" %(exe, i, i, ti_prmtop, j, i, i))

    #2. Total disappearing---------------------------------------------------
    if rev == 1:
        for i in range(1, ti_windows+1):

            lamada = lamadas[i-1]
            print_pmemd_1step_mdf(str(i)+'_bwd.in', 'ti', ti_window_steps, lamada, -1)
            print_message('tot', ti_window_steps, 'backward', i)

            if i == 1:
                os.system("%s -O -i %d_bwd.in -o %d_bwd.out -p %s -c %d_fwd.rst -r %d_bwd.rst -x %d_bwd.netcdf" %(exe, i, i, ti_prmtop, ti_windows, i, i))
            else:
                j = i - 1
                os.system("%s -O -i %d_bwd.in -o %d_bwd.out -p %s -c %d_bwd.rst -r %d_bwd.rst -x %d_bwd.netcdf" %(exe, i, i, ti_prmtop, j, i, i))

    #--------------------------------------------------------------------------
    #                             Obtain the Results
    #--------------------------------------------------------------------------

    dG1 = get_results(ti_windows, '_fwd.out', ti_sample_steps)
    if rev == 1:
        dG2 = get_results(ti_windows, '_bwd.out', ti_sample_steps)
        dG = (dG1 - dG2)/2
        dG = round(dG, 2)
    else:
        dG = dG1

    return dG

#------------------------------------------------------------------------------
#                               1 step TI using sander
#------------------------------------------------------------------------------

def print_sander_1step_mdf(file_name, nxt, window_steps, lamada, ifc4):

    print_md_inputf(file_name, nxt, window_steps, lamada)

    ti_mdf = open(file_name, 'a')
    print("  ifsc=1,", file=ti_mdf)
    print("  scmask = \":1\",", file=ti_mdf)
    if ifc4 == 1:
        print("  lj1264=1,", file=ti_mdf)
    print("/", file=ti_mdf)
    ti_mdf.close()

def OneStep_sTI(ti_windows, ti_window_steps, ti_sample_steps, minexe, exe,
                md_prmtop, md_inpcrd, md0_prmtop, md0_inpcrd,
                ti_min_steps, ti_nvt_steps, ti_npt_steps, rev, ifc4):

    prog = exe.split()[-1]

    print("****Perform TI calculation using one step method and %s program" %prog)

    # Get lamadas
    lamadas = get_lamadas(ti_windows)

    print("TOT lamadas are: ", lamadas)

    # 1. Total appearing--------------------------------------------------------
    for i in range(1, ti_windows+1):
        lamada = lamadas[i-1]

        print_sander_1step_mdf(str(i)+'_v1.in', 'ti', ti_window_steps, lamada, ifc4)
        os.system("cp %d_v1.in %d_v2.in" %(i, i))

        if i == 1:

            # 1.1 Minimization
            print_message('min', ti_min_steps, 'preparation')
 
            # Input file and group file
            print_sander_1step_mdf('ti_min1.in', 'min', ti_min_steps, lamada, ifc4)
            os.system("cp ti_min1.in ti_min2.in")
            groupf = open('ti_min.group', 'w')
            print("-O -i ti_min1.in -o ti_min1.out -p %s -c %s -inf ti_min1.info -x ti_min1.netcdf -r ti_min1.rst" %(md0_prmtop, md0_inpcrd), file=groupf)
            print("-O -i ti_min2.in -o ti_min2.out -p %s -c %s -inf ti_min2.info -x ti_min2.netcdf -r ti_min2.rst" %(md_prmtop, md_inpcrd), file=groupf)
            groupf.close()

            # Run the job
            os.system("%s -ng 2 -groupfile ti_min.group" %minexe)

            # 1.2 NVT heating
            print_message('heating', ti_nvt_steps, 'preparation')

            # Input file and group file
            print_sander_1step_mdf('ti_nvt1.in', 'nvt', ti_nvt_steps, lamada, ifc4)
            os.system("cp ti_nvt1.in ti_nvt2.in")
            groupf = open('ti_nvt.group', 'w')
            print("-O -i ti_nvt1.in -o ti_nvt1.out -p %s -c ti_min1.rst -inf ti_nvt1.info -x ti_nvt1.netcdf -r ti_nvt1.rst" %md0_prmtop, file=groupf)
            print("-O -i ti_nvt2.in -o ti_nvt2.out -p %s -c ti_min2.rst -inf ti_nvt2.info -x ti_nvt2.netcdf -r ti_nvt2.rst" %md_prmtop, file=groupf)
            groupf.close()

            # Run the job
            os.system("%s -ng 2 -groupfile ti_nvt.group" %exe)

            # 1.3 NPT equil
            print_message('equ', ti_npt_steps, 'preparation')

            # Input file and group file
            print_sander_1step_mdf('ti_npt1.in', 'npt', ti_npt_steps, lamada, ifc4)
            os.system("cp ti_npt1.in ti_npt2.in")
            groupf = open('ti_npt.group', 'w')
            print("-O -i ti_npt1.in -o ti_npt1.out -p %s -c ti_nvt1.rst -inf ti_npt1.info -x ti_npt1.netcdf -r ti_npt1.rst" %md0_prmtop, file=groupf)
            print("-O -i ti_npt2.in -o ti_npt2.out -p %s -c ti_nvt2.rst -inf ti_npt2.info -x ti_npt2.netcdf -r ti_npt2.rst" %md_prmtop, file=groupf)
            groupf.close()

            # Run the job
            os.system("%s -ng 2 -groupfile ti_npt.group" %exe)

            # 1.4 TI
            print_message('tot', ti_window_steps, 'forward', i)
            # Group file
            groupf = open(str(i)+'_fwd.group', 'w')
            print("-O -i %d_v1.in -o %d_fwd1.out -p %s -c ti_npt1.rst -inf %d_fwd1.info -x %d_fwd1.netcdf -r %d_fwd1.rst" %(i, i, md0_prmtop, i, i, i), file=groupf)
            print("-O -i %d_v2.in -o %d_fwd2.out -p %s -c ti_npt2.rst -inf %d_fwd2.info -x %d_fwd2.netcdf -r %d_fwd2.rst" %(i, i, md_prmtop, i, i, i), file=groupf)
            groupf.close()

            # Run the job
            os.system("%s -ng 2 -groupfile %d_fwd.group" %(exe, i))

        else:
            print_message('tot', ti_window_steps, 'forward', i)
            j = i - 1

            # Group file
            groupf = open(str(i)+'_fwd.group', 'w')
            print("-O -i %d_v1.in -o %d_fwd1.out -p %s -c %d_fwd1.rst -inf %d_fwd1.info -x %d_fwd1.netcdf -r %d_fwd1.rst" %(i, i, md0_prmtop, j, i, i, i), file=groupf)
            print("-O -i %d_v2.in -o %d_fwd2.out -p %s -c %d_fwd2.rst -inf %d_fwd2.info -x %d_fwd2.netcdf -r %d_fwd2.rst" %(i, i, md_prmtop, j, i, i, i), file=groupf)
            groupf.close()

            # Run the job
            os.system("%s -ng 2 -groupfile %d_fwd.group" %(exe, i))

    #2. Total disappearing---------------------------------------------------
    if rev == 1:
        for i in range(1, ti_windows+1):
            print_message('tot', ti_window_steps, 'backward', i)

            if i == 1:
                # Group file
                groupf = open(str(i)+'_bwd.group', 'w')
                print("-O -i %d_v2.in -o %d_bwd2.out -p %s -c %d_fwd2.rst -inf %d_bwd2.info -x %d_bwd2.netcdf -r %d_bwd2.rst" %(i, i, md_prmtop, ti_windows, i, i, i), file=groupf)
                print("-O -i %d_v1.in -o %d_bwd1.out -p %s -c %d_fwd1.rst -inf %d_bwd1.info -x %d_bwd1.netcdf -r %d_bwd1.rst" %(i, i, md0_prmtop, ti_windows, i, i, i), file=groupf)
                groupf.close()

                # Run the job
                os.system("%s -ng 2 -groupfile %d_bwd.group" %(exe, i))
            else:
                j = i - 1
                # Group file
                groupf = open(str(i)+'_bwd.group', 'w')
                print("-O -i %d_v2.in -o %d_bwd2.out -p %s -c %d_bwd2.rst -inf %d_bwd2.info -x %d_bwd2.netcdf -r %d_bwd2.rst" %(i, i, md_prmtop, j, i, i, i), file=groupf)
                print("-O -i %d_v1.in -o %d_bwd1.out -p %s -c %d_bwd1.rst -inf %d_bwd1.info -x %d_bwd1.netcdf -r %d_bwd1.rst" %(i, i, md0_prmtop, j, i, i, i), file=groupf)
                groupf.close()

                # Run the job
                os.system("%s -ng 2 -groupfile %d_bwd.group" %(exe, i))


    #--------------------------------------------------------------------------
    #                             Obtain the Results
    #--------------------------------------------------------------------------

    dG1 = get_results(ti_windows, '_fwd1.out', ti_sample_steps)
    if rev == 1:
        dG2 = get_results(ti_windows, '_bwd1.out', ti_sample_steps)
        dG = (dG1 - dG2)/2
        dG = round(dG, 2)
    else:
        dG = dG1

    return dG

#------------------------------------------------------------------------------
#                               2 step TI using pmemd
#------------------------------------------------------------------------------

def print_pmemd_2step_mdf(file_name, nxt, window_steps, lamada, fwdorbwd, force, ntxo=2):

    print_md_inputf(file_name, nxt, window_steps, lamada)

    ti_mdf = open(file_name, 'a')
    if fwdorbwd == 1:
        if force == 'v':
            print("  timask1 = \":1\",", file=ti_mdf)
            print("  timask2 = \":2\",", file=ti_mdf)
            print("  ifsc=1,", file=ti_mdf)
            print("  scmask1 = \":1\",", file=ti_mdf)
            print("  scmask2 = \":2\",", file=ti_mdf)
            print("  crgmask = \":1|:2\"", file=ti_mdf)
        elif force == 'c':
            print("  timask1 = \":1\",", file=ti_mdf)
            print("  timask2 = \":2\",", file=ti_mdf)
    elif fwdorbwd == -1:
        if force == 'v':
            print("  timask1 = \":2\",", file=ti_mdf)
            print("  timask2 = \":1\",", file=ti_mdf)
            print("  ifsc=1,", file=ti_mdf)
            print("  scmask1 = \":2\",", file=ti_mdf)
            print("  scmask2 = \":1\",", file=ti_mdf)
            print("  crgmask = \":1|:2\"", file=ti_mdf)
        elif force == 'c':
            print("  timask1 = \":2\",", file=ti_mdf)
            print("  timask2 = \":1\",", file=ti_mdf)

    if ntxo == 1:
        print("  ntxo = 1,", file=ti_mdf)

    print("/", file=ti_mdf)
    ti_mdf.close()

def TwoStep_pTI(ti_vdw_windows, ti_chg_windows, vdw_window_steps,
               chg_window_steps, vdw_sample_steps, chg_sample_steps,
               exe, vdw_prmtop, ti_prmtop, ti_inpcrd, ti_min_steps,
               ti_nvt_steps, ti_npt_steps, rev):

    prog = exe.split()[-1]

    print("****Perform TI calculation using two step method and %s program" %prog)

    #Get lamadas
    vdw_lamadas = get_lamadas(ti_vdw_windows)
    chg_lamadas = get_lamadas(ti_chg_windows)

    print("VDW lamadas are: ", vdw_lamadas)
    print("CHG lamadas are: ", chg_lamadas)

    #1. VDW appearing-------------------------------------------------------
    for i in range(1, ti_vdw_windows+1):
        vdw_lamada = vdw_lamadas[i-1]
        print_pmemd_2step_mdf(str(i)+'_vdw_fwd.in', 'ti', vdw_window_steps, vdw_lamada, 1, 'v')

        if i == 1:
            #1.1 Minimization
            print_message('min', ti_min_steps, 'preparation')
            print_pmemd_2step_mdf('ti_min.in', 'min', ti_min_steps, vdw_lamada, 1, 'v')
            os.system("%s -O -i ti_min.in -o ti_min.out -p %s -c %s -r ti_min.rst -x ti_min.netcdf" %(exe, ti_prmtop, ti_inpcrd))

            #1.2 NVT heating
            print_message('heating', ti_nvt_steps, 'preparation')
            print_pmemd_2step_mdf('ti_nvt.in', 'nvt', ti_nvt_steps, vdw_lamada, 1, 'v')
            os.system("%s -O -i ti_nvt.in -o ti_nvt.out -p %s -c ti_min.rst -r ti_nvt.rst -x ti_nvt.netcdf" %(exe, ti_prmtop))

            #1.3 NPT equil
            print_message('equ', ti_npt_steps, 'preparation')
            print_pmemd_2step_mdf('ti_npt.in', 'npt', ti_npt_steps, vdw_lamada, 1, 'v')
            os.system("%s -O -i ti_npt.in -o ti_npt.out -p %s -c ti_nvt.rst -r ti_npt.rst -x ti_npt.netcdf" %(exe, ti_prmtop))

            #1.4 TI
            print_message('vdw', vdw_window_steps, 'forward', i)
            os.system("%s -O -i %d_vdw_fwd.in -o %d_vdw_fwd.out -p %s -c ti_npt.rst -r %d_vdw_fwd.rst -x %d_vdw_fwd.netcdf" %(exe, i, i, ti_prmtop, i, i))
        else:
            j = i - 1
            print_message('vdw', vdw_window_steps, 'forward', i)
            os.system("%s -O -i %d_vdw_fwd.in -o %d_vdw_fwd.out -p %s -c %d_vdw_fwd.rst -r %d_vdw_fwd.rst -x %d_vdw_fwd.netcdf" %(exe, i, i, ti_prmtop, j, i, i))

    #VDW TI
    k = ti_vdw_windows+1
    print_message('vdw', vdw_window_steps, 'forward and equ', k)
    print_pmemd_2step_mdf(str(k)+'_vdw_fwd.in', 'ti', vdw_window_steps, 0.98000, 1, 'v', 1)
    os.system("%s -O -i %d_vdw_fwd.in -o %d_vdw_fwd.out -p %s -c %d_vdw_fwd.rst -r %d_vdw_fwd.rst -x %d_vdw_fwd.netcdf" %(exe,
               k, k, ti_prmtop, ti_vdw_windows, k, k))

    #transfer the rst file
    print("Transfer the rst file...") # To have the metal ions (the first two atoms) have the same coordinates
    os.system("awk 'NR<=2' %d_vdw_fwd.rst > %d_vdw_fwd_merge.rst" %(k, k))
    os.system("awk 'NR==3' %d_vdw_fwd.rst | awk '{printf( \"%%12.7f%%12.7f%%12.7f%%12.7f%%12.7f%%12.7f\\n\", $1, $2, $3, $1, $2, $3)}' >> %d_vdw_fwd_merge.rst" %(k, k))
    os.system("awk 'NR>3' %d_vdw_fwd.rst >> %d_vdw_fwd_merge.rst" %(k, k))
    #os.system('cp %d_vdw_fwd.rst %d_vdw_fwd_merge.rst' %(k, k)) # This line is for test

    #2. Charge appearing-------------------------------------------------------
    for i in range(1, ti_chg_windows+1):
        chg_lamada = chg_lamadas[i-1]
        print_pmemd_2step_mdf(str(i)+'_chg_fwd.in', 'ti', chg_window_steps, chg_lamada, 1, 'c')
        print_message('chg', chg_window_steps, 'forward', i)

        if i == 1:
            os.system("%s -O -i %d_chg_fwd.in -o %d_chg_fwd.out -p %s -c %d_vdw_fwd_merge.rst -r %d_chg_fwd.rst -x %d_chg_fwd.netcdf" %(exe, i, i, vdw_prmtop, k, i, i))
        else:
            j = i - 1
            os.system("%s -O -i %d_chg_fwd.in -o %d_chg_fwd.out -p %s -c %d_chg_fwd.rst -r %d_chg_fwd.rst -x %d_chg_fwd.netcdf" %(exe, i, i, vdw_prmtop, j, i, i))

    if rev == 1:
        #3. Charge disappearing----------------------------------------------------
        for i in range(1, ti_chg_windows+1):
            chg_lamada = chg_lamadas[i-1]

            if i == ti_chg_windows:
                print_pmemd_2step_mdf(str(i)+'_chg_bwd.in', 'ti', chg_window_steps, chg_lamada, -1 , 'c', 1)
            else:
                print_pmemd_2step_mdf(str(i)+'_chg_bwd.in', 'ti', chg_window_steps, chg_lamada, -1 , 'c')

            print_message('chg', chg_window_steps, 'backward', i)

            if i == 1:
                os.system("%s -O -i %d_chg_bwd.in -o %d_chg_bwd.out -p %s -c %d_chg_fwd.rst -r %d_chg_bwd.rst -x %d_chg_bwd.netcdf" %(exe, i, i, vdw_prmtop, ti_chg_windows, i, i))
            else:
                j = i - 1
                os.system("%s -O -i %d_chg_bwd.in -o %d_chg_bwd.out -p %s -c %d_chg_bwd.rst -r %d_chg_bwd.rst -x %d_chg_bwd.netcdf" %(exe, i, i, vdw_prmtop, j, i, i))

        #transfer the rst file
        print("Transfer the RST file...")
        os.system("awk 'NR<=2' %d_chg_bwd.rst > %d_chg_bwd_merge.rst" %(ti_chg_windows, ti_chg_windows))
        os.system("awk 'NR==3' %d_chg_bwd.rst | awk '{printf( \"%%12.7f%%12.7f%%12.7f%%12.7f%%12.7f%%12.7f\\n\", $1, $2, $3, $1, $2, $3)}' >> %d_chg_bwd_merge.rst" %(ti_chg_windows, ti_chg_windows))
        os.system("awk 'NR>3' %d_chg_bwd.rst >> %d_chg_bwd_merge.rst" %(ti_chg_windows, ti_chg_windows))

        #4. VDW disappearing----------------------------------------------------
        for i in range(1, ti_vdw_windows+1):
            vdw_lamada = vdw_lamadas[i-1]
            print_pmemd_2step_mdf(str(i)+'_vdw_bwd.in', 'ti', vdw_window_steps, vdw_lamada, -1, 'v')
            print_message('vdw', vdw_window_steps, 'backward', i)

            if i == 1:
                os.system("%s -O -i %d_vdw_bwd.in -o %d_vdw_bwd.out -p %s -c %d_chg_bwd_merge.rst -r %d_vdw_bwd.rst -x %d_vdw_bwd.netcdf" %(exe, i, i, ti_prmtop, ti_chg_windows, i, i))
            else:
                j = i - 1
                os.system("%s -O -i %d_vdw_bwd.in -o %d_vdw_bwd.out -p %s -c %d_vdw_bwd.rst -r %d_vdw_bwd.rst -x %d_vdw_bwd.netcdf" %(exe, i, i, ti_prmtop, j, i, i))

    #--------------------------------------------------------------------------
    #                             Obtain the Results
    #--------------------------------------------------------------------------

    dG1 = get_results(ti_vdw_windows, '_vdw_fwd.out', vdw_sample_steps)
    dG2 = get_results(ti_chg_windows, '_chg_fwd.out', chg_sample_steps)

    if rev == 1:
        dG3 = get_results(ti_chg_windows, '_chg_bwd.out', chg_sample_steps)
        dG4 = get_results(ti_vdw_windows, '_vdw_bwd.out', vdw_sample_steps)
        dG = (dG1 + dG2 - dG3 - dG4)/2
        dG = round(dG, 2)
    else:
        dG = dG1 + dG2
        dG = round(dG, 2)

    return dG

#------------------------------------------------------------------------------
#                               2 step TI using sander
#------------------------------------------------------------------------------

def print_sander_2step_mdf(file_name, nxt, window_steps, lamada, force, ifc4):

    print_md_inputf(file_name, nxt, window_steps, lamada)

    ti_mdf = open(file_name, 'a')
    if force == 'v':
        print("  scmask=\":1\",", file=ti_mdf)
        print("  crgmask=\":1\",", file=ti_mdf)
    elif force == 'c':
        if ifc4 == 1:
            print("  lj1264=1,", file=ti_mdf)
    print("/", file=ti_mdf)
    ti_mdf.close()

def TwoStep_sTI(ti_vdw_windows, ti_chg_windows, vdw_window_steps,
            chg_window_steps, vdw_sample_steps, chg_sample_steps,
            minexe, exe, md0_prmtop, md0_inpcrd, mdv_prmtop, mdv_inpcrd,
            md_prmtop, ti_min_steps, ti_nvt_steps, ti_npt_steps, rev, ifc4):

    prog = exe.split()[-1]

    print("****Perform TI calculation using two step method and %s program" %prog)

    #Get lamadas
    vdw_lamadas = get_lamadas(ti_vdw_windows)
    chg_lamadas = get_lamadas(ti_chg_windows)

    print("VDW lamadas are: ", vdw_lamadas)
    print("CHG lamadas are: ", chg_lamadas)

    #1. VDW appearing-------------------------------------------------------
    for i in range(1, ti_vdw_windows+1):
        vdw_lamada = vdw_lamadas[i-1]
        print_sander_2step_mdf(str(i)+'_vdw1.in', 'ti', vdw_window_steps, vdw_lamada, 'v', 0)
        os.system("cp %d_vdw1.in %d_vdw2.in" %(i, i))

        if i == 1:
            #1.1 Minimization
            print_message('min', ti_min_steps, 'preparation')
            print_sander_2step_mdf('ti_min1.in', 'min', ti_min_steps, vdw_lamada, 'v', 0)
            os.system("cp ti_min1.in ti_min2.in")

            #Group file
            groupf = open('ti_min.group', 'w')
            print("-O -i ti_min1.in -o ti_min1.out -p %s -c %s -inf ti_min1.info -x ti_min1.netcdf -r ti_min1.rst" %(md0_prmtop, md0_inpcrd), file=groupf)
            print("-O -i ti_min2.in -o ti_min2.out -p %s -c %s -inf ti_min2.info -x ti_min2.netcdf -r ti_min2.rst" %(mdv_prmtop, mdv_inpcrd), file=groupf)
            groupf.close()
            os.system("%s -ng 2 -groupfile ti_min.group" %minexe)

            #1.2 NVT heating
            print_message('heating', ti_nvt_steps, 'preparation')
            print_sander_2step_mdf('ti_nvt1.in', 'nvt', ti_nvt_steps, vdw_lamada, 'v', 0)
            os.system("cp ti_nvt1.in ti_nvt2.in")

            #Group file
            groupf = open('ti_nvt.group', 'w')
            print("-O -i ti_nvt1.in -o ti_nvt1.out -p %s -c ti_min1.rst -inf ti_nvt1.info -x ti_nvt1.netcdf -r ti_nvt1.rst" %md0_prmtop, file=groupf)
            print("-O -i ti_nvt2.in -o ti_nvt2.out -p %s -c ti_min2.rst -inf ti_nvt2.info -x ti_nvt2.netcdf -r ti_nvt2.rst" %mdv_prmtop, file=groupf)
            groupf.close()
            os.system("%s -ng 2 -groupfile ti_nvt.group" %exe)

            #1.3 NPT equil
            print_message('equ', ti_npt_steps, 'preparation')
            print_sander_2step_mdf('ti_npt1.in', 'npt', ti_npt_steps, vdw_lamada, 'v', 0)
            os.system("cp ti_npt1.in ti_npt2.in")

            #Group file
            groupf = open('ti_npt.group', 'w')
            print("-O -i ti_npt1.in -o ti_npt1.out -p %s -c ti_nvt1.rst -inf ti_npt1.info -x ti_npt1.netcdf -r ti_npt1.rst" %md0_prmtop, file=groupf)
            print("-O -i ti_npt2.in -o ti_npt2.out -p %s -c ti_nvt2.rst -inf ti_npt2.info -x ti_npt2.netcdf -r ti_npt2.rst" %mdv_prmtop, file=groupf)
            groupf.close()
            os.system("%s -ng 2 -groupfile ti_npt.group" %exe)

            #1.4 TI
            print_message('vdw', vdw_window_steps, 'forward', i)
            #Group file
            groupf = open(str(i)+'_vdw_fwd.group', 'w')
            print("-O -i %d_vdw1.in -o %d_vdw_fwd1.out -p %s -c ti_npt1.rst -inf %d_vdw_fwd1.info -x %d_vdw_fwd1.netcdf -r %d_vdw_fwd1.rst" %(i, i, md0_prmtop, i, i, i), file=groupf)
            print("-O -i %d_vdw2.in -o %d_vdw_fwd2.out -p %s -c ti_npt2.rst -inf %d_vdw_fwd2.info -x %d_vdw_fwd2.netcdf -r %d_vdw_fwd2.rst" %(i, i, mdv_prmtop, i, i, i), file=groupf)
            groupf.close()
            os.system("%s -ng 2 -groupfile %d_vdw_fwd.group" %(exe, i))
        else:
            j = i - 1
            print_message('vdw', vdw_window_steps, 'forward', i)
            #Group file
            groupf = open(str(i)+'_vdw_fwd.group', 'w')
            print("-O -i %d_vdw1.in -o %d_vdw_fwd1.out -p %s -c %d_vdw_fwd1.rst -inf %d_vdw_fwd1.info -x %d_vdw_fwd1.netcdf -r %d_vdw_fwd1.rst" %(i, i, md0_prmtop, j, i, i, i), file=groupf)
            print("-O -i %d_vdw2.in -o %d_vdw_fwd2.out -p %s -c %d_vdw_fwd2.rst -inf %d_vdw_fwd2.info -x %d_vdw_fwd2.netcdf -r %d_vdw_fwd2.rst" %(i, i, mdv_prmtop, j, i, i, i), file=groupf)
            groupf.close()
            os.system("%s -ng 2 -groupfile %d_vdw_fwd.group" %(exe, i))

    #One more step to equ the structure
    k = ti_vdw_windows + 1
    print_message('vdw', vdw_window_steps, 'forward and equ', k)
    print_sander_2step_mdf(str(k)+'_vdw1.in', 'ti', vdw_window_steps, 0.98000, 'v', 0)
    os.system("cp %d_vdw1.in %d_vdw2.in" %(k, k))

    #Group file
    j = k - 1
    groupf = open(str(k)+'_vdw_fwd.group', 'w')
    print("-O -i %d_vdw1.in -o %d_vdw_fwd1.out -p %s -c %d_vdw_fwd1.rst -inf %d_vdw_fwd1.info -x %d_vdw_fwd1.netcdf -r %d_vdw_fwd1.rst" %(k, k, md0_prmtop, j, k, k, k), file=groupf)
    print("-O -i %d_vdw2.in -o %d_vdw_fwd2.out -p %s -c %d_vdw_fwd2.rst -inf %d_vdw_fwd2.info -x %d_vdw_fwd2.netcdf -r %d_vdw_fwd2.rst" %(k, k, mdv_prmtop, j, k, k, k), file=groupf)
    groupf.close()
    os.system("%s -ng 2 -groupfile %d_vdw_fwd.group" %(exe, k))

    #2. Charge appearing-------------------------------------------------------
    for i in range(1, ti_chg_windows+1):
        chg_lamada = chg_lamadas[i-1]
        print_sander_2step_mdf(str(i)+'_chg1.in', 'ti', chg_window_steps, chg_lamada, 'c', ifc4)
        os.system("cp %d_chg1.in %d_chg2.in" %(i, i))
        print_message('chg', chg_window_steps, 'forward', i)

        if i == 1:
            #Group file
            groupf = open(str(i)+'_chg_fwd.group', 'w')
            print("-O -i %d_chg1.in -o %d_chg_fwd1.out -p %s -c %d_vdw_fwd1.rst -inf %d_chg_fwd1.info -x %d_chg_fwd1.netcdf -r %d_chg_fwd1.rst" %(i, i, mdv_prmtop, k, i, i, i), file=groupf)
            print("-O -i %d_chg2.in -o %d_chg_fwd2.out -p %s -c %d_vdw_fwd2.rst -inf %d_chg_fwd2.info -x %d_chg_fwd2.netcdf -r %d_chg_fwd2.rst" %(i, i, md_prmtop, k, i, i, i), file=groupf)
            groupf.close()
            os.system("%s -ng 2 -groupfile %d_chg_fwd.group" %(exe, i))
        else:
            j = i - 1
            #Group file
            groupf = open(str(i)+'_chg_fwd.group', 'w')
            print("-O -i %d_chg1.in -o %d_chg_fwd1.out -p %s -c %d_chg_fwd1.rst -inf %d_chg_fwd1.info -x %d_chg_fwd1.netcdf -r %d_chg_fwd1.rst" %(i, i, mdv_prmtop, j, i, i, i), file=groupf)
            print("-O -i %d_chg2.in -o %d_chg_fwd2.out -p %s -c %d_chg_fwd2.rst -inf %d_chg_fwd2.info -x %d_chg_fwd2.netcdf -r %d_chg_fwd2.rst" %(i, i, md_prmtop, j, i, i, i), file=groupf)
            groupf.close()
            os.system("%s -ng 2 -groupfile %d_chg_fwd.group" %(exe, i))

    if rev == 1:
        #3. Charge disappearing----------------------------------------------------
        for i in range(1, ti_chg_windows+1):
            print_message('chg', chg_window_steps, 'backward', i)
            if i == 1:
                #Group file
                groupf = open(str(i)+'_chg_bwd.group', 'w')
                print("-O -i %d_chg2.in -o %d_chg_bwd2.out -p %s -c %d_chg_fwd2.rst -inf %d_chg_bwd2.info -x %d_chg_bwd2.netcdf -r %d_chg_bwd2.rst" %(i, i, md_prmtop, ti_chg_windows, i, i, i), file=groupf)
                print("-O -i %d_chg1.in -o %d_chg_bwd1.out -p %s -c %d_chg_fwd1.rst -inf %d_chg_bwd1.info -x %d_chg_bwd1.netcdf -r %d_chg_bwd1.rst" %(i, i, mdv_prmtop, ti_chg_windows, i, i, i), file=groupf)
                groupf.close()
                os.system("%s -ng 2 -groupfile %d_chg_bwd.group" %(exe, i))
            else:
                j = i - 1
                #Group file
                groupf = open(str(i)+'_chg_bwd.group', 'w')
                print("-O -i %d_chg2.in -o %d_chg_bwd2.out -p %s -c %d_chg_bwd2.rst -inf %d_chg_bwd2.info -x %d_chg_bwd2.netcdf -r %d_chg_bwd2.rst" %(i, i, md_prmtop, j, i, i, i), file=groupf)
                print("-O -i %d_chg1.in -o %d_chg_bwd1.out -p %s -c %d_chg_bwd1.rst -inf %d_chg_bwd1.info -x %d_chg_bwd1.netcdf -r %d_chg_bwd1.rst" %(i, i, mdv_prmtop, j, i, i, i), file=groupf)
                groupf.close()
                os.system("%s -ng 2 -groupfile %d_chg_bwd.group" %(exe, i))

        #4. VDW disappearing----------------------------------------------------
        for i in range(1, ti_vdw_windows+1):
            print_message('vdw', vdw_window_steps, 'backward', i)
            if i == 1:
                #Group file
                groupf = open(str(i)+'_vdw_bwd.group', 'w')
                print("-O -i %d_vdw2.in -o %d_vdw_bwd2.out -p %s -c %d_chg_bwd2.rst -inf %d_vdw_bwd2.info -x %d_vdw_bwd2.netcdf -r %d_vdw_bwd2.rst" %(i, i, mdv_prmtop, ti_chg_windows, i, i, i), file=groupf)
                print("-O -i %d_vdw1.in -o %d_vdw_bwd1.out -p %s -c %d_chg_bwd1.rst -inf %d_vdw_bwd1.info -x %d_vdw_bwd1.netcdf -r %d_vdw_bwd1.rst" %(i, i, md0_prmtop, ti_chg_windows, i, i, i), file=groupf)
                groupf.close()
                os.system("%s -ng 2 -groupfile %d_vdw_bwd.group" %(exe, i))
            else:
                j = i - 1
                #Group file
                groupf = open(str(i)+'_vdw_bwd.group', 'w')
                print("-O -i %d_vdw2.in -o %d_vdw_bwd2.out -p %s -c %d_vdw_bwd2.rst -inf %d_vdw_bwd2.info -x %d_vdw_bwd2.netcdf -r %d_vdw_bwd2.rst" %(i, i, mdv_prmtop, j, i, i, i), file=groupf)
                print("-O -i %d_vdw1.in -o %d_vdw_bwd1.out -p %s -c %d_vdw_bwd1.rst -inf %d_vdw_bwd1.info -x %d_vdw_bwd1.netcdf -r %d_vdw_bwd1.rst" %(i, i, md0_prmtop, j, i, i, i), file=groupf)
                groupf.close()
                os.system("%s -ng 2 -groupfile %d_vdw_bwd.group" %(exe, i))

    #--------------------------------------------------------------------------
    #                             Obtain the Results
    #--------------------------------------------------------------------------

    dG1 = get_results(ti_vdw_windows, '_vdw_fwd1.out', vdw_sample_steps)
    dG2 = get_results(ti_chg_windows, '_chg_fwd1.out', chg_sample_steps)

    if rev == 1:
        dG3 = get_results(ti_chg_windows, '_chg_bwd1.out', chg_sample_steps)
        dG4 = get_results(ti_vdw_windows, '_vdw_bwd1.out', vdw_sample_steps)
        dG = (dG1 + dG2 - dG3 - dG4)/2
        dG = round(dG, 2)
    else:
        dG = dG1 + dG2
        dG = round(dG, 2)

    return dG

