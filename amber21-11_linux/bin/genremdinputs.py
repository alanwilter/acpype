#!/home/conda/feedstock_root/build_artifacts/ambertools_1635195478443/_h_env_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_/bin/python


from itertools import product
import os, sys, re, random

__author__ = 'Vinicius Wilian D Cruzeiro'
__email__ = 'vwcruzeiro@ufl.edu'
__version__ = '1.0'

try:
    from argparse import ArgumentParser
except ImportError:
    if sys.version_info < (2, 7):
        raise ImportError('%s requires Python 2.7 or later' %
                          (os.path.split(sys.argv[0])[1]))

# We start by setting up the text to be printed when the user types the --help

parser = ArgumentParser(epilog='''This program generates the input files for any REMD simulations (including
                                  MultiD-REMD). It generates: the groupfile, mdin files, and the -remd-file.''', usage='%(prog)s [Options]')
parser.add_argument('-v', '--version', action='version', version='%s: %s' %
                    (parser.prog, __version__), help='''show the program's version and exit''')
parser.add_argument('--author', action='version', version='%s author: %s (E-mail: %s )' %
                    (parser.prog, __author__,__email__), help='''show the program's author name and exit''')
parser.add_argument('-O', '--overwrite', dest='overwrite', action='store_const',
                   help='''Allow existing outputs to be overwritten. Default: False''',
                   const=True, default=False)
group = parser.add_argument_group('Required Arguments')
group.add_argument('-inputs', dest='inputs', metavar='FILE', nargs='*',
                   help='''Input files containing pH, Redox Potential, Temperature, or Hamiltonian values. Each file must state
                   the type of exchange on the first line (same as in the exch_type flag of the remd-file for M-REMD
                   simulations), a description in the second line, and all variable values on the following lines (one value per
                   line). As the number of replicas on each REMD dimension needs to be even, the number of values needs to be even.''',
                   default=None)
group.add_argument('-groupfile', dest='groupfiles', metavar='FILE', nargs='*',
                   help='''Reference groupfiles. Each reference groupfile must contain only a single block referring to
                   a single replica. In this block the replica number must be replaced by REPNUM (the program will replace
                   this flag later in order to create the whole groupfile file). If doing a REMD simulation with the Hamiltonian
                   dimension, the prmtop file name may be replaced by the same flag entered in the first line of the hamiltonian
                   file given in the -inputs flag. The reference groupfiles must be given in the same order as their corresponding
                   reference mdin files.''',
                   default=None)
group.add_argument('-i', dest='mdins', metavar='FILE', nargs='*',
                   help='''Reference mdin files. Each reference mdin file must contain the variable(s) being exchanged replaced by
                   the same flag entered in the first line of the file given in the -inputs flag. Examples: solvph=PH,
                   solve=REDOX, temp0=TEMPERATURE. In order to insure each replica has a different random seed (recommended) you must place
                   ig=RANDOMNUM. The reference mdin
                   files must be given in the same order as their corresponding reference groupfiles.''',
                   default=None)
group = parser.add_argument_group('Non-required Arguments')
group.add_argument('-disang', dest='disangs', metavar='FILE', nargs='*',
help='''Reference DISANG files. If doing Umbrella Sampling in a Hamiltonian dimension, the target distance or angle must be
replaced by the same flag entered in the first line of the file given in the -inputs flag. Example: r1=-99.0, r2=HAMILTONIAN,
r3=HAMILTONIAN, r4=99.0. The reference disang files must be given in the same order as their corresponding reference mdin files.''',
default=None)
group.add_argument('-remd-file-out', dest='remdfileout', metavar='FILE',
                   help='''Name of the -remd-file output file. Default: remd.dim.[REMD dimensions types]remd''',
                   type=str, default=None)
group.add_argument('-randomseed', dest='randomseed', metavar='INTEGER',
                   help='''Seed for the random number generator. Default: 10''',
                   type=int, default=10)
group.add_argument('-nosort', dest='nosort', action='store_const',
                   help='''If stated, the replica ordering per dimension will not be sorted. If not stated, sorting will be done
                   if the input values are float or integer.''',
                   const=True, default=False)
group.add_argument('-verbose','--verbose', dest='verbose', action='store_const',
                   help='''If stated, prints information on the screen while the program is executed.''',
                   const=True, default=False)

def ErrorExit(text):
    """
    This function prints an error and kills the execution of the program
    """

    print(('\nERROR: %s' % text))
    print ('       The execution of genremdinputs.py stopped')
    sys.exit(0)

def isfloat(value):
    """
    This function checks if a value is a float number
    """

    try:
        float(value)
        return True
    except ValueError:
        return False

def isint(value):
    """
    This function checks if a value is an integer number
    """

    try:
        int(value)
        return True
    except ValueError:
        return False

def main(opt):
    """
    This is the main function to execute the program
    """

    # Number of REMD dimensions
    if (not opt.inputs):
        ErrorExit('Please provide the input file(s) (using the flag -inputs). For help type: genremdinputs.py --help')
    else:
        ndims = len(opt.inputs)
    if (opt.verbose): print("The number of REMD dimensions is %d"%(ndims))

    # Checking required arguments
    if (not opt.groupfiles):
        ErrorExit('Please provide the groupfile(s) (using the flag -groupfile). For help type: genremdinputs.py --help')
    else:
        for i in range(len(opt.groupfiles)):
            if (not os.path.exists(opt.groupfiles[i])): ErrorExit('Groupfile %s does not exist!'%(opt.groupfiles[i]))
    if (not opt.mdins):
        ErrorExit('Please provide the mdin file(s) (using the flag -i). For help type: genremdinputs.py --help')
    else:
        for i in range(len(opt.mdins)):
            if (not os.path.exists(opt.mdins[i])): ErrorExit('Mdin file %s does not exist!'%(opt.mdins[i]))
    if (len(opt.groupfiles) != len(opt.mdins)): ErrorExit('The number of group files provided (%d) do not match the number of mdin files provided (%d)'%(len(opt.groupfiles),len(opt.mdins)))
    if (opt.disangs):
        for i in range(len(opt.disangs)):
            if (not os.path.exists(opt.disangs[i])): ErrorExit('DISANG file %s does not exist!'%(opt.disangs[i]))
        if (len(opt.disangs) != len(opt.mdins)): ErrorExit('The number of DISANG files provided (%d) do not match the number of mdin files provided (%d)'%(len(opt.disangs),len(opt.mdins)))
    # Number of reference groupfiles
    nrefgrpfls = len(opt.groupfiles)

    # Reading files from -inputs
    remd_txttypes = []
    remd_descs = []
    lists = []
    DoingHamiltonian = False
    hdims = []
    # lists will contain all values of all REMD variables
    for i in range(ndims):
        infile = open(opt.inputs[i], "r")
        cnt = 0
        vector = []
        while True:
            line = infile.readline()
            if (line):
                if (cnt == 0):
                    txt = line.strip().lower()
                    # Removing numbers from the string
                    txt = ''.join([l for l in txt if not l.isdigit()])
                    if (txt != "hamiltonian" and txt != "hremd" and txt != "temperature" and txt != "temp" and txt != "ph" and txt != "redox"):
                        ErrorExit('The input file %s does not contain a valid REMD type. \'%s\' is not a valid REMD type' % (opt.inputs[i],txt.upper()))
                    if (txt == "hamiltonian" or txt == "hremd"):
                        DoingHamiltonian = True
                        hdims.append(i)
                    remd_txttypes.append(line.strip())
                    if(opt.verbose): print("  REMD Dimension #%d: %s"%(i+1,remd_txttypes[i].upper()))
                elif (cnt == 1):
                    remd_descs.append(line.strip())
                else:
                    if(len(line.strip())>0):
                        vector.append(line.strip().split()[0])
            else:
                break
            cnt += 1
        if(len(vector)%2 != 0): ErrorExit('The number of %s values (%d) in the file %s is not even! The number of replicas on each REMD dimension needs to even'%(remd_txttypes[i],len(vector),opt.inputs[i]))
        if (opt.nosort):
            lists.append(vector)
        else:
            if (isfloat(vector[0]) or isint(vector[0])):
                lists.append(sorted(vector, key=float))
            else:
                lists.append(vector)
        if(opt.verbose):
            for cnt2 in range(len(vector)):
                print("    Value #%d: %s"%(cnt2+1,lists[i][cnt2]))
        infile.close()
    del infile, vector, cnt
    if(opt.verbose): print("")

    # Setting -remd-file default name, if user did not pick a name
    if (not opt.remdfileout):
        txt = ""
        for dim in range(ndims):
            txt2 = remd_txttypes[dim].lower()
            # Removing numbers from the txt2 string
            txt2 = ''.join([l for l in txt2 if not l.isdigit()])
            if (txt2 == "temperature" or txt2 == "temp"):
                txt += "t"
            elif (txt2 == "hamiltonian" or txt2 == "hremd"):
                txt += "h"
            elif (txt2 == "ph"):
                txt += "ph"
            elif (txt2 == "redox"):
                txt += "e"
        suffix= txt+"remd"
        del txt2
        opt.remdfileout = "remd.dim.%s"%(suffix)
        del suffix

    # Checking on the groupfile(s) the number of blocks and if the required pointer are there
    FoundHamiltonian = [False for dim in range(len(hdims))]
    grpfllines = []
    for i in range(nrefgrpfls):
        group_file = open(opt.groupfiles[i], "r")
        grpfllines.append(group_file.readlines())
        numblocks = 0
        FoundFlag = False
        for line in grpfllines[i]:
            # Is this a comment line or a blank line?
            if (line.lower().replace(" ","")[0] != "#" and len(line.replace(" ","")) > 1):
                if re.search("REPNUM", line): FoundFlag = True
                numblocks += 1
                if (DoingHamiltonian):
                    for j in range(len(hdims)):
                        if re.search(remd_txttypes[hdims[j]], line): FoundHamiltonian[j] = True
        if (numblocks>1): ErrorExit('The groupfile %s must contain only a single block referring to a single replica'%(opt.groupfiles[i]))
        if (not FoundFlag): ErrorExit('Pointer REPNUM not found on groupfile %s'%(opt.groupfiles[i]))
        group_file.close()
    del group_file,numblocks,FoundFlag

    # Checking on the DISANG file(s) if the required pointer are there
    if (opt.disangs):
        disanglines = []
        for i in range(nrefgrpfls):
            disang_file = open(opt.disangs[i], "r")
            disanglines.append(disang_file.readlines())
            for line in disanglines[i]:
                if (DoingHamiltonian):
                    for j in range(len(hdims)):
                        if re.search(remd_txttypes[hdims[j]], line): FoundHamiltonian[j] = True
            disang_file.close()
        del disang_file

    # Checking on the mdin file(s) if the required pointer are there
    mdinlines = []
    for i in range(nrefgrpfls):
        mdin_file = open(opt.mdins[i], "r")
        mdinlines.append(mdin_file.readlines())
        FoundFlag = [False for dim in range(ndims+1)]
        for line in mdinlines[i]:
            # Is this a comment line or a blank line?
            if (line.lower().replace(" ","")[0] != "#" and len(line.replace(" ","")) > 1):
                for dim in range(ndims):
                    flag = remd_txttypes[dim]+","
                    if re.search(flag, line.replace(" ","")):
                        FoundFlag[dim] = True
        for dim in range(ndims):
            flag = remd_txttypes[dim]
            # Removing numbers from the flag string
            txt = ''.join([l for l in flag if not l.isdigit()])
            if (not FoundFlag[dim] and txt.lower() != "hamiltonian" and txt.lower() != "hremd"):
                ErrorExit('Pointer %s not found on mdin file %s'%(flag,opt.mdins[i]))
            elif (txt.lower() == "hamiltonian" or txt.lower() == "hremd"):
                if (not FoundHamiltonian[hdims[hdims.index(dim)]] and not FoundFlag[dim]):
                    if (not opt.disangs):
                        ErrorExit('Pointer %s not found on groupfile %s or on mdin file %s'%(flag,opt.groupfiles[i],opt.mdins[i]))
                    else:
                        ErrorExit('Pointer %s not found on groupfile %s, on mdin file %s, or on DISANG file %s'%(flag,opt.groupfiles[i],opt.mdins[i],opt.disangs[i]))
        mdin_file.close()
    del mdin_file, flag

    # Number of replicas per group for each REMD dimension
    nrepspgrp = []
    for dim in range(ndims):
        nrepspgrp.append(len(lists[dim]))

    # Here is where the list of REMD replicas is made and consequently the replica numbers are defined
    prodlist = list(product(*lists))

    # Number of REMD replicas
    nreps = len(prodlist)

    # Generating the groupfiles and mdin files
    if(opt.verbose): print("New groupfiles")
    # Opening new groupfiles
    grpfiles = []
    for i in range(nrefgrpfls):
        if re.search('.ref',opt.groupfiles[i]):
            txt = opt.groupfiles[i].replace('.ref','')
        else:
            txt = opt.groupfiles[i]+".final"
        filename = txt
        if (os.path.exists(filename) and not opt.overwrite): ErrorExit('The output groupfile %s already exists. Use the flag -O if you want it to be overwritten'%(filename))
        grpfiles.append(open(filename, "w"))
        if(opt.verbose): print("  %s"%(filename))
    # Opening new DISANG files
    if (opt.disangs):
        disangfiles = [["" for rep in range(nreps)] for i in range(nrefgrpfls)]
        for i in range(nrefgrpfls):
            if re.search('.ref',opt.disangs[i]):
                txt = opt.disangs[i].replace('.ref','')
            else:
                txt = opt.disangs[i]
            for rep in range(nreps):
                filename = txt+".rep.%03d"%(rep+1)
                if (os.path.exists(filename) and not opt.overwrite): ErrorExit('The output DISANG file %s already exists. Use the flag -O if you want it to be overwritten'%(filename))
                disangfiles[i][rep] = open(filename, "w")
    # Opening new mdin files
    mdinfiles = [["" for rep in range(nreps)] for i in range(nrefgrpfls)]
    for i in range(nrefgrpfls):
        if re.search('.ref',opt.mdins[i]):
            txt = opt.mdins[i].replace('.ref','')
        else:
            txt = opt.mdins[i]
        for rep in range(nreps):
            filename = txt+".rep.%03d"%(rep+1)
            if (os.path.exists(filename) and not opt.overwrite): ErrorExit('The output mdin file %s already exists. Use the flag -O if you want it to be overwritten'%(filename))
            mdinfiles[i][rep] = open(filename, "w")
    del filename
    # Printing some info on the screen, if requested
    if(opt.verbose):
        print("")
        print("Information printed to the new groupfiles")
        dashes = "-------"
        txt    = " Rep # "
        for dim in range(ndims):
            dashes += "--------------"
            txt    += " %12s "%(remd_txttypes[dim].center(12))
        print(dashes)
        print(txt)
        print(dashes)
    # Setting random seed
    random.seed(opt.randomseed)
    # Writting files
    for rep in range(nreps):
        randomnum = int(random.uniform(1.0,10000.0))
        for i in range(nrefgrpfls):
            # Writting to groupfile
            for line in grpfllines[i]:
                line = line.replace('REPNUM','%03d'%(rep+1))
                if (DoingHamiltonian):
                    for j in range(len(hdims)):
                        line = line.replace(remd_txttypes[hdims[j]],prodlist[rep][hdims[j]])
                grpfiles[i].write(line)
            grpfiles[i].write("\n")
            # Writting to DISANG file
            if (opt.disangs):
                for line in disanglines[i]:
                    if (DoingHamiltonian):
                        for j in range(len(hdims)):
                            line = line.replace(remd_txttypes[hdims[j]],prodlist[rep][hdims[j]])
                    disangfiles[i][rep].write(line)
            # Writting to mdin
            for line in mdinlines[i]:
                for dim in range(ndims+2):
                    if (dim == ndims+1):
                        flag = "REPNUM"
                        val = str('%03d'%(rep+1))
                    elif (dim == ndims):
                        flag = "RANDOMNUM"
                        val = str(randomnum)
                    else:
                        flag = remd_txttypes[dim]
                        val = prodlist[rep][dim]
                    if re.search(flag, line):
                        line = line.replace(flag,val)
                mdinfiles[i][rep].write(line)
        if(opt.verbose):
            txt = "%6d "%(rep+1)
            for dim in range(ndims):
                txt += " %12s "%(prodlist[rep][dim].center(12))
            print(txt)
    del flag, val, randomnum
    if(opt.verbose):
        print(dashes)
        print("")
        del dashes

    # Opening the -remd-file output file
    if (os.path.exists(opt.remdfileout) and not opt.overwrite): ErrorExit('The -remd-file output file %s already exists. Use the flag -O if you want it to be overwritten'%(opt.remdfileout))
    remdfl = open(opt.remdfileout, "w")
    # Writting the -remd-file output file
    # Obs.: This was the hardest part to program! :0
    if(opt.verbose): print("Information printed to the REMD file %s"%(opt.remdfileout))
    delta = [0 for dim in range(ndims)]
    for dim in range(ndims):
        # Writting some info to the -remd-file
        remdfl.write('&multirem\n')
        remdfl.write("exch_type='%s'\n"%(''.join([l for l in remd_txttypes[dim] if not l.isdigit()])))
        # Verbose
        if(opt.verbose): print("  Groups for dimension #%d (%s)"%(dim+1,remd_txttypes[dim]))
        # Number of groups in this dimension
        ngrps = int(round(nreps/len(lists[dim])))
        delta[dim] = nreps
        for dim2 in range(dim+1):
            delta[dim] = int(round(delta[dim]/nrepspgrp[dim2]))
        # delta[dim] = numbering distance between two replicas of the same group
        grp = 1
        cnt = 0
        for i in range(int(round(ngrps/delta[dim]))):
            if (dim == 0):
                num = 1
            else:
                num = delta[dim-1]
            frep = 1+cnt*num
            # frep = First replica of a group
            for j in range(delta[dim]):
                txt = "group(%d,:)="%(grp)
                rep = frep
                for k in range(nrepspgrp[dim]):
                    txt +="%d,"%(rep)
                    rep += delta[dim]
                # Writting group to the -remd-file
                remdfl.write(txt+"\n")
                # Verbose
                if(opt.verbose): print("    ",txt)
                frep += 1
                grp += 1
            cnt += 1
        # Writting some info to the -remd-file
        remdfl.write("desc='%s'\n"%(remd_descs[dim]))
        remdfl.write('/\n')
    remdfl.close()
    del remdfl

if __name__ == '__main__':
    opt = parser.parse_args()

    # Go ahead and execute the program.
    main(opt)
    sys.exit(0)
