#!/opt/anaconda1anaconda2anaconda3/bin/python
from __future__ import division

import os, sys, subprocess, re, math, signal

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

parser = ArgumentParser(epilog='''This program will perform constant pH or
                        constant redox potential simulations in order to find
                        the value of Delta G reference (DELTAGREF) that gives around 50% fraction
                        of protonated or reduced species for a given residue at a given target pH
                        or redox potential. In order to run the program, you need to replace at
                        least one of the values at the STATENE flag by DELTAGREF on your cpin or
                        cein file.''', usage='%(prog)s [Options]')
parser.add_argument('-v', '--version', action='version', version='%s: %s' %
                    (parser.prog, __version__), help='''show the program's version and exit''')
parser.add_argument('--author', action='version', version='%s author: %s (E-mail: %s )' %
                    (parser.prog, __author__,__email__), help='''show the program's author information and exit''')
group = parser.add_argument_group('Required Arguments')
group.add_argument('-mdexec', dest='mdexec', metavar='FILE',
                   help='Path to the AMBER executable file. Exemple: $AMBERHOME/bin/pmemd',
                   type=str, default='mdexec')
group = parser.add_argument_group('Required Arguments - With Replica Exchange')
group.add_argument('-target', dest='target', metavar='FLOAT',
                   help='''Value of pH or Redox Potential (in Volts) that we expect to obtain a converged fraction of protonated or reduced species close to 50%%.
                           This is the target value of the pKa or Standard Redox Potential (Eo) of the system at the end of the execution.
                           Default: None''',
                   type=float, default=None)
group = parser.add_argument_group('Not-required Arguments')
group.add_argument('-do_parallel', dest='do_parallel', metavar='STRING',
                   help='Command preciding mdexec for parallel execution. Used only with Replica Exchange. Default: mpirun -np [-ng]',
                   type=str, default='do_parallel')
group.add_argument('-log', dest='log', metavar='FILE',
                   help='''When set, prints the log of the program execution to an external file (-log FILENAME). If not set, print it at the screen.
                           Default: None''',
                   type=str, default=None)
group.add_argument('-resnum', dest='resnum', metavar='INT',
                   help='''Number of the residue in which the fraction of protonated or reduced species will be monitored.
                           (REQUIRED if the number of pH or Redox titratable residues is larger than 1)''',
                   type=int, default=None)
group.add_argument('-dgrefest', dest='dgrefest', metavar='FLOAT',
                   help='''Estimated value of Delta G reference. When this flag is given, the program starts in the last phase of
                           execution, that is, on the phase of making more accurate estimatives of Delta G reference. Note: if the
                           value of -dgrefest is not close enough to the true value of Delta G reference, the execution will fail.
                           Default: None''',
                   type=float, default=None)
group.add_argument('-dgrefrange', dest='dgrefrange', metavar='FLOAT', nargs=2,
                   help='''Range of values for Delta G reference. The desired Delta G reference value has to be
                           inside this range. If -dgrefest and -dgrefrange are not given, the program will try to find a range automatically.
                           Suggestion: choose one value in which the fraction of protonated or
                           reduced species is ~ 0 and the other value in which it is ~ 1.
                           Default: None''',
                   type=float, default=None)
group.add_argument('-dginterval', dest='dginterval', metavar='FLOAT',
                   help='''When the values of the argument -dgrefrange are to be found automatically, dginterval is the interval
                           of trial values. Default: 100.0 kcal/mol''',
                   type=float, default=100.0)
group.add_argument('-maxsteps', dest='maxsteps', metavar='INT',
                   help='''Maximum number of AMBER executions. Default: 100''',
                   type=int, default=100)
group.add_argument('-fracthreshold', dest='fracthreshold', metavar='FLOAT',
                   help='''Fraction threshold. The fraction convergence criterium is: 0.5-fracthreshold/2 >= frac >=
                           0.5+fracthreshold/2. Default: 0.03''',
                   type=float, default=0.03)
group.add_argument('-noequi', dest='noequi', action='store_const',
                   help='''If stated, the equilibration simulation for a new DELTAGREF value will not be performed. Equilibration
                           runs for 10%% the number of steps of the production simulation. Default: False''',
                   const=True, default=False)
group.add_argument('-rmouts', dest='rmouts', action='store_const',
                   help='''If stated, at the end of the execution of the program, erases all output files generated by AMBER
                           (all files not stated as REQUIRED at "AMBER Arguments" below). Default: False''',
                   const=True, default=False)
group.add_argument('-bin-path', dest='binpath', metavar='FILE', required=False,
                   help='Path to the AMBER bin directory. Used to locate cphstats, cestats or fitpkaeo.py (Example: $AMBERHOME/bin ; Default: not set).',
                   type=str, default='')
group = parser.add_argument_group('AMBER Arguments - Without Replica Exchange', 'These are the arguments to be executed together with mdexec.')
group.add_argument('-i', dest='mdin', metavar='FILE', required=False,
                   help='AMBER mdin file (REQUIRED)',
                   type=str, default='mdin')
group.add_argument('-p', dest='prmtop', metavar='FILE', required=False,
                   help='AMBER parmtop file (REQUIRED)',
                   type=str, default='prmtop')
group.add_argument('-c', dest='inpcrd', metavar='FILE', required=False,
                   help='AMBER inpcrd (input coordinates) file (REQUIRED)',
                   type=str, default='inpcrd')
group.add_argument('-x', dest='mdcrd', metavar='FILE', required=False,
                   help='AMBER mdcrd (output coordinates) file',
                   type=str, default='mdcrd')
group.add_argument('-inf', dest='mdinfo', metavar='FILE', required=False,
                   help='AMBER mdinfo file',
                   type=str, default='mdinfo')
group.add_argument('-o', dest='mdout', metavar='FILE', required=False,
                   help='AMBER mdout (log) file',
                   type=str, default='mdout')
group.add_argument('-r', dest='restrt', metavar='FILE', required=False,
                   help='AMBER mdout file',
                   type=str, default='restrt')
group.add_argument('-cpin', dest='cpin', metavar='FILE', required=False,
                   help='AMBER cpin file (REQUIRED if cein file is not given)',
                   type=str, default='cpin')
group.add_argument('-cpout', dest='cpout', metavar='FILE', required=False,
                   help='AMBER cpout file',
                   type=str, default='cpout')
group.add_argument('-cprestrt', dest='cprestrt', metavar='FILE', required=False,
                   help='AMBER cprestrt file',
                   type=str, default='cprestrt')
group.add_argument('-cein', dest='cein', metavar='FILE', required=False,
                   help='AMBER cein file (REQUIRED if cpin file is not given)',
                   type=str, default='cein')
group.add_argument('-ceout', dest='ceout', metavar='FILE', required=False,
                   help='AMBER ceout file',
                   type=str, default='ceout')
group.add_argument('-cerestrt', dest='cerestrt', metavar='FILE', required=False,
                   help='AMBER cerestrt file',
                   type=str, default='cerestrt')
group.add_argument('-ref', dest='refc', metavar='FILE', required=False,
                   help='AMBER ref file',
                   type=str, default=None)
group = parser.add_argument_group('AMBER Arguments - With Replica Exchange', 'These are the arguments to be executed together with do_parallel and mdexec.')
group.add_argument('-ng', dest='ng', metavar='INT',
                   help='''Number of groups/replicas (REQUIRED)''',
                   type=int, default=None)
group.add_argument('-groupfile', dest='groupfile', metavar='FILE', required=False,
                   help='AMBER groupfile file (REQUIRED)',
                   type=str, default=None)
    
def PrintLog(opt,log,text):
    """
    This function prints a text on the screen or at the log file
    """
    if (not opt.log):
        print ('%s' % text)
    else:
        log.write('%s\n' % text)
        log.flush()

def ErrorExit(opt,log,text):
    """
    This function prints an error and kills the execution of the program
    """
    if (not opt.log):
        print ('\nERROR: %s' % text)
        print ('       The execution of finddgref.py stopped')
    else:
        log.write('\nERROR: %s\n' % text)
        log.write('       The execution of finddgref.py stopped\n')
        log.flush()
    sys.exit(0)

def RunAMBER(opt,log,data,equi):
    """
    This function execute AMBER for DELTAGREF = dg and returns the fraction of protonated or reduced species
    The total number of steps to be executed is: data.num for serial or data.num * data.nstlim for REMD
    """
    # Running total simulation with no equilibration (equi == 0)
    if (equi == 0):
        if (data.remd):
            PrintLog(opt,log," AMBER execution #%d: running %d MD steps (%d replica exchange attempts) for DELTAGREF = %.6f kcal/mol" % (data.step,data.num*data.nstlim,data.num,data.dg))
        else:
            PrintLog(opt,log," AMBER execution #%d: running %d MD steps for DELTAGREF = %.6f kcal/mol" % (data.step,data.num,data.dg))
    # Running equilibration, prior production (equi == 1)
    elif (equi == 1):
        if (data.remd):
            PrintLog(opt,log," AMBER execution #%d: running %d  MD steps (%d  replica exchange attempts) of equilibration for DELTAGREF = %.6f kcal/mol" % (data.step,int(round(float(data.num*data.nstlim)/10.0)),int(round(float(data.num)/10.0)),data.dg))
        else:
            PrintLog(opt,log," AMBER execution #%d: running %d  MD steps of equilibration for DELTAGREF = %.6f kcal/mol" % (data.step,int(round(float(data.num)/10.0)),data.dg))
    # Running production, after equilibration (equi == 2)
    elif (equi == 2):
        if (data.remd):
            PrintLog(opt,log," AMBER execution #%d: running %d MD steps (%d replica exchange attempts) of production    for DELTAGREF = %.6f kcal/mol" % (data.step,data.num*data.nstlim,data.num,data.dg))
        else:
            PrintLog(opt,log," AMBER execution #%d: running %d MD steps of production    for DELTAGREF = %.6f kcal/mol" % (data.step,data.num,data.dg))
    # Prepare temporary mdin files
    for i in range(0,opt.ng):
        mdin_file = open(data.mdins[i], 'r')
        tempfile = open(data.mdins[i]+".temporary", 'w')
        cnt = 0
        for line in mdin_file:
            if (re.search("nstlim", line) and not data.remd):
                cnt += 1
                lines = line.split(",")
                txt=""
                for j in range(len(lines)):
                    if (re.search("nstlim", lines[j])):
                        if (equi == 1):
                            txt+="nstlim=%d," % (int(round(float(data.num)/10.0)))
                        else:    
                            txt+="nstlim=%d," % (data.num)
                    elif (j<len(lines)-1):
                        txt+=lines[j]+","
                    else:
                        txt+=lines[j]
                tempfile.write(txt)
            elif (re.search("numexchg", line) and data.remd):
                cnt += 1
                lines = line.split(",")
                txt=""
                for j in range(len(lines)):
                    if (re.search("numexchg", lines[j])):
                        if (equi == 1):
                            txt+="numexchg=%d," % (int(round(float(data.num)/10.0)))
                        else:    
                            txt+="numexchg=%d," % (data.num)
                    elif (j<len(lines)-1):
                        txt+=lines[j]+","
                    else:
                        txt+=lines[j]
                tempfile.write(txt)
            else:
                tempfile.write(line)
        # In case nstlim was not given in the mdin file, we add it
        if (cnt == 0 and not data.remd):
            tempfile.close()
            mdin_file.close()
            mdin_file = open(data.mdins[i], 'r')
            tempfile = open(data.mdins[i]+".temporary", 'w')
            var = 0
            for line in mdin_file:
                if (var == 0):
                    tempfile.write(line)
                else:
                    if (equi == 1):
                        tempfile.write("nstlim=%d,\n" % (int(round(float(data.num)/10.0))))
                    else:    
                        tempfile.write("nstlim=%d,\n" % (data.num))
                    tempfile.write(line)
                    var = 0
                if re.search("&cntrl", line):
                    var = 1
        tempfile.close()
        mdin_file.close()
    
    # Prepare temporary cpin or cein files
    if (data.modeopt == 1):
        cin = opt.cpin
    else:
        cin = opt.cein
    cin_file = open(cin, 'r')
    tempfile = open(cin+".temporary", 'w')
    for line in cin_file:
        tempfile.write(line.lower().replace("deltagref","%f"%data.dg).upper())
    tempfile.close()
    cin_file.close()
    
    # Generating temporary groupfile if Replica Exchange is going to be used
    if (data.remd):
        group_file = open(opt.groupfile, "r")
        tempfile = open(opt.groupfile+".temporary", 'w')
        cnt = 0
        for line in group_file:
            if (line.lower().replace(" ","")[0] != "#" and len(line.replace(" ","")) > 1):
                lines = line.split()
                cnt1=0
                cnt2=0
                for i in range(0,len(lines)):
                    if (lines[i] == "-i"):
                        lines[i+1] += ".temporary"
                    if (equi == 0):
                        if (data.modeopt == 1):
                            if (lines[i] == "-cpin"):
                                lines[i+1] += ".temporary"
                        else:
                            if (lines[i] == "-cein"):
                                lines[i+1] += ".temporary"
                    elif (equi == 1):
                        if (data.modeopt == 1):
                            if (lines[i] == "-cpin"):
                                lines[i+1] += ".temporary"
                        else:
                            if (lines[i] == "-cein"):
                                lines[i+1] += ".temporary"
                        if (lines[i] == "-cprestrt"):
                            lines[i+1] = data.cprestrts[cnt]+".temporary"
                            cnt1+=1
                        if (lines[i] == "-cerestrt"):
                            lines[i+1] = data.cerestrts[cnt]+".temporary"
                            cnt2+=1
                    elif (equi == 2):
                        if (lines[i] == "-cpin"):
                            lines[i+1] = data.cprestrts[cnt]+".temporary"
                        if (lines[i] == "-cein"):
                            lines[i+1] = data.cerestrts[cnt]+".temporary"
                        if (lines[i] == "-cprestrt"):
                            lines[i+1] = data.cprestrts[cnt]
                            cnt1+=1
                        if (lines[i] == "-cerestrt"):
                            lines[i+1] = data.cerestrts[cnt]
                            cnt2+=1
                if (equi>0):
                    if (data.modeopt == 1 and cnt1==0):
                        if (equi==1):
                            lines.append(" -cprestrt "+data.cprestrts[cnt]+".temporary")
                        elif (equi==2):
                            lines.append(" -cprestrt "+data.cprestrts[cnt])
                    elif (data.modeopt == 2 and cnt2==0):
                        if (equi==1):
                            lines.append(" -cprestrt "+data.cprestrts[cnt]+".temporary")
                        elif (equi==2):
                            lines.append(" -cprestrt "+data.cprestrts[cnt])
                line = ""
                for i in range(0,len(lines)):
                    line += lines[i]+" "
                line += "\n"
                tempfile.write(line)
                cnt += 1
        tempfile.close()
        group_file.close()
        del cnt1, cnt2
    
    # Writting the command to execute AMBER
    if (data.remd):
        cmd = "%s %s -ng %d -groupfile %s.temporary" % (opt.do_parallel,opt.mdexec,opt.ng,opt.groupfile)
    else:
        if (equi == 0):
            crd_md = opt.inpcrd
            if (data.modeopt == 1):
                cpin_md = opt.cpin+".temporary"
                cein_md = opt.cein
            else:
                cpin_md = opt.cpin
                cein_md = opt.cein+".temporary"
            cprst_md = opt.cprestrt
            cerst_md = opt.cerestrt
        elif (equi == 1):
            crd_md = opt.inpcrd
            if (data.modeopt == 1):
                cpin_md = opt.cpin+".temporary"
                cein_md = opt.cein
            else:
                cpin_md = opt.cpin
                cein_md = opt.cein+".temporary"
            cprst_md = opt.cprestrt+".temporary"
            cerst_md = opt.cerestrt+".temporary"
        elif (equi == 2):
            crd_md = opt.restrt
            cpin_md = opt.cprestrt+".temporary"
            cein_md = opt.cerestrt+".temporary"
            cprst_md = opt.cprestrt
            cerst_md = opt.cerestrt
        cmd = "%s -O -i %s.temporary -p %s -c %s -x %s -inf %s -o %s -r %s" % (opt.mdexec,opt.mdin,opt.prmtop,crd_md,opt.mdcrd,opt.mdinfo,opt.mdout,opt.restrt)
        if (data.modeopt == 1):
            cmd += " -cpin %s -cpout %s -cprestrt %s" % (cpin_md,opt.cpout,cprst_md)
            if (os.path.exists(opt.cein)):
                cmd += " -cein %s -ceout %s -cerestrt %s" % (cein_md,opt.ceout,cerst_md)
        else:
            cmd += " -cein %s -ceout %s -cerestrt %s" % (cein_md,opt.ceout,cerst_md)
            if (os.path.exists(opt.cpin)):
                cmd += " -cpin %s -cpout %s -cprestrt %s" % (cpin_md,opt.cpout,cprst_md)
        if (opt.refc):
            cmd += " -ref %s" % opt.refc
        
    # In order to insure the program is going to run, we erase any previous mdout files
    for i in range(0,opt.ng):
        if (os.path.exists(data.mdouts[i])): os.remove(data.mdouts[i])
        
    # Executing AMBER
    FNULL = open(os.devnull, 'w')
    pipe = subprocess.Popen([cmd], stdout=FNULL, stderr=subprocess.PIPE, shell=True)
    FNULL.close()
    stderr = pipe.stderr.readline()
    
    # Check the mdout file(s) was(were) generated
    for i in range(0,opt.ng):
        if(not os.path.exists(data.mdouts[i])):
            ErrorExit(opt,log,
	    	  'The execution of AMBER using the binary %s failed. The execution returned the following STDERR:\n\n'
                      '       %s\n\n       The command executed was:\n\n       %s'
                      '\n\n       Check the mdout file (%s) to see what the error might be. Also check if all libraries '
                      'necessary to run AMBER are properly set.' % (opt.mdexec,stderr,cmd,data.mdouts[i]))
        
    # Opening mdout file(s) to check if the execution ended with success
    for i in range(0,opt.ng):
        mdout_file = open(data.mdouts[i], 'r')
        pos, lines = 21, []
        while len(lines) <= 20:
            try:
                mdout_file.seek(-pos, 2)
            except IOError:
                mdout_file.seek(0)
                break
            finally:
                lines = list(mdout_file)
            pos *= 2
        mdout_file.close()
        cnt = 0
        for line in lines:
            if re.search("Final Performance Info:", line):
                cnt += 1
                break
        if(cnt == 0):
            ErrorExit(opt,log,
	    	  'The execution of AMBER using the binary %s failed. The execution returned the following STDERR:\n\n'
                      '       %s\n\n       The command executed was:\n\n       %s'
                      '\n\n       Check the mdout file (%s) to see what the error might be. Also check if all libraries '
                      'necessary to run AMBER are properly set.' % (opt.mdexec,stderr,cmd,data.mdouts[i]))
            
    # If equi == 1, then we don't need to compute the fraction.
    if (equi != 1):
        # If Replica Exchange is going to be done, reorder the cpout or ceout files
        if (data.remd):
            if (data.modeopt == 1):
                cmd = "%scphstats -O --fix-remd reordered.cpout" % (opt.binpath.strip())
            else:
                cmd = "%scestats -O --fix-remd reordered.ceout" % (opt.binpath.strip())
            for i in range(0,opt.ng):
                if (data.modeopt == 1):
                    cmd += " "+data.cpouts[i]
                else:
                    cmd += " "+data.ceouts[i]
            pipe = subprocess.Popen([cmd], stdout=subprocess.PIPE, shell=True)
            # The next line forces the program to wait until the execution of the subprocess finishes
            lines = pipe.stdout.readlines()
                
        # Computing the fraction of protonated or reduced species
        fracs = [0.0 for rep in range(0,opt.ng)]
        for rep in range(0,opt.ng):
            if (data.remd):
                if (data.modeopt == 1):
                    cmd = "%scphstats -i %s.temporary reordered.cpout.pH_%.2f" % (opt.binpath.strip(),opt.cpin,data.solvphe[rep])
                else:
                    cmd = "%scestats -i %s.temporary reordered.ceout.E_%.5f" % (opt.binpath.strip(),opt.cein,data.solvphe[rep])
            else:
                if (data.modeopt == 1):
                    cmd = "%scphstats -i %s.temporary %s" % (opt.binpath.strip(),opt.cpin,opt.cpout)
                else:
                    cmd = "%scestats -i %s.temporary %s" % (opt.binpath.strip(),opt.cein,opt.ceout)
            pipe = subprocess.Popen([cmd], stdout=subprocess.PIPE, shell=True)
            # Counting how many fractions cphstats or cestats returns
            lines = pipe.stdout.readlines()
            cnt = 0
            for line in lines:
                if re.search("Offset", line):
                    cnt += 1
        
            if (data.modeopt == 1):
                txt1 = 'pH'
                txt2 = 'cphstats'
            else:
                txt1 = 'Redox'
                txt2 = 'cestats'
            # If cnt is zero this means we couldn't read the cphstats or cestats output
            if (cnt == 0):
                ErrorExit(opt,log,'Could not read the fraction of protonated or reduced species from the output of cphstats or cestats. '
                                  'The command executed was:\n\n       %s\n' % cmd)
            # If the argument -resnum is not given, compute the fraction only if the number of residues on cphstats or cestats is one
            elif (cnt > 1 and not opt.resnum):
                txt =  'The number of '+txt1+' titratable residues on the output of '+txt2+' is larger than 1. '
                txt += 'Please run the program again informing on the argument -resnum the number of the residue to be monitored.\n'
                txt += '       According to '+txt2+' the possible options are:'
                for line in lines:
                    if re.search("Offset", line):
                        txt += "\n         Residue '%s': to monitore this residue add the argument -resnum %d" % (line.split(":")[0].strip(),int(line.split(":")[0].strip().split(" ")[1]))
                ErrorExit(opt,log,txt)
            elif (cnt == 1):
                for line in lines:
                    if re.search("Offset", line):
                        break
                if (opt.resnum and not opt.resnum==int(line.split(":")[0].strip().split(" ")[1])):
                    ErrorExit(opt,log,'The residue number informed on -resnum (Res #%d) does not match with the residue number printed by '
                                'cphstats or cestats (Res #%d).' % (opt.resnum,int(line.split(":")[0].strip().split(" ")[1])))
                resname = line.split(":")[0].strip()
                # Reading the fraction
                if (data.modeopt == 1):
                    fracs[rep] = float(line.split("Frac Prot")[1].lstrip().split(" ")[0])
                else:
                    fracs[rep] = float(line.split("Frac Redu")[1].lstrip().split(" ")[0])
            # If the number of residues on cphstats or cestats is larger than one, compute the fraction for the residue provided at -resnum
            else:
                i = 0
                for line in lines:
                    if (re.search("Offset", line) and int(line.strip().split(" ")[1])==opt.resnum):
                        i += 1
                        break
                if (i == 0):
                    txt =  'The residue number informed on -resnum (Res #%d) is not a ' % opt.resnum
                    txt += txt1+' valid residue. Please run the program again informing correctly on the argument -resnum the number of the residue to be monitored.\n'
                    txt += '       According to '+txt2+' the possible options are:'
                    for line in lines:
                        if re.search("Offset", line):
                            txt += "\n         Residue '%s': to monitore this residue set the argument as -resnum %d" % (line.split(":")[0].strip(),int(line.split(":")[0].strip().split(" ")[1]))
                    ErrorExit(opt,log,txt)
                else:
                    resname = line.split(":")[0].strip()
                    # Reading the fraction
                    if (data.modeopt == 1):
                        fracs[rep] = float(line.split("Frac Prot")[1].lstrip().split(" ")[0])
                    else:
                        fracs[rep] = float(line.split("Frac Redu")[1].lstrip().split(" ")[0])

        # If Replica Exchange was used, compute the fraction at the target pH or Redox potential using fitpkaeo.py
        if (data.remd):
            # Testing if all fractions are equal
            AllEqual = True
            if ((max(fracs)-min(fracs))>0.03): AllEqual = False
            if (AllEqual):
                data.frac = 0.0
                for rep in range(0,opt.ng):
                    data.frac += fracs[rep]
                data.frac = data.frac/opt.ng
            else:
                if (data.modeopt == 1):
                    cmd  ="for file in reordered.cpout*; do %scphstats -i %s.temporary $file; done | %sfitpkaeo.py" % (opt.binpath.strip(),opt.cpin,opt.binpath.strip())
                else:
                    cmd  ="for file in reordered.ceout*; do %scestats -i %s.temporary $file; done | %sfitpkaeo.py" % (opt.binpath.strip(),opt.cein,opt.binpath.strip())
                pipe = subprocess.Popen([cmd], stdout=subprocess.PIPE, shell=True)
                lines = pipe.stdout.readlines()
                cnt = 0
                for line in lines:
                    if re.search(resname, line):
                        data.pkaeo = float(line.split()[4])
                        data.hillc = float(line.split()[12])
                        cnt += 1
                        break
                if (cnt == 0): ErrorExit(opt,log,'The execution of fitpkaeo.py failed. The command executed was:\n\n       %s' % (cmd))
                if (data.modeopt == 1):
                    data.frac = 1.0/(1.0+math.exp(data.hillc*math.log(10.0)*(opt.target-data.pkaeo)))
                else:
                    data.frac = 1.0/(1.0+math.exp(data.hillc*data.FARADAY*(opt.target-data.pkaeo)/(data.KB*data.temp0)))
        else:
            data.frac = fracs[0]
        
        # Deleting temporary files
        for rep in range(0,opt.ng):
            os.remove(data.mdins[rep]+".temporary")
            if (equi == 2):
                if (data.remd):
                    if (len(data.cprestrts)>0):
                        if (os.path.exists(data.cprestrts[rep]+".temporary")): os.remove(data.cprestrts[rep]+".temporary")
                    if (len(data.cerestrts)>0):
                        if (os.path.exists(data.cerestrts[rep]+".temporary")): os.remove(data.cerestrts[rep]+".temporary")
                else:
                    if (os.path.exists(opt.cprestrt+".temporary")): os.remove(opt.cprestrt+".temporary")
                    if (os.path.exists(opt.cerestrt+".temporary")): os.remove(opt.cerestrt+".temporary")
        os.remove(cin+".temporary")
        if (data.remd): os.remove(opt.groupfile+".temporary")
        
        if (data.remd):
            for rep in range(0,opt.ng):
                if (data.modeopt == 1):
                    PrintLog(opt,log,"   The fraction of protonated species for pH = %7.3f is %6.2f%% for the Residue '%s'" % (data.solvphe[rep],fracs[rep]*100, resname))
                else:
                    PrintLog(opt,log,"   The fraction of reduced species for E = %10.5f V is %6.2f%% for the Residue '%s'" % (data.solvphe[rep],fracs[rep]*100, resname))
            if (data.modeopt == 1):
                if (not AllEqual):
                    PrintLog(opt,log,"   Fitted values for Residue '%s': pKa = %7.3f and Hill coefficient = %.3f" % (resname,data.pkaeo,data.hillc))
                    PrintLog(opt,log,"   The computed fraction of protonated species at the target pH = %7.3f is %6.2f%% for the Residue '%s'" % (opt.target,data.frac*100, resname))
            else:
                if (not AllEqual):
                    PrintLog(opt,log,"   Fitted values for Residue '%s': Eo = %10.5f V and Hill coefficient = %.3f" % (resname,data.pkaeo,data.hillc))
                    PrintLog(opt,log,"   The computed fraction of reduced species at the target E = %10.5f V is %6.2f%% for the Residue '%s'" % (opt.target,data.frac*100, resname))
        else:
            if (data.modeopt == 1):
                PrintLog(opt,log,"   The fraction of protonated species is %6.2f%% for the Residue '%s'" % (data.frac*100, resname))
            else:
                PrintLog(opt,log,"   The fraction of reduced species is %6.2f%% for the Residue '%s'" % (data.frac*100, resname))
                
def EstimateDG(opt,log,data):
    """
    When this routine is activated, the program starts trying to finding the value of DELTAGREF
    using estimatives and keeps increasing the number of steps on the simulation
    """
    
    # Tell the user at which point of the execution we are
    if (data.step == 1):
        PrintLog(opt,log," We start with accurate estimatives of DELTAGREF values starting from the value given at the argument -dgrefest (%.6f kcal/mol), and we will keep increasing the number of steps." % data.dg)
    else:
        PrintLog(opt,log,"\n From this step on, we will make more accurate estimatives of DELTAGREF values and we will keep increasing the number of steps.")
    if (data.modeopt == 1):
        if (data.remd):
            txt = "target pH = "+str(opt.target)
        else:
            txt = "pH = "+str(data.solvphe[0])
    else:
        if (data.remd):
            txt = "target E = "+str(opt.target)+" V"
        else:
            txt = "E = "+str(data.solvphe[0])+" V"
    PrintLog(opt,log," The convergence criterium is: %.2f%% >= fraction >= %.2f%% for %s.\n" % ((0.5-opt.fracthreshold/2.0)*100.0,(0.5+opt.fracthreshold/2.0)*100.0,txt))
    
    if (data.remd):
        data.num = data.numexchg
    else:
        data.num = data.nstlim
    while (data.step <= opt.maxsteps):
        # Correcting dg with an estimative of the true value
        data.dg += data.KB*data.temp0*data.factor*math.log(1.0/data.frac-1.0)
        # Running AMBER, if no equi equilibration is required
        if (opt.noequi):
            RunAMBER(opt,log,data,0)
        else:
            RunAMBER(opt,log,data,1)
            RunAMBER(opt,log,data,2)
        # Testing if fraction was inside the desired range for 3 consecutive steps
        if(data.frac >= (0.5-opt.fracthreshold/2.0) and data.frac <= (0.5+opt.fracthreshold/2.0)):
            # We found a converged value of DELTAGREF
            data.FinalPoint = True
            break
        # Increasing the execution step
        data.step += 1

def main(opt):
    """
    This is the main function to execute the program
    """
    # Creating a simple class to be able to store some variables in a simple manner
    class GroupVariables(object):
        pass
    
    # Creating a simple object that will contain several different variables
    data = GroupVariables()
    
    # Constants to be used
    data.KB = 0.001987192639 # kcal/(mol.K)
    data.FARADAY = 23.06054801 # (kcal/mol)/V
    
    # Opening the log file, if set:
    if (opt.log):
        log = open(opt.log, "w")
    else:
        log = None
    
    # Checking if Replica Exchange will be performed or not
    if (opt.groupfile):
        data.remd = True
    else:
        data.remd = False
    
    # Check the required arguments
    if not os.path.exists(opt.mdexec):
        ErrorExit(opt,log,'mdexec file (%s) is missing' % opt.mdexec)
    if (data.remd):
        if not os.path.exists(opt.groupfile):
            ErrorExit(opt,log,'groupfile file (%s) is missing' % opt.groupfile)
        if (opt.ng < 2):
            ErrorExit(opt,log,'The value of ng needs to be higher than 1.')
        if (not opt.target):
            ErrorExit(opt,log,'The value of target pKa or Eo was not set. Please set it using the -target flag and run the program again.')
        if (opt.do_parallel == "do_parallel"):
            opt.do_parallel = "mpirun -np "+str(opt.ng)
    else:
        if not os.path.exists(opt.mdin):
            ErrorExit(opt,log,'mdin file (%s) is missing' % opt.mdin)
        if not os.path.exists(opt.prmtop):
            ErrorExit(opt,log,'prmtop file (%s) is missing' % opt.prmtop)
        if not os.path.exists(opt.inpcrd):
            ErrorExit(opt,log,'inpcrd file (%s) is missing' % opt.inpcrd)
        if (not os.path.exists(opt.cpin) and not os.path.exists(opt.cein)):
            ErrorExit(opt,log,'cpin file (%s) or cein file (%s) is missing' % (opt.cpin,opt.cein))
    if (opt.maxsteps <= 2):
        ErrorExit(opt,log,'The value of maxsteps needs to be higher than 2.')
    if (opt.dginterval == 0.0 and not opt.dgrefrange and not opt.dgrefest):
        ErrorExit(opt,log,'The value of dginterval cannot be 0.')
    if (opt.dginterval < 0.0):
        opt.dginterval = -opt.dginterval
    if (opt.fracthreshold > 0.3):
        ErrorExit(opt,log,'The value of fracthreshold cannot be higher than 0.3.')
    if (opt.binpath.strip()!="" and opt.binpath.strip()[-1]!="/"):
        opt.binpath+="/"

    # If dgrefrange was given, check if values are equal
    if (opt.dgrefrange):
        if(opt.dgrefrange[0] == opt.dgrefrange[1]):
            ErrorExit(opt,log,'The values you provided on the -dgrefrange argument are equal. Please provide different values.')
            
    data.mdins = []
    data.mdouts = []
    # If Replica Exchange is going to be performed, read the mdin, mdout, cpin and/or cein file names from the groupfile
    if (data.remd):
        # Defining some extra variables
        data.remlogs = []
        data.mdcrds = []
        data.mdinfos = []
        data.restrts = []
        data.cpouts = []
        data.cprestrts = []
        data.ceouts = []
        data.cerestrts = []
        cnt = 0
        group_file = open(opt.groupfile, "r")
        for line in group_file:
            if (line.lower().replace(" ","")[0] != "#" and len(line.replace(" ","")) > 1):
                for i in range(0,len(line.split())):
                    if (line.split()[i] == "-remlog"):
                        data.remlogs.append(line.split()[i+1])
                    elif (line.split()[i] == "-i"):
                        data.mdins.append(line.split()[i+1])
                    elif (line.split()[i] == "-o"):
                        data.mdouts.append(line.split()[i+1])
                    elif (line.split()[i] == "-x"):
                        data.mdcrds.append(line.split()[i+1])
                    elif (line.split()[i] == "-inf"):
                        data.mdinfos.append(line.split()[i+1])
                    elif (line.split()[i] == "-r"):
                        data.restrts.append(line.split()[i+1])
                    elif (line.split()[i] == "-cpin"):
                        opt.cpin = line.split()[i+1]
                    elif (line.split()[i] == "-cpout"):
                        data.cpouts.append(line.split()[i+1])
                    elif (line.split()[i] == "-cprestrt"):
                        data.cprestrts.append(line.split()[i+1])
                    elif (line.split()[i] == "-cein"):
                        opt.cein = line.split()[i+1]
                    elif (line.split()[i] == "-ceout"):
                        data.ceouts.append(line.split()[i+1])
                    elif (line.split()[i] == "-cerestrt"):
                        data.cerestrts.append(line.split()[i+1])
                cnt += 1
        group_file.close()
        # Checking if groupfile is valid
        if (cnt==0):
            ErrorExit(opt,log,'The groupfile %s is not valid.'%opt.groupfile)
        if ((len(data.cpouts)>0 and len(data.cpouts)<cnt) or (len(data.ceouts)>0 and len(data.ceouts)<cnt)):
            ErrorExit(opt,log,'The cpout or ceout files are not properly set inside the groupfile %s .'%opt.groupfile)
    # If Replica Exchange is not going to be done, we set some variables for simplicity
    else:
        data.mdins.append(opt.mdin)
        data.mdouts.append(opt.mdout)
        opt.ng = 1
            
    # If Replica Exchange will be done, check if the number of groups match with what is on the group file
    if (data.remd and opt.ng != cnt):
        ErrorExit(opt,log,'The number of groups on the groupfile do not match the number provided at the ng flag. Please check it and run the program again.')
        
    # If Replica Exchange will be done, check if the mdin files contain the keyword numexchg
    if (data.remd):
        for i in range(0,opt.ng):
            # Checking first if the mdin file exists
            if (not os.path.exists(data.mdins[i])): ErrorExit(opt,log,'The mdin file %s does not exist!'%data.mdins[i])
            cnt = 0
            mdin_file = open(data.mdins[i], "r")
            for line in mdin_file:
                if re.search("numexchg", line.lower().replace(" ","")):
                    cnt += 1
                    break
            mdin_file.close()
            if (cnt ==0):
                ErrorExit(opt,log,'The mdin file (%s) provided is not valid for Replica Exchange simulations because the keyword numexchg was not provided. Please check it and run the program again.'%(data.mdins[i]))
    
    # Check if the cpin or cein file has DELTAGREF written at the STATENE flag.
    # If everything is okay, we decide wheater we are going to find DELTAGREF for
    # a pH or a Redox titratable residue
    cpincnt = 0
    if (os.path.exists(opt.cpin)):
        cpin_file = open(opt.cpin, "r")
        for line in cpin_file:
            if re.search("deltagref", line.lower()):
                 cpincnt += 1
        cpin_file.close()
    ceincnt = 0
    if (os.path.exists(opt.cein)):
        cein_file = open(opt.cein, "r")
        for line in cein_file:
            if re.search("deltagref", line.lower()):
                 ceincnt += 1
        cein_file.close()
    
    PrintLog(opt,log," Checking cpin file and/or cein file.")
    if (cpincnt > 0 and ceincnt > 0):
        ErrorExit(opt,log,
		  'Both the cpin file (%s) and the cein file (%s) have '
                  'DELTAGREF written on the STATENE flag. finddgref.py '
                  'can only run for one of them at the same time. Let DELTAGREF '
                  'in only one of the files and run the program again.' % (opt.cpin,opt.cein))
    elif (cpincnt == 0 and ceincnt == 0):
        if (not os.path.exists(opt.cpin)):
            txt = ('The cein file (%s) does not have ' % opt.cein)
        elif (not os.path.exists(opt.cein)):
            txt = ('The cpin file (%s) does not have ' % opt.cpin)
        else:
            txt = ('Neither the cpin file (%s) nor the cein file (%s) have ' % (opt.cpin,opt.cein))
        txt += ('DELTAGREF written on the STATENE flag. Please substitute at '
                'least one of the values on STATENE by DELTAGREF (not all of them) and run the program again.')
        ErrorExit(opt,log,txt)
    elif (cpincnt > 1 or ceincnt > 1):
        ErrorExit(opt,log,
		  'DELTAGREF is written in more than one line on your cpin file (%s)'
                  'or the cein file (%s). Please write DELTAGREF on only one STATENE line'
                  'and run the program again.' % (opt.cpin,opt.cein))
    
    # Telling the user the type of residue
    txt = " We are going to find DELTAGREF for a "
    if (cpincnt > 0):
        data.modeopt = 1
        txt += "pH"
    else:
        data.modeopt = 2
        txt += "Redox"
    txt += " titratable residue"
    if (data.remd):
        PrintLog(opt,log,txt+" using Replica Exchange.")
    else:
        PrintLog(opt,log,txt+" without using Replica Exchange.")
    del cpincnt, ceincnt
    
    # If doing Replica Exchange, check if the cpout/cprestrt or ceout/cerestrt files are defined inside the groupfile
    if (data.remd):
        # cpout or ceout files
        if (data.modeopt == 1 and len(data.cpouts)==0):
            for i in range(opt.ng):
                data.cpouts.append("cpout.%03d"%i)
        elif (data.modeopt == 2 and len(data.ceouts)==0):
            for i in range(opt.ng):
                data.ceouts.append("ceout.%03d"%i)
        # cprestrt or cerestrt files
        if (data.modeopt == 1 and len(data.cprestrts)==0):
            for i in range(opt.ng):
                data.cprestrts.append("cprestrt.%03d"%i)
        elif (data.modeopt == 2 and len(data.cerestrts)==0):
            for i in range(opt.ng):
                data.cerestrts.append("cerestrt.%03d"%i)            
    
    # If Replica Exchange is to be done, check if all cpin or cein files provided on the groupfile are the same
    # Checking also if the -rem is present in all lines
    if (data.remd):
        cnt = 0
        group_file = open(opt.groupfile, "r")
        for line in group_file:
            if (line.lower().replace(" ","")[0] != "#" and len(line.replace(" ","")) > 1):
                for i in range(0,len(line.split())):
                    if (line.split()[i] == "-cpin"):
                        if (line.split()[i+1] != opt.cpin):
                            ErrorExit(opt,log,'The name of the cpin file is not the same for all replicas in the groupfile. Please check it and run the program again.')
                    elif (line.split()[i] == "-cein"):
                        if (line.split()[i+1] != opt.cein):
                            ErrorExit(opt,log,'The name of the cein file is not the same for all replicas in the groupfile. Please check it and run the program again.')
                    elif (line.split()[i] == "-rem"):
                        cnt += 1
        if (cnt == 0):
            ErrorExit(opt,log,'The groupfile (%s) provided is not valid for Replica Exchange simulations because the flag -rem was not provided. Please check it and run the program again.'%(opt.groupfile))
        elif (cnt != opt.ng):
            ErrorExit(opt,log,'The groupfile (%s) provided is not valid for Replica Exchange simulations because the flag -rem was not provided for all replicas. Please check it and run the program again.'%(opt.groupfile))
        group_file.close()
    
    # Get the value(s) of the solvent pH or solvent redox potential from the mdin file(s)
    if (data.modeopt == 1):
        txt = "solvph="
    else:
        txt = "solve="
    data.solvphe = []
    for i in range(0,opt.ng):
        cnt = 0
        mdin_file = open(data.mdins[i], "r")
        for line in mdin_file:
            if re.search(txt, line.lower().replace(" ","")):
                cnt += 1
                break
        mdin_file.close()
        if (data.remd):
            txt2 = " of replica %4d"%(i+1)
        else:
            txt2 = ""
        # If the solvph or solve flag is not present the the mdin file, we print an error
        if (cnt == 0):
            if (data.modeopt == 1):
                ErrorExit(opt,log,"The solvent pH value"+txt2+" could not be read from the mdin file %s. The flag solvph is not present in this file." % (data.solvphe[i],data.mdins[i]))
            else:
                ErrorExit(opt,log,"The solvent redox potential value"+txt2+" could not be read from the mdin file %s. The flag solve is not present in this file." % (data.solvphe[i],data.mdins[i]))
        # If it is present, we read it.
        else:
            lines = line.lower().replace(" ","").split(txt)
            lines = lines[1].split(",")
            data.solvphe.append(float(lines[0]))
            if (data.modeopt == 1):
                PrintLog(opt,log," The solvent pH value"+txt2+" is %7.3f and was loaded from the mdin file (%s)." % (data.solvphe[i],data.mdins[i]))
            else:
                PrintLog(opt,log," The solvent redox potential value"+txt2+" is %10.5f V and was loaded from the mdin file (%s)." % (data.solvphe[i],data.mdins[i]))
        del txt2
        
    # If Replica Exchange is to be done, print the opt.target value
    if (data.remd):
        if (data.modeopt == 1):
            PrintLog(opt,log," The target pKa is %10.5f." % (opt.target))
        else:
            PrintLog(opt,log," The target Standard Redox Potential (Eo) is %10.5f V." % (opt.target))

    # Checking if all values on STATENE are equal to DELTAGREF
    if (data.modeopt == 1):
        cin_file = open(opt.cpin, "r")
    else:
        cin_file = open(opt.cein, "r")
    for line in cin_file:
        if re.search("deltagref", line.lower()):
            break
    cin_file.close()
    lines = line.lower().replace(" ","").rstrip().split("statene=")[1][:-1].split(",")
    AllEqual = True
    for i in range(1,len(lines)):
        if (not lines[i] == lines[i-1]):
            AllEqual = False
    if (AllEqual):
        ErrorExit(opt,log,
                  'All values of STATENE are equal to DELTAGREF on your cpin or cein file. '
                  'At least one of the values must be numeric. Please fix that and run the '
                  'program again.')

    # We don't know at which position of the STATENE flag the string DELTAGREF was written.
    # We need to check the values of PROTCNT or ELECCNT in order to decide how we are going to increase
    # or decrease the value of DELTAGREF on the Delta G ref estimative phase (last phase of execution).
    if (data.modeopt == 1):
        cin_file = open(opt.cpin, "r")
        txt = "protcnt="
    else:
        cin_file = open(opt.cein, "r")
        txt = "eleccnt="
    cnt = 0
    for line in cin_file:
        if re.search("statene=", line.lower()):
            line1 = line
            cnt += 1
        if re.search(txt, line.lower()):
            line2 = line
            cnt += 1
        if (cnt == 2):
            break
    cin_file.close()
    lines1 = line1.lower().replace(" ","").rstrip().split("statene=")[1][:-1].split(",")
    lines2 = line2.lower().replace(" ","").rstrip().split(txt)[1][:-1].split(",")
    AllEqual = True
    cnt = 0
    for i in range(0,len(lines)):
        if (lines1[i]=="deltagref"):
            if (cnt == 0):
                line = lines2[i]
                cnt = 1
            else:
                if (lines2[i] != line):
                    AllEqual = False
    if (not AllEqual):
        ErrorExit(opt,log,
                  'The positions that have DELTAGREF on the STATENE flag must have the same '
                  'value on the corresponding positions at the PROTCNT or ELECCNT flag on your cpin or cein file. '
                  'Please fix that and run the program again.')
    AmImin = True
    for i in range(0,len(lines)):
        if (int(lines2[i]) > int(line)):
            AmImin = False
            break
    if (AmImin):
        data.factor = 1
    else:
        data.factor = -1
    del line1, lines1, line2, lines2, AllEqual, AmImin

    # Get the value of temp0 from the mdin file
    cnt = 0
    mdin_file = open(data.mdins[0], "r")
    for line in mdin_file:
        if re.search("temp0=", line.lower().replace(" ","")):
            cnt += 1
            break
    mdin_file.close()
    # If the temp0 flag is not present the the mdin file, this means the default
    # temperature (300 K) was used.
    if (cnt == 0):
        data.temp0 = 300.0
        PrintLog(opt,log," The temperature is %8.2f K. Using the default value because temp0 was not found at the mdin file (%s)." % (data.temp0,data.mdins[0]))
    # If it is present, we read it.
    else:
        lines = line.lower().replace(" ","").split("temp0=")
        lines = lines[1].split(",")
        data.temp0 = float(lines[0])
        PrintLog(opt,log," The temperature is %8.2f K and was loaded from the mdin file (%s)." % (data.temp0,data.mdins[0]))
    
    # Get the value of ntcnstph or ntcnste from the mdin file
    mdin_file = open(data.mdins[0], "r")
    if (data.modeopt == 1):
        txt="ntcnstph"
    else:
        txt="ntcnste"
    cnt = 0
    for line in mdin_file:
        if re.search(txt, line):
            cnt += 1
            break
    if (cnt == 0):
        data.ntcnst = 500
        PrintLog(opt,log," The value of %s is not present at the mdin file (%s). Please set it if you want finddgref.py to run more efficiently." % (txt,data.mdins[0]))
    else:
        lines = line.replace(" ","").split(txt+"=")
        lines = lines[1].split(",")
        data.ntcnst = int(lines[0])
        mdin_file.close()
        PrintLog(opt,log," According to the mdin file (%s), the value of %s is %d." % (data.mdins[0],txt,data.ntcnst))
    
    # Get the value of nstlim from the mdin file
    mdin_file = open(data.mdins[0], "r")
    cnt = 0
    for line in mdin_file:
        if re.search("nstlim", line):
            cnt += 1
            break
    if (cnt == 0):
        data.nstlim = 1
        PrintLog(opt,log," The value of nstlim is not set at the mdin file (%s). Please set it if you want finddgref.py to run more efficiently." % (data.mdins[0]))
    else:
        lines = line.replace(" ","").split("nstlim=")
        lines = lines[1].split(",")
        data.nstlim = int(lines[0])
        mdin_file.close()
        PrintLog(opt,log," According to the mdin file (%s), the value of nstlim is %d." % (data.mdins[0],data.nstlim))
    # Get the value of numexchg from the mdin file, if doing Replica Exchange
    if (data.remd):
        mdin_file = open(data.mdins[0], "r")
        cnt = 0
        for line in mdin_file:
            if re.search("numexchg", line):
                cnt += 1
                break
        if (cnt == 0):
            ErrorExit(opt,log," The value of numexchg is not set at the mdin file (%s)." % (data.mdins[0]))
        else:
            lines = line.replace(" ","").split("numexchg=")
            lines = lines[1].split(",")
            data.numexchg = int(lines[0])
            mdin_file.close()
            PrintLog(opt,log," According to the mdin file (%s), the value of numexchg is %d." % (data.mdins[0],data.numexchg))
    PrintLog(opt,log,"")
        
    # Check if cphstats or cestats are set
    if (data.modeopt == 1):
        try:
            job = subprocess.Popen(["%scphstats"%opt.binpath.strip()], stdout=subprocess.PIPE)
        except:
            amberhome = os.getenv('AMBERHOME') or '$AMBERHOME'
            if (opt.binpath.strip()==""):
                ErrorExit(opt,log,
		          'Could not execute cphstats. Please make sure '
                          'you have sourced %s/AMBER.sh (if you are using sh/ksh/'
                          'bash/zsh) or %s/AMBER.csh (if you are using csh/tcsh). '
                          'You may also provide the path to cphstats rerunning the '
                          'program again setting -bin-path' %
                          (amberhome, amberhome))
            else:
                ErrorExit(opt,log,
		          'Could not execute cphstats at %s. Please set -bin-path '
                          'to a folder where cphstats can be found. Otherwise, '
                          'do not state -bin-path and source %s/AMBER.sh (if you '
                          'are using sh/ksh/bash/zsh) or %s/AMBER.csh (if you are '
                          'using csh/tcsh)' % (opt.binpath.strip(), amberhome, amberhome))
    else:
        try:
            job = subprocess.Popen(["%scestats"%opt.binpath.strip()], stdout=subprocess.PIPE)
        except:
            amberhome = os.getenv('AMBERHOME') or '$AMBERHOME'
            if (opt.binpath.strip()==""):
                ErrorExit(opt,log,
		          'Could not execute cestats. Please make sure '
                          'you have sourced %s/AMBER.sh (if you are using sh/ksh/'
                          'bash/zsh) or %s/AMBER.csh (if you are using csh/tcsh). '
                          'You may also provide the path to cestats rerunning the '
                          'program again setting -bin-path' %
                          (amberhome, amberhome))
            else:
                ErrorExit(opt,log,
		          'Could not execute cestats at %s. Please set -bin-path '
                          'to a folder where cestats can be found. Otherwise, '
                          'do not state -bin-path and source %s/AMBER.sh (if you '
                          'are using sh/ksh/bash/zsh) or %s/AMBER.csh (if you are '
                          'using csh/tcsh)' % (opt.binpath.strip(), amberhome, amberhome))
    job.send_signal(signal.SIGINT)
    
    # If Replica Exchange will be done, check if fitpkaeo.py is set
    if (data.remd):
        try:
            job = subprocess.Popen(["%sfitpkaeo.py"%opt.binpath.strip()], stdout=subprocess.PIPE)
        except:
            amberhome = os.getenv('AMBERHOME') or '$AMBERHOME'
            if (opt.binpath.strip()==""):
                ErrorExit(opt,log,
		          'Could not execute fitpkaeo.py. Please make sure '
                          'you have sourced %s/AMBER.sh (if you are using sh/ksh/'
                          'bash/zsh) or %s/AMBER.csh (if you are using csh/tcsh). '
                          'You may also provide the path to fitpkaeo.py rerunning the '
                          'program again setting -bin-path' %
                          (amberhome, amberhome))
            else:
                ErrorExit(opt,log,
		          'Could not execute fitpkaeo.py at %s. Please set -bin-path '
                          'to a folder where fitpkaeo.py can be found. Otherwise, '
                          'do not state -bin-path and source %s/AMBER.sh (if you '
                          'are using sh/ksh/bash/zsh) or %s/AMBER.csh (if you are '
                          'using csh/tcsh)' % (opt.binpath.strip(), amberhome, amberhome))
    job.send_signal(signal.SIGINT)

    # We now execute AMBER to load the fraction values for the dgrefrange values given as input
    
    # Did we have success on finding one value of DELTAGREF where 0.2 > data.frac > 0.8? Not yet.
    if(not opt.dgrefest):
        FoundPoint = False
    # If -dgrefest is given, we go straight to the accurate estimatives phase of the program. That is, we force FoundPoint to be true
    else:
        data.dg = opt.dgrefest
        data.frac = 0.5
        FoundPoint = True
    
    # Initiating the number of steps
    data.step = 1
    
    if (data.modeopt == 1):
        txt1 = "protonated"
    else:
        txt1 = "reduced"
    # If the argument -dgrefest and -dgrefrange are not given, we try to find the range of values for DELTAGREF automatically
    if (not opt.dgrefrange and not FoundPoint):
        PrintLog(opt,log," The program will try to find a range of values for DELTAGREF automatically, as the argument -dgrefrange was not given\n")
        while (data.step <= opt.maxsteps):
            if ((data.step-1)%2 == 0):
                data.dg = opt.dginterval*float(data.step-1)/2.0
            else:
                data.dg = -opt.dginterval*float(data.step)/2.0
            data.num = 10*data.ntcnst
            RunAMBER(opt,log,data,0)
            if(data.frac > 0.2 and data.frac < 0.8):
                # We found a dg value in which an initial estimative of the desirable DELTAGREF value can be done
                FoundPoint = True
                data.step += 1
                break
            else:
                if (data.step == 1):
                    frac_old = data.frac
                else:
                    if (abs(data.frac-frac_old)>0.2):
                        if (data.step == 2):
                            if(data.frac >= 0.8):
                                dgmin = 0.0
                                dgmax = data.dg
                            else:
                                dgmin = data.dg
                                dgmax = 0.0
                        else:
                            if(data.frac >= 0.8):
                                if (data.dg>0.0):
                                    dgmin = data.dg-opt.dginterval
                                else:
                                    dgmin = data.dg+opt.dginterval
                                dgmax = data.dg
                            else:
                                dgmin = data.dg
                                if (data.dg>0.0):
                                    dgmax = data.dg-opt.dginterval
                                else:
                                    dgmax = data.dg+opt.dginterval
                        data.step += 1
                        break
            data.step += 1
    # If -dgrefest is not given, compute the fractions for the values given for dgrefrange if they were given
    elif(not FoundPoint):
        PrintLog(opt,log," The range of values for DELTAGREF is %.6f to %.6f kcal/mol (given on the argument -dgrefrange)\n" % (opt.dgrefrange[0],opt.dgrefrange[1]))
        data.dg = opt.dgrefrange[0]
        data.num = 10*data.ntcnst
        RunAMBER(opt,log,data,0)
        frac1 = data.frac
        data.step += 1
        if(frac1 > 0.2 and frac1 < 0.8):
            # We found a dg value in which an initial estimative of the desirable DELTAGREF value can be done
            FoundPoint = True
        else:
            data.dg = opt.dgrefrange[1]
            data.num = 10*data.ntcnst
            RunAMBER(opt,log,data,0)
            frac2 = data.frac
            data.step += 1
            if(frac2 > 0.2 and frac2 < 0.8):
                # We found a dg value in which an initial estimative of the desirable DELTAGREF value can be done
                FoundPoint = True
            else:
                # Checking if the range given on dgrefrange is acceptable
                if(abs(frac1-frac2)<=0.2):
                    txt  = 'The values entered for dgrefrange (%.6f and %.6f kcal/mol) are not valid because the ' % (opt.dgrefrange[0],opt.dgrefrange[1])
                    txt += 'fraction of '+txt1+' species is the same (or almost the same) for both values. '
                    txt += 'Please try to run the program again with a valid range of values for -dgrefrange, or '
                    txt += 'let the program find this range automatically by not providing the argument -dgrefrange.'
                    ErrorExit(opt,log,txt)
                else:
                    if(frac1 >= 0.8):
                        dgmax = opt.dgrefrange[0]
                        dgmin = opt.dgrefrange[1]
                    else:
                        dgmax = opt.dgrefrange[1]
                        dgmin = opt.dgrefrange[0]
    
    # Did we find the final value of DELTAGREF? Not yet
    data.FinalPoint = False
    
    # If we were lucky to find one value of DELTAGREF where 0.2 < data.frac < 0.8 already, start doing estimatives
    if (FoundPoint and not data.step == (opt.maxsteps+1)):
        EstimateDG(opt,log,data)
    # Else, we try to find one value of DELTAGREF where 0.2 > data.frac > 0.8
    elif(not data.step == (opt.maxsteps+1)):
        if (data.modeopt == 1):
            if (data.remd):
                str1="target solvent pH = "+str(opt.target)
            else:
                str1="solvent pH = "+str(data.solvphe[0])
        else:
            if (data.remd):
                str1="target solvent redox potential = "+str(opt.target)+" V"
            else:
                str1="solvent redox potential = "+str(data.solvphe[0])+" V"
        PrintLog(opt,log,"\n The value of DELTAGREF in which 20% > fraction of "+txt1+" species > 80% for "+str1+" is between %.6f and %.6f kcal/mol.\n" % (dgmin,dgmax))
        # Try to find a value of DELTAGREF in which 20% > fraction > 80%
        data.num = 40*data.ntcnst
        while (data.step <= opt.maxsteps):
            data.dg = (dgmin+dgmax)/2.0
            RunAMBER(opt,log,data,0)
            if(data.frac > 0.2 and data.frac < 0.8):
                # We found a dg value in which an initial estimative of the desirable DELTAGREF value can be done
                FoundPoint = True
                data.step += 1
                break
            else:
                if(data.frac >= 0.2):
                    dgmax = data.dg
                else:
                    dgmin = data.dg
            data.step += 1

    # At this point the program either exited due the maximum number of steps reached, or FoundPoint is true
    if (FoundPoint and not data.FinalPoint):
        EstimateDG(opt,log,data)

    # Removing all AMBER output files, if the user so desire.
    if (opt.rmouts):
        txt = '\n Erasing all output files '
        if (data.remd): txt += 'present on the groupfile (%s) '%(opt.groupfile)
        txt += 'generated by AMBER'
        PrintLog(opt,log,txt)
        if (data.remd):
            for rep in range(0,opt.ng):
                if (len(data.remlogs) > rep):
                    if (os.path.exists(data.remlogs[rep])): os.remove(data.remlogs[rep])
                if (len(data.mdcrds) > rep):
                    if (os.path.exists(data.mdcrds[rep])): os.remove(data.mdcrds[rep])
                if (len(data.mdinfos) > rep):
                    if (os.path.exists(data.mdinfos[rep])): os.remove(data.mdinfos[rep])
                if (len(data.mdouts) > rep):
                    if (os.path.exists(data.mdouts[rep])): os.remove(data.mdouts[rep])
                if (len(data.restrts) > rep):
                    if (os.path.exists(data.restrts[rep])): os.remove(data.restrts[rep])
                if (len(data.cpouts) > rep):
                    if (os.path.exists(data.cpouts[rep])): os.remove(data.cpouts[rep])
                if (len(data.ceouts) > rep):
                    if (os.path.exists(data.ceouts[rep])): os.remove(data.ceouts[rep])
                if (len(data.cprestrts) > rep):
                    if (os.path.exists(data.cprestrts[rep])): os.remove(data.cprestrts[rep])
                if (len(data.cerestrts) > rep):
                    if (os.path.exists(data.cerestrts[rep])): os.remove(data.cerestrts[rep])
                if (data.modeopt == 1):
                    if (os.path.exists("reordered.cpout.pH_%.2f"%(data.solvphe[rep]))): os.remove("reordered.cpout.pH_%.2f"%(data.solvphe[rep]))
                else:
                    if (os.path.exists("reordered.ceout.E_%.5f"%(data.solvphe[rep]))): os.remove("reordered.ceout.E_%.5f"%(data.solvphe[rep]))
        else:
            os.remove(opt.mdcrd)
            os.remove(opt.mdinfo)
            os.remove(opt.mdout)
            os.remove(opt.restrt)
            if (os.path.exists(opt.cpin)):
                if (os.path.exists(opt.cpout)): os.remove(opt.cpout)
                if (os.path.exists(opt.cprestrt)):os.remove(opt.cprestrt)
            if (os.path.exists(opt.cein)):
                if (os.path.exists(opt.ceout)): os.remove(opt.ceout)
                if (os.path.exists(opt.cerestrt)):os.remove(opt.cerestrt)
        
    # Now we write the final line on the program saying whether we had sucess or not
    if (data.modeopt == 1):
        if (data.remd):
            str1="target solvent pH = %7.3f"%(opt.target)
        else:
            str1="solvent pH = %7.3f"%(data.solvphe[0])
    else:
        if (data.remd):
            str1="target solvent redox potential = %10.5f V"%(opt.target)
        else:
            str1="solvent redox potential = %10.5f V"%(data.solvphe[0])
    if (data.step == (opt.maxsteps+1) and not FoundPoint):
        txt  = 'After %d AMBER executions we could not find a value of DELTAGREF in which 20% > fraction of ' % opt.maxsteps
        txt += txt1+' species > 80% for '+str1+'. Try to run the program again with a larger value on the -maxsteps argument (> %d), ' % opt.maxsteps
        txt += 'or a larger value on the -dginterval argument (> %.2f kcal/mol).' % opt.dginterval
        ErrorExit(opt,log,txt)
    elif (data.step == (opt.maxsteps+1) and not data.FinalPoint):
        txt  = 'After %d AMBER executions we could not find a value of DELTAGREF that gives a converged fraction of ' % opt.maxsteps
        txt += txt1+' species close to 50% for '+str1+'. Try to run the program again with a larger value on the -maxsteps argument (> %d).' % opt.maxsteps
        ErrorExit(opt,log,txt)
    else:
        if (data.remd):
            txt = "for %d MD steps"%(data.num*data.nstlim)
        else:
            txt = "for %d MD steps"%(data.num)
        PrintLog(opt,log,"\n The value of DELTAGREF that gives a converged fraction of "+txt1+" species "+txt+" and for "+str1+" equals to %.2f%% is: DELTAGREF = %.6f kcal/mol" % (data.frac*100,data.dg))
        PrintLog(opt,log," The execution of finddgref.py ended with success.")
    
    # Close the log file, if it exists
    if (opt.log):
        log.close()
      
if __name__ == '__main__':
    opt = parser.parse_args()

    # Go ahead and execute the program.
    main(opt)
    sys.exit(0)
