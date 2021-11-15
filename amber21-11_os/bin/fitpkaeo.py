#!/Users/runner/miniforge3/conda-bld/ambertools_1635195607525/_h_env_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehol/bin/python


import os, sys
import numpy as np
import scipy.optimize as optimization

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

parser = ArgumentParser(epilog='''This program will perform the fitting of the pKa (based on the Henderson-Hasselbalch equation) or of
                        the standard Redox potential (Nernst equation) values and of the Hill coefficient. USAGE: pass the output of
                        cphstats or cestats for all pH or Redox potential values through pipe.''', usage='%(prog)s [Options]')
parser.add_argument('-v', '--version', action='version', version='%s: %s' %
                    (parser.prog, __version__), help='''show the program's version and exit''')
parser.add_argument('--author', action='version', version='%s author: %s (E-mail: %s )' %
                    (parser.prog, __author__,__email__), help='''show the program's author name and exit''')
group = parser.add_argument_group('Not-required Arguments')
group.add_argument('--usage', dest='usage', action='store_const',
                   help='''If stated, prints detailed information on how to execute the program.''',
                   const=True, default=False)
group.add_argument('--verbose', dest='verbose', action='store_const',
                   help='''If stated, prints verbose showing initial and final values of X2. Default: False''',
                   const=True, default=False)
group.add_argument('--graph', dest='graph', action='store_const',
                   help='''If stated, generates a plot of the fitted curve. Default: False''',
                   const=True, default=False)
group.add_argument('--graph-path', dest='graphpath', metavar='FILE', required=False,
                   help='Path to save the graph. Default: ./graph.png',
                   type=str, default='./graph.png')

def ErrorExit(text):
    """
    This function prints an error and kills the execution of the program
    """

    print(('\nERROR: %s' % text))
    print ('       The execution of fitpkaeo.py stopped')
    sys.exit(0)

def f(x,xc,k):
    """
    Compute the fitting function f
    """

    return 1.0/(1.0+np.exp(k*(x-xc)))

def X2(xc,k,xvals,yvals):
    """
    Compute the X2 factor
    """

    val = 0.0

    for n in range(len(xvals)):
        val+=(yvals[n]-f(xvals[n],xc,k))**2

    return val

def usage():
    print("The input for fitpkaeo.py needs to be passed through pipe.")
    print(" 1) In order to fit the pKa (using the Henderson-Hasselbalch equation) execute the following:")
    print("      for cpout_file in [list all cpout files]; do cphstats -i [cpin file] $cpout_file; done | fitpkaeo.py [additonal non-required arguments]")
    print("")
    print(" 2) In order to fit the standard Redox potential (using the Nernst equation) execute the following:")
    print("      for ceout_file in [list all ceout files]; do cestats -i [cein file] $ceout_file; done | fitpkaeo.py [additonal non-required arguments]")
    print("")
    sys.exit(0)

def main(opt):
    """
    This is the main function to execute the program
    """

    # Check if pipe input was given
    if (not sys.stdin.isatty()):
        input_stream = sys.stdin.readlines()
    elif (not opt.usage):
        print("ERROR: No input detected.\n       The execution of fitpkaeo.py stopped\n")
        usage()

    # If user requires more information on how to run the program
    if (opt.usage): usage()

    # Checking if we are going to perform a pKa or a Eo fitting
    words = input_stream[0].split()
    if (words[0]=="Solvent" and words[1]=="pH"):
        mode = 0
    elif (words[0]=="Redox" and words[1]=="potential"):
        mode = 1
    else:
        ErrorExit('The input provided is invalid. For more details about usage execute: fitpkaeo.py --usage')

    # Reading all pH or all Redox potential values
    cntph = 0
    cnteo = 0
    xvals = []
    temp0 = -1.0
    for line in input_stream:
        if(len(line.strip())>1):
            words = line.strip().split()
            if (words[0]=="Solvent" and words[1]=="pH"):
                xvals.append(float(words[3]))
                for j in range(cntph):
                    if (xvals[j]==xvals[cntph]): ErrorExit('Duplicate values of pH (%s) found!'%(words[3]))
                cntph+=1
            elif (words[0]=="Redox" and words[1]=="potential"):
                xvals.append(float(words[3]))
                val = float(words[-2])
                if (temp0 == -1.0):
                    temp0 = val
                elif (val != temp0):
                    ErrorExit('Different temperature values (%.2f and %.2f) found!'%(temp0,val))
                for j in range(cnteo):
                    if (xvals[j]==xvals[cnteo]): ErrorExit('Duplicate values of Redox potential (%s) found!'%(words[3]))
                cnteo+=1
    if (cntph > 0 and cnteo > 0): ErrorExit('The input provided is invalid, fitpkaeo.py cannot work for both pKa and Eo fitting at the same time. For more details about usage execute: fitpkaeo.py --usage')
    if (mode == 0):
        nxlines = cntph
    else:
        nxlines = cnteo
    del cntph, cnteo

    # Counting the number of residues
    cnt=0
    while True:
        words = input_stream[cnt+1].strip().split()
        if (len(words)>0):
            if (not words[3]=="Offset"):
                break
            else:
                cnt+=1
        else:
            break
    nres = cnt

    # Reading all fractions and names for each residue
    allyvals = np.zeros((nxlines,nres))
    resnames = []
    for i in range(nxlines):
        for res in range(nres):
            if (mode==0):
                allyvals[i,res] = float(input_stream[i*(nres+3)+1+res].strip().split()[9])
            else:
                allyvals[i,res] = float(input_stream[i*(nres+3)+1+res].strip().split()[12])
            if (i==0):
                resnames.append(input_stream[1+res].split(":")[0].strip())

    # Setting the initial parameter values
    xc0 = np.mean(xvals)
    KB = 0.001987192639 # kcal/(mol.K)
    FARADAY = 23.06054801 # (kcal/mol)/V
    k0_redox = FARADAY/(KB*temp0)
    k0_ph    = 2.302585093 # = ln(10)
    if (mode == 1):
        k0 = k0_redox
    else:
        k0 = k0_ph

    # Printing initial guess, if verbose is on
    if(opt.verbose):
        for res in range(nres):
            if (mode == 1):
                print(("Initial guess for residue %8s: Eo = %12.6f V; Hill coefficient = %6.3f; X2 factor = %10.6f" % (resnames[res],xc0,k0/k0_redox,X2(xc0,k0,xvals,allyvals[:,res]))))
            else:
                print(("Initial guess for residue %8s: pKa = %9.5f ; Hill coefficient = %6.3f; X2 factor = %10.6f" % (resnames[res],xc0,k0/k0_ph,X2(xc0,k0,xvals,allyvals[:,res]))))

    # Obtaining fitted values for each residue
    x0 = np.array([xc0,k0])
    xc = np.zeros(nres)
    k = np.zeros(nres)
    xc_err = np.zeros(nres)
    k_err = np.zeros(nres)
    for res in range(nres):
        # Checking if all y values are equal
        if (np.array_equal(allyvals[:,res],np.full((nxlines),allyvals[0,res]))):
            if(opt.verbose):
                if (mode == 1):
                    print(("Final results for residue %8s: Eo =     ALL FRACTIONS EQUAL      ; Hill coefficient =ALL FRACTIONS EQUAL; X2 factor = ALL FRAC. EQU.; Temperature = %.2f K" % (resnames[res],temp0)))
                else:
                    print(("Final results for residue %8s: pKa =   ALL FRACTIONS EQUAL   ; Hill coefficient =ALL FRACTIONS EQUAL; X2 factor = ALL FRACTIONS EQUAL" % (resnames[res])))
            else:
                if (mode == 1):
                    print(("%8s: Eo =     ALL FRACTIONS EQUAL      ; Hill coefficient =ALL FRACTIONS EQUAL; Temperature = %.2f K" % (resnames[res],temp0)))
                else:
                    print(("%8s: pKa =    ALL FRACTIONS EQUAL  ; Hill coefficient = ALL FRACTIONS EQUAL" % (resnames[res])))
        # If all y values are not equal, proceed
        else:
            try:
                vec, mat = optimization.curve_fit(f, xvals, allyvals[:,res], x0)
            except RuntimeError:
                if(opt.verbose):
                    if (mode == 1):
                        print(("Final results for residue %8s: Eo =        FITTING FAILED        ; Hill coefficient =   FITTING FAILED  ; X2 factor = FITTING FAILED; Temperature = %.2f K" % (resnames[res],temp0)))
                    else:
                        print(("Final results for residue %8s: pKa =      FITTING FAILED     ; Hill coefficient =   FITTING FAILED  ; X2 factor = FITTING FAILED" % (resnames[res])))
                else:
                    if (mode == 1):
                        print(("%8s: Eo =        FITTING FAILED        ; Hill coefficient =   FITTING FAILED  ; Temperature = %.2f K" % (resnames[res],temp0)))
                    else:
                        print(("%8s: pKa =      FITTING FAILED     ; Hill coefficient =   FITTING FAILED" % (resnames[res])))
            else:
                xc[res] = vec[0]
                k[res] = vec[1]
                # Computing standard error
                err = np.sqrt(np.diag(mat))
                xc_err[res] = err[0]
                k_err[res] = err[1]
                # Printing final results
                if(opt.verbose):
                    if (mode == 1):
                        print(("Final results for residue %8s: Eo = %12.6f (+- %.6f ) V; Hill coefficient = %6.3f (+- %.3f ); X2 factor = %14.6f; Temperature = %.2f K" % (resnames[res],xc[res],xc_err[res],k[res]/k0_redox,k_err[res]/k0_redox,X2(xc[res],k[res],xvals,allyvals[:,res]),temp0)))
                    else:
                        print(("Final results for residue %8s: pKa = %9.5f (+- %.5f ) ; Hill coefficient = %6.3f (+- %.3f ); X2 factor = %10.6f" % (resnames[res],xc[res],xc_err[res],k[res]/k0_ph,k_err[res]/k0_ph,X2(xc[res],k[res],xvals,allyvals[:,res]))))
                else:
                    if (mode == 1):
                        print(("%8s: Eo = %12.6f (+- %.6f ) V; Hill coefficient = %6.3f (+- %.3f ); Temperature = %.2f K" % (resnames[res],xc[res],xc_err[res],k[res]/k0_redox,k_err[res]/k0_redox,temp0)))
                    else:
                        print(("%8s: pKa = %9.5f (+- %.5f ) ; Hill coefficient = %6.3f (+- %.3f )" % (resnames[res],xc[res],xc_err[res],k[res]/k0_ph,k_err[res]/k0_ph)))

    # If requested, print a graph with the final results
    if (opt.graph):
        import matplotlib.pylab as plt
        import matplotlib.cm as cm
        colors = iter(cm.rainbow(np.linspace(0, 1, nres)))
        X = np.linspace(min(xvals), max(xvals), 500)
        fig = plt.figure()
        yvals = np.zeros(nxlines)
        for res in range(nres):
            for i in range(nxlines):
                yvals[i] = allyvals[i,res]
            c = next(colors)
            plt.plot(X, f(X,xc[res],k[res]), color=c, label="Fitted curve:%s"%resnames[res])
            plt.scatter(xvals, yvals, color=c, label="Simulation points:%s"%resnames[res])
        plt.ylim([-0.05,1.05])
        if (mode == 1):
            plt.xlabel('Redox Potential (V)')
            plt.ylabel('Fraction of Reduced Species')
        else:
            plt.xlabel('Solvent pH')
            plt.ylabel('Fraction of Protonated Species')
        plt.legend(loc='best')
        fig.savefig(opt.graphpath)
        print("Graph generated to ",opt.graphpath)

if __name__ == '__main__':
    opt = parser.parse_args()

    # Go ahead and execute the program.
    main(opt)
    sys.exit(0)
