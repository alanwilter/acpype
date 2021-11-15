#!/Users/runner/miniforge3/conda-bld/ambertools_1635195607525/_h_env_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehol/bin/python


import os, sys

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

parser = ArgumentParser(epilog='''This program will reorder Replica Exchange CPOUT and/or CEOUT files.
                                  It can be used even when pH or Redox Potential REMD are not used, for example:
                                  to reconstruct CPOUT files per temperature on a T-REMD simulation with constant
                                  pH on. This tool can also be used with Multidimentional REMD CPOUT and/or CEOUT
                                  files.''', usage='%(prog)s [Options]')
parser.add_argument('-v', '--version', action='version', version='%s: %s' %
                    (parser.prog, __version__), help='''show the program's version and exit''')
parser.add_argument('--author', action='version', version='%s author: %s (E-mail: %s )' %
                    (parser.prog, __author__,__email__), help='''show the program's author name and exit''')
parser.add_argument('-O', '--overwrite', dest='overwrite', action='store_const',
                   help='''Allow existing outputs to be overwritten. Default: False''',
                   const=True, default=False)
group = parser.add_argument_group('Required Arguments')
group.add_argument('-couts', dest='couts', metavar='FILE', nargs='*',
                   help='AMBER CPOUT and/or CEOUT files',
                   default=None)
group = parser.add_argument_group('Non-required Arguments')
group.add_argument('-prefix', dest='prefix', metavar='STRING',
                   help='Prefix of the reordered file names. Default: reordered',
                   type=str, default='reordered')

def ErrorExit(text):
    """
    This function prints an error and kills the execution of the program
    """
    
    print(('\nERROR: %s' % text))
    print ('       The execution of fixremdcouts.py stopped')
    sys.exit(0)
    
def main(opt):
    """
    This is the main function to execute the program
    """
    
    # Number of CPOUT and/or CEOUT files
    if (not opt.couts):
        ErrorExit('Please provide the CPOUT and/or CEOUT files (using the flag -couts)')
    else:
        numcouts = len(opt.couts)
        for i in range(numcouts):
            if (not os.path.exists(opt.couts[i])): ErrorExit('File %d does not exist!'%(opt.couts[i]))
    
    # Opening all CPOUT and/or CEOUT files and checking if their are valid
    couts_type = [0 for i in range(numcouts)]
    coutfiles = [open(opt.couts[i], "r") for i in range(numcouts)]
    numcpouts = 0
    numceouts = 0
    for i in range(numcouts):
        words = coutfiles[i].readline().lower().split()
        if (words[0] == "solvent" and words[1] == "ph:"):
            couts_type[i] = 4
            numcpouts += 1
        elif (words[0] == "redox" and words[1] == "potential:"):
            couts_type[i] = 5
            numceouts += 1
        else:
            ErrorExit('The file %s is not a valid CPOUT or CEOUT file' % (opt.couts[i]))
    
    # Closing all CPOUT and/or CEOUT files
    for i in range(numcouts):
        coutfiles[i].close()
    
    # We will first reconstruct CPOUT files (if provided), and then repeat the same process for CEOUT files (if provided)
    counter = 0
    while True:
        if (counter == 0):
            if (numcpouts>0):
                # We will work with CPOUT files
                typ = 4
            elif (numceouts>0):
                # We will work with CEOUT files
                typ = 5
            else:
                ErrorExit('Something is wrong! The number of CPOUT and CEOUT files provided is zero.')
        elif (counter == 1 and numcpouts>0 and numceouts>0):
            # We will work with CEOUT files
            typ = 5
        else:
            break
        
        remd_types = []
        remd_prefixes = []
        remd_values = []
        if (typ == 4):
           prefix = opt.prefix+".cpout"
        else:
           prefix = opt.prefix+".ceout"
        # Opening only CPOUT or CEOUT files
        coutfls = []
        for i in range(numcouts):
            if (couts_type[i] == typ):
                coutfls.append(open(opt.couts[i], "r"))
        numcoutfls = len(coutfls)
        
        FirstFile = True 
        for i in range(numcoutfls):
            # Reading how many REMD dimensions we have and what are the types
            # from the first line that has residue information
            while True:
                line = coutfls[i].readline()
                if (line):
                    words = line.split()
                    if (words[0] == "Residue" and words[2] == "State:"):
                        if (len(remd_types) != 0): FirstFile = False
                        cnt = 0
                        vector = []
                        for j in range(3,len(words)-1):
                            if (FirstFile):
                                if (words[j] == "T:"):
                                    remd_types.append(1)
                                    remd_prefixes.append(".T_")
                                    vector.append(words[j+1])
                                elif (words[j] == "H:"):
                                    remd_types.append(3)
                                    remd_prefixes.append(".H_")
                                    vector.append(words[j+1])
                                elif (words[j] == "pH:"):
                                    remd_types.append(4)
                                    remd_prefixes.append(".pH_")
                                    vector.append(words[j+1])
                                elif (words[j] == "E:"):
                                    remd_types.append(5)
                                    remd_prefixes.append(".E_")
                                    vector.append(words[j+1])
                                elif (words[j].find(":")!=-1):
                                    ErrorExit("Unknown Replica Exchange dimension '%s' stated on file %s."%(words[j],coutfls[i].name))
                            else:
                                HadError = False
                                if (words[j] == "T:"):
                                    vector.append(words[j+1])
                                    if (remd_types[cnt] != 1): HadError = True
                                    cnt += 1
                                elif (words[j] == "H:"):
                                    vector.append(words[j+1])
                                    if (remd_types[cnt] != 3): HadError = True
                                    cnt += 1
                                elif (words[j] == "pH:"):
                                    vector.append(words[j+1])
                                    if (remd_types[cnt] != 4): HadError = True
                                    cnt += 1
                                elif (words[j] == "E:"):
                                    vector.append(words[j+1])
                                    if (remd_types[cnt] != 5): HadError = True
                                    cnt += 1
                                if (HadError): ErrorExit('The REMD dimension types on the file %s do not match with the dimensions on the file %s' % (opt.couts[i],opt.couts[0]))
                                del HadError
                        if (FirstFile and len(remd_types)==0): ErrorExit('The file %s does not come from a Replica Exchange simulation.'%(coutfls[i].name))
                        for j in range(len(remd_values)):
                            if (remd_values[j] == vector): ErrorExit('Duplicate Replica Exchange values on the files %s and %s !'%(coutfls[j].name,coutfls[i].name))
                        remd_values.append(vector)
                        del vector
                        # Defining the characther positions from which we will read values of pH, Redox Potential, Temperature, etc
                        if (FirstFile):
                            charpos = [[0,0] for j in range(len(remd_types))]
                            cnt = 0
                            # Searching for ":" characters positions
                            for j in range(22,len(line)):
                                if (line[j]==":"):
                                    charpos[cnt][0] = j+1
                                    cnt+=1
                            # Searching for ending characters positions
                            for j in range(cnt):
                                BlankChar = True
                                for k in range(charpos[j][0],len(line)):
                                    if (not line[k]==" "): BlankChar = False
                                    if (not BlankChar and (line[k]==" " or k == (len(line)-1))):
                                        charpos[j][1] = k
                                        break
                            del BlankChar
                        break
                else:
                    ErrorExit('The file %s does not contain any residue information.' % (opt.couts[i]))
            coutfls[i].close()
        del FirstFile
        
        # Number of REMD dimensions
        ndims = len(remd_types)
        
        # Opening the new blank sorted CPOUT or CEOUT files. This step is necessary to insure any existing file will be overwritten
        for i in range(numcoutfls):
            filename = prefix
            for j in range(ndims):
                filename += remd_prefixes[j]+remd_values[i][j]
            if (os.path.exists(filename) and not opt.overwrite):
                ErrorExit('The file %s already exists! Use the flag -O if you want it to be overwritten.'%(filename))
            else:
                tempscoutfl=open(filename, "w")
                tempscoutfl.close()
                del tempscoutfl
        
        # Opening all the new blank sorted CPOUT or CEOUT files as appended
        scoutfl = []
        scoutfl_filename = []
        for i in range(numcoutfls):
            filename = prefix
            for j in range(ndims):
                filename += remd_prefixes[j]+remd_values[i][j]
            scoutfl.append(open(filename, "a"))
            scoutfl_filename.append(filename)
        del remd_values
        
        # Reopening CPOUT or CEOUT files so we read them starting from the first line again
        coutfls = []
        for i in range(numcouts):
            if (couts_type[i] == typ): coutfls.append(open(opt.couts[i], "r"))
        numcoutfls = len(coutfls)
        
        # Setting some important variables
        AllFinished = False
        AmIFinished = [0 for i in range(numcoutfls)]
        
        # Writting the sorted CPOUT or CEOUT files
        while True:
            for i in range(numcoutfls):
                # Read file until we reach a line that has REMD information
                final_line = ""
                while True:
                    line = coutfls[i].readline()
                    # If we didn't reach the end of the file
                    if (line):
                        if (len(line)==1):
                            final_line += line
                        elif (line[:7] == "Residue"):
                            # Setting file name
                            filename = prefix
                            for j in range(ndims):
                                filename+=remd_prefixes[j]+line[charpos[j][0]:charpos[j][1]].lstrip()
                            # Chunking the line to be printted
                            final_line += line[0:22]+"\n"
                            # Printing lines to sorted file
                            cnt = scoutfl_filename.index(filename)
                            scoutfl[cnt].write(final_line)
                            break
                        else:
                            # Accumulating lines
                            final_line += line
                    # If the end of the file was reached
                    else:
                        AmIFinished[i] = 1
                        # Checking if all files were closed
                        if (AmIFinished.count(0) == 0): AllFinished = True
                        break
                if (AllFinished): break
            if (AllFinished): break
        
        # Closing all REMD CPOUT or CEOUT files
        for i in range(numcoutfls):
            coutfls[i].close()
            
        # Increasing the counter to let the program know we finished to reconstruct either CPOUT or CEOUT files
        counter += 1
      
if __name__ == '__main__':
    opt = parser.parse_args()

    # Go ahead and execute the program.
    main(opt)
    sys.exit(0)
