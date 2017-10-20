#!/usr/bin/env python

# usage: ./run_test_acpype_db_ligands.py (opt: n=10 or "['001','Rca']")

import sys, os, time

from subprocess import Popen

numCpu = 20

def elapsedTime(seconds, suffixes = ['y', 'w', 'd', 'h', 'm', 's'], add_s = False, separator = ' '):
    """
    Takes an amount of seconds and turns it into a human-readable amount of time.
    """
    # the formatted time string to be returned
    if seconds == 0:
        return '0s'
    time = []

    # the pieces of time to iterate over (days, hours, minutes, etc)
    # - the first piece in each tuple is the suffix (d, h, w)
    # - the second piece is the length in seconds (a day is 60s * 60m * 24h)
    parts = [(suffixes[0], 60 * 60 * 24 * 7 * 52),
          (suffixes[1], 60 * 60 * 24 * 7),
          (suffixes[2], 60 * 60 * 24),
          (suffixes[3], 60 * 60),
          (suffixes[4], 60),
          (suffixes[5], 1)]

    # for each time piece, grab the value and remaining seconds, and add it to
    # the time string
    for suffix, length in parts:
        value = seconds / length
        if value > 0:
            seconds = seconds % length
            time.append('%s%s' % (str(value),
                           (suffix, (suffix, suffix + 's')[value > 1])[add_s]))
        if seconds < 1:
            break

    return separator.join(time)

def runConversionJobs(chemCompVarFiles, scriptName):

    _timeStamp = time.strftime("%Y_%m_%d_%H_%M_%S")

    startCode = 'start'

    currentChemCompVarFile = None
    currentProcesses = {startCode: None}
    currentJobOut = {}
    currentIndex = -1
    endChemCompVarFile = chemCompVarFiles[-1]

    outputHandle = sys.__stdout__

    while (currentProcesses):

        if startCode in currentProcesses.keys():
            del(currentProcesses[startCode])

        if len(currentProcesses.keys()) < numCpu:

            tempIndex = currentIndex + 1
            for _i in range(currentIndex, currentIndex + numCpu - len(currentProcesses.keys())):
                # Don't start a job if it's at the end!
                if currentChemCompVarFile != endChemCompVarFile:

                    chemCompVarFile = chemCompVarFiles[tempIndex]

                    #
                    # TODO: ensure that stdout and stdin OK for this job!! Might want to reroute!!
                    #

                    currentOutFile = chemCompVarFile.replace('.mol2', '.out')
                    jobOut = open(currentOutFile, 'w')

                    varDir, varFile = os.path.split(chemCompVarFile)
                    os.chdir(varDir)
                    process = Popen(['nice', '-19', scriptName, '-i', varFile, '-d'], stdout = jobOut, stderr = jobOut)

                    currentJobOut[chemCompVarFile] = jobOut
                    currentProcesses[chemCompVarFile] = process
                    currentChemCompVarFile = chemCompVarFile

                    outputHandle.write("\n*** Job %s started ***\n\n" % chemCompVarFile)

                    tempIndex += 1

                    # Allow jobs to start one by one... nicer for output
                    time.sleep(1)

            currentIndex = tempIndex - 1

        time.sleep(3)

        #
        # Check if finished
        #

        for chemCompVarFile in currentProcesses.keys():

            # Finished...
            if currentProcesses[chemCompVarFile].poll() != None:
                del(currentProcesses[chemCompVarFile])
                currentJobOut[chemCompVarFile].close()
                del(currentJobOut[chemCompVarFile])
                outputHandle.write("\n *** Job %s finished ***\n" % chemCompVarFile)
                if not currentProcesses:
                    currentProcesses = {startCode: None}

        outputHandle.flush()

if __name__ == '__main__':

    t0 = time.time()

    chemCompVarFiles = []
    curDir = os.getcwd()
    ccpCodes = os.listdir('other')
    ccpCodes.sort()

    # Only run this on 'other'!
    if len(sys.argv) > 1:
        args = sys.argv[1:]
        if args[0].startswith('n='):
            num = int(args[0][2:])
            ccpCodes = ccpCodes[:num]
        else:
            args = list(set(eval(args[0])))
            args.sort()
            ccpCodes = args

    for ccpCode in ccpCodes:
        # HACK to restart after powercut
        #if ccpCode < 'NA':
        #  continue

        ccvNames = os.listdir(os.path.join('other', ccpCode))

        for ccvName in  ccvNames:
            if ccvName[-5:] == '.mol2' and not ccvName.count("bcc_gaff"):
                chemCompVarFiles.append(os.path.join(curDir, 'other', ccpCode, ccvName))

    runConversionJobs(chemCompVarFiles, 'acpype')
    execTime = int(round(time.time() - t0))
    msg = elapsedTime(execTime)
    print("Total time of execution: %s" % msg)
    print("ALL DONE")

# nohup ./run_test_acpype_db_ligands.py &
# grep started nohup.out | wc -l ; grep finished nohup.out | wc -l # 17405 17405
# grep -r "FAILED: semi-QM taking" other/*
# grep -r "invalid" other/*
