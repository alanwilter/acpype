#!/home/conda/feedstock_root/build_artifacts/ambertools_1635195478443/_h_env_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_/bin/python

import argparse as ap
import numpy as np
import sys
import f90nml
import math
from scipy.linalg import lu_factor, lu_solve

###### parameters ######
maxq = 8000
maxlgr = 900

###### iostuf ######
inopt,ioutopt,iqopt,iunits = 0,0,1,0

###### files ######
Input,output,qin,qout,punch,espot,qwts,esout = [None]*8
input_file,output_file,qin_file,qout_file,punch_file,espot_file,qwts_file,esout_file = [None]*8

###### infoa ######
nat,iuniq,nesp,natpl1,ihfree,irstrnt = 0,0,0,0,1,1

###### runlab ######
title = None

###### espcom ######
apot,awt,bpot,bwt = [None]*4
ssvpot,chipot,vavrg = 0.0,0.0,0.0

###### calcul ###### 
qcal,a,b,qwtval,iqcntr = [None]*5

###### lagrng ###### 
grpchg = np.zeros((maxlgr))
lgrcnt = np.zeros((maxlgr,maxq), dtype=int)
nlgrng = 0

###### orig ###### 
q0, crd, izan, ivary = [None]*4
qwt = 0.0005

###### worker ###### 
awork,bwork = [None]*2

###### propty ###### 
#cmas = np.ndarray((3))
dipol = np.zeros((3))
quad = np.zeros((6))
dipmom = 0.0

###### mltmol ###### 
wtmol,ibeg,iend = [None]*3
nmol = 1

def file_in():
	parser = ap.ArgumentParser(usage='py_resp [-O] -i input -o output -p punch -q qin -t qout -e espot -w qwts -s esout')
	parser.add_argument("-O", action='store_true', help="Overwrite output files if they exist")
	parser.add_argument("-i", "--input", required=True, help="type: input, required; description: input options")
	parser.add_argument("-e", "--espot", required=True, help="type: input, required; description: input of MEP and coordinates")
	parser.add_argument("-q", "--qin", help="type: input, optional; description: replacement charges")
	parser.add_argument("-w", "--qwts", help="type: input, optional; description: input of new weight factors")
	parser.add_argument("-o", "--output", required=True, help="type: output, always produced; description: output of results")
	parser.add_argument("-p", "--punch", required=True, help="type: output, always produced; description: synopsis of results")
	parser.add_argument("-t", "--qout", required=True, help="type: output, always produced; description: output of current charges")
	parser.add_argument("-s", "--esout", help="type: output, optional; description: generated MEP values for new charges")
	
	args = parser.parse_args()
	Input,espot,qin,qwts = args.input,args.espot,args.qin,args.qwts
	output,punch,qout,esout = args.output,args.punch,args.qout,args.esout

	return Input,output,qin,qout,punch,espot,qwts,esout

def read_in():
	global inopt,ioutopt,iuniq,nmol,iqopt,irstrnt,ihfree,qwt
	global title,wtmol,ibeg,iend,izan,ivary,nlgrng,q0

	# start of molecule input
	output_file.write('\n -----------------------------------------------')
	output_file.write('\n      Py_RESP Alpha Version  ')
	output_file.write('\n -----------------------------------------------')
	output_file.write('\n '+input_file.readline())
	output_file.write(' -----------------------------------------------\n\n')

	# read in charge, number of charge centers, and control parameters
	if 'cntrl' in f90nml.read(Input):
		cntrl_nml = f90nml.read(Input)['cntrl']
		inopt = cntrl_nml['inopt'] if 'inopt' in cntrl_nml else 0 
		ioutopt = cntrl_nml['ioutopt'] if 'ioutopt' in cntrl_nml else 0
		iuniq = cntrl_nml['iuniq'] if 'iuniq' in cntrl_nml else 0
		nmol = cntrl_nml['nmol'] if 'nmol' in cntrl_nml else 1
		iqopt = cntrl_nml['iqopt'] if 'iqopt' in cntrl_nml else 1
		irstrnt = cntrl_nml['irstrnt'] if 'irstrnt' in cntrl_nml else 1
		ihfree = cntrl_nml['ihfree'] if 'ihree' in cntrl_nml else 1
		qwt = cntrl_nml['qwt'] if 'qwt' in cntrl_nml else 0.0005
	else:
		output_file.write('Sorry, you must use namelist input\n')
		sys.exit()

	output_file.write('\n inopt       = %d   ioutopt     = %d'%(inopt, ioutopt))
	output_file.write('\n nmol        = %d   iqopt       = %d'%(nmol, iqopt))
	output_file.write('\n ihfree      = %d   irstrnt     = %d'%(ihfree, irstrnt))
	output_file.write('\n iunits      = %d   qwt         = %.8f'%(iunits, qwt))

	for line in input_file:
		if "&end" in line:
			break;

	# if nmol > 1, this is a multiple molecule run
	# and so we should leave this function now because function mult_mol
	# is responsible for the rest of the multiple-molecule reading in.
	if nmol > 1:
		return

	# read in fitting weight for q0 and esp point weighting
	wtmol = np.ndarray((1))
	wtmol[0] = float(input_file.readline())
	output_file.write("\n wtmol[0] = %12.6f"%wtmol[0])
	title = input_file.readline()
	output_file.write('\n subtitle: \n\n%s'%title)
	line = input_file.readline().split()
	ich, iuniq = int(line[0]), int(line[1])
	output_file.write('ich = %3d  iuniq = %3d'%(ich, iuniq))

	ibeg = np.ndarray((1), dtype=int)
	iend = np.ndarray((1), dtype=int)
	ibeg[0] = 0
	iend[0] = iuniq-1
	wtmol[0] = 1.0

	# read in atomic number izan[i] and ivary[i]
	izan = np.ndarray((iuniq), dtype=int)
	ivary = np.ndarray((iuniq), dtype=int)
	for i in range(iuniq):
		line = input_file.readline().split()
		izan[i], ivary[i] = int(line[0]), int(line[1])
		output_file.write("\n%5d%5d%5d"%(i+1, izan[i], ivary[i]))

	# read in lagrange (charge) constraints
	lagrange(ich, ibeg[0], iend[0], nmol)

	# replacement initial charges q0 from unit iuniq if iqopt=2
	q0 = np.ndarray((iuniq))
	if iqopt > 1:
		qin_file = open(qin, 'r')
		output_file.write("\n new q0 values to be read %4d"%iuniq)

		q0_index = 0
		for line in qin_file:
			q0_list = line.split()
			for q0_val in q0_list:
				q0[q0_index] = float(q0_val)
				q0_index += 1

		qin_file.close()

		# now is a simple trap for when not enough replacement q0 are given.
		# tactic: just keep going with old q0 (assume iqopt=2 was not intended)
		if q0_index < iuniq-1:
			output_file.write("\n  not enough (possibly none) q0 are given in file")
			output_file.write("\n  ESP.Q0, so the remaining old ones will be use.")
	# end of "replacement initial charge" section.

	# set initial charges to 0; done if iqopt=1
	if iqopt == 1:
		output_file.write("\n iqopt=1, all q0 values will be set to 0")
		q0.fill(0.0)

	output_file.write("\n\n  ----------------------------------------------------------------------------")
	output_file.write("\n     ATOM                       COORDINATES                         CHARGE")
	output_file.write("\n                       X             Y                Z")
	output_file.write("\n  ----------------------------------------------------------------------------")
	output_file.write("\n  ----------------------------------------------------------------------------\n")

	output_file.write("\n\n Charge on the molecule(ich) =%5d\
		\n Total number of atoms (iuniq)     =%5d\
		\n Weight factor on initial charge restraints(qwt)=%16.5f\n"%(ich, iuniq, qwt))

	# charge constraint info
	output_file.write("\n\n  there are %3d charge constraints:\n\n"%nlgrng)
	for i in range(iuniq):
		output_file.write("\n %5d     "%(i+1))
		for j in range(nlgrng):
			output_file.write("%3d"%(lgrcnt[j,i]))

	if iuniq > maxq:
		output_file.write("\nNumber of atoms exceeds program dimensions")
		sys.exit()

def mult_mol():
	######################################################################
	# this function reads in multiple molecule input. In function readin
 	# it has already read the control variable input for the run, namely:
 	# ich, iuniq, inopt, iqopt, ihfree, irstrnt, nlgrng, nmol
 	# where nmol > 1 caused this function to be called.
 	# The input form for the other molecules is that their entire control
 	# decks are appended to the initial 2 control lines just read in,
 	# each control deck separated by a blank line, and then comes the
 	# multiple-molecule specific input, which is
 	#  - equivalencing of centers between molecules in the series
 	#  - the lagrange constraints to be applied between molecules
 	#
 	# the control characters read in the individual job decks are ignored
 	# except for ich, and icntrs. The lagrange (charge) constraints
 	# contained in the individual-molecule inputs ARE included.
 	#
 	# NOTE: the following are NOT implemented:
 	#       - output QM/calculated esp's (ioutOPT=1)
 	######################################################################
	global nlgrng,iuniq,title,wtmol,ibeg,iend,izan,ivary,q0

	imol = 0                # local variable
	iuniq = 0

	output_file.write("\n\n %%RESP-I-MULT_MOL,  multiple-molecule run of %3d molecules"%nmol)

	wtmol = np.ndarray((nmol))
	ibeg = np.ndarray((nmol), dtype=int)
	iend = np.ndarray((nmol), dtype=int)
	izan = np.ndarray((1), dtype=int)
	ivary = np.ndarray((1), dtype=int)

	for imol in range(nmol):
		# read in the molecule weight and control variables
		# the lagrange constraints to be applied between molecules
		wtmol[imol] = float(input_file.readline())
		output_file.write("\n\n Reading input for molecule %3d weight:%10.3f\n"%(imol+1,wtmol[imol]))
		title = input_file.readline()
		output_file.write(title)

		# read in charge, number of charge centers, and control parameters
		line = input_file.readline().split()
		ich, icntrs = int(line[0]), int(line[1])
		output_file.write('\n Total charge (ich):%3d'%ich)
		output_file.write('\n Number of centers:%3d'%icntrs)

		# read in fitting weight for q0 and esp point weighting

		# now some book-keeping: iuniq is the global variable for the total
		# number of centers over all molecules. The first center of this
		# mol therefore starts in iuniq and goes to iuniq+icntrs-1.
		ibeg[imol] = iuniq
		iend[imol] = iuniq+icntrs-1

		# trap for having too many centers
		if iend[imol]+1 > maxq:
			output_file.write('\n ERROR: more than %5d centers'%maxq)
			sys.exit()

		# Read in atomic number izan[i] and ivary[i]
		# Since ivary[i] is supposed to correspond to a center-number in the
		# same molecule, this has to be adjusted to ivary[i]+ibeg[imol]
		# convert angstroms to bohrs if necessary
		izan.resize((iuniq+icntrs), refcheck=False)
		ivary.resize((iuniq+icntrs), refcheck=False)
		for i in range(ibeg[imol], iend[imol]+1):
			line = input_file.readline().split()
			izan[i], ivary[i] = int(line[0]), int(line[1])
			output_file.write("\n%5d%5d%5d"%(i+1, izan[i], ivary[i]))
			if ivary[i] > 0:
				ivary[i] += ibeg[imol]

		# now set iuniq to iuniq+icntrs
		iuniq += icntrs

		# now read in the lagrange (charge) constraints for this molecule
		lagrange(ich, ibeg[imol], iend[imol], imol)

	# end of molecule input, now do other preparation stuff

	# read past a blank line after the final molecule job deck and then
	# read in inter-molecule lagrange (charge) constraints (mmlgr).
	# The "-99" for the total charge tells lgrange to drop the total charge
	# constraint
	lagrange(-99, 0, iuniq-1, imol)

	# this is for different types of charge resetting according to iqopt :

	# replacement initial charges q0 from qin if iqopt=1
	q0 = np.ndarray((iuniq))
	if iqopt > 1:
		qin_file = open(qin, 'r')
		output_file.write("\n since IQOPT=1, %4d new q0 values"%iuniq)
		output_file.write("\nwill be read in from file ESP.Q0 (unit 3)")

		# now read in replacement charges
		q0_index = 0
		for line in qin_file:
			q0_list = line.split()
			for q0_val in q0_list:
				q0[q0_index] = float(q0_val)
				q0_index += 1

		qin_file.close()

		# now is a simple trap for when not enough replacement q0 are given.
		# tactic: just keep going with old q0 (assume iqopt=2 was not intended)
		if q0_index < iuniq-1:
			output_file.write("\n  not enough (possibly none) q0 are given in file")
			output_file.write("\n  ESP.Q0, so the remaining old ones will be use.")
	# end of "replacement initial charge" section.

	# begin section: setting initial charges to 0; done if iqopt=1
	if iqopt == 1:
		output_file.write("\n iqopt=1, all q0 values will be set to 0")
		q0.fill(0.0)

	# now carry out the inter-molecule equivalencing.  This is done by
	#
	# First : read the cards saying how many centers will be read in in the
	#         next card. a zero means we have finished input
	#
	# Second: read the first-occurrence-in-each-molecule of the centers to
	#        be equivalenced.
	#
	#        The specifcations MUST be in ascending order.
	#
	#        The expanding of the centers within each
	#        molecule is based on the ivary values for the individual mol.
	#
	#        if ivary for mol 2+ is zero it is replaced with the atom number
	#        of mol 1.
	output_file.write("\n --------------------------------")
	output_file.write("\n reading mult_mol constraint info")
	output_file.write("\n --------------------------------")

	line = input_file.readline().split()
	while line:  # line is not empty
		ntmp1 = int(line[0])
		itmp = np.ndarray((ntmp1), dtype=int)
		imoll = np.ndarray((ntmp1), dtype=int)

		cnt = ntmp1
		for row in range(math.ceil(ntmp1/8)):
			line = input_file.readline().split()
			output_file.write("\n")
			for j in range(8):
				index = row*8+j
				imoll[index] = int(line[2*j])
				itmp[index] = int(line[2*j+1])
				output_file.write("%5d%5d"%(imoll[index], itmp[index]))
				cnt -= 1
				if cnt == 0:
					break

		for i in range(ntmp1):
			jmol = imoll[i] - 1  # -1 since python is 0 indexed
			icntrs = ibeg[jmol]
			itmp[i] += icntrs

		for i in range(1, ntmp1):
			ivary[itmp[i]-1] = itmp[0]

		line = input_file.readline().split()

	for i in range(iuniq):
		ntmp1 = ivary[i]
		if ntmp1 > 0:
			ntmp2 = ivary[ntmp1-1]
			if ntmp2 > 0:
				ivary[i] = ntmp2

	output_file.write("\n\n  --------------------")
	output_file.write("\n     Atom   Ivary")
	output_file.write("\n  --------------------")
	icnt = 0
	jcnt = 0
	for iat in range(iuniq):
		output_file.write("\n%5d%5d"%(izan[iat], ivary[iat]))
		jcnt += 1
		if (jcnt > iend[icnt]):
			output_file.write("\n")
			icnt += 1

	output_file.write("\n  ----------------------------------------------------------------------------\n")
	output_file.write("\n\n Total number of atoms      =%5d"%iuniq)
	output_file.write("\n Weight factor on initial charge restraints=%10.6f\n"%qwt)

	# charge constraint info
	output_file.write("\n\n There are%3d charge constraints"%nlgrng)

def lagrange(ncharge, ifirst, ilast, imol):
	###################################################
	# read in and assign lagrange constraint pointers #
	# called from "readin" and "mult_mol"             #
	###################################################
	global nlgrng

	line = input_file.readline().split()
	if line:  # line is not empty
		while line:
			ntmp, gtemp = int(line[0]), float(line[1])
			nlgrng += 1
			if nlgrng > maxlgr:
				output_file.write("\n Too many charge-group constraints: %3d"%nlgrng)
				output_file.write("\n Maximum allowed: %3d"%maxlgr)
				sys.end()

			grpchg[nlgrng-1] = gtemp
			itmp = np.ndarray((ntmp), dtype=int)
			imoll = np.ndarray((ntmp), dtype=int)

			cnt = ntmp
			for row in range(math.ceil(ntmp/8)):
				line = input_file.readline().split()
				output_file.write("\n")
				for j in range(8):
					index = row*8+j
					imoll[index] = int(line[2*j])
					itmp[index] = int(line[2*j+1])
					output_file.write("%5d%5d"%(imoll[index], itmp[index]))
					cnt -= 1
					if cnt == 0:
						break

			for i in range(ntmp):
				jmol = imoll[i] - 1  # -1 since python is 0 indexed
				icntrs = ibeg[jmol]
				itmp[i] += icntrs

			for j in range(ntmp):
				if itmp[j] > 0:
					itmp[j] += ifirst
					lgrcnt[nlgrng-1][itmp[j]-1] = 1 # minus 1 since python is 0 indexed
				elif itmp[j] < 0:
					itmp[j] -= ifirst
					lgrcnt[nlgrng-1][itmp[j]-1] = -1 # minus 1 since python is 0 indexed
			line = input_file.readline().split()

	# as long as ncharge is not -99, implement the "total charge" constraint
	if ncharge > -99:
		nlgrng += 1
		if nlgrng > maxlgr:
			output_file.write("\n Too many charge-group constraints: %3d"%nlgrng)
			output_file.write("\n Maximum allowed: %3d"%maxlgr)
			sys.end()

		grpchg[nlgrng-1] = float(ncharge)
		for j in range(ifirst, ilast+1):
			lgrcnt[nlgrng-1][j] = 1

def matpot():
	# read in the electrostatic potential points used in the fitting,
	# building up as we go the matrices for LU decomposition
	#
	# called from Main
	global apot, awt, bpot, bwt, crd, ssvpot, vavrg, nesp
	
	apot = np.zeros((iuniq,iuniq))
	awt = np.zeros((iuniq,iuniq))
	bpot = np.zeros((iuniq))
	bwt = np.zeros((iuniq))

	espot_file = open(espot, 'r')

	vavtmp = 0.0 # local variable
	ioff = 0     # local variable

	if nmol > 0:
		inmol = nmol
	else:
		inmol = 1
		ibeg[0] = 0
		iend[0] = iuniq
		wtmol[0] = 1.0

	crd = np.zeros((3,iuniq))

	for imol in range(inmol):
		line = espot_file.readline().split()
		inat, nesp = int(line[0]), int(line[1])
		output_file.write("\n\n Reading esp's for molecule %3d"%(imol+1))
		output_file.write("\n total number of atoms      = %5d"%inat)
		output_file.write("\n total number of esp points = %5d"%nesp)

		output_file.write("\n\n\n center     X       Y       Z")
		for i in range(inat):
			line = espot_file.readline().split()
			crd[0][ioff], crd[1][ioff], crd[2][ioff] = float(line[0]),float(line[1]),float(line[2])
			output_file.write("\n {:4d}{:16.7E}{:16.7E}{:16.7E}".format(i+1, crd[0][ioff], crd[1][ioff], crd[2][ioff]))
			ioff += 1

		# build up matrix elements Ajk according to (SUMi 1/Rik SUMj 1/Rij)
		for i in range(nesp):
			line = espot_file.readline().split()
			if not line:
				output_file.write("     premature end of potential file")
				sys.end()
			espi, xi, yi, zi = float(line[0]),float(line[1]),float(line[2]),float(line[3])
			wt = wtmol[imol]
			wt2 = wt*wt
			vavtmp += wt*espi
			ssvpot += wt2*espi*espi
			vavrg += vavtmp/nesp
			for k in range(ibeg[imol], iend[imol]+1):
				rik = 1/math.sqrt((xi-crd[0][k])**2 + (yi-crd[1][k])**2 + (zi-crd[2][k])**2)
				bpot[k] += espi*rik
				bwt[k] += wt2*espi*rik
				apot[k][k] += rik*rik
				awt[k][k] += wt2*rik*rik
				for j in range(k+1, iend[imol]+1):
					rij = 1/math.sqrt((xi-crd[0][j])**2 + (yi-crd[1][j])**2 + (zi-crd[2][j])**2)
					apot[j][k] += rij*rik
					awt[j][k] += wt2*rij*rik

	# symmetrize the potenitial and weighted potential matrices
	for k in range(iuniq):
		for j in range(k+1, iuniq):
			awt[k][j] = awt[j][k]
			apot[k][j] = apot[j][k]

	espot_file.close()
	output_file.write("\n Initial ssvpot = %10.3f"%ssvpot)

def data_prep():
	# setup pointers for groups of charges  based on "ivary" info
	#
	# called from Main

	#************************************************************************
	# begin section: set lists for combined and frozen charges
	#
	# ivary[i] = 0, it is a new charge center to be fitted
	# ivary[i] =+n, it is a charge center to be fitted with center n
	#                    (center n must be a previous center entered with
	#                    ivary[n] = 0
	# ivary[i] =-n, it is a frozen charge center to be kept at q0[i]
	#************************************************************************
	global iqcntr, nat, natpl1

	iqcntr = np.zeros((iuniq+nlgrng), dtype=int)
	for i in range(iuniq):
		if ivary[i] == 0:
			nat += 1
			iqcntr[i] = nat
		elif ivary[i] > 0:
			iqcntr[i] = iqcntr[ivary[i]-1]

			if iqcntr[i] > nat:
				output_file.write("\n data_prep: charge vary input is screwy")
				sys.end()
		else:
			iqcntr[i] = -1

	output_file.write("\n\n\n Number of unique UNfrozen centers=%5d"%nat)
	if nat == 0:
		output_file.write("\n ALL charges are frozen!!!")

	# finish off list with Lagrange constraints
	for i in range(nlgrng):
		iqcntr[iuniq+i] = nat+i+1

	# set natpl1 to the total number of row elements
	# ( charges to be independantly fit + constraints )
	# in fitting matrix
	natpl1 = nat + nlgrng

	# done adding Lagrange constraints to elements list

	# read in charges must now be averaged
	# a posteriori averaging of replacement charges according
	# to current ivary charge-combining pointers
	if iqopt == 3:
		for i in range(iuniq-1):
			qcntrs = q0[i]
			tmpctr = 1.0
			for j in range(i+1, iuniq):
				if ivary[j] == i:
					qcntrs += q0[j]
					tmpctr += 1.0

			if tmpctr > 0.99:
				qcntrs /= tmpctr
				q0[i] = qcntrs
				for j in range(i+1,iuniq):
					if ivary[j] == i:
						q0[j] = qcntrs

def cycle(icycle):
	nqwt = None
	# if inopt=1, reads values new of qwt from unit 4 and
	# writes a summary of the resulting fit,
	# as well as a list of the fit charges.
	if inopt < 1:
		return 0

	qwts_file = open(qwts, 'r')
	line = qwts_file.readline().split()

	if icycle == -1:
		nqwt = 0

		# first pass
		nqwt = int(line[0])
		if nqwt == 0:
			output_file.write("\n subroutine icycle: INPUT ERROR...")
			output_file.write("\n INOPT=1 so qwt cycling expected, but reading non-zero nqwt failed")
			return
		qwt = float(line[1])
		icycle = 1
		output_file.write("\n cycle   1: weighting factor=%10.4f"%qwt)
	else:
		# rest of the time
		nqwt = int(line[0])
		icycle += 1
		output_file.write("\n cycle%4d: weighting factor=%10.4f"%(icycle,qwt))
		if icycle >= nqwt:
			icycle = 0
			qwts_file.close()

	return icycle

def charge_opt():
	# driver for the charge determinization/optimizaton
	#
	# called from Main
	global irstrnt, awork, bwork

	qold = np.zeros((iuniq)) # local variable

	irsave = 0 # local variable
	nitern = 0 # local variable

	# qtol & maxit are criteria for convergence & maximum iterations for
	# the non-linear optimizations.
	qtol = 0.000001
	maxit = 24

	# only on first pass through this function (indicated by nitern= 0),
	# if irstrnt > 0, transfer irstrnt to irsave and reset irstrnt to 0,
	# in order to get an initial guess using a harmonic constraint.  This is
	# done so restraint function rstran() will use a harmonic restraint.
	if irstrnt > 0:
		irsave = irstrnt
		irstrnt = 0
		output_file.write("\n\n Non-linear optimization requested.")

	# now go do a "harmonic restraint" run, restraint= qwt(qcal[i]-q0[i])**2
	# -- loop to convergence

	while nitern < maxit:
		matbld()
	
		# solve (Ax = b) where A and b are input, x is solution
		#              awork x = bwork
		# 
		# the solution "x" is returned in "b" (bwork)
		# 
		# -- condition the matrix diagonal to avoid DGETRF() detectingNN (wrapped in lu_factor)
		#    singularity
		for jn in range(natpl1):
			if abs(awork[jn][jn]) < 1.0E-10:
				awork[jn][jn] = 1.0E-10
	
		awork, piv = lu_factor(awork, overwrite_a=True)
		bwork = lu_solve((awork, piv), bwork, overwrite_b=True)
	
		# -- copy solution vector "bwork" to 'calculated charges' vector
		#    qcal
		for k in range(iuniq):
			icntr = iqcntr[k]
			if icntr >= 1:
				# -- new charge
				qcal[k] = bwork[icntr-1]
			else:
				# -- frozen charge
				qcal[k] = q0[k]
	
		# -- a quick check from rstrn: if irstrnt is now negative,
		#    there are no restraints because no qwtval(i) > 0.1e-10,
		#    so reset irsave= 0
		if irstrnt < 0:
			irsave = 0
			output_file.write("\n\n WARNING: Restraints were requested, but the restraint weights were all zero")
	
		# -- we're finished if it's only a "harmonic restraint" run,
		#    but if it's a non-linear optimization (irsave>0)...
		#    we've only just begun (i.e. we have our initial guess)
		if irsave <= 0:
			return
		else:
			# -- it's a non-linear optimization: reset irstrnt (to now
			#    calculate the proper non-linear restraint derivatives
			#    in routine rstran)
			irstrnt = irsave
	
		# -- begin iterative optimization loop with comparison of
		#    old & new charges; calculate the convergence and replace
		#    the old charges with the new
		qchnge = 0.0
		for i in range(iuniq):
			qdiff = qcal[i] - qold[i]
			qchnge += qdiff*qdiff
			qold[i] = qcal[i]
	
		qchnge = math.sqrt(qchnge)/iuniq
		output_file.write("\n qchnge ={:20.10E}".format(qchnge))
	
		# -- if this is less than qtol then we're done
		if qchnge < qtol and nitern > 1:
			output_file.write("\n\n Convergence in%5d iterations"%nitern)
			return

		# loop again
		nitern += 1

	output_file.write("\n after %5d iterations, no convergence!"%maxit)

def matbld():
	# called from "chgopt"
	#
	# build up matrices for LU decomposition:
	#
	#   stage 1: copy weighted matrices awt and bwt to work arrays awork and bwork
	#            (which are destroyed in the LU decomp & back subst)
	#
	#   stage 2: if charge restraints are to be included,
	#            then modify awork and bwork appropriately
	global a, b, awork, bwork

	a = np.zeros((iuniq+nlgrng,iuniq+nlgrng))
	b = np.zeros((iuniq+nlgrng))

	for k in range(iuniq):
		b[k] = bwt[k]
		for j in range(iuniq):
			a[j][k] = awt[j][k]

	# fill in the final columns & rows of A with the Lagrange
	# constraints which keep the charge on groups of atoms to a
	# constant
	#
	# note index counters!
	for i in range(nlgrng):
		b[iuniq+i] = grpchg[i]
		for j in range(iuniq+nlgrng):
			a[iuniq+i][j] = float(lgrcnt[i][j])
			a[j][iuniq+i] = float(lgrcnt[i][j])

	# add restraint to initial charge q0[i]:
	rstran()

	# build awork and bwork based on "combined and frozen centers" info:
	#
	# 1) frozen centers do not appear in the matrix of fitted charges
	# 2) combined centers appear as one single charge center for fitting
	#
	# first, since we accumulate values, zero out awork & bwork up to natpl1
	# (the independant + contraint number):
	awork = np.zeros((natpl1,natpl1))
	bwork = np.zeros((natpl1))

	# loop over all centers, building awork & bwork from A and
	# B based on iqcntr: for each center, iqcntr[i] dictates which of
	# the fitted charges it is and therefore where it goes in the matrices.
	# If iqcntr[j] < 1, this center is a frozen charge and it is skipped as
	# far as forming a row in awork, and its esp contribution is subtracted
	# from bwork to take care of it's awork jth column-element for each i.
	for i in range(iuniq+nlgrng):
		icntr = iqcntr[i]
		if icntr >= 1:
			# i is "active"
			bwork[icntr-1] += b[i]
			for j in range(iuniq+nlgrng):
				jcntr = iqcntr[j]
				if jcntr >= 1:
					# j i is active
					awork[icntr-1][jcntr-1] += a[i][j]
				else:
					# j is a frozen charge
					bwork[icntr-1] -= q0[j]*a[i][j]

def rstran():
	# routine to assign the retraint weights
	# to the diagonal of A and to B
	#
	# called from "matbld"
	global irstrnt, qwtval

	#----------------------------------------------------------------------
	# two kinds of restraint are available:
	#
	# a) a harmonic restraint to the initial charge.  Fine as long as there
	#  aren't any large charges that SHOULD be large... these really feel a
	#  strong force if they are restrained to a low value.
	#
	# b) a hyperbolic restraint to a charge of 0.  This gets asymptotic at
	#  "large" values, so "large" charges aren't pulled down any stronger
	#  than some (reasonable) limiting force.  This is a non-linear
	#  weighting function, so the fit procedure is iterative.
	#
	# other options for restraints to initial charge q0[i]:
	# if requested, restrain the charges by modifying the sum-of-squares
	# cost function derivative.  The scheme for doing this is as follows:
	#
	# if control variable ihfree > 0, let hydrogen charges float free
	#                                   (i.e. reset their qwtval to 0.0).
	#
	#-----------------------------------------------------------------------
	qwtval = np.ndarray((iuniq))
	qwtval.fill(qwt)

	for i in range(iuniq):
		if ihfree > 0 and izan[i] == 1:
			qwtval[i] = 0.0

		if irstrnt == 0:
			a[i][i] += qwtval[i]

			# q0 has the initial and/or frozen charge
			b[i] += qwtval[i]*q0[i]

		elif irstrnt > 0 and qwtval[i] > 0.1E-10:
			# use analytic gradient of the hyperbola

			# qcal has the current (calculated) charge
			qwtval[i] = qwt/math.sqrt(qcal[i]*qcal[i] + 0.01)
			a[i][i] += qwtval[i]

	# if all qwtval[i] are 0.0, no restraints so set irstrnt= -1
	ihit = 0
	for i in range(iuniq):
		if qwtval[i] > 0.1E-10:
			ihit = 1
	if ihit == 0:
		irstrnt = -1

def evlchi():
	# called from Main
	#
	# Evaluate chi-square for linear function  yclci= sum-j(paramj*termij),
	# where j= number of terms, paramj is the coefficient to termij, and
	# chi-square is the merit function: chi-square= sum-i((yi-yclci)**2),
	# where i is the number of data points for which y is known.
	#
	# To avoid going through all i data points every time, the chi-square
	# function is expanded into: chi-square= sum-i(yi**2 - 2*yi*yclci + yclci**2),
	# and re-expressed as sum-i(yi**2) - 2*sum-i(yi*yclci) + sum-i(yclci**2).
	# The first term is calculated once as the data is read in (sum-i(yi**2)== ssy).
	# The second and third term depend upon the parameters: for each parameter
	# paramj, sum-i(yi*termij) is the jth element in xprod, and sum-i(yclci**2) is
	# a row in apot where element ajk= sum-i(termij*termik).  These elements
	# are also built up once as the (possibly large) set of data is read in.
	global chipot

	cross = 0.0   # local variable
	ssycle = 0.0  # local variable

	for j in range(iuniq):
		cross += qcal[j]*bpot[j]
		for k in range(iuniq):
			ssycle += qcal[j]*qcal[k]*apot[j][k]

	chipot = ssvpot - 2.0*cross + ssycle

def elec_mom():
	# routine to calculate the dipole and quadrupole moments
	global dipmom

	bohr = 0.52917725
	debye = 2.541765

	# calculate dipole moment
	for i in range(iuniq):
		for j in range(3):
			dipol[j] += qcal[i]*crd[j][i]

	dipmom = math.sqrt(dipol[0]**2 + dipol[1]**2 + dipol[2]**2)

	# calculate the quadrapole moment
	quadm()

	# convert dipoles from a.u. to debyes, and quadrupoles to debye*angstroms
	for i in range(3):
		dipol[i] *= debye
		quad[i] *= debye*bohr
		quad[i+3] *= debye*bohr
	dipmom *= debye

def quadm():
	# called from "elec_mom"
	#
	# this function calculates the components of the quadropole moment
	for i in range(iuniq):
		x2, y2, z2 = crd[0][i]**2, crd[1][i]**2, crd[2][i]**2
		rq = x2 + y2 + z2
		quad[0] += qcal[i]*(3.0*crd[0][i]*crd[0][i]-rq)
		quad[1] += qcal[i]*(3.0*crd[1][i]*crd[1][i]-rq)
		quad[2] += qcal[i]*(3.0*crd[2][i]*crd[2][i]-rq)
		quad[3] += qcal[i]*(3.0*crd[0][i]*crd[1][i]-rq)
		quad[4] += qcal[i]*(3.0*crd[0][i]*crd[2][i]-rq)
		quad[5] += qcal[i]*(3.0*crd[1][i]*crd[2][i]-rq)

def pun_sum():
	global punch_file
	# called from Main

	qcrtrn = math.sqrt(chipot/ssvpot)  # local variable

	# punch name, one-line summary, and charges
	punch_file = open(punch, 'w')
	punch_file.write("\n%s"%title)
	punch_file.write("\n\niqopt   irstrnt  ihfree     qwt")
	punch_file.write("\n%3d    %3d    %3d    %12.6f"%(iqopt, irstrnt, ihfree, qwt))
	punch_file.write("\n\n rel.rms   dipole mom       Qxx      Qyy      Qzz")
	punch_file.write("\n %10.5f%10.5f%10.5f%10.5f%10.5f"%(qcrtrn, dipmom, quad[0], quad[1], quad[2]))
	punch_file.write("\n\n          Point charges before & after optimization")
	punch_file.write("\n    NO   At.No.    q0           q(opt)   IVARY  d(rstr)/dq")

	icnt, jcnt = 0, 0
	for j in range(iuniq):
		punch_file.write("\n  %4d  %4d %10.6f     %10.6f%5d%12.6f"%(j+1,izan[j],q0[j],qcal[j],ivary[j],qwtval[j]))
		jcnt += 1
		if (jcnt > iend[icnt]):
			output_file.write("\n")
			icnt += 1

	qout_file = open(qout, 'w')

	counter = iuniq
	for i in range(math.ceil(iuniq/8)):
		for j in range(min(8, counter)):
			qout_file.write("%10.6f"%qcal[i*8+j])
		counter -= 8
		qout_file.write("\n")

	qout_file.close()

def pot_out():
	# called from Main

	# ---- print the optimized charges and coordinates ----
	output_file.write("\n\n%s"%title)

	# print the charges
	output_file.write("\n          Point Charges Before & After Optimization")
	output_file.write("\n\n    no.  At.no.    q(init)       q(opt)     ivary    d(rstr)/dq")

	icnt, jcnt = 0, 0
	chge = 0.0
	for j in range(iuniq):
		output_file.write("\n %4d%4d     %10.6f     %10.6f%7d%15.6f"%(j+1,izan[j],q0[j],qcal[j],ivary[j],qwtval[j]))
		chge += qcal[j]
		jcnt += 1
		if (jcnt > iend[icnt]):
			output_file.write("\n")
			icnt += 1
	output_file.write("\n Sum over the calculated charges: %10.3f"%chge)

	qcrtrn = math.sqrt(chipot/ssvpot)  # local variable

	# calculate standard error of estimate
	sigma = math.sqrt(chipot/nesp)

	# now write all these stuff out
	output_file.write("\n\n        Statistics of the fitting:")
	output_file.write("\n  The initial sum of squares (ssvpot)      %15.3f"%ssvpot)
	output_file.write("\n  The residual sum of squares (chipot)     %15.3f"%chipot)
	output_file.write("\n  The std err of estimate (sqrt(chipot/N)) %15.5f"%sigma)
	output_file.write("\n  ESP relative RMS (SQRT(chipot/ssvpot))   %15.5f"%qcrtrn)

	punch_file.write("\n\n        Statistics of the fitting:")
	punch_file.write("\n  The initial sum of squares (ssvpot)      %15.3f"%ssvpot)
	punch_file.write("\n  The residual sum of squares (chipot)     %15.3f"%chipot)
	punch_file.write("\n  The std err of estimate (sqrt(chipot/N)) %15.5f"%sigma)
	punch_file.write("\n  ESP relative RMS (SQRT(chipot/ssvpot))   %15.5f"%qcrtrn)

	# ----- print the dipole, quadrupole and center of mass ----
	if nmol == 1:
		#output_file.write("\n\n Center of Mass (Angst.):")
		#output_file.write("\n\n X  =%10.5f      Y  =%10.5f      Z  =%10.5f"%(cmas[0], cmas[1], cmas[2]))
		output_file.write("\n\n Dipole (Debye):")
		output_file.write("\n\n X  =%10.5f      Y  =%10.5f      Z  =%10.5f"%(dipol[0], dipol[1], dipol[2]))
		output_file.write("\n\n Dipole Moment (Debye)=%10.5f"%dipmom)
		punch_file.write("\n\n Dipole Moment (Debye)=%10.5f"%dipmom)
		output_file.write("\n\n Quadrupole (Debye*Angst.):")
		output_file.write("\n\n QXX =%10.5f     QYY =%10.5f     QZZ =%10.5f"%(quad[0], quad[1], quad[2]))
		output_file.write("\n QXY =%10.5f     QXZ =%10.5f     QYZ =%10.5f"%(quad[3], quad[4], quad[5]))

def wrt_pot():
	# called from Main
	#
	# read in the electrostatic potential points used in the fitting,
	# calculate esp using existing charges, and write out both esp's & residual

	# open the file containing the qm esp points & read in the no. of points
	au2cal = 627.5095
	bohr = 0.52917725
	espot_file = open(espot, 'r')
	esout_file = open(esout, 'w')

	line = espot_file.readline().split()
	idum, nesp = int(line[0]), int(line[1])
	ssvkcl = ssvpot*au2cal*au2cal
	chikcl = chipot*au2cal*au2cal
	esout_file.write("%5d%6d    %20.10f%20.10f"%(nesp, izan[0], ssvkcl, chikcl))

	for i in range(idum):
		line = espot_file.readline().split()

	for i in range(nesp):
		line = espot_file.readline().split()
		if not line:
			output_file.write("\n unexpected eof in %s"%espot)
			sys.exit()
		espqmi, xi, yi, zi = float(line[0]),float(line[1]),float(line[2]),float(line[3])
		espclc = 0.0
		for k in range(iuniq):
			xik = xi - crd[0][k]
			yik = yi - crd[1][k]
			zik = zi - crd[2][k]
			espclc += qcal[k]/math.sqrt(xik**2 + yik**2 + zik**2)
		vresid = espqmi - espclc
		xa = xi*bohr
		ya = yi*bohr
		za = zi*bohr
		espclc *= au2cal
		espqmi *= au2cal
		vresid *= au2cal
		esout_file.write("\n%10.5f%10.5f%10.5f%12.5f%12.5f%12.5f"%(xa,ya,za,espqmi,espclc,vresid))

	espot_file.close()
	esout_file.close()


###### get the file names ######
Input,output,qin,qout,punch,espot,qwts,esout = file_in()

input_file = open(Input, 'r')
output_file = open(output, 'w')

###### read the atomic centers and q0's, then read the potential inf ######
read_in()

###### if its a multiple molecule run (nmol>0), do mult. mol. input ######
if nmol > 1:
	mult_mol()

###### center & reorient molecule once in preparation for dipole & quadrupole ######
#if nmol == 1:
#	reornt()   ### NOT USED

###### read in the qm esp, forming the matrices apot(awt) and bpot(bwt)
matpot()

###### process the input (freezing, equivalencing charges)
data_prep()

# set up cycle control structure: if icycle .gt. 0, come back to stmt 10
# subroutine cycle and read new "qwt"

# icycle is initially set in readin
# subsequently decremented in cycle
icycle = -1

qcal = np.zeros((iuniq))

while icycle != 0:
	# call cycle to see if we are supposed to cycle. reset icycle as needed
	icycle = cycle(icycle)
	
	if irstrnt == 2:
		# if irstrnt= 2 then we just want to compare esp's to q0's
		for k in range(iuniq):
			qcal[k] = q0[k]
		qwt = 0.0
	else:
		# do the charge fitting
		charge_opt()
	
	# calculate residuals sum-of-squares (chi-square) for the esp's
	evlchi()
	
	# now calculate and print dipole and quadrupole moments
	if nmol == 1:
		elec_mom()
	
	# now punch short summary of charges & important evaluation criteria
	pun_sum()

# now calculate and print sum-of-squares, sigma, and rms
pot_out()

if ioutopt == 1:
	# write the coords, old esps, new esps and residuals
	wrt_pot()

input_file.close()
output_file.close()
punch_file.close()


