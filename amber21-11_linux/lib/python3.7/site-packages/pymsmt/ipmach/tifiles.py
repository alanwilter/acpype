# This module for Generating the TI files

import os

Mass = {'H':   1.008,  'C':  12.01,  'N':  14.01,  'O':  16.00,  'S':  32.06,
        'P':   30.97,  'F':  19.00, 'Cl':  35.45, 'Br':  79.90,  'I':  126.9,
        'Li':   6.94, 'Na':  22.99,  'K':  39.10, 'Rb':  85.47, 'Cs': 132.91,
        'Be':   9.01, 'Cu':  63.55, 'Ni':  58.69, 'Pt': 195.08, 'Zn':   65.4,
        'Co':  58.93, 'Pd': 106.42, 'Ag': 107.87, 'Cr':  52.00, 'Fe':  55.85,
        'Mg': 24.305,  'V':  50.94, 'Mn':  54.94, 'Hg': 200.59, 'Cd': 112.41,
        'Yb': 173.05, 'Ca':  40.08, 'Sn': 118.71, 'Pb':  207.2, 'Eu': 151.96,
        'Sr':  87.62, 'Sm': 150.36, 'Ba': 137.33, 'Ra': 226.03, 'Al':  26.98,
	    'Fe':  55.85, 'Cr':  52.00, 'In': 114.82, 'Tl': 204.38, 'Y' :  88.91,
	    'La': 138.91, 'Ce': 140.12, 'Pr': 140.91, 'Nd': 144.24, 'Sm': 150.36,
	    'Eu': 151.96, 'Gd': 157.25, 'Tb': 158.93, 'Dy':  162.5, 'Er': 167.26,
	    'Tm': 168.93, 'Lu': 174.97, 'As':  74.92, 'Ru': 101.07, 'Hf': 178.49,
	    'Zr':  91.22, 'Ce': 140.12, 'U' : 238.03, 'Pu': 244.06, 'Th': 232.04,
	    'Mo':  95.96
        }

def write_frcmod(ion, fname):
    mass = Mass[ion.element]
    frcmodf = open(fname, 'w')
    print("# remark goes here", file=frcmodf)
    print("MASS", file=frcmodf)
    print("%s  %5.2f" %(ion.attype, mass), file=frcmodf)
    print(" ", file=frcmodf)
    print("BOND", file=frcmodf)
    print(" ", file=frcmodf)
    print("ANGL", file=frcmodf)
    print(" ", file=frcmodf)
    print("DIHE", file=frcmodf)
    print(" ", file=frcmodf)
    print("NONB", file=frcmodf)
    print("%s    %5.3f    %10.8f" %(ion.attype, ion.rmin, ion.ep), file=frcmodf)
    print(" ", file=frcmodf)
    frcmodf.close()

def write_cmd(ion):
    cmdf = open(ion.resname+'.cmd', 'w')
    print("i = createAtom %s %s %3.1f" %(ion.atname, ion.attype, ion.charge), file=cmdf)
    print("set i element \"%s\"" %ion.element, file=cmdf)
    print("set i position { 0 0 0 }", file=cmdf)
    print("r = createResidue %s" %ion.resname, file=cmdf)
    print("add r i", file=cmdf)
    print("%s = createUnit %s" %(ion.resname, ion.resname), file=cmdf)
    print("add %s r" %ion.resname, file=cmdf)
    print("saveOff %s ./%s.lib" %(ion.resname, ion.resname), file=cmdf)
    print("quit", file=cmdf)
    cmdf.close()

def write_leapin(ion0, ion1, watermodel, distance):

    leapf = open('tleap.in', 'w')

    print("source leaprc.protein.ff14SB", file=leapf)
    print("loadOff solvents.lib", file=leapf)

    print("loadOff %s.lib" %ion0.resname, file=leapf)
    print("loadOff %s.lib" %ion1.resname, file=leapf)

    # 1st is the metal ion without VDW or CHG, 2nd is the one with VDW and CHG
    print("ion = loadpdb ION.pdb", file=leapf)

    print("loadamberparams %s.frcmod" %ion0.attype, file=leapf)  # VDW=0
    print("loadamberparams %s.frcmod" %ion1.attype, file=leapf)  # VDW!=0

    if watermodel == 'OPC':
        print("loadamberparams frcmod.opc", file=leapf)
        print("solvatebox ion OPCBOX %5.1f" %distance, file=leapf)
    elif watermodel == 'OPC3':
        print("loadamberparams frcmod.opc3", file=leapf)
        print("solvatebox ion OPC3BOX %5.1f" %distance, file=leapf)
    elif watermodel == 'FB3':
        print("loadamberparams frcmod.tip3pfb", file=leapf)
        print("solvatebox ion FB3BOX %5.1f" %distance, file=leapf)
    elif watermodel == 'FB4':
        print("loadamberparams frcmod.tip4pfb", file=leapf)
        print("solvatebox ion FB4BOX %5.1f" %distance, file=leapf)
    elif watermodel == 'SPC':
        print("solvatebox ion SPCBOX %5.1f" %distance, file=leapf)
    elif watermodel == 'SPCE':
        print("loadamberparams frcmod.spce", file=leapf)
        print("solvatebox ion SPCBOX %5.1f" %distance, file=leapf)
    elif watermodel == 'TIP3P':
        print("solvatebox ion TIP3PBOX %5.1f" %distance, file=leapf)
    elif watermodel == 'TIP4P':
        print("loadamberparams frcmod.tip4p", file=leapf)
        print("solvatebox ion TIP4PBOX %5.1f" %distance, file=leapf)
    elif watermodel == 'TIP4PEW':
        print("loadamberparams frcmod.tip4pew", file=leapf)
        print("solvatebox ion TIP4PEWBOX %5.1f" %distance, file=leapf)
    elif watermodel == 'TIP5P':
        print("loadamberparams frcmod.tip5p", file=leapf)
        print("solvatebox ion TIP5PBOX %5.1f" %distance, file=leapf)

    print("ion2 = copy ion", file=leapf)

    # Only keep the ion without VDW or CHG
    print("remove ion2 ion2.2", file=leapf)

    # For sander, prmtop without VDW or CHG
    print("saveamberparm ion2 %s_wat_s0.prmtop %s_wat_s0.inpcrd" %(ion0.element, ion0.element), file=leapf)

    # For TI one-step method using pmemd, from the ion 0 to ion 1
    print("saveamberparm ion %s_wat_pti.prmtop %s_wat_pti.inpcrd" %(ion0.element, ion0.element), file=leapf)

    # VDW != 0
    print("loadamberparams %s_vdw.frcmod" %ion0.attype, file=leapf)

    # For sander, prmtop with VDW but not CHG
    print("saveamberparm ion2 %s_wat_sv.prmtop %s_wat_sv.inpcrd" %(ion0.element, ion0.element), file=leapf)

    # For TI two-step method using pmemd: from ion 0 (with VDW) to ion 1
    print("saveamberparm ion %s_wat_pvdw.prmtop %s_wat_pvdw.inpcrd" %(ion0.element, ion0.element), file=leapf)

    # For sander/MD, prmtop with VDW or CHG
    print("remove ion ion.1", file=leapf)
    print("saveamberparm ion %s_wat_svc.prmtop %s_wat_svc.inpcrd" %(ion0.element, ion0.element), file=leapf)

    print("quit", file=leapf)
    leapf.close()

def addc4(ion0, ion1, c4v):

    # Only for sander
    print("Add C4 parameters...")

    #Without charge and VDW
    c4f0 = open('c4_0.txt', 'w')
    print("%s %f" %(ion0.element+str(int(ion0.charge)), 0.0), file=c4f0)
    c4f0.close()

    #Without VDW or charge
    parmf = open('addc4_0.in', 'w')
    print("loadRestrt %s_wat_s0.inpcrd" %ion0.element, file=parmf)
    print("setOverwrite True", file=parmf)
    print("add12_6_4 :1@%s c4file c4_0.txt" %ion0.atname, file=parmf)
    print("outparm %s_wat_s0.prmtop %s_wat_s0.inpcrd" %(ion0.element, ion0.element), file=parmf)
    parmf.close()

    os.system("parmed -i addc4_0.in -p %s_wat_s0.prmtop | tail -n 3" %ion0.element)

    #With VDW but not charge
    parmf = open('addc4_v.in', 'w')
    print("loadRestrt %s_wat_sv.inpcrd" %ion0.element, file=parmf)
    print("setOverwrite True", file=parmf)
    print("add12_6_4 :1@%s c4file c4_0.txt" %(ion0.atname), file=parmf)
    print("outparm %s_wat_sv.prmtop %s_wat_sv.inpcrd" %(ion0.element, ion0.element), file=parmf)
    parmf.close()
    os.system("parmed -i addc4_v.in -p %s_wat_sv.prmtop | tail -n 3" %ion0.element)

    #With Charge and VDW
    c4f1 = open('c4_vc.txt', 'w')
    print("%s %f" %(ion1.element+str(int(ion1.charge)), c4v), file=c4f1)
    c4f1.close()

    parmf = open('addc4_vc.in', 'w')
    print("loadRestrt %s_wat_svc.inpcrd" %ion1.element, file=parmf)
    print("setOverwrite True", file=parmf)
    print("add12_6_4 :1@%s c4file c4_vc.txt" %(ion1.atname), file=parmf)
    print("outparm %s_wat_svc.prmtop %s_wat_svc.inpcrd" %(ion1.element, ion1.element), file=parmf)
    parmf.close()

    os.system("parmed -i addc4_vc.in -p %s_wat_svc.prmtop | tail -n 3" %ion1.element)

def gene_topcrd(ion0, ion1, watermodel, distance, ifc4=0, c4v=0.0):

    print("Generate the Topology and Coordinate Files...")

    #lib file
    write_cmd(ion0)
    write_cmd(ion1)
    os.system("tleap -s -f %s.cmd > %s.out" %(ion0.resname, ion0.resname))
    os.system("tleap -s -f %s.cmd > %s.out" %(ion1.resname, ion1.resname))

    #pdb file
    pdbf = open('ION.pdb', 'w')
    print("HETATM 2032 %2s    %2s A   1       0.000   0.000   0.000  1.00  0.00" %(ion0.atname, ion0.resname), file=pdbf)
    print("TER", file=pdbf)
    print("HETATM 2032 %2s    %2s A   1       0.000   0.000   0.000  1.00  0.00" %(ion1.atname, ion1.resname), file=pdbf)
    pdbf.close()

    #frcmod file
    write_frcmod(ion0, ion0.attype + '.frcmod')
    write_frcmod(ion1, ion1.attype + '.frcmod')

    ion0.rmin = ion1.rmin
    ion0.ep = ion1.ep

    write_frcmod(ion0, ion0.attype + '_vdw.frcmod')

    #tleap file
    write_leapin(ion0, ion1, watermodel, distance)
    os.system("tleap -s -f tleap.in > tleap.out")

    #add c4 parameters
    if ifc4 == 1:
        addc4(ion0, ion1, c4v)


