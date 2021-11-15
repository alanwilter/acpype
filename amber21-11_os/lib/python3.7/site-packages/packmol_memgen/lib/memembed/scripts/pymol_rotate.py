#!/usr/bin/python

import __main__
__main__.pymol_argv = [ 'pymol', '-qc'] # Quiet and no GUI

import pymol
import sys, time, os, random

pymol.finish_launching()
InStructurePath = os.path.abspath(sys.argv[1])
OutStructurePath = os.path.abspath(sys.argv[2])
x_rot = random.uniform(0,360)
y_rot = random.uniform(0,360)
z_rot = random.uniform(0,360)
z_trans = random.uniform(-20,20)
pymol.cmd.load(InStructurePath)
pymol.cmd.rotate("x", x_rot)
pymol.cmd.rotate("y", y_rot)
pymol.cmd.rotate("z", z_rot)
pymol.cmd.translate([0,0,z_trans])
pymol.cmd.save(OutStructurePath)
pymol.cmd.quit()

