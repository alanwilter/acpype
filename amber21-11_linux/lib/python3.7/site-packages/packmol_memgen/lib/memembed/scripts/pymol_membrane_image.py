#!/usr/bin/python

import __main__
__main__.pymol_argv = [ 'pymol', '-qc'] # Quiet and no GUI

import pymol
import sys, time, os

pymol.finish_launching()
InStructurePath = os.path.abspath(sys.argv[1])
OutPngPath = os.path.abspath(sys.argv[2])
pymol.cmd.set("max_threads",2)
pymol.cmd.load(InStructurePath)
view = pymol.cmd.get_view()
camera = list(view)
camera[0] = -1
camera[1] = 0
camera[2] = 0
camera[3] = 0
camera[4] = 0
camera[5] = 1
camera[6] = 0
camera[7] = 1
camera[8] = 0
pymol.cmd.set_view(camera)
pymol.cmd.hide("lines","all")
pymol.cmd.show("cartoon","all")
pymol.cmd.ray(1280,720)
pymol.cmd.png(OutPngPath, 300, 0 )
pymol.cmd.quit()

