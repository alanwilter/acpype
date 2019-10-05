#!/opt/anaconda1anaconda2anaconda3/bin/python

from tkFileDialog import askopenfilenames
from Tkinter import Tk, BOTH
import matplotlib
matplotlib.use('TkAgg')
try:
    from mdoutanalyzer import __version__, __author__, __date__
    from mdoutanalyzer.graphproperties import GraphControlWindow
    from mdoutanalyzer.mdout import AmberMdout, OpenMMMdout
    from mdoutanalyzer.toplevel_app import MdoutAnalyzerApp
except ImportError:
    import os
    amberhome = os.getenv('AMBERHOME') or '$AMBERHOME'
    raise ImportError('Could not import Amber Python modules. Please make sure '
                      'you have sourced %s/amber.sh (if you are using sh/ksh/'
                      'bash/zsh) or %s/amber.csh (if you are using csh/tcsh)' %
                      (amberhome, amberhome))
from optparse import OptionParser
import re
import sys

geore = re.compile(r'(\d+)x(\d+)\+(-?\d+)\+(-?\d+)')
geoformat = '%dx%d+%d+%d'

verstring = """
   %%prog : An AMBER MD output file parser and graphing utility

                              Version %s
                             %s

   Written by %s
""" % (__version__, __date__, __author__)

parser = OptionParser(usage='%prog [Options] [mdout1] [mdout2] ... [mdoutN]',
                      version=verstring)
parser.add_option('--openmm', dest='openmm', default=False,
      action='store_true', help='Treat the input files as comma-delimited '
      'output from the Amber/OpenMM-interface in ParmEd.')

opt, arg = parser.parse_args()

if opt.openmm:
    MdoutClass = OpenMMMdout
else:
    MdoutClass = AmberMdout

root = Tk()
root.title('Mdout Analyzer')
if not arg:
    arg = askopenfilenames(title='Select Mdout File(s)', parent=root,
                           filetypes=[('Mdout Files', '*.mdout'),
                                      ('All Files', '*')])

if not arg:
    print ('No mdout files chosen!')
    sys.exit(1)

for f in arg:
    try:
        mdout += MdoutClass(f)
    except NameError:
        mdout = MdoutClass(f)

app = MdoutAnalyzerApp(root, mdout)
app.pack(fill=BOTH)
# Update idle tasks here to make sure the whole app is filled in before making
# our window non-resizable.  In some instances, this can chop off the second
# frame
app.update_idletasks()
# Now make our window non-resizable
root.resizable(False, False)
root.update_idletasks()
# Load up the graph options window and move it to the right of the root window
rootgeo = [int(i) for i in geore.match(root.geometry()).groups()]
graphprops = GraphControlWindow(root, app.graph_props)
root.update_idletasks()
ggeo = [int(i) for i in geore.match(graphprops.geometry()).groups()]
graphprops.geometry(geoformat % (ggeo[0], ggeo[1], rootgeo[2] + rootgeo[0] + 5,
                                 rootgeo[3])
                   )

# For some reason it seems like mainloop is not exited properly if we simply
# call root.destroy(). Therefore, we'll replace WM_DELTE_WINDOW with root.quit()
# to bust out of the mainloop.
root.protocol('WM_DELETE_WINDOW', root.quit)

# Enter our mainloop
root.mainloop()
