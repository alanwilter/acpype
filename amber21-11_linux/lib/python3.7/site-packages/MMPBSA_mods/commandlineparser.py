"""
This module contains classes and such that are responsible for parsing
command-line arguments for MMPBSA.py.  All of the files specified for use
in MMPBSA.py will be assigned as attributes to the returned class.

                          GPL LICENSE INFO                             

  Copyright (C) 2009  Dwight McGee, Billy Miller III, and Jason Swails

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
   
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place - Suite 330,
  Boston, MA 02111-1307, USA.
"""

import os
from MMPBSA_mods import __version__

class OptionList(object):
   """
   Just a container to hold the command-line options. Necessary when reading in
   a MMPBSA.py info file to have a container to load the results from the
   parser.
   """
   pass

# Set up the MM/PBSA parser here. It clutters up the MMPBSA_App to do it there
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
parser = ArgumentParser(epilog='''This program will calculate binding free
                        energies using end-state free energy methods on an
                        ensemble of snapshots using a variety of implicit
                        solvent models''',
                        formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('-v', '--version', action='version',
                    version='%%(prog)s %s' % __version__)
parser.add_argument('--input-file-help', dest='infilehelp', action='store_true',
                    help='Print all available options in the input file.',
                    default=False)
group = parser.add_argument_group('Miscellaneous Options')
group.add_argument('-O', '--overwrite', default=False, action='store_true',
                  help='Allow output files to be overwritten', dest='overwrite')
group.add_argument('-prefix', dest='prefix', default='_MMPBSA_',
                  metavar='<file prefix>',
                  help='Prefix for intermediate files.')
group = parser.add_argument_group('Input and Output Files', '''These options
                        specify the input files and optional output files.''')
group.add_argument('-i', dest='input_file', metavar='FILE',
                   help='MM/PBSA input file.')
group.add_argument('-xvvfile', dest='xvvfile', help='XVV file for 3D-RISM.',
                  default=os.path.join(os.getenv('AMBERHOME'), 'dat', 'mmpbsa',
                                       'spc.xvv'))
group.add_argument('-o', dest='output_file', default='FINAL_RESULTS_MMPBSA.dat',
                  metavar='FILE', help='Output file with MM/PBSA statistics.')
group.add_argument('-do', dest='decompout', metavar='FILE',
                   default='FINAL_DECOMP_MMPBSA.dat',
                   help='Output file for decomposition statistics summary.')
group.add_argument('-eo', dest='energyout', metavar='FILE',
                  help='''CSV-format output of all energy terms for every frame
                  in every calculation. File name forced to end in [.csv].
                  This file is only written when specified on the
                  command-line.''')
group.add_argument('-deo', dest='dec_energies', metavar='FILE',
                  help='''CSV-format output of all energy terms for each printed
                  residue in decomposition calculations. File name forced to end
                  in [.csv]. This file is only written when specified on the
                  command-line.''')
group.add_argument_group('Input Topology Files', '''These topology files must be
                      consistent with each other. ante-MMPBSA.py can be used to
                      generate consistent topology files if necessary.''')
group.add_argument('-sp', dest='solvated_prmtop', metavar='<Topology File>',
                  help='''Topology file of a fully solvated system. If provided,
                  the atoms specified by <strip_mask> will be stripped from the
                  trajectory file. The complex topology file (-cp) must be
                  consistent with this stripped trajectory''')
group.add_argument('-cp', dest='complex_prmtop', metavar='<Topology File>',
                  default='complex_prmtop', help='''Topology file of the
                  bound complex (or the single system for 'stability'
                  calculations)''')
group.add_argument('-rp', dest='receptor_prmtop', metavar='<Topology File>',
                  help='''Topology file of the unbound receptor. If omitted (and
                  -lp is omitted, too), a stability calculation with just the
                  complex will be performed.''')
group.add_argument('-lp', dest='ligand_prmtop', metavar='<Topology File>',
                  help='''Topology file of the unbound ligand. If omitted (and
                  -rp is omitted, too), a stability calculation with just the
                  complex will be performed.''')
group.add_argument('-mc', dest='mutant_complex_prmtop',
                  metavar='<Topology File>', help='''Complex topology file of
                  the mutant complex in which one residue has been mutated to
                  either a glycine or alanine to perform computational alanine
                  (or glycine) scanning.''')
group.add_argument('-mr', dest='mutant_receptor_prmtop',
                  metavar='<Topology File>', help='''Receptor topology file of
                  the mutant receptor (see -mc above). If omitted, the mutation
                  is assumed to be in the ligand.''')
group.add_argument('-ml', dest='mutant_ligand_prmtop',
                  metavar='<Topology File>', help='''Ligand topology file of the
                  mutant receptor (see -mc above). If omitted, the mutation is
                  assumed to be in the receptor.''')
group.add_argument('-srp', dest='solvated_receptor_prmtop',
                  metavar='<Topology File>', help='''Receptor ligand topology
                  file. For use with multiple-trajectory simulations when the
                  receptor trajectory is solvated. This will trigger the atoms
                  specified by 'strip_mask' to be removed from the receptor
                  trajectory''')
group.add_argument('-slp', dest='solvated_ligand_prmtop',
                  metavar='<Topology File>', help='''Solvated ligand topology
                  file. See -srp description above.''')
group = parser.add_argument_group('Input Trajectory Files', '''These files
                  contain the snapshots analyzed by MM/PBSA-type
                  calculations.''')
group.add_argument('-y', dest='mdcrd', nargs='*', default=['mdcrd'],
                  help='''Input trajectories of the (maybe solvated) complex.
                  (specify as many as you'd like).''', metavar='MDCRD')
group.add_argument('-yr', dest='receptor_mdcrd', help='''Receptor trajectory
                  file for multiple trajectory approach''', nargs='*',
                  metavar='MDCRD')
group.add_argument('-yl', dest='ligand_mdcrd', metavar='MDCRD', nargs='*',
                  help='''Ligand trajectory file for multiple trajectory
                  approach.''')
group = parser.add_argument_group('Miscellaneous Actions')
group.add_argument('-make-mdins', dest='make_mdins', default=False,
                  action='store_true', help='''Create the input files for each
                  calculation and quit. This allows you to modify them and
                  re-run using -use-mdins''')
group.add_argument('-use-mdins', dest='use_mdins', default=False,
                  action='store_true', help='''Use existing input files for each
                  calculation. If they do not exist with the appropriate names,
                  %(prog)s will quit in error.''')
group.add_argument('-rewrite-output', dest='rewrite_output', default=False,
                  action='store_true', help='''Do not re-run any calculations,
                  just parse the output files from the previous calculation and
                  rewrite the output files.''')
group.add_argument('--clean', dest='clean', action='store_true', default=False,
                  help='''Clean temporary files and quit.''')
