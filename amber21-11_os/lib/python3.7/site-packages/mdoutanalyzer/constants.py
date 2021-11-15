"""
This module holds some constant data to aid in the operation of MdoutAnalyzer.py
"""
import re

class LabelDesc(dict):
   """
   A dict of label descriptions, but suppressing KeyError exceptions
   """
   
   scfre = re.compile(r'([A-Z0-9]*)(?:-)?ESCF')

   def __getitem__(self, key):
      """
      Suppress key errors, and pretend the key name is the best alternative
      description
      """
      try:
         return dict.__getitem__(self, key)
      except KeyError:
         # Special-case the SCF energies
         if self.scfre.match(key):
            return '%s SCF energy' % self.scfre.match(key).groups()[0]
         return key

LABEL_DESC = LabelDesc(
       {
        'NSTEP' : 'MD or minimization Step number',
        'TIME(PS)' : 'Time in MD trajectory (picoseconds)',
        'TEMP(K)' : '\'Instantaneous\' Temperature (Kelvin)',
        'Etot' : 'Total energy',
        'EKtot' : 'Total kinetic energy',
        'EPtot' : 'Total potential energy',
        'BOND' : 'Bond energy',
        'ANGLE' : 'Angle energy',
        'DIHED' : 'Torsion energy',
        'UB' : 'Urey-Bradley energy (CHARMM FF)',
        'IMP' : 'Improper torsion energy (CHARMM FF)',
        'CMAP' : 'Correction map (coupled-torsion) energy',
        '1-4 NB' : '1-4 van der Waals energy',
        '1-4 EEL' : '1-4 electrostatic energy',
        'EELEC' : 'Electrostatic energy',
        'VDWAALS' : 'van der Waals Energy',
        'EHBOND' : '10-12 potential energy (H-bond term)',
        'RESTRAINT' : 'Restraint potential energy',
        'EGB' : 'Generalized Born polar solvation energy',
        'EPB' : 'Poisson-Boltzmann polar solvation energy',
        'ERISM' : '3D-RISM solvation energy',
        'TEMP0' : 'Target (thermostat) temperature',
        'SGFT' : 'Self-guided Langevin Guiding Factor',
        'TEMPSG' : 'Self-guided Langevin temperature',
        'STAGE' : 'Self-guided Langevin Stage ID',
        'EMAP' : 'Energy map restraint energy',
        'REPNUM' : 'Replica number in REMD runs',
        'EXCHANGE#' : 'The REMD exchange number',
        'SOLVPH' : 'Solution pH',
        'REMD_DIMENSION' : 'multi-dimensional REMD dimension',
        'SURFTEN' : 'Non-polar solvation free energy',
        'EKCMT' : 'Kinetic energy of system COM translation',
        'ECAVITY' : 'Cavitation non-polar solvation free energy',
        'EDISPER' : 'Attractive dispersion non-polar solvation free energy',
        'VIRIAL' : 'System Virial',
        'DV/DL' : 'Derivative of potential with respect to lambda',
        'PUPESCF' : 'PUPIL QM SCF energy',
        'ENERGY' : 'Total potential energy', # minimization
        'RMS' : 'Minimization Root-mean-squared deviation',
        'GMAX' : 'Maximum minimization gradient',
        'VOLUME' : 'System volume (cubic Angstroms)',
        'PRESS' : 'Instantaneous system pressure (bar)'
       }
                      )
