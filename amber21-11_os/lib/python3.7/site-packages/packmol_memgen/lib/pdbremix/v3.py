# encoding: utf-8

__doc__ = """

Interface 3D vector geometry library.

This is an interface to choose between the pure Python version
or the numpy version. This depends on numpy availability.
"""

try:
  import numpy
  from .v3numpy import *
except:
  from .v3array import *

