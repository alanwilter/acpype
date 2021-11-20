""""""
from __future__ import absolute_import
from glob import glob
from . import c_analysis

analysislist = []
for act in c_analysis.__dict__.keys():
    if 'Analysis' in act:
        analysislist.append(act)

__all__ = analysislist
