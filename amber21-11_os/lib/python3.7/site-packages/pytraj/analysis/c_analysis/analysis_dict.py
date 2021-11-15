from __future__ import print_function
from pytraj.analysis.c_analysis import c_analysis as allanalyses

ADICT = {}

for key in allanalyses.__dict__.keys():
    if "Analysis_" in key:
        act = key.split('Analysis_')[1]
        # create dict of analysis objects
        ADICT[act.lower()] = allanalyses.__dict__[key]

# make another dict to convert something like 'MolSurf' to 'molsurf'


class AnalysisDict:
    def __init__(self):
        self.analysisdict = ADICT

    def __getitem__(self, key):
        # return Action object
        return self.analysisdict[key]()

    def keys(self):
        return self.analysisdict.keys()
