"""
Exceptions which are used in the pyMSMT package
"""

from sys import stderr

class pymsmtError(Exception):
    def __init__(self, info='pyMSMT Error'):
        self.info = info
    def __str__(self):
        return self.info

class pymsmtWarning(Warning, pymsmtError):
    def __init__(self, info='pyMSMT Warning'):
        self.info = info
    def __str__(self):
        return self.info

class BadResidueNumError(pymsmtError):
    """Bad Resdiue Number"""
    pass
