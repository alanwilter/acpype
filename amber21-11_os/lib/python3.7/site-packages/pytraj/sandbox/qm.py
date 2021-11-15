def read_orca_trj(fname):
    """return numpy 2D array
    """
    # http://stackoverflow.com/questions/14645789/
    # numpy-reading-file-with-filtering-lines-on-the-fly
    import numpy as np
    regexp = r'\s+\w+' + r'\s+([-.0-9]+)' * 3 + r'\s*\n'
    return np.fromregex(fname, regexp, dtype='f')
