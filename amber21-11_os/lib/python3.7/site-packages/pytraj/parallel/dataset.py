from __future__ import absolute_import
import numpy as np
from collections import OrderedDict
from pytraj import matrix
from pytraj import mean_structure
from pytraj import volmap
from pytraj import Frame
from pytraj import ired_vector_and_matrix
from pytraj import rotation_matrix
from pytraj.utils.tools import concat_dict

from .base import concat_hbond


class PmapDataset(object):
    '''Dataset handlder for both pmap and pmap_mpi. This class is for internal use.

    Parameters
    ----------
    data_collection : List[Tuple[OrdereDict[key, array], n_frames]] or numpy
    array based on returning type
    func : callable
        a function that produces data_collection. We use this to determine data type
    traj : Trajectory-like, original input
    kwargs : Dict
        original kwargs parameters, used to determine return type
    '''

    def __init__(self, data_collection, func=None, traj=None, kwargs=None):
        self.data = data_collection
        self.func = func
        self.traj = traj
        self.kwargs = kwargs

    def process(self):
        # val : Tuple[OrdereDict, n_frames]

        if self.func in [matrix.dist, matrix.idea, volmap]:
            mat = np.sum((val[0] * val[1]
                          for val in self.data)) / self.traj.n_frames
            return mat
        elif self.func in [
                ired_vector_and_matrix,
        ]:
            # val : Tuple[(vecs, mat), n_frames]
            mat = np.sum((val[0][1] * val[1]
                          for val in self.data)) / self.traj.n_frames
            vecs = np.column_stack(val[0][0] for val in self.data)
            return (vecs, mat)
        elif self.func in [
                rotation_matrix,
        ]:
            if 'with_rmsd' in self.kwargs.keys() and self.kwargs['with_rmsd']:
                # val : Tuple[(mat, rmsd), n_frames]
                mat = np.row_stack(val[0][0] for val in self.data)
                rmsd_ = np.hstack(val[0][1] for val in self.data)
                return OrderedDict(out=(mat, rmsd_))
            else:
                # val : Tuple[mat, n_frames]
                mat = np.row_stack(val[0] for val in self.data)
                return OrderedDict(mat=mat)
        elif self.func == mean_structure:
            xyz = np.sum((x[1] * x[0].xyz
                          for x in self.data)) / self.traj.n_frames
            frame = Frame(xyz.shape[0])
            frame.xyz[:] = xyz
            return frame
        elif 'hbond' in self.func.__name__:
            return concat_hbond(self.data)
        else:
            return concat_dict((x[0] for x in self.data))
