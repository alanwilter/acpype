# distutils: language = c++
from libcpp.vector cimport vector
from libcpp.string cimport string
from ..frame cimport _Frame, Frame
from ...core.c_core cimport _AtomMask, AtomMask, _ArgList, ArgList
from ...core.coordinfo cimport _CoordinateInfo, CoordinateInfo
from ...core.c_core cimport AtomMask
from ...core.box cimport Box
from ...topology.topology cimport _Topology, Topology
from ...datasets.c_datasets cimport _DatasetCoords
from ...datasets.c_datasetlist cimport DatasetList as CpptrajDatasetList, _DatasetList as _CpptrajDatasetList
from ...analysis.c_action.actionlist cimport _ActionList, ActionList


cdef extern from "DataSet_Coords_TRJ.h": 
    cdef cppclass _TrajectoryCpptraj "DataSet_Coords_TRJ" (_DatasetCoords):
        _TrajectoryCpptraj() 
        int AddSingleTrajin(const string&, _ArgList&, _Topology *)
        #size_t Size() const 
        void GetFrame(int idx, _Frame& fIn)
        #void GetFrame(int idx, _Frame& fIn, _AtomMask& mIn)
        #void CoordsSetup(const _Topology&, const CoordinateInfo &)
        #const _Topology& Top() const 
        #const CoordinateInfo& CoordsInfo()
        #_Frame AllocateFrame()
        _Topology* TopPtr()

cdef class TrajectoryCpptraj:
    cdef Topology _top
    cdef _TrajectoryCpptraj* thisptr
    cdef object tmpfarray
    cdef list _filelist
    cdef public _base
    cdef public bint _own_memory
    cdef public ActionList _actionlist 
    cdef public CpptrajDatasetList _cdslist
    cdef public bint _being_transformed
    cdef public bint _being_superposed
    cdef public list _transform_commands
    cdef public int _max_count_to_reset
