# distutil: language = c++

from libcpp.string cimport string
from ...core.c_core cimport _DispatchObject, DispatchObject, DispatchAllocatorType, FunctPtr
from ...core.box cimport Box
from ...core.c_core cimport _ArgList, ArgList, _AtomMask, AtomMask
from ...datafiles.datafiles cimport _DataFileList, DataFileList
from ...topology.topology cimport _Topology, Topology
from ...trajectory.frame cimport _Frame, Frame
from ...datasets.c_datasetlist cimport _DatasetList, DatasetList
from .c_action cimport _Action, Action, _ActionInit, _ActionSetup, _ActionFrame, CoordinateInfo

cdef extern from "ActionList.h":
    cdef cppclass _ActionList "ActionList":
        _ActionList()
        void Clear()
        void SetDebug(int)
        int Debug()
        int AddAction(_Action*, _ArgList&,
                      _ActionInit&,)
        int SetupActions(_ActionSetup, bint exit_on_error)
        bint DoActions(int, _ActionFrame)
        void PrintActions()
        void List()
        bint Empty()
        int Naction()
        const string& CmdString(int)
        DispatchAllocatorType ActionAlloc(int i)

cdef class ActionList:
    cdef _ActionList* thisptr

    # alias for TopologyList (self.process(top))
    cdef object top

    # check if self.process is already called or not
    cdef public bint is_setup
    cdef public object _dslist
    cdef public object _dflist
    cdef public object _crdinfo
    cdef public unsigned int n_frames
