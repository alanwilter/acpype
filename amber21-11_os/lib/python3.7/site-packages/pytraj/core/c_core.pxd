# distutil: language = c++

from libcpp.string cimport string
from posix.unistd cimport off_t
from ..cython_extra_header.cpp_vector cimport vector

from .box cimport _Box, Box
from ..datafiles.datafiles cimport _DataFileList, DataFileList, _DataFile, DataFile
from ..datasets.c_datasetlist cimport _DatasetList, DatasetList

ctypedef _BaseIOtype* (*AllocatorType)()
ctypedef void (*HelpType)()

cdef extern from "CpptrajState.h": 
    ctypedef enum RetType "CpptrajState::RetType":
        pass
    cdef cppclass _CpptrajState "CpptrajState":
        _CpptrajState()
        _DatasetList& DSL()
        _DataFileList& DFL()
        int AddTrajin(_ArgList &, bint)
        int AddTrajin(const string&)
        int RunAnalyses()
        inline int AddTrajout "AddOutputTrajectory" (const _ArgList&)
        inline int AddTrajout "AddOutputTrajectory" (const string&)
        int AddReference(const string&, _ArgList &)
        inline int AddReference(const string&)
        inline int AddAction(DispatchAllocatorType, _ArgList &)
        inline int AddAnalysis(DispatchAllocatorType, _ArgList &)
        int TrajLength(const string&, const vector[string]&)
        int Run()
        void MasterDataFileWrite()
        bint EmptyState()

cdef class CpptrajState:
    cdef _CpptrajState* thisptr
    cdef public DataFileList datafilelist
    cdef public DatasetList _datasetlist

cdef extern from "Command.h": 
    cdef cppclass _Command "Command":
        @staticmethod
        void Init()
        @staticmethod
        void Free()
        @staticmethod
        RetType ProcessInput(_CpptrajState&, const string&)
        @staticmethod
        RetType Dispatch(_CpptrajState&, const string&)

cdef extern from "BaseIOtype.h":
    #ctypedef _BaseIOtype* (*AllocatorType)()
    #ctypedef void (*HelpType)()
    cdef cppclass _BaseIOtype "BaseIOtype":
        pass

cdef class BaseIOtype:
    cdef _BaseIOtype* baseptr0

ctypedef _DispatchObject* (*DispatchAllocatorType)()
cdef extern from "DispatchObject.h":
    cdef cppclass _DispatchObject "DispatchObject":
        pass

cdef class DispatchObject:
    cdef _DispatchObject* thisptr

# dummy class to hold function pointer
cdef class FunctPtr:
    cdef DispatchAllocatorType ptr
    # used for BaseIOtype
    cdef AllocatorType allocptr


cdef extern from "AtomMask.h": 
    cdef cppclass _AtomMask "AtomMask":
        _AtomMask()
        _AtomMask(const string&)
        _AtomMask(int, int)
        _AtomMask(int)
        _AtomMask(vector[int], int)
        _AtomMask(const _AtomMask &)
        #_AtomMask & operator =(const _AtomMask &)
        const vector [int]& Selected()const 
        vector[int].const_iterator begin()const 
        vector[int].const_iterator end()const 
        int back()const 
        int Nselected()const 
        const int & index_opr "operator[]"(int idx)const 
        const char * MaskString()const 
        const string& MaskExpression()const 
        bint MaskStringSet()const 
        bint None()const 
        bint IsCharMask()const 
        void ResetMask()
        void ClearSelected()
        void InvertMask() except +
        int NumAtomsInCommon(const _AtomMask&)
        void AddSelectedAtom(int i)
        void AddAtom(int)
        void AddAtoms(const vector [int]&)
        void AddAtomRange(int, int)
        void AddMaskAtPosition(const _AtomMask&, int)
        void PrintMaskAtoms(const char *)const 
        #int SetMaskString(const char *)
        int SetMaskString(const string&)
        void SetupIntMask(const char *, int, int)
        void SetupCharMask(const char *, int, int)
        bint AtomInCharMask(int)const 
        bint AtomsInCharMask(int, int)const 
        void SetNatom(int a)
        int ConvertToCharMask()
        int ConvertToIntMask()
        void MaskInfo()const 
        void BriefMaskInfo()const 
        #inline token_iterator begintoken()const 
        #inline token_iterator endtoken()const 

#ctypedef fused charstring:
#    char*
#    string

cdef class AtomMask:
    cdef _AtomMask* thisptr


cdef extern from "FileName.h":
    cdef cppclass _FileName "FileName":
        _FileName()
        _FileName(_FileName)
        int SetFileName(string)
        int SetFileNameWithExpansion(string)
        int SetFileName(string, bool)
        void clear()
        bint MatchFullOrBase(string)
        string Full()
        string Base()
        char * full()
        char * base()
        string Ext()
        string Compress()
        string DirPrefix()
        bint empty()

cdef class FileName:
    cdef _FileName* thisptr

cdef extern from "NameType.h":
    cdef cppclass _NameType "NameType":
        _NameType() 
        _NameType(const _NameType&)
        _NameType(const char *)
        _NameType(const string&)
        #_NameType& operator =(const _NameType&)
        void ToBuffer(char *) const 
        bint Match(const _NameType&) const 
        bint operator ==(const _NameType&) const 
        bint operator ==(const char *) const 
        #bint opr_ne "operator !="(const _NameType&) const 
        bint operator !=(const _NameType&) const 
        bint operator !=(const char *) const 
        const char* opr_star "operator*" () const 
        char opr_idx "operator[]"(int) const 
        string Truncated() const 


cdef class NameType:
        cdef _NameType* thisptr


cdef extern from "ArgList.h": 
    cdef cppclass _ArgList "ArgList":
        _ArgList() 
        _ArgList(const char *)
        _ArgList(const string&)
        _ArgList(const string&, const char *)
        _ArgList(const _ArgList&)
        #_ArgList& operator =(const _ArgList&)
        const string& operator[](int) const 
        const vector[string]& List() const 
        vector[string].const_iterator begin() const 
        vector[string].const_iterator end() const 
        int Nargs() const 
        bint empty() const 
        const char * ArgLine() const 
        void ClearList() 
        int SetList(const string&, const char *)
        _ArgList RemainingArgs() 
        void AddArg(const string&)
        void MarkArg(int)
        bint CheckForMoreArgs() const 
        void PrintList() const 
        void PrintDebug() const 
        void RemoveFirstArg() 
        const char * Command() const 
        bint CommandIs(const char *) const 
        const string& GetStringNext() 
        const string& GetMaskNext() 
        const string& getNextTag() 
        bint ValidInteger(int) const 
        int IntegerAt(int) const 
        bint ValidDouble(int) const 
        int getNextInteger(int)
        double getNextDouble(double)
        const string& GetStringKey(const char *)
        int getKeyInt(const char *, int)
        double getKeyDouble(const char *, double)
        bint hasKey(const char *)
        bint Contains(const char *) const 

cdef class ArgList:
    cdef _ArgList* thisptr
