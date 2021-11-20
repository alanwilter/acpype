# distutils: language = c++
from libcpp.string cimport string

# seriously I need to use absolute import here
from pytraj.datasets.c_datasets cimport _Dataset, Dataset, DataType, _MetaData
from ..cython_extra_header.cpp_vector cimport vector as cppvector

ctypedef cppvector[_Dataset*] DataListType
ctypedef cppvector[_Dataset*].const_iterator const_iterator

cdef extern from "DataSetList.h": 
    cdef cppclass _DatasetList "DataSetList":
        _DatasetList() 
        #~_DatasetList() 
        void Clear() 
        void ClearTop() 
        void ClearRef() 
        _DatasetList& addequal "operator +="(const _DatasetList&)
        const_iterator begin() const 
        const_iterator end() const 
        bint empty() const 
        size_t size() const 
        int EnsembleNum() const 
        void RemoveSet(const_iterator)
        void RemoveSet(_Dataset *)
        _Dataset * index_opr "operator[]"(int didx)
        void SetDebug(int)
        void SetEnsembleNum(int i)
        void AllocateSets(long int)
        void SetPrecisionOfDataSets(const string&, int, int)
        void SynchronizeData()
        _Dataset * GetSet(const string&, int, const string&) const 
        _Dataset * GetDataSet(const string&) const 
        _DatasetList GetMultipleSets(const string&) const 
        string GenerateDefaultName(const char *) const 
        _Dataset * AddSet(DataType, const string&, const char *)
        _Dataset * AddSetIdx(DataType, const string&, int)
        _Dataset * AddSetAspect(DataType, const string&, const string&)
        _Dataset * AddSetIdxAspect(DataType, const string&, int, const string&)
        _Dataset * AddSetIdxAspect(DataType, const string&, int, const string&, const string&)
        void AddCopyOfSet(_Dataset *)
        void AddSet(_Dataset *)
        void List() const 
        void SynchronizeData() 
        _Dataset * FindSetOfType(const string&, DataType) const 
        _Dataset * FindCoordsSet(const string&)
        _Dataset* GetReferenceFrame(string name_tag)
        #ReferenceFrame GetReferenceFrame(ArgList&) const;
        #void ListReferenceFrames() const;


cdef class DatasetList:
    cdef _DatasetList* thisptr
    cdef bint _own_memory
    cdef list _parent_lists
