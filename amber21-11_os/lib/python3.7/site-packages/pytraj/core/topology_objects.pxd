# distutils: language = c++
from libcpp.vector cimport vector
from libcpp.string cimport string
#from libcpp.set cimport set
from .c_core cimport _NameType, NameType
from ..cython_extra_header.cpp_vector cimport vector as cppvector


ctypedef cppvector[int].const_iterator bond_iterator
ctypedef cppvector[int].const_iterator excluded_iterator

cdef extern from "Atom.h": 
    ctypedef enum AtomicElementType "Atom::AtomicElementType":
        pass
    cdef cppclass _Atom "Atom":
        _Atom() 
        #virtual ~_Atom() 
        _Atom(const _NameType&, char, const char *)
        _Atom(const _NameType&, const _NameType&, double)
        _Atom(const _NameType&, double, double, const _NameType&)
        _Atom(const _NameType&, double, double, int, double, int, const _NameType&, double, double)
        _Atom(const _Atom&)
        void swap(_Atom&, _Atom&)
        #_Atom& operator =(_Atom)
        inline bond_iterator bondbegin() const 
        inline bond_iterator bondend() const 
        inline excluded_iterator excludedbegin() const 
        inline excluded_iterator excludedend() const 
        void SetResNum(int resnumIn)
        void SetMol(int molIn)
        void SetCharge(double qin)
        void SetGBradius(double rin)
        inline bint NoMol() const 
        inline const char * c_str() const 
        inline int ResNum() const 
        inline AtomicElementType Element() const 
        inline int AtomicNumber() const 
        inline const char * ElementName() const 
        inline const _NameType& Name() const 
        inline const _NameType& Type() const 
        inline int TypeIndex() const 
        inline int MolNum() const 
        inline char ChainID() const 
        inline int Nbonds() const 
        inline int Nexcluded() const 
        inline double Mass() const 
        inline double Charge() const 
        inline double Polar() const 
        inline double GBRadius() const 
        inline double Screen() const 
        void AddBond(int)
        void ClearBonds() 
        void SortBonds() 
        bint IsBondedTo(int)
        #void AddExclusionList(const set[int]&)
        @staticmethod
        double GetBondLength(AtomicElementType, AtomicElementType)


cdef class Atom:
    cdef _Atom* thisptr
    cdef int _index
    cdef public object resname
    cdef bint own_memory

# distutils: language = c++
from pytraj.core.c_core cimport _NameType, NameType

cdef extern from "Residue.h": 
    cdef cppclass _Residue "Residue":
        _Residue()
        _Residue(int onum, const _NameType& resname, int first_AtomIn)
        _Residue(_NameType& n, int r, char ic, char cid)
        inline void SetLastAtom(int i)
        inline void SetOriginalNum(int i)
        inline int FirstAtom() const 
        inline int LastAtom() const 
        inline int OriginalResNum() const 
        inline const char * c_str() const 
        inline const _NameType& Name() const 
        inline int NumAtoms() const 
        inline bint NameIsSolvent() const 

cdef class Residue:
    cdef _Residue* thisptr

#void  distutils: language = c++


cdef extern from "Molecule.h": 
    cdef cppclass _Molecule "Molecule":
        _Molecule()
        _Molecule(int begin, int end)
        void SetFirst(int begin)
        void SetLast(int last)
        void SetSolvent() 
        void SetNoSolvent() 
        inline int BeginAtom() const 
        inline int EndAtom() const 
        inline bint IsSolvent() const 
        inline int NumAtoms() const 

cdef class Molecule:
    cdef _Molecule* thisptr

