# distutils: language = c++

from ..cython_extra_header.cpp_vector cimport vector as cppvector
from ..core.topology_objects cimport _Atom, Atom, _Residue, Residue, _Molecule, Molecule
from ..core.box cimport _Box, Box
from ..core.parameter_types cimport *
from ..core.c_core cimport (_FileName, FileName, _NameType, NameType)
from ..core.c_core cimport _AtomMask, AtomMask
from ..trajectory.frame cimport _Frame, Frame
from libcpp.string cimport string
from ..core.c_core cimport _FileName, FileName, _ArgList, ArgList


ctypedef cppvector[_Atom].const_iterator atom_iterator
ctypedef cppvector[_Residue].const_iterator res_iterator
ctypedef cppvector[_Molecule].const_iterator mol_iterator
ctypedef cppvector[int].const_iterator bond_iterator

cdef extern from "CoordinateInfo.h": 
    cdef cppclass _CoordinateInfo "CoordinateInfo" nogil:
        _CoordinateInfo() 
        _CoordinateInfo(const _Box& b, bint v, bint t, bint m)
        bint HasBox() const 
        const _Box& TrajBox() const 
        bint HasVel() const 
        bint HasTemp() const 
        bint HasTime() const 
        bint HasForce() const 
        bint HasReplicaDims() const 
        void SetTime(bint m)
        void SetTemperature(bint t)
        void SetVelocity(bint v)
        void SetBox(const _Box& b)

cdef extern from "Topology.h": 
    cdef cppclass _Topology "Topology" nogil:
        _Topology() 
        int DetermineMolecules()
        void SetOffset(double oIn)
        void SetDebug(int dIn)
        void SetIpol(int iIn)
        void SetPindex(int pIn)
        void Increase_Frames(int fIn)
        void SetTag(const string& t)
        void SetVelInfo(bint v)
        void SetNrepDim(int n)
        void SetGBradiiSet(const string& s)
        void SetParmName(const string&, const _FileName&)
        void SetDistMaskRef(_Frame)
        _Atom& GetAtomView "SetAtom" (int idx)
        const string& Tag() const 
        int Ipol() const 
        int Pindex() const 
        int Natom() const 
        int Nres() const 
        int Nmol() const 
        int Nsolvent() const 
        int Nframes() const 
        int NextraPts() const 
        bint HasVelInfo() const 
        int NrepDims "NrepDim"() const 
        const string& ParmName() const 
        const _FileName& OriginalFilename() const 
        const string& GBradiiSet() const 
        bint NoRefCoords()
        int FinalSoluteRes()
        const char * c_str()
        atom_iterator begin()
        atom_iterator end()
        const _Atom& index_opr "operator[]"(int idx)
        const vector[_Atom]& Atoms()
        inline res_iterator ResStart()
        inline res_iterator ResEnd()
        const _Residue& Res(int idx)
        _Residue& SetRes(int idx)
        inline mol_iterator MolStart() const 
        inline mol_iterator MolEnd() const 
        const _Molecule& Mol(int idx) const 
        void ClearMoleculeInfo() 
        const BondArray& Bonds() const 
        const BondArray& BondsH() const 
        const BondParmArray& BondParm() const 
        void AddBond(int, int)
        int SetBondInfo(const BondArray&, const BondArray&, const BondParmArray&)
        const AngleArray& Angles() const 
        const AngleArray& AnglesH() const 
        const AngleParmArray& AngleParm() const 
        int SetAngleInfo(const AngleArray&, const AngleArray&, const AngleParmArray&)
        const DihedralArray& Dihedrals() const 
        const DihedralArray& DihedralsH() const 
        const DihedralParmArray& DihedralParm() const 
        int SetDihedralInfo(const DihedralArray&, const DihedralArray&, const DihedralParmArray&)
        const _NonbondParmType& Nonbond() const 
        int SetNonbondInfo(const _NonbondParmType&)
        inline const _NonbondType& GetLJparam(int, int) const 
        const _CapParmType& Cap() const 
        void SetCap(const _CapParmType& c)
        const _LES_ParmType& LES() const 
        void SetLES(const _LES_ParmType& l)
        const _ChamberParmType& Chamber() const 
        void SetChamber(const _ChamberParmType& c)
        inline const vector[double]& Solty() const 
        inline const vector[_NameType]& Itree() const 
        inline const vector[int]& Join() const 
        inline const vector[int]& Irotat() const 
        string TruncResAtomName(int) const 
        string AtomMaskName(int atom) const 
        string TruncResNameNum(int) const 
        int FindAtomInResidue(int, const _NameType&) const 
        #int FindResidueMaxNatom() const 
        int SoluteAtoms() const 
        int SetSolvent(const string&)
        void Summary() const 
        void Brief(const char *) const 
        void PrintAtomInfo(const string&) const 
        void PrintBondInfo(const string&) const 
        void PrintAngleInfo(const string&) const 
        void PrintDihedralInfo(const string&) const 
        void PrintMoleculeInfo(const string&) const 
        void PrintResidueInfo(const string&) const 
        int PrintChargeMassInfo(const string&, int) const 
        void PrintBonds(const BondArray&, _AtomMask&, int&) const
        void PrintAngles(const AngleArray&, const _AtomMask&, int&) const
        void PrintDihedrals(const DihedralArray&, const _AtomMask&, int&) const
        inline const _Box& ParmBox() const 
        void SetParmBox(_Box& bIn)
        int AddTopAtom(_Atom&, _Residue&)
        void AddAngle(int, int, int)
        void AddDihedral(int, int, int, int)
        void StartNewMol() 
        int CommonSetup(bint)
        int SetAmberExtra(const vector[double]&, const vector[_NameType]&, const vector[int]&, const vector[int]&)
        bint SetupIntegerMask(_AtomMask&) const 
        bint SetupCharMask(_AtomMask&) const 
        bint SetupIntegerMask(_AtomMask&, const _Frame&) const 
        bint SetupCharMask(_AtomMask&, const _Frame&) const 
        void ScaleDihedralK(double)
        _Topology* partialModifyStateByMask(const _AtomMask& m) const 
        _Topology* modifyStateByMask(const _AtomMask& m) const 
        _Topology* ModifyByMap(const vector[int]& m) const 
        int AppendTop(const _Topology &)
        # add more
        _CoordinateInfo& ParmCoordInfo() const
        double GetVDWradius(int) except +

cdef class Topology:
    cdef _Topology* thisptr
    cdef public bint _own_memory
    cdef cppvector[int] _get_atom_bond_indices(self, _Atom)

cdef extern from "ParmFile.h": 
    ctypedef enum ParmFormatType "ParmFile::ParmFormatType":
        pass
        UNKNOWN_PARM "ParmFile::UNKNOWN_PARM"
    cdef cppclass _ParmFile "ParmFile" nogil:
        @staticmethod
        void ReadOptions() 
        @staticmethod
        void WriteOptions() 
        _ParmFile() 
        int ReadTopology(_Topology&, const string&, const _ArgList&, int)
        int ReadTopology(_Topology& t, const string& n, int d)
        int WritePrefixTopology(const _Topology&, const string&, ParmFormatType, int)
        int WriteTopology(const _Topology&, const string&, const _ArgList&, ParmFormatType, int)
        int WriteTopology(const _Topology& t, const string& n, ParmFormatType f, int d)
        const _FileName ParmFilename() 


cdef class ParmFile:
    cdef _ParmFile* thisptr
