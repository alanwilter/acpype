# distutils: language = c++
from __future__ import absolute_import
from libcpp.map cimport map as cmap
from libcpp.vector cimport vector
from libcpp.string cimport string 
from ..math.cpp_math cimport _Grid, _Vec3, Vec3, _Matrix_3x3, Matrix_3x3, _Matrix
from ..trajectory.frame cimport _Frame, Frame
from ..topology.topology cimport _Topology, Topology
from ..core.c_core cimport _ArgList, ArgList, _AtomMask, AtomMask
from ..core.coordinfo cimport _CoordinateInfo


ctypedef vector[size_t] SizeArray

cdef extern from "MetaData.h": 
    ctypedef enum scalarType "MetaData::scalarType":
        pass
    cdef cppclass _MetaData "MetaData":
        string& Name()
        string& Aspect()
        string& Legend()
        void SetName(string & n)
        void SetAspect(string& a)
        void SetLegend(string& l)
        void SetScalarType(scalarType s)
        scalarType ScalarType()

cdef extern from "DataSet.h": 
    ctypedef enum DataType "DataSet::DataType":
        pass
    cdef cppclass _Dataset "DataSet":
        _Dataset() 
        _Dataset(DataType, int, int, int)
        _Dataset(const _Dataset&)
        #_Dataset& operator =(const Dataset&)
        const string& Name() const 
        void SetLegend(string)
        DataType Type() const 
        size_t Ndim() const 
        size_t Size()
        _MetaData& Meta()
        int SetMeta(_MetaData&)
        double Coord(unsigned int d, size_t p)

cdef class Dataset:
    cdef _Dataset* baseptr0
    cdef public object _base


cdef extern from "DataSet_1D.h": 
    cdef cppclass _Dataset1D "DataSet_1D" (_Dataset):
        _Dataset1D() 
        _Dataset1D(Dataset)
        # virtual methods
        #virtual ~_Dataset1D() 
        int Allocate (SizeArray)
        double Dval(size_t) const
        double Xcrd(size_t) const
        # end virtual methods
        inline bint IsTorsionArray() const 
        double Avg() const 
        double Avg(double& sd) const 
        double Min() const 
        double Max() const 
        int CrossCorr(const _Dataset1D&, _Dataset1D&, int, bint, bint) const 
        double CorrCoeff(const _Dataset1D&) const 


cdef class Dataset1D (Dataset):
    # baseptr0 is from Dataset
    cdef _Dataset1D* baseptr_1
# distutils: language = c++
cdef extern from "DataSet_double.h": 
    cdef cppclass _DatasetDouble "DataSet_double" (_Dataset1D):
        _DatasetDouble() 
        @staticmethod
        _Dataset * Alloc() 
        double& operator[](size_t idx)
        double& index_opr "operator[]"(size_t idx)
        const vector[double]& Data() const 
        void assign_opr "operator =" (const vector[double]& rhs)
        void AddElement(double d)
        void Resize(size_t sizeIn)
        size_t Size()
        int Allocate(SizeArray)
        void Add(size_t, const void *)
        double Dval(size_t idx) const 
        double Xcrd(size_t idx) const 
        void Append(const _DatasetDouble&)
        void SetNOE(double b, double bh, double r)
        double NOE_bound() const 
        double NOE_boundH() const 
        double NOE_rexp() const 
        void ShiftTorsions(double, double)


cdef class DatasetDouble (Dataset1D):
    cdef _DatasetDouble* thisptr
    cdef bint _own_memory 

cdef extern from "DataSet_float.h": 
    cdef cppclass _DatasetFloat "DataSet_float" (_Dataset1D):
        _DatasetFloat() 
        @staticmethod
        _Dataset * Alloc() 
        void AddElement(float d)
        float& operator[](size_t idx)
        float& index_opr "operator[]"(size_t idx)
        int Size()
        void Resize(size_t)

cdef class DatasetFloat (Dataset1D):
    cdef _DatasetFloat* thisptr
    cdef bint _own_memory 

cdef extern from "DataSet_integer_mem.h": 
    cdef cppclass _DatasetIntegerMem "DataSet_integer_mem" (_Dataset1D):
        _DatasetIntegerMem() 
        @staticmethod
        _Dataset * Alloc() 
        void SetElement(size_t idx, int val)
        int index_opr "operator[]"(size_t idx)
        void AddElement(int i)
        int Size()
        void Resize(size_t)
        void Add( size_t, const void* )

cdef class DatasetInteger (Dataset1D):
    cdef _DatasetIntegerMem* thisptr
    cdef bint _own_memory 

cdef extern from "DataSet_string.h": 
    cdef cppclass _DatasetString "DataSet_string" (_Dataset1D):
        _DatasetString()
        _Dataset * Alloc() 
        string& index_opr "operator[]"(size_t idx)
        void AddElement(const string& s)
        void Resize(size_t sizeIn)
        int Size()

cdef class DatasetString(Dataset1D):
    cdef _DatasetString* thisptr
    cdef bint _own_memory

cdef extern from "DataSet_Vector.h": 
    cdef cppclass _DatasetVector "DataSet_Vector" (_Dataset):
        _DatasetVector() 
        _Dataset * Alloc() 
        void Resize(size_t s)
        void Resize(size_t s, const _Vec3& v)
        bint Empty() const 
        bint HasOrigins() const
        _Vec3& index_opr "operator[]" (int i)
        const _Vec3& VXYZ(int i) const 
        const _Vec3& OXYZ(int i) const 
        void ReserveVecs(size_t n)
        void AddVxyz(const _Vec3& v)
        void AddVxyz(const _Vec3& v, const _Vec3& c)

cdef class DatasetVector (Dataset):
    cdef _DatasetVector* thisptr
    cdef bint _own_memory

cdef extern from "DataSet_2D.h": 
    # DataSet_2D.h
    ctypedef enum MatrixType "DataSet_2D::MatrixType":
        pass
    ctypedef enum MatrixKind "DataSet_2D::MatrixKind":
        pass
    cdef cppclass _Dataset2D "DataSet_2D" (_Dataset):
        _Dataset2D() 
        _Dataset2D(DataType tIn, int wIn, int pIn)
        # virtual methods
        int Allocate(SizeArray)
        int AllocateHalf(size_t) 
        int AllocateTriangle(size_t) 
        int Allocate2D(size_t, size_t)
        double GetElement(size_t, size_t) const  
        void UpdateElement(size_t, size_t, double) const  
        size_t Nrows() const  
        size_t Ncols() const  
        double * MatrixArray() const  
        MatrixKind Kind "MatrixKind"() const  
        # end virtual methods

        void Add(size_t, const void *)
        const char * MatrixTypeString(MatrixType m)
        const char * MatrixOutputString(MatrixType m)

cdef class Dataset2D (Dataset):
    cdef _Dataset2D* baseptr_1


#ctypedef Matrix[double].iterator iterator
ctypedef vector[double] Darray

cdef extern from "DataSet_MatrixDbl.h": 
    cdef cppclass _DatasetMatrixDouble "DataSet_MatrixDbl" (_Dataset2D):
        _DatasetMatrixDouble() 
        double& index_opr "operator[]"(size_t idx)
        @staticmethod
        _Dataset * Alloc() 
        size_t Size() const 
        int Allocate(SizeArray)
        int AllocateHalf(size_t x)
        int AllocateTriangle(size_t x)
        double GetElement(size_t x, size_t y) const 
        size_t Nrows() const 
        size_t Ncols() const 
        #double * MatrixArray() const # not implemented
        MatrixKind Kind() const 
        # make alias to avoid naming conflict with Dataset (DataType)
        MatrixType matType "Type"() const 
        unsigned int Nsnapshots() const 
        void IncrementSnapshots() 
        double& Element(size_t x, size_t y)
        int AddElement(double d)
        void SetElement(size_t x, size_t y, double d)
        #iterator begin() 
        #iterator end() 
        const Darray& Vect() const 
        Darray& V1() 
        void AllocateVector(size_t vsize)
        #Darray.iterator v1begin() 
        #Darray.iterator v1end() 
        void SetTypeAndKind(MatrixType tIn, MatrixKind kIn)
        void StoreMass(const Darray& mIn)
        const Darray& Mass() const 


cdef class DatasetMatrixDouble (Dataset2D):
    cdef _DatasetMatrixDouble* thisptr
    cdef public bint _own_memory


cdef extern from "DataSet_MatrixFlt.h": 
    cdef cppclass _DatasetMatrixFloat  "DataSet_MatrixFlt" (_Dataset2D):
        _DatasetMatrixFlt() 
        float& index_opr "operator[]" (size_t idx)
        @staticmethod
        _Dataset * Alloc() 


cdef class DatasetMatrixFloat(Dataset2D):
    cdef _DatasetMatrixFloat * thisptr
    cdef bint _own_memory


cdef extern from "DataSet_3D.h": 
    cdef cppclass _Dataset3D "DataSet_3D" (_Dataset):
        _Dataset3D() 
        _Dataset3D(DataType tIn, int wIn, int pIn)
        void Add(size_t, const void *)
        inline bint CalcBins(double, double, double, int&, int&, int&) const 
        inline double DX() const 
        inline double DY() const 
        inline double DZ() const 
        inline double OX() const 
        inline double OY() const 
        inline double OZ() const 
        inline double MX() const 
        inline double MY() const 
        inline double MZ() const 

cdef class Dataset3D (Dataset):
    cdef _Dataset3D* baseptr_1

cdef extern from "DataSet_GridFlt.h": 
    cdef cppclass _DatasetGridFloat "DataSet_GridFlt" (_Dataset3D):
        _DatasetGridFloat()
        float& index_opr "operator[]"(size_t idx)
        _Dataset * Alloc() 
        const _Grid[float]& InternalGrid() const 
        size_t Size() const 
        int Allocate3D(size_t x, size_t y, size_t z)
        double GetElement(int x, int y, int z) const 
        void SetElement(int x, int y, int z, float v)
        double operator[](size_t idx) const 
        size_t NX() const 
        size_t NY() const 
        size_t NZ() const 
        float GridVal(int x, int y, int z) const 
        long int CalcIndex(int i, int j, int k) const 

cdef class DatasetGridFloat (Dataset3D):
    cdef _DatasetGridFloat* thisptr
    cdef public bint _own_memory

cdef extern from "DataSet_GridDbl.h": 
    cdef cppclass _DatasetGridDouble "DataSet_GridDbl" (_Dataset3D):
        _DatasetGridDouble()
        double& index_opr "operator[]"(size_t idx)
        _Dataset * Alloc() 
        size_t Size() const 
        int Allocate3D(size_t x, size_t y, size_t z)
        double GetElement(int x, int y, int z) const 
        void SetElement(int x, int y, int z, float v)
        double operator[](size_t idx) const 
        size_t NX() const 
        size_t NY() const 
        size_t NZ() const 
        double GridVal(int x, int y, int z) const 
        long int CalcIndex(int i, int j, int k) const 

cdef class DatasetGridDouble (Dataset3D):
    cdef _DatasetGridDouble* thisptr
    cdef public bint _own_memory

cdef extern from "DataSet_Modes.h": 
    cdef cppclass _DatasetModes "DataSet_Modes" (_Dataset):
        _DatasetModes() 
        _Dataset * Alloc() 
        size_t Size() const 
        void Add(size_t, const void *)
        const Darray& AvgCrd() const 
        const Darray& Mass() const 
        int NavgCrd() const 
        double * AvgFramePtr() 
        void AllocateAvgCoords(int n)
        void SetAvgCoords(const _Dataset2D&)
        int SetModes(bint, int, int, const double *, const double *)
        int CalcEigen(const _Dataset2D&, int)
        void PrintModes() 
        int EigvalToFreq(double)
        int MassWtEigvect()
        int ReduceVectors() 
        double Eigenvalue(int i)
        double * Eigenvectors()
        double * Eigenvector(int i)
        int Nmodes() const 
        int VectorSize() const 
        #MatrixType Type() const 
        bint IsReduced() const 

cdef class DatasetModes (Dataset):
    cdef _DatasetModes* thisptr
    cdef public bint _own_memory

ctypedef cmap[double, int] TcmapType
cdef extern from "DataSet_RemLog.h": 
    cdef cppclass _DatasetRemLog "DataSet_RemLog":
        _DatasetRemLog() 
        _Dataset * Alloc() 
        void AllocateReplicas(int)
        void AddRepFrame(int rep, const _ReplicaFrame& frm)
        const _ReplicaFrame& RepFrame(int exch, int rep) const 
        int NumExchange() const 
        bint ValidEnsemble() const 
        void TrimLastExchange() 
        size_t Size() const 
        void Add(size_t, const void *)


    cdef cppclass _ReplicaFrame "DataSet_RemLog::ReplicaFrame":
        _Replica_Frame() 
        int SetTremdFrame(const char *, const TcmapType&)
        int SetHremdFrame(const char *, const vector[int]&)
        int ReplicaIdx() const 
        int PartnerIdx() const 
        int CoordsIdx() const 
        bint Success() const 
        double Temp0() const 
        double PE_X1() const 
        double PE_X2() const 


cdef class DatasetRemLog:
    cdef _DatasetRemLog* thisptr

cdef class ReplicaFrame:
    cdef _ReplicaFrame* thisptr



cdef extern from "DataSet_Mat3x3.h": 
    ctypedef vector[_Matrix_3x3].iterator mat_iterator
    cdef cppclass _DatasetMatrix3x3 "DataSet_Mat3x3" (_Dataset):
        _DatasetMatrix3x3()
        @staticmethod
        _Dataset * Alloc() 
        bint Empty()
        void AddMat3x3(_Matrix_3x3)
        mat_iterator begin()
        mat_iterator end()
        _Matrix_3x3& operator[](int i)


cdef class DatasetMatrix3x3(Dataset):
    cdef _DatasetMatrix3x3* thisptr
    cdef bint _own_memory 


cdef extern from "DataSet_Mesh.h": 
    cdef cppclass _DatasetMesh "DataSet_Mesh" (_Dataset1D):
        _DatasetMesh()
        _DatasetMesh(int, double, double)
        _Dataset * Alloc() 
        size_t Size() const 
        int Allocate(SizeArray)
        void Add(size_t, const void *)
        inline void AddXY(double, double)
        double X(int i)
        double Y(int i)
        void CalculateMeshX(int, double, double)
        int SetMeshXY(const _Dataset1D&)
        double Integrate_Trapezoid(_DatasetMesh&) const 
        double Integrate_Trapezoid() const 
        int SetSplinedMeshY(const vector[double]&, const vector[double]&)
        int SetSplinedMesh(const _Dataset1D&)
        int LinearRegression(double&, double&, double&, bint) const 

cdef class DatasetMesh(Dataset1D):
    cdef _DatasetMesh* thisptr
    cdef public bint _own_memory


cdef extern from "DataSet_Coords.h": 
    cdef cppclass _DatasetCoords "DataSet_Coords" (_Dataset):
        _DatasetCoords() 
        _DatasetCoords(DataType)
        #virtual ~_DatasetCoords() 
        _Frame AllocateFrame() const 
        
        # virtual methods
        void AddFrame(const _Frame&) 
        void SetCRD(int, const _Frame&) 
        void GetFrame(int, _Frame&) 
        void GetFrame(int, _Frame&, const _AtomMask&) 
        # end virtual methods

        #void SetTopology(const _Topology&)
        inline _Topology& Top()
        void CoordsSetup(const _Topology&, const _CoordinateInfo &)
        const _CoordinateInfo& CoordsInfo()


cdef class DatasetCoords (Dataset):
    # Dataset has baseptr0
    cdef _DatasetCoords* baseptr_1
    cdef Topology _top
    cdef bint _own_memory

    # use tmpfarray object to hold Frame or Trajectory 
    # (if we want to use dset[0][0] correctly)
    cdef object tmpfarray
# distutils: language = c++

cdef extern from "DataSet_Coords_CRD.h": 
    cdef cppclass _DatasetCoordsCRD "DataSet_Coords_CRD" (_DatasetCoords):
        _DatasetCoordsCRD() 
        @staticmethod
        _Dataset * Alloc() 
        size_t Size() const 
        int Allocate1D(size_t)
        void Add(size_t, const void *)
        double Dval(size_t)const 
        double Xcrd(size_t idx)const 
        inline void AddFrame(const _Frame& fIn)
        inline void GetFrame(int idx, _Frame & fIn)
        inline void GetFrame(int idx, _Frame & fIn, const _AtomMask& mIn)
        inline void SetCRD(int idx, const _Frame& fIn)


cdef class DatasetCoordsCRD (DatasetCoords):
    cdef _DatasetCoordsCRD* thisptr

cdef extern from "DataSet_Coords_REF.h": 
    cdef cppclass _DatasetCoordsRef "DataSet_Coords_REF" (_DatasetCoords):
        _DatasetCoordsRef() 

        # turn off those methods since they are in parent class
        @staticmethod
        _Dataset * Alloc() 
        size_t Size() const 

        int LoadRef(const string&, const _Topology&, int)
        int SetupRef_Frame(const string&, const string&, const _Topology&, _ArgList&, int)
        int SetupRef_Frame(_DatasetCoords *, const string&, int, int)
        int StripRef(const string&)
        int StripRef(const _AtomMask&)
        const _Frame& RefFrame() const 
        int RefIndex() const 
        #void SetCRD(int idx, _Frame& fIn)

cdef class DatasetCoordsRef (DatasetCoords):
    cdef _DatasetCoordsRef* thisptr

cdef extern from "DataSet_Topology.h":
    cdef cppclass _DatasetTopology "DataSet_Topology" (_Dataset):
        _DatasetTopology()
        _Dataset * Alloc()
        size_t Size() const
        int LoadTopFromFile(const _ArgList&, int)
        void SetTop(const _Topology& t)
        int StripTop(const string&)
        void SetPindex(int p)
        _Topology * TopPtr()
        _Topology Top()


cdef class DatasetTopology (Dataset):
    cdef _DatasetTopology* thisptr
    cdef bint _own_memory
