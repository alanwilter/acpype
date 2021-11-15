# distutil: language = c++
from libcpp.vector cimport vector

cdef extern from "Grid.h":
    #ctypedef GridIterator "Grid::ArrayIterator[T]" iterator
    cdef cppclass _Grid "Grid" [T]:
        _Grid()
        #~Grid() 
        _Grid(const _Grid &)
        #not supported yet>
        #Grid & operator =(const Grid &)
        #T & operator [](size_t idx)
        T& index_opr "operator[]"(size_t idx)
        T& operator[](size_t idx)
        size_t size() const 
        int resize(size_t, size_t, size_t)
        size_t NX() const 
        size_t NY() const 
        size_t NZ() const 
        long int incrementBy(int, int, int, const T &)
        void setGrid(int, int, int, const T &)
        const T& element(int, int, int)const 
        long int CalcIndex(int x, int y, int z)const 
        #iterator begin() 
        #iterator end() 

cdef class Grid:
    pass
    cdef _Grid[float]* thisptr

cdef extern from "Matrix.h":
    cdef cppclass _Matrix "Matrix" [T]:
        _Matrix()
        _Matrix(const _Matrix&)
        T& operator[](size_t idx)
        const T& operator[] (size_t)
        size_t size()
        int resize(size_t, size_t)
        const T& element(int, int) const
        T& element(int, int)
        size_t Nrows()
        size_t Ncols()
        int addElement(const T&)
        void setElement(int, int, const T&)
        const T* Ptr()
        T* Ptr()
        size_t CalcIndex(int, int)


cdef extern from "Matrix_3x3.h": 
    cdef cppclass _Matrix_3x3 "Matrix_3x3":
        _Matrix_3x3() 
        _Matrix_3x3(const _Matrix_3x3&)
        _Matrix_3x3(const double *)
        _Matrix_3x3(double)
        _Matrix_3x3(double, double, double)
        #_Matrix_3x3& operator =(const _Matrix_3x3&)
        double operator[](int idx) const 
        double& operator[](int idx)
        _Vec3 Row1() 
        _Vec3 Row2() 
        _Vec3 Row3() 
        _Vec3 Col1() 
        _Vec3 Col2() 
        _Vec3 Col3() 
        void Zero() 
        void Print(const char *) const 
        int Diagonalize(_Vec3&)
        int Diagonalize_Sort(_Vec3&)
        int Diagonalize_Sort_Chirality(_Vec3&, int)
        void Transpose() 
        _Matrix_3x3& star_equal "operator *=" (const _Matrix_3x3&)
        void RotationAroundZ(double, double)
        void RotationAroundY(double, double)
        void CalcRotationMatrix(const _Vec3&, double)
        void CalcRotationMatrix(double, double, double)
        double RotationAngle() 
        _Vec3 AxisOfRotation(double)
        _Vec3 operator *(const _Vec3& rhs) const 
        #_Vec3 TransposeMult(const _Vec3& rhs) const  #not yet implemented in cpptraj?
        _Matrix_3x3 operator *(const _Matrix_3x3&) const 
        _Matrix_3x3 TransposeMult(const _Matrix_3x3&) const 
        #const double * Dptr() const 
        double * Dptr() 

cdef class Matrix_3x3:
    cdef _Matrix_3x3* thisptr


cdef extern from "Vec3.h": 
    cdef cppclass _Vec3 "Vec3":
        _Vec3() 
        _Vec3(const _Vec3& rhs)
        _Vec3(double vx, double vy, double vz)
        _Vec3(double vxyz)
        _Vec3(const double * XYZ)
        _Vec3(const float * XYZ)
        _Vec3(const int * XYZ)
        #_Vec3& operator =(const _Vec3& rhs)
        void Assign(const double * XYZ)
        void divequal "operator /=" (double xIn)
        _Vec3 operator /(double xIn) const 
        void mulequal "operator *=" (double xIn)
        _Vec3 operator *(double xIn) const 
        void addequal "operator +=" (double xIn)
        _Vec3 operator +(double xIn) const 
        void subequal "operator -=" (const _Vec3& rhs)
        _Vec3 operator -(const _Vec3& rhs) const 
        void addequal "operator +=" (const _Vec3& rhs)
        _Vec3 operator +(const _Vec3& rhs) const 
        double operator *(const _Vec3& rhs) const 
        _Vec3 Cross(const _Vec3& rhs) const 
        #double operator[](int idx) const 
        #double& index_opr "operator[]"(int idx)
        double& index_opr "operator[]"(int idx) const
        double Magnitude2() const 
        void Zero() 
        bint IsZero() const 
        void Neg() 
        void SetVec(double vx, double vy, double vz)
        double Normalize() 
        void Print(const char *) const 
        double Angle(const _Vec3&) const 
        double SignedAngle(const _Vec3&, const _Vec3&) const 
        #const double * Dptr() const 
        double * Dptr() 

cdef class Vec3:
    cdef _Vec3* thisptr
    cdef bint _own_memory
