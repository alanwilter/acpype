# distutils: language = c++
from libcpp.string cimport string
from ...core.c_core cimport (_DispatchObject, DispatchObject,  FunctPtr)
from ...datafiles.datafiles cimport  _DataFileList, DataFileList
from ...topology.topology cimport _Topology, Topology
from ...core.c_core cimport _ArgList, ArgList
from ...datasets.c_datasetlist cimport _DatasetList, DatasetList
from ...trajectory.frame cimport _Frame, Frame

cdef extern from "ActionState.h":
    cdef cppclass _ActionInit "ActionInit":
        _ActionInit()
        _ActionInit(_DatasetList& dslIn, _DataFileList& dflIn)
        _DatasetList& DSL()
        const _DatasetList& DSL()
        _DatasetList * DslPtr()
        _DataFileList& DFL()
        const _DataFileList& DFL()
        _DatasetList * DSL_Ptr()

ctypedef _ActionInit _AnalysisSetup

cdef extern from "Analysis.h": 
    ctypedef enum RetType "Analysis::RetType":
        OKANALYSIS "Analysis::OK"
        ERRANALYSIS "Analysis::ERR"

    cdef cppclass _Analysis "Analysis" nogil:
        RetType Setup(_ArgList&, _AnalysisSetup&, int)
        RetType Analyze()


cdef class Analysis:
    cdef _Analysis* baseptr

cdef extern from "Analysis_AmdBias.h": 
    cdef cppclass _Analysis_AmdBias "Analysis_AmdBias" (_Analysis) nogil:
        _Analysis_AmdBias() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_AmdBias (Analysis):
    cdef _Analysis_AmdBias* thisptr



cdef extern from "Analysis_AutoCorr.h": 
    cdef cppclass _Analysis_AutoCorr "Analysis_AutoCorr" (_Analysis) nogil:
        _Analysis_AutoCorr() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_AutoCorr (Analysis):
    cdef _Analysis_AutoCorr* thisptr



cdef extern from "Analysis_Average.h": 
    cdef cppclass _Analysis_Average "Analysis_Average" (_Analysis) nogil:
        _Analysis_Average() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_Average (Analysis):
    cdef _Analysis_Average* thisptr



cdef extern from "Analysis_Clustering.h": 
    cdef cppclass _Analysis_Clustering "Analysis_Clustering" (_Analysis) nogil:
        _Analysis_Clustering() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_Clustering (Analysis):
    cdef _Analysis_Clustering* thisptr



cdef extern from "Analysis_Corr.h": 
    cdef cppclass _Analysis_Corr "Analysis_Corr" (_Analysis) nogil:
        _Analysis_Corr() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_Corr (Analysis):
    cdef _Analysis_Corr* thisptr



cdef extern from "Analysis_CrankShaft.h": 
    cdef cppclass _Analysis_CrankShaft "Analysis_CrankShaft" (_Analysis) nogil:
        _Analysis_CrankShaft() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_CrankShaft (Analysis):
    cdef _Analysis_CrankShaft* thisptr


cdef extern from "Analysis_CrdFluct.h": 
    cdef cppclass _Analysis_CrdFluct "Analysis_CrdFluct" (_Analysis) nogil:
        _Analysis_CrdFluct() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_CrdFluct (Analysis):
    cdef _Analysis_CrdFluct* thisptr


cdef extern from "Analysis_CrossCorr.h": 
    cdef cppclass _Analysis_CrossCorr "Analysis_CrossCorr" (_Analysis) nogil:
        _Analysis_CrossCorr() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_CrossCorr (Analysis):
    cdef _Analysis_CrossCorr* thisptr


cdef extern from "Analysis_Divergence.h": 
    cdef cppclass _Analysis_Divergence "Analysis_Divergence" (_Analysis) nogil:
        _Analysis_Divergence() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_Divergence (Analysis):
    cdef _Analysis_Divergence* thisptr



cdef extern from "Analysis_FFT.h": 
    cdef cppclass _Analysis_FFT "Analysis_FFT" (_Analysis) nogil:
        _Analysis_FFT() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_FFT (Analysis):
    cdef _Analysis_FFT* thisptr



cdef extern from "Analysis_Hist.h": 
    cdef cppclass _Analysis_Hist "Analysis_Hist" (_Analysis) nogil:
        _Analysis_Hist() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_Hist (Analysis):
    cdef _Analysis_Hist* thisptr



cdef extern from "Analysis_IRED.h": 
    cdef cppclass _Analysis_IRED "Analysis_IRED" (_Analysis) nogil:
        _Analysis_IRED() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_IRED (Analysis):
    cdef _Analysis_IRED* thisptr



cdef extern from "Analysis_Integrate.h": 
    cdef cppclass _Analysis_Integrate "Analysis_Integrate" (_Analysis) nogil:
        _Analysis_Integrate() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_Integrate (Analysis):
    cdef _Analysis_Integrate* thisptr



cdef extern from "Analysis_KDE.h": 
    cdef cppclass _Analysis_KDE "Analysis_KDE" (_Analysis) nogil:
        _Analysis_KDE() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_KDE (Analysis):
    cdef _Analysis_KDE* thisptr



cdef extern from "Analysis_Lifetime.h": 
    cdef cppclass _Analysis_Lifetime "Analysis_Lifetime" (_Analysis) nogil:
        _Analysis_Lifetime() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_Lifetime (Analysis):
    cdef _Analysis_Lifetime* thisptr



cdef extern from "Analysis_Matrix.h": 
    cdef cppclass _Analysis_Matrix "Analysis_Matrix" (_Analysis) nogil:
        _Analysis_Matrix() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_Matrix (Analysis):
    cdef _Analysis_Matrix* thisptr



cdef extern from "Analysis_MeltCurve.h": 
    cdef cppclass _Analysis_MeltCurve "Analysis_MeltCurve" (_Analysis) nogil:
        _Analysis_MeltCurve() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_MeltCurve (Analysis):
    cdef _Analysis_MeltCurve* thisptr



cdef extern from "Analysis_Modes.h": 
    cdef cppclass _Analysis_Modes "Analysis_Modes" (_Analysis) nogil:
        _Analysis_Modes() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_Modes (Analysis):
    cdef _Analysis_Modes* thisptr



cdef extern from "Analysis_MultiHist.h": 
    cdef cppclass _Analysis_MultiHist "Analysis_MultiHist" (_Analysis) nogil:
        _Analysis_MultiHist() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_MultiHist (Analysis):
    cdef _Analysis_MultiHist* thisptr



cdef extern from "Analysis_Overlap.h": 
    cdef cppclass _Analysis_Overlap "Analysis_Overlap" (_Analysis) nogil:
        _Analysis_Overlap() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_Overlap (Analysis):
    cdef _Analysis_Overlap* thisptr



cdef extern from "Analysis_Regression.h": 
    cdef cppclass _Analysis_Regression "Analysis_Regression" (_Analysis) nogil:
        _Analysis_Regression() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_Regression (Analysis):
    cdef _Analysis_Regression* thisptr



cdef extern from "Analysis_RemLog.h": 
    cdef cppclass _Analysis_RemLog "Analysis_RemLog" (_Analysis) nogil:
        _Analysis_RemLog() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_RemLog (Analysis):
    cdef _Analysis_RemLog* thisptr



cdef extern from "Analysis_Rms2d.h": 
    cdef cppclass _Analysis_Rms2d "Analysis_Rms2d" (_Analysis) nogil:
        _Analysis_Rms2d() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_Rms2d (Analysis):
    cdef _Analysis_Rms2d* thisptr



cdef extern from "Analysis_RmsAvgCorr.h": 
    cdef cppclass _Analysis_RmsAvgCorr "Analysis_RmsAvgCorr" (_Analysis) nogil:
        _Analysis_RmsAvgCorr() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_RmsAvgCorr (Analysis):
    cdef _Analysis_RmsAvgCorr* thisptr



cdef extern from "Analysis_RunningAvg.h": 
    cdef cppclass _Analysis_RunningAvg "Analysis_RunningAvg" (_Analysis) nogil:
        _Analysis_RunningAvg() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_RunningAvg (Analysis):
    cdef _Analysis_RunningAvg* thisptr


cdef extern from "Analysis_Spline.h": 
    cdef cppclass _Analysis_Spline "Analysis_Spline" (_Analysis) nogil:
        _Analysis_Spline() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_Spline (Analysis):
    cdef _Analysis_Spline* thisptr



cdef extern from "Analysis_Statistics.h": 
    cdef cppclass _Analysis_Statistics "Analysis_Statistics" (_Analysis) nogil:
        _Analysis_Statistics() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_Statistics (Analysis):
    cdef _Analysis_Statistics* thisptr


cdef extern from "Analysis_Timecorr.h": 
    cdef cppclass _Analysis_Timecorr "Analysis_Timecorr" (_Analysis) nogil:
        _Analysis_Timecorr() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_Timecorr (Analysis):
    cdef _Analysis_Timecorr* thisptr


cdef extern from "Analysis_VectorMath.h": 
    cdef cppclass _Analysis_VectorMath "Analysis_VectorMath" (_Analysis) nogil:
        _Analysis_VectorMath() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_VectorMath (Analysis):
    cdef _Analysis_VectorMath* thisptr


cdef extern from "Analysis_Rotdif.h": 
    cdef cppclass _Analysis_Rotdif "Analysis_Rotdif" (_Analysis) nogil:
        _Analysis_Rotdif() 
        _DispatchObject * Alloc() 
        void Help()

cdef class Analysis_Rotdif(Analysis):
    cdef _Analysis_Rotdif* thisptr

cdef extern from "Analysis_LowestCurve.h": 
    cdef cppclass _Analysis_LowestCurve  "Analysis_LowestCurve" (_Analysis) nogil:
        _Analysis_LowestCurve() 
        _DispatchObject * Alloc() 
        void Help()

cdef class Analysis_LowestCurve(Analysis):
    cdef _Analysis_LowestCurve* thisptr


cdef extern from "Analysis_PhiPsi.h":
    cdef cppclass _Analysis_PhiPsi "Analysis_PhiPsi" (_Analysis) nogil:
        _Analysis_PhiPsi()
        _DispatchObject * Alloc()
        void Help()


cdef class Analysis_PhiPsi (Analysis):
    cdef _Analysis_PhiPsi* thisptr


cdef extern from "Analysis_TI.h":
    cdef cppclass _Analysis_TI "Analysis_TI" (_Analysis) nogil:
        _Analysis_TI()
        _DispatchObject * Alloc()
        void Help()


cdef class Analysis_TI (Analysis):
    cdef _Analysis_TI* thisptr


cdef extern from "Analysis_Wavelet.h":
    cdef cppclass _Analysis_Wavelet "Analysis_Wavelet" (_Analysis) nogil:
        _Analysis_Wavelet()
        _DispatchObject * Alloc()
        void Help()


cdef class Analysis_Wavelet (Analysis):
    cdef _Analysis_Wavelet* thisptr


cdef extern from "Analysis_CurveFit.h":
    cdef cppclass _Analysis_CurveFit "Analysis_CurveFit" (_Analysis) nogil:
        _Analysis_CurveFit()
        _DispatchObject * Alloc()
        void Help()


cdef class Analysis_CurveFit (Analysis):
    cdef _Analysis_CurveFit* thisptr


cdef extern from "Analysis_Multicurve.h":
    cdef cppclass _Analysis_Multicurve "Analysis_Multicurve" (_Analysis) nogil:
        _Analysis_Multicurve()
        _DispatchObject * Alloc()
        void Help()


cdef class Analysis_Multicurve (Analysis):
    cdef _Analysis_Multicurve* thisptr


cdef extern from "Analysis_State.h":
    cdef cppclass _Analysis_State "Analysis_State" (_Analysis) nogil:
        _Analysis_State()
        _DispatchObject * Alloc()
        void Help()


cdef class Analysis_State (Analysis):
    cdef _Analysis_State* thisptr


cdef extern from "Analysis_HausdorffDistance.h":
    cdef cppclass _Analysis_Hausdorff "Analysis_HausdorffDistance" (_Analysis) nogil:
        _Analysis_Hausdorff()
        _DispatchObject * Alloc()
        void Help()

cdef class Analysis_Hausdorff (Analysis):
    cdef _Analysis_Hausdorff* thisptr
