# distutils: language = c++
from libcpp.string cimport string
from ...core.c_core cimport (_DispatchObject, DispatchObject,  FunctPtr)
from ...core.box cimport Box
from ...datafiles.datafiles cimport  _DataFileList, DataFileList
from ...topology.topology cimport _Topology, Topology
from ...core.c_core cimport _ArgList, ArgList
from ...datasets.c_datasetlist cimport _DatasetList, DatasetList
from ...trajectory.frame cimport _Frame, Frame
from ...core.coordinfo cimport _CoordinateInfo, CoordinateInfo


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


    cdef cppclass _ActionSetup "ActionSetup":
        _ActionSetup()
        void Set(_Topology * p, const _CoordinateInfo& c, int n)
        _ActionSetup(_Topology * topIn, const _CoordinateInfo& cInfoIn, int nIn)
        const _Topology& Top() const
        _Topology * TopAddress()
        const _CoordinateInfo& CoordInfo() const
        int Nframes() const
        void SetTopology(_Topology * p)
        void SetCoordInfo(_CoordinateInfo * c)


    cdef cppclass _ActionFrame "ActionFrame":
        _ActionFrame()
        _ActionFrame(_Frame * fIn, int trajout_index)
        const _Frame& Frm() const
        _Frame& ModifyFrm()
        _Frame * _FramePtr()
        void SetFrame(_Frame * f)


cdef extern from "Action.h": 
    # Action.h
    ctypedef enum RetType "Action::RetType":
        OK "Action::OK"
        ERR "Action::ERR"
        USE_ORIGINAL_FRAME "Action::USE_ORIGINAL_FRAME"
        SUPPRESS_COORD_OUTPUT "Action::SUPPRESS_COORD_OUTPUT"
        SKIP "Action::SKIP"
        MODIFY_TOPOLOGY "Action::MODIFY_TOPOLOGY"
        MODIFY_COORDS "Action::MODIFY_COORDS"
    cdef cppclass _Action "Action" nogil:
        RetType Init(_ArgList&, _ActionInit&, int)
        RetType Setup(_ActionSetup&)
        RetType DoAction(int, _ActionFrame&)
        void Print()


cdef class Action:
    cdef _Action* baseptr
    cdef public int n_frames
    cdef bint top_is_processed
    cdef object top
    cdef public object _dslist
    cdef public object _dflist
    cdef public object _command
    # create pointer to pass to ActionList
    cdef bint own_memory # Stop mark for generated script


cdef extern from "Action_Align.h":
    cdef cppclass _Action_Align "Action_Align" (_Action) nogil:
        _Action_Align()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_Align(Action):
    cdef _Action_Align* thisptr


cdef extern from "Action_Angle.h":
    cdef cppclass _Action_Angle "Action_Angle" (_Action) nogil:
        _Action_Angle()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_Angle(Action):
    cdef _Action_Angle* thisptr


cdef extern from "Action_AreaPerMol.h":
    cdef cppclass _Action_AreaPerMol "Action_AreaPerMol" (_Action) nogil:
        _Action_AreaPerMol()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_AreaPerMol(Action):
    cdef _Action_AreaPerMol* thisptr


cdef extern from "Action_AtomMap.h":
    cdef cppclass _Action_AtomMap "Action_AtomMap" (_Action) nogil:
        _Action_AtomMap()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_AtomMap(Action):
    cdef _Action_AtomMap* thisptr


cdef extern from "Action_AtomicCorr.h":
    cdef cppclass _Action_AtomicCorr "Action_AtomicCorr" (_Action) nogil:
        _Action_AtomicCorr()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_AtomicCorr(Action):
    cdef _Action_AtomicCorr* thisptr


cdef extern from "Action_AtomicFluct.h":
    cdef cppclass _Action_AtomicFluct "Action_AtomicFluct" (_Action) nogil:
        _Action_AtomicFluct()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_AtomicFluct(Action):
    cdef _Action_AtomicFluct* thisptr


cdef extern from "Action_AutoImage.h":
    cdef cppclass _Action_AutoImage "Action_AutoImage" (_Action) nogil:
        _Action_AutoImage()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_AutoImage(Action):
    cdef _Action_AutoImage* thisptr


cdef extern from "Action_Average.h":
    cdef cppclass _Action_Average "Action_Average" (_Action) nogil:
        _Action_Average()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_Average(Action):
    cdef _Action_Average* thisptr


cdef extern from "Action_Bounds.h":
    cdef cppclass _Action_Bounds "Action_Bounds" (_Action) nogil:
        _Action_Bounds()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_Bounds(Action):
    cdef _Action_Bounds* thisptr


cdef extern from "Action_Box.h":
    cdef cppclass _Action_Box "Action_Box" (_Action) nogil:
        _Action_Box()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_Box(Action):
    cdef _Action_Box* thisptr


cdef extern from "Action_Center.h":
    cdef cppclass _Action_Center "Action_Center" (_Action) nogil:
        _Action_Center()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_Center(Action):
    cdef _Action_Center* thisptr


cdef extern from "Action_Channel.h":
    cdef cppclass _Action_Channel "Action_Channel" (_Action) nogil:
        _Action_Channel()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_Channel(Action):
    cdef _Action_Channel* thisptr


cdef extern from "Action_CheckChirality.h":
    cdef cppclass _Action_CheckChirality "Action_CheckChirality" (_Action) nogil:
        _Action_CheckChirality()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_CheckChirality(Action):
    cdef _Action_CheckChirality* thisptr


cdef extern from "Action_CheckStructure.h":
    cdef cppclass _Action_CheckStructure "Action_CheckStructure" (_Action) nogil:
        _Action_CheckStructure()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_CheckStructure(Action):
    cdef _Action_CheckStructure* thisptr


cdef extern from "Action_Closest.h":
    cdef cppclass _Action_Closest "Action_Closest" (_Action) nogil:
        _Action_Closest()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_Closest(Action):
    cdef _Action_Closest* thisptr


cdef extern from "Action_ClusterDihedral.h":
    cdef cppclass _Action_ClusterDihedral "Action_ClusterDihedral" (_Action) nogil:
        _Action_ClusterDihedral()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_ClusterDihedral(Action):
    cdef _Action_ClusterDihedral* thisptr


cdef extern from "Action_Contacts.h":
    cdef cppclass _Action_Contacts "Action_Contacts" (_Action) nogil:
        _Action_Contacts()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_Contacts(Action):
    cdef _Action_Contacts* thisptr


cdef extern from "Action_CreateCrd.h":
    cdef cppclass _Action_CreateCrd "Action_CreateCrd" (_Action) nogil:
        _Action_CreateCrd()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_CreateCrd(Action):
    cdef _Action_CreateCrd* thisptr


cdef extern from "Action_DNAionTracker.h":
    cdef cppclass _Action_DNAionTracker "Action_DNAionTracker" (_Action) nogil:
        _Action_DNAionTracker()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_DNAionTracker(Action):
    cdef _Action_DNAionTracker* thisptr


cdef extern from "Action_DSSP.h":
    cdef cppclass _Action_DSSP "Action_DSSP" (_Action) nogil:
        _Action_DSSP()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_DSSP(Action):
    cdef _Action_DSSP* thisptr


cdef extern from "Action_Density.h":
    cdef cppclass _Action_Density "Action_Density" (_Action) nogil:
        _Action_Density()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_Density(Action):
    cdef _Action_Density* thisptr


cdef extern from "Action_Diffusion.h":
    cdef cppclass _Action_Diffusion "Action_Diffusion" (_Action) nogil:
        _Action_Diffusion()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_Diffusion(Action):
    cdef _Action_Diffusion* thisptr


cdef extern from "Action_Dihedral.h":
    cdef cppclass _Action_Dihedral "Action_Dihedral" (_Action) nogil:
        _Action_Dihedral()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_Dihedral(Action):
    cdef _Action_Dihedral* thisptr


cdef extern from "Action_Dipole.h":
    cdef cppclass _Action_Dipole "Action_Dipole" (_Action) nogil:
        _Action_Dipole()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_Dipole(Action):
    cdef _Action_Dipole* thisptr


cdef extern from "Action_DistRmsd.h":
    cdef cppclass _Action_DistRmsd "Action_DistRmsd" (_Action) nogil:
        _Action_DistRmsd()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_DistRmsd(Action):
    cdef _Action_DistRmsd* thisptr


cdef extern from "Action_Distance.h":
    cdef cppclass _Action_Distance "Action_Distance" (_Action) nogil:
        _Action_Distance()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_Distance(Action):
    cdef _Action_Distance* thisptr


cdef extern from "Action_Energy.h":
    cdef cppclass _Action_Energy "Action_Energy" (_Action) nogil:
        _Action_Energy()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_Energy(Action):
    cdef _Action_Energy* thisptr


cdef extern from "Action_Esander.h":
    cdef cppclass _Action_Esander "Action_Esander" (_Action) nogil:
        _Action_Esander()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_Esander(Action):
    cdef _Action_Esander* thisptr


cdef extern from "Action_FilterByData.h":
    cdef cppclass _Action_FilterByData "Action_FilterByData" (_Action) nogil:
        _Action_FilterByData()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_FilterByData(Action):
    cdef _Action_FilterByData* thisptr


cdef extern from "Action_FixAtomOrder.h":
    cdef cppclass _Action_FixAtomOrder "Action_FixAtomOrder" (_Action) nogil:
        _Action_FixAtomOrder()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_FixAtomOrder(Action):
    cdef _Action_FixAtomOrder* thisptr


cdef extern from "Action_FixImagedBonds.h":
    cdef cppclass _Action_FixImagedBonds "Action_FixImagedBonds" (_Action) nogil:
        _Action_FixImagedBonds()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_FixImagedBonds(Action):
    cdef _Action_FixImagedBonds* thisptr


cdef extern from "Action_GIST.h":
    cdef cppclass _Action_GIST "Action_GIST" (_Action) nogil:
        _Action_GIST()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_GIST(Action):
    cdef _Action_GIST* thisptr


cdef extern from "Action_Grid.h":
    cdef cppclass _Action_Grid "Action_Grid" (_Action) nogil:
        _Action_Grid()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_Grid(Action):
    cdef _Action_Grid* thisptr


cdef extern from "Action_GridFreeEnergy.h":
    cdef cppclass _Action_GridFreeEnergy "Action_GridFreeEnergy" (_Action) nogil:
        _Action_GridFreeEnergy()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_GridFreeEnergy(Action):
    cdef _Action_GridFreeEnergy* thisptr


cdef extern from "Action_HydrogenBond.h":
    cdef cppclass _Action_HydrogenBond "Action_HydrogenBond" (_Action) nogil:
        _Action_HydrogenBond()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_HydrogenBond(Action):
    cdef _Action_HydrogenBond* thisptr


cdef extern from "Action_Image.h":
    cdef cppclass _Action_Image "Action_Image" (_Action) nogil:
        _Action_Image()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_Image(Action):
    cdef _Action_Image* thisptr


cdef extern from "Action_Jcoupling.h":
    cdef cppclass _Action_Jcoupling "Action_Jcoupling" (_Action) nogil:
        _Action_Jcoupling()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_Jcoupling(Action):
    cdef _Action_Jcoupling* thisptr


cdef extern from "Action_LESsplit.h":
    cdef cppclass _Action_LESsplit "Action_LESsplit" (_Action) nogil:
        _Action_LESsplit()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_LESsplit(Action):
    cdef _Action_LESsplit* thisptr


cdef extern from "Action_LIE.h":
    cdef cppclass _Action_LIE "Action_LIE" (_Action) nogil:
        _Action_LIE()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_LIE(Action):
    cdef _Action_LIE* thisptr


cdef extern from "Action_LipidOrder.h":
    cdef cppclass _Action_LipidOrder "Action_LipidOrder" (_Action) nogil:
        _Action_LipidOrder()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_LipidOrder(Action):
    cdef _Action_LipidOrder* thisptr


cdef extern from "Action_MakeStructure.h":
    cdef cppclass _Action_MakeStructure "Action_MakeStructure" (_Action) nogil:
        _Action_MakeStructure()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_MakeStructure(Action):
    cdef _Action_MakeStructure* thisptr


cdef extern from "Action_Mask.h":
    cdef cppclass _Action_Mask "Action_Mask" (_Action) nogil:
        _Action_Mask()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_Mask(Action):
    cdef _Action_Mask* thisptr


cdef extern from "Action_Matrix.h":
    cdef cppclass _Action_Matrix "Action_Matrix" (_Action) nogil:
        _Action_Matrix()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_Matrix(Action):
    cdef _Action_Matrix* thisptr


cdef extern from "Action_MinImage.h":
    cdef cppclass _Action_MinImage "Action_MinImage" (_Action) nogil:
        _Action_MinImage()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_MinImage(Action):
    cdef _Action_MinImage* thisptr


cdef extern from "Action_Molsurf.h":
    cdef cppclass _Action_Molsurf "Action_Molsurf" (_Action) nogil:
        _Action_Molsurf()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_Molsurf(Action):
    cdef _Action_Molsurf* thisptr


cdef extern from "Action_MultiDihedral.h":
    cdef cppclass _Action_MultiDihedral "Action_MultiDihedral" (_Action) nogil:
        _Action_MultiDihedral()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_MultiDihedral(Action):
    cdef _Action_MultiDihedral* thisptr


cdef extern from "Action_MultiVector.h":
    cdef cppclass _Action_MultiVector "Action_MultiVector" (_Action) nogil:
        _Action_MultiVector()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_MultiVector(Action):
    cdef _Action_MultiVector* thisptr


cdef extern from "Action_NAstruct.h":
    cdef cppclass _Action_NAstruct "Action_NAstruct" (_Action) nogil:
        _Action_NAstruct()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_NAstruct(Action):
    cdef _Action_NAstruct* thisptr


cdef extern from "Action_NMRrst.h":
    cdef cppclass _Action_NMRrst "Action_NMRrst" (_Action) nogil:
        _Action_NMRrst()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_NMRrst(Action):
    cdef _Action_NMRrst* thisptr


cdef extern from "Action_NativeContacts.h":
    cdef cppclass _Action_NativeContacts "Action_NativeContacts" (_Action) nogil:
        _Action_NativeContacts()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_NativeContacts(Action):
    cdef _Action_NativeContacts* thisptr


cdef extern from "Action_OrderParameter.h":
    cdef cppclass _Action_OrderParameter "Action_OrderParameter" (_Action) nogil:
        _Action_OrderParameter()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_OrderParameter(Action):
    cdef _Action_OrderParameter* thisptr


cdef extern from "Action_Outtraj.h":
    cdef cppclass _Action_Outtraj "Action_Outtraj" (_Action) nogil:
        _Action_Outtraj()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_Outtraj(Action):
    cdef _Action_Outtraj* thisptr


cdef extern from "Action_PairDist.h":
    cdef cppclass _Action_PairDist "Action_PairDist" (_Action) nogil:
        _Action_PairDist()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_PairDist(Action):
    cdef _Action_PairDist* thisptr


cdef extern from "Action_Pairwise.h":
    cdef cppclass _Action_Pairwise "Action_Pairwise" (_Action) nogil:
        _Action_Pairwise()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_Pairwise(Action):
    cdef _Action_Pairwise* thisptr


cdef extern from "Action_Principal.h":
    cdef cppclass _Action_Principal "Action_Principal" (_Action) nogil:
        _Action_Principal()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_Principal(Action):
    cdef _Action_Principal* thisptr


cdef extern from "Action_Projection.h":
    cdef cppclass _Action_Projection "Action_Projection" (_Action) nogil:
        _Action_Projection()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_Projection(Action):
    cdef _Action_Projection* thisptr


cdef extern from "Action_Pucker.h":
    cdef cppclass _Action_Pucker "Action_Pucker" (_Action) nogil:
        _Action_Pucker()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_Pucker(Action):
    cdef _Action_Pucker* thisptr


cdef extern from "Action_Radgyr.h":
    cdef cppclass _Action_Radgyr "Action_Radgyr" (_Action) nogil:
        _Action_Radgyr()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_Radgyr(Action):
    cdef _Action_Radgyr* thisptr


cdef extern from "Action_Radial.h":
    cdef cppclass _Action_Radial "Action_Radial" (_Action) nogil:
        _Action_Radial()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_Radial(Action):
    cdef _Action_Radial* thisptr


cdef extern from "Action_RandomizeIons.h":
    cdef cppclass _Action_RandomizeIons "Action_RandomizeIons" (_Action) nogil:
        _Action_RandomizeIons()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_RandomizeIons(Action):
    cdef _Action_RandomizeIons* thisptr


cdef extern from "Action_Remap.h":
    cdef cppclass _Action_Remap "Action_Remap" (_Action) nogil:
        _Action_Remap()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_Remap(Action):
    cdef _Action_Remap* thisptr


cdef extern from "Action_ReplicateCell.h":
    cdef cppclass _Action_ReplicateCell "Action_ReplicateCell" (_Action) nogil:
        _Action_ReplicateCell()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_ReplicateCell(Action):
    cdef _Action_ReplicateCell* thisptr


cdef extern from "Action_Rmsd.h":
    cdef cppclass _Action_Rmsd "Action_Rmsd" (_Action) nogil:
        _Action_Rmsd()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_Rmsd(Action):
    cdef _Action_Rmsd* thisptr


cdef extern from "Action_Rotate.h":
    cdef cppclass _Action_Rotate "Action_Rotate" (_Action) nogil:
        _Action_Rotate()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_Rotate(Action):
    cdef _Action_Rotate* thisptr


cdef extern from "Action_RunningAvg.h":
    cdef cppclass _Action_RunningAvg "Action_RunningAvg" (_Action) nogil:
        _Action_RunningAvg()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_RunningAvg(Action):
    cdef _Action_RunningAvg* thisptr


cdef extern from "Action_STFC_Diffusion.h":
    cdef cppclass _Action_STFC_Diffusion "Action_STFC_Diffusion" (_Action) nogil:
        _Action_STFC_Diffusion()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_STFC_Diffusion(Action):
    cdef _Action_STFC_Diffusion* thisptr


cdef extern from "Action_Scale.h":
    cdef cppclass _Action_Scale "Action_Scale" (_Action) nogil:
        _Action_Scale()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_Scale(Action):
    cdef _Action_Scale* thisptr


cdef extern from "Action_SetVelocity.h":
    cdef cppclass _Action_SetVelocity "Action_SetVelocity" (_Action) nogil:
        _Action_SetVelocity()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_SetVelocity(Action):
    cdef _Action_SetVelocity* thisptr


cdef extern from "Action_Spam.h":
    cdef cppclass _Action_Spam "Action_Spam" (_Action) nogil:
        _Action_Spam()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_Spam(Action):
    cdef _Action_Spam* thisptr


cdef extern from "Action_Strip.h":
    cdef cppclass _Action_Strip "Action_Strip" (_Action) nogil:
        _Action_Strip()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_Strip(Action):
    cdef _Action_Strip* thisptr


cdef extern from "Action_Surf.h":
    cdef cppclass _Action_Surf "Action_Surf" (_Action) nogil:
        _Action_Surf()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_Surf(Action):
    cdef _Action_Surf* thisptr


cdef extern from "Action_SymmetricRmsd.h":
    cdef cppclass _Action_SymmetricRmsd "Action_SymmetricRmsd" (_Action) nogil:
        _Action_SymmetricRmsd()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_SymmetricRmsd(Action):
    cdef _Action_SymmetricRmsd* thisptr


cdef extern from "Action_Temperature.h":
    cdef cppclass _Action_Temperature "Action_Temperature" (_Action) nogil:
        _Action_Temperature()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_Temperature(Action):
    cdef _Action_Temperature* thisptr


cdef extern from "Action_Translate.h":
    cdef cppclass _Action_Translate "Action_Translate" (_Action) nogil:
        _Action_Translate()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_Translate(Action):
    cdef _Action_Translate* thisptr


cdef extern from "Action_Unstrip.h":
    cdef cppclass _Action_Unstrip "Action_Unstrip" (_Action) nogil:
        _Action_Unstrip()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_Unstrip(Action):
    cdef _Action_Unstrip* thisptr


cdef extern from "Action_Unwrap.h":
    cdef cppclass _Action_Unwrap "Action_Unwrap" (_Action) nogil:
        _Action_Unwrap()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_Unwrap(Action):
    cdef _Action_Unwrap* thisptr


cdef extern from "Action_Vector.h":
    cdef cppclass _Action_Vector "Action_Vector" (_Action) nogil:
        _Action_Vector()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_Vector(Action):
    cdef _Action_Vector* thisptr


cdef extern from "Action_VelocityAutoCorr.h":
    cdef cppclass _Action_VelocityAutoCorr "Action_VelocityAutoCorr" (_Action) nogil:
        _Action_VelocityAutoCorr()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_VelocityAutoCorr(Action):
    cdef _Action_VelocityAutoCorr* thisptr


cdef extern from "Action_Volmap.h":
    cdef cppclass _Action_Volmap "Action_Volmap" (_Action) nogil:
        _Action_Volmap()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_Volmap(Action):
    cdef _Action_Volmap* thisptr


cdef extern from "Action_Volume.h":
    cdef cppclass _Action_Volume "Action_Volume" (_Action) nogil:
        _Action_Volume()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_Volume(Action):
    cdef _Action_Volume* thisptr


cdef extern from "Action_Watershell.h":
    cdef cppclass _Action_Watershell "Action_Watershell" (_Action) nogil:
        _Action_Watershell()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_Watershell(Action):
    cdef _Action_Watershell* thisptr


cdef extern from "Action_XtalSymm.h":
    cdef cppclass _Action_XtalSymm "Action_XtalSymm" (_Action) nogil:
        _Action_XtalSymm()
        _DispatchObject * Alloc()
        void Help()


cdef class Action_XtalSymm(Action):
    cdef _Action_XtalSymm* thisptr
