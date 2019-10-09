#ifndef MDGX_H
#define MDGX_H

//-----------------------------------------------------------------------------
// This is the Application Programmer Interface as a C header file to the mdgx
// library libmdgx.  This is a pruned version of the mdgxhcat.sh output;
// unneeded stuff was removed according to Phenix getmdgxfrc usage.
// See the tail for functions not currently used but possibly needed.
// For problem cases try using the header produced by mdgxhcat.sh.
// Potential future work includes a logical reorganization, merging
// of mdgx/Wrappers.c and xtalutil/Phenix/getmdgxfrc.c, and other
// techniques to employ information hiding and not expose of guts of mdgx.
// Scott Brozell and Pawel Janowski June 9, 2014
//-----------------------------------------------------------------------------

#include <sys/time.h>
#ifdef MPI
#include <mpi.h>
#endif
#include "fftw3.h"

// String length variables
#define MAXNAME 512

typedef struct BSplineCoeff bcof;

struct IMatrix {
  int row;
  int col;
  int* data;
  int** map;
};
typedef struct IMatrix imat;

struct DMatrix {
  int row;
  int col;
  int pfft;
  double* data;
  double** map;
  fftw_complex* fdata;
  fftw_complex** fmap;
};
typedef struct DMatrix dmat;

struct CMatrix {
  int row;
  int col;
  char* data;
  char** map;
};
typedef struct CMatrix cmat;

struct FloatBook {
  int pag;                  // Number of pages (bins in x-dimension)
  int row;                  // Number of rows (bins in y-dimension)
  int col;                  // Number of columns (bins in z-dimension)
  int isortho;              // Flag to indicate that the grid is orthogonal
  double orig[3];           // The Cartesian origin of the grid
  dmat U;                   // The transformation matrix that puts everything
                            //   in units of the grid length vectors
  dmat invU;                // The inverse of U
  dmat L;                   // The transformation matrix that puts everything
                            //   in units of the grid BIN length vectors
  dmat invL;                // The inverse of L
  float* data;              // Grid information
  float*** map;             // Three-dimensional map into data
};
typedef struct FloatBook fbook;

struct DoubleBook {
  int pag;
  int row;
  int col;
  int pfft;
  double* data;
  fftw_complex* fdata;
  double*** map;
  fftw_complex*** fmap;
};
typedef struct DoubleBook dbook;

//-----------------------------------------------------------------------------
// PROCESSCELLGROUP: information about the size and composition of a group of
//                   cells controlled by one process which must communicate
//                   with another process.                  
//-----------------------------------------------------------------------------
struct ProcessCellGroup {
  int ncell;       // The number of cells in this group
  int partner;     // The origin or destination process with which this group
                   //   communicates as referenced in MPI_COMM_WORLD
  int BaseID;      // The message base, as determined by the rank of the
                   //   sending process in the cell grid communicator and the
                   //   direction of the message.  The message base is further
                   //   offset by increments for each function that relies on
                   //   such messages.
  int* cellpt;     // Pointer to the segment of a list of all cells controlled
                   //   by this process; the segment lists cells in this group
};
typedef struct ProcessCellGroup pcgrp;

//-----------------------------------------------------------------------------
// ATOMSHAREPLAN: details the cell-to-cell communication that must occur in
//                order to complete the ShareCoordinates function in the
//                CellManip library.                                
//-----------------------------------------------------------------------------
struct AtomSharePlan {
  int nsend;       // The number of processes to which this process will
                   //   send information, including itself (this process will
                   //   send ncomm-1 messages over MPI while executing this
                   //   plan)
  int nrecv;       // The number of receives to expect / post
  int* slist;      // The list of all cells controlled by this process (the
                   //   "send from" list)
  int* rlist;      // The list of all cell not controlled by this process but
                   //   which send information to this process (the "receive
                   //   from" list)
  pcgrp* send;     // Detailed list of sends from this process
  pcgrp* recv;     // List of receives for this process
  pcgrp selfrecv;  // The special "self receive" of cells that this plan
                   //   delivered to the same process
};
typedef struct AtomSharePlan ashr;

//-----------------------------------------------------------------------------
// DIRECTCOMMPLAN: holds all plans needed for direct-space communication by
//                 this process.                                     
//-----------------------------------------------------------------------------
struct DirectCommPlan {
  ashr mvshr[3];   // Atom sharing communication
  ashr frcmg[3];   // Force merger communication
};
typedef struct DirectCommPlan dcplan;

//-----------------------------------------------------------------------------
// GRIDSPLICE: details the information present in a grid struct imported from
//             another process.                                    
//-----------------------------------------------------------------------------
struct GridSplice {
  int npag;        // The number of pages controlled, one X value and all Y and
                   //   Z coordinates at that X value
  int ncol;        // The number of pencils / columns controlled, one X and one
                   //   Y value and all Z values
  int npc;         // The number of partial columns controlled, one X and one Y
                   //   value, Z values from Zo to Zf
  int* pagl;       // The list of pages controlled by a process
  int* coll;       // The list of columns controlled by a process, format
                   //   [X1][Y1][X2][Y2] ...
  int* pcl;        // The list of column pieces controlled by a process, format
                   //   [X1][Y1][Z1o][Z1f][X2][Y2][Z2o][Z2f] ...
};
typedef struct GridSplice gsplc;

struct AtomInCell {
  int id;
  int lj;
  double q;
  double loc[3];
  double frc[3];
};
typedef struct AtomInCell atomc;

typedef struct AtomBuffer atomb;

typedef struct AtomBufferPlusVelocity atombv;

typedef struct AtomBufferPlusAllInfo atombx;

typedef struct RangeBuffer rngbuff;

struct CellBlock {
  int nexp;            // The number of atoms the cell is exporting
  int nimp;            // The number of atoms the cell has imported
  int maxatom;         // The maximum number of atoms that this cell can hold
  int CGRank;          // The rank of the process responsible for this cell
                       //   within the cell grid private communicator
  int pmordr[3];       // Reciprocal space particle -> mesh interpolation order
                       //   (stored in the cell for convenience; otherwise it
                       //   becomes a headache to keep referencing that
                       //   information from the reccon struct)
  int gbin[4];         // The location of the cell in the cell grid (X, Y, and
                       //   Z position followed by absolute number in the
                       //   range 0... nX*nY*nZ)
  int nFscr;           // The number of scoring grids maintained by this cell
                       //   (initialized to zero by default, and only set to
                       //   nonzero values by restraint module functions)
  int* nr;             // Total number of atoms in each sector
  int* nsr;            // Total number of atoms in each sub-sector list
  int* ljIDbuff;       // Buffer array for CELL ID numbers of atoms that have
                       //   passed all the necessary cutoff tests; corresponds
                       //   to r2 values stored in the ljr2buff array
  int* qIDbuff;        // Buffer array for CELL ID numbers of atoms that have
                       //   passed the electrostatic cutoff test; corresponds
                       //   to r2 values stored in the qr2buff array
  int* GPSptr;         // Pointer into parent cell grid's AtmGPS 2D array
  rngbuff* ljr2buff;   // Buffer array of r2 values for LJ interactions
  rngbuff* qr2buff;    // Buffer array of r2 values for Q-Q interactions
  imat ordr;           // Ordering array for direct space interactions
  imat supordr;        // "Super-ordering" array for direct space interactions,
                       //   holds additional sorting after ordr is computed
  double orig[3];      // Coordinate origin of the cell's primary sector
  double midp[3];      // Coordinate origin of the cell's midpoint

  // Auxiliary arrays to store data that is not always needed
  bcof* xcof;           // xcof, ycof, and zcof store B-spline coefficients for
  bcof* ycof;           //   mapping atoms in this cell to a mesh.
  bcof* zcof;           //
  atomc* atmscr;        // Scratch space for merge sort

  // Atom information used during nonbonded loop.  The elements of the 
  // map array point to regions of data, similar to a matrix structure.
  atomc* data;
  atomc** map;

  // Restraint grid information; a series of 
  fbook* Fscr;

  // Buffer arrays
  atomb* import;
  atomb* pexport;
  atombv* Vimport;
  atombv* Vexport;
  atombx* Ximport;
  atombx* Xexport;
};
typedef struct CellBlock cell;

struct CellBlockGrid {
  int nsend;            // The maximum number of sends or receives for which
  int nrecv;            //   this cell grid is prepared
  int MyCellCount;      // Number of cells controlled by a particular process
  int sysID;            // The number of the system to which this grid pertains
  int tid;              // The rank of this particular process in the cell
                        //   grid's private communicator  
  int nthreads;         // The number of (MPI) threads assigned to this cell
                        //   grid, or 1 if no MPI
  int MasterHalfLoad;   // Half the master process's cell load
  int* nexp;            // The number of atoms that this cell grid holds in its
                        //   its export buffers
  int* maxexp;          // The maximum number of atoms that may be sent or
  int* maximp;          //   received in any of the expected messages
  int* MyCellDomain;    // Cells controlled by this particular process
  int* CrdPoolSize;     // 
  atomb** import;       // Pooled import and export buffers for all cells,
  atomb** pexport;       //   process-to-process communication involving atomb,
  atombv** Vimport;     //   atombv, and atombx types, similar to their
  atombv** Vexport;     //   counterparts in the individual cell structs
  atombx** Ximport;     //
  atombx** Xexport;     //
  atombx** CrdPool;     // Buffer for pooling coordinates to the master process
  dcplan DirCommPlan;   // Direct space communication plan
  int maxatom;          // The maximum number of atoms that any one cell may
                        //   hold
  int ncell;            // The number of cells in this grid
  int ng[3];            // The dimensions of the cell block grid
  double dbng[3];       // Dimensions of the cell block grid recast as doubles
  double celldim[7];    // celldim contains the lengths, in Angstroms, of the
                        //   edges of each cell, then the length of the
                        //   maximum direct space cutoff (Mcut copied over from
                        //   a dircon struct found in pmeDirectDS.h), and
                        //   finally (in the last three elements) the location
                        //   of the origin of the cell's central region in
                        //   units of the cell lengths
  cell* data;           // The linear array of cells in this grid
  cell*** map;          // Map to cells by x / y / z indices



  gsplc* MeshCommPlan;  // Mesh communication plans (allocates nthreads
                        //   on the master process, 1 on every other
                        //   process)
  imat AtmGPS;          // Lists of cell data indices where atoms ordered
                        //   according to the master topology file can be found
};
typedef struct CellBlockGrid cellgrid;

//-----------------------------------------------------------------------------
// BNDANGDIHE: a group of three integers, bonds, angles, and dihedrals. 
//-----------------------------------------------------------------------------
struct bndangdihe {
  int nbond;
  int nangl;
  int ndihe;
};
typedef struct bndangdihe bah;

typedef struct bondidx bond;

typedef struct bondcommand bondcomm;

typedef struct atombondlist bondlist;

typedef struct BondDef bonddef;

typedef struct anglidx angle;

typedef struct anglecommand anglcomm;

typedef struct atomanglelist angllist;

typedef struct AngleDef angldef;

typedef struct diheidx dihedral;

typedef struct dihedralcommand dihecomm;

typedef struct dihedrallist dihelist;

typedef struct DihedralDef dihedef;

typedef struct ConstraintCommand cnstcomm;

//-----------------------------------------------------------------------------
// SETTLEPARAMETERS: three parameters that can be conveniently pre-computed
//                   for SETTLE rigid water constraints.       
//-----------------------------------------------------------------------------
struct SettleParameters {
  double ra;
  double rb;
  double rc;
  double mO;
  double mH;
};
typedef struct SettleParameters settleparm;

typedef struct ExtraPoint expt;

typedef struct ExtraPointRule eprule;

typedef struct ConnectorMatrix map1234;

typedef struct ListOfGroup lgrp;

typedef struct NB14Pair nixpr;

typedef struct AuxiliaryElimination auxelim;

//-----------------------------------------------------------------------------
// AMBPRMTOP: a structure for AMBER 7 topology file data.               
//-----------------------------------------------------------------------------
struct ambprmtop {

  // Integers
  int natom;    // The number of atoms
  int nres;     // The number of residues
  int ntypes;   // Total number of distinct atom types
  int nhparm;   // Currently not used
  int nparm;    // Currently not used
  int tnexcl;   // The number of excluded atoms
  int natyp;    // The number of atom types in the parameter file (see Solty)
  int nphb;     // The number of distinct 10-12 interaction types
  int ifpert;   // Set to 1 if perturbation information is to be read
  int ifbox;    // Set to 1 if standard periodic box, 2 if truncated octahedron
  int nmxrs;    // Number of atoms in the largest residue
  int ifcap;    // Set to 1 if the CAP option from edit was specified
  int iptres;   // The final residues that is considered part of the solute
  int iptatom;  // The final atom that is considered part of the solute
  int nspm;     // The total number of molecules
  int nspsol;   // The first solvent molecule
  int natcap;   // The last atom before the start of cap waters placed by edit
  int blank;    // Not sure what this is
  int rattle;   // Flag to activate RATTLE
  int settle;   // Flag to activate SETTLE
  int nwat;     // The number of rigid three point water molecules
  int ncnst;    // The number of constraints in the system
  int ndf;      // The number of degrees of freedom in the system
  int nprtcl;   // The number of particles in the system (each constrained
                //   group of atoms counts as one particle, and each atom not
                //   in a constrained group, counts as one particle)
  int ljbuck;   // Flag to signal whether a Lennard-Jones 12-6 or Buckingham
                //   exp-6 potential is to be used
  int qshape;   // Flag to signal whether point charges or spherical Gaussian-
                //   distributed charges are to be used
  int neprule;  // The number of extra point rules, obtained from an auxiliary
                //   file (the name of the file is stored in the eprulesource
                //   field)
  int nxtrapt;  // The number of extra points, as counted by the extra point
                //   rule implementation routines
  int maxtrapt; // The maximum number of extra points that can be stored; kept
                //   to indicate when the extra points array must be expanded
  int numextra; // The number of extra points, as stated in the preamble to the
                //   topology file
  int ncopy;    // The number of Path Integral Molecular Dynamics slices /
                //   number of beads
  int ngrp;     // The number of bonded atom groups

  int EPInserted;    // Flag to indicate that extra points have been added to
                     //   this topology
  int norigatom;     // The original number of atoms, relevant if new extra
                     //   points are inserted
  int nclingrule;    // The number of cling rules defining how atoms react to
                     //   a grid-based restraint potential
  int RattleGrpMax;  // The maximum number of atoms that could be in any of the
                     //   RATTLE-constrained groups in this topology
  int AtomsNoElec;   // Flags to indicate that this topology contains atoms
                     //   which have no electrostatic or Lennard-Jones
  int AtomsNoLJ;     //   properties, respectively
  int ExclMarked;    // Flag to indicate that exclusions must be eliminated
                     //   from non-bonded interactions apart from what bond,
                     //   angle, and dihedral routines would already exclude

  // Reals
  double cutcap;    // The distance from the center of the cap to the outside
  double xcap;      // X coordinate for the center of the cap
  double ycap;      // Y coordinate for the center of the cap
  double zcap;      // Z coordinate for the center of the cap
  double lj14fac;   // The van-der Waals 1:4 interaction adjustment factor
  double elec14fac; // The electrostatic 1:4 interaction adjustment factor
  double TotalMass; // The total mass of the system (atomic mass units) 
  double initq;     // The initial total charge on the system (for output
                    //   record keeping--system charges are typically rounded
                    //   to the nearest integer, often zero, by adding a small
                    //   amount of charge to every atom)

  // In what follows, (b, a, h) means "(bonds, angles, and dihedrals)"
  bah withH;    // Number of (b, a, h) with hydrogen
  bah woH;      // Number of (b, a, h) without hydrogen
  bah woHC;     // Number of (b, a, h) without hydrogen plus constraints
  bah nBAH;     // Number of unique (b, a, h) types
  bah pert;     // Number of (b, a, h) to be perturbed
  bah wpert;    // Number of (b, a, h) with all atoms in the perturbed group

  // Arrays of elemental types
  int* LJIdx;        // The atom type index for each atom
  int* NExcl;        // The number of excluded atoms for each atom
  int* ConExcl;      // Start positions in the excluded atom list
  int* ExclList;     // The excluded atom list
  int* NBParmIdx;    // Nonbonded parameter index, ntypes*ntypes elements
  int* ResLims;      // The residue limits
  int* Join;         // Tree joining information, used by ancient programs
  int* Rotat;        // No longer in use
  int* Nsp;          // The total number of atoms in each molecule
  int* OldAtomNum;   // Old atom number sequence (initially not allocated, only
                     //   allocated if EPInserted is 1); this array stores the
                     //   number of the atom in the topology as initially read
                     //   from a file, before any extra points have been
                     //   inserted, or -1 for extra points that were not part
                     //   of the original topology
  int* MobileAtoms;  // Identifies mobile atoms in the system
  double* Charges;   // The atomic charges
  double* Masses;    // The atomic masses (atomic mass units)
  double* InvMasses; // The square root inverse atomic masses
  double* BondK;     // The bond force constants (kcal/mol-A^2)
  double* BondEq;    // The equilibrium bond lengths (Angstroms)
  double* AnglK;     // The angle force constants (kcal/mol-rad^2)
  double* AnglEq;    // The equilibrium angle values (radians)
  double* DiheK;     // The dihedral force constants (kcal/mol)
  double* DiheN;     // The dihedral periodicities
  double* DihePhi;   // The dihedral phase
  double* scee;      // The 1-4 electrostatic adjustments for all dihedrals
  double* scnb;      // The 1-4 van-der Waals adjustments for all dihedrals
  double* LJA;       // The Lennard-Jones A parameters
  double* LJB;       // The Lennard-Jones B parameters
  double* LJC;       // The "Lennard-Jones" C parameters (if these are given,
                     //   it means that the van-der Waals model is a modified
                     //   Buckingham potential, and thus the Lennard-Jones A
                     //   and B parameters are really Buckingham A and B
                     //   parameters)
  double* lVDWc;     // Long-Range van-der Waals corrections, stored by type,
                     //   to be divided by the instantaneous unit cell volume
                     //   when implemented
  double* SolA;      // The HBond A parameters
  double* SolB;      // The HBond B parameters
  double* HBCut;     // No longer in use
  double* Radii;     // Atomic radii
  double* Screen;    // Atomic screening terms (?)
  char RadSet[128];  // Radius set name
  char vstamp[128];  // Unknown         
  char WaterName[8]; // The name of water molecules (for implementing SETTLE)
  char* AtomNames;   // The atom names
  char* ResNames;    // The residue names
  char* AtomTypes;   // The atom type names
  char* TreeSymbols; // The atom tree chain symbols
  char* rattlemask;  // Mask to specially request RATTLE bond constraints
  char* norattlemask;// Mask to specially omit RATTLE bond constraints
  double gdim[6];    // The box dimensions

  // Arrays of structures
  bond* BIncH;       // Bonds including Hydrogen
  bond* BNoH;        // Bonds not including Hydrogen
  angle* AIncH;      // Angles including Hydrogen
  angle* ANoH;       // Angles not including Hydrogen
  dihedral* HIncH;   // Dihedrals including Hydrogen
  dihedral* HNoH;    // Dihedrals not including Hydrogen
  bondlist* BLC;     // The list of bonds controlled by each atom
  angllist* ALC;     // The list of angles controlled by each atom
  dihelist* HLC;     // The list of dihedral angles controlled by each atom
  bonddef* BParam;   // The re-organized, unified bond parameter array
  angldef* AParam;   // The re-organized, unified angle parameter array
  dihedef* HParam;   // The re-organized, unified dihedral parameter array
  cnstcomm* SHL;     // Constraint commands for each individual atom (SETTLE is
                     //   performed starting with the oxygens)
  eprule* eprules;   // List of rules for extra points
  expt* xtrapts;     // List of extra points (field nxtrapts holds the size)
  lgrp* lgrps;       // List of atom groups in this topology, defining
                     //   separate molecules within the simulation
  map1234* nb1234;   // A connectivity array detailing 1:1, 1:2, 1:3, and 1:4
                     //   interactions which is reflexive and symmetric
  auxelim* ElimPair; // A list of 1:1, 1:2, and 1:3 bonded or 1:4 nonbonded
                     //   interactions not involving dihedrals that must be
                     //   calculated based on added extra points
  lgrp* FR1Idx;      // Index into extra point array xtrapts (each atom has 
                     //   is own element of FR1Idx, but it points to no 
                     //   elements of the xtrapts array unless the atom is
                     //   frame atom 1 of one or more extra points).  When
                     //   atoms go into cells, extra points are ultimately
                     //   assigned after their frame atoms, and frame atom 1,
                     //   to which the extra point is attached, is the first
                     //   thing to look for.

  // Lennard-Jones pre-computed parameter tables
  dmat LJftab;                 // Forces
  dmat LJutab;                 // Energies

  // SETTLE pre-computed parameters
  settleparm FWtab;            

  // Source information
  char source[MAXNAME];        // Topology source file
  char eprulesource[MAXNAME];  // Extra points rules source file
};
typedef struct ambprmtop prmtop;

typedef struct UniqueBondList ublist;

//-----------------------------------------------------------------------------
// PrmtopCorrespondence: stores the correspondence between two topology
//                       structs, in terms of the atom-to-atom lineup and 
//                       global matching style.                     
//-----------------------------------------------------------------------------
struct PrmtopCorrespondence {
  int relate;           // How these topologies correspond generally...
                        //   0: atom N of topology A corresponds to atom N of
                        //      topology B
                        //   1: if atoms N and N+1 of topology A correspond to
                        //      anything in topology B, they will correspond to
                        //      atoms K and K+P of topology B, where P > 0
                        //   2: atoms of topology A may or may not correspond
                        //      atoms of topology B, but there is no guarantee
                        //      of ordering
  int vdw;              // Flag to indicate whether any van-der Waals
                        //   interactions change between the corresponding
                        //   atoms of topologies A and B
  int vdw14;            // Flag to indicate whether any 1:4 van-der Waals
                        //   interactions change between the corresponding
                        //   atoms of topologies A and B
  int elec;             // Flag to indicate whether any charges change between
                        //   the corresponding atoms of topologies A and B
  int elec14;           // Flag to indicate whether any 1:4 electrostatics
                        //   change between corresponding atoms of topologies
                        //   A and B
  int bond;             // Flag to indicate whether any bond terms change
                        //   between corresponding atoms of topologies A and B
                        //   (also set to 1 if RATTLE is toggled)
  int angl;             // Flag to indicate whether any angle terms change
                        //   between corresponding atoms of topologies A and B
  int dihe;             // Flag to indicate whether any dihedral terms change
                        //   between corresponding atoms of topologies A and B
  int appear;           // Flag to indicate that topology B has atoms that
                        //   have appeared amidst the contents of topology A
  int disappear;        // Flag to indicate that topology B has atoms that
                        //   have appeared from the contents of topology A
  int nxdf;             // The number of degrees of freedom in the extended
                        //   system which is the union of particles in A and B
  int nxprtcl;          // The number of particles in the extended system  
  int uniA;             // Number of unique atoms in system A
  int uniB;             // Number of unique atoms in system B
  int comAB;            // Number of common atoms between systems A and B
  int* matchA;          // The correspondences of atoms in topology A to B
  int* matchB;          // The correspondences of atoms in topology B to A
  int* corrA;           // The correspondences of atoms in topology A to B
  int* corrB;           // The correspondences of atoms in topology B to A
  double* dQ;
  ublist* SubBondA;     // Bonds that must be subtracted from topology A to
                        //   get topology B (bonds unique to topology A)
  ublist* AddBondB;     // Bonds unique to topology B
  ublist* SubAnglA;     // Angles that must be subtracted from topology A to
                        //   get topology B (angles unique to topology A)
  ublist* AddAnglB;     // Angles unique to topology B
  ublist* SubDiheA;     // Dihedrals that must be subtracted from topology A to
                        //   get topology B (dihedrals unique to topology A)
  ublist* AddDiheB;     // Dihedrals unique to topology B
  prmtop *tpA;          // Pointers to topologies A and B are stored for
  prmtop *tpB;          //   reference, so things like the number of atoms in
                        //   each topology is known
};
typedef struct PrmtopCorrespondence prmcorr;

typedef struct RestraintData nail;

typedef struct CubicSpline CSpln;

//-----------------------------------------------------------------------------
// ForceTable: this structure holds a force (and energy) piecewise cubic
//             spline table.                                            
//-----------------------------------------------------------------------------
struct ForceTable {

  // Attributes needed for force and energy evaluation
  int nbin;
  double rmax;
  double dr;
  double ivdr;
  CSpln* SD;
  CSpln* dSD;

  // Error analysis
  int FitType;            // The type of spline interpolation
                          //   0: best fit at intervals 0, 1/3, 2/3, 1, ...
                          //   1: fit with continuous derivatives at 0, 1, ...
  double fmaxerr;         // Maximum absolute error in the derivative (force)
  double fmaxerrloc;      // Location of the maximum error in the derivative
  double fmaxrelerr;      // Maximum relative error in the derivative
  double fmaxrelerrloc;   // Location of maximum relative error...
  double umaxerr;         // Maximum absolute error in the function (potential)
  double umaxerrloc;      // Location of the maximum error in the function
  double umaxrelerr;      // Maximum relative error in the function
  double umaxrelerrloc;   // Location of maximum relative error in the function
};
typedef struct ForceTable FrcTab;

struct Coordinates {
  int natom;        // The number of atoms in the system
  int isortho;      // Flag to indicate whether the unit cell is orthorhombic
  int* atmid;       // The ID numbers of atoms in the system; starts out as a
                    //   series of numbers from 0 to natom-1, but may be useful
                    //   for tracking atoms if coord structs are copied and
                    //   rearranged
  double* loc;      // The (current) locations of the atoms
  double* prvloc;   // The previous locations of the atoms
  double* scrloc;   // Scratch space for locations of the atoms
  double* vel;      // The (current) velocities of the atoms
  double* prvvel;   // The previous velocities of the atoms
  double* frc;      // The (current) forces on atoms
  double* prvfrc;   // The previous forces on atoms
  double* scrfrc;   // Scratch array of forces on atoms
  double gdim[6];   // The simulation cell dimensions
  double hgdim[3];  // The simulation cell vector half lengths, stored for
                    //   convenience but crucial to update
  dmat U;           // The transformation matrix for taking real-space
                    //   coordinates into fractional coordinate space
  dmat invU;        // The inverse transformation matrix, invU = U^-1
  dmat fcnorm;      // Normals to each of the pair interaction interfaces; rows
                    //   of this matrix contain the unit normal vectors for
                    //   easy memory access
};
typedef struct Coordinates coord;

typedef struct GridToGridMap g2gmap;

typedef struct EquivAtomGroup eagrp;

typedef struct ExtendedAtomDef xatomdef;

typedef struct ExtendedBondDef xbonddef;

typedef struct ExtendedHBondDef xhb1012def;

typedef struct ExtendedAngleDef xangldef;

typedef struct TorsionTerm torterm;

typedef struct BondIndex bidx;

typedef struct AnglIndex aidx;

typedef struct TorsionIndex hidx;

struct BondMap {
  int nbond;        // The number of bonds in this system
  bidx* id;         // Numbers of atoms participating in each bond, and indices
                    //   into the master list of bonds
  double* val;      // Values of the underlying bond length in the fitting set
  double* Ukernel;  // The contribution to the fitting matrix
  double* Ucontrib; // The kernel x stiffness constant = contribution to energy
  double* Fkernel;  // Same as Ukernel and Ucontrib above, but for forces.
  double* Fcontrib; //   These are significantly larger arrays, giving x, y,
                    //   and z contributions for each atom in the bond
};
typedef struct BondMap bondmap;

struct AngleMap {
  int nangl;        // The number of angles in this system
  aidx* id;         // Indices into the master list of angles
  double* val;      // Values of the underlying angle in the fitting set
  double* Ukernel;  // The contribution to the fitting matrix
  double* Ucontrib; // The kernel x stiffness constant = contribution to energy
  double* Fkernel;  // Same as Ukernel and Ucontrib above, but for forces.
  double* Fcontrib; //   These are significantly larger arrays, giving x, y,
                    //   and z contributions for each atom in the bond
};
typedef struct AngleMap anglmap;

struct TorsionMap {
  int ntterm;       // The number of torsional terms in this system
  hidx* id;         // Indices into the master list of torsion terms
  double* val;      // Values of the underlying angle in the fitting set
  double* Ukernel;  // The contribution to the fitting matrix
  double* Ucontrib; // The kernel x stiffness constant = contribution to energy
  double* Fkernel;  // Same as Ukernel and Ucontrib above, but for forces.
  double* Fcontrib; //   These are significantly larger arrays, giving x, y,
                    //   and z contributions for each atom in the bond
};
typedef struct TorsionMap tormap;

typedef struct AtomTypeSwitch typeswitch;

typedef struct AtomTypeBranch typebranch;

typedef struct MMSystem mmsys;

struct pmeDirectControlData {
  int LRvdw;         // Flag to activate long-ranged van-der Waals correction
  double Ecut;       // The electrostatic non-bonded cutoff
  double Vcut;       // The van-der Waals non-bonded cutoff
  double Mcut;       // The maximum of Ecut and Vcut
  double MaxDens;    // The maximum expected density of atoms (determines the
                     //   storage size in direct space decomposition cells)
  double invMcut;    // The inverse of Mcut
  double invEcut;    // The inverse of Ecut
  double ewcoeff;    // The Ewald coefficient, 0.5/sigma
  double sigma;      // The Gaussian width for spreading charges in preparation
                     //   for Ewald calculations
  double Dtol;       // The direct sum tolerance, the point at which the
                     //   difference in the interactions of Gaussian and point
                     //   charges is so small as to be deemed negligible
  double lkpspc;     // The spacing of the direct-space lookup table
};
typedef struct pmeDirectControlData dircon;

struct pmeRecipControlData {

  // Smooth particle Mesh Ewald
  int ordr[3];
  int* ng;
  double S;

  // Multi-Level (Smooth Particle Mesh) Ewald
  int nlev;
  int nslab;
  int nstrip;
  int ggordr;
  int PadYZ[4];
  double cfac[4];
  g2gmap* SPrv;
  g2gmap* SPcv;
  dbook* Urec;
  dbook* QL;
  fftw_plan* forwplan;
  fftw_plan* backplan;
};
typedef struct pmeRecipControlData reccon;

struct BCMeshKit {
  int plans;             // Flag to indicate that FFT plans have been made
                         //   especially for this struct (the plans may be
                         //   borrowed from elsewhere)
  double SelfEcorr;      // Self energy correction for this mesh setup
  double* Bx;            // B mesh prefactors in X, Y, and Z
  double* By;            //
  double* Bz;            //
  double* mx;            // M value (reciprocal space vector index) in X, Y,
  double* my;            //   and Z
  double* mz;            //
  double* mxs;           // Shifted M values in X, Y, or Z
  double* mys;           //
  double* mzs;           //
  double* Ex;            // Exponential tables for X, Y, and Z (used for both
  double* Ey;            //   orthorhombic and non-orthorhombic unit cells)
  double* Ez;            //
  fftw_plan forwplan;    // Forward FFT plan
  fftw_plan backplan;    // Backward FFT plan
};
typedef struct BCMeshKit bckit;

typedef struct ElectronDensityGaussian edgauss;

struct RestraintControls {
  int active;          // Flag to indicate that any restraints are active
  int usegrid;         // Flag to indicate a grid-based restraint is in use
  int usebelly;        // Flag to activate a belly mask and freeze all other
                       //   atoms
  int XpandGrid;       // Number to indicate what grid promotion should be used
                       //   in order to get smoother lookup tables (default of
                       //   1 means no expansion, expand by a factor of 1)
  double GridScale;    // Scaling factor to make the grid on file uniformly
                       //   stronger or weaker at run time
  char* GridFile;      // Name of the grid restraint file
  char* GridDefsFile;  // Name of file describing how atoms of the simulation
                       //   interact with the grid
  char* BellyMask;     // The belly mask string (ambmask format)
  char* FrozenMask;    // The belly mask string (ambmask format)
  fbook Rgrd;          // The restraint grid
};
typedef struct RestraintControls rstrcon;

struct ExecutionControl {
  struct timeval t0;
  struct timeval tC;
  struct timeval tti;
  struct timeval ttf;
  double bonds;
  double cellcomm;
  double nbInt;
  double nbDirAll;
  double nbBsp;
  double nbPtM;
  double nbFFT;
  double nbCnv;
  double nbMtP;
  double nbMtM;
  double nbRecAll;
  double Setup;
  double Integ;
  double Write;
  double Constraints;
  double Barostat;
  double Thermostat;

};
typedef struct ExecutionControl execon;

//-----------------------------------------------------------------------------
// CoupledHoover: a structure for storing all the necessary parameters to
//                manage a coupled Hoover thermo-barostat.           
//-----------------------------------------------------------------------------
struct CoupledHoover {
  double pmass;     // The mass of the barostat
  double qmass;     // The mass of the thermostat
  double chi;       // The friction coefficient of the thermostat
  double eta;       // The friction coefficient of the barostat
  double sigma;     // Computed from the number of degrees of freedom
  double TauT;      // Time constant for temperature fluctuations
  double TauP;      // Time constant for pressure fluctuations
};
typedef struct CoupledHoover choov;

//-----------------------------------------------------------------------------
// LangevinThermostat: a structure for storing constants related to a Langevin
//                     thermostat.                             
//-----------------------------------------------------------------------------
struct LangevinThermostat {
  double gamma_ln;
  double c_implic;
  double c_explic;
  double sdfac;
};
typedef struct LangevinThermostat lnbath;
  
//-----------------------------------------------------------------------------
// TrajectoryControlData: a structure for storing all parameters for managing  
//                        a trajectory.                         
//-----------------------------------------------------------------------------
struct TrajectoryControlData {

  // Molecular dynamics / minimization run parameters
  int mode;                // The mode for this run (0 = MD, 1 = minimization,
                           //   2 = force and energy evaluation for a single
                           //   snapshot, 3 = charge fitting)
  int nsys;                // The number of systems in play (1 for standard
                           //   MD, minimization, or force evaluation, 2 for
                           //   thermodynamic integration, more for replica
                           //   exchange)
  int MySystemCount;       // The number of systems tended by this process
  int ntop;                // The number of topologies (1 for standard MD,
                           //   minimization, or force evaluation, 2 for
                           //   thermodynamic integration or replica exchange)
  long long int nstep;     // The total number of steps to perform
  long long int nfistep;   // The total number of steps to include in each
                           //   file; additional files may be written if
                           //   nstep > nfistep
  long long int currstep;  // The current step number
  int currfi;              // The current file number
  int RemoveMomentum;      // Remove net system momentum at this interval   
  int irest;               // Flag to tell whether run is a restart or begins
                           //   from initial coordiantes only (no velocities)
  int ntwr;                // The restart file writing frequency
  int ntwx;                // The coordinate trajectory file writing frequency
  int ntwv;                // The velocity trajectory file writing frequency
  int ntwf;                // The force trajectory file writing frequency
  int ntpr;                // The output diagnostics file writing frequency
  int ioutfm;              // Format for writing trajectory coordinates
  int OverwriteOutput;     // Flag to activate overwriting of trajectory and
                           //   other output information
  int Reckless;            // Flag to deactivate case inhibition; default 0, do
                           //   not deactivate inhibition of potentially buggy
                           //   or unphysical cases
  int igseed;              // The random number generator seed
  int SyncRNG;             // Flag to indicate that, in MPI implementations,
                           //   separate processes should (set to 1) or should
                           //   not (default, set to 0) be synchronized
  int MaxRattleIter;       // The maximum number of RATTLE iterations
  int topchk;              // Flag to activate topology checking
  int ntt;                 // Thermostat type (default 0, no thermostating; see
                           //   Manual.c for all options)
  int ntp;                 // Coordinate rescaling type (default 0, no
                           //   coordinate rescaling; set to 1 for isotropic
                           //   and 2 for anisotropic rescaling)
  int barostat;            // Choice of barostat (default 1 for Berendsen, 2
                           //   for Monte Carlo)
  int vrand;               // Interval for random velocity reset
  int MCBarostatFreq;      // The frequency at which the Monte-Carlo barostat
                           //   is applied.  By default, every 100 steps.
  int TI;                  // Flag to activate thermodynamic integration;
                           //   default 0 but set to 1 to activate T.I. along
                           //   a linear mixing path
  int mxorder;             // The order to which the mixing coefficient lambda
                           //   is raised in order to mix forces between the
                           //   initial and final states
  int nsynch;              // In TI calculations, coordinates of the the two
                           //   trajectories will be explicitly synchronized 
                           //   every nsynch steps.  Because forces between
                           //   corresponding particles are synchronized at
                           //   every step and the order of operations hitting
                           //   each particle location should be the same, the
                           //   coordinates should never differ in theory.
                           //   However, this occasional check is just to be on
                           //   the safe side.
  long int rndcon;         // The random number counter

  double mxA;              // The mixing factors for states A and B in
  double mxB;              //   thermodynamic integration calculations
  double dmxA;             // Derivatives of the mixing factors for states A
  double dmxB;             //   and B in thermodynamic integration calculations
  double starttime;        // The starting simulation time
  double currtime;         // The current simulation time
  double dt;               // The time step
  double rattletol;        // The RATTLE tolerance
  double Ttarget;          // The target temperature
  double Tinit;            // The initial temperature
  double BerendsenTCoupl;  // Time constant for Berendsen temperature coupling,
                           //   units of ps^-1, default 0.4 ps^-1
  double BerendsenPTime;   // Time constant for Berendsen pressure coupling,
                           //   units of ps^-1, default 1.0 ps^-1
  double BerendsenPCoupl;  // Pressure constant for Berendsen pressure
                           //   coupling, units of bar^-1, default 44.6 x 10^-6
  double MCBarostatFac[3]; // The rescaling constant multipliers for random
                           //   moves using a Monte-Carlo barostat.  Up to
                           //   three constants may be specified.  Initially,
                           //   only the first has a non-negative value and is
                           //   set to the default of 0.0001.  The default is
                           //   to perform isotropic rescaling, but specifying
                           //   non-negative values for other variables will
                           //   cause the system to rescale anisotropically in
                           //   the other dimensions.
  double mcdVmax;          // Monte-Carlo barostat volume rescaling margin
  double Ptarget;          // The target pressure, units of bar, default 1.0
  double lambda;           // The value of the mixing coefficient in
                           //   thermodynamic integration (TI) calculations
  double EMinStep0;        // The initial energy minimization step size
                           //   (default 0.01 Angstroms)
  double EMinStep;         // The instantaneous energy minimization step size
  double EMinTol;          // The energy minimization convergence criterion
  choov npth;              // Nose-Hoover thermobarostat parameters
  lnbath lnth;             // Langevin thermostat constants
  rstrcon Leash;           // Restraint control structure

  // Topology correlations
  prmcorr prc;             // This is the topology correlation map, created for
                           //   the purposes of running TI

  // Force report parameters
  int DMPcrd;        // Activate coordinate dumping
  int DMPbond;       // Activate bond force dumping
  int DMPangl;       // Activate angle force dumping
  int DMPdihe;       // Activate dihedral force dumping
  int DMPdelec;      // Activate direct sum electrostatic force dumping
  int DMPrelec;      // Activate reciprocal sum electrostatic force dumping
  int DMPvdw;        // Activate van-der Waals force dumping
  int DMPall;        // Activate summed (total) force dumping
  char DMPvar[32];   // Variable name for Matlab output

  // File names
  char inpname[MAXNAME];   // Input command file
  cmat ipcname;            // Input coordinates file(s)
  char dumpname[MAXNAME];  // The force dump file (for printing a comprehensive
                           //   report on all forces and energies)
  char rsrptname[MAXNAME]; // The residue report file (for printing a human-
                           //   readable description of all residues)
  char parmfile[MAXNAME];  // The force field parameter file (for fitting
                           //   torsions and other terms in a new model)
  char fmodfile[MAXNAME];  // The force field parameter file (for fitting
                           //   torsions and other terms in a new model)

  // File base names (these serve as the file names in case the
  // corresponding suffixes are NULL strings                   
  cmat rstbase;           // Restart
  cmat trjbase;           // Coordinates trajectory
  cmat velbase;           // Velocity trajectory
  cmat frcbase;           // Force trajectory
  char outbase[MAXNAME];  // Output diagnostics

  // File suffixes
  cmat rstsuff;           // Restart
  cmat trjsuff;           // Coordinates trajectory
  cmat velsuff;           // Velocity trajectory
  cmat frcsuff;           // Force trajectory
  char outsuff[32];       // Output diagnostics

  // Input file text, verbatim
  cmat inptext;
  char* inpline;

  // Parallel control data
  int tid;                   // Thread rank in MPI_COMM_WORLD, if MPI is
                             //   defined, or zero otherwise
  int nthreads;              // The total number of threads in the parallel run

  int nCPUcluster;           // Number of tightly connected CPU clusters
  lgrp* CPUcluster;          // Clusters of tightly connected CPUs
  imat SystemCPUs;           // List of CPUs devoted to each system, as ranked
                             //   in MPI_COMM_WORLD if MPI is defined or simply
                             //   zeros if the implementation is serial
  int* MySystemDomain;       // List of systems tended by this process
};
typedef struct TrajectoryControlData trajcon;

//-----------------------------------------------------------------------------
// EnergyTracker: a structure for tracking instantaneous energy and related
//                quantities such as pressure and temperature.  
//-----------------------------------------------------------------------------
struct EnergyTracker {

  // Components of energy
  double delec;     // Direct space electrostatic energy
  double relec;     // Reciprocal space electrostatic energy
  double vdw12;     // Total repulsive vdW energy
  double vdw6;      // Total attractive (London dispersion) vdW energy
  double bond;      // Energy due to bonded terms
  double angl;      // Energy due to angle terms
  double dihe;      // Energy due to dihedral terms
  double kine;      // Total kinetic energy
  double P;         // The current pressure
  double V;         // The current volume
  double T;         // The current temperature
  double dVdL;      // Derivative of the mixed potential energy function with
                    //   respect to the mixing parameter lambda

  // Sums of other energy quantities
  double elec;      // Total electrostatic energy
  double eptot;     // Total potential energy
  double etot;      // Total energy

  // Running averages
  double AVEdelec;
  double AVErelec;
  double AVEvdw;
  double AVEbond;
  double AVEangl;
  double AVEdihe;
  double AVEkine;
  double AVEelec;
  double AVEeptot;
  double AVEetot;
  double AVEP;
  double AVEV;
  double AVET;
  double AVEVir;
  double AVEdVdL;

  // Running standard deviations
  double RMSdelec;
  double RMSrelec;
  double RMSvdw;
  double RMSbond;
  double RMSangl;
  double RMSdihe;
  double RMSkine;
  double RMSelec;
  double RMSeptot;
  double RMSetot;
  double RMSP;
  double RMSV;
  double RMST;
  double RMSVir;
  double RMSdVdL;

  // Virial tensor
  double Vir[9];

  // Energy decomposition
  double* BondUdc;  // Quantities for recalculating bond energy

  // Directive flags
  int updateU;      // Flag to activate system energy update
  int updateV;      // Flag to activate system virial update
  int nUdc;         // Number of types in the bond energy decomposition
  int Esummed;      // Flag to indicate whether the energy has been
                    //   summed since the last initialization
};
typedef struct EnergyTracker Energy;

//-----------------------------------------------------------------------------
// uform: all the potential function structs wrapped up in one type.    
//-----------------------------------------------------------------------------
struct PotentialFunction {
  prmtop tp;        // Topology
  dircon dcinp;     // Direct-space controls
  FrcTab Etab;      // Direct-space standard (coarse) lookup table
  FrcTab EHtab;     // Direct-space high-resolution lookkup table
  reccon rcinp;     // Reciprocal space controls
  bckit PPk;        // Convolution support data
};
typedef struct PotentialFunction uform;

//-----------------------------------------------------------------------------
// mdsys: all structs for a molecular dynamics trajectory in one type.  
//-----------------------------------------------------------------------------
struct MolecularDynamicsSystem {
  coord crd;        // Coordinates (phone book organization of the system)
  cellgrid CG;      // Cell grid (spatial reorganization of the system)
  Energy sysUV;     // System energy decomposition and virial
  execon etimers;   // Timing data (for this system only)
};
typedef struct MolecularDynamicsSystem mdsys;


void FreeTopology(prmtop *tp);

coord ReadRst(prmtop *tp, char* source);

void DestroyCoord(coord *crd);

void FreeFrcTab(FrcTab *F);

cellgrid CreateCellGrid(coord *crd, dircon *dcinp, reccon *rcinp, prmtop *tp,
                        trajcon *tj, int sysnum);
#ifdef MPI
void LinkCellGrid(cellgrid *CG, coord *crd, reccon *rcinp);
#else
void LinkCellGrid(cellgrid *CG, reccon *rcinp);
#endif

void PrepPME(cellgrid *CG, reccon *rcinp, coord *crd);

bckit CreateBCKit(reccon *rcinp, dbook *Q, coord *crd, prmtop *tp, int DoPlan);

void InitHistory(coord *crd);

void InitBasicTrajcon(trajcon *tj);

void InitExecon(execon *tm);

uform InitPotential(char* topsrc, double NBcut, trajcon *tj);

void MMForceEnergy(uform *U, mdsys *MD, trajcon *tj);

void AtomsToCells(coord *crd, cellgrid *CG, prmtop *tp);

void ImageBondedGroups(coord *crd, prmtop *tp);

void DestroyTrajCon(trajcon *tj);

void DestroyCellGrid(cellgrid *CG);

void DestroyAdvancedRecCon(reccon *rcinp, cellgrid *CG);

void DestroyBCKit(bckit *PPk);

void DestroyEnergyTracker(Energy *sysUV);

//-----------------------------------------------------------------------------
// The following functions are not currently used by any program calling mdgx.
// They are included for completeness and possible future use.
//-----------------------------------------------------------------------------
// void GetPrmTop(prmtop *tp, trajcon *tj, int adjbnd);
// 
// void LongRangeVDW(prmtop *tp, dircon *dcinp);
// 
// FrcTab DirectSpaceR2(double range, double spc, double ewcoeff,
//                      int contderiv);
// 
// void DestroyRecCon(reccon *rcinp, cellgrid *CG);
// 
// void AtomForces(coord *crd, cellgrid *CG, prmtop *tp, dircon *dcinp,
//                  FrcTab *Etab, FrcTab *EHtab, reccon *rcinp, bckit *PPk,
//                  Energy *sysUV, execon *etimers, trajcon *tj);
// 
// double KineticEnergy(cellgrid *CG, coord *crd, prmtop *tp, trajcon *tj);
// 
// double SystemTemperature(cellgrid *CG, coord *crd, prmtop *tp,
//                          Energy *sysUV, trajcon *tj, int updateKE);
// 
// #ifdef MPI
// void SumTotalEnergy(cellgrid *CG, Energy *sysUV);
// #else
// void SumTotalEnergy(Energy *sysUV);
// #endif
// 
// void InitBasicDircon(dircon *dcinp, double NBcut);
// 
// void InitBasicReccon(reccon *rcinp, dircon *dcinp);
// 
// void PrimeBasicTopology(prmtop *tp);
// 
// void InitializeEnergy(Energy *sysUV, trajcon *tj, prmtop *tp,
//                        int allocBondUdc);
// 
// mdsys LoadCoordToGrid(char* crdname, uform *U, trajcon *tj);
// 
// void DefineMPITypes(trajcon *tj);
// 
// void FreeMPITypes(trajcon *tj);
// 
// void TransCrd(double* crds, int natom, double* tvec, double step);
// 
// void RotateCrd(double* crds, int natom, dmat U);
// 
// void FindCoordCenter(double* C, double* mass, int usem, int n,
//                      double* cofm);
// 
// void QuatAlign(double* frameI, double* frameII, int natom, double* mass,
//                int m, dmat *U);
// 
// imat CreateImat(int N, int M);
// 
// void DestroyImat(imat *A);
// 
// imat ReallocImat(imat *A, int M, int N);
// 
// dmat CreateDmat(int N, int M, int prepFFT);
// 
// void DestroyDmat(dmat *A);
// 
// dmat ReallocDmat(dmat *A, int M, int N);
// 
// void CopyDmat(dmat *Ac, dmat *A, int Acex);
// 
// void DestroyCmat(cmat *A);

#endif 
