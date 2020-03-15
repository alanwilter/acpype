#ifndef SFF_H
#define SFF_H

#include <stdio.h>

#ifndef	FALSE
#define	FALSE	0
#endif

#ifndef	TRUE
#define	TRUE	1
#endif

#define	UNDEF	(-1)

	/* Fundamental sff types:	*/

typedef	int	INT_T;
typedef	size_t	SIZE_T;

#define NAB_DOUBLE_PRECISION 1
#define REAL_T double

	/* other nab types:	*/
typedef	char	STRING_T;
typedef	FILE	FILE_T;

// The NAME parameters control how the mask parser expects strings to look.
// Originally, the parameter file assumes the atom, residue, symbol, etc. 
// names to be four characters. When stored as a string, the NULL character
// is required, requiring a size of 5. It has been increased to 6 in cpptraj
// to accomodate the slightly larger names that can be found in MOL2 files.
// NAMESIZE: Default size for atom and residue names, 5 + NULL.
// Amber atom/residue names are 4, but some mol2 atom types are larger.
#define NAMESIZE 6
#define NAME_DEFAULT "     "
// Default type for atom names, res names, atom types etc
typedef char NAME_T[NAMESIZE];

typedef struct parm {
	char	ititl[81];
	INT_T 	IfBox, Nmxrs, IfCap,
		 Natom,  Ntypes,  Nbonh,  Mbona,  Ntheth,  Mtheta, 
		 Nphih,  Mphia,  Nhparm, Nparm, Nnb, Nres,
		 Nbona,  Ntheta,  Nphia,  Numbnd,  Numang,  Nptra,
		 Natyp,  Nphb, Nat3, Ntype2d, Nttyp, Nspm, Iptres, Nspsol,
		 Ipatm, Natcap, Numextra;
	NAME_T *AtomNames, *ResNames, *AtomType, *AtomTree;
	REAL_T	*Charges, *Masses, *Rk, *Req, *Tk, *Teq, *Pk, *Pn, *Phase,
		 *Solty, *Cn1, *Cn2, *HB12, *HB10, *Rborn, *Fs;
	REAL_T	Box[4], Cutcap, Xcap, Ycap, Zcap, Fsmax;
	INT_T 	*Iac, *Iblo, *Cno, *Ipres, *ExclAt, *TreeJoin, 
		*AtomRes, *BondHAt1, *BondHAt2, *BondHNum, *BondAt1, *BondAt2, 
		*BondNum, *AngleHAt1, *AngleHAt2, *AngleHAt3, *AngleHNum, 
		*AngleAt1, *AngleAt2, *AngleAt3, *AngleNum, *DihHAt1, 
		*DihHAt2, *DihHAt3, *DihHAt4, *DihHNum, *DihAt1, *DihAt2, 
		*DihAt3, *DihAt4, *DihNum, *Boundary, *AtomicNum;
	INT_T	*N14pairs, *N14pairlist;
	REAL_T *Gvdw, *Scee, *Scnb, *N14sceelist, *N14scnblist;
	INT_T	Nstrand, Ncomplex;       /* HCP: number of strands/complexes */
	INT_T	*Ipstrand, *Ipcomplex;   /* HCP: index into residue/strand list */
} PARMSTRUCT_T;

typedef struct xmin_opt{ 
    INT_T  mol_struct_opt;
    INT_T  maxiter;
    REAL_T grms_tol;
    INT_T  method;
    INT_T  numdiff;
    INT_T  m_lbfgs;
    INT_T  iter;
    REAL_T xmin_time;
    INT_T  ls_method;
    INT_T  ls_maxiter;
    REAL_T ls_maxatmov;
    REAL_T beta_armijo;
    REAL_T c_armijo;
    REAL_T mu_armijo;
    REAL_T ftol_wolfe;
    REAL_T gtol_wolfe;
    INT_T  ls_iter;
    INT_T  print_level;
    INT_T  error_flag; 
} XMIN_OPT_T;

typedef struct lmod_opt {
    INT_T niter;
    INT_T nmod;
    INT_T kmod;
    INT_T nrotran_dof;
    INT_T nconf;
    REAL_T minim_grms;
    REAL_T energy_window;
    REAL_T conf_separation_rms;
    INT_T eig_recalc;
    INT_T ndim_arnoldi;
    INT_T lmod_restart;
    INT_T n_best_struct;
    INT_T mc_option;
    REAL_T rtemp;
    REAL_T lmod_step_size_min;
    REAL_T lmod_step_size_max;
    INT_T nof_lmod_steps;
    REAL_T lmod_relax_grms;
    INT_T nlig;
    INT_T apply_rigdock;
    INT_T nof_poses_to_try;
    INT_T random_seed;
    INT_T print_level;
    REAL_T lmod_time;
    REAL_T aux_time;
    INT_T error_flag;
} LMOD_OPT_T;

   /* the output for all non error emissions and some error ones too */
extern FILE	*nabout;

   /*  signatures for the public routines:  */

#ifdef __cplusplus
extern "C" {
#endif

INT_T		conjgrad( REAL_T*, INT_T*, REAL_T*,
			  REAL_T ( *func )( REAL_T*, REAL_T*, INT_T* ),
			  REAL_T*, REAL_T*, INT_T* );
REAL_T      gauss( REAL_T*, REAL_T* );
INT_T		getxv( STRING_T*, INT_T, REAL_T, REAL_T*, REAL_T* );
INT_T       getxyz( STRING_T**, INT_T*, REAL_T* );
REAL_T      lmod(INT_T*natm, REAL_T *xyz, REAL_T *grad, REAL_T *energy,
            REAL_T *conflib, REAL_T *lmod_trajectory,
            INT_T *lig_start, INT_T *lig_end, INT_T *lig_cent,
            REAL_T *tr_min, REAL_T *tr_max, REAL_T *rot_min, REAL_T *rot_max,
            struct xmin_opt *xo, struct lmod_opt *lo );
INT_T       lmod_opt_init( struct lmod_opt *lo, struct xmin_opt *xo );
INT_T		md( INT_T, INT_T, REAL_T*, REAL_T*, REAL_T*,
		    REAL_T ( *func )( REAL_T*, REAL_T*, INT_T* ) );
INT_T		mdrat( INT_T, INT_T, REAL_T*, REAL_T*, REAL_T*, REAL_T*, REAL_T ( *mme )() );
INT_T		mm_options( STRING_T* );
void            mm_set_checkpoint( STRING_T** );
REAL_T		mme( REAL_T*, REAL_T*, INT_T* );
REAL_T		mme2( REAL_T*, REAL_T*, REAL_T*, REAL_T*, 
		      REAL_T*, INT_T*, INT_T*, INT_T*,
		      INT_T *, INT_T *, INT_T *, INT_T *,
		      INT_T*, INT_T* , char* );
REAL_T		mme4( REAL_T*, REAL_T*, INT_T* );
REAL_T		mme_rattle( REAL_T*, REAL_T*, INT_T* );
INT_T		mme_init_sff( PARMSTRUCT_T*, INT_T*, INT_T*, REAL_T*, STRING_T* );
INT_T		newton( REAL_T*, INT_T*, REAL_T*,
			REAL_T ( *func1 )( REAL_T*, REAL_T*, INT_T* ),
			REAL_T ( *func2 )( REAL_T*, REAL_T*, REAL_T*, REAL_T*,
					   REAL_T*, INT_T*, INT_T*, INT_T*,
					   INT_T*, INT_T*, INT_T*, INT_T*,
					   INT_T*, INT_T*, char* ),
			REAL_T*, REAL_T*, INT_T* );
INT_T		nmode( REAL_T*, INT_T,
			REAL_T ( *func)( REAL_T*, REAL_T*, REAL_T*, REAL_T*,
					 REAL_T*, INT_T*, INT_T*, INT_T*,
					 INT_T*, INT_T*, INT_T*, INT_T*,
					 INT_T*, INT_T*, char*  ),
		       INT_T, INT_T, REAL_T, REAL_T, INT_T);
INT_T		putxv( STRING_T*, STRING_T*, INT_T, REAL_T, REAL_T*, REAL_T* );
INT_T           putxyz( STRING_T**, INT_T*, REAL_T* );
REAL_T		rand2( void );
INT_T		sasad( REAL_T*, REAL_T*, REAL_T*, INT_T, REAL_T );
INT_T   xmin_opt_init( struct xmin_opt* );
REAL_T	xmin( REAL_T ( *func )( REAL_T*, REAL_T*, INT_T* ),
             INT_T*, REAL_T*, REAL_T*, REAL_T*, REAL_T*, struct xmin_opt* );
REAL_T          xminC(INT_T*, INT_T*, INT_T*, REAL_T*, INT_T*, INT_T*,
                      INT_T*, REAL_T*, REAL_T*, REAL_T*, REAL_T*,
                      INT_T*, REAL_T*, INT_T*, INT_T*, INT_T*, INT_T*,
                      REAL_T*, REAL_T*, REAL_T*, REAL_T*,
                      REAL_T*, REAL_T*, INT_T*, INT_T*);
INT_T            mme_rism_max_memory();

   /*  prmtop routine public interfaces: */
INT_T    free_prm( PARMSTRUCT_T* prm );
PARMSTRUCT_T*  rdparm( STRING_T* name );
INT_T    wrparm( PARMSTRUCT_T* prm, STRING_T* name );

	/* chiral centers:	*/

typedef	struct	chiral_t	{
	INT_T	c_anum[ 4 ];
	REAL_T	c_dist;
} CHIRAL_T;
void	chirvol( int, int, int, int, int, REAL_T *, REAL_T *, REAL_T * );

	/* Distance geometry "bounds" matrix	*/

typedef	struct	bounds_t	{
	PARMSTRUCT_T	*b_prm;
	INT_T		b_nres;
	INT_T		b_natoms;
	INT_T		*b_aoff;
    void        **b_aindex;   /* was ATOM_T; just an opaque pointer for now */
	INT_T		*b_nconnect;
	INT_T		**b_connect;
	REAL_T		**b_bmat;
	INT_T		**b_boundList;
	REAL_T		*b_lowerBoundPathLength;
	REAL_T		*b_upperBoundPathLength;
	INT_T		b_schiral;
	INT_T		b_nchiral;
	CHIRAL_T	*b_chiral;
} BOUNDS_T;

  /*  rism routine interfaces: */
  /* 3D RISM section  */
  /*  N.B.: must match the rismprm_t struct in amber_rism_interface.F90 */
  typedef struct {
    REAL_T solvcut;
    REAL_T buffer;
    REAL_T grdspc[3];
    REAL_T solvbox[3];
    REAL_T mdiis_del;
    REAL_T mdiis_restart;
    REAL_T fcecut;
    REAL_T fcenormsw;
    REAL_T uccoeff[4];
    REAL_T biasPotential;
    REAL_T treeDCFMAC;
    REAL_T treeTCFMAC;
    REAL_T treeCoulombMAC;
    REAL_T asympKSpaceTolerance;
    REAL_T ljTolerance;
    REAL_T chargeSmear;
    INT_T closureOrder;
    INT_T ng3[3];
    INT_T rism;      /* non-zero if RISM is turned on */
    INT_T asympCorr;
    INT_T mdiis_nvec;
    INT_T mdiis_method;
    INT_T maxstep;
    INT_T npropagate;
    INT_T centering;
    INT_T zerofrc;
    INT_T apply_rism_force;
    INT_T polarDecomp;
    INT_T entropicDecomp;
    INT_T gfCorrection;
    INT_T pcplusCorrection;
    INT_T rismnrespa;
    INT_T fcestride;
    INT_T fcenbasis;
    INT_T fcenbase;
    INT_T fcecrd;
    INT_T fceweigh;
    INT_T fcetrans;
    INT_T fcesort;
    INT_T fceifreq;
    INT_T fcentfrcor;
    INT_T fcewrite;
    INT_T fceread;
    INT_T saveprogress;
    INT_T ntwrism;
    INT_T verbose;
    INT_T progress;
    INT_T write_thermo; 
    INT_T treeDCFOrder;
    INT_T treeTCFOrder;
    INT_T treeCoulombOrder;
    INT_T treeDCFN0;
    INT_T treeTCFN0;
    INT_T treeCoulombN0;
    INT_T selftest;
    INT_T treeDCF;
    INT_T treeTCF;
    INT_T treeCoulomb;
    /*This is an unused variable that aligns
      the type on eight byte boundaries*/
    //INT_T padding; 
  } RismData;

#ifdef RISMSFF
  void rism_force_( REAL_T*, REAL_T*, REAL_T*, INT_T*, INT_T* );
  void rism_setparam_( RismData*, INT_T*, REAL_T*,
                       INT_T*, INT_T*, STRING_T[10][8] ,
                       INT_T*, STRING_T*,  INT_T*, STRING_T*,  INT_T*, STRING_T*,
                       INT_T*, STRING_T*,  INT_T*, STRING_T*,  INT_T*, STRING_T*,
                       INT_T*, STRING_T*,  INT_T*, STRING_T*,  INT_T*, STRING_T*,
                       INT_T*, STRING_T*,  INT_T*, STRING_T*,  INT_T*, STRING_T*,
                       INT_T*, STRING_T*,  INT_T*, STRING_T*,  INT_T*, STRING_T*,
                       INT_T*, STRING_T*,  INT_T*, STRING_T*,  INT_T*, STRING_T*,
                       INT_T*, STRING_T*,  INT_T*, STRING_T*,  INT_T*, STRING_T*,
                       INT_T*, STRING_T*,  INT_T*, STRING_T*,  INT_T*, STRING_T*,
                       INT_T*, STRING_T*,
                       INT_T*, INT_T*, INT_T*,
                       REAL_T*, REAL_T*, REAL_T*, REAL_T*,
                       INT_T*, INT_T*);
  void rism_init_( INT_T*);
  void rism_list_param_();
  void rism_writesolvdist_(INT_T*);
  void rism_solvdist_thermo_calc_(INT_T*,INT_T*);
  void rism_thermo_print_(INT_T*,REAL_T*);
  void rism_printtimer_();
  void rism_max_memory_();
#endif /*RISMSFF*/

#ifdef __cplusplus
}
#endif

/*  from arpack:  */
void arsecond_( double * );

#ifndef INC_AMBER_NETCDF_H
#define INC_AMBER_NETCDF_H
/*! \file AmberNetcdf.h 
  * \author Daniel R. Roe
  * \date 2010-12-07
  * \brief A C implementation of routines for reading and writing the Amber 
  * Netcdf trajectory and restart formats.
  *
  * Based on Cpptraj implementation.
  * Original implementation of netcdf in Amber by Jon Mongan.
  */

// NOTE: It would be better to allocate memory for single-precision coords
//       upon loading netcdf traj, but since NAB does not really handle pointers
//       the memory must be allocated during every read/write.
// NOTE: Any changes made to the structure below must also be made to 
//       nab_netcdf.h, however due to nab not recognizing the double
//       type all double vars here must be float vars there. NAB float is
//       equivalent to double (see nab/defreal.h etc), but this is why there 
//       have to be two separate definitions of this structure.
/// Hold info for Amber Netcdf trajectory or restart
struct AmberNetcdf {
  double temp0;       // Temperature of current frame (if TempVID!=-1)
  double restartTime; // Simulation time if Amber restart
  int isNCrestart;    // 0 if trajectory, 1 if restart
  int ncid;           // Netcdf ID of the file when open
  int frameDID;       // ID of frame dimension
  int ncframe;        // Number of frames in the file
  int currentFrame;   // Current frame number 
  int atomDID;        // ID of atom dimension
  int ncatom;         // Number of atoms
  int ncatom3;        // Number of coordinates (ncatom * 3)
  int coordVID;       // ID of coordinates variable
  int velocityVID;    // ID of velocities variable
  int cellAngleVID;   // ID of box angle variable
  int cellLengthVID;  // ID of box length variable
  int spatialDID;
  int labelDID;
  int cell_spatialDID;
  int cell_angularDID;
  int spatialVID;
  int timeVID;
  int cell_spatialVID;
  int cell_angularVID;
  int TempVID;
};

// NOTE: To be NAB-useable, any functions added below must also be referenced in
//       nab/nabcode.h and nab/symbol.c
int netcdfDebug(struct AmberNetcdf*);
int netcdfLoad(struct AmberNetcdf*,char *);
int netcdfClose(struct AmberNetcdf*);
int netcdfWriteRestart(char *, int, double *, double *, double *, double, double);
int netcdfCreate(struct AmberNetcdf*,char *, int, int);
int netcdfGetVelocity(struct AmberNetcdf*, int, double *);
int netcdfGetFrame(struct AmberNetcdf*, int, double*, double*);
int netcdfGetNextFrame(struct AmberNetcdf*,double*,double*);
int netcdfWriteFrame(struct AmberNetcdf*, int, double *, double *);
int netcdfWriteNextFrame(struct AmberNetcdf*, double *, double *);
int netcdfInfo(struct AmberNetcdf*);

#endif

/// ptrajmask: The enhanced atom selection mask parser from ptraj.
/// Originally written by Viktor Hornak, Stony Brook University.
/// Adapted as standalone code by Dan Roe, NIST.
#ifdef __cplusplus
extern "C" {
#endif
// The NAME parameters control how the mask parser expects strings to look.
// Originally, the parameter file assumes the atom, residue, symbol, etc. 
// names to be four characters. When stored as a string, the NULL character
// is required, requiring a size of 5. It has been increased to 6 in cpptraj
// to accomodate the slightly larger names that can be found in MOL2 files.
// NAMESIZE: Default size for atom and residue names, 5 + NULL.
// Amber atom/residue names are 4, but some mol2 atom types are larger.
#define NAMESIZE 6
#define NAME_DEFAULT "     "
// Default type for atom names, res names, atom types etc
typedef char NAME[NAMESIZE];
// More defines
#define  MAXSELE 1000
#define  ALL      0
#define  NUMLIST  1
#define  NAMELIST 2
#define  TYPELIST 3
#define  ELEMLIST 4
/* parseMaskString()
 * The main interface to the mask parser. Takes a mask expression and some
 * information from a parameter file (stored in a PARMSTRUCT_T struct),
 * atomic coords in X0Y0Z0X1Y1Z1... format, and a debug value controlling how
 * much debug information is printed (internally the global int prnlev).
 * It returns a integer mask array imask[i]=1|0, i=0,atoms-1
 * which contains the resulting atom selection.
 */
int *parseMaskString(char*,PARMSTRUCT_T*,double*,int);
#ifdef __cplusplus
}
#endif
#endif
