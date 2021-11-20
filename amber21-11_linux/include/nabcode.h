#ifndef	NABCODE_H
#define	NABCODE_H

	/* "stub" types:	*/

#include "defreal.h"
#include "nabtypes.h"

typedef	char		HASH_T;
typedef	struct	curhash_t	{
	int	index;
	char	*pointer;
} CURHASH_T;

	/* nab builtins (but no libc or libm calls):	*/

#define		I2R(i)	((REAL_T)(i))
INT_T		addresidue( MOLECULE_T*, STRING_T*, RESIDUE_T* );
INT_T		addstrand( MOLECULE_T*, STRING_T* );
INT_T		alignframe( MOLECULE_T*, MOLECULE_T* );
INT_T		andbounds( BOUNDS_T*, MOLECULE_T*, STRING_T*, STRING_T*, REAL_T, REAL_T );
REAL_T		angle( MOLECULE_T*, STRING_T*, STRING_T*, STRING_T* );
REAL_T		anglep( POINT_T, POINT_T, POINT_T );
INT_T		axis2frame( MOLECULE_T*, POINT_T, POINT_T );
MOLECULE_T	*bdna( STRING_T** );
INT_T		bonded_atoms( ATOM_T*, ATOM_T** );
INT_T		cap( MOLECULE_T*, STRING_T*, INT_T, INT_T );
INT_T		circle( REAL_T*, REAL_T*, REAL_T*, REAL_T* );
INT_T		conjgrad( REAL_T*, INT_T*, REAL_T*,
			  REAL_T ( *func )( REAL_T*, REAL_T*, INT_T* ),
			  REAL_T*, REAL_T*, INT_T* );
INT_T		connectres( MOLECULE_T*, STRING_T*, INT_T, STRING_T*, INT_T, STRING_T* );
MOLECULE_T	*copymolecule( MOLECULE_T* );
INT_T		countmolatoms( MOLECULE_T*, STRING_T* );
INT_T		countmolres( MOLECULE_T*, STRING_T* );
INT_T		countmolstrands( MOLECULE_T*, STRING_T* );
INT_T		countstrandresidues( MOLECULE_T*, INT_T );
STRING_T	*date();
REAL_T		db_viol( REAL_T*, REAL_T*, INT_T* );
REAL_T		db_viol3( REAL_T*, REAL_T*, INT_T* );
MOLECULE_T	*dg_helix( STRING_T**, STRING_T**, STRING_T**, STRING_T**, STRING_T**, STRING_T**, 
		REAL_T*, REAL_T*, REAL_T*, REAL_T*, STRING_T** );
INT_T		dg_options( BOUNDS_T*, STRING_T* );
REAL_T		dist( MOLECULE_T*, STRING_T*, STRING_T* );
REAL_T		distp( POINT_T, POINT_T );
BOUNDS_T	*dt_to_bmat( MOLECULE_T**, STRING_T**, STRING_T** );
INT_T		dt_to_prmtop( STRING_T**, STRING_T**, STRING_T**, STRING_T**, STRING_T** );
INT_T		dumpatom( FILE_T*, RESIDUE_T*, INT_T, INT_T );
INT_T		dumpbounds( FILE_T*, BOUNDS_T*, INT_T );
REAL_T		dumpboundsviolations( FILE_T*, BOUNDS_T*, REAL_T );
REAL_T		dumpchiviolations( FILE_T*, BOUNDS_T*, REAL_T);
INT_T		dumpmatrix( FILE_T*, MATRIX_T );
INT_T		dumpmolecule( FILE_T*, MOLECULE_T*, INT_T, INT_T, INT_T);
INT_T		dumpresidue( FILE_T*, RESIDUE_T*, INT_T, INT_T );
INT_T		embed( BOUNDS_T*, REAL_T* );
MOLECULE_T     *fd_helix( STRING_T**, STRING_T**, STRING_T** );
INT_T		freemolecule( MOLECULE_T* );
INT_T		freeparm( MOLECULE_T* );
INT_T		freeresidue ( RESIDUE_T* );
#ifndef WIN32
STRING_T	*ftime( STRING_T* );
#endif
REAL_T      gauss( REAL_T*, REAL_T* );
INT_T		geodesics( BOUNDS_T* );
REAL_T		getchivol( MOLECULE_T**, STRING_T**, STRING_T**, STRING_T**, STRING_T** );
REAL_T		getchivolp( POINT_T, POINT_T, POINT_T, POINT_T );
MOLECULE_T	*getcif( STRING_T*, STRING_T* );
MOLECULE_T	*getcompcif( STRING_T*, STRING_T* );
STRING_T	*NAB_getline( FILE_T* );
REF_MATRIX_T	getmatrix( STRING_T* );
INT_T		getseq_from_pdb( STRING_T**, INT_T*, STRING_T**, STRING_T**, STRING_T** );
MOLECULE_T	*getpdb( STRING_T*, STRING_T* );
MOLECULE_T	*getpdb_prm( STRING_T**, STRING_T**, STRING_T**, INT_T* );
RESIDUE_T	*getresidue( STRING_T*, STRING_T* );
STRING_T	*getreslibkind( STRING_T* );
STRING_T	*getresname( RESIDUE_T* );
INT_T		getxv( STRING_T*, INT_T, REAL_T, REAL_T*, REAL_T* );
INT_T           getxyz( STRING_T**, INT_T*, REAL_T* );
INT_T		getxyz_from_pdb( STRING_T**, MOLECULE_T**, STRING_T**, INT_T* );
INT_T		helixanal( MOLECULE_T** );
INT_T		length( STRING_T*** );
MOLECULE_T 	*linkprot( STRING_T**, STRING_T**, STRING_T** );
MOLECULE_T 	*link_na( STRING_T**, STRING_T**, STRING_T**, STRING_T**, STRING_T** );
REAL_T		lmodC(INT_T*, INT_T*, INT_T*, INT_T*, INT_T*, REAL_T*, REAL_T*,
                 REAL_T*, INT_T*, REAL_T*, REAL_T*, REAL_T*, INT_T*, INT_T*,
                 INT_T*, INT_T*, INT_T*, REAL_T*, REAL_T*, REAL_T*, INT_T*,
                 INT_T*, INT_T*, INT_T*, INT_T*, INT_T*, REAL_T*, REAL_T*,
                 INT_T*, REAL_T*, REAL_T*, INT_T*, INT_T*, REAL_T*, REAL_T*,
                 INT_T*, INT_T* );
/*
REAL_T lmod(INT_T*, REAL_T*, REAL_T*, REAL_T*, REAL_T*, REAL_T*, 
            INT_T*, INT_T*, INT_T*,
            REAL_T*, REAL_T*, REAL_T*, REAL_T*, 
            struct xmin_opt*, struct lmod_opt* );
*/
REAL_T lmod();   /* don't test arguments for now....*/
INT_T  lmod_opt_init();   /* don't test arguments for now....*/
            
INT_T		md( INT_T, INT_T, REAL_T*, REAL_T*, REAL_T*,
		    REAL_T ( *func )( REAL_T*, REAL_T*, INT_T* ) );
INT_T		mdrat( INT_T, INT_T, REAL_T*, REAL_T*, REAL_T*, REAL_T*, REAL_T ( *mme )() );
INT_T		mergestr( MOLECULE_T*, STRING_T*, STRING_T*, MOLECULE_T*, STRING_T*, STRING_T* );
INT_T		metrize( BOUNDS_T*, INT_T );
INT_T		mm_options( STRING_T* );
void            mm_set_checkpoint( STRING_T** );
REAL_T		mme( REAL_T*, REAL_T*, INT_T* );
REAL_T		mme2( REAL_T*, REAL_T*, REAL_T*, REAL_T*, 
		      REAL_T*, INT_T*, INT_T*, INT_T*,
		      INT_T *, INT_T *, INT_T *, INT_T *,
		      INT_T*, INT_T* , char* );
INT_T           mme2_timer();
REAL_T		mme4( REAL_T*, REAL_T*, INT_T* );
REAL_T		mme_rattle( REAL_T*, REAL_T*, INT_T* );
INT_T		mme_init( MOLECULE_T*, STRING_T*, STRING_T*, REAL_T*, STRING_T* );
INT_T		mme_timer();
void            mme_rism_max_memory();
REAL_T		molsurf( MOLECULE_T**, STRING_T**, REAL_T* );
INT_T           mpierror( INT_T );
INT_T           mpifinalize( void );
INT_T           mpiinit( INT_T*, STRING_T**, INT_T*, INT_T* );
BOUNDS_T	*newbounds( MOLECULE_T*, STRING_T* );
MOLECULE_T	*newmolecule();
INT_T		newton( REAL_T*, INT_T*, REAL_T*,
			REAL_T ( *func1 )( REAL_T*, REAL_T*, INT_T* ),
			REAL_T ( *func2 )( REAL_T*, REAL_T*, REAL_T*, REAL_T*,
					   REAL_T*, INT_T*, INT_T*, INT_T*,
					   INT_T*, INT_T*, INT_T*, INT_T*,
					   INT_T*, INT_T*, char* ),
			REAL_T*, REAL_T*, INT_T* );
REF_MATRIX_T	newtransform( REAL_T, REAL_T, REAL_T, REAL_T, REAL_T, REAL_T );
INT_T		nmode( REAL_T*, INT_T,
			REAL_T ( *func)( REAL_T*, REAL_T*, REAL_T*, REAL_T*,
					 REAL_T*, INT_T*, INT_T*, INT_T*,
					 INT_T*, INT_T*, INT_T*, INT_T*,
					 INT_T*, INT_T*, char*  ),
		       INT_T, INT_T, REAL_T, REAL_T , INT_T);
INT_T           nm_timer();
INT_T		orbounds( BOUNDS_T*, MOLECULE_T*, STRING_T*, STRING_T*, REAL_T, REAL_T );
REAL_T		pair_ener( STRING_T*, INT_T );
INT_T		plane( MOLECULE_T**, STRING_T**, REAL_T*, REAL_T*, REAL_T* );
INT_T		putarc( STRING_T**, MOLECULE_T** );
INT_T		putbnd( STRING_T*, MOLECULE_T* );
INT_T		putcif( STRING_T*, STRING_T*, MOLECULE_T* );
INT_T		putdist( STRING_T*, MOLECULE_T*, STRING_T*, STRING_T* );
INT_T		putmatrix( STRING_T*, MATRIX_T );
INT_T		putpdb( STRING_T*, MOLECULE_T*, STRING_T* );
INT_T		putx( STRING_T**, MOLECULE_T** );
INT_T		putxv( STRING_T*, STRING_T*, INT_T, REAL_T, REAL_T*, REAL_T* );
INT_T           putxyz( STRING_T**, INT_T*, REAL_T* );
REAL_T		rand2( void );
INT_T		readbinposhdr( FILE* );
INT_T		readbinposfrm( INT_T, REAL_T*, FILE* );
INT_T		readparm( MOLECULE_T*, STRING_T* );
INT_T		rmsd( MOLECULE_T**, STRING_T**, MOLECULE_T**, STRING_T**, REAL_T*);
REF_MATRIX_T	rot4( MOLECULE_T*, STRING_T*, STRING_T*, REAL_T );
REF_MATRIX_T	rot4p( POINT_T, POINT_T, REAL_T );
INT_T       rseed( void );
FILE_T		*safe_fopen( STRING_T*, STRING_T* );
INT_T		sasad( REAL_T*, REAL_T*, REAL_T*, INT_T, REAL_T );
INT_T		setbounds( BOUNDS_T*, MOLECULE_T*, STRING_T*, STRING_T*, REAL_T, REAL_T );
INT_T		setboundsfromdb( BOUNDS_T**, MOLECULE_T**, STRING_T**, STRING_T**, 
		STRING_T**, REAL_T* );
INT_T		setchiplane( BOUNDS_T**, MOLECULE_T**, STRING_T** );
INT_T		setchivol( BOUNDS_T*, MOLECULE_T*, STRING_T*, STRING_T*, STRING_T*, STRING_T*, 
		REAL_T );
INT_T		setframe( INT_T, MOLECULE_T*, STRING_T*, STRING_T*, STRING_T*, STRING_T*, STRING_T* );
INT_T		setframep( INT_T, MOLECULE_T*, POINT_T, POINT_T, POINT_T, POINT_T, POINT_T );
INT_T		setmol_from_xyz( MOLECULE_T**, STRING_T**, REAL_T* );
INT_T		setmol_from_xyzw( MOLECULE_T**, STRING_T**, REAL_T* );
INT_T		setpoint( MOLECULE_T*, STRING_T*, POINT_T );
INT_T		setreskind( MOLECULE_T*, STRING_T*, STRING_T* );
INT_T		setreslibkind( STRING_T*, STRING_T* );
INT_T       setseed( INT_T* );
INT_T		setxyz_from_mol( MOLECULE_T**, STRING_T**, REAL_T* );
INT_T		setxyzw_from_mol( MOLECULE_T**, STRING_T**, REAL_T* );
INT_T		showbounds( BOUNDS_T*, MOLECULE_T*, STRING_T*, STRING_T* );
INT_T		split( STRING_T*, STRING_T**, STRING_T* );
REAL_T		step_ener( STRING_T*, INT_T );
STRING_T	*substr( STRING_T*, INT_T, INT_T );
INT_T		sugarpuckeranal( MOLECULE_T**, INT_T*, INT_T*, INT_T* );
REF_MATRIX_T		superimpose( MOLECULE_T*, STRING_T*, MOLECULE_T*, STRING_T* );
STRING_T	*timeofday();
REF_MATRIX_T	trans4( MOLECULE_T*, STRING_T*, STRING_T*, REAL_T );
REF_MATRIX_T	trans4p( POINT_T, POINT_T, REAL_T );
REAL_T		torsion( MOLECULE_T*, STRING_T*, STRING_T*, STRING_T*, STRING_T* );
REAL_T		torsionp( POINT_T, POINT_T, POINT_T, POINT_T );
INT_T		transformmol( MATRIX_T, MOLECULE_T*, STRING_T* );
INT_T		transformpts( MATRIX_T, POINT_T, INT_T );
RESIDUE_T	*transformres( MATRIX_T, RESIDUE_T*, STRING_T* );
INT_T		tsmooth( BOUNDS_T*, REAL_T );
INT_T		useboundsfrom( BOUNDS_T*, MOLECULE_T*, STRING_T*, MOLECULE_T*, STRING_T*, REAL_T );
INT_T		usemodeldist( BOUNDS_T*, MOLECULE_T*, STRING_T*, STRING_T* );
REF_MATRIX_T	updtransform( MATRIX_T, MATRIX_T );
MOLECULE_T	*wc_basepair( RESIDUE_T**, RESIDUE_T** );
STRING_T	*wc_complement( STRING_T**, STRING_T**, STRING_T** );
MOLECULE_T	*wc_helix( STRING_T**, STRING_T**, STRING_T**, STRING_T**, STRING_T**, STRING_T**,
		REAL_T*, REAL_T*, REAL_T*, REAL_T*, STRING_T** );
INT_T		writebinposhdr( FILE* );
INT_T		writebinposfrm( INT_T, REAL_T*, FILE* );
INT_T		writeparm( MOLECULE_T*, STRING_T* );
/*
REAL_T		xmin( REAL_T ( *func )( REAL_T*, REAL_T*, INT_T* ),
                      INT_T*, REAL_T*, REAL_T*, REAL_T*, REAL_T*, 
                      struct xmin_opt* );
*/
REAL_T      xmin();  /*  don't test arguments for now.... */
INT_T       xmin_opt_init();  /*  don't test arguments for now.... */
REAL_T      xminC(INT_T*, INT_T*, INT_T*, REAL_T*, INT_T*, INT_T*,
            INT_T*, REAL_T*, REAL_T*, REAL_T*, REAL_T*,
            INT_T*, REAL_T*, INT_T*, INT_T*, INT_T*, INT_T*,
            REAL_T*, REAL_T*, REAL_T*, REAL_T*,
            REAL_T*, REAL_T*, INT_T*, INT_T*);

/*  AmberNetcdf routines: (no argument checking) */
INT_T netcdfDebug();
INT_T netcdfLoad();
INT_T netcdfClose();
INT_T netcdfWriteRestart();
INT_T netcdfCreate();
INT_T netcdfGetVelocity();
INT_T netcdfGetFrame();
INT_T netcdfGetNextFrame();
INT_T netcdfWriteFrame();
INT_T netcdfWriteNextFrame();
INT_T netcdfInfo();

	/* output from molecular mechanics routines goes to nabout  */
extern FILE_T*  nabout;

	/* defines for hash table references:	*/

typedef	struct	{
	int	v_type;
	union	{
		int	v_ival;
		SIZE_T	v_szval;
		REAL_T	v_fval;
		char	*v_cval;
		REAL_T	v_ptval[ 3 ];
		char	*v_matval;
		FILE	*v_fpval;
		char	*v_atomval;
		char	*v_molval;
		char	*v_resval;
		char	*v_bval;
		char	*v_uval;
	} v_value;
} VALUE_T;

VALUE_T		*NAB_href();
STRING_T	*NAB_hin();
STRING_T	*NAB_hfirst();
STRING_T	*NAB_hnext();

#define	HRI(h,k,t)	NAB_href((h),(k),(t),0)->v_value.v_ival
#define	HRSZ(h,k,t)	NAB_href((h),(k),(t),0)->v_value.v_szval
#define	HRF(h,k,t)	NAB_href((h),(k),(t),0)->v_value.v_fval
#define	HRC(h,k,t)	NAB_href((h),(k),(t),0)->v_value.v_cval
#define	HRPT(h,k,t)	NAB_href((h),(k),(t),0)->v_value.v_ptval
#define	HRMAT(h,k,t)	NAB_href((h),(k),(t),0)->v_value.v_matval
#define	HRFP(h,k,t)	NAB_href((h),(k),(t),0)->v_value.v_fpval
#define	HRATOM(h,k,t)	NAB_href((h),(k),(t),0)->v_value.v_atomval
#define	HRRES(h,k,t)	NAB_href((h),(k),(t),0)->v_value.v_resval
#define	HRMOL(h,k,t)	NAB_href((h),(k),(t),0)->v_value.v_molval
#define	HRB(h,k,t)	NAB_href((h),(k),(t),0)->v_value.v_bval
#define	HRU(h,k,t,s,c)	(*((c)NAB_href((h),(k),(t),(s))->v_value.v_uval))

	/* defines & declares for points	*/

#define	RECIP(f)	( 1./(f) )

#define	PTX(p)	((p)[0])
#define	PTY(p)	((p)[1])
#define	PTZ(p)	((p)[2])

#define	PTEQ(p,q)	\
	((p)[0]==(q)[0]&&(p)[1]==(q)[1]&&(p)[2]==(q)[2])
#define	PTNE(p,q)	\
	((p)[0]!=(q)[0]||(p)[1]!=(q)[1]||(p)[2]!=(q)[2])
#define	PT_ISTRUE(p)	((p)[0]!=0.0||(p)[1]!=0.0||(p)[2]!=0.0)

POINT_T	*NAB_ptcpy();
POINT_T	*NAB_ptadd();
POINT_T	*NAB_ptsub();
POINT_T	*NAB_ptscl();
POINT_T	*NAB_ptcrs();
REAL_T	NAB_ptdot();

	/* functions for accessing parts of atoms & residues */

INT_T		*NAB_ari();
REAL_T		*NAB_arf();
STRING_T	**NAB_arc();
POINT_T		*NAB_arp( ATOM_T *ap, char key[] );
INT_T		*NAB_rri();
STRING_T	**NAB_rrc();
INT_T		*NAB_mri();

	/* functions for for( a in m ) etc	*/

ATOM_T		*NAB_mnext (MOLECULE_T *mol, ATOM_T *cap);
ATOM_T		*NAB_anext();
RESIDUE_T	*NAB_rnext();

	/* defines for string compares:	*/

int	NAB_strcmp( char *, char * );
#define	LT(a,b)	(NAB_strcmp((a),(b))<0)
#define	LE(a,b)	(NAB_strcmp((a),(b))<=0)
#define	EQ(a,b)	(NAB_strcmp((a),(b))==0)
#define	NE(a,b)	(NAB_strcmp((a),(b))!=0)
#define	GE(a,b)	(NAB_strcmp((a),(b))>=0)
#define	GT(a,b)	(NAB_strcmp((a),(b))>0)

	/* String stuff	*/

#define	length(s)	strlen(s)

char	*NAB_readstring();
char	*NAB_strcpy();
char	*NAB_strcat();
int	NAB_index( char [], char [] );
int	NAB_rematch( char [], char [] );
int NAB_gsub(int, char **, char **, char **);

	/*  Other NAB declarations:  */

int	NAB_matcpy();
int	unlink();
int	NAB_aematch( ATOM_T *ap, char aex[] );

	/* defines for assigning then using temp vars in exprs	*/

#define	ITEMP(t,e)	((t)=(e),&(t))
#define	SZTEMP(t,e)	((t)=(e),&(t))
#define	FTEMP(t,e)	((t)=(e),&(t))
#define	STEMP(t,e)	(NAB_strcpy(&(t),(e)),&(t))
#define	FPTEMP(t,e)	((t)=(e),&(t))
#define	SPRINTF(s,t)	((s),(t))

	/* trig functions in degrees:	*/

#define	R2D	57.29577951308232090712
#define	D2R	 0.01745329251994329576
#define	ACOS(c)	(R2D*acos(c))
#define	ASIN(s)	(R2D*asin(s))
#define	ATAN(t)	(R2D*atan(t))
#define	ATAN2(y,x)	(R2D*atan2((y),(x)))
#define	COS(a)	cos(D2R*(a))
#define	SIN(a)	sin(D2R*(a))
#define	TAN(a)	tan(D2R*(a))

	/* dynamic array allocation macro:	*/

#define DA_ALLOC(a,fn,an)	\
	if(!(a)){fprintf( stderr,	\
	"%s: can't allocate space for dynamic array \"%s\"\n",	\
		(fn),(an));exit(1);}

	/* local string array zero macro:	*/

#define	NAB_AINIT(a)	memset((a),0,sizeof(a))
#define	NAB_SINIT(s)	memset(&(s),0,sizeof(s))

#endif
