#ifndef	NABTYPES_H
#define	NABTYPES_H

#include <stdio.h>
#include "defreal.h"

	/* Fundamental nab types:	*/

typedef	int	INT_T;
typedef	size_t	SIZE_T;

	/* geometric types: Frames are not directly accessible 	*/
	/* at the nab level.					*/

typedef	REAL_T	POINT_T[ 3 ];		/* 0 == x, 1 == y, 2 == z */
typedef	REAL_T	FRAME_T[ 4 ][ 3 ];	/* org, x, y, z axes	*/
typedef	REAL_T	MATRIX_T[ 4 ][ 4 ];	
/* If MATRIX_T m; then REF_MATRIX_T is the type of the array name m. */
/* Functions returning, in principle, MATRIX_T should have a return value */
/* type of REF_MATRIX_T.  Note that MATRIX_T* is not REF_MATRIX_T. */
typedef	REAL_T	( *REF_MATRIX_T )[ 4 ];

	/* other nab types:	*/
typedef	char	STRING_T;
typedef	FILE	FILE_T;

	/* atom, residue & molecule types, includes types	*/
	/* INTBOND_T, EXTBOND_T, CHIRAL_T & STRAND_T that have	*/
	/* no nab equivalent					*/

#define	A_CONNECT_SIZE	8
typedef	struct	atom_t	{
	STRING_T *a_atomname;
	STRING_T *a_atomtype;
	INT_T	a_attr;
	INT_T	a_nconnect;
	INT_T	a_connect[ A_CONNECT_SIZE ];
	struct	residue_t	*a_residue;
	REAL_T	a_charge;
	REAL_T	a_radius;
	REAL_T	a_bfact;
	REAL_T	a_occ;
	STRING_T *a_element;	/* from SD (mol) files		*/
	INT_T	a_int1;		/* user property		*/
	REAL_T	a_float1;	/* user property		*/
	REAL_T	a_float2;	/* user property		*/
	INT_T	a_tatomnum;
	INT_T	a_atomnum;
	STRING_T *a_fullname;
	POINT_T	a_pos;
	REAL_T	a_w;		/* 4th dimension		*/
} ATOM_T;

typedef	struct	extbond_t	{
	struct	extbond_t	*eb_next;
	INT_T	eb_anum;	/* atom in current residue	*/
	INT_T	eb_rnum;	/* other residue number		*/
	INT_T	eb_ranum;	/* atom in other residue	*/
} EXTBOND_T;

typedef	INT_T	INTBOND_T[ 2 ];

	/* chiral centers:	*/

typedef	struct	chiral_t	{
	INT_T	c_anum[ 4 ];
	REAL_T	c_dist;
} CHIRAL_T;

typedef	struct	residue_t	{
	struct	residue_t	*r_next;
	INT_T	r_num;
	INT_T	r_tresnum;	/* set by NAB_rri( a, "tresnum" )  */
	INT_T	r_resnum;	/* set by NAB_rri( a, "resnum" )  */
	STRING_T *r_resname;	/* set by NAB_rrc( a, "resname" ) */
	STRING_T *r_resid;	/* set by NAB_rrc( a, "resid" ) */
	INT_T	r_attr;
	INT_T	r_kind;
	INT_T	r_atomkind;
	struct	strand_t	*r_strand;
	EXTBOND_T	*r_extbonds;	
	INT_T	r_nintbonds;	/* INTERNAL bonds		*/
	INTBOND_T	*r_intbonds;
	INT_T		r_nchiral;
	CHIRAL_T	*r_chiral;
	INT_T	r_natoms;
	INT_T	*r_aindex;
	ATOM_T	*r_atoms;
} RESIDUE_T;

typedef	struct	strand_t	{	/* not visible in nab 	*/
	STRING_T *s_strandname;
	INT_T	s_strandnum;
	INT_T	s_attr;
	struct	molecule_t	*s_molecule;
	struct	strand_t	*s_next;
	INT_T	s_nresidues;
	INT_T	s_res_size;
	RESIDUE_T	**s_residues;
} STRAND_T;

typedef struct parm {
	char	ititl[81];
	INT_T 	IfBox, Nmxrs, IfCap,
		 Natom,  Ntypes,  Nbonh,  Mbona,  Ntheth,  Mtheta, 
		 Nphih,  Mphia,  Nhparm, Nparm, Nnb, Nres,
		 Nbona,  Ntheta,  Nphia,  Numbnd,  Numang,  Nptra,
		 Natyp,  Nphb, Nat3, Ntype2d, Nttyp, Nspm, Iptres, Nspsol,
		 Ipatm, Natcap, Numextra;
	STRING_T *AtomNames, *ResNames, *AtomSym, *AtomTree;
	REAL_T	*Charges, *Masses, *Rk, *Req, *Tk, *Teq, *Pk, *Pn, *Phase,
		 *Solty, *Cn1, *Cn2, *HB12, *HB10, *Rborn, *Fs;
	REAL_T	Box[4], Cutcap, Xcap, Ycap, Zcap, Fsmax;
	INT_T 	*Iac, *Iblo, *Cno, *Ipres, *ExclAt, *TreeJoin, 
		*AtomRes, *BondHAt1, *BondHAt2, *BondHNum, *BondAt1, *BondAt2, 
		*BondNum, *AngleHAt1, *AngleHAt2, *AngleHAt3, *AngleHNum, 
		*AngleAt1, *AngleAt2, *AngleAt3, *AngleNum, *DihHAt1, 
		*DihHAt2, *DihHAt3, *DihHAt4, *DihHNum, *DihAt1, *DihAt2, 
		*DihAt3, *DihAt4, *DihNum, *Boundary;
	INT_T	*N14pairs, *N14pairlist;
	REAL_T *Gvdw;
} PARMSTRUCT_T;


typedef	struct	molecule_t	{
	FRAME_T	m_frame;
	INT_T	m_nstrands;
	STRAND_T	*m_strands;
	INT_T	m_nresidues;
	INT_T	m_natoms;
	INT_T	m_nvalid;	/* "numbers" valid	*/
	PARMSTRUCT_T	*m_prm;	/* points forcefield stuff */
} MOLECULE_T;

	/* Distance geometry "bounds" matrix	*/

typedef	struct	bounds_t	{
	MOLECULE_T	*b_mol;
	INT_T		b_nres;
	INT_T		b_natoms;
	INT_T		*b_aoff;
	ATOM_T		**b_aindex;
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

#endif
