#ifndef __NAB__
#define	__NAB__

#ifndef	FALSE
#define	FALSE	0
#endif

#ifndef	TRUE
#define	TRUE	1
#endif

#define	UNDEF	(-1)

	/* Run time stuff:	*/

#include "defreal.h"
#include "nabtypes.h"

	/* value type for expressions	*/
	/* also used by compiler	*/

typedef	struct	value_t	{
	int	v_type;
	union	{
		char	*v_aval;
		int	v_ival;
		size_t	v_szval;
		REAL_T	v_fval;
		char	*v_cval;
		POINT_T	v_ptval;
		MATRIX_T *v_matval;
		FILE	*v_fpval;
		ATOM_T	*v_atomval;
		RESIDUE_T *v_resval;
		MOLECULE_T *v_molval;
		BOUNDS_T *v_bval;
		char	*v_uval;
	} v_value;
} VALUE_T;

	/* "hashed" or associative arrays:	*/

typedef	struct	h_entry_t	{
	struct	h_entry_t	*he_next;
	char	*he_key;
	VALUE_T	he_val;
} H_ENTRY_T;

typedef	struct	hash_t	{
	int	h_type;
	int	h_llen;		/* length of the longest list	*/
	int	h_n_entries;
	int	h_he_size;
	H_ENTRY_T	**h_entries;
} HASH_T;

typedef	struct	curhash_t	{
	int		c_index;
	H_ENTRY_T	*c_entry;
} CURHASH_T;

	/* attributes of atoms, residues and molecules.	*/
	/* AT_SELECT is set by select_atoms() when it 	*/
	/* evaluates a regular expression.  AT_SELECTED	*/
	/* is used to store the previous select value	*/
	/* and is used in the construction of bounds	*/
	/* and dist matrices				*/
	
#define	AT_SELECT	0001
#define	AT_SELECTED	0002
#define	AT_WORK		0200

	/* Residue & Atom "kinds"	*/

#define	RT_UNDEF	0
#define	RT_DNA		1
#define	RT_RNA		2
#define	RT_AA		3

#define	RAT_UNDEF	0
#define	RAT_UNITED	1
#define	RAT_ALLATOM	2

	/* Compile time stuff:			*/
	/* nab types (used by both r/t & c/t:	*/
	/* user types go from T_UNDEF = 0, to T_ERROR = ?	*/

#define	T_UNDEF		0
#define	T_INT		1
#define	T_SIZE_T	2
#define	T_FLOAT		3
#define	T_STRING	4
#define	T_POINT		5
#define	T_MATRIX	6
#define	T_FILE		7
#define	T_ATOM		8
#define	T_RESIDUE	9
#define	T_MOLECULE	10
#define	T_BOUNDS	11
#define	T_NULL		12
#define	T_USER		13
#define	T_ERROR		14	/* stops printing of spur. e-msgs	*/
#define	T_HASH		15	/* never input */
#define	T_CURHASH	16	/* never input */
#define	T_LDABOUND	17	/* local dyn array bnd. never input	*/
#define	T_GDABOUND	18	/* global dyn array bnd. never input	*/
#define	N_TYPES		19
#define	T_LASTUSER	T_ERROR	/* last type that appears in nab input	*/

	/* From here to end used only at c/t	*/
	/* nab classes:	*/

#define	C_UNDEF		0
#define	C_LIT		1
#define	C_STRTAG	2	/* struct tag	*/
#define	C_VAR		3
#define	C_FUNC		4
#define	C_DEFINE	5
#define	C_EXPR		6
#define	C_NULL		7
#define	C_ERROR		8

	/* expression element kinds:	*/

#define	K_UNDEF		0
#define	K_SCALAR	1
#define	K_HASHEL	2
#define	K_DARRAYEL	3
#define	K_ARRAY		4
#define	K_DARRAY	5
#define	K_HASHED	6
#define	K_ERROR		7

	/* attribute access:	*/

#define	A_UNDEF		0
#define	A_DIRECT	1
#define	A_FUNCTION	2
#define	A_STRUCT	3

	/* parse tree:	*/

typedef	struct	node_t	{
	int	n_lineno;
	int	n_sym;
	int	n_type;
	int	n_class;
	int	n_kind;
	VALUE_T	n_val;
	struct	node_t	*n_left;
	struct	node_t	*n_right;
} NODE_T;

	/* nab symbol scopes:	*/

#define	S_UNDEF		0
#define	S_SYSTEM	1
#define	S_GLOBAL	2
#define	S_LOCAL		3
#define	S_USER		4	/* decl's in struct's have user scope	*/

	/* nab calling conventions */

#define	CC_UNDEF	0
#define	CC_NAB		1
#define	CC_CC		2
#define	CC_FORTRAN	3	
#define	CC_IO		4

	/* Symbols:	*/

typedef	struct	symrec_t	{
	char	*s_name;
	int	s_type;
	int	s_class;
	int	s_kind;
	int	s_scope;
	int	s_cconv;	/* calling convention			*/
	int	s_isdyn;	/* symbol is a dynamic array 	*/
	int	s_isparm;	/* symbol is a formal			*/
	int	s_init;		/* struct w/strings, needs init	*/
	int	s_pcount;	/* components, func parms, array dims	*/
	NODE_T	**s_parts;
	struct	symrec_t	*s_left;
	struct	symrec_t	*s_right;
	char	*s_uname;	/* the symbol tag	*/
	struct	symrec_t	*s_user;	/* hold decl's for structs	*/
	int	s_nusyms;
} SYMREC_T;

	/* code generator debug stuff:	*/

#define	CGD_NONE	000
#define	CGD_VARDECL	001
#define	CGD_CHKEXPR	002
#define	CGD_FIXEXPR	004
#define	CGD_EXPRCODE	010
#define	CGD_PARSER	020
#define	CGD_ALL		( CGD_VARDECL |\
			  CGD_CHKEXPR | CGD_FIXEXPR | CGD_EXPRCODE |\
			  CGD_PARSER )

	/* AVS stuff:	*/

#define	AI_UNDEF	0
#define	AI_PARM		1
#define	AI_FUNC		2
#define	AI_PORT		3
#define	AI_FREE		4
#define	AI_SEND		5
typedef	struct	avsinfo_t	{
	struct	avsinfo_t	*a_next;
	int	a_type;
	char	*a_name;
	char	*a_data;
} AVSINFO_T;

NODE_T	*node( int, VALUE_T *, NODE_T *, NODE_T * );
NODE_T	*copynode( NODE_T * );

	/*  Commonly used functions:    */

int	atom_in_aexpr( ATOM_T *, char [] );
char *  compile( char *, char *, char *, int );
int step( char *,  char * );

	/* Size of the buffer used by snprintf()	*/
	/* the buffer itself is in stringutil.c		*/

#define	NAB_RSBUF_SIZE	10000

	/* the output for all non error emissions */

extern FILE	*nabout;

#endif
