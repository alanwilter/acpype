/** This header file defines the C/C++ interface into the SANDER API that can be
 * used to compute energies and forces for a given topology file and arbitrary
 * coordinates.  It supports a standard MM potential as well as a hybrid QM-MM
 * potential as well.
 *
 * Input options are handled through structs, and this module provides an
 * interface to interact with the Fortran subroutines defined in libsander.so.
 * Make sure that the structs here match up correctly with the compound types
 * defined in the sander interface (and they have to be in the same order!!!)
 *
 * To get the LES API, compile with -DLES and link libsanderles.so instead of
 * libsander.so
 */
#ifndef SANDER_H
#define SANDER_H
// Symbol translation table (calling Fortran in C)

// These are the raw functions that we expose unaltered from the F90 API
#define pme_sander_input ext_pme_sander_input_
#define qm_sander_input ext_qm_sander_input_
#define energy_forces ext_energy_forces_
#define set_positions ext_set_positions_
#define sander_cleanup ext_sander_cleanup_
#define get_positions ext_get_positions_
#define get_box ext_get_box_

// These are the functions whose C-interface we want to change slightly from the
// F90 interface, so we name-mangle to keep them internal. These are _NOT_ part
// of the public C/C++ API
#define __internal_gas_sander_input ext_gas_sander_input_
#define __internal_set_box ext_set_box_
#define __internal_sander_setup ext_sander_setup_
#define __internal_sander_setup2 ext_sander_setup2_
#define __internal_sander_natom ext_sander_natom_
#define __internal_read_inpcrd_file ext_read_inpcrd_file_
#define __internal_get_inpcrd_natom ext_get_inpcrd_natom_
#define __internal_is_setup ext_is_setup_

#ifdef __cplusplus
#   include <cstdlib>
#   include <cstring>
#   include <string>
#else
#   include <stdlib.h>
#   include <string.h>
#endif

// Data types
typedef struct {
    // Floating point input parameters
    double extdiel;
    double intdiel;
    double rgbmax;
    double saltcon;
    double cut;
    double dielc;
    double rdt;
    double fswitch;
    double restraint_wt;

    // Integer choice flags
    int igb;
    int alpb;
    int gbsa;
    int lj1264;
    int ipb;
    int inp;
    int vdwmeth;
    int ew_type;
    int ntb;
    int ifqnt;
    int jfastw;
    int ntf;
    int ntc;
    int ntr;
    int ibelly;
    int mask_from_ref;

    // Strings
    char restraintmask[256];
    char bellymask[256];
    char refc[256];
} sander_input;

 /* Make sure this is the same as the number defined in lib/constants.F90!!! */
#define MAX_QUANTUM_ATOMS 10000
typedef struct {
    // Floats (input parameters)
    double qmcut;
    double lnk_dis;
    double scfconv;
    double errconv;
    double dftb_telec;
    double dftb_telec_step;
    double fockp_d1;
    double fockp_d2;
    double fockp_d3;
    double fockp_d4;
    double damp;
    double vshift;
    double kappa;
    double pseudo_diag_criteria;
    double min_heavy_mass;
    double r_switch_hi;
    double r_switch_lo;
    // Integers
    int iqmatoms[MAX_QUANTUM_ATOMS];
    int qmgb;
    int lnk_atomic_no;
    int ndiis_matrices;
    int ndiis_attempts;
    int lnk_method;
    int qmcharge;
    int corecharge;
    int buffercharge;
    int spin;
    int qmqmdx;
    int verbosity;
    int printcharges;
    int printdipole;
    int print_eigenvalues;
    int peptide_corr;
    int itrmax;
    int printbondorders;
    int qmshake;
    int qmmmrij_incore;
    int qmqm_erep_incore;
    int pseudo_diag;
    int qm_ewald;
    int qm_pme;
    int kmaxqx;
    int kmaxqy;
    int kmaxqz;
    int ksqmaxq;
    int qmmm_int;
    int adjust_q;
    int tight_p_conv;
    int diag_routine;
    int density_predict;
    int fock_predict;
    int vsolv;
    int dftb_maxiter;
    int dftb_disper;
    int dftb_chg;
    int abfqmmm;
    int hot_spot;
    int qmmm_switch;
    int core_iqmatoms[MAX_QUANTUM_ATOMS];
    int buffer_iqmatoms[MAX_QUANTUM_ATOMS];
    // Strings
    char qmmask[8192];
    char coremask[8192];
    char buffermask[8192];
    char centermask[8192];
    char dftb_3rd_order[256];
    char dftb_slko_path[256];
    char qm_theory[12];
} qmmm_input_options;

typedef struct {
    double tot;
    double vdw;
    double elec;
    double gb;
    double bond;
    double angle;
    double dihedral;
    double vdw_14;
    double elec_14;
    double constraint;
    double polar;
    double hbond;
    double surf;
    double scf;
    double disp;
    double dvdl;
    double angle_ub;
    double imp;
    double cmap;
    double emap;
    double les;
    double noe;
    double pb;
    double rism;
    double ct;
    double amd_boost;
} pot_ene;

typedef struct {

    // Integer pointers

    int natom, ntypes, nbonh,  mbona, ntheth, mtheta, nphih,
        mphia, nhparm, nparm,  nnb  , nres  , nbona , ntheta,
        nphia, numbnd, numang, nptra, natyp , nphb  , ifpert,
        nbper, ngper , ndper , mbper, mgper , mdper , ifbox,
        nmxrs, ifcap,
        numextra, ncopy;

    int iptres, nspm, nspsol; // solvent pointers

    int natcap, nlesty;

    // Bool flags

    int is_chamber;
    int ipol;

    // Other single-number parameters

    double cutcap, xcap, ycap, zcap;

    // CHARMM-specific pointers

    int cmap_term_count, cmap_type_count, charmm_nub, charmm_nubtypes,
        charmm_nimphi, charmm_nimprtyp;

    // Data arrays

    char title[80];
    char *atom_name;
    char *residue_label;
    char *amber_atom_type;
    char *tree_chain_classification;

    int *atomic_number;
    int *atom_type_index;
    int *number_excluded_atoms;
    int *nonbonded_parm_index;
    int *residue_pointer;
    int *bonds_inc_hydrogen;
    int *bonds_without_hydrogen;
    int *angles_inc_hydrogen;
    int *angles_without_hydrogen;
    int *dihedrals_inc_hydrogen;
    int *dihedrals_without_hydrogen;
    int *excluded_atoms_list;
    int *join_array;
    int *irotat;
    int *atoms_per_molecule;
    int *les_type;
    int *les_fac;
    int *les_cnum;
    int *les_id;
    int *charmm_cmap_resolution;
    int *charmm_cmap_index;
    int *charmm_impropers;
    int *charmm_urey_bradley;

    double *charge;
    double *mass;
    double *bond_force_constant;
    double *bond_equil_value;
    double *angle_force_constant;
    double *angle_equil_value;
    double *dihedral_force_constant;
    double *dihedral_periodicity;
    double *dihedral_phase;
    double *scee_scale_factor;
    double *scnb_scale_factor;
    double *solty;
    double *lennard_jones_acoef;
    double *lennard_jones_bcoef;
    double *lennard_jones_ccoef;
    double *expvdwmodel_beta;
    double *expvdwmodel_a;
    double *expvdwmodel_b;
    double *hbond_acoef;
    double *hbond_bcoef;
    double *hbcut;
    double *box_dimensions;
    double *radii;
    double *screen;
    double *polarizability;
    double *dipole_damp_factor;
    double *ti_mass;
    double *charmm_cmap_parameter;
    double *charmm_improper_force_constant;
    double *charmm_improper_phase;
    double *lennard_jones_14_acoef;
    double *lennard_jones_14_bcoef;
    double *charmm_urey_bradley_equil_value;
    double *charmm_urey_bradley_force_constant;

} prmtop_struct;

#ifdef __cplusplus
extern "C" {
#endif
/// Prepare a sander input struct for PME electrostatics
void pme_sander_input(sander_input*);

/// Prepare a QM input struct with default values
void qm_sander_input(qmmm_input_options*);

/* I've found that you really need to fix strings to the same number of
 * characters when you want to pass them from C to Fortran or vice-versa. As a
 * result, all file names passed into here will be copied into a container of
 * the same length as the file name holders in the Fortran routines. As a
 * result, you MUST (!!) make sure that __MAX_FN_LEN here stays consistent with
 * MAX_FN_LEN as defined in file_io_dat.F90
 */
#define __MAX_FN_LEN 256

/* Define the function prototype for internal functions, these should not be
 * called by programs using the API
 */
void __internal_sander_setup(const char[__MAX_FN_LEN], double*, double*,
                             sander_input*, qmmm_input_options*, int*);
void __internal_sander_natom(int *);
void __internal_read_inpcrd_file(const char[__MAX_FN_LEN], double*, double*, int*);
void __internal_get_inpcrd_natom(const char[__MAX_FN_LEN], int*);
void __internal_set_box(double *a, double *b, double *c,
                        double *alpha, double *beta, double *gamma);
void __internal_sander_setup2(prmtop_struct*, double*, double*, sander_input*,
                              qmmm_input_options*, int*);
void __internal_is_setup(int*);
void __internal_gas_sander_input(sander_input*, int*);

static inline void gas_sander_input(sander_input *inp, const int gb) {
    int dum = gb;
    __internal_gas_sander_input(inp, &dum);
}

/** Set up an energy calculation
 * \param prmname Name of the prmtop file defining the force field
 * \param crdname Name of the input coordinate file defining the initial
 *                positions (and box dimensions, if applicable)
 * \param input_options struct of input options for MM terms
 * \param qmmm_options struct of input options for QM part
 * \returns 0 for success, 1 for failure
 */
static inline int sander_setup(const char *prmname, double *coords, double *box,
                        sander_input *input_options, qmmm_input_options *qmmm_options) {
    int ierr;
    char *prmtop;
    prmtop = (char*)malloc(__MAX_FN_LEN*sizeof(char));
    strncpy(prmtop, prmname, __MAX_FN_LEN);
    __internal_sander_setup(prmtop, coords, box, input_options, qmmm_options, &ierr);
    free(prmtop);
    return ierr;
}
static inline int sander_setup_mm(const char *prmname, double *coords, double *box,
                           sander_input *input_options) {
    int ierr;
    char *prmtop;
    qmmm_input_options dummy;
    prmtop = (char*)malloc(__MAX_FN_LEN*sizeof(char));
    strncpy(prmtop, prmname, __MAX_FN_LEN);
    __internal_sander_setup(prmtop, coords, box, input_options, &dummy, &ierr);
    free(prmtop);
    return ierr;
}

static inline int sander_setup2(prmtop_struct *parm, double *coords, double *box,
                         sander_input *inp, qmmm_input_options *qm_inp) {
   int ierr;
   __internal_sander_setup2(parm, coords, box, inp, qm_inp, &ierr);
   return ierr;
}

static inline int sander_setup2_mm(prmtop_struct *parm, double *coords, double *box,
                            sander_input *input_options) {
    qmmm_input_options dummy;
    return sander_setup2(parm, coords, box, input_options, &dummy);
}

/** Sets the particle positions
 * \param positions Vector of doubles with 3N cartesian coordinates
                    (X1, Y1, Z1, * X2, Y2, Z2, ...)
 */
void set_positions(double *positions);

/** Sets the periodic box vectors
 * \param a Length of the first unit cell vector
 * \param b Length of the second unit cell vector
 * \param c Length of the third unit cell vector
 * \param alpha Angle between the second and third vectors
 * \param beta Angle between the first and third vectors
 * \param gamma Angle between the first and second vectors
 */
static inline void
set_box(double a, double b, double c, double alpha, double beta, double gamma) {
   __internal_set_box(&a, &b, &c, &alpha, &beta, &gamma);
}

/** Gets the periodic box vectors
 * \param a Length of the first unit cell vector
 * \param b Length of the second unit cell vector
 * \param c Length of the third unit cell vector
 * \param alpha Angle between the second and third vectors
 * \param beta Angle between the first and third vectors
 * \param gamma Angle between the first and second vectors
 */
void
get_box(double *a, double *b, double *c,
        double *alpha, double *beta, double *gamma);

/** Returns 1 if sander has been set up and 0 otherwise
 */
static inline int is_setup() {
   int i;
   __internal_is_setup(&i);
   return i;
}

/// Gets the current positions of the currently set up system
void get_positions(double *coordinates);

/** Computes energies and forces
 * \param energy pot_ene struct where the energy components will be stored
 * \param forces Double array of length 3*natom where forces will be stored
 */
void energy_forces(pot_ene *energy, double *forces);

/// Cleans up and frees allocated data.
void sander_cleanup(void);

/// Returns the number of atoms defined in the system
static inline int sander_natom(void) {
    int natom;
    __internal_sander_natom(&natom);
    return natom;
}

/// Fills the coordinates and box. Returns 0 for success, 1 for error
static inline int
read_inpcrd_file(const char *filename, double *coordinates, double *box) {
    int inerr;
    char *fname;

    fname = (char *)malloc(__MAX_FN_LEN*sizeof(char));
    strncpy(fname, filename, __MAX_FN_LEN);
    __internal_read_inpcrd_file(fname, coordinates, box, &inerr);
    free(fname);
    return inerr;
}

static inline int get_inpcrd_natom(const char *filename) {
    int natom;
    char *fname;

    fname = (char *)malloc(__MAX_FN_LEN*sizeof(char));
    strncpy(fname, filename, __MAX_FN_LEN);
    __internal_get_inpcrd_natom(fname, &natom);
    free(fname);
    return natom;
}

#define SAFE_FREE(var) if (parm->var) free(parm->var)
static inline void destroy_prmtop_struct(prmtop_struct *parm) {
    SAFE_FREE(atom_name);
    SAFE_FREE(residue_label);
    SAFE_FREE(amber_atom_type);
    SAFE_FREE(tree_chain_classification);

    SAFE_FREE(atomic_number);
    SAFE_FREE(atom_type_index);
    SAFE_FREE(number_excluded_atoms);
    SAFE_FREE(nonbonded_parm_index);
    SAFE_FREE(residue_pointer);
    SAFE_FREE(bonds_inc_hydrogen);
    SAFE_FREE(bonds_without_hydrogen);
    SAFE_FREE(angles_inc_hydrogen);
    SAFE_FREE(angles_without_hydrogen);
    SAFE_FREE(dihedrals_inc_hydrogen);
    SAFE_FREE(dihedrals_without_hydrogen);
    SAFE_FREE(excluded_atoms_list);
    SAFE_FREE(join_array);
    SAFE_FREE(irotat);
    SAFE_FREE(atoms_per_molecule);
    SAFE_FREE(les_type);
    SAFE_FREE(les_fac);
    SAFE_FREE(les_cnum);
    SAFE_FREE(les_id);
    SAFE_FREE(charmm_cmap_resolution);
    SAFE_FREE(charmm_cmap_index);
    SAFE_FREE(charmm_impropers);
    SAFE_FREE(charmm_urey_bradley);

    SAFE_FREE(charge);
    SAFE_FREE(mass);
    SAFE_FREE(bond_force_constant);
    SAFE_FREE(bond_equil_value);
    SAFE_FREE(angle_force_constant);
    SAFE_FREE(angle_equil_value);
    SAFE_FREE(dihedral_force_constant);
    SAFE_FREE(dihedral_periodicity);
    SAFE_FREE(dihedral_phase);
    SAFE_FREE(scee_scale_factor);
    SAFE_FREE(scnb_scale_factor);
    SAFE_FREE(solty);
    SAFE_FREE(lennard_jones_acoef);
    SAFE_FREE(lennard_jones_bcoef);
    SAFE_FREE(lennard_jones_ccoef);
    SAFE_FREE(expvdwmodel_beta);
    SAFE_FREE(expvdwmodel_a);
    SAFE_FREE(expvdwmodel_b);
    SAFE_FREE(hbond_acoef);
    SAFE_FREE(hbond_bcoef);
    SAFE_FREE(hbcut);
    SAFE_FREE(box_dimensions);
    SAFE_FREE(radii);
    SAFE_FREE(screen);
    SAFE_FREE(polarizability);
    SAFE_FREE(dipole_damp_factor);
    SAFE_FREE(ti_mass);
    SAFE_FREE(charmm_cmap_parameter);
    SAFE_FREE(charmm_improper_force_constant);
    SAFE_FREE(charmm_improper_phase);
    SAFE_FREE(lennard_jones_14_acoef);
    SAFE_FREE(lennard_jones_14_bcoef);
    SAFE_FREE(charmm_urey_bradley_equil_value);
    SAFE_FREE(charmm_urey_bradley_force_constant);
}
#undef SAFE_FREE
#ifdef __cplusplus
} /* extern "C" */

// Some C++-overloaded functions to accept strings in addition to char*
static inline int get_incprd_natom(const std::string& filename) {
   return get_inpcrd_natom(filename.c_str());
}

static inline int
read_inpcrd_file(const std::string& filename, double *coordinates, double *box)
{
   return read_inpcrd_file(filename.c_str(), coordinates, box);
}
#endif

#endif /* SANDER_H */
