#ifndef INC_XMIN_OPT_H
#define INC_XMIN_OPT_H

//  options structure for xmin()

struct xmin_opt {
	int mol_struct_opt;
	int maxiter;
	float grms_tol;
	int method;
	int numdiff;
	int m_lbfgs;
	int iter;        // output
	float xmin_time; // output
        int ls_method;
        int ls_maxiter;
        float ls_maxatmov;
        float beta_armijo;
        float c_armijo;
        float mu_armijo;
        float ftol_wolfe;
        float gtol_wolfe;
        int ls_iter;     // output
	int print_level;
	int error_flag;  // output
};

#endif
