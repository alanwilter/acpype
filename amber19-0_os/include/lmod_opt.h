//  options structure for lmod()

struct lmod_opt {
	int niter; 
	int nmod; 
	int kmod;
	int nrotran_dof; 
	int nconf; 
	float minim_grms;
	float energy_window;
	float conf_separation_rms;
	int eig_recalc;
	int ndim_arnoldi; 
	int lmod_restart; 
	int n_best_struct;
	int mc_option; 
	float rtemp; 
	float lmod_step_size_min;
	float lmod_step_size_max; 
	int nof_lmod_steps;
	float lmod_relax_grms; 
	int nlig; 
	int apply_rigdock;
	int nof_poses_to_try;
	int random_seed;
	int print_level; 
	float lmod_time; 
	float aux_time;
	int error_flag;
};
