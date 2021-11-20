#ifndef PDB_HPP
#define PDB_HPP

#include <vector>
#include <string>
#include <map>

struct Point3d{	
	int res_, num_;
	double x_, y_, z_;
	Point3d(int r, int n, double x, double y, double z) : res_(r), num_(n), x_(x), y_(y), z_(z) {}              
};  

class PDB{

public:

	PDB() : force_span(false), polar_flg(false), output_all_chains(false), beta(false), flip(false), pre_pos(false), pre_flip(false), max_x(-10000.0), max_y(-10000.0), max_z(-10000.0), max_c_dist(0.0), x_rot_opt(0.0), y_rot_opt(0.0), z_trans_opt(0.0), y_rot_pre(0.0), z_rot_pre(0.0), z_trans_pre(0.0), best_extra_shift(0.0), best_cyto_shift(0.0) {};
	~PDB(){};
	int parse_pdb(std::string&, std::string&, std::vector<int>&);
	void write_pdb(double = 0.0);
	int parse_potential(std::string&);	
	double orientate(double, double, double);	
	int calc_helix_tilt(std::vector<int>&);
	void calc_thickness(double, double, double);
	int pre_position(std::string&, std::vector<int>&);
	void origin_shift();
	void transform_atom(Point3d&, double, double, double);
	void set_polar(bool t){polar_flg = t;};
	void set_output(const char* t){output = t;};
	void set_forcespan(bool t){force_span = t;};
	void set_allchains(bool t){output_all_chains = t;};
	void set_chain(std::string t){chains.push_back(t);};
	void set_optparams(double x, double y, double z){x_rot_opt = x; y_rot_opt = y; z_trans_opt = z;};
	void set_beta(bool t);
	void set_flip(bool t){flip = t;};
	std::string get_nterm(double, double, double);
	double get_maxcdist(){return max_c_dist;};
	double get_optx(){return x_rot_opt;};
	double get_opty(){return y_rot_opt;};
	double get_optz(){return z_trans_opt;};

private:
	int get_slice_index(double&);
	std::vector<std::string> chains;
	std::vector<std::string> all_atoms;
	std::vector<Point3d> backbone;
	std::string target, output;
	std::map <int, int> pdb_res_nums;
	bool force_span, polar_flg, output_all_chains ,beta, flip, pre_pos, pre_flip;
	double max_x, max_y, max_z, max_c_dist, x_rot_opt, y_rot_opt, z_trans_opt, y_rot_pre, z_rot_pre, z_trans_pre, best_extra_shift, best_cyto_shift;

};

#endif
