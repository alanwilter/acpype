#ifndef GRID_HPP
#define GRID_HPP

#include <pdb.hpp>	
#include <boost/thread.hpp>  

class Grid{
public:
	Grid() : threads(1), ncalls(0), prevbest(1000000.0) {};
	~Grid(){};
	void set_target(PDB* t){protein = t;};
	void set_threads(int t){threads = t;};
	void set_zrange(double t){z_range = t;};
	int run();

private:
	void exhaustive_search(double, double, int);
	int threads, ncalls;
	double z_range, prevbest;
	std::string n;
	boost::thread_group g;
	boost::mutex result_mutex;
	PDB* protein;
};

#endif
