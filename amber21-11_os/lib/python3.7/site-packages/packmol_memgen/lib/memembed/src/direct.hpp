#ifndef DIRECT_HPP
#define DIRECT_HPP

#include <pdb.hpp>	

class Pattern{
public:
	Pattern() : funevals(0) {};
	~Pattern(){};
	void set_target(PDB* t){protein = t;};
	int hooke(int, double*, double*, double, double, int);

private:
	double best_nearby(double*, double*, double, int);
	double func(double*);
	int funevals;
	PDB* protein;
};

#endif
