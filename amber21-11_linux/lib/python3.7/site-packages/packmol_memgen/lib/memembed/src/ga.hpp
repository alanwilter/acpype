#ifndef GA_HPP
#define GA_HPP

#include <pdb.hpp>	
#include <boost/thread.hpp>  

struct Schema{

	Schema() : perfval(0.0), selval(0.0), evalflg(true) {}
	std::vector<double> genome;
	double perfval, selval;
	bool evalflg;
	
};

class GA{

public:

	GA(std::vector<double>&, std::vector<double>&);
	~GA();
	std::vector<double> run_ga();
	void set_threads(int);
	void set_poolsize(int);
	void set_maxcalls(int);
	void set_target(PDB*);
	int get_calls();

private:	
	void randomise();
	void init();
	void statistics(Schema*);
	void threaded_eval(int, int, int, Schema*);
	void crossovr();
	void mutate(double);
	void gaselect();
	void sortpool(Schema*);
	double JKISS();
	double gaussrnd();
	double eval(std::vector<double>&, int);
	static bool schcmp(const Schema&, const Schema&);
	std::vector<int> samparr;
	std::vector<double> minparam, maxparam, results;
	int	ncalls, threads, genlen, poolsize, prevbest, vbest, max_ga_calls, besti, totalcalls;
	unsigned int w, x, y, z, c;
	double mutrate, crosrate, mutscfac, worst, best, avc_perf;
	Schema *curpool, *newpool;		
	boost::mutex result_mutex;
	PDB* protein;
};

#endif
