#include <pdb.hpp>	
#include <ga.hpp>	
#include <random>
#include <iostream>

using namespace std;

GA::GA(vector<double>& minp, vector<double>& maxp){

	threads = 1;
	minparam = minp;
	maxparam = maxp;
	genlen = minparam.size();
	max_ga_calls = 1000000;
	poolsize = 10000;
	mutrate = 0.1;
	crosrate = 0.9;
	mutscfac = 1;
	prevbest = 100000;
	vbest = 100000;
	totalcalls = 0;
	ncalls = 0;
	for(int i = 0; i < poolsize; i++){
		samparr.push_back(0.0);
	}
	randomise();
	init();
}

GA::~GA(){

    	delete [] curpool;
    	delete [] newpool;

}

void GA::set_threads(int t){

	threads = t;

}	

void GA::set_poolsize(int t){

	poolsize = t;

}	

void GA::set_maxcalls(int t){

	max_ga_calls = t;

}	

void GA::set_target(PDB* t){

	protein = t;

}	

int GA::get_calls(){

	return(ncalls);

}	

void GA::randomise(){

	std::random_device rd("/dev/urandom");
	x = rd();
	while (!(y = rd()));
	z = rd();
	c = rd() % 698769068 + 1;

}

void GA::init(){

	curpool = new Schema[poolsize];
	newpool = new Schema[poolsize];
	for (int i = 0; i < poolsize; i++){
		for (int j = 0; j < genlen; j++){
			curpool[i].genome.push_back(JKISS()*(maxparam[j]-minparam[j])+minparam[j]);
		}
	}
}

vector<double> GA::run_ga(){

	int gen, prevgen = 0;
    	Schema *temp;
    	statistics(curpool);

    	for (gen = 1;; gen++){
		gaselect();
		crossovr();
		mutate(mutrate);	    
		statistics(newpool);
		if (best < prevbest){
	    	prevbest = best;
	    	prevgen = gen;
		}
		if(gen - prevgen > 50 || worst == best || ncalls > max_ga_calls){
	    	cout << "*** Convergence detected!" << endl;
	    	return(results);
		}	
		temp = newpool;
		newpool = curpool;
		curpool = temp;
    	}
}

void GA::statistics(Schema * pool){

	// Split up pool - assign ppt to each thread
	double ppt = 1+(int)poolsize/threads;
	boost::thread_group *g = new boost::thread_group();

	int first_pool_index = 0;
	for(int t = 0; t < threads; t++){
		int last_pool_index = (int)ppt*(t+1);
		if(last_pool_index > poolsize){
			last_pool_index = poolsize;
		}	
		boost::thread *t1 = new boost::thread(boost::bind(&GA::threaded_eval,this,t,first_pool_index,last_pool_index,pool));
		g->add_thread(t1);
		first_pool_index = (int)ppt*(t+1)+1;
	}	
	// Wait for threads to join
	g->join_all();
	delete g;

	// When thread_group g is deleted,
	// all threads should get deleted.
	// Valgrind complains they're not but 
	// this appears to be harmless, e.g.:
	// http://www.boost.org/doc/libs/1_54_0/doc/html/thread/thread_management.html#thread.thread%5Fmanagement.threadgroup.destructor
	// http://boost.2283326.n4.nabble.com/thread-Memory-leak-in-Boost-Thread-td2648030.html

    	avc_perf = best = worst = pool[0].perfval;

    	for (int i = 1; i < poolsize; i++){
		avc_perf += pool[i].perfval;
		if (pool[i].perfval > worst){
	   		worst = pool[i].perfval;
		}
		if (pool[i].perfval < best){
	    	best = pool[i].perfval;
		}
    	}
    	avc_perf /= (double)poolsize;
}

void GA::threaded_eval(int t, int start, int stop, Schema * pool){

	for (int i = start; i < stop; i++){
		if(pool[i].evalflg){
			pool[i].perfval = eval(pool[i].genome,t);
			pool[i].evalflg = false;
		}
	} 	
}

/* Randomly crossover pool */
void GA::crossovr(){

	vector<double> old1, old2;
    	for (int i = 0; i < poolsize - 1; i += 2){
		if(JKISS() < crosrate){

	    		/* Multi point crossover */
			old1 = newpool[i].genome;
			old2 = newpool[i+1].genome;

	    		for (int p1 = 0; p1 < genlen; p1++){
				if (JKISS() < 0.5){
			    		newpool[i].genome[p1] = newpool[i+1].genome[p1];
			    		newpool[i+1].genome[p1] = old1[p1];
				}
			}	
			if(newpool[i].genome != old1){
				newpool[i].evalflg = true;
			}
			if(newpool[i+1].genome != old2){
				newpool[i+1].evalflg = true;
			}
		}		
	}	
}

void GA::mutate(double prob){

	double delta;
	static double scale = 0.25;

    	if (prob > 0.0){
		for (int i = 0; i < poolsize; i++){
			if (i != besti){
				for (int j = 0; j < genlen; j++){
		   			if (JKISS() < prob){
						delta = scale * gaussrnd() * (maxparam[j] - minparam[j]);
						newpool[i].genome[j] += delta;
						if (newpool[i].genome[j] > maxparam[j]){
				    			newpool[i].genome[j] = maxparam[j];
						}
						if (newpool[i].genome[j] < minparam[j]){
				    			newpool[i].genome[j] = minparam[j];
				    		}	
				    		newpool[i].evalflg = true;
			    		}
				}
			}	
		}
	}	    
    	/* Apply mutation scaling factor */
    	scale *= mutscfac;
}

/* Select new population from old */
void GA::gaselect(){

	double ptr; /* determines fractional selection */
	double sum; /* control for selection loop */
	double fitsum; /* sum of fitness values */
	int i, k;

	sortpool(curpool);

	/* denominator for ordinal selection probabilities */
    	for (fitsum = i = 0; i < poolsize; i++){
		curpool[i].selval = poolsize - i;
		fitsum += curpool[i].selval;
		if (curpool[i].perfval == best){
	    	besti = i;
		}
	}


	for (i = 0; i < poolsize; i++){
		sum = 0.0;
		k = -1; /* index of next Selected structure */
		ptr = fitsum * JKISS();	/* spin the wheel one time */
		do{
	    	k++;
	    	sum += curpool[k].selval;
		}while (sum < ptr && k < poolsize - 1);
		samparr[i] = k;
	}

    	/* Form the new population */
    	for (i = 0; i < poolsize; i++){
		k = samparr[i];
		newpool[i].genome = curpool[k].genome;
		newpool[i].perfval = curpool[k].perfval;
		newpool[i].evalflg = false;
    	}
}

bool GA::schcmp (const Schema& sch1, const Schema& sch2){

	return(sch1.perfval < sch2.perfval);

}

void GA::sortpool(Schema * pool){

	sort(pool, pool+poolsize, schcmp);

}

/* Generate gaussian deviate with mean 0 and stdev 1 */
double GA::gaussrnd(){

	double x1, x2, w;
 
	do{
		x1 = 2.0 * JKISS() - 1.0;
		x2 = 2.0 * JKISS() - 1.0;
		w = x1 * x1 + x2 * x2;
    	}while (w >= 1.0);

    	w = sqrt((-2.0 * log(w)) / w);
   	return x1 * w;
}

/* Evaluate given set of parameter values */
double GA::eval(vector<double>& t, int thread){

	double v = protein->orientate(t[0],t[1],t[2]);
	boost::mutex::scoped_lock lock(result_mutex);	

	ncalls++;
	if (v < vbest){
		if(threads > 1){
			printf("Lowest energy %f after %d function evaluations (thread %d).\n", v, ncalls, thread+1);
		}else{
			printf("Lowest energy %f after %d function evaluations.\n", v, ncalls);
		}
		results = t;
		vbest = v;
	}
    	return v;
}

double GA::JKISS(){
	
	unsigned long long t;
	x = 314527869 * x + 1234567;
	y ^= y << 5; y ^= y >> 7; y ^= y << 22;
	t = 4294584393ULL * z + c; c = t >> 32; z = t;
	return (double)(x + y + z)/4294967296.0;
	
}
