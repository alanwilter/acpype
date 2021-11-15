/*
	Knowledge-based membrane protein positioning
	(c) Timothy Nugent 2013
*/
#include <pdb.hpp>	
#include <ga.hpp>
#include <grid.hpp>
#include <direct.hpp>
#include <stdio.h>
#include <math.h>
#include <string>
#include <iostream>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>

using namespace std;
using namespace boost;
using namespace boost::accumulators;

#define VERSION 1.15
#define _USE_MATH_DEFINES

void usage(const char* progname){
	cout << endl << "Knowledge-based membrane protein positioning" << endl;
	cout << "(c) Timothy Nugent 2012, Version " << VERSION << endl << endl;
	cout << "Usage : " << progname << " [options] <input pdb file>" << endl << endl;
	cout << "Options:" << endl << endl;
	cout << "-o <output>\tOutput pdb file. Default <input pdb file>_EMBED.pdb" << endl;
	cout << "-c <chains>\tComma separated list of transmembrane chains. Default all" << endl;
	cout << "-x\t\tOutput all chains regardless of chains passed by -c parameter." << endl;	
	cout << "-m <file>\tUse alternative potential file." << endl;
	cout << "-t <string>\tComma seperated list of helix boundaries used to calculate helix tilt angles." << endl;
	cout << "-f <string>\tComma seperated list of residues to skip." << endl;
	cout << "-r <string>\tChain that topology refers to. Default 'A'" << endl;
	cout << "-v <int>\tMaximum Cb-Cb distance to allow." << endl;		
	cout << "-a <int>\tThreads to use." << endl;	
	cout << "-s <int>\tSearch type. 0 = Genetic algorithm, 1 = Grid, 2 = Direct, 3 = GA repeated 5 times. Default 0." << endl;	
	cout << "-q <int>\tOptimise membrane thickness. 0 = Do not optimise, 1 = After orientation, 2 = Do not orientate. Default 0." << endl;	
	cout << "-n <in|out>\tLocation of N-terminal (first residue of first chain)." << endl;
	cout << "-b\t\tTarget is a beta-barrel." << endl;
	cout << "-l\t\tForce target to span membrane." << endl;
	cout << "-p\t\tDraw lines representing polar head groups." << endl;	
	cout << "-e\t\tJust compute energy." << endl;
	cout << "-z\t\tJust compute helix tilt angles." << endl;
	cout << "-h\t\tDisplay usage." << endl << endl;
}	

int main(int argc, const char* argv[]){

	int threads = 1, i = 1, e = 0, s = 0, calc_tilt_only = 0, calc_hydro_thickness = 0, maxcc = 5000;
	string target, mempotfile, chain_string, n_term, topology_string, residue_skip_string;
	string top_chain = "A";
	vector<int> topology,skip_residues;
	vector<double> x_rot_values,y_rot_values,z_trans_values;

	PDB* protein = new PDB;

	if (argc < 2){
		usage(argv[0]);
		delete protein;
		return(1);
	}
	while(i < argc){
  		if( argv[i][0] == '-'){
    			i++;
    			switch(argv[i-1][1]){
				case 'h' : {usage(argv[0]);delete protein;return(1);}
				case 'e' : {e = 1; i--; break;}	
				case 'n' : {n_term = argv[i];; break;}	   
				case 'a' : {threads=atoi(argv[i]); break;}
				case 'v' : {maxcc=atoi(argv[i]); break;}
				case 's' : {s = atoi(argv[i]); break;}
				case 'b' : {protein->set_beta(true); i--; break;}
				case 't' : {topology_string=argv[i]; break;}
				case 'f' : {residue_skip_string=argv[i]; break;}
				case 'p' : {protein->set_polar(true); i--; break;}
				case 'o' : {protein->set_output(argv[i]); break;}	
				case 'm' : {mempotfile = argv[i]; break;}
				case 'x' : {protein->set_allchains(true); i--; break;}
				case 'l' : {protein->set_forcespan(true); i--; break;}
				case 'z' : {calc_tilt_only = 1; i--; break;}
				case 'q' : {calc_hydro_thickness = atoi(argv[i]); break;}
				case 'r' : {top_chain = argv[i]; break;}
				case 'c' : {chain_string = argv[i]; break;}
   				default  : {usage(argv[0]);delete protein;return(1);}
			}   
   		}
    		i++;
  	}
	target = argv[i-1];

	cout << "Membrane protein orientation using a knowledge-based statistical potential" << endl;
	cout << "(c) Tim Nugent 2013" << endl;
	printf("Version %1.2f\n\n",VERSION);

	if(!e && chain_string.length()){
		cout << "Processing chain(s) " << chain_string << endl;
		for(unsigned int l = 0; l < chain_string.length(); l++){
			if(chain_string.substr(l,1) != ","){
				protein->set_chain(chain_string.substr(l,1));
			}
		}
	}
	if(!e && topology_string.length()){
  		typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
 		boost::char_separator<char> sep(",-");
  		tokenizer tokens(topology_string, sep);
  		for(tokenizer::iterator tok_iter = tokens.begin(); tok_iter != tokens.end(); ++tok_iter){
  			int boundary = boost::lexical_cast<int>(*tok_iter);
  			topology.push_back(boundary);
  		}
	}
	if(!e && residue_skip_string.length()){
  		typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
 		boost::char_separator<char> sep(",-");
  		tokenizer tokens(residue_skip_string, sep);
  		for(tokenizer::iterator tok_iter = tokens.begin(); tok_iter != tokens.end(); ++tok_iter){
  			int res = boost::lexical_cast<int>(*tok_iter);
  			skip_residues.push_back(res);
  		}
	}
	if(target.length() < 5){
		cout << "Doesn't look like a PDB file." << endl;
		delete protein;
		return(1);	
	}		

	// Parse potential file
	if(!protein->parse_potential(mempotfile)){
		delete protein;
		return(1);
	}

	// Parse PDB file
	if(!protein->parse_pdb(target,top_chain,skip_residues)){
		delete protein;
		return(1);
	}

	if(calc_tilt_only){
		if(topology.size()){
			protein->calc_helix_tilt(topology);
		}else{
			cout << "No topology provided." << endl;
		}
		delete protein;
		return(1);	
	}
	
	if(calc_hydro_thickness == 2){
		protein->calc_thickness(0.0,0.0,0.0);
		protein->write_pdb();
		delete protein;
		return(0);	
	}

	// Fit chains into membrane
	if(!e){

		double initial_e = protein->orientate(0.0,0.0,0.0);
		protein->origin_shift();
		if(protein->get_maxcdist() > maxcc){
			cout << "Max C-C distance is > " << maxcc << " - aborting." << endl << endl;
			delete protein;
			return(1);
		}

		// Set search parameters
		vector<double> minparam(3), maxparam(3);

		minparam[0] = 0.0;
		maxparam[0] = 2*M_PI;
		minparam[1] = 0.0;
		maxparam[1] = 2*M_PI;
		minparam[2] = -15.0;
		maxparam[2] = protein->get_maxcdist()+15;

		if(topology.size() && n_term != "" && !calc_tilt_only){
			cout << "Pre-positioning using topology..." << endl;
			if(!protein->pre_position(n_term,topology)){
				cout << "Couldn't find helix boundaries for chain " << top_chain << endl;
				delete protein;
				return(1);
			}
		}

		int ncalls = 0;

		// Grid search
		if(s == 1){
			if(threads > 1){
				cout << "Performing grid search using " << threads << " threads...\n";
			}else{	
				cout << "Performing grid search...\n";
			}
			Grid* gs = new Grid();
			gs->set_target(protein);
			gs->set_threads(threads);
			gs->set_zrange((maxparam[2]-minparam[2])/threads);
			ncalls = gs->run();
			delete gs;

		// Use Hooke pattern search	
		}else if(s == 2){			
			cout << "Performing direct search...\n";			
		   	double startpt[3], endpt[3];
		   	int itermax = 5000;
		   	double rho = 0.9;
		   	double epsilon = 1E-6;
		   	double prevbest = 1000000;
		   	// Quick grid search to get starting values
			for(double x = 0; x < maxparam[0]; x+= 0.104719755){
				for(double y = 0; y < maxparam[0]; y+= 0.1047197555){
					for(double z = minparam[2]; z < maxparam[2]; z+= 5){
						double lowest_e = protein->orientate(x,y,z);
						if(lowest_e < prevbest){
							prevbest = lowest_e;
						    startpt[0] = x;
						    startpt[1] = y;
						    startpt[2] = z;
						}	
					}
				}
			}
			Pattern* pa = new Pattern();
			pa->set_target(protein);
		    ncalls = pa->hooke(3, startpt, endpt, rho, epsilon, itermax);
		    protein->set_optparams(endpt[0],endpt[1],endpt[2]);
		    delete pa;

		// Genetic algorithm    
		}else{

			int repeats = 1;
			if(s == 3) repeats = 5;
			for(int r = 0; r < repeats; r++){

				if(threads > 1){
					cout << "Performing GA search using " << threads << " threads...\n";
				}else{	
					cout << "Performing GA search search...\n";
				}

				GA* ga = new GA(minparam, maxparam);
				ga->set_threads(threads);
				ga->set_poolsize(10000);
				ga->set_maxcalls(1000000);
				ga->set_target(protein);
				vector<double> results = ga->run_ga();
				ncalls += ga->get_calls();
				delete ga;

				if(s == 3) cout << endl;

				x_rot_values.push_back(results[0]);
				y_rot_values.push_back(results[1]);
				z_trans_values.push_back(results[2]);
				protein->set_optparams(results[0],results[1],results[2]);

			}	

			double e_big = 1000000;
			for(int r = 0; r < repeats; r++){

				double lowest_rep_e = protein->orientate(x_rot_values[r],y_rot_values[r],z_trans_values[r]);
				if(lowest_rep_e < e_big){					
					e_big = lowest_rep_e;
					protein->set_optparams(x_rot_values[r],y_rot_values[r],z_trans_values[r]);
				}

			}	
			accumulator_set<double, stats<tag::variance> > accx, accy, accz; 
			for(int r = 0; r < repeats; r++){

				x_rot_values[r] = x_rot_values[r]*(180/M_PI);
				y_rot_values[r] = y_rot_values[r]*(180/M_PI);
				if(x_rot_values[r] > 180){
					x_rot_values[r] -= 180;
				}
				if(y_rot_values[r] > 180){
					y_rot_values[r] -= 180;
				}
				accx(x_rot_values[r]);
				accy(y_rot_values[r]);
				accz(z_trans_values[r]);

			}	
			if(s == 3){
    			
    			cout << "X-axis rotation from starting position (degrees):" << endl;  
    			double t =  protein->get_optx() *(180/M_PI);
    			if(t > 180){
    				t -= 180;
    			}
    			cout << "Lowest energy:\t" << t << endl;     
    			cout << "Mean:\t\t" << mean(accx) << endl; 
    			cout << "St. error:\t" << sqrt(variance(accx)) << endl;  
    			cout << "Y-axis rotation from starting position (degrees):" << endl; 
    			t =  protein->get_opty() *(180/M_PI);
    			if(t > 180){
    				t -= 180;
    			}
    			cout << "Lowest energy:\t" << t << endl;          
    			cout << "Mean:\t\t" << mean(accy) << endl; 
    			cout << "St. error:\t" << sqrt(variance(accy)) << endl;    
    			cout << "Z-axis translation from starting position (angstroms):" << endl;
    			cout << "Lowest energy:\t" << protein->get_optz() << endl;            
    			cout << "Mean:\t\t" << mean(accz) << endl; 
    			cout << "St. error:\t" << sqrt(variance(accz)) << endl << endl;     			
			}
		}	

		double lowest_e = protein->orientate(protein->get_optx(),protein->get_opty(),protein->get_optz());
		if(initial_e <= lowest_e){
			cout << "Couldn't find lower energy than starting orientation." << endl;
			lowest_e = initial_e;
		}

		string n = protein->get_nterm(protein->get_optx(),protein->get_opty(),protein->get_optz());
		if(n_term != "" && n != n_term){
			cout << "Inverting to satisfy N-terminal constraint..." << endl;
			protein->set_flip(true);
		}

		cout << "Total function calls:\t" << ncalls << endl;	
		cout << "Final Energy:\t\t" << setprecision(10) << lowest_e << endl;

		// Calculate membrane thickness
		if(calc_hydro_thickness){
			protein->calc_thickness(protein->get_optx(),protein->get_opty(),protein->get_optz());
		}

		protein->write_pdb(lowest_e);

	// Just compute energy
	}else{

		double lowest_e = protein->orientate(0.0,0.0,0.0);
		cout << setprecision(10) << lowest_e << endl;
	}

	delete protein;
	return(0);
}
