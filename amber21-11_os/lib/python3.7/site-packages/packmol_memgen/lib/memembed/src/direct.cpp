#include <pdb.hpp>	
#include <direct.hpp>
#include <math.h> 

/* Calculate energy of structure */
double Pattern::func(double *x){

	double e = protein->orientate(x[0],x[1],x[2]);
	funevals++;
    return e;

}

/* given a point, look for a better one nearby, one coord at a time */
double Pattern::best_nearby(double *delta, double *point, double prevbest, int nvars){

	   double* z = new double[nvars];
	   double minf, ftmp;
	   int i;

	   minf = prevbest;
	   for (i = 0; i < nvars; i++){
		   z[i] = point[i];
	   }
	   for (i = 0; i < nvars; i++){
		   z[i] = point[i] + delta[i];
		   ftmp = func(z);
		   if (ftmp < minf){
			   minf = ftmp;
		   }else{
			   delta[i] = 0.0 - delta[i];
			   z[i] = point[i] + delta[i];
			   ftmp = func(z);
			   if (ftmp < minf){
				   minf = ftmp;
			   }else{
				   z[i] = point[i];
			   }
		   }
	   }
	   for (i = 0; i < nvars; i++){
		   point[i] = z[i];
	   }

	   delete [] z;

	   return (minf);
}


int Pattern::hooke(int nvars, double *startpt, double *endpt, double rho, double epsilon, int itermax){
    
    double newf, fbefore, steplength, tmp, last_e = 1000000;
    double* delta = new double[nvars];
    double* xbefore = new double[nvars];
    double* newx = new double[nvars];
    int i, keep, iters, iadj;
    
    for (i = 0; i < nvars; i++){
		newx[i] = xbefore[i] = startpt[i];
		delta[i] = fabs(startpt[i] * rho);
		if (delta[i] == 0.0){
		    delta[i] = rho;
		}
    }

    iadj = 0;
    steplength = rho;
    iters = 0;
    fbefore = func(newx);
    newf = fbefore;
    while ((iters < itermax) && (steplength > epsilon)){
		iters++;
		iadj++;

		if(fbefore < last_e){
			printf("Lowest energy %f after %d function evaluations.\n", fbefore, funevals);
		}
		last_e = fbefore;
		/* find best new point, one coord at a time */
		for (i = 0; i < nvars; i++){
		    newx[i] = xbefore[i];
		}
		newf = best_nearby(delta, newx, fbefore, nvars);
		/* if we made some improvements, pursue that direction */
		keep = 1;
		while ((newf < fbefore) && (keep == 1)) {
		    iadj = 0;
		    for (i = 0; i < nvars; i++) {
				/* firstly, arrange the sign of delta[] */
				if (newx[i] <= xbefore[i])
				    delta[i] = 0.0 - fabs(delta[i]);
				else
				    delta[i] = fabs(delta[i]);
				/* now, move further in this direction */
				tmp = xbefore[i];
				xbefore[i] = newx[i];
				newx[i] = newx[i] + newx[i] - tmp;
		    }
		    fbefore = newf;
		    newf = best_nearby(delta, newx, fbefore, nvars);
		    /* if the further (optimistic) move was bad.... */
		    if (newf >= fbefore)
			break;
		    /* make sure that the differences between the new */
		    /* and the old points are due to actual */
		    /* displacements; beware of roundoff errors that */
		    /* might cause newf < fbefore */
		    keep = 0;
		    for(i = 0; i < nvars; i++){
				keep = 1;
				if(fabs(newx[i] - xbefore[i]) > (0.5 * fabs(delta[i]))){
				    break;
				}else{
				    keep = 0;
				}
		    }
		}
		if((steplength >= epsilon) && (newf >= fbefore)){
		    steplength = steplength * rho;
		    for (i = 0; i < nvars; i++){
				delta[i] *= rho;
		    }
		}
    }
    for (i = 0; i < nvars; i++){
		endpt[i] = xbefore[i];
    }

    delete [] delta;
    delete [] xbefore;
    delete [] newx;

    return(funevals);
}
