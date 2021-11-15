#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "n_gaussian_raw.cpp"

#define element tl1plus_6
static char *label = "Tl+";
static int nelec = 80;

int main(){


int i,j;
double r, rho, rhoint;
double ap[6], bp[6];

for( j=0; j<6; j++){

    ap[j] = element[ 2*j ] * pow( 12.5663706 / element[ 2*j + 1 ], 1.5 );
    bp[j] = 157.91367 / element[ 2*j + 1 ];
    
}

printf( "# %s density, from 6-gaussian atomic form factor from cctbx\n", 
       label );
printf( "# total number of electrons:\n%d\n", nelec );
printf( "#r(A) e-density (A^-3)\n" );

r = 0.; rhoint = 0.;

for( i=0; i<160; i++ ){
    r += 0.025;
    rho = 0.;
    for( j=0; j<6; j++ ){
        rho += ap[j]* exp( -bp[j] * r * r / 4. );
    }
    rhoint += rho * r * r;
    printf( "%5.3f  %15.6g\n", r, rho );
}

rhoint *= 12.566370 * 0.025;
fprintf( stderr, "rough estimate of integrated density: %10.5f\n", rhoint);

}
