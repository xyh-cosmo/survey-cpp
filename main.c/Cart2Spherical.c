#include "SurveySim.h"

int main( int argc, char* argv[] ){
    if( argc != 3 ){
        printf("usage:\n\t%s ra_in dec_in\n",argv[0]);
        exit(0);
    }

    double ra_in = atof(argv[1]);
    double dec_in= atof(argv[2]);
 
    double x[3];
    x[0] = cos(dec_in*M_PI/180)*cos(ra_in*M_PI/180);
    x[1] = cos(dec_in*M_PI/180)*sin(ra_in*M_PI/180);
    x[2] = sin(dec_in*M_PI/180);

    double radec_out[2];

    Cartesian2Equatorial(x, radec_out);

    printf("Input:  ra = %.18g, dec = %.18g\n",ra_in,dec_in);
    printf("Output: ra = %.18g, dec = %.18g\n",radec_out[0],radec_out[1]);
    printf("diff_ra  = %.25g\ndiff_dec = %.25g\n",radec_out[0]-ra_in,radec_out[1]-dec_in);

    return 0;
}