#include "SurveySim.h"

int main( int argc, char *argv[] ){

	if( argc < 2 ){
		printf("usage: %s angle (in degree)\n",argv[0]);
		exit(0);
	}

#if defined(_SLEW_TIME_GSL_INTERP_)
    init_TransTime();
#endif

	double ang = atof(argv[1]);
	double t = calculateTransTime(ang);
	printf("Transition time for %8.5f deg rotation is %8.5f secs\n",ang,t);

#if defined(_SLEW_TIME_GSL_INTERP_)
    free_TransTime();
#endif

	return 0;
}