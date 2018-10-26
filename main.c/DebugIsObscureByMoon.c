#include "SurveySim.h"

int main( int argc, char* argv[] ){
	if( argc < 7 ){
		printf("\nusage: %s moon_x moon_y moon_z sky_x sky_y sky_z\n",argv[0]);
		exit(0);
	}

	double moon_x = atof(argv[1]);
	double moon_y = atof(argv[2]);
	double moon_z = atof(argv[3]);

	double px = atof(argv[4]);
	double py = atof(argv[5]);
	double pz = atof(argv[6]);

	double moon[3] = {moon_x, moon_y, moon_z};

	IsObscureByMoon(moon,px,py,pz);
}