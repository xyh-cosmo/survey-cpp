//	Test angle calculation between panel normal vector and the position vector of Sun.
//
//	Method: place the Sun at (1,0,0), and choose several specific sky-directions for the telescope;
//		for these choose disrections, the angles can be found by hand.

#include "SurveySim.h"

int main( int argc, char* argv[] ){

	if( argc < 3 ){
		printf("usage:\n\t %s ra dec\n",argv[0]);
		printf("-------------------------------\n");
		printf("Note: sun is located at (1,0,0)\n");
		printf("-------------------------------\n");
		exit(0);
	}
	

	double ra = atof(argv[1]);
	double dec= atof(argv[2]);

	printf(" (ra, dec) = (%g, %g)\n",ra,dec);

//	Sun's position and its distance
	double sun[3] = {1,0,0};
	double dist_sun = 1.0;

//	transform ra,dec to the unit vector that points at it
	double phi = ra*PI_180;
	double theta = dec*PI_180;
	
	double x = cos(theta)*cos(phi);
	double y = cos(theta)*sin(phi);
	double z = sin(theta);
	double p[3] = {x,y,z};

//	get the panel rotation axis A
	double ez[3] = {0,0,1};
	
	if( fabs(dec) > 89.99 ){
		ez[0] = cos(phi);
		ez[1] = sin(phi);
		ez[2] = 0;
	}
	
	double n[3];
	CrossProduct(p,ez,n);

	int i=0;
	double n_norm = VecNorm(n);
	for(i=0;i<3;i++){
		n[i] /= n_norm;
	}

//	call TestPanelAngle(***)
	double cosval=-999;
	TestPanelAngle(sun,dist_sun,x,y,z,n[0],n[1],n[2],&cosval,25,0);

	double ang = acos(cosval)*180/M_PI;

	printf("Angle = %g\n",ang);

	return 0;
}