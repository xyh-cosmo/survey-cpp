#include "SurveySim.h"


int main(int argc, char* argv[]){

	if( argc < 5 ){
		printf("usage: %s ra_old dec_old ra_new dec_new\n",argv[0]);
		exit(0);
	}

	double ra_old = atof(argv[1]);
	double dec_old= atof(argv[2]);
	double ra_new = atof(argv[3]);
	double dec_new= atof(argv[4]);
	double angle_rot1,angle_rot2;

	Get_RotationAngle_faster( ra_old,dec_old,
							   ra_new,dec_new,&angle_rot1);

	Get_RotationAngle_faster2( ra_old,dec_old,
							   ra_new,dec_new,&angle_rot2,0);

	printf("========================================================\n");
	printf("ra_old = %10.6f dec_old = %10.6f\nra_new = %10.6f dec_new = %10.6f\n",
			ra_old,dec_old,ra_new,dec_new);

	printf("x0 = %10.6f y0 = %10.6f z0 = %10.6f\n",
					sin(dec_old*PI_180)*cos(ra_old*PI_180),
					sin(dec_old*PI_180)*sin(ra_old*PI_180),
					cos(dec_old*PI_180));
	printf("x1 = %10.6f y1 = %10.6f z1 = %10.6f\n",
					sin(dec_new*PI_180)*cos(ra_new*PI_180),
					sin(dec_new*PI_180)*sin(ra_new*PI_180),
					cos(dec_new*PI_180));

	printf("rotation angle1 = %10.6f\n",angle_rot1);
	printf("rotation angle2 = %10.6f\n",angle_rot2);

	return 0;
}