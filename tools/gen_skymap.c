#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#define NUM 1000000
#define PI_180 (M_PI/180.0)

void CoordinateSpin_x(double coorO[3], double coorR[3], double angle ) {

    double arcAngle = angle * PI_180;
    coorR[0] = coorO[0];
    coorR[1] = coorO[1] * cos(arcAngle) + coorO[2] * sin(arcAngle);
    coorR[2] = coorO[1] * (-sin(arcAngle)) + coorO[2] * cos(arcAngle);
}

void Cartesian2Equatorial(double* carCoor, double* eCoor) {
    double x1 = carCoor[0], x2 = carCoor[1], x3 = carCoor[2];
    double r = sqrt(x1*x1+x2*x2+x3*x3);
    double theta = asin(x3/r);

    *(eCoor+1) = theta*180*M_1_PI;
    *(eCoor+0) = atan(x2/(r*cos(theta)+x1)) *360*M_1_PI;

	*eCoor += (*eCoor < 0)*360;
}

/////////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] ){

	gsl_rng* rand_seed = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(rand_seed,time(NULL));

	int i,j;
	double u;
	double phi, theta;
	double x, y, z;
	double b;

    double beta, delta_beta;
    double bgal, delta_bgal;
    
    if( argc != 5 ){
        printf("usage: %s beta delta_beta bgal delta_bgal\n",argv[0]);
        exit(0);
    }
    
    beta = atof(argv[1]);
    delta_beta = atof(argv[2]);
    bgal = atof(argv[3]);
    delta_bgal = atof(argv[4]);    
/*    scanf("%lf %lf %lf %lf",&beta,&delta_beta,&bgal,&delta_bgal);*/
    printf("beta = %g, delta_beta = %g\n",beta,delta_beta);
    printf("bgal = %g, delta_bgal = %g\n",bgal,delta_bgal);

	FILE *fp1 = fopen("area_new.txt","w");
	FILE *fp2 = fopen("area_normal.txt","w");

    int n_good=0;
    int cnt=0;
	while( cnt < NUM ){
		phi 	= gsl_ran_flat(rand_seed, 0.0, 2.0*M_PI);
		u 		= gsl_ran_flat(rand_seed, -1.0, 1.0);
		theta 	= acos(u) - M_PI_2;

		double t_coor[3] = {cos(theta)*cos(phi), cos(theta)*sin(phi), sin(theta)};
		double r_coor[3], ll_coor[3];
		CoordinateSpin_x(t_coor, r_coor, -23.4522);
		Cartesian2Equatorial(r_coor, ll_coor);

		b = asin(-0.8676660 * cos(ll_coor[0] * PI_180) * cos(ll_coor[1] * PI_180)
				- 0.1980764 * sin(ll_coor[0] * PI_180) * cos(ll_coor[1] * PI_180)
				+ 0.4559840 * sin(ll_coor[1] * PI_180) ) / PI_180;

		theta *= (180*M_1_PI);
		int flag1 = (int)((fabs(theta) > beta && fabs(theta)<75) && (fabs(b) > bgal));
		int flag2 = (int)( (fabs(theta) < beta+delta_beta) && (fabs(b) < bgal+delta_bgal) );
		int flag = (int)(flag1 || flag2);
		
		if( flag == 1 ){
		    if( flag2 == 1 ){
		        fprintf(fp1,"%f %f\n",phi*180*M_1_PI,theta);
		    }
		    else if( flag1 == 1 ){
		        fprintf(fp2,"%f %f\n",phi*180*M_1_PI,theta);
		    }
		    n_good++;
		}

		cnt++;
	}

    double area = 4*M_PI*pow(180/M_PI,2)*n_good/NUM;
    printf("survey area: %f\n",area);

	fclose(fp1);
	fclose(fp2);

	gsl_rng_free(rand_seed);
	return 0;
}
