#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#define NUM 100000000
#define PI_180 (M_PI/180.0)
#define AREA_MIN 17500

#define B_MIN 10
#define B_MAX 20
#define DELTA_B 1

#define BETA_MIN 7.5
#define BETA_MAX 15.5
#define DELTA_BETA 1

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


double fun_area(double beta, double bgal, gsl_rng* rand_seed){
	FILE* fp = fopen("xxxx.txt","w");
	double theta,b;
	double cnt = 0, n_good = 0.0, n_bad = 0.0;
    while( cnt < NUM ){
	    double phi 	= gsl_ran_flat(rand_seed, 0.0, 2.0*M_PI);
	    double u 	= gsl_ran_flat(rand_seed, -1.0, 1.0);
	    double theta= acos(u) - M_PI_2;

	    double t_coor[3] = {cos(theta)*cos(phi), cos(theta)*sin(phi), sin(theta)};//  黄道坐标系里的xyz
	    double r_coor[3], ll_coor[3];
	    CoordinateSpin_x(t_coor, r_coor, -23.4522);
	    Cartesian2Equatorial(r_coor, ll_coor);

	    b = asin( -0.8676660 * cos(ll_coor[0] * PI_180) * cos(ll_coor[1] * PI_180)
			    - 0.1980764 * sin(ll_coor[0] * PI_180) * cos(ll_coor[1] * PI_180)
			    + 0.4559840 * sin(ll_coor[1] * PI_180) );

		b *= (180*M_1_PI);
	    theta *= (180*M_1_PI);
	    int flag = (int)( (fabs(theta) > beta && fabs(theta) < 75 ) && (fabs(b) > bgal));
	    switch ( flag ) {
		    case 1:
			    n_good++;
			    fprintf(fp,"%6.3f %6.3f\n",phi*180/M_PI,theta);
			    break;
		    case 0:
			    n_bad++;
			    break;
	    }

	    cnt++;
    }

    double area = 4*M_PI*pow(180/M_PI,2)*n_good/NUM;
    
    fclose(fp);
    
    return area;
}

double fun_area2(double beta, double bgal, gsl_rng* rand_seed, double *s1, double *s2, double *s3, double *s4){
	double theta,b;
	double cnt = 0, n_good = 0.0, n_bad = 0.0;
	double n_good1=0, n_good2=0, n_good3=0, n_good4=0;
    while( cnt < NUM ){
	    double phi 	= gsl_ran_flat(rand_seed, 0.0, 2.0*M_PI);
	    double u 	= gsl_ran_flat(rand_seed, -1.0, 1.0);
	    double theta= acos(u) - M_PI_2;

	    double t_coor[3] = {cos(theta)*cos(phi), cos(theta)*sin(phi), sin(theta)};  //	榛勯亾鍧愭爣绯婚噷鐨剎yz
	    double r_coor[3], ll_coor[3];
	    CoordinateSpin_x(t_coor, r_coor, -23.4522);
	    Cartesian2Equatorial(r_coor, ll_coor);

	    b = asin( -0.8676660 * cos(ll_coor[0] * PI_180) * cos(ll_coor[1] * PI_180)
			    - 0.1980764 * sin(ll_coor[0] * PI_180) * cos(ll_coor[1] * PI_180)
			    + 0.4559840 * sin(ll_coor[1] * PI_180) );

		b *= (180*M_1_PI);
	    theta *= (180*M_1_PI);
	    int flag = (int)( (fabs(theta) > beta ) && (fabs(b) > bgal));
	    switch ( flag ) {
		    case 1:
			    n_good++;
			    if(theta > beta && b > bgal) n_good1++;
			    if(theta <-beta && b > bgal) n_good2++;
			    if(theta > beta && b <-bgal) n_good3++;
			    if(theta <-beta && b <-bgal) n_good4++;
			    
			    break;
		    case 0:
			    n_bad++;
			    break;
	    }

	    cnt++;
    }

    double area = 4*M_PI*pow(180/M_PI,2)*n_good/NUM;
    *s1 = 4*M_PI*pow(180/M_PI,2)*n_good1/NUM;
   	*s2 = 4*M_PI*pow(180/M_PI,2)*n_good2/NUM;
    *s3 = 4*M_PI*pow(180/M_PI,2)*n_good3/NUM;
    *s4 = 4*M_PI*pow(180/M_PI,2)*n_good4/NUM;
    
    return area;
}

/////////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] ){

	gsl_rng* rand_seed = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(rand_seed,time(NULL));

    if( argc == 1) {

	    FILE *fp = fopen("area.txt","w");

	    double B = B_MIN;

	    int i,j;
	    double bgal_max = 20.5, bgal_min = 15;
	    double beta_max = 20.5, beta_min = 15;
	    
	    double area;
	    double bgal = bgal_min;
	    while( bgal <= bgal_max ){
	    
		    printf("galactic-disk intersection angle b = %g\n",bgal);
		    fprintf(fp,"galactic-disk intersection angle b = %g\n",bgal);

	        double beta = beta_min;
		    while( beta <= beta_max ){
			    area = fun_area(beta,bgal,rand_seed);
			    printf("\tbeta = +/-%5.2f,  observable area is %g\n", beta, area);
			    fprintf(fp,"\tbeta = +/-%5.2f,  observable area is %g\n", beta, area);
			    beta += 0.5;
		    }

		    fprintf(fp,"\n");
		    bgal += 0.5;
	    }

	    fclose(fp);
    
    }
    else if( argc==3 ){
        double beta = atoi(argv[1]);
        double bgal = atoi(argv[2]);
        double s1,s2,s3,s4;
        double area = fun_area2(beta,bgal,rand_seed,&s1,&s2,&s3,&s4);
	    printf("> observable area is %g\n", area);
	    printf("> area1 = %g\n",s1);
	    printf("> area2 = %g\n",s2);
	    printf("> area3 = %g\n",s3);
	    printf("> area4 = %g\n",s4);
    }


	gsl_rng_free(rand_seed);
	return 0;
}
