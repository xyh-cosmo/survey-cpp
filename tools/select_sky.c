#include <stdio.h>
#include <stdlib.h>
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
int main(){

	gsl_rng* rand_seed = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(rand_seed,time(NULL));

	int i,j;
	double u;
	double phi, theta;
	double x, y, z;
	double b;

	double beta = 25, BETA_MIN=20;
	double bgal = 25, BGAL_MIN=20;
	double cnt, n_good, n_bad;
	
	int N_TRY = 20;
	FILE *fp = fopen("beta_bgal_area.txt","w");

	while( beta >= BETA_MIN ){
		bgal = 25;
		while( bgal >= BGAL_MIN ){
			double area=0;
			int n_try=0;
			double Atmp[N_TRY];
			while( n_try < N_TRY ){
				cnt = 0;
				n_good = 0;
				n_bad = 0;
				while( cnt < NUM ){
					phi 	= gsl_ran_flat(rand_seed, 0.0, 2.0*M_PI);
					u 		= gsl_ran_flat(rand_seed, -1.0, 1.0);
					theta 	= acos(u) - M_PI_2;

					double t_coor[3] = {cos(theta)*cos(phi), cos(theta)*sin(phi), sin(theta)};  //	黄道坐标系里的xyz
					double r_coor[3], ll_coor[3];
					CoordinateSpin_x(t_coor, r_coor, -23.4522);
					Cartesian2Equatorial(r_coor, ll_coor);

					b = asin(-0.8676660 * cos(ll_coor[0] * PI_180) * cos(ll_coor[1] * PI_180)
							- 0.1980764 * sin(ll_coor[0] * PI_180) * cos(ll_coor[1] * PI_180)
							+ 0.4559840 * sin(ll_coor[1] * PI_180) ) / PI_180;

					theta *= (180*M_1_PI);
					int flag1 = (int)((fabs(theta) > beta) && (fabs(b) > bgal));
					int flag2 = (int)( (fabs(theta) < beta+0.25) && (fabs(b) < bgal+0.25) );
					int flag = (int)(flag1 || flag2);
					switch ( flag ) {
						case 1:
							n_good++;
							break;
						case 0:
							n_bad++;
							break;
					}

					cnt++;
				}

				double area_tmp = 4*M_PI*pow(180/M_PI,2)*n_good/NUM;
				// printf("A = %10.6f\t",area_tmp);
				area += area_tmp;
				Atmp[n_try] = area_tmp;
				n_try+=1;
			}

			area /= n_try;

			int i=0;
			double area_std=0;
			for(i=0; i<n_try; i++){
				area_std += (Atmp[i]-area)*(Atmp[i]-area);
			}

			area_std = pow(area_std,0.5);

			printf("%6.3f %6.3f area = %10.6f, std(area) = %10.6f\n",beta,bgal,area,area_std);
			fprintf(fp,"%6.3f %6.3f %8.4f %6.3f\n",beta,bgal,area,area_std);

			bgal -= 0.25;
		}

		beta -= 0.25;
	}

	fclose(fp);

	gsl_rng_free(rand_seed);
	return 0;
}
