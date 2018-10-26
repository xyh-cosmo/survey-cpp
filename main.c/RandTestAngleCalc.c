//	利用随机生成的指向来测试旋转角度的计算，最终的结果还可以进行统计
//	可以预期的是，当“新”、“旧”指向的经度相差不超过90度的时候，新版本和旧版本给出
//	的旋转角度值是一致的，只有当经度差值大于90度的时候旋转角度才会出现差异。

#include "SurveySim.h"

int main(int argc, char* argv[]){

	if( argc < 3 ){
		printf("usage: %s test_num result.txt\n",argv[0]);
		exit(0);
	}

	gsl_rng* r_seed = gsl_rng_alloc(gsl_rng_mt19937);

    srand(time(NULL));
    gsl_rng_set(r_seed,rand());
    char filename[1024];

	int num = atof(argv[1]);
	sprintf(filename,"%s",argv[2]);

	double ra_old,dec_old;
	double ra_new,dec_new;
	double angle_rot1,angle_rot2,angle_rot3;

	FILE *fp = fopen(filename,"w");

	int cnt=0;
	while( cnt < num ){

		// ra_old = gsl_ran_flat(r_seed, 0, 360);
  //       ra_new = gsl_ran_flat(r_seed, 0, 360);
  //       dec_old = gsl_ran_flat(r_seed, 0, 180);
  //       dec_new = gsl_ran_flat(r_seed, 0, 180);

		ra_old 	= gsl_ran_flat(r_seed, 0.0, 360);
		ra_new 	= gsl_ran_flat(r_seed, 0.0, 360);
	    dec_old = asin(gsl_ran_flat(r_seed, -1.0, 1.0))*180/M_PI;  // astronomical convention
	    dec_new = asin(gsl_ran_flat(r_seed, -1.0, 1.0))*180/M_PI;  // astronomical convention

		Get_RotationAngle_faster( ra_old,90-dec_old,
								   ra_new,90-dec_new,&angle_rot1);

		Get_RotationAngle_faster2( ra_old,90-dec_old,
								   ra_new,90-dec_new,&angle_rot2,0);

		Get_RotationAngle_Zhang( ra_old,90-dec_old,
								   ra_new,90-dec_new,&angle_rot3);
		
		cnt++;

		fprintf(fp,"%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n",
				ra_old,dec_old,ra_new,dec_new,angle_rot1,angle_rot2,angle_rot3);
	}

	fclose(fp);
	gsl_rng_free(r_seed);

	return 0;
}