//	利用随机生成的指向来测试旋转角度的计算，最终的结果还可以进行统计
//	可以预期的是，当“新”、“旧”指向的经度相差不超过90度的时候，新版本和旧版本给出
//	的旋转角度值是一致的，只有当经度差值大于90度的时候旋转角度才会出现差异。

#include "SurveySim.h"


int main(int argc, char* argv[]){

	if( argc < 5 ){
		printf("Purpose: calculation rotation angle for given 'old' and 'new' directions\n");
		printf("Note the difference between astronomical and physical/mathematical defintions of 'ra' and 'dec' !\n");
		printf("usage:\n\t%s ra_old dec_old ra_new dec_new\n",argv[0]);
		exit(0);
	}

	double ra_old  = atof(argv[1]);
	double dec_old = atof(argv[2]);
	double ra_new  = atof(argv[3]);
	double dec_new = atof(argv[4]);
	double angle_zx,angle_xyh;

	Get_RotationAngle_Zhang( ra_old, 90-dec_old, ra_new, 90-dec_new, &angle_zx);
	Get_RotationAngle_faster2( ra_old, 90-dec_old, ra_new, 90-dec_new, &angle_xyh,0);

	double cmg_use_zx = get_cmg_use(angle_zx);
	double cmg_use_xyh= get_cmg_use(angle_xyh);
	double tTime_zx   = calculateTransTime( angle_zx );
	double tTime_xyh  = calculateTransTime( angle_xyh );
	
	printf("Angle_zx  : %8.5f, cmg_use = %5.3f, rot_time = %8.5f\n", angle_zx, cmg_use_zx, tTime_zx);
	printf("Angle_xyh : %8.5f, cmg_use = %5.3f, rot_time = %8.5f\n", angle_xyh, cmg_use_xyh, tTime_xyh);

	return 0;
}
