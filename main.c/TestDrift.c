//	利用随机生成的指向来测试旋转角度的计算，最终的结果还可以进行统计
//	可以预期的是，当“新”、“旧”指向的经度相差不超过90度的时候，新版本和旧版本给出
//	的旋转角度值是一致的，只有当经度差值大于90度的时候旋转角度才会出现差异。

#include "SurveySim.h"

int main(int argc, char* argv[]){

    const gsl_rng_type *T;
	gsl_rng *r;

    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

	SatelliteOrbits *sat_orbits = malloc(sizeof(SatelliteOrbits));

    SatelliteOrbits_init(sat_orbits,"orbit20160925");

    double t_start = 2459766;
    double dist;

    double sat_old[3],sat_new[3];
    double t = t_start;
	int cnt=0, cnt_good=0;

    while( (t < t_start + 365.25*10) ){
        int status1 = get_satellite_position(t,sat_old,&dist,sat_orbits);
		int status2 = get_satellite_position(t+90./86400,sat_new,&dist,sat_orbits);

		if( status1 == 0 && status2 == 0 ){
			double sat_v[3];
			get_tel_velocity(sat_old,sat_new,sat_v);
			double ra_old = gsl_ran_flat(r,0,360);
			double dec_old = gsl_ran_flat(r,-89.99,89.99);
			double ra_tmp = ra_old, dec_tmp = dec_old;

			double drift_angle;
			Drift_by_5deg( &ra_tmp, &dec_tmp, sat_old, sat_v, &drift_angle, 0);

			double ra_new  = ra_tmp;
			double dec_new = dec_tmp;
			double angle_zx,angle_xyh;

			// Get_RotationAngle_Zhang( ra_old, 90-dec_old, ra_new, 90-dec_new, &angle_zx);
			Get_RotationAngle_faster2( ra_old, 90-dec_old, ra_new, 90-dec_new, &angle_xyh,0);

			// double cmg_use_zx = get_cmg_use(angle_zx);
            // double tTime_zx   = calculateTransTime( angle_zx );
			double cmg_use_xyh= get_cmg_use(angle_xyh);
			double tTime_xyh  = calculateTransTime( angle_xyh );
			
			// printf("Angle_zx  : %8.5f, cmg_use = %5.3f, rot_time = %8.5f\n", angle_zx, cmg_use_zx, tTime_zx);
			// printf("Yr = %6.3g, Angle_xyh : %8.5f, cmg_use = %5.3f, rot_time = %8.5f\n", 
			// 			(t-2459766)/365.25, angle_xyh, cmg_use_xyh, tTime_xyh);

			if( angle_xyh <= 5.0 ){
				cnt_good++;
			}

			cnt++;
		}

        t += 100./86400;
    }

	printf("cnt_total = %10d\n",cnt);
	printf("cnt_good  = %10d\n",cnt_good);


    SatelliteOrbits_free(sat_orbits);

    gsl_rng_free(r);

	return 0;
}
