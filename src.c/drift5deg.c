//	分析模拟结果时发现，会出现很多长达2000秒甚至3500秒的间隔期，在此期间望远镜没有进行
//	任何观测. 为了充分利用这些间隔期,我们让望远镜在搜索可观测天区失败的时候把指向转动5度.
//	选择转动5度是因为这样的转动对CMG的影响可以忽略不计. 
//
//	具体转动的方式有两种:
//	1) 望远镜的指向与天顶位置不同,那么我们将望远镜朝着天顶方向转动5度;
//	2) 望远镜的指向与天顶十分接近,那么我们将望远镜朝着它运动的方向转动5度;
//
//	这里需要详细说明第二种情况. 具体转动的指向只是估算的结果,我们在望远镜所在(天球上的)
//	位置附近,将其速度方向沿着经纬圈分解. 如果将速度向量理解为一个普通的向量,那么其末端
//	就定义了一个新的位置,对应着一个新的天球坐标. 接下来要做的就是给望远镜的指向增加一个
//	"增量". 这个增量的计算是将望远镜的速度放在(ra=0,dec=0)的位置,然后以末端对应的
//	ra,dec作为"增量".
//
//	注意:第二种情况里的做法在数学上很不严格,但我们的主要目的是让望远镜尽可能地朝着天顶
//	的方向运动,而天顶又严格跟随望远镜运动,在天球上扫过,因此总的来讲可以满足我们设定的目标.
//	另外,出现第二种情况的概率其实非常低,所以即使这里的建模不合理也不会有明显的影响.
//
//	2018-09-18: 发现原先设想的5度角度转动矩阵并不能保证“不产生像旋”，其效果其实久等价于
//				张新最初采取的转动方式。现在通过暴力搜索的办法，找到一个合适的新指向，是得
//				“不产生像旋的”转动角度小于5度（目前看来，大部分这类转动的角度都是略小于5度

#include "SurveySim.h"

void get_tel_velocity( double sat_old[], double sat_new[], double sat_vel[]){
	double dp = sqrt( (sat_new[0] - sat_old[0])*(sat_new[0] - sat_old[0]) 
					+ (sat_new[1] - sat_old[1])*(sat_new[1] - sat_old[1]) 
					+ (sat_new[2] - sat_old[2])*(sat_new[2] - sat_old[2]));

	sat_vel[0] = (sat_new[0] - sat_old[0])/dp;
	sat_vel[1] = (sat_new[1] - sat_old[1])/dp;
	sat_vel[2] = (sat_new[2] - sat_old[2])/dp;
}

void Drift_by_5deg( double *cur_ra_tel, 
					double *cur_dec_tel, 
					double sat[],
					double velocity_tel[],
					double *drift_angle,
					int p_rank ){

	// if( p_rank == 0)
	// 	printf("@ drifting ...\n");

//  第一步：根据望远镜指向和天顶指向或者望远镜飞行的速度，得到旋转矩阵

	double ra_old = *cur_ra_tel;
	double dec_old= *cur_dec_tel;
	double ra_new, dec_new;
	double radec_tmp[2];

	//  根据当前望远镜的位置获取天顶位置的黄道坐标
	Cartesian2Equatorial(sat,radec_tmp);
	double cur_ra_zenith = radec_tmp[0];
	double cur_dec_zenith= radec_tmp[1];

	if( fabs(ra_old-cur_ra_zenith) > 1e-2 || fabs(dec_old-cur_dec_zenith) > 1e-2 ){
	//  望远镜指向与天顶不重合
		ra_new = cur_ra_zenith;
		dec_new= cur_dec_zenith;
	}
	else{
	//  此时望远镜指向与天定方向可以认为是重合的，所以改为往望远镜运动的方向转动5度。
	//  具体转动的方向直接按照望远镜当时的速度方向来估算，做法是将估算速度（实际上时位移矢量）的
	//  起始点移动到坐标原点，得到末点对应的黄道坐标，分别加到ra_old,dec_old 上去作为要转向的指向。
		Cartesian2Equatorial(velocity_tel,radec_tmp);
		ra_new  = ra_old+radec_tmp[0];
		dec_new = dec_old+radec_tmp[1];

		if( ra_new > 360 ) ra_new -= 360;
		if( ra_new < 0 ) ra_new += 360;
		if( dec_new > 90 ) dec_new = 90 - (dec_new-90);
		if( dec_new < -90 ) dec_new = -90 + (-dec_new-90);
	}

	int status = 0;
	double rot_axis[3];
	
	// double R5deg[9] = {0};
	// status = GenRotationMatrix5deg( ra_old, dec_old, ra_new, dec_new, R5deg, rot_axis, p_rank );
//	直接确定出“最近”的旋转轴
	double p_old[3] = { cos(dec_old*M_PI/180)*cos(ra_old*M_PI/180),
						cos(dec_old*M_PI/180)*sin(ra_old*M_PI/180),
						sin(dec_old*M_PI/180) };
	double p_new[3] = { cos(dec_new*M_PI/180)*cos(ra_new*M_PI/180),
						cos(dec_new*M_PI/180)*sin(ra_new*M_PI/180),
						sin(dec_new*M_PI/180) };
	
	if( VecDiffNorm(p_old,p_new) < 1e-15 ){
		status = 1;
		*drift_angle = 0.0;
		return;
	}

	CrossProduct(p_old,p_new,rot_axis);
	double rot_axis_norm = VecNorm(rot_axis);
	rot_axis[0] /= rot_axis_norm;
	rot_axis[1] /= rot_axis_norm;
	rot_axis[2] /= rot_axis_norm;

//  第二步：根据旋转轴，绕着这个轴转动5度。 （其实有可能只需转动很小的角度就可以将望远镜转动到当时的天顶方向，暂时忽略这种可能）
	double ra_tmp, dec_tmp;

	if( status == 0 ){

		double rotation_matrix[9] = {0};
		double theta=0.25;
		double tAngle=0;
		int search_status=1;

		while( theta <= 10.0 ){
			GenRotationMatrix( rot_axis, theta, rotation_matrix, p_rank );
			Rotate_about_AxisP_by_5deg( rotation_matrix, *cur_ra_tel, *cur_dec_tel, &ra_tmp, &dec_tmp, p_rank );
			
			Get_RotationAngle_faster2( *cur_ra_tel, 90 - *cur_dec_tel, 
										ra_tmp, 90-dec_tmp, &tAngle, p_rank);

			if( fabs(tAngle-4.9995) <= 0.05 ){
				*drift_angle = tAngle;
				search_status = 0;
				break;
			}

			theta += 0.025;
		}

		if( search_status == 0 ){
			*cur_ra_tel = ra_tmp;
			*cur_dec_tel= dec_tmp;
		} else {
			// 保持原状
			*drift_angle = 0;
		}
	}

}

void Write_drift_state( FILE* outf,
						int p_rank,
						double jdate, double dec, double ra,
						double satx, double saty, double satz,
						double sunx, double suny, double sunz,
						double moonx, double moony, double moonz,
						double area1, double area2, double deepflag,
						int indisk,
						double exptime,
						double tAngle,
						int insunside,
						double cmgtotal,
						double battery,
						double panel_sun_angle,
						double saa_time,
						int skyid ){
	if( p_rank == 0 ){
	//  time, lat, lon, (1,2,3)
	    fprintf(outf, "%15.8f %8.4f %8.4f ", jdate, dec, ra);
	//  satx, saty, satz, (4,5,6)
	    fprintf(outf, "%10.4f %10.4f %10.4f ", satx, saty, satz);
	//  sunx, suny, sunz, (7,8,9)
	    fprintf(outf, "%10.4f %10.4f %10.4f ", sunx, suny, sunz);
	//  moonx, moony, moonz, (10,11,12)
	    fprintf(outf, "%10.4f %10.4f %10.4f ", moonx, moony, moonz);
	//  survey_area, survey_area2, isdeep, 大面积巡天，极深度巡天 (13,14,15)
	    fprintf(outf, "%12.6f %12.6f %10.4f ", area1,area2,deepflag);
	//  isInGalaxydisk, (16)
	    fprintf(outf, "%2d ", indisk);
	//  exposuretime, (17)
	    fprintf(outf, "%10.4f ", exptime);
	//  transAngle, 转动的角度 (18)
	    fprintf(outf, "%10.5f ", tAngle);
	//  isInSunsid (19)
	    fprintf(outf, "%2d ", insunside);
	//  CMG ... (20)
	    fprintf(outf, "%10.4f ", cmgtotal);
	//  battery level, (21)
	    fprintf(outf, "%10.4f ", battery);
	//  太阳和帆板法线的夹角 (22)
	    fprintf(outf, "%10.4f ", panel_sun_angle);
	//  SAA 时间 (23)
	    fprintf(outf, "%10.4f ", saa_time);
	//  天区的ID (24)
	    fprintf(outf, "%8d ", skyid);
	//  switch to next line
	    fprintf(outf, "\n");
	}
}
