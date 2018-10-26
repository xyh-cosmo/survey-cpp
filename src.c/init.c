#include "SurveySim.h"


void init_default(){

//	==========================================
//	停靠维护的时间表
	stopTime[0][0] = 2460131.56;
	stopTime[0][1] = 2460205.53;

	stopTime[1][0] = 2460887.37;
	stopTime[1][1] = 2460997.22;

	stopTime[2][0] = 2461688.66;
	stopTime[2][1] = 2461796.99;

	stopTime[3][0] = 2462429.30;
	stopTime[3][1] = 2462538.41;

	stopTime[4][0] = 2463166.49;
	stopTime[4][1] = 2463275.57;

//	==========================================

    CCD_X = 0.5195;   //degree
    CCD_Y = 0.1024;   //degree
    CCD_X_DEEP = 0.1509; // degree

    CCD_OVERLAP_X = 0.0027778;  //0.0027778 //degree, 10 arcsec
    CCD_OVERLAP_Y = 0.0027778;
    CCD_OVERLAP_X_DEEP = 0.0027778;

    COS_SUN_POINT_ANGLE = 0.64279;
    COS_MOON_POINT_ANGLE = 0.76604;
    COS_SUN_PLANE_ANGLE = 0.90631;
    COS_POINT_ZENITH_ANGLE_LIGHT_Max = 0.866025;  // earth light limb 80 deg
    COS_POINT_ZENITH_ANGLE_LIGHT_Min = 0.866025;  // earth light limb 80 deg
    COS_POINT_ZENITH_ANGLE_DARK = 0.173648;

    EXTIME_DEEP = 200.0;
    EXTIME = 150.0;
    EXTIME_G_E = 100.0;
    EXTIME_SPEC = 200.0;

    RANGE_DEC_N = 90;
    RANGE_DEC_S = -90;

    ECLIPTIC_LAT_LIMIT_N = 21.5;
    ECLIPTIC_LAT_LIMIT_S = -21.5;

    Low_Galaxy_Img = 20;
    Ecliptic_Lat_Sec_Low = 15;

    Galaxy_B_Sec_Low = 15;

    F_STATUS = 20000;

    DEEP_AREA = 0.2;
    ULTRALDEPP_AREA = 0.1;
    //SPEC_AREA = 5;

    startTime = 0;
    endTime = 0;
    PANEL_TRANSE_ANGLE=0;
    strcpy(jpl_fileName,"jpl.405");
    strcpy(result_fileName,"SS_TEST.dat");
    strcpy(status_img_fileName,"status_img.dat");
    //strcpy(status_spec_fileName,"status_spec.dat");

    sky_id_tracker = NULL;
    sky_num_remained = 0;

    HIGH_LATITUDE_PRIOR_TIME = 5; // 默认前5年内采用高纬度优先的策略
    use_high_altitude_prior=0;
    
    CONTINOUS_OBS_START_TIME = 0; // 默认从最开始的时候就使用连续观测策略

    jump_time = 60.0;   // in unit of seconds
    DEC60_PRIOR_TIME = 0;

#if defined(_SLEW_TIME_GSL_INTERP_)
    init_TransTime();
#endif
}
