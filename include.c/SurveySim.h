#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <stdbool.h>

#include "ephcom.h"
#include <mpi.h>

#include "macros.h"

//	some user defined functions (making use of GSL lib)
#include "GSL_funcs.h"

#include "cmg.h"
#include "sky.h"
#include "transform.h"
#include "CoverAreaStat.h"
#include "satellite.h"

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#ifndef _SURVEY_SIM_H_
#define _SURVEY_SIM_H_

#ifdef DEBUG
	#define _CRTDBG_MAP_ALLOC
	#define NEW   new( _NORMAL_BLOCK, __FILE__, __LINE__)
	#include <stdlib.h>
	#include <crtdbg.h>
#else
	#define NEW   new
#endif


// radius of the Earth
#define EARTH_RADIUS 6371.0 //unit, km

//	distance from the Earth to the Sun
#define DIST_SUN_TO_EARTH 149597870.7 // km

//	User defined maximum value, just for conveience, no particular meaning(s)
#define MAX_VALUE 100000000.0
#define VECTOR_VALUE  1.0

//	diameter of the deep survey (circular) areas
//#define DEEP_DIAM 7.95//7.7 //degree, the diameter of the area of deep survey
//#define DEEP_DIAM 8.0
#define DEEP_DIAM 8.9//7.7 //degree, the diameter of the area of deep survey

#define TARGET_COVER_NUM_LARGE 1
#define MAX_COVER_NUM_LARGE 1

#define TARGET_COVER_NUM_DEEP 4
#define MAX_COVER_NUM_DEEP 4

#define ULTRALDEEP_EX_NUM 8

#define Orbit_File_Num 50
#define Stop_Time_Num 5

#ifndef _STOP_TIME_TABLE_
#define _STOP_TIME_TABLE_
double stopTime[Stop_Time_Num][2];	//停靠维护的时间表
#endif


#define MAXLINE_STRING 1024

// 电池的最高、最低电量（这是张鑫设置的数值）
// #define BATTERY_MAX   5530.0
// #define BATTERY_LOW   -104.5

//	这是我重新设置的数值，依据文档里给出的一些设计参数计算得到
#define BATTERY_MAX   	97200.0
#define BATTERY_LOW		19440.0   // 19440.0=97200*0.2


// 归一化后的总体功率消耗（见张鑫所写的文档中的相关讨论）
// #define POWER_CONSUMPTION 1.0
#define POWER_CONSUMPTION 7.5   // kW

void init_default();
void init_ccd_pos(double ccd_pos_in_focus[18][2]);

#ifndef _CCD_X
	#define _CCD_X
	double CCD_X;
#endif
#ifndef _CCD_Y
	#define _CCD_Y
	double CCD_Y;
#endif
#ifndef _CCD_X_DEEP
	#define _CCD_X_DEEP
	double CCD_X_DEEP;
#endif
#ifndef _CCD_OVERLAP_X
	#define _CCD_OVERLAP_X
	double CCD_OVERLAP_X;
#endif
#ifndef _CCD_OVERLAP_Y
	#define _CCD_OVERLAP_Y
	double CCD_OVERLAP_Y;
#endif
#ifndef _CCD_OVERLAP_X_DEEP
	#define _CCD_OVERLAP_X_DEEP
	double CCD_OVERLAP_X_DEEP;
#endif
#ifndef _COS_SUN_POINT_ANGLE
	#define _COS_SUN_POINT_ANGLE
	double COS_SUN_POINT_ANGLE;
#endif
#ifndef _COS_MOON_POINT_ANGLE
	#define _COS_MOON_POINT_ANGLE
	double COS_MOON_POINT_ANGLE;
#endif
#ifndef _COS_SUN_PLANE_ANGLE
	#define _COS_SUN_PLANE_ANGLE
	double COS_SUN_PLANE_ANGLE;
#endif
#ifndef _COS_POINT_ZENITH_ANGLE_LIGHT_Max
	#define _COS_POINT_ZENITH_ANGLE_LIGHT_Max
	double COS_POINT_ZENITH_ANGLE_LIGHT_Max;
#endif
#ifndef _COS_POINT_ZENITH_ANGLE_LIGHT_Min
	#define _COS_POINT_ZENITH_ANGLE_LIGHT_Min
	double COS_POINT_ZENITH_ANGLE_LIGHT_Min;
#endif
#ifndef _COS_POINT_ZENITH_ANGLE_DARK
	#define _COS_POINT_ZENITH_ANGLE_DARK
	double COS_POINT_ZENITH_ANGLE_DARK;
#endif
#ifndef _EXTIME_DEEP
	#define _EXTIME_DEEP
	double EXTIME_DEEP;
#endif
#ifndef _EXTIME
	#define _EXTIME
	double EXTIME;
#endif
#ifndef _EXTIME_G_E
	#define _EXTIME_G_E
	double EXTIME_G_E;
#endif
#ifndef _EXTIME_SPEC
	#define _EXTIME_SPEC
	double EXTIME_SPEC;
#endif
#ifndef _RANGE_DEC_N
	#define _RANGE_DEC_N
	double RANGE_DEC_N;
#endif
#ifndef _RANGE_DEC_S
	#define _RANGE_DEC_S
	double RANGE_DEC_S;
#endif
#ifndef _ECLIPTIC_LAT_LIMIT_N
	#define _ECLIPTIC_LAT_LIMIT_N
	double ECLIPTIC_LAT_LIMIT_N;
#endif
#ifndef _ECLIPTIC_LAT_LIMIT_S
	#define _ECLIPTIC_LAT_LIMIT_S
	double ECLIPTIC_LAT_LIMIT_S;
#endif
//#ifndef _ECLIPTIC_LAT_LIMIT_N_SPEC
	//#define _ECLIPTIC_LAT_LIMIT_N_SPEC
	//double ECLIPTIC_LAT_LIMIT_N_SPEC;
//#endif
//#ifndef _ECLIPTIC_LAT_LIMIT_S_SPEC
	//#define _ECLIPTIC_LAT_LIMIT_S_SPEC
	//double ECLIPTIC_LAT_LIMIT_S_SPEC;
//#endif
#ifndef _Low_Galaxy_Img
	#define _Low_Galaxy_Img
	double Low_Galaxy_Img;
#endif
#ifndef _Ecliptic_Lat_Sec_Low
	#define _Ecliptic_Lat_Sec_Low
	double Ecliptic_Lat_Sec_Low;
#endif
#ifndef _Galaxy_B_Sec_Low
	#define _Galaxy_B_Sec_Low
	double Galaxy_B_Sec_Low;
	//double Low_Galaxy_Spec;
#endif
#ifndef _F_STATUS
	#define _F_STATUS
	int F_STATUS;
#endif
#ifndef _DEEP_AREA
	#define _DEEP_AREA
	double DEEP_AREA;
#endif
#ifndef _ULTRALDEPP_AREA
	#define _ULTRALDEPP_AREA
	double ULTRALDEPP_AREA;
#endif
//#ifndef _SPEC_AREA
	//#define _SPEC_AREA
	//double SPEC_AREA;
//#endif
#ifndef _START_TIME
	#define _START_TIME
	double startTime;
#endif
#ifndef _END_TIME
	#define _END_TIME
	double endTime;
#endif
#ifndef _JPL_FILE
	#define _JPL_FILE
	char jpl_fileName[100];
#endif
#ifndef _RESUALT_FILE
	#define _RESUALT_FILE
	char result_fileName[100];
#endif
#ifndef _STUTAS_FILE_IMG
	#define _STUTAS_FILE_IMG
	char status_img_fileName[100];
#endif
//#ifndef _STUTAS_FILE_SPEC
	//#define _STUTAS_FILE_SPEC
	//char status_spec_fileName[100];
//#endif

#ifndef _PANEL_TRANSE_ANGLE
     #define _PANEL_TRANSE_ANGLE
     double PANEL_TRANSE_ANGLE;
#endif

#ifndef _ORBIT_TIME_POINT
     #define _ORBIT_TIME_POINT
     double timePoint[Orbit_File_Num][2];
#endif

#ifndef JUMP_TIME
   #define JUMP_TIME
   double jump_time;
#endif



//	电池容量。每一步计算中都会更新电池容量的数值
#ifndef BATTERY_CAPACITY
   #define BATTERY_CAPACITY
   double battery_q;
#endif

#ifndef SHADOW_TIME_START
   #define SHADOW_TIME_START
   double shadow_st;
#endif

#ifndef SHADOW_TIME_END
   #define SHADOW_TIME_END
   double shadow_et;
#endif

#ifndef FILE_JPL
   #define FILE_JPL
   FILE *infp;
#endif


//	轨道数据
#ifndef DATA_ORBIT
   #define DATA_ORBIT
   double** orbitData;
#endif

#ifndef DATA_TIME
   #define DATA_TIME
   double* orbitTime[Orbit_File_Num];
#endif

#ifndef DATA_ORBIT_NUMS
   #define DATA_ORBIT_NUMS
   int orbitDataNum[Orbit_File_Num];
#endif

//==================================================================
//	The following are used to search satellite positions

#ifndef _GSL_INTERP_ACC_
#define _GSL_INTERP_ACC_
	gsl_interp_accel *locate_acc[Orbit_File_Num];
#endif

#ifndef _ORBIT_INDEX_
#define _ORBIT_INDEX_
	int orbit_index;
	int current_orbit_file_idx;
	int last_orbit_file_idx;
	int next_orbit_file_idx;
	int at_orbit_file_end;
#endif

int update_orbit_info( 	double 	time,
						int* 	ofile_idx,
						int*	last_ofile_idx,
						int*	next_ofile_idx,
						int*	at_end );

#ifndef _ORBIT_TIME_POINT
     #define _ORBIT_TIME_POINT
     double timePoint[Orbit_File_Num][2];
#endif

#ifndef DATA_ORBIT
   #define DATA_ORBIT
   double** orbitData;
#endif

#ifndef DATA_ORBIT_NUMS
   #define DATA_ORBIT_NUMS
   int orbitDataNum[Orbit_File_Num];
#endif

//	BETA, 太阳与望远镜轨道面的夹角
#ifndef _BETA_ANGLE_
	#define _BETA_ANGLE_
	double BETA_ANGLE;
#endif

//	在前几年内采用高纬度优先的策略，之后取消
#ifndef _HIGH_LATITUDE_PRIOR_TIME
	#define _HIGH_LATITUDE_PRIOR_TIME
	double HIGH_LATITUDE_PRIOR_TIME;
	int use_high_altitude_prior;
#endif

#ifndef _CONTINOUS_OBS_START_TIME
	#define _CONTINOUS_OBS_START_TIME
	double CONTINOUS_OBS_START_TIME;
#endif

#ifndef _DEC60_PRIOR_
#define _DEC60_PRIOR_
double DEC60_PRIOR_TIME;
#endif

//==================================================================

// #ifndef CMG_LIST
//    #define CMG_LIST
//    struct CMG_List* cmg_list;
// #endif


// LocateSatellite.c
void crossMultiple(double*, double*, double*);
double* readOrbitData(char*, int*);
void readAllorbitsFile(double**, int*);
void freeAllOrbitData(double**, int);
double hermite3(double ,double , double , double , double , double ,double);
//void matrixMultiple(double**, double**, double**);

// void locateSat1(double* , double, double**, int*, int);
double locateSat1(double* , double, double**, int*, int); // 直接返回卫星到地心的距离，简化判断卫星的位置是否获取成功，如果小于0则表示失败

void locateSat(double*, double, double, double, double, double, double);
void locateSun(FILE*, double, double*);
void locateMoon(FILE*, double, double*);
double getSun2OrbitAngle1(FILE*, double,double**, int*, int);
void calculateSatEarthPoint(double , double , const double , double* , double*, FILE *);

// //transformTools.c
// void CoordinateSpin(double*, double*, double, int);
// void CoordinateSpin_x(double coorO[3], double coorR[3], double angle );
// void CoordinateSpinEquatorial2Ecliptic(double*, double*);
// void Cartesian2Equatorial(double*, double*);
// double getAngle132(double, double, double, double, double, double, double, double, double);
// double calculateAngle(double , double , double , double );
// void get_solar_normal_30( double px, double py, double pz, 
// 						  double nx, double ny, double nz, 
// 						  double* n30, double* n30_);

//SkyAreaSplit.c


//SurveyConditionLimit.c
double isInDeepSurveyArea(double, double, double);
int IsInSunSide(double*, double*);
int IsInSunSideByTime(double, double , double );
// int IsObscureByEarth_1(double*, double*, double, double, double, double*);
int IsObscureByEarth_1_old( double sat[3], 
							double sun[3], 
							double p_x, 
							double p_y, 
							double p_z, 
							double* angleValue );
int IsObscureByEarth_1( double sat[3], 
                        double sun[3], 
                        double p_x, 
                        double p_y, 
                        double p_z, 
                        double* angleValue );

bool debug_IsObscureByEarth( double sat[3], 
							 double sun[3], 
							 double p_x, 
							 double p_y, 
							 double p_z );

double IsObscureBySun(double*, double, double, double);
double IsObscureByMoon(double*, double, double, double);
int IsInSAA(double*);
// double calculateTransTime(double); // moved into cmg.h
double getExposureTime_Zodical(double, double);
double getPlaneEnergy(double, double);
void aquireShadowTime(double*, double*, double);

//programStatusInfo.c
void parseStatusFile(char* , double* , int* ,double* ,double* ,double* ,double* ,double* ,double* ,double* ,double* ,int* , SKY_Coord* );
void parseConfigFile(char*);
void printConfig();


void setDoubleArrayZero(double* array, int num);

/*
 * @curTelPoint:望远镜指向，curTelPoint[0]--RA, curTelPoint[1]--DEC
 */
void FindNewPointByAllSearch(	double* 	currentTime,
								int* 		id,
								SKY_Coord* 	skyMap,
								int 		all_sky_num,
								FILE*		infp,
								double 		greenStartTime,
								double		initalGreenWithRa,
								double*		outputresult,
								double*		outputresult_fail,
								double**	orbitData,
								int*		orbitDataNum,
								int			p_rank,
								int 		p_size);

//  稍微新一些的版本，主要修改了权重和转动、能源等的计算
void FindNewPointByAllSearch2(  double*     currentTime,
                                int*        id,
                                SKY_Coord*  skyMap,
                                int         all_sky_num,
                                FILE*       infp,
                                double      greenStartTime,
                                double      initalGreenWithRa,
                                double*     outputresult,
                                double*     outputresult_fail,
                                double**    orbitData,
                                int*        orbitDataNum,
                                int         p_rank,
                                int         p_size);

//	新版本的搜索函数增加了两个输入变量:cur_ra_tel,cur_dec_tel，用于记录当前望远镜的指向。
//	旧版本中，望远镜当前的某一时刻的指向总是指向上一次观测天区的指向。
void FindNewPointByAllSearch3(	double* 	currentTime,
                              	int* 		id,
                              	double*		cur_ra_tel,
                              	double*		cur_dec_tel,
                              	SKY_Coord* 	skyMap,
                              	int 		all_sky_num,
                              	FILE*		infp,
                              	double 		greenStartTime,
                              	double		initalGreenWithRa,
                              	double*		outputresult,
                              	double*		outputresult_fail,
                              	double**	orbitData,
                              	int*		orbitDataNum,
                              	int			p_rank,
                              	int         p_size,
                                FILE*       fp_drift );

void FindNewPointByAllSearch4(	double* 	currentTime,
                                double*     lastObsEndTime,
                                int*        trysteps,
                              	int* 		id,
                              	// double*		cur_ra_tel,
                              	// double*		cur_dec_tel,
                              	SKY_Coord* 	skyMap,
                              	int 		all_sky_num,
                              	FILE*		infp,
                              	double 		greenStartTime,
                              	double		initalGreenWithRa,
                              	double*		outputresult,
                              	double*		outputresult_fail,
                              	double**	orbitData,
                              	int*		orbitDataNum,
                              	int			p_rank,
                              	int         p_size );

void TestCondition( double* 	currentTime,
					int* 		id,
					SKY_Coord* 	skyMap,
					int 		all_sky_num,
					FILE*		infp,
					double 		greenStartTime,
					double		initalGreenWithRa,
					double*		outputresult,
					double*		outputresult_fail,
					double**	orbitData,
					int*		orbitDataNum,
					int			p_rank,
					int         p_size );

#endif // _SURVEY_SIM_H_
