//  This header file is added by XYH @ 2018-05-09

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

//	GSL header files
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_linalg.h>


#ifndef _GSL_FUNCS_H_
#define _GSL_FUNCS_H_

#define _newline_ {printf("\n");}


//  some mathematical constants
#define E 2.718281828459045
#define PI  M_PI  // M_PI is taken from GSL
#define PI_180 (M_PI/180)

#define radian2degree (180*M_1_PI)
#define degree2radian (M_PI/180)

#define COS_25deg 0.9063077870366499	//	cos(25deg)
#define COS_50deg 0.6427876096865394  	// cos(50deg)

void ErrorInfo( char* info );
void print_vec(char* vecname, double* vec);
void print_vec_diff(char* vecname1, double* vec1, char* vecname2, double* vec2);
void print_mat(char* matname, double* mat);
void print_rotation_axis(char* matname, double* mat);
void print_rotation_angle(char* matname, double* mat);
double DotProduct( double *u, double *v );
double VecNorm( double *v );
double VecDiffNorm( double *u, double *v );
double GetMatrixTrace( double *M );
double GetRotationAngleFromMatrix(double *M);

int CrossProduct( double *u, double *v, double *w );
int MatrixVecProduct( double *M, double *V, double *W );
int MatrixMultiplication( double* R, double* S, double* T );
int GenRotationMatrix( double* u, double angle_deg, double* R, int rank );



int Get_RotationAngle2( double ra_old, double dec_old,
                        double ra_new, double dec_new,
                        double* angle_rot );

int Get_RotationAngle_faster( double ra_old, double dec_old,
                              double ra_new, double dec_new,
                              double* angle_rot );

//  新版本的获取转动角度的函数
int Get_RotationAngle_faster2( double ra_old, double dec_old,
                               double ra_new, double dec_new,
                               double* angle_rot,
                               int p_rank );

int Get_RotationAngle_Axis( double ra_old, double dec_old,
                            double ra_new, double dec_new,
                            double* angle_rot, double* axis );

//  张鑫采用的转动角计算方式
int Get_RotationAngle_Zhang( double ra_old, double dec_old, double ra_new, double dec_new, double *angle );


//  按照给定的新旧指向，生成一个将旧指向往新指向转动5度的旋转矩阵
int GenRotationMatrix5deg(  double ra_old, double dec_old, 
                            double ra_new, double dec_new,
                            double R5deg[],
                            double rot_axis[],
                            int p_rank );


//  将望远镜的指向按照旋转矩阵R5deg[] 旋转5度
int Rotate_about_AxisP_by_5deg( double R5deg[], 
                                double ra_old, 
                                double dec_old, 
                                double *ra_new, 
                                double *dec_new,
                                int p_rank );

//  估算（归一化的）望远镜运动速度
void get_tel_velocity( double sat_old[], double sat_new[], double sat_vel[]);

void Drift_by_5deg( double *cur_ra_tel, 
                    double *cur_dec_tel, 
                    double sat[],
                    double velocity_tel[],
                    double *drift_angle,
                    int p_rank );

void Write_drift_state( FILE* fp,
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
                        int skyid );

int init_gsl_interp_accel();
int free_gsl_interp_accel();

int Test_OrbitFile( double t, int idx );
int Find_OrbitFile_Index( double t, int* idx );


//  检验帆板转动后能否满足“25度夹角”的要求
int TestPanelAngle(	double *sun, 
					double dist_sun,
					double x,
					double y,
					double z,
					double nx,
					double ny,
					double nz,
					double *cosval,
          double rotate_angle_max,
          int p_rank );

double locateSat_GSL( double    t,
                      double    *result );
#endif
