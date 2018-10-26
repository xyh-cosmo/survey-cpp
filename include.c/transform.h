#ifndef _TRANSFORM_H_
#define _TRANSFORM_H_

//transformTools.c
void CoordinateSpin(double*, double*, double, int);
void CoordinateSpin_x(double coorO[3], double coorR[3], double angle );
void CoordinateSpinEquatorial2Ecliptic(double*, double*);
void Cartesian2Equatorial(double*, double*);
double getAngle132(double, double, double, double, double, double, double, double, double);
double calculateAngle(double , double , double , double );
void get_solar_normal_30( double px, double py, double pz, 
						  double nx, double ny, double nz, 
						  double* n30, double* n30_);

#endif