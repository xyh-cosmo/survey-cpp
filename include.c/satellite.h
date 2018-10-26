#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "GSL_funcs.h"

#ifndef _SATELLITE_H_
#define _SATELLITE_H_


struct SatelliteOrbits_{
	double t_head, t_tail;
	double *t[50];
	double *x[50];
	double *y[50];
	double *z[50];
	double *t_start, *t_end;
	double *x_last, *y_last, *z_last;
	int *size;
	int orbit_num;
	gsl_interp_accel *accel[50];
};

typedef struct SatelliteOrbits_ SatelliteOrbits;

void SatelliteOrbits_init( SatelliteOrbits* orbits, char *orbit_dir );
void SatelliteOrbits_free( SatelliteOrbits* orbits );

int get_satellite_position( double jdate, 
							double x[], 
							double *dist, 
							SatelliteOrbits* orbits );

#endif //_SATELLITE_H_