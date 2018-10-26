#include "satellite.h"

void SatelliteOrbits_init( SatelliteOrbits* orbits, char *orbit_dir ){

	if( orbits == NULL ){
		printf("orbits is NULL, please allocate memory for it!\n");
		exit(0);
	}

	int i,j;

	orbits->orbit_num = 50;
	orbits->size    = malloc(sizeof(int)*orbits->orbit_num);
	orbits->t_start = malloc(sizeof(double)*orbits->orbit_num);
	orbits->t_end   = malloc(sizeof(double)*orbits->orbit_num);
	orbits->x_last  = malloc(sizeof(double)*orbits->orbit_num);
	orbits->y_last  = malloc(sizeof(double)*orbits->orbit_num);
	orbits->z_last  = malloc(sizeof(double)*orbits->orbit_num);

	for( i=0; i<orbits->orbit_num; i++ ){
		orbits->t[i] = NULL;
		orbits->x[i] = NULL;
		orbits->y[i] = NULL;
		orbits->z[i] = NULL;
		orbits->accel[i] = NULL;
	}

	char orbitfile[1024];

	int i_start = 0;

	for( i=i_start; i<50; i++ ){

		sprintf(orbitfile,"%s/%d.txt",orbit_dir,i+1);
		printf("debug --> loading orbit data from : %s\n",orbitfile);

		FILE *fp = NULL;
		fp = fopen(orbitfile,"r");

		if( fp == NULL ){
			printf("failed to open orbit data file: %s\n",orbitfile);
			exit(0);
		}

		char line[1024];

		//	count the number of lines
		int line_cnt = 0;
		while( fgets(line,1024,fp) != NULL ){
			if( line[0] != '#' )
				line_cnt++;
		}

		//	read data into *_tmp
		double *t_tmp = malloc(line_cnt*sizeof(double));
		double *x_tmp = malloc(line_cnt*sizeof(double));
		double *y_tmp = malloc(line_cnt*sizeof(double));
		double *z_tmp = malloc(line_cnt*sizeof(double));

		rewind(fp);
		int line_cnt2=0;
		while( fgets(line,1024,fp) != NULL ){
			if( line[0] == '#' )
				continue;

			char *p = strtok(line," \n");
			double tmp[7];
			int cnt=0;
			while (p != NULL) {
	            tmp[cnt] = atof(p);
	            cnt++;
	            p = strtok(NULL, " \n");
        	}

        	t_tmp[line_cnt2] = tmp[0];
        	x_tmp[line_cnt2] = tmp[1];
        	y_tmp[line_cnt2] = tmp[2];
        	z_tmp[line_cnt2] = tmp[3];
        	line_cnt2++;
		}

		if( line_cnt2 != line_cnt ){
			printf("&&&& line_cnt2 != line_cnt !!!!\n");
			exit(0);
		}

		int add_last_end = 0;
		if( i > i_start ){
			double t_diff = fabs(t_tmp[0]-orbits->t_end[i-1])*86400;
			// printf("@@ t_diff = %.10g\n",t_diff);
			if( t_diff < 250 ){
				add_last_end = 1;
			}
		}

		if( add_last_end == 0 ){

			orbits->size[i] = line_cnt;
			orbits->t[i] = malloc(line_cnt*sizeof(double));
			orbits->x[i] = malloc(line_cnt*sizeof(double));
			orbits->y[i] = malloc(line_cnt*sizeof(double));
			orbits->z[i] = malloc(line_cnt*sizeof(double));

			if( orbits->t[i] == NULL || 
				orbits->x[i] == NULL ||
				orbits->y[i] == NULL ||
				orbits->z[i] == NULL ){
				printf("Error in allocating memory for t or x or y or z!\n");
				exit(0);
			}

			for( j=0; j<line_cnt; j++ ){
				orbits->t[i][j] = t_tmp[j];
				orbits->x[i][j] = x_tmp[j];
				orbits->y[i][j] = y_tmp[j];
				orbits->z[i][j] = z_tmp[j];
			}

			orbits->t_start[i] = orbits->t[i][0];
			orbits->t_end[i]   = orbits->t[i][line_cnt-1];
			orbits->x_last[i]  = orbits->x[i][line_cnt-1];
			orbits->y_last[i]  = orbits->y[i][line_cnt-1];
			orbits->z_last[i]  = orbits->z[i][line_cnt-1];
		}
		else if( add_last_end == 1 ){

			orbits->size[i] = line_cnt+1;
			orbits->t[i] = malloc((line_cnt+1)*sizeof(double));
			orbits->x[i] = malloc((line_cnt+1)*sizeof(double));
			orbits->y[i] = malloc((line_cnt+1)*sizeof(double));
			orbits->z[i] = malloc((line_cnt+1)*sizeof(double));	

			if( orbits->t[i] == NULL || 
				orbits->x[i] == NULL ||
				orbits->y[i] == NULL ||
				orbits->z[i] == NULL ){
				printf("Error in allocating memory for t or x or y or z!\n");
				exit(0);
			}

			orbits->t[i][0] = orbits->t_end[i-1];
			orbits->x[i][0] = orbits->x_last[i-1];
			orbits->y[i][0] = orbits->y_last[i-1];
			orbits->z[i][0] = orbits->z_last[i-1];

			for( j=1; j<=line_cnt; j++ ){
				orbits->t[i][j] = t_tmp[j-1];
				orbits->x[i][j] = x_tmp[j-1];
				orbits->y[i][j] = y_tmp[j-1];
				orbits->z[i][j] = z_tmp[j-1];	
			}

			orbits->t_start[i] = orbits->t[i][0];
			orbits->t_end[i]   = orbits->t[i][line_cnt];
			orbits->x_last[i]  = orbits->x[i][line_cnt];
			orbits->y_last[i]  = orbits->y[i][line_cnt];
			orbits->z_last[i]  = orbits->z[i][line_cnt];
		}

		fclose(fp);
		free(t_tmp);
		free(x_tmp);
		free(y_tmp);
		free(z_tmp);

		orbits->accel[i] = gsl_interp_accel_alloc();
	}

	orbits->t_head = orbits->t_start[0];
	orbits->t_tail = orbits->t_end[orbits->orbit_num-1];

}

void SatelliteOrbits_free( SatelliteOrbits* orbits ){

	int i;
	for( i=0; i<orbits->orbit_num; i++ ){
		free(orbits->t[i]);
		free(orbits->x[i]);
		free(orbits->y[i]);
		free(orbits->z[i]);
		gsl_interp_accel_free(orbits->accel[i]);
	}

	free(orbits->size);
	free(orbits->t_start);
	free(orbits->t_end);
	free(orbits->x_last);
	free(orbits->y_last);
	free(orbits->z_last);
	
	orbits = NULL;
}

int get_satellite_position( 	double jdate, 
								double x[], 
								double *dist, 
								SatelliteOrbits* orbits ){
	int i,idx=-1;

	for( i=0; i<orbits->orbit_num; i++ ){
		if( jdate >= orbits->t_start[i] && jdate <= (orbits->t_end[i]-1e-15) ){
			idx = i;
			break;
		}
	}

	double *ptr_t=NULL, *ptr_x=NULL, *ptr_y=NULL, *ptr_z=NULL;

	if( idx != -1 ){
		ptr_t = orbits->t[idx];
		ptr_x = orbits->x[idx];
		ptr_y = orbits->y[idx];
		ptr_z = orbits->z[idx];
	}
	else if( idx == -1 ){
		x[0] = -999;
		x[1] = -999;
		x[2] = -999;
		*dist = -999;
		return 1;
	}

	int t_idx = gsl_interp_accel_find(	orbits->accel[idx],
										orbits->t[idx],
										orbits->size[idx],
										jdate );

#if defined(_USE_LAGRANGE_INTERP_)
	int idx_start = -1, idx_end = -1;

	if( t_idx < 3 ){
		idx_start = 0;
		idx_end   = idx_start + 8;
	}
	else if( t_idx > orbits->size[idx]-5 ){
		idx_end   = orbits->size[idx]-1;
		idx_start = idx_end - 8;
	}
	else {
		idx_start = t_idx - 3;
		idx_end   = idx_start + 8;
	}

	double tmp_x=0, tmp_y=0, tmp_z=0;

	for( i=idx_start; i<=idx_end; i++ ){
		double li = 1.0;
		int j;
		for( j=idx_start; j<=idx_end; j++ ){
			if( j != i ){
				li *= ( jdate-ptr_t[j] ) / ( ptr_t[i] - ptr_t[j] );
			}
		}
		tmp_x += ptr_x[i]*li;
		tmp_y += ptr_y[i]*li;
		tmp_z += ptr_z[i]*li;
	}

	x[0] = tmp_x;
	x[1] = tmp_y;
	x[2] = tmp_z;
	*dist = sqrt(tmp_x*tmp_x+tmp_y*tmp_y+tmp_z*tmp_z);
#else

	double t1,t2;
	double x1,y1,z1,x2,y2,z2;
	
	t1 = ptr_t[t_idx];
	t2 = ptr_t[t_idx+1];

	x1 = ptr_x[t_idx];
	y1 = ptr_y[t_idx];
	z1 = ptr_z[t_idx];
	x2 = ptr_x[t_idx+1];
	y2 = ptr_y[t_idx+1];
	z2 = ptr_z[t_idx+1];

	double l1 = sqrt(x1*x1 + y1*y1 + z1*z1);
	double l2 = sqrt(x2*x2 + y2*y2 + z2*z2);
	double theta = acos((x1*x2 + y1*y2 + z1*z2)/(l1*l2));
	double theta1 = (jdate-t1)/(t2-t1)*theta;
	double theta2 = theta-theta1;
	double l = (t2-jdate)/(t2-t1)*l1 + (jdate-t1)/(t2-t1)*l2;

	double sin_theta1 = sin(theta1);
	double sin_theta2 = sin(theta2);

	double x0 = sin_theta2*x1/l1 + sin_theta1*x2/l2;
	double y0 = sin_theta2*y1/l1 + sin_theta1*y2/l2;
	double z0 = sin_theta2*z1/l1 + sin_theta1*z2/l2;
	double l_ = sqrt(x0*x0 + y0*y0 + z0*z0);

	x[0] = x0*l/l_;
	x[1] = y0*l/l_;
	x[2] = z0*l/l_;
	*dist = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
#endif

	return 0;
}