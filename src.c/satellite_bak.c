#include "satellite.h"

void SatelliteOrbits_init( SatelliteOrbits* orbits, char orbit_dir[] ){

	if( orbits == NULL ){
		printf("orbits is NULL, please allocate memory for it!\n");
		exit(0);
	}

	orbits->orbit_num = 50;
	orbits->orbit_data = malloc(orbits->orbit_num*sizeof(OrbitData));

	int i,j;
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
			double t_diff = fabs(t_tmp[0]-orbits->orbit_data[i-1].t_end)*86400;
			// printf("@@ t_diff = %.10g\n",t_diff);
			if( t_diff < 250 ){
				add_last_end = 1;
			}
		}

// printf("debug -- add_last_end = %d\n",add_last_end);

		if( add_last_end == 0 ){

			orbits->orbit_data[i].size = line_cnt;
			orbits->orbit_data[i].t = malloc(line_cnt*sizeof(double));
			orbits->orbit_data[i].x = malloc(line_cnt*sizeof(double));
			orbits->orbit_data[i].y = malloc(line_cnt*sizeof(double));
			orbits->orbit_data[i].z = malloc(line_cnt*sizeof(double));

			for( j=0; j<line_cnt; j++ ){
				orbits->orbit_data[i].t[j] = t_tmp[j];
				orbits->orbit_data[i].x[j] = x_tmp[j];
				orbits->orbit_data[i].y[j] = y_tmp[j];
				orbits->orbit_data[i].z[j] = z_tmp[j];
			}

			orbits->orbit_data[i].t_start = orbits->orbit_data[i].t[0];
			orbits->orbit_data[i].t_end   = orbits->orbit_data[i].t[line_cnt-1];
			orbits->orbit_data[i].x_last  = orbits->orbit_data[i].x[line_cnt-1];
			orbits->orbit_data[i].y_last  = orbits->orbit_data[i].y[line_cnt-1];
			orbits->orbit_data[i].z_last  = orbits->orbit_data[i].z[line_cnt-1];
		}
		else if( add_last_end == 1 ){

			orbits->orbit_data[i].size = line_cnt+1;
			orbits->orbit_data[i].t = malloc((line_cnt+1)*sizeof(double));
			orbits->orbit_data[i].x = malloc((line_cnt+1)*sizeof(double));
			orbits->orbit_data[i].y = malloc((line_cnt+1)*sizeof(double));
			orbits->orbit_data[i].z = malloc((line_cnt+1)*sizeof(double));	

			orbits->orbit_data[i].t[0] = orbits->orbit_data[i-1].t_end;
			orbits->orbit_data[i].x[0] = orbits->orbit_data[i-1].x_last;
			orbits->orbit_data[i].y[0] = orbits->orbit_data[i-1].y_last;
			orbits->orbit_data[i].z[0] = orbits->orbit_data[i-1].z_last;

			for( j=1; j<=line_cnt; j++ ){
				orbits->orbit_data[i].t[j] = t_tmp[j-1];
				orbits->orbit_data[i].x[j] = x_tmp[j-1];
				orbits->orbit_data[i].y[j] = y_tmp[j-1];
				orbits->orbit_data[i].z[j] = z_tmp[j-1];	
			}

			orbits->orbit_data[i].t_start = orbits->orbit_data[i].t[0];
			orbits->orbit_data[i].t_end   = orbits->orbit_data[i].t[line_cnt];
			orbits->orbit_data[i].x_last  = orbits->orbit_data[i].x[line_cnt];
			orbits->orbit_data[i].y_last  = orbits->orbit_data[i].y[line_cnt];
			orbits->orbit_data[i].z_last  = orbits->orbit_data[i].z[line_cnt];
		}

		fclose(fp);
		free(t_tmp);
		free(x_tmp);
		free(y_tmp);
		free(z_tmp);

	}

	orbits->t_head = orbits->orbit_data[0].t[0];
	orbits->t_tail = orbits->orbit_data[orbits->orbit_num-1].t_end;
	orbits->accel = gsl_interp_accel_alloc();
}

void SatelliteOrbits_free( SatelliteOrbits* orbits ){

	int i;
	for( i=0; i<orbits->orbit_num; i++ ){
		free(orbits->orbit_data[i].t);
		free(orbits->orbit_data[i].x);
		free(orbits->orbit_data[i].y);
		free(orbits->orbit_data[i].z);
	}

	free(orbits->orbit_data);
	gsl_interp_accel_free(orbits->accel);
	
	orbits = NULL;
}

int get_satellite_position( 	double jdate, 
								double x[], 
								double *dist, 
								SatelliteOrbits* orbits ){
	int i,j,idx=-1;

	for( i=0; i<orbits->orbit_num; i++ ){
		if( jdate >= orbits->orbit_data[i].t_start && jdate <= (orbits->orbit_data[i].t_end-1e-15) ){
			idx = i;
			break;
		}
	}

	if( idx == -1 ){
		x[0] = -999;
		x[1] = -999;
		x[2] = -999;
		*dist = -999;
		return 1;
	}

	int t_idx = gsl_interp_accel_find(	orbits->accel,
										orbits->orbit_data[idx].t,
										orbits->orbit_data[idx].size,
										jdate );

	int idx_start = -1, idx_end = -1;

	if( t_idx < 3 ){
		idx_start = 0;
		idx_end   = idx_start + 8;
	}
	else if( t_idx > orbits->orbit_data[idx].size-5 ){
		idx_end   = orbits->orbit_data[idx].size-1;
		idx_start = idx_end - 8;
	}
	else {
		idx_start = t_idx - 3;
		idx_end   = idx_start + 8;
	}

	double tmp_x=0, tmp_y=0, tmp_z=0;
	for( i=idx_start; i<=idx_end; i++ ){
		double li = 1.0;
		for( j=idx_start; j<=idx_end; j++ ){
			if( j != i ){
				li *= ( jdate-orbits->orbit_data[idx].t[j] ) / ( orbits->orbit_data[idx].t[i] - orbits->orbit_data[idx].t[j] );
			}
		}
		tmp_x += orbits->orbit_data[idx].x[i]*li;
		tmp_y += orbits->orbit_data[idx].y[i]*li;
		tmp_z += orbits->orbit_data[idx].z[i]*li;
	}

	x[0] = tmp_x;
	x[1] = tmp_y;
	x[2] = tmp_z;

	*dist = sqrt(tmp_x*tmp_x+tmp_y*tmp_y+tmp_z*tmp_z);

	return 0;
}