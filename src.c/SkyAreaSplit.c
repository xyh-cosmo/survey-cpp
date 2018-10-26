#include "SurveySim.h"


int* readData(char* fileName, int elemNum, int* datNum) {
    FILE *skyData;
    skyData = fopen(fileName, "r");

    if (skyData == NULL) {
        printf("No such file1 !!!\n");
        exit(1);
    }
    int MAX_ROW = 1500000;
    int* tDat = (int *)malloc(MAX_ROW*elemNum*sizeof(int));
    int num = 0;
    char line[MAXLINE_STRING];
    char cmp_char = '#';
    while (fgets(line, MAXLINE_STRING, skyData) != NULL) {
        if(line[0]==cmp_char) {
            continue;
        }
        char *p;
        p = strtok(line, " \t\n");
        int cNum = 0;
        while (p != NULL && cNum < elemNum) {
            double cc = atof(p);
            *(tDat + num * elemNum + cNum) = (int)cc;
            cNum++;
            p = strtok(NULL, " \t\n");
        }
        num++;
    }
    *datNum = num;
    int i;
    int* dat = (int *)malloc(num*elemNum*sizeof(int));
    for (i = 0; i < num; i ++) {
        dat[i] = tDat[i];
    }
    free(tDat);
    fclose(skyData);
    return dat;
}

void init_SKY_Coord( SKY_Coord* sky ){
    sky->flag               = 0;
    sky->id                 = -1;
    sky->IsInSunSideFlag    = -1;
    sky->isObserve          = -1;
    sky->targetCoverNum     = TARGET_COVER_NUM_LARGE;
    sky->maxCoverNum        = MAX_COVER_NUM_LARGE;

    sky->left_neighbour     = -1;
    sky->right_neighbour    = -1;
    sky->down_neighbour     = -1;
    sky->up_neighbour       = -1;

    int i;
    for( i=0; i<10; i++)
        sky->filter[i] = 0;

    sky->weight = MAX_VALUE;
    sky->wFactor = 1.0;

    sky->ra     = MAX_VALUE;
    sky->dec    = MAX_VALUE;
    sky->inDeepFlag = -1;
    sky->elat = 0;
    sky->gb = 0;
    
    sky->x = 0;
    sky->y = 0;
    sky->z = 0;
    sky->dist2orig = 0;
    
    sky->nx = 0;
    sky->ny = 0;
    sky->nz = 0;
    
    sky->nx30 = 0;
    sky->ny30 = 0;
    sky->nz30 = 0;
    
    sky->nx30_ = 0;
    sky->ny30_ = 0;
    sky->nz30_ = 0;
    
    sky->nnx = 0;
    sky->nny = 0;
    sky->nnz = 0;

    sky->interval_lon   = 0.0;
    sky->interval_lat   = 0.0;
}

void copySkyMapData(SKY_Coord* from, SKY_Coord* to) {
    to->weight = from->weight;
    to->wFactor = from->wFactor;
    to->flag = from->flag;
    to->dec = from->dec;
    to->ra = from->ra;
    to->id = from->id;
    to->inDeepFlag = from->inDeepFlag;
    to->IsInSunSideFlag = from->IsInSunSideFlag;
    to->elat = from->elat; // ecliptic latitude
    to->gb = from->gb; // galaxy latitude

    //Cartesian coordinate ,related with equator coordinate
    to->x = from->x;
    to->y = from->y;
    to->z = from->z;
    to->dist2orig = from->dist2orig;

    to->nx = from->nx;
    to->ny = from->ny;
    to->nz = from->nz;
    to->nx30 = from->nx30;
    to->ny30 = from->ny30;
    to->nz30 = from->nz30;
    to->nx30_ = from->nx30_;
    to->ny30_ = from->ny30_;
    to->nz30_ = from->nz30_;
    to->nnx = from->nnx;
    to->nny = from->nny;
    to->nnz = from->nnz;

    to->interval_lon = from->interval_lon;
    to->interval_lat = from->interval_lat;

    to->isObserve = from->isObserve;

    to->targetCoverNum = from->targetCoverNum;
    to->maxCoverNum = from->maxCoverNum;

    to->up_neighbour = from->up_neighbour;
    to->down_neighbour = from->down_neighbour;
    to->left_neighbour = from->left_neighbour;
    to->right_neighbour = from->right_neighbour;
}

/**
 * size[0] = x_length
 * size[1] = y_length
 */
void calculateSkySize(  double size[2],
                        double x_size, 
                        double y_size,
                        double x_overLap, 
                        double y_overlap ) {
    int xLen, yLen;
    double xDelt, yDelt;
    xDelt = x_size - x_overLap;
    yDelt = y_size - y_overlap;

    yLen = (int)(ceil(180/yDelt) + 1);
    xLen = (int)(ceil(360/xDelt) + 1);

    size[0] = xLen;
    size[1] = yLen;
}

/*
* 划分天区
*/
int SplitSkyArea( int yLen, 
                  int xLen, 
                  SKY_Coord* skyMap,
                  double x_size, 
                  double y_size ) {
    int i,j;
    int cnt=0;  //for debug
    int skyAreaNum = 0;
    int deepSkyNum = 0;
    
    double xDelt = x_size - CCD_OVERLAP_X;
    double yDelt = y_size - CCD_OVERLAP_Y;
    double deltY = 180.0/(yLen - 1);

    for( i = 0; i < yLen; i++ ) {

        struct Coordinate_SM coor;
        init_SKY_Coord( &coor );

        coor.dec = -90 + i * deltY;
        if(coor.dec>90) {
            coor.dec = 90;
        }

        double new_x_delt;
        int x_len = 0;

        if ( fabs(coor.dec) == 90 ) {
            new_x_delt = 90;
            x_len = 4;
        } else {
            // double borderAngle = (fabs(coor.dec) -  yDelt) / 57.29578;
            // new_x_delt = asin (sin(0.5 * xDelt / 57.29578)/cos(borderAngle) ) * 57.29578*2;
            double borderAngle = (fabs(coor.dec) -  yDelt) * PI_180;
            new_x_delt = asin( sin(0.5*xDelt*PI_180) / cos(borderAngle) ) * 360*M_1_PI;
            x_len = (int)(ceil(360/new_x_delt));
            new_x_delt = 360/(ceil(360/new_x_delt));
        }

        for (j = 0; j < x_len; j++) {
            coor.ra = j * new_x_delt;
            // if(coor.ra > 360) cnt++;
            coor.z = VECTOR_VALUE * sin(coor.dec * PI_180);
            coor.x = VECTOR_VALUE * cos(coor.dec * PI_180) * cos(coor.ra * PI_180);
            coor.y = VECTOR_VALUE * cos(coor.dec * PI_180) * sin(coor.ra * PI_180);
            coor.dist2orig = VECTOR_VALUE;

            double temp_z = VECTOR_VALUE;
            double temp_x = 0;
            double temp_y = 0;
            double t_a[3],t_b[3],n_v[3];
            t_a[0] = coor.x;
            t_a[1] = coor.y;
            t_a[2] = coor.z;
            t_b[0] = temp_x;
            t_b[1] = temp_y;
            t_b[2] = temp_z;
            crossMultiple(t_a, t_b, n_v);

            double n_v_norm = sqrt( n_v[0]*n_v[0] + n_v[1]*n_v[1] + n_v[2]*n_v[2] );

            coor.nx = n_v[0]/n_v_norm;
            coor.ny = n_v[1]/n_v_norm;
            coor.nz = n_v[2]/n_v_norm;

            double n30[3],n30_[3];
            get_solar_normal_30(coor.x,coor.y,coor.z,coor.nx,coor.ny,coor.nz,n30,n30_);
            coor.nx30 = n30[0];
            coor.ny30 = n30[1];
            coor.nz30 = n30[2];
            coor.nx30_ = n30_[0];
            coor.ny30_ = n30_[1];
            coor.nz30_ = n30_[2];
            double nn[3];
            crossMultiple(t_a, n_v, nn);
            double m_nn = sqrt( nn[0]*nn[0] + nn[1]*nn[1] + nn[2]*nn[2] );
            coor.nnx = nn[0]/m_nn;
            coor.nny = nn[1]/m_nn;
            coor.nnz = nn[2]/m_nn;

            double t_coor[3] = { coor.x, coor.y, coor.z };//    黄道坐标系坐标
            double r_coor[3];
            // CoordinateSpin(t_coor, r_coor, -23.4522, 3);
            CoordinateSpin_x(t_coor, r_coor, -23.4522);//转动到赤道坐标系
            double ll_coor[2];
            Cartesian2Equatorial(r_coor, ll_coor);

            double bn = 57.2957795
                        * asin(-0.8676660 * cos(ll_coor[0] * PI_180) * cos(ll_coor[1] * PI_180)
                              - 0.1980764 * sin(ll_coor[0] * PI_180) * cos(ll_coor[1] * PI_180)
                              + 0.4559840 * sin(ll_coor[1] * PI_180) );

            coor.gb = bn;
            coor.elat = coor.dec;

            coor.inDeepFlag = isInDeepSurveyArea(coor.ra, coor.dec, DEEP_DIAM);
            if( coor.inDeepFlag >= 0 ) {
                coor.isObserve = 1;
                coor.targetCoverNum = TARGET_COVER_NUM_DEEP;
                coor.maxCoverNum = MAX_COVER_NUM_DEEP;
                skyAreaNum++;
                deepSkyNum++;
            } else {
                coor.targetCoverNum = TARGET_COVER_NUM_LARGE;
                coor.maxCoverNum = MAX_COVER_NUM_LARGE;
//////////////////////////////////////////////////////////////////////////////////////
// 修改此处的天区划分方式：

#if defined(_USE_OLD_SKY_SPLIT_)
                if (   (coor.dec >= ECLIPTIC_LAT_LIMIT_N && coor.dec <= RANGE_DEC_N)
                    || (coor.dec >= RANGE_DEC_S && coor.dec <= ECLIPTIC_LAT_LIMIT_S) ) {
                    if( fabs(bn) < Low_Galaxy_Img ) {
                        coor.isObserve = 0;
                    } else {
                        coor.isObserve = 1;
                        skyAreaNum ++;
                    }
                } else {
                    coor.isObserve = 0;
                }
#elif defined(_USE_SKY_SPLIT_2_)
                coor.isObserve = 0;
            //  第一：按照原来的划分方式来判断
                if (   (coor.dec >= ECLIPTIC_LAT_LIMIT_N && coor.dec <= RANGE_DEC_N)
                    || (coor.dec >= RANGE_DEC_S && coor.dec <= ECLIPTIC_LAT_LIMIT_S) ) {
                    if( fabs(bn) < Low_Galaxy_Img ) {
                        coor.isObserve = 0;
                    } else {
                        coor.isObserve = 1;
                        // skyAreaNum ++;
                    }
                }

            //  第二：如果银纬、黄纬同时小于某个范围
                if( fabs(bn) < 10 ){
                    if( fabs(coor.dec) < 10 ){
						coor.isObserve = 1;
//						skyAreaNum++;
					}
                }

                if( coor.isObserve == 1 )
                    skyAreaNum++;
#else
                coor.isObserve = 0;
            //  第一：按照原来的划分方式来判断
                if (   (coor.dec >= ECLIPTIC_LAT_LIMIT_N && coor.dec <= RANGE_DEC_N)
                    || (coor.dec >= RANGE_DEC_S && coor.dec <= ECLIPTIC_LAT_LIMIT_S) ) {
                    if( fabs(bn) < Low_Galaxy_Img ) {
                        coor.isObserve = 0;
                    } else {
                        coor.isObserve = 1;
                        // skyAreaNum ++;
                    }
                }

            //  第二：如果银纬、黄纬同时小于某个范围
                if( fabs(bn) < Low_Galaxy_Img ){
                    if( (coor.dec >= ECLIPTIC_LAT_LIMIT_S-1.5) && (coor.dec <= ECLIPTIC_LAT_LIMIT_N+1.5) )
                    coor.isObserve = 1;
                }

                if( coor.isObserve == 1 )
                    skyAreaNum++;
#endif
//////////////////////////////////////////////////////////////////////////////////////
            }

            copySkyMapData( &coor, &skyMap[i * xLen + j] );

        }
    }

//  DEBUG
    // printf("--> finised debuging the condition test coord.ra > 360!\n");
    // printf("%d of %d have ra > 360\n",cnt,skyAreaNum);
    // exit(0);

    //printf("sky number is %d, deep num is %d \n",skyAreaNum,deepSkyNum);
    return skyAreaNum;
}


/*
*在划分天区的基础上添加相邻关系
*/
SKY_Coord* produceSkyArea(int* normal_sky_num) {
    //printf("%f      %f      %f \n",CCD_X,CCD_Y,CCD_X_DEEP);
    int i = 0, j = 0;
    int m;
    int i_sky = 0;
    
    double size[2];
    calculateSkySize( size, 
                      CCD_X, 
                      CCD_Y,
                      CCD_OVERLAP_X,
                      CCD_OVERLAP_Y );
    
    int xlen_normal = size[0];
    int ylen_normal = size[1];
    SKY_Coord* skyMap_norm = (SKY_Coord*)malloc( sizeof(CoordSM) * xlen_normal * ylen_normal );

    *normal_sky_num = SplitSkyArea(ylen_normal, xlen_normal, skyMap_norm, CCD_X, CCD_Y);

    SKY_Coord* skyMap = (SKY_Coord*)malloc( sizeof(CoordSM) * (*normal_sky_num) );

    for ( i = 0; i < xlen_normal*ylen_normal; i++ ) {
        if( skyMap_norm[i].isObserve == 1 ) {
            copySkyMapData( &skyMap_norm[i], &skyMap[i_sky] );
            skyMap[i_sky].id = i_sky;
            skyMap_norm[i].id = i_sky;
            i_sky++;
        }
    }

    double mdist = 0.26;
    for( i = 0; i < ylen_normal; i++ ) {

        int first_id = skyMap_norm[i * xlen_normal].id;
        int last_id = first_id;

        for ( j = 0; j < xlen_normal; j++ ) {

            int n_i = i * xlen_normal + j;

            if( skyMap_norm[n_i].isObserve == 1 ) {

                int c_id = skyMap_norm[n_i].id;
                if( first_id == -1 ) {
                    first_id = c_id;
                }

                if( i-1 >= 0 ) {

                    int down_id=0;
                    double c_dist = 10000;

                    for ( m=-10; m <=10; m++ ) {
                        if ( j+m >= 0 && j+m < xlen_normal ) {
                            int candi_id = (i-1) * xlen_normal + j+m;
                            double candi_dist = calculateAngle( skyMap_norm[n_i].ra,
                                                                skyMap_norm[n_i].dec,
                                                                skyMap_norm[candi_id].ra,
                                                                skyMap_norm[candi_id].dec );

                            if( candi_dist < c_dist ) {
                                c_dist = candi_dist;
                                down_id = candi_id;
                            }

                        }
                    }

                    if( c_dist < mdist ) {
                        if( skyMap_norm[down_id].isObserve == 1 ) {
                            skyMap[c_id].down_neighbour = skyMap_norm[down_id].id;
                        } else {
                            skyMap[c_id].down_neighbour = -1;
                        }
                    } else {
                        skyMap[c_id].down_neighbour = -1;
                    }
                }

                if ( j-1 >= 0 ) {
                    int left_id = i *xlen_normal + j - 1;
                    if(skyMap_norm[left_id].isObserve == 1) {
                        skyMap[c_id].left_neighbour = skyMap_norm[left_id].id;
                    } else {
                        skyMap[c_id].left_neighbour = -1;
                    }
                }

                if ( i+1 < ylen_normal ) {
                    int up_id=0;
                    double c_dist = 10000;
                    for( m=-10; m <=10; m++ ) {
                        if ( j+m >= 0 && j+m < xlen_normal ) {
                            int candi_id = (i + 1) *xlen_normal + j+m;
                            double candi_dist = calculateAngle( skyMap_norm[n_i].ra,
                                                                skyMap_norm[n_i].dec,
                                                                skyMap_norm[candi_id].ra,
                                                                skyMap_norm[candi_id].dec );
                            if( candi_dist < c_dist ) {
                                c_dist = candi_dist;
                                up_id = candi_id;
                            }

                        }
                    }

                    if( c_dist < mdist ) {
                        if( skyMap_norm[up_id].isObserve == 1 ) {
                            skyMap[c_id].up_neighbour = skyMap_norm[up_id].id;
                        } else {
                            skyMap[c_id].up_neighbour = -1;
                        }
                    } else {
                        skyMap[c_id].up_neighbour = -1;
                    }

                }

                if ( j+1 < xlen_normal ) {
                    int right_id = i *xlen_normal + j + 1;
                    if( skyMap_norm[right_id].isObserve == 1 ) {
                        skyMap[c_id].right_neighbour = skyMap_norm[right_id].id;
                    } else {
                        skyMap[c_id].right_neighbour = -1;
                    }

                }

                last_id = c_id;
            }
        }

        if( first_id != -1 && last_id != -1 ) {
            if( skyMap[first_id].ra+360-skyMap[last_id].ra < CCD_X ) {
                skyMap[first_id].left_neighbour = last_id;
                skyMap[last_id].right_neighbour = first_id;
            }
        }
    }

    // printf("%f      %f      %f \n",CCD_X,CCD_Y,CCD_X_DEEP);
    free(skyMap_norm);

    return skyMap;
}

int count_sky_num( SKY_Coord* sky, int sky_num, int rank ){
    int num=0,i=0;
    for( i=0; i<sky_num; i++ ){
        if( sky[i].isObserve == 1
        &&  sky[i].flag < sky[i].maxCoverNum )
            num++;
    }

    return num;
}

int* update_sky_id_tracker( SKY_Coord* sky, int sky_num, int rank, int* remained_num ){
    // printf("DEBUG:: calling update_sky_id_tracker(***)\n");
    int i,j;
    int cnt = 0;
    for( i=0; i<sky_num; i++ ){
        if( sky[i].isObserve == 1
        &&  sky[i].flag < sky[i].maxCoverNum )
            cnt++;
    }

    *remained_num = cnt;

    int* tracker = (int*)malloc(sizeof(int)*cnt);
    j = 0;
    for( i=0; i<sky_num; i++ ){
        if( sky[i].isObserve == 1
        &&  sky[i].flag < sky[i].maxCoverNum ){
            tracker[j] = sky[i].id;
            j++;
        }
    }

    if( j != cnt ){
        if(rank == 0) {
            printf("Error happened when reducing the # of observable sky!!!\n");
            printf("--> cnt = %d\n", cnt);
            printf("--> j   = %d\n", j);
        }
        MPI_Finalize();
        exit(0);
    }
    
    return tracker;
}
