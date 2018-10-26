#include "SurveySim.h"

void parseConfigFile(char* fileName) {
    FILE *configFile;
    configFile = fopen(fileName,"r");

    if(configFile == NULL) {
        printf("No such config file: %s !!!\n",fileName);
        exit(1);
    }

    char line[MAXLINE_STRING];
    char *p;
    char start[3] = "#";

    while (fgets(line, MAXLINE_STRING, configFile) != NULL) {

        p = strtok(line, " \n\t");

        if (p==NULL) {
            continue;
        }
        if( strcmp(p,"")==0 ) {
            continue;
        }
        if(p[0] == start[0]) {
            continue;
        }

        if(strcmp((char*)"START_TIME",p)==0) {
            p = strtok(NULL," \n\t");
            if(p == NULL) {
                continue;
            }
            if (p[0] == start[0]) {
                continue;
            }
            startTime = atof(p);
        } else if(strcmp((char*)"END_TIME",p)==0) {
            p = strtok(NULL," \n\t");
            if(p == NULL) {
                continue;
            }
            if (p[0] == start[0]) {
                continue;
            }
            endTime = atof(p);
        } else if(strcmp((char*)"JPL_FILE",p)==0) {
            p = strtok(NULL," \n\t");
            if(p == NULL) {
                continue;
            }
            if (p[0] == start[0]) {
                continue;
            }
            strcpy(jpl_fileName,p);
        } else if(strcmp((char*)"RESULT_FILE",p)==0) {
            p = strtok(NULL," \n\t");
            if(p == NULL) {
                continue;
            }
            if (p[0] == start[0]) {
                continue;
            }
            strcpy(result_fileName,p);
        } else if(strcmp((char*)"STATUS_FILE_IMG",p)==0) {
            p = strtok(NULL," \n\t");
            if(p == NULL) {
                continue;
            }
            if (p[0] == start[0]) {
                continue;
            }
            strcpy(status_img_fileName,p);
        } else if(strcmp((char*)"DEEP_AREA",p)==0) {
            p = strtok(NULL," \n\t");
            if(p == NULL) {
                continue;
            }
            if (p[0] == start[0]) {
                continue;
            }
            DEEP_AREA = atof(p);
        } else if(strcmp((char*)"ULTRALDEPP_AREA",p)==0) {
            p = strtok(NULL," \n\t");
            if(p == NULL) {
                continue;
            }
            if (p[0] == start[0]) {
                continue;
            }
            ULTRALDEPP_AREA = atof(p);
        } else if(strcmp((char*)"STATUS_INTVAL",p)==0) {
            p = strtok(NULL," \n\t");
            if(p == NULL) {
                continue;
            }
            if (p[0] == start[0]) {
                continue;
            }
            F_STATUS = atoi(p);
        } else if(strcmp((char*)"CCD_X",p)==0) {
            p = strtok(NULL," \n\t");
            if(p == NULL) {
                continue;
            }
            if (p[0] == start[0]) {
                continue;
            }
            CCD_X = atof(p);
        } else if(strcmp((char*)"CCD_Y",p)==0) {
            p = strtok(NULL," \n\t");
            if(p == NULL) {
                continue;
            }
            if (p[0] == start[0]) {
                continue;
            }
            CCD_Y = atof(p);
        } else if(strcmp((char*)"CCD_X_DEEP",p)==0) {
            p = strtok(NULL," \n\t");
            if(p == NULL) {
                continue;
            }
            if (p[0] == start[0]) {
                continue;
            }
            CCD_X_DEEP = atof(p);
        } else if(strcmp((char*)"CCD_OVERLAP_X",p)==0) {
            p = strtok(NULL," \n\t");
            if(p == NULL) {
                continue;
            }
            if (p[0] == start[0]) {
                continue;
            }
            CCD_OVERLAP_X = atof(p);
        } else if(strcmp((char*)"CCD_OVERLAP_Y",p)==0) {
            p = strtok(NULL," \n\t");
            if(p == NULL) {
                continue;
            }
            if (p[0] == start[0]) {
                continue;
            }
            CCD_OVERLAP_Y = atof(p);
        } else if(strcmp((char*)"CCD_OVERLAP_X_DEEP",p)==0) {
            p = strtok(NULL," \n\t");
            if(p == NULL) {
                continue;
            }
            if (p[0] == start[0]) {
                continue;
            }
            CCD_OVERLAP_X_DEEP = atof(p);
        } else if(strcmp((char*)"SUN_POINT_ANGLE",p)==0) {
            p = strtok(NULL," \n\t");
            if(p == NULL) {
                continue;
            }
            if (p[0] == start[0]) {
                continue;
            }
            double sun_point_angle = atof(p);
            COS_SUN_POINT_ANGLE = cos(sun_point_angle * PI_180);
        } else if(strcmp((char*)"MOON_POINT_ANGLE",p)==0) {
            p = strtok(NULL," \n\t");
            if(p == NULL) {
                continue;
            }
            if (p[0] == start[0]) {
                continue;
            }
            double moon_point_angle = atof(p);
            COS_MOON_POINT_ANGLE = cos(moon_point_angle * PI_180);
        } else if(strcmp((char*)"SUN_PLANE_ANGLE",p)==0) {
            p = strtok(NULL," \n\t");
            if(p == NULL) {
                continue;
            }
            if (p[0] == start[0]) {
                continue;
            }
            double sun_plane_angle = atof(p);
            COS_SUN_PLANE_ANGLE = cos(sun_plane_angle * PI_180);
        } else if(strcmp((char*)"POINT_EARTH_LLIMB_ANGLE_MAX",p)==0) {
            p = strtok(NULL," \n\t");
            if(p == NULL) {
                continue;
            }
            if (p[0] == start[0]) {
                continue;
            }
            double p_e_angle = atof(p);
            COS_POINT_ZENITH_ANGLE_LIGHT_Max = cos((110-p_e_angle) * PI_180);
        } else if(strcmp((char*)"POINT_EARTH_LLIMB_ANGLE_MIN",p)==0) {
            p = strtok(NULL," \n\t");
            if(p == NULL) {
                continue;
            }
            if (p[0] == start[0]) {
                continue;
            }
            double p_e_angle = atof(p);
            COS_POINT_ZENITH_ANGLE_LIGHT_Min = cos((110-p_e_angle) * PI_180);
        } else if(strcmp((char*)"POINT_EARTH_DLIMB_ANGLE",p)==0) {
            p = strtok(NULL," \n\t");
            if(p == NULL) {
                continue;
            }
            if (p[0] == start[0]) {
                continue;
            }
            double p_e_angle = atof(p);
            COS_POINT_ZENITH_ANGLE_DARK = cos((110-p_e_angle) * PI_180);
        } else if(strcmp((char*)"EXTIME_ULTRAL_DEEP",p)==0) {
            p = strtok(NULL," \n\t");
            if(p == NULL) {
                continue;
            }
            if (p[0] == start[0]) {
                continue;
            }
            EXTIME_DEEP = atof(p);
        } else if(strcmp((char*)"EXTIME_DEEP",p)==0) {
            p = strtok(NULL," \n\t");
            if(p == NULL) {
                continue;
            }
            if (p[0] == start[0]) {
                continue;
            }
            EXTIME = atof(p);
        } else if(strcmp((char*)"EXTIME_SPEC",p)==0) {
            p = strtok(NULL," \n\t");
            if(p == NULL) {
                continue;
            }
            if (p[0] == start[0]) {
                continue;
            }
            EXTIME_SPEC = atof(p);
        } else if(strcmp((char*)"EXTIME_MIN",p)==0) {
            p = strtok(NULL," \n\t");
            if(p == NULL) {
                continue;
            }
            if (p[0] == start[0]) {
                continue;
            }
            EXTIME_G_E = atof(p);
        } else if(strcmp((char*)"RANGE_NORTH_SEMISPHERE_N",p)==0) {
            p = strtok(NULL," \n\t");
            if(p == NULL) {
                continue;
            }
            if (p[0] == start[0]) {
                continue;
            }
            RANGE_DEC_N = atof(p);
        } else if(strcmp((char*)"RANGE_NORTH_SEMISPHERE_S_IMG",p)==0) {
            p = strtok(NULL," \n\t");
            if(p == NULL) {
                continue;
            }
            if (p[0] == start[0]) {
                continue;
            }
            ECLIPTIC_LAT_LIMIT_N = atof(p);
        } else if(strcmp((char*)"RANGE_SOUTH_SEMISPHERE_S",p)==0) {
            p = strtok(NULL," \n\t");
            if(p == NULL) {
                continue;
            }
            if (p[0] == start[0]) {
                continue;
            }
            RANGE_DEC_S = atof(p);
        } else if(strcmp((char*)"RANGE_SOUTH_SEMISPHERE_N_IMG",p)==0) {
            p = strtok(NULL," \n\t");
            if(p == NULL) {
                continue;
            }
            if (p[0] == start[0]) {
                continue;
            }
            ECLIPTIC_LAT_LIMIT_S = atof(p);
        } else if(strcmp((char*)"Low_Galaxy_Img",p)==0) {
            p = strtok(NULL," \n\t");
            if(p == NULL) {
                continue;
            }
            if (p[0] == start[0]) {
                continue;
            }
            Low_Galaxy_Img = atof(p);
        } else if(strcmp((char*)"Galaxy_B_Sec_Low",p)==0) {
            p = strtok(NULL," \n\t");
            if(p == NULL) {
                continue;
            }
            if (p[0] == start[0]) {
                continue;
            }
            Galaxy_B_Sec_Low = atof(p);
        } else if(strcmp((char*)"Ecliptic_Lat_Sec_Low",p)==0) {
            p = strtok(NULL," \n\t");
            if(p == NULL) {
                continue;
            }
            if (p[0] == start[0]) {
                continue;
            }
            Ecliptic_Lat_Sec_Low = atof(p);
        } else if(strcmp((char*)"PANEL_TRANSE_ANGLE",p)==0) {
            p = strtok(NULL," \n\t");
            if(p == NULL) {
                continue;
            }
            if (p[0] == start[0]) {
                continue;
            }
            double panel_trans_angle = atof(p);
            PANEL_TRANSE_ANGLE = panel_trans_angle;
        }
    //  The following are added by Youhua Xu
         else if( strcmp((char*)"BETA_ANGLE",p) == 0 ){
            p = strtok(NULL," \n\t");
            if(p == NULL) {
                continue;
            }
            if (p[0] == start[0]) {
                continue;
            }
            double beta_angle = atof(p);
            BETA_ANGLE = beta_angle;
        } else if( strcmp((char*)"HIGH_LATITUDE_PRIOR_TIME",p) == 0 ){
            p = strtok(NULL," \n\t");
            if(p == NULL) {
                continue;
            }
            if (p[0] == start[0]) {
                continue;
            }
            double high_latitude_prior_time = atof(p);
            HIGH_LATITUDE_PRIOR_TIME = high_latitude_prior_time;
        } else if( strcmp((char*)"CONTINOUS_OBS_START_TIME",p) == 0 ){
            p = strtok(NULL," \n\t");
            if(p == NULL) {
                continue;
            }
            if (p[0] == start[0]) {
                continue;
            }
            CONTINOUS_OBS_START_TIME = atof(p);
        }
        else if( strcmp((char*)"JUMP_TIME",p) == 0 ){
            p = strtok(NULL," \n\t");
            if(p == NULL) {
                continue;
            }
            if (p[0] == start[0]) {
                continue;
            }
            jump_time = atof(p);
            // printf("*** JUMP_TIME = %g ***\n",jump_time);
        }
        else if( strcmp((char*)"DEC60_PRIOR_TIME",p) == 0 ){
            p = strtok(NULL," \n\t");
            if(p == NULL) {
                continue;
            }
            if (p[0] == start[0]) {
                continue;
            }
            DEC60_PRIOR_TIME = atof(p);
        }
    }
}

void printConfig() {
    printf("start time: %f\n",startTime);
    printf("end time: %f\n",endTime);
    printf("jpl file name: %s\n",jpl_fileName);
    printf("result file name: %s\n",result_fileName);
    printf("status img file name: %s\n",status_img_fileName);
    //printf("status spec file name: %s\n",status_spec_fileName);
    printf("DEEP_AREA: %f\n",DEEP_AREA);
    printf("ULTRALDEPP_AREA: %f\n",ULTRALDEPP_AREA);
    //printf("SPEC_AREA: %f\n",SPEC_AREA);
    printf("STATUS_INTVAL: %d\n",F_STATUS);
    printf("CCD_X: %f\n",CCD_X);
    printf("CCD_Y: %f\n",CCD_Y);
    printf("CCD_X_DEEP: %f\n",CCD_X_DEEP);
    printf("CCD_OVERLAP_X: %f\n",CCD_OVERLAP_X);
    printf("CCD_OVERLAP_Y: %f\n",CCD_OVERLAP_Y);
    printf("CCD_OVERLAP_X_DEEP: %f\n",CCD_OVERLAP_X_DEEP);
    printf("COS_SUN_POINT_ANGLE: %f\n",COS_SUN_POINT_ANGLE);
    printf("COS_MOON_POINT_ANGLE: %f\n",COS_MOON_POINT_ANGLE);
    printf("COS_SUN_PLANE_ANGLE: %f\n",COS_SUN_PLANE_ANGLE);
    printf("COS_POINT_ZENITH_ANGLE_LIGHT_Max: %f\n",COS_POINT_ZENITH_ANGLE_LIGHT_Max);
    printf("COS_POINT_ZENITH_ANGLE_LIGHT_Min: %f\n",COS_POINT_ZENITH_ANGLE_LIGHT_Min);
    printf("COS_POINT_ZENITH_ANGLE_DARK: %f\n",COS_POINT_ZENITH_ANGLE_DARK);
    printf("EXTIME_DEEP: %f\n",EXTIME_DEEP);
    printf("EXTIME: %f\n",EXTIME);
    printf("EXTIME_SPEC: %f\n",EXTIME_SPEC);
    printf("EXTIME_G_E: %f\n",EXTIME_G_E);
    printf("RANGE_NORTH_SEMISPHERE_N: %f\n",RANGE_DEC_N);
    printf("RANGE_NORTH_SEMISPHERE_S_IMG: %f\n",ECLIPTIC_LAT_LIMIT_N);
    printf("RANGE_SOUTH_SEMISPHERE_S: %f\n",RANGE_DEC_S);
    printf("RANGE_SOUTH_SEMISPHERE_N_IMG: %f\n",ECLIPTIC_LAT_LIMIT_S);
    printf("Low_Galaxy_Img: %f\n",Low_Galaxy_Img);
    printf("Galaxy_B_Sec_Low: %f\n",Galaxy_B_Sec_Low);
    printf("Low_EcliEcliptic_Lat_Sec_Lowptic: %f\n",Ecliptic_Lat_Sec_Low);
    printf("PANEL_TRANSE_ANGLE: %f\n",PANEL_TRANSE_ANGLE);

    printf("\n============================================================\n");
    printf("=== the following variables are newly added by Youhua Xu ===\n");
    printf("============================================================\n");
    printf("BETA_ANGLE: %f\n", BETA_ANGLE);
    printf("HIGH_LATITUDE_PRIOR_TIME: %f\n", HIGH_LATITUDE_PRIOR_TIME);
    printf("CONTINOUS_OBS_START_TIME: %f\n", CONTINOUS_OBS_START_TIME);
    printf("DEC60_PRIOR_TIME: %f\n",DEC60_PRIOR_TIME);
    printf("jump_time: %f\n",jump_time);

    printf("==> End of printfConfig()\n\n");
}

void parseStatusFile( char* fileName, 
                      double* curTime, 
                      int* point_id,
                      double *surveyArea,
                      double* surveyArea2,
                      double* surveyArea3,
                      double* surveyArea4,
                      double* surveyAreaG,
                      double* surveyAreaG2,
                      double* deepSurveyArea,
                      double* toearthTime,
                      int* excuteNum, 
                      SKY_Coord* skyMap ) {
    FILE *statusFile;
    statusFile = fopen(fileName,"r");

    if(statusFile == NULL) {
        printf("No such file1 !!!\n");
        exit(1);
    }
    // int num = 0;
    char line[MAXLINE_STRING];
    int elemNum = 0, i = 0;  //, j = 0;
    char *p;
    char start[3] = "#";
    while (fgets(line, MAXLINE_STRING, statusFile) != NULL) {
        p = strtok(line, " \n");
        if(p[0] != start[0]) {
            break;
        }
    }

    if( p == NULL) {
        printf("error in read status file \n");
        exit(1);
    }
    *curTime = atof(p) ;
    p = strtok(NULL, " \n");
    if( p == NULL) {
        printf("error in read status file \n");
        exit(1);
    }
    printf(" time is %f \n",*curTime);
    *point_id = atoi(p) ;
    p = strtok(NULL, " \n");
    if( p == NULL) {
        printf("error in read status file \n");
        exit(1);
    }
    *surveyArea = atof(p) ;
    p = strtok(NULL, " \n");
    if( p == NULL) {
        printf("error in read status file \n");
        exit(1);
    }
    *surveyArea2 = atof(p) ;
    p = strtok(NULL, " \n");
    if( p == NULL) {
        printf("error in read status file \n");
        exit(1);
    }
    *surveyArea3 = atof(p) ;
    p = strtok(NULL, " \n");
    if( p == NULL) {
        printf("error in read status file \n");
        exit(1);
    }
    *surveyArea4 = atof(p) ;
    p = strtok(NULL, " \n");
    if( p == NULL) {
        printf("error in read status file \n");
        exit(1);
    }
    *surveyAreaG = atof(p) ;
    p = strtok(NULL, " \n");
    if( p == NULL) {
        printf("error in read status file \n");
        exit(1);
    }
    *surveyAreaG2 = atof(p) ;
    p = strtok(NULL, " \n");
    if( p == NULL) {
        printf("error in read status file \n");
        exit(1);
    }
    *deepSurveyArea = atof(p) ;
    // p = strtok(NULL, " \n");
    // if( p == NULL){
    // 		printf("error in read status file \n");
    // 		exit(1);
    // }
    // *lapTime = atof(p) ;
    p = strtok(NULL, " \n");
    if( p == NULL) {
        printf("error in read status file \n");
        exit(1);
    }
    *toearthTime = atof(p) ;
    p = strtok(NULL, " \n");
    if( p == NULL) {
        printf("error in read status file \n");
        exit(1);
    }
    *excuteNum = atoi(p) ;
    while (fgets(line, MAXLINE_STRING, statusFile) != NULL) {
        p = strtok(line, " \n");
        if(p[0] != start[0]) {
            break;
        }
    }
    if( p == NULL) {
        printf("error in read status file \n");
        exit(1);
    }
    elemNum = atoi(p);

    for(i = 0; i<elemNum; i ++) {
        fgets(line, MAXLINE_STRING, statusFile);
        if(line != NULL) {
            int id = 0;
            int flag = 0;
            p = strtok(line, " \n");
            if( p == NULL) {
                printf("error in read status file \n");
                exit(1);
            }
            id = atoi(p);
            p = strtok(NULL, " \n");
            if( p == NULL) {
                printf("error in read status file \n");
                exit(1);
            }
            flag = atoi(p);
            skyMap[id].flag = flag;
        }
    }
    fclose(statusFile);
}
