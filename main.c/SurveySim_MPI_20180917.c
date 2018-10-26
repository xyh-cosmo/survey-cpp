#include "SurveySim.h"

int main(int argc, char *argv[]) {

    int p_rank;
    int p_size = 0;
    MPI_Init(NULL, NULL );
    MPI_Comm_rank(MPI_COMM_WORLD, &p_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p_size);

    init_default();

    int i = 0;
    // int isReadStatus = 0;
    double curTime = 0;

    FILE *outf = NULL;
    FILE *outf_fail = NULL;
    // FILE *statusFile = NULL;


    if (argc < 2) {
        fprintf(stderr, "\nFormat:\n\n         %s config-file\n\n", argv[0]);
        MPI_Finalize();
        return 0;
    }

    parseConfigFile(argv[1]);

    if(p_rank == 0) {
        printConfig();
        printf("interval time (JUMP_TIME) is %f \n",jump_time);
    }

    //  读取星表数据
    infp = NULL;
    if ((infp = fopen(jpl_fileName, "r")) == NULL ) {
        fprintf(stderr, "\nERROR: Can't open ephemeris file %s for input.\n\n", argv[1]);
        MPI_Finalize();
        return 0;
    }

    int outfile_status1 = 0, outfile_status2 = 0;

    double earthStartTime, initalGreenWithRa = 0;
    if (p_rank == 0) {
        //计算某一时刻地区本初子午线与赤道交点午夜时在天球赤道坐标系中的位置
        earthStartTime = floor(startTime - 0.5) + 0.5; // midnight time (unit day)
        double sunMid[3];
        locateSun(infp, earthStartTime, sunMid); // 太阳的笛卡尔坐标，以地球为原点(J2000)

        double reEarth[2]; // (ra,dec) of Sun
        Cartesian2Equatorial(sunMid, reEarth);  // get (ra,dec) of Sun

        if (reEarth[0] >= 180) {
            initalGreenWithRa = reEarth[0] - 180;
        } else {
            initalGreenWithRa = reEarth[0] + 180;
        }

        char r_file_name[100] = "";
        strcat(r_file_name,result_fileName);
        strcat(r_file_name,".dat");

        char r_fail_file_name[100] = "";
        strcat(r_fail_file_name,result_fileName);
        strcat(r_fail_file_name,"_fail_stat.dat");

        if ((outf = fopen(r_file_name, "r")) == NULL ) { /* File doesn't exist */
            if ((outf = fopen(r_file_name, "w")) == NULL ) {
                fprintf(stderr, "\nERROR: Can't open %s for output.\n\n", r_file_name);
                outfile_status1 = 999;
            }
        } else {
            fprintf(stderr, "\nERROR: output file %s already exists.\n", r_file_name);
            outfile_status1 = 999;
        }

        if ((outf_fail = fopen(r_fail_file_name, "r")) == NULL ) { /* File doesn't exist */
            if ((outf_fail = fopen(r_fail_file_name, "w")) == NULL ) {
                fprintf(stderr, "\nERROR: Can't open %s for output.\n\n", r_fail_file_name);
                outfile_status2 = 999;
            }
        } else {
            fprintf(stderr, "\nERROR: output file %s already exists.\n", r_fail_file_name);
            outfile_status2 = 999;
        }
    }

    MPI_Bcast(&outfile_status1, 1, MPI_INT, 0, MPI_COMM_WORLD );
    MPI_Bcast(&outfile_status2, 1, MPI_INT, 0, MPI_COMM_WORLD );

    if( outfile_status1 != 0 ){
        MPI_Finalize();
        return 0;
    }

    if( outfile_status2 != 0 ){
        MPI_Finalize();
        return 0;
    }

    MPI_Bcast(&initalGreenWithRa, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast(&earthStartTime, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    // MPI_Bcast(&isReadStatus, 1, MPI_INT, 0, MPI_COMM_WORLD );

//  读取所有的轨道数据
    orbitData = (double**)malloc(sizeof(double*) * Orbit_File_Num);
    readAllorbitsFile(orbitData, orbitDataNum);

    if(p_rank == 0) {
        printf("start time is %f \n", startTime);
        printf("rank size is %d \n", p_size);
    }

//  对天区进行划分
    int normal_sky_num = 0;
    SKY_Coord* skyMap = produceSkyArea(&normal_sky_num);

    if(p_rank == 0) {
        printf("survey sky num is %d\n",normal_sky_num);

#if defined(_ENABLE_SAVE_SKYMAP_)
		FILE* skyF = NULL;
		char skyF_fname[1024];
		sprintf(skyF_fname,"%s_skymap.dat",result_fileName);
		skyF = fopen(skyF_fname,"w");

		int ii;
		for( ii=0; ii< normal_sky_num; ii++ ){
			fprintf(skyF,"%d %f %f %f %d %d %d %d\n",
					skyMap[ii].id,
					skyMap[ii].dec,
					skyMap[ii].ra,
					skyMap[ii].inDeepFlag,
					skyMap[ii].left_neighbour,
					skyMap[ii].right_neighbour,
					skyMap[ii].up_neighbour,
					skyMap[ii].down_neighbour );
		}

        printf("==> finished writing skymap ...\n");
		fclose(skyF);
#endif
    }

	MPI_Barrier(MPI_COMM_WORLD);

//  更新可观测天区数目的信息
    if( sky_id_tracker!=NULL )
        free(sky_id_tracker);
        
    sky_id_tracker = update_sky_id_tracker( skyMap, 
                                            normal_sky_num, 
                                            p_rank,
                                            &sky_num_remained );

    if( sky_num_remained != normal_sky_num ){
        if( p_rank == 0 ){
            printf("Error happened when initializing sky_id_tracker!!!\n");
            printf("--> sky_num_remained = %d\n",sky_num_remained);
            printf("--> normal_sky_num = %d\n", normal_sky_num);
        }
        MPI_Finalize();
        exit(0);
    }

    curTime = startTime;

    double coor_e[3]; // initial position of satellite
    double dist_sat = locateSat1(coor_e, curTime, orbitData, orbitDataNum,Orbit_File_Num);
    while( dist_sat < 0 ){
        curTime += jump_time / 86400.0;
        dist_sat = locateSat1(coor_e,curTime,orbitData,orbitDataNum,Orbit_File_Num);
    }

    double coor_s[3];
    CoordinateSpinEquatorial2Ecliptic(coor_e, coor_s);  //从赤道坐标系转动到黄道坐标系
    double p_init[2];
    Cartesian2Equatorial(coor_s, p_init);
    double min_init_point = MAX_VALUE;
    int init_id = 0;
    for(i = 0; i < (normal_sky_num); i ++) {// 搜寻里初始指向最接近的天区
        double dist_now = calculateAngle(skyMap[i].ra,skyMap[i].dec,p_init[0], p_init[1]);
        if(dist_now < min_init_point) {
            min_init_point=dist_now;
            init_id = i;
        }
    }

    int point_id = init_id;
    int excuteNum = 0;
    int noExcuteNum = 0;
    int m = 0; //, n = 0;
    double toearthTime = 0;

    time_t stime, etime;
    time(&stime);

//==============================================================================

    cmg_list = (struct CMG_List*)malloc(sizeof(struct CMG_List));
    initCMGList(cmg_list);

    battery_q = BATTERY_MAX;    // sun plane, energy balance（全局变量）
    shadow_st = startTime;      // 进入阴影区的时刻（全局变量）
    shadow_et = startTime;      // 离开阴影区的时刻（全局变量）

    double area1 = 0; // 正常成像面积
    double area2 = 0; // 极深度巡天面积

#if defined(_ENABLE_MARK_OBS_POINT_)
    double deltArea;
    double ccd_pos_in_focus[18][2];
    struct FourTree* Trees12 = NULL;

    if( p_rank == 0 ) {
        Trees12 = (struct FourTree*)malloc(sizeof(struct FourTree)*12);
        int layer = 10;
        produceHealPixeTrees(Trees12, layer);
        deltArea = 41253.0/(12.0*pow(4,layer));
        init_ccd_pos(ccd_pos_in_focus);
    }
#endif
    MPI_Barrier(MPI_COMM_WORLD);

//  ====================================================
	int NMAX = 1000;
	int cnt_good = 0;
	int cnt_fail = 0;
//    double* outputresult = (double*)malloc(sizeof(double)*16);
    double* outputresult[NMAX];
    double dec[NMAX];
    double ra[NMAX];
    double deepflag[NMAX];
    int sunside[NMAX];
    double cmgtotal[NMAX];
    int posflag_[NMAX];
    double battery_q_[NMAX];
    int id_[NMAX];

    double area1_[NMAX];
    double area2_[NMAX];

#if defined(_ENABLE_NEW_SEARCH4_)
    int try_steps[NMAX];
#endif

	// output_fail 用于存储和统计搜索可观测天区失败的原因
    // 0--- time
    // 1--- sun
    // 2--- moon
    // 3--- earth
    // 4--- energy
    // 5--- cmg
    // 6--- sun_plane
    // 7--- suppass cover number
    // 8--- SAA
    // 9--- is In sun side
    
//    double* output_fail = (double*)malloc(sizeof(double)*10);
	double* output_fail[NMAX];

	for( m=0; m<NMAX; m++ ){
		outputresult[m] = (double*)malloc(sizeof(double)*16);
		output_fail[m] = (double*)malloc(sizeof(double)*10);
	}

    // MPI_Finalize();
    // exit(0);

//////////////////////////////////////////////////////
// 这里的while循环才是最为耗时的部分!所以在测试的时候从这个地方
// 开始计时,以便得到更精确的测试运行所需要的时间.
    time(&stime);

#if defined(_ENABLE_NEW_SEARCH3_)
    double cur_ra = skyMap[init_id].ra;
    double cur_dec= skyMap[init_id].dec;
#endif

#if defined(_ENABLE_NEW_SEARCH4_)
double lastObsEndTime;
#endif

    int sky_num_remained_tmp = sky_num_remained;
    int sky_id_tracker_cnt = 0;
    while ( curTime < endTime ) {

        sky_id_tracker_cnt++;  //放在这里可以确保这个计数器在每一次步进时都会增加1

        // if( sky_num_remained < 1 ){
        if( sky_num_remained_tmp < 1 ){
            if( p_rank == 0 ){
                printf("*** There is no enough observable sky area left, so let's stop here! ***\n");
            }
            break;
        }

        // 检查是否需要停靠
        for( m = 0; m < Stop_Time_Num; m ++ ) {
            if( curTime > stopTime[m][1] ) {
                continue;
            } else {
                if( curTime < stopTime[m][0] ) {
                    break;
                } else {
                    while( curTime <= stopTime[m][1] ) {
                        curTime = curTime + jump_time / 86400.0;
#if defined(_ENABLE_NEW_SEARCH4_)
                        lastObsEndTime = curTime;
#endif
                        battery_q = BATTERY_MAX;
                        free(cmg_list);
                        cmg_list = (struct CMG_List*)malloc(sizeof(struct CMG_List));
                        initCMGList(cmg_list);
                    }
                    break;
                }
            }
        }


//      ========================================================================
#if defined(_ENABLE_BETA_ANGLE_)
        double sunAngle = getSun2OrbitAngle1(infp, curTime,orbitData, orbitDataNum,Orbit_File_Num);

    //  当轨道面与太阳夹角足够小的时候,直接开启高纬度优先的策略
        // use_high_altitude_prior = 0;
        // if( fabs(sunAngle) < BETA_ANGLE + 20 ){
        //     use_high_altitude_prior = 1;
        // }

        while ( fabs(sunAngle) < BETA_ANGLE ) {
            curTime = curTime + jump_time / 86400.0;
#if defined(_ENABLE_NEW_SEARCH4_)
            lastObsEndTime = curTime;
#endif
            toearthTime += jump_time;
            sunAngle = getSun2OrbitAngle1(infp, curTime,orbitData, orbitDataNum,Orbit_File_Num);
            battery_q = BATTERY_MAX;
       
            free(cmg_list);
            cmg_list = (struct CMG_List*)malloc(sizeof(struct CMG_List));
            initCMGList(cmg_list);
       
            if(curTime >= endTime) {
                break;
            }
        }

        for(m = 0; m < Stop_Time_Num; m ++) {
            if(curTime > stopTime[m][1]) {
                continue;
            } else {
                if(curTime < stopTime[m][0]) {
                    break;
                } else {
                    while(curTime<=stopTime[m][1]) {
                        curTime = curTime + jump_time / 86400.0;
#if defined(_ENABLE_NEW_SEARCH4_)
                        lastObsEndTime = curTime;
#endif
                        battery_q = BATTERY_MAX;
                        free(cmg_list);
                        cmg_list = (struct CMG_List*)malloc(sizeof(struct CMG_List));
                        initCMGList(cmg_list);
                    }
                    break;
                }
            }
        }
#endif
//      ========================================================================

        if ( curTime >= endTime ) {
            break;
        }

        if ( curTime > shadow_et ) {
        //  如果当前处于阴影区以外了，那就获取下一次进入和离开阴影区的时刻。
            aquireShadowTime( &shadow_st, 
                              &shadow_et, 
                              curTime );
        }

        setDoubleArrayZero( output_fail[cnt_fail], 10 );


#if defined(_ENABLE_NEW_SEARCH_)
        FindNewPointByAllSearch( &curTime,
                                 &point_id,
                                 skyMap,
                                 normal_sky_num,
                                 infp,
                                 earthStartTime,
                                 initalGreenWithRa,
                                 outputresult[cnt_good],
                                 output_fail[cnt_fail],
                                 orbitData,
                                 orbitDataNum,
                                 p_rank,
                                 p_size );
#elif defined(_ENABLE_NEW_SEARCH2_)     
        FindNewPointByAllSearch2(&curTime,
                                 &point_id,
                                 skyMap,
                                 normal_sky_num,
                                 infp,
                                 earthStartTime,
                                 initalGreenWithRa,
                                 outputresult[cnt_good],
                                 output_fail[cnt_fail],
                                 orbitData,
                                 orbitDataNum,
                                 p_rank,
                                 p_size );
#elif defined(_ENABLE_NEW_SEARCH3_)
        FindNewPointByAllSearch3(&curTime,
                                 &point_id,
                                 &cur_ra,
                                 &cur_dec,
                                 skyMap,
                                 normal_sky_num,
                                 infp,
                                 earthStartTime,
                                 initalGreenWithRa,
                                 outputresult[cnt_good],
                                 output_fail[cnt_fail],
                                 orbitData,
                                 orbitDataNum,
                                 p_rank,
                                 p_size );

#elif defined(_ENABLE_NEW_SEARCH4_)
        FindNewPointByAllSearch4(&curTime,
                                 &lastObsEndTime,
                                 &try_steps[cnt_good],
                                 &point_id,
                                //  &cur_ra,
                                //  &cur_dec,
                                 skyMap,
                                 normal_sky_num,
                                 infp,
                                 earthStartTime,
                                 initalGreenWithRa,
                                 outputresult[cnt_good],
                                 output_fail[cnt_fail],
                                 orbitData,
                                 orbitDataNum,
                                 p_rank,
                                 p_size );
#endif

        int resultId = (int)outputresult[cnt_good][0];
        
        if ( resultId >= 0 ) {

            point_id = resultId;
            double b = skyMap[resultId].gb;
            excuteNum++;

            int posflag = 0;
            if ( fabs(b) <= Galaxy_B_Sec_Low )
                posflag = 1;

#if defined(_ENABLE_MARK_OBS_POINT_)
            if ( p_rank == 0 ) {
                MarkObervePoint( skyMap[resultId].dec,
                                 skyMap[resultId].ra,
                                 Trees12,
                                 ccd_pos_in_focus,
                                 &area1,
                                 &area2,
                                 deltArea );
            }

            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Bcast(&area1, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
            MPI_Bcast(&area2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
#endif

            area1_[cnt_good]    = area1;
            area2_[cnt_good]    = area2;

            dec[cnt_good]       = skyMap[resultId].dec;
            ra[cnt_good]        = skyMap[resultId].ra;
            deepflag[cnt_good]  = skyMap[point_id].inDeepFlag;
            sunside[cnt_good]   = skyMap[point_id].IsInSunSideFlag;
            cmgtotal[cnt_good]  = cmg_list->cmg_total;

            posflag_[cnt_good]  = posflag;
            battery_q_[cnt_good]= battery_q;
            id_[cnt_good]       = skyMap[resultId].id;

            cnt_good = cnt_good + 1;

            int idx;
            if(p_rank==0 && cnt_good >= NMAX) {
                for( idx=0; idx<NMAX; idx++ ){

                //  time, lat, lon, (1,2,3)
                    fprintf(outf, "%15.8f %8.4f %8.4f ", outputresult[idx][1], dec[idx], ra[idx]);
                //  satx, saty, satz, (4,5,6)
                    fprintf(outf, "%10.4f %10.4f %10.4f ", outputresult[idx][2], outputresult[idx][3], outputresult[idx][4]);
                //  sunx, suny, sunz, (7,8,9)
                    fprintf(outf, "%10.4f %10.4f %10.4f ", outputresult[idx][5], outputresult[idx][6], outputresult[idx][7]);
                //  moonx, moony, moonz, (10,11,12)
                    fprintf(outf, "%10.4f %10.4f %10.4f ", outputresult[idx][8], outputresult[idx][9], outputresult[idx][10]);
                //  survey_area, survey_area2, isdeep, 大面积巡天，极深度巡天 (13,14,15)
                    fprintf(outf, "%12.6f %12.6f %10.4f ", area1_[idx],area2_[idx], deepflag[idx]);
                //  isInGalaxydisk, (16)
                    fprintf(outf, "%2d ", posflag_[idx]);
                //  exposuretime, (17)
                    fprintf(outf, "%10.4f ", outputresult[idx][11]);
                //  transAngle, 转动的角度 (18)
                    fprintf(outf, "%10.5f ", outputresult[idx][12]);
                //  isInSunsid (19)
                    fprintf(outf, "%2d ", sunside[idx]);
                //  CMG ... (20)
                    fprintf(outf, "%10.4f ", cmgtotal[idx]);
                //  battery level, (21)
                    fprintf(outf, "%10.4f ", battery_q_[idx]);
                //  太阳和帆板法线的夹角 (22)
                    fprintf(outf, "%10.4f ", outputresult[idx][14]);
                //  SAA 时间 (23)
                    fprintf(outf, "%10.4f ", outputresult[idx][15]);
                //  天区的ID (24)
                    fprintf(outf, "%8d ", id_[idx]);
#if defined(_ENABLE_NEW_SEARCH4_)
                //  5deg try steps
                    fprintf(outf, "%3d", try_steps[idx]);
#endif
                //  switch to next line
                    fprintf(outf, "\n");
                }
                fflush(outf);
                cnt_good = 0;
            }

            MPI_Bcast(&cnt_good, 1, MPI_INT, 0, MPI_COMM_WORLD );
        } else {
			cnt_fail = cnt_fail + 1;

            if( p_rank == 0 && cnt_fail == NMAX ) {
                int idx;

                for( idx=0; idx<NMAX; idx++ ){
                    int o_i = 0;
                    for (o_i = 0; o_i < 10; o_i ++) {
                        fprintf(outf_fail,"%15.8f ",output_fail[idx][o_i]);
                    }
                    fprintf(outf_fail,"\n");
                }
                fflush(outf_fail);
                cnt_fail = 0;
            }
            
            MPI_Bcast(&cnt_fail, 1, MPI_INT, 0, MPI_COMM_WORLD );
            noExcuteNum++;
        }

        //  更新可观测天区的信息

        if( sky_id_tracker_cnt >= 10000 ){

#if defined(_UPDATE_ONLY_SKY_NUM_)
            sky_num_remained_tmp = count_sky_num(skyMap,normal_sky_num,p_rank);
            if( p_rank == 0 ){
                double yr = (curTime-2459766.0)/365.25;
                printf("> t = %8.4f yr,",yr);
                printf(" %8d patches to be observed\n", sky_num_remained_tmp);
            }
#else

            if( sky_id_tracker!=NULL )
                free(sky_id_tracker);
    
            sky_id_tracker = update_sky_id_tracker( skyMap, 
                                                    normal_sky_num, 
                                                    p_rank,
                                                    &sky_num_remained );
            if( p_rank == 0 ){
                double yr = (curTime-2459766.0)/365.25;
                printf("> t = %8.4f yr,",yr);
                printf(" %8d patches to be observed\n", sky_num_remained);
            }
#endif

            sky_id_tracker_cnt = 0;
        }

        // if( p_rank == 0 ){
        //     printf("--> curTime = %15.8f\n",curTime);
        // }

        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    if(p_rank==0 && cnt_good > 0 ) {
        int idx;

        for( idx=0; idx<cnt_good; idx++ ){
        //  time, lat, lon,
            fprintf(outf, "%15.8f %8.4f %8.4f ", outputresult[idx][1], dec[idx], ra[idx]);
        //  satx, saty, satz,
            fprintf(outf, "%10.4f %10.4f %10.4f ", outputresult[idx][2], outputresult[idx][3], outputresult[idx][4]);
        //  sunx, suny, sunz,
            fprintf(outf, "%10.4f %10.4f %10.4f ", outputresult[idx][5], outputresult[idx][6], outputresult[idx][7]);
        //  moonx, moony, moonz,
            fprintf(outf, "%10.4f %10.4f %10.4f ", outputresult[idx][8], outputresult[idx][9], outputresult[idx][10]);
        //  survey_area, survey_area2, isdeep, 大面积巡天，极深度巡天
            fprintf(outf, "%12.6f %12.6f %10.4f ", area1_[idx],area2_[idx], deepflag[idx]);
        //  isInGalaxydisk,
            fprintf(outf, "%2d ", posflag_[idx]);
        //  exposuretime, 
            fprintf(outf, "%10.4f ", outputresult[idx][11]);
        //  transAngle, 转动的角度
            fprintf(outf, "%10.5f ", outputresult[idx][12]);
        //  isInSunsid
            fprintf(outf, "%2d ", sunside[idx]);
        //  CMG ...
            fprintf(outf, "%10.4f ", cmgtotal[idx]);
        //  battery level,
            fprintf(outf, "%10.4f ", battery_q_[idx]);
        //  太阳和帆板法线的夹角
            fprintf(outf, "%10.4f ", outputresult[idx][14]);
        //  SAA 时间
            fprintf(outf, "%10.4f ", outputresult[idx][15]);
        //  天区的ID
            fprintf(outf, "%8d ", id_[idx]);
#if defined(_ENABLE_NEW_SEARCH4_)
        //  5deg try steps
            fprintf(outf, "%3d", try_steps[idx]);
#endif
        //  switch to next line
            fprintf(outf, "\n");
        }
    }

    if( p_rank==0 && cnt_fail > 0 ) {
        int idx;

        for( idx=0; idx<cnt_fail; idx++ ){
            int o_i = 0;
            for (o_i = 0; o_i < 10; o_i ++) {
                fprintf(outf_fail,"%15.8f ",output_fail[idx][o_i]);
            }
            fprintf(outf_fail,"\n");
        }
    }


    for( m=0; m<NMAX; m++ ){
		free(outputresult[m]);
		free(output_fail[m]);
	}


    time(&etime);
    if ( p_rank == 0 ) {
        fprintf(outf, "# num of executions  %d,  # num of un-executions %d\n", excuteNum, noExcuteNum);
        fprintf(outf, "# num of sun angle with orbit_angle > 15 is %f\n", toearthTime);
        fprintf(outf, "# total simulation time is %ld seconds\n", etime - stime);
        printf("执行的次数为：  %d    未执行的次数为:    %d   \n", excuteNum, noExcuteNum);
        printf("use time is %ld s\n", etime - stime);

#if defined(_UPDATE_ONLY_SKY_NUM_)
        sky_num_remained = count_sky_num(skyMap,normal_sky_num,p_rank);
#endif
        printf("# of un-observed sky patches: %d\n",sky_num_remained);
    }
    

#if defined(_ENABLE_MARK_OBS_POINT_)
    if( p_rank == 0)
        freeTrees(Trees12);
#endif

    if( sky_id_tracker != NULL )
        free(sky_id_tracker);
    free(skyMap);
    free(orbitData);
    MPI_Finalize();
    return 0;
}
