//  这个源文件中的主函数由SurveySim_MPI.c修改而来.
//  读入生成的指向序列,检测相邻两次指向是否为"正常的连续观测",也就是
//  "第一次观测结束后的时间" + "转动到下一次观测开始时所需要的时间" = "下一次观测开始的时间"
//  如果这个等式不成立,那么就需要首先排除一些特殊的情况:
//      1) 处于停靠时期; 
//      2) 处于beta角小于15°或者10.25°; 
//  排除了以上这些情况之后,从上一次曝光结束的时间开始,以60秒的间隔开始尝试寻找可观测天区,并记录下
//  所有的限制因素(注:在被某一个因素限制之后,不会通过continue直接跳过随后的判断)

#include "SurveySim.h"


bool in_stopTime(double t){
    int i=0;
    for(i=0; i<Stop_Time_Num; i++){
        if( t>= stopTime[i][0] && t<= stopTime[i][1]){
            return true;
        }
    }

    return false;
}

int read_simulation_result( char* result_file,
                            int nobs,
                            double* obs_time,
                            double* dec,
                            double* ra,
                            double* satx,
                            double* saty,
                            double* satz,
                            double* sunx,
                            double* suny,
                            double* sunz,
                            double* moonx,
                            double* moony,
                            double* moonz,
                            double* area_norm,
                            double* area_deep,
                            double* is_deep,
                            double* in_gdisk,
                            double* exp_time,
                            double* tAnge,
                            double* isinsunside,
                            double* cmgtotal,
                            double* battery_level,
                            double* solar_panel_angle,
                            double* saa_time,
                            int* sky_id ){

    printf("reading simulation result from %s\n",result_file);
    int i=0;
    FILE* fp = fopen(result_file,"r");
    while( fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d",
                    obs_time+i,
                    dec+i,
                    ra+i,
                    satx+i,
                    saty+i,
                    satz+i,
                    sunx+i,
                    suny+i,
                    sunz+i,
                    moonx+i,
                    moony+i,
                    moonz+i,
                    area_norm+i,
                    area_deep+i,
                    is_deep+i,
                    in_gdisk+i,
                    exp_time+i,
                    tAnge+i,
                    isinsunside+i,
                    cmgtotal+i,
                    battery_level+i,
                    solar_panel_angle+i,
                    saa_time+i,
                    sky_id+i) > 0 ){

        i++;

        // printf("t = %15.8f, ra = %8.4f, dec = %8.4f, sunx = %15.8f, suny = %15.8f, sunz = %15.8f, sky_id = %d\n",
        //         *(obs_time+i),*(ra+i),*(dec+i),
        //         *(sunx+i),*(suny+i),*(sunz+i),
        //         *(sky_id+i));

        // if( i > 5 ){
        //     exit(0);
        // }

        if( i>=nobs )
            break;
    }

    fclose(fp);
    printf("finished reading simulation result...\n");
    
    if( i != nobs ){
        printf("Error in reading %s\n",result_file);
        exit(0);
    }

    return 0;
}

int main(int argc, char *argv[]) {

    int p_rank;
    int p_size = 0;
    
    MPI_Init(NULL, NULL );
    MPI_Comm_rank(MPI_COMM_WORLD, &p_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p_size);

    init_default();

    int i = 0; //, j = 0;
    int isReadStatus = 0;
    double curTime = 0;

    infp = NULL;
    double earthStartTime, initalGreenWithRa = 0;

    if (argc < 4) {
        fprintf(stderr, "usage:   %s config-file result_root nobs\n", argv[0]);
        MPI_Finalize();
        return 0;
    }

    char infile[1024];
    char outfile[1024];

    sprintf(infile,"%s.dat",argv[2]);
    sprintf(outfile,"%s.dT",argv[2]);

    int nobs = atoi(argv[3]);

    parseConfigFile(argv[1]);

    if(p_rank == 0) {
        printConfig();
        printf("interval time is %f \n",jump_time);
    }

    if ((infp = fopen(jpl_fileName, "r")) == NULL ) {
        fprintf(stderr, "\nERROR: Can't open ephemeris file %s for input.\n\n", argv[1]);
        MPI_Finalize();
        return 0;
    }

    if (p_rank == 0) {
        //计算某一时刻地区本初子午线与赤道交点午夜时在天球赤道坐标系中的位置
        earthStartTime = floor(startTime - 0.5) + 0.5; // midnight ut time (unit day)
        double sunMid[3];
        locateSun(infp, earthStartTime, sunMid); // 太阳的笛卡尔坐标，以地球为原点(J2000)

        double reEarth[2]; // (ra,dec) of Sun
        Cartesian2Equatorial(sunMid, reEarth);  // get (ra,dec) of Sun

        if (reEarth[0] >= 180) {
            initalGreenWithRa = reEarth[0] - 180;
        } else {
            initalGreenWithRa = reEarth[0] + 180;
        }
    }

    MPI_Bcast(&initalGreenWithRa, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast(&earthStartTime, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast(&isReadStatus, 1, MPI_INT, 0, MPI_COMM_WORLD );

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

    time_t stime, etime;
    time(&stime);

//==============================================================================

    cmg_list = (struct CMG_List*)malloc(sizeof(struct CMG_List));
    initCMGList(cmg_list);

    battery_q = BATTERY_MAX;    // sun plane, energy balance（全局变量）
    shadow_st = startTime;      // 进入阴影区的时刻（全局变量）
    shadow_et = startTime;      // 离开阴影区的时刻（全局变量）

    
    double* outputresult = (double*)malloc(sizeof(double)*16);

	// output_fail 用于存储一些输出变量,和限制因素导致无法观测的次数有关
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
    
    double* output_fail = (double*)malloc(sizeof(double)*10);

//  =============================================================================
    double *obs_time = calloc(nobs,sizeof(double));
    double *dec = calloc(nobs,sizeof(double));
    double *ra = calloc(nobs,sizeof(double));
    double *satx = calloc(nobs,sizeof(double));
    double *saty = calloc(nobs,sizeof(double));
    double *satz = calloc(nobs,sizeof(double));
    double *sunx = calloc(nobs,sizeof(double));
    double *suny = calloc(nobs,sizeof(double));
    double *sunz = calloc(nobs,sizeof(double));
    double *moonx = calloc(nobs,sizeof(double));
    double *moony = calloc(nobs,sizeof(double));
    double *moonz = calloc(nobs,sizeof(double));
    double *area_norm = calloc(nobs,sizeof(double));
    double *area_deep = calloc(nobs,sizeof(double));
    double *is_deep = calloc(nobs,sizeof(double));
    double *in_gdisk= calloc(nobs,sizeof(double));
    double *exp_time= calloc(nobs,sizeof(double));
    double *tAngle = calloc(nobs,sizeof(double));
    double *isinsunside = calloc(nobs,sizeof(double));
    double *cmgtotal = calloc(nobs,sizeof(double));
    double *battery_level = calloc(nobs,sizeof(double));
    double *solar_panel_angle = calloc(nobs,sizeof(double));
    double *saa_time = calloc(nobs,sizeof(double));
    int *sky_id = calloc(nobs,sizeof(int));

    read_simulation_result( infile,
                            nobs,
                            obs_time,
                            dec,
                            ra,
                            satx,
                            saty,
                            satz,
                            sunx,
                            suny,
                            sunz,
                            moonx,
                            moony,
                            moonz,
                            area_norm,
                            area_deep,
                            is_deep,
                            in_gdisk,
                            exp_time,
                            tAngle,
                            isinsunside,
                            cmgtotal,
                            battery_level,
                            solar_panel_angle,
                            saa_time,
                            sky_id );

//  ==================================================================================================
    int cnt=1;

    int err_cnt=0;
    double wasted_time = 0;

    double dt=0.5;

    double dtime_max = 0;
    int idx_max=0;
    int beta_cnt=0;

    // 用于记录哪些指向与前一次指向之间的时间差超出了预计值.
    int *is_err = calloc(nobs,sizeof(int));

    FILE *fp=NULL;
    if ( p_rank == 0 )
        fp = fopen(outfile,"w");

    while ( cnt < nobs ) {

        double dtime = ( obs_time[cnt] - obs_time[cnt-1] )*86400;   //转换为秒
        double trans_time = 0;

        if( tAngle[cnt] == 0 ) // 可能会遇到转动角为0的情况，之前没有将这种情况考虑进来。
            trans_time = 0;
        else
            trans_time = calculateTransTime(tAngle[cnt]);

        int beta_flag=0;
        int saa_flag=0;
        int sat_flag=0;

        if( dtime/86400 < 2 ){

            double t0=obs_time[cnt-1]+(exp_time[cnt-1]+trans_time)/86400;
            double t1=obs_time[cnt];
            double tt=t0;

            // printf("tt = %g\tt1 = %g\n",tt,t0);

            double sunAngle = getSun2OrbitAngle1(infp, tt, orbitData, orbitDataNum,Orbit_File_Num);
            while( tt < t1 ) { 

                if( fabs(sunAngle) < BETA_ANGLE ){
                    // printf("cnt = %8d, T = %15.8f, found beta-angle < %g\n",cnt,tt,BETA_ANGLE);
                    beta_flag += 1;
                }

/*                double tSat[3], ttSat[3], tSatE[3], ttSatE[3];*/
				double ttSatE[3];

                if( locateSat1(ttSatE, tt, orbitData, orbitDataNum,Orbit_File_Num) < 0.0 ){
                    printf("failed to obtain position of Sat ...\n");
                    sat_flag++;
                }

                // CoordinateSpinEquatorial2Ecliptic(ttSatE, tSatE);
            //  计算星下点的位置
                double unSPoint[2];
                calculateSatEarthPoint( earthStartTime, 
                                        initalGreenWithRa, 
                                        tt, 
                                        ttSatE, 
                                        unSPoint, 
                                        infp );

                // printf("usSPoint[0] = %10.5f, usSPoint[1] = %10.5f\n",unSPoint[0],unSPoint[1]);

                if ( IsInSAA(unSPoint) == 1 ) {
                    // printf("find positions in SAA ...\n");
                    saa_flag += 1;
                }

                tt += 20./86400;
                sunAngle = getSun2OrbitAngle1(infp, tt, orbitData, orbitDataNum,Orbit_File_Num);

            }
        } else {
            // double sunAngle1 = getSun2OrbitAngle1(infp, obs_time[cnt-1]+12*3600/86400, orbitData, orbitDataNum,Orbit_File_Num);
            // double sunAngle2 = getSun2OrbitAngle1(infp, obs_time[cnt]-12*3600/86400, orbitData, orbitDataNum,Orbit_File_Num);
            // if( fabs(sunAngle1) < 10.25 || fabs(sunAngle2) < 10.25 ){
            // 超过两天的时间间隔,应该都是由于"停靠"或者是"轨道面与太阳的夹角"小于"BETA_ANGLE"!
            // printf("sunAngle1 = %8.3f, sunAngle2 = %8.3f\n",sunAngle1,sunAngle2);
            cnt++;
            continue;
        }


        if( dtime > (exp_time[cnt-1]+trans_time+dt)
            && beta_flag == 0
            && saa_flag == 0
            && sat_flag == 0 ){
            
            err_cnt++;
            
            wasted_time += ( dtime - (trans_time+exp_time[cnt-1]) );

            double ddtime = dtime - (exp_time[cnt-1]+trans_time);
            if ( p_rank == 0 ){
                fprintf(fp,"%15.8f  %15.8f %15.8f %15.8f %15.8f %d\n",
                            (obs_time[cnt-1]-2459766)/365.25,
                            ddtime,
                            cmgtotal[cnt-1],
                            cmgtotal[cnt],
                            tAngle[cnt],
                            cnt);
            }

            if( dtime > dtime_max ){
                dtime_max = dtime;
                idx_max = cnt;
            }

            is_err[cnt] = 1;
        }

        if( beta_flag != 0 ){
            beta_cnt++;
        }

        cnt++;
    }

    if ( p_rank == 0 )
        fclose(fp);

    if ( p_rank == 0 ) {
        time(&etime);
        printf("--> finished search for unused observational time intervals\n");
        printf("--> use time is %ld s\n", etime - stime);
    }

    printf("-->  %d unused observational time interval are found!\n", err_cnt);
    printf("--> total wasted time is %f days\n",wasted_time/86400);
    printf("--> max time interval = %g days\n",dtime_max/86400);
    printf("--> corresponding idx : %d\n",idx_max);
    printf("--> BETA_CNT = %d\n",beta_cnt);

//  ==================================================================================================
//  开始对前面选择出来的 is_err 进行循环.

    // if ( p_rank == 0 )
    //     fp = fopen("Error_conditions.txt","w");

    // cnt=1;
    // while( cnt < 100 ){
    //     if( is_err[cnt] == 1 ){
    //     //  t0: 上一次曝光介绍后的时间
    //     //  t1: 本次曝光开始的时间
    //         double trans_time = calculateTransTime(tAngle[cnt]);
    //         double t0 = obs_time[cnt-1]+(exp_time[cnt-1]+trans_time)/86400;
    //         double t1 = obs_time[cnt];
            
    //         double tt = t0;
    //         double loc_time=0;
    //         double loc_sun=0;
    //         double loc_moon=0;
    //         double loc_earth=0;
    //         double loc_energy=0;
    //         double loc_cmg = 0;
    //         double loc_panel=0;
    //         double loc_observed=0;
    //         double loc_saa=0;
    //         double in_sunside=0;

    //         point_id = sky_id[cnt-1];
    //         battery_q = battery_level[cnt-1];

    //         int cnt_tmp=0;
    //         while( tt < t1 ){

    //             if ( tt > shadow_et ) {
    //             //  如果当前处于阴影区以外了，那就获取下一次进入和离开阴影区的时刻。
    //                 aquireShadowTime( &shadow_st, &shadow_et, tt );
    //             }

    //             setDoubleArrayZero( output_fail, 10 );

    //             double tt_tmp=tt;

    //             // printf(".... tt = %15.8f  ",tt_tmp);

    //             TestCondition(  &tt_tmp,
    //                             &point_id,
    //                             skyMap,
    //                             normal_sky_num,
    //                             infp,
    //                             earthStartTime,
    //                             initalGreenWithRa,
    //                             outputresult,   // 这个其实不需要
    //                             output_fail,    // 这里记录了每次无法观测时的各个限制因素所占的次数
    //                             orbitData,
    //                             orbitDataNum,
    //                             p_rank,
    //                             p_size );

    //             printf("tt_tmp - tt = %15.8f\n",tt_tmp-tt);

    //         //  将output_fail中的结果临时保存,最后进行平均
    //             if( outputresult[0] < 0 ){
    //                 loc_time += output_fail[0];
    //                 loc_sun += output_fail[1];
    //                 loc_moon += output_fail[2];
    //                 loc_earth += output_fail[3];
    //                 loc_energy += output_fail[4];
    //                 loc_cmg += output_fail[5];
    //                 loc_panel += output_fail[6];
    //                 loc_observed += output_fail[7];
    //                 loc_saa += output_fail[8];
    //                 in_sunside += output_fail[9];

    //                 printf("%15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f\n",
    //                     loc_time,
    //                     loc_sun,
    //                     loc_moon,
    //                     loc_earth,
    //                     loc_energy,
    //                     loc_cmg,
    //                     loc_panel,
    //                     loc_observed,
    //                     loc_saa,
    //                     in_sunside);
    //             }

    //             tt += 60./86400;   //以60秒为间隔来进行各项限制条件的判断
    //             cnt_tmp++;
    //         }
            
    //         if( cnt_tmp >= 1 && p_rank == 0 ){
    //             // printf("--->  cnt_tmp = %d\n",cnt_tmp);
    //             fprintf(fp,"%15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f\n",
    //                     loc_time/cnt_tmp,
    //                     loc_sun/cnt_tmp,
    //                     loc_moon/cnt_tmp,
    //                     loc_earth/cnt_tmp,
    //                     loc_energy/cnt_tmp,
    //                     loc_cmg/cnt_tmp,
    //                     loc_panel/cnt_tmp,
    //                     loc_observed/cnt_tmp,
    //                     loc_saa/cnt_tmp,
    //                     in_sunside/cnt_tmp);

    //             printf("%15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f\n",
    //                     loc_time/cnt_tmp,
    //                     loc_sun/cnt_tmp,
    //                     loc_moon/cnt_tmp,
    //                     loc_earth/cnt_tmp,
    //                     loc_energy/cnt_tmp,
    //                     loc_cmg/cnt_tmp,
    //                     loc_panel/cnt_tmp,
    //                     loc_observed/cnt_tmp,
    //                     loc_saa/cnt_tmp,
    //                     in_sunside/cnt_tmp);
    //         }
    //         // else{
    //         //     printf("---> failed to re-test the conditions ...\n");
    //         // }
    //     }

    //     cnt++;
    // }

    // if( p_rank == 0 )
    //     fclose(fp);

//  ==================================================================================================

    if ( p_rank == 0 ) {
        time(&etime);
        printf("--> finished re-testing the conditions\n");
        printf("--> use time is %ld s\n", etime - stime);
    }
    
    free(outputresult);
    free(output_fail);
    free(skyMap);
    free(orbitData);
        
    if( sky_id_tracker != NULL )
        free(sky_id_tracker);

    free(obs_time);
    free(dec);
    free(ra);
    free(satx);
    free(saty);
    free(satz);
    free(sunx);
    free(suny);
    free(sunz);
    free(moonx);
    free(moony);
    free(moonz);
    free(area_norm);
    free(area_deep);
    free(is_deep);
    free(in_gdisk);
    free(exp_time);
    free(tAngle);
    free(isinsunside);
    free(cmgtotal);
    free(battery_level);
    free(solar_panel_angle);
    free(saa_time);
    free(sky_id);

    MPI_Finalize();
    return 0;
}
