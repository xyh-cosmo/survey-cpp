
#include "SurveySim.h"

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

    if (argc < 2) {
        fprintf(stderr, "usage:   %s config-file\n", argv[0]);
        MPI_Finalize();
        return 0;
    }

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

//==============================================================================

    double t_start = 2459767;
    double t = t_start;
    double dist;

    time(&stime);
    while( t < t_start + 365.25*10 ){
        dist = locateSat1(coor_e,t,orbitData,orbitDataNum,Orbit_File_Num);
        t += 1./86400;
    }
    time(&etime);

    if ( p_rank == 0 ) {
        printf("--> finished search for unused observational time intervals\n");
        printf("--> use time is %ld s\n", etime - stime);
    }

    MPI_Finalize();
    return 0;
}
