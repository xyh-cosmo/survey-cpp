#include "SurveySim.h"
#include "CoverAreaStat.h"

// double stopTime[Stop_Time_Num][2]= {
//     {2460131.56,2460205.53},
//     {2460887.37,2460997.22},
//     {2461688.66,2461796.99},
//     {2462429.30,2462538.41},
//     {2463166.49,2463275.57},
// };

int main(int argc, char *argv[]) {

    int p_rank;
    int p_size = 0;
    MPI_Init(NULL, NULL );
    MPI_Comm_rank(MPI_COMM_WORLD, &p_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p_size);

    init_default();

    int isReadStatus = 0;
    double curTime = 0;

    double earthStartTime, initalGreenWithRa = 0;

    if (argc < 2) {
        fprintf(stderr, "\nFormat:\n\n         %s config-file\n\n", argv[0]);
        MPI_Finalize();
        return 0;
    }

    parseConfigFile(argv[1]);

    jump_time = 60.0;   // in unit of seconds

    if(p_rank == 0) {
        printConfig();
        printf("interval time is %f \n",jump_time);
    }

    if ((infp = fopen(jpl_fileName, "r")) == NULL ) {
        fprintf(stderr, "\nERROR: Can't open ephemeris file %s for input.\n\n", argv[1]);
        MPI_Finalize();
        return 0;
    }

    int outfile_status1 = 0;
    int outfile_status2 = 0;

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

    MPI_Bcast(&outfile_status1, 1, MPI_INT, 0, MPI_COMM_WORLD );
    if( outfile_status1 != 0 ){
        MPI_Finalize();
        return 0;
    }

    MPI_Bcast(&outfile_status2, 1, MPI_INT, 0, MPI_COMM_WORLD );
    if( outfile_status2 != 0 ){
        MPI_Finalize();
        return 0;
    }

    MPI_Bcast(&initalGreenWithRa, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast(&earthStartTime, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast(&isReadStatus, 1, MPI_INT, 0, MPI_COMM_WORLD );

//  读取所有的轨道数据
    orbitData = (double**)malloc(sizeof(double*) * Orbit_File_Num);
    readAllorbitsFile(orbitData, orbitDataNum);

    int i,j;
    double beta_time_i, beta_time_tot = 0;
    double sunAngle;


    for( i=0; i<50; i++ ){
        printf("proceesing %3d th orbit data file:\n",i+1);
        beta_time_i = 0;
        
        int in_beta = 0;
        double time_tmp = 0;
        for( j=0; j<orbitDataNum[i]-1; j++ ){
            curTime = 0.5*(orbitData[i][j*7]+orbitData[i][(j+1)*7]);
            // double sunAngle = getSun2OrbitAngle1(infp, curTime,orbitData, orbitDataNum, Orbit_File_Num);

            double a[3],b[3],normalVect[3],sun[3];
            a[0] = orbitData[i][j*7+1];
            a[1] = orbitData[i][j*7+2];
            a[2] = orbitData[i][j*7+3];
            b[0] = orbitData[i][(j+1)*7+1];
            b[1] = orbitData[i][(j+1)*7+2];
            b[2] = orbitData[i][(j+1)*7+3];

            crossMultiple(a,b,normalVect);
            locateSun(infp, curTime, sun);

            double pointMul = sun[0]*normalVect[0] + sun[1]*normalVect[1]+ sun[2]*normalVect[2];
            double modSun = sqrt(sun[0]*sun[0] + sun[1]*sun[1] + sun[2]*sun[2]);
            double modNormal = sqrt(normalVect[0]*normalVect[0] + normalVect[1]*normalVect[1] + normalVect[2]*normalVect[2]);
            sunAngle = acos(pointMul/(modSun*modNormal))*57.29577951;
            sunAngle = 90 - sunAngle;

            if( fabs(sunAngle) < BETA_ANGLE ){
                beta_time_i += 120./86400;
                in_beta = 1;
                time_tmp += 120;
            }
            else{
                in_beta = 0;
            }
            
            if( in_beta == 0 && time_tmp > 0 ){
                printf("beta_time_interval = %10.5f seconds\n",time_tmp);
                time_tmp = 0;
            }

        }

        beta_time_tot += beta_time_i;

        printf("--> beta_time[%d] = %10.5f days, beta_time_tot = %10.5f days\n", i+1,beta_time_i, beta_time_tot);
    }
   
    free(orbitData);

    MPI_Finalize();
    return 0;
}
