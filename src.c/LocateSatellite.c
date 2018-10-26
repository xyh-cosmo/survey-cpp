#include "SurveySim.h"

#define MAXLINE 10000

#define elemNum 7

void crossMultiple(double* a, double* b, double* result) {
    result[0] = a[1]*b[2] - a[2]*b[1];//a2b3-a3b2
    result[1] = a[2]*b[0] - a[0]*b[2];//a3b1-a1b3
    result[2] = a[0]*b[1] - a[1]*b[0];//a1b2-a2b1
}

// Read orbit data
double* readOrbitData(char* fileName, int* static_num) {
    FILE *orbitFile = fopen(fileName,"r");
    *static_num = 0;
    if(orbitFile == NULL) {
        printf("No such file:'%s'\n",fileName);
        exit(1);
    }

    int num = 0;
    char line[MAXLINE];

    char* ta = "#";
    while (fgets(line, MAXLINE, orbitFile) != NULL) {
        if(line[0] == ta[0]) {
            continue;
        }
        *static_num = *static_num + 1;
    }
    double* data = (double*) malloc(sizeof(double)*(*static_num)*elemNum);
    fclose(orbitFile);

    orbitFile = fopen(fileName,"r");
    if(orbitFile == NULL) {
        printf("No such file:'%s'\n",fileName);
        exit(1);
    }
    while (fgets(line, MAXLINE, orbitFile) != NULL) {
        if(line[0] == ta[0])
            continue;

        char *p;
        p = strtok(line, " \n");
        int cNum = 0;
        while (p != NULL) {
            double cc = atof(p);
            *(data + num * elemNum + cNum) = cc;
            cNum++;
            p = strtok(NULL, " \n");
        }
        num ++;
    }
    fclose(orbitFile);
    return data;
}

//Read all orbit data
void readAllorbitsFile(double** OrbitData, int* dataNum) {
    int i;
    for(i = 0; i < Orbit_File_Num; i ++) {
        char orbitfileName[200] = "";
        strcat(orbitfileName,"orbit20160925/");
        char idn[10];
        sprintf(idn,"%d",i+1);
        strcat(orbitfileName,idn);
        strcat(orbitfileName,".txt");
        //printf("%s\n",orbitfileName);
        OrbitData[i] = readOrbitData(orbitfileName, &dataNum[i]);
    }

    for(i = 0; i < Orbit_File_Num; i ++) {
        timePoint[i][0] = *(OrbitData[i]);
        timePoint[i][1] = *(OrbitData[i] + (dataNum[i]-1) * elemNum);
    }
}

void freeAllOrbitData(double** OrbitData, int FileNum) {
    int i;
    for(i = 0; i< FileNum; i ++) {
        free(OrbitData[i]);
    }
    free(OrbitData);
}

double hermite3(double x,double x0, double y0, double x1, double y1, double d0,double d1) {
    double y;
    double h0 = (1+2*((x-x0)/(x1-x0)))*((x-x1)/(x0-x1))*((x-x1)/(x0-x1));
    double h1 = (1+2*((x-x1)/(x0-x1)))*((x-x0)/(x1-x0))*((x-x0)/(x1-x0));
    double wh0 = (x - x0)*((x-x1)/(x0-x1))*((x-x1)/(x0-x1));
    double wh1 = (x - x1)*((x-x0)/(x1-x0))*((x-x0)/(x1-x0));
    y = y0*h0 +y1*h1+d0*wh0+d1*wh1;
    return y;
}

void matrixMultiple(double a[3][3], double b[3][3], double result[3][3]) {

    if (a == NULL || b == NULL ) {
        return;
    }

    int i, j, k;
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            result[i][j] = 0;
            for (k = 0; k < 3; k++) {
                result[i][j] += a[i][k] * b[k][j];
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////
int init_gsl_interp_accel(){
    int i;
    for( i=0; i<Orbit_File_Num; i++ ){
        locate_acc[i] = gsl_interp_accel_alloc();
        gsl_interp_accel_reset(locate_acc[i]);
    }
    return 0;
}

int free_gsl_interp_accel(){
    int i;
    for( i=0; i<Orbit_File_Num; i++ ){
        gsl_interp_accel_free(locate_acc[i]);
    }
    return 0;
}

int update_orbit_info( 	double 	time,
						int* 	ofile_idx,
						int*	last_ofile_idx,
						int*	next_ofile_idx,
						int*	at_end ){
    int i, status=0;
    *ofile_idx = -1;
    for( i=0; i<Orbit_File_Num; i++ ){
        int data_size = orbitDataNum[i];
        if( time >= timePoint[i][0] && time <= timePoint[i][1] ){
            *ofile_idx = i;
            *last_ofile_idx = i - 1;
            *next_ofile_idx = i + 1;
            if( time > orbitTime[i][data_size-2] ){
                *at_end = 1;
            }
            break;
        }
    }

    if( *ofile_idx == -1 ){
        status = 1;
    }

    return status;
}

//////////////////////////////////////////////////////////////////////

// get satellite position from the orbit data
// void locateSat1(double *result, double time, double** OrbitData, int* dataNum, int FileNum) {
double locateSat1(  double*     result, 
                    double      time, 
                    double**    OrbitData, 
                    int*        dataNum, 
                    int         FileNum ) {
    double dist = -9999.9;
    result[0] = 0.0;
    result[1] = 0.0;
    result[2] = 0.0;

//  这里的搜索太费事了，重复了许多不必要的搜索
    int i;
    // int isIn = 0;

    for (i = 0; i < FileNum; i ++) {
        if( time >= timePoint[i][0] && time <= timePoint[i][1]) {
            double *dataSeg = OrbitData[i];
            int segNum = dataNum[i];
            int s_id = 0;
            int e_id = segNum - 1;
            int mid_id = (s_id + e_id)/2;

            while( e_id - s_id > 1 ) {
                if( time <= dataSeg[mid_id*elemNum] ){
                    e_id = mid_id;
                } else {
                    s_id = mid_id;
                }
                mid_id = (s_id + e_id)/2;
            }

            double t1 = dataSeg[s_id*elemNum];
            double t2 = dataSeg[e_id*elemNum];

            // if( (t2-t1) > 130.0/86400.0 ) { // should be 120s, for error, set 120+10s
            if( (t2-t1) > 0.0015046296296296296 ){ // should be 120s, for error, set 120+10s
                break;
            }

            // isIn = 1;

            double x1 = dataSeg[s_id*elemNum + 1];
            double y1 = dataSeg[s_id*elemNum + 2];
            double z1 = dataSeg[s_id*elemNum + 3];

            double x2 = dataSeg[e_id*elemNum + 1];
            double y2 = dataSeg[e_id*elemNum + 2];
            double z2 = dataSeg[e_id*elemNum + 3];

            double l1 = sqrt(x1*x1 + y1*y1 + z1*z1);
            double l2 = sqrt(x2*x2 + y2*y2 + z2*z2);
            double theta = acos((x1*x2 + y1*y2 + z1*z2)/(l1*l2));
            double theta1 = (time-t1)/(t2-t1)*theta;
            double theta2 = theta-theta1;
            double l = (t2-time)/(t2-t1)*l1 + (time-t1)/(t2-t1)*l2;

            double sin_theta1 = sin(theta1);
            double sin_theta2 = sin(theta2);

            double x0 = sin_theta2*x1/l1 + sin_theta1*x2/l2;
            double y0 = sin_theta2*y1/l1 + sin_theta1*y2/l2;
            double z0 = sin_theta2*z1/l1 + sin_theta1*z2/l2;
            double l_ = sqrt(x0*x0 + y0*y0 + z0*z0);

            result[0] = x0*l/l_;
            result[1] = y0*l/l_;
            result[2] = z0*l/l_;

            dist = sqrt( result[0]*result[0] + result[1]*result[1] + result[2]*result[2] );
            break;
        } else if( time<timePoint[i][0] ) {
            if(i>0) {
                double *dataSeg1 = OrbitData[i-1];
                double *dataSeg2 = OrbitData[i];
                int s_id = dataNum[i-1]-1;
                int e_id = 0;

                double t1 = dataSeg1[s_id*elemNum];
                double t2 = dataSeg2[e_id*elemNum];

                // if((t2-t1)>130.0/86400.0) { // should be 120s, for error, set 120+10s
                if((t2-t1)>0.0015046296296296296) { // should be 120s, for error, set 120+10s
                    break;
                }

                // isIn = 1;
                double x1 = dataSeg1[s_id*elemNum + 1];
                double y1 = dataSeg1[s_id*elemNum + 2];
                double z1 = dataSeg1[s_id*elemNum + 3];
                double x2 = dataSeg2[e_id*elemNum + 1];
                double y2 = dataSeg2[e_id*elemNum + 2];
                double z2 = dataSeg2[e_id*elemNum + 3];

                double l1 = sqrt(x1*x1 + y1*y1 + z1*z1);
                double l2 = sqrt(x2*x2 + y2*y2 + z2*z2);
                double theta = acos((x1*x2 + y1*y2 + z1*z2)/(l1*l2));
                double theta1 = (time-t1)/(t2-t1)*theta;
                double theta2 = theta-theta1;
                double l = (t2-time)/(t2-t1)*l1 + (time-t1)/(t2-t1)*l2;
                double ef = sin(theta2)/sin(theta1);

                double x0 = ef*x1/l1 + x2/l2;
                double y0 = ef*y1/l1 + y2/l2;
                double z0 = ef*z1/l1 + z2/l2;
                double l_ = sqrt(x0*x0 + y0*y0 + z0*z0);

                result[0] = x0*l/l_;
                result[1] = y0*l/l_;
                result[2] = z0*l/l_;

                dist = sqrt( result[0]*result[0] + result[1]*result[1] + result[2]*result[2] );
                break;
            } else {
                break;
            }
        }
    }

    // if(isIn == 0) {
    //     result[0] = 0.0;
    //     result[1] = 0.0;
    //     result[2] = 0.0;

    //     dist = -9999.0; // -9999.0 is arbitrarily chosen
    // }

    return dist;
}

//  这个是用来干嘛的？？？？
void locateSat( double *result,
                double raan,  // 升交点赤经
                double altitude,
                double inclination,
                double obPeriod,
                double presession,
                double time ) {
//	printf("raan = %f,   altitude = %f,    inclination = %f, obPeriod=%f,   presession = %f, time=%f\n", raan, altitude, inclination,obPeriod,presession,time);
    // double cos(double), sin(double), fmod(double, double);
    double obritAV = 2. * PI / obPeriod;  // obrit angular velocity
    double presessionAV = 0.;             // presession angular velocity
    if (presession != 0) {
        presessionAV = 2. * PI / presession;
    }
    double f = fmod(obritAV * time, 2. * PI);	  // 轨道面与x轴的夹角
    double pa = fmod(presessionAV * time, 2. * PI);	  // 轨道面绕地轴转动角度，从北极看顺时针转动
//  obrit planet coodinate
    double obX = altitude * cos(f);
    double obY = altitude * sin(f);
    double obZ = 0.;

    // double xaa, yaa, zaa;
    double xaa, zaa;
    double arcInclination = inclination * 2. * PI / 360.;

//  变换到赤道惯性坐标系的变换矩阵
    xaa = arcInclination; // = 0
    zaa = fmod(pa + raan * 2. * PI / 360., 2.*PI); // pa=0

    double tranMatrix4[3][3] = { { 1.0, 0.0, 0.0 },
                                 { 0.0, cos(-xaa), sin(-xaa) },
                                 { 0.0,-sin(-xaa), cos(-xaa) } };

    double tranMatrix5[3][3] = { { cos(-zaa), sin(-zaa), 0.0 },
                                 {-sin(-zaa), cos(-zaa), 0.0 },
                                 { 0.0, 0.0, 1.0 } };

    double eMatrix[3][3];
    matrixMultiple(tranMatrix5, tranMatrix4, eMatrix);
    double intertiaX = eMatrix[0][0] * obX + eMatrix[0][1] * obY + eMatrix[0][2] * obZ;
    double intertiaY = eMatrix[1][0] * obX + eMatrix[1][1] * obY + eMatrix[1][2] * obZ;
    double intertiaZ = eMatrix[2][0] * obX + eMatrix[2][1] * obY + eMatrix[2][2] * obZ;

    *result = intertiaX;
    *(result + 1) = intertiaY;
    *(result + 2) = intertiaZ;
}

// void init_sun_infp(){
    
// }

//get sun position
void locateSun(FILE *infp, double time, double *result) {
    struct ephcom_Header header;
    struct ephcom_Coords coords;
    double *datablock; /* Will hold coefficients from a data block */
    // int datablocknum;
    double testr[6];
    ephcom_readbinary_header(infp, &header);
    // ephcom_writeascii_header(outfp, &header1);
    /*
     Done with header.  Now we'll read and write data blocks.
     */
    datablock = (double *) malloc((int) header.ncoeff * sizeof(double));
    // datablocknum = 0;
    coords.km = 1; /* not AU, use kilometers */
    coords.seconds = 0; /* Timescale is days, not seconds */
    coords.bary = 1; /* Center is Solar System Barycenter */
    coords.coordtype = 0; /* No correction for light travel time or
	 relativistic effects from Sun */

    coords.et2[0] = time; /* Good enough precision for test dates */
    coords.et2[1] = 0.0;
    if (ephcom_get_coords(infp, &header, &coords, datablock) == 0) {
        ephcom_pleph(&coords, EPHCOM_SUN, EPHCOM_EARTH, testr);
        *result = testr[0];
        *(result + 1) = testr[1];
        *(result + 2) = testr[2];
    } else {
        printf("Coordinates not found for %f \n",time);
        *result = 0;
        *(result + 1) = 0;
        *(result + 2) = 0;
    }

    free(datablock);
}


void init_moon_infp(){

}

//get moon position
void locateMoon(FILE *infp, double time, double *result) {
    struct ephcom_Header header;
    struct ephcom_Coords coords;
    double *datablock; /* Will hold coefficients from a data block */
    // int datablocknum;
    double testr[6];
    ephcom_readbinary_header(infp, &header);
    // ephcom_writeascii_header(outfp, &header1);
    /*
     Done with header.  Now we'll read and write data blocks.
     */
    datablock = (double *) malloc(header.ncoeff * sizeof(double));
    // datablocknum = 0;
    coords.km = 1; /* not AU, use kilometers */
    coords.seconds = 0; /* Timescale is days, not seconds */
    coords.bary = 1; /* Center is Solar System Barycenter */
    coords.coordtype = 0; /* No correction for light travel time or
	 relativistic effects from Sun */

    coords.et2[0] = time; /* Good enough precision for test dates */
    coords.et2[1] = 0.0;
    if (ephcom_get_coords(infp, &header, &coords, datablock) == 0) {
        ephcom_pleph(&coords, EPHCOM_MOON, EPHCOM_EARTH, testr);
        *result = testr[0];
        *(result + 1) = testr[1];
        *(result + 2) = testr[2];
    } else {
        printf("Coordinates not found for %f\n",time);
        *result = 0;
        *(result + 1) = 0;
        *(result + 2) = 0;
    }

    free(datablock);
}


//  get the angle between sun vector with orbit plane
//  轨道平面指的是什么？
double getSun2OrbitAngle1(FILE *infp, double time,double** OrbitData, int* dataNum, int FileNum) {
    double angle = -1;
    double satPos[6];
    int i;
    int isIn = 0;

    for (i = 0; i < FileNum; i ++) {
        if(time>=timePoint[i][0] && time <=timePoint[i][1]) {
            isIn = 1;
            
            double *dataSeg = OrbitData[i];
            int segNum = dataNum[i];
            int s_id = 0;
            int e_id = segNum - 1;
            int mid_id = (s_id + e_id)/2;
            while(e_id - s_id > 1) {
                if(time <= dataSeg[mid_id*elemNum]) {
                    e_id = mid_id;
                } else {
                    s_id = mid_id;
                }
                mid_id = (s_id + e_id)/2;
            }
            satPos[0] = dataSeg[s_id*elemNum + 1];
            satPos[1] = dataSeg[s_id*elemNum + 2];
            satPos[2] = dataSeg[s_id*elemNum + 3];
            satPos[3] = dataSeg[e_id*elemNum + 1];
            satPos[4] = dataSeg[e_id*elemNum + 2];
            satPos[5] = dataSeg[e_id*elemNum + 3];

            break;
        } else if(time<timePoint[i][0]) {
            if(i>0) {
                isIn = 1;
                double *dataSeg1 = OrbitData[i-1];
                double *dataSeg2 = OrbitData[i];
                int s_id = dataNum[i-1]-1;
                int e_id = 0;

                satPos[0] = dataSeg1[s_id*elemNum + 1];
                satPos[1] = dataSeg1[s_id*elemNum + 2];
                satPos[2] = dataSeg1[s_id*elemNum + 3];
                satPos[3] = dataSeg2[e_id*elemNum + 1];
                satPos[4] = dataSeg2[e_id*elemNum + 2];
                satPos[5] = dataSeg2[e_id*elemNum + 3];
                break;
            } else {
                break;
            }
        }
    }

    if(isIn == 0) {
        angle = -1;  // 这种情况是否会出现？？？ 如果不会出现，干脆就去掉这一层的判断。

        // 在单独统计beta角的时候,有可能出现无法获取到卫星位置的情况,因此将返回值设为999可以有效排除
        //  一部分错误的统计
        angle = 999;
    } else {
        double a[3],b[3],normalVect[3],sun[3];

        a[0] = satPos[0];
        a[1] = satPos[1];
        a[2] = satPos[2];
        b[0] = satPos[3];
        b[1] = satPos[4];
        b[2] = satPos[5];

        crossMultiple(a,b,normalVect);
        locateSun(infp, time, sun);

        double pointMul = sun[0]*normalVect[0] + sun[1]*normalVect[1]+ sun[2]*normalVect[2];
        double modSun = sqrt(sun[0]*sun[0] + sun[1]*sun[1] + sun[2]*sun[2]);
        double modNormal = sqrt(normalVect[0]*normalVect[0] + normalVect[1]*normalVect[1] + normalVect[2]*normalVect[2]);
        angle = acos(pointMul/(modSun*modNormal))*57.29577951;
        angle = 90 - angle;
        //printf("%f    %f    %f    %f    %f    %f     %f    %f    %f \n",sun[0],sun[1],sun[2],satPos[0],satPos[1],satPos[2],satPos[3],satPos[4],satPos[5]);
    }
    //printf("time:%f, angle:  %f\n",time,angle);
    return angle;
}

// Calculate the point under the star in earth （计算星下点）
void calculateSatEarthPoint(double        greenStartTime,
                            double        initalGreenWithRa,
                            const double  currentTime,
                            double        sat[3],
                            double        result[2],
                            FILE          *infp ) {

    double greenWithDurTime = (currentTime - greenStartTime) * 86400.0;
    double greeWithPoint[3];
    double greenWith[2];
    double sat_point[2];

    locateSat(greeWithPoint, initalGreenWithRa, EARTH_RADIUS, 0.0, 86400.0, 0.0, greenWithDurTime);
    Cartesian2Equatorial(greeWithPoint, greenWith);
    Cartesian2Equatorial(sat, sat_point);

    result[0] = sat_point[0] - greenWith[0];
    result[1] = sat_point[1];

    if (result[0] > 180.) {
        result[0] = result[0] - 360.;
    }
    if (result[0] < -180.) {
        result[0] = result[0] + 360.;
    }

}
