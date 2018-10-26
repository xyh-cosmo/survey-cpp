#include "SurveySim.h"

int main(int argc, char *argv[]) {

	init_default();
	parseConfigFile(argv[1]);

	double jdate_origin = 2459766;
	double jdate_start = atof(argv[2]);
	double jdate_end   = atof(argv[3]);

    infp = NULL;
    double earthStartTime, initalGreenWithRa = 0;

    if ((infp = fopen("jpl.405", "r")) == NULL ) {
        fprintf(stderr, "\nERROR: Can't open ephemeris file %s for input.\n\n", argv[1]);
        return 0;
    }

    //计算某一时刻地区本初子午线与赤道交点午夜时在天球赤道坐标系中的位置
    earthStartTime = floor(jdate_origin - 0.5) + 0.5; // midnight ut time (unit day)
    double sunMid[3];
    locateSun(infp, earthStartTime, sunMid); // 太阳的笛卡尔坐标，以地球为原点(J2000)

    double reEarth[2]; // (ra,dec) of Sun
    Cartesian2Equatorial(sunMid, reEarth);  // get (ra,dec) of Sun

    if (reEarth[0] >= 180) {
        initalGreenWithRa = reEarth[0] - 180;
    } else {
        initalGreenWithRa = reEarth[0] + 180;
    }

//  读取所有的轨道数据
    orbitData = (double**)malloc(sizeof(double*) * Orbit_File_Num);
    readAllorbitsFile(orbitData, orbitDataNum);

	int num=10;
	double jt = jdate_start;
	while( jt <= jdate_end ){
		double sat[3];
		locateSat1(sat,jt,orbitData,orbitDataNum,Orbit_File_Num);
		double underStarPoint[2];
        calculateSatEarthPoint(earthStartTime, initalGreenWithRa, jt, sat, underStarPoint, infp);
        
        printf("Jdate: %14.10g, ra = %15.10g dec = %15.10g\n",jt,underStarPoint[0],underStarPoint[1]);
		jt += (jdate_end-jdate_start)/num;
	}

    free(orbitData);

    return 0;
}
