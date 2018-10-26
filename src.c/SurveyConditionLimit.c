#include "SurveySim.h"


// 通过判断天区方位角与所选定的极深度天区（圆形）的中心的夹角的余弦来反推视向夹角，
// 如果夹角小于某个阈值，则该指向所对应的天区为极深度区域

#ifndef  _N_DEEP_
#define _N_DEEP_ 8  // 每次修改极深场区域的数目以后需要修改这个"_N_DEEP_"的数值！！！
#endif

double isInDeepSurveyArea(double ra, double dec, double diam) {
    double isInTheArea = -1;
    double r = 0.5*diam;   //单位是度，不是弧度

    // double areaCenter[_N_DEEP_][2] = { { 266.85, -5.59 },
    //                                     { 86.84, 5.52368 },
    //                                     {148.328147,57.301580},
    //                                     {41.125056,-45.196693},
    //                                     {151.417734,-9.375990},
    //                                     {179.937094,60.262282},
    //                                     {30.301973,-17.884744},
    //                                     {230,50},
    //                                     {310,-60},
    //                                     {350,45} };

/*
    double areaCenter[_N_DEEP_][2] = {  {266.85,     -5.59},
                                        {86.84,     5.52},
                                        {135.33,    40.30},
                                        {60,    -35.19},
                                        {151.42,    -9.37},
                                        {179.94,    40.26},
                                        {30.3,    -17.88},
                                        {230,      30},
                                        {310,      -30},
                                        {350,      45}};
*/

    // double areaCenter[_N_DEEP_][2] = { {320,    -25},
    //                                     {120,    23},
    //                                     {148.33,  57.30},
    //                                     {41.13,  -45.20},
    //                                     {180,     -25},
    //                                     {179.94,  60.26},
    //                                     {30.30,  -17.88},
    //                                     {230,     50},
    //                                     {310,    -60},
    //                                     {350,     45}};

// double areaCenter[_N_DEEP_][2] = {  {320,    -25},
//                                     {105,     55},
//                                     {155,     30},
//                                     {70,      -45},
//                                     {180,    -25},
//                                     {265,     75},
//                                     {30,      -22},
//                                     {250,     35},
//                                     {280,    -40},
//                                     {350,     40}};


    double areaCenter[_N_DEEP_][2] = { 
                                 //  {266.85,  -5.59},  // galactic center
                                 //  { 77.50,   6.09},  // anti-galatic center
                                  {345.97, -43.18},  // RLAIS S1
                                  { 31.04, -17.90},  // XMM-LSS
                                  { 40.29, -45.47},  // Extended Chandra Deep Filed-South (ECDFS)
                                  {150.70,  -9.39},  // COSMOS
                                  {148.33,  57.30},  // GOODS-N
                                  // { 41.13, -45.19},  // GOODS-S (very close to XMM-LSS)
                                  {179.94,  60.26},  // Extended Groth Strip (EGS)
                                  {230.00,  50.00},  // Arbitrarily chosen
                                  {310.00, -60.00}   // Arbitrarily chosen
                                   };

    int i = 0;
    double cr;
    for (i = 0; i < _N_DEEP_; i++) {
        cr = calculateAngle(ra, dec, areaCenter[i][0], areaCenter[i][1]);
        if (cr <= r) {
            isInTheArea = cr;
            break;
        }
    }

    return isInTheArea;
}

/*
 * 根据卫星和太阳的位置来判断是否在阳照区，在阳照区为1，否则为0
 */
int IsInSunSide(double sat[3], double sun[3]) {
    int in_sun_side = -999;
//	double angle = getAngle132(sat[0], sat[1], sat[2], sun[0], sun[1], sun[2], 0, 0, 0);
    double innerM_sat_sun = sat[0] * sun[0] + sat[1] * sun[1] + sat[2] * sun[2];
    double modSat = sat[0]*sat[0] + sat[1]*sat[1] + sat[2]*sat[2];
    double modSun = sun[0]*sun[0] + sun[1]*sun[1] + sun[2]*sun[2];
    double cosAngle = innerM_sat_sun / sqrt(modSat * modSun);
/*
    if (cosAngle >= -0.3385737) { // cos109.79
        in_sun_side = 1;
    }
*/
	in_sun_side = (int)(cosAngle >= -0.3385737);

	return in_sun_side;
}

/*
 * 根据进入和离开阴影区的时间来判断是否在阳照区，在阳照区为1，否则为0
 */
int IsInSunSideByTime(double shadow_start, double shadow_end, double curTime) {

    // if ( curTime >= shadow_start && curTime <= shadow_end ) {
    //     return 0;
    // } else {
    //     return 1;
    // }

    return (int)(!(curTime >= shadow_start && curTime <= shadow_end));
}

// 判断是否被地球遮挡
// (px,py,pz) 表示所指向的目标天区
int IsObscureByEarth_1( double sat[3], 
                        double sun[3], 
                        double p_x, 
                        double p_y, 
                        double p_z, 
                        double* angleValue ) {
    int isObscure = 1;
    int isInSunSide = 0;
    double modSat = sqrt(sat[0]*sat[0] + sat[1]*sat[1]+sat[2]*sat[2]);
    double modPoint = sqrt(p_x*p_x + p_y*p_y+p_z*p_z);
    double withLocalZenithAngle = (p_x*sat[0] + p_y*sat[1] + p_z*sat[2]) / (modPoint*modSat); //与卫星所在地天顶方向的夹角

    double innerM_sat_sun = sat[0] * sun[0] + sat[1] * sun[1] + sat[2] * sun[2];
    double cosAngle = innerM_sat_sun / (modSat * DIST_SUN_TO_EARTH); // 粗略估算？

    double sat0, sat1, sat2;
    double sat_sp_sat, sat_sp_point;
    double v1[3];
    double v1mod;
    double tan20=0.359825;
    double v2[3];
    double sun_sp_v2;
    
    if ( cosAngle < -0.3385737 ) { // cos109.79
        isInSunSide = -1; //在这个区域所看到的地球全部是阴影区
    } else if ( cosAngle >= -0.3385737 && cosAngle <= 0.3385737 ) {
        isInSunSide = 0;  //在这个区域可以看到被太阳照亮的部分，也可以看到阴影的部分
    } else {
        isInSunSide = 1;  //在这个区域只能看到被太阳照亮的部分
    }

#if defined(_USE_SWITCH_)

    switch (isInSunSide) {
        case 1://只能看到阴影区
            *angleValue = withLocalZenithAngle;
            // isObscure = (int)(withLocalZenithAngle <  COS_POINT_ZENITH_ANGLE_LIGHT_Min);
            isObscure = withLocalZenithAngle < COS_POINT_ZENITH_ANGLE_LIGHT_Min;
            break;
        case -1://只能看到阳照区
            *angleValue = MAX_VALUE;
            // isObscure = (int)(withLocalZenithAngle < COS_POINT_ZENITH_ANGLE_DARK);
            isObscure = withLocalZenithAngle < COS_POINT_ZENITH_ANGLE_DARK;
            break;
        case 0://同时看到阳照区和阴影区
            sat0 = sat[0] / modSat;
            sat1 = sat[1] / modSat;
            sat2 = sat[2] / modSat;
            sat_sp_sat = sat0 * sat0 + sat1 * sat1 + sat2 * sat2;
            sat_sp_point = sat0 * p_x + sat1 * p_y + sat2 * p_z;
            
            v1[0] = p_x * sat_sp_sat - sat0 * sat_sp_point;
            v1[1] = p_y * sat_sp_sat - sat1 * sat_sp_point;
            v1[2] = p_z * sat_sp_sat - sat2 * sat_sp_point;

            v1mod = sqrt( v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2] );

            v1[0] /= v1mod;
            v1[1] /= v1mod;
            v1[2] /= v1mod;

            v2[0] = tan20 * v1[0] + sat0;
            v2[1] = tan20 * v1[1] + sat1;
            v2[2] = tan20 * v1[2] + sat2;
            sun_sp_v2 = sun[0] * v2[0] + sun[1] * v2[1] + sun[2] * v2[2];

            // if ( sun_sp_v2 >= 0 ) {
            //     *angleValue = withLocalZenithAngle;
            //     // isObscure = (int)(withLocalZenithAngle <  COS_POINT_ZENITH_ANGLE_LIGHT_Min);
            //     isObscure = withLocalZenithAngle < COS_POINT_ZENITH_ANGLE_LIGHT_Min;
            // } else {
            //     *angleValue = MAX_VALUE;
            //     // isObscure = (int)(withLocalZenithAngle < COS_POINT_ZENITH_ANGLE_DARK);
            //     isObscure = withLocalZenithAngle < COS_POINT_ZENITH_ANGLE_DARK;
            // }

            switch ( (int)(sun_sp_v2 >= 0) ) {
                case 1:
                    *angleValue = withLocalZenithAngle;
                    isObscure = withLocalZenithAngle < COS_POINT_ZENITH_ANGLE_LIGHT_Min;
                    break;
                case 0:
                    *angleValue = MAX_VALUE;
                    isObscure = withLocalZenithAngle < COS_POINT_ZENITH_ANGLE_DARK;
                    break;
            }

            break;
        default:
            printf("failed to get correct isInSunSide!\n");
    }

#else
    if ( isInSunSide == 1 ) {
        *angleValue = withLocalZenithAngle;
        if ( withLocalZenithAngle >=  COS_POINT_ZENITH_ANGLE_LIGHT_Min ) {
            isObscure = 0;
        } else {
            isObscure = 1;
        }
    } else if ( isInSunSide == -1 ) {
        *angleValue = MAX_VALUE;//这个值有什么特殊意义么？
        if ( withLocalZenithAngle >= COS_POINT_ZENITH_ANGLE_DARK ) {
            isObscure = 0;
        } else {
            isObscure = 1;
        }
    } else {
        sat0 = sat[0] / modSat;
        sat1 = sat[1] / modSat;
        sat2 = sat[2] / modSat;

        sat_sp_sat = sat0 * sat0 + sat1 * sat1 + sat2 * sat2;
        sat_sp_point = sat0 * p_x + sat1 * p_y + sat2 * p_z;

        v1[0] = p_x * sat_sp_sat - sat0 * sat_sp_point;
        v1[1] = p_y * sat_sp_sat - sat1 * sat_sp_point;
        v1[2] = p_z * sat_sp_sat - sat2 * sat_sp_point;
        
        v1mod = sqrt(v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2]);
        v1[0] = v1[0]/v1mod;
        v1[1] = v1[1]/v1mod;
        v1[2] = v1[2]/v1mod;

        v2[0] = tan20 * v1[0] + sat0;
        v2[1] = tan20 * v1[1] + sat1;
        v2[2] = tan20 * v1[2] + sat2;

        sun_sp_v2 = sun[0] * v2[0] + sun[1] * v2[1] + sun[2] * v2[2];

        if ( sun_sp_v2 >= 0 ) {
            *angleValue = withLocalZenithAngle;
            if ( withLocalZenithAngle >=  COS_POINT_ZENITH_ANGLE_LIGHT_Min ) {
                isObscure = 0;
            } else {
                isObscure = 1;
            }
        } else {
            *angleValue = MAX_VALUE;
            if ( withLocalZenithAngle >= COS_POINT_ZENITH_ANGLE_DARK ) {
                isObscure = 0;
            } else {
                isObscure = 1;
            }
        }
    }
#endif
    return isObscure;
}

int IsObscureByEarth_1_old( double sat[3], 
                            double sun[3], 
                            double p_x, 
                            double p_y, 
                            double p_z, 
                            double* angleValue ) {
    int isObscure = 1;
    int isInSunSide = 0;
    double modSat = sqrt(sat[0]*sat[0] + sat[1]*sat[1]+sat[2]*sat[2]);
    double modPoint = sqrt(p_x*p_x + p_y*p_y+p_z*p_z);
    double withLocalZenithAngle = (p_x*sat[0] + p_y*sat[1] + p_z*sat[2]) / (modPoint*modSat); //与卫星所在地天顶方向的夹角

    double innerM_sat_sun = sat[0] * sun[0] + sat[1] * sun[1] + sat[2] * sun[2];
    double cosAngle = innerM_sat_sun / (modSat * DIST_SUN_TO_EARTH); // 粗略估算？

//  这一部分需要继续改进
    if ( cosAngle < -0.3385737 ) { // cos109.79
        isInSunSide = -1; //在这个区域所看到的地球全部是阴影区
    } else if ( cosAngle >= -0.3385737 && cosAngle <= 0.3385737 ) {
        isInSunSide = 0;  //在这个区域可以看到被太阳照亮的部分，也可以看到阴影的部分
    } else {
        isInSunSide = 1;  //在这个区域只能看到被太阳照亮的部分
    }

    if ( isInSunSide == 1 ) {
        *angleValue = withLocalZenithAngle;
        if ( withLocalZenithAngle >=  COS_POINT_ZENITH_ANGLE_LIGHT_Min ) {
            isObscure = 0;
        } else {
            isObscure = 1;
        }
    } else if ( isInSunSide == -1 ) {
        *angleValue = MAX_VALUE;//这个值有什么特殊意义么？
        if ( withLocalZenithAngle >= COS_POINT_ZENITH_ANGLE_DARK ) {
            isObscure = 0;
        } else {
            isObscure = 1;
        }
    } else {
        double sat0 = sat[0] / modSat;
        double sat1 = sat[1] / modSat;
        double sat2 = sat[2] / modSat;

        double sat_sp_sat = sat0 * sat0 + sat1 * sat1 + sat2 * sat2;
        double sat_sp_point = sat0 * p_x + sat1 * p_y + sat2 * p_z;
        double v1[3] = { p_x * sat_sp_sat - sat0 * sat_sp_point,
                         p_y * sat_sp_sat - sat1 * sat_sp_point,
                         p_z * sat_sp_sat - sat2 * sat_sp_point };
        double v1mod = sqrt(v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2]);
        v1[0] = v1[0]/v1mod;
        v1[1] = v1[1]/v1mod;
        v1[2] = v1[2]/v1mod;

        double tan20 = 0.359825; //tan 19.79
        double v2[3] = { tan20 * v1[0] + sat0, tan20 * v1[1] + sat1, tan20 * v1[2] + sat2 };
        double sun_sp_v2 = sun[0] * v2[0] + sun[1] * v2[1] + sun[2] * v2[2];

        if ( sun_sp_v2 >= 0 ) {
            *angleValue = withLocalZenithAngle;
            if ( withLocalZenithAngle >=  COS_POINT_ZENITH_ANGLE_LIGHT_Min ) {
                isObscure = 0;
            } else {
                isObscure = 1;
            }
        } else {
            *angleValue = MAX_VALUE;
            if ( withLocalZenithAngle >= COS_POINT_ZENITH_ANGLE_DARK ) {
                isObscure = 0;
            } else {
                isObscure = 1;
            }
        }
    }

    return isObscure;
}

bool debug_IsObscureByEarth( double sat[3], 
							 double sun[3], 
							 double p_x, 
							 double p_y, 
							 double p_z ){

    int test_result_old, test_result_new;
    double angleValue;

    test_result_old = IsObscureByEarth_1_old(sat,sun,p_x,p_y,p_z,&angleValue);
    test_result_new = IsObscureByEarth_1(sat,sun,p_x,p_y,p_z,&angleValue);

    // bool cmp_result = false;
    // if( test_result_new == test_result_old ){
    //      cmp_result = true;
    // }

    // return cmp_result;

    return (test_result_new == test_result_old);
}

/**
 * obscure return 2
 * or return cos(sun_point_angle);
 **/
double IsObscureBySun(  double sun[3], 
                        double px, 
                        double py, 
                        double pz ) {
    double sunObscurFlag = 0;
    double modP = px*px + py*py + pz*pz;
    double modSun = sun[0]*sun[0] + sun[1]*sun[1] + sun[2]*sun[2];
    double cosA = (sun[0] * px + sun[1] * py + sun[2] * pz) / sqrt(modSun * modP);

    // sunObscurFlag = cosA;
    // if (cosA > COS_SUN_POINT_ANGLE) {
    //     sunObscurFlag = 2;
    // }
    // else {
    //     sunObscurFlag = cosA;
    // }

    sunObscurFlag = (cosA > COS_SUN_POINT_ANGLE) ? 2 : cosA;

    return sunObscurFlag;
}

/**
 * obscure return 2
 * or return cos(moon_point_angle);
 **/
double IsObscureByMoon( double moon[3], 
                        double px, 
                        double py, 
                        double pz) {
    double moonObscureFlag = 0;
    double modP = px*px + py*py + pz*pz;
    double modMoon = moon[0]*moon[0] + moon[1]*moon[1] + moon[2]*moon[2];
    double cosA = (moon[0] * px + moon[1] * py + moon[2] * pz) / sqrt(modMoon * modP);

    // moonObscureFlag = cosA;
    // if (cosA > COS_MOON_POINT_ANGLE) {
    //     moonObscureFlag = 2;
    // }
    // else {
    //     moonObscureFlag = cosA;
    // }

    moonObscureFlag =  (cosA > COS_MOON_POINT_ANGLE) ? 2 : cosA;


    // COS_MOON_POINT_ANGLE = cos(40*M_PI/180);
    // printf("@debug modP    = %10.8f\n",modP);
    // printf("@debug modMoon = %10.8f\n",modMoon);

    // printf("@debug COS_MOON_POINT_ANGLE = %10.8f\n",COS_MOON_POINT_ANGLE);
    // printf("@debug cosA = %10.8f\n", cosA);
    // printf("@debug angle = %10.8f\n",acos(cosA)*180/M_PI);

    return moonObscureFlag;
}

/**
 * in saa = 1;
 * no in saa = 0;
 */
int IsInSAA( double uderStarPoint[2] ) {

    double rangePoint[5][2] = { { -13, -42.5 },
                                { -72, -42.5 },
                                { -81, -25 },
                                { -54, -15 },
                                { 5, -32 } };
    int newFlag = 1;
    int i, next;
    double directFlag[5];
    double a1,a2,b1,b2;

    for ( i = 0; i < 5; i++ ) {
        a1 = uderStarPoint[0] - rangePoint[i][0];
        a2 = uderStarPoint[1] - rangePoint[i][1];
        next = i + 1;

        if (next > 4) {
            next = next - 5;
        }

        b1 = rangePoint[next][0] - rangePoint[i][0];
        b2 = rangePoint[next][1] - rangePoint[i][1];
        directFlag[i] = a1 * b2 - a2 * b1;
    }

    for ( i = 0; i < 5; i++ ) {
        next = i + 1;
        if ( next > 4 ) {
            next = next - 5;
        }
        if ( directFlag[i] * directFlag[next] < 0 ) {
            newFlag = 0;
            break;
        }
    }

    return newFlag;
}

//  calculate rotate time of tel point
double calculateTransTime(double transAngle) {

#if defined(_SLEW_TIME_GSL_INTERP_)
    return gsl_spline_eval(spline_trans,transAngle,acc_trans);
#else

    //double data[9][2] = { {0.5, 70},{1, 80}, {5, 95},
    //{ 10, 105 }, { 15, 115 }, { 20, 120 }, { 30, 135 }, { 45, 150 },{180,200} };

    //double data[4][2] = { {0.1, 70},{1, 80}, {45, 161},{180,200} };
    // double data[4][2] = { {1, 80}, {20,127},{45, 196},{180,581} };
    double data[4][2] = { {1, 45}, {20,92},{45, 196},{180,581} };  // 减少稳定时间
    // double data[4][2] = { {1, 76}, {20,123},{45, 192},{180,577} };
    // double data[3][2] = { {1, 80}, {45,170},{180,445} };
    int i = 0;
    double tTime = 0;
    
    // if(transAngle > 180) {
    //     printf("%f \n",transAngle);
    // }

    if(transAngle < 1) {
        tTime = 70;
    } else if(transAngle == 1) {
        tTime = 80;
    } else {
        for(i = 0; i < 3 ; i ++) {
            if(transAngle>data[i][0] && transAngle <= data[i + 1][0] ) {
                tTime = data[i][1] * ((transAngle - data[i+1][0])) / (((data[i][0] - data[i+1][0]))) 
                      + data[i+1][1] * ((transAngle - data[i][0])) / (((data[i+1][0] - data[i][0])));
                break;
            }
        }
    }
    return tTime;
#endif

}

double getExposureTime_Zodical(double p_lat, double p_lon2sun) {

    return EXTIME; // 简化处理

//  以下部分是为了根据黄道光强度来判断可以进行曝光的时间。
    double  table_times[4][4]= { {298,168,146,138},
                                 {164,146,138,138},
                                 {160,146,138,138},
                                 {173,150,140,138} };
    double lon[4] = {50,110,145,180};
    double lat[4] = {0,30,60,90};

    int lati, loni,i;

    if (p_lon2sun<lon[0]) {
        p_lon2sun = lon[0];
    }

    lati = 0;

    if (p_lat == lat[0]) {
        lati = 0;
    }

    for(i = 0; i < 3; i ++) {
        if (p_lat > lat[i] && p_lat <= lat[i+1]) {
            lati = i;
            break;
        }
    }

    loni = 0;
    if (p_lon2sun == lon[0]) {
        loni = 0;
    }

    for(i = 0; i < 3; i ++) {
        if (p_lon2sun > lon[i] && p_lon2sun <= lon[i+1]) {
            loni = i;
            break;
        }
    }

    double k = (table_times[loni+1][lati]-table_times[loni][lati])/(lon[loni+1]-lon[loni]);
    double y1 = k*(p_lon2sun-lon[loni])+table_times[loni][lati];

    k = (table_times[loni+1][lati+1]-table_times[loni][lati+1])/(lon[loni+1]-lon[loni]);
    double y2 = k*(p_lon2sun-lon[loni])+table_times[loni][lati+1];

    k = (y2-y1)/(lat[lati+1]-lat[lati]);

    // printf("%d    %d \n",lati,loni);

    double y = k*(p_lat-lat[lati])+y1;
    return y;
}


/*
cosAngle: cos value of angle, normal of plane to sun
*/
// double getPlaneEnergy(double cosAngle, double yr) {
// //  计算发电功率。
// //  值得注意的是，这个函数依赖于帆板法线与太阳方位的夹角
//     double yr_s = 0.0;
//     double yr_e = 10.0;
//     double energy_s[5] = {1.6478530148909550,1.6353357198900544,1.5981641665618982,1.5374677937446495,1.4550908300125129};
//     double energy_e[5] = {1.6332290528848927,1.6208228433053604,1.5839811709669922,1.5238234514410043,1.4421775465939182};
//     double angleCos[5] = {0.9961947,0.9848078,0.9659258,0.9396926,0.0};
//     double slop, energy_val = 0.;
//
//     int i;
//     for (i = 0; i < 5; i ++) {
//         if(fabs(cosAngle)>=angleCos[i]) {
//             slop = (energy_e[i] - energy_s[i])/(yr_e - yr_s);
//             energy_val = slop * (yr - yr_s)+energy_s[i];
//             break;
//         }
//     }
//
//     return energy_val;
// }


//  新版本是依据文档中提供的参数所设计
double getPlaneEnergy(double cos_value, double yr) {
//   return 19.11799422*pow(cos_value,3)*(1.-0.007*yr);
	 return 63.7266474*pow(cos_value,4)*(0.3-0.007*yr);
    // return 63.7266474*pow(cos_value,3)*(0.3-0.007*yr);
}

//get the time in earth shadow and the time out earth shadow
// 这个函数可以优化改进一下，使得获取的时间精度更高
void aquireShadowTime(double* shadow_start, double* shadow_end, double curTime) {

    int inSunflage_sdw = 0;
    double sdw_sat[3],sdw_sun[3];
    double shadow_n_time = curTime;

    if( locateSat1(sdw_sat, shadow_n_time, orbitData, orbitDataNum,Orbit_File_Num) < 0 ){
        return;
    }

    locateSun(infp, shadow_n_time, sdw_sun);
    inSunflage_sdw = IsInSunSide(sdw_sat, sdw_sun);

    if (inSunflage_sdw == 0) {
        *shadow_start = shadow_n_time;
    } else {
        while(inSunflage_sdw==1) {
            shadow_n_time +=10.0/86400.0;
            // locateSat1(sdw_sat, shadow_n_time, orbitData, orbitDataNum,Orbit_File_Num);
            // if(sdw_sat == NULL) {
            if( locateSat1(sdw_sat, shadow_n_time, orbitData, orbitDataNum,Orbit_File_Num) < 0 ) {
                *shadow_start = shadow_n_time;
                *shadow_end = shadow_n_time;
                return;
            }
            locateSun(infp, shadow_n_time, sdw_sun);
            inSunflage_sdw = IsInSunSide(sdw_sat, sdw_sun);
        }
        *shadow_start = shadow_n_time;
    }
    
    while(inSunflage_sdw==0) {
        shadow_n_time +=10.0/86400.0;
        if( locateSat1(sdw_sat, shadow_n_time, orbitData, orbitDataNum,Orbit_File_Num) < 0 ){
            break;
        }
// #if defined(_FULL_DEBUG_)
        // if(sdw_sat == NULL) {// 这个判断是否必要？？？
        //     break;
        // }
// #endif
        locateSun(infp, shadow_n_time, sdw_sun);
        inSunflage_sdw = IsInSunSide(sdw_sat, sdw_sun);
    }
    *shadow_end = shadow_n_time;
    // printf("%f    \n",(*shadow_end-*shadow_start)*1440);
}
