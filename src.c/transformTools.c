#include "SurveySim.h"

/*
*太阳帆板可旋转，计算一定旋转角度范围内，帆板法线的范围
*/
void get_solar_normal_30(double px, double py, double pz, double nx, double ny, double nz,double* n30, double* n30_) {

    double angle=PANEL_TRANSE_ANGLE/57.29577951;
    double a = cos(angle);
    double b = sin(angle);
    double p = sqrt(px*px+py*py+pz*pz);
    double n = sqrt(nx*nx+ny*ny+nz*nz);
    double t_nx30 = a*nx*p/n + b*px;
    double t_ny30 = a*ny*p/n + b*py;
    double t_nz30 = a*nz*p/n + b*pz;

    double t_nx30_ = a*nx*p/n - b*px;
    double t_ny30_ = a*ny*p/n - b*py;
    double t_nz30_ = a*nz*p/n - b*pz;
    n30[0] = t_nx30;
    n30[1] = t_ny30;
    n30[2] = t_nz30;

    n30_[0] = t_nx30_;
    n30_[1] = t_ny30_;
    n30_[2] = t_nz30_;
}

/**
 * 坐标变化沿某一平面旋转angle，angle单位为 度
 * coorO表示原始坐标
 * coorR存储变换后的结果
 * flag 标识 1：x-y 平面 2： x-z平面 3：y-z 平面
 * 转动方向为：例如x-y平面，逆时针转动
 */
void CoordinateSpin(double coorO[3], double coorR[3], double angle, int flag) {

    double arcAngle = angle * PI_180;

    if (flag > 3 || flag < 1) { // after debug, this should be commented out!
        return;
    }

    if (flag == 1) {
        coorR[0] = coorO[0] * cos(arcAngle) + coorO[1] * sin(arcAngle);
        coorR[1] = coorO[0] * (-sin(arcAngle)) + coorO[1] * cos(arcAngle);
        coorR[2] = coorO[2];
    }
    if (flag == 2) {
        coorR[0] = coorO[0] * cos(arcAngle) + coorO[2] * sin(arcAngle);
        coorR[1] = coorO[1];
        coorR[2] = coorO[0] * (-sin(arcAngle)) + coorO[2] * cos(arcAngle);
    }
    if (flag == 3) {
        coorR[0] = coorO[0];
        coorR[1] = coorO[1] * cos(arcAngle) + coorO[2] * sin(arcAngle);
        coorR[2] = coorO[1] * (-sin(arcAngle)) + coorO[2] * cos(arcAngle);
    }
}

// 绕着x轴逆时针转动
void CoordinateSpin_x(double coorO[3], double coorR[3], double angle ) {

    double arcAngle = angle * PI_180;
    coorR[0] = coorO[0];
    coorR[1] = coorO[1] * cos(arcAngle) + coorO[2] * sin(arcAngle);
    coorR[2] = coorO[1] * (-sin(arcAngle)) + coorO[2] * cos(arcAngle);
}

// 从赤道坐标系变换到黄道坐标系
void CoordinateSpinEquatorial2Ecliptic(double coorO[3], double coorR[3]) {
    // 沿着x轴旋转
    //angle = 23.4522
    //cos(angle) = 0.917392
    //sin(angle) = 0.394984

    double cosa = 0.917392;
	double sina = 0.397983;

    coorR[0] = coorO[0];
    coorR[1] = coorO[1] * cosa + coorO[2] * sina;
    coorR[2] =-coorO[1] * sina + coorO[2] * cosa;
}

/**
 * *eCoor = ra, *eCoor+1 = dec
 */
// 从笛卡尔坐标系变换到（天文的）球坐标系
// void Cartesian2Equatorial(double* carCoor, double* eCoor) {
//     double x1 = carCoor[0], x2 = carCoor[1], x3 = carCoor[2];
//     double r = sqrt(x1*x1+x2*x2+x3*x3);
//     double theta = asin(x3/r);

//     double tmp = r*cos(theta)+x1;
//     if( fabs(tmp) < 1e-15 ) // 避免数值问题导致的NAN
//         tmp = 1e-15;

//     *(eCoor+1) = theta*180./PI;
//     *(eCoor+0) = atan(x2/tmp) *360./PI;

// //  ra只能是正数！
//     if( *eCoor < 0 )
//         *eCoor += 360.;
// }

//  新版本的误差控制更好！
void Cartesian2Equatorial(double* carCoor, double* eCoor) {
    double x1 = carCoor[0], x2 = carCoor[1], x3 = carCoor[2];
    double r = sqrt(x1*x1+x2*x2+x3*x3);
    double theta = asin(x3/r);

    double phi=0;
    if( fabs(x1/r) > 1e-8 && fabs(x2/r) > 1e-8 ){
        phi = 2.0*atan(x2/(r*cos(theta)+x1));
    }
    else{
        double x = x1*1e8;
        double y = x2*1e8;
        phi = 2.0*atan(y/(sqrt(x*x+y*y)+x));
    }

    *(eCoor+0) = phi*180.0/PI;
    *(eCoor+1) = theta*180.0/PI;

//  ra只能是正数！
    if( *eCoor < 0 )
        *eCoor += 360.;
}


double getAngle132(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3) {

    double cosValue = 0;
    double angle = 0;

    double x11 = x1-x3;
    double y11 = y1-y3;
    double z11 = z1-z3;

    double x22 = x2-x3;
    double y22 = y2-y3;
    double z22 = z2-z3;

    double tt = sqrt((x11*x11 + y11*y11 + z11* z11) * (x22*x22 + y22*y22 + z22*z22));

#if defined(_FULL_DEBUG_)
    if(tt==0) {//实际应用的过程中不可能出现tt=0的情况，因为（）和（）不可能有至少一个出现在原点。
        return 0;
    }
#endif

    cosValue = (x11*x22+y11*y22+z11*z22)/tt;

    if (cosValue > 1) {
        cosValue = 1;
    }
    if (cosValue < -1) {
        cosValue = -1;
    }
    angle = acos(cosValue);
    return angle * 180 * M_1_PI;
}

double calculateAngle(double ra1, double dec1, double ra2, double dec2) {

    double x1, y1, z1, x2, y2, z2, angle, cosValue;
    x1 = cos(dec1 * PI_180) * cos(ra1 * PI_180);
    y1 = cos(dec1 * PI_180) * sin(ra1 * PI_180);
    z1 = sin(dec1 * PI_180);
    x2 = cos(dec2 * PI_180) * cos(ra2 * PI_180);
    y2 = cos(dec2 * PI_180) * sin(ra2 * PI_180);
    z2 = sin(dec2 * PI_180);

    cosValue = x1*x2+y1*y2+z1*z2;

    if( cosValue < -1.0 ){
      return 180.0;
    }else if( cosValue > 1.0 ){
      return 0.0;
    }

    angle = acos(x1*x2+y1*y2+z1*z2);
    return angle * 180 * M_1_PI;  //返回值为角度（以度数为单位，不是以弧度为单位）
}
