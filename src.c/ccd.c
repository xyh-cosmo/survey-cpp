#include "SurveySim.h"

void init_ccd_pos(double ccd_pos_in_focus[18][2]){
    //----------------NUV  4----------------------//
    ccd_pos_in_focus[0][0] =   0.0030/28.;  //单位球上的坐标？
    ccd_pos_in_focus[0][1] =  -0.1475/28.;
    ccd_pos_in_focus[1][0] =   0.0030/28.;
    ccd_pos_in_focus[1][1] =  -0.0462/28.;
    ccd_pos_in_focus[2][0] =  -0.0952/28.;
    ccd_pos_in_focus[2][1] =  -0.0462/28.;
    ccd_pos_in_focus[3][0] =  -0.0952/28.;
    ccd_pos_in_focus[3][1] =   0.0551/28.;

    //--------------------u  2---------------------//
    ccd_pos_in_focus[4][0] = -0.0952/28.;
    ccd_pos_in_focus[4][1] = -0.1475/28.;

    ccd_pos_in_focus[5][0] =  0.0030/28.;
    ccd_pos_in_focus[5][1] =  0.0551/28.;

    //--------------------g  2---------------------//
    ccd_pos_in_focus[6][0] =  0.1012/28.;
    ccd_pos_in_focus[6][1] = -0.0462/28.;

    ccd_pos_in_focus[7][0] = -0.1934/28.;
    ccd_pos_in_focus[7][1] = -0.0462/28.;

    //--------------------r  2---------------------//
    ccd_pos_in_focus[8][0] = -0.1934/28.;
    ccd_pos_in_focus[8][1] = -0.1475/28.;

    ccd_pos_in_focus[9][0] =  0.1012/28.;
    ccd_pos_in_focus[9][1] =  0.0551/28.;

    //--------------------i  2---------------------//
    ccd_pos_in_focus[10][0] =  0.1012/28.;
    ccd_pos_in_focus[10][1] = -0.1475/28.;

    ccd_pos_in_focus[11][0] = -0.1934/28.;
    ccd_pos_in_focus[11][1] =  0.0551/28.;

    //--------------------z  2---------------------//
    ccd_pos_in_focus[12][0] = -0.0952/28.;
    ccd_pos_in_focus[12][1] = -0.2488/28.;

    ccd_pos_in_focus[13][0] =  0.0030/28.;
    ccd_pos_in_focus[13][1] =  0.1564/28.;

    //--------------------Y  4---------------------//
    ccd_pos_in_focus[14][0] =  0.1012/28.;
    ccd_pos_in_focus[14][1] = -0.2488/28.;

    ccd_pos_in_focus[15][0] =  0.0030/28.;
    ccd_pos_in_focus[15][1] = -0.2488/28.;

    ccd_pos_in_focus[16][0] = -0.0952/28.;
    ccd_pos_in_focus[16][1] =  0.1564/28.;

    ccd_pos_in_focus[17][0] = -0.1934/28.;
    ccd_pos_in_focus[17][1] =  0.1564/28.;
}