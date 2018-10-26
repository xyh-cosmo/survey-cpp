// 这个版本基于SearchNewPoints_new2.c,但是策略部分有更大的改动：
// 当搜索观测天区失败时，将望远镜的指向往天顶方向或者望远镜的运动方向旋转5度。
// 转动5度所需的时间为90秒。
//
// 注意：目前的版本里，5度的调整仅仅在经过SAA区域和搜寻天区失败的时候才进行，而在停靠和
// beta角限制的时候不进行调整。

#include "SurveySim.h"

//  新版本的搜索函数增加了两个输入变量:cur_ra_tel,cur_dec_tel，用于记录当前望远镜的指向。
//  旧版本中，望远镜当前的某一时刻的指向总是指向上一次观测天区的指向。
void FindNewPointByAllSearch3(double* 	    currentTime,
                              int* 		    id,
                              double*       cur_ra_tel,
                              double*       cur_dec_tel,
                              SKY_Coord* 	skyMap,
                              int 		    all_sky_num,
                              FILE*		    infp,
                              double 		greenStartTime,
                              double		initalGreenWithRa,
                              double*		outputresult,
                              double*		outputresult_fail,
                              double**	    orbitData,
                              int*		    orbitDataNum,
                              int			p_rank,
                              int           p_size,
                              FILE*         fp_drift ) {
    int i; //, j;
    int isInSunSide = 0;
    double returnTime = 0;
    double extraDeep = 5;
    double extraUDeep = 12;
    double curTime = 0;
    double saaTime = 0;
    double sat[3], sun[3], moon[3],sun_lon_lat[2];
    double oldCurTime = *currentTime;

    curTime = *currentTime;
    isInSunSide = IsInSunSideByTime( shadow_st, shadow_et, curTime );

#if defined(_TURN_ON_5DEG_DRIFT_)
    double velocity_tel[3]; //望远镜飞行速度，实际上只是两个时间点之间的位移，可以用来指引望远镜指向的改变。
#endif

    if( locateSat1(sat, curTime, orbitData, orbitDataNum,Orbit_File_Num) > 0 ){

        double underStarPoint[2];
        calculateSatEarthPoint(greenStartTime, initalGreenWithRa, curTime, sat, underStarPoint, infp);

#if defined(_TURN_ON_5DEG_DRIFT_)
        double sat_pos_now[3],sat_pos_next[3];
        CoordinateSpinEquatorial2Ecliptic(sat, sat_pos_now);
#endif
        while ( IsInSAA(underStarPoint) == 1) {
        //  如果处在SAA区域，那么就逐步将时间往后按照步长jump_time（单位是秒）递进，直到离开SAA区域为止。
            double lastSAATime = curTime;
            curTime = curTime + jump_time / 86400.0;

            if( locateSat1(sat, curTime, orbitData, orbitDataNum,Orbit_File_Num) < 0 ){
                continue;
            }

#if defined(_TURN_ON_5DEG_DRIFT_)
            CoordinateSpinEquatorial2Ecliptic(sat, sat_pos_next);
#endif

            double sun_front_time = 0; // 1 orbit, 一轨结束前剩下的阳照时间（单位是天）
            double sun_back_time = 0;  // 1 orbit, 一轨结束前剩下的阴影时间（单位是天）

            // 判断剩下的阳照时间以及阴影时间，用于计算电池电量
            if( curTime > shadow_et && lastSAATime < shadow_et ){
                sun_front_time = curTime - shadow_et;
                sun_back_time = shadow_et - lastSAATime;
            } else if( curTime > shadow_st && lastSAATime < shadow_st ){
                sun_front_time = shadow_st - lastSAATime;
                sun_back_time = curTime - shadow_st;
            } else if( curTime <= shadow_st ){
                sun_front_time = curTime - lastSAATime;
                sun_back_time = 0;
            } else if( lastSAATime >= shadow_st ){
                sun_front_time = 0;
                sun_back_time = curTime - lastSAATime;
            }

            //  计算电池电量
            double yr = (curTime-2459766.0)/365.25;

            battery_q = battery_q
                      + sun_front_time*86400.0*getPlaneEnergy(1.0,yr)
                      - jump_time*POWER_CONSUMPTION;

            if( battery_q > BATTERY_MAX ){  // 电池电量存在一个最大值
                battery_q = BATTERY_MAX;
            }

            if( battery_q < BATTERY_LOW ){  // 50这个数字是怎么来的？？？
                printf("Energy not balance, Position 3,battery is %f, sun time is %f, shadow time is %f\n",battery_q,sun_front_time,sun_back_time);
                printf("shadow s is %f , e is %f, curTime is %f, old time is %f \n",shadow_st,shadow_et,curTime,oldCurTime);
                MPI_Finalize();
                exit(0);
            }

            // 如果已经离开阴影区，那么需要获取下一次进入和离开阴影区的时间
            if( curTime > shadow_et ){
                aquireShadowTime( &shadow_st, &shadow_et, curTime);
            }

            calculateSatEarthPoint(greenStartTime, initalGreenWithRa, curTime, sat, underStarPoint, infp);

#if defined(_TURN_ON_5DEG_DRIFT_)
            get_tel_velocity( sat_pos_next, sat_pos_now, velocity_tel );

            double drift_angle;
            Drift_by_5deg( cur_ra_tel, 
                           cur_dec_tel, 
                           sat_pos_next,
                           velocity_tel,
                           &drift_angle,
                           p_rank );

#if defined(_USE_CMG_ONE_ORBIT_STATE_)
            CMG_one_orbit_state_update(cmg_state, curTime, 90, 0, 0, p_rank);
#endif

            Write_drift_state(  fp_drift,
                                p_rank,
                                curTime, *cur_dec_tel, *cur_ra_tel,
                                -999, -999, -999,
                                -999, -999, -999,
                                -999, -999, -999,
                                0, 0, -999,
                                -999,
                                -999,
                                drift_angle,
                                -999,
                                0,
                                battery_q,
                                -999,
                                -999,
                                -999 );
#endif
        }

        if( sat[0]==0.0 && sat[1]==0.0 && sat[2]==0.0 ) {
            outputresult[0] = -1;
            returnTime = curTime + jump_time / 86400.0;
            *currentTime = returnTime;
            return;
        }

        isInSunSide = IsInSunSideByTime(shadow_st, shadow_et,curTime);
        saaTime = curTime-oldCurTime;
        locateSun(infp, curTime, sun);
        locateMoon(infp, curTime, moon);

        double sun_t[3], moon_t[3], sat_t[3];

//      从赤道坐标系转换到黄道坐标系
        CoordinateSpinEquatorial2Ecliptic(sat, sat_t);
        CoordinateSpinEquatorial2Ecliptic(sun, sun_t);
        CoordinateSpinEquatorial2Ecliptic(moon, moon_t);

        sat[0] = sat_t[0];
        sat[1] = sat_t[1];
        sat[2] = sat_t[2];
        sun[0] = sun_t[0];
        sun[1] = sun_t[1];
        sun[2] = sun_t[2];
        moon[0] = moon_t[0];
        moon[1] = moon_t[1];
        moon[2] = moon_t[2];

        Cartesian2Equatorial(sun,sun_lon_lat);
    }

    if( sat[0]==0.0 && sat[1]==0.0 && sat[2]==0.0 ) {
        outputresult[0] = -1;
        returnTime = curTime + jump_time / 86400.0;
        *currentTime = returnTime;
        return;
    }

    int final_index = 0;
    double final_use_time = 0;
    double final_exp_time = 0;
    double final_energy_time=0;
    double final_min_value = MAX_VALUE;
    double final_Angle = 0;
    double final_cos_sun_plane_angle=2;
    double final_cos_panel_point_angle = 2;
    double final_cmg_value = 0;

    // double final_left_1o_energy = 0;
    // double final_temp_sun = 0;
    // double final_temp_shadow = 0;

    struct {
        double val_weight;
        int n_rank;
    } rank_result, final_result;
    rank_result.val_weight = MAX_VALUE;
    rank_result.n_rank=p_rank;

    double local_fail_out[10];
    setDoubleArrayZero(local_fail_out,10);

//  计算sky_id_tracker[]的指标范围
    int idx_start, idx_end, idx_size;
    int sky_num_per_rank = sky_num_remained / p_size;
    int idx_remained     = sky_num_remained % p_size;

    if( sky_num_per_rank >= 1 ){
        idx_size = sky_num_per_rank;
        if( p_rank < idx_remained ){
            idx_size    += 1;
            idx_start   = p_rank*idx_size;
            idx_end     = idx_start + idx_size;
        }
        else{
            idx_start   = (p_rank-idx_remained)*idx_size + idx_remained*(idx_size+1);
            idx_end     = idx_start + idx_size;
        }
    }
    else{
        idx_size = 1;
        if( p_rank < idx_remained ){
            idx_start = p_rank;
            idx_end   = idx_start + idx_size;
        }
        else{
            idx_start = 0;
            idx_end   = 1;
        }
    }

	// 测试指标的计算结果是否正确
    // printf("rank = %3d, idx_start = %8d, idx_end = %8d\n", p_rank, idx_start, idx_end);

    // MPI_Barrier(MPI_COMM_WORLD);
    // MPI_Finalize();
    // exit(0);

    double curTime_yr = (curTime-2459766)/365.25;   //将当前的儒略日转换为【0~10】中的年份
    double dist_sun = sqrt( sun[0]*sun[0] + sun[1]*sun[1] + sun[2]*sun[2] );

    if ( isInSunSide == 1 ) {// 在阳照区时

        for( i=idx_start; i<idx_end; i++ ){

            int s_id = sky_id_tracker[i];

            // 如果该天区对应的黄经大于360，那么直接跳过本次循环（不过这个情况有可能发生吗？？？）
            // 已经测试过，这种情况不会发生
            // if (skyMap[s_id].ra > 360) {
            //     continue;
            // }

            // 如果该天区不可观测，那么直接跳过本次循环
            // if (skyMap[s_id].isObserve == 0) {
            //     continue;
            // }

            
            //  ================================================================================
            //  将下面三个判断放到此处可以加速模拟速度
            if ( skyMap[s_id].flag >= skyMap[s_id].maxCoverNum ) {
                local_fail_out[7] += 1.0;
                continue;
            }

            if ( IsObscureBySun(sun,skyMap[s_id].x,skyMap[s_id].y, skyMap[s_id].z) > 1 ) {
                local_fail_out[1] += 1.0;
                continue;
            }
            
            if ( IsObscureByMoon(moon, skyMap[s_id].x,skyMap[s_id].y, skyMap[s_id].z) > 1 ) {
                local_fail_out[2] += 1.0;
                continue;
            }
            //  ================================================================================

            double value_panel_point_angle = 0;  // NOT USED!!!
            double value_sun_angle_to_nomal= 0;
            int panel_is_ok = -1;

#if defined(_ENABLE_PANEL_ROTATION_)
            panel_is_ok = TestPanelAngle( sun,
                                          dist_sun,
                                          skyMap[s_id].x,
                                          skyMap[s_id].y,
                                          skyMap[s_id].z,
                                          skyMap[s_id].nx,
                                          skyMap[s_id].ny,
                                          skyMap[s_id].nz,
                                          &value_sun_angle_to_nomal,
                                          25,
                                          p_rank );
#else
            // 此处简化处理，假设帆板不转动
            double thres_for_sun_angle = COS_SUN_PLANE_ANGLE;   // 帆板法线与太阳矢量夹角
            value_sun_angle_to_nomal = ( sun[0]*skyMap[s_id].nx 
                                        + sun[1]*skyMap[s_id].ny 
                                        + sun[2]*skyMap[s_id].nz ) / dist_sun;

            value_sun_angle_to_nomal = fabs(value_sun_angle_to_nomal);
            panel_is_ok = (int)(value_sun_angle_to_nomal >= thres_for_sun_angle);
#endif

            if( panel_is_ok == 1 ) {

                // // 这一步不可以注释掉，因为不是每一步都进行sky_id_tracker的更新，所以当前仍然有可能碰到
                // // 下面条件已经满足了情况!!!
                // if ( skyMap[s_id].flag >= skyMap[s_id].maxCoverNum ) {
                //     local_fail_out[7] += 1.0;
                //     continue;
                // }

                // if ( IsObscureBySun(sun,skyMap[s_id].x,skyMap[s_id].y, skyMap[s_id].z) > 1 ) {
                //     local_fail_out[1] += 1.0;
                //     continue;
                // }
                
                // if ( IsObscureByMoon(moon, skyMap[s_id].x,skyMap[s_id].y, skyMap[s_id].z) > 1 ) {
                //     local_fail_out[2] += 1.0;
                //     continue;
                // }

            //=========================================================================
            //  计算CMG的使用情况
            //  计算从当前指向到下一个指向所需要转过的角度，由此角度可计算CMG功耗以及转动所需的时间。
            //  p1和p2分别表示当前望远镜的指向和下一个指向
                double tAngle, cmg_use, tTime;
                // double ra_old = skyMap[*id].ra;
                // double dec_old = skyMap[*id].dec;
                double ra_old = *cur_ra_tel;
                double dec_old= *cur_dec_tel;
                double ra_new = skyMap[s_id].ra;
                double dec_new = skyMap[s_id].dec;

#if defined(_USE_OLD_ROTATION_)
                Get_RotationAngle_Zhang(ra_old, 90-dec_old, ra_new, 90-dec_new, &tAngle);
#else
                Get_RotationAngle_faster2(ra_old, 90-dec_old, ra_new, 90-dec_new, &tAngle, p_rank);
#endif

                cmg_use = get_cmg_use( tAngle );
                tTime = calculateTransTime( tAngle ); // 计算转动时间

                double tSat[3], ttSat[3], tSatE[3], ttSatE[3];

                if( locateSat1(ttSat, curTime+tTime/86400.0, orbitData, orbitDataNum,Orbit_File_Num) < 0.0 ) {
                    continue;
                }

                CoordinateSpinEquatorial2Ecliptic(ttSat, tSat);
                double angleValue;

            //  判断是否被地球所遮挡
                if( IsObscureByEarth_1( tSat, 
                                        sun, 
                                        skyMap[s_id].x, 
                                        skyMap[s_id].y, 
                                        skyMap[s_id].z, 
                                        &angleValue ) == 1 ){
                    local_fail_out[3] += 1.0;
                    continue;
                }

// =============================================================================
// 确定曝光时间
                double exTime = 0;  // 曝光时间
                if(angleValue > 1) { // 跟暗边的夹角，不需要考虑增加曝光时间（以达到信噪比）
                    if (skyMap[s_id].inDeepFlag < 0) {
                        double p2sun_angle = fabs(sun_lon_lat[0] - skyMap[s_id].ra);
                        p2sun_angle = 360-p2sun_angle*(p2sun_angle>180);
                        exTime = getExposureTime_Zodical(fabs(skyMap[s_id].dec),p2sun_angle);
                    }
                    else{
                        exTime = EXTIME_DEEP;
                    }
                } else {
                    double cos70 = 0.34202; //  跟亮边夹角的余弦值
                    if (skyMap[s_id].inDeepFlag < 0) {
                        if(angleValue < COS_POINT_ZENITH_ANGLE_LIGHT_Max) {// angle with earth light edge is 70 deg
                            if(angleValue >= cos70)
                                exTime = EXTIME + (COS_POINT_ZENITH_ANGLE_LIGHT_Max-angleValue)*extraDeep/(COS_POINT_ZENITH_ANGLE_LIGHT_Max - cos70);
                            else
                                exTime = EXTIME + extraDeep +  (cos70-angleValue)*(100-extraDeep)/(cos70 - COS_POINT_ZENITH_ANGLE_LIGHT_Min);
                        } else
                            exTime = EXTIME;
                    } else {
                        if(angleValue < COS_POINT_ZENITH_ANGLE_LIGHT_Max) {// angle with earth light edge is 70 deg
                            if(angleValue >= cos70)
                                exTime = EXTIME_DEEP + (COS_POINT_ZENITH_ANGLE_LIGHT_Max-angleValue)*extraUDeep/(COS_POINT_ZENITH_ANGLE_LIGHT_Max - cos70);
                            else
                                exTime = EXTIME_DEEP + extraUDeep + (cos70-angleValue)*(50-extraUDeep)/(cos70 - COS_POINT_ZENITH_ANGLE_LIGHT_Min);
                        } else{
                            exTime = EXTIME_DEEP;
                        }
                    }
                }


//  将CMG条件判断提前
#if defined(_USE_CMG_ONE_ORBIT_STATE_)
                if( CMG_one_orbit_state_check(cmg_state, curTime, tTime, exTime, cmg_use, p_rank) != 0 ){
                    local_fail_out[5] += 1.0;
                    continue;
                }
#else
                double cmg_total_use = Test_CMGConstraint( cmg_list, 
                                                           curTime, 
                                                           tTime, 
                                                           exTime, 
                                                           cmg_use );
                if( cmg_total_use > CMG_THRES ){
                    local_fail_out[5] += 1.0;
                    continue;
                }
#endif
// =============================================================================
                double exEndTime = curTime+(tTime+exTime)/86400.0;

                if( locateSat1(ttSatE, exEndTime, orbitData, orbitDataNum,Orbit_File_Num) < 0.0 ){
                    continue;
                }

                double energyT1 = battery_q;
                double e_t1 = 0; // 1 orbit, sun time left
                double e_t2 = 0; // 1 orbit, shadow time left
                double yr = (curTime-2459766.0)/365.25;

// =============================================================================
                if (exEndTime > shadow_st && curTime <= shadow_st) {
                    energyT1 = energyT1
                             + (getPlaneEnergy(value_sun_angle_to_nomal,yr)-POWER_CONSUMPTION)*(shadow_st - curTime)*86400.0
                             - POWER_CONSUMPTION*(exEndTime-shadow_st)*86400.0;
                    e_t1 = 0;
                    e_t2 = shadow_et-exEndTime;
                } else {
                    energyT1 = energyT1
                             + (getPlaneEnergy(value_sun_angle_to_nomal,yr)-POWER_CONSUMPTION)*(exEndTime - curTime)*86400.0;
                    e_t1 = shadow_st - exEndTime;
                    e_t2 = shadow_et - shadow_st;
                }

                // 离开阴影区之前所需要消耗的能量（包括可能的太阳能发电部分）
                double left_1o_energy = -getPlaneEnergy(1.0,yr)*e_t1*86400.0+POWER_CONSUMPTION*(e_t1+e_t2)*86400.0;

                // 防止由于电池能量过低，导致望远镜在阴影区的时候死翘翘
                if (energyT1 - left_1o_energy < BATTERY_LOW || energyT1 < BATTERY_LOW) {
                    local_fail_out[4] += 1.0;
                    continue;
                }
// =============================================================================

                // 从赤道坐标系转换到黄道坐标系
                CoordinateSpinEquatorial2Ecliptic(ttSatE, tSatE);
                double unSPoint[2];
                calculateSatEarthPoint( greenStartTime, 
                                        initalGreenWithRa, 
                                        curTime + (tTime + exTime) / 86400.0, 
                                        ttSatE, 
                                        unSPoint, 
                                        infp );

                if( IsInSAA(unSPoint) == 1 ){
                    local_fail_out[8] += 1.0;
                    continue;
                }

                if( IsObscureByEarth_1( tSatE, 
                                        sun, 
                                        skyMap[s_id].x, 
                                        skyMap[s_id].y, 
                                        skyMap[s_id].z, 
                                        &angleValue ) == 1 ){
                    local_fail_out[3] += 1.0;
                    continue;
                }

            /////////////////
            //  CMG 条件判断
            /////////////////
                // double cmg_times = 0;  // 当前队列对应的时间长度（从一轨的
                // if( cmg_list->head != NULL ){
                //     cmg_times = curTime - cmg_list->head->s_time + (tTime+exTime)/86400.0;
                // }

                // double cmg_temp_value = cmg_list->cmg_total + cmg_use;
                // struct CMG_Node* temp_node = cmg_list->head;

                // // 这一段代码的逻辑暂时没有弄明白
                // while(   cmg_times > ORBIT_TIME
                //       && temp_node != NULL
                //       && temp_node->next != NULL ){
                //     cmg_temp_value = cmg_temp_value - temp_node->cmg_value;
                //     double ct1 = temp_node->s_time;
                //     temp_node = temp_node->next;
                //     double ct2 = temp_node->s_time;
                //     cmg_times = cmg_times-(ct2-ct1);
                // }

                // if( cmg_temp_value > CMG_THRES ){
                //     local_fail_out[5] += 1.0;
                //     continue;
                // }

////////////////////////////////////////////////////////////////////
//              以上都是硬性条件的判断，下面开始进行策略的权重分配
////////////////////////////////////////////////////////////////////
                skyMap[s_id].weight = 10000; // 初始权重
                double weight = skyMap[s_id].weight;

// ================================================================
// 帆板的法线必须满足一定的条件，与太阳的夹角小于等于25度
// 观测尽量寻找在“移动的”区域（被太阳和帆板所限定）
                // double lon1 = sun_lon_lat[0] - 90.;
                // double lon2 = sun_lon_lat[0] + 90.;

                // if ( lon1 < 0. )
                //     lon1 = lon1 + 360.;

                // if ( lon2 >= 360. )
                //     lon2 = lon2 - 360.;

                // double lon1_l = lon1 - 25;
                // double lon1_r = lon1 + 25;
                // if ( lon1_l < 0 )
                //     lon1_l = lon1_l + 360;

                // if ( lon1_r >= 360 )
                //     lon1_r = lon1_r - 360;

                // if (lon1_l < lon1_r) {
                //     if( !(skyMap[s_id].ra >= lon1_l && skyMap[s_id].ra <= lon1_r) )
                //         weight = weight + 500;
                // } else {
                //     if( skyMap[s_id].ra > lon1_r && skyMap[s_id].ra < lon1_l)
                //         weight = weight + 500;
                // }

                // double lon2_l = lon2 - 25;
                // if( lon2_l < 0 ) {
                //     lon2_l = lon2_l + 360;
                // }
                // double lon2_r = lon2 + 25;
                // if ( lon2_r >= 360 ) {
                //     lon2_r = lon2_r - 360;
                // }
                // if ( lon2_l < lon2_r ) {
                //     if ( !(skyMap[s_id].ra >= lon2_l && skyMap[s_id].ra <= lon2_r) ) {
                //         weight = weight + 500;
                //     }
                // } else {
                //     if ( skyMap[s_id].ra > lon2_r && skyMap[s_id].ra < lon2_l ) {
                //         weight = weight + 500;
                //     }
                // }
//              ================================================================

//              根据太阳和帆板法线的夹角来修正权重，尽量观测小角度的天区
                // if ( value_sun_angle_to_nomal >= 0.9961947 ) {
                //     weight = weight + 0;
                // } else if ( value_sun_angle_to_nomal >= 0.9848078 ) {
                //     weight = weight + 1000;
                // } else if ( value_sun_angle_to_nomal >= 0.9659258 ) {
                //     weight = weight + 2000;
                // } else if ( value_sun_angle_to_nomal >= 0.9396926 ) {
                //     weight = weight + 3000;
                // } else if ( value_sun_angle_to_nomal >= 0.9063078 ) {
                //     weight = weight + 4000;
                // }

                if ( value_sun_angle_to_nomal >= 0.9961947 ) {
                    weight = weight + 0;
                } else if ( value_sun_angle_to_nomal >= 0.9848078 ) {
                    weight = weight + 1000*0.5;
                } else if ( value_sun_angle_to_nomal >= 0.9659258 ) {
                    weight = weight + 2000*0.5;
                } else if ( value_sun_angle_to_nomal >= 0.9396926 ) {
                    weight = weight + 3000*0.5;
                } else if ( value_sun_angle_to_nomal >= 0.9063078 ) {
                    weight = weight + 4000*0.5;
                }

                // ===========
                // 相邻优先条件
                // 尽可能使得观测天区连续地被覆盖
                // ===========
#if defined(_ENABLE_CONTINOUS_OBS_)
                if( curTime_yr > CONTINOUS_OBS_START_TIME ){
                    int n_up_id,n_down_id,n_left_id,n_right_id;
                    n_up_id = skyMap[s_id].up_neighbour;
                    n_down_id = skyMap[s_id].down_neighbour;
                    n_left_id = skyMap[s_id].left_neighbour;
                    n_right_id = skyMap[s_id].right_neighbour;

                    int neighbourFlag = 4000;
                    if (skyMap[s_id].flag == 0) { // 优先观测还未被观测过的天区
                        if ( (n_up_id >= 0 && skyMap[n_up_id].flag > 0) || n_up_id == -1 ) {
                            neighbourFlag -= 1000;
                        }

                        if ( (n_down_id >= 0 && skyMap[n_down_id].flag > 0) || n_down_id == -1 ) {
                            neighbourFlag -= 1000;
                        }

                        if ( (n_left_id >= 0 && skyMap[n_left_id].flag > 0) || n_left_id == -1 ) {
                            neighbourFlag -= 1000;
                        }

                        if ( (n_right_id >= 0 && skyMap[n_right_id].flag > 0) || n_right_id == -1 ) {
                            neighbourFlag -= 1000;
                        }
                    }

                    weight = weight + neighbourFlag;
                }
#endif

#if defined(_ENABLE_HIGH_LATITUDE_PRIOR_)
                // if( curTime_yr > 3 && curTime_yr < 8 ) {
                if( curTime_yr <  HIGH_LATITUDE_PRIOR_TIME ) {
                    // weight += MAX_VALUE/3.0*( fabs(skyMap[s_id].gb) <= Galaxy_B_Sec_Low );
                    // weight += MAX_VALUE/4.0*( fabs(skyMap[s_id].dec) <= Ecliptic_Lat_Sec_Low );
                    weight += 1500*( fabs(skyMap[s_id].gb) <= Galaxy_B_Sec_Low );
                    weight += 1500*( fabs(skyMap[s_id].dec) <= Ecliptic_Lat_Sec_Low );
                }
#endif

#if defined(_ENABLE_GRI_WEIGHT_)                    
                // gri 三波段覆盖权重
                int griIds[6],t_i;
                for ( t_i = 0; t_i < 6; t_i++ ) {
                    griIds[t_i]=-1;
                }

                // 先找到当前指向对应的gri三个波段在这个天区（CCD）上的位置
                int tem_n_id = s_id;
                if ( skyMap[tem_n_id].right_neighbour > 0 ) {
                    tem_n_id=skyMap[tem_n_id].right_neighbour;
                    if ( skyMap[tem_n_id].down_neighbour > 0 ) {
                        tem_n_id=skyMap[tem_n_id].down_neighbour;
                        griIds[0]=skyMap[tem_n_id].filter[4];
                        if( skyMap[tem_n_id].right_neighbour > 0 ) {
                            tem_n_id=skyMap[tem_n_id].right_neighbour;
                            griIds[1]=skyMap[tem_n_id].filter[2];
                            if( skyMap[tem_n_id].right_neighbour > 0 ) {
                                tem_n_id=skyMap[tem_n_id].right_neighbour;
                                griIds[2]=skyMap[tem_n_id].filter[3];
                            }
                        }
                    }
                }

                tem_n_id = s_id;
                int f_row = 0;
                while( tem_n_id > 0 && f_row < 5 ) {
                    tem_n_id = skyMap[tem_n_id].down_neighbour;
                    f_row++;
                }

                if ( tem_n_id > 0 ) {
                    tem_n_id=skyMap[tem_n_id].right_neighbour;
                    if( tem_n_id > 0 ) {
                        griIds[3]=skyMap[tem_n_id].filter[3];;
                        tem_n_id=skyMap[tem_n_id].right_neighbour;
                        if( tem_n_id > 0 ) {
                            griIds[4]=skyMap[tem_n_id].filter[2];;
                            tem_n_id=skyMap[tem_n_id].right_neighbour;
                            if( tem_n_id > 0 ) {
                                griIds[5]=skyMap[tem_n_id].filter[4];
                            }
                        }
                    }
                }

                int griFlag = 1;
                for ( t_i = 0; t_i < 6; t_i++ ) {
                    if ( griIds[t_i] != -1 && griIds[t_i] < 2 ) {
                    // 只要有一个没有足够的观测次数
                        griFlag = 0;
                        break;
                    } else if ( griIds[t_i] != -1 && griIds[t_i] >= 2 ) {
                    // 全部波段都观测过
                        griFlag = 1;
                    }
                }

                if ( griFlag == 1 ) {
                    weight += 10000;
                }
#endif

//  8.5年后开启极深场优先,即使需要大角度转动
                // if( curTime_yr > 8.5 && skyMap[s_id].inDeepFlag > 0 ){
                //     weight -= 3000;
                // }

/************************************************************************/
                // 根据CMG转动角度修正权重

                double angles[8] = {5,10,20,35,45,75,90,180};
                int ii, angles_idx = -1;
                for( ii=0; ii<8; ii++ ){
                    if( tAngle <= angles[ii] ){
                        angles_idx = ii;
                        break;
                    }
                }

                switch ( angles_idx ) {
                    case 0:
                        weight = weight/(-5*tAngle+60);
                        if (skyMap[s_id].inDeepFlag > 0){
                            weight = weight - 200*(5-skyMap[s_id].flag);
                            // weight = weight - 100*(5-skyMap[s_id].flag)*( *id == s_id );
                        }
                        
                        if( curTime_yr <= CONTINOUS_OBS_START_TIME )
                            weight = weight - 50.*exp(-(tAngle-5)*(tAngle-5)/5);
                        break;
                    case 1:
                        weight = weight/(-2*tAngle+45);
                        if (skyMap[s_id].inDeepFlag > 0){
                            weight = weight - 150*(5-skyMap[s_id].flag);
                        }

                        if( curTime_yr <= CONTINOUS_OBS_START_TIME )
                            weight = weight - 50.*exp(-(tAngle-5)*(tAngle-5)/5);
                        break;
                    case 2:
                        weight = weight/(-tAngle+35.0);
                        if (skyMap[s_id].inDeepFlag > 0){
                            weight = weight - 100*(5-skyMap[s_id].flag);
                        }
                        break;
                    case 3:
                        weight = weight/(-0.333333*tAngle+21.666667);
                        break;
                    case 4:
                        weight = weight/(-0.2*tAngle+17.0);
                        break;
                    case 5:
                        weight = weight/(-0.1*tAngle+12.5);
                        break;
                    case 6:
                        weight = weight/(-0.0666667*tAngle+10.0);
                        break;
                    case 7:
                        weight = weight/(-0.0222222*tAngle+6.0);
                        break;
                    default:
                        printf("Error happened when correcting weights from CMG rotation angle!\n");
                        break;
                }


                //  ==================================================================
                //  根据望远镜的运动方向来分配权重，尽可能观测与望远镜运动方向一致的区域
                double newPx, newPy, newPz, lastPx, lastPy, lastPz; // pointAngle;
                newPx = skyMap[s_id].x;
                newPy = skyMap[s_id].y;
                newPz = skyMap[s_id].z;

                // lastPx = skyMap[*id].x;
                // lastPy = skyMap[*id].y;
                // lastPz = skyMap[*id].z;

                lastPx = cos(dec_old*PI/180.)*cos(ra_old*PI/180);
                lastPy = cos(dec_old*PI/180.)*sin(ra_old*PI/180);
                lastPz = sin(dec_old*PI/180.);

                // skyMap[s_id].wFactor = 1;

                // tSatE 是下一次曝光开始时望远镜的位置，sat是望远镜在当前时刻的位置

                double dp1 = sqrt( (tSatE[0] - sat[0])*(tSatE[0] - sat[0]) 
                                 + (tSatE[1] - sat[1])*(tSatE[1] - sat[1]) 
                                 + (tSatE[2] - sat[2])*(tSatE[2] - sat[2]));
                double dp2 = sqrt( (newPx - lastPx)*(newPx - lastPx) 
                                 + (newPy - lastPy)*(newPy - lastPy)
                                 + (newPz - lastPz)*(newPz - lastPz) );
                double dp = (tSatE[0] - sat[0]) * (newPx - lastPx)
                          + (tSatE[1] - sat[1]) * (newPy - lastPy)
                          + (tSatE[2] - sat[2]) * (newPz - lastPz);

                dp = dp/dp1/dp2;
                // 如果可以观测同一极深场区域，那么这里的权重就不改变。
                weight -= 1000*dp;

                // double lat_dir = tSatE[2] - sat[2];
                // double p_lat_dir = newPz - lastPz;

                // // 这个判断是什么意思？
                // // dp>0 使得天区选择的方向与卫星运动的方向同向
                // // "lat_dir * p_lat_dir > 0" 使得尽量在南半球看南半球天区，在北半球看北半球天区

                // if ( lat_dir * p_lat_dir > 0 )
                //     weight -= 500;

#if defined(_TURN_ON_5DEG_DRIFT_)
                velocity_tel[0] = (tSatE[0] - sat[0])/dp1;
                velocity_tel[1] = (tSatE[1] - sat[1])/dp1;
                velocity_tel[2] = (tSatE[2] - sat[2])/dp1;
#endif

#if defined(_TURN_DEC60_PRIOR_)
                if( curTime_yr < DEC60_PRIOR_TIME ){
                    // weight -= 500*exp(-(dec_new-60)*(dec_new-60)/15);
                    // weight -= 500*exp(-(dec_new+60)*(dec_new+60)/15);

                    double delta_ra = fmin(fabs(ra_new-0),fabs(ra_new-360));
                    weight -= (500+(10-curTime_yr)*50)*exp(-(dec_new-60)*(dec_new-60)/800)*exp(-(ra_new-180)*(ra_new-180)/1000);
                    // weight -= (500+(10-curTime_yr)*50)*exp(-(dec_new-60)*(dec_new-60)/800)*exp(-(ra_new-150)*(ra_new-150)/1000);
                    weight -= (500+(10-curTime_yr)*50)*exp(-(dec_new+60)*(dec_new+60)/800)*exp(-delta_ra*delta_ra/1000);
                }
#endif

/************************************************************************/

                if ( weight < final_min_value ) {
                    final_min_value = weight;
                    final_index = s_id;
                    final_use_time = tTime + exTime;
                    final_cmg_value = cmg_use;
                    final_exp_time = exTime;
                    final_energy_time = energyT1;
                    final_Angle = tAngle;
                    final_cos_sun_plane_angle = value_sun_angle_to_nomal;
                    final_cos_panel_point_angle = value_panel_point_angle;
                    rank_result.n_rank = p_rank;
                    rank_result.val_weight = weight;

                    // final_left_1o_energy = left_1o_energy;
                    // final_temp_sun = e_t1;
                    // final_temp_shadow = e_t2;
                }

            } else {
                local_fail_out[6] += 1.0;
            }
        }
    } else {
        //  在阴影区时

        double next_shadow_s=0;
        double next_shadow_e=0;
        aquireShadowTime(&next_shadow_s, &next_shadow_e, shadow_et+10/86400.);

        int temp_i = 0;
        for( temp_i=idx_start; temp_i<idx_end; temp_i++ ){

            i = sky_id_tracker[temp_i];

            // 已经测试过，这种情况不会发生
            // if ( skyMap[i].ra > 360 ) {
            //     continue;
            // }

            // if (skyMap[i].isObserve == 0) {
            //     continue;
            // }
                
            if ( skyMap[i].flag >= skyMap[i].maxCoverNum ) {
                local_fail_out[7] += 1.0;
                continue;
            }

            if ( IsObscureBySun(sun,skyMap[i].x,skyMap[i].y, skyMap[i].z) > 1 ) {
                local_fail_out[1] += 1.0;
                continue;
            }
            
            if ( IsObscureByMoon(moon, skyMap[i].x,skyMap[i].y, skyMap[i].z) > 1 ) {
                local_fail_out[2] += 1.0;
                continue;
            }

        //  CMG
            double tAngle, cmg_use, tTime;
            // double ra_old = skyMap[*id].ra;
            // double dec_old = skyMap[*id].dec;
            double ra_old = *cur_ra_tel;
            double dec_old= *cur_dec_tel;
            double ra_new = skyMap[i].ra;
            double dec_new = skyMap[i].dec;
            
#if defined(_USE_OLD_ROTATION_)
            Get_RotationAngle_Zhang(ra_old, 90-dec_old, ra_new, 90-dec_new, &tAngle); 
#else
            Get_RotationAngle_faster2(ra_old, 90-dec_old, ra_new, 90-dec_new, &tAngle, p_rank);
#endif

            cmg_use = get_cmg_use( tAngle );
            tTime = calculateTransTime( tAngle );

            double tSat[3], ttSat[3], tSatE[3], ttSatE[3];

            int flag_sat = locateSat1(ttSat, curTime + tTime/86400.0, orbitData, orbitDataNum,Orbit_File_Num) > 0.0;
            if( flag_sat != 1 ){
                continue;
            }

            CoordinateSpinEquatorial2Ecliptic(ttSat, tSat);

//  ================================================================================
#if defined(_ENABLE_DARKSIDE_PANEL_SUN_ANGLE_CHECK_)
            //  旧版本中，这里不检查望远镜从地影跑出来时帆板的情况。之前检查模拟结果时发现，望远镜
            //  如果在曝光前刚好出了地影，那么会可能出现帆板法线与太阳夹角大于25度的情况。
            int panel_is_ok = -1;
            double cosval_tmp=1.0;
            if ( flag_sat == 1 && IsInSunSide( tSat, sun ) == 1 ) {

#if defined(_ENABLE_PANEL_ROTATION_)
                panel_is_ok = TestPanelAngle( sun,
                                              dist_sun,
                                              skyMap[i].x,
                                              skyMap[i].y,
                                              skyMap[i].z,
                                              skyMap[i].nx,
                                              skyMap[i].ny,
                                              skyMap[i].nz,
                                              &cosval_tmp,
                                              25,
                                              p_rank );
#else
            //  以下是张鑫原来的做法
                double thres_for_sun_angle = COS_SUN_PLANE_ANGLE;   // 帆板法线与太阳矢量夹角
                cosval_tmp = fabs( sun[0]*skyMap[i].nx 
                                 + sun[1]*skyMap[i].ny 
                                 + sun[2]*skyMap[i].nz ) / dist_sun;
                panel_is_ok = (int)(cosval_tmp >= thres_for_sun_angle);
#endif
                
                if( panel_is_ok != 1  ){
                    local_fail_out[6] += 1.0;
                    continue;
                }

            }
#endif
//  ================================================================================
            
            double angleValue;

            if( IsObscureByEarth_1( tSat, 
                                    sun, 
                                    skyMap[i].x, 
                                    skyMap[i].y, 
                                    skyMap[i].z, 
                                    &angleValue ) == 1 ) {
                local_fail_out[3] += 1.0;
                continue;
            }

            double exTime = 0;
            if( angleValue > 1 ) {
                if ( skyMap[i].inDeepFlag < 0 ) {
                    double p2sun_angle = fabs(sun_lon_lat[0] - skyMap[i].ra);
                    p2sun_angle = 360-p2sun_angle*( p2sun_angle > 180 );
                    exTime = getExposureTime_Zodical(fabs(skyMap[i].dec),p2sun_angle);
                } else {
                    exTime = EXTIME_DEEP;
                }
            } else {
                double cos70 = 0.34202;
                if ( skyMap[i].inDeepFlag < 0 ) {
                    if( angleValue < COS_POINT_ZENITH_ANGLE_LIGHT_Max ) {// angle with earth light edge is 70 deg
                        if( angleValue >= cos70 ) {
                            exTime = EXTIME + (COS_POINT_ZENITH_ANGLE_LIGHT_Max-angleValue)*extraDeep/(COS_POINT_ZENITH_ANGLE_LIGHT_Max - cos70);
                        } else {
                            exTime = EXTIME + extraDeep +  (cos70-angleValue)*(100-extraDeep)/(cos70 - COS_POINT_ZENITH_ANGLE_LIGHT_Min);
                        }
                    } else {
                        exTime = EXTIME;
                    }
                } else {
                    if( angleValue < COS_POINT_ZENITH_ANGLE_LIGHT_Max ) {// angle with earth light edge is 70 deg
                        if( angleValue >= cos70 ) {
                            exTime = EXTIME_DEEP + (COS_POINT_ZENITH_ANGLE_LIGHT_Max-angleValue)*extraUDeep/(COS_POINT_ZENITH_ANGLE_LIGHT_Max - cos70);
                        } else {
                            exTime = EXTIME_DEEP + extraUDeep + (cos70-angleValue)*(50-extraUDeep)/(cos70 - COS_POINT_ZENITH_ANGLE_LIGHT_Min);
                        }
                    } else {
                        exTime = EXTIME_DEEP;
                    }
                }
            }

//  =========================================================================
//  将CMG条件判断提前
#if defined(_USE_CMG_ONE_ORBIT_STATE_)
            if( CMG_one_orbit_state_check(cmg_state, curTime, tTime, exTime, cmg_use, p_rank) != 0 ){
                local_fail_out[5] += 1.0;
                continue;
            }
#else
            double cmg_total_use = Test_CMGConstraint( cmg_list, 
                                                       curTime, 
                                                       tTime, 
                                                       exTime, 
                                                       cmg_use );
            if( cmg_total_use > CMG_THRES ){
                local_fail_out[5] += 1.0;
                continue;
            }
#endif
//  =========================================================================

            double exEndTime = curTime+(tTime+exTime)/86400.0;
            flag_sat = locateSat1(ttSatE, exEndTime, orbitData, orbitDataNum,Orbit_File_Num) > 0;
            if( flag_sat != 1 ){
                continue;
            }

            CoordinateSpinEquatorial2Ecliptic(ttSatE, tSatE);

//  ================================================================================
#if defined(_ENABLE_DARKSIDE_PANEL_SUN_ANGLE_CHECK_)
            //  旧版本中，这里不检查望远镜从地影跑出来时帆板的情况。之前检查模拟结果时发现，望远镜
            //  如果在曝光前刚好出了地影，那么会可能出现帆板法线与太阳夹角大于25度的情况。
            if ( flag_sat == 1 && IsInSunSide( tSatE, sun ) == 1 ) {

#if defined(_ENABLE_PANEL_ROTATION_)
                panel_is_ok = TestPanelAngle( sun,
                                              dist_sun,
                                              skyMap[i].x,
                                              skyMap[i].y,
                                              skyMap[i].z,
                                              skyMap[i].nx,
                                              skyMap[i].ny,
                                              skyMap[i].nz,
                                              &cosval_tmp,
                                              25,
                                              p_rank );
#else
            //  以下是张鑫原来的做法
                double thres_for_sun_angle = COS_SUN_PLANE_ANGLE;   // 帆板法线与太阳矢量夹角
                cosval_tmp = fabs( sun[0]*skyMap[i].nx 
                                 + sun[1]*skyMap[i].ny 
                                 + sun[2]*skyMap[i].nz ) / dist_sun;
                panel_is_ok = (int)(cosval_tmp >= thres_for_sun_angle);
#endif
                
                if( panel_is_ok != 1  ){
                    local_fail_out[6] += 1.0;
                    continue;
                }
            
            }
#endif
//  ================================================================================

            double energyT1 = battery_q;
            double e_t1 = 0; // 1 orbit, sun time left
            double e_t2 = 0; // 1 orbit , shadow time left

            double yr = (curTime-2459766.0)/365.25;  // 此处时间以“年”为单位

            if ( exEndTime > shadow_et && curTime <= shadow_et ) {
                energyT1 = energyT1 + (getPlaneEnergy(0.0,yr)-POWER_CONSUMPTION)*(exEndTime-shadow_et)*86400.0-POWER_CONSUMPTION*(shadow_et-curTime)*86400.0;
                e_t1 = next_shadow_s-exEndTime;
                e_t2 = next_shadow_e- next_shadow_s;
            } else {
                energyT1 = energyT1 - (POWER_CONSUMPTION)*(exEndTime - curTime)*86400.0;
                e_t1 = 0;
                e_t2 = shadow_et - exEndTime;
            }

            double left_1o_energy = -getPlaneEnergy(1.0,yr)*e_t1*86400.0+POWER_CONSUMPTION*(e_t1+e_t2)*86400.0;

            if ( energyT1 - left_1o_energy < BATTERY_LOW || energyT1 < BATTERY_LOW ) {
                local_fail_out[4] += 1.0;
                continue;
            }

            // CoordinateSpinEquatorial2Ecliptic(ttSatE, tSatE);
            double unSPoint[2];
            calculateSatEarthPoint( greenStartTime, 
                                    initalGreenWithRa, 
                                    curTime + (tTime + exTime)/86400.0, 
                                    ttSatE, 
                                    unSPoint, 
                                    infp );
            
            if ( IsInSAA(unSPoint) == 1 ) {
                local_fail_out[8] += 1.0;
                continue;
            }

            if( IsObscureByEarth_1( tSatE, 
                                    sun, 
                                    skyMap[i].x, 
                                    skyMap[i].y, 
                                    skyMap[i].z, 
                                    &angleValue) == 1 ) {
                local_fail_out[3] += 1.0;
                continue;
            }

        ////////////////////////////////////////////////////////////////////
        //  CMG 条件判断
        ////////////////////////////////////////////////////////////////////
            // double cmg_times = 0;
            // if ( cmg_list->head != NULL ) {
            //     cmg_times = curTime - cmg_list->head->s_time + (tTime+exTime)/86400.0;
            // }

            // double cmg_temp_value = cmg_list->cmg_total + cmg_use;
            // struct CMG_Node* temp_node = cmg_list->head;
            
            // while( cmg_times > ORBIT_TIME
            //     && temp_node != NULL
            //     && temp_node->next != NULL ) {
            //     cmg_temp_value = cmg_temp_value-temp_node->cmg_value;
            //     double ct1 = temp_node->s_time;
            //     temp_node = temp_node->next;
            //     double ct2 = temp_node->s_time;
            //     cmg_times = cmg_times-(ct2-ct1);
            // }

            // if( cmg_temp_value > CMG_THRES ) {
            //     local_fail_out[5] += 1.0;
            //     continue;
            // }

        ////////////////////////////////////////////////////////////////////
        //  开始分配权重
            skyMap[i].weight=10000;
            double weight = skyMap[i].weight;

#if defined(_ENABLE_CONTINOUS_OBS_)
            // ===========
            // 相邻优先条件
            // ===========
            if( curTime_yr > CONTINOUS_OBS_START_TIME ){
                int n_up_id,n_down_id,n_left_id,n_right_id;
                n_up_id = skyMap[i].up_neighbour;
                n_down_id = skyMap[i].down_neighbour;
                n_left_id = skyMap[i].left_neighbour;
                n_right_id = skyMap[i].right_neighbour;

                int neighbourFlag = 4000;
                if(skyMap[i].flag == 0) { // 优先观测还未被观测过的天区
                    if( (n_up_id>=0 && skyMap[n_up_id].flag > 0) || n_up_id==-1 ) {
                        neighbourFlag -= 1000;
                    }

                    if( (n_down_id>=0 && skyMap[n_down_id].flag > 0) || n_down_id==-1 ) {
                        neighbourFlag -= 1000;
                    }

                    if( (n_left_id>=0 && skyMap[n_left_id].flag > 0) || n_left_id==-1 ) {
                        neighbourFlag -= 1000;
                    }

                    if( (n_right_id>=0 && skyMap[n_right_id].flag > 0) || n_right_id==-1 ) {
                        neighbourFlag -= 1000;
                    }
                }
                weight = weight + neighbourFlag;
            }
#endif

#if defined(_ENABLE_HIGH_LATITUDE_PRIOR_)
            // if( curTime_yr > 3 && curTime_yr < 8 ) {
            if( curTime_yr <  HIGH_LATITUDE_PRIOR_TIME ) {
                // weight += MAX_VALUE/3.0*( fabs(skyMap[i].gb) < Galaxy_B_Sec_Low );
                // weight += MAX_VALUE/4.0*( fabs(skyMap[i].dec) <= Ecliptic_Lat_Sec_Low );
                weight += 1500*( fabs(skyMap[i].gb) <= Galaxy_B_Sec_Low );
                weight += 1500*( fabs(skyMap[i].dec) <= Ecliptic_Lat_Sec_Low );
            }
#endif

#if defined(_ENABLE_GRI_WEIGHT_)  
            int griIds[6],t_i;
            for ( t_i = 0; t_i < 6; t_i++ ) {
                griIds[t_i] = -1;
            }

            int tem_n_id = i;
            if ( skyMap[tem_n_id].right_neighbour > 0 ) {
                tem_n_id=skyMap[tem_n_id].right_neighbour;
                if ( skyMap[tem_n_id].down_neighbour > 0 ) {
                    tem_n_id=skyMap[tem_n_id].down_neighbour;
                    griIds[0]=skyMap[tem_n_id].filter[4];
                    if ( skyMap[tem_n_id].right_neighbour > 0 ) {
                        tem_n_id=skyMap[tem_n_id].right_neighbour;
                        griIds[1]=skyMap[tem_n_id].filter[2];
                        if ( skyMap[tem_n_id].right_neighbour > 0 ) {
                            tem_n_id=skyMap[tem_n_id].right_neighbour;
                            griIds[2]=skyMap[tem_n_id].filter[3];
                        }
                    }
                }
            }

            tem_n_id = i;
            int f_row = 0;
            while ( tem_n_id > 0 && f_row<5 ) {
                tem_n_id = skyMap[tem_n_id].down_neighbour;
                f_row++;
            }

            if ( tem_n_id > 0 ) {
                tem_n_id=skyMap[tem_n_id].right_neighbour;
                if ( tem_n_id > 0 ) {
                    griIds[3]=skyMap[tem_n_id].filter[3];;
                    tem_n_id=skyMap[tem_n_id].right_neighbour;
                    if ( tem_n_id > 0 ) {
                        griIds[4]=skyMap[tem_n_id].filter[2];;
                        tem_n_id=skyMap[tem_n_id].right_neighbour;
                        if( tem_n_id > 0 ) {
                            griIds[5]=skyMap[tem_n_id].filter[4];
                        }
                    }
                }
            }
            int griFlag = 1;
            for ( t_i = 0; t_i < 6; t_i++ ) {
                if (griIds[t_i]!=-1 && griIds[t_i]<2) {
                    griFlag = griFlag*0;
                } else if (griIds[t_i]!=-1 && griIds[t_i]>=2) {
                    griFlag = griFlag*1;
                }
            }
            if ( griFlag == 1 ) {
                weight += 10000;
            }
#endif

//  8.5年后开启极深场优先,即使需要大角度转动
            // if( curTime_yr > 8.5 && skyMap[i].inDeepFlag > 0 ){
            //     weight -= 3000;
            // }

/************************************************************************/
            // 根据CMG转动角度修正权重

            double angles[8] = {5,10,20,35,45,75,90,180};
            int ii, angles_idx = -1;
            for( ii=0; ii<8; ii++ ){
                if( tAngle <= angles[ii] ){
                    angles_idx = ii;
                    break;
                }
            }

            switch ( angles_idx ) {
                case 0:
                    weight = weight/(-5*tAngle+60);
                    // if ( curTime_yr > 2.5 && skyMap[i].inDeepFlag > 0){
                    if ( skyMap[i].inDeepFlag > 0){
                            weight = weight - 200*(5-skyMap[i].flag);
                            // weight = weight - 100*(5-skyMap[i].flag)*(*id == i);
                    }

                    if( curTime_yr <= CONTINOUS_OBS_START_TIME )
                        weight = weight - 250.*exp(-(tAngle-5)*(tAngle-5)/5);
                    break;
                case 1:
                    weight = weight/(-2*tAngle+45);
                    if (skyMap[i].inDeepFlag > 0){
                            weight = weight - 150*(5-skyMap[i].flag);
                    }

                    if( curTime_yr <= CONTINOUS_OBS_START_TIME )
                        weight = weight - 250.*exp(-(tAngle-5)*(tAngle-5)/5);
                    break;
                case 2:
                    weight = weight/(-tAngle+35.0);
                    if (skyMap[i].inDeepFlag > 0){
                            weight = weight - 100*(5-skyMap[i].flag);
                    }
                    break;
                case 3:
                    weight = weight/(-0.333333*tAngle+21.666667);
                    break;
                case 4:
                    weight = weight/(-0.2*tAngle+17.0);
                    break;
                case 5:
                    weight = weight/(-0.1*tAngle+12.5);
                    break;
                case 6:
                    weight = weight/(-0.0666667*tAngle+10.0);
                    break;
                case 7:
                    weight = weight/(-0.0222222*tAngle+6.0);
                    break;
                default:
                    printf("Error happened when correcting weights from CMG rotation angle!\n");
                    break;
            }

        //  引导望远镜观测那些沿着它运动方向的天区
            double newPx, newPy, newPz, lastPx, lastPy, lastPz; // pointAngle;
            newPx = skyMap[i].x;
            newPy = skyMap[i].y;
            newPz = skyMap[i].z;

            // lastPx = skyMap[*id].x;
            // lastPy = skyMap[*id].y;
            // lastPz = skyMap[*id].z;

            lastPx = sin(dec_old*PI/180)*sin(ra_old*PI/180);
            lastPy = sin(dec_old*PI/180)*sin(ra_old*PI/180);
            lastPz = cos(dec_old*PI/180);


            double dp1 = sqrt( (tSatE[0] - sat[0])*(tSatE[0] - sat[0]) 
                             + (tSatE[1] - sat[1])*(tSatE[1] - sat[1]) 
                             + (tSatE[2] - sat[2])*(tSatE[2] - sat[2]));
            double dp2 = sqrt( (newPx - lastPx)*(newPx - lastPx) 
                             + (newPy - lastPy)*(newPy - lastPy)
                             + (newPz - lastPz)*(newPz - lastPz) );

            double dp = (tSatE[0] - sat[0]) * (newPx - lastPx)
                        + (tSatE[1] - sat[1]) * (newPy - lastPy)
                        + (tSatE[2] - sat[2]) * (newPz - lastPz);

            dp = dp/dp1/dp2;
            
            weight -= 1000*dp;

#if defined(_TURN_ON_5DEG_DRIFT_)
            velocity_tel[0] = (tSatE[0] - sat[0])/dp1;
            velocity_tel[1] = (tSatE[1] - sat[1])/dp1;
            velocity_tel[2] = (tSatE[2] - sat[2])/dp1;
#endif

#if defined(_TURN_DEC60_PRIOR_)
            if( curTime_yr < DEC60_PRIOR_TIME ){
                double delta_ra = fmin(fabs(ra_new-0),fabs(ra_new-360));
                weight -= (500+(10-curTime_yr)*50)*exp(-(dec_new-60)*(dec_new-60)/800)*exp(-(ra_new-180)*(ra_new-180)/1000);
                // weight -= (500+(10-curTime_yr)*50)*exp(-(dec_new-60)*(dec_new-60)/800)*exp(-(ra_new-150)*(ra_new-150)/1000);
                weight -= (500+(10-curTime_yr)*50)*exp(-(dec_new+60)*(dec_new+60)/800)*exp(-delta_ra*delta_ra/1000);
            }
#endif
/************************************************************************/

            if ( weight < final_min_value ) {
                final_min_value = weight;
                final_index = i;
                final_use_time = tTime + exTime;
                final_cmg_value = cmg_use;
                final_energy_time = energyT1;
                final_exp_time = exTime;
                final_Angle = tAngle;
                rank_result.n_rank = p_rank;
                rank_result.val_weight = weight;

                // final_left_1o_energy = left_1o_energy;
                // final_temp_sun = e_t1;
                // final_temp_shadow = e_t2;
            }
            
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Allreduce(&rank_result,&final_result,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);

    MPI_Bcast(&final_min_value, 1, MPI_DOUBLE, final_result.n_rank, MPI_COMM_WORLD );
    MPI_Bcast(&final_index, 1, MPI_INT, final_result.n_rank, MPI_COMM_WORLD );
    MPI_Bcast(&final_use_time, 1, MPI_DOUBLE, final_result.n_rank, MPI_COMM_WORLD );
    MPI_Bcast(&final_cmg_value, 1, MPI_DOUBLE, final_result.n_rank, MPI_COMM_WORLD );
    MPI_Bcast(&final_exp_time, 1, MPI_DOUBLE, final_result.n_rank, MPI_COMM_WORLD );
    MPI_Bcast(&final_energy_time, 1, MPI_DOUBLE, final_result.n_rank, MPI_COMM_WORLD );
    MPI_Bcast(&final_Angle, 1, MPI_DOUBLE, final_result.n_rank, MPI_COMM_WORLD );
    MPI_Bcast(&final_cos_sun_plane_angle, 1, MPI_DOUBLE, final_result.n_rank, MPI_COMM_WORLD );
    MPI_Bcast(&final_cos_panel_point_angle, 1, MPI_DOUBLE, final_result.n_rank, MPI_COMM_WORLD );

    // if( p_rank == 0 ){
    //     printf("final_min_value = %20.5f,  MAX_VALUE = %g\n",final_min_value,MAX_VALUE);
    // }


    if ( final_min_value >= MAX_VALUE ) {

        int local_i;
        for (local_i = 1; local_i < 9; local_i ++) {
            MPI_Allreduce(local_fail_out+local_i,outputresult_fail+local_i,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        }

        outputresult_fail[0] = curTime;
        outputresult_fail[9] = isInSunSide;
        outputresult[0] = -1;
        //sun_status=oSunStatus;

        returnTime = curTime + jump_time / 86400.0;

        double e_t1 = 0; // 1 orbit, sun time left
        // double e_t2 = 0;// 1 orbit , shadow time left

        if ( returnTime > shadow_et && curTime < shadow_et ) {
            e_t1 = returnTime - shadow_et;
            // e_t2 = shadow_et - curTime;
        } else if ( returnTime> shadow_st && curTime < shadow_st ) {
            e_t1 = shadow_st - curTime;
            // e_t2 = returnTime - shadow_st;
        } else if ( returnTime <= shadow_st ) {
            e_t1 = returnTime-curTime;
            // e_t2 = 0;
        } else if ( curTime >= shadow_st ) {
            e_t1 = 0;
            // e_t2 = returnTime-curTime;
        }

        double aBatt = battery_q;

        double yr = (curTime-2459766.0)/365.25;

        battery_q = battery_q + e_t1*86400.0*getPlaneEnergy(1.0,yr) - jump_time*POWER_CONSUMPTION;

        if ( battery_q > BATTERY_MAX ) {
            battery_q = BATTERY_MAX;
        }

        if ( battery_q < BATTERY_LOW ) {
            battery_q = BATTERY_LOW;
            printf("Energy not balance, Position 1, battery is %f, before battery is %f, shadow st is %f, shadow et is %f, curTime is %f\n",
                    battery_q,aBatt,shadow_st,shadow_et,curTime);
            MPI_Finalize();
            exit(0);
        }


#if defined(_TURN_ON_5DEG_DRIFT_)
    //  借助这段不观测时间，将望远镜的指向逐步转向天顶或望远镜飞行的方向
    double drift_angle;
    Drift_by_5deg(  cur_ra_tel, 
                    cur_dec_tel, 
                    sat,
                    velocity_tel,
                    &drift_angle,
                    p_rank );
    
#if defined(_USE_CMG_ONE_ORBIT_STATE_)
    CMG_one_orbit_state_update(cmg_state, curTime, 90, 0, 0, p_rank);
#endif

    Write_drift_state(  fp_drift,
                        p_rank,
                        returnTime, *cur_dec_tel, *cur_ra_tel,
                        -999, -999, -999,
                        -999, -999, -999,
                        -999, -999, -999,
                        0, 0, -999,
                        -999,
                        -999,
                        drift_angle,
                        -999,
                        0,
                        battery_q,
                        -999,
                        -999,
                        -999 );
#endif

    } else {
        returnTime = curTime + (final_use_time) / 86400.0;
        skyMap[final_index].flag++;

        //  不要忘记在此处更新望远镜的指向
        *cur_ra_tel = skyMap[final_index].ra;
        *cur_dec_tel= skyMap[final_index].dec;

        double otime = curTime+(final_use_time-final_exp_time)/86400.0; // 下一次曝光开始的时间
        double sun_r1[3], sun_r1_ecl[3];
        double moon_r1[3],sat_r1[3];

        locateSun(infp, otime, sun_r1);
        locateMoon(infp, otime, moon_r1);
        locateSat1(sat_r1, otime, orbitData, orbitDataNum,Orbit_File_Num);

    //  有这样一种可能，在搜寻目标观测天区的时候，望远镜还处在阴影区，但是下一次曝光开始的时候，望远镜已经
    //  跑到了阳照区，于是这里重新判断一下是否在阳照区，并且重新计算帆板法线与太阳的夹角
        int isInSunSide_tmp = IsInSunSide( sat_r1, sun_r1);

        if ( isInSunSide_tmp == 1 ) {
            skyMap[final_index].IsInSunSideFlag = 1;
        
            //  计算帆板法线与太阳夹角时需要将太阳的坐标转换到黄道坐标系
            CoordinateSpinEquatorial2Ecliptic(sun_r1, sun_r1_ecl);
            
            //  final_cos_sun_plane_angle 需要根据搜寻到的天区指向重新计算该余弦值

            double dist_sun = sqrt( sun_r1_ecl[0]*sun_r1_ecl[0] 
                                  + sun_r1_ecl[1]*sun_r1_ecl[1] 
                                  + sun_r1_ecl[2]*sun_r1_ecl[2] );

#if defined(_ENABLE_PANEL_ROTATION_)
        //  这里是我重新实现的算法，与唐怀金写的一个matlab小程序进行了对比，结果完全一致
            TestPanelAngle( sun_r1_ecl,
                            dist_sun,
                            skyMap[final_index].x,
                            skyMap[final_index].y,
                            skyMap[final_index].z,
                            skyMap[final_index].nx,
                            skyMap[final_index].ny,
                            skyMap[final_index].nz,
                            &final_cos_sun_plane_angle,
                            25,
                            p_rank );
            outputresult[14] = acos(final_cos_sun_plane_angle)*180/PI;
#else
        //  以下是张鑫原来的做法
            final_cos_sun_plane_angle = fabs( sun_r1_ecl[0]*skyMap[final_index].nx 
                                            + sun_r1_ecl[1]*skyMap[final_index].ny 
                                            + sun_r1_ecl[2]*skyMap[final_index].nz ) / dist_sun;

            outputresult[14] = acos(final_cos_sun_plane_angle)*180/PI;
            if ( outputresult[14] > 90 ) { // 这种情况不应该出现
                outputresult[14] = 180 - outputresult[14];
            }
#endif
        } else {
            skyMap[final_index].IsInSunSideFlag = 0;
            outputresult[14] = 1000;
        }


        outputresult[0] = skyMap[final_index].id;
        outputresult[1] = otime;    // 下一次曝光开始的时间
        //outputresult[1] = curTime;
        outputresult[2] = sat_r1[0];
        outputresult[3] = sat_r1[1];
        outputresult[4] = sat_r1[2];
        outputresult[5] = sun_r1[0];
        outputresult[6] = sun_r1[1];
        outputresult[7] = sun_r1[2];
        outputresult[8] = moon_r1[0];
        outputresult[9] = moon_r1[1];
        outputresult[10] = moon_r1[2];
        outputresult[11] = final_exp_time;
        outputresult[12] = final_Angle; //  实际上是从上次指向转动到当前指向所需的转动角度！！！！

        // outputresult[13] is NOT used !!!

#if defined(_USE_CMG_ONE_ORBIT_STATE_)
        CMG_one_orbit_state_update( cmg_state, 
                                    curTime, 
                                    final_use_time-final_exp_time, 
                                    final_exp_time, 
                                    final_cmg_value,
                                    p_rank );

        // if( p_rank == 0 ){
        //     printf("\t\t--> cmg_total_time = %.10g, cmg_total = %.10g\n", 
        //                     (cmg_state->time_tail-cmg_state->time_head)*24*60,
        //                     cmg_state->cmg_total);
        // }
#else
        struct CMG_Node* cmg_node = NULL;
        cmg_node = (struct CMG_Node*) malloc(sizeof(struct CMG_Node));

        initCMGNode( cmg_node,
                     otime,
                     final_cmg_value,
                     final_use_time/86400.,
                     NULL );

        addCMG_Node( cmg_list,
                     cmg_node );

        double total_dura_time = cmg_list->time_dura;
        while ( total_dura_time > ORBIT_TIME ) {
            struct CMG_Node* del_node = PopCMGNode(cmg_list);
            if(del_node != NULL) {
                freeCMG_Node(del_node);
            } else {
                cmg_list->time_dura = 0;
            }
            total_dura_time = cmg_list->time_dura ;
        }

        if ( cmg_list->cmg_total > CMG_THRES ) {
            printf("CMG have ill !!!!!!!!!!!!!!!!!!!!");
        }
#endif

        battery_q = final_energy_time;
        if ( battery_q > BATTERY_MAX ) {
            battery_q = BATTERY_MAX;
        }

        if ( battery_q < BATTERY_LOW ) {
            printf("Energy not balance, Position 2, battery is %f\n",battery_q);
            MPI_Finalize();
            exit(0);
        }

        outputresult[15] = saaTime;
    }

    *currentTime = returnTime;
}
