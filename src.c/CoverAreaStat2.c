#include "CoverAreaStat.h"
#include "GSL_funcs.h"


//  统计指向一次所增加的面积
void MarkObervePoint_A( double lat,
                        double lon,
                        double isdeep, 
                        struct FourTree* Trees12,
                        double ccd_pos_in_focus[18][2],
                        double *area1,
                        double *area2,
                        double deltArea ) {
    int i, j;
    // int bandNum = 7;
    // int obsNum1 = 2;
    // int obsNum2 = 8;

    struct searchNodes* r_nodes = (struct searchNodes*)malloc(sizeof(struct searchNodes));
    r_nodes->s_len = 0;
    for(i = 0; i < 12; i ++) {
        findSearchRange(lat, lon, (Trees12+i)->head, 0, (Trees12+i)->layer, r_nodes);
    }

    for( i = 0; i < 4; i++ ) {
        struct Rectangle* ccdpos = (struct Rectangle*)malloc(sizeof(struct Rectangle));
        computCCDPos( lat, lon, ccd_pos_in_focus[i][0], ccd_pos_in_focus[i][1], ccdpos );

        for(j =0; j < r_nodes->s_len; j ++) {
            if(IsInRect(r_nodes->s_nodes[j]->midPoint,ccdpos) == 1) {
                r_nodes->s_nodes[j]->cnt_A++;
                if( r_nodes->s_nodes[j]->cnt_A == 1 ){
                    *area1 = *area1 + deltArea;
                    if( isdeep > 0 )
                        *area2 = *area2 + deltArea;
                }
            }
        }

        free(ccdpos);
    }

    for(i = 4; i < 6; i ++) {
        struct Rectangle* ccdpos = (struct Rectangle*)malloc(sizeof(struct Rectangle));
        computCCDPos(lat, lon,ccd_pos_in_focus[i][0], ccd_pos_in_focus[i][1], ccdpos);

        for(j =0; j < r_nodes->s_len; j ++) {
            if(IsInRect(r_nodes->s_nodes[j]->midPoint,ccdpos) == 1) {
                r_nodes->s_nodes[j]->cnt_A++;
                if( r_nodes->s_nodes[j]->cnt_A == 1 ){
                    *area1 = *area1 + deltArea;
                    if( isdeep > 0 )
                        *area2 = *area2 + deltArea;
                }
            }
        }

        free(ccdpos);
    }

    for(i = 6; i < 8; i ++) {
        struct Rectangle* ccdpos = (struct Rectangle*)malloc(sizeof(struct Rectangle));
        computCCDPos(lat, lon,ccd_pos_in_focus[i][0], ccd_pos_in_focus[i][1], ccdpos);

        for(j =0; j < r_nodes->s_len; j ++) {
            if(IsInRect(r_nodes->s_nodes[j]->midPoint,ccdpos) == 1) {
                r_nodes->s_nodes[j]->cnt_A++;
                if( r_nodes->s_nodes[j]->cnt_A == 1 ){
                    *area1 = *area1 + deltArea;
                    if( isdeep > 0 )
                        *area2 = *area2 + deltArea;
                }
            }
        }

        free(ccdpos);
    }

    for(i = 8; i < 10; i ++) {
        struct Rectangle* ccdpos = (struct Rectangle*)malloc(sizeof(struct Rectangle));
        computCCDPos(lat, lon,ccd_pos_in_focus[i][0], ccd_pos_in_focus[i][1], ccdpos);

        for(j =0; j < r_nodes->s_len; j ++) {
            if(IsInRect(r_nodes->s_nodes[j]->midPoint,ccdpos) == 1) {
                r_nodes->s_nodes[j]->cnt_A++;
                if( r_nodes->s_nodes[j]->cnt_A == 1 ){
                    *area1 = *area1 + deltArea;
                    if( isdeep > 0 )
                        *area2 = *area2 + deltArea;
                }
            }
        }

        free(ccdpos);
    }

    for(i = 10; i < 12; i ++) {
        struct Rectangle* ccdpos = (struct Rectangle*)malloc(sizeof(struct Rectangle));
        computCCDPos(lat, lon,ccd_pos_in_focus[i][0], ccd_pos_in_focus[i][1], ccdpos);

        for(j =0; j < r_nodes->s_len; j ++) {
            if(IsInRect(r_nodes->s_nodes[j]->midPoint,ccdpos) == 1) {
                r_nodes->s_nodes[j]->cnt_A++;
                if( r_nodes->s_nodes[j]->cnt_A == 1 ){
                    *area1 = *area1 + deltArea;
                    if( isdeep > 0 )
                        *area2 = *area2 + deltArea;
                }
            }
        }

        free(ccdpos);
    }

    for(i = 12; i < 14; i ++) {
        struct Rectangle* ccdpos = (struct Rectangle*)malloc(sizeof(struct Rectangle));
        computCCDPos(lat, lon,ccd_pos_in_focus[i][0], ccd_pos_in_focus[i][1], ccdpos);

        for(j =0; j < r_nodes->s_len; j ++) {
            if(IsInRect(r_nodes->s_nodes[j]->midPoint,ccdpos) == 1) {
                r_nodes->s_nodes[j]->cnt_A++;
                if( r_nodes->s_nodes[j]->cnt_A == 1 ){
                    *area1 = *area1 + deltArea;
                    if( isdeep > 0 )
                        *area2 = *area2 + deltArea;
                }
            }
        }

        free(ccdpos);
    }
    
    for(i = 14; i < 18; i ++) {
        struct Rectangle* ccdpos = (struct Rectangle*)malloc(sizeof(struct Rectangle));
        computCCDPos(lat, lon,ccd_pos_in_focus[i][0], ccd_pos_in_focus[i][1], ccdpos);
        
        for(j =0; j < r_nodes->s_len; j ++) {
            if(IsInRect(r_nodes->s_nodes[j]->midPoint,ccdpos) == 1) {
                r_nodes->s_nodes[j]->cnt_A++;
                if( r_nodes->s_nodes[j]->cnt_A == 1 ){
                    *area1 = *area1 + deltArea;
                    if( isdeep > 0 )
                        *area2 = *area2 + deltArea;
                }
            }
        }

        free(ccdpos);
    }

    free(r_nodes);
}

//  统计7个波段都覆盖一次以上所增加的面积
void MarkObervePoint_B( double lat,
                        double lon,
                        double isdeep, 
                        struct FourTree* Trees12,
                        double ccd_pos_in_focus[18][2],
                        double *area1,
                        double *area2,
                        double deltArea ) {
    int i, j;
    // int bandNum = 7;
    // int obsNum1 = 2;
    // int obsNum2 = 8;

    struct searchNodes* r_nodes = (struct searchNodes*)malloc(sizeof(struct searchNodes));
    r_nodes->s_len = 0;
    for(i = 0; i < 12; i ++) {
        findSearchRange(lat, lon, (Trees12+i)->head, 0, (Trees12+i)->layer, r_nodes);
    }

    for( i = 0; i < 4; i++ ) {
        struct Rectangle* ccdpos = (struct Rectangle*)malloc(sizeof(struct Rectangle));
        computCCDPos( lat, lon, ccd_pos_in_focus[i][0], ccd_pos_in_focus[i][1], ccdpos );

        for(j =0; j < r_nodes->s_len; j ++) {
            if(IsInRect(r_nodes->s_nodes[j]->midPoint,ccdpos) == 1) {
                r_nodes->s_nodes[j]->band[0]++;
                
                if(    r_nodes->s_nodes[j]->band[0]>=1
                    && r_nodes->s_nodes[j]->band[1]>=1
                    && r_nodes->s_nodes[j]->band[2]>=1
                    && r_nodes->s_nodes[j]->band[3]>=1
                    && r_nodes->s_nodes[j]->band[4]>=1
                    && r_nodes->s_nodes[j]->band[5]>=1
                    && r_nodes->s_nodes[j]->band[6]>=1 ) {
                    // *area1 = *area1 + deltArea;
                    r_nodes->s_nodes[j]->cnt_B++;
                }

                if( r_nodes->s_nodes[j]->cnt_B == 1 ){
                    *area1 = *area1 + deltArea;
                    if( isdeep > 0 )
                        *area2 = *area2 + deltArea;
                }
            }
        }

        free(ccdpos);
    }

    for(i = 4; i < 6; i ++) {
        struct Rectangle* ccdpos = (struct Rectangle*)malloc(sizeof(struct Rectangle));
        computCCDPos(lat, lon,ccd_pos_in_focus[i][0], ccd_pos_in_focus[i][1], ccdpos);

        for(j =0; j < r_nodes->s_len; j ++) {
            if(IsInRect(r_nodes->s_nodes[j]->midPoint,ccdpos) == 1) {
                r_nodes->s_nodes[j]->band[1]++;

                if(    r_nodes->s_nodes[j]->band[0]>=1
                    && r_nodes->s_nodes[j]->band[1]>=1
                    && r_nodes->s_nodes[j]->band[2]>=1
                    && r_nodes->s_nodes[j]->band[3]>=1
                    && r_nodes->s_nodes[j]->band[4]>=1
                    && r_nodes->s_nodes[j]->band[5]>=1
                    && r_nodes->s_nodes[j]->band[6]>=1 ) {
                    // *area1 = *area1 + deltArea;
                    r_nodes->s_nodes[j]->cnt_B++;
                }

                if( r_nodes->s_nodes[j]->cnt_B == 1 ){
                    *area1 = *area1 + deltArea;
                    if( isdeep > 0 )
                        *area2 = *area2 + deltArea;
                }
            }
        }

        free(ccdpos);
    }

    for(i = 6; i < 8; i ++) {
        struct Rectangle* ccdpos = (struct Rectangle*)malloc(sizeof(struct Rectangle));
        computCCDPos(lat, lon,ccd_pos_in_focus[i][0], ccd_pos_in_focus[i][1], ccdpos);

        for(j =0; j < r_nodes->s_len; j ++) {
            if(IsInRect(r_nodes->s_nodes[j]->midPoint,ccdpos) == 1) {
                r_nodes->s_nodes[j]->band[2]++;
                
                if(    r_nodes->s_nodes[j]->band[0]>=1
                    && r_nodes->s_nodes[j]->band[1]>=1
                    && r_nodes->s_nodes[j]->band[2]>=1
                    && r_nodes->s_nodes[j]->band[3]>=1
                    && r_nodes->s_nodes[j]->band[4]>=1
                    && r_nodes->s_nodes[j]->band[5]>=1
                    && r_nodes->s_nodes[j]->band[6]>=1 ) {
                    // *area1 = *area1 + deltArea;
                    r_nodes->s_nodes[j]->cnt_B++;
                }

                if( r_nodes->s_nodes[j]->cnt_B == 1 ){
                    *area1 = *area1 + deltArea;
                    if( isdeep > 0 )
                        *area2 = *area2 + deltArea;
                }
            }
        }

        free(ccdpos);
    }

    for(i = 8; i < 10; i ++) {
        struct Rectangle* ccdpos = (struct Rectangle*)malloc(sizeof(struct Rectangle));
        computCCDPos(lat, lon,ccd_pos_in_focus[i][0], ccd_pos_in_focus[i][1], ccdpos);

        for(j =0; j < r_nodes->s_len; j ++) {
            if(IsInRect(r_nodes->s_nodes[j]->midPoint,ccdpos) == 1) {
                r_nodes->s_nodes[j]->band[3]++;
                
                if(    r_nodes->s_nodes[j]->band[0]>=1
                    && r_nodes->s_nodes[j]->band[1]>=1
                    && r_nodes->s_nodes[j]->band[2]>=1
                    && r_nodes->s_nodes[j]->band[3]>=1
                    && r_nodes->s_nodes[j]->band[4]>=1
                    && r_nodes->s_nodes[j]->band[5]>=1
                    && r_nodes->s_nodes[j]->band[6]>=1 ) {
                    // *area1 = *area1 + deltArea;
                    r_nodes->s_nodes[j]->cnt_B++;
                }

                if( r_nodes->s_nodes[j]->cnt_B == 1 ){
                    *area1 = *area1 + deltArea;
                    if( isdeep > 0 )
                        *area2 = *area2 + deltArea;
                }
            }
        }

        free(ccdpos);
    }

    for(i = 10; i < 12; i ++) {
        struct Rectangle* ccdpos = (struct Rectangle*)malloc(sizeof(struct Rectangle));
        computCCDPos(lat, lon,ccd_pos_in_focus[i][0], ccd_pos_in_focus[i][1], ccdpos);

        for(j =0; j < r_nodes->s_len; j ++) {
            if(IsInRect(r_nodes->s_nodes[j]->midPoint,ccdpos) == 1) {
                r_nodes->s_nodes[j]->band[4]++;
                
                if(    r_nodes->s_nodes[j]->band[0]>=1
                    && r_nodes->s_nodes[j]->band[1]>=1
                    && r_nodes->s_nodes[j]->band[2]>=1
                    && r_nodes->s_nodes[j]->band[3]>=1
                    && r_nodes->s_nodes[j]->band[4]>=1
                    && r_nodes->s_nodes[j]->band[5]>=1
                    && r_nodes->s_nodes[j]->band[6]>=1 ) {
                    // *area1 = *area1 + deltArea;
                    r_nodes->s_nodes[j]->cnt_B++;
                }

                if( r_nodes->s_nodes[j]->cnt_B == 1 ){
                    *area1 = *area1 + deltArea;
                    if( isdeep > 0 )
                        *area2 = *area2 + deltArea;
                }
            }
        }

        free(ccdpos);
    }

    for(i = 12; i < 14; i ++) {
        struct Rectangle* ccdpos = (struct Rectangle*)malloc(sizeof(struct Rectangle));
        computCCDPos(lat, lon,ccd_pos_in_focus[i][0], ccd_pos_in_focus[i][1], ccdpos);

        for(j =0; j < r_nodes->s_len; j ++) {
            if(IsInRect(r_nodes->s_nodes[j]->midPoint,ccdpos) == 1) {
                r_nodes->s_nodes[j]->band[5]++;
                
                if(    r_nodes->s_nodes[j]->band[0]>=1
                    && r_nodes->s_nodes[j]->band[1]>=1
                    && r_nodes->s_nodes[j]->band[2]>=1
                    && r_nodes->s_nodes[j]->band[3]>=1
                    && r_nodes->s_nodes[j]->band[4]>=1
                    && r_nodes->s_nodes[j]->band[5]>=1
                    && r_nodes->s_nodes[j]->band[6]>=1 ) {
                    // *area1 = *area1 + deltArea;
                    r_nodes->s_nodes[j]->cnt_B++;
                }

                if( r_nodes->s_nodes[j]->cnt_B == 1 ){
                    *area1 = *area1 + deltArea;
                    if( isdeep > 0 )
                        *area2 = *area2 + deltArea;
                }
            }
        }

        free(ccdpos);
    }
    
    for(i = 14; i < 18; i ++) {
        struct Rectangle* ccdpos = (struct Rectangle*)malloc(sizeof(struct Rectangle));
        computCCDPos(lat, lon,ccd_pos_in_focus[i][0], ccd_pos_in_focus[i][1], ccdpos);
        
        for(j =0; j < r_nodes->s_len; j ++) {
            if(IsInRect(r_nodes->s_nodes[j]->midPoint,ccdpos) == 1) {
                r_nodes->s_nodes[j]->band[6]++;
                
                if(    r_nodes->s_nodes[j]->band[0]>=1
                    && r_nodes->s_nodes[j]->band[1]>=1
                    && r_nodes->s_nodes[j]->band[2]>=1
                    && r_nodes->s_nodes[j]->band[3]>=1
                    && r_nodes->s_nodes[j]->band[4]>=1
                    && r_nodes->s_nodes[j]->band[5]>=1
                    && r_nodes->s_nodes[j]->band[6]>=1 ) {
                    // *area1 = *area1 + deltArea;
                    r_nodes->s_nodes[j]->cnt_B++;
                }

                if( r_nodes->s_nodes[j]->cnt_B == 1 ){
                    *area1 = *area1 + deltArea;

                    if( isdeep > 0 )
                        *area2 = *area2 + deltArea;
                }
            }
        }

        free(ccdpos);
    }

    free(r_nodes);
}

//  ========================================
//  统计所有7个波段都覆盖两次以上所增加的面积
//  (与原先的面积统计函数没有任何的区别)
void MarkObervePoint_C( double lat,
                        double lon,
                        double isdeep, 
                        struct FourTree* Trees12,
                        double ccd_pos_in_focus[18][2],
                        double *area1,
                        double *area2,
                        double deltArea ) {
    int i, j;
    // int bandNum = 7;
    int obsNum1 = 2;
    int obsNum2 = 8;

    struct searchNodes* r_nodes = (struct searchNodes*)malloc(sizeof(struct searchNodes));
    r_nodes->s_len = 0;
    for(i = 0; i < 12; i ++) {
        findSearchRange(lat, lon, (Trees12+i)->head, 0, (Trees12+i)->layer, r_nodes);
    }

    for( i = 0; i < 4; i++ ) {
        struct Rectangle* ccdpos = (struct Rectangle*)malloc(sizeof(struct Rectangle));
        computCCDPos( lat, lon, ccd_pos_in_focus[i][0], ccd_pos_in_focus[i][1], ccdpos );

        for(j =0; j < r_nodes->s_len; j ++) {
            if(IsInRect(r_nodes->s_nodes[j]->midPoint,ccdpos) == 1) {
                r_nodes->s_nodes[j]->band[0]++;
                if(r_nodes->s_nodes[j]->band[0]==obsNum1) {
                    if(    r_nodes->s_nodes[j]->band[1]>=obsNum1
                        && r_nodes->s_nodes[j]->band[2]>=obsNum1
                        && r_nodes->s_nodes[j]->band[3]>=obsNum1
                        && r_nodes->s_nodes[j]->band[4]>=obsNum1
                        && r_nodes->s_nodes[j]->band[5]>=obsNum1
                        && r_nodes->s_nodes[j]->band[6]>=obsNum1 ) {
                        *area1 = *area1 + deltArea;
                    }
                }
                if(r_nodes->s_nodes[j]->band[0]==obsNum2) {
                    if(    r_nodes->s_nodes[j]->band[1]>=obsNum2
                        && r_nodes->s_nodes[j]->band[2]>=obsNum2
                        && r_nodes->s_nodes[j]->band[3]>=obsNum2
                        && r_nodes->s_nodes[j]->band[4]>=obsNum2
                        && r_nodes->s_nodes[j]->band[5]>=obsNum2
                        && r_nodes->s_nodes[j]->band[6]>=obsNum2 ) {
                        *area2 = *area2 + deltArea;
                    }
                }
            }
        }

        free(ccdpos);
    }

    for(i = 4; i < 6; i ++) {
        struct Rectangle* ccdpos = (struct Rectangle*)malloc(sizeof(struct Rectangle));
        computCCDPos(lat, lon,ccd_pos_in_focus[i][0], ccd_pos_in_focus[i][1], ccdpos);

        for(j =0; j < r_nodes->s_len; j ++) {
            if(IsInRect(r_nodes->s_nodes[j]->midPoint,ccdpos) == 1) {
                r_nodes->s_nodes[j]->band[1]++;
                if(r_nodes->s_nodes[j]->band[1]==obsNum1) {
                    if(    r_nodes->s_nodes[j]->band[0]>=obsNum1
                        && r_nodes->s_nodes[j]->band[2]>=obsNum1
                        && r_nodes->s_nodes[j]->band[3]>=obsNum1
                        && r_nodes->s_nodes[j]->band[4]>=obsNum1
                        && r_nodes->s_nodes[j]->band[5]>=obsNum1
                        && r_nodes->s_nodes[j]->band[6]>=obsNum1 ) {
                        *area1 = *area1 + deltArea;
                    }
                }
                if(r_nodes->s_nodes[j]->band[1]==obsNum2) {
                    if(    r_nodes->s_nodes[j]->band[0]>=obsNum2
                        && r_nodes->s_nodes[j]->band[2]>=obsNum2
                        && r_nodes->s_nodes[j]->band[3]>=obsNum2
                        && r_nodes->s_nodes[j]->band[4]>=obsNum2
                        && r_nodes->s_nodes[j]->band[5]>=obsNum2
                        && r_nodes->s_nodes[j]->band[6]>=obsNum2 ) {
                        *area2 = *area2 + deltArea;
                    }
                }
            }
        }

        free(ccdpos);
    }

    for(i = 6; i < 8; i ++) {
        struct Rectangle* ccdpos = (struct Rectangle*)malloc(sizeof(struct Rectangle));
        computCCDPos(lat, lon,ccd_pos_in_focus[i][0], ccd_pos_in_focus[i][1], ccdpos);

        for(j =0; j < r_nodes->s_len; j ++) {
            if(IsInRect(r_nodes->s_nodes[j]->midPoint,ccdpos) == 1) {
                r_nodes->s_nodes[j]->band[2]++;
                if(r_nodes->s_nodes[j]->band[2]==obsNum1) {
                    if(    r_nodes->s_nodes[j]->band[0]>=obsNum1
                        && r_nodes->s_nodes[j]->band[1]>=obsNum1
                        && r_nodes->s_nodes[j]->band[3]>=obsNum1
                        && r_nodes->s_nodes[j]->band[4]>=obsNum1
                        && r_nodes->s_nodes[j]->band[5]>=obsNum1
                        && r_nodes->s_nodes[j]->band[6]>=obsNum1 ) {
                        *area1 = *area1 + deltArea;
                    }
                }
                if(r_nodes->s_nodes[j]->band[2]==obsNum2) {
                    if(    r_nodes->s_nodes[j]->band[0]>=obsNum2
                        && r_nodes->s_nodes[j]->band[1]>=obsNum2
                        && r_nodes->s_nodes[j]->band[3]>=obsNum2
                        && r_nodes->s_nodes[j]->band[4]>=obsNum2
                        && r_nodes->s_nodes[j]->band[5]>=obsNum2
                        && r_nodes->s_nodes[j]->band[6]>=obsNum2 ) {
                        *area2 = *area2 + deltArea;
                    }
                }
            }
        }

        free(ccdpos);
    }

    for(i = 8; i < 10; i ++) {
        struct Rectangle* ccdpos = (struct Rectangle*)malloc(sizeof(struct Rectangle));
        computCCDPos(lat, lon,ccd_pos_in_focus[i][0], ccd_pos_in_focus[i][1], ccdpos);

        for(j =0; j < r_nodes->s_len; j ++) {
            if(IsInRect(r_nodes->s_nodes[j]->midPoint,ccdpos) == 1) {
                r_nodes->s_nodes[j]->band[3]++;
                if(r_nodes->s_nodes[j]->band[3]==obsNum1) {
                    if(    r_nodes->s_nodes[j]->band[0]>=obsNum1
                        && r_nodes->s_nodes[j]->band[1]>=obsNum1
                        && r_nodes->s_nodes[j]->band[2]>=obsNum1
                        && r_nodes->s_nodes[j]->band[4]>=obsNum1
                        && r_nodes->s_nodes[j]->band[5]>=obsNum1
                        && r_nodes->s_nodes[j]->band[6]>=obsNum1 ) {
                        *area1 = *area1 + deltArea;
                    }
                }
                if(r_nodes->s_nodes[j]->band[3]==obsNum2) {
                    if(    r_nodes->s_nodes[j]->band[0]>=obsNum2
                        && r_nodes->s_nodes[j]->band[1]>=obsNum2
                        && r_nodes->s_nodes[j]->band[2]>=obsNum2
                        && r_nodes->s_nodes[j]->band[4]>=obsNum2
                        && r_nodes->s_nodes[j]->band[5]>=obsNum2
                        && r_nodes->s_nodes[j]->band[6]>=obsNum2 ) {
                        *area2 = *area2 + deltArea;
                    }
                }
            }
        }

        free(ccdpos);
    }

    for(i = 10; i < 12; i ++) {
        struct Rectangle* ccdpos = (struct Rectangle*)malloc(sizeof(struct Rectangle));
        computCCDPos(lat, lon,ccd_pos_in_focus[i][0], ccd_pos_in_focus[i][1], ccdpos);

        for(j =0; j < r_nodes->s_len; j ++) {
            if(IsInRect(r_nodes->s_nodes[j]->midPoint,ccdpos) == 1) {
                r_nodes->s_nodes[j]->band[4]++;
                if(r_nodes->s_nodes[j]->band[4]==obsNum1) {
                    if(    r_nodes->s_nodes[j]->band[0]>=obsNum1
                        && r_nodes->s_nodes[j]->band[1]>=obsNum1
                        && r_nodes->s_nodes[j]->band[2]>=obsNum1
                        && r_nodes->s_nodes[j]->band[3]>=obsNum1
                        && r_nodes->s_nodes[j]->band[5]>=obsNum1
                        && r_nodes->s_nodes[j]->band[6]>=obsNum1 ) {
                        *area1 = *area1 + deltArea;
                    }
                }
                if(r_nodes->s_nodes[j]->band[4]==obsNum2) {
                    if(    r_nodes->s_nodes[j]->band[0]>=obsNum2
                        && r_nodes->s_nodes[j]->band[1]>=obsNum2
                        && r_nodes->s_nodes[j]->band[2]>=obsNum2
                        && r_nodes->s_nodes[j]->band[3]>=obsNum2
                        && r_nodes->s_nodes[j]->band[5]>=obsNum2
                        && r_nodes->s_nodes[j]->band[6]>=obsNum2 ) {
                        *area2 = *area2 + deltArea;
                    }
                }
            }
        }

        free(ccdpos);
    }

    for(i = 12; i < 14; i ++) {
        struct Rectangle* ccdpos = (struct Rectangle*)malloc(sizeof(struct Rectangle));
        computCCDPos(lat, lon,ccd_pos_in_focus[i][0], ccd_pos_in_focus[i][1], ccdpos);

        for(j =0; j < r_nodes->s_len; j ++) {
            if(IsInRect(r_nodes->s_nodes[j]->midPoint,ccdpos) == 1) {
                r_nodes->s_nodes[j]->band[5]++;
                if(r_nodes->s_nodes[j]->band[5]==obsNum1) {
                    if(    r_nodes->s_nodes[j]->band[0]>=obsNum1
                        && r_nodes->s_nodes[j]->band[1]>=obsNum1
                        && r_nodes->s_nodes[j]->band[2]>=obsNum1
                        && r_nodes->s_nodes[j]->band[3]>=obsNum1
                        && r_nodes->s_nodes[j]->band[4]>=obsNum1
                        && r_nodes->s_nodes[j]->band[6]>=obsNum1 ) {
                        *area1 = *area1 + deltArea;
                    }
                }
                if(r_nodes->s_nodes[j]->band[5]==obsNum2) {
                    if(    r_nodes->s_nodes[j]->band[0]>=obsNum2
                        && r_nodes->s_nodes[j]->band[1]>=obsNum2
                        && r_nodes->s_nodes[j]->band[2]>=obsNum2
                        && r_nodes->s_nodes[j]->band[3]>=obsNum2
                        && r_nodes->s_nodes[j]->band[4]>=obsNum2
                        && r_nodes->s_nodes[j]->band[6]>=obsNum2 ) {
                        *area2 = *area2 + deltArea;
                    }
                }
            }
        }

        free(ccdpos);
    }
    
    for(i = 14; i < 18; i ++) {
        struct Rectangle* ccdpos = (struct Rectangle*)malloc(sizeof(struct Rectangle));
        computCCDPos(lat, lon,ccd_pos_in_focus[i][0], ccd_pos_in_focus[i][1], ccdpos);
        
        for(j =0; j < r_nodes->s_len; j ++) {
            if(IsInRect(r_nodes->s_nodes[j]->midPoint,ccdpos) == 1) {
                r_nodes->s_nodes[j]->band[6]++;
                if(r_nodes->s_nodes[j]->band[6]==obsNum1) {
                    if(    r_nodes->s_nodes[j]->band[0]>=obsNum1
                        && r_nodes->s_nodes[j]->band[1]>=obsNum1
                        && r_nodes->s_nodes[j]->band[2]>=obsNum1
                        && r_nodes->s_nodes[j]->band[3]>=obsNum1
                        && r_nodes->s_nodes[j]->band[4]>=obsNum1
                        && r_nodes->s_nodes[j]->band[5]>=obsNum1 ) {
                        *area1 = *area1 + deltArea;
                    }
                }
                if(r_nodes->s_nodes[j]->band[6]==obsNum2) {
                    if(    r_nodes->s_nodes[j]->band[0]>=obsNum2
                        && r_nodes->s_nodes[j]->band[1]>=obsNum2
                        && r_nodes->s_nodes[j]->band[2]>=obsNum2
                        && r_nodes->s_nodes[j]->band[3]>=obsNum2
                        && r_nodes->s_nodes[j]->band[4]>=obsNum2
                        && r_nodes->s_nodes[j]->band[5]>=obsNum2 ) {
                        *area2 = *area2 + deltArea;
                    }
                }
            }
        }

        free(ccdpos);
    }

    free(r_nodes);
}
