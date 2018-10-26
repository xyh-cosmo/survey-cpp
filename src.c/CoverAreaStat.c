#include "CoverAreaStat.h"
#include "GSL_funcs.h"

double getDistPP(double x1, double y1, double z1, double x2, double y2, double z2) {
    double dx = x1-x2;
    double dy = y1-y2;
    double dz = z1-z2;
    return sqrt(dx*dx + dy*dy + dz*dz);
}

double getAngleInSphere(double x1, double y1, double z1, double x2, double y2, double z2) {
    double tt = sqrt((x1*x1 + y1*y1 + z1* z1) * (x2*x2 + y2*y2 + z2*z2));
#if defined(_FULL_DEBUG_)
    if(tt==0) {
        return -1;
    }
#endif
    double cosValue = (x1*x2+y1*y2+z1*z2)/tt;

    if (cosValue > 1) {
        cosValue = 1;
    }
    if (cosValue < -1) {
        cosValue = -1;
    }
    double angle = acos(cosValue);
    return angle;
}

void findmidArc( double x1,double y1,double z1,
                 double x2,double y2,double z2, 
                 double* midP ) {

    double a_sphere = getAngleInSphere(x1,y1,z1,x2,y2,z2);
    double ratio = 0;
    if (fabs(a_sphere-M_PI)>0.00001) {
        ratio = 1/cos(a_sphere*0.5);
    }

    double x_m = 0.5*(x1+x2);
    double y_m = 0.5*(y1+y2);
    double z_m = 0.5*(z1+z2);

    midP[0] = x_m*ratio;
    midP[1] = y_m*ratio;
    midP[2] = z_m*ratio;
}

void getCenterRect(struct Rectangle rect,double* center) {
    double mid1[3],mid2[3],mid3[3],mid4[3];
    double mid13[3],mid24[3];

    findmidArc(rect.p1[0],rect.p1[1],rect.p1[2],rect.p2[0],rect.p2[1],rect.p2[2],mid1);
    findmidArc(rect.p2[0],rect.p2[1],rect.p2[2],rect.p3[0],rect.p3[1],rect.p3[2],mid2);
    findmidArc(rect.p3[0],rect.p3[1],rect.p3[2],rect.p4[0],rect.p4[1],rect.p4[2],mid3);
    findmidArc(rect.p4[0],rect.p4[1],rect.p4[2],rect.p1[0],rect.p1[1],rect.p1[2],mid4);

    findmidArc(mid1[0],mid1[1],mid1[2],mid3[0],mid3[1],mid3[2],mid13);
    findmidArc(mid2[0],mid2[1],mid2[2],mid4[0],mid4[1],mid4[2],mid24);
    center[0] = 0.5*(mid13[0]+mid24[0]);
    center[1] = 0.5*(mid13[1]+mid24[1]);
    center[2] = 0.5*(mid13[2]+mid24[2]);
}

// 将一个矩形划分为2×2的4个小矩形
void devideRect( struct Rectangle rect, 
                 struct Rectangle *rect1, 
                 struct Rectangle *rect2, 
                 struct Rectangle *rect3, 
                 struct Rectangle *rect4) {

    double mid1[3],mid2[3],mid3[3],mid4[3];
    double mid13[3],mid24[3];
    double mid[3];

    findmidArc(rect.p1[0],rect.p1[1],rect.p1[2],rect.p2[0],rect.p2[1],rect.p2[2],mid1);
    findmidArc(rect.p2[0],rect.p2[1],rect.p2[2],rect.p3[0],rect.p3[1],rect.p3[2],mid2);
    findmidArc(rect.p3[0],rect.p3[1],rect.p3[2],rect.p4[0],rect.p4[1],rect.p4[2],mid3);
    findmidArc(rect.p4[0],rect.p4[1],rect.p4[2],rect.p1[0],rect.p1[1],rect.p1[2],mid4);

    findmidArc(mid1[0],mid1[1],mid1[2],mid3[0],mid3[1],mid3[2],mid13);
    findmidArc(mid2[0],mid2[1],mid2[2],mid4[0],mid4[1],mid4[2],mid24);

    mid[0] = 0.5*(mid13[0]+mid24[0]);
    mid[1] = 0.5*(mid13[1]+mid24[1]);
    mid[2] = 0.5*(mid13[2]+mid24[2]);

///////////////////////////////////////////////////
//  这里有个疑问，每一个矩形的1234个角的顺序是什么样子的？
///////////////////////////////////////////////////

//  first
    rect1->p1[0] = rect.p1[0];
    rect1->p1[1] = rect.p1[1];
    rect1->p1[2] = rect.p1[2];

    rect1->p2[0] = mid1[0];
    rect1->p2[1] = mid1[1];
    rect1->p2[2] = mid1[2];

    rect1->p3[0] = mid[0];
    rect1->p3[1] = mid[1];
    rect1->p3[2] = mid[2];

    rect1->p4[0] = mid4[0];
    rect1->p4[1] = mid4[1];
    rect1->p4[2] = mid4[2];

//  second
    rect2->p1[0] = rect.p2[0];
    rect2->p1[1] = rect.p2[1];
    rect2->p1[2] = rect.p2[2];

    rect2->p2[0] = mid2[0];
    rect2->p2[1] = mid2[1];
    rect2->p2[2] = mid2[2];

    rect2->p3[0] = mid[0];
    rect2->p3[1] = mid[1];
    rect2->p3[2] = mid[2];

    rect2->p4[0] = mid1[0];
    rect2->p4[1] = mid1[1];
    rect2->p4[2] = mid1[2];

//  third
    rect3->p1[0] = rect.p3[0];
    rect3->p1[1] = rect.p3[1];
    rect3->p1[2] = rect.p3[2];

    rect3->p2[0] = mid3[0];
    rect3->p2[1] = mid3[1];
    rect3->p2[2] = mid3[2];

    rect3->p3[0] = mid[0];
    rect3->p3[1] = mid[1];
    rect3->p3[2] = mid[2];

    rect3->p4[0] = mid2[0];
    rect3->p4[1] = mid2[1];
    rect3->p4[2] = mid2[2];

//  fourth
    rect4->p1[0] = rect.p4[0];
    rect4->p1[1] = rect.p4[1];
    rect4->p1[2] = rect.p4[2];

    rect4->p2[0] = mid4[0];
    rect4->p2[1] = mid4[1];
    rect4->p2[2] = mid4[2];

    rect4->p3[0] = mid[0];
    rect4->p3[1] = mid[1];
    rect4->p3[2] = mid[2];

    rect4->p4[0] = mid3[0];
    rect4->p4[1] = mid3[1];
    rect4->p4[2] = mid3[2];
}


void divideIterate( struct Rectangle rect, 
                    int pn,
                    int layer, 
                    struct treeNode *tn ) {

    double mid[3];

    getCenterRect(rect,mid);
    tn->midPoint[0] = mid[0];
    tn->midPoint[1] = mid[1];
    tn->midPoint[2] = mid[2];
    
    int i = 0;
    for(i=0; i<7; i ++) {
        tn->band[i] = 0;
    }

    double angle1 = getAngleInSphere(mid[0],mid[1],mid[2],rect.p1[0],rect.p1[1],rect.p1[2]);
    double angle2 = getAngleInSphere(mid[0],mid[1],mid[2],rect.p2[0],rect.p2[1],rect.p2[2]);
    double angle3 = getAngleInSphere(mid[0],mid[1],mid[2],rect.p3[0],rect.p3[1],rect.p3[2]);
    double angle4 = getAngleInSphere(mid[0],mid[1],mid[2],rect.p4[0],rect.p4[1],rect.p4[2]);

    double anglef = angle1;
    if(anglef<angle2) {
        anglef = angle2;
    }
    if(anglef<angle3) {
        anglef = angle3;
    }
    if(anglef<angle4) {
        anglef = angle4;
    }

    ///////////////////////
    //  这两步是在做什么？？？
    ///////////////////////
    anglef += 1/180.*M_PI;
    tn->radius = sqrt(2-2*cos(anglef));

    //printf("%f\n",tn->radius);

    if (pn == layer) {
        tn->child1 = NULL;
        tn->child2 = NULL;
        tn->child3 = NULL;
        tn->child4 = NULL;
        return;
    }

    struct Rectangle *rect1 = (struct Rectangle*)malloc(sizeof(struct Rectangle)*1);
    struct Rectangle *rect2 = (struct Rectangle*)malloc(sizeof(struct Rectangle)*1);
    struct Rectangle *rect3 = (struct Rectangle*)malloc(sizeof(struct Rectangle)*1);
    struct Rectangle *rect4 = (struct Rectangle*)malloc(sizeof(struct Rectangle)*1);

    devideRect(rect,rect1,rect2,rect3,rect4);

    pn ++;
    struct treeNode* cNode1 = (struct treeNode*) malloc(sizeof(struct treeNode));
    cNode1->parent = NULL;
    cNode1->child1 = NULL;
    cNode1->child2 = NULL;
    cNode1->child3 = NULL;
    cNode1->child4 = NULL;
    cNode1->parent = tn;
    tn->child1 = cNode1;
    tn->cnt_A = 0; // added by XYH
    tn->cnt_B = 0; // added by XYH
    divideIterate(*rect1, pn,layer, cNode1);  // 递归调用

    struct treeNode* cNode2 = (struct treeNode*) malloc(sizeof(struct treeNode));
    cNode2->parent = NULL;
    cNode2->child1 = NULL;
    cNode2->child2 = NULL;
    cNode2->child3 = NULL;
    cNode2->child4 = NULL;
    cNode2->parent = tn;
    tn->child2 = cNode2;
    tn->cnt_A = 0; // added by XYH
    tn->cnt_B = 0; // added by XYH
    divideIterate(*rect2, pn,layer, cNode2);

    struct treeNode* cNode3 = (struct treeNode*) malloc(sizeof(struct treeNode));
    cNode3->parent = NULL;
    cNode3->child1 = NULL;
    cNode3->child2 = NULL;
    cNode3->child3 = NULL;
    cNode3->child4 = NULL;
    cNode3->parent = tn;
    tn->child3 = cNode3;
    tn->cnt_A = 0; // added by XYH
    tn->cnt_B = 0; // added by XYH
    divideIterate(*rect3, pn,layer, cNode3);

    struct treeNode* cNode4 = (struct treeNode*) malloc(sizeof(struct treeNode));
    cNode4->parent = NULL;
    cNode4->child1 = NULL;
    cNode4->child2 = NULL;
    cNode4->child3 = NULL;
    cNode4->child4 = NULL;
    cNode4->parent = tn;
    tn->child4 = cNode4;
    tn->cnt_A = 0; // added by XYH
    tn->cnt_B = 0; // added by XYH
    divideIterate(*rect4, pn,layer, cNode4);

    free(rect1);
    free(rect2);
    free(rect3);
    free(rect4);
}

void searchTrees( struct treeNode* node, 
                  int nLayer, 
                  int searchLayer, 
                  FILE* fn ) {
///////////////////////////////////////////
//  nLayer与searchLayer之间的大小关系是啥？
///////////////////////////////////////////

    if (nLayer == searchLayer) {
        double eqCoor[2];
        if(    node->band[0]==0
            && node->band[1]==0
            && node->band[2]==0
            && node->band[3]==0
            && node->band[4]==0
            && node->band[5]==0
            && node->band[6]==0 ) {
//  这里其实啥也不干。。。
        } else {
            //Cartesian2Equatorial(node->midPoint, eqCoor);
            fprintf(fn, "%f    %f    %f    %d    %d    %d    %d    %d    %d    %d\n",
                        node->midPoint[0],
                        node->midPoint[1],
                        node->midPoint[2],
                        node->band[0], 
                        node->band[1], 
                        node->band[2], 
                        node->band[3], 
                        node->band[4], 
                        node->band[5], 
                        node->band[6]);
        }

        //fprintf(fn,"%f    %f    %f\n",node->midPoint[0], node->midPoint[1], node->midPoint[2]);

        return;
    }

    nLayer ++;
    struct treeNode* ch1 = node->child1;
    if(ch1 == NULL) {
        return;
    } else {
        searchTrees(ch1, nLayer, searchLayer, fn);
    }

    struct treeNode* ch2 = node->child2;
    if(ch2 == NULL) {
        return;
    } else {
        searchTrees(ch2, nLayer, searchLayer, fn);
    }

    struct treeNode* ch3 = node->child3;
    if(ch3 == NULL) {
        return;
    } else {
        searchTrees(ch3, nLayer, searchLayer, fn);
    }

    struct treeNode* ch4 = node->child4;
    if(ch4 == NULL) {
        return;
    } else {
        searchTrees(ch4, nLayer, searchLayer, fn);
    }
}

void freeNode(struct treeNode* node) {
    //if(node->parent != NULL){
    //freeNode(node->parent);
    //}
    if(node->child1 != NULL) {
        freeNode(node->child1);
    }
    if(node->child2 != NULL) {
        freeNode(node->child2);
    }
    if(node->child3 != NULL) {
        freeNode(node->child3);
    }
    if(node->child4 != NULL) {
        freeNode(node->child4);
    }
    free(node);
    node = NULL;
}

void freeTrees(struct FourTree * trees) {
    if (trees->head != NULL) {
        freeNode(trees->head);
    }
    free(trees);
}

void findSearchRange( double lat, 
                      double lon, 
                      struct treeNode *node, 
                      int pn, 
                      int layer, 
                      struct searchNodes* r_nodes ) {

    double px = cos(lat *M_PI/180) * cos(lon *M_PI/180);
    double py = cos(lat *M_PI/180) * sin(lon *M_PI/180);
    double pz = sin(lat *M_PI/180);
    double dist = getDistPP(px,py,pz,node->midPoint[0],node->midPoint[1],node->midPoint[2]);

    if (dist > node->radius) {
        return;
    }

    if (pn == layer) {
        r_nodes-> s_nodes[r_nodes->s_len] = node;
        r_nodes-> s_len ++;
        //printf("%f    %f    %f\n",node->midPoint[0],node->midPoint[1],node->midPoint[2]);
        return;
    }

    pn ++;
    
    struct treeNode *ch1 = node->child1;
    if (ch1 == NULL) {
        return;
    }
    findSearchRange(lat, lon, ch1, pn, layer,r_nodes);

    struct treeNode *ch2 = node->child2;
    if (ch2 == NULL) {
        return;
    }
    findSearchRange(lat, lon, ch2, pn, layer, r_nodes);

    struct treeNode *ch3 = node->child3;
    if (ch3 == NULL) {
        return;
    }
    findSearchRange(lat, lon, ch3, pn, layer, r_nodes);

    struct treeNode *ch4 = node->child4;
    if (ch4 == NULL) {
        return;
    }
    findSearchRange(lat, lon, ch4, pn, layer, r_nodes);
}

void rotateByAxisZ( double* coor,
                    double  theta,
                    double* coor_new ) {
    coor_new[0] = coor[0]*cos(theta)+coor[1]*sin(theta);
    coor_new[1] = -sin(theta)*coor[0]+cos(theta)*coor[1];
    coor_new[2] = coor[2];
}

void rotateByAxisY( double* coor,
                    double  theta, 
                    double* coor_new ) {
    coor_new[0] = cos(theta)*coor[0]-sin(theta)*coor[2];
    coor_new[1] = coor[1];
    coor_new[2] = sin(theta)*coor[0]+cos(theta)*coor[2];
}

void rotateByAxisX( double* coor,
                    double  theta,
                    double* coor_new ) {
    coor_new[0] = coor[0];
    coor_new[1] = cos(theta)*coor[1]+sin(theta)*coor[2];
    coor_new[2] = -sin(theta)*coor[1]+cos(theta)*coor[2];
}


// ???????????????????????
// ccd position
// y,z ----distance from center(lat,lon),left-down
void computCCDPos(  double lat, 
                    double lon,
                    double y, 
                    double z, 
                    struct Rectangle* ccdpos ) {

    double ccd_z = 0.0924/28.0;
    double ccd_y = 0.0922/28.0;
    double p1[3],p2[3],p3[3],p4[3],p11[3],p22[3],p33[3],p44[3];

    double py = y;
    double pz = z;
    p1[0]=sqrt(1 - py*py-pz*pz);
    p1[1]=py;
    p1[2]=pz;
    rotateByAxisY(p1,lat *M_PI/180, p11);
    rotateByAxisZ(p11,-lon *M_PI/180,ccdpos->p1);

    py = y+ccd_y;
    pz = z;
    p2[0]=sqrt(1 - py*py-pz*pz);
    p2[1]=py;
    p2[2]=pz;
    rotateByAxisY(p2,lat *M_PI/180, p22);
    rotateByAxisZ(p22,-lon *M_PI/180,ccdpos->p2);


    py = y+ccd_y;
    pz = z+ccd_z;
    p3[0]=sqrt(1 - py*py-pz*pz);
    p3[1]=py;
    p3[2]=pz;
    rotateByAxisY(p3,lat *M_PI/180, p33);
    rotateByAxisZ(p33,-lon *M_PI/180,ccdpos->p3);

    py = y;
    pz = z+ccd_z;
    p4[0]=sqrt(1 - py*py-pz*pz);
    p4[1]=py;
    p4[2]=pz;
    rotateByAxisY(p4,lat *M_PI/180, p44);
    rotateByAxisZ(p44,-lon *M_PI/180,ccdpos->p4);
}

int IsInRect(double *point,struct Rectangle* rect) {

    //printf("%f    %f    %f\n",rect->p1[0]-point[0],rect->p1[0]-point[1],rect->p1[0]-point[2]);

    double a1,a2,a3;
    double b1,b2,b3;
    double vec[4][3];

    //p1
    a1 = rect->p1[0]-point[0];
    a2 = rect->p1[1]-point[1];
    a3 = rect->p1[2]-point[2];

    b1 = rect->p2[0]-point[0];
    b2 = rect->p2[1]-point[1];
    b3 = rect->p2[2]-point[2];

    vec[0][0]=a2*b3-a3*b2;
    vec[0][1]=a3*b1-a1*b3;
    vec[0][2]=a1*b2-a2*b1;

    //p2
    a1 = rect->p2[0]-point[0];
    a2 = rect->p2[1]-point[1];
    a3 = rect->p2[2]-point[2];

    b1 = rect->p3[0]-point[0];
    b2 = rect->p3[1]-point[1];
    b3 = rect->p3[2]-point[2];

    vec[1][0]=a2*b3-a3*b2;
    vec[1][1]=a3*b1-a1*b3;
    vec[1][2]=a1*b2-a2*b1;

    //p3
    a1 = rect->p3[0]-point[0];
    a2 = rect->p3[1]-point[1];
    a3 = rect->p3[2]-point[2];

    b1 = rect->p4[0]-point[0];
    b2 = rect->p4[1]-point[1];
    b3 = rect->p4[2]-point[2];

    vec[2][0]=a2*b3-a3*b2;
    vec[2][1]=a3*b1-a1*b3;
    vec[2][2]=a1*b2-a2*b1;

    //p4
    a1 = rect->p4[0]-point[0];
    a2 = rect->p4[1]-point[1];
    a3 = rect->p4[2]-point[2];

    b1 = rect->p1[0]-point[0];
    b2 = rect->p1[1]-point[1];
    b3 = rect->p1[2]-point[2];

    vec[3][0]=a2*b3-a3*b2;
    vec[3][1]=a3*b1-a1*b3;
    vec[3][2]=a1*b2-a2*b1;

    int i = 0;
    int reFlag = 1;

    //double flag[4];
    double scale = 100000000;
    double flag1 = (vec[0][0]*vec[1][0]+vec[0][1]*vec[1][1]+vec[0][2]*vec[1][2])*scale;
    //printf("flag %f \n",flag1);

    for(i=1; i<4; i++) {
        double flagn = (vec[i][0]*vec[(i+1)%4][0]+vec[i][1]*vec[(i+1)%4][1]+vec[i][2]*vec[(i+1)%4][2])*scale;
        //printf("flagn %f \n",flagn);
        if(flag1*flagn < 0) {
            reFlag = 0;
            break;
        }
    }

    return reFlag;
}

//--------------------------divide sphere use HealPix----------------------------------------//
void produceHealPixeTrees(struct FourTree* Trees12, int layer) {
/*---------------------1--------------------------*/
    int pn = 0;

    struct Rectangle rect;
    rect.p1[0]=sqrt(2)/2.;
    rect.p1[1]=sqrt(2)/2.;
    rect.p1[2]=0;

    rect.p2[0]=0;
    rect.p2[1]=cos(0.61547970867038737);
    rect.p2[2]=sin(0.61547970867038737);

    rect.p3[0]=-sqrt(2)/2.;
    rect.p3[1]=sqrt(2)/2.;
    rect.p3[2]=0;

    rect.p4[0]=0;
    rect.p4[1]=cos(0.61547970867038737);
    rect.p4[2]=-sin(0.61547970867038737);

    struct treeNode *head1 = (struct treeNode*)malloc(sizeof(struct treeNode));
    head1->parent = NULL;
    head1->child1 = NULL;
    head1->child2 = NULL;
    head1->child3 = NULL;
    head1->child4 = NULL;
    (Trees12+0)->head = head1;
    (Trees12+0)->layer = layer;
    divideIterate(rect,pn,layer, head1);

/*---------------------2--------------------------*/
    pn = 0;
    rect.p1[0]=-sqrt(2)/2.;
    rect.p1[1]=sqrt(2)/2.;
    rect.p1[2]=0;

    rect.p2[0]=-cos(0.61547970867038737);
    rect.p2[1]=0;
    rect.p2[2]=sin(0.61547970867038737);

    rect.p3[0]=-sqrt(2)/2.;
    rect.p3[1]=-sqrt(2)/2.;
    rect.p3[2]=0;

    rect.p4[0]=-cos(0.61547970867038737);
    rect.p4[1]=0;
    rect.p4[2]=-sin(0.61547970867038737);

    struct treeNode *head2 = (struct treeNode*)malloc(sizeof(struct treeNode));
    head2->parent = NULL;
    head2->child1 = NULL;
    head2->child2 = NULL;
    head2->child3 = NULL;
    head2->child4 = NULL;
    (Trees12+1)->head = head2;
    (Trees12+1)->layer = layer;
    divideIterate(rect,pn,layer, head2);

/*---------------------3--------------------------*/
    pn = 0;
    rect.p1[0]=-sqrt(2)/2.;
    rect.p1[1]=-sqrt(2)/2.;
    rect.p1[2]=0;

    rect.p2[0]=0;
    rect.p2[1]=-cos(0.61547970867038737);
    rect.p2[2]=sin(0.61547970867038737);

    rect.p3[0]=sqrt(2)/2.;
    rect.p3[1]=-sqrt(2)/2.;
    rect.p3[2]=0;

    rect.p4[0]=0;
    rect.p4[1]=-cos(0.61547970867038737);
    rect.p4[2]=-sin(0.61547970867038737);

    struct treeNode *head3 = (struct treeNode*)malloc(sizeof(struct treeNode));
    head3->parent = NULL;
    head3->child1 = NULL;
    head3->child2 = NULL;
    head3->child3 = NULL;
    head3->child4 = NULL;
    (Trees12+2)->head = head3;
    (Trees12+2)->layer = layer;
    divideIterate(rect,pn,layer, head3);

/*---------------------4--------------------------*/
    pn = 0;
    rect.p1[0]=sqrt(2)/2.;
    rect.p1[1]=-sqrt(2)/2.;
    rect.p1[2]=0;

    rect.p2[0]=cos(0.61547970867038737);
    rect.p2[1]=0;
    rect.p2[2]=sin(0.61547970867038737);

    rect.p3[0]=sqrt(2)/2.;
    rect.p3[1]=sqrt(2)/2.;
    rect.p3[2]=0;

    rect.p4[0]=cos(0.61547970867038737);
    rect.p4[1]=0;
    rect.p4[2]=-sin(0.61547970867038737);

    struct treeNode *head4 = (struct treeNode*)malloc(sizeof(struct treeNode));
    head4->parent = NULL;
    head4->child1 = NULL;
    head4->child2 = NULL;
    head4->child3 = NULL;
    head4->child4 = NULL;
    (Trees12+3)->head = head4;
    (Trees12+3)->layer = layer;
    divideIterate(rect,pn,layer, head4);

/*---------------------5--------------------------*/
    pn = 0;
    rect.p1[0]=sqrt(2)/2.;
    rect.p1[1]=sqrt(2)/2.;
    rect.p1[2]=0;

    rect.p2[0]=0;
    rect.p2[1]=cos(0.61547970867038737);
    rect.p2[2]=sin(0.61547970867038737);

    rect.p3[0]=0;
    rect.p3[1]=0;
    rect.p3[2]=1;

    rect.p4[0]=cos(0.61547970867038737);
    rect.p4[1]=0;
    rect.p4[2]=sin(0.61547970867038737);

    struct treeNode *head5 = (struct treeNode*)malloc(sizeof(struct treeNode));
    head5->parent = NULL;
    head5->child1 = NULL;
    head5->child2 = NULL;
    head5->child3 = NULL;
    head5->child4 = NULL;
    (Trees12+4)->head = head5;
    (Trees12+4)->layer = layer;
    divideIterate(rect,pn,layer, head5);

/*---------------------6--------------------------*/
    pn = 0;
    rect.p1[0]=-sqrt(2)/2.;
    rect.p1[1]=sqrt(2)/2.;
    rect.p1[2]=0;

    rect.p2[0]=-cos(0.61547970867038737);
    rect.p2[1]=0;
    rect.p2[2]=sin(0.61547970867038737);

    rect.p3[0]=0;
    rect.p3[1]=0;
    rect.p3[2]=1;

    rect.p4[0]=0;
    rect.p4[1]=cos(0.61547970867038737);
    rect.p4[2]=sin(0.61547970867038737);

    struct treeNode *head6 = (struct treeNode*)malloc(sizeof(struct treeNode));
    head6->parent = NULL;
    head6->child1 = NULL;
    head6->child2 = NULL;
    head6->child3 = NULL;
    head6->child4 = NULL;

    (Trees12+5)->head = head6;
    (Trees12+5)->layer = layer;
    divideIterate(rect,pn,layer, head6);

/*---------------------7--------------------------*/
    pn = 0;
    rect.p1[0]=-sqrt(2)/2.;
    rect.p1[1]=-sqrt(2)/2.;
    rect.p1[2]=0;

    rect.p2[0]=0;
    rect.p2[1]=-cos(0.61547970867038737);
    rect.p2[2]=sin(0.61547970867038737);

    rect.p3[0]=0;
    rect.p3[1]=0;
    rect.p3[2]=1;

    rect.p4[0]=-cos(0.61547970867038737);
    rect.p4[1]=0;
    rect.p4[2]=sin(0.61547970867038737);

    struct treeNode *head7 = (struct treeNode*)malloc(sizeof(struct treeNode));
    head7->parent = NULL;
    head7->child1 = NULL;
    head7->child2 = NULL;
    head7->child3 = NULL;
    head7->child4 = NULL;
    (Trees12+6)->head = head7;
    (Trees12+6)->layer = layer;
    divideIterate(rect,pn,layer, head7);

/*---------------------8--------------------------*/
    pn = 0;
    rect.p1[0]=sqrt(2)/2.;
    rect.p1[1]=-sqrt(2)/2.;
    rect.p1[2]=0;

    rect.p2[0]=cos(0.61547970867038737);
    rect.p2[1]=0;
    rect.p2[2]=sin(0.61547970867038737);

    rect.p3[0]=0;
    rect.p3[1]=0;
    rect.p3[2]=1;

    rect.p4[0]=0;
    rect.p4[1]=-cos(0.61547970867038737);
    rect.p4[2]=sin(0.61547970867038737);

    struct treeNode *head8 = (struct treeNode*)malloc(sizeof(struct treeNode));
    head8->parent = NULL;
    head8->child1 = NULL;
    head8->child2 = NULL;
    head8->child3 = NULL;
    head8->child4 = NULL;
    (Trees12+7)->head = head8;
    (Trees12+7)->layer = layer;
    divideIterate(rect,pn,layer, head8);

/*---------------------9--------------------------*/
    pn = 0;
    rect.p1[0]=sqrt(2)/2.;
    rect.p1[1]=sqrt(2)/2.;
    rect.p1[2]=0;

    rect.p2[0]=0;
    rect.p2[1]=cos(0.61547970867038737);
    rect.p2[2]=-sin(0.61547970867038737);

    rect.p3[0]=0;
    rect.p3[1]=0;
    rect.p3[2]=-1;

    rect.p4[0]=cos(0.61547970867038737);
    rect.p4[1]=0;
    rect.p4[2]=-sin(0.61547970867038737);

    struct treeNode *head9 = (struct treeNode*)malloc(sizeof(struct treeNode));
    head9->parent = NULL;
    head9->child1 = NULL;
    head9->child2 = NULL;
    head9->child3 = NULL;
    head9->child4 = NULL;
    (Trees12+8)->head = head9;
    (Trees12+8)->layer = layer;
    divideIterate(rect,pn,layer, head9);

/*---------------------10--------------------------*/
    pn = 0;
    rect.p1[0]=-sqrt(2)/2.;
    rect.p1[1]=sqrt(2)/2.;
    rect.p1[2]=0;

    rect.p2[0]=-cos(0.61547970867038737);
    rect.p2[1]=0;
    rect.p2[2]=-sin(0.61547970867038737);

    rect.p3[0]=0;
    rect.p3[1]=0;
    rect.p3[2]=-1;

    rect.p4[0]=0;
    rect.p4[1]=cos(0.61547970867038737);
    rect.p4[2]=-sin(0.61547970867038737);

    struct treeNode *head10 = (struct treeNode*)malloc(sizeof(struct treeNode));
    head10->parent = NULL;
    head10->child1 = NULL;
    head10->child2 = NULL;
    head10->child3 = NULL;
    head10->child4 = NULL;
    (Trees12+9)->head = head10;
    (Trees12+9)->layer = layer;
    divideIterate(rect,pn,layer, head10);

/*---------------------11--------------------------*/
    pn = 0;
    rect.p1[0]=-sqrt(2)/2.;
    rect.p1[1]=-sqrt(2)/2.;
    rect.p1[2]=0;

    rect.p2[0]=0;
    rect.p2[1]=-cos(0.61547970867038737);
    rect.p2[2]=-sin(0.61547970867038737);

    rect.p3[0]=0;
    rect.p3[1]=0;
    rect.p3[2]=-1;

    rect.p4[0]=-cos(0.61547970867038737);
    rect.p4[1]=0;
    rect.p4[2]=-sin(0.61547970867038737);

    struct treeNode *head11 = (struct treeNode*)malloc(sizeof(struct treeNode));
    head11->parent = NULL;
    head11->child1 = NULL;
    head11->child2 = NULL;
    head11->child3 = NULL;
    head11->child4 = NULL;
    (Trees12+10)->head = head11;
    (Trees12+10)->layer = layer;
    divideIterate(rect,pn,layer, head11);

/*---------------------12--------------------------*/
    pn = 0;
    rect.p1[0]=sqrt(2)/2.;
    rect.p1[1]=-sqrt(2)/2.;
    rect.p1[2]=0;

    rect.p2[0]=cos(0.61547970867038737);
    rect.p2[1]=0;
    rect.p2[2]=-sin(0.61547970867038737);

    rect.p3[0]=0;
    rect.p3[1]=0;
    rect.p3[2]=-1;

    rect.p4[0]=0;
    rect.p4[1]=-cos(0.61547970867038737);
    rect.p4[2]=-sin(0.61547970867038737);

    struct treeNode *head12 = (struct treeNode*)malloc(sizeof(struct treeNode));
    head12->parent = NULL;
    head12->child1 = NULL;
    head12->child2 = NULL;
    head12->child3 = NULL;
    head12->child4 = NULL;
    (Trees12+11)->head = head12;
    (Trees12+11)->layer = layer;
    divideIterate(rect,pn,layer, head12);

//------------------------division end-------------------------------------------//
}

void MarkObervePoint( double lat,
                      double lon,
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

////////////////////////////////////
//  不太明白这里具体是怎么进行面积统计的
//  是否有可能循环20000次以后还是没有找到ccd所在的位置，从而导致面积的统计偏小？？？？
////////////////////////////////////
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