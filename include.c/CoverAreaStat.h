#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <string.h>

#ifndef _COVER_AREA_STAT_H_
#define _COVER_AREA_STAT_H_

struct treeNode{
	struct treeNode* parent;
	struct treeNode* child1;
	struct treeNode* child2;
	struct treeNode* child3;
	struct treeNode* child4;

	double midPoint[3];
	double radius;
	
	int band[7];
	int cnt_A;
	int cnt_B;
};

struct FourTree{
	struct treeNode *head;
	int layer;
};


struct searchNodes{
	struct treeNode* s_nodes[20000];
	int s_len;
};

//	？？？ 需要搞清楚1234的对应位置
struct Rectangle{
	double p1[3];
	double p2[3];
	double p3[3];
	double p4[3];
};

double getDistPP(double x1, double y1, double z1, double x2, double y2, double z2);

double getAngleInSphere(double x1, double y1, double z1, double x2, double y2, double z2);

void findmidArc(double x1,double y1,double z1,double x2,double y2,double z2, double* midP);

void getCenterRect(struct Rectangle rect,double* center);

void devideRect(struct Rectangle rect, struct Rectangle *rect1, struct Rectangle *rect2, struct Rectangle *rect3, struct Rectangle *rect4);

void divideIterate(struct Rectangle rect, int pn,int layer, struct treeNode *tn);

void searchTrees(struct treeNode* node, int nLayer, int searchLayer, FILE* fn);

void freeNode(struct treeNode* node);

void freeTrees(struct FourTree * trees);

void findSearchRange(double lat, double lon, struct treeNode *node, int pn, int layer, struct searchNodes* r_nodes);

void rotateByAxisZ(double* coor,double theta,double* coor_new);
void rotateByAxisY(double *coor,double theta,double *coor_new);
void rotateByAxisX(double* coor,double theta,double* coor_new);

// ccd position
// x,y----distance from center(lat,lon),left-down
void computCCDPos(double lat, double lon,double y, double z, struct Rectangle* ccdpos);

int IsInRect(double *point,struct Rectangle* rect);

void produceHealPixeTrees(struct FourTree* Trees12, int layer);

void MarkObervePoint(double lat, double lon, struct FourTree* Trees12, double ccd_pos_in_focus[18][2], double *area1, double *area2, double deltArea);
void MarkObervePoint_A(double lat, double lon, double isdeep, struct FourTree* Trees12, double ccd_pos_in_focus[18][2], double *area1, double *area2, double deltArea);
void MarkObervePoint_B(double lat, double lon, double isdeep, struct FourTree* Trees12, double ccd_pos_in_focus[18][2], double *area1, double *area2, double deltArea);
void MarkObervePoint_C(double lat, double lon, double isdeep, struct FourTree* Trees12, double ccd_pos_in_focus[18][2], double *area1, double *area2, double deltArea);
#endif //_COVER_AREA_STAT_H_
