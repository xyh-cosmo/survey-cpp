
#include "GSL_funcs.h"

#ifndef _CMG_H_
#define _CMG_H_


//  compute slew time. added @2018-10-25
#if defined(_SLEW_TIME_GSL_INTERP_)
gsl_interp_accel *acc_trans;
gsl_spline *spline_trans;
void init_TransTime();
void free_TransTime();
#endif

double calculateTransTime(double);


#define CMG_THRES 1.0
#define ORBIT_TIME 0.0638

struct CMG_Node{
	double s_time;		// 队列起始时间
	double use_time; 	// exposure time + trans time
	double cmg_value;	// 权重。。。
	struct CMG_Node *next;
};

struct CMG_List {
	struct CMG_Node* head;
	struct CMG_Node* tail;
	double cmg_total;	// 总的权重
	double time_dura;	// 队列的总时间
};

// #ifndef CMG_LIST
//    #define CMG_LIST
   struct CMG_List* cmg_list;
// #endif

struct CMG_Node* PopCMGNode(struct CMG_List*);
void initCMGNode(struct CMG_Node*,double, double, double, struct CMG_Node*);
void initCMGList(struct CMG_List*);
void addCMG_Node(struct CMG_List*,struct CMG_Node*);
void freeCMG_Node(struct CMG_Node*);
void freeCMG_List(struct CMG_List*);

double get_cmg_use(double tAngle);

void reset_cmg(struct CMG_List* cmg_list);

double Test_CMGConstraint(  struct CMG_List* cmg_list, 
                            double curTime, 
                            double tTime, 
                            double exTime,
                            double cmg_use );

void Update_CMGList( struct CMG_List *cmg_list, 
                     double otime, 
                     double final_cmg_value,
                     double final_use_time );

//////////////////////////////////
//  以下是用于一轨、一轨地测试CMG条件的
struct CMG_one_orbit_state_{
    double one_orbit_time;
    double time_head, time_tail;
    double cmg_total;
};

typedef struct CMG_one_orbit_state_ CMG_one_orbit_state;

int CMG_one_orbit_state_init( CMG_one_orbit_state *cmg, 
                              double jt_start );
int CMG_one_orbit_state_reset( CMG_one_orbit_state *cmg,
                                double jt_new,
                                double cmg_use, 
                                int p_rank );
int CMG_one_orbit_state_check( CMG_one_orbit_state *cmg, 
                               double curTime,
                               double tTime,
                               double exTime,
                               double cmg_use,
                               int p_rank );
void CMG_one_orbit_state_update( CMG_one_orbit_state *cmg, 
                                 double curTime,
                                 double tTime,
                                 double exTime,
                                 double cmg_use,
                                 int p_rank );


CMG_one_orbit_state *cmg_state;


//////////////////////////////////
//  根据一轨内的总的机动时间来判断CMG限制条件
#endif