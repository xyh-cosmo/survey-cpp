#include "SurveySim.h"

#if defined(_SLEW_TIME_GSL_INTERP_)

void init_TransTime(){
    double trans_angle[12] = {0,1,5,10,20,35,45,75,90,110,150,180};
    double trans_time[12]  = {41,41,63,80,96,139,168,253,296,382,468,553};
    acc_trans = gsl_interp_accel_alloc();
    spline_trans = gsl_spline_alloc(gsl_interp_linear,12);
    gsl_spline_init(spline_trans,trans_angle,trans_time,12);
}

void free_TransTime(){
    gsl_spline_free(spline_trans);
    gsl_interp_accel_free(acc_trans);
}

#endif

void initCMGNode( struct CMG_Node* node,
                  double st,
                  double cv,
                  double ut,
                  struct CMG_Node* next ) {
    node->s_time=st;
    node->cmg_value = cv;
    node->use_time = ut;
    node->next = next;
}

void initCMGList(struct CMG_List* list) {
    list->head = NULL;
    list->tail = NULL;
    list->cmg_total = 0;
    list->time_dura = 0;
}

void addCMG_Node( struct CMG_List* list,
                  struct CMG_Node* node ) {
    if (list->head == NULL) {
        list->head = node;
        list->tail = node;
        list->cmg_total = list->cmg_total + node->cmg_value;
        list->time_dura = node->use_time;
    } else {
        list->tail->next = node;
        list->tail = list->tail->next;
        list->tail->next = NULL;
        list->cmg_total = list->cmg_total + node->cmg_value;
        double t1 = list->head->s_time;
        double t2 = list->tail->s_time;
        list->time_dura = t2-t1+list->tail->use_time;
    }
}


struct CMG_Node* PopCMGNode(struct CMG_List* list) {
    struct CMG_Node* rNode = list->head;
    if(list->head->next != NULL) {
        double t1 = list->head->s_time;
        double cmg_v = list->head->cmg_value;

        list->head = list->head->next;

        double t2 = list->head->s_time;
        list->time_dura = list->time_dura - (t2-t1);
        list->cmg_total = list->cmg_total - cmg_v;
    } else {
        list->head = NULL;
        list->tail = NULL;
        list->time_dura = 0;
        list->cmg_total = 0;
    }
    return rNode;
}

void freeCMG_Node(struct CMG_Node* node) {
    node->next = NULL;
    free(node);
}

void freeCMG_List(struct CMG_List* list) {
    while(list->head != NULL) {
        struct CMG_Node* node = list->head;
        list->head = list->head->next;
        freeCMG_Node(node);
    }
    free(list);
}

double get_cmg_use(double tAngle){

    double cmg_use=0.0;
    
#if defined(_OLD_CMG_USE_)

    double tm_k;
    if( tAngle <= 5.0 ) {
        cmg_use = 0;
    } else if( tAngle <= 10.0 ) {
        tm_k = (1/29.0-0)/(10.0-5.0);
        cmg_use = tm_k*(tAngle-5);
    } else if( tAngle <= 20.0 ) {
        tm_k = (1/19.0-1/29.0)/(20.0-10.0);
        cmg_use = tm_k*(tAngle-10)+1/29.0;
    } else if( tAngle <= 35.0 ) {
        tm_k = (1/13.0-1/19.0)/(35.0-20.0);
        cmg_use = tm_k*(tAngle-20)+1/19.0;
    } else if( tAngle <= 45.0 ) {
        tm_k = (1/10.0-1/13.0)/(45.0-35.0);
        cmg_use = tm_k*(tAngle-35)+1/13.0;
    } else if( tAngle <= 75.0 ) {
        tm_k = (1/6.0-1/10.0)/(75.0-45.0);
        cmg_use = tm_k*(tAngle-45)+1/10.0;
    } else if( tAngle <= 90.0 ) {
        tm_k = (1/5.0-1/6.0)/(90.0-75.0);
        cmg_use = tm_k*(tAngle-75)+1/6.0;
    } else if( tAngle <= 135.0 ) {
        tm_k = (1/3.0-1/5.0)/(135.0-90.0);
        cmg_use = tm_k*(tAngle-90)+1/5.0;
    } else if( tAngle <= 180.0 ) {
        tm_k = (1/2.0-1/3.0)/(180.0-135.0);
        cmg_use = tm_k*(tAngle-135)+1/3.0;
    }

#else

    if( tAngle < 5.0 ) {
        cmg_use = 0;
    } else if( tAngle >=  5.0 && tAngle < 10.0 ) {
        cmg_use = 1./29.;
    } else if( tAngle >= 10.0 && tAngle < 20.0 ) {
        cmg_use = 1./19.;
    } else if( tAngle >= 20.0 && tAngle < 35.0 ) {
        cmg_use = 1./13.;
    } else if( tAngle >= 35.0 && tAngle < 45.0 ) {
        cmg_use = 1./10.;
    } else if( tAngle >= 45.0 && tAngle < 75.0 ) {
        cmg_use = 1./6.;
    } else if( tAngle >= 75.0 && tAngle < 90.0 ) {
        cmg_use = 1./5.;
    } else if( tAngle >= 90.0 && tAngle < 135.0 ) {
        cmg_use = 1./3.;
    } else if( tAngle >= 135.0 && tAngle <= 180.0 ) {
        cmg_use = 1./2.;
    }

#endif

    return cmg_use;
}


void reset_cmg(struct CMG_List* cmg_list){
    free(cmg_list);
    cmg_list = (struct CMG_List*)malloc(sizeof(struct CMG_List));
    initCMGList(cmg_list);
}


double Test_CMGConstraint( struct CMG_List* cmg_list, 
                           double curTime, 
                           double tTime, 
                           double exTime,
                           double cmg_use ){
    double cmg_times = 0;  // 用于获取当前队列对应的时间长度

    if( cmg_list->head != NULL ){
        cmg_times = curTime - cmg_list->head->s_time + (tTime+exTime)/86400.0; // ????? 这里是否应该包含exTime ?????
    }

    double cmg_total_use = cmg_list->cmg_total + cmg_use;
    struct CMG_Node* temp_node = cmg_list->head;

    //  如果当前队列的时间长度超过了一轨的时间，就将队列头部的元素给pop掉，并重新计算CMG的总消耗量
    while(    cmg_times > ORBIT_TIME
           && temp_node != NULL
           && temp_node->next != NULL ){
        cmg_total_use = cmg_total_use - temp_node->cmg_value;
        double ct1 = temp_node->s_time;
        temp_node = temp_node->next;
        double ct2 = temp_node->s_time;
        cmg_times = cmg_times-(ct2-ct1);
    }

    return cmg_total_use;
}



void Update_CMGList( struct CMG_List *cmg_list, 
                     double otime, 
                     double final_cmg_value,
                     double final_use_time ){
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
}

// double SmartSearch( double curTime, 
//                     double lastObsEndTime,
//                     double tTime,
//                     double tAngle ){
//     double cmg_use = 0;
//     double t_total = (curTime - *lastObsEndTime)*86400+tTime;
//     int try_max_steps = (int)( t_total/90 );
//     if( try_max_steps > 0 ){
//         int steps = 1;
//         while( steps <= try_max_steps ){
//             double t_left = t_total - steps*90.0;
//             double angle_left = tAngle - steps*5.0;
//             double tTime_tmp = calculateTransTime(angle_left);
//             if( tTime_tmp <= t_left ){
//                 try5degs_works = true;
//                 cmg_use = get_cmg_use(angle_left);
//                 tTime = steps*90 + tTime_tmp;
//                 break;
//             }
//             steps += 1;
//         }
//     }

//     return cmg_use;
// }


//////////////////////////////////////////////////////////////////////


int CMG_one_orbit_state_init( CMG_one_orbit_state *cmg,
                              double jt_start ){
    cmg->one_orbit_time = ORBIT_TIME;
    cmg->time_head = jt_start;
    cmg->time_tail = cmg->time_head;
    cmg->cmg_total = 0.0;
    return 0;
}

int CMG_one_orbit_state_reset( CMG_one_orbit_state *cmg,
                                double jt_new,
                                double cmg_use,
                                int p_rank ){
    // if( p_rank == 0 ){
    //     printf("### resetting CMG state ....\n");
    // }
    cmg->time_head = jt_new;
    cmg->time_tail = cmg->time_head;
    cmg->cmg_total = cmg_use;
    return 0;
}

int CMG_one_orbit_state_check( CMG_one_orbit_state *cmg, 
                               double curTime,
                               double tTime,
                               double exTime,
                               double cmg_use,
                               int p_rank ){
    // if( p_rank == 0 ){
    //     printf("### checking CMG state ....\n");
    // }
    int status = 0;
    double cmg_total_use = cmg->cmg_total;
    double cmg_total_time= cmg->time_tail - cmg->time_head;
    
    if( cmg_total_time + tTime/86400 < cmg->one_orbit_time ){
        if( cmg_total_use + cmg_use > CMG_THRES )
            status = 1;
    } else {
        double frac = (cmg->one_orbit_time - cmg_total_time)/(tTime/86400);
        if( cmg_total_use + cmg_use*frac > CMG_THRES )
            status = 1;
    }

    return status;
}

void CMG_one_orbit_state_update( CMG_one_orbit_state *cmg, 
                                 double curTime,
                                 double tTime,
                                 double exTime,
                                 double cmg_use,
                                 int p_rank ){
    // if( p_rank == 0 ){
    //     printf("### updating CMG state ....\n");
    // }
    double cmg_total_time = cmg->time_tail - cmg->time_head;

    if( cmg_total_time + (tTime+exTime)/86400 < cmg->one_orbit_time ){
        cmg->time_tail += (tTime+exTime)/86400;
        cmg->cmg_total += cmg_use;
    } else {
        //  此处重新设置时间起点的时候要小心
        //  时间终点以曝光结束为基准
        if( cmg_total_time + tTime/86400 < cmg->one_orbit_time ){
            double dt = cmg->one_orbit_time - (cmg_total_time + tTime/86400);
            cmg->time_head = curTime + tTime/86400 + dt;
            cmg->time_tail = curTime + (tTime+exTime)/86400;
            cmg->cmg_total = 0;
        } else if ( cmg_total_time + tTime/86400 > cmg->one_orbit_time ) {
            double frac = (cmg->one_orbit_time - cmg_total_time)/(tTime/86400);
            cmg->time_head = curTime + (cmg->one_orbit_time - cmg_total_time);
            cmg->time_tail = curTime + (tTime+exTime)/86400;
            cmg->cmg_total = (1.0-frac)*cmg_use;
        }
    }
}


// int main(int argc, char* argv[]){

// 	struct CMG_List * list = NULL;
// 	list = (struct CMG_List*)malloc(sizeof(struct CMG_List));
// 	initCMGList(list);

// 	if (list ==NULL){
// 		printf("Error !!\n");
// 	}
// 	int i = 0;
// 	for(i = 0;i < 100; i ++){

// 		struct CMG_Node* node = NULL;
// 		node = (struct CMG_Node*) malloc(sizeof(struct CMG_Node));
// 		initCMGNode(node,i+10,0.01,2.0,NULL);
// 		addCMG_Node(list,node);

// 	}

// 	struct CMG_Node* i_node = list->head;
// 	while(i_node !=NULL){
// 		printf("%f    %f\n",list->cmg_total,list->time_dura);
// 		i_node = i_node->next;
// 	}


// 	freeCMG_List(list);

// 	return 1;
// }
