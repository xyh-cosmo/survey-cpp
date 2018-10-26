#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multiroots.h>

#define R_EARTH 6400.0
#define H_SAT	400.0

struct SatSun{
	double x_sat;
	double y_sat;
	double z_sat;

	double x_sun;
	double y_sun;
	double z_sun;

	double R;
};

void init_satsun( struct SatSun* ss,
				  double theta_sat,
				  double phi_sat,
				  double theta_sun,
				  double phi_sun ){
	ss->R = (R_EARTH);
	ss->x_sat = (R_EARTH+H_SAT)*sin(theta_sat*M_PI/180)*cos(phi_sat*M_PI/180);
	ss->y_sat = (R_EARTH+H_SAT)*sin(theta_sat*M_PI/180)*sin(phi_sat*M_PI/180);
	ss->z_sat = (R_EARTH+H_SAT)*cos(theta_sat*M_PI/180);

	ss->x_sun = sin(theta_sun*M_PI/180)*cos(phi_sun*M_PI/180);
	ss->y_sun = sin(theta_sun*M_PI/180)*sin(phi_sun*M_PI/180);
	ss->z_sun = cos(theta_sun*M_PI/180);
	
	printf("--> position of Sat:\n");
	printf("(%15.8f, %15.8f, %15.8f)\n",ss->x_sat,ss->y_sat,ss->z_sat);
	
	printf("--> position of Sun:\n");
	printf("(%15.8f, %15.8f, %15.8f)\n",ss->x_sun,ss->y_sun,ss->z_sun);
}

int target_function1( const gsl_vector* 	X,
					  void* 				params,
					  gsl_vector* 			fvals ){
	
	struct SatSun* SS = (struct SatSun*)(params);
	double d = SS->x_sat;
	double e = SS->y_sat;
	double f = SS->z_sat;
	double g = SS->x_sun;
	double h = SS->y_sun;
	double i = SS->z_sun;
	
	double R = SS->R;

	double x = gsl_vector_get(X, 0);
	double y = gsl_vector_get(X, 1);
	double z = gsl_vector_get(X, 2);

/*	double df0 = z - sqrt(R*R - x*x - y*y);*/
	double df0 = z*z + x*x + y*y - R*R;
	double df1 = d*x+e*y+f*z - R*R;
	double df2 = g*x+h*y+i*z;
	
	gsl_vector_set(fvals, 0, df0);
	gsl_vector_set(fvals, 1, df1);
	gsl_vector_set(fvals, 2, df2);

	return GSL_SUCCESS;
}

int target_function2(	const gsl_vector*	X,
						void*				params,
						gsl_vector*			fvals ){
	struct SatSun* SS = (struct SatSun*)(params);
	double d = SS->x_sat;
	double e = SS->y_sat;
	double f = SS->z_sat;
	double g = SS->x_sun;
	double h = SS->y_sun;
	double i = SS->z_sun;

	double R = SS->R;

	double x = gsl_vector_get(X,0);
	double y = gsl_vector_get(X,1);
	double z = gsl_vector_get(X,2);

	double df0 = z + sqrt(R*R-x*x-y*y);
	double df1 = d*x+e*y+f*z - R*R;
	double df2 = g*x+h*y+i*z;

	gsl_vector_set(fvals,0,df0);
	gsl_vector_set(fvals,1,df1);
	gsl_vector_set(fvals,2,df2);

	return GSL_SUCCESS;
}

int print_state( size_t iter, gsl_multiroot_fsolver *s ){
    printf("iter = %3lu x = %15.8f %15.8f %15.8f f(x) = % .3e % .3e % .3e\n",
            iter,
            gsl_vector_get(s->x,0),
            gsl_vector_get(s->x,1),
            gsl_vector_get(s->x,2),
            gsl_vector_get(s->f,0),
            gsl_vector_get(s->f,1),
            gsl_vector_get(s->f,2));
}

int main(){
	const gsl_multiroot_fsolver_type *T;
	gsl_multiroot_fsolver *s1, *s2;

	int status;
	size_t i, iter;

	const size_t n = 3;
	struct SatSun sat_sun;

	init_satsun( &sat_sun,
				 0.0,
				 0.0,
				 90.0,
				 90.0 );

    gsl_multiroot_function fun1 = {&target_function1, n, &sat_sun};
    gsl_multiroot_function fun2 = {&target_function2, n, &sat_sun};
    
/*    double x_init[3] = {(R_EARTH+H_SAT),0,0};*/
    double x_init[3] = {-(R_EARTH+H_SAT),0,0};
    
    gsl_vector* x = gsl_vector_alloc(n);
    
    T = gsl_multiroot_fsolver_hybrids;
    s1 = gsl_multiroot_fsolver_alloc(T,n);
    s2 = gsl_multiroot_fsolver_alloc(T,n);
    
    //  solve first set of eqns
    gsl_vector_set(x,0,x_init[0]-1);
    gsl_vector_set(x,1,x_init[1]+2);
    gsl_vector_set(x,2,x_init[2]-400);    
    gsl_multiroot_fsolver_set(s1,&fun1,x);
    
    iter = 0;
    do{
        iter++;
        status = gsl_multiroot_fsolver_iterate(s1);
        
        print_state(iter,s1);

        if( status )
            break;

        status = gsl_multiroot_test_residual(s1->f, 1e-10);
    } while( status == GSL_CONTINUE && iter < 1000 );
    
    printf("\nsolutions for the 1st set of eqns:\n");
    printf("x = %20.10f\n", gsl_vector_get(s1->x,0));
    printf("y = %20.10f\n", gsl_vector_get(s1->x,1));
    printf("z = %20.10f\n\n", gsl_vector_get(s1->x,2));
   
    //  solve second set of eqns
    gsl_vector_set(x,0,x_init[0]+1);
    gsl_vector_set(x,1,x_init[1]-2);
    gsl_vector_set(x,2,x_init[2]+400);    
    gsl_multiroot_fsolver_set(s2,&fun2,x);

/*    iter = 0;*/
/*    do{*/
/*        iter++;*/
/*        status = gsl_multiroot_fsolver_iterate(s2);*/
/*        */
/*        print_state(iter,s2);*/
/*        */
/*        if( status )*/
/*            break;*/

/*        status = gsl_multiroot_test_residual(s2->f, 1e-10);*/
/*    } while( status == GSL_CONTINUE && iter < 1000 );*/

/*    printf("\nsolutions for the 2st set of eqns:\n");*/
/*    printf("x = %20.10f\n", gsl_vector_get(s2->x,0));*/
/*    printf("y = %20.10f\n", gsl_vector_get(s2->x,1));*/
/*    printf("z = %20.10f\n\n", gsl_vector_get(s2->x,2));*/
    
    gsl_multiroot_fsolver_free(s1);
    gsl_multiroot_fsolver_free(s2);
    gsl_vector_free(x);
}
