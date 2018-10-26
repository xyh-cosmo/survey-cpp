#include "GSL_funcs.h"

int main(){

    gsl_rng* r_seed = gsl_rng_alloc(gsl_rng_mt19937);

    srand(time(NULL));
    gsl_rng_set(r_seed,rand());

    double theta_old, phi_old;
    double theta_new, phi_new;

    double angle1, angle2;

    int cnt = 0;
    while( cnt < 10 ){

        theta_old = gsl_ran_flat(r_seed, 0, 2*M_PI);
        theta_new = gsl_ran_flat(r_seed, 0, 2*M_PI);
        phi_old = gsl_ran_flat(r_seed, 0, M_PI);
        phi_new = gsl_ran_flat(r_seed, 0, M_PI);

        // Get_RotationAngle(phi_old, theta_old, phi_new, theta_new, &angle1);
        // Get_RotationAngle_faster(phi_old, theta_old, phi_new, theta_new, &angle2);

        printf("angle results: angle1 = %18.15f angle2 = %18.15f\n", angle1, angle2);

        cnt++;
    }

    gsl_rng_free(r_seed);

}