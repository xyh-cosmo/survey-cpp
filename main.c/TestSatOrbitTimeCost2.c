
#include "SurveySim.h"

int main(int argc, char *argv[]) {

    time_t stime, etime;

    SatelliteOrbits *sat_orbits = malloc(sizeof(SatelliteOrbits));

    SatelliteOrbits_init(sat_orbits,"orbit20160925");

    double t_start = 2459766;
    double dist;

    time(&stime);

    double x[3];
    double t = t_start;
    while( t < t_start + 365.25*10 ){
        if( get_satellite_position(t,x,&dist,sat_orbits) == 0 ){
            // printf(" Yr = %12.6f, x = %12.6f, y = %12.6f, z = %12.6f, dist = %2.6f\n",
            //         (t-2459766)/365.25,x[0],x[1],x[2],dist);
        }

        t += 1./86400;
        // t += 10;
    }

    time(&etime);

    printf("--> finished search for unused observational time intervals\n");
    printf("--> use time is %ld s\n", etime - stime);

    SatelliteOrbits_free(sat_orbits);

    return 0;
}
