/*
   vtransit - program to print positions of Sun, Venus, and Earth
              during the Venus Transit of 8 June 2004 at 08:19:44
*/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h> /* needed because we call sqrt() */

#include "ephcom.h"


main(int argc, char* argv[], int pacth) {

    struct ephcom_Header header1;
    struct ephcom_Coords coords;
    double *datablock;  /* Will hold coefficients from a data block */
    int datablocknum;
    int i;
    double dist, sqrt(double); /* To get Geocentric distance */
    /*
           deltat = 32.184 000 sec + (TAI-UTC) - (UT1-UTC)

       From IERS Bulletin A, vol. XVII no. 019, 13 May 2004:

                    32.184 000 sec
          TAI-UTC = 32.000 000 sec
       -(UT1-UTC) =  0.469 74  sec
                    --------------
          deltat  = 64.653 74  sec
    */
    double deltat = 64.65374; /* value on 8 June 2004 */
    /*
       Test parameter variables.
    */
    double testjd; /* Test Julian Day */
    int testntarg, testnctr; /* Target body and center body */
    int testncoord; /* Coordinate number in testr to compare */
    double testr[6]; /* To hold x, xdot, y, ydot, z, zdot for all bodies */
    int idate[6];    /* To hold year, month, day, hour, minute, second */
    double ephcom_cal2jd(int *, int);
    /*
       Output file parameters.
    */
    FILE *infp, *outfp=stdout, *fopen();

    if (argc < 2) {
        fprintf(stderr,
                "\nFormat:\n\n         %s ephemeris-file\n\n",
                argv[0]);
        exit(1);
    }

    if ((infp = fopen(argv[1],"r")) == NULL) {
        fprintf(stderr,"\nERROR: Can't open ephemeris file %s for input.\n\n", argv[1]);
        exit(1);
    }

    ephcom_readbinary_header(infp, &header1);
// ephcom_writeascii_header(outfp, &header1);
    /*
       Done with header.  Now we'll read and write data blocks.
    */
    datablock = (double *)malloc(header1.ncoeff * sizeof(double));
    datablocknum = 0;
    coords.km = 0;        /* Use AU, not kilometers */
    coords.seconds = 0;   /* Timescale is days, not seconds */
    coords.bary = 1;      /* Center is Solar System Barycenter */
    coords.coordtype=0;   /* No correction for light travel time or
                         relativistic effects from Sun */
    /*
       2004.June.08 00:00:00
    */
    idate[0] = 2004;
    idate[1] =  6;
    idate[2] =  8;
    idate[3] = idate[4] = idate[5] = 0;

    testjd = ephcom_cal2jd(idate, 0);

    fprintf(outfp,"ephemeris.com vtransit program.  ");
    fprintf(outfp,"Last modified May 2004.\n\n");
    fprintf(outfp,"vtransit was written to commemorate the historic");
    fprintf(outfp," Venus Transit of 8 June 2004,\n");
    fprintf(outfp,"which reaches a Geocentric mid-transit at");
    fprintf(outfp," approximately 08:19:44 UTC.\n\n");

    /*
       Make sure time is in range that ephemeris covers.
    */
    if (testjd >= header1.ss[0] && testjd <= header1.ss[1]) {
        /* If time is in range, process data. */
        coords.et2[0] = testjd;  /* Good enough precision for test dates */
        coords.et2[1] = 0.0;
        if (ephcom_get_coords(infp, &header1, &coords, datablock) == 0) {
            ephcom_pleph(&coords, EPHCOM_SUN, EPHCOM_EARTH, testr);
            fprintf(outfp, "Check of program operation at 2004.06.08 00:00:00\n\n");
            fprintf(outfp, "   Julian Day is %15.7f\n\n", testjd);
            fprintf(outfp, "           Astronomical        Calculated\n");
            fprintf(outfp, "           Almanac 2004          Value\n");
            fprintf(outfp, "              (in AU)           (in AU)\n");
            fprintf(outfp, "   Sun x    +0.2199082   %20.17f\n", testr[0]);
            fprintf(outfp, "   Sun y    +0.9091691   %20.17f\n", testr[1]);
            fprintf(outfp, "   Sun z    +0.3941587   %20.17f\n", testr[2]);
            dist = sqrt(testr[0]*testr[0] + testr[1]*testr[1] + testr[2]*testr[2]);
            fprintf(outfp, "   Sun d     1.0150415   %20.17f\n", dist);
            ephcom_pleph(&coords, EPHCOM_VENUS, EPHCOM_EARTH, testr);
            dist = sqrt(testr[0]*testr[0] + testr[1]*testr[1] + testr[2]*testr[2]);
            fprintf(outfp, "   Venus d   0.2888953   %20.17f\n", dist);
            fprintf(outfp, "\n   d = sqrt(x^2 + y^2 + z^2)\n\n");
        }
        /*
           Now adjust for the time we want, 08:19:44 + Delta T
        */
        idate[3] =  8;
        idate[4] = 19;
        idate[5] = 44;
        testjd = ephcom_cal2jd(idate, 0) + deltat / (60 * 60 * 24);
        coords.et2[0] = testjd;  /* Good enough precision for test dates */
        coords.et2[1] = 0.0;
        if (ephcom_get_coords(infp, &header1, &coords, datablock) == 0) {
            ephcom_pleph(&coords, EPHCOM_SUN, EPHCOM_EARTH, testr);
            fprintf(outfp,
                    "Calculated positions (in AU) and velocities (in AU/day)");
            fprintf(outfp,
                    " at 2004.06.08 08:19:44\n\n");
            fprintf(outfp, "   Julian Day is %15.7f (Delta T is %9.5f seconds)\n\n",
                    testjd, deltat);
            for (i=0; i<3; i++) {
                fprintf(outfp, "   %-18s %c: %20.17f", "Geocentric Sun",
                        'x'+i, testr[i]);
                fprintf(outfp, "   %c dot: %20.17f\n", 'x'+i, testr[3+i]);
            }
            fprintf(outfp, "\n");
            ephcom_pleph(&coords, EPHCOM_VENUS, EPHCOM_EARTH, testr);
            for (i=0; i<3; i++) {
                fprintf(outfp, "   %-18s %c: %20.17f", "Geocentric Venus",
                        'x'+i, testr[i]);
                fprintf(outfp, "   %c dot: %20.17f\n", 'x'+i, testr[3+i]);
            }
            fprintf(outfp, "\n");
        }
        else {
            fprintf(outfp, "Julian Day %15.2f not found.\n", testjd);
        }
    }
    else {
        fprintf(outfp,
                "\nERROR: Ephemeris file doesn't cover Julian Day %15.2f\n",
                testjd);
        exit(1);
    }

    exit(0);
}

