/*
 testeph - program to interpolate positions and velocieties
 from a binary ephemeris file, and compare results
 with test data from JPL.
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "ephcom.h"

main(int argc, char* argv[], int pacth) {

    struct ephcom_Header header1;
    struct ephcom_Coords coords;
    double *datablock; /* Will hold coefficients from a data block */
    int datablocknum;
    int i;
    /*
     Test parameter variables.
     */
    char testline[EPHCOM_MAXLINE + 1]; /* To read line from TEST file */
    int testnumde; /* ephemeris number; make sure it matches ephemeris */
    char testdate[16]; /* To hold "yyyy.mm.dd" with possible negative year */
    double testjd; /* Test Julian Day */
    int testntarg, testnctr; /* Target body and center body */
    int testncoord; /* Coordinate number in testr to compare */
    double testxi; /* Pre-calculated value to compare with testr */
    int outline; /* Line we've read in so far */
    double testr[6]; /* To hold x, xdot, y, ydot, z, zdot for all bodies */
    double testdel; /* Difference between TEST file result and calculation */
    /*
     Output file parameters.
     */
    FILE *infp, *outfp, *testfp, *fopen();

    if (argc < 3) {
        fprintf(stderr,
                "\nFormat:\n\n         %s ephemeris-file testpo-file [output-file]\n\n",
                argv[0]);
        exit(1);
    }

    if (strcmp(argv[1], "-") == 0) {
        fprintf(stderr, "\nERROR: Can't open ephemeris file on stdin.\n\n");
        exit(1);
    }
    if ((infp = fopen(argv[1], "r")) == NULL) {
        fprintf(stderr, "\nERROR: Can't open ephemeris file %s for input.\n\n",
                argv[1]);
        exit(1);
    }
    if (strcmp(argv[2], "-") == 0) {
        fprintf(stderr, "\nERROR: Can't read TESTPO file on stdin.\n\n");
        exit(1);
    } else if ((testfp = fopen(argv[2], "r")) == NULL) {
        fprintf(stderr, "\nERROR: Can't open testpo file %s for input.\n\n",
                argv[2]);
        exit(1);
    }
    if (argc >= 4) {
        if (strcmp(argv[3], "-") == 0) {
            outfp = stdout;
        } else if ((outfp = fopen(argv[3], "r")) == NULL) { /* File doesn't exist */
            if ((outfp = fopen(argv[3], "w")) == NULL) {
                fprintf(stderr, "\nERROR: Can't open %s for output.\n\n",
                        argv[3]);
                exit(1);
            }
        } else {
            fprintf(stderr, "\nERROR: output file %s already exists.\n",
                    argv[3]);
            exit(1);
        }
    } else
        outfp = stdout;

    ephcom_readbinary_header(infp, &header1);
// ephcom_writeascii_header(outfp, &header1);
    /*
     Done with header.  Now we'll read and write data blocks.
     */
    datablock = (double *) malloc(header1.ncoeff * sizeof(double));
    datablocknum = 0;
    coords.km = 0; /* Use AU, not kilometers */
    coords.seconds = 0; /* Timescale is days, not seconds */
    coords.bary = 1; /* Center is Solar System Barycenter */
    coords.coordtype = 0; /* No correction for light travel time or
	 relativistic effects from Sun */

    testjd = header1.ss[0]; /* Initialize to be within range */
    outline = 0;

    fprintf(outfp, "ephemeris.com testeph program.  ");
    fprintf(outfp, "Last modified May 2004.\n\n");
    fprintf(outfp, "  line -- jed -- t# c# x# ");
    fprintf(outfp, "-- jpl value ---  -- user value --  -- difference --\n");

// print header
    int j = 0;
    i = 0;
    fprintf(outfp, "ksize=%d  ncoeff=%d  ncon=%d    nval=%d    au=%f \n",
            header1.ksize, header1.ncoeff, header1.ncon, header1.nval,
            header1.au);
    fprintf(outfp, "emrat=%f   clight=%f  numde=%d  numle=%d    maxcheby=%d \n",
            header1.emrat, header1.clight, header1.numde, header1.numle,
            header1.maxcheby);

    for (i = 0; i < 12; i++) {
        for (j = 0; j < 3; j++) {
            fprintf(outfp, "cnam[%d][%d]=%d    ", i, j, header1.ipt[i][j]);
        }
        fprintf(outfp, "\n\n");
    }
    for (j = 0; j < 3; j++) {
        fprintf(outfp, "ss[%d]=%f    ", j, header1.ss[j]);
    }
    fprintf(outfp, "\n\n");

    fprintf(outfp, "\n\n");

    fprintf(outfp, "\n\n");
    i = 0;

    coords.et2[0] = 2453371.5;
    coords.et2[1] = 0;
    ephcom_get_coords(infp, &header1, &coords, datablock);
    ephcom_pleph(&coords, 11, 13, testr);
    printf("***********************************************************************************************\n");
    printf("%f    %f     %f     %f      %f     %f \n",testr[0],testr[1],testr[2],testr[3],testr[4],testr[5]);
    printf("***********************************************************************************************\n");

    while (testjd <= header1.ss[1]
            && fgets(testline, EPHCOM_MAXLINE, testfp) != NULL ) {
        if (isdigit(testline[0])) { /* read data */
            sscanf(testline, "%d %s %lf %d %d %d %lf", &testnumde, testdate,
                   &testjd, &testntarg, &testnctr, &testncoord, &testxi);
            /* Make sure time is in range that ephemeris covers. */
            if (testjd >= header1.ss[0] && testjd <= header1.ss[1]) {
//       fprintf(outfp, "%s\n", testline); */
//       fprintf(outfp, "%d\t%s\t%25.17E\t%d\t%d\t%d\t%25.17E\n",
//               testnumde, testdate, testjd, testntarg, testnctr,
//               testncoord, testxi);
                /* If time is in range, process data. */
                coords.et2[0] = testjd; /* Good enough precision for test dates */
                coords.et2[1] = 0.0;
                if (ephcom_get_coords(infp, &header1, &coords, datablock)
                        == 0) {
                    ephcom_pleph(&coords, testntarg, testnctr, testr);
                    testdel = abs(testr[testncoord - 1] - testxi);
                    if (testntarg == 15 && testncoord == 3)
                        testdel /= 0.23 * (testjd - 2451545.0);
//          fprintf(outfp, "diff=%g, r[%d]=%g\n",
//                 testdel, testncoord-1, testr[testncoord-1]);
//          fprintf(outfp, "\n\n");
                    outline++;
                    if ((outline % 100) == 0) {
//             fprintf(outfp, "%s", testline);
                        fprintf(outfp,
                                "%6d %9.1f %2d %2d %1d %17.13f %17.13f %17.13f\n",
                                outline, testjd, testntarg, testnctr,
                                testncoord, testxi, testr[testncoord - 1],
                                testdel);
                    }
                    if (testdel >= 1.0e-13) { /* error is too large */
                        fprintf(stderr,
                                "**** WARNING : miscompare at line %d, JD %9.1f\n",
                                outline, testjd);
                        fprintf(outfp,
                                "**** WARNING : next difference >= 1.e-13 ****\n");
                        fprintf(outfp,
                                "%6d %9.1f %2d %2d %1d %16.13f %16.13f %16.13f\n",
                                outline, testjd, testntarg, testnctr,
                                testncoord, testxi, testr[testncoord - 1],
                                testdel);
                    }
                } else {
                    fprintf(outfp, "Coordinates not found for\n");
                    fprintf(outfp, "%6d %9.1f %2d %2d %1d %16.13f\n", outline,
                            testjd, testntarg, testnctr, testncoord, testxi);
                }
            }
        }
    }

    fclose(infp);
    fclose(testfp);
    if (outfp != stdout)
        fclose(outfp);

    exit(0);
}

