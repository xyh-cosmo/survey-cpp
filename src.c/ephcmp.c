/*
   ephcmp - program to compare all header constants and data block coefficients
            in two JPL ephemeris files.

*/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "ephcom.h"


main(int argc, char* argv[], int pacth) {

    struct ephcom_Header header1, header2;
    double *datablock1, *datablock2;  /* Will hold coefficients from a data block */
    int nblocks1, nblocks2;
    int i;
    int contotal = 0;       /* Number of constants in header */
    double conerror = 0.0;  /* Error from constants in header */
    int ntotal = 0;         /* Number of coefficients in data blocks */
    int miscompares=0;      /* Number of times we found a difference */
    double xerror = 0.0;    /* Error from coefficients in data blocks */
    double x;
    double fabs(double);
    /*
       Output file parameters.
    */
    FILE *infp1, *infp2, *outfp, *fopen();


    if (argc < 3) {
        fprintf(stderr,
                "\nFormat:\n\n         %s first-ephemeris second-ephemeris [output-file]\n\n",
                argv[0]);
        exit(1);
    }

    if (strcmp(argv[1], "-") == 0 || strcmp(argv[2],"-") == 0) {
        fprintf(stderr,"\nERROR: Can't read ephemeris file from stdin.\n\n");
        exit(1);
    }
    if ((infp1 = fopen(argv[1],"r")) == NULL) {
        fprintf(stderr,"\nERROR: Can't open binary ephemeris file %s for input.\n\n", argv[1]);
        exit(1);
    }
    if ((infp2 = fopen(argv[2],"r")) == NULL) {
        fprintf(stderr,"\nERROR: Can't open binary ephemeris file %s for input.\n\n", argv[2]);
        exit(1);
    }
    if (argc >= 4) {
        if (strcmp(argv[3], "-") == 0) {
            outfp = stdout;
        }
        else if ((outfp = fopen(argv[3],"r")) == NULL) { /* File doesn't exist */
            if ((outfp = fopen(argv[3],"w")) == NULL) {
                fprintf(stderr,"\nERROR: Can't open %s for output.\n\n", argv[3]);
                exit(1);
            }
        }
        else {
            fprintf(stderr,"\nERROR: Output file %s already exists.\n\n", argv[3]);
            exit(1);
        }
    }
    else outfp = stdout;


    ephcom_readbinary_header(infp1, &header1);
    ephcom_readbinary_header(infp2, &header2);
    printf("Binary Ephemeris Comparison.\n");
    printf("   Header1: DE%03d/LE%03d, %d constants and %d coefficients per block.\n",
           header1.numde, header1.numle, header1.ncon, header1.ncoeff);
    printf("   Header2: DE%03d/LE%03d, %d constants and %d coefficients per block.\n",
           header2.numde, header2.numle, header2.ncon, header2.ncoeff);
    if (header1.ncoeff != header2.ncoeff) {
        fprintf(stderr,"Different ephemeris files; different coefficients.\n");
        exit(1);
    }
    for (i=0; i<header1.ncon; i++) {
        if (fabs(header1.cval[i]) != 0.0) {
            x = fabs((header1.cval[i] - header2.cval[i]) / header1.cval[i]);
            conerror += x;
        }
        if (header1.cval[i] != header2.cval[i]) miscompares++;
        contotal++;
    }

    /*
       Done with header.  Now we'll read and write data blocks.
    */
    datablock1 = (double *)malloc(header1.ncoeff * sizeof(double));
    datablock2 = (double *)malloc(header2.ncoeff * sizeof(double));
    nblocks1 = 0;
    nblocks2 = 0;
    while (ephcom_readbinary_block(infp1, &header1, nblocks1, datablock1) > 0 &&
            ephcom_readbinary_block(infp2, &header2, nblocks2, datablock2) > 0) {
        for (i=0; i<header1.ncoeff; i++) {
            if (fabs(datablock1[i]) != 0.0) {
                x = fabs((datablock2[i] - datablock1[i]) / datablock1[i]);
                xerror += x;
            }
            if (datablock1[i] != datablock2[i]) miscompares++;
            ntotal++;
        }
        nblocks1++;
        nblocks2++;
    }
    fclose(infp1);
    fclose(infp2);
    printf("   Total of %d coefficients compared in %d blocks\n", ntotal, nblocks1);
    printf("   Average fractional constant error: %23.18e\n", conerror / contotal);
    printf("   Average fractional coefficient error: %23.18e\n", xerror / ntotal);
    printf("   Total miscompares: %d\n", miscompares);

    if (outfp != stdout) fclose(outfp);

    exit(0);
}

