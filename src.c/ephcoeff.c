/*
   ephcoeff - program to print the parsed coefficients from one data block
              in an ephemeris that contains the desired Julian Day.
*/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "ephcom.h"


main(int argc, const char* argv[], int pacth) {

    struct ephcom_Header header;
    double *datablock;  /* Will hold coefficients from a data block */
    double *prevblock;  /* Will hold coefficients from a data block */
    int nblocks;
    int i;
    double testjd;
    double atof(const char *);
    FILE *infp, *outfp=stdout, *fopen();

    if (argc < 3) {
        fprintf(stderr,
                "\nFormat:\n\n         %s ephemeris-file Julian-Day\n\n",
                argv[0]);
        exit(1);
    }

    if (strcmp(argv[1], "-") == 0) {
        fprintf(stderr,"\nERROR: can't read binary ephemeris file from stdin.\n\n");
        exit(1);
    }
    else if ((infp = fopen(argv[1],"r")) == NULL) {
        fprintf(stderr,"\nERROR: Can't open %s for input.\n\n", argv[1]);
        exit(1);
    }
    testjd = atof(argv[2]);

    /*
       Make sure ephemeris is within the desired range.
    */
    ephcom_readbinary_header(infp, &header);
    if (header.ss[0] > testjd) {
        fprintf(stderr,
                "\nERROR: ephemeris begins (%15.2f) after desired date (%15.2f).\n\n",
                header.ss[0], testjd);
        exit(1);
    }
    if (header.ss[1] < testjd) {
        fprintf(stderr,
                "\nERROR: ephemeris ends (%15.2f) before desired date (%15.2f).\n\n",
                header.ss[1], testjd);
        exit(1);
    }

    /*
       Done with header.  Now we'll read desired data block.
    */
    datablock = (double *)malloc(header.ncoeff * sizeof(double));
    nblocks = (testjd - header.ss[0]) / header.ss[2];
    fprintf(outfp,"JD is %f\n", testjd);
    fprintf(outfp,"Reading data block %d\n", nblocks);
    if (ephcom_readbinary_block(infp, &header, nblocks, datablock) > 0) {
        ephcom_parse_block(outfp, &header, datablock);
    }
    else {
        fprintf(stderr,"\nERROR: couldn't read block %d from ephemeris file.\n\n",
                nblocks+2);
        exit(1);
    }
    fclose(infp);
    if (outfp != stdout) fclose(outfp);

    exit(0);
}
