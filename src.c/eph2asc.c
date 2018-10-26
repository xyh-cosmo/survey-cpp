/*
   eph2asc - program to convert JPL binary ephemerides to ASCII format.
*/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "ephcom.h"


main(int argc, char* argv[], int pacth) {

    struct ephcom_Header header1;
    double *datablock;  /* Will hold coefficients from a data block */
    double *prevblock;  /* Will hold coefficients from a data block */
    int inblocks;       /* Data block number we're on in input file */
    int outblocks;      /* Data block number we're on in output file */
    int i;
    int same;
    /*
       Output file parameters.
    */
    FILE *infp, *outfp, *fopen();
    double startjd, stopjd;

    if (argc < 4) {
        fprintf(stderr,
                "\nFormat:\n\n         %s binary-file ascii-header ascii-data [startJD [stopJD]]\n\n",
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
    if (strcmp(argv[2],"-") == 0 || strcmp(argv[3],"-") == 0) {
        fprintf(stderr,"\nERROR: can't write ASCII ephemeris to stdout.\n\n");
        exit(1);
    }
    if ((outfp = fopen(argv[2], "r")) != NULL) {
        fprintf(stderr,"\nERROR: output header file %s already exists.\n\n",
                argv[2]);
        exit(1);
    }
    if ((outfp = fopen(argv[3], "r")) != NULL) {
        fprintf(stderr,"\nERROR: output data file %s already exists.\n\n",
                argv[3]);
        exit(1);
    }

    if (argc > 4) {
        startjd = atof(argv[4]);
        if (argc > 5)
            stopjd = atof(argv[5]);
        else stopjd = EPHCOM_MAXJD;
    }
    else {
        startjd = EPHCOM_MINJD;
        stopjd = EPHCOM_MAXJD;
    }

    /*
       Make sure ephemeris is within the desired range.
    */
    ephcom_readbinary_header(infp, &header1);
    if (header1.ss[0] > stopjd) {
        fprintf(stderr,
                "\nERROR: ephemeris begins (%15.2f) after desired end date (%15.2f).\n\n",
                header1.ss[0], stopjd);
        exit(1);
    }
    if (header1.ss[1] < startjd) {
        fprintf(stderr,
                "\nERROR: ephemeris ends (%15.2f) before desired start date (%15.2f).\n\n",
                header1.ss[1], startjd);
        fprintf(stderr, "\nERROR: StopJD is before ephemeris begins.\n\n");
        exit(1);
    }

    if ((outfp = fopen(argv[3],"w")) == NULL) {
        fprintf(stderr,"Can't open %s for output.\n\n", argv[3]);
        exit(1);
    }
// fprintf(stderr, "%d coefficients per block.\n", header1.ncoeff);
    /*
       Done with header.  Now we'll read and write data blocks.
    */
    datablock = (double *)malloc(header1.ncoeff * sizeof(double));
    prevblock = (double *)malloc(header1.ncoeff * sizeof(double));
    outblocks = inblocks = 0;
    while (ephcom_readbinary_block(infp, &header1, inblocks, datablock) > 0 &&
            datablock[0] <= stopjd) {
        same = 1; /* Assume this and previous block are the same */
        for (i=0; same && i<header1.ncoeff; i++)
            if (prevblock[i] != datablock[i]) same = 0; /* Found difference */
        if (same) {
            fprintf(stderr,
                    "WARNING: Ignoring data block %d - identical to block %d.\n",
                    inblocks+1, inblocks);
        }
        else if (datablock[1] >= startjd) {
            ephcom_writeascii_block(outfp, &header1, outblocks, datablock);
            if (outblocks == 0) header1.ss[0] = datablock[0]; /* Get real start */
            header1.ss[1] = datablock[1];
            outblocks++;
        }
        memcpy(prevblock, datablock, sizeof(double) * header1.ncoeff);
// for (i=0; i<header1.ncoeff; i++) prevblock[i] = datablock[i];
        inblocks++; /* Update even if block is duplicate, so we read next block */
    }
    fclose(infp);
    fclose(outfp);
    if ((outfp = fopen(argv[2],"w")) == NULL) {
        fprintf(stderr,"\nERROR: Can't open %s for output.\n\n", argv[2]);
        exit(1);
    }
    ephcom_writeascii_header(outfp, &header1);

    exit(0);
}

