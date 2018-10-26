/*
   eph2eph - program to read part or all of two binary JPL ephemeris files
             into an output binary ephemeris file.  Also allows the two
             input ephemerides to be the same, in which case it will
             clean up the header portion and write a new binary file.

         As a safety precaution, eph2eph will not write to the output file
         if it already exists.
*/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "ephcom.h"


main(int argc, char* argv[], int pacth) {

    struct ephcom_Header header[3];
    double *datablock;  /* Will hold coefficients from a data block */
    int inblock;        /* Data block number in input file, 0.. */
    int outblock;       /* Data block number in output file, 0.. */
    int totalblocks;    /* Number of blocks in current ephemeris input file */
    int bytesread;      /* Number of bytes read in binary ephemeris data block */
    int i;
    int init;
    int size1, size2;   /* Check for bugs by computing file size two ways */
    /*
       Output file parameters.
    */
    FILE *infp[2], *outfp, *readfp, *fopen();
    int currentfile; /* Holds 0 or 1, index into infp[] */
    double startjd, stopjd, prevjd;
    void *memcpy();

    if (argc < 4) {
        fprintf(stderr,
                "\nFormat:\n\n         %s input-file-1 input-file-2 output-file [startJD [stopJD]]\n\n",
                argv[0]);
        exit(1);
    }

    /*
       If output file exists, stop without destroying it.
    */
    if (strcmp(argv[1], "-") == 0 || strcmp(argv[2], "-") == 0) {
        fprintf(stderr,"\nERROR: Can't read binary ephemeris from stdin.\n\n");
        exit(1);
    }
    if (strcmp(argv[3], "-") == 0) {
        fprintf(stderr,"\nERROR: Can't write binary ephemeris to stdout.\n\n");
        exit(1);
    }
    if ((infp[0] = fopen(argv[1],"r")) == NULL) {
        fprintf(stderr,"\nERROR: Can't open %s for input.\n", argv[1]);
        fprintf(stderr,"Can't continue.\n\n");
        exit(1);
    }
    /*
       If the file names are the same, user just wanted to clean up the
       data in the header portion of the file.  It's okay.
    */
    if (strcmp(argv[1], argv[2]) == 0) infp[1] = infp[0];
    else if ((infp[1] = fopen(argv[2],"r")) == NULL) {
        fprintf(stderr,"\nERROR: Can't open %s for input.\n", argv[2]);
        fprintf(stderr,"Can't continue.\n\n");
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
       Make sure first ephemeris is within the desired range.
    */
    ephcom_readbinary_header(infp[0], &header[0]);
    if (header[0].ss[0] > stopjd) {
        fprintf(stderr,
                "\nERROR: ephemeris %s begins (%15.2f)\n", argv[1], header[0].ss[0]);
        fprintf(stderr,
                "       after desired end date (%15.2f).\n\n", stopjd);
        exit(1);
    }
    if (header[0].ss[1] < startjd) {
        fprintf(stderr,
                "\nERROR: ephemeris %s ends (%15.2f)\n", argv[1], header[0].ss[1]);
        fprintf(stderr,
                "       before desired start date (%15.2f).\n\n", startjd);
        exit(1);
    }

    /*
       Make sure second ephemeris is within the desired range and that format
       is compatible with first file (same DE number).
    */
    if (infp[1] == infp[0]) {
        memcpy(&header[1], &header[0], sizeof(struct ephcom_Header));
    }
    else {
        ephcom_readbinary_header(infp[1], &header[1]);
        if (header[1].ss[0] > stopjd) {
            fprintf(stderr,
                    "\nERROR: ephemeris %s begins (%15.2f)\n", argv[1], header[1].ss[0]);
            fprintf(stderr,
                    "       after desired end date (%15.2f).\n\n", stopjd);
            exit(1);
        }
        if (header[1].ss[1] < startjd) {
            fprintf(stderr,
                    "\nERROR: ephemeris %s ends (%15.2f)\n", argv[1], header[1].ss[1]);
            fprintf(stderr,
                    "       before desired start date (%15.2f).\n\n", startjd);
            exit(1);
        }

// printf("File 1: %d coefficients; File 2: %d coefficients.\n",
//        header[0].ncoeff, header[1].ncoeff);
        if (header[0].numde != header[1].numde) {
            fprintf(stderr, "\nERROR:\n");
            fprintf(stderr, "      %s is in ephemeris DE%d format.\n",
                    argv[1], header[0].numde);
            fprintf(stderr, "      %s is in ephemeris DE%d format.\n",
                    argv[2], header[1].numde);
            fprintf(stderr, "Can't continue.\n\n");
            exit(1);
        }
        if (header[0].ncoeff != header[1].ncoeff) {
            fprintf(stderr,
                    "\nERROR: input files have different number of coefficients.\n\n");
            exit(1);
        }
        if (header[0].emrat != header[1].emrat) {
            fprintf(stderr,
                    "\nERROR: input files have different Earth-Moon mass ratios.\n\n");
            exit(1);
        }
        if (header[0].au != header[1].au) {
            fprintf(stderr,
                    "\nERROR: input files have different values for an AU.\n\n");
            exit(1);
        }
    }
    /*
       If second input file begins before first one, swap all values.
    */
    if (header[1].ss[0] < header[0].ss[0]) {
        memcpy(&header[2], &header[0], sizeof(struct ephcom_Header));
        memcpy(&header[0], &header[1], sizeof(struct ephcom_Header));
        memcpy(&header[1], &header[2], sizeof(struct ephcom_Header));
        readfp = infp[0];
        infp[0] = infp[1];
        infp[1] = readfp;
    }
    /*
       If the first file spans a greater range than the second file,
       ignore the second file (not the usual case, but handle it).
    */
    if (header[1].ss[1] < header[0].ss[1]) {
        memcpy(&header[1], &header[0], sizeof(struct ephcom_Header));
        fclose(infp[1]);
        infp[1] = infp[0];
    }
    /*
       If the two files are different and they don't overlap in time,
       we can't merge them.  Stop.
    */
    if (header[0].ss[1] < header[1].ss[0]) {
        fprintf(stderr,"\nERROR: There is a time gap between %s and %s --\n",
                argv[1], argv[2]);
        fprintf(stderr,"       Can't merge them.\n\n");
        exit(1);
    }
    /*
       The header files are compatible.  We can start to write output.
    */
    outfp = fopen(argv[3], "r"); /* File should not exist; this should fail. */
    if (outfp == NULL) {
        outfp = fopen(argv[3], "w");  /* What we want to happen */
    }
    else {
        fclose(outfp);
        fprintf(stderr, "\nERROR: Output file %s already exists; can't continue.\n\n",
                argv[3]);
        exit(1);
    }
    /*
       Initialize output parameters in header3 struct.
    */
    memcpy(&header[2], &header[0], sizeof(struct ephcom_Header));
    header[2].ss[1] = header[1].ss[1]; /* later end date */
    strncpy(header[2].ttl[2], header[1].ttl[2], 84); /* later end date */

    ephcom_writebinary_header(outfp, &header[2]);

// fprintf(stderr, "%d coefficients per block.\n", header[1].ncoeff);

    /*
       Done with header.  Now we'll read and write data blocks.
    */
    datablock = (double *)malloc(header[0].ncoeff * sizeof(double));
    currentfile = 0;
    /*
       Determine start date in output file, and start block in first file.
    */
    inblock = (startjd - header[0].ss[0]) / header[0].ss[2] + 1;
    if (inblock < 0) inblock = 0; /* If first file starts before startjd */
    header[2].ss[0] = header[0].ss[0] + inblock*header[0].ss[2];

    outblock = 0;

// printf("Allocated %d bytes for datablock, for %d coefficients.\n",
//        header[0].ncoeff * sizeof(double), header[0].ncoeff);
// printf("Starting in first file at block %d\n", inblock);
    prevjd = header[2].ss[0] - header[2].ss[2];
    for (currentfile = 0; currentfile < 2; currentfile++) {
// printf("Reading block %d...", inblock, currentfile+1);
        readfp = infp[currentfile];
        inblock = (prevjd - header[currentfile].ss[0])/header[currentfile].ss[2];
        if (inblock < 0) {
            if (currentfile == 0) inblock = 0;
            else {
                fprintf(stderr,"\nERROR: time gap in data.  Can't continue.\n\n");
                exit(1);
            }
        }
        totalblocks  = (header[currentfile].ss[1] - header[currentfile].ss[0]) /
                       header[currentfile].ss[2];
        bytesread = 0;
        do {
            if (inblock < totalblocks) {
                bytesread = ephcom_readbinary_block(readfp, &header[currentfile],
                                                    inblock, datablock);
                if (bytesread > 0 && datablock[0] < stopjd) {

//          printf("file%d: inblock %d, outblock %d\n",
//                 currentfile, inblock, outblock);
                    ephcom_writebinary_block(outfp, &header[0], outblock, datablock);
                    if (currentfile == 0 && outblock == 0)
                        header[2].ss[0] = datablock[0]; /* Guarantee start date */
                    header[2].ss[1] = datablock[1];
                    inblock++;
                    outblock++;
                }
            }
            prevjd = datablock[1];
//       printf("bytesread=%d, inblock=%d, totalblocks=%d, %12.1f < %12.1f\n",
//              bytesread, inblock, totalblocks, datablock[0], stopjd);
        } while (bytesread > 0 && inblock < totalblocks &&
                 datablock[0] < stopjd);
        if (infp[currentfile] != stdin) fclose(infp[currentfile]);
        if (infp[0] == infp[1]) currentfile++; /* Only loop once if same file */
    }
    /*
       Now that we know the real start and stop Julian Days that we wrote
       into the data blocks, update the header with the new dates.
    */
    ephcom_writebinary_header(outfp, &header[2]);
    /*
       N.B.: the following null terminations in header[].ttl[] corrupts
       the ttl lines by inserting a null.  Don't do this if you will later
       write these lines to an ASCII or binary ephemeris file directly.

       ephcom_writeascii_header() and ephcom_writebinary_header() both
       regenerate these lines before writing them, by using the Julian Days
       in header[].ss[0] and header[].ss[1].
    */
    header[2].ttl[1][49] = header[2].ttl[2][49] = '\0';

    size1 = (2 + (header[2].ss[1] - header[2].ss[0]) / header[2].ss[2]) *
            header[2].ncoeff * 8;
    size2 = (2 + outblock) * header[2].ncoeff * 8;
    printf("Days per block: %5.0f\n%s\n%s\n",
           header[2].ss[2], header[2].ttl[1], header[2].ttl[2]);
    printf("Output ephemeris %s should be exactly %d bytes\n",
           argv[3], size1);
    printf("   = (2 header blocks + %d data blocks) * %d coeff/block *",
           outblock, header[2].ncoeff);
    printf(" 8 bytes/coeff\n");
    if (size1 != size2) {
        fprintf(stderr,"\nERROR: output ephemeris should have %d bytes\n", size1);
        fprintf(stderr,"       but instead has %d bytes.  This is a bug.\n", size2);
        exit(1);
    }
    if (outfp != stdout) fclose(outfp);

    exit(0);
}
