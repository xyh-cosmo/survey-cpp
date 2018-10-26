/*
   headcomp - program to compare headers of two JPL ephemeris files

*/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "ephcom.h"


main(int argc, char* argv[], int pacth) {

    struct ephcom_Header header1, header2;
    int nblocks1, nblocks2;
    int i, j;
    unsigned char inchar1[8], inchar2[8]; /* Two strings, for compare_ascii() */
    double x, y;                          /* Two doubles, for compare_double() */
    int int1, int2;                       /* Two ints, for compare_int() */
    char outstring[86];                   /* To write label for line */
    double ephcom_indouble(FILE *);
    int ephcom_inint(FILE *);
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

    if (strcmp(argv[1], "-") == 0 || strcmp(argv[2], "-") == 0) {
        fprintf(stderr,"\nERROR: Can't read ephemeris file from stdin.\n\n");
        exit(1);
    }
    if ((infp1 = fopen(argv[1],"r")) == NULL) {
        fprintf(stderr,"\nERROR: Can't open ephemeris file %s for input.\n\n", argv[1]);
        exit(1);
    }
    if ((infp2 = fopen(argv[2],"r")) == NULL) {
        fprintf(stderr,"\nERROR: Can't open ephemeris file %s for input.\n\n", argv[2]);
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
            fprintf(stderr,"\nERROR: output file %s alredy exists.\n\n", argv[3]);
            exit(1);
        }
    }
    else outfp = stdout;


    ephcom_readbinary_header(infp1, &header1);
    ephcom_readbinary_header(infp2, &header2);
    fprintf(outfp, "Header1: %d coefficients per block.\n", header1.ncoeff);
    fprintf(outfp, "Header2: %d coefficients per block.\n", header2.ncoeff);
    fprintf(outfp, "\n\n");
    if (header1.ncoeff != header2.ncoeff) {
        fprintf(stderr,"Different ephemeris files; different coefficients.\n");
        exit(1);
    }

    rewind(infp1);
    rewind(infp2);

    /*
       Go through first two blocks of data (the header).
       One block is 8*ncoeff bytes.
    */
    strcpy(outstring, "  ");
    for (i=0; i < 2648 /* 16*header1.ncoeff */ ; i += 8) {
        for (j=0; j<8; j++) inchar1[j] = fgetc(infp1);
        for (j=0; j<8; j++) inchar2[j] = fgetc(infp2);
        compare_ascii(outfp, i, outstring, inchar1, inchar2, 8);
    }
    for (j=0; j<4; j++) inchar1[j] = fgetc(infp1);
    for (j=0; j<4; j++) inchar2[j] = fgetc(infp2);
    strcpy(outstring, "                  ");
    compare_ascii(outfp, i, outstring, inchar1, inchar2, 4);
    i += 4;
    x = ephcom_indouble(infp1);
    y = ephcom_indouble(infp2);
    compare_double(outfp, i, "StartJD   ", x, y);
    i += 8;
    x = ephcom_indouble(infp1);
    y = ephcom_indouble(infp2);
    compare_double(outfp, i, "EndJD     ", x, y);
    i += 8;
    x = ephcom_indouble(infp1);
    y = ephcom_indouble(infp2);
    compare_double(outfp, i, "BlockSpan ", x, y);
    i += 8;
    int1 = ephcom_inint(infp1);
    int2 = ephcom_inint(infp2);
    compare_int(outfp, i, "ncon      ", int1, int2);
    i += 4;
    x = ephcom_indouble(infp1);
    y = ephcom_indouble(infp2);
    compare_double(outfp, i, "au        ", x, y);
    i += 8;
    x = ephcom_indouble(infp1);
    y = ephcom_indouble(infp2);
    compare_double(outfp, i, "emrat     ", x, y);
    i += 8;
    for (j=0; j<36; j++) { /* ipt values */
        int1 = ephcom_inint(infp1);
        int2 = ephcom_inint(infp2);
        sprintf(outstring, "ipt[%02d]   ", j);
        compare_int(outfp, i, outstring, int1, int2);
        i += 4;
    }
    int1 = ephcom_inint(infp1);
    int2 = ephcom_inint(infp2);
    compare_int(outfp, i, "numde     ", int1, int2);
    i += 4;
    for (j=0; j<3; j++) { /* ipt values */
        int1 = ephcom_inint(infp1);
        int2 = ephcom_inint(infp2);
        sprintf(outstring, "lpt[%02d]   ", j);
        compare_int(outfp, i, outstring, int1, int2);
        i += 4;
    }
    /*
       End of first block can just be blank (but it's not in JPL DE406).
    */
    strcpy(outstring, "# ");
    for ( ; i < 8 * header1.ncoeff - 7; i += 8) {
        for (j=0; j<8; j++) inchar1[j] = fgetc(infp1);
        for (j=0; j<8; j++) inchar2[j] = fgetc(infp2);
        compare_ascii(outfp, i, outstring, inchar1, inchar2, 8);
    }
    for (j=0; j<400; j++) {
        x = ephcom_indouble(infp1);
        y = ephcom_indouble(infp2);
        sprintf(outstring, "val[%03d]  ", j);
        compare_double(outfp, i, outstring, x, y);
        i += 8;
    }
    /*
       End of second block can just be blank (but it's not in JPL DE406).
    */
    strcpy(outstring, "# ");
    for ( ; i < 16 * header1.ncoeff; i += 8) { /* print to end of 2nd block */
        for (j=0; j<8; j++) inchar1[j] = fgetc(infp1);
        for (j=0; j<8; j++) inchar2[j] = fgetc(infp2);
        compare_ascii(outfp, i, outstring, inchar1, inchar2, 8);
    }

    fclose(infp1);
    fclose(infp2);
    if (outfp != stdout) fclose(outfp);

    exit(0);
}

/*
   Compare two ASCII strings (not null terminated).
*/
compare_ascii(FILE *outfp, int position, char *name,
              unsigned char *c1, unsigned char *c2, int length) {

    int i;
    int same;

    fprintf(outfp, "%06d", position);
    same = 1; /* Assume two strings are the same */
    for (i=0; i<length && c1[i] == c2[i]; i++);
    if (i<8 && c1[i] != c2[i]) same=0;
    if (same)
        fputc(' ', outfp);
    else
        fputc('~', outfp);

    fprintf(outfp, "%s ", name); /* gcc wouldn't print this in previous printf! */
    for (i=0; i<length; i++) fputc(isprint(c1[i]) ? c1[i] : ' ', outfp);
    for (i=0; i<length; i++) fprintf(outfp, " %02x", c1[i]);
    if (same)
        fprintf(outfp, "  :  ");
    else
        fprintf(outfp, " <:> ");
    for (i=0; i<length; i++) fputc(isprint(c2[i]) ? c2[i] : ' ', outfp);
    for (i=0; i<length; i++) fprintf(outfp, " %02x", c2[i]);
    fputc('\n', outfp);
    return(0);
}


/*
   Compare two integer values.
*/
compare_int(FILE *outfp, int position, char *name, int int1, int int2) {

    fprintf(outfp, "%06d%c%s ", position, int1 == int2 ? ' ' : '~', name);
    fprintf(outfp, "%24d", int1);
    if (int1 == int2)
        fprintf(outfp, "  :  ");
    else
        fprintf(outfp, " <:> ");
    fprintf(outfp, "%24d\n", int2);
    return(0);
}

/*
   Compare two double precision floating point values.
*/
compare_double(FILE *outfp, int position, char *name, double d1, double d2) {

    fprintf(outfp, "%06d%c%s ", position, d1 == d2 ? ' ' : '~', name);
    fprintf(outfp, "%24.17e", d1);
    if (d1 == d2)
        fprintf(outfp, "  :  ");
    else
        fprintf(outfp, " <:> ");
    fprintf(outfp, "%24.17e\n", d2);
    return(0);
}
