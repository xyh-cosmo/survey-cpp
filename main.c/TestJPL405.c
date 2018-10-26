//  测试获取太阳、月球的位置

#include "SurveySim.h"

int main(int argc, char *argv[]) {

    infp = NULL;

    if (argc < 3) {
        fprintf(stderr, "usage:   %s jpl405_file julian_date\n", argv[0]);
        return 0;
    }

    char jpl_fileName[1024];
    sprintf(jpl_fileName,"%s",argv[1]);

    if ((infp = fopen(jpl_fileName, "r")) == NULL ) {
        fprintf(stderr, "\nERROR: Can't open ephemeris file %s for input.\n\n", argv[1]);
    }

    double jdate = atof(argv[2]);

    double sun[3],moon[3];

    locateSun(infp, jdate, sun);
    locateMoon(infp, jdate, moon);

    int i;

    printf("Sun:\n");
    for( i=0; i<3; i++ )
        printf("%18.15g ",sun[i]);
    printf("\n");

    printf("Moon:\n");
    for( i=0; i<3; i++ )
        printf("%18.15g ",moon[i]);
    printf("\n");

    if( infp != NULL )
        free(infp);

    return 0;
}
