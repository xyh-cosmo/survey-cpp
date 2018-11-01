#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "../src/iniparser.h"

int main( int argc, char *argv[] ){

    if( argc < 2 ){
        printf("usage: %s *.ini\n",argv[0]);
        return 0;
    }

    dictionary *ini;

    ini = iniparser_load(argv[1]);
    if( ini == NULL ){
        fprintf(stderr, "cannot parse file: %s\n", argv[1]);
        return -1;
    }

    FILE *fp_dump = fopen("xxxx.ini","w");

    iniparser_dump_ini(ini, fp_dump);
    fclose(fp_dump);

    const char *ch;
    ch = iniparser_getstring(ini,"CMG:usage_max","temp");
    double cmg_usage_max = iniparser_getdouble(ini,"CMG:usage_max",-1);
    double panel_angle_max = iniparser_getdouble(ini,"Solar_Panel:max_solar_angle",-99);

    printf("cmg_usage_max (str) : %s\n", ch);
    printf("cmg_usage_max : %g\n", cmg_usage_max);
    printf("panel_angle_max : %g\n", panel_angle_max);

    iniparser_freedict(ini);

    return 0;
}
