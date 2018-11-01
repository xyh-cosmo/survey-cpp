#include <iostream>
#include <stdlib.h>
#include "../src/iniparser.h"

using namespace std;

int main( int argc, char *argv[] ){

    if( argc < 2 ){
        cout << "usage: " << argv[0] << "  *.ini\n";
        return 0;
    }

    dictionary *ini;

    ini = iniparser_load(argv[1]);
    if( ini == NULL ){
        cout << "cannot parse file: " << argv[1] << "\n";
        return 0;
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
