#include "SurveySim.h"
#define MAXLINE 400


void readData(char* fileName, int lineNum,int elemNum, double *da) {
    FILE *f;
    f = fopen(fileName,"r");

    if(f == NULL) {
        printf("No such file !!!\n");
        exit(1);
    }
    int num = 0;
    char line[MAXLINE];
    while (fgets(line, MAXLINE, f) != NULL && num < lineNum) {
        char *p;
        p = strtok(line, " \n");
        int cNum = 0;
        while (p != NULL) {
            double cc = atof(p);
            *(da + num * elemNum + cNum) = cc;
            cNum++;
            p = strtok(NULL, " \n");
        }
        num ++;
    }
}

int main(int argc, char *argv[]) {
    char* fileName = "SS_plan1.dat";
    char* infile = "jpl.405";
    char* fileResult = "spectral";

    int cNum = 19;
    int rNum = 226526;//1184175;
    double *data = (double*)malloc(sizeof(double)*rNum*cNum);
    readData(fileName,rNum,  cNum, data);
    FILE *infp = fopen(infile,"r");
    if(infp == NULL) {
        printf("no this ephemeris file \n");
        exit(1);
    }
    FILE *rf = fopen(fileResult,"w");
    if(rf == NULL) {
        printf("no result file \n");
        exit(1);
    }
    double startTime = 2459671.5;
    int i = 0;
    for(i = 0; i < rNum; i ++) {
        double time = data[cNum*i];
        //printf("time is %f \n",time);
        double sat[3], sun[3], moon[3];
        double altitude = EARTH_RADIUS + SAT_ATITUDE;
        double obPeriod = ORBIT_PERIOD * 60;
        double presession = PRESESSION_PERIOD * 24 * 60 * 60;
        double durTime = 0; //, curTime = currentTime;
        durTime = (time - startTime) * 86400.0;
        locateSat(sat, RAAN, altitude, INCLINATION, obPeriod, presession, durTime);
        locateSun(infp, time, sun);
        locateMoon(infp, time, moon);
        fprintf(rf,"%f    %f    %f    %f    %f    %f    %f    %f    %f    %f\n",time,sun[0],sun[1],sun[2],moon[0],moon[1],moon[2],sat[0],sat[1],sat[2]);
    }
    fclose(infp);
    fclose(rf);

}
