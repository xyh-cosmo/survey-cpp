#include "SurveySim.h"
#include "CoverAreaStat.h"

int main(int argc, char *argv[]) {

    double ccd_pos_in_focus[18][2];
    init_ccd_pos(ccd_pos_in_focus);
    
    int layer = 10;
    struct FourTree* Trees12 = (struct FourTree*)malloc(sizeof(struct FourTree)*12);
    double deltArea = 41253.0/(12.0*pow(4,layer));

    printf("start growing the tree ...\n");
    produceHealPixeTrees(Trees12, layer);
    printf("the tree is ready!\n");
    
    double t,ra,dec,isdeep;
    char line[1024];
    char start[3] = "#";

    FILE *fp_in  = fopen("debug.dat","r");
    FILE *fp_out = fopen("debug.txt","w");

    double wide_area1=0, wide_area2=0, wide_area3=0; // 正常成像面积
    double deep_area1=0, deep_area2=0, deep_area3=0; // 极深度巡天面积
    while( fgets(line,1024,fp_in) != NULL ){
        char *p;

        p = strtok(line," \n\t");
        if( p==NULL ){
            continue;
        }

        if( strcmp(p,"")==0 ){
            continue;
        }

        if( p[0]==start[0] ){
            continue;
        }

        int cNum=0;
        while( p!=NULL && cNum < 15 ){
            if( cNum == 0 )
                t = atof(p);
            if( cNum == 1 )
                dec = atof(p);
            if( cNum == 2 )
                ra = atof(p);
            if( cNum == 14 ){
                isdeep = atof(p);
            }
            cNum++;
            p = strtok(NULL," \t\n");
        }

        MarkObervePoint_A(  dec,
                            ra,
                            isdeep,
                            Trees12,
                            ccd_pos_in_focus,
                            &wide_area1,
                            &deep_area1,
                            deltArea );

        MarkObervePoint_B(  dec,
                            ra,
                            isdeep,
                            Trees12,
                            ccd_pos_in_focus,
                            &wide_area2,
                            &deep_area2,
                            deltArea );

        MarkObervePoint_C(  dec,
                            ra,
                            isdeep,
                            Trees12,
                            ccd_pos_in_focus,
                            &wide_area3,
                            &deep_area3,
                            deltArea );

        fprintf(fp_out,"%15.5f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %f\n",
                t,ra,dec,
                wide_area1,deep_area1,
                wide_area2,deep_area2,
                wide_area3,deep_area3,
                isdeep);
    }
    
    fclose(fp_in);
    fclose(fp_out);

    freeTrees(Trees12);
    return 0;
}
