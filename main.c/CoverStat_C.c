#include "SurveySim.h"
#include "CoverAreaStat.h"

void CoordinateSpin_x0(double coorO[3], double coorR[3], double angle ) {

    double arcAngle = angle * PI_180;
    coorR[0] = coorO[0];
    coorR[1] = coorO[1] * cos(arcAngle) + coorO[2] * sin(arcAngle);
    coorR[2] = coorO[1] * (-sin(arcAngle)) + coorO[2] * cos(arcAngle);
}

void Cartesian2Equatorial0(double* carCoor, double* eCoor) {
    double x1 = carCoor[0], x2 = carCoor[1], x3 = carCoor[2];
    double r = sqrt(x1*x1+x2*x2+x3*x3);
    double theta = asin(x3/r);

    *(eCoor+1) = theta*180*M_1_PI;
    *(eCoor+0) = atan(x2/(r*cos(theta)+x1)) *360*M_1_PI;

	*eCoor += (*eCoor < 0)*360;
}


double galactic_b(double ra, double dec){
    double phi  = ra*M_PI/180;
    double theta= dec*M_PI/180;
    double t_coor[3] = {cos(theta)*cos(phi), cos(theta)*sin(phi), sin(theta)};  //	黄道坐标系里的xyz
    double r_coor[3], ll_coor[3];
    CoordinateSpin_x0(t_coor, r_coor, -23.4522);
    Cartesian2Equatorial0(r_coor, ll_coor);

    double b = asin( -0.8676660 * cos(ll_coor[0] * PI_180) * cos(ll_coor[1] * PI_180)
                    - 0.1980764 * sin(ll_coor[0] * PI_180) * cos(ll_coor[1] * PI_180)
                    + 0.4559840 * sin(ll_coor[1] * PI_180) ) / PI_180;

    return b;
}


int main(int argc, char *argv[]) {
    
    char file_in[1024];
    char file_out[1024];
    char file_out_deep[1024];
    char file_out_norm[1024];
    char file_out_low[1024];

    if( argc != 5 ){
        printf("usage: %s survey_result.dat out_root E_dec_min G_dec_min\n",argv[0]);
        exit(0);
    }

    sprintf(file_in,"%s",argv[1]);
    printf("\n==> processing %s (all filter twice)\n",argv[1]);

	sprintf(file_out,"%s_area_C.txt",argv[2]);
    sprintf(file_out_deep,"%s_deep_C.txt",argv[2]);
    sprintf(file_out_norm,"%s_norm_C.txt",argv[2]);
    sprintf(file_out_low,"%s_low_C.txt",argv[2]);


    double E_dec_min = atof(argv[3]);
    double G_dec_min = atof(argv[4]);

    double ccd_pos_in_focus[18][2];
    init_ccd_pos(ccd_pos_in_focus);
    
    int layer = 10;
    struct FourTree* Trees12 = (struct FourTree*)malloc(sizeof(struct FourTree)*12);
    double deltArea = 41253.0/(12.0*pow(4,layer));

    printf("==> start growing the tree \n");
    produceHealPixeTrees(Trees12, layer);
    printf("==> tree is ready! \n");
    
    double t=-999,ra=-999,dec=-999,isdeep=-999;
    char line[1024];
    char start[3] = "#";

    FILE *fp_in  = fopen(file_in,"r");
    FILE *fp_out = fopen(file_out,"w");
    FILE *fp_out_deep = fopen(file_out_deep,"w");
    FILE *fp_out_norm = fopen(file_out_norm,"w");
    FILE *fp_out_low = fopen(file_out_low,"w");

    double area1=0, area2=0;
    double b;
    double area1_old = area1;
    double area2_old = area2;
    int in_deep=-1;
    int in_low=-1;

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

        area1_old = area1;
        area2_old = area2;

        MarkObervePoint_C(  dec,
                            ra,
                            isdeep,
                            Trees12,
                            ccd_pos_in_focus,
                            &area1,
                            &area2,
                            deltArea );

        fprintf(fp_out,"%15.5f %10.6f %10.6f %10.6f %10.6f %f\n",
                t,ra,dec,area1,area2,isdeep);

        b = galactic_b(ra,dec);

        in_low = -1;
        if( fabs(b)<G_dec_min || fabs(dec)<E_dec_min ){
            in_low = 1;
        }

        in_deep = (int)(isdeep > 0);

        if( in_deep != 1){
            if( area1 > area1_old ){
                if( in_low != 1 )
                    fprintf(fp_out_norm,"%15.8f %15.8f %15.8f %15.8f\n",t, ra, dec, area1);
                else
                    fprintf(fp_out_low,"%15.8f %15.8f %15.8f %15.8f\n",t, ra,dec,area1);
            }
        }
        else{
            if( area2 > area2_old ){
                fprintf(fp_out_deep,"%15.8f %15.8f %15.8f %15.8f\n",t, ra,dec,area2);
            }
        }
    }

    printf("Current Time: %8.6f yr\n",(t-2459766.)/365.25);
    printf("Normal area: %8.2f\n",area1);
    printf("U-Deep area: %8.2f\n",area2);
    
    fclose(fp_in);
    fclose(fp_out);
    fclose(fp_out_deep);
    fclose(fp_out_norm);
    fclose(fp_out_low);

    freeTrees(Trees12);
    return 0;
}
