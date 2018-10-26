#include "SurveySim.h"

int main( int argc, char* argv[] ){
	
	if( argc < 3 ){
		printf("usage: %s simulation_result.txt out.txt\n",argv[0]);
		exit(0);
	}
	
	char file_in[1024];
	char file_out[1024];
	
	sprintf(file_in,"%s",argv[1]);
	sprintf(file_out,"%s",argv[2]);

    FILE *fp_in  = fopen(file_in,"r");
    FILE *fp_out = fopen(file_out,"w");

    double tAngle=-999, cmg_use, cmg_total=0;
    double t,t_min=-999;
    int cnt_orbit = 0;
    
    char line[1024];
    char start[3] = "#";
    
/*    int cnt_debug=0;*/
    int cnt_err=0;
    int cnt_empty_orbit=0;

    double first_cmg_use=-999, last_cmg_use=0;
    
    while( fgets(line,1024,fp_in) != NULL ){

        char *p;

        p = strtok(line," \n\t");
        if( p==NULL ) continue;
        if( strcmp(p,"")==0 ) continue;
        if( p[0]==start[0] ) continue;

        int cNum=0;
        while( p!=NULL && cNum < 18 ){
            if( cNum == 0 ){
                t = (atof(p)-2459766);

                if( t_min <= -999 )
	                t_min = t;
            }
            
            if( cNum == 17 ){
                tAngle = atof(p);
            }
            
            cNum++;
            p = strtok(NULL," \t\n");
        }

        
        
        if( t <= t_min + (cnt_orbit+1)*0.0638 ){
        	cmg_use = get_cmg_use(tAngle);

            if( first_cmg_use < 0 ){
                first_cmg_use = cmg_use;
                // first_cmg_use = tAngle;
            }

            last_cmg_use = cmg_use;
            
            cmg_total += cmg_use;
        	if( cmg_total > 1.0 ){
        		cnt_err++;
        	}
        }
        else{
        	if( cmg_total > 0 && first_cmg_use >= 0){
	        	fprintf(fp_out, "%12d  %10.6f  %10.6f  %10.6f\n", 
                        cnt_orbit+1, 
                        cmg_total,
                        first_cmg_use,
                        last_cmg_use);
                first_cmg_use = -1;
                last_cmg_use = 0;
            }
	        else
	        	cnt_empty_orbit++;
        	
        	cnt_orbit++;
        	// cmg_total = 0;
			cmg_total = get_cmg_use(tAngle);
        }
        
/*        if( cnt_debug > 100 ) break;*/
/*        cnt_debug++;*/
	}
	
	printf("--> found %5d CMG errors in total %5d orbits!\n",cnt_err,cnt_orbit);
	printf("--> found %5d empty orbit.\n",cnt_empty_orbit);
	
	fclose(fp_in);
	fclose(fp_out);
}
