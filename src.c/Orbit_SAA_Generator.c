/**************************************************************
*此程序主要用于判断SAA的区域以及输出坐标点的辐射质子能谱 
*此程序主要应用于需要计算流量的数据点较多且已确定的情况(统计辐射计量)
* 输入量:
      太阳活动状态：极大期为1，极小期为0  c1_sunstate=1（可由时间判断，判断函数未写，标准自定）.
      坐标值 :c2_altitude(360~420km，已经扣除6371.20km),
             c3_long(0~360deg),c4_lat(-89~89deg) .
      地球表面高度设置为6371.20km


           
*输出量：
         质子能谱en_flue[10] (10～100Mev),间隔10Mev

*注：目前使用的空间高度间隔为20km，辐射场的空间分辨率为0.1deg
***************************************************************/
#include <stdlib.h>
#include<stdio.h>
#include<math.h>     
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

#define FILE_NUM 10890

// max:360,380,400,420
//     360: offset- 0
//     380: offset- 90*121
//     400: offset- 90*121*2
//     420: offset- 90*121*3
// min is same as max
// aquire flue >-10MeV data

void aquireFlueData(float* flue_max, float* flue_min){
  int hf,en;
  double  c1_B,c2_L,c3_FLUE;
  for (hf=0;hf<=3;hf++){
    // char str1[1024]="SAA_data/AP_MAX_10.txt";
    char str1[1024]="SAA_data/AP_MAX_";
    char str2[10];
    int height_SAA=20*hf+360;
    sprintf(str2,"%d",height_SAA);
    char str3[]="KM_10.txt";
    strcat(str1,str2);
    strcat(str1,str3);
   
    char str1s[1024]="SAA_data/AP_MIN_";
    char str2s[10];
    int height_SAAs=20*hf+360;
    sprintf(str2s,"%d",height_SAAs);
    char str3s[]="KM_10.txt";
    strcat(str1s,str2s);
    strcat(str1s,str3s);
    
    //数据读取
    FILE *fp1;
    char buf[1024];                                
    if( (fp1=fopen(str1,"r"))==NULL ){
        printf("打开输入文件失败，可能文件%s还没有创建!\n",str1);
        exit(0);
    }

    int nn=0;
    while(!feof(fp1)){
      fscanf(fp1,"%le,%le,%le",&c1_B,&c2_L,&c3_FLUE);
      memset(buf, 0, 1024); 
      fgets(buf,1024,fp1);
      if(c3_FLUE<0){
        c3_FLUE=0;
      }
      // long_grid=(nn%121);
      // lat_grid=floor(nn/121);

      flue_max[FILE_NUM*hf + nn] =c3_FLUE;

      nn ++; 

      if(nn >= FILE_NUM){
        break;
      }
    }
    fclose(fp1);

    FILE *fp1s;
    char bufs[1024];
    if( (fp1s=fopen(str1s,"r"))==NULL ){
      printf("打开输入文件失败，可能文件%s还没有创建!\n",str1s);
      exit(0);
    }
    int nns=0;
    while(!feof(fp1s)){
      fscanf(fp1s,"%le,%le,%le",&c1_B,&c2_L,&c3_FLUE);
      memset(bufs, 0, 1024); 
      fgets(bufs,1024,fp1s);
      if(c3_FLUE<0){
        c3_FLUE=0;
      }

      flue_min[FILE_NUM*hf + nn] = c3_FLUE;

      nns ++;

      if(nns >= FILE_NUM){
        break;
      } 
    }
    fclose(fp1s);
  }
}

double getPositionFlueMax(float dist, float lat, float lon, float* flue){
  int i,j,nL,nH;
  float dL,dH;

  if(dist>=360 && dist <=380){
    nL = 0;
    nH = 1;
    dL = 360.0;
    dH = 380.0;
  }else if(dist>380 && dist <= 400){
    nL = 1;
    nH = 2;
    dL = 380.0;
    dH = 400.0;
  }else if(dist>400 && dist <=420){
    nL = 2;
    nH = 3;
    dL = 400.0;
    dH = 420.0;
  }

  const size_t n_lon = 121;
  const size_t n_lat = 90;
  const gsl_interp2d_type *T_flue =gsl_interp2d_bilinear;
  double long_grid2[n_lon], lat_grid2[n_lat];

  for(i=0;i< n_lon;i++){
    long_grid2[i] = 3*i-180;
  }
  for (i = 0; i < n_lat; i ++){
    lat_grid2[i] = 2*i-89;
  }

  printf("%d    %d\n",nL,nH );

  double *za_flue_L = (double*)malloc(n_lat * n_lon * sizeof(double));
  double *za_flue_H = (double*)malloc(n_lat * n_lon * sizeof(double));

  gsl_spline2d *spline_flue_L = gsl_spline2d_alloc(T_flue,n_lat, n_lon);
  gsl_spline2d *spline_flue_H = gsl_spline2d_alloc(T_flue,n_lat, n_lon);

  gsl_interp_accel *lat_grid2acc = gsl_interp_accel_alloc();
  gsl_interp_accel *long_grid2acc = gsl_interp_accel_alloc();

  for(i = 0; i < n_lat; i ++){
    for(j = 0; j < n_lon; j ++){
      gsl_spline2d_set(spline_flue_L,za_flue_L,i,j,flue[FILE_NUM*nL + i*n_lon + j]);
      // printf("%f \n" ,flue_max[FILE_NUM*nL + i*n_lon + j]);
      gsl_spline2d_set(spline_flue_H,za_flue_H,i,j,flue[FILE_NUM*nH + i*n_lon + j]);
      // printf("%f \n" ,flue_max[FILE_NUM*nH + i*n_lon + j]);
    }
  }


  gsl_spline2d_init(spline_flue_L,  lat_grid2 , long_grid2, za_flue_L,n_lat, n_lon);
  gsl_spline2d_init(spline_flue_H,  lat_grid2 , long_grid2, za_flue_H,n_lat, n_lon);

  double val_L = gsl_spline2d_eval(spline_flue_L, lat, lon,lat_grid2acc,long_grid2acc);
  double val_H = gsl_spline2d_eval(spline_flue_H, lat, lon,lat_grid2acc,long_grid2acc);

  printf("%f    %f \n",val_L, val_H);

  double slop = (val_H - val_L)/(dH - dL);

  gsl_spline2d_free(spline_flue);
  gsl_interp_accel_free(lat_grid2acc);
  gsl_interp_accel_free(long_grid2acc);

  return slop*(dist-dL)+val_L;
}

double* readOrbitData(char* fileName, int* static_num){
    FILE *orbitFile = fopen(fileName,"r");
    *static_num = 0;   
    if(orbitFile == NULL){
        printf("No such file:'%s'\n",fileName);
        exit(1);
    }

    int num = 0;
    char line[MAXLINE];
    
    while (fgets(line, MAXLINE, orbitFile) != NULL){
    //int cmp_res = strcmp(line[0],"#");
    //if(*static_num < 1){
      //printf("%d \n",cmp_res);
    //}
    char* ta = "#";
    if(line[0] == ta[0]){
      
      continue;
    }

    *static_num = *static_num + 1;
  }
  double* data = (double*) malloc(sizeof(double)*(*static_num)*elemNum);
    fclose(orbitFile);
    orbitFile = fopen(fileName,"r");
    if(orbitFile == NULL){
        printf("No such file:'%s'\n",fileName);
        exit(1);
    }
    while (fgets(line, MAXLINE, orbitFile) != NULL) {
    char* ta = "#";
    if(line[0] == ta[0]){
      
      continue;
    }
        char *p;
        p = strtok(line, " \n");
        int cNum = 0;
        while (p != NULL) {
            double cc = atof(p);
            *(data + num * elemNum + cNum) = cc;
            cNum++;
            p = strtok(NULL, " \n");
        }
        num ++;
    }
    fclose(orbitFile);
    return data;
}

void readAllorbitsFile(double** OrbitData, int* dataNum){
  
  int i;
  for(i = 0; i < Orbit_File_Num; i ++){
    char orbitfileName[200] = "";
    strcat(orbitfileName,"orbit20160925/");
    char idn[10];
    sprintf(idn,"%d",i+1);
    strcat(orbitfileName,idn);
    strcat(orbitfileName,".txt");
    //printf("%s\n",orbitfileName);
    OrbitData[i] = readOrbitData(orbitfileName, &dataNum[i]);
  }
  
  for(i = 0; i < Orbit_File_Num; i ++){
    timePoint[i][0] = *(OrbitData[i]);
    timePoint[i][1] = *(OrbitData[i] + (dataNum[i]-1) * elemNum);
  }
}

int main(void) {

  float* flue_max = (float*)malloc(sizeof(float)*FILE_NUM*4);
  float* flue_min = (float*)malloc(sizeof(float)*FILE_NUM*4);
  aquireFlueData(flue_max,flue_min);

  double flue = getPositionFlueMax(370, -35, -30, flue_max);

  printf("%f\n", flue);


return 0;
}



