from pylab import *
import os,sys
import numpy as np


def getAngle132( x1,  y1,  z1,  x2,  y2,  z2,  x3,  y3,  z3):
    
    cosValue = 0;
    angle = 0;
    
    x11 = x1-x3;
    y11 = y1-y3;
    z11 = z1-z3;
    
    x22 = x2-x3;
    y22 = y2-y3;
    z22 = z2-z3;
    
    tt = sqrt((x11*x11 + y11*y11 + z11* z11) * (x22*x22 + y22*y22 + z22*z22));   

    cosValue = (x11*x22+y11*y22+z11*z22)/tt;

    if (cosValue > 1):
        cosValue = 1;
    
    if (cosValue < -1):
        cosValue = -1;
    angle = np.arccos(cosValue);
    return angle * 360 / (2 * pi);



PANEL_TRANSE_ANGLE = 25.0
sun = zeros(3);
sun[0]=1.0;
sun[1]=0.0;
sun[2]=0.0;


ra = np.float(sys.argv[1]);
dec = np.float(sys.argv[2]);

p = zeros(3)
p[0] = cos(dec * 2 * pi / 360) * cos(ra * 2 * pi / 360);
p[1] = cos(dec * 2 * pi / 360) * sin(ra * 2 * pi / 360);
p[2] = sin(dec * 2 * pi / 360);

p_z = zeros(3);
p_z[0] = 0.0;
p_z[1] = 0.0;
p_z[2] = 1.0;

p_n = np.cross(p,p_z);#法线

p_n = p_n/linalg.norm(p_n)
p_nn = np.cross(p,p_n);#法线的法线

p_nn = p_nn/linalg.norm(p_nn)





t_tran_angle = PANEL_TRANSE_ANGLE*np.pi/180; 
t_cos= cos(t_tran_angle);
t_sin = sin(t_tran_angle);

normal_m = p_n[0]*p_n[0] + p_n[1]*p_n[1] + p_n[2]*p_n[2];
p_m = p[0]*p[0] + p[1]*p[1] + p[2]*p[2];

p_n1 = zeros(3);
p_n2 = zeros(3);

p_n1[0] = t_cos*p_n[0]*p_m/normal_m + t_sin*p[0];
p_n1[1] = t_cos*p_n[1]*p_m/normal_m + t_sin*p[1];
p_n1[2] = t_cos*p_n[2]*p_m/normal_m + t_sin*p[2];

p_n1 = p_n1/linalg.norm(p_n1);

p_n2[0] = t_cos*p_n[0]*p_m/normal_m - t_sin*p[0];
p_n2[1] = t_cos*p_n[1]*p_m/normal_m - t_sin*p[1];
p_n2[2] = t_cos*p_n[2]*p_m/normal_m - t_sin*p[2];

p_n2 = p_n2/linalg.norm(p_n2)


Sun_sqr = sun[0]*sun[0] + sun[1]*sun[1] + sun[2]*sun[2];
			# // double thres_for_sun_angle = COS_SUN_PLANE_ANGLE;	// 帆板法线与太阳矢量夹角
			
h_val = sun[0]*p_nn[0] + sun[1]*p_nn[1] + sun[2]*p_nn[2];
proj_sun = zeros(3);

proj_sun[0] = sun[0]-h_val*p_nn[0];
proj_sun[1] = sun[1]-h_val*p_nn[1];
proj_sun[2] = sun[2]-h_val*p_nn[2];

proj_sun = proj_sun/linalg.norm(proj_sun)

proj_sun_m = proj_sun[0]*proj_sun[0] + proj_sun[1]*proj_sun[1] + proj_sun[2]*proj_sun[2];
cos_proj_nor = (proj_sun[0]*p_n[0]+proj_sun[1]*p_n[1]+proj_sun[2]*p_n[2])/sqrt(normal_m*proj_sun_m);

cos_v = 0.0;
if(fabs(cos_proj_nor)>=t_cos):
    cos_v = (sun[0]*proj_sun[0] + sun[1]*proj_sun[1]+ sun[2]*proj_sun[2])/sqrt(Sun_sqr*proj_sun_m);
else:
    pn1_m = p_n1[0]*p_n1[0] + p_n1[1]*p_n1[1] + p_n1[2]*p_n1[2];
    pn2_m = p_n2[0]*p_n2[0] + p_n2[1]*p_n2[1] + p_n2[2]*p_n2[2];
    cos_v1 = (sun[0]*p_n1[0] + sun[1]*p_n1[1]+ sun[2]*p_n1[2])/sqrt(Sun_sqr*pn1_m);
    cos_v2 = (sun[0]*p_n2[0] + sun[1]*p_n2[1]+ sun[2]*p_n2[2])/sqrt(Sun_sqr*pn2_m);
    if cos_v1>cos_v2:
        cos_v = cos_v1;
    else:
        cos_v = cos_v2;
rangle = np.arccos(cos_v)*180.0/pi;
print(rangle)