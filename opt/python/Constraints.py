#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import os,sys
import numpy as np
import matplotlib.pylab as plt
from mpl_toolkits.basemap import Basemap

from Transformation import RaDec_to_XYZ,XYZ_to_RaDec,RotateAroundP, \
EquatorialToEcliptic, EclipticToEquatorial, GalacticToEcliptic, \
EclipticToGalactic

class Constraints:
    def __init__(self,
                sun_angle_min = 50,
                moon_angle_min = 40,
                earth_angle_light = 70,
                earth_angle_dark = 30,
                panel_rotate_angle=25,
                panel_rotate=False):
        # 光轴指向与太阳之间的最小夹角
        self._sun_angle_min = sun_angle_min
        # 光轴指向与月球之间的最小夹角
        self._moon_angle_min = moon_angle_min
        # 光轴指向与地球亮边的最小夹角
        self._earth_angle_light = earth_angle_light
        # 光轴指向与地球暗边的最小夹角
        self._earth_angle_dark = earth_angle_dark
        # 帆板可以转动的（最大）角度
        self._panel_rotate_angle = panel_rotate_angle
        # 帆板是否可以转动
        self._panel_rotate = panel_rotate
        
    def makeCircle(self,ra=None,dec=None,angle=None,size=100):
        """
        生成一系列点，绕着给定的天区坐标形成一个圆：每一个点对应的单位向量与(ra,dec)对应的单位
        向量的点积=cos(angle); 这些点的数目为size
        """
        # print('ra = %.10g, dec = %.10g'%(ra,dec))
        phi, theta = ra*np.pi/180, dec*np.pi/180
        X_sun = np.array([np.cos(theta)*np.cos(phi), np.cos(theta)*np.sin(phi), np.sin(theta)])

        # 第一个点的位置为(ra,dec)所在子午圈内的纬度值最大的那个点，即（ra,dec+angle);随后将这个点逆时针
        # 绕着(ra,dec)对应的单位向量旋转一周，共size-1次。
        phi, theta = ra*np.pi/180, (dec+angle)*np.pi/180
        x0 = np.array([np.cos(theta)*np.cos(phi), np.cos(theta)*np.sin(phi), np.sin(theta)])
        ra,dec = [],[]
        cnt = 0
        while cnt < size:
            angle = 360./(size-1)*cnt
            x = RotateAroundP(x0,X_sun,-angle)
            ra_tmp,dec_tmp = XYZ_to_RaDec(x)
            ra.append(ra_tmp)
            dec.append(dec_tmp)
            cnt += 1

        return np.array(ra), np.array(dec)
        
    def getSunConstraints(self,ra_sun=None,dec_sun=None,angle=50,size=100):
        """
        给出太阳对可观测天区的限制
        """
        return self.makeCircle(ra_sun,dec_sun,angle=angle,size=size)
    
    
    def getMoonConstraints(self,ra_moon=None,dec_moon=None,angle=40,size=100):
        """
        给出月球对可观测天区的限制
        """
        return self.makeCircle(ra_moon,dec_moon,angle=angle,size=size)
    
    
    def getEarthConstraints(self):
        """
        给出在地球地平线与地气光的限制下，还有哪些天区可以进行观测
        """
        pass

    def getPanelConstraints(self,sun,force_no_rotate=False):
        """
        给出帆板转动对可观测天区的限制
        """

    def getLargeSurveyBoundary(self,Bmin=20,bmin=17,size=200):
        """
        给出大面积巡天的边界
        """
        # step 1: create the boundary defined by B>=Bmin
        boundary_Bmin_up = []
        boundary_Bmin_down = []
        ra = 0
        while ra <= 360:
            dec = Bmin
            l,b = EclipticToGalactic(ra,dec)
            if np.abs(b)>bmin:
                boundary_Bmin_up.append([ra,dec])
            l,b = EclipticToGalactic(ra,-dec)
            if np.abs(b)>bmin:
                boundary_Bmin_down.append([ra,-dec])
            ra += 360./(size-1)


        # step 2: create the boundary defined by B<=Bmax
        boundary_bmin_up = []
        boundary_bmin_down = []
        l = 0
        while l <= 360:
            b = bmin
            ra,dec = GalacticToEcliptic(l,b)
            if np.abs(dec) > Bmin:
                boundary_bmin_up.append([ra,dec])
            ra,dec = GalacticToEcliptic(l,-b)
            if np.abs(dec) > Bmin:
                boundary_bmin_down.append([ra,dec])
            l += 360./(size-1)

        return np.array(boundary_Bmin_up),np.array(boundary_Bmin_down),\
                np.array(boundary_bmin_up),np.array(boundary_bmin_down)
        
    def getDeepSurveyBoundary(self,deep_sites,size=100):
        """
        给出极深场巡天区域的边界。目前假设每一个极深场区域都是圆形，且面积相同。
        """
        if type(deep_sites) is np.ndarray:
            if len(deep_sites.shape) != 2:
                print('Error: deep_sites must be a 2D array!')
                sys.exit(0)
            circles = []
            for i in range(deep_sites.shape[0]):
                # print('## debug: deep_sites[%d] = '%i)
                # print(deep_sites[i])
                ra,dec,ang = deep_sites[i][0],deep_sites[i][1],deep_sites[i][2]
                ra_tmp,dec_tmp = self.makeCircle(ra,dec,angle=ang*0.5,size=size)
                circles.append([ra_tmp,dec_tmp])
            return circles
        else:
            print('Error: deep_sites must be np.ndarray of dim nx3!')
            sys.exit(0)

        
#################################################################
# run test
if __name__ == '__main__':

    C = Constraints()

    fig = plt.figure(figsize=(8,6))
    
    m =Basemap(projection='hammer',lat_0=0,lon_0=180,resolution='l')
    m.drawparallels(np.arange(-90.,90.,30.),labels=[1,0,1,1],fontsize=14)
    m.drawmeridians(np.arange(-180.,180.,60.))
    
#    m = fig.add_subplot(111)

    # plot Sun at (ra=180,dec=0) and add a circle surrounding it
#    ra_sun,dec_sun = 90,0
#    m.scatter(ra_sun,dec_sun,50,marker='o',color='y',latlon=True,alpha=1)

    # generate the date used to plot the Sun
#    ra,dec = C.getSunConstraints(ra_sun=ra_sun,dec_sun=dec_sun,angle=50,size=360)
#    m.scatter(ra,dec,1,marker='.',color='y',latlon=True,alpha=1)

    # plot Sun at (ra=180,dec=0) and add a circle surrounding it
#    ra_moon,dec_moon = 120,0
#    m.scatter(ra_moon,dec_moon,25,marker='o',color='g',latlon=True,alpha=1)
    
    # generate the date used to plot the Moon
#    ra,dec = C.getMoonConstraints(ra_moon=ra_moon,dec_moon=dec_moon,angle=40,size=360)
#    m.scatter(ra,dec,1,marker='.',color='g',latlon=True,alpha=1)

    # 画出大面积巡天的边界
    B1,B2,b1,b2 = C.getLargeSurveyBoundary(Bmin=20,bmin=17,size=500)
    m.scatter(B1[:,0],B1[:,1],1,marker='.',color='b',latlon=True,alpha=1)
    m.scatter(B2[:,0],B2[:,1],1,marker='.',color='b',latlon=True,alpha=1)
    m.scatter(b1[:,0],b1[:,1],1,marker='.',color='b',latlon=True,alpha=1)
    m.scatter(b2[:,0],b2[:,1],1,marker='.',color='b',latlon=True,alpha=1)
    
#    m.scatter(B1[:,0],B1[:,1],1,marker='.',color='b',alpha=1)
#    m.scatter(B2[:,0],B2[:,1],1,marker='.',color='b',alpha=1)
#    m.scatter(b1[:,0],b1[:,1],1,marker='.',color='b',alpha=1)
#    m.scatter(b2[:,0],b2[:,1],1,marker='.',color='b',alpha=1)
    
    # 画出极深场的边界
    deep_sites = np.loadtxt('deep_survey_bak3.txt')
    deep_site_circles = C.getDeepSurveyBoundary(deep_sites,size=50)
    for i in range(len(deep_site_circles)):
        ra,dec = deep_site_circles[i][0],deep_site_circles[i][1]
        m.scatter(ra,dec,1,marker='.',color='r',latlon=True)
#        m.scatter(ra,dec,1,marker='.',color='r')
    
    x = np.linspace(150,210,200)
    y = np.linspace(-30,30,200)
    X,Y = np.meshgrid(x,y)
    Z = np.ones(X.shape)
    m.contourf(X,Y,Z,cmap='BuGn',latlon=True,alpha=0.5)
    
    plt.show()
