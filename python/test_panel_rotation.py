#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  8 18:01:04 2018

@author: xyh
"""

import os,sys
import numpy as np
import matplotlib.pylab as plt

def GetOriginalPanelNorm(ra,dec):
    phi = ra*np.pi/180
    theta = dec*np.pi/180
    p = np.array([np.cos(theta)*np.cos(phi),
                  np.cos(theta)*np.sin(phi),
                  np.sin(theta)])
    
    tmp = None
    if np.abs(dec) < 89:# use z-axis
        z = np.array([0,0,1])
        tmp = np.cross(p,z)
    else:# use the unit vector lay in the ecliptic plane
        zz = np.array([np.cos(phi),np.sin(phi),0])
        tmp = np.cross(p,zz)
        
    # normalize tmp
    p_norm = tmp/np.linalg.norm(tmp)

    return p_norm

def GetPanelAxis(ra,dec):
    phi = ra*np.pi/180
    theta = dec*np.pi/180
    p = np.array([np.cos(theta)*np.cos(phi),
                  np.cos(theta)*np.sin(phi),
                  np.sin(theta)])
    
    p_norm = GetOriginalPanelNorm(ra,dec)
    
    p_axis = np.cross(p,p_norm)
    p_axis /= np.linalg.norm(p_axis)

    return p_axis
    

def GenRotationMatrix(axis,theta):
    """
    Generate a rotation matrix, which rotates a vector or a point around the
    axes 'axis', by angle theta.
    Inputs:
        axis: rotation axis
        theta: rotation angle (in degree)
    """
    if len(axis) == 3:
        n = np.array([axis[0],axis[1],axis[2]])
        n_norm = np.linalg.norm(n)
        n = n/n_norm
        nx,ny,nz = n[0],n[1],n[2]
        
        cos_theta = np.cos(theta*np.pi/180)
        sin_theta = np.sin(theta*np.pi/180)
        
        R11 = cos_theta + nx*nx*(1-cos_theta)
        R12 = nx*ny*(1-cos_theta) - nz*sin_theta
        R13 = nx*nz*(1-cos_theta) + ny*sin_theta
        
        R21 = ny*nx*(1-cos_theta) + nz*sin_theta
        R22 = cos_theta + ny*ny*(1-cos_theta)
        R23 = ny*nz*(1-cos_theta) - nx*sin_theta
        
        R31 = nz*nx*(1-cos_theta) - ny*sin_theta
        R32 = nz*ny*(1-cos_theta) + nx*sin_theta
        R33 = cos_theta + nz*nz*(1-cos_theta)
    else:
        print('==> axes must have 3 elements!')
        sys.exit(0)

    
    return np.array([[R11,R12,R13],[R21,R22,R23],[R31,R32,R33]])


def GetAlphaRelation(ra,dec):
    x_sun = np.array([1,0,0]) # assume Sun is located at (1,0,0)
    p_norm0 = GetOriginalPanelNorm(ra,dec)
	
    if np.dot(p_norm0,x_sun) < 0:
        p_norm0 *= -1

    p_axis = GetPanelAxis(ra,dec)

    cosval = []
    alpha = np.linspace(-25,25,51)
    for i in range(len(alpha)):
        R = GenRotationMatrix(p_axis,alpha[i])
        p_norm = np.matmul(R,p_norm0)
        # 此处输出法线矢量，用于检验
#        if abs(abs(alpha[i])-0)<1e-5:
#            print('p_norm = ')
#            print(p_norm)
        cosval.append(np.dot(p_norm,x_sun))
    
    return alpha, np.array(cosval)

if __name__ == '__main__':
	if len(sys.argv) < 3:
		print('usage: %s ra dec (ra,dec in degrees)'%(sys.argv[0]))
		sys.exit(0)
#
#	plt.figure(figsize=(11,9))
#	plt.subplot(2,1,1)
#	ra = 30
#	dec = -0
#	while ra < 120:
#		ang,cosval = GetAlphaRelation(ra,dec)
#		plt.plot(ang,cosval,'-',label='ra = '+str(ra)+', dec = '+str(dec))
#		ra += 5
#	
#	plt.hlines(y=np.cos(25*np.pi/180),xmin=-25,xmax=25,linestyle='dashed',label='cos(25deg)')
#	plt.legend()
#	
#	plt.subplot(2,1,2)
#	ra = 80
#	dec = -45
#	while dec < 45:
#		ang,cosval = GetAlphaRelation(ra,dec)
#		plt.plot(ang,cosval,'-',label='ra = '+str(ra)+', dec = '+str(dec))
#		dec += 10
#	    
#	plt.hlines(y=np.cos(25*np.pi/180),xmin=-25,xmax=25,linestyle='dashed',label='cos(25deg)')
#	plt.legend()
#	
#	plt.show()
#	    
	ra = np.float(sys.argv[1])
	dec= np.float(sys.argv[2])
	alpha, cosval = GetAlphaRelation(ra,dec)
	ang = np.arccos(cosval.max())*180/np.pi
	
	print('ra = %g, dec = %g'%(ra,dec))
	print('angle = %g'%(ang))
	
#	radec = np.loadtxt('rand_radec.txt',delimiter=',')
#	fp = open('xxx.txt','w')
#	for i in range(len(radec)):
#		alpha,cosval = GetAlphaRelation(radec[i,0],radec[i,1])
#		ang = np.arccos(cosval.max())*180/np.pi
#		fp.write('%15.10f'%ang)
#	fp.close()