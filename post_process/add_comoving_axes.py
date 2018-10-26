# -*- coding:utf-8 -*-
#
# 目的：给巡天编排模拟的结果增加一个固定在望远镜上的局部的直角坐标系。
# 该局部坐标系由其三个轴的单位矢量来确定。这三个轴分别是：
#	1）望远镜光轴的指向
#	2）朝着太阳的帆板面的法线
#	3）由以上两个矢量叉乘得到的第三个单位矢量，也就是帆板的旋转轴
#	
#	这三个矢量分别以3维笛卡尔坐标的形式给出（黄道坐标系），并按照顺序
#	添加到原来模拟结果的每一行的最后，共计9列。

import os,sys
from pylab import *
import numpy as np

if len(sys.argv) < 3:
	print('usage: %s input_root out_root'%(sys.argv[0]))
	sys.exit(0)


################################################################
#	将赤道坐标转换成黄道坐标
def transform_equatorial_to_ecliptic(x_eq,y_eq,z_eq):
	cosa = 0.917392
	sina = 0.397983
	x_ec = x_eq
	y_ec = y_eq*cosa + z_eq*sina
	z_ec =-y_eq*sina + z_eq*cosa
	return x_ec,y_ec,z_ec

################################################################
def add_axes(sunx,suny,sunz,ra,dec):
	tmp = None
	p_n = None # 朝着太阳的帆板面的法线
	p_r = None # 帆板的旋转轴
	
	axes_t = '' # 望远镜的光轴
	axes_n = '' # 法线
	axes_r = '' # 帆板的旋转轴

	# get the unit vector pointing to the sky patch with corrdinate (ra,dec)
	phi = ra*np.pi/180
	theta = dec*np.pi/180
	p_t = np.array([np.cos(theta)*np.cos(phi),
					np.cos(theta)*np.sin(phi),
					np.sin(theta) ])

	# get the unit vector pointing to the Sun
	sunx1, suny1, sunz1 = transform_equatorial_to_ecliptic(sunx,suny,sunz)
	x_sun = np.array([sunx1,suny1,sunz1])
	x_sun /= np.linalg.norm(x_sun)

	# 根据天区指向与z轴（或者是处在黄赤道面上的一个指向）的叉乘来得到帆板的法线方向（只是候选，后面需要根据与太阳的关系来确定最终的结果）
	if abs(dec) < 89.999:
		z = np.array([0,0,1])
		tmp = np.cross(p_t,z)
	else:
		z = np.array([np.cos(phi), np.sin(phi), 0])
		tmp = np.cross(p_t,z)

	# 依据与太阳方位矢量的余弦值的符号来判断法线的方向。
	tmp /= np.linalg.norm(tmp)  #将矢量归一化
	if np.dot(x_sun,tmp) >= 0:
		p_n = tmp
	else:
		p_n = -1.0*tmp

	# 跟据p_t,p_n来得到p_r
	p_r = np.cross(p_t,p_n) # 此处无需再进行归一化

	# 将数值结果转换为字符串
	axes_t += '%8.5f '%(p_t[0])
	axes_t += '%8.5f '%(p_t[1])
	axes_t += '%8.5f '%(p_t[2])

	axes_n += '%8.5f '%(p_n[0])
	axes_n += '%8.5f '%(p_n[1])
	axes_n += '%8.5f '%(p_n[2])
	
	axes_r += '%8.5f '%(p_r[0])
	axes_r += '%8.5f '%(p_r[1])
	axes_r += '%8.5f '%(p_r[2])

	return ' ' + axes_t + axes_n + axes_r + '\n'

################################################################

ifile = sys.argv[1]+'.dat'
ofile = sys.argv[2]+'_with_local_axes.dat'

sr = loadtxt(ifile)

fp0= open(ifile,'r')
fp = open(ofile,'w')

for i in range(sr.shape[0]):
	content = fp0.readline()
	content = content[0:(len(content)-1)] # remove the newline char
	sunx,suny,sunz = sr[i,6],sr[i,7],sr[i,8]
	ra,dec = sr[i,2],sr[i,1]
	content += add_axes(sunx,suny,sunz,ra,dec)
	fp.write(content)


fp0.close()
fp.close()
