# -*- coding: utf-8 -*-
import os,sys
import numpy as np

#	先定义一些常用的坐标转换函数

def RaDec_to_XYZ(ra=None,dec=None):
	if ra is not None and dec is not None:
		phi,theta = ra*np.pi/180,dec*np.pi/180
		x = np.array([np.cos(theta)*np.cos(phi), np.cos(theta)*np.sin(phi), np.sin(theta)])
		return x
	else:
		if ra is None:
			print('Error: ra is None, give it a meaningful value!')
			sys.exit(0)
		if dec is None:
			print('Error: dec is None, give it a meaningful value!')
			sys.exit(0)

def XYZ_to_RaDec(X=None):
	"""
	Convert Cartesian coordinates (X,Y,Z) to spherical coordinates (ra,dec)
	"""
	if len(X) == 3:
		x,y,z=X[0],X[1],X[2]
		r = (x**2+y**2+z**2)**0.5
		theta = np.arcsin(z/r)
		ra = np.arctan(y/(r*np.cos(theta)+x+1e-15))*360/np.pi # add 1e-15 to avoid dividing by ZERO
		dec = theta*180/np.pi
		if ra < 0: ra += 360
		return ra,dec
	else:
		print('Error: len(x) should be 3.')
		print('-- you input is:')
		print(X)
		sys.exit(0)


def EquatorialToEcliptic(lon_eq,lat_eq):
	phi,theta = lon_eq*np.pi/180,lat_eq*np.pi/180
	x_old = np.array([np.cos(theta)*np.cos(phi),np.cos(theta)*np.sin(phi),np.sin(theta)])
	# 绕x轴逆时针旋转23.4522度
	cosa = 0.917392
	sina = 0.397983
	R = np.array([  [1.0, 0.0, 0.0],
					[0.0, cosa, sina ],
					[0.0,-sina, cosa]])
	x_new = np.matmul(R,x_old)
	return XYZ_to_RaDec(x_new)

def EclipticToEquatorial(lon_ep,lat_ep):
	phi,theta = lon_ep*np.pi/180,lat_ep*np.pi/180
	x_old = np.array([np.cos(theta)*np.cos(phi),np.cos(theta)*np.sin(phi),np.sin(theta)])
	cosa = 0.917392
	sina = 0.397983
	R = np.array([  [1.0, 0.0, 0.0],
					[0.0, cosa,-sina ],
					[0.0, sina, cosa]])
	x_new = np.matmul(R,x_old)
	return XYZ_to_RaDec(x_new)

def EquatorialToGalactic(lon_eq,lat_eq):
	phi,theta = lon_eq*np.pi/180,lat_eq*np.pi/180
	x_old = np.array([np.cos(theta)*np.cos(phi),np.cos(theta)*np.sin(phi),np.sin(theta)])
	R = np.array([  [-0.05488, -0.87344, -0.48384],
					[ 0.49411, -0.44483,  0.74698],
					[-0.86767, -0.19808,  0.45598]])

	x_new = np.matmul(R,x_old)
	return XYZ_to_RaDec(x_new)


def GalacticToEquatorial(l,b):
	phi,theta = l*np.pi/180,b*np.pi/180
	x_old = np.array([np.cos(theta)*np.cos(phi),np.cos(theta)*np.sin(phi),np.sin(theta)])
	R = np.array([  [-0.05488, -0.87344, -0.48384],
					[ 0.49411, -0.44483,  0.74698],
					[-0.86767, -0.19808,  0.45598]])
	invR = np.linalg.inv(R)
	x_new = np.matmul(invR,x_old)
	ra,dec= XYZ_to_RaDec(x_new)
	return ra,dec

def GalacticToEcliptic(l,b):
	lon_eq,lat_eq = GalacticToEquatorial(l,b)
	return EquatorialToEcliptic(lon_eq,lat_eq)

def EclipticToGalactic(lon_ep,lat_ep):
	lon_eq,lat_eq = EclipticToEquatorial(lon_ep,lat_ep)
	return EquatorialToGalactic(lon_eq,lat_eq)


######################################################################
def GetRotationAngle(ra_old,dec_old,ra_new,dec_new):
	"""
	Get rotation angle for two given directions.
	"""
	phi_old, theta_old = ra_old*np.pi/180, (90-dec_old)*np.pi/180
	phi_new, theta_new = ra_new*np.pi/180, (90-dec_new)*np.pi/180
	x_old = np.array([  np.sin(theta_old)*np.cos(phi_old), 
						np.sin(theta_old)*np.sin(phi_old), 
						np.cos(theta_old)])
	x_new = np.array([  np.sin(theta_new)*np.cos(phi_new), 
						np.sin(theta_new)*np.sin(phi_new), 
						np.cos(theta_new)])
	cosval= np.dot(x_old,x_new)
	if cosval > 1: cosval = 1.0
	if cosval <-1: cosval = -1.0
	return np.arccos(cosval)*180/np.pi

def GetRotationAngleFromMatrix(M):
	"""
	Get rotation angle from given rotation matrix.
	Reture value has unit of degree
	"""
	if type(M) is not np.ndarray:
		print('Error: M must be a 3x3 np.ndarray object')
		sys.exit(0)
	trM = np.trace(M)
	cosval = 0.5*(trM-1)
	if cosval > 1: cosval = 1
	if cosval <-1: cosval =-1
	return np.arccos(cosval)*180/np.pi


def GenRotationMatrix(P,theta=0):
	"""
	Generate a rotation matrix that rotate a vector or a point around a given axis P,
	by angle theta. The axis P is defined by the unit vector n=(nx,ny,nz).
	Note that theta MUST be in unit of degree!
	"""
	if len(P) == 3:
		n = np.array([P[0],P[1],P[2]])
		n /= np.linalg.norm(n)
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

		return np.array([[R11,R12,R13],[R21,R22,R23],[R31,R32,R33]])
	else:
		print('==> P must have 3 elements!')
		sys.exit(0)



def RotateAroundP(X,P,theta):
	"""
	Rotate a vector X around the given direction P=(px,py,pz);
	The vector X has components (x,y,z) on the unit sphere
	"""
	rotation_matrix = GenRotationMatrix(P,theta)
	if type(X) is np.ndarray:
		tmp = []
		if len(X.shape) == 1 and len(X)==3:
			Xnew = np.matmul(rotation_matrix,X)
			tmp.append(Xnew[0])
			tmp.append(Xnew[1])
			tmp.append(Xnew[2])
		else:
			if X.shape[1] == 3:
				for i in range(X.shape[0]):
					Xnew = np.matmul(rotation_matrix,X[i,:])
					tmp.append([Xnew[0],Xnew[1],Xnew[2]])
			else:
				print('Error: X.shape[1] must be 3!')
				sys.exit(0)
		return np.array(tmp)

	elif type(X) is list:
		if len(X) ==3:
			x = np.array([X[0],X[1],X[2]])
			return np.matmul(rotation_matrix,x)
		else:
			print('Error: list X must have 3 lements!')
			sys.exit(0)



#################################################################
# run test
if __name__ == '__main__':
	print('---> start running test')

	# # test rotation matrix generation
	# n = [0,0,1]
	# m = GenRotationMatrix(n,45)
	# print('rotation matrix is:')
	# print(m)
	# print('m is of type:')
	# print(type(m))
	# print('rotation angle is:')
	# print(GetRotationAngleFromMatrix(m))

	# x_old = np.array([1,0,0])
	# x_new = RotateAroundP(x_old,n,45)
	# print('x_old:')
	# print(x_old)
	# print('x_new:')
	# print(x_new)

	print('EclipticToGalactic(0,0) = ')
	print(EclipticToGalactic(0,0))

	print('GalacticToEcliptic(0,0) = ')
	print(GalacticToEcliptic(0,0))

	print('GalacticToEcliptic(180,0) = ')
	print(GalacticToEcliptic(180,0))
