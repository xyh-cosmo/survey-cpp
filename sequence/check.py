from pylab import *

# def transTime(angle):
# 	if angle < 

raw_data = loadtxt('temp.txt')

T = raw_data[:,0]
dT = T[1:]-T[0:(len(T)-1)]
dt = dT*24.*3600.
t_exp = raw_data[:,16]
rot_angle = raw_data[:,17]

plot(dT,'.')


show()