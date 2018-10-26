import os,sys
from pylab import *

if len(sys.argv) < 2:
	print('usage: %s rotation_angle.txt\n'%sys.argv[0])
	sys.exit(0)

d = loadtxt(sys.argv[1])

angle_old = d[:,6]
angle_new1= d[:,4]
angle_new2= d[:,5]

figure(figsize=(6,7))

subplot(3,1,1)
plot(angle_new1-angle_old,'.',label='new1 - old')
# xlabel('1000 random trials')
ylabel(r'$\Delta\theta$',fontsize=13)
legend()

subplot(3,1,2)
plot(angle_new2-angle_old,'.',label='new2 - old')
ylabel(r'$\Delta\theta$',fontsize=13)
# xlabel('1000 random trials')
legend()

subplot(3,1,3)
plot(angle_new2-angle_new1,'.',label='new2 - new1')
ylabel(r'$\Delta\theta$',fontsize=13)
xlabel('1000 random trials')
legend()

subplots_adjust(top=0.975,bottom=0.075)

####
figure()
plot(angle_new2,angle_new2-angle_old,'.')

show()