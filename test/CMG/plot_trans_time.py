from pylab import *

def GetTransTime(angle):
	data = array([[1,80],[20,127],[45,196],[180,581]])
	if angle > 180.:
		print('trans_angle is larger than 180.')
		sys.exit(0)
	time = 0.
	if angle < 1.:
		time = 70.0
	elif angle == 1.:
		time = 80.0
	else:
		for i in range(3):
			if angle > data[i,0] and angle <= data[i+1,0]:
				time = data[i,1] \
					+ (angle-data[i,0])*(data[i+1,1]-data[i,1])/(data[i+1,0]-data[i,0])
				break

	return time


angle = linspace(0,180,50)
T = []

for i in range(len(angle)):
	T.append(GetTransTime(angle[i]))

T = array(T)

plot(angle,T,'-o',label='Translation_Time vs Angle')
legend()

show()
