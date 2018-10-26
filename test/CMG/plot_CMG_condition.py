from pylab import *

def CMG_allowed_rot_times(angle):
	lens = len(angle)
	results = []

	def fun(x):
		if x < 5:
			return 0
		elif x < 10:
			return 1./29
		elif x < 20:
			return 1./19
		elif x < 35:
			return 1./13
		elif x < 45:
			return 1./10
		elif x < 75:
			return 1./6
		elif x < 90:
			return 1./5
		elif x < 135:
			return 1./3
		elif x <= 180:
			return 1./2

	for n in range(lens):
		results.append(fun(angle[n]))

	return array(results)

angle = linspace(0,180,50)

plot(angle, CMG_allowed_rot_times(angle),'-o')
xlabel('Angle (deg)')
ylabel('CMG cost')

show()
