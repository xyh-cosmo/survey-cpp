# -*- coding: utf-8 -*-

from pylab import *

class SolarPanel:
	def __init__(self,
				 solar_constant=1.361E3,
				 panel_area=75,
				 cover_coefficient=0.87,
				 efficiency_at_25deg=0.3,
				 lost_ratio1=0.08,
				 lost_ratio2=0.22,
				 efficiency_decrease=0.07,
				 n_order=2):
		# 帆板的面积（平方米）
		self._solar_constant		= solar_constant
		# 帆板的面积（平方米）
		self._panel_aera			= panel_area
		# 布片系数（无量纲？）
		self._cover_coefficient		= cover_coefficient
		# 地面实验室25摄氏度时的转化效率（无量纲）
		self._efficiency_at_25deg	= efficiency_at_25deg
		# 发电能力损失8%
		self._lost_ratio1			= lost_ratio1
		# 在轨工作温度差导致的功率损失
		self._lost_ratio2			= lost_ratio2
		# 由于空间辐照导致的性能衰减(10年后)
		self._efficiency_decrease	= efficiency_decrease
		# cos的幂次
		self._n_order 				= n_order

	def EnergyGainPerSec(self,theta,T_yr=0.0):
		return 	self._solar_constant * \
				self._panel_aera * \
				self._cover_coefficient * \
				self._efficiency_at_25deg * \
				(1.0 - self._lost_ratio1) * \
				(1.0 - self._lost_ratio2) * \
				(1.0 - 0.1*self._efficiency_decrease*T_yr)*cos(theta)**(self._n_order)



#######################################

SP = SolarPanel(n_order=4)

# print('Energy gained per sec: %g'%(SP.EnergyGainPerSec(0,0)))
# print('Energy gained per sec: %g'%(SP.EnergyGainPerSec(0,5)))
# print('Energy gained per sec: %g'%(SP.EnergyGainPerSec(0,10)))

# print('total E-cell: ', 3*90*100*3600*1e-3)
# print('total E-cost: ', 7.5*92*60)
# print('total E-gain: ', 57*60*75*1.361*0.3*(1-0.08)*(1-0.22)*0.87*0.93*cos(25/180*pi)**4)


Q_ini = 3*90*100*3600*1e-3 # kJ
Q_cost1 = 7.5*35*60  # kJ
# E_charge = 1.35*75*0.87*0.3*(1-0.08)*(1-0.22)*cos(25./180*pi)**3.85*(1-0.07) - 7.5  # kW
E_charge = 1.361*75*0.87*(0.3-0.07)*(1-0.08)*(1-0.22)*cos(20./180*pi)**3.5 - 7.5  # kW

Q_tmp = Q_ini
cnt = 0
while( Q_tmp > 0.8*Q_ini ):
	Q_tmp = Q_tmp - Q_cost1 + E_charge*57*60
	cnt += 1
	print(Q_tmp)

print('cnt = ',cnt)