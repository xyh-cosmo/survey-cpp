#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  5 18:01:04 2018

@author: xyh
"""

# class for holding the survey simulation results and time-interval data

import os,sys
import numpy as np
import matplotlib.pylab as plt
from mpl_toolkits.basemap import Basemap

from Transformation import GetRotationAngle

class Debugger:
	def __init__(self,
				simulation_result=None,
				t_interval=None,
				num_of_check=10,
				back_search_steps=20,
				start_time=0.0,
				end_time=10.3,
				outdir=None):
	#	save file names
		self._simulation_result = simulation_result  # simulation result file
		self._t_interval = t_interval	 # time interval file
		self._back_search_steps=back_search_steps # number of steps to be searched in revrse order
		
	#	output dir
		self._outdir = outdir
		if self._outdir is None:
			print('==> outdir must be provided!')
			sys.exit(0)
		else:
			if os.path.isdir(self._outdir):
				print('--> fodler %s already exists ...'%(self._outdir))
			else:
				tmp=os.path.join('.',self._outdir)
				os.system('mkdir '+tmp)
			

	#	check input 'start_time' and 'end_time'
		if end_time <= start_time:
			print('==> end_time MUST be larger than start_time!')
			sys.exit(0)
		self._start_time = start_time	# start time after which the time intervals will be used
		self._end_time = end_time		# end time before which the time intervals will be used

	#	the number of checks to be made
		self._num_of_check = num_of_check

	#	simulation result part:
		self._obs_time = None
		self._dec = None
		self._ra = None
		self._satx = None
		self._saty = None
		self._satz = None
		self._sunx = None
		self._suny = None
		self._sunz = None
		self._moonx = None
		self._moony = None
		self._moonz = None
		self._area_norm = None
		self._area_deep = None
		self._is_deep = None
		self._in_disk = None
		self._exp_time = None
		self._tAngle = None
		self._insunside = None
		self._cmg = None
		self._battery = None
		self._solar_panel_angle = None
		self._saa_time = None
		self._sky_id = None
	
	#	time interval part
		self._dt_obs_time = None
		self._dtime = None
		self._dt_idx = None
		self._largest_dt_idx = None
		
	#	back search part
		self._bs_outfilename_root = None

		if self._simulation_result is None:
			print('--> Error: simulation_result is None, exit ...')
			sys.exit(0)
		else:
			if os.path.isfile(self._simulation_result):
				print('--> loading simulation results from %s'%(self._simulation_result))
				self._load_results()
			else:
				print('==> cannot find %s'%(self._simulation_result))
		if self._t_interval is None:
			print('--> Error: t_interval is None, exit ...')
			sys.exit(0)
		else:
			if os.path.isfile(self._t_interval):
				print('--> loading time intervals from %s'%(self._t_interval))
				self._load_time_interval()
			else:
				print('==> cannot find %s'%(self._t_interval))
				sys.exit(0)

	def _load_results(self):
		result = np.loadtxt(self._simulation_result,comments='#')
		self._obs_time = result[:,0]
		self._dec = result[:,1]
		self._ra = result[:,2]
		self._satx = result[:,3]
		self._saty = result[:,4]
		self._satz = result[:,5]
		self._sunx = result[:,6]
		self._suny = result[:,7]
		self._sunz = result[:,8]
		self._moonx = result[:,9]
		self._moony = result[:,10]
		self._moonz = result[:,11]
		self._area_norm = None
		self._area_deep = None
		self._is_deep = result[:,14]
		self._in_disk = result[:,15]
		self._exp_time = result[:,16]
		self._tAngle = result[:,17]
		self._insunside = result[:,18]
		self._cmg = result[:,19]
		self._battery = result[:,20]
		self._solar_panel_angle = result[:,21]
		self._saa_time = result[:,22]
		tmp = []
		for i in range(len(result[:,23])):
			tmp.append(int(result[i,23]))
		self._sky_id = np.array(tmp)
		print('--> finished loading results, %d lines are loaded'%(len(self._sky_id)))

	def _load_time_interval(self):
		# added the functionality of add only time inverals that are within [start_time, end_time]
		dt = np.loadtxt(self._t_interval)
		tmp1 = []
		tmp2 = []
		tmp3 = []
		for i in range(len(dt[:,5])):
			# t = self.J2yr(dt[i,0])
			t = dt[i,0] # no need to make transformation from Julian date to time in year
			if t >= self._start_time and t <= self._end_time:
				# print('--> adding time interval ...')
				tmp1.append(dt[i,0])
				tmp2.append(dt[i,1])
				tmp3.append(int(dt[i,5]))
		self._dt_obs_time = np.array(tmp1)
		self._dtime = np.array(tmp2)
		self._dt_idx = np.array(tmp3)

		# print('==> len(self._dt_obs_time) = %d'%len(self._dt_obs_time))
		# print('==> len(self._dtime) = %d'%len(self._dtime))
		# print('==> len(self._dt_idx) = %d'%len(self._dt_idx))

		print('--> finished loading time intervals, %d lines are loaded'%(len(self._dt_idx)))
		print('--> start time: %g yr'%(self._start_time))
		print('--> end_time: %g yr'%(self._end_time))
		self._select_largest_interval()


	def J2yr(self,t):
    	# translate Julian time into survey year
		return (t-2459766)/365.25

	def _select_largest_interval(self):
		print('--> finding the largest %d time interval and the corresponding result line #'%(self._num_of_check))
		idx_sort = np.argsort(self._dtime)
		tmp = []
		for i in range(self._num_of_check):
			tmp.append(int(self._dt_idx[idx_sort[len(idx_sort)-1-i]]))
		self._largest_dt_idx = np.array(tmp)
		print('--> selected largest time intervals are:')
		for i in range(len(self._largest_dt_idx)):
			print('dt[%2d] = %g'%(i,self._dtime[idx_sort[len(idx_sort)-1-i]]))
    	

	def print_sky(self,idx=-1):
		if idx >= 0:
			print('----------------------------------------')
			print('[ %7d-th sky patch\'s info ]'%(self._sky_id[idx]))
			print('%15s = %15.8f [Jyr]'%('obs-time',self._obs_time[idx]))
			print('%15s = %15.8f [secs]'%('exp-time',self._exp_time[idx]))
			print('%15s = %g'%('ra',self._ra[idx]))
			print('%15s = %g'%('dec',self._dec[idx]))
			print('%15s = %g'%('cmg value',self._cmg[idx]))
			print('%15s = %g [min = 19440, max = 97200]'%('battery',self._battery[idx]))
			print('%15s = %g'%('tAngle',self._tAngle[idx]))
			print('%15s = %g'%('SAA time',self._saa_time[idx]))
			isdeep = 'False'
			if self._is_deep[idx] > 0: isdeep = 'True'
			print('%15s = %s'%('is_deep',isdeep))

	def _get_one_direction(self,idx,ref_time):
		content = ''
		content += "{:<15}".format(str(idx))
		content += "{:<15}".format(str(round((self._obs_time[idx]-2459766)/365.25,6)) + ' ')
		dtime = (self._obs_time[idx]-ref_time)*86400
		content += "{:<15}".format(str(round(dtime,3)) + ' ')
		content += "{:<10}".format(str(self._exp_time[idx]) + ' ')
		content += "{:<12}".format(str(self._ra[idx]) + ' ')
		content += "{:<12}".format(str(self._dec[idx]) + ' ')
		content += "{:<12}".format(str(self._battery[idx]) + ' ')
		content += "{:<10}".format(str(self._cmg[idx]) + ' ')
		content += "{:<12}".format(str(self._tAngle[idx]) + ' ')
		content += "{:<8}".format(str(self._is_deep[idx]) + ' ')
		content += "{:<10}".format(str(self._insunside[idx]) + ' ')
		content += "{:<15}".format(str(self._solar_panel_angle[idx]) + ' ')
		content += "{:<10}".format(str(round(self._saa_time[idx]*86400,3)) + ' ')
		content += "{:<15}".format(str(self._satx[idx]) + ' ')
		content += "{:<15}".format(str(self._saty[idx]) + ' ')
		content += "{:<15}".format(str(self._satz[idx]) + ' ')
		content += "{:<15}".format(str(self._sunx[idx]) + ' ')
		content += "{:<15}".format(str(self._suny[idx]) + ' ')
		content += "{:<15}".format(str(self._sunz[idx]) + ' ')
		content += "{:<15}".format(str(self._moonx[idx]) + ' ')
		content += "{:<15}".format(str(self._moony[idx]) + ' ')
		content += "{:<15}".format(str(self._moonz[idx]) + ' ')
		content += '\n'
		return content

	def back_search(self,outfile_root=None):
		print('--> back-searching ...')
		self._bs_outfilename_root = outfile_root
		outfilename_root = outfile_root

		if outfile_root is None:
			print('--> no root name for the back-search results are provided, use test for simplicity')
			outfilename_root = 'test'
			self._bs_outfilename_root = 'test'

		if self._largest_dt_idx is None:
			print('==> self._largest_dt_idx is None, check it!')
			sys.exit(0)

		for i in range(len(self._largest_dt_idx)):
			idx = self._largest_dt_idx[i]
			fname = os.path.join(self._outdir,outfilename_root+'_'+str(i+1)+'.txt')
			fp = open(fname,'w')
			
			header = ''
			header += "{:<15}".format('#result-idx ')
			header += "{:<15}".format('time[yr] ')
			header += "{:<15}".format('time[secs] ')
			header += "{:<10}".format('exp-time ')
			header += "{:<12}".format('ra ')
			header += "{:<12}".format('dec ')
			header += "{:<12}".format('battery ')
			header += "{:<10}".format('cmg ')
			header += "{:<12}".format('tAnge ')
			header += "{:<8}".format('isdeep ')
			header += "{:<10}".format('sunside ')
			header += "{:<15}".format('panel-angle ')
			header += "{:<10}".format('saa-time ')
			header += "{:<15}".format('sat-x ')
			header += "{:<15}".format('sat-y ')
			header += "{:<15}".format('sat-z ')
			header += "{:<15}".format('sun-x ')
			header += "{:<15}".format('sun-y ')
			header += "{:<15}".format('sun-z ')
			header += "{:<15}".format('moon-x ')
			header += "{:<15}".format('moon-y ')
			header += "{:<15}".format('moon-z')
			header += '\n'
			fp.write(header)
			
			ref_time = self._obs_time[idx]
							
			for j in range(self._back_search_steps):
				if idx-j >= 0:
					content = self._get_one_direction(idx-j,ref_time)
					ref_time = self._obs_time[idx-j]
					fp.write(content)
			fp.close()
			print('--> back-searched results are saved into: %s'%(fname))

	def transform_equatorial_to_ecliptic(self,x_eq,y_eq,z_eq):
		cosa = 0.917392
		sina = 0.397983
		x_ec = x_eq
		y_ec = y_eq*cosa + z_eq*sina
		z_ec =-y_eq*sina + z_eq*cosa
		return x_ec,y_ec,z_ec

	def get_ra_dec(self,x,y,z):
		r = (x**2+y**2+z**2)**0.5
		theta = np.arcsin(z/r)
		ra = np.arctan(y/(r*np.cos(theta)+x))*360/np.pi
		dec = theta*180/np.pi
		if dec < 0: dec += 360
		return ra,dec

	def make_cone(self,ra_ref,dec_ref,deg,size=10000):
		phi0, theta0 = np.pi/180*ra_ref, np.pi/180*dec_ref
		x_ref = np.array([np.cos(theta0)*np.cos(phi0), np.cos(theta0)*np.sin(phi0), np.sin(theta0)])
		ra = []
		dec= []
		cosval_min = np.cos(deg*np.pi/180)
		cnt=0
		while cnt < size:
			ra_i = np.random.rand()*360
			u = 2*np.random.rand()-1
			dec_i= (np.arccos(u)-np.pi/2)*180
#			dec_i= (np.random.rand()-0.5)*180
			phi,theta = np.pi/180*ra_i, np.pi/180*dec_i
			x = np.array([np.cos(theta)*np.cos(phi), np.cos(theta)*np.sin(phi), np.sin(theta)])
			if np.dot(x_ref,x) >= cosval_min:
				ra.append(ra_i)
				dec.append(dec_i)
				cnt += 1
				
		return np.array(ra), np.array(dec)


	def add_panel_constrants(self,sunx,suny,sunz,deg=25,size=10000):
		norm_sun = np.sqrt(sunx*sunx + suny*suny + sunz*sunz)
		x_sun = np.array([sunx/norm_sun,suny/norm_sun,sunz/norm_sun])
		ra = []
		dec= []
		cosval_min = np.cos(deg*np.pi/180)
		cnt=0
		while cnt < size:
			ra_i = np.random.rand()*360
			u = 2*np.random.rand()-1
			dec_i= (np.arccos(u)-np.pi/2)*180
			phi,theta = np.pi/180*ra_i, np.pi/180*dec_i
			x_north = np.array([0,0,1])
			x0 = np.array([np.cos(theta)*np.cos(phi), np.cos(theta)*np.sin(phi), np.sin(theta)])
			x1 = np.cross(x_north,x0)
			x1 = x1/np.linalg.norm(x1)	
			if np.abs(np.dot(x_sun,x1)) >= cosval_min:
				ra.append(ra_i)
				dec.append(dec_i)
				cnt += 1
				
		return np.array(ra), np.array(dec)
		

	def make_one_plot(self, searched_file_root,idx,out_dir):
		# read info from the searched file
		sf = np.loadtxt(os.path.join(out_dir,searched_file_root+'_'+str(idx)+'.txt'))
		idx_tmp = int(sf[:,0].max())
		
		fig_dir = os.path.join(self._outdir,'png_'+str(idx))
		if os.path.isdir(fig_dir):
			print('--> fig folder %s already exists...'%(fig_dir))
		else:
			cmd = 'mkdir '+fig_dir
			os.system(cmd)

		# step-2: plot the back-serached directions
		for i in range(sf.shape[0]):
			fig = plt.figure(figsize=(12,6))
			m =Basemap(projection='hammer',lat_0=0,lon_0=-180.5,resolution='l')
			m.drawparallels(np.arange(-90.,90.,30.),labels=[1,0,1,1],fontsize=18)
			m.drawmeridians(np.arange(-180.,180.,60.))
			
			# step-0: plot ra=0,90,180,270
#			m.scatter(180,0,10,marker='o',color='y',latlon=True,alpha=1)

			ax = plt.gca()
			
			# step-1: plot the covered parts
			ra = self._ra[0:idx_tmp]
			dec= self._dec[0:idx_tmp]
			m.scatter(ra,dec,1,marker='.',color='c',latlon=True,alpha=0.15)

			# step-2: plot positions of the sun
			sunx_eq,suny_eq,sunz_eq = sf[i,16],sf[i,17],sf[i,18]
			sunx_ec,suny_ec,sunz_ec = self.transform_equatorial_to_ecliptic(sunx_eq,suny_eq,sunz_eq)
			ra,dec = self.get_ra_dec(sunx_ec,suny_ec,sunz_ec)
			m.scatter(ra,dec,10,marker='o',color='r',latlon=True,alpha=1)
			ra50,dec50 = self.make_cone(ra,dec,50)
			m.scatter(ra50,dec50,2,marker='.',color='r',latlon=True,alpha=0.3)
			x_sun,y_sun=m(ra,dec)

			# step-3: plot positions of the moon
			moonx_eq,moony_eq,moonz_eq = sf[i,19],sf[i,20],sf[i,21]
			moonx_ec,moony_ec,moonz_ec = self.transform_equatorial_to_ecliptic(moonx_eq,moony_eq,moonz_eq)
			ra,dec = self.get_ra_dec(moonx_ec,moony_ec,moonz_ec)
			m.scatter(ra,dec,10,marker='o',color='magenta',latlon=True,alpha=1)
			ra40,dec40 = self.make_cone(ra,dec,40)
			m.scatter(ra40,dec40,2,marker='.',color='m',latlon=True,alpha=0.3)
			x_moon,y_moon=m(ra,dec)

			# step-4: plot regions constrained by panel
			ra,dec=self.add_panel_constrants(sunx_ec,suny_ec,sunz_ec,25)
			m.scatter(ra,dec,2,marker='.',color='g',latlon=True,alpha=0.5)

			# step-5: plot the target observational directions
			ra,dec=sf[i,4],sf[i,5]
			m.scatter(ra,dec,10,marker='D',color='k',latlon=True,alpha=1)
			x_sky,y_sky=m(ra,dec)

		 	# step-6: plot positions of the telescope
			satx_eq,saty_eq,satz_eq = sf[i,13],sf[i,14],sf[i,15]
			satx_ec,saty_ec,satz_ec = self.transform_equatorial_to_ecliptic(satx_eq,saty_eq,satz_eq)
			ra,dec = self.get_ra_dec(satx_ec,saty_ec,satz_ec)
			m.scatter(ra,dec,10,marker='s',color='k',latlon=True,alpha=1)
			ra20,dec20 = self.make_cone(ra,dec,60)
			m.scatter(ra20,dec20,2,marker='.',color='b',latlon=True,alpha=0.4)
			x_sat,y_sat=m(ra,dec)

			ax.text(x_sun,y_sun,'SUN',color='r',fontsize=10)
			ax.text(x_moon,y_moon,'MOON',color='magenta',fontsize=10)
			ax.text(x_sky,y_sky,'obj-sky',color='k',fontsize=10)
			ax.text(x_sat,y_sat,'SAT',color='k',fontsize=10)
			if sf[i,10] == 1: # in sunside
				ax.set_title('in Sun side')

			figname = os.path.join(self._outdir,'png_'+str(idx)+'/'+str(i)+'.png')
			fig.savefig(figname,dvi=1000)
			print('--> successfully saved figure: %s'%(figname))
			plt.close(fig)
		

	def make_plots(self):
		for idx in range(1,self._num_of_check+1):
			self.make_one_plot(self._bs_outfilename_root,idx,self._outdir)

    		
