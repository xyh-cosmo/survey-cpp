#coding=utf-8

#############################################################################################
# 作者：许优华
# 日期：2018年1月31日
#
# 一些必要的说明：
# 1）平台机动包括“转动指向到下一个目标天区”和“稳像”
#	平台的转动在数据开始读取后0.1秒启动，稳像在平台机动结束前16秒开始，并假定持续
#   时间为15.5秒。
# 2）CCD焦面的刷新安排在平台机动结束前2秒开始，仅刷新一次，持续时间为1秒。
#############################################################################################

#############################################################################################
# 相机事件的定义：（每一个事件用一个数字来表示）
# 注意：以下10个事件是按照发生时间的先后顺序所排列的！
#############################################################################################

# 打开快门：1 （打开快门的时刻定义为一个拍摄周期的开始）
# 快门全开：2 （快门完全打开的时刻定义为曝光的开始）
# 关闭快门：3 （快门开始关闭的时刻定义为曝光的结束）
# 快门全关：4  (快门完全关闭后立即进行图像数据的读取，假定读取持续的时间为40秒)
# 平台机动：5 （假定平台机动发生在数据读取开始后0.1秒）
# 结束读取：6
# 开始稳像：7 （假定平台机动结束前16秒开始稳像）
# 开始刷新：8 （假定刷新在稳像结束前两秒开始，并只刷新一次，所需时间为1秒）
# 结束刷新：9
# 结束稳像：10（假定稳像所需时间固定为15.5秒,稳像结束后剩余0.5秒到达本次拍摄周期的末尾）


#############################################################################################
# 组部件标号及状态说明：
# 注意：巡天观测时所有部件都加电工作！
# 主焦面前端电箱		A（0:不工作；1:工作，具体功耗暂时未定）
# 主焦面数传电箱		B（0:不工作；1:工作，具体功耗暂时未定）
# 主焦面电源箱			C（0:不工作；1:工作，具体功耗暂时未定）
# 主焦面控制箱			D（0:不工作；1:工作，具体功耗暂时未定）
# 短波红外电箱			E（0:不工作；1:工作，10W；2:工作，25W）
# 短波红外制冷			F（0:不工作；1:工作，210W；2:工作，222W）
# 主焦面制冷			G（0:不工作；1:工作，633W；2:工作，660W；3:工作，6663W；4:工作，690W；5:工作750W）
# 热控				H（0:不工作；1:工作，28W；2:工作，31W）
# 定标光源组件			I（0:不工作；1:前期，0W，2:中后期，0.6W）
# 快门组件			J（0:完全关闭状态，75.6W；1:正在打开，75.6W；2:完全打开，75.6W；3:正在关闭，75.6W）
# 像旋微调			K（0:不工作；1:工作，35W）
# 主电控箱			L（0:不工作；1:工作，34.4W)
#
# 目前仅“快门组件”的状态完全清楚
#############################################################################################

# -------------------------------------------------------------------------------------------
import os,sys
from pylab import *

sec_per_day = 24.0*3600.


# -------------------------------------------------------------------------------------------
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

# -------------------------------------------------------------------------------------------
# 从笛卡尔坐标系到（天文球坐标系的变换）
def Cart2Spherical_deg(x,y,z):
	r = (x**2+y**2+z**2)**0.5
	theta = arcsin(z/r)
	phi = 2*arctan(y/(r*cos(theta)+x))
	
	ra  = phi*180/pi
	dec = theta*180/pi

	if ra<0:
		ra += 360

	return ra,dec

# -------------------------------------------------------------------------------------------
# 获取旋转矩阵，同时保证不产生像旋
def get_rotation_axes(	theta_old, 
						phi_old, 
						theta_new, 
						phi_new ):
	
	sin_theta_old	= sin(theta_old)
	cos_theta_old	= cos(theta_old)
	sin_phi_old		= sin(phi_old)
	cos_phi_old		= cos(phi_old)
	
	sin_theta_new	= sin(theta_new)
	cos_theta_new	= cos(theta_new)
	sin_phi_new		= sin(phi_new)
	cos_phi_new		= cos(phi_new)
	
	R3_	= [ [ cos_phi_old, sin_phi_old, 0.0],
			[-sin_phi_old, cos_phi_old,  0.0],
			[ 0.0,         0.0,         1.0] ]
	
	R2_	= [ [ cos_theta_old, 0.0,-sin_theta_old],
			[ 0.0,           1.0, 0.0],
			[ sin_theta_old, 0.0, cos_theta_old] ] 
	
	R2 	= [ [ cos_theta_new, 0.0, sin_theta_new],
			[ 0.0,           1.0, 0.0],
			[-sin_theta_new, 0.0, cos_theta_new] ] 

	R3	= [ [ cos_phi_new,-sin_phi_new, 0.0],
			[ sin_phi_new, cos_phi_new,  0.0],
			[ 0.0,         0.0,         1.0] ]
			
	R3_ = array(R3_)
	R2_ = array(R2_)
	R2  = array(R2)
	R3  = array(R3)
	
	R = matmul(R3, matmul(R2, matmul(R2_,R3_)))
	
	axes_x = R[2,1]-R[1,2]
	axes_y = R[0,2]-R[2,0]
	axes_z = R[1,0]-R[0,1]
	
	return axes_x, axes_y, axes_z
	
def get_rotation_angle(	theta_old, 
						phi_old, 
						theta_new, 
						phi_new ):
	sin_theta_old	= sin(theta_old)
	cos_theta_old	= cos(theta_old)
	sin_phi_old		= sin(phi_old)
	cos_phi_old		= cos(phi_old)
	
	sin_theta_new	= sin(theta_new)
	cos_theta_new	= cos(theta_new)
	sin_phi_new		= sin(phi_new)
	cos_phi_new		= cos(phi_new)
	
	Trace_of_R =  cos_phi_old*cos_theta_old*cos_theta_new*cos_phi_new \
				+ sin_phi_old*sin_phi_new \
				+ cos_phi_old*sin_theta_old*sin_theta_new*cos_phi_new \
				+ sin_phi_old*cos_theta_old*cos_theta_new*sin_phi_new \
				+ cos_phi_new*cos_phi_old \
				+ sin_phi_old*sin_theta_old*sin_theta_new*sin_phi_new \
				+ sin_theta_old*sin_theta_new \
				+ cos_theta_old*cos_theta_new

	cosval = 0.5*(Trace_of_R -1.)
	if cosval < -1: cosval = -1
	if cosval > 1.: cosval = 1
	
	return arccos(cosval)

# -------------------------------------------------------------------------------------------
# 根据当前和下一次指向（以及转动时间和转动角度）来计算中间任意时刻的指向
# 假设转动是匀速的
def get_tmp_direction(	ra_current,
						dec_current,
						ra_next,
						dec_next,
						trans_angle,
						trans_time_total,
						trans_time_passed):
	"""
	(ra_current,dec_current) 当前指向
	(ra_next,dec_next) 		 下一个指向
	trans_angle, 			 总的转动角度
	trans_time_total, 		 总的转动时间 (seconds)
	trans_time_passed, 		 已经花费的转动时间 (seconds)
	"""
	theta_now 	= 0.5*pi-dec_current*pi/180
	phi_now   	= ra_current*pi/180
	theta_next 	= 0.5*pi-dec_next*pi/180
	phi_next   	= ra_next*pi/180

	# 当前指向的笛卡尔坐标
	a = [sin(theta_now)*cos(phi_now), sin(theta_now)*sin(phi_now), cos(theta_now)]
	a = array(a)
	a = a/norm(a)
#	b = [sin(theta_next)*cos(phi_next), sin(theta_next)*sin(phi_next), cos(theta_next)]
	
	#get rotation angle for debug
#	anglexx = get_rotation_angle(theta_now,phi_now,theta_next,phi_next)
#	if abs(anglexx*180/pi-trans_angle) > 1.5:
#		print('anglexx = %10.6f,  tAngle = %10.6f'%(anglexx*180/pi, trans_angle))
	
	# rotation axes
	# get rotation axes
	r_x, r_y, r_z = get_rotation_axes(theta_now, phi_now, theta_next, phi_next)
	u = array([r_x,r_y,r_z])
	u = u/norm(u)

	# rotation matrix
	angle = trans_angle*trans_time_passed/trans_time_total*pi/180
	cosA = cos(angle)
	sinA = sin(angle)
	R = [ [cosA+u[0]*u[0]*(1-cosA),      u[0]*u[1]*(1-cosA)-u[2]*sinA, u[0]*u[2]*(1-cosA)+u[1]*sinA],
		  [u[1]*u[0]*(1-cosA)+u[2]*sinA, cosA+u[1]*u[1]*(1-cosA),      u[1]*u[2]*(1-cosA)-u[0]*sinA],
		  [u[2]*u[0]*(1-cosA)-u[1]*sinA, u[2]*u[1]*(1-cosA)+u[0]*sinA, cosA+u[2]*u[2]*(1-cosA)] ]
		  
#	print('R:')
#	print R

	R = array(R)
	c = matmul(R,a)
	
#	print a
#	print c
	ra, dec = Cart2Spherical_deg(c[0],c[1],c[2])
	
	return ra, dec


# -------------------------------------------------------------------------------------------
def check_time_interval(time_current,
						time_next,
						shutter_time,
						exposure_time,
						trans_angle
						):
	transTime = GetTransTime(trans_angle)
	time_interval = (time_next - time_current)*sec_per_day
	time_diff = time_interval - (2.*shutter_time+exposure_time+transTime)
	if time_diff > 0:
		return True
	else:
		print('\n--- Error ---')
		print('time_exposure = %g'%(exposure_time))
		print('time_interval = %g'%(time_interval))
		print('angle         = %g'%(trans_angle))
		print('time_trans    = %g'%(transTime))
		print('time_needed   = %g'%(2.*shutter_time+exposure_time+transTime))
		print('time diff     = %g'%(time_diff))
		return False

# ----------------------------------------------------------------------------------------------------------

def print_state_head():
#	state_head = "{:<20}".format('time cost (seconds)')+'\t '
	state_head = "{:<20}".format('事件持续时间  (seconds)')+'\t '
	state_head += 'sA\t '
	state_head += 'sB\t '
	state_head += 'sC\t '
	state_head += 'sD\t '
	state_head += 'sE\t '
	state_head += 'sF\t '
	state_head += 'sG\t '
	state_head += 'sH\t '
	state_head += 'sI\t '
	state_head += 'sJ\t '
	state_head += 'sK\t '
	state_head += 'sL'

	return state_head

def print_state(time_cost=-1,
				sA=0,sB=0,sC=0,sD=0,
				sE=0,sF=0,sG=0,sH=0,
				sI=0,sJ=0,sK=0,sL=0):
	state = "{:<20}".format(str(time_cost))+'\t '
	state += str(sA)+'\t '
	state += str(sB)+'\t '
	state += str(sC)+'\t '
	state += str(sD)+'\t '
	state += str(sE)+'\t '
	state += str(sF)+'\t '
	state += str(sG)+'\t '
	state += str(sH)+'\t '
	state += str(sI)+'\t '
	state += str(sJ)+'\t '
	state += str(sK)+'\t '
	state += str(sL)

	return state

def produce_one_sequence(time_current,
						 time_next,
						 ra_current,
						 dec_current,
						 ra_next,
						 dec_next,
						 shutter_time,
						 exposure_time,
						 trans_angle):
	#
	# 现在增加了指向信息的输出：指向在整个序列中会有所变化，按照如下的顺利：
	# 1）指向当前目标天区
	# 2）曝光结束后机动开始，指向开始逐步移向下一个目标天区
	# 3）机动结束，指向下一个目标天区
	# 备注：目前假设整个转动过程是匀速的，不考虑加速和减速的过程

	sequence = ""
	
	# step 1) 快门开始打开
	event_time = time_current
	state = print_state(time_cost=1.5,
						sA=1,sB=0,sC=1,sD=1,
						sE=2,sF=2,sG=1,sH=1,
						sI=0,sJ=1,sK=1,sL=1)
	temp = "{:<20}".format(str(event_time)) + '\t1\t\t' + state + '\t\t'
	temp += "{:<20}".format(str(round(ra_current,4))) # ra
	temp += "{:<20}".format(str(round(dec_current,4))) # dec
	sequence += temp + '\n'


	# step 2) 快门完全打开，开始曝光
	event_time = time_current + shutter_time/sec_per_day
	state = print_state(time_cost=round(exposure_time,5),
						sA=1,sB=0,sC=1,sD=1,
						sE=2,sF=2,sG=3,sH=1,
						sI=0,sJ=2,sK=1,sL=1)
	temp = "{:<20}".format(str(event_time)) + '\t2\t\t' + state + '\t\t'
	temp += "{:<20}".format(str(round(ra_current,4))) # ra
	temp += "{:<20}".format(str(round(dec_current,4))) # dec
	sequence += temp + '\n'

	# step 3) 曝光结束，快门开始关闭
	event_time = time_current + (shutter_time+exposure_time)/sec_per_day
	state = print_state(time_cost=1.5,
						sA=1,sB=0,sC=1,sD=1,
						sE=2,sF=2,sG=3,sH=1,
						sI=0,sJ=3,sK=1,sL=1)
	temp = "{:<20}".format(str(event_time)) + '\t3\t\t' + state + '\t\t'
	temp += "{:<20}".format(str(round(ra_current,4))) # ra
	temp += "{:<20}".format(str(round(dec_current,4))) # dec
	sequence += temp + '\n'

	# step 4) 快门完全关闭，开始读取 （假定读取数据所需时间为40秒）
	event_time = time_current + (shutter_time*2+exposure_time)/sec_per_day
	state = print_state(time_cost=40,
						sA=1,sB=1,sC=1,sD=1,
						sE=2,sF=2,sG=4,sH=2,
						sI=0,sJ=0,sK=0,sL=1)
	temp = "{:<20}".format(str(event_time)) + '\t4\t\t' + state + '\t\t'
	temp += "{:<20}".format(str(round(ra_current,4))) # ra
	temp += "{:<20}".format(str(round(dec_current,4))) # dec
	sequence += temp + '\n'

	# step 5) 平台机动开始，指向下一个目标天区, 此时刻的指向依然为(ra_current, dec_current)
	event_time = time_current + (shutter_time*2+exposure_time+0.1)/sec_per_day
	transTime = GetTransTime(trans_angle) - 0.1 #减去0.1秒是为了补偿前面增加的0.1秒
	state = print_state(time_cost=round(transTime,5),
						sA=1,sB=1,sC=1,sD=1,
						sE=1,sF=2,sG=5,sH=2,
						sI=0,sJ=0,sK=0,sL=1)
	temp = "{:<20}".format(str(event_time)) + '\t5\t\t' + state + '\t\t'
	temp += "{:<20}".format(str(round(ra_current,4))) # ra
	temp += "{:<20}".format(str(round(dec_current,4))) # dec
	sequence += temp + '\n'

	# step 6) 平台机动结束前16秒开始稳像

	# 计算此时刻的指向
	ra_tmp,dec_tmp = get_tmp_direction(	ra_current,
										dec_current,
										ra_next,
										dec_next,
										trans_angle,
										transTime,
										transTime - 16 )

	event_time = time_current + (shutter_time*2+exposure_time + 0.1 + (transTime-16))/sec_per_day
	state = print_state(time_cost=15,
						sA=0,sB=0,sC=1,sD=1,
						sE=1,sF=1,sG=2,sH=2,
						sI=0,sJ=0,sK=0,sL=1)
	temp = "{:<20}".format(str(event_time)) + '\t6\t\t' + state + '\t\t'
	temp += "{:<20}".format(str(round(ra_tmp,4))) # ra
	temp += "{:<20}".format(str(round(dec_tmp,4))) # dec
	sequence += temp + '\n'

	# step 7) 平台机动结束前2秒开始刷新
	ra_tmp,dec_tmp = get_tmp_direction(	ra_current,
										dec_current,
										ra_next,
										dec_next,
										trans_angle,
										transTime,
										transTime - 2 )

	event_time = time_current + (shutter_time*2+exposure_time + 0.1 + (transTime-2))/sec_per_day
	state = print_state(time_cost=1,
						sA=0,sB=0,sC=1,sD=1,
						sE=1,sF=1,sG=3,sH=2,
						sI=0,sJ=0,sK=1,sL=1)
	temp = "{:<20}".format(str(event_time)) + '\t7\t\t' + state + '\t\t'
	temp += "{:<20}".format(str(round(ra_tmp,4))) # ra
	temp += "{:<20}".format(str(round(dec_tmp,4))) # dec
	sequence += temp + '\n'

	# step 8) 平台机动结束前1秒开始刷新结束
	ra_tmp,dec_tmp = get_tmp_direction(	ra_current,
										dec_current,
										ra_next,
										dec_next,
										trans_angle,
										transTime,
										transTime - 1 )

	event_time = time_current + (shutter_time*2+exposure_time + 0.1 + (transTime-1))/sec_per_day
	state = print_state(time_cost=1,
						sA=0,sB=0,sC=1,sD=1,
						sE=1,sF=1,sG=3,sH=2,
						sI=0,sJ=0,sK=1,sL=1)
	temp = "{:<20}".format(str(event_time)) + '\t8\t\t' + state + '\t\t'
	temp += "{:<20}".format(str(round(ra_tmp,4))) # ra
	temp += "{:<20}".format(str(round(dec_tmp,4))) # dec
	sequence += (temp + '\n')
	
	# debug
#	anglexx = get_rotation_angle(ra_current*pi/180,dec_current*pi/180,ra_tmp*pi/180,dec_tmp*pi/180)
#	if abs(anglexx*180/pi-trans_angle*(transTime-1)/transTime) > 1e-1 :
#		print('WARNING: anglexx = %g, angle00 = %g'%(anglexx*180/pi,trans_angle*(transTime-1)/transTime))

	# step 9) 平台机动结束、稳像结束，进入下一个拍摄周期.当前指向已经到达下一个目标指向(ra_next, dec_next)
	ra_tmp,dec_tmp = get_tmp_direction(	ra_current,
										dec_current,
										ra_next,
										dec_next,
										trans_angle,
										transTime,
										transTime )
	# for debug only, count possible errors
	state_is_ok = 1
	diff_ra = ra_tmp-ra_next
	diff_ra0 = (ra_tmp-ra_next-360)
	diff_dec= dec_tmp-dec_next
	if (abs(diff_ra) > 0.01 and abs(diff_ra0) > 0.01) or abs(diff_dec) > 0.01:
		state_is_ok = 0
		print('ra_old = %10.6f, dec_old = %10.6f'%(ra_current, dec_current))
		print('ra_new = %10.6f, dec_new = %10.6f'%(ra_next, dec_next))
		print('ra_tmp = %10.6f, dec_tmp = %10.6f'%(ra_tmp, dec_tmp))
		print('ra_diff= %10.6f, dec_diff= %10.6f'%(diff_ra, diff_dec))
		print('tAngle = %g\n'%(trans_angle))
		# ra_err.append(diff_ra)
		# dec_err.append(diff_dec)
		# sys.exit(0)

	event_time = time_current + (shutter_time*2+exposure_time + 0.1 + (transTime-0.5))/sec_per_day
	state = print_state(time_cost=0.5,
						sA=0,sB=0,sC=1,sD=1,
						sE=1,sF=1,sG=3,sH=2,
						sI=0,sJ=0,sK=1,sL=1)
	temp = "{:<20}".format(str(event_time)) + '\t9\t\t' + state + '\t\t'
	temp += "{:<20}".format(str(round(ra_next,4))) # ra
	temp += "{:<20}".format(str(round(dec_next,))) # dec
	sequence += (temp + '\n')

	sequence += '\n'

	return sequence, state_is_ok

def generate_sequence(time_start,
					  time_exposure,
					  trans_angle,
					  sky_ra,
					  sky_dec,
					  shutter_time=1.5,
					  seq_filename='action_seq.txt',
					  operation_time_in_yr=1.5):
	fp = open(seq_filename,'w')
#	print >> fp, '# event time (JD)       event   ' + print_state_head()
	# print >> fp, '# 事件触发时间 (JD)     事件   ' + print_state_head()
	fp.write('# 事件触发时间 (JD)     事件   ' + print_state_head()+ '    指向 (ra, dec) ' + '\n')
	
	time_head = time_start.min()
	seq_len = len(time_start)-1

	cnt = 0
	cnt_err = 0
	while cnt < seq_len:
#	while cnt < 10:
#		print('counting : %7d'%(cnt))
		time_current	= time_start[cnt]
		time_next	 	= time_start[cnt+1]
		exposure_time	= time_exposure[cnt] # 这里将快门时间包含在曝光时间内
		angle 			= trans_angle[cnt+1]
		ra_current 		= sky_ra[cnt]
		ra_next 		= sky_ra[cnt+1]
		dec_current		= sky_dec[cnt]
		dec_next		= sky_dec[cnt+1]
		shutter_time	= shutter_time

		# 修正快门时间，从150秒或250秒减去快门开闭的时间（3秒）
		exposure_time -= 3.0

		# if check_time_interval(time_current,time_next,shutter_time,exposure_time,angle) == True:

		seq, state = produce_one_sequence( time_current,
											time_next,
											ra_current,
											dec_current,
											ra_next,
											dec_next,
											shutter_time,
											exposure_time,
											angle)
		# print >> fp, '%s'%(seq)
		fp.write('%s'%(seq))
		if state != 1:
			cnt_err += 1

		# 判断整个运行序列是否满足所需要的时间长度
		if time_current > time_head+operation_time_in_yr*365.25:
			print('total time %g yr has reached, stop here.'%(operation_time_in_yr))
			break
		cnt += 1

	fp.close()
	print('cnt = %d, cnt_err = %d'%(cnt,cnt_err))


####################################################################################################
#raw_data = loadtxt('temp.txt')
#raw_data = loadtxt('new_sky_4_no_gal_center.dat')
#raw_data = loadtxt('0603-A.dat')
raw_data = loadtxt('0603-A_beta_of.dat')

time_start 		= raw_data[:,0]
time_exposure	= raw_data[:,16]
trans_angle		= raw_data[:,17]
sky_position_ra	= raw_data[:,2]
sky_position_dec= raw_data[:,1]

#新增的指向(ra,dec)
shutter_time	= 1.5

generate_sequence(	time_start,
					time_exposure,
					trans_angle,
					sky_position_ra,
					sky_position_dec,
					shutter_time,
					operation_time_in_yr=1.5 )

#ra_old, dec_old = 0, 0
#ra_new, dec_new = 45, 0

#ra_tmp, dec_tmp = get_tmp_direction(ra_old,dec_old,ra_new, dec_new,45,1,0.5)

#print('ra_tmp = %g, dec_tmp = %g'%(ra_tmp,dec_tmp))
