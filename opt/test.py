import os, sys

sys.path.append('python')

from SurveyDebugger import Debugger

result_file = '../0710/E20_b17_beta_15_with_high_altitude_prior.dat'
t_interval  = '../dT.txt'

db = Debugger(	result_file,
				t_interval,
				num_of_check=20,
				back_search_steps=40,
				start_time=2,
				end_time=7.5)

db.back_search('xxx')
db.make_plots()
