#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  8 18:01:04 2018

@author: xyh
"""

import os, sys

sys.path.append('python')

from SurveyDebugger import Debugger

# result_file = '../tmp/E20_b17_beta_15_deep8.5yr_highL8.5yr_cont8.5yr_panel_rotates_Y.dat'
# t_interval  = '../tmp/E20_b17_beta_15_deep8.5yr_highL8.5yr_cont8.5yr_panel_rotates_Y.dT'

result_file = '../results_debug/0808/E20_b17_beta_15_deep8.5yr_highL8.5yr_cont8.5yr_panel_rotates_Y.dat'
t_interval  = '../results_debug/0808/E20_b17_beta_15_deep8.5yr_highL8.5yr_cont8.5yr_panel_rotates_Y.dT'

db = Debugger(	result_file,
				t_interval,
				num_of_check=50,
				back_search_steps=50,
				start_time=0,
				end_time=8.5,
				outdir='0808')

db.back_search('xxx')
db.make_plots()
