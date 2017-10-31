# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 13:46:18 2017

@author: SunJingyu
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from Constants import *
from Data import *
from workbench import workbench
from compare_exp_sim import compare_exp_sim
from plot_sim_all import plot_sim_all
from plot_sim_pv import plot_sim_pv
from plot_exp_vs_sim_pv import plot_exp_vs_sim_pv
from plot_exp_pv import plot_exp_pv

experiment_log = ExperimentLog(ExperimentLogFile)

#for name in ['7047']: #TC-IP
for name in ['7005']: #TC-IP
#for name in ['7018','7025','7017']: #NPR-IP
    experiment_log = ExperimentLog(ExperimentLogFile)
    experiment_log.output(name)
    regular = r'.*'
    load_type = experiment_log.obtainItem(name,'load_type',regular)[0]
    regular = r'\d+\.?\d*'
    temperature_mode = experiment_log.obtainItem(name,'temperature_mode',regular)
    d_out = float(experiment_log.obtainItem(name,'d_out',regular)[0])
    gauge_length = float(experiment_log.obtainItem(name,'gauge_length',regular)[0])
    axial_strain = float(experiment_log.obtainItem(name,'axial_strain',regular)[0])
    angel_strain = float(experiment_log.obtainItem(name,'angel_strain',regular)[0])
    equivalent_strain = float(experiment_log.obtainItem(name,'equivalent_strain',regular)[0])
    period = float(experiment_log.obtainItem(name,'period',regular)[0])
    axial_temperature_phase = float(experiment_log.obtainItem(name,'axial_temperature_phase',regular)[0])
    life = float(experiment_log.obtainItem(name,'comments',regular)[0])
    
    workbench(name,loading_cycles=10,copy=True)

#    compare_exp_sim(name,15,'axial_strain','axial_stress')
#    plot_sim_all(name,begin_cycle=1,end_cycle=9,xitem='axial_strain',yitem='axial_stress')

#    compare_exp_sim(name,5,'axial_stress','shear_stress')
#    plot_sim_all(name,begin_cycle=1,end_cycle=9,xitem='axial_stress',yitem='shear_stress')
    
#    plot_sim_pv(name,item='axial_stress')
#    plot_exp_vs_sim_pv(name,item='axial_stress')

#==============================================================================
# OUTPUT
#==============================================================================
#=========================  1796.8/  1800.0========== 99.82%
# PACC  4.283127271598782E-002
# TEMP   462.454128668044     
# TRIAXIALITY  0.235359841543725     
# YD  5.960464388721221E-006