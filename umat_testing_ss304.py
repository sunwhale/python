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


#for name in experiment_type_dict['TC']:
for name in ['3112']:
    exp_filename = ExperimentDirectory + name + '.csv'
    experiment = ExperimentData(exp_filename)
#    workbench(name,loading_cycles=experiment.total_axial_count,copy=True)
#    workbench(name,loading_cycles=5,copy=True)
#    
#    compare_exp_sim(name,5,'axial_strain','axial_stress')
#    plot_sim_all(name,begin_cycle=1,end_cycle=100,xitem='axial_strain',yitem='axial_stress')
    
    compare_exp_sim(name,4,'axial_stress','shear_stress')
#    plot_sim_all(name,begin_cycle=1,end_cycle=9,xitem='axial_stress',yitem='shear_stress')
    
#    plot_exp_vs_sim_pv(name,item='axial_stress')

#for name in experiment_type_dict['TC']:
#    exp_filename = ExperimentDirectory + name + '.csv'
#    experiment = ExperimentData(exp_filename)
#    print name,experiment.axial_count[-1],experiment.runing_time[-1],experiment.runing_time[-1]/experiment.axial_count[-1]
#    print name,experiment.total_axial_count

#plot_exp_pv(experiment_type_dict['TC'],item='axial_stress')

#==============================================================================
# OUTPUT
#==============================================================================
#=========================   119.5/   120.0========== 99.58%
# PACC  4.099855867681935E-003
# TEMP  0.000000000000000E+000
# TRIAXIALITY  0.000000000000000E+000
# YD  0.391979424862558  