# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 13:46:18 2017

@author: SunJingyu
"""

from Constants import *
from Data import *
from workbench import workbench

for name in experiment_type_dict['BIAXIAL']:
    
    
name = '3001'
exp_filename = ExperimentDirectory + name + '.csv'
experiment = ExperimentData(exp_filename)
workbench(name,loading_cycles=experiment.total_axial_count,copy=True)