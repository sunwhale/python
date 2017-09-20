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
from workbench_convection import workbench
from compare_exp_sim import compare_exp_sim
from plot_sim_all import plot_sim_all
from plot_sim_pv import plot_sim_pv
from plot_exp_vs_sim_pv import plot_exp_vs_sim_pv
from plot_exp_pv import plot_exp_pv


for name in ['0000']:
    exp_filename = ExperimentDirectory + name + '.csv'
    experiment = ExperimentData(exp_filename)
#    workbench(name,loading_cycles=experiment.total_axial_count,copy=True)
    workbench(name,loading_cycles=1,copy=True)