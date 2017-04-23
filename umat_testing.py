# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 13:46:18 2017

@author: SunJingyu
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from Work import UMAT
from Data import ExperimentData,ExperimentLog
from Functions import obtain_kinematic_hardening_parameters,calculate_elastic_by_temperature_in718
from Constants import *
from workbench import workbench
from compare_exp_sim import compare_exp_sim

for name in ['7047']:
    workbench(name,loading_cycles=10,copy=False)
#    compare_exp_sim(name,2,'axial_stress','shear_stress')
    compare_exp_sim(name,5,'axial_strain','axial_stress')