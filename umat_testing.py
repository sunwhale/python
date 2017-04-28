# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 13:46:18 2017

@author: SunJingyu
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from Constants import *
from workbench import workbench
from compare_exp_sim import compare_exp_sim

for name in ['7047']: #TC-IP
#for name in ['7046']: #NPR-IP
    workbench(name,loading_cycles=10,copy=False)
#    compare_exp_sim(name,2,'axial_stress','shear_stress')
#    compare_exp_sim(name,5,'axial_strain','axial_stress')


#==============================================================================
# OUTPUT
#==============================================================================
#=========================  1796.8/  1800.0========== 99.82%
# PACC  4.283127271598782E-002
# TEMP   462.454128668044     
# TRIAXIALITY  0.235359841543725     
# YD  5.960464388721221E-006