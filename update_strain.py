# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 13:46:18 2017

@author: SunJingyu
"""

import numpy as np
from Constants import *
from Data import *

#for name in experiment_type_dict['TC-IP-TGMF']+experiment_type_dict['TC-OP-TGMF']:
for name in ['14']:
    experiment_log = ExperimentLog(ExperimentLogFile)
    experiment_log.output(name)
    regular = r'.*'
    load_type = experiment_log.obtainItem(name,'load_type',regular)[0]
    regular = r'\d+\.?\d*'
    temperature_mode = experiment_log.obtainItem(name,'temperature_mode',regular)
    if len(temperature_mode) == 1:
        temperature_list = [float(temperature_mode[0]), float(temperature_mode[0])]
    if len(temperature_mode) == 2:
        temperature_list = [float(temperature_mode[0]), float(temperature_mode[1])]
    d_out = float(experiment_log.obtainItem(name,'d_out',regular)[0])
    gauge_length = float(experiment_log.obtainItem(name,'gauge_length',regular)[0])
    axial_strain = float(experiment_log.obtainItem(name,'axial_strain',regular)[0])
    angel_strain = float(experiment_log.obtainItem(name,'angel_strain',regular)[0])
    equivalent_strain = float(experiment_log.obtainItem(name,'equivalent_strain',regular)[0])
    period = float(experiment_log.obtainItem(name,'period',regular)[0])
    axial_temperature_phase = float(experiment_log.obtainItem(name,'axial_temperature_phase',regular)[0])
    life = float(experiment_log.obtainItem(name,'comments',regular)[0])
    
    exp_filename = ExperimentDirectory + name + '.csv'
    experiment = ExperimentData(exp_filename)
    experiment.updateStrain(period,[0.01*axial_strain*1.0,-0.01*axial_strain*1.0])