# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 13:46:18 2017

@author: SunJingyu
"""

import numpy as np
from Constants import *
from Data import *
from Work import *
from Plot import *
from multicolored_lines import *
from plot_format import plot_format

name = '7037'

#==============================================================================
# experiment_log
#==============================================================================
experiment_log= ExperimentLog('F:\\Database\\IN718\\Timed\\Inconel718_test_log.csv')
experiment_log.output(name)
#number_list = experiment_log.keyFilter("(self.load_type == 'cyclic diamond path')")
#print number_list
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
period = float(experiment_log.obtainItem(name,'period',regular)[0])
#==============================================================================
# experiment
#==============================================================================
exp = ExperimentData('F:\\Database\\IN718\\Timed\\' + name + '.csv')
#==============================================================================
# load
#==============================================================================
temperature_mean = (temperature_list[0] + temperature_list[1])/2.0
temperature_min = min(temperature_list)
temperature_max = max(temperature_list)
load = Load(runing_time=[0], temprature=[temperature_mean], axial_strain=[0], shear_strain=[0], first_cycle_shift=1)
axial_strain = axial_strain/100.0
shear_strain = np.deg2rad(angel_strain)*d_out/2.0/gauge_length
#==============================================================================
# Diamond path TMF IP
#==============================================================================
#load.setLoadBiaxial(exp.total_axial_count,
#                    [0,period/4.0,period/2.0,period/4.0*3.0,period],
#                    [temperature_mean,temperature_max,temperature_mean,temperature_min,temperature_mean],
#                    [0,axial_strain,0,-1*axial_strain,0],
#                    [-1*shear_strain,0,shear_strain,0,-1*shear_strain])
#==============================================================================
# Proportional path TMF IP
#==============================================================================
#load.setLoadBiaxial(exp.total_axial_count,
#                    [0,period/4.0,period/2.0,period/4.0*3.0,period],
#                    [temperature_mean,temperature_max,temperature_mean,temperature_min,temperature_mean],
#                    [0,axial_strain,0,-1*axial_strain,0],
#                    [0,shear_strain,0,-1*shear_strain,0])
#==============================================================================
# Uniaxial TMF IP
#==============================================================================
#load.setLoadBiaxial(exp.total_axial_count,
#                    [0,period/4.0,period/2.0,period/4.0*3.0,period],
#                    [temperature_mean,temperature_max,temperature_mean,temperature_min,temperature_mean],
#                    [0,axial_strain,0,-1*axial_strain,0],
#                    [0,shear_strain,0,-1*shear_strain,0])
#==============================================================================
# Uniaxial TMF IP
#==============================================================================
load.setLoadBiaxial(exp.total_axial_count,
                    [0,period/4.0,period/2.0,period/4.0*3.0,period],
                    [temperature_mean,temperature_min,temperature_mean,temperature_max,temperature_mean],
                    [0,axial_strain,0,-1*axial_strain,0],
                    [0,shear_strain,0,-1*shear_strain,0])
#==============================================================================
# Step
#==============================================================================
#step = Step(predefined_temperature = int(exp.initial_temperature), 
#          time_period = int(load.total_runing_time), initial_inc = 0.005, 
#          min_inc = 0.0001, max_inc = 5, nonlinear = 'ON')
#step = Step(predefined_temperature = temperature_mean, 
#          time_period = int(load.total_runing_time), initial_inc = 0.005, 
#          min_inc = 0.0001, max_inc = 5, nonlinear = 'ON')
step = Step(predefined_temperature = temperature_mean, 
          time_period = int(load.total_runing_time), initial_inc = 0.005, 
          min_inc = 0.0001, max_inc = period/40.0, nonlinear = 'ON')
#==============================================================================
# UMAT
#==============================================================================
umat = UMAT(UMATDirectory = 'F:\\UMAT\\CurrentVersion\\', UMATMainFile = 'MAIN_IN718.for')
#==============================================================================
# Job
#==============================================================================
job = Job(JobName=name, UMAT=umat, Step=step, Load=load)
#==============================================================================
# plot
#==============================================================================
sim = SimulationData(job.CSVFullName,period)



nth = 10
xitem = 'axial_stress'
yitem = 'shear_stress'
zitem = 'temperature'
xitem = 'axial_strain'
yitem = 'axial_stress'
#fig, ax = plt.subplots()
x1 = sim.obtainNthCycle(xitem,nth)
y1 = sim.obtainNthCycle(yitem,nth)
z1 = sim.obtainNthCycle(zitem,nth)

x2 = exp.obtainNthCycle(xitem,nth)
y2 = exp.obtainNthCycle(yitem,nth)
z2 = exp.obtainNthCycle(zitem,nth)

plot_format()

multicolored_lines(x1,y1,z1)
multicolored_lines(x2,y2,z2)


#plt.gca().set_aspect('auto')
plt.xlabel(xylabels[xitem])
plt.ylabel(xylabels[yitem])
plt.show()