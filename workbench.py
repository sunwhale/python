# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 13:46:18 2017

@author: SunJingyu
"""

import numpy as np
import os
import matplotlib.pyplot as plt
from Constants import *
from Data import SimulationData,ExperimentData,ExperimentLog
from Material import Material
from Work import Step,UMAT,Load,Job

def workbench(name,loading_cycles=None,copy=True):
    """
    某试件对应的边界条件下的数值模拟。
    """
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
#==============================================================================
# material
#==============================================================================
    material = Material()
    material.setName(name='IN718')
    material.setTemperature(temperature=650.0)
    material.setMonotonic(youngs_modulus=167100.0,poisson_ratio=0.2886,yield_stress=1064.0)
    material.setCyclicAxial(sigma_f=1034.0,b=-0.04486,epsilon_f=0.11499,c=-0.52436)
    material.setCyclicTorsion(tau_f=1034.0/np.sqrt(3),b0=-0.04486,gamma_f=0.11499*np.sqrt(3),c0=-0.52436)
    predicted_life = material.calculateMansonCoffinLife(equivalent_strain/100.0)
#==============================================================================
# experiment
#==============================================================================
    exp_full_name = ExperimentDirectory + name + '.csv'
    if os.path.exists(exp_full_name):
        exp = ExperimentData(exp_full_name)
        experimental_life = exp.total_axial_count
    else:
        print ('%s is not existed' % exp_full_name)
#==============================================================================
# load
#==============================================================================
    temperature_mean = (temperature_list[0] + temperature_list[1])/2.0
    temperature_min = min(temperature_list)
    temperature_max = max(temperature_list)
    load = Load(runing_time=[0], temperature=[temperature_mean], axial_strain=[0], shear_strain=[0], first_cycle_shift=1)
    axial_strain = axial_strain/100.0
    shear_strain = np.deg2rad(angel_strain)*d_out/2.0/gauge_length
    if loading_cycles == None:
        loading_cycles = int(predicted_life/4.0)
    use_exp_data = True
#==============================================================================
# Diamond path TMF IP
#==============================================================================
    if load_type == 'cyclic diamond path' and axial_temperature_phase == 0.0:
        use_exp_data = False
        load.setLoadBiaxial(int(loading_cycles),
                            [0,period/4.0,period/2.0,period/4.0*3.0,period],
                            [temperature_mean,temperature_max,temperature_mean,temperature_min,temperature_mean],
                            [0,axial_strain,0,-1*axial_strain,0],
                            [-1*shear_strain,0,shear_strain,0,-1*shear_strain])
#==============================================================================
# Proportional path TMF IP
#==============================================================================
    if load_type == 'cyclic proportional path' and axial_temperature_phase == 0.0:
        use_exp_data = False
        load.setLoadBiaxial(int(loading_cycles),
                            [0,period/4.0,period/2.0,period/4.0*3.0,period],
                            [temperature_mean,temperature_max,temperature_mean,temperature_min,temperature_mean],
                            [0,axial_strain,0,-1*axial_strain,0],
                            [0,shear_strain,0,-1*shear_strain,0])
#==============================================================================
# Uniaxial TMF IP
#==============================================================================
    if load_type == 'cyclic tension compression' and axial_temperature_phase == 0.0:
        use_exp_data = False
        load.setLoadBiaxial(int(loading_cycles),
                            [0,period/4.0,period/2.0,period/4.0*3.0,period],
                            [temperature_mean,temperature_max,temperature_mean,temperature_min,temperature_mean],
                            [0,axial_strain,0,-1*axial_strain,0],
                            [0,shear_strain,0,-1*shear_strain,0])
#==============================================================================
# Uniaxial TMF OP
#==============================================================================
    if load_type == 'cyclic tension compression' and axial_temperature_phase == 180.0:
        use_exp_data = False
        load.setLoadBiaxial(int(loading_cycles),
                            [0,period/4.0,period/2.0,period/4.0*3.0,period],
                            [temperature_mean,temperature_min,temperature_mean,temperature_max,temperature_mean],
                            [0,axial_strain,0,-1*axial_strain,0],
                            [0,shear_strain,0,-1*shear_strain,0])
#==============================================================================
# Uniaxial TMF 90
#==============================================================================
    if load_type == 'cyclic tension compression' and axial_temperature_phase == 90.0:
        use_exp_data = False
        load.setLoadBiaxial(int(loading_cycles),
                            [0,period/4.0,period/2.0,period/4.0*3.0,period],
                            [temperature_min,temperature_mean,temperature_max,temperature_mean,temperature_min],
                            [0,axial_strain,0,-1*axial_strain,0],
                            [0,shear_strain,0,-1*shear_strain,0])
#==============================================================================
# load from experiment data
#==============================================================================
    if use_exp_data:
        if loading_cycles == None:
            load.setLoadFromExperiment(exp,runing_time=None)
        else:
            load.setLoadFromExperiment(exp,runing_time=period*loading_cycles)
#==============================================================================
# Step
#==============================================================================
#    step = Step(predefined_temperature = int(exp.initial_temperature), 
#              time_period = int(load.total_runing_time), initial_inc = 0.005, 
#              min_inc = 0.0001, max_inc = 5, nonlinear = 'ON')
    step = Step(predefined_temperature = temperature_mean, 
              time_period = int(load.total_runing_time), initial_inc = 0.005, 
              min_inc = 0.0001, max_inc = period/40.0, nonlinear = 'OFF')
#==============================================================================
# UMAT
#==============================================================================
#    umat = UMAT(UMATDirectory = 'F:\\GitHub\\umat\\', 
#                UMATMainFile = 'MAIN_IN718.for', 
#                ParameterFortranFile = 'PARAMETERS_IN718_TMF.for',
#                OutputFortranFile = 'OUTPUT.for',
#                OutputTextFile = name + '_output.txt')
                
#    umat = UMAT(UMATDirectory = 'F:\\GitHub\\umat\\', 
#                UMATMainFile = 'MAIN_IN718.for', 
#                ParameterFortranFile = 'PARAMETERS_SS304.for',
#                OutputFortranFile = 'OUTPUT.for',
#                OutputTextFile = name + '_output.txt')
                
    umat = UMAT(UMATDirectory = 'F:\\UMAT\\Fangjie\\', 
                UMATMainFile = 'AO304Tanakav6VariDamageMyisokin2.for', 
                ParameterFortranFile = 'PARAMETERS_SS304.for',
                OutputFortranFile = 'OUTPUT.for',
                OutputTextFile = name + '_output.txt')
                
#    umat = UMAT(UMATDirectory = 'F:\\UMAT\\SS304\\', 
#                UMATMainFile = 'MAIN_SS304.for', 
#                ParameterFortranFile = 'PARAMETERS_SS304.for',
#                OutputFortranFile = 'OUTPUT.for',
#                OutputTextFile = name + '_output.txt')
#==============================================================================
# Job
#==============================================================================
    job = Job(JobName=name, UMAT=umat, Step=step, Load=load, copy=copy)
    job.allProc()
#    job.createDirectory()
#    job.copyFiles()
#    job.creatBatchFile()
#    job.createAbaqusCAE()
#    job.createAbaqusInput()
#    job.run()
#    job.autoPostProc()

#==============================================================================
# experiment_log
#==============================================================================
#experiment_log = ExperimentLog(ExperimentLogFile)
#number_list = experiment_log.keyFilter("(self.temperature_mode == '300-650') & (self.calculate == '1') & (self.load_type == 'cyclic diamond path')")
#number_list = experiment_log.keyFilter("(self.calculate == '1') & (self.load_type == 'cyclic diamond path')")
#number_list += experiment_log.keyFilter("(self.calculate == '1') & (self.load_type == 'cyclic proportional path')")
#number_list += experiment_log.keyFilter("(self.calculate == '1') & (self.load_type == 'cyclic tension compression')")
#print number_list

#for name in number_list:
#    work(name)
#workbench('7116',loading_cycles=1000)