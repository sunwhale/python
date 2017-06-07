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

def plot_exp_vs_sim_pv(name,item='axial_stress',save_types=['.png','.pdf']):
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

    sim_filename = SimulationDirectory + name + '.csv'
    simulation = SimulationData(sim_filename,period)
    cycle,peak,valley = simulation.obtainPeakValley(item)
    plt.plot(cycle,peak,label='sim')
    plt.plot(cycle,valley,label='sim')
        
    exp_filename = ExperimentDirectory + name + '.csv'
    experiment = ExperimentData(exp_filename)
    cycle,peak,valley = experiment.obtainPeakValley(item)
    plt.plot(cycle,peak,label='exp')
    plt.plot(cycle,valley,label='exp')
    
    plt.xscale('log')
#    plt.yscale('log')
    
    plt.legend(loc=0)
    
    figure_path = ArticleFigureDirectory
    figure_name = name + '_' + item + '_pv'
    if not os.path.isdir(ArticleFigureDirectory):
        os.makedirs(ArticleFigureDirectory)
        print 'Create new directory:',ArticleFigureDirectory
    if figure_path <> None and figure_name<> None:
        for save_type in save_types:
            plt.savefig(figure_path + figure_name + save_type, dpi=150)
            print 'save as', figure_path + figure_name + save_type
            
#    plt.show()
    plt.close()