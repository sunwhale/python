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
    film_coefficient = 0.3
    sink_temperature = 20.0
    temperature_list = [650.0,650.0]
    workbench(name,loading_cycles=1,copy=True,
              film_coefficient=film_coefficient,
              sink_temperature=sink_temperature,
              temperature_list=temperature_list)
              
    sim_filename = SimulationDirectory + name + '.csv'
    simulation = SimulationData(sim_filename,period=1)
    print simulation.temperature
    print simulation.heat_flux_1