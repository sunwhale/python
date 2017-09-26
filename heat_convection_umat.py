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

air_temperature = 20.0 + 273.15
pressure = 6.5e5
dynamic_viscosity_0 = 1.7894e-5
rin = 0.0065/2.0
d = rin*2

#outfile = open('heat_convection.csv', 'w')

for name in ['7201']:
    for volume_flow in [40]:
        for out_face_temperature in [650]:
            velocity = volume_flow/60.0/1000.0/(np.pi*rin**2)
            density = pressure/287.058/air_temperature
            thermal_conductivity = -0.00037+0.000103*air_temperature-4.657E-08*air_temperature**2
            specific_heat = 1070.3-0.564*air_temperature+0.001507*air_temperature**2-0.000001102*air_temperature**3-0.000000014*air_temperature**4
            dynamic_viscosity = dynamic_viscosity_0 * (air_temperature/288.15)**(1.5)*(288.15+110.4)/(air_temperature+110.4)
            kinematic_viscosity = dynamic_viscosity/density
            Re = density * velocity * d / dynamic_viscosity
            Pr = specific_heat*dynamic_viscosity/thermal_conductivity
            #Gnielinski correlation
            f = (0.79*np.log(Re)-1.64)**(-2.0)
            Nu = (f/8.0)*(Re-1000.0)*Pr/(1.0+12.7*(f/8.0)**0.5*(Pr**(2.0/3.0)-1))
            #Dittus-Boelter equation
            #Nu = 0.023*Re**(4.0/5.0)*Pr**0.4           
        
            film_coefficient = thermal_conductivity/d*Nu/1000.0
            sink_temperature = air_temperature - 273.15
            temperature_list = [650,300]
            workbench(name,loading_cycles=10,copy=True,
                      film_coefficient=film_coefficient,
                      sink_temperature=sink_temperature,
                      temperature_list=temperature_list)
                      
            sim_filename = SimulationDirectory + name + '.csv'
            simulation = SimulationData(sim_filename,period=240)
            power = simulation.heat_flux_1[-1]*np.pi*6.5*80.0/1000.0
#            print >>outfile, '%s,%s,%s,%s,%s,%s,%s' % (out_face_temperature,volume_flow,Nu,film_coefficient,simulation.temperature[-1],simulation.heat_flux_1[-1],power)
            
#outfile.close()