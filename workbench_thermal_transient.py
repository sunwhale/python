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
from Work_thermal_transient import Step,UMAT,Load,Job
from Functions import copy_file,read_json_file,write_json_file
from ambient_temperature import *

def workbench(name,loading_cycles=None,copy=True,film_coefficient=0.0,sink_temperature=0.0,
              temperature_list=[],thermal_strain_list=[0.0,-0.0],heat_flux=0.0,
              film_coefficient_outer=0.0,film_coefficient_inner=0.0,emissivity=0.80,
              sink_temperature_inner=293.15,sink_temperature_outer=293.15,outer_temperature=293.15):
    """
    某试件对应的边界条件下的数值模拟。
    """
    experiment_log = ExperimentLog(ExperimentLogFile)
    experiment_log.output(name)
    regular = r'.*'
    load_type = experiment_log.obtainItem(name,'load_type',regular)[0]
    regular = r'\d+\.?\d*'
    temperature_mode = experiment_log.obtainItem(name,'temperature_mode',regular)
    if temperature_list == []:
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
        loading_cycles = min(int(predicted_life/4.0),5000)
    use_exp_data = True
    thermal_strain_min = min(thermal_strain_list)/100.0
    thermal_strain_max = max(thermal_strain_list)/100.0
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
                            [0,axial_strain*1,0,-1*axial_strain*1,0],
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
# load of convection
#==============================================================================
    load.setConvection(film_coefficient,sink_temperature)
#==============================================================================
# load of thermal
#==============================================================================
    load.setThermal(heat_flux,film_coefficient_outer,film_coefficient_inner,emissivity,sink_temperature_inner,sink_temperature_outer,outer_temperature)
#==============================================================================
# Step
#==============================================================================
#    step = Step(predefined_temperature = int(exp.initial_temperature), 
#              time_period = int(load.total_runing_time), initial_inc = 0.005, 
#              min_inc = 0.0001, max_inc = 5, nonlinear = 'ON')
    step = Step(predefined_temperature = 20.0, 
              time_period = int(load.total_runing_time), initial_inc = 0.00001, 
              min_inc = 0.00000001, max_inc = period/40.0, nonlinear = 'OFF')
#==============================================================================
# UMAT
#==============================================================================
#    umat = UMAT(UMATDirectory = 'F:\\GitHub\\umat\\', 
#                UMATMainFile = 'MAIN_IN718.for', 
#                ParameterFortranFile = 'PARAMETERS_IN718_TMF.for',
#                OutputFortranFile = 'OUTPUT.for',
#                OutputTextFile = name + '_output.txt')
                
    umat = UMAT(UMATDirectory = 'F:\\GitHub\\umat\\', 
                UMATMainFile = 'MAIN_IN718.for', 
                ParameterFortranFile = 'PARAMETERS_SS304.for',
                OutputFortranFile = 'OUTPUT.for',
                OutputTextFile = name + '_output.txt')
                
#    umat = UMAT(UMATDirectory = 'F:\\UMAT\\Fangjie\\', 
#                UMATMainFile = 'AO304Tanakav6VariDamageMyisokin3.for', 
#                ParameterFortranFile = 'PARAMETERS_SS304.for',
#                OutputFortranFile = 'OUTPUT.for',
#                OutputTextFile = name + '_output.txt')
                
#    umat = UMAT(UMATDirectory = 'F:\\UMAT\\SS304\\', 
#                UMATMainFile = 'MAIN_SS304.for', 
#                ParameterFortranFile = 'PARAMETERS_SS304.for',
#                OutputFortranFile = 'OUTPUT.for',
#                OutputTextFile = name + '_output.txt')
#==============================================================================
# Job
#==============================================================================
    job = Job(JobName=name, UMAT=umat, Step=step, Load=load, copy=copy)

#    job.allProc()
    job.createDirectory()
    job.copyFiles()
    job.creatBatchFile()
    job.createAbaqusCAE()
    job.createAbaqusInput()
    job.run()
    job.autoPostProc()


def calculate_film_coefficient(volume_flow):
    if volume_flow == 0.0:
        return 0.02
    air_temperature = 20.0 + 273.15
    pressure = 6.5e5
    dynamic_viscosity_0 = 1.7894e-5
    rin = 0.0065/2.0
    d = rin*2
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
    film_coefficient_inner = thermal_conductivity/d*Nu/1000.0
    return film_coefficient_inner
        
absolute_zero = -273.15
predefined_temperature = 20
d_out = 8.5e-3
height = 80.0e-3
total_power = 6400.0
reflect = 0.2
emissivity = 0.80
power_percent_exp_list = [0.3,0.5,0.7,0.9]
exp_list = []
outer_temperature_list = [
215,
338,
495,
654,
]
exp_list.append([power_percent_exp_list,outer_temperature_list,67])
outer_temperature_list = [
277,
435,
603,
757.5,
]
exp_list.append([power_percent_exp_list,outer_temperature_list,50])
outer_temperature_list = [
400,
575,
730,
870,
]
exp_list.append([power_percent_exp_list,outer_temperature_list,25])
outer_temperature_list = [
440,
602.5,
759.5,
899,
]
exp_list.append([power_percent_exp_list,outer_temperature_list,0])

print exp_list
plot_list = []
out_list = []
threshold = 0.20
name = '0001'
AbaqusWorkDirectory = AbaqusTempDirectory + name + '\\'

for exp in exp_list[:]:
    result_list = []
    for i in range(4):
        outer_temperature = exp[1][i]
        power_percent_exp = exp[0][i]
        if power_percent_exp == 0.3:
            write_json_file(AbaqusWorkDirectory+'ambient_temperature.txt',ambient_temperature_30)
        if power_percent_exp == 0.5:
            write_json_file(AbaqusWorkDirectory+'ambient_temperature.txt',ambient_temperature_50)
        if power_percent_exp == 0.7:
            write_json_file(AbaqusWorkDirectory+'ambient_temperature.txt',ambient_temperature_70)
        if power_percent_exp == 0.9:
            write_json_file(AbaqusWorkDirectory+'ambient_temperature.txt',ambient_temperature_90)
        volume_flow = exp[2]
        film_coefficient_inner = calculate_film_coefficient(volume_flow) * 1.0
        film_coefficient_outer = 0.005
#        sink_temperature_inner = 50/0.6*(power_percent_exp-threshold)*(1.0)/(1.0-threshold)+20+273.15
        sink_temperature_inner = 20
#        sink_temperature_outer = 500/0.9*(power_percent_exp-threshold)*(1.0)/(1.0-threshold)+20+273.15
        sink_temperature_outer = 20
        heat_flux = emissivity*reflect*total_power*(power_percent_exp-threshold)*(1.0)/(1.0-threshold)/np.pi/d_out/height*1.0e-3
        workbench(name,loading_cycles=1,
                  heat_flux=heat_flux,
                  film_coefficient_outer=film_coefficient_outer,
                  film_coefficient_inner=film_coefficient_inner,
                  emissivity=emissivity,
                  sink_temperature_inner=sink_temperature_inner,
                  sink_temperature_outer=sink_temperature_outer,
                  outer_temperature=outer_temperature)
        sim_filename = SimulationDirectory + name + '.csv'
        out_filename = 'F:\\Cloud\\Simulation\\Temperature\\%s_%s_.csv' % (int(power_percent_exp*100),volume_flow)
        copy_file(sim_filename,out_filename)
        simulation = SimulationData(sim_filename,1)
        temperature = simulation.temperature[-1] + 273.15
        heat_flux_1 = simulation.heat_flux_1[-1]
        heat_flux_2 = simulation.heat_flux_2[-1]
        power_percent_cal = -1.0*heat_flux_1*np.pi*8.5e-3*80e-3/6.4
        efficiency = power_percent_cal/power_percent_exp
        result_list.append([power_percent_exp,film_coefficient_inner,temperature,heat_flux_1])
        out_list.append([exp[2],power_percent_exp,film_coefficient_inner,temperature,heat_flux_1,power_percent_cal,efficiency])
    plot_list.append([exp[2],result_list])

#AbaqusWorkDirectory = AbaqusTempDirectory + name + '\\'
#
#out_name = '%s_%s_%s_%s' % (film_coefficient_inner,sink_temperature_inner,film_coefficient_outer,sink_temperature_outer)
#
#for plot in plot_list:
#    p = np.array(plot[1])
#    plt.plot(p[:,0],p[:,3],label=plot[0],marker='s')
#plt.legend(loc=0)
#plt.savefig(AbaqusWorkDirectory + out_name + '.png', dpi=150, transparent=False)
#plt.show()
#
#outfile = open(AbaqusWorkDirectory + out_name + '.csv', 'w')
#for r in out_list:
#    print >>outfile, '%s,%s,%s,%s,%s,%s,%s' % (r[0],r[1],r[2],r[3],r[4],r[5],r[6])
#outfile.close()