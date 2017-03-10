# -*- coding: utf-8 -*-
"""
Created on Thu Jan 05 11:50:47 2017

@author: j.Sun
"""

from Node import Node
from Material import Material,material_in718
from Data import SimulationData,ExperimentData,ExperimentLog
from Constants import *
from Work import *

#==============================================================================
# job
#==============================================================================
def calculate_fatigue_life(name,material=material_in718()):
    material.show()
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
    
    exp = ExperimentData()
    sim = SimulationData(AbaqusTempDirectory+name+'//'+name+'.csv',period)
    nodelabel = sim.node_label[0]
    nth = sim.axial_count_index_list[-2]
    time = sim.obtainNthCycle('runing_time',nth)
    temperature = sim.obtainNthCycle('temperature',nth)
    length = len(time)
    s11 = sim.obtainNthCycle('axial_stress',nth)
    s22 = np.zeros(length)
    s33 = np.zeros(length)
    s12 = sim.obtainNthCycle('shear_stress',nth)
    s13 = np.zeros(length)
    s23 = np.zeros(length)
    e11 = sim.obtainNthCycle('axial_strain',nth)
    e22 = e11*-1.0*material.poisson_ratio
    e33 = e11*-1.0*material.poisson_ratio
    e12 = sim.obtainNthCycle('shear_strain',nth)/2.0
    e13 = np.zeros(length)
    e23 = np.zeros(length)
    
    stress=[]
    strain=[]
    for i in range(len(time)):
        stress.append([[s11[i],s12[i],s13[i]],[s12[i],s22[i],s23[i]],[s13[i],s23[i],s33[i]]])
        strain.append([[e11[i],e12[i],e13[i]],[e12[i],e22[i],e23[i]],[e13[i],e23[i],e33[i]]])
    
    node = Node(nodelabel=nodelabel, dimension=2, time=time, coordinate=[], 
                displacement=[], stress=stress, strain=strain, temperature=temperature)

    fatigue_life = node.fatigueLifeFSModel(material)
    
#    resultfile = open('out0.dat', 'w')
    line = ''
    line += '%s    ' % (node.nodelabel)
    line += '%s'     % (fatigue_life)
#    print >>resultfile, line
#    resultfile.close()
    
#for csv_files in all_files[:]:
#    resultfile = open(OutputDirectiory+csv_files[0]+'.csv', 'w')
#    Headers = 'Number of Cycles to Failure N\-(f),Mises Equivalent Strain Amplitude \i(\g(De))\-(eq)/2,Stress Amplitude e \i(\g(Ds))/2,Specimen,Critical Plane,sigma_n_max,delta_sigma,delta_epsilon,tau_n_max,delta_tau,delta_gamma,Predicted Fatigue Lifetime N\-(p),Fatigue Coefficient,Temperature'
#    Units = 'cycles,mm/mm,Mpa,-,deg,Mpa,Mpa,mm/mm,Mpa,Mpa,mm/mm,cycles,-,C'
#    print >>resultfile, Headers
#    print >>resultfile, Units
name = '7037'
calculate_fatigue_life(name)