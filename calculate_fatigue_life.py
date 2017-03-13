# -*- coding: utf-8 -*-
"""
Created on Thu Jan 05 11:50:47 2017

@author: j.Sun
"""

import xlsxwriter
from Node import Node
from Material import Material,material_in718
from Data import SimulationData,ExperimentData,ExperimentLog
from Constants import *
from Work import *

#==============================================================================
# calculate_data_fatigue_life
#==============================================================================
def calculate_data_fatigue_life(data,material):
#    nodelabel = data.node_label[0]
    nodelabel = 36
    nth = data.axial_count_index_list[-2]
    time = data.obtainNthCycle('runing_time',nth)
    temperature = data.obtainNthCycle('temperature',nth)
    length = len(time)
    s11 = data.obtainNthCycle('axial_stress',nth)
    s22 = np.zeros(length)
    s33 = np.zeros(length)
    s12 = data.obtainNthCycle('shear_stress',nth)
    s13 = np.zeros(length)
    s23 = np.zeros(length)
    e11 = data.obtainNthCycle('axial_strain',nth)
    e22 = e11*-1.0*material.poisson_ratio
    e33 = e11*-1.0*material.poisson_ratio
    e12 = data.obtainNthCycle('shear_strain',nth)/2.0
    e13 = np.zeros(length)
    e23 = np.zeros(length)
    stress=[]
    strain=[]
    for i in range(len(time)):
        stress.append([[s11[i],s12[i],s13[i]],[s12[i],s22[i],s23[i]],[s13[i],s23[i],s33[i]]])
        strain.append([[e11[i],e12[i],e13[i]],[e12[i],e22[i],e23[i]],[e13[i],e23[i],e33[i]]])
    node = Node(nodelabel=nodelabel, dimension=2, time=time, coordinate=[], 
                displacement=[], stress=stress, strain=strain, temperature=temperature)
    fatigue_data = node.fatigueLifeFSModel(material)
    return fatigue_data
#==============================================================================
# calculate_exp_fatigue_life
#==============================================================================
def calculate_exp_fatigue_life(exp,material):
    return calculate_data_fatigue_life(exp,material)
#==============================================================================
# calculate_exp_fatigue_life
#==============================================================================
def calculate_sim_fatigue_life(sim,material):
    return calculate_data_fatigue_life(sim,material)
#==============================================================================
# calculate_fatigue_life
#==============================================================================
def calculate_fatigue_life(material=material_in718()):
    material.show()
    experiment_log = ExperimentLog(ExperimentLogFile)
    
    OutputDirectiory = 'F:\\Database\\IN718\\Fatigue\\'

    headers = 'Number of Cycles to Failure N\-(f),Mises Equivalent Strain Amplitude \i(\g(De))\-(eq)/2,Stress Amplitude e \i(\g(Ds))/2,Specimen,Critical Plane,sigma_n_max,delta_sigma,delta_epsilon,tau_n_max,delta_tau,delta_gamma,Predicted Fatigue Lifetime N\-(p),Fatigue Coefficient,Temperature'
    units = 'cycles,mm/mm,Mpa,-,deg,Mpa,Mpa,mm/mm,Mpa,Mpa,mm/mm,cycles,-,C'   
    
    workbook = xlsxwriter.Workbook(OutputDirectiory + 'FS' + '.xlsx') # write to excel
    
    for tmf_test in tmf_tests:
        resultfile = open(OutputDirectiory + tmf_test[0] + '.csv', 'w') # write to csv
        print >>resultfile, headers # write to csv
        print >>resultfile, units # write to csv
        
        worksheet = workbook.add_worksheet(tmf_test[0]) # write to excel
        bold = workbook.add_format({'bold': 1}) # write to excel
        row_number = 1 # write to excel
        header_list = headers.split(',') # write to excel
        unit_list = units.split(',') # write to excel
        comment_list = [ tmf_test[0] for i in range(len(header_list))] # write to excel
        worksheet.write_row('A'+str(row_number), header_list, bold); row_number += 1 # write to excel
        worksheet.write_row('A'+str(row_number), unit_list, bold); row_number += 1 # write to excel
        worksheet.write_row('A'+str(row_number), comment_list, bold); row_number += 1 # write to excel
        
        for name in tmf_test[1]:
            regular = r'\d+\.?\d*'
            period = float(experiment_log.obtainItem(name,'period',regular)[0])
            expriment_life = int(experiment_log.obtainItem(name,'comments',regular)[0])
            equivalent_strain = float(experiment_log.obtainItem(name,'equivalent_strain',regular)[0])
            sim = SimulationData(AbaqusTempDirectory+name+'//'+name+'.csv',period)
            exp = ExperimentData(ExperimentDirectiory+name+'.csv')
            data = calculate_sim_fatigue_life(sim,material)
            
            line = '' # write to csv
            line += '%s,' % (expriment_life) # write to csv
            line += '%s,' % (equivalent_strain) # write to csv
            line += '%s,' % (0) # write to csv
            line += '%s,' % (name) # write to csv
            for d in data: # write to csv
                line += '%s,' % (d) # write to csv
            print >>resultfile, line # write to csv
            
            data_list = [expriment_life,equivalent_strain,0,name] + data # write to excel
            worksheet.write_row('A'+str(row_number), data_list); row_number += 1 # write to excel
            
        resultfile.close() # write to csv
        
    workbook.close() # write to excel

calculate_fatigue_life()