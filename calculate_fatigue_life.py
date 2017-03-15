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
def calculate_data_fatigue_life(data,material,fatigue_model):
    nodelabel = data.node_label[0]
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
    if fatigue_model == 'FS':
        fatigue_data = node.fatigueLifeFSModel(material,k=1.0)
    if fatigue_model == 'SWT':
        fatigue_data = node.fatigueLifeSWTModel(material)
    if fatigue_model == 'BM':
        fatigue_data = node.fatigueLifeBMModel(material,S=0.36)
    if fatigue_model == 'Liu1':
        fatigue_data = node.fatigueLifeLiu1Model(material)
    if fatigue_model == 'Liu2':
        fatigue_data = node.fatigueLifeLiu2Model(material)
    if fatigue_model == 'Chu':
        fatigue_data = node.fatigueLifeChuModel(material)
    return fatigue_data

#==============================================================================
# calculate_fatigue_life
#==============================================================================
def calculate_fatigue_life(fatigue_model,material=material_in718()):
    material.show()
    experiment_log = ExperimentLog(ExperimentLogFile)
    
    OutputDirectiory = 'F:\\Database\\Fatigue\\%s\\' % fatigue_model

    if not os.path.isdir(OutputDirectiory):
        os.makedirs(OutputDirectiory)
        print 'Create new directory:',OutputDirectiory
            
    headers = 'Number of Cycles to Failure N\-(f),Mises Equivalent Strain Amplitude \i(\g(De))\-(eq)/2,Stress Amplitude e \i(\g(Ds))/2,Specimen,Critical Plane,sigma_n_max,delta_sigma,delta_epsilon,tau_n_max,delta_tau,delta_gamma,Predicted Fatigue Lifetime N\-(p),Fatigue Coefficient,Temperature'
    units = 'cycles,mm/mm,MPa,-,deg,MPa,MPa,mm/mm,MPa,MPa,mm/mm,cycles,-,C'
    
    workbook = xlsxwriter.Workbook(OutputDirectiory + fatigue_model + '.xlsx') # write to excel
    
    allresultfile = open(OutputDirectiory + fatigue_model + '.csv', 'w') # write to csv all
    print >>allresultfile, headers + ',Load Type' # write to csv all
    print >>allresultfile, units + ',-' # write to csv all
        
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
#            sim = SimulationData(AbaqusTempDirectory+name+'//'+name+'.csv',period)
            sim = SimulationData(SimulationDirectiory+name+'.csv',period)
            exp = ExperimentData(ExperimentDirectiory+name+'.csv')
            data = calculate_data_fatigue_life(sim,material,fatigue_model)
            
            line = '' # write to csv
            line += '%s,' % (expriment_life) # write to csv
            line += '%s,' % (equivalent_strain) # write to csv
            line += '%s,' % (0) # write to csv
            line += '%s,' % (name) # write to csv
            for d in data: # write to csv
                line += '%s,' % (d) # write to csv
            print >>resultfile, line[:-1] # write to csv, ignore the last comma
            
            line += '%s' % (tmf_test[0]) # write to csv all
            print >>allresultfile, line # write to csv all
            
            data_list = [expriment_life,equivalent_strain,0,name] + data # write to excel
            worksheet.write_row('A'+str(row_number), data_list); row_number += 1 # write to excel
            
        resultfile.close() # write to csv
        
    allresultfile.close() # write to csv all
    workbook.close() # write to excel

fatigue_model_list = ['BM','FS','SWT','Liu1','Liu2','Chu']
for fatigue_model in fatigue_model_list:
    calculate_fatigue_life(fatigue_model)