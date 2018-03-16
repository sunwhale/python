# -*- coding: utf-8 -*-
"""
Created on Thu Jan 05 11:50:47 2017

@author: j.Sun
"""

import xlsxwriter
from Node import Node
from Material import Material,material_in718,material_cmsx4
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
    heat_flux_1 = data.obtainNthCycle('heat_flux_1',nth)
    heat_flux_2 = data.obtainNthCycle('heat_flux_2',nth)
    heat_flux_3 = data.obtainNthCycle('heat_flux_3',nth)
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
    if heat_flux_1 == []:
        heat_flux_1 = np.zeros(length)
    if heat_flux_2 == []:
        heat_flux_2 = np.zeros(length)
    if heat_flux_3 == []:
        heat_flux_3 = np.zeros(length)
    stress = []
    strain = []
    heatflux = []
    for i in range(len(time)):
        stress.append([[s11[i],s12[i],s13[i]],[s12[i],s22[i],s23[i]],[s13[i],s23[i],s33[i]]])
        strain.append([[e11[i],e12[i],e13[i]],[e12[i],e22[i],e23[i]],[e13[i],e23[i],e33[i]]])
        heatflux.append([heat_flux_1[i],heat_flux_2[i],heat_flux_3[i]])
    node = Node(nodelabel=nodelabel, dimension=2, time=time, coordinate=[], 
                displacement=[], stress=stress, strain=strain,
                temperature=temperature, heatflux=heatflux)
#    if fatigue_model == 'FS':
#        fatigue_data = node.fatigueLifeFSModel(material,k=1.0)
    if fatigue_model == 'FS':
        life = nth*2.0
        k = ((material.tau_f/material.shear_modulus*(2.0*life)**material.b0 + material.gamma_f*(2.0*life)**material.c0)/((1.0+material.poisson_ratio)*material.sigma_f/material.youngs_modulus*(2.0*life)**material.b+(1.0+0.5)*material.epsilon_f*(2.0*life)**material.c)-1.0)*material.yield_stress/(material.sigma_f*(2.0*life)**material.b)
        print k
        fatigue_data = node.fatigueLifeFSModel(material,k=k)
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
    if fatigue_model == 'Our':
        fatigue_data = node.fatigueLifeOurModel(material)
    return fatigue_data, node

#==============================================================================
# calculate_fatigue_life
#==============================================================================
def calculate_fatigue_life(fatigue_model,material=material_in718()):
    material.show()
    experiment_log = ExperimentLog(ExperimentLogFile)
    
#    FatigueDirectory = 'F:\\Database\\Fatigue\\%s\\' % fatigue_model
#    FatigueDirectory = 'F:\\Database\\Fatigue\\'

    if not os.path.isdir(FatigueDirectory):
        os.makedirs(FatigueDirectory)
        print 'Create new directory:',FatigueDirectory
    
    headers = ''
    units = ''
    headers += 'Number of Cycles to Failure N\-(f),';                   units +='cycle,'
    headers += 'Mises Equivalent Strain Amplitude \i(\g(De))\-(eq)/2,'; units +='mm/mm,'
    headers += 'Stress Amplitude e \i(\g(Ds))/2,';                      units +='MPa,'
    headers += 'Specimen,';                                             units +='-,'
    headers += 'Critical Plane,';                                       units +='deg,'
    headers += 'sigma_n_max,';                                          units +='MPa,'
    headers += 'delta_sigma,';                                          units +='MPa,'
    headers += 'delta_epsilon,';                                        units +='mm/mm,'
    headers += 'tau_n_max,';                                            units +='MPa,'
    headers += 'delta_tau,';                                            units +='MPa,'
    headers += 'delta_gamma,';                                          units +='-,'
    headers += 'Predicted Fatigue Lifetime N\-(p),';                    units +='cycle,'
    headers += 'Fatigue Coefficient,';                                  units +='-,'
    headers += 'Temperature,';                                          units +='C,'
    headers += 'Temperature Gradient 1,';                               units +='C/mm,'
    headers += 'Temperature Gradient 2,';                               units +='C/mm,'
    headers += 'Load Type,';                                            units +='-,'
    
    workbook = xlsxwriter.Workbook(FatigueDirectory + fatigue_model + '.xlsx') # write to excel
    
    allresultfile = open(FatigueDirectory + fatigue_model + '.csv', 'w') # write to csv all
    print >>allresultfile, headers[:-1] # write to csv all
    print >>allresultfile, units[:-1] # write to csv all
        
    for experiment_type in experiment_type_list:
#        resultfile = open(FatigueDirectory + experiment_type[0] + '.csv', 'w') # write to csv
#        print >>resultfile, headers # write to csv
#        print >>resultfile, units # write to csv
        
        worksheet = workbook.add_worksheet(experiment_type[0]) # write to excel
        bold = workbook.add_format({'bold': 1}) # write to excel
        row_number = 1 # write to excel
        header_list = headers.split(',') # write to excel
        unit_list = units.split(',') # write to excel
        comment_list = [ experiment_type[0] for i in range(len(header_list))] # write to excel
        worksheet.write_row('A'+str(row_number), header_list, bold); row_number += 1 # write to excel
        worksheet.write_row('A'+str(row_number), unit_list, bold); row_number += 1 # write to excel
        worksheet.write_row('A'+str(row_number), comment_list, bold); row_number += 1 # write to excel
        
        for name in experiment_type[1]:
            print name
            regular = r'\d+\.?\d*'
            period = float(experiment_log.obtainItem(name,'period',regular)[0])
            expriment_life = int(experiment_log.obtainItem(name,'comments',regular)[0])
            equivalent_strain = float(experiment_log.obtainItem(name,'equivalent_strain',regular)[0])
            """
            使用计算模拟结果。
            """
#            sim = SimulationData(SimulationDirectory+name+'.csv',period)
#            data, node = calculate_data_fatigue_life(sim,material,fatigue_model)
            """
            使用试验结果。
            """
            exp = ExperimentData(ExperimentDirectory+name+'.csv')
            data, node = calculate_data_fatigue_life(exp,material,fatigue_model)
            
            line = '' # write to csv
            line += '%s,' % (expriment_life) # write to csv
            line += '%s,' % (equivalent_strain) # write to csv
            line += '%s,' % (0) # write to csv
            line += '%s,' % (name) # write to csv
            for d in data[1:]: # write to csv
                line += '%s,' % (d) # write to csv
#            print >>resultfile, line[:-1] # write to csv, ignore the last comma
            
            line += '%s' % (experiment_type[0]) # write to csv all
            print >>allresultfile, line # write to csv all
            
            data_list = [expriment_life,equivalent_strain,0,name] + data # write to excel
            worksheet.write_row('A'+str(row_number), data_list); row_number += 1 # write to excel
            
#        resultfile.close() # write to csv
        
    allresultfile.close() # write to csv all
    workbook.close() # write to excel

if __name__ == '__main__':
#    fatigue_model_list = ['BM','FS','SWT','Liu1','Liu2','Chu','Our']
    fatigue_model_list = ['BM']
    for fatigue_model in fatigue_model_list:
#        calculate_fatigue_life(fatigue_model,material=material_in718())
        calculate_fatigue_life(fatigue_model,material=material_cmsx4())