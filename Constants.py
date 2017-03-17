# -*- coding: utf-8 -*-
"""
Created on Thu Jan 05 11:50:47 2017

@author: j.Sun
"""

WorkbenchDirectory = 'F:\\UMAT\\'

UMATDirectory = 'F:\\UMAT\\CurrentVersion\\'

AbaqusTempDirectory = 'F:\\Temp\\IN7183\\'

PythonDirectiory = 'F:\\GitHub\\python\\'

InputDirectiory = WorkbenchDirectory + 'Input\\'

ExperimentDirectiory = 'F:\\Database\\IN718\\Timed\\'

SimulationDirectiory = 'F:\\Temp\\IN718_Sim\\'

ExperimentLogFile = 'Inconel718_test_log.csv'

FatigueDirectiory = 'F:\\Database\\Fatigue\\'

InputTemplate = 'Tension3DTemplate.inp'

xylabels = {'axial_count':'N, cycle',\
'runing_time':'Running Time, s',\
'temprature':'Temperature, $^{\circ}$C',\
'axial_disp':'Axial Displacement, mm',\
'axial_force':'Axial Force, N',\
'axial_strain':'Axial Strain $\\varepsilon$, %',\
'axial_stress':'Axial Stress $\sigma$, MPa',\
'rotation':'Rotation, $^{\circ}$',\
'torque':'Torque, N$\cdot$m',\
'shear_strain':'Engineering Shear Strain $\gamma$, %',\
'shear_stress':'Shear Stress $\\tau$, MPa',\
'eqpl_strain':'Equvialent Plastic Strain $p$',\
'thermal_strain':'Thermal Strain, %',\
'total_strain':'Total Strain, %',\
'shear_stress_eq':'Equvialent Shear Stress $\sqrt{3}\\tau$, MPa',\
'shear_strain_eq':'Equvialent Shear Strain $\\gamma/\sqrt{3}$, %',\
'axial_log_strain':'Axial Logarithmic Strain $\\varepsilon$, %',\
'axial_true_stress':'Axial True Stress $\sigma$, MPa',\
}

experiment_type_list = []
experiment_type_list.append(['TC-IP',['7031','7047','7030','7018']])
experiment_type_list.append(['TC-OP',['7033','7048','7032','7017']])
experiment_type_list.append(['PRO-IP',['7040','7029','7039','7038']])
experiment_type_list.append(['NPR-IP',['7036','7034','7045','7046','7028','7037']])
experiment_type_list.append(['TC-90',['7025']])
experiment_type_list.append(['TC-IF',['7110','7111','7112','7113','7114','7115']])

#xylabels = {'axial_count':'N [cycle]',\
#'runing_time':'Running Time [s]',\
#'temprature':'Temperature [$^{\circ}$C]',\
#'axial_disp':'Axial Displacement [mm]',\
#'axial_force':'Axial Force [N]',\
#'axial_strain':'Axial Strain $\\varepsilon$ [mm/mm]',\
#'axial_stress':'Axial Stress $\sigma$ [MPa]',\
#'rotation':'Rotation [$^{\circ}$]',\
#'torque':'Torque [N$\cdot$m]',\
#'shear_strain':'Engineering Shear Strain $\gamma$ [-]',\
#'shear_stress':'Shear Stress $\\tau$ [MPa]',\
#'eqpl_strain':'Equvialent Plastic Strain $p$',\
#'thermal_strain':'Thermal Strain [mm/mm]',\
#'total_strain':'Total Strain [mm/mm]',\
#'shear_stress_eq':'Equvialent Shear Stress $\sqrt{3}\\tau$ [MPa]',\
#'shear_strain_eq':'Equvialent Shear Strain $\\gamma/\sqrt{3}$',\
#'axial_log_strain':'Axial Logarithmic Strain $\\varepsilon$',\
#'axial_true_stress':'Axial True Stress $\sigma$ [MPa]',\
#}