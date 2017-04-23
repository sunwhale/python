# -*- coding: utf-8 -*-
"""
Created on Thu Jan 05 11:50:47 2017

@author: j.Sun
"""

WorkbenchDirectory = 'F:\\UMAT\\'
UMATDirectory = 'F:\\GitHub\\umat\\'
UMATDirectory = 'F:\\UMAT\\CurrentVersion\\'
PythonDirectory = 'F:\\GitHub\\python\\'
InputDirectory = WorkbenchDirectory + 'Input\\'
InputTemplate = 'Tension3DTemplate.inp'

#==============================================================================
# Inconel 718
#==============================================================================
AbaqusTempDirectory = 'F:\\Temp\\IN7183\\'
ExperimentDirectory = 'F:\\Database\\IN718\\Timed\\'
SimulationDirectory = 'F:\\Temp\\IN718_Sim\\'
ExperimentLogFile = PythonDirectory + 'Inconel718_test_log.csv'
FatigueDirectory = 'F:\\Database\\Fatigue\\'
ArticleFigureDirectory = 'F:\\Cloud\\Article\\Fatigue\\Figs\\'

experiment_type_list = []
experiment_type_list.append(['TC-IP',['7031','7047','7030','7018']])
experiment_type_list.append(['TC-OP',['7033','7048','7032','7017']])
experiment_type_list.append(['PRO-IP',['7040','7029','7039','7038']])
experiment_type_list.append(['NPR-IP',['7036','7034','7045','7046','7028','7037']])
experiment_type_list.append(['TC-90',['7025']])
experiment_type_list.append(['TC-IF',['7110','7111','7112','7113','7114','7115','7116']])

experiment_type_dict = {}
experiment_type_dict['TC-IP']=['7031','7047','7030','7018']
experiment_type_dict['TC-OP']=['7033','7048','7032','7017']
experiment_type_dict['PRO-IP']=['7040','7029','7039','7038']
experiment_type_dict['NPR-IP']=['7036','7034','7045','7046','7028','7037']
experiment_type_dict['TC-90']=['7025']
experiment_type_dict['TC-IF']=['7110','7111','7112','7113','7114','7115','7116']

#==============================================================================
# CMSX4
#==============================================================================
#ExperimentDirectory = 'F:\\Database\\CMSX4\\Timed\\'
#ExperimentLogFile = PythonDirectory + 'CMSX4_test_log.csv'
#ArticleFigureDirectory = 'F:\\Database\\CMSX4\\'
#experiment_type_dict = {}
#experiment_type_dict['TC-IP']=['08','09','10']
#experiment_type_dict['TC-OP']=['11','12']
#experiment_type_dict['TC-IF']=['03','04','05','06','07']

#==============================================================================
# dictionary
#==============================================================================
xylabels = {'axial_count':'N [cycle]',
            'runing_time':'Running Time [s]',
            'temperature':'Temperature [$^{\circ}$C]',
            'axial_disp':'Axial Displacement [mm]',
            'axial_force':'Axial Force [N]',
            'axial_strain':'Axial Strain $\\varepsilon$ [%]',
            'axial_stress':'Axial Stress $\sigma$ [MPa]',
            'rotation':'Rotation [$^{\circ}$]',
            'torque':'Torque [N$\cdot$m]',
            'shear_strain':'Engineering Shear Strain $\gamma$ [%]',
            'shear_stress':'Shear Stress $\\tau$ [MPa]',
            'eqpl_strain':'Equvialent Plastic Strain $p$]',
            'thermal_strain':'Thermal Strain [%]',
            'total_strain':'Total Strain [%]',
            'shear_stress_eq':'Equvialent Shear Stress $\sqrt{3}\\tau$ [MPa]',
            'shear_strain_eq':'Equvialent Shear Strain $\\gamma/\sqrt{3}$ [%]',
            'axial_log_strain':'Axial Logarithmic Strain $\\varepsilon$ [%]',
            'axial_true_stress':'Axial True Stress $\sigma$ [MPa]',
            'axial_stress_amplitude':'Axial Stress Amplitude $\Delta\sigma/2$ [MPa]',
            'axial_strain_amplitude':'Axial Strain Amplitude $\Delta\\varepsilon/2$ [%]',
            'equivalent_strain_amplitude':'Equivalent Strain Amplitude $\Delta\\varepsilon_{\\rm{ep}}/2$ [%]',
            'experimental_life':'Experimental Fatigue Lifetime $N_{\\rm{f}}$',
            'predicted_life':'Predicted Fatigue Lifetime $N_{\\rm{p}}$',
            }