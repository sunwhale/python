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
#AbaqusTempDirectory = 'F:\\Temp\\IN7183\\'
#ExperimentDirectory = 'F:\\Database\\IN718\\Timed\\'
#SimulationDirectory = 'F:\\Temp\\IN718_Sim\\'
#SimulationDirectory = 'F:\\Temp\\IN718_Sim_2\\'
#ExperimentLogFile = PythonDirectory + 'Inconel718_test_log.csv'
#FatigueDirectory = 'F:\\Database\\Fatigue\\'
#ArticleFigureDirectory = 'F:\\Cloud\\Article\\Fatigue\\Figs\\'
#
#experiment_type_list = []
#experiment_type_list.append(['TC-IP',['7031','7047','7030','7018']])
#experiment_type_list.append(['TC-OP',['7033','7048','7032','7017']])
#experiment_type_list.append(['PRO-IP',['7040','7029','7039','7038']])
#experiment_type_list.append(['NPR-IP',['7036','7034','7045','7046','7028','7037']])
#experiment_type_list.append(['TC-90',['7025']])
#experiment_type_list.append(['TC-IF',['7110','7111','7112','7113','7114','7115','7116']])
#
#experiment_type_dict = {}
#experiment_type_dict['TC-IP']=['7031','7047','7030','7018']
#experiment_type_dict['TC-OP']=['7033','7048','7032','7017']
#experiment_type_dict['PRO-IP']=['7040','7029','7039','7038']
#experiment_type_dict['NPR-IP']=['7036','7034','7045','7046','7028','7037']
#experiment_type_dict['TC-90']=['7025']
#experiment_type_dict['TC-IF']=['7110','7111','7112','7113','7114','7115','7116']

#==============================================================================
# SS304
#==============================================================================
AbaqusTempDirectory = 'F:\\Temp\\SS304\\'
ExperimentDirectory = 'F:\\Database\\SS304\\Timed\\'
#SimulationDirectory = 'F:\\Temp\\SS304_Sim\\'
SimulationDirectory = 'F:\\Temp\\SS304_Sim_fang_20170607\\'
ExperimentLogFile = PythonDirectory + 'SS304_test_log.csv'
FatigueDirectory = 'F:\\Database\\Fatigue\\'
ArticleFigureDirectory = 'F:\\GitHub\\python\\output\\'
experiment_type_dict = {}
experiment_type_dict['Tensile']=['3000']
experiment_type_dict['TC']=[
'3001', #ramp:+-0.15
'3002', #ramp:+-0.20
'3003', #ramp:+-0.25
'3004', #ramp:+-0.30
'3005', #ramp:+-0.35
'3006', #ramp:+-0.40
'3007', #ramp:+-0.45
'3008', #ramp:+-0.50
'3009', #ramp:+-0.55
'3010', #ramp:+-0.60
'3011', #ramp:+-0.70
'3012', #ramp:+-0.80
'3013', #ramp:+-1.00
'3014', #ramp:+-1.00
'3015', #ramp:+-1.00
'3016', #ramp:+-1.00
]
experiment_type_dict['BIAXIAL']=[
'3101', #cyclic cross path	ramp:+-0.10
'3102', #cyclic cross path	ramp:+-0.20
'3103', #cyclic left triangle path	ramp:+-0.15
'3104', #cyclic left triangle path	ramp:+-0.30
'3105', #cyclic circle path	sin:+-0.10
'3106', #cyclic circle path	sin:+-0.15
'3107', #cyclic circle path	sin:+-0.30
'3108', #cyclic proportional path	ramp:+-0.80
'3109', #cyclic proportional path	ramp:+-0.40
'3110', #cyclic circle path	sin:+-0.80
'3111', #cyclic circle path	sin:+-0.40
'3112', #cyclic circle path	sin:+-0.60
]
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