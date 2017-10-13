# -*- coding: utf-8 -*-
"""
Created on Thu Jan 05 11:50:47 2017

@author: j.Sun
"""

from Data import *
from Constants import *

"""
导入word，页边距1mm，分两栏，宋体7号，固定行距10.9
"""
exp_all=[
'7005',
'7006',
'7007',
'7008',
'7009',
'7010',
'7011',
'7012',
'7013',
'7014',
'7101',
'7102',
'7103',
'7211',
'7212',
'7213',
'7200',
'7018',
'7017',
'7025',
'7026',
'7028',
'7029',
'7030',
'7031',
'7032',
'7033',
'7034',
'7035',
'7036',
'7037',
'7038',
'7039',
'7040',
'7041',
'7042',
'7043',
'7044',
'7045',
'7046',
'7047',
'7048',
'7049',
'7110',
'7111',
'7112',
'7113',
'7114',
'7115',
'7116',
'7201',
'7202',
'7203',
'7204',
'7205',
'7206',
'7207',
'7208',
'7209',
'7301'
]

exp_sem=[
'7033',
'7036',
'7040',
'7046',
'7047',
'7111',
'7112',
'7206',
'7207',
'7209',
'7301',
]

experiment_log = ExperimentLog(ExperimentLogFile)

#for name in exp_all:
for name in exp_sem:
    experiment_log.output(name)
    print
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
    axial_displacement = float(experiment_log.obtainItem(name,'axial_displacement',regular)[0])
    angel_strain = float(experiment_log.obtainItem(name,'angel_strain',regular)[0])
    equivalent_strain = float(experiment_log.obtainItem(name,'equivalent_strain',regular)[0])
    period = float(experiment_log.obtainItem(name,'period',regular)[0])
    axial_temperature_phase = float(experiment_log.obtainItem(name,'axial_temperature_phase',regular)[0])
    life = float(experiment_log.obtainItem(name,'comments',regular)[0])