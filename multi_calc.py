# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 13:46:18 2017

@author: SunJingyu
"""

from Constants import *
from Data import *

#for name in experiment_type_dict['BIAXIAL']:
for name in experiment_type_dict['TC-IP']+experiment_type_dict['TC-OP']+experiment_type_dict['TC-90']:
    filename = 'umat_cal_in718_%s.cal' % name
    outfile = open(filename, 'w')
    print >>outfile,"""
# -*- coding: utf-8 -*-
from Constants import *
from Data import *
from workbench_convection import workbench
name = '%s'
exp_filename = ExperimentDirectory + name + '.csv'
experiment = ExperimentData(exp_filename)
workbench(name,loading_cycles=experiment.total_axial_count,copy=True)
""" % name
    outfile.close()

outfile = open('multi_calc.bat', 'w')
#for name in experiment_type_dict['BIAXIAL']:
for name in experiment_type_dict['TC-IP']+experiment_type_dict['TC-OP']+experiment_type_dict['TC-90']:
    print >>outfile,'start python umat_cal_in718_%s.cal' % name
outfile.close()