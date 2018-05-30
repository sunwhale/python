# -*- coding: utf-8 -*-
"""
Created on Wed Nov 21 14:03:59 2012

@author: Sun
"""

from abaqus import *
from abaqusConstants import *
import numpy as np
import os
from Constants_reflection import *

#===============================================================================
# OPTPUT TEMP
#===============================================================================
filename = AbaqusTempDirectory + 'OnePlaneTrn.odb'
o1 = session.openOdb(name=filename)
odb = session.odbs[filename]
frame = odb.steps['Heat'].frames[-1]
temperature_field = frame.fieldOutputs['NT11']
nset_region = odb.rootAssembly.instances['PART-1-1'].nodeSets['_HORIZONTALLINE']
temperature_field_comp = temperature_field.getSubset(region=nset_region)
nset_val = map(lambda x,y:[x.label,x.coordinates[0],x.coordinates[1],x.coordinates[2], y.data], nset_region.nodes, temperature_field_comp.values)
#===============================================================================
# WRITE
#===============================================================================
result_reaction_table_filename = SimulationDirectory + 'HORIZONTALLINE.txt'
file_disp_result = open(result_reaction_table_filename, 'w')
for val in nset_val:
    print >>file_disp_result, '%s,  %s,  %s,  %s,  %s,' % (val[0], val[1], val[2], val[3], val[4])
file_disp_result.close()

#===============================================================================
# VerticalLine
#===============================================================================
nset_region = odb.rootAssembly.instances['PART-1-1'].nodeSets['_VERTICALLINE']
temperature_field_comp = temperature_field.getSubset(region=nset_region)
nset_val = map(lambda x,y:[x.label,x.coordinates[0],x.coordinates[1],x.coordinates[2], y.data], nset_region.nodes, temperature_field_comp.values)
#===============================================================================
# WRITE
#===============================================================================
result_reaction_table_filename = SimulationDirectory + 'VERTICALLINE.txt'
file_disp_result = open(result_reaction_table_filename, 'w')
for val in nset_val:
    print >>file_disp_result, '%s,  %s,  %s,  %s,  %s,' % (val[0], val[1], val[2], val[3], val[4])
file_disp_result.close()