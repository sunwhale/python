# -*- coding: utf-8 -*-
"""
Created on Wed Nov 21 14:03:59 2012

@author: Sun
"""

from abaqus import * 
from abaqusConstants import *
from odbAccess import *
import numpy as np
import os
from Constants import *
from PreprocABAQUSParameters import *

os.chdir(AbaqusWorkDirectory)
filenames = os.listdir(AbaqusWorkDirectory)

H_field_name = 'HFL'

if nonlinear == OFF:
    U_field_name = 'U'
    S_field_name = 'S'
    E_field_name = 'E'
    T_field_name = 'NT11'
    
if nonlinear == ON:
    U_field_name = 'U'
    S_field_name = 'S'
    E_field_name = 'LE'
    T_field_name = 'NT11'
    
for filename in filenames:
    if filename.split('.')[-1] == 'odb':
        odb = session.openOdb(name=filename)
        frames = odb.steps['Step-1'].frames
#        for i in range(len(odb.rootAssembly.instances['PART-1-1'].nodes)):
#            resultfile = open('node_' + str(i+1) + '.csv', 'w')
        for i in [7]:
            resultfile = open(filename.split('.')[0] + '.csv', 'w')
            header = ''
            header += 'NodeLabel,'
#            header += 'x,'
#            header += 'y,'
#            header += 'z,'
            header += 'Frame,'
            header += 'Time,'
            header += 'Temperature,'
#            header += U_field_name+'1,'
#            header += U_field_name+'2,'
#            header += (E_field_name+'11,')
            header += (E_field_name+'22,')
#            header += (E_field_name+'33,')
#            header += (E_field_name+'12,')
#            header += (E_field_name+'13,')
            header += (E_field_name+'23,')
#            header += (S_field_name+'11,')
            header += (S_field_name+'22,')
#            header += (S_field_name+'33,')
#            header += (S_field_name+'12,')
#            header += (S_field_name+'13,')
            header += (S_field_name+'23,')
            header += 'Mises,'
            header += (H_field_name+'1,')
            header += (H_field_name+'2,')
            print >>resultfile, header[:-1]
            
            for frame in frames:
                x=odb.rootAssembly.instances['PART-1-1'].nodes[i].coordinates[0]
                y=odb.rootAssembly.instances['PART-1-1'].nodes[i].coordinates[1]
                z=odb.rootAssembly.instances['PART-1-1'].nodes[i].coordinates[2]
#                U_field_comp = frame.fieldOutputs[U_field_name].getSubset(position=NODAL,region=odb.rootAssembly.instances['PART-1-1'].nodes[i])
                T_field_comp = frame.fieldOutputs[T_field_name].getSubset(position=NODAL,region=odb.rootAssembly.instances['PART-1-1'].nodes[i])                
                H_field_comp = frame.fieldOutputs[H_field_name].getSubset(position=ELEMENT_NODAL,region=odb.rootAssembly.instances['PART-1-1'].nodes[i])                
#                S_field_comp = frame.fieldOutputs[S_field_name].getSubset(position=ELEMENT_NODAL,region=odb.rootAssembly.instances['PART-1-1'].nodes[i])
#                E_field_comp = frame.fieldOutputs[E_field_name].getSubset(position=ELEMENT_NODAL,region=odb.rootAssembly.instances['PART-1-1'].nodes[i])

                nset_val_T = T_field_comp.values[0]
                
                line = ''
                line += '%s,' % (nset_val_T.nodeLabel)
#                line += '%s,' % (x)
#                line += '%s,' % (y)
#                line += '%s,' % (z)
                line += '%s,' % (frame.frameId)
                line += '%s,' % (frame.frameValue)
                line += '%s,' % (T_field_comp.values[0].data + 273.15)
#                line += '%s,' % (nset_val_U.data[0])
#                line += '%s,' % (nset_val_U.data[1])
#                line += '%s,' % (E_field_comp.getScalarField(componentLabel=E_field_name+'11').values[0].data)
                line += '%s,' % (0)
#                line += '%s,' % (E_field_comp.getScalarField(componentLabel=E_field_name+'33').values[0].data)
#                line += '%s,'  % (E_field_comp.getScalarField(componentLabel=E_field_name+'12').values[0].data)
#                line += '%s,'  % (E_field_comp.getScalarField(componentLabel=E_field_name+'13').values[0].data)
                line += '%s,'  % (0)
#                line += '%s,' % (S_field_comp.getScalarField(componentLabel=S_field_name+'11').values[0].data)
                line += '%s,' % (0)
#                line += '%s,' % (S_field_comp.getScalarField(componentLabel=S_field_name+'33').values[0].data)
#                line += '%s,' % (S_field_comp.getScalarField(componentLabel=S_field_name+'12').values[0].data)
#                line += '%s,' % (S_field_comp.getScalarField(componentLabel=S_field_name+'13').values[0].data)
                line += '%s,' % (0)
                line += '%s,' % (0)
                line += '%s,' % (H_field_comp.getScalarField(componentLabel=H_field_name+'1').values[0].data)
                line += '%s,' % (H_field_comp.getScalarField(componentLabel=H_field_name+'2').values[0].data)
                print >>resultfile, line[:-1]
            resultfile.close()