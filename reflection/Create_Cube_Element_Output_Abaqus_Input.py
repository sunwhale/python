# -*- coding: utf-8 -*-
"""
Created on Tue Nov 27 16:45:05 2012

@author: sunwhale
"""
import matplotlib.pyplot as plt
import numpy as np
import os
from Constants_reflection import *

#===============================================================================
# Shape function of C3D8
#===============================================================================
def N1(xi,eta,zeta):
    return 1.0/8.0*(1-xi)*(1-eta)*(1-zeta)
    
def N2(xi,eta,zeta):
    return 1.0/8.0*(1+xi)*(1-eta)*(1-zeta)
    
def N3(xi,eta,zeta):
    return 1.0/8.0*(1+xi)*(1+eta)*(1-zeta)
    
def N4(xi,eta,zeta):
    return 1.0/8.0*(1-xi)*(1+eta)*(1-zeta)
    
def N5(xi,eta,zeta):
    return 1.0/8.0*(1-xi)*(1-eta)*(1+zeta)
    
def N6(xi,eta,zeta):
    return 1.0/8.0*(1+xi)*(1-eta)*(1+zeta)
    
def N7(xi,eta,zeta):
    return 1.0/8.0*(1+xi)*(1+eta)*(1+zeta)
    
def N8(xi,eta,zeta):
    return 1.0/8.0*(1-xi)*(1+eta)*(1+zeta)
#===============================================================================
# Calculate element middle location
#===============================================================================
def Element_mid_location(node1,node2,node3,node4,node5,node6,node7,node8):
    xi=0.0
    eta=0.0
    zeta=0.0
    x=node1[0]*N1(xi,eta,zeta)+node2[0]*N2(xi,eta,zeta)+node3[0]*N3(xi,eta,zeta)+node4[0]*N4(xi,eta,zeta)+node5[0]*N5(xi,eta,zeta)+node6[0]*N6(xi,eta,zeta)+node7[0]*N7(xi,eta,zeta)+node8[0]*N8(xi,eta,zeta)
    y=node1[1]*N1(xi,eta,zeta)+node2[1]*N2(xi,eta,zeta)+node3[1]*N3(xi,eta,zeta)+node4[1]*N4(xi,eta,zeta)+node5[1]*N5(xi,eta,zeta)+node6[1]*N6(xi,eta,zeta)+node7[1]*N7(xi,eta,zeta)+node8[1]*N8(xi,eta,zeta)
    z=node1[2]*N1(xi,eta,zeta)+node2[2]*N2(xi,eta,zeta)+node3[2]*N3(xi,eta,zeta)+node4[2]*N4(xi,eta,zeta)+node5[2]*N5(xi,eta,zeta)+node6[2]*N6(xi,eta,zeta)+node7[2]*N7(xi,eta,zeta)+node8[2]*N8(xi,eta,zeta)
    mid_node=np.array([x,y,z])    
    return mid_node
#===============================================================================
# Calculate node number with i,j,k
#===============================================================================
def Node_number(i,j,k,p,q,r):
    return k*len(p)*len(q)+j*len(p)+i+1
#===============================================================================
# Calculate element number with i,j,k
#===============================================================================
def Element_number(i,j,k,P,Q,R):
    return k*P*Q+j*P+i+1
#===============================================================================
# Define cubic length,width,height
#===============================================================================
length = 100.0e-3
width = 100.0e-3
height = 2.35e-3
div_in_length = 128
div_in_width = 128
div_in_height = 4
Area = length * width / div_in_length / div_in_width
node_in_x = div_in_length + 1
node_in_y = div_in_length + 1
node_in_level = node_in_x * node_in_y
#===============================================================================
# p,q,r are nodes lists in x,y,z directions
#===============================================================================
p=[-length/2 + i*length/(div_in_length-1) for i in range(0,div_in_length+1)]
q=[-width/2 + i*width/(div_in_width-1) for i in range(0,div_in_width+1)]
r=[-height/2 + i*height/(div_in_height-1) for i in range(0,div_in_height+1)]
#===============================================================================
# P,Q,R are elements number in x,y,z directions
#===============================================================================
P = len(p)-1
Q = len(q)-1
R = len(r)-1
#===============================================================================
# nodes[0] is original reference node
#===============================================================================
nodes = []
nodes.append(np.array([0,0,0]))
#===============================================================================
# Generate node list
#===============================================================================
for rr in r:
    for qq in q:
        for pp in p:
            nodes.append(np.array([nodes[0][0] + pp ,nodes[0][1] + qq, nodes[0][2] + rr]))
#===============================================================================
# create command file and component
#===============================================================================
abaqus_input_filename = AbaqusTempDirectory + 'OnePlaneTrn.inp'
abaqus_input_file = open(abaqus_input_filename, 'w')
abaqus_input_file.writelines('*Node' + '\n')
#===============================================================================
# create nodess
#===============================================================================
node_number = 0
for node in nodes[1:]:
    node_number += 1
    abaqus_input_file.writelines('      ' + str(node_number) + ',' + str(node[0]) + ',' + str(node[1]) + ',' + str(node[2]) + '\n')
#===============================================================================
# create elememts
#===============================================================================
abaqus_input_file.writelines('*Element, type=DC3D8' + '\n')
element_number = 0
remove_element = []
for k in range(0, R):
    for j in range(0, Q):
        for i in range(0, P):
            element_number += 1
            node_number1 = Node_number(i,j,k,p,q,r)
            node_number2 = Node_number(i+1,j,k,p,q,r)
            node_number3 = Node_number(i+1,j+1,k,p,q,r)
            node_number4 = Node_number(i,j+1,k,p,q,r)
            node_number5 = Node_number(i,j,k+1,p,q,r)
            node_number6 = Node_number(i+1,j,k+1,p,q,r)
            node_number7 = Node_number(i+1,j+1,k+1,p,q,r)
            node_number8 = Node_number(i,j+1,k+1,p,q,r)
            node_list =  str(node_number1) + ',' + str(node_number2) + ',' + str(node_number3) + ',' + str(node_number4) + ',' + str(node_number5) + ',' + str(node_number6) + ',' + str(node_number7) + ',' + str(node_number8)         
            abaqus_input_file.writelines('      '+ str(element_number) + ',' + node_list + '\n')
            if (q[j]<-40e-3 or q[j]>40e-3):
                remove_element.append(element_number)
#===============================================================================
# Set
#===============================================================================
abaqus_input_file.writelines('*Elset, elset=ElementAll, generate' + '\n')
abaqus_input_file.writelines('1,' + str(element_number) + ',1' + '\n')
abaqus_input_file.writelines('*Nset, nset=NodeAll, generate' + '\n')
abaqus_input_file.writelines('1,' + str(node_number) + ',1' + '\n')
abaqus_input_file.writelines('*Nset, nset=_HorizontalLine, generate' + '\n')
abaqus_input_file.writelines(str(node_in_x * int(node_in_y/2) + 1 + node_in_level * div_in_height/2 ) + ',' + str(node_in_x * (int(node_in_y/2) + 1) + node_in_level * div_in_height/2 ) + ',1' + '\n')
abaqus_input_file.writelines('*Nset, nset=_VerticalLine, generate' + '\n')
abaqus_input_file.writelines(str(int(node_in_y/2) + 1 + node_in_level * div_in_height/2) + ',' + str(node_in_x * (node_in_y-1) + (int(node_in_y/2) + 1) + node_in_level * div_in_height/2 ) + ',' + str(node_in_x) + '\n')
#===============================================================================
# revemove element set
#===============================================================================
abaqus_input_file.writelines('*Elset, elset=ElementRemove' + '\n')
count = 0
for s in remove_element:
    count += 1
    abaqus_input_file.writelines( str(s) )
    if ( count%10==0 or count==len(remove_element) ):
        abaqus_input_file.writelines( '\n' )
    else:
        abaqus_input_file.writelines( ',' )
#===============================================================================
# section
#===============================================================================
abaqus_input_file.writelines('*Solid Section, elset=ElementAll, material="Strainless Steel"' + '\n')
#===============================================================================
# Surface
#===============================================================================
for i in range(div_in_length*div_in_width):
    abaqus_input_file.writelines('*Surface, type=ELEMENT, name=_Surface' + str(i+1) + '\n')
    abaqus_input_file.writelines(str(i+1) + ',S1' + '\n')
#===============================================================================
# convection surface
#===============================================================================
abaqus_input_file.writelines('*Elset, elset=_ConvectionS1, generate' + str(i+1) + '\n')
abaqus_input_file.writelines(str(1) + ',' + str(div_in_length*div_in_width) + ',1' + '\n')
abaqus_input_file.writelines('*Elset, elset=_ConvectionS2, generate' + str(i+1) + '\n')
abaqus_input_file.writelines(str(div_in_length*div_in_width*(div_in_height-1) + 1) + ',' + str(div_in_length*div_in_width*div_in_height) + ',1' + '\n')

abaqus_input_file.writelines('*Surface, type=ELEMENT, name=_SurfaceConvectionS1' + '\n')
abaqus_input_file.writelines('_ConvectionS1,S1' + '\n')
abaqus_input_file.writelines('*Surface, type=ELEMENT, name=_SurfaceConvectionS2' + '\n')
abaqus_input_file.writelines('_ConvectionS2,S2' + '\n')
#===============================================================================
# Material
#===============================================================================
abaqus_input_file.writelines('*Material, name="Strainless Steel"' + '\n')
abaqus_input_file.writelines('*Conductivity' + '\n')
abaqus_input_file.writelines(' 16.221, 400.' + '\n')
abaqus_input_file.writelines(' 18.158, 500.' + '\n')
abaqus_input_file.writelines('  20.09, 600.' + '\n')
abaqus_input_file.writelines('  22.03, 700.' + '\n')
abaqus_input_file.writelines('  23.97, 800.' + '\n')
abaqus_input_file.writelines('   25.9, 900.' + '\n')
abaqus_input_file.writelines('  26.87, 950.' + '\n')
abaqus_input_file.writelines(' 27.845,1000.' + '\n')
abaqus_input_file.writelines('  28.81,1050.' + '\n')
abaqus_input_file.writelines('  29.78,1100.' + '\n')
abaqus_input_file.writelines('*Density' + '\n')
abaqus_input_file.writelines('8400.,' + '\n')
abaqus_input_file.writelines('*Specific Heat' + '\n')
abaqus_input_file.writelines('490.,' + '\n')
abaqus_input_file.writelines('*Material, name=Wolfram' + '\n')
abaqus_input_file.writelines('*Conductivity' + '\n')
abaqus_input_file.writelines('170.,' + '\n')
abaqus_input_file.writelines('*Density' + '\n')
abaqus_input_file.writelines('19300.,' + '\n')
abaqus_input_file.writelines('*Specific Heat' + '\n')
abaqus_input_file.writelines('138.,' + '\n')
abaqus_input_file.writelines('*Physical Constants, absolute zero=-273.15, stefan boltzmann=5.67e-08' + '\n')
abaqus_input_file.writelines('*Initial Conditions, type=TEMPERATURE' + '\n')
abaqus_input_file.writelines('NodeAll, 25.' + '\n')
#===============================================================================
# Step
#===============================================================================
abaqus_input_file.writelines('*Step, name=Heat, inc=1000' + '\n')
abaqus_input_file.writelines('*Heat Transfer, end=PERIOD, deltmx=50.' + '\n')
abaqus_input_file.writelines('0.1, 200., 1e-07, 20., ' + '\n')
#===============================================================================
# Add boundary conditions
#===============================================================================
TracePro_file = TraceProDirectory + 'Heat128Pixels_1mm_1mm.txt'
emissivity = 0.95
efficiency = 0.95
list1 = open(TracePro_file,"r").readlines()	
list2 = [i.strip() for i in list1]
list3 = [j.split('\t') for j in list2]
list4 = [[float(j) for j in i] for i in list3]
for i in list4:
    for j in i:
        s = s + j
HeatFlux = [[j * 400 * 0.94 * efficiency * 0.277 * emissivity / s / Area  for j in i] for i in list4]
s = 0
for i in HeatFlux:
    for j in i:
        s = s + j*Area

print 'Total Heat Power is ' , s
for i in range(div_in_length):
    for j in range(div_in_length):
        abaqus_input_file.writelines('*Dsflux' + '\n')
        abaqus_input_file.writelines('_Surface' + str(Element_number(i,j,0,P,Q,R)) + ',S,' + str(HeatFlux[i][j]) + '\n')
#===============================================================================
# Remove element
#===============================================================================
abaqus_input_file.writelines('*Model Change, remove' + '\n')
abaqus_input_file.writelines('ElementRemove, ' + '\n')
#===============================================================================
# Convection        
#===============================================================================
abaqus_input_file.writelines('*Sfilm' + '\n')
abaqus_input_file.writelines('_SurfaceConvectionS1, F, 25., 5.5' + '\n')
abaqus_input_file.writelines('*Sfilm' + '\n')
abaqus_input_file.writelines('_SurfaceConvectionS2, F, 25., 5.5' + '\n')
#===============================================================================
# Radiation
#===============================================================================
abaqus_input_file.writelines('*Sradiate' + '\n')
abaqus_input_file.writelines('_SurfaceConvectionS1, R, 25., 0.95' + '\n')
abaqus_input_file.writelines('*Sradiate' + '\n')
abaqus_input_file.writelines('_SurfaceConvectionS2, R, 25., 0.95' + '\n')
#===============================================================================
# Output
#===============================================================================
abaqus_input_file.writelines('*Restart, write, frequency=0' + '\n')
abaqus_input_file.writelines('*Output, field, frequency=1' + '\n')
abaqus_input_file.writelines('*Node Output' + '\n')
abaqus_input_file.writelines('NT, ' + '\n')
abaqus_input_file.writelines('*Element Output, directions=YES' + '\n')
abaqus_input_file.writelines('HFL, TEMP' + '\n')
abaqus_input_file.writelines('*Output, history, frequency=0' + '\n')
abaqus_input_file.writelines('*End Step' + '\n')
#===============================================================================
# copy input file to work directionary
#===============================================================================
abaqus_input_file.close()

os.chdir(AbaqusTempDirectory)
cmd1 = 'del OnePlaneTrn.odb'
os.system(cmd1)
cmd3 = 'abaqus job=OnePlaneTrn cpus=4 int'
os.system(cmd3)
cmd4 = 'abaqus viewer noGUI="%sPostproc_Output_TEMP_and_HFL.py"' % (PythonDirectory)
os.system(cmd4)