# -*- coding: utf-8 -*-
"""
Created on Tue Nov 27 16:45:05 2012

@author: sunwhale
"""
import matplotlib.pyplot as plt
import numpy as np
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
# Define cubic length,width,height
#===============================================================================
length = 100.0
width = 80.0
height = 2.0
#===============================================================================
# p,q,r are nodes lists in x,y,z directions
#===============================================================================
p=[-length/2 + i*length/255.0 for i in range(0,10)]
q=[-width/2 + i*width/255.0 for i in range(0,10)]
r=[-height/2 + i*height/255.0 for i in range(0,2)]
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
filename = 'nodes.cmf'
outfile = open(filename, 'w')
outfile.writelines('*collectorcreateonly(components,"solid","",4)' + '\n')
#===============================================================================
# create nodess
#===============================================================================
for node in nodes[1:]:
    outfile.writelines('*createnode(' + str(node[0])  +',' + str(node[1]) + ',' + str(node[2]) + ',0,0,0)' + '\n')
#===============================================================================
# create elememts
#===============================================================================
for k in range(0, R):
    for j in range(0, Q):
        for i in range(0, P):
            node_number1 = Node_number(i,j,k,p,q,r)
            node_number2 = Node_number(i+1,j,k,p,q,r)
            node_number3 = Node_number(i+1,j+1,k,p,q,r)
            node_number4 = Node_number(i,j+1,k,p,q,r)
            node_number5 = Node_number(i,j,k+1,p,q,r)
            node_number6 = Node_number(i+1,j,k+1,p,q,r)
            node_number7 = Node_number(i+1,j+1,k+1,p,q,r)
            node_number8 = Node_number(i,j+1,k+1,p,q,r)
            node_list =  str(node_number1) + ' ' + str(node_number2) + ' ' + str(node_number3) + ' ' + str(node_number4) + ' ' + str(node_number5) + ' ' + str(node_number6) + ' ' + str(node_number7) + ' ' + str(node_number8)         
            outfile.writelines('*createlist(node,1)' + node_list + '\n')            
            outfile.writelines('*createelement(208,5,1,1)' + '\n')
#===============================================================================
# create component and nodes for creating lines
#===============================================================================
outfile.writelines('*nodecleartempmark()' + '\n')

#outfile.writelines('*morphstoremorphvolumes(3)' + '\n')
#outfile.writelines('*createmark(nodes,1) "all"' + '\n')
#outfile.writelines('*createmark(nodes,2)' + '\n')
#outfile.writelines('*createlist(lines,1) 4' + '\n')
#outfile.writelines('*createlist(nodes,1)' + '\n')
#outfile.writelines('*createlist(lines,2) 1 2 3' + '\n')
#outfile.writelines('*createlist(nodes,2)' + '\n')
#outfile.writelines('*createplane(1,1.0000,0.0000,0.0000,0.0000,0.0000,0.0000)' + '\n')
#outfile.writelines('*morphmapdifference(nodes,1,nodes,2,1,1,2,2,1,1,2,0,1,0,1,1)' + '\n')

outfile.writelines('*writefile("output.hm",1)' + '\n')

bat_outfile = open('CREAT_HM.bat', 'w')
bat_outfile.writelines('call hmbatch -cnodes.cmf \n')