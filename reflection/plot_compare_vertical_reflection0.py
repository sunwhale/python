# -*- coding: utf-8 -*-
"""
Created on Tue Dec 18 13:49:06 2012

@author: Sun
"""

import matplotlib.pyplot as plt

#===============================================================================
# node position
#===============================================================================
node = 16*[0]
node[0] = [0,0]
node[1] = [25e-3,0]
node[2] = [15e-3,0]
node[3] = [10e-3,0]
node[4] = [5e-3,0]
node[5] = [0e-3,0]
node[6] = [-5e-3,0]
node[7] = [-10e-3,0]
node[8] = [-15e-3,0]
node[9] = [-25e-3,0]
node[10] = [0,10e-3]
node[11] = [0,20e-3]
node[12] = [0,30e-3]
node[13] = [0,-10e-3]
node[14] = [0,-20e-3]
node[15] = [0,-30e-3]

node_x = [ n[0] for n in node]
node_y = [ n[1] for n in node]

#begin_time = '11:16:58'
begin_time = '19:29:43'
num_node = 10
txt="VERTICALLINE.txt"
Liste1 = open(txt,"r").readlines()	
Liste2 = [i.strip() for i in Liste1]
Liste3 = [j.split(',') for j in Liste2]
y_coordinate = [ float(k[2])  for k in Liste3 ]
temp_sim = [ float(k[4])  for k in Liste3 ]


#txt="..\\Experiment\\ExperimentalData\\7\\2013_08_21_B11h16m58s_E11h18m28s_800W_NI9213.txt"
txt="..\\Experiment\\ExperimentalData\\4\\2013_08_20_B19h29m43s_E19h34m13s_400W_NI9213.txt"
Liste1 = open(txt,"r").readlines()	
Liste2 = [i.strip() for i in Liste1 if len(i)>50]
Liste3 = [j.split('\t') for j in Liste2]
system_time =  [ k[0]  for k in Liste3 ]
num = [ k[0]  for k in Liste3 ].index(begin_time)
time_exp = [ float(k[1])  for k in Liste3 ]
time_shift = time_exp[num]
time_exp_shift = [ t - time_shift  for t in time_exp ]
temp_exp = [ float(k)  for k in Liste3[num+200][2:] ]

#===============================================================================
# Plot in x-y plane
#===============================================================================
title=''
xlabel='x'
ylabel='y'
plt.figure(figsize=(8,6))
plt.title(title, fontsize=16)
plt.xlabel(xlabel)
plt.ylabel(ylabel)
#plt.xlim(0,1000)
#plt.ylim(0,200)

plt.grid(True)
#plt.gca().set_aspect('equal')

for i in range(len(y_coordinate)):
    print y_coordinate[i],temp_sim[i]
for i in range(len(node_y[10:16])):
    print node_y[10:16][i],temp_exp[10:16][i]
    
plt.plot(y_coordinate, temp_sim, label='ABAQUS',linewidth=3,ls='',marker='s')
plt.plot(node_y[5], temp_exp[5], label='TK_MID',linewidth=3,ls='',marker='s')
plt.plot(node_y[10:16], temp_exp[10:16], label='TK',linewidth=3,ls='',marker='s') 
plt.legend(loc=0)
plt.show()