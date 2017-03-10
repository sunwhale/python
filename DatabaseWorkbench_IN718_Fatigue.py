# -*- coding: utf-8 -*-
"""
Workbench
"""
#==============================================================================
# Import lib
#==============================================================================
import sys
import numpy as np
import os
import xlsxwriter
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from subprocess import Popen, PIPE
from shutil import copy

#==============================================================================
# Our Experiments
#==============================================================================
sigma_f=1034.0
b=-0.04486
epsilon_f=0.11499
c=-0.52436
#==============================================================================

TEMP=650.0
EMOD=206308.7426+(-51.20306)*TEMP+0.01109*TEMP**2+(-3.84391E-05)*TEMP**3
ENU=2.901300E-01+(1.457750E-05)*TEMP+(-2.067420E-07)*TEMP**2+(2.780300E-10)*TEMP**3
EG=0.5*EMOD/(1.0+ENU)
print "EMOD: ",EMOD,"MPa"
print "ENU: ",ENU
print "EG: ",EG,"MPa"

tau_f=sigma_f/np.sqrt(3)
b0=b*1.0
gamma_f=epsilon_f*np.sqrt(3)
c0=c
shear_modulus=EG
k=1.0
S=0.36
sigma_yield=1064.0

WorkbenchDirectory = 'F:\\Database\\'
#==============================================================================
# IO
#==============================================================================
InputDirectiory = WorkbenchDirectory + 'IN718\\Fatigue\\FS'
#InputDirectiory = WorkbenchDirectory + 'IN718\\Fatigue\\BM'
#InputDirectiory = WorkbenchDirectory + 'IN718\\Fatigue\\SWT'
#InputDirectiory = WorkbenchDirectory + 'IN718\\Fatigue\\Liu'
#InputDirectiory = WorkbenchDirectory + 'IN718\\Fatigue\\Liu2'
#InputDirectiory = WorkbenchDirectory + 'IN718\\Fatigue\\Chu'
OutputDirectiory= InputDirectiory
if not os.path.isdir(OutputDirectiory):
    os.makedirs(OutputDirectiory)
    print 'Create new directory:',OutputDirectiory
#==============================================================================
# Read
#==============================================================================
#==============================================================================
# TMF tests
#==============================================================================
#sigma_f=1651
#b=-0.1017
#epsilon_f=0.83332
#c=-0.9123
#==============================================================================
# http://www.sciencedirect.com/science/article/pii/S0142112313002247
# Effect of hot corrosion on low cycle fatigue behavior of superalloy IN718 650C
#==============================================================================
#sigma_f=929.4
#b=-0.03917
#epsilon_f=0.21357
#c=-0.53731

#TEMP=650
#EMOD=206308.7426+(-51.20306)*TEMP+0.01109*TEMP**2+(-3.84391E-05)*TEMP**3
#ENU=2.901300E-01+(1.457750E-05)*TEMP+(-2.067420E-07)*TEMP**2+(2.780300E-10)*TEMP**3
#EG=0.5*EMOD/(1.0+ENU)
#print "EMOD: ",EMOD,"MPa"
#print "ENU: ",ENU
#print "EG: ",EG,"MPa"
#
#tau_f=sigma_f/np.sqrt(3)
#b0=b*1.0
#gamma_f=epsilon_f*np.sqrt(3)
#c0=c
#shear_modulus=EG
#k=1.0
#sigma_yield=1064.0
    
xyarray = []

os.chdir(InputDirectiory)
filenames = ['IP Diamond.csv','IP Proportional.csv','IP Uniaxial.csv','OP Uniaxial.csv']
workbook = xlsxwriter.Workbook('Excel' + '.xlsx')
for filename in filenames:
    print 'Reading:',filename
    infile = open(filename,'r')
    list1 = infile.readlines()
    list2 = [i.strip() for i in list1]
    list3 = [i.split(',') for i in list2]
    for i in list3:
        for j in range(len(i)):
            if i[j]=='--' or i[j]=='':
                i[j]='0'
    list4 = [[float(i[0]),float(i[1]),float(i[2]),float(i[3]),float(i[4]),float(i[5]),float(i[6]),float(i[7]),float(i[8]),float(i[9]),float(i[10]),float(i[11]),float(i[12]),float(i[13])] for i in list3[2:]]
    cycles_of_failure_experiment_list= [i[0] for i in list4]
    equvialent_strain_amplitude_list = [i[1] for i in list4]
    equvialent_stress_amplitude_list = [i[2] for i in list4]
    specimen_number_list = [i[3] for i in list4]
    critical_plane_list = [i[4] for i in list4]
    sigma_nmax_list = [i[5] for i in list4]
    delta_sigma_list = [i[6] for i in list4]
    delta_epsilon_list = [i[7] for i in list4]
    tau_nmax_list = [i[8] for i in list4]
    delta_tau_list = [i[9] for i in list4]
    delta_gamma_list = [i[10] for i in list4]
    cycles_of_failure_predicted_list_list = [i[11] for i in list4]
    fatigue_coefficient_list = [i[12] for i in list4]
    temperature_list = [i[13] for i in list4]
    
    fatemi_socie_list  = [i[12] for i in list4]
    brown_miller_list  = [i[12] for i in list4]
    smith_watson_topper_list = [i[12] for i in list4]
    liu_list = [i[12] for i in list4]
#==============================================================================
# Write excel
#==============================================================================
    worksheet = workbook.add_worksheet(filename[:-4])
    bold = workbook.add_format({'bold': 1})
    
    row_number = 1
    
#    headers = ['Number of Cycles to Failure N\-(f)','Mises Equivalent Strain Amplitude \i(\g(De))\-(eq)/2','Stress Amplitude e \i(\g(Ds))/2','Specimen','Critical Plane','sigma_n_max','delta_sigma','delta_epsilon','tau_n_max','delta_tau','delta_gamma','Predicted Fatigue Lifetime N\-(p)','Fatigue Coefficient','Temperature']
#    units = ['cycles','mm/mm','Mpa','-','deg','Mpa','Mpa','mm/mm','Mpa','Mpa','mm/mm','cycles','-','C']
    headers = list3[0]
    units = list3[1]
    comments = [ filename[:-4] for i in range(len(list3[0]))]

    # Add the worksheet data that the charts will refer to.
    worksheet.write_row('A'+str(row_number), headers, bold); row_number += 1
    worksheet.write_row('A'+str(row_number), units, bold); row_number += 1
    worksheet.write_row('A'+str(row_number), comments, bold); row_number += 1
    
    for l in list3[2:]:
        worksheet.write_row('A'+str(row_number), l); row_number += 1
#==============================================================================
# FS
#==============================================================================
    if InputDirectiory.split('\\')[-1] == 'FS':
        for i in range(len(cycles_of_failure_experiment_list)):
            fatemi_socie_list[i] = delta_gamma_list[i] * (1.0+k*sigma_nmax_list[i]/sigma_yield)
        
            def f(x):
                x0 = float(x[0])
                return [
                    tau_f/shear_modulus*(2*x0)**b0 + gamma_f*(2*x0)**c0-fatemi_socie_list[i]
                ]
            
            result = fsolve(f, [1])
            cycles_of_failure_predicted_list_list[i] = result
        
        xyarray.append([cycles_of_failure_experiment_list,cycles_of_failure_predicted_list_list,filename])
#==============================================================================
# BM
#==============================================================================
    if InputDirectiory.split('\\')[-1] == 'BM':
        for i in range(len(cycles_of_failure_experiment_list)):
            S = 0.36
            
            brown_miller_list[i] = delta_gamma_list[i] + S*delta_epsilon_list[i]
            sigma_n_mean_max = sigma_nmax_list[i] - delta_sigma_list[i]/2.0
            
            def f(x):
                x0 = float(x[0])
                return [
                    (1.3+0.7*S)*(sigma_f-2.0*sigma_n_mean_max)/EMOD*(2*x0)**b + (1.5+0.5*S)*epsilon_f*(2*x0)**c-brown_miller_list[i]
                ]
            
            result = fsolve(f, [1])
            cycles_of_failure_predicted_list_list[i] = result
      
        xyarray.append([cycles_of_failure_experiment_list,cycles_of_failure_predicted_list_list,filename])
#==============================================================================
# SWT
#==============================================================================
    if InputDirectiory.split('\\')[-1] == 'SWT':
        for i in range(len(cycles_of_failure_experiment_list)):
            smith_watson_topper_list[i] = sigma_nmax_list[i]*delta_epsilon_list[i]/2.0
    
            def f(x):
                x0 = float(x[0])
                return [
                    sigma_f*sigma_f/EMOD*(2*x0)**(2*b) + sigma_f*epsilon_f*(2*x0)**(b+c)-smith_watson_topper_list[i]
                ]
            
            result = fsolve(f, [1])
            cycles_of_failure_predicted_list_list[i] = result
        
        xyarray.append([cycles_of_failure_experiment_list,cycles_of_failure_predicted_list_list,filename])
#==============================================================================
# Liu
#==============================================================================
    if InputDirectiory.split('\\')[-1] == 'Liu':
        for i in range(len(cycles_of_failure_experiment_list)):
            delta_w1_max = delta_sigma_list[i]*delta_epsilon_list[i] + delta_tau_list[i]*delta_gamma_list[i]*2.0
            stress_ratio = (sigma_nmax_list[i]-delta_sigma_list[i])/sigma_nmax_list[i]
            delta_w1_max *= 2.0/(1.0-stress_ratio)
            
            def f(x):
                x0 = float(x[0])
                return [
                    4*sigma_f*sigma_f/EMOD*(2*x0)**(2*b) + 4*sigma_f*epsilon_f*(2*x0)**(b+c)-delta_w1_max
                ]
            
            result = fsolve(f, [1])
            cycles_of_failure_predicted_list_list[i] = result
        
        xyarray.append([cycles_of_failure_experiment_list,cycles_of_failure_predicted_list_list,filename])  
#==============================================================================
# Liu2
#==============================================================================
    if InputDirectiory.split('\\')[-1] == 'Liu2':
        for i in range(len(cycles_of_failure_experiment_list)):
            delta_w2_max = delta_sigma_list[i]*delta_epsilon_list[i] + delta_tau_list[i]*delta_gamma_list[i]*2.0
            sigma_n_mean_max = sigma_nmax_list[i] - delta_sigma_list[i]/2.0
            delta_w2_max *= sigma_f/(sigma_f-sigma_n_mean_max)

            def f(x):
                x0 = float(x[0])
                return [
                    4*tau_f*tau_f/EG*(2*x0)**(2*b0) + 4*tau_f*gamma_f*(2*x0)**(b0+c0)-delta_w2_max
                ]
            
            result = fsolve(f, [1])
            cycles_of_failure_predicted_list_list[i] = result
        
        xyarray.append([cycles_of_failure_experiment_list,cycles_of_failure_predicted_list_list,filename])  
#==============================================================================
# Chu
#==============================================================================
    if InputDirectiory.split('\\')[-1] == 'Chu':
        for i in range(len(cycles_of_failure_experiment_list)):
            delta_w_max = sigma_nmax_list[i]*delta_epsilon_list[i]/2.0 + tau_nmax_list[i]*delta_gamma_list[i]
            
            def f(x):
                x0 = float(x[0])
                return [
                    1.02*sigma_f*sigma_f/EMOD*(2*x0)**(2*b) + 1.04*sigma_f*epsilon_f*(2*x0)**(b+c)-delta_w_max
                ]
            
            result = fsolve(f, [1])
            cycles_of_failure_predicted_list_list[i] = result
        
        xyarray.append([cycles_of_failure_experiment_list,cycles_of_failure_predicted_list_list,filename])  

workbook.close()
#==============================================================================
# title
#==============================================================================
title=''
figurename='TMF 90 +-1% 300C-650C'
#xitem = 'axial_stress'
#yitem = 'shear_stress'
##labels = ['+-0.6%','+-0.726%','+-1.0%']
labels = []
save = False
save = True

#==============================================================================
# figure size
# http://matplotlib.org/users/customizing.html?highlight=rcparams
#==============================================================================
mpl.rcParams['axes.linewidth'] = 1.5

mpl.rcParams['xtick.major.size'] = 8
mpl.rcParams['xtick.major.width'] = 1.5
mpl.rcParams['xtick.minor.size'] = 4
mpl.rcParams['xtick.minor.width'] = 1.5
mpl.rcParams['xtick.labelsize'] = 18

mpl.rcParams['ytick.major.size'] = 8
mpl.rcParams['ytick.major.width'] = 1.5
mpl.rcParams['ytick.minor.size'] = 4
mpl.rcParams['ytick.minor.width'] = 1.5
mpl.rcParams['ytick.labelsize'] = 18

mpl.rcParams['figure.figsize'] = (8,6)
mpl.rcParams['axes.labelsize'] = 22  # fontsize of the x any y labels
mpl.rcParams['figure.subplot.left'] = 0.125
mpl.rcParams['figure.subplot.right'] = 0.95
mpl.rcParams['figure.subplot.bottom'] = 0.125
mpl.rcParams['figure.subplot.top'] = 0.95

mpl.rcParams['legend.fontsize'] = 18
mpl.rcParams['legend.frameon'] = False
mpl.rcParams['legend.numpoints'] = 1
mpl.rcParams['legend.scatterpoints'] = 1

fig, ax = plt.subplots()
#==============================================================================
# grid set up
#==============================================================================
#plt.grid(True, which='major',linestyle='-')
#plt.grid(True, which='minor',linestyle='-')
#plt.grid(True, which='major')
#plt.grid(True, which='minor')
#==============================================================================
# print title
#==============================================================================
plt.title(title, fontsize=16)
#==============================================================================
# x,y limite
#==============================================================================
plt.xlim(1E1,1E4)
plt.ylim(1E1,1E4)

plt.xlim(10,10000)
plt.ylim(10,10000)

#==============================================================================
# x,y label
#==============================================================================
plt.xlabel(r'Experimental Fatigue Lifetime $N_{\rm{f}}$')
plt.ylabel(r'Predicted Fatigue Lifetime $N_{\rm{p}}$')
#==============================================================================
# xy axial equal
#==============================================================================
plt.gca().set_aspect('equal')
plt.gca().set_aspect('auto')
#==============================================================================
# xy log scale
#==============================================================================
plt.xscale('log')
plt.yscale('log')
#==============================================================================
# print arrow
#==============================================================================
#plt.annotate('local max', xy = (100, 400), xytext = (100, 100),arrowprops = dict(facecolor = 'black', shrink = 0.1))
#plt.arrow(100,180,0,270,head_width=10)
#==============================================================================
# plot lines
#==============================================================================
i = 0
marker_list = ['D','^','s','o']
labels = ['NPR-IP','PRO-IP','TC-IP','TC-OP']
for xy in xyarray:
#==============================================================================
# plot with filename
#==============================================================================
    if labels == []:
        plt.plot(xy[0], xy[1], label=xy[2], linewidth=2, linestyle='',
                 marker=marker_list[i], markersize=12)
#    if labels == []:
#        plt.plot(xy[0], xy[1], label=xy[2], linewidth=2, linestyle='', \
#                 marker=marker_list[i], markersize=9, markerfacecolor='black')
#==============================================================================
# plot with given label
#==============================================================================
    else:
        plt.plot(xy[0], xy[1], label=labels[i], linewidth=2, linestyle='',
                 marker=marker_list[i], markersize=12)
#        plt.plot(xy[0], xy[1], label=labels[i], linewidth=2, linestyle='', \
#                 marker='', markersize=5, markerfacecolor='none')
    i = i + 1
#==============================================================================
# plot 2x lines
#==============================================================================
linewidth = 0.5
plt.plot([20,1e4],[10,5e3],color='black',linewidth=linewidth)
plt.plot([10,5e3],[20,1e4],color='black',linewidth=linewidth)
#==============================================================================
# plot 5x lines
#==============================================================================
plt.plot([50,1e4],[10,2e3],color='black',linewidth=linewidth)
plt.plot([10,2e3],[50,1e4],color='black',linewidth=linewidth)
#==============================================================================
# show legend
#==============================================================================
plt.legend(loc=2)
#plt.legend(loc=0,fontsize='small',frameon=True,numpoints=1,title='Temperature')
#plt.legend(loc=0,fontsize='small',frameon=True,numpoints=1,title='$\gamma/\sqrt{3}$')
#==============================================================================
# save figures
#==============================================================================
if save == True:
    plt.savefig(OutputDirectiory + figurename + '.png', dpi=150)
    print 'save as',OutputDirectiory + figurename + '.png'
    plt.savefig(OutputDirectiory + figurename + '.eps', dpi=150)
    print 'save as',OutputDirectiory + figurename + '.eps'
    plt.savefig(OutputDirectiory + figurename + '.pdf', dpi=150)
    print 'save as',OutputDirectiory + figurename + '.pdf'
#==============================================================================
# show
#==============================================================================
plt.show()
#==============================================================================
# end
#==============================================================================
print 'Finished'