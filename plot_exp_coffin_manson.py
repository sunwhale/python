# -*- coding: utf-8 -*-
"""
Created on Thu Jan 05 11:50:47 2017

@author: j.Sun
"""

import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator,ScalarFormatter,FormatStrFormatter 
from Data import FatigueData
from Constants import *
from plot_format import plot_format
from Material import material_in718,material_in718_NASA,material_in718_BHU

life_NASA = [33,35,200,515,830,1072,1785,1850,2868,4323,8903,13332,13762]
strain_amplitude_NASA = [1.75,1.75,0.85,0.85,0.575,0.5,0.507,0.575,0.398,0.402,0.329,0.323,0.328] # %
stess_amplitude_NASA = [873.5,873.5,726.6,726.6,636.8,601,605,636.8,538,541,480,474,479] # MPa

life_BHU = [100000,48745,50900,35940,12980,6300,2220,678,515]
strain_amplitude_BHU = [0.4,0.42,0.43,0.45,0.47,0.5,0.6,0.8,1.0] # %
stess_amplitude_BHU = [621,638,607,606,653,646,749,756,779] # MPa

def plot_exp_coffin_manson(fatigue_data,figure_name=None,figure_path=None,save_types=[]):
#==============================================================================
# title
#==============================================================================
    title=''
#==============================================================================
# figure format
# http://matplotlib.org/users/customizing.html?highlight=rcparams
#==============================================================================
    plot_format()   
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
    plt.xlim(1E1,2E6)
    plt.ylim(0.2,2.0)
#==============================================================================
# x,y label
#==============================================================================
    plt.xlabel(r'Experimental Fatigue Lifetime $N_{\rm{f}}$')
    plt.ylabel(r'Equivalent Strain Amplitude $\Delta\varepsilon_{\rm{ep}}/2$, %')
#==============================================================================
# xy log scale
#==============================================================================
    plt.xscale('log')
    plt.yscale('log')
#==============================================================================
# xy axial equal
#==============================================================================
    ax = plt.gca()
    ax.set_aspect('equal')
    ax.set_aspect('auto')
#==============================================================================
# http://stackoverflow.com/questions/21920233/matplotlib-log-scale-tick-label-number-formatting
#==============================================================================
    ax.yaxis.set_major_locator(MultipleLocator(0.5))
    ax.yaxis.set_minor_locator(MultipleLocator(0.05))
    ax.yaxis.set_major_formatter(ScalarFormatter())
#==============================================================================
# plot lines
#==============================================================================
    i = 0
    marker_list = ['s','o','^','D','<']
    for load_type in ['TC-IF']:
        experimental_life = fatigue_data.loadTypeFilter(load_type,'experimental_life')
        equivalent_strain_amplitude = fatigue_data.loadTypeFilter(load_type,'equivalent_strain_amplitude')
        plt.plot(experimental_life, equivalent_strain_amplitude, label='Tsinghua', linewidth=2, linestyle='',
                     marker=marker_list[i], markersize=12)
        i += 1
    plt.plot(life_NASA, strain_amplitude_NASA, label='NASA', linewidth=2, linestyle='',
                     marker=marker_list[i], markersize=12)
    i += 1
    plt.plot(life_BHU, strain_amplitude_BHU, label='G.S. Mahobia 2014', linewidth=2, linestyle='',
                     marker=marker_list[i], markersize=12)
                     
    material = material_in718()
    life,epsilon_amplitude = material.plotMansonCoffinAxial()
    plt.plot(life,epsilon_amplitude*100,linewidth=1.5,color='blue')
    
    material = material_in718_NASA()
    life,epsilon_amplitude = material.plotMansonCoffinAxial()
    plt.plot(life,epsilon_amplitude*100,linewidth=1.5,color='green')
    
    material = material_in718_BHU()
    life,epsilon_amplitude = material.plotMansonCoffinAxial()
    plt.plot(life,epsilon_amplitude*100,linewidth=1.5,color='red')
#==============================================================================
# show legend
#==============================================================================
    plt.legend(loc=0)
#    plt.legend(loc=0,fontsize='small',frameon=True,numpoints=1,title='Temperature')
#==============================================================================
# save figures
#==============================================================================
    if figure_path <> None and figure_name<> None:
        for save_type in save_types:
            plt.savefig(figure_path + figure_name + save_type, dpi=150)
            print 'save as', figure_path + figure_name + save_type
    plt.show()
    plt.close()

fatigue_file = '%s%s.csv' % (FatigueDirectory,'BM')
fatigue_data = FatigueData(fatigue_file)
plot_exp_coffin_manson(fatigue_data,figure_path=ArticleFigureDirectory,figure_name='plot_exp_coffin_manson',save_types=['.pdf'])