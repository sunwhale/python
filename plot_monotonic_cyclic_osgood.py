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
    plt.xlim(0,1.5)
    plt.ylim(0,1200.0)
#==============================================================================
# x,y label
#==============================================================================
    plt.xlabel(xylabels['axial_strain_amplitude'])
    plt.ylabel(xylabels['axial_stress_amplitude'])
#==============================================================================
# xy log scale
#==============================================================================
#    plt.xscale('log')
#    plt.yscale('log')
#==============================================================================
# xy axial equal
#==============================================================================
    ax = plt.gca()
    ax.set_aspect('equal')
    ax.set_aspect('auto')
#==============================================================================
# http://stackoverflow.com/questions/21920233/matplotlib-log-scale-tick-label-number-formatting
#==============================================================================
#    ax.yaxis.set_major_locator(MultipleLocator(0.5))
#    ax.yaxis.set_minor_locator(MultipleLocator(0.05))
#    ax.yaxis.set_major_formatter(ScalarFormatter())
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
                     
    material = material_in718()
    epsilon,sigma = material.plotStrengthCoefficient()
    plt.plot(epsilon*100,sigma,linewidth=1.5,color='blue')

    epsilon,sigma = material.plotCyclicStrengthCoefficient()
    plt.plot(epsilon*100,sigma,linewidth=1.5,color='red')
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

fatigue_file = '%s%s.csv' % (FatigueDirectiory,'BM')
fatigue_data = FatigueData(fatigue_file)
plot_exp_coffin_manson(fatigue_data,figure_path=ArticleFigureDirectiory,figure_name='plot_exp_coffin_manson',save_types=['.pdf'])