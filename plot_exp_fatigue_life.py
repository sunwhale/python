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
from Material import material_in718


def plot_exp_fatigue_life(fatigue_data,figure_name=None,figure_path=None,save_types=[]):
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
    plt.xlim(1E1,1E4)
    plt.ylim(0.4,1.1)
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
    ax.yaxis.set_major_locator(MultipleLocator(0.1))
    ax.yaxis.set_minor_locator(MultipleLocator(0.01))
    ax.yaxis.set_major_formatter(ScalarFormatter())
#==============================================================================
# print arrow
#==============================================================================
#plt.annotate('local max', xy = (100, 400), xytext = (100, 100),arrowprops = dict(facecolor = 'black', shrink = 0.1))
#plt.arrow(100,180,0,270,head_width=10)
#==============================================================================
# plot lines
#==============================================================================
    i = 0
    marker_list = ['s','o','^','D']
    for load_type in ['TC-IP','TC-OP','PRO-IP','NPR-IP']:
        experimental_life = fatigue_data.loadTypeFilter(load_type,'experimental_life')
        equivalent_strain_amplitude = fatigue_data.loadTypeFilter(load_type,'equivalent_strain_amplitude')
        plt.plot(experimental_life, equivalent_strain_amplitude, label=load_type, linewidth=2, linestyle='',
                     marker=marker_list[i], markersize=12)
        i += 1
    material = material_in718()
    life,epsilon_amplitude = material.plotMansonCoffinAxial()
    plt.plot(life,epsilon_amplitude*100,label='TC-IF',linewidth=1.5,color='black')
#==============================================================================
# show legend
#==============================================================================
    plt.legend(loc=0)
#    plt.legend(loc=0,fontsize='small',frameon=True,numpoints=1,title='Temperature')
#    plt.legend(loc=0,fontsize='small',frameon=True,numpoints=1,title='$\gamma/\sqrt{3}$')
#==============================================================================
# save figures
#==============================================================================
    if figure_path <> None and figure_name<> None:
        for save_type in save_types:
            plt.savefig(figure_path + figure_name + save_type, dpi=150)
            print 'save as', figure_path + figure_name + save_type

    plt.show()
    plt.close()

ArticleFigureDirectiory = 'F:\\Articles\\Fatigue\\Figs\\'
fatigue_file = '%s%s.csv' % (FatigueDirectiory,'BM')
fatigue_data = FatigueData(fatigue_file)
plot_exp_fatigue_life(fatigue_data,figure_path=ArticleFigureDirectiory,figure_name='Test',save_types=['.pdf'])