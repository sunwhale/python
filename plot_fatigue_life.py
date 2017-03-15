# -*- coding: utf-8 -*-
"""
Created on Thu Jan 05 11:50:47 2017

@author: j.Sun
"""

import numpy as np
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
from Data import FatigueData

def plot_fatigue_life(fatigue_data,figure_name=None,figure_path=None):
#==============================================================================
# title
#==============================================================================
    title=''
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
    marker_list = ['s','o','^','D']
    for load_type in ['TC-IP','TC-OP','PRO-IP','NPR-IP']:
        experimental_life,predicted_life = fatigue_data.loadTypeFilter(load_type)
        plt.plot(experimental_life, predicted_life, label=load_type, linewidth=2, linestyle='',
                     marker=marker_list[i], markersize=12)
        i += 1
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
#    plt.legend(loc=0,fontsize='small',frameon=True,numpoints=1,title='Temperature')
#    plt.legend(loc=0,fontsize='small',frameon=True,numpoints=1,title='$\gamma/\sqrt{3}$')
#==============================================================================
# save figures
#==============================================================================
    if figure_path <> None and figure_name<> None:
        save_types = ['.png','.pdf']
        for save_type in save_types:
            plt.savefig(figure_path + figure_name + save_type, dpi=150)
            print 'save as', figure_path + figure_name + save_type

#    plt.show()
    plt.close()

FatigueDirectiory = 'F:\\Database\\Fatigue\\'
fatigue_model_list = ['BM','FS','SWT','Liu1','Liu2','Chu']

for fatigue_model in fatigue_model_list:
    fatigue_file = '%s%s\\%s.csv' % (FatigueDirectiory,fatigue_model,fatigue_model)
    fatigue_data = FatigueData(fatigue_file)
    plot_fatigue_life(fatigue_data,figure_path=FatigueDirectiory,figure_name=fatigue_model)
    
