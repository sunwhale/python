# -*- coding: utf-8 -*-
"""
Created on Thu Jan 05 11:50:47 2017

@author: j.Sun
"""

import matplotlib.pyplot as plt
import shutil
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator,ScalarFormatter,FormatStrFormatter 
from Data import FatigueData,PlotData,ExperimentData
from Constants import *
from plot_format import plot_format
from Material import material_in718,material_in718_NASA,material_in718_BHU

def create_plot_data_exp_half_life_cycle(figure_path=None,figure_name=None):
#==============================================================================
# x,y label
#==============================================================================
    xlabel = xylabels['axial_strain']
    ylabel = xylabels['axial_stress']
#==============================================================================
# plot lines
#==============================================================================
    i = 0
    marker_list = ['s','o','^','D']
    plot_data = PlotData()
    for name in ['7116','7111','7115','7114']:
        filename = ExperimentDirectory + name + '.csv'
        experiment = ExperimentData(filename)
        strain = experiment.obtainNthCycle('axial_strain',experiment.half_life_cycle)
        stress = experiment.obtainNthCycle('axial_stress',experiment.half_life_cycle)
        plot_data.addLine(strain*100,
                          stress,
                          xlabel=xlabel,
                          ylabel=ylabel,
                          linelabel='',
                          linewidth=2,
                          linestyle='-',
                          marker=None,
                          markersize=12,
                          color='auto')
        i += 1
    
    material = material_in718()
    epsilon,sigma = material.plotCyclicStrengthCoefficient()
    plot_data.addLine(epsilon*100,
                      sigma,
                      xlabel=xlabel,
                      ylabel=ylabel,
                      linelabel='Ramberg-Osgood',
                      linewidth=4,
                      linestyle='--',
                      marker=None,
                      markersize=12,
                      color='black')
    
    plot_data.writeToFile(figure_path,figure_name)
    
def plot_exp_half_life_cycle(figure_path=None,figure_name=None,save_types=[]):
#==============================================================================
# title
#==============================================================================
    title=''
#==============================================================================
# figure format
# http://matplotlib.org/users/customizing.html?highlight=rcparams
#==============================================================================
    plot_format()
    mpl.rcParams['figure.subplot.left'] = 0.15
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
    plt.xlim(-1,1)
    plt.ylim(-1000,1000)
#==============================================================================
# xy log scale
#==============================================================================
#    plt.xscale('log')
#    plt.yscale('log')
#==============================================================================
# xy axial equal
#==============================================================================
    ax = plt.gca()
#    ax.set_aspect('equal')
    ax.set_aspect('auto')
#==============================================================================
# plot lines
#==============================================================================
    plot_data = PlotData()
    plot_data.readFromFile(figure_path,figure_name)
    plot_data.plot()
#==============================================================================
# http://stackoverflow.com/questions/21920233/matplotlib-log-scale-tick-label-number-formatting
#==============================================================================
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.yaxis.set_major_locator(MultipleLocator(500))
    ax.yaxis.set_minor_locator(MultipleLocator(100))
    ax.yaxis.set_major_formatter(ScalarFormatter())
#==============================================================================
# show legend
#==============================================================================
    plt.legend(loc=0)
#==============================================================================
# save figures
#==============================================================================
    if figure_path <> None and figure_name<> None:
        for save_type in save_types:
            plt.savefig(figure_path + figure_name + save_type, dpi=150, transparent=True)
            print 'save as', figure_path + figure_name + save_type
    plt.show()
    plt.close()

figure_path=ArticleFigureDirectory
figure_name='plot_exp_half_life_cycle-TC-IF'
create_plot_data_exp_half_life_cycle(figure_path,figure_name)
plot_exp_half_life_cycle(figure_path,figure_name,save_types=['.pdf'])

shutil.copy(__file__,ArticleFigureDirectory)