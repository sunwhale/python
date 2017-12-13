# -*- coding: utf-8 -*-
"""
Created on Thu Jan 05 11:50:47 2017

@author: j.Sun
"""

import matplotlib.pyplot as plt
import shutil
from matplotlib.ticker import MultipleLocator,ScalarFormatter,FormatStrFormatter 
from Data import FatigueData,PlotData,ExperimentData
from Constants import *
from plot_format import plot_format
from Material import material_in718,material_in718_NASA,material_in718_BHU

def create_plot_data_monotonic_cyclic_osgood(figure_path=None,figure_name=None):
#==============================================================================
# x,y label
#==============================================================================
    xlabel = xylabels['axial_strain_amplitude']
    ylabel = xylabels['axial_stress_amplitude']
#==============================================================================
# plot lines
#==============================================================================
    i = 0
    marker_list = ['s','o','^','D']
    plot_data = PlotData()
    cyclic_strain = []
    cyclic_stress = []
    for name in experiment_type_dict['TC-IF']:
#        print name
        filename = ExperimentDirectory + name + '.csv'
        experiment = ExperimentData(filename)
#        print experiment.half_life_cycle
        stress = experiment.obtainNthCycle('axial_stress',experiment.half_life_cycle)
        strain = experiment.obtainNthCycle('axial_strain',experiment.half_life_cycle)
        stress_amplitude = (max(stress) - min(stress))/2.0
        strain_amplitude = (max(strain) - min(strain))/2.0*100
        cyclic_strain.append(strain_amplitude)
        cyclic_stress.append(stress_amplitude)
    plot_data.addLine(cyclic_strain,
                      cyclic_stress,
                      xlabel=xlabel,
                      ylabel=ylabel,
                      linelabel='Experiment',
                      linewidth=2,
                      linestyle='',
                      marker=marker_list[i],
                      markersize=12,
                      color='red')
    i += 1

    name = '7102' # 650C monotonic tension
    filename = ExperimentDirectory + name + '.csv'
    experiment = ExperimentData(filename)
    monotonic_stress = list(experiment.axial_stress)[::20]
    monotonic_strain = list(experiment.axial_strain*100.0)[::20]
    plot_data.addLine(monotonic_strain,
                      monotonic_stress,
                      xlabel=xlabel,
                      ylabel=ylabel,
                      linelabel='Experiment',
                      linewidth=2,
                      linestyle='',
                      marker=marker_list[i],
                      markersize=12,
                      color='blue')
    i += 1
        
    material = material_in718()
    epsilon,sigma = material.plotStrengthCoefficient()
    plot_data.addLine(epsilon*100,
                      sigma,
                      xlabel=xlabel,
                      ylabel=ylabel,
                      linelabel='Monotonic',
                      linewidth=2,
                      linestyle='-',
                      marker=None,
                      markersize=12,
                      color='blue')
                      
    epsilon,sigma = material.plotCyclicStrengthCoefficient()
    plot_data.addLine(epsilon*100,
                      sigma,
                      xlabel=xlabel,
                      ylabel=ylabel,
                      linelabel='Cyclic',
                      linewidth=2,
                      linestyle='-',
                      marker=None,
                      markersize=12,
                      color='red')
    
    plot_data.writeToFile(figure_path,figure_name)
    
def plot_monotonic_cyclic_osgood(figure_path=None,figure_name=None,save_types=[]):
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
    ax.yaxis.set_major_locator(MultipleLocator(200))
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
figure_name='plot_monotonic_cyclic_osgood'
create_plot_data_monotonic_cyclic_osgood(figure_path,figure_name)
plot_monotonic_cyclic_osgood(figure_path,figure_name,save_types=['.pdf'])

shutil.copy(__file__,ArticleFigureDirectory)