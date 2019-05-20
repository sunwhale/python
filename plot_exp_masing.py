# -*- coding: utf-8 -*-
"""
Created on Thu Jan 05 11:50:47 2017

@author: j.Sun
"""

import matplotlib.pyplot as plt
import shutil
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator,ScalarFormatter,FormatStrFormatter 
from Data import FatigueData,PlotData,ExperimentData,ExperimentLog
from Constants import *
from plot_format import plot_format
from Material import material_in718,material_in718_NASA,material_in718_BHU
from Functions import obtain_masing_curve, obtain_youngs_modulus

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
    experiment_log = ExperimentLog(ExperimentLogFile)
#    for name in ['7110','7111','7112','7113','7114','7115','7116']:
    for name in ['7045']:
        experiment_log.output(name)
        regular = r'.*'
        load_type = experiment_log.obtainItem(name,'load_type',regular)[0]
        regular = r'\d+\.?\d*'
        temperature_mode = experiment_log.obtainItem(name,'temperature_mode',regular)
        if len(temperature_mode) == 1:
            temperature_list = [float(temperature_mode[0]), float(temperature_mode[0])]
        if len(temperature_mode) == 2:
            temperature_list = [float(temperature_mode[0]), float(temperature_mode[1])]
        d_out = float(experiment_log.obtainItem(name,'d_out',regular)[0])
        gauge_length = float(experiment_log.obtainItem(name,'gauge_length',regular)[0])
        axial_strain = float(experiment_log.obtainItem(name,'axial_strain',regular)[0])
        angel_strain = float(experiment_log.obtainItem(name,'angel_strain',regular)[0])
        equivalent_strain = float(experiment_log.obtainItem(name,'equivalent_strain',regular)[0])
        period = float(experiment_log.obtainItem(name,'period',regular)[0])
        axial_temperature_phase = float(experiment_log.obtainItem(name,'axial_temperature_phase',regular)[0])
        
        
        filename = ExperimentDirectory + name + '.csv'
        experiment = ExperimentData(filename)
        print experiment.axial_count_index_list
        e_list = []
        count = []
        for i in experiment.axial_count_index_list[2:200:50]:
            strain = experiment.obtainNthCycle('axial_strain',i)
            stress = experiment.obtainNthCycle('axial_stress',i)
            masing = obtain_masing_curve([strain,stress])
            e = obtain_youngs_modulus(masing,0.002)
            count.append(i)
            e_list.append(e)

            plot_data.addLine(masing[0]*100,
                              masing[1],
                              xlabel=xlabel,
                              ylabel=ylabel,
                              linelabel=str(int(temperature_mode[0])+273)+'K',
                              linewidth=2,
                              linestyle='-',
                              marker=None,
                              markersize=12,
                              color='auto')
    
#        plot_data.addLine(count,
#                          e_list,
#                          xlabel=xlabel,
#                          ylabel=ylabel,
#                          linelabel='',
#                          linewidth=2,
#                          linestyle='',
#                          marker='o',
#                          markersize=12,
#                          color='auto')
    
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
#    plt.xlim(0,1)
#    plt.ylim(0,800)
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
#    ax.xaxis.set_major_locator(MultipleLocator(0.5))
#    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
#    ax.xaxis.set_major_formatter(ScalarFormatter())
#    ax.yaxis.set_major_locator(MultipleLocator(200))
#    ax.yaxis.set_minor_locator(MultipleLocator(100))
#    ax.yaxis.set_major_formatter(ScalarFormatter())
#==============================================================================
# show legend
#==============================================================================
#    plt.legend(loc=0)
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
figure_name='plot_exp_masing'
create_plot_data_exp_half_life_cycle(figure_path,figure_name)
plot_exp_half_life_cycle(figure_path,figure_name,save_types=['.pdf'])

shutil.copy(__file__,ArticleFigureDirectory)