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

def create_plot_data_exp_temperature(figure_path=None,figure_name=None):
#==============================================================================
# x,y label
#==============================================================================
    xlabel = xylabels['runing_time']
    ylabel = xylabels['temperature']
#==============================================================================
# plot lines
#==============================================================================
    i = 0
    marker_list = ['s','o','^','D']
    plot_data = PlotData()
#    for name in ['30_0_','50_0_','70_0_','90_0_']:
#    for name in ['30_25_','50_25_','70_25_','90_25_']:
#    for name in ['30_50_','50_50_','70_50_','90_50_']:
    for name in ['30_67_','50_67_','70_67_','90_67_']:
#    for name in ['90_0_','90_25_','90_50_','90_67_']:
#    for name in ['70_0_','70_25_','70_50_','70_67_']:
#    for name in ['50_0_','50_25_','50_50_','50_67_']:
#    for name in ['30_0_','30_25_','30_50_','30_67_']:
        ExperimentDirectory = 'F:\\Database\\Temperature2\\'
        filename = ExperimentDirectory + name + '.csv'
        experiment = ExperimentData(filename)
        time = experiment.runing_time
        temperature = experiment.temperature
        plot_data.addLine(time,
                          temperature,
                          xlabel=xlabel,
                          ylabel=ylabel,
                          linelabel=name.split('_')[0] + '%',
                          linewidth=2,
                          linestyle='-',
                          marker=None,
                          markersize=12,
                          color='auto')
        i += 1
    
    plot_data.writeToFile(figure_path,figure_name)
    
def plot_exp_temperature(figure_path=None,figure_name=None,save_types=[]):
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
    plt.xlim(1,1e4)
    plt.ylim(0,1100)
#==============================================================================
# xy log scale
#==============================================================================
    plt.xscale('log')
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
#    ax.yaxis.set_major_locator(MultipleLocator(500))
#    ax.yaxis.set_minor_locator(MultipleLocator(100))
#    ax.yaxis.set_major_formatter(ScalarFormatter())
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

ArticleFigureDirectory = 'F:\\Cloud\\GitHub\\furnace\\figs\\'
figure_path=ArticleFigureDirectory
figure_name='plot_exp_temperature_67'
create_plot_data_exp_temperature(figure_path,figure_name)
plot_exp_temperature(figure_path,figure_name,save_types=['.pdf','.png'])

shutil.copy(__file__,ArticleFigureDirectory)