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

def create_plot_data_exp_nth_cycle(name,nth_cycle,figure_path=None,figure_name=None):
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

    filename = ExperimentDirectory + name + '.csv'
    experiment = ExperimentData(filename)
    strain = experiment.obtainNthCycle('axial_strain',nth_cycle)
    stress = experiment.obtainNthCycle('axial_stress',nth_cycle)
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
    
    plot_data.writeToFile(figure_path,figure_name)
    
def plot_exp_nth_cycle(figure_title='',figure_path=None,figure_name=None,save_types=[]):
#==============================================================================
# title
#==============================================================================
    title=figure_title
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
#    plt.show()
    plt.close()

for name in experiment_type_dict['TC-OP']+experiment_type_dict['TC-IP']:
    print name
#    filename = ExperimentDirectory + name + '.csv'
#    experiment = ExperimentData(filename)
#    nth_cycle = experiment.half_life_cycle
    nth_cycle = 1
    figure_name = name + '_' + '1st'
    figure_title = '#%s 1st cycle' % (name)
    figure_path=ArticleFigureDirectory
    create_plot_data_exp_nth_cycle(name,nth_cycle,figure_path=figure_path,figure_name=figure_name)
    plot_exp_nth_cycle(figure_title=figure_title,figure_path=figure_path,figure_name=figure_name,save_types=['.png'])
    
    shutil.copy(__file__,ArticleFigureDirectory)