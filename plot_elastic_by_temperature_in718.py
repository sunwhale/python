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
from Functions import calculate_elastic_by_temperature_in718
from Material import material_in718,material_in718_NASA,material_in718_BHU

def create_plot_data_elastic_by_temperature_in718(figure_path=None,figure_name=None):
#==============================================================================
# x,y label
#==============================================================================
    xlabel = xylabels['axial_strain']
    ylabel = xylabels['axial_stress']
#==============================================================================
# plot lines
#==============================================================================
    youngs_modulus = EMOD
    poisson_ratio = ENU
    shear_modulus = EG
    return 

    temperature = [16,100,200,300,400,500,600,700]
    youngs_modulus = [208209.8593,205706.8624,198567.8644,193263.6415,185010.342,177691.2534,170836.6662,157353.0804]
    shear_modulus = [77312.0387,75575.14351,73848.08965,72482.74175,70214.82371,68174.05477,65158.86763,60376.24277]
    poisson_ratio = [0.34656,0.36094,0.34443,0.33317,0.31746,0.30322,0.31092,0.3031]
    
    plot_data = PlotData()
    marker_list = ['s','o','^','D']
    
    i = 0
    plot_data.addLine(temperature,
                      youngs_modulus,
                      xlabel=xlabel,
                      ylabel=ylabel,
                      linelabel='',
                      linewidth=2,
                      linestyle='',
                      marker=marker_list[i],
                      markersize=12,
                      color='auto')
    i += 1
    plot_data.addLine(temperature,
                      shear_modulus,
                      xlabel=xlabel,
                      ylabel=ylabel,
                      linelabel='',
                      linewidth=2,
                      linestyle='',
                      marker=marker_list[i],
                      markersize=12,
                      color='auto')
    i += 1
    plot_data.addLine(temperature,
                      poisson_ratio,
                      xlabel=xlabel,
                      ylabel=ylabel,
                      linelabel='',
                      linewidth=2,
                      linestyle='',
                      marker=marker_list[i],
                      markersize=12,
                      color='auto')
    
    youngs_modulus,poisson_ratio,shear_modulus
    
    plot_data.addLine(temperature,
                      youngs_modulus,
                      xlabel=xlabel,
                      ylabel=ylabel,
                      linelabel='',
                      linewidth=2,
                      linestyle='-',
                      marker=None,
                      markersize=12,
                      color='auto')
    
    plot_data.addLine(temperature,
                      shear_modulus,
                      xlabel=xlabel,
                      ylabel=ylabel,
                      linelabel='',
                      linewidth=2,
                      linestyle='-',
                      marker=None,
                      markersize=12,
                      color='auto')
    
    plot_data.addLine(temperature,
                      poisson_ratio,
                      xlabel=xlabel,
                      ylabel=ylabel,
                      linelabel='',
                      linewidth=2,
                      linestyle='-',
                      marker=None,
                      markersize=12,
                      color='auto')
    plot_data.writeToFile(figure_path,figure_name)
    
def plot_elastic_by_temperature_in718(figure_path=None,figure_name=None,save_types=[]):
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
            plt.savefig(figure_path + figure_name + save_type, dpi=150)
            print 'save as', figure_path + figure_name + save_type
    plt.show()
    plt.close()

figure_path=ArticleFigureDirectory
figure_name='plot_exp_half_life_cycle'
#create_plot_data_exp_half_life_cycle(figure_path,figure_name)
plot_exp_half_life_cycle(figure_path,figure_name,save_types=['.pdf'])

shutil.copy(__file__,ArticleFigureDirectory)