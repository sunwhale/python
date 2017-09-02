# -*- coding: utf-8 -*-
"""
Created on Thu Jan 05 11:50:47 2017

@author: j.Sun
"""

import numpy as np
import matplotlib.pyplot as plt
import shutil
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator,ScalarFormatter,FormatStrFormatter 
from Data import FatigueData,PlotData,ExperimentData
from Constants import *
from plot_format import plot_format
from Functions import calculate_elastic_by_temperature_in718
from Material import material_in718,material_in718_NASA,material_in718_BHU

def create_plot_data_elastic_by_temperature_in718(plot_data):
#==============================================================================
# x,y label
#==============================================================================
    xlabel = xylabels['temperature']
    ylabel = 'Young\'s and Shear Modulus [GPa]'
#==============================================================================
# plot lines
#==============================================================================
    temperature = [16,100,200,300,400,500,600,700]
    youngs_modulus = [208209.8593,205706.8624,198567.8644,193263.6415,185010.342,177691.2534,170836.6662,157353.0804]
    shear_modulus = [80786.07201,79741.37215,77069.32555,75333.41595,72269.49729,69436.2758,66462.31754,60762.26995]
    poisson_ratio = [0.28865,0.28984,0.28824,0.28272,0.28,0.27953,0.28521,0.29483]
    
    temperature = np.array(temperature)
    youngs_modulus = np.array(youngs_modulus)
    shear_modulus = np.array(shear_modulus)
    poisson_ratio = np.array(poisson_ratio)
    
    marker_list = ['s','o','^','D']
    
    i = 0
    plot_data.addLine(temperature,
                      youngs_modulus/1000.0,
                      xlabel=xlabel,
                      ylabel=ylabel,
                      linelabel='Young\'s Modulus',
                      linewidth=2,
                      linestyle='',
                      marker=marker_list[i],
                      markersize=12,
                      color='green')
    i += 1
    plot_data.addLine(temperature,
                      shear_modulus/1000.0,
                      xlabel=xlabel,
                      ylabel=ylabel,
                      linelabel='Shear Modulus',
                      linewidth=2,
                      linestyle='',
                      marker=marker_list[i],
                      markersize=12,
                      color='blue')
    i += 1
    plot_data.addLine(temperature,
                      poisson_ratio,
                      xlabel=xlabel,
                      ylabel='Poisson\'s Ratio',
                      linelabel='Poisson\'s Ratio',
                      linewidth=2,
                      linestyle='',
                      marker=marker_list[i],
                      markersize=12,
                      color='red')
    
    temperature = range(0,800,1)
    youngs_modulus = []
    shear_modulus = []
    poisson_ratio = []
    for t in temperature:
        y,p,s,ys = calculate_elastic_by_temperature_in718(t)
        youngs_modulus.append(y)
        shear_modulus.append(s)
        poisson_ratio.append(p)
        
    temperature = np.array(temperature)
    youngs_modulus = np.array(youngs_modulus)
    shear_modulus = np.array(shear_modulus)
    poisson_ratio = np.array(poisson_ratio)
    
    plot_data.addLine(temperature,
                      youngs_modulus/1000.0,
                      xlabel=xlabel,
                      ylabel=ylabel,
                      linelabel='',
                      linewidth=2,
                      linestyle='-',
                      marker=None,
                      markersize=12,
                      color='green')
    
    plot_data.addLine(temperature,
                      shear_modulus/1000.0,
                      xlabel=xlabel,
                      ylabel=ylabel,
                      linelabel='',
                      linewidth=2,
                      linestyle='-',
                      marker=None,
                      markersize=12,
                      color='blue')
    
    plot_data.addLine(temperature,
                      poisson_ratio,
                      xlabel=xlabel,
                      ylabel='Poisson\'s Ratio',
                      linelabel='',
                      linewidth=2,
                      linestyle='-',
                      marker=None,
                      markersize=12,
                      color='red')
                      
    plot_data.writeToFile()

def plot_elastic_by_temperature_in718(plot_data):
#==============================================================================
# title
#==============================================================================
    title=''
#==============================================================================
# figure format
# http://matplotlib.org/users/customizing.html?highlight=rcparams
#==============================================================================
    plot_format()
    mpl.rcParams['figure.subplot.right'] = 0.9
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
    plt.xlim(0,750)
    plt.ylim(0,220)
#==============================================================================
# xy log scale
#==============================================================================
#    plt.xscale('log')
#    plt.yscale('log')
#==============================================================================
# xy axial equal
#==============================================================================
    ax1 = plt.gca()
    ax1.set_aspect('auto')
#==============================================================================
# plot lines
#==============================================================================
    plot_data.readFromFile()
    plot_data.plot([0,1,3,4])
    ax1.legend(loc=3)
    
    ax2 = ax1.twinx()  # this is the important function
    plot_data.plot([2,5])
    ax2.set_xlim(0,750)
    ax2.set_ylim(0.1,0.4)
    ax2.legend(loc=4)
#==============================================================================
# http://stackoverflow.com/questions/21920233/matplotlib-log-scale-tick-label-number-formatting
#==============================================================================
    ax1.xaxis.set_major_locator(MultipleLocator(100))
    ax1.xaxis.set_minor_locator(MultipleLocator(50))
    ax1.xaxis.set_major_formatter(ScalarFormatter())
    ax1.yaxis.set_major_locator(MultipleLocator(100))
    ax1.yaxis.set_minor_locator(MultipleLocator(10))
    ax1.yaxis.set_major_formatter(ScalarFormatter())
    ax2.yaxis.set_major_locator(MultipleLocator(0.1))
    ax2.yaxis.set_minor_locator(MultipleLocator(0.01))
    ax2.yaxis.set_major_formatter(ScalarFormatter())
#==============================================================================
# show legend
#==============================================================================
    
#==============================================================================
# save figures
#==============================================================================
    plot_data.saveFigure()
    plt.show()
    plt.close()

figure_path = ArticleFigureDirectory
figure_name = 'plot_elastic_by_temperature_in718'
plot_data = PlotData(figure_path,figure_name,save_types=['.pdf'])

create_plot_data_elastic_by_temperature_in718(plot_data)
plot_elastic_by_temperature_in718(plot_data)

shutil.copy(__file__,ArticleFigureDirectory)