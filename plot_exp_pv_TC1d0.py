# -*- coding: utf-8 -*-
"""
Created on Thu Jan 05 11:50:47 2017

@author: j.Sun
"""

import matplotlib.pyplot as plt
import shutil
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator,ScalarFormatter,FormatStrFormatter 
from Data import *
from Constants import *
from plot_format import plot_format
from Material import material_in718,material_in718_NASA,material_in718_BHU

def create_plot_data(figure_path=None,figure_name=None):
#==============================================================================
# x,y label
#==============================================================================
    xlabel = xylabels['axial_count']
    ylabel = xylabels['axial_stress']
#==============================================================================
# plot lines
#==============================================================================
    i = 0
    marker_list = ['s','o','^','D']
    label_list = ['TC-IF','TC-IP','TC-OP','TC-90','PRO-IP','NPR-IP']
#    color_list = ['blue','red','black','green']
    plot_data = PlotData()
    experiment_log = ExperimentLog(ExperimentLogFile)
#    for name in ['7114','7018','7017','7025','7038','7037']:
    for name in ['7114','7018']:
        print name

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
        cycle,peak,valley = experiment.obtainPeakValley('axial_stress')
        mean = (np.array(peak) + np.array(valley))/2
        plot_data.addLine(cycle,
                          peak,
                          xlabel=xlabel,
                          ylabel=ylabel,
                          linelabel=label_list[i],
                          linewidth=2,
                          linestyle='-',
                          marker=None,
                          markersize=12,
                          color=color_list[i])
#        plot_data.addLine(cycle,
#                          mean,
#                          xlabel=xlabel,
#                          ylabel=ylabel,
#                          linelabel='',
#                          linewidth=2,
#                          linestyle='-',
#                          marker=None,
#                          markersize=12,
#                          color=color_list[i])
#        plot_data.addLine(cycle,
#                          valley,
#                          xlabel=xlabel,
#                          ylabel=ylabel,
#                          linelabel='',
#                          linewidth=2,
#                          linestyle='-',
#                          marker=None,
#                          markersize=12,
#                          color=color_list[i])
        i += 1
    
    plot_data.writeToFile(figure_path,figure_name)
    
def plot(figure_path=None,figure_name=None,save_types=[]):
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
#    plt.xlim(-1,1)
    plt.ylim(600,1200)
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
    ax.yaxis.set_minor_locator(MultipleLocator(100))
#    ax.yaxis.set_major_formatter(ScalarFormatter())
#==============================================================================
# show legend
#==============================================================================
    lg = plt.legend(loc=1)
    title = lg.get_title()
    title.set_fontsize(16)
#==============================================================================
# annotate
#==============================================================================
#    plt.text(10,1000,r'Peak stress',fontsize=14)
#    plt.text(10,0,r'Mean stress',fontsize=14)
#    plt.text(10,-900,r'Valley stress',fontsize=14)
    
#    plt.annotate(r'Peak stresses',xy=(1,1000),xytext=(100,1000),fontsize=14,color='black',arrowprops=dict(arrowstyle='->',color='black'))
#    plt.annotate(r'Mean stresses',xy=(1,0),xytext=(100,0),fontsize=14,color='black',arrowprops=dict(arrowstyle='->',color='black'))
#    plt.annotate(r'Valley stresses',xy=(1,-1000),xytext=(100,-1000),fontsize=14,color='black',arrowprops=dict(arrowstyle='->',color='black'))
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
figure_name='plot_exp_pv_TC1d0'
create_plot_data(figure_path,figure_name)
plot(figure_path,figure_name,save_types=['.pdf'])

shutil.copy(__file__,ArticleFigureDirectory)