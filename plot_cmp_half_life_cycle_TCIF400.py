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
    xlabel = xylabels['axial_strain']
    ylabel = xylabels['axial_stress']
#==============================================================================
# plot lines
#==============================================================================
    i = 0
#    marker_list = ['s','o','^','D']
#    color_list = ['blue','red','black','green','yellow','orange','magenta','cyan']
    plot_data = PlotData()
    experiment_log = ExperimentLog(ExperimentLogFile)
#    for name in experiment_type_dict['TC-IP-TGMF']+experiment_type_dict['TC-OP-TGMF']+['7301']:
    for name in ['7042']:
#    for name in ['7208']:
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
        axial_displacement = float(experiment_log.obtainItem(name,'axial_displacement',regular)[0])
        angel_strain = float(experiment_log.obtainItem(name,'angel_strain',regular)[0])
        equivalent_strain = float(experiment_log.obtainItem(name,'equivalent_strain',regular)[0])
        period = float(experiment_log.obtainItem(name,'period',regular)[0])
        axial_temperature_phase = float(experiment_log.obtainItem(name,'axial_temperature_phase',regular)[0])
        life = float(experiment_log.obtainItem(name,'comments',regular)[0])
        
        exp_filename = ExperimentDirectory + name + '.csv'
        experiment = ExperimentData(exp_filename)
        strain = experiment.obtainNthCycle('axial_strain',experiment.half_life_cycle)
        stress = experiment.obtainNthCycle('axial_stress',experiment.half_life_cycle)
    
        plot_data.addLine(strain*100,
                          stress,
                          xlabel=xlabel,
                          ylabel=ylabel,
                          linelabel='Exp.',
                          linewidth=2,
                          linestyle='',
                          marker=marker_list[i],
                          markersize=12,
                          color=color_list[i],
                          skip=10,
                          log_skip=1)
                          
        sim_filename = SimulationDirectory + name + '.csv'
        simulation = SimulationData(sim_filename,period)
        strain = simulation.obtainNthCycle('axial_strain',simulation.total_axial_count-1)
        stress = simulation.obtainNthCycle('axial_stress',simulation.total_axial_count-1)
        
        plot_data.addLine(strain*100,
                          stress,
                          xlabel=xlabel,
                          ylabel=ylabel,
                          linelabel='Sim.',
                          linewidth=2,
                          linestyle='-',
                          marker=None,
                          markersize=12,
                          color=color_list[i])
        
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
    plt.xlim(-1,1)
    plt.ylim(-1200,1200)
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
    ax.yaxis.set_major_locator(MultipleLocator(400))
    ax.yaxis.set_minor_locator(MultipleLocator(100))
    ax.yaxis.set_major_formatter(ScalarFormatter())
#==============================================================================
# show legend
#==============================================================================
    lg = plt.legend(title='',loc=0)
    title = lg.get_title()
    title.set_fontsize(16)
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
figure_name='plot_cmp_half_life_cycle_TCIF400'
create_plot_data(figure_path,figure_name)
plot(figure_path,figure_name,save_types=['.pdf'])

shutil.copy(__file__,ArticleFigureDirectory)