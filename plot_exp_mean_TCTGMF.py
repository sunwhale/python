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
#    color_list = ['blue','red','black','green','yellow','orange','magenta','cyan']
    plot_data = PlotData()
    experiment_log = ExperimentLog(ExperimentLogFile)
#    for name in experiment_type_dict['TC-IP-TGMF']+experiment_type_dict['TC-OP-TGMF']+['7301']:
    for name in experiment_type_dict['TC-IP-TGMF'] + experiment_type_dict['TC-OP-TGMF']:
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
        
        filename = ExperimentDirectory + name + '.csv'
        experiment = ExperimentData(filename)
        cycle,peak,valley = experiment.obtainPeakValley('axial_stress')
#        cycle,peak,valley = experiment.obtainPeakValley('total_strain')
        mean = (np.array(peak) + np.array(valley))/2
        ylabel = xylabels['mean_stress']
#        strain_amp = (np.array(peak) - np.array(valley))
#        
#        plt.plot(strain_amp[:10])
#        plt.show()
#        material = material_in718()
#        print name,axial_displacement,material.calcStrainAmplitude(peak[int(life/2.0)])*2.0,peak[int(life/2.0)],valley[int(life/2.0)],life
#        print name,axial_displacement,material.calcStrain(-1.0*valley[0]),(peak[0]-valley[0])/2.0,life
        
#        plot_data.addLine(cycle,
#                          peak,
#                          xlabel=xlabel,
#                          ylabel=ylabel,
#                          linelabel=str(axial_strain) + '%',
#                          linewidth=2,
#                          linestyle='-',
#                          marker=None,
#                          markersize=12,
#                          color=color_list[i])
        plot_data.addLine(cycle,
                          mean,
                          xlabel=xlabel,
                          ylabel=ylabel,
                          linelabel=str(axial_strain) + '%',
                          linewidth=2,
                          linestyle='-',
                          marker=None,
                          markersize=12,
                          color=color_list[i])
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
    mpl.rcParams['figure.figsize'] = (12,9)
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
    plt.xlim(1,10000)
    plt.ylim(-150,150)
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
    ax.yaxis.set_major_locator(MultipleLocator(50))
    ax.yaxis.set_minor_locator(MultipleLocator(10))
    ax.yaxis.set_major_formatter(ScalarFormatter())
#==============================================================================
# annotate
#==============================================================================
    plt.text(5,-100,r'TC-IP-TGMF',fontsize=14)
    plt.text(5,100,r'TC-OP-TGMF',fontsize=14)
#    plt.text(10,-700,r'Valley stress',fontsize=14)
    
#    plt.annotate(r'Peak stresses',xy=(1,1000),xytext=(100,1000),fontsize=14,color='black',arrowprops=dict(arrowstyle='->',color='black'))
#    plt.annotate(r'Mean stresses',xy=(1,0),xytext=(100,0),fontsize=14,color='black',arrowprops=dict(arrowstyle='->',color='black'))
#    plt.annotate(r'Valley stresses',xy=(1,-1000),xytext=(100,-1000),fontsize=14,color='black',arrowprops=dict(arrowstyle='->',color='black'))
#==============================================================================
# show legend
#==============================================================================
    lg = plt.legend(title='$\Delta\\varepsilon/2$',loc=4,frameon=1)
    lg.get_frame().set_color('white')
    lg.get_frame().set_alpha(0.9)
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
figure_name='plot_exp_mean_TCTGMF'
#create_plot_data(figure_path,figure_name)
plot(figure_path,figure_name,save_types=['.pdf'])

shutil.copy(__file__,ArticleFigureDirectory)