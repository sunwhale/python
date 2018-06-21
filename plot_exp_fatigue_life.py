# -*- coding: utf-8 -*-
"""
Created on Thu Jan 05 11:50:47 2017

@author: j.Sun
"""

import shutil
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator,ScalarFormatter,FormatStrFormatter 
from Data import FatigueData,PlotData
from Constants import *
from plot_format import plot_format
from Material import material_in718

def create_plot_data_exp_fatigue_life(fatigue_data,figure_path=None,figure_name=None):
#==============================================================================
# x,y label
#==============================================================================
    xlabel = xylabels['experimental_life']
#    ylabel = xylabels['equivalent_strain_amplitude']
    ylabel = xylabels['axial_strain_amplitude']
    label_dict = {'TC-IP':'TMF-IP',
                  'TC-OP':'TMF-OP',
                  'TC-IP-TGMF':'TGMF-IP',
                  'TC-OP-TGMF':'TGMF-OP',
                  'TC-IF':'IF',
                  'TC-IP-TGMF-TBC':'TGMF-IP-TBC',
                  }
#==============================================================================
# plot lines
#==============================================================================
    i = 0
    marker_list = ['s','o','^','D']
    plot_data = PlotData()
    i = 0
    marker_list = ['s','o','^','D','<','^','>']
#    for load_type in ['TC-IP','TC-OP','PRO-IP','NPR-IP','TC-90','TC-IP-TGMF','TC-OP-TGMF']:
    for load_type in ['TC-IP','TC-OP','TC-IP-TGMF','TC-OP-TGMF','TC-IP-TGMF-TBC']:
#    for load_type in ['TC-IP-TGMF','TC-IP-TGMF-TBC']:
#    for load_type in ['TC-IP','TC-IP-TGMF']:
#    for load_type in ['TC-OP','TC-OP-TGMF']:
        experimental_life = fatigue_data.loadTypeFilter(load_type,'experimental_life')
        equivalent_strain_amplitude = fatigue_data.loadTypeFilter(load_type,'equivalent_strain_amplitude')
        plot_data.addLine(experimental_life,
                          equivalent_strain_amplitude,
                          xlabel=xlabel,
                          ylabel=ylabel,
#                          linelabel=label_dict[load_type],
                          linelabel=load_type,
                          linewidth=2,
                          linestyle='',
                          marker=marker_list[i],
                          markersize=12,
                          color='auto')
        i += 1

#        material = material_in718()
#        material.calculateMansonCoffinAxial(np.array(equivalent_strain_amplitude)/100.0,np.array(experimental_life)*2)
#        life,epsilon_amplitude = material.plotMansonCoffinAxial()
#        plot_data.addLine(life,
#                          epsilon_amplitude*100,
#                          xlabel=xlabel,
#                          ylabel=ylabel,
#                          linelabel='',
#                          linewidth=2,
#                          linestyle='-',
#                          marker=None,
#                          markersize=12,
#                          color='black')
                      
    material = material_in718()
#    material.calculateMansonCoffinAxial(np.array(equivalent_strain_amplitude)/100.0,np.array(experimental_life)*2)
    life,epsilon_amplitude = material.plotMansonCoffinAxial()
    plot_data.addLine(life,
                      epsilon_amplitude*100,
                      xlabel=xlabel,
                      ylabel=ylabel,
                      linelabel='TC-IF',
                      linewidth=2,
                      linestyle='-',
                      marker=None,
                      markersize=12,
                      color='black')
    
    plot_data.writeToFile(figure_path,figure_name)
    
def plot_exp_fatigue_life(figure_path=None,figure_name=None,save_types=[]):
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
    plt.xlim(1E1,1E4)
    plt.ylim(0.35,1.1)
#==============================================================================
# xy log scale
#==============================================================================
    plt.xscale('log')
    plt.yscale('log')
#==============================================================================
# xy axial equal
#==============================================================================
    ax = plt.gca()
#    ax.set_aspect('equal')
    ax.set_aspect('auto')
#==============================================================================
# http://stackoverflow.com/questions/21920233/matplotlib-log-scale-tick-label-number-formatting
#==============================================================================
    ax.yaxis.set_major_locator(MultipleLocator(0.1))
    ax.yaxis.set_minor_locator(MultipleLocator(0.01))
    ax.yaxis.set_major_formatter(ScalarFormatter())
#==============================================================================
# plot lines
#==============================================================================
    plot_data = PlotData()
    plot_data.readFromFile(figure_path,figure_name)
    plot_data.plot()
#==============================================================================
# show legend
#==============================================================================
    plt.legend(loc=0)
#    plt.legend(loc=0,fontsize='small',frameon=True,numpoints=1,title='Temperature')
#==============================================================================
# save figures
#==============================================================================
    if figure_path <> None and figure_name<> None:
        for save_type in save_types:
            plt.savefig(figure_path + figure_name + save_type, dpi=150, transparent=True)
            print 'save as', figure_path + figure_name + save_type
    plt.show()
    plt.close()

fatigue_file = '%s%s.csv' % (FatigueDirectory,'Study2')
fatigue_data = FatigueData(fatigue_file)
figure_path = ArticleFigureDirectory
#figure_name = 'plot_exp_fatigue_life'
figure_name = 'plot_exp_fatigue_life_TGMF'
#figure_name = 'plot_exp_fatigue_life_TGMF-TBC'
#figure_name = 'plot_exp_fatigue_life_IP'
#figure_name = 'plot_exp_fatigue_life_OP'
create_plot_data_exp_fatigue_life(fatigue_data,figure_path,figure_name)
plot_exp_fatigue_life(figure_path,figure_name,save_types=['.pdf','.png'])

shutil.copy(__file__,ArticleFigureDirectory)