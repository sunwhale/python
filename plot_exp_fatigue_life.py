# -*- coding: utf-8 -*-
"""
Created on Thu Jan 05 11:50:47 2017

@author: j.Sun
"""

import shutil
import matplotlib.pyplot as plt
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
    ylabel = xylabels['equivalent_strain_amplitude']
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
        experimental_life = fatigue_data.loadTypeFilter(load_type,'experimental_life')
        equivalent_strain_amplitude = fatigue_data.loadTypeFilter(load_type,'equivalent_strain_amplitude')
        plot_data.addLine(experimental_life,
                          equivalent_strain_amplitude,
                          xlabel=xlabel,
                          ylabel=ylabel,
                          linelabel=load_type,
                          linewidth=2,
                          linestyle='',
                          marker=marker_list[i],
                          markersize=12,
                          color='auto')
        i += 1

    material = material_in718()
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
            plt.savefig(figure_path + figure_name + save_type, dpi=150)
            print 'save as', figure_path + figure_name + save_type
    plt.show()
    plt.close()

fatigue_file = '%s%s.csv' % (FatigueDirectory,'BM')
fatigue_data = FatigueData(fatigue_file)
figure_path = ArticleFigureDirectory
figure_name = 'plot_exp_fatigue_life'
create_plot_data_exp_fatigue_life(fatigue_data,figure_path,figure_name)
plot_exp_fatigue_life(figure_path,figure_name,save_types=['.pdf','.png'])

shutil.copy(__file__,ArticleFigureDirectory)