# -*- coding: utf-8 -*-
"""
Created on Thu Jan 05 11:50:47 2017

@author: j.Sun
"""

import matplotlib.pyplot as plt
import shutil
from Data import FatigueData,PlotData
from Constants import xylabels,FatigueDirectory,ArticleFigureDirectory
from plot_format import plot_format

def create_plot_data_fatigue_life(fatigue_data,figure_path=None,figure_name=None):
#==============================================================================
# x,y label
#==============================================================================
    xlabel = xylabels['experimental_life']
    ylabel = xylabels['predicted_life']
#==============================================================================
# plot lines
#==============================================================================
    i = 0
    marker_list = ['s','o','^','<','D','>','*']
    plot_data = PlotData()
#    for load_type in ['TC-IP','TC-OP','TC-90','PRO-IP','NPR-IP']:
    for load_type in ['TC-IP','TC-OP','TC-90','PRO-IP','NPR-IP','TC-IP-TGMF','TC-OP-TGMF']:
        experimental_life = fatigue_data.loadTypeFilter(load_type,'experimental_life')
        predicted_life = fatigue_data.loadTypeFilter(load_type,'predicted_life')
        plot_data.addLine(experimental_life,
                          predicted_life,
                          xlabel=xlabel,
                          ylabel=ylabel,
                          linelabel=load_type,
                          linewidth=2,
                          linestyle='',
                          marker=marker_list[i],
                          markersize=12)
        i += 1
    plot_data.writeToFile(figure_path,figure_name)
    
def plot_fatigue_life(figure_path=None,figure_name=None,save_types=[]):
#==============================================================================
# title
#==============================================================================
    title=''
#==============================================================================
# figure format
# http://matplotlib.org/users/customizing.html?highlight=rcparams
#==============================================================================
    plot_format()
    fig, ax = plt.subplots()
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
    plt.xlim(1E1,1E5)
    plt.ylim(1E1,1E5)
#==============================================================================
# xy axial equal
#==============================================================================
#    plt.gca().set_aspect('equal')
    plt.gca().set_aspect('auto')
#==============================================================================
# xy log scale
#==============================================================================
    plt.xscale('log')
    plt.yscale('log')
#==============================================================================
# plot lines
#==============================================================================
    plot_data = PlotData()
    plot_data.readFromFile(figure_path,figure_name)
    plot_data.plot()
#==============================================================================
# plot 1x lines
#==============================================================================
    linewidth = 1.0
    plt.plot([10,1e5],[10,1e5],color='black',linewidth=linewidth)
#==============================================================================
# plot 2x lines
#==============================================================================
    linewidth = 0.5
    plt.plot([20,1e5],[10,5e4],color='black',linewidth=linewidth)
    plt.plot([10,5e4],[20,1e5],color='black',linewidth=linewidth)
#==============================================================================
# plot 5x lines
#==============================================================================
    plt.plot([50,1e5],[10,2e4],color='black',linewidth=linewidth)
    plt.plot([10,2e4],[50,1e5],color='black',linewidth=linewidth)
#==============================================================================
# show legend
#==============================================================================
    plt.legend(loc=2)
    plt.show()
#==============================================================================
# save figures
#==============================================================================
    if figure_path <> None and figure_name<> None:
        for save_type in save_types:
            plt.savefig(figure_path + figure_name + save_type, dpi=150)
            print 'save as', figure_path + figure_name + save_type
    plt.close()

#fatigue_model_list = ['BM','FS','SWT','Liu1','Liu2','Chu']
fatigue_model_list = ['BM']

for fatigue_model in fatigue_model_list:
    fatigue_file = '%s%s.csv' % (FatigueDirectory,fatigue_model)
    fatigue_data = FatigueData(fatigue_file)
    figure_path = ArticleFigureDirectory
    figure_name = 'NF-NP-'+fatigue_model
    create_plot_data_fatigue_life(fatigue_data,figure_path,figure_name)
    
for fatigue_model in fatigue_model_list:
    figure_path = ArticleFigureDirectory
    figure_name = 'NF-NP-'+fatigue_model
    plot_fatigue_life(figure_path,figure_name,save_types=['.pdf','.png'])
    
shutil.copy(__file__,ArticleFigureDirectory)