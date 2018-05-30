# -*- coding: utf-8 -*-
"""
Created on Tue Dec 18 13:49:06 2012

@author: Sun
"""

import sys
sys.path.append('F:\\GitHub\\python')

import matplotlib.pyplot as plt
import shutil
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator,ScalarFormatter,FormatStrFormatter 
from Data import *
from plot_format import plot_format
from Constants_reflection import *

def create_plot_data(figure_path=None,figure_name=None):
#==============================================================================
# x,y label
#==============================================================================
    xlabel = '$y$ coordinate [mm]'
    ylabel = 'Temperature [K]'
#==============================================================================
# plot lines
#==============================================================================
    plot_data = PlotData()

    node = 16*[0]
    node[0] = [0,0]
    node[1] = [25e-3,0]
    node[2] = [15e-3,0]
    node[3] = [10e-3,0]
    node[4] = [5e-3,0]
    node[5] = [0e-3,0]
    node[6] = [-5e-3,0]
    node[7] = [-10e-3,0]
    node[8] = [-15e-3,0]
    node[9] = [-25e-3,0]
    node[10] = [0,10e-3]
    node[11] = [0,20e-3]
    node[12] = [0,30e-3]
    node[13] = [0,-10e-3]
    node[14] = [0,-20e-3]
    node[15] = [0,-30e-3]
    
    node_x = [ n[0] for n in node]
    node_y = [ n[1] for n in node]
    
#    begin_time = '11:16:58'
    begin_time = '19:29:43'
    num_node = 10

    txt = SimulationDirectory + "VERTICALLINE.txt"
    list1 = open(txt,"r").readlines()	
    list2 = [i.strip() for i in list1]
    list3 = [j.split(',') for j in list2]
    y_coordinate = [ float(k[2])  for k in list3 ]
    temp_sim = [ float(k[4])  for k in list3 ]
    
    txt = SimulationDirectory + "VERTICALLINE_400W_200s_0.80.txt"
    list1 = open(txt,"r").readlines()	
    list2 = [i.strip() for i in list1]
    list3 = [j.split(',') for j in list2]
    y_coordinate = [ float(k[2])  for k in list3 ]
    temp_sim_80 = [ float(k[4])  for k in list3 ]
    
    txt = SimulationDirectory + "VERTICALLINE_400W_200s_0.90.txt"
    list1 = open(txt,"r").readlines()	
    list2 = [i.strip() for i in list1]
    list3 = [j.split(',') for j in list2]
    y_coordinate = [ float(k[2])  for k in list3 ]
    temp_sim_90 = [ float(k[4])  for k in list3 ]

    txt = SimulationDirectory + "VERTICALLINE_400W_200s_1.00.txt"
    list1 = open(txt,"r").readlines()	
    list2 = [i.strip() for i in list1]
    list3 = [j.split(',') for j in list2]
    y_coordinate = [ float(k[2])  for k in list3 ]
    temp_sim_100 = [ float(k[4])  for k in list3 ]

#    txt = ExperimentDirectory + "ExperimentalData\\7\\2013_08_21_B11h16m58s_E11h18m28s_800W_NI9213.txt"
    txt = ExperimentDirectory + "ExperimentalData\\4\\2013_08_20_B19h29m43s_E19h34m13s_400W_NI9213.txt"
    list1 = open(txt,"r").readlines()	
    list2 = [i.strip() for i in list1 if len(i)>50]
    list3 = [j.split('\t') for j in list2]
    system_time =  [ k[0]  for k in list3 ]
    num = [ k[0]  for k in list3 ].index(begin_time)
    time_exp = [ float(k[1])  for k in list3 ]
    time_shift = time_exp[num]
    time_exp_shift = [ t - time_shift  for t in time_exp ]
    temp_exp = [ float(k)  for k in list3[num+200][2:] ]

    i = 0

    plot_data.addLine(np.array(y_coordinate)*1000,
                      np.array(temp_sim)+273,
                      xlabel=xlabel,
                      ylabel=ylabel,
                      linelabel='Sim. 80',
                      linewidth=2,
                      linestyle='-',
                      marker=None,
                      markersize=12,
                      color='blue')
                      
    plot_data.addLine(np.array(y_coordinate)*1000,
                      np.array(temp_sim_80)+273,
                      xlabel=xlabel,
                      ylabel=ylabel,
                      linelabel='Sim. 80',
                      linewidth=2,
                      linestyle='-',
                      marker=None,
                      markersize=12,
                      color='black')
                      
    plot_data.addLine(np.array(y_coordinate)*1000,
                      np.array(temp_sim_90)+273,
                      xlabel=xlabel,
                      ylabel=ylabel,
                      linelabel='Sim. 90',
                      linewidth=2,
                      linestyle='-',
                      marker=None,
                      markersize=12,
                      color='black')

    plot_data.addLine(np.array(y_coordinate)*1000,
                      np.array(temp_sim_100)+273,
                      xlabel=xlabel,
                      ylabel=ylabel,
                      linelabel='Sim. 100',
                      linewidth=2,
                      linestyle='-',
                      marker=None,
                      markersize=12,
                      color='black')
                      
    plot_data.addLine(np.array(node_y[10:16]+node_y[5:6])*1000,
                      np.array(temp_exp[10:16]+temp_exp[5:6])+273,
                      xlabel=xlabel,
                      ylabel=ylabel,
                      linelabel='Exp.',
                      linewidth=2,
                      linestyle='',
                      marker='o',
                      markersize=12,
                      color='red')
    
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
    plt.xlim(-40,40)
    plt.ylim(273,350+273)
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
    print figure_path,figure_name
    plot_data.readFromFile(figure_path,figure_name)
    
    plot_data.plot()
#==============================================================================
# http://stackoverflow.com/questions/21920233/matplotlib-log-scale-tick-label-number-formatting
#==============================================================================
#    ax.xaxis.set_major_locator(MultipleLocator(500))
#    ax.xaxis.set_minor_locator(MultipleLocator(100))
#    ax.xaxis.set_major_formatter(ScalarFormatter())
#    ax.yaxis.set_major_locator(MultipleLocator(500))
#    ax.yaxis.set_minor_locator(MultipleLocator(100))
#    ax.yaxis.set_major_formatter(ScalarFormatter())
#==============================================================================
# show legend
#==============================================================================
    lg = plt.legend(title='',loc=0, frameon = 1)
    lg.get_frame().set_color('white')
    lg.get_frame().set_alpha(0.8)
    title = lg.get_title()
    title.set_fontsize(16)
    plt.plot([0,0],[0,1000],lw=1.5,color='black')
#==============================================================================
# save figures
#==============================================================================
    if figure_path <> None and figure_name<> None:
        for save_type in save_types:
            plt.savefig(figure_path + figure_name + save_type, dpi=150, transparent=True)
            print 'save as', figure_path + figure_name + save_type
    plt.show()
    plt.close()

if __name__ == '__main__':
    ArticleFigureDirectory = 'F:\\Cloud\\GitHub\\furnace\\figs\\'
    figure_path=ArticleFigureDirectory
    figure_name='plot_compare_horizontal_reflection'
    create_plot_data(figure_path,figure_name)
    plot(figure_path,figure_name,save_types=['.pdf','.png'])
    
    shutil.copy(__file__,ArticleFigureDirectory)