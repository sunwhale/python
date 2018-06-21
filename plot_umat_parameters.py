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
    xlabel = xylabels['temperature']
    ylabel = '$r_0^k$ [MPa]'
#    ylabel = '$r_{\Delta \\rm{s}}^k$ [MPa]'
#    ylabel = '$a_1^k$'
#    ylabel = '$b_1^k$'
#==============================================================================
# plot lines
#==============================================================================
    i = 0
    plot_data = PlotData()

#    temperature_list = [300,550,650]
    temperature_list = [573,823,923]
    temperature = np.array(temperature_list)
    r0_list = [[292.65,60.11,74.7,89.45,82.86,91.27,113.42,135.81,140.01],
    [242.66,57.06,71.17,85.54,79.54,87.9,109.64,131.78,136.93],
    [227.08,54.94,68.47,82.22,76.38,84.34,105.11,126.21,130.89]]
    rdelta_list = [[-37.83,-10.89,-14.3,-18.06,-17.2,-20.2,-26.21,-32.72,-32.38],
    [-100.14,-13.56,-17.29,-21.23,-20.15,-22.68,-28.84,-35.34,-33.02],
    [-139,-15.16,-18.98,-22.89,-21.36,-23.68,-29.65,-35.76,-32.75]]
    a1_list = [0.4,0.41,0.42]
    b1_list = [1.27,0.87,0.49]
    b2_list = [9.45,11.8,9.18]
    r0 = np.array(r0_list)
    rdelta = np.array(rdelta_list)

    for i in range(len(r0[0])):
        plot_data.addLine(temperature,
                          r0[:,i],
                          xlabel=xlabel,
                          ylabel=ylabel,
                          linelabel=('$k=%s$' % str(i+1)),
                          linewidth=2,
                          linestyle='-',
                          marker=marker_list[i],
                          markersize=12,
                          color=color_list[i])

#    for i in range(len(r0[0])):
#        plot_data.addLine(temperature,
#                          rdelta[:,i],
#                          xlabel=xlabel,
#                          ylabel=ylabel,
#                          linelabel=('$k=%s$' % str(i+1)),
#                          linewidth=2,
#                          linestyle='-',
#                          marker=marker_list[i],
#                          markersize=12,
#                          color=color_list[i])
                          
#    plot_data.addLine(temperature,
#                      a1_list,
#                      xlabel=xlabel,
#                      ylabel=ylabel,
#                      linelabel='$a_1^k$',
#                      linewidth=2,
#                      linestyle='-',
#                      marker=marker_list[i],
#                      markersize=12,
#                      color=color_list[i])
#    i += 1
    
#    plot_data.addLine(temperature,
#                      b1_list,
#                      xlabel=xlabel,
#                      ylabel=ylabel,
#                      linelabel='$b_1^k$',
#                      linewidth=2,
#                      linestyle='-',
#                      marker=marker_list[i],
#                      markersize=12,
#                      color=color_list[i])
#    i += 1
#                      
#    plot_data.addLine(temperature,
#                      b2_list,
#                      xlabel=xlabel,
#                      ylabel=ylabel,
#                      linelabel='$b_2^k$',
#                      linewidth=2,
#                      linestyle='-',
#                      marker=marker_list[i],
#                      markersize=12,
#                      color=color_list[i])
#    i += 1                  
                          
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
    plt.xlim(400,1000)
#    plt.ylim(-150,0)
    plt.ylim(0,300)
#    plt.ylim(0.2,0.6)
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
#    plot_data.lines[0]['ylabel'] = '$b_1^k,b_2^k$'
#    plot_data.lines[1]['ylabel'] = '$b_1^k,b_2^k$'
    plot_data.plot()
#==============================================================================
# http://stackoverflow.com/questions/21920233/matplotlib-log-scale-tick-label-number-formatting
#==============================================================================
#    ax.xaxis.set_major_locator(MultipleLocator(0.5))
#    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
#    ax.xaxis.set_major_formatter(ScalarFormatter())
#    ax.yaxis.set_major_locator(MultipleLocator(50))
#    ax.yaxis.set_minor_locator(MultipleLocator(10))
#    ax.yaxis.set_major_formatter(ScalarFormatter())
#==============================================================================
# show legend
#==============================================================================
    lg = plt.legend(title='',loc=2)
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
figure_name='plot_umat_parameters_r0k'
#figure_name='plot_umat_parameters_rdeltak'
#figure_name='plot_umat_parameters_a1'
#figure_name='plot_umat_parameters_b1b2'
create_plot_data(figure_path,figure_name)
plot(figure_path,figure_name,save_types=['.pdf'])

shutil.copy(__file__,ArticleFigureDirectory)