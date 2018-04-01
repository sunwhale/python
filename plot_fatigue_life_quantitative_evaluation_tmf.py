# -*- coding: utf-8 -*-
"""
Created on Thu Jan 05 11:50:47 2017

@author: j.Sun
"""

import matplotlib.pyplot as plt
import matplotlib as mpl
import shutil
import json

from Data import FatigueData,PlotData
from Constants import xylabels,FatigueDirectory,ArticleFigureDirectory
from matplotlib.ticker import MultipleLocator,ScalarFormatter,FormatStrFormatter 
from plot_format import plot_format


"""
## Writing JSON data
with open(r'F:\Temp\test.txt', 'w') as f:
    json.dump(data, f)

## Reading data back
with open(r'F:\Temp\test.txt', 'r') as f:
    data = json.load(f)
"""

def read_file(file_name):
    """
    Writing JSON data to file.
    """
    with open(file_name,'r') as data_file:
        return json.loads(data_file.read())

def write_file(file_name, data):
    """
    Reading JSON data to file.
    """
    with open(file_name,'w') as data_file:
        return json.dump(data, data_file)
    
def saveFigure(figure_path=None,figure_name=None,save_types=[]):
    if figure_path <> None and figure_name<> None:
        for save_type in save_types:
            plt.savefig(figure_path + figure_name + save_type, dpi=150, transparent=True)
            print 'save as', figure_path + figure_name + save_type
                
#fatigue_model_list = ['BM','FS','SWT','Liu1','Liu2','Chu']
fatigue_model_list = ['BM','FS','SWT','Chu','Liu1','Study']
#fatigue_model_list = ['BM']

TN_list = []
TRMS_list = []

for fatigue_model in fatigue_model_list:
    fatigue_file = '%s%s.csv' % (FatigueDirectory,fatigue_model)
    fatigue_data = FatigueData(fatigue_file)
    TN_list.append(fatigue_data.quantitativeTN())
    TRMS_list.append(fatigue_data.quantitativeTRMS())

write_file(ArticleFigureDirectory + 'TN_list.txt', TN_list)
write_file(ArticleFigureDirectory + 'TRMS_list.txt', TRMS_list)

TN_list = read_file(ArticleFigureDirectory + 'TN_list.txt')
TRMS_list = read_file(ArticleFigureDirectory + 'TRMS_list.txt')

plot_format()
mpl.rcParams['xtick.major.size'] = 0
mpl.rcParams['xtick.major.width'] = 0
#==============================================================================
# x,y limite
#==============================================================================
plt.xlim(-0.5,6.5)
plt.ylim(0,2)
#==============================================================================
# http://blog.csdn.net/lanchunhui/article/details/52931883 关闭坐标刻度
#==============================================================================
plt.xticks([])
#plt.yticks([])
plt.xlabel('')
plt.ylabel('')
#==============================================================================
# xy log scale
#==============================================================================
#plt.xscale('log')
#plt.yscale('log')
#==============================================================================
# http://stackoverflow.com/questions/21920233/matplotlib-log-scale-tick-label-number-formatting
#==============================================================================
ax = plt.gca()
#ax.xaxis.set_major_locator(MultipleLocator(100))
#ax.xaxis.set_minor_locator(MultipleLocator(10))
#ax.xaxis.set_major_formatter(ScalarFormatter())
#ax.yaxis.set_major_locator(MultipleLocator(1))
#ax.yaxis.set_minor_locator(MultipleLocator(1))
#ax.yaxis.set_major_formatter(ScalarFormatter())
fatigue_model_list = ['BM','FS','SWT','Chu','Liu I','Study']
#==============================================================================
# annotate
#==============================================================================
plt.text(0.0,1.3,r'conservative',fontsize=20,color='red')
plt.text(0.0,0.7,r'non-conservative',fontsize=20,color='red')

n = range(len(TN_list))
plt.plot([-10,10],[1,1],ls='--',lw='4',color='red')
plt.bar(n,TN_list)
    
#n = range(len(TRMS_list))
#plt.bar(n,TRMS_list)

for t in TN_list:
    print t
    
for t in TRMS_list:
    print t
    
n = [i+0.5 for i in n]
plt.xticks(n, fatigue_model_list)

figure_path = ArticleFigureDirectory
#figure_name = 'plot_fatigue_life_quantitative_evaluation_tmf_TRMS'
#figure_name = 'plot_fatigue_life_quantitative_evaluation_tmf_TN'
#saveFigure(figure_path,figure_name,save_types=['.pdf','.png'])
#    
#shutil.copy(__file__,ArticleFigureDirectory)
#
#plt.show()