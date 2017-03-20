# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 13:46:18 2017

@author: SunJingyu
"""

import numpy as np
import matplotlib.pyplot as plt
import re
import os
from read_exp_files import read_exp_files
from Functions import (obtain_masing_curve,obtain_youngs_modulus,
                      obtain_kinematic_hardening_parameters,linear_interpolation,
                      osgood_interpolation,obtain_plastic_strain_by_stress)
from Constants import *

#==============================================================================
# ExperimentData
#==============================================================================
class ExperimentData:
    """
    Read and analyze experiment data.
    """
    def __init__(self, filename):
        """
        读取csv文件，表头和单位如下：
        
        ['Axial Segment Count' 'Running Time' 'Temperature' 'Axial Displacement'
         'Axial Force' 'Axial Strain' 'Axial Stress' 'Rotation' 'Torque'
         'Angle Strain' 'Shear Stress' 'Equivalent Plastic Strain' 'Thermal Strain'
         'Axial Total Strain']
         
        ['cycles' 'sec' 'C' 'mm' 'N' 'mm/mm' 'Mpa' 'deg' 'N*m' '-' 'Mpa' '-'
         'mm/mm' 'mm/mm']
         
        文件大于200MB的时候，32位系统中，np.genfromtxt()函数会造成内存溢出。
        这时候使用自定义的read_exp_files()函数
        """
        if os.path.getsize(filename) < 150*1024*1024: # 如果文件小于150MB
            data = np.genfromtxt(filename, delimiter=',', skip_header=2, dtype=float)
        else:
            data = read_exp_files(filename, n=2)
            
        self.filename = filename
        self.nodelabel = [36] # 试验对应的有限元模型节点为36号
        
        self.axial_count = data[:,0]
        self.runing_time = data[:,1]
        self.temperature = data[:,2]
        self.axial_disp = data[:,3]
        self.axial_force = data[:,4]
        self.axial_strain = data[:,5]
        self.axial_stress = data[:,6]
        self.rotation = data[:,7]
        self.torque = data[:,8]
        self.shear_strain = data[:,9]
        self.shear_stress = data[:,10]
        self.eqpl_strain = data[:,11]
        self.thermal_strain = data[:,12]
        self.total_strain = data[:,13]
        del data
        
        self.shear_stress_eq = self.shear_stress*np.sqrt(3.0)
        self.shear_strain_eq = self.shear_strain/np.sqrt(3.0)
        self.axial_log_strain = np.log(1.0+self.axial_strain)
        self.axial_true_stress = self.axial_stress*(1.0+self.axial_strain)

        # 数组长度
        self.length = len(self.axial_count)
        # 总试验时间
        self.total_runing_time = self.runing_time[-1]
        # 总循环数
        self.total_axial_count = int(self.axial_count[-1])
        # 初始半疲劳寿命
        self.half_life_cycle = 0
        # 初始温度
        self.initial_temperature = self.temperature[0]
        
        self.axial_count_index_list = []
        self.axial_count_begin_index = {}
        self.axial_count_end_index = {}
        
        self.obtainCountIndex()
        
    def obtainCountIndex(self):
        """
        建立每个循环开始位置的索引表：axial_count_begin_index。
        建立每个循环结束位置的索引表：axial_count_end_index。
        计算半疲劳寿命：half_life_cycle。
        """
        self.axial_count_begin_index = {}
        self.axial_count_end_index = {}
        self.axial_count_begin_index[0] = 0
        self.axial_count_begin_index[int(self.axial_count[0])] = 0
        for i in range(self.length-1):
            if int(self.axial_count[i]) <> int(self.axial_count[i+1]): # 如果当前行计数与下一行计数不同
                self.axial_count_end_index[int(self.axial_count[i])] = i # 当前循环结尾为i-1,这里用i,因为[a:b]中截取a到b-1
                self.axial_count_begin_index[int(self.axial_count[i+1])] = i # 下一循环的开头为i
        self.axial_count_index_list = list(self.axial_count_begin_index.keys()) # 所有存在的循环列表
        for c in self.axial_count_index_list:
            if abs(c-int(self.total_axial_count/2.0))<=abs(self.half_life_cycle-int(self.total_axial_count/2.0)):
                self.half_life_cycle = c
                
    def obtainNthCycle(self,item,begin_cycle,end_cycle=None):
        """
        截取指定数组的第begin_cycle到end_cycle循环中的数据。
        """
        if end_cycle is None:
            end_cycle = begin_cycle
        if self.total_axial_count == []:
            print 'self.total_axial_count is empty.'
            return []
        if self.axial_count_begin_index.has_key(begin_cycle) and self.axial_count_end_index.has_key(end_cycle):
            return eval('self.'+item)[self.axial_count_begin_index[begin_cycle]:self.axial_count_end_index[end_cycle]]
        else:
            return []
            
    def obtainPeakValley(self,item):
        """
        返回指定数组的峰谷值。
        """
        cycle = []
        peak = []
        valley = []
        for n in self.axial_count_index_list:
            data_nth_cycle = self.obtainNthCycle(item,n)
            if data_nth_cycle <> []:
                maxitem = max(data_nth_cycle)
                minitem = min(data_nth_cycle)
                cycle.append(n)
                peak.append(maxitem)
                valley.append(minitem)
        return cycle,peak,valley
        
#==============================================================================
# SimulationData
#==============================================================================
class SimulationData:
    """
    Read and analyze simulation data.
    """
    def __init__(self, filename, period):
        data = np.genfromtxt(filename, delimiter=',', skip_header=0, dtype=None)
        self.header = [i for i in data[0]]
        del data
        
        data = np.genfromtxt(filename, delimiter=',', skip_header=1, dtype=float)
        self.filename = filename
        self.axial_count=[]
        self.frame = []
        self.runing_time = []
        self.node_label = []
        self.temperature = []
        self.axial_strain = []
        self.axial_stress = []
        self.shear_strain = []
        self.shear_stress = []
        self.mises_stress = []
        self.axial_count_index = {}
        self.axial_count_index_list = []
        self.total_axial_count = []
        self.item = ['axial_count', 'axial_strain', 'axial_stress', 'frame', 'mises_stress', 'node_label', 'runing_time', 'shear_strain', 'shear_stress', 'temperature']
        
        for h in self.header:
            if h in ['Frame']:
                self.frame = data[:,self.header.index(h)]
            if h in ['Time']:
                self.runing_time = data[:,self.header.index(h)]
            if h in ['NodeLabel']:
                self.node_label = data[:,self.header.index(h)]
            if h in ['Temperature']:
                self.temperature = data[:,self.header.index(h)]
            if h in ['LE22','E22','E33']:
                self.axial_strain = data[:,self.header.index(h)]
            if h in ['S22','S33']:
                self.axial_stress = data[:,self.header.index(h)]
            if h in ['LE23','E23']:
                self.shear_strain = data[:,self.header.index(h)]
            if h in ['S23']:
                self.shear_stress = data[:,self.header.index(h)]
            if h in ['Mises']:
                self.mises_stress = data[:,self.header.index(h)]
        del data
        
        self.length = len(self.runing_time)
        self.axial_count = np.zeros(self.length)
        
        self.obtainCountIndex(period)
        
    def obtainCountIndex(self,period):
        """
        建立每个循环开始位置的索引表：axial_count_index。
        """
        self.total_axial_count = int(self.runing_time[-1]/period)
        self.axial_count_index = {}
        self.axial_count_index[0] = 0
        self.axial_count_index[1] = 0
        count = 1
        for i in range(self.length):
            if self.runing_time[i] <= count * period and self.runing_time[i] >= (count-1) * period:
                self.axial_count[i] = count
            if self.runing_time[i] > count * period:
                count = int(self.runing_time[i]/period) + 1
                self.axial_count_index[count] = i
        self.axial_count_index_list = list(self.axial_count_index.keys())

    def obtainNthCycle(self,item,begin_cycle,end_cycle=None):
        """
        截取指定数组的第begin_cycle到end_cycle循环中的数据。
        """
        if end_cycle is None:
            end_cycle = begin_cycle
        if self.total_axial_count == []:
            print 'self.total_axial_count is empty.'
        if self.axial_count_index.has_key(begin_cycle) and self.axial_count_index.has_key(end_cycle+1):
            return eval('self.'+item)[self.axial_count_index[begin_cycle]:self.axial_count_index[end_cycle+1]]
            
    def obtainPeakValley(self,item):
        """
        返回指定数组的峰谷值。
        """
        cycle = []
        peak = []
        valley = []
        for n in range(1,self.total_axial_count+1):
            data_nth_cycle = self.obtainNthCycle(item,n)
            if data_nth_cycle <> []:
                maxitem = max(data_nth_cycle)
                minitem = min(data_nth_cycle)
                cycle.append(n)
                peak.append(maxitem)
                valley.append(minitem)
        return cycle,peak,valley
        
#==============================================================================
# ExperimentLog
#==============================================================================
class ExperimentLog:
    """
    Read experiment log.
    """
    def __init__(self, filename):
        data = np.genfromtxt(filename, delimiter=',', skip_header=0, dtype=None)
        self.header = [i for i in data[0]]
        del data

        data = np.genfromtxt(filename, delimiter=',', skip_header=1, dtype=str)
        """
        ['number', 'din mm', 'dout mm', 'gauge length mm', 'load type', 'axial control mode', 
        'rotational control mode', 'temperature mode', 'axial displacement', 'axial strain %', 
        'axial force', 'rotation', 'angel strain deg', 'torque', 'equivalent strain %', 
        'axial rotational phase deg', 'axial temperature phase deg', 'Period', 
        'load rate or frequence', 'test date', 'comments', 'calculate']
        """

        self.number = data[:,0]
        self.d_in = data[:,1]
        self.d_out = data[:,2]
        self.gauge_length = data[:,3]
        self.load_type = data[:,4]
        self.axial_control_mode = data[:,5]
        self.rotational_control_mode = data[:,6]
        self.temperature_mode = data[:,7]
        self.axial_displacement = data[:,8]
        self.axial_strain = data[:,9]
        self.axial_force = data[:,10]
        self.rotation = data[:,11]
        self.angel_strain = data[:,12]
        self.torque = data[:,13]
        self.equivalent_strain = data[:,14]
        self.axial_rotational_phase = data[:,15]
        self.axial_temperature_phase = data[:,16]
        self.period = data[:,17]
        self.load_rate = data[:,18]
        self.test_date = data[:,19]
        self.comments = data[:,20]
        self.calculate = data[:,21]
        self.data = data
        self.length = len(self.number)

    def output(self, name):
        number = self.find(name)
        print '=========================Experiment========================='
        for i in range(len(self.header)):
            print '%-40s%-20s' % (self.header[i],self.data[number][i])
        
    def find(self, name):
        number = int(np.where(self.number == name)[0][0])
        return number
        
    def keyFilter(self, key):
        numbers = np.where(eval(key))[0]
        for number in numbers:
            self.find(self.number[number])
        number_list = []
        for number in numbers:
            number_list.append(self.number[number])
        return number_list

    def obtainItem(self, name, item, regular):
        """
        总结
        ^ 匹配字符串的开始。
        $ 匹配字符串的结尾。
        \b 匹配一个单词的边界。
        \d 匹配任意数字。
        \D 匹配任意非数字字符。
        x? 匹配一个可选的 x 字符 (换言之，它匹配 1 次或者 0 次 x 字符)。
        x* 匹配0次或者多次 x 字符。
        x+ 匹配1次或者多次 x 字符。
        x{n,m} 匹配 x 字符，至少 n 次，至多 m 次。
        (a|b|c) 要么匹配 a，要么匹配 b，要么匹配 c。
        (x) 一般情况下表示一个记忆组 (remembered group)。你可以利用 re.search 函数返回对象的 groups() 函数获取它的值。
        正则表达式中的点号通常意味着 “匹配任意单字符”
        """        
        number = self.find(name)
        string = eval('self.%s[number]' % item)
        result_list = re.findall(regular, string)
        if len(result_list) == 0:
            return [0.0]
        elif len(result_list) >= 1:
            return result_list

#==============================================================================
# FatigueData
#==============================================================================
class FatigueData:
    """
    Read and analyze fatigue data.
    """
    def __init__(self, filename):
        """
        ['Number of Cycles to Failure N\\-(f)', 'Mises Equivalent Strain Amplitude \\i(\\g(De))\\-(eq)/2', 
        'Stress Amplitude e \\i(\\g(Ds))/2', 'Specimen', 'Critical Plane', 'sigma_n_max', 'delta_sigma', 
        'delta_epsilon', 'tau_n_max', 'delta_tau', 'delta_gamma', 'Predicted Fatigue Lifetime N\\-(p)', 
        'Fatigue Coefficient', 'Temperature', 'Load Type']
        """
        data = np.genfromtxt(filename, delimiter=',', skip_header=0, dtype=None)
        self.header = [i for i in data[0]]
        self.unit = [i for i in data[1]]
        del data
        
        data = np.genfromtxt(filename, delimiter=',', skip_header=2, dtype=float)
        self.experimental_life = data[:,0]
        self.equivalent_strain_amplitude = data[:,1]
        self.stress_amplitude = data[:,2]
        self.specimen = data[:,3]
        self.critical_plane = data[:,4]
        self.sigma_n_max = data[:,5]
        self.delta_sigma = data[:,6]
        self.delta_epsilon = data[:,7]
        self.tau_n_max = data[:,8]
        self.delta_tau = data[:,9]
        self.delta_gamma = data[:,10]
        self.predicted_life = data[:,11]
        self.fatigue_coefficient = data[:,12]
        self.temperature = data[:,13]
        del data
        
        data = np.genfromtxt(filename, delimiter=',', skip_header=2, dtype=str)
        self.load_type = data[:,14]
        del data

    def loadTypeFilter(self,load_type,item):
        numbers = np.where(self.load_type == load_type)[0]
        result = []
        for number in numbers:
            result.append(eval('self.%s[number]' % item))
        return result
        
#==============================================================================
# PlotData
#==============================================================================
class PlotData:
    """
    Define a data structure for plot.
    
    lines = [{line1},{line2},...]
    
    Example:
    
    life_NASA = [33,35,200,515,830,1072,1785,1850,2868,4323,8903,13332,13762]
    strain_amplitude_NASA = [1.75,1.75,0.85,0.85,0.575,0.5,0.507,0.575,0.398,0.402,0.329,0.323,0.328] # %
    stess_amplitude_NASA = [873.5,873.5,726.6,726.6,636.8,601,605,636.8,538,541,480,474,479] # MPa
    
    life_BHU = [100000,48745,50900,35940,12980,6300,2220,678,515]
    strain_amplitude_BHU = [0.4,0.42,0.43,0.45,0.47,0.5,0.6,0.8,1.0] # %
    stess_amplitude_BHU = [621,638,607,606,653,646,749,756,779] # MPa
    
    plot_data = PlotData()
    plot_data.addLine(life_NASA,strain_amplitude_NASA)
    plot_data.addLine(life_BHU,strain_amplitude_BHU)
    """
    def __init__(self):
        self.lines = []
        self.number_of_header_lines = 8
        
    def addLine(self,
                x,
                y,
                xlabel='xlabel',
                ylabel='ylabel',
                linelabel='',
                linestyle='',
                linewidth=2,
                marker=None,
                markersize=12,
                color='auto'):
        
        self.lines.append({'x':x,
                           'y':y,
                           'xlabel':xlabel,
                           'ylabel':ylabel,
                           'linelabel':linelabel,
                           'linestyle':linestyle,
                           'linewidth':linewidth,
                           'marker':marker,
                           'markersize':markersize,
                           'color':color})
        
    def plot(self):
        for line in self.lines:
            plt.xlabel(line['xlabel'])
            plt.ylabel(line['ylabel'])
            if line['color'] == 'auto':
                plt.plot(line['x'],
                         line['y'],
                         label=line['linelabel'],
                         linestyle=line['linestyle'],
                         linewidth=line['linewidth'],
                         marker=line['marker'],
                         markersize=line['markersize'])
            else:
                plt.plot(line['x'],
                         line['y'],
                         label=line['linelabel'],
                         linestyle=line['linestyle'],
                         linewidth=line['linewidth'],
                         marker=line['marker'],
                         markersize=line['markersize'],
                         color=line['color'])
            
    def writeToFile(self,directory,filename):
        resultfile = open(directory + filename + '.csv', 'w') # write to csv
        data_list = []
        for line in self.lines:
            xarray = []
            xarray.append('x')
            xarray.append(line['xlabel'])
            xarray.append(line['linelabel'])
            xarray.append(line['linestyle'])
            xarray.append(line['linewidth'])
            xarray.append(line['marker'])
            xarray.append(line['markersize'])
            xarray.append(line['color'])
            for x in line['x']:
                xarray.append(x)
            data_list.append(xarray)
            
            yarray = []
            yarray.append('y')
            yarray.append(line['ylabel'])
            yarray.append(line['linelabel'])
            yarray.append(line['linestyle'])
            yarray.append(line['linewidth'])
            yarray.append(line['marker'])
            yarray.append(line['markersize'])
            yarray.append(line['color'])
            for y in line['y']:
                yarray.append(y)
            data_list.append(yarray)
        
        data_length_list = []
        for data in data_list:
            data_length_list.append(len(data))
        data_length_max = max(data_length_list)
        data_height = len(data_list)

        for i in range(data_length_max):
            l = ''
            for j in range(data_height):
                if i < data_length_list[j]:
                    l += '%s,' % str(data_list[j][i])
                else:
                    l += '%s,' % ''
            print >>resultfile, l[:-1]
            
        resultfile.close()
        print 'save as', directory + filename + '.csv'
            
    def readFromFile(self,directory,filename):
        self.lines = []
        data = np.genfromtxt(directory + filename + '.csv', delimiter=',', skip_header=0, dtype=str) # read from csv
        data = data.transpose()
        data_item_list = []
        for d in data:
            data_item_list.append(d[0])
        data_item = np.array(data_item_list)
        data_item_index = np.where(data_item == 'x')[0]
        for i in data_item_index:
            x_index = i
            y_index = i + 1
            x_data = data[x_index]
            y_data = data[y_index]
            x = [float(x) for x in x_data[self.number_of_header_lines:] if x <> '']
            y = [float(y) for y in y_data[self.number_of_header_lines:] if y <> '']
            xlabel = x_data[1]
            ylabel = y_data[1]
            linelabel = x_data[2]
            linestyle = x_data[3]
            linewidth = float(x_data[4])
            marker = x_data[5]
            markersize = int(x_data[6])
            color = x_data[7]
            self.lines.append({'x':x,
                               'y':y,
                               'xlabel':xlabel,
                               'ylabel':ylabel,
                               'linelabel':linelabel,
                               'linestyle':linestyle,
                               'linewidth':linewidth,
                               'marker':marker,
                               'markersize':markersize,
                               'color':color})