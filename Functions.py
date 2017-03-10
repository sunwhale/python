# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 13:46:18 2017

@author: SunJingyu
"""

import numpy as np
import matplotlib.pyplot as plt
import re
import os
import shutil
from scipy.optimize import leastsq

#==============================================================================
# some constants
#==============================================================================
plastic_strain_list=[0.00005,0.0001,0.0002,0.0005,0.001,0.002,0.004,0.01,0.02]
plastic_strain_list=[0.001,0.002,0.003,0.004,0.005,0.01,0.02]
#==============================================================================
# get_files
#==============================================================================
def get_files(path):
    fns=[]
    for root, dirs, files in os.walk( path ):
        for fn in files:
            fns.append( [ root , fn ] )
    return fns
#==============================================================================
# is_text_file
#==============================================================================
def is_text_file(filename):
    b = False
    suffixs = ['txt','dat','log']
    for suffix in suffixs:
        if filename.split('.')[-1] == suffix:
            b = True
    return b
#==============================================================================
# is_csv_file
#==============================================================================
def is_csv_file(filename):
    b = False
    suffixs = ['csv']
    for suffix in suffixs:
        if filename.split('.')[-1] == suffix:
            b = True
    return b
    
AbaqusTempDirectory = 'F:\\Temp\\IN7183\\'
AbaqusOutputDirectory = 'F:\\Temp\\IN718_Sim\\'
if not os.path.isdir(AbaqusOutputDirectory):
    os.makedirs(AbaqusOutputDirectory)
    print 'Create new directory:',AbaqusOutputDirectory
            
filenames = get_files(AbaqusTempDirectory)
for filename in filenames:
    if is_csv_file(filename[1]):
        print filename
        fullname = filename[0] + '\\' + filename[1]
        shutil.copy(fullname, AbaqusOutputDirectory)
#==============================================================================
# obtain_masing_curve
#==============================================================================
def obtain_masing_curve(strain_stress_curve):
    """
    输入应力应变关系[strain,stress]，输出[masing_strain, masing_stress]
    """
    strain = strain_stress_curve[0]
    stress = strain_stress_curve[1]
    
    if len(strain) == len(stress):
        strain_list = [i for i in strain]
        stress_list = [i for i in stress]
        
        stress_max_index = stress_list.index(max(stress_list))
        stress_min_index = stress_list.index(min(stress_list))
        
        strain_max_index = strain_list.index(max(strain_list))
        strain_min_index = strain_list.index(min(strain_list))
        
#        if stress_min_index<stress_max_index:
#            tmp_strain = np.array(strain_list[stress_min_index:stress_max_index])
#            tmp_stress = np.array(stress_list[stress_min_index:stress_max_index])
#        else:
#            tmp_strain = np.array(strain_list[stress_min_index:]+strain_list[:stress_max_index])
#            tmp_stress = np.array(stress_list[stress_min_index:]+stress_list[:stress_max_index])
        
        if strain_min_index<strain_max_index:
            tmp_strain = np.array(strain_list[strain_min_index:strain_max_index])
            tmp_stress = np.array(stress_list[strain_min_index:strain_max_index])
        else:
            tmp_strain = np.array(strain_list[strain_min_index:]+strain_list[:strain_max_index])
            tmp_stress = np.array(stress_list[strain_min_index:]+stress_list[:strain_max_index])
            
        masing_strain = (tmp_strain + abs(tmp_strain[0]))/2
        masing_stress = (tmp_stress + abs(tmp_stress[0]))/2
        
        return [masing_strain, masing_stress]
    else:
        return []
#==============================================================================
# obtain_youngs_modulus
#==============================================================================
def obtain_youngs_modulus(strain_stress_curve,elastic_limit=0.002):
    """
    根据应力应变曲线计算杨氏模量，设置应变范围大于零且小于elastic_limit的数据进行拟合。
    输入应力应变关系[strain,stress]，输出youngs_modulus。
    """
    strain = strain_stress_curve[0]
    stress = strain_stress_curve[1]
    if len(strain) == len(stress):
        strain_list = []
        stress_list = []
        for i in range(len(strain)):
            if strain[i]>=0.0 and strain[i]<=elastic_limit:
                strain_list.append(strain[i])
                stress_list.append(stress[i])
        x = np.array(strain_list)
        y = np.array(stress_list)
        z = np.polyfit(x,y,1)
        p = np.poly1d(z)
        return z[0]
    else:
        return 0.0
#==============================================================================
# linear_interpolation
#==============================================================================
def linear_interpolation(xy,position_x_list=plastic_strain_list):
    """
    线性差值。
    """
    position_y_list = []
    for position_x in position_x_list:
        xx=xy[0]
        yy=xy[1]
        for i in range(len(xx)):
            if xx[i]>position_x:
                lower_index=i-1
                upper_index=i
                break
        x=np.array([xx[lower_index],xx[upper_index]])
        y=np.array([yy[lower_index],yy[upper_index]])
        z=np.polyfit(x,y,1)
        p = np.poly1d(z)
        position_y_list.append(p(position_x))
    if len(position_x_list) == len(position_y_list):
        return [np.array(position_x_list),np.array(position_y_list)]
    else:
        print 'linear interpolation error.'
        return []
#==============================================================================
# osgood_interpolation
#==============================================================================
def osgood_interpolation(xy,position_x_list=plastic_strain_list):
    """
    指数差值。
    """
    def func(x,p):
        a,b = p
        return a*x**b
    def residuals(p,y,x):
        return y-func(x,p)
        
    position_y_list = []
    xx=xy[0]
    yy=xy[1]
    x=[]
    y=[]
    for i in range(len(xx)):
        if xx[i]>0.001 and xx[i]<0.03: # 300C 980, 550C 850, 650C 780   
            x.append(xx[i])
            y.append(yy[i])
    p0 = [1000,1]
    plsq = leastsq(residuals, p0, args=(y, x))
    
    for position_x in position_x_list:
        position_y_list.append(func(position_x,plsq[0]))
        
    if len(position_x_list) == len(position_y_list):
        return [np.array(position_x_list),np.array(position_y_list)]
    else:
        print 'osgood interpolation error.'
        return []
#==============================================================================
# 
#==============================================================================
def obtain_plastic_strain_by_stress(strain_stress_curve,youngs_modulus):
    strain = strain_stress_curve[0]
    stress = strain_stress_curve[1]
    plastic_strain = strain - stress/youngs_modulus
    return [plastic_strain,stress]
#==============================================================================
# obtain_kinematic_hardening_parameters
#==============================================================================
def obtain_kinematic_hardening_parameters(monotonic_curve,cyclic_curve,
                                          position_x_list=plastic_strain_list,yield_stress=300.0,
                                          youngs_modulus=180000.0,show=False):
    """
    给入应变分区、单拉曲线和循环滞回曲线，计算UMAT模型所需要的参数。
    """    
    position_x = np.array(position_x_list)
    zeta = 1.0/position_x
    seg = len(position_x_list)
    r0 = []
    ri = []
    
    m = osgood_interpolation(obtain_plastic_strain_by_stress(monotonic_curve,youngs_modulus),position_x_list)
    H = []
    H.append( (m[1][0]-yield_stress)/(position_x[0]-0.0) )
    for i in range(1,seg):
        H.append( (m[1][i]-m[1][i-1])/(position_x[i]-position_x[i-1]) )
    for i in range(0,seg-1):
        r0.append( (H[i]-H[i+1])/zeta[i] )
    r0.append( (H[-1]-0.0)/zeta[i] )
    
    cyclic_marsing_curve = obtain_masing_curve(cyclic_curve)
    c = osgood_interpolation(obtain_plastic_strain_by_stress(cyclic_marsing_curve,youngs_modulus),position_x_list)
    H = []
    H.append( (c[1][0]-yield_stress)/(position_x[0]-0.0) )
    for i in range(1,seg):
        H.append( (c[1][i]-c[1][i-1])/(position_x[i]-position_x[i-1]) )
    for i in range(0,seg-1):
        ri.append( (H[i]-H[i+1])/zeta[i] - r0[i] )
    ri.append( (H[-1]-0.0)/zeta[i] - r0[i] )
    
    """
    绘制差值对比图。
    """
    
    if show:
        monotonic_curve = obtain_plastic_strain_by_stress(monotonic_curve,youngs_modulus)
        plt.plot(monotonic_curve[0],monotonic_curve[1])
        monotonic_curve = osgood_interpolation(monotonic_curve)
        plt.plot(monotonic_curve[0],monotonic_curve[1])
        
        cyclic_curve = obtain_masing_curve(cyclic_curve)
        cyclic_curve = obtain_plastic_strain_by_stress(cyclic_curve,youngs_modulus)
        plt.plot(cyclic_curve[0],cyclic_curve[1])
        cyclic_curve = osgood_interpolation(cyclic_curve)
        plt.plot(cyclic_curve[0],cyclic_curve[1])
        plt.show()

    return seg,zeta,r0,ri
#==============================================================================
# calculate_elastic_by_temperature_in718
#==============================================================================
def calculate_elastic_by_temperature_in718(temperature):
    TEMP = temperature
    EMOD = 206308.7426+(-51.20306)*TEMP+0.01109*TEMP**2+(-3.84391E-05)*TEMP**3
    ENU = 2.901300E-01+(1.457750E-05)*TEMP+(-2.067420E-07)*TEMP**2+(2.780300E-10)*TEMP**3
    EG = 0.5*EMOD/(1.0+ENU)
    youngs_modulus = EMOD
    poisson_ratio = ENU
    shear_modulus = EG
    return youngs_modulus,poisson_ratio,shear_modulus
#==============================================================================
# write_umat_output_file
#==============================================================================
def write_umat_output_file(fortran_file_name,text_file_name,total_time):
    outfile = open(fortran_file_name, 'w')
    print >>outfile, """
      OPEN(UNIT=101, FILE='%s')
200   FORMAT('=========================',F8.1,'/',F8.1,'==========',F6.2,'%%')
      TOTALTIME=%s
      IF ( (NOEL .EQ. 25) .AND. (NPT .EQ. 1) ) THEN
      WRITE(101,200) TIME(2),TOTALTIME,TIME(2)/TOTALTIME*100.0
      WRITE(101,*) 'PACC',PACC
      WRITE(101,*) 'TEMP',TEMP
      WRITE(101,*) 'TRIAXIALITY',TRIAXIALITY
      WRITE(101,*) 'YD',YD
C       IF ( MOD(INT(TIME(2)/TOTALTIME*100.0),1) .EQ. 0 ) THEN
      WRITE(*,200) TIME(2),TOTALTIME,TIME(2)/TOTALTIME*100.0
      WRITE(*,*) 'PACC',PACC
      WRITE(*,*) 'TEMP',TEMP
      WRITE(*,*) 'TRIAXIALITY',TRIAXIALITY
      WRITE(*,*) 'YD',YD
C       ENDIF
      ENDIF""" % (text_file_name,total_time)
    outfile.close()