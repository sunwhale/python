# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 13:46:18 2017

@author: SunJingyu
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from Work import UMAT
from Data import ExperimentData,ExperimentLog
from Functions import obtain_kinematic_hardening_parameters,calculate_elastic_by_temperature_in718
from Constants import *
from workbench import workbench
from compare_exp_sim import compare_exp_sim

def calculate_umat_parameters_in718(name='_output.txt'):
    parameters = []
    yield_stress = 500.0
    show = False
#    plastic_strain_list=[0.001,0.002,0.005,0.01,0.02]
#    plastic_strain_list=[0.002,0.005,0.01,0.02,0.05]
    plastic_strain_list=[0.00005,0.0001,0.0002,0.0005,0.001,0.002,0.004,0.01,0.02]
    umat = UMAT(UMATDirectory = UMATDirectory, 
                UMATMainFile = 'MAIN_IN718.for',
                ParameterFortranFile = 'PARAMETERS_IN718_TMF.for',
                OutputFortranFile = 'OUTPUT.for',
                OutputTextFile = name + '_output.txt')
    outfile = open(umat.ParameterFortranFileFullName, 'w')
    print outfile.name
#==============================================================================
# 300C
#==============================================================================
    temperature = 300.0
#    youngs_modulus, poisson_ratio, shear_modulus = calculate_elastic_by_temperature_in718(temperature)
    youngs_modulus, poisson_ratio, shear_modulus, y = calculate_elastic_by_temperature_in718(temperature)
    
    name = '7101'
    exp_full_name = ExperimentDirectory + name + '.csv'
    exp_monotonic = ExperimentData(exp_full_name)
    
    name = '7211'
    exp_full_name = ExperimentDirectory + name + '.csv'
    exp_cyclic = ExperimentData(exp_full_name)
    
    monotonic_curve = [exp_monotonic.axial_log_strain,exp_monotonic.axial_true_stress]
    cyclic_curve = [exp_cyclic.obtainNthCycle('axial_strain',190),exp_cyclic.obtainNthCycle('axial_stress',190)]
    seg,zeta,r0,ri = obtain_kinematic_hardening_parameters(monotonic_curve,cyclic_curve,
                                                           position_x_list=plastic_strain_list,
                                                           yield_stress=yield_stress,
                                                           youngs_modulus=youngs_modulus,show=show)
    parameters.append([temperature,seg,zeta,r0,ri])
#==============================================================================
# 550C
#==============================================================================
    temperature = 550.0
#    youngs_modulus, poisson_ratio, shear_modulus = calculate_elastic_by_temperature_in718(temperature)
    youngs_modulus, poisson_ratio, shear_modulus, y = calculate_elastic_by_temperature_in718(temperature)
    
    name = '7103'
    exp_full_name = ExperimentDirectory + name + '.csv'
    exp_monotonic = ExperimentData(exp_full_name)
    
    name = '7212'
    exp_full_name = ExperimentDirectory + name + '.csv'
    exp_cyclic = ExperimentData(exp_full_name)
    
    monotonic_curve = [exp_monotonic.axial_log_strain,exp_monotonic.axial_true_stress]
    cyclic_curve = [exp_cyclic.obtainNthCycle('axial_strain',190),exp_cyclic.obtainNthCycle('axial_stress',190)]
    seg,zeta,r0,ri = obtain_kinematic_hardening_parameters(monotonic_curve,cyclic_curve,
                                                           position_x_list=plastic_strain_list,
                                                           yield_stress=yield_stress,
                                                           youngs_modulus=youngs_modulus,show=show)
    parameters.append([temperature,seg,zeta,r0,ri])
#==============================================================================
# 650C
#==============================================================================
    temperature = 650.0
#    youngs_modulus, poisson_ratio, shear_modulus = calculate_elastic_by_temperature_in718(temperature)
    youngs_modulus, poisson_ratio, shear_modulus, y = calculate_elastic_by_temperature_in718(temperature)
    
    name = '7102'
    exp_full_name = ExperimentDirectory + name + '.csv'
    exp_monotonic = ExperimentData(exp_full_name)
    
    name = '7213'
    exp_full_name = ExperimentDirectory + name + '.csv'
    exp_cyclic = ExperimentData(exp_full_name)
    
    monotonic_curve = [exp_monotonic.axial_log_strain,exp_monotonic.axial_true_stress]
    cyclic_curve = [exp_cyclic.obtainNthCycle('axial_strain',190),exp_cyclic.obtainNthCycle('axial_stress',190)]
    seg,zeta,r0,ri = obtain_kinematic_hardening_parameters(monotonic_curve,cyclic_curve,
                                                           position_x_list=plastic_strain_list,
                                                           yield_stress=yield_stress,
                                                           youngs_modulus=youngs_modulus,show=show)
    parameters.append([temperature,seg,zeta,r0,ri])    
    print parameters
#==============================================================================
# output elastic
#==============================================================================
    print >>outfile, """C <ELASTKIC CONSTANTS>
C <Young's modulus 3rd order fitting>
      EMOD=206308.7426+(-51.20306)*TEMP+0.01109*TEMP**2
     1     +(-3.84391E-05)*TEMP**3
C <Poisson's Ratio 3rd order fitting>
      ENU=2.901300E-01+(1.457750E-05)*TEMP
     1     +(-2.067420E-07)*TEMP**2+(2.780300E-10)*TEMP**3
      EG=0.5D0*EMOD/(1.0D0+ENU)
      EBULK=EMOD/(3.0D0*(1.0D0-2.0D0*ENU))
      ELAM=EBULK-EG*2.0D0/3.0D0
      ELAMK0=0.0D0-EG*2.0D0/3.0D0
C <YIELD STRESS>
      SY=%sD0
C <KINEMATIC HARDENING PARAMETERS>
      MU=0.2D0
      TQ=3.0D1""" % yield_stress
#    print >>outfile, """C <ELASTKIC CONSTANTS>
#C <Young's modulus 3rd order fitting>
#      EMOD=206308.7426+(-51.20306)*TEMP+0.01109*TEMP**2
#     1     +(-3.84391E-05)*TEMP**3
#C <Poisson's Ratio 3rd order fitting>
#      ENU=2.901300E-01+(1.457750E-05)*TEMP
#     1     +(-2.067420E-07)*TEMP**2+(2.780300E-10)*TEMP**3
#      EG=0.5D0*EMOD/(1.0D0+ENU)
#      EBULK=EMOD/(3.0D0*(1.0D0-2.0D0*ENU))
#      ELAM=EBULK-EG*2.0D0/3.0D0
#      ELAMK0=0.0D0-EG*2.0D0/3.0D0
#C <YIELD STRESS>
#      SY=1.15203105E+03+(-1.02927853E-01)*TEMP
#     1     +7.30935695E-05*TEMP**2+(-2.15967142E-07)*TEMP**3
#C <KINEMATIC HARDENING PARAMETERS>
#      MU=0.2D0
#      TQ=3.0D1"""
#==============================================================================
# output zeta
#==============================================================================
    for i in range(seg):
        lines='      ZKI(%u)=(%.2f)' %(i+1,parameters[0][2][i])
        print >>outfile, lines
#==============================================================================
# output r0
#==============================================================================
    for i in range(seg):
        x=np.array([parameters[0][0],parameters[1][0],parameters[2][0]])
        y=np.array([parameters[0][3][i],parameters[1][3][i],parameters[2][3][i]])
        z=np.polyfit(x,y,2)
        p = np.poly1d(z)
        lines='      R0KI(%u)=(%g)*TEMP**2+(%g)*TEMP+(%g)' %(i+1,z[0],z[1],z[2])
        print >>outfile, lines
#==============================================================================
# output r_delta
#==============================================================================
    for i in range(seg):
        x=np.array([parameters[0][0],parameters[1][0],parameters[2][0]])
        y=np.array([parameters[0][4][i],parameters[1][4][i],parameters[2][4][i]])
        z=np.polyfit(x,y,2)
        p = np.poly1d(z)
        lines='      RDSSKI(%u)=(%g)*TEMP**2+(%g)*TEMP+(%g)' %(i+1,z[0],z[1],z[2])
        print >>outfile, lines
#==============================================================================
# others
#==============================================================================
    print >>outfile, """
C <300C 200CYCLES>
C      CA1=0.4061D0
C      CB1=1.2746D0
C      CA2=1.0D0-CA1
C      CB2=9.4506D0
C      CB3=0.0D0
C <300C 800CYCLES>
C      CA1=0.3073D0
C      CB1=0.42549D0
C      CA2=1.0D0-CA1
C      CB2=5.78033D0
C      CB3=0.0D0
C <550C 200CYCLES>
C      CA1=0.4007D0
C      CB1=0.8672D0
C      CA2=1.0D0-CA1
C      CB2=11.81D0
C      CB3=0.0D0
C <650C 200CYCLES>
      CA1=0.4163D0
      CB1=0.4952D0
      CA2=1.0D0-CA1
      CB2=9.108D0
      CB3=0.0D0
C <ISOTROPIC PARAMETERS>
      YDNONPSS=100.0D0
      GAMAP=5.0D1
      GAMAQ=5.0D0
      YDPTS=0.0D0
      T1=20.0D0
      M1=2.0D0
C <NONPROPORTIONAL PARAMETERS>
      CCOEF=50.0D0
C <DAMAGE PARAMETERS>
      PD=1000.0D0"""
    outfile.close()

if __name__ == '__main__':
    calculate_umat_parameters_in718()

##for name in experiment_type_dict['NPR-IP']:
##for name in ['7211']:
#for name in ['7028']:    
#    calculate_umat_parameters_in718()
#    workbench(name,copy=False)
##    compare_exp_sim(name,2,'axial_stress','shear_stress')
#    compare_exp_sim(name,2,'axial_strain','axial_stress')