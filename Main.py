# -*- coding: utf-8 -*-
"""
Created on Thu Jan 05 11:50:47 2017

@author: j.Sun
"""

from Node import *
from Material import *
from Data import *
from Constants import *
from Work import *
from Plot import *

#==============================================================================
# material
#==============================================================================
material = Material()
material.setName(name='IN718')
material.setTemperature(temperature=650.0)
material.setMonotonic(youngs_modulus=167100.0,poisson_ratio=0.2886,yield_stress=1064.0)
material.setCyclicAxial(sigma_f=1034.0,b=-0.04486,epsilon_f=0.11499,c=-0.52436)
material.setCyclicTorsion(tau_f=1034.0/np.sqrt(3),b0=-0.04486,gamma_f=0.11499*np.sqrt(3),c0=-0.52436)
material.show()
#==============================================================================
# job
#==============================================================================
name = '7037'

experiment_log = ExperimentLog(ExperimentDirectory + 'Inconel718_test_log.csv')
experiment_log.find(name)

sim = SimulationData(r'F:\Temp\IN7182\7037\7037.csv')
sim.obtainCountIndex(180)
nth = 20
time = sim.obtainNthCycle('runing_time',nth)
length = len(time)
s11 = sim.obtainNthCycle('axial_stress',nth)
s22 = np.zeros(length)
s33 = np.zeros(length)
s12 = sim.obtainNthCycle('shear_stress',nth)
s13 = np.zeros(length)
s23 = np.zeros(length)
e11 = sim.obtainNthCycle('axial_strain',nth)
e22 = e11*-1.0*material.poisson_ratio
e33 = e11*-1.0*material.poisson_ratio
e12 = sim.obtainNthCycle('shear_strain',nth)/2.0
e13 = np.zeros(length)
e23 = np.zeros(length)

#exp = ExperimentData('F:\\Database\\IN718\\Timed\\7036.csv')
#exp.obtainCountIndex(180)
#nth = 10
#time = exp.obtainNthCycle('runing_time',nth)[::2]
#length = len(time)
#s11 = exp.obtainNthCycle('axial_stress',nth)[::2]
#s22 = np.zeros(length)
#s33 = np.zeros(length)
#s12 = exp.obtainNthCycle('shear_stress',nth)[::2]
#s13 = np.zeros(length)
#s23 = np.zeros(length)
#e11 = exp.obtainNthCycle('axial_strain',nth)[::2]
#e22 = e11*-1.0*material.poisson_ratio
#e33 = e11*-1.0*material.poisson_ratio
#e12 = exp.obtainNthCycle('shear_strain',nth)[::2]/2.0
#e13 = np.zeros(length)
#e23 = np.zeros(length)
#==============================================================================
# plot
#==============================================================================
#xitem = 'axial_stress'
#yitem = 'shear_stress'
#x1 = sim.obtainNthCycle(xitem,nth)
#y1 = sim.obtainNthCycle(yitem,nth)
#x2 = exp.obtainNthCycle(xitem,nth)
#y2 = exp.obtainNthCycle(yitem,nth)
#plt.xlabel(xylabels[xitem])
#plt.ylabel(xylabels[yitem])
#plt.plot(x1, y1, label='1',linewidth=1,marker='',ls='-')
#plt.plot(x2, y2, label='1',linewidth=1,marker='',ls='-')
#plt.gca().set_aspect('auto')
#plt.show()
#==============================================================================
# calculate life
#==============================================================================
stress=[]
strain=[]
for i in range(len(time)):
    stress.append([[s11[i],s12[i],s13[i]],[s12[i],s22[i],s23[i]],[s13[i],s23[i],s33[i]]])
    strain.append([[e11[i],e12[i],e13[i]],[e12[i],e22[i],e23[i]],[e13[i],e23[i],e33[i]]])

node = Node(nodelabel=1, dimension=2, time=time, coordinate=[], displacement=[], stress=stress, strain=strain)
#print node.stress_list[0]
#print node.strain_list[0]
fatigue_life = node.fatigueLifeFSModel(material)
print 'fatigue_life',fatigue_life

resultfile = open('out0.dat', 'w')
line = ''
line += '%s    ' % (node.nodelabel)
line += '%s'     % (fatigue_life)
print >>resultfile, line
resultfile.close()            