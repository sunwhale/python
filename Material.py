# -*- coding: utf-8 -*-
"""
Created on Thu Jan 05 11:50:47 2017

@author: j.Sun
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.optimize import leastsq
from Data import ExperimentData
from Constants import ExperimentDirectory
from Functions import calculate_elastic_by_temperature_in718

class Material:
    def __init__(self):
        self.name = None
        self.temperature = None
        self.youngs_modulus = None
        self.poisson_ratio = None
        self.shear_modulus = None
        self.yield_stress = None
        self.K = None
        self.n = None
        self.K_cyclic = None
        self.n_cyclic = None
        self.sigma_f = None
        self.b = None
        self.epsilon_f = None
        self.c = None
        self.tau_f = None
        self.b0 = None
        self.gamma_f = None
        self.c0 = None
        
    def show(self):
        line_format = '%-40s'
        line_format_string = line_format + '%-20s'
        line_format_temperature = line_format + '%.1fC'
        line_format_strain = line_format + '%.4f%%'
        line_format_stress = line_format + '%.2fMPa'
        line_format_coefficient = line_format + '%.6f'
        line_format_life = line_format + '%d'
        print '=========================Material========================='
        print line_format_string % ('Material:',self.name)
        print line_format_temperature % ('Temperature:',self.temperature)
        print line_format_stress % ('Young\'s modulus:',self.youngs_modulus)
        print line_format_coefficient % ('Poisson ratio:',self.poisson_ratio)
        print line_format_stress % ('Shear modulus:',self.shear_modulus)
        if self.K <> None:
            print line_format_stress % ('Strength coefficient:',self.K)
        if self.n <> None:
            print line_format_coefficient % ('Strain hardening exponent:',self.n)
        if self.sigma_f <> None:
            print line_format_stress % ('Fatigue strength coefficient:',self.sigma_f)
        if self.b <> None:
            print line_format_coefficient % ('Fatigue strength exponent:',self.b)
        if self.epsilon_f <> None:
            print line_format_stress % ('Fatigue ductility coefficient:',self.epsilon_f)
        if self.c <> None:
            print line_format_coefficient % ('Fatigue ductility exponent:',self.c)
        if self.K_cyclic <> None:
            print line_format_stress % ('Cyclic strength coefficient:',self.K_cyclic)
        if self.n_cyclic <> None:
            print line_format_coefficient % ('Cyclic strain hardening exponent:',self.n_cyclic)
        if self.tau_f <> None:
            print line_format_stress % ('Shear fatigue strength coefficient:',self.tau_f)
        if self.b0 <> None:
            print line_format_coefficient % ('Shear fatigue strength exponent:',self.b0)
        if self.gamma_f <> None:
            print line_format_stress % ('Shear fatigue ductility coefficient:',self.gamma_f)
        if self.c0 <> None:
            print line_format_coefficient % ('Shear fatigue ductility exponent:',self.c0)
    def setName(self,name=None):
        self.name = name
        
    def setTemperature(self,temperature=None):
        self.temperature = temperature
        
    def setMonotonic(self,youngs_modulus=None,poisson_ratio=None,yield_stress=None,K=None,n=None):
        self.youngs_modulus = youngs_modulus
        self.poisson_ratio = poisson_ratio
        self.shear_modulus = 0.5*youngs_modulus/(1.0+poisson_ratio)
        self.yield_stress = yield_stress
        self.K = K
        self.n = n
        
    def setCyclicAxial(self,sigma_f=None,b=None,epsilon_f=None,c=None,K_cyclic=None,n_cyclic=None):
        self.sigma_f = sigma_f
        self.b = b
        self.epsilon_f = epsilon_f
        self.c = c
        self.K_cyclic = K_cyclic
        self.n_cyclic = n_cyclic  
        
    def setCyclicTorsion(self,tau_f=None,b0=None,gamma_f=None,c0=None):
        self.tau_f = tau_f
        self.b0 = b0
        self.gamma_f = gamma_f
        self.c0 = c0
        if tau_f == None:
            self.tau_f = self.sigma_f/np.sqrt(3)
            self.b0 = self.b
            self.gamma_f = self.epsilon_f*np.sqrt(3)
            self.c0 = self.c

    def calculateRambergOsgood(self,x,y):
        def func(x,p):
            a,b = p
            return a*x**b
            
        def residuals(p,y,x):
            return y-func(x,p)
        
        p0 = [1000,1]
        plsq = leastsq(residuals, p0, args=(y, x))
        
        return plsq[0]
        
    def calculateCyclicStrengthCoefficient(self,epsilon,sigma):
        epsilon_p = epsilon - sigma/self.youngs_modulus
        epsilon_p_list = []
        sigma_list = []
        for i in range(len(epsilon_p)):
            if epsilon_p[i]>0.0005:
                epsilon_p_list.append(epsilon_p[i])
                sigma_list.append(sigma[i])
        epsilon_p = np.array(epsilon_p_list)
        sigma = np.array(sigma_list)
        
        self.K_cyclic = self.calculateRambergOsgood(epsilon_p,sigma)[0]
        self.n_cyclic = self.calculateRambergOsgood(epsilon_p,sigma)[1]
               
#        plt.plot(epsilon_p,sigma)
#        plt.plot(epsilon_p,self.K_cyclic*epsilon_p**self.n_cyclic)
#        plt.show()

    def calculateStrengthCoefficient(self,epsilon,sigma):
        epsilon_p = epsilon - sigma/self.youngs_modulus
        epsilon_p_list = []
        sigma_list = []
        for i in range(len(epsilon_p)):
            if epsilon_p[i]>0.0005 and epsilon_p[i]<0.015:
                epsilon_p_list.append(epsilon_p[i])
                sigma_list.append(sigma[i])
        epsilon_p = np.array(epsilon_p_list)
        sigma = np.array(sigma_list)
        self.K = self.calculateRambergOsgood(epsilon_p,sigma)[0]
        self.n = self.calculateRambergOsgood(epsilon_p,sigma)[1]

    def plotStrengthCoefficient(self):
        epsilon_p = np.arange(0,0.2,0.00001)
        sigma = self.K*epsilon_p**self.n
        epsilon = sigma/self.youngs_modulus + epsilon_p
        plt.plot(epsilon,sigma)
        plt.show()
        return epsilon,sigma
    
    def calcStrain(self,stress):
        return stress/self.youngs_modulus + (stress/self.K)**(1.0/self.n)

    def calcStrainAmplitude(self,stress_amplitude):
        return stress_amplitude/self.youngs_modulus + (stress_amplitude/self.K_cyclic)**(1.0/self.n_cyclic)
        
    def plotCyclicStrengthCoefficient(self):
        epsilon_p = np.arange(0,0.1,0.00001)
        sigma = self.K_cyclic*epsilon_p**self.n_cyclic
        epsilon = sigma/self.youngs_modulus + epsilon_p
        plt.plot(epsilon,sigma)
        plt.show()
        return epsilon,sigma
        
    def calculateMansonCoffinAxial(self,epsilon,life):
        sigma = []
        for e in epsilon:
            def f(x):
                x0 = float(x[0])
                return [
                    x0/self.youngs_modulus + (x0/self.K_cyclic)**(1.0/self.n_cyclic) - e
                ]
            
            result = fsolve(f, [1])
        
            sigma.append(result[0])
            
        epsilon_e = np.array(sigma)/self.youngs_modulus
        epsilon_p = epsilon - epsilon_e
        
        self.sigma_f = self.calculateRambergOsgood(life, sigma)[0]
        self.b = self.calculateRambergOsgood(life, sigma)[1]
        self.epsilon_f = self.calculateRambergOsgood(life, epsilon_p)[0]
        self.c = self.calculateRambergOsgood(life, epsilon_p)[1]
        
    def calculateMansonCoffinTorsion(self,gamma,life):
        sigma = []
        epsilon = gamma/np.sqrt(3.0)

        for e in epsilon:
            def f(x):
                x0 = float(x[0])
                return [
                    x0/self.youngs_modulus + (x0/self.K_cyclic)**(1.0/self.n_cyclic) - e
                ]
            
            result = fsolve(f, [1])
        
            sigma.append(result[0])

        tau = np.array(sigma)/np.sqrt(3.0)
        
        gamma_e = tau/self.shear_modulus
        gamma_p = gamma - gamma_e
        
        self.tau_f = self.calculateRambergOsgood(life, tau)[0]
        self.b0 = self.calculateRambergOsgood(life, tau)[1]
        self.gamma_f = self.calculateRambergOsgood(life, gamma_p)[0]
        self.c0 = self.calculateRambergOsgood(life, gamma_p)[1]
    
    def calculateMansonCoffinLife(self,epsilon):
        def f(x):
            x0 = float(x[0])
            return [
                self.sigma_f/self.youngs_modulus*(2*x0)**self.b + self.epsilon_f*(2*x0)**self.c-epsilon
            ]   
        result = fsolve(f, [1])
        return int(result[0])
        
    def plotMansonCoffinAxial(self):
        life_list = []
        seg = 10
        for i in range(1,8):
            for j in range(1,seg):
                life_list.append(float(j)/float(seg)*10**i)
        life = np.array(life_list)
        epsilon_amplitude = self.sigma_f/self.youngs_modulus*(2*life)**self.b + self.epsilon_f*(2*life)**self.c
#        plt.plot(life,epsilon)
#        plt.show()
        return life,epsilon_amplitude
        
#==============================================================================
# material inconel718 at 650
#==============================================================================
def material_in718():
    material = Material()
    material.setName(name='IN718')
    material.setTemperature(temperature=650.0)
    material.setMonotonic(youngs_modulus=167100.0,poisson_ratio=0.2886,yield_stress=1064.0,K=1433.0,n=0.048958)
    material.setCyclicAxial(sigma_f=1034.0,b=-0.04486,epsilon_f=0.11499,c=-0.52436,K_cyclic=1406.0,n_cyclic=0.10527)
    material.setCyclicTorsion()
    return material
    
def material_in718_NASA():
    material = Material()
    material.setName(name='IN718')
    material.setTemperature(temperature=650.0)
    material.setMonotonic(youngs_modulus=162600.0,poisson_ratio=0.2886,yield_stress=1064.0,K=1260.0,n=0.05030)
    material.setCyclicAxial(sigma_f=1348.0,b=-0.10052,epsilon_f=0.12445,c=-0.55218,K_cyclic=1827.0,n_cyclic=0.16723)
    material.setCyclicTorsion()
    return material
    
def material_in718_BHU():
    material = Material()
    material.setName(name='IN718')
    material.setTemperature(temperature=650.0)
    material.setMonotonic(youngs_modulus=177200.0,poisson_ratio=0.2886,yield_stress=1064.0)
    material.setCyclicAxial(sigma_f=985.0,b=-0.03917,epsilon_f=0.24721,c=-0.55682,K_cyclic=1420.0,n_cyclic=0.11332)
    material.setCyclicTorsion()
    return material
#==============================================================================
# calculate material inconel718
#==============================================================================\
def calc_material_in718():
    material = Material()
    material.setName(name='IN718')
    
    name = '7101' #300C
    temperature = 300.0
    material.setTemperature(temperature=temperature)
    youngs_modulus,poisson_ratio,shear_modulus = calculate_elastic_by_temperature_in718(temperature)
    material.setMonotonic(youngs_modulus=youngs_modulus,poisson_ratio=poisson_ratio,yield_stress=1121.9,K=1449.4,n=0.041214)
    
    name = '7103' #550C
    temperature = 550.0
    material.setTemperature(temperature=temperature)
    youngs_modulus,poisson_ratio,shear_modulus = calculate_elastic_by_temperature_in718(temperature)
    material.setMonotonic(youngs_modulus=youngs_modulus,poisson_ratio=poisson_ratio,yield_stress=1081.6,K=1398.7,n=0.041366)
    
    name = '7102' #650C
    temperature = 650.0
    material.setTemperature(temperature=temperature)
    youngs_modulus,poisson_ratio,shear_modulus = calculate_elastic_by_temperature_in718(temperature)
    material.setMonotonic(youngs_modulus=youngs_modulus,poisson_ratio=poisson_ratio,yield_stress=1056.7,K=1433.8,n=0.049103)
     
    exp_full_name = ExperimentDirectory + name + '.csv'
    exp = ExperimentData(exp_full_name)
    epsilon = exp.axial_strain
    sigma = exp.axial_stress
    plt.plot(epsilon,sigma)
    material.calculateStrengthCoefficient(epsilon,sigma)
    material.show()
    material.plotStrengthCoefficient()
    print material.K*0.002**material.n
    return material
    
#==============================================================================
# material ss304 at 20
#==============================================================================
def material_ss304():
    material = Material()
    material.setName(name='SS304')
    material.setTemperature(temperature=20.0)
    material.setMonotonic(youngs_modulus=198000.0,poisson_ratio=0.3,yield_stress=220.0)
    material.setCyclicAxial(sigma_f=671.22,b=-0.0842,epsilon_f=0.0931,c=-0.3792)
    material.setCyclicTorsion(tau_f=315.12,b0=-0.0556,gamma_f=0.1043,c0=-0.2649)
    
    epsilon, sigma = np.genfromtxt(r'F:\Work\2017-01-07_Exercise\CyclicData.txt', unpack=True)
    material.calculateCyclicStrengthCoefficient(epsilon,sigma)
    
    life_x2, epsilon_amp  = np.genfromtxt(r'F:\Work\2017-01-07_Exercise\epsilonN.txt', unpack=True)
    material.calculateMansonCoffinAxial(epsilon_amp/100.0,life_x2)
    
    life_x2, gamma_amp  = np.genfromtxt(r'F:\Work\2017-01-07_Exercise\gammaN.txt', unpack=True)
    material.calculateMansonCoffinTorsion(gamma_amp/100.0,life_x2)
    return material
    
#==============================================================================
# material cmsx-4 at 1000
#==============================================================================
def material_cmsx4():
    material = Material()
    material.setName(name='CMSX-4')
    material.setTemperature(temperature=1000.0)
    material.setMonotonic(youngs_modulus=77120.0,poisson_ratio=0.3,yield_stress=500.0)

    epsilon, sigma = np.genfromtxt(r'F:\Cloud\Database\CMSX4\MonotonicData.txt', unpack=True)
    material.calculateCyclicStrengthCoefficient(epsilon,sigma)
    
    life_x2 = np.array([1300,476,249])*2.0
    epsilon_amp = np.array([0.5,0.6,0.7])
    material.calculateMansonCoffinAxial(epsilon_amp/100.0,life_x2)
    
    material.setCyclicAxial(sigma_f=1700.22,b=-0.189252,epsilon_f=0.0,c=-0.0)
    material.setCyclicTorsion()
    return material
    
if __name__ == '__main__':
    material = material_cmsx4()
    material.show()