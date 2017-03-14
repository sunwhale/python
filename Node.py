# -*- coding: utf-8 -*-
"""
Created on Thu Jan 05 11:50:47 2017

@author: j.Sun
"""

import numpy as np
from scipy.optimize import fsolve

class Node:
    def __init__(self, nodelabel=1, dimension=3, time=[], coordinate=[], 
                 displacement=[], stress=[], strain=[], temperature=[]):
        self.nodelabel = nodelabel
        self.dimension = dimension
        self.time_list = [t for t in time]
        self.temperature_list = [t for t in temperature]
        self.coordinate_list = [[i for i in c[:dimension]] for c in coordinate]
        self.displacement_list = [[i for i in d[:dimension]] for d in displacement]
        self.strain_list = [[i[:dimension] for i in s[:dimension]] for s in strain]
        self.stress_list = [[i[:dimension] for i in s[:dimension]] for s in stress]
        
    def __str__(self):
        return 'node %s at time %s' % (str(int(self.nodelabel)),str(self.time))
        
    def stressTransform(self, transformation):
        stress_list = []
        for stress in self.stress_list:
            stress_list.append(np.dot(np.dot(np.array(transformation),np.array(stress)),np.array(transformation).T))
        return stress_list
        
    def strainTransform(self, transformation):
        strain_list = []
        for strain in self.strain_list:
            strain_list.append(np.dot(np.dot(np.array(transformation),np.array(strain)),np.array(transformation).T))
        return strain_list
        
    def normalStress(self, transformation):
        transform_stress_list = self.stressTransform(transformation)
        normalstress_list = []
        for transform_stress in transform_stress_list:
            normalstress_list.append(transform_stress[0][0])
        return normalstress_list
        
    def normalStrain(self, transformation):
        transform_strain_list = self.strainTransform(transformation)
        normalstrain_list = []
        for transform_strain in transform_strain_list:
            normalstrain_list.append(transform_strain[0][0])
        return normalstrain_list
        
    def shearStressEqv(self, transformation):
        transform_stress_list = self.stressTransform(transformation)
        shearstress_list = []
        for transform_stress in transform_stress_list:
            if self.dimension == 3:
                shearstress_list.append(np.sqrt(transform_stress[0][1]**2+transform_stress[0][2]**2))
            if self.dimension == 2:
                shearstress_list.append(transform_stress[0][1])
        return shearstress_list

    def shearStrainEqv(self, transformation):
        transform_strain_list = self.strainTransform(transformation)
        shearstrain_list = []
        for transform_strain in transform_strain_list:
            if self.dimension == 3:
                shearstrain_list.append(2.0*np.sqrt(transform_strain[0][1]**2+transform_strain[0][2]**2))
            if self.dimension == 2:
                shearstrain_list.append(2.0*transform_strain[0][1])
        return shearstrain_list
        
    def shearStress(self, transformation):
        transform_stress_list = self.stressTransform(transformation)
        shearstress_list = []
        for transform_stress in transform_stress_list:
            if self.dimension == 3:
                shearstress_list.append([transform_stress[0][1],transform_stress[0][2]])
            if self.dimension == 2:
                shearstress_list.append([transform_stress[0][1],0.0])
        return shearstress_list
        
    def shearStrain(self, transformation):
        transform_strain_list = self.strainTransform(transformation)
        shearstrain_list = []
        for transform_strain in transform_strain_list:
            if self.dimension == 3:
                shearstrain_list.append([transform_strain[0][1],transform_strain[0][2]])
            if self.dimension == 2:
                shearstrain_list.append([transform_strain[0][1],0.0])
        return shearstrain_list
        
    def deltaGamma(self, transformation):
        gamma_list = np.array(self.shearStrain(transformation))
        if self.dimension == 3:
            delta_gamma_max = 0.0
            for i in range(len(gamma_list)):
                for j in range(i,len(gamma_list)):
                    delta_gamma = 2.0*np.linalg.norm(gamma_list[i]-gamma_list[j])
                    if delta_gamma > delta_gamma_max:
                        delta_gamma_max = delta_gamma
            return delta_gamma_max
        if self.dimension == 2:
            return 2*(max(gamma_list[:,0])-min(gamma_list[:,0]))

    def deltaTau(self, transformation):
        tau_list = np.array(self.shearStress(transformation))
        if self.dimension == 3:
            delta_tau_max = 0.0
            for i in range(len(tau_list)):
                for j in range(i,len(tau_list)):
                    delta_tau = np.linalg.norm(tau_list[i]-tau_list[j])
                    if delta_tau > delta_tau_max:
                        delta_tau_max = delta_tau               
            return delta_tau_max
        if self.dimension == 2:
            return max(tau_list[:,0])-min(tau_list[:,0])
            
    def deltaEpsilon(self, transformation):
        return max(self.normalStrain(transformation))-min(self.normalStrain(transformation))
        
    def deltaSigma(self, transformation):
        return max(self.normalStress(transformation))-min(self.normalStress(transformation))
        
    def sigmaNmax(self, transformation):
        return max(self.normalStress(transformation))
    
    def temperatureAtSigmaNmax(self, transformation):
        normalstress_list = self.normalStress(transformation)
        max_normalstress = max(normalstress_list)
        max_normalstress_index = normalstress_list.index(max_normalstress)
        temperature_at_sigma_nmax = self.temperature_list[max_normalstress_index]
        return temperature_at_sigma_nmax
        
    def tauNmax(self, transformation):
        return max(self.shearStressEqv(transformation))
    
    def initialValues(self):
        delta_gamma_max = 0.0
        theta_critical_plane = 0.0
        phi_critical_plane = 0.0
        transformation_critical_plane = []

        if self.dimension == 3:
            theta_deg_list = range(0,180,5)
            phi_deg_list = range(0,180,5)
        if self.dimension == 2:
            theta_deg_list = [90]
            phi_deg_list = range(0,180,1)
            
        return (delta_gamma_max,theta_critical_plane,phi_critical_plane,
                transformation_critical_plane,theta_deg_list,phi_deg_list)
                
    def calcTransformation(self,theta,phi):
        transformation=[
        [np.sin(theta)*np.cos(phi) ,np.sin(theta)*np.sin(phi) ,np.cos(theta)],
        [-np.sin(phi)              ,np.cos(phi)               ,0            ],
        [-np.cos(theta)*np.cos(phi),-np.cos(theta)*np.sin(phi),np.sin(theta)]
        ]
        transformation=[t[:self.dimension] for t in transformation[:self.dimension]]
        return transformation
        
    def fatigueLifeFSModel(self,Material):
        print '=========================FS model ========================='
        k=1.0

        (delta_gamma_max,theta_critical_plane,phi_critical_plane,
         transformation_critical_plane,theta_deg_list,phi_deg_list) = self.initialValues()
        
        for theta_deg in theta_deg_list:
            theta=np.radians(theta_deg)
            for phi_deg in phi_deg_list:
                phi=np.radians(phi_deg)
                transformation = self.calcTransformation(theta,phi)
                delta_gamma_theta_phi = self.deltaGamma(transformation)
                if delta_gamma_theta_phi > delta_gamma_max:
                    delta_gamma_max = delta_gamma_theta_phi
                    theta_critical_plane = theta_deg
                    phi_critical_plane = phi_deg
                    transformation_critical_plane = transformation
                    
        if self.dimension == 3:
            print 'theta_critical_plane',theta_critical_plane
            print 'phi_critical_plane',phi_critical_plane
        if self.dimension == 2:
            print 'phi_critical_plane',phi_critical_plane
        
        sigma_nmax_critical_plane = self.sigmaNmax(transformation_critical_plane)
        delta_sigma_critical_plane = self.deltaSigma(transformation_critical_plane)
        delta_epsilon_critical_plane = self.deltaEpsilon(transformation_critical_plane)
        tau_nmax_critical_plane = self.tauNmax(transformation_critical_plane)
        delta_tau_critical_plane = self.deltaTau(transformation_critical_plane)
        temperature_at_sigma_nmax_critical_plane = self.temperatureAtSigmaNmax(transformation_critical_plane)
        
        fs_coefficient=delta_gamma_max/2.0*(1.0+k*sigma_nmax_critical_plane/Material.yield_stress)
                
        def f(x):
            x0 = float(x[0])
            return [
                Material.tau_f/Material.shear_modulus*(2*x0)**Material.b0 + Material.gamma_f*(2*x0)**Material.c0-fs_coefficient
            ]
        
        fatigue_coefficient = fs_coefficient
        result = fsolve(f, [1])
        fatigue_life = int(result[0])
        
        self.outputFatigueLife(delta_gamma_max,sigma_nmax_critical_plane,
                          delta_sigma_critical_plane,delta_epsilon_critical_plane,
                          tau_nmax_critical_plane,delta_tau_critical_plane,
                          temperature_at_sigma_nmax_critical_plane,fatigue_coefficient,fatigue_life)
        
        return [phi_critical_plane,sigma_nmax_critical_plane,delta_sigma_critical_plane,
                delta_epsilon_critical_plane,tau_nmax_critical_plane,
                delta_tau_critical_plane,delta_gamma_max,fatigue_life,
                fatigue_coefficient,temperature_at_sigma_nmax_critical_plane]
                
    def fatigueLifeSWTModel(self,Material):
        print '=========================SWT model ========================='

        (delta_gamma_max,theta_critical_plane,phi_critical_plane,
         transformation_critical_plane,theta_deg_list,phi_deg_list) = self.initialValues()
        
        for theta_deg in theta_deg_list:
            theta=np.radians(theta_deg)
            for phi_deg in phi_deg_list:
                phi=np.radians(phi_deg)
                transformation = self.calcTransformation(theta,phi)
                delta_gamma_theta_phi = self.deltaGamma(transformation)
                if delta_gamma_theta_phi > delta_gamma_max:
                    delta_gamma_max = delta_gamma_theta_phi
                    theta_critical_plane = theta_deg
                    phi_critical_plane = phi_deg
                    transformation_critical_plane = transformation
                    
        if self.dimension == 3:
            print 'theta_critical_plane',theta_critical_plane
            print 'phi_critical_plane',phi_critical_plane
        if self.dimension == 2:
            print 'phi_critical_plane',phi_critical_plane
        
        sigma_nmax_critical_plane = self.sigmaNmax(transformation_critical_plane)
        delta_sigma_critical_plane = self.deltaSigma(transformation_critical_plane)
        delta_epsilon_critical_plane = self.deltaEpsilon(transformation_critical_plane)
        tau_nmax_critical_plane = self.tauNmax(transformation_critical_plane)
        delta_tau_critical_plane = self.deltaTau(transformation_critical_plane)
        temperature_at_sigma_nmax_critical_plane = self.temperatureAtSigmaNmax(transformation_critical_plane)
        
        swt_coefficient = sigma_nmax_critical_plane*delta_epsilon_critical_plane/2.0
    
        def f(x):
            x0 = float(x[0])
            return [
                Material.sigma_f*Material.sigma_f/Material.youngs_modulus*(2*x0)**(2*Material.b) + Material.sigma_f*Material.epsilon_f*(2*x0)**(Material.b+Material.c)-swt_coefficient
            ]
                
        fatigue_coefficient = swt_coefficient
        result = fsolve(f, [1])
        fatigue_life = int(result[0])
        
        self.outputFatigueLife(delta_gamma_max,sigma_nmax_critical_plane,
                          delta_sigma_critical_plane,delta_epsilon_critical_plane,
                          tau_nmax_critical_plane,delta_tau_critical_plane,
                          temperature_at_sigma_nmax_critical_plane,fatigue_coefficient,fatigue_life)
        
        return [phi_critical_plane,sigma_nmax_critical_plane,delta_sigma_critical_plane,
                delta_epsilon_critical_plane,tau_nmax_critical_plane,
                delta_tau_critical_plane,delta_gamma_max,fatigue_life,
                fatigue_coefficient,temperature_at_sigma_nmax_critical_plane]
                
    def outputFatigueLife(self,delta_gamma_max,sigma_nmax_critical_plane,
                          delta_sigma_critical_plane,delta_epsilon_critical_plane,
                          tau_nmax_critical_plane,delta_tau_critical_plane,
                          temperature_at_sigma_nmax_critical_plane,fs_coefficient,fatigue_life):
        line_format = '%-40s'
        line_format_strain = line_format + '%.4f%%'
        line_format_stress = line_format + '%.2fMPa'
        line_format_coefficient = line_format + '%.6f'
        line_format_life = line_format + '%d'
        line_format_temperature = line_format + '%.1fC'
        print line_format_strain % ('delta_gamma_max',delta_gamma_max*100)
        print line_format_stress % ('sigma_nmax_critical_plane',sigma_nmax_critical_plane)
        print line_format_stress % ('delta_sigma_critical_plane',delta_sigma_critical_plane)
        print line_format_strain % ('delta_epsilon_critical_plane',delta_epsilon_critical_plane*100)
        print line_format_stress % ('tau_nmax_critical_plane',tau_nmax_critical_plane)
        print line_format_stress % ('delta_tau_critical_plane',delta_tau_critical_plane)
        print line_format_temperature % ('temperature_at_sigma_nmax_critical_plane',temperature_at_sigma_nmax_critical_plane)
        print line_format_coefficient % ('fs_coefficient',fs_coefficient)
        print line_format_life % ('fatigue_life',fatigue_life)
        
    def test(self):
        theta_deg=90.0
        phi_deg=0.0
        theta=np.radians(theta_deg)
        phi=np.radians(phi_deg)
        
        stress=np.array([[2.0,1.0,0.0],[1.0,0.0,0.0],[0.0,0.0,0.0]])
        strain=np.array([[2.0,1.0,0.0],[1.0,0.0,0.0],[0.0,0.0,0.0]])

#        print stress
#        print strain
        transformation=np.array([
        [np.sin(theta)*np.cos(phi) ,np.sin(theta)*np.sin(phi) ,np.cos(theta)],
        [-np.sin(phi)              ,np.cos(phi)               ,0            ],
        [-np.cos(theta)*np.cos(phi),-np.cos(theta)*np.sin(phi),np.sin(theta)]
        ])
        print transformation
        
#        print np.dot(np.dot(transformation,stress),transformation.T)
#        print np.dot(np.dot(transformation,strain),transformation.T)
        
        s11 = stress[0][0]
        s22 = stress[1][1]
        s12 = stress[0][1]
        
        e11 = strain[0][0]
        e22 = strain[1][1]
        e12 = strain[0][1]*2
        
        sigma_phi      =s11/2+s11/2*np.cos(2*phi)+s12*np.sin(2*phi)
        tau_phi        =s11/2*np.sin(2*phi)-s12*np.cos(2*phi)
#        sigma_phi      =(s11+s22)/2+(s11-s22)/2*np.cos(2*phi)+s12*np.sin(2*phi)
#        tau_phi        =(s11-s22)/2*np.sin(2*phi)-s12*np.cos(2*phi)
        epsilon_phi    =(e11+e22)/2+(e11-e22)/2*np.cos(2*phi)+e12/2*np.sin(2*phi)
        gamma_phi_half =(e11-e22)/2*np.sin(2*phi)-e12/2*np.cos(2*phi)
        
#        print sigma_phi, tau_phi
#        print epsilon_phi, gamma_phi_half
        
#n = Node()
#n.test()