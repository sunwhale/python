# -*- coding: utf-8 -*-
"""
Created on Thu Jan 05 11:50:47 2017

@author: j.Sun
"""

import numpy as np
from scipy.optimize import fsolve
from Material import Material

class Node:
    def __init__(self, nodelabel=1, dimension=2, time=[], coordinate=[], 
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
        normal_stress_list = []
        for transform_stress in transform_stress_list:
            normal_stress_list.append(transform_stress[0][0])
        return normal_stress_list
        
    def normalStrain(self, transformation):
        transform_strain_list = self.strainTransform(transformation)
        normal_strain_list = []
        for transform_strain in transform_strain_list:
            normal_strain_list.append(transform_strain[0][0])
        return normal_strain_list
        
    def shearStress(self, transformation):
        transform_stress_list = self.stressTransform(transformation)
        shear_stress_list = []
        for transform_stress in transform_stress_list:
            if self.dimension == 3:
                shear_stress_list.append([transform_stress[0][1],transform_stress[0][2]])
            if self.dimension == 2:
                shear_stress_list.append([transform_stress[0][1],0.0])
        return shear_stress_list
        
    def shearStrain(self, transformation):
        transform_strain_list = self.strainTransform(transformation)
        shear_strain_list = []
        for transform_strain in transform_strain_list:
            if self.dimension == 3:
                shear_strain_list.append([2.0*transform_strain[0][1],2.0*transform_strain[0][2]])
            if self.dimension == 2:
                shear_strain_list.append([2.0*transform_strain[0][1],0.0])
        return shear_strain_list
        
    def shearStressEqv(self, transformation):
        shear_stress_list = self.shearStress(transformation)
        shear_stress_eqv_list = []
        for shear_stress in shear_stress_list:
            shear_stress_eqv_list.append(np.sqrt(shear_stress[0]**2+shear_stress[1]**2))
        return shear_stress_eqv_list

    def shearStrainEqv(self, transformation):
        shear_strain_list = self.shearStrain(transformation)
        shear_strain_eqv_list = []
        for shear_strain in shear_strain_list:
            shear_strain_eqv_list.append(np.sqrt(shear_strain[0]**2+shear_strain[1]**2))
        return shear_strain_eqv_list
        
    def deltaGamma(self, transformation):
        gamma_list = np.array(self.shearStrain(transformation))
        if self.dimension == 3:
            delta_gamma_max = 0.0
            for i in range(len(gamma_list)):
                for j in range(i,len(gamma_list)):
                    delta_gamma = np.linalg.norm(gamma_list[i]-gamma_list[j])
                    if delta_gamma > delta_gamma_max:
                        delta_gamma_max = delta_gamma
            return delta_gamma_max
        if self.dimension == 2:
            return (max(gamma_list[:,0])-min(gamma_list[:,0]))

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
        
    def sigmaNMax(self, transformation):
        return max(self.normalStress(transformation))
        
    def tauNMax(self, transformation):
        return max(self.shearStressEqv(transformation))
        
    def tensionEnergy(self, transformation):
        return self.deltaSigma(transformation)*self.deltaEpsilon(transformation)

    def shearEnergy(self, transformation):
        return self.deltaTau(transformation)*self.deltaGamma(transformation)
        
    def temperatureAtsigmaNMax(self, transformation):
        normal_stress_list = self.normalStress(transformation)
        max_normal_stress = max(normal_stress_list)
        max_normal_stress_index = normal_stress_list.index(max_normal_stress)
        temperature_at_sigma_nmax = self.temperature_list[max_normal_stress_index]
        return temperature_at_sigma_nmax
        
    def fatigueLifeFSModel(self,Material,k=0.3):
        print '=========================FS model========================='

        (delta_gamma_max,
         delta_epsilon_max,
         tension_energy_max,
         shear_energy_max,
         total_energy_max,
         theta_critical_plane,
         phi_critical_plane,
         transformation_critical_plane,
         theta_deg_list,
         phi_deg_list) = self.initialValues()
        
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
                    
        (sigma_nmax_critical_plane,
         delta_sigma_critical_plane,
         delta_epsilon_critical_plane,
         tau_nmax_critical_plane,
         delta_tau_critical_plane,
         delta_gamma_critical_plane,
         temperature_at_sigma_nmax_critical_plane) = self.calcValuesAtCriticalPlane(transformation_critical_plane)
        
        fs_coefficient=delta_gamma_max/2.0*(1.0+k*sigma_nmax_critical_plane/Material.yield_stress)
                
        def f(x):
            x0 = float(x[0])
            return [
                Material.tau_f/Material.shear_modulus*(2*x0)**Material.b0 + Material.gamma_f*(2*x0)**Material.c0-fs_coefficient
            ]
        
        fatigue_coefficient = fs_coefficient
        result = fsolve(f, [1])
        fatigue_life = int(result[0])
        
        return self.outputFatigueLife(theta_critical_plane,
                                      phi_critical_plane,
                                      delta_gamma_critical_plane,
                                      sigma_nmax_critical_plane,
                                      delta_sigma_critical_plane,
                                      delta_epsilon_critical_plane,
                                      tau_nmax_critical_plane,
                                      delta_tau_critical_plane,
                                      temperature_at_sigma_nmax_critical_plane,
                                      fatigue_coefficient,
                                      fatigue_life)
                
    def fatigueLifeSWTModel(self,Material):
        print '=========================SWT model========================='

        (delta_gamma_max,
         delta_epsilon_max,
         tension_energy_max,
         shear_energy_max,
         total_energy_max,
         theta_critical_plane,
         phi_critical_plane,
         transformation_critical_plane,
         theta_deg_list,
         phi_deg_list) = self.initialValues()
        
        for theta_deg in theta_deg_list:
            theta=np.radians(theta_deg)
            for phi_deg in phi_deg_list:
                phi=np.radians(phi_deg)
                transformation = self.calcTransformation(theta,phi)
                delta_epsilon_theta_phi = self.deltaEpsilon(transformation)
                if delta_epsilon_theta_phi > delta_epsilon_max:
                    delta_epsilon_max = delta_epsilon_theta_phi
                    theta_critical_plane = theta_deg
                    phi_critical_plane = phi_deg
                    transformation_critical_plane = transformation
                    
        (sigma_nmax_critical_plane,
         delta_sigma_critical_plane,
         delta_epsilon_critical_plane,
         tau_nmax_critical_plane,
         delta_tau_critical_plane,
         delta_gamma_critical_plane,
         temperature_at_sigma_nmax_critical_plane) = self.calcValuesAtCriticalPlane(transformation_critical_plane)
        
        swt_coefficient = sigma_nmax_critical_plane*delta_epsilon_critical_plane/2.0
    
        def f(x):
            x0 = float(x[0])
            return [
                Material.sigma_f*Material.sigma_f/Material.youngs_modulus*(2*x0)**(2*Material.b) + Material.sigma_f*Material.epsilon_f*(2*x0)**(Material.b+Material.c)-swt_coefficient
            ]
                
        fatigue_coefficient = swt_coefficient
        result = fsolve(f, [1])
        fatigue_life = int(result[0])
        
        return self.outputFatigueLife(theta_critical_plane,
                                      phi_critical_plane,
                                      delta_gamma_critical_plane,
                                      sigma_nmax_critical_plane,
                                      delta_sigma_critical_plane,
                                      delta_epsilon_critical_plane,
                                      tau_nmax_critical_plane,
                                      delta_tau_critical_plane,
                                      temperature_at_sigma_nmax_critical_plane,
                                      fatigue_coefficient,
                                      fatigue_life)
                                      
    def fatigueLifeBMModel(self,Material,S=0.45):
        print '=========================BM model========================='

        (delta_gamma_max,
         delta_epsilon_max,
         tension_energy_max,
         shear_energy_max,
         total_energy_max,
         theta_critical_plane,
         phi_critical_plane,
         transformation_critical_plane,
         theta_deg_list,
         phi_deg_list) = self.initialValues()
        
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
                    
        (sigma_nmax_critical_plane,
         delta_sigma_critical_plane,
         delta_epsilon_critical_plane,
         tau_nmax_critical_plane,
         delta_tau_critical_plane,
         delta_gamma_critical_plane,
         temperature_at_sigma_nmax_critical_plane) = self.calcValuesAtCriticalPlane(transformation_critical_plane)
        
        bm_coefficient = delta_gamma_critical_plane/2.0 + S*delta_epsilon_critical_plane
        sigma_nmean_critical_plane = sigma_nmax_critical_plane - delta_sigma_critical_plane/2.0
                
        def f(x):
            x0 = float(x[0])
            return [
                (1.3+0.7*S)*(Material.sigma_f-2.0*sigma_nmean_critical_plane)/Material.youngs_modulus*(2*x0)**Material.b + (1.5+0.5*S)*Material.epsilon_f*(2*x0)**Material.c-bm_coefficient
            ]
                
        fatigue_coefficient = bm_coefficient
        result = fsolve(f, [1])
        fatigue_life = int(result[0])
        
        return self.outputFatigueLife(theta_critical_plane,
                                      phi_critical_plane,
                                      delta_gamma_critical_plane,
                                      sigma_nmax_critical_plane,
                                      delta_sigma_critical_plane,
                                      delta_epsilon_critical_plane,
                                      tau_nmax_critical_plane,
                                      delta_tau_critical_plane,
                                      temperature_at_sigma_nmax_critical_plane,
                                      fatigue_coefficient,
                                      fatigue_life)

    def fatigueLifeLiu1Model(self,Material):
        print '=========================Liu1 model========================='

        (delta_gamma_max,
         delta_epsilon_max,
         tension_energy_max,
         shear_energy_max,
         total_energy_max,
         theta_critical_plane,
         phi_critical_plane,
         transformation_critical_plane,
         theta_deg_list,
         phi_deg_list) = self.initialValues()
        
        for theta_deg in theta_deg_list:
            theta=np.radians(theta_deg)
            for phi_deg in phi_deg_list:
                phi=np.radians(phi_deg)
                transformation = self.calcTransformation(theta,phi)
                tension_energy_theta_phi = self.tensionEnergy(transformation)
                if tension_energy_theta_phi > tension_energy_max:
                    tension_energy_max = tension_energy_theta_phi
                    theta_critical_plane = theta_deg
                    phi_critical_plane = phi_deg
                    transformation_critical_plane = transformation
                    
        (sigma_nmax_critical_plane,
         delta_sigma_critical_plane,
         delta_epsilon_critical_plane,
         tau_nmax_critical_plane,
         delta_tau_critical_plane,
         delta_gamma_critical_plane,
         temperature_at_sigma_nmax_critical_plane) = self.calcValuesAtCriticalPlane(transformation_critical_plane)
        
        delta_w1_max = self.tensionEnergy(transformation_critical_plane) + self.shearEnergy(transformation_critical_plane)
        stress_ratio = (sigma_nmax_critical_plane-delta_sigma_critical_plane)/sigma_nmax_critical_plane
        delta_w1_max *= 2.0/(1.0-stress_ratio)
        
        def f(x):
            x0 = float(x[0])
            return [
                4*Material.sigma_f*Material.sigma_f/Material.youngs_modulus*(2*x0)**(2*Material.b) + 4*Material.sigma_f*Material.epsilon_f*(2*x0)**(Material.b+Material.c)-delta_w1_max
            ]
                
        fatigue_coefficient = delta_w1_max
        result = fsolve(f, [1])
        fatigue_life = int(result[0])
        
        return self.outputFatigueLife(theta_critical_plane,
                                      phi_critical_plane,
                                      delta_gamma_critical_plane,
                                      sigma_nmax_critical_plane,
                                      delta_sigma_critical_plane,
                                      delta_epsilon_critical_plane,
                                      tau_nmax_critical_plane,
                                      delta_tau_critical_plane,
                                      temperature_at_sigma_nmax_critical_plane,
                                      fatigue_coefficient,
                                      fatigue_life)
                                      
    def fatigueLifeLiu2Model(self,Material):
        print '=========================Liu2 model========================='

        (delta_gamma_max,
         delta_epsilon_max,
         tension_energy_max,
         shear_energy_max,
         total_energy_max,
         theta_critical_plane,
         phi_critical_plane,
         transformation_critical_plane,
         theta_deg_list,
         phi_deg_list) = self.initialValues()
        
        for theta_deg in theta_deg_list:
            theta=np.radians(theta_deg)
            for phi_deg in phi_deg_list:
                phi=np.radians(phi_deg)
                transformation = self.calcTransformation(theta,phi)
                shear_energy_theta_phi = self.shearEnergy(transformation)
                if shear_energy_theta_phi > shear_energy_max:
                    shear_energy_max = shear_energy_theta_phi
                    theta_critical_plane = theta_deg
                    phi_critical_plane = phi_deg
                    transformation_critical_plane = transformation
                    
        (sigma_nmax_critical_plane,
         delta_sigma_critical_plane,
         delta_epsilon_critical_plane,
         tau_nmax_critical_plane,
         delta_tau_critical_plane,
         delta_gamma_critical_plane,
         temperature_at_sigma_nmax_critical_plane) = self.calcValuesAtCriticalPlane(transformation_critical_plane)
        
        delta_w2_max = self.tensionEnergy(transformation_critical_plane) + self.shearEnergy(transformation_critical_plane)
        sigma_nmean_critical_plane = sigma_nmax_critical_plane - delta_sigma_critical_plane/2.0
        delta_w2_max *= Material.sigma_f/(Material.sigma_f-sigma_nmean_critical_plane)
        
        def f(x):
            x0 = float(x[0])
            return [
                4*Material.tau_f*Material.tau_f/Material.shear_modulus*(2*x0)**(2*Material.b0) + 4*Material.tau_f*Material.gamma_f*(2*x0)**(Material.b0+Material.c0)-delta_w2_max
            ]
                
        fatigue_coefficient = delta_w2_max
        result = fsolve(f, [1])
        fatigue_life = int(result[0])
        
        return self.outputFatigueLife(theta_critical_plane,
                                      phi_critical_plane,
                                      delta_gamma_critical_plane,
                                      sigma_nmax_critical_plane,
                                      delta_sigma_critical_plane,
                                      delta_epsilon_critical_plane,
                                      tau_nmax_critical_plane,
                                      delta_tau_critical_plane,
                                      temperature_at_sigma_nmax_critical_plane,
                                      fatigue_coefficient,
                                      fatigue_life)
                                      
    def fatigueLifeChuModel(self,Material):
        print '=========================Chu model========================='

        (delta_gamma_max,
         delta_epsilon_max,
         tension_energy_max,
         shear_energy_max,
         total_energy_max,
         theta_critical_plane,
         phi_critical_plane,
         transformation_critical_plane,
         theta_deg_list,
         phi_deg_list) = self.initialValues()
        
        for theta_deg in theta_deg_list:
            theta=np.radians(theta_deg)
            for phi_deg in phi_deg_list:
                phi=np.radians(phi_deg)
                transformation = self.calcTransformation(theta,phi)
                total_energy_theta_phi = self.sigmaNMax(transformation)*self.deltaEpsilon(transformation)/2.0+self.tauNMax(transformation)*self.deltaGamma(transformation)/2.0
                if total_energy_theta_phi > total_energy_max:
                    total_energy_max = total_energy_theta_phi
                    theta_critical_plane = theta_deg
                    phi_critical_plane = phi_deg
                    transformation_critical_plane = transformation
                    
        (sigma_nmax_critical_plane,
         delta_sigma_critical_plane,
         delta_epsilon_critical_plane,
         tau_nmax_critical_plane,
         delta_tau_critical_plane,
         delta_gamma_critical_plane,
         temperature_at_sigma_nmax_critical_plane) = self.calcValuesAtCriticalPlane(transformation_critical_plane)
        
        delta_w_max = sigma_nmax_critical_plane*delta_epsilon_critical_plane/2.0 + tau_nmax_critical_plane*delta_gamma_critical_plane/2.0
                           
        def f(x):
            x0 = float(x[0])
            return [
                1.02*Material.sigma_f*Material.sigma_f/Material.youngs_modulus*(2*x0)**(2*Material.b) + 1.04*Material.sigma_f*Material.epsilon_f*(2*x0)**(Material.b+Material.c)-delta_w_max
            ]
                
        fatigue_coefficient = delta_w_max
        result = fsolve(f, [1])
        fatigue_life = int(result[0])
        
        return self.outputFatigueLife(theta_critical_plane,
                                      phi_critical_plane,
                                      delta_gamma_critical_plane,
                                      sigma_nmax_critical_plane,
                                      delta_sigma_critical_plane,
                                      delta_epsilon_critical_plane,
                                      tau_nmax_critical_plane,
                                      delta_tau_critical_plane,
                                      temperature_at_sigma_nmax_critical_plane,
                                      fatigue_coefficient,
                                      fatigue_life)
                                      
    def initialValues(self):
        delta_gamma_max = 0.0
        delta_epsilon_max = 0.0
        tension_energy_max = 0.0
        shear_energy_max = 0.0
        total_energy_max = 0.0
        theta_critical_plane = 0.0
        phi_critical_plane = 0.0
        transformation_critical_plane = []

        if self.dimension == 3:
            theta_deg_list = range(0,180,1)
            phi_deg_list = range(0,180,1)
        if self.dimension == 2:
            theta_deg_list = [90]
            phi_deg_list = range(0,180,1)
            
        return (delta_gamma_max,
                delta_epsilon_max,
                tension_energy_max,
                shear_energy_max,
                total_energy_max,
                theta_critical_plane,
                phi_critical_plane,
                transformation_critical_plane,
                theta_deg_list,
                phi_deg_list)
                
    def calcTransformation(self,theta,phi):
        transformation=[
        [np.sin(theta)*np.cos(phi) ,np.sin(theta)*np.sin(phi) ,np.cos(theta)],
        [-np.sin(phi)              ,np.cos(phi)               ,0            ],
        [-np.cos(theta)*np.cos(phi),-np.cos(theta)*np.sin(phi),np.sin(theta)]
        ]
        transformation=[t[:self.dimension] for t in transformation[:self.dimension]]
        return transformation
    
    def calcValuesAtCriticalPlane(self,transformation_critical_plane):
        sigma_nmax_critical_plane = self.sigmaNMax(transformation_critical_plane)
        delta_sigma_critical_plane = self.deltaSigma(transformation_critical_plane)
        delta_epsilon_critical_plane = self.deltaEpsilon(transformation_critical_plane)
        tau_nmax_critical_plane = self.tauNMax(transformation_critical_plane)
        delta_tau_critical_plane = self.deltaTau(transformation_critical_plane)
        delta_gamma_critical_plane = self.deltaGamma(transformation_critical_plane)
        temperature_at_sigma_nmax_critical_plane = self.temperatureAtsigmaNMax(transformation_critical_plane)
        return (sigma_nmax_critical_plane,
                delta_sigma_critical_plane,
                delta_epsilon_critical_plane,
                tau_nmax_critical_plane,
                delta_tau_critical_plane,
                delta_gamma_critical_plane,
                temperature_at_sigma_nmax_critical_plane)
                
    def outputFatigueLife(self,
                          theta_critical_plane,
                          phi_critical_plane,
                          delta_gamma_critical_plane,
                          sigma_nmax_critical_plane,
                          delta_sigma_critical_plane,
                          delta_epsilon_critical_plane,
                          tau_nmax_critical_plane,
                          delta_tau_critical_plane,
                          temperature_at_sigma_nmax_critical_plane,
                          fatigue_coefficient,
                          fatigue_life):
        
        if self.dimension == 3:
            print 'theta_critical_plane',theta_critical_plane
            print 'phi_critical_plane',phi_critical_plane
        if self.dimension == 2:
            print 'phi_critical_plane',phi_critical_plane
            
        line_format = '%-40s'
        line_format_strain = line_format + '%.4f%%'
        line_format_stress = line_format + '%.2fMPa'
        line_format_coefficient = line_format + '%.6f'
        line_format_life = line_format + '%d'
        line_format_temperature = line_format + '%.1fC'
        print line_format_strain % ('delta_gamma_critical_plane',delta_gamma_critical_plane*100)
        print line_format_stress % ('sigma_nmax_critical_plane',sigma_nmax_critical_plane)
        print line_format_stress % ('delta_sigma_critical_plane',delta_sigma_critical_plane)
        print line_format_strain % ('delta_epsilon_critical_plane',delta_epsilon_critical_plane*100)
        print line_format_stress % ('tau_nmax_critical_plane',tau_nmax_critical_plane)
        print line_format_stress % ('delta_tau_critical_plane',delta_tau_critical_plane)
        print line_format_temperature % ('temperature_at_sigma_nmax_critical_plane',temperature_at_sigma_nmax_critical_plane)
        print line_format_coefficient % ('fatigue_coefficient',fatigue_coefficient)
        print line_format_life % ('fatigue_life',fatigue_life)
        
        if self.dimension == 2:
            return [phi_critical_plane,
                    sigma_nmax_critical_plane,
                    delta_sigma_critical_plane,
                    delta_epsilon_critical_plane,
                    tau_nmax_critical_plane,
                    delta_tau_critical_plane,
                    delta_gamma_critical_plane,
                    fatigue_life,
                    fatigue_coefficient,
                    temperature_at_sigma_nmax_critical_plane]
                
    def mathematicsTest(self):
        theta_deg=90.0
        phi_deg=0.0
        theta=np.radians(theta_deg)
        phi=np.radians(phi_deg)
        stress_list=[[2.0,1.0,0.0],[1.0,0.0,0.0],[0.0,0.0,0.0]]
        strain_list=[[2.0,1.0,0.0],[1.0,-1.0,0.0],[0.0,0.0,-1.0]]
        transformation_list=[
        [np.sin(theta)*np.cos(phi) ,np.sin(theta)*np.sin(phi) ,np.cos(theta)],
        [-np.sin(phi)              ,np.cos(phi)               ,0            ],
        [-np.cos(theta)*np.cos(phi),-np.cos(theta)*np.sin(phi),np.sin(theta)]
        ]
        
        dimension = 3
        stress_list = [s[:dimension] for s in stress_list[:dimension]]
        strain_list = [s[:dimension] for s in strain_list[:dimension]]
        transformation_list = [t[:dimension] for t in transformation_list[:dimension]]
        stress=np.array(stress_list)
        strain=np.array(strain_list)
        transformation=np.array(transformation_list)
        print stress
        print strain
        print 'theta = %s deg' % theta_deg
        print 'phi = %s deg' % phi_deg
        print transformation
        print np.dot(np.dot(transformation,stress),transformation.T)
        print np.dot(np.dot(transformation,strain),transformation.T)
        
        dimension = 2
        stress_list = [s[:dimension] for s in stress_list[:dimension]]
        strain_list = [s[:dimension] for s in strain_list[:dimension]]
        transformation_list = [t[:dimension] for t in transformation_list[:dimension]]
        stress=np.array(strain_list)
        strain=np.array(stress_list)
        transformation=np.array(transformation_list)
        print stress
        print strain
        print transformation
        print np.dot(np.dot(transformation,stress),transformation.T)
        print np.dot(np.dot(transformation,strain),transformation.T)
        
        s11 = stress[0][0]
        s22 = stress[1][1]
        s12 = stress[0][1]
        e11 = strain[0][0]
        e22 = strain[1][1]
        e12 = strain[0][1]*2
        sigma_phi      =s11/2+s11/2*np.cos(2*phi)+s12*np.sin(2*phi)
        tau_phi        =s11/2*np.sin(2*phi)-s12*np.cos(2*phi)
        sigma_phi      =(s11+s22)/2+(s11-s22)/2*np.cos(2*phi)+s12*np.sin(2*phi)
        tau_phi        =(s11-s22)/2*np.sin(2*phi)-s12*np.cos(2*phi)
        epsilon_phi    =(e11+e22)/2+(e11-e22)/2*np.cos(2*phi)+e12/2*np.sin(2*phi)
        gamma_phi_half =(e11-e22)/2*np.sin(2*phi)-e12/2*np.cos(2*phi)
        print sigma_phi, tau_phi
        print epsilon_phi, gamma_phi_half

        
    def lifeTest(self):
        youngs_modulus = 210000.0
        shear_modulus = 80800.0
        poisson_ratio = youngs_modulus/shear_modulus/2.0-1.0
        material = Material()
        material.setName(name='Test')
        material.setTemperature(temperature=20.0)
        material.setMonotonic(youngs_modulus=210000.0,poisson_ratio=poisson_ratio,yield_stress=573.0,K=1550.0,n=0.16)
        material.setCyclicAxial(sigma_f=1323.0,b=-0.097,epsilon_f=0.375,c=-0.60)
        material.setCyclicTorsion(tau_f=703.0,b0=-0.087,gamma_f=0.775,c0=-0.61)
        sigma_x = [-411,-281,-152,-23,106,222,304,376,437,486,521,521,514,473,430,389,374,245,116,-14,-123,-199,-263,-314,-357,-395,-411]
        epsilon_x = [0,0.00062,0.00123,0.00185,0.00246,0.00308,0.00369,0.00431,0.00492,0.00554,0.006,0.006,0.006,0.006,0.006,0.006,0.006,0.00538,0.00477,0.00415,0.00354,0.00292,0.00231,0.00169,0.00108,0.00046,0]
        epsilon_y = [-0.0004,-0.00058,-0.00076,-0.00094,-0.00113,-0.00133,-0.00156,-0.0018,-0.00205,-0.00231,-0.0025,-0.0025,-0.00251,-0.00255,-0.00259,-0.00263,-0.00264,-0.00246,-0.00227,-0.00209,-0.00189,-0.00165,-0.0014,-0.00115,-0.00088,-0.00061,-0.0004]
        tau = [20,41,66,91,116,132,131,132,132,131,130,41,-47,-117,-172,-215,-229,-204,-179,-154,-123,-87,-56,-31,-9,9,20]
        gamma = [0,0.00031,0.00062,0.00092,0.00123,0.00154,0.00185,0.00215,0.00246,0.00277,0.003,0.00189,0.00078,-0.00033,-0.00144,-0.00255,-0.003,-0.00269,-0.00238,-0.00208,-0.00177,-0.00146,-0.00115,-0.00085,-0.00054,-0.00023,0]
        length = len(sigma_x)
        time = range(1,length+1)
        temperature = range(1,length+1)
        s11 = np.array(sigma_x)
        s22 = np.zeros(length)
        s33 = np.zeros(length)
        s12 = np.array(tau)
        s13 = np.zeros(length)
        s23 = np.zeros(length)
        e11 = np.array(epsilon_x)
        e22 = np.array(epsilon_y)
        e33 = e11*-1.0*material.poisson_ratio
        e12 = np.array(gamma)/2.0
        e13 = np.zeros(length)
        e23 = np.zeros(length)
        stress=[]
        strain=[]
        for i in range(length):
            stress.append([[s11[i],s12[i],s13[i]],[s12[i],s22[i],s23[i]],[s13[i],s23[i],s33[i]]])
            strain.append([[e11[i],e12[i],e13[i]],[e12[i],e22[i],e23[i]],[e13[i],e23[i],e33[i]]])
        dimension = self.dimension
        self.time_list = [t for t in time]
        self.temperature_list = [t for t in temperature]
        self.strain_list = [[i[:dimension] for i in s[:dimension]] for s in strain]
        self.stress_list = [[i[:dimension] for i in s[:dimension]] for s in stress]
        
        self.fatigueLifeBMModel(material)
        self.fatigueLifeFSModel(material)
        self.fatigueLifeSWTModel(material)
        self.fatigueLifeLiu1Model(material)
        self.fatigueLifeLiu2Model(material)
        self.fatigueLifeChuModel(material)
        
n = Node(dimension=3)
n.lifeTest()
#n.mathematicsTest()