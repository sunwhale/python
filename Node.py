# -*- coding: utf-8 -*-
"""
Created on Thu Jan 05 11:50:47 2017

@author: j.Sun
"""

import numpy as np
from scipy.optimize import fsolve
from Material import Material,material_in718
from Data import ExperimentData
from Functions import calculate_conductivity_by_temperature_in718,calculate_elastic_by_temperature_in718
from Constants import *

class Node:
    def __init__(self, nodelabel=1, dimension=2, time=[], coordinate=[], 
                 displacement=[], stress=[], strain=[], temperature=[], heatflux=[],
                 theta_interval=5, phi_interval=5):
        self.nodelabel = nodelabel
        self.dimension = dimension
        self.time_list = [t for t in time]
        self.temperature_list = [t for t in temperature]
        self.heatflux_list = [[i for i in h[:dimension]] for h in heatflux]
        self.coordinate_list = [[i for i in c[:dimension]] for c in coordinate]
        self.displacement_list = [[i for i in d[:dimension]] for d in displacement]
        self.strain_list = [[i[:dimension] for i in s[:dimension]] for s in strain]
        self.stress_list = [[i[:dimension] for i in s[:dimension]] for s in stress]
        self.theta_interval = theta_interval
        self.phi_interval = phi_interval
        self.transformation_critical_plane = self.calcTransformation(0,0)
        self.TMFCefficient()
        
    def __str__(self):
        return 'node %s at time %s' % (str(int(self.nodelabel)),str(self.time))

    def outputJson(self):
        return {'nodelabel':self.nodelabel,
                'dimension':self.dimension,
                'time_list':self.time_list,
                'temperature_list':self.temperature_list,
                'heatflux_list':self.heatflux_list,
                'coordinate_list':self.coordinate_list,
                'displacement_list':self.displacement_list,
                'strain_list':self.strain_list,
                'stress_list':self.stress_list,
                'theta_interval':self.theta_interval,
                'phi_interval':self.phi_interval
                }
                
    def inputJson(self,node,interval=1):
        interval = int(interval)
        self.nodelabel = node['nodelabel']
        self.dimension = node['dimension']
        self.time_list = node['time_list'][::interval]
        self.temperature_list = node['temperature_list'][::interval]
        self.heatflux_list = node['heatflux_list'][::interval]
        self.coordinate_list = node['coordinate_list'][::interval]
        self.displacement_list = node['displacement_list'][::interval]
        self.strain_list = node['strain_list'][::interval]
        self.stress_list = node['stress_list'][::interval]
        self.theta_interval = node['theta_interval']
        self.phi_interval = node['phi_interval']
        
    def stressTriaxiality(self):
        if self.dimension == 2:
            stress_tensor = [[i[:self.dimension]+[0.0] for i in s[:self.dimension]] + [[0.0,0.0,0.0]] for s in self.stress_list]
        if self.dimension == 3:
            stress_tensor = [[i[:self.dimension] for i in s[:self.dimension]] for s in self.stress_list]
        stress_tensor_deviatoric = []
        stress_hydrostatic = []
        stress_mises = []
        stress_triaxiality = []
#        stress_tensor = [[[1,0,0],[0,0,0],[0,0,0]]]
        for s in stress_tensor:
            hydrostatic = np.trace(s)/3.0
            s[0][0] -= hydrostatic
            s[1][1] -= hydrostatic
            s[2][2] -= hydrostatic
            stress_hydrostatic.append(hydrostatic)
            stress_tensor_deviatoric.append(s)
            J2 = 0.0
            for i in range(3):
                for j in range(3):
                    J2 += s[i][j]**2/2.0
            mises = np.sqrt(3.0*J2)
            stress_mises.append(mises)
            stress_triaxiality.append(hydrostatic/mises)
        return stress_hydrostatic,stress_mises,stress_triaxiality

    def strainMises(self):
        if self.dimension == 2:
            strain_tensor = [[i[:self.dimension]+[0.0] for i in s[:self.dimension]] + [[0.0,0.0,0.0]] for s in self.strain_list]
        if self.dimension == 3:
            strain_tensor = [[i[:self.dimension] for i in s[:self.dimension]] for s in self.strain_list]
#        strain_tensor = [[[1,0,0],[0,-0.5,0],[0,0,-0.5]]]
        strain_mises = []
        for s in strain_tensor:
            mises = np.sqrt(s[0][0]**2 + (2*s[0][1])**2/3.0)
            strain_mises.append(mises)
        return strain_mises
        
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
        
    def tensionIntegralEnergy(self,transformation):
        epsilon = self.normalStrain(transformation)
        sigma = self.normalStress(transformation)
        integral = 0.0
        for i in range(len(epsilon)-1):
            integral += (sigma[i]+sigma[i+1])/2*(epsilon[i+1]-epsilon[i])
        temperature = 475.0
        youngs_modulus,poisson_ratio,shear_modulus,yield_stress = calculate_elastic_by_temperature_in718(temperature)
        (max(sigma)-min(sigma))**2/youngs_modulus
        return integral + (max(sigma)-min(sigma))**2/youngs_modulus*2.0
        
    def shearIntegralEnergy(self,transformation):
        gamma = self.shearStrain(transformation)
        tau = self.shearStress(transformation)
        integral = 0.0
        for i in range(len(gamma)-1):
            integral += (tau[i][0]+tau[i+1][0])/2*(gamma[i+1][0]-gamma[i][0]) + (tau[i][1]+tau[i+1][1])/2*(gamma[i+1][1]-gamma[i][1])
        temperature = 475.0
        youngs_modulus,poisson_ratio,shear_modulus,yield_stress = calculate_elastic_by_temperature_in718(temperature)
        tau = np.array(tau)
        return integral + (max(tau[:,0])-min(tau[:,0]))**2/shear_modulus*2.0 + (max(tau[:,1])-min(tau[:,1]))**2/shear_modulus*2.0
        
    def totalIntegralEnergy(self):
        transformation = self.calcTransformation(np.radians(90),np.radians(0))
        transform_stress_list = self.stressTransform(transformation)
        transform_strain_list = self.strainTransform(transformation)
        integral = 0.0
        for i in range(len(transform_stress_list)-1):
            for j in range(self.dimension):
                for k in range(self.dimension):
                    integral += abs((transform_stress_list[i][j][k]+transform_stress_list[i+1][j][k])/2*(transform_strain_list[i+1][j][k]-transform_strain_list[i][j][k]))

        return integral

    def temperatureAtsigmaNMax(self, transformation):
        normal_stress_list = self.normalStress(transformation)
        max_normal_stress = max(normal_stress_list)
        max_normal_stress_index = normal_stress_list.index(max_normal_stress)
        temperature_at_sigma_nmax = self.temperature_list[max_normal_stress_index]
        return temperature_at_sigma_nmax
        
    def heatfluxAtsigmaNMax(self, transformation):
        normal_stress_list = self.normalStress(transformation)
        max_normal_stress = max(normal_stress_list)
        max_normal_stress_index = normal_stress_list.index(max_normal_stress)
        heatflux_at_sigma_nmax = self.heatflux_list[max_normal_stress_index]
        return heatflux_at_sigma_nmax
    
    def outputValuesAtCriticalPlane(self, transformation):
        normal_stress_list = self.normalStress(transformation)
        return [self.time_list,normal_stress_list,self.temperature_list,self.heatflux_list]
    
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
         temperature_at_sigma_nmax_critical_plane,
         heatflux_at_sigma_nmax_critical_plane) = self.calcValuesAtCriticalPlane(transformation_critical_plane)
        
        thermal_corrected_term = self.thermalCorrectedTerm(heatflux_at_sigma_nmax_critical_plane,temperature_at_sigma_nmax_critical_plane)
        thermal_corrected_term = 1.0
        tmf_corrected_term = self.TMFCorrectedTerm(transformation_critical_plane)
        tmf_corrected_term = 1.0
        
        fs_coefficient=delta_gamma_max/2.0*(1.0+k*sigma_nmax_critical_plane*tmf_corrected_term/Material.yield_stress)
        
        fs_coefficient *= thermal_corrected_term
        
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
                                      heatflux_at_sigma_nmax_critical_plane,
                                      fatigue_coefficient,
                                      fatigue_life)

    def fatigueLifeShowModel(self,Material,k=0.3):
        print '=========================Show model========================='

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
        for i in range(0,len(self.time_list),10):
            print self.time_list[i],self.strain_list[i][0][0],self.strain_list[i][0][1],self.strain_list[i][1][1],self.stress_list[i][0][0],self.stress_list[i][0][1],self.temperature_list[i]
        for theta_deg in theta_deg_list:
            theta=np.radians(theta_deg)
            for phi_deg in phi_deg_list:
                phi=np.radians(phi_deg)
                transformation = self.calcTransformation(theta,phi)

                sigma_nmax_theta_phi = self.sigmaNMax(transformation)
                delta_sigma_theta_phi = self.deltaSigma(transformation)
                delta_epsilon_theta_phi = self.deltaEpsilon(transformation)
                tau_nmax_theta_phi = self.tauNMax(transformation)
                delta_tau_theta_phi = self.deltaTau(transformation)
                delta_gamma_theta_phi = self.deltaGamma(transformation)
                temperature_at_sigma_nmax_theta_phi = self.temperatureAtsigmaNMax(transformation)
                heatflux_at_sigma_nmax_theta_phi = self.heatfluxAtsigmaNMax(transformation)

                print phi_deg,sigma_nmax_theta_phi,delta_sigma_theta_phi,delta_epsilon_theta_phi,tau_nmax_theta_phi,delta_tau_theta_phi,delta_gamma_theta_phi,temperature_at_sigma_nmax_theta_phi
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
         temperature_at_sigma_nmax_critical_plane,
         heatflux_at_sigma_nmax_critical_plane) = self.calcValuesAtCriticalPlane(transformation_critical_plane)
        
        thermal_corrected_term = self.thermalCorrectedTerm(heatflux_at_sigma_nmax_critical_plane,temperature_at_sigma_nmax_critical_plane)
        thermal_corrected_term = 1.0
        tmf_corrected_term = self.TMFCorrectedTerm(transformation_critical_plane)
        tmf_corrected_term = 1.0
        
        fs_coefficient=delta_gamma_max/2.0*(1.0+k*sigma_nmax_critical_plane*tmf_corrected_term/Material.yield_stress)
        
        fs_coefficient *= thermal_corrected_term
        
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
                                      heatflux_at_sigma_nmax_critical_plane,
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
         temperature_at_sigma_nmax_critical_plane,
         heatflux_at_sigma_nmax_critical_plane) = self.calcValuesAtCriticalPlane(transformation_critical_plane)
        
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
                                      heatflux_at_sigma_nmax_critical_plane,
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
         temperature_at_sigma_nmax_critical_plane,
         heatflux_at_sigma_nmax_critical_plane) = self.calcValuesAtCriticalPlane(transformation_critical_plane)
        
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
                                      heatflux_at_sigma_nmax_critical_plane,
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
         temperature_at_sigma_nmax_critical_plane,
         heatflux_at_sigma_nmax_critical_plane) = self.calcValuesAtCriticalPlane(transformation_critical_plane)
        
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
                                      heatflux_at_sigma_nmax_critical_plane,
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
         temperature_at_sigma_nmax_critical_plane,
         heatflux_at_sigma_nmax_critical_plane) = self.calcValuesAtCriticalPlane(transformation_critical_plane)
        
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
                                      heatflux_at_sigma_nmax_critical_plane,
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
         temperature_at_sigma_nmax_critical_plane,
         heatflux_at_sigma_nmax_critical_plane) = self.calcValuesAtCriticalPlane(transformation_critical_plane)
        
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
                                      heatflux_at_sigma_nmax_critical_plane,
                                      fatigue_coefficient,
                                      fatigue_life)

    def fatigueLifeZamrikModel(self,Material):
        print '=========================Zamrik model========================='

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
            theta=np.radians(90)
            for phi_deg in phi_deg_list:
                phi=np.radians(0)
                transformation = self.calcTransformation(theta,phi)
                theta_critical_plane = theta_deg
                phi_critical_plane = phi_deg
                transformation_critical_plane = transformation
                    
        (sigma_nmax_critical_plane,
         delta_sigma_critical_plane,
         delta_epsilon_critical_plane,
         tau_nmax_critical_plane,
         delta_tau_critical_plane,
         delta_gamma_critical_plane,
         temperature_at_sigma_nmax_critical_plane,
         heatflux_at_sigma_nmax_critical_plane) = self.calcValuesAtCriticalPlane(transformation_critical_plane)
         
        A = 10**(-3.40627845)
        B = -3.68266722
        epsilon_fracture = 0.234
        sigma_ultimate = 1305.0
        Q = 240.0
        R = 8.31e-3
        T_0 = 300 + 273.15
        
        stress_hydrostatic,stress_mises,stress_triaxiality = self.stressTriaxiality()
        sigma_max = max(stress_mises)
        print sigma_max
        strain_mises = self.strainMises()
        epsilon_tension = (max(strain_mises) - min(strain_mises))/2
        print epsilon_tension
        T_max = max(self.temperature_list) + 273.15
        T_max = 650.0 + 273.15
        
        fatigue_life = int(A*((epsilon_tension/epsilon_fracture)*(sigma_max/sigma_ultimate))**B)
        
        print fatigue_life
#        print T_max
#        print np.exp(-1.0*Q/R/(T_max-T_0))
        
        fatigue_coefficient = (epsilon_tension/epsilon_fracture)*(sigma_max/sigma_ultimate)
        
        return self.outputFatigueLife(theta_critical_plane,
                                      phi_critical_plane,
                                      delta_gamma_critical_plane,
                                      sigma_nmax_critical_plane,
                                      delta_sigma_critical_plane,
                                      delta_epsilon_critical_plane,
                                      tau_nmax_critical_plane,
                                      delta_tau_critical_plane,
                                      temperature_at_sigma_nmax_critical_plane,
                                      heatflux_at_sigma_nmax_critical_plane,
                                      fatigue_coefficient,
                                      fatigue_life)
                                      
    def fatigueLifeVoseModel(self,Material):
        print '=========================Vose model========================='

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
            theta=np.radians(90)
            for phi_deg in phi_deg_list:
                phi=np.radians(0)
                transformation = self.calcTransformation(theta,phi)
                theta_critical_plane = theta_deg
                phi_critical_plane = phi_deg
                transformation_critical_plane = transformation
                    
        (sigma_nmax_critical_plane,
         delta_sigma_critical_plane,
         delta_epsilon_critical_plane,
         tau_nmax_critical_plane,
         delta_tau_critical_plane,
         delta_gamma_critical_plane,
         temperature_at_sigma_nmax_critical_plane,
         heatflux_at_sigma_nmax_critical_plane) = self.calcValuesAtCriticalPlane(transformation_critical_plane)

        Q0=240.0
        upsilon0=3.50e-4
        sigma_ult=1305.0
        R=8.31e-3
        C5=1.28710938e+12
        k=1.2
        C1,C2,C3,C4 = [  9.94816905e-02,  -2.30395385e-01,   1.79726262e-07,   7.85958724e-01]
#        C1,C2,C3,C4 = [ 0.54659446, -0.04706091, -0.4971955,  -0.04087346]
        stress_ratio = (sigma_nmax_critical_plane-delta_sigma_critical_plane)/sigma_nmax_critical_plane
        m = 1

        integral = 0.0
        for i in range(1,len(self.stress_list)-1):
            sigma_alt = delta_sigma_critical_plane/2.0
            T = self.temperature_list[i] + 273.15
            integral += np.exp(-1.0*(Q0-upsilon0*sigma_alt*(1.0-0.5*sigma_alt/sigma_ult))/(R*T)) * (self.time_list[i+1]-self.time_list[i])
        
        fatigue_coefficient = delta_epsilon_critical_plane * 2.0**(1-m)*(1.0-stress_ratio)**(m-1) * (1.0+C5*integral)**k

        def f(x):
            x0 = float(x[0])
            return [
#                Material.sigma_f/Material.youngs_modulus*(2*x0)**(Material.b) + Material.epsilon_f*(2*x0)**(Material.c)-fatigue_coefficient
                C1*(2*x0)**(C2) + C3*(2*x0)**(C4)-fatigue_coefficient
            ]
        
        result = fsolve(f, [1])
        fatigue_life = int(result[0])
        
#        C1 = 10**(-3.40627845)
#        C2 = -3.68266722
                
#        fatigue_life = int(C1*fatigue_coefficient**C2)
        
        print fatigue_life   
        
        return self.outputFatigueLife(theta_critical_plane,
                                      phi_critical_plane,
                                      delta_gamma_critical_plane,
                                      sigma_nmax_critical_plane,
                                      delta_sigma_critical_plane,
                                      delta_epsilon_critical_plane,
                                      tau_nmax_critical_plane,
                                      delta_tau_critical_plane,
                                      temperature_at_sigma_nmax_critical_plane,
                                      heatflux_at_sigma_nmax_critical_plane,
                                      fatigue_coefficient,
                                      fatigue_life)

    def fatigueLifeStudyModel(self,Material):
        print '=========================Study model========================='

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
         temperature_at_sigma_nmax_critical_plane,
         heatflux_at_sigma_nmax_critical_plane) = self.calcValuesAtCriticalPlane(transformation_critical_plane)
        
        tmf_corrected_term = self.TMFCorrectedTerm(transformation_critical_plane)
#        delta_w1_max = self.tensionIntegralEnergy(transformation_critical_plane)*tmf_corrected_term + self.shearIntegralEnergy(transformation_critical_plane)
#        delta_w1_max *= 1.732
        delta_w1_max = self.tensionEnergy(transformation_critical_plane)*tmf_corrected_term + self.shearEnergy(transformation_critical_plane)*0.5
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
                                      heatflux_at_sigma_nmax_critical_plane,
                                      fatigue_coefficient,
                                      fatigue_life)

    def fatigueLifeStudy2Model(self,Material):
        print '=========================Study2 model========================='

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
         temperature_at_sigma_nmax_critical_plane,
         heatflux_at_sigma_nmax_critical_plane) = self.calcValuesAtCriticalPlane(transformation_critical_plane)
        
        tmf_corrected_term = self.TMFCorrectedTerm(transformation_critical_plane)
#        tmf_corrected_term = 1.0
        tgmf_corrected_term = self.TGMFCorrectedTerm(heatflux_at_sigma_nmax_critical_plane,temperature_at_sigma_nmax_critical_plane)
#        delta_w1_max = self.tensionIntegralEnergy(transformation_critical_plane)*tmf_corrected_term + self.shearIntegralEnergy(transformation_critical_plane)
#        delta_w1_max *= 1.732
        delta_w1_max = self.tensionEnergy(transformation_critical_plane)*tmf_corrected_term + self.shearEnergy(transformation_critical_plane)*0.5
        stress_ratio = (sigma_nmax_critical_plane-delta_sigma_critical_plane)/sigma_nmax_critical_plane
        m = 0
        delta_w1_max *= 2.0**(1-m)*(1.0-stress_ratio)**(m-1)
        print 2.0/(1.0-stress_ratio)
        print tmf_corrected_term
        print tgmf_corrected_term
        delta_w1_max *= tgmf_corrected_term
        
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
                                      heatflux_at_sigma_nmax_critical_plane,
                                      fatigue_coefficient,
                                      fatigue_life)
                                      
    def fatigueLifeOurModel(self,Material):
        print '=========================Our model========================='

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
         temperature_at_sigma_nmax_critical_plane,
         heatflux_at_sigma_nmax_critical_plane) = self.calcValuesAtCriticalPlane(transformation_critical_plane)

        thermal_corrected_term = self.thermalCorrectedTerm(heatflux_at_sigma_nmax_critical_plane,temperature_at_sigma_nmax_critical_plane)
        
        delta_w_max = sigma_nmax_critical_plane*delta_epsilon_critical_plane/2.0 + tau_nmax_critical_plane*delta_gamma_critical_plane/2.0
        
        delta_w_max *= thermal_corrected_term
        
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
                                      heatflux_at_sigma_nmax_critical_plane,
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
            theta_deg_list = range(0,180,self.theta_interval)
            phi_deg_list = range(0,180,self.phi_interval)
        if self.dimension == 2:
            theta_deg_list = [90]
            phi_deg_list = range(0,180,self.phi_interval)
            
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
        
    def thermalCorrectedTerm(self,heatflux_at_sigma_nmax_critical_plane,temperature_at_sigma_nmax_critical_plane):
        hf_norm = 0.0
        for hf in heatflux_at_sigma_nmax_critical_plane:
            hf_norm += (hf/calculate_conductivity_by_temperature_in718(temperature_at_sigma_nmax_critical_plane))**2
        thermal_corrected_term = 1.0 + hf_norm/(750.0-temperature_at_sigma_nmax_critical_plane)
        return thermal_corrected_term

    def TGMFCorrectedTerm(self,heatflux_at_sigma_nmax_critical_plane,temperature_at_sigma_nmax_critical_plane):
        hf_norm = 0.0
        for hf in heatflux_at_sigma_nmax_critical_plane:
            hf_norm += (hf/calculate_conductivity_by_temperature_in718(temperature_at_sigma_nmax_critical_plane))**2
        tgmf_corrected_term = 1.0 + 1.0*hf_norm/(850.0-temperature_at_sigma_nmax_critical_plane)
        return tgmf_corrected_term
        
    def TMFCefficient(self,Q0=240,upsilon0=3.50e-4,sigma_ult=1305.0,R=8.31e-3,C=1.28710938e+12,k=1.2):
        self.Q0 = Q0
        self.upsilon0 = upsilon0
        self.sigma_ult = sigma_ult
        self.R = R
        self.C = C
        self.C = 1.0e+12
        self.k = k
        
    def TMFCorrectedTerm(self,transformation):
        Q0 = self.Q0
        upsilon0 = self.upsilon0
        sigma_ult = self.sigma_ult
        R = self.R
        C = self.C
        k = self.k
        
        stress_hydrostatic,stress_mises,stress_triaxiality = self.stressTriaxiality()
        
        integral = 0.0
        for i in range(1,len(self.stress_list)-1):
            sigma_alt = self.normalStress(transformation)[i] - self.shearStress(transformation)[i][0]
            sigma_alt = abs(sigma_alt)
            triaxiality = abs(stress_triaxiality[i])
            triaxiality = stress_triaxiality[i]
#            triaxiality = 1.0
            T = self.temperature_list[i] + 273.15
            integral += triaxiality * np.exp(-1.0*(Q0-upsilon0*sigma_alt*(1.0-0.5*sigma_alt/sigma_ult))/(R*T)) * (self.time_list[i+1]-self.time_list[i])
#            if sigma_alt >= 0:
#                integral += np.exp(-1.0*(Q0-nu0*sigma_alt*abs(triaxiality)*(1.0-0.5*sigma_alt/sigma_ult))/(R*T)) * (self.time_list[i+1]-self.time_list[i])
#            if sigma_alt <= 0:
#                sigma_alt = abs(sigma_alt)
#                integral -= np.exp(-1.0*(Q0-nu0*sigma_alt*abs(triaxiality)*(1.0-0.5*sigma_alt/sigma_ult))/(R*T)) * (self.time_list[i+1]-self.time_list[i])

#        print (1.0+C*integral)**k
        return (1.0+C*integral)**k
        
    def calcValuesAtCriticalPlane(self,transformation_critical_plane):
        self.transformation_critical_plane = transformation_critical_plane
        sigma_nmax_critical_plane = self.sigmaNMax(transformation_critical_plane)
        delta_sigma_critical_plane = self.deltaSigma(transformation_critical_plane)
        delta_epsilon_critical_plane = self.deltaEpsilon(transformation_critical_plane)
        tau_nmax_critical_plane = self.tauNMax(transformation_critical_plane)
        delta_tau_critical_plane = self.deltaTau(transformation_critical_plane)
        delta_gamma_critical_plane = self.deltaGamma(transformation_critical_plane)
        temperature_at_sigma_nmax_critical_plane = self.temperatureAtsigmaNMax(transformation_critical_plane)
        heatflux_at_sigma_nmax_critical_plane = self.heatfluxAtsigmaNMax(transformation_critical_plane)
        return (sigma_nmax_critical_plane,
                delta_sigma_critical_plane,
                delta_epsilon_critical_plane,
                tau_nmax_critical_plane,
                delta_tau_critical_plane,
                delta_gamma_critical_plane,
                temperature_at_sigma_nmax_critical_plane,
                heatflux_at_sigma_nmax_critical_plane)
                
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
                          heatflux_at_sigma_nmax_critical_plane,
                          fatigue_coefficient,
                          fatigue_life):
        
#        if self.dimension == 3:
#            print 'theta_critical_plane',theta_critical_plane
#        print 'phi_critical_plane',phi_critical_plane
#            
#        line_format = '%-40s'
#        line_format_strain = line_format + '%.4f%%'
#        line_format_stress = line_format + '%.2fMPa'
#        line_format_coefficient = line_format + '%.6f'
#        line_format_life = line_format + '%d'
#        line_format_temperature = line_format + '%.1fC'
#        print line_format_strain % ('delta_gamma_critical_plane',delta_gamma_critical_plane*100)
#        print line_format_stress % ('sigma_nmax_critical_plane',sigma_nmax_critical_plane)
#        print line_format_stress % ('delta_sigma_critical_plane',delta_sigma_critical_plane)
#        print line_format_strain % ('delta_epsilon_critical_plane',delta_epsilon_critical_plane*100)
#        print line_format_stress % ('tau_nmax_critical_plane',tau_nmax_critical_plane)
#        print line_format_stress % ('delta_tau_critical_plane',delta_tau_critical_plane)
#        print line_format_temperature % ('temperature_at_sigma_nmax_critical_plane',temperature_at_sigma_nmax_critical_plane)
#        print line_format_coefficient % ('fatigue_coefficient',fatigue_coefficient)
#        print line_format_life % ('fatigue_life',fatigue_life)
        
        temperature_gradient_at_sigma_nmax_critical_plane_x = heatflux_at_sigma_nmax_critical_plane[0]/calculate_conductivity_by_temperature_in718(temperature_at_sigma_nmax_critical_plane)
        temperature_gradient_at_sigma_nmax_critical_plane_y = heatflux_at_sigma_nmax_critical_plane[1]/calculate_conductivity_by_temperature_in718(temperature_at_sigma_nmax_critical_plane)
        
        if self.dimension == 3:
            temperature_gradient_at_sigma_nmax_critical_plane_z = heatflux_at_sigma_nmax_critical_plane[2]/calculate_conductivity_by_temperature_in718(temperature_at_sigma_nmax_critical_plane)                    

        if self.dimension == 3:
            return [theta_critical_plane,
                    phi_critical_plane,
                    sigma_nmax_critical_plane,
                    delta_sigma_critical_plane,
                    delta_epsilon_critical_plane,
                    tau_nmax_critical_plane,
                    delta_tau_critical_plane,
                    delta_gamma_critical_plane,
                    fatigue_life,
                    fatigue_coefficient,
                    temperature_at_sigma_nmax_critical_plane,
                    temperature_gradient_at_sigma_nmax_critical_plane_x,
                    temperature_gradient_at_sigma_nmax_critical_plane_y,
                    temperature_gradient_at_sigma_nmax_critical_plane_z]
                    
        return [theta_critical_plane,
                phi_critical_plane,
                sigma_nmax_critical_plane,
                delta_sigma_critical_plane,
                delta_epsilon_critical_plane,
                tau_nmax_critical_plane,
                delta_tau_critical_plane,
                delta_gamma_critical_plane,
                fatigue_life,
                fatigue_coefficient,
                temperature_at_sigma_nmax_critical_plane,
                temperature_gradient_at_sigma_nmax_critical_plane_x,
                temperature_gradient_at_sigma_nmax_critical_plane_y]
                
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
        print 'stress'
        print stress
        print 'strain'
        print strain
        print 'theta = %s deg' % theta_deg
        print 'phi = %s deg' % phi_deg
        print 'transformation'
        print transformation
#        print np.dot(np.dot(transformation,stress),transformation.T)
#        print np.dot(np.dot(transformation,strain),transformation.T)
        
        dimension = 2
        stress_list = [s[:dimension] for s in stress_list[:dimension]]
        strain_list = [s[:dimension] for s in strain_list[:dimension]]
        transformation_list = [t[:dimension] for t in transformation_list[:dimension]]
        stress=np.array(strain_list)
        strain=np.array(stress_list)
        transformation=np.array(transformation_list)
        print 'stress'
        print stress
        print 'strain'
        print strain
        print 'transformation'
        print transformation
        print 'stress_transformation'
        print np.dot(np.dot(transformation,stress),transformation.T)
        print 'strain_transformation'
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
        heat_flux_1 = np.zeros(length)
        heat_flux_2 = np.zeros(length)
        self.heatflux_list = []
        for i in range(length):
            self.heatflux_list.append([heat_flux_1[i],heat_flux_2[i]])
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
        
#        self.fatigueLifeBMModel(material)
#        self.fatigueLifeFSModel(material)
#        self.fatigueLifeSWTModel(material)
#        self.fatigueLifeLiu1Model(material)
#        self.fatigueLifeLiu2Model(material)
#        self.fatigueLifeChuModel(material)
        
        transformation = self.calcTransformation(np.radians(90),np.radians(0))
#        print self.tensionIntegralEnergy(transformation)
#        print self.shearIntegralEnergy(transformation)
#        print self.tensionIntegralEnergy(transformation) + self.shearIntegralEnergy(transformation)
#        print self.totalIntegralEnergy()
        self.stressTriaxiality()
        self.strainEquivalent()

def calculate_data_fatigue_life(data,material):
    nodelabel = data.node_label[0]
    nth = data.axial_count_index_list[-2]
    nth = data.half_life_cycle
    time = data.obtainNthCycle('runing_time',nth)
    temperature = data.obtainNthCycle('temperature',nth)
    heat_flux_1 = data.obtainNthCycle('heat_flux_1',nth)
    heat_flux_2 = data.obtainNthCycle('heat_flux_2',nth)
    heat_flux_3 = data.obtainNthCycle('heat_flux_3',nth)
    length = len(time)
    s11 = data.obtainNthCycle('axial_stress',nth)
    s22 = np.zeros(length)
    s33 = np.zeros(length)
    s12 = data.obtainNthCycle('shear_stress',nth)
    s13 = np.zeros(length)
    s23 = np.zeros(length)
    e11 = data.obtainNthCycle('axial_strain',nth)
    e22 = e11*-1.0*material.poisson_ratio
    e33 = e11*-1.0*material.poisson_ratio
    e12 = data.obtainNthCycle('shear_strain',nth)/2.0
    e13 = np.zeros(length)
    e23 = np.zeros(length)
    if heat_flux_1 == []:
        heat_flux_1 = np.zeros(length)
    if heat_flux_2 == []:
        heat_flux_2 = np.zeros(length)
    if heat_flux_3 == []:
        heat_flux_3 = np.zeros(length)
    stress = []
    strain = []
    heatflux = []
    for i in range(len(time)):
        stress.append([[s11[i],s12[i],s13[i]],[s12[i],s22[i],s23[i]],[s13[i],s23[i],s33[i]]])
        strain.append([[e11[i],e12[i],e13[i]],[e12[i],e22[i],e23[i]],[e13[i],e23[i],e33[i]]])
        heatflux.append([heat_flux_1[i],heat_flux_2[i],heat_flux_3[i]])
    node = Node(nodelabel=nodelabel, dimension=2, time=time, coordinate=[], 
                displacement=[], stress=stress, strain=strain,
                temperature=temperature, heatflux=heatflux)
    node.TMFCorrectedTerm(node.calcTransformation(np.radians(90),0))
                
if __name__ == '__main__':
    n = Node(dimension=2)
#    n.lifeTest()
#
#    for name in ['7114','7018','7017','7025','7037']:
#    for name in ['7114']: