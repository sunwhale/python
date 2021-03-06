# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 13:46:18 2017

@author: SunJingyu
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import subprocess
import shutil
from Constants import *
from Data import SimulationData,ExperimentData,ExperimentLog
from Functions import write_umat_output_file,copy_suffix_files

#==============================================================================
# class Step
#==============================================================================
class Step:
    """
    Step类，定义求解时的一些参数。
    """
    def __init__(self, predefined_temperature = 300.0, time_period = 1000.0, 
                 initial_inc = 0.005, min_inc = 0.0001, max_inc = 5.5, nonlinear = 'ON'):
        self.predefined_temperature = predefined_temperature
        self.time_period = time_period
        self.initial_inc = initial_inc
        self.min_inc = min_inc
        self.max_inc = max_inc
        self.nonlinear = nonlinear

#==============================================================================
# class UMAT  
#==============================================================================
class UMAT:
    """
    UMAT类，定义UMAT文件的路径，名称。
    """
    def __init__(self, UMATDirectory = 'F:\\GitHub\\umat\\', 
                 UMATMainFile = 'MAIN_IN718.for', 
                 ParameterFortranFile = 'PARAMETERS_IN718_TMF.for',
                 OutputFortranFile = 'OUTPUT.for',
                 OutputTextFile = 'output.txt'):
        self.UMATDirectory = UMATDirectory
        self.UMATMainFile = UMATMainFile
        self.ParameterFortranFile = ParameterFortranFile
        self.OutputFortranFile = OutputFortranFile
        self.OutputTextFile = OutputTextFile
        self.UMATMainFileFullName = UMATDirectory + UMATMainFile
        self.ParameterFortranFileFullName = UMATDirectory + ParameterFortranFile
        self.OutputFortranFileFullName = UMATDirectory + OutputFortranFile
        
#==============================================================================
# class Load
#==============================================================================
class Load:
    """
    Load类，定义载荷。
    """
    def __init__(self, runing_time=[], temperature=[], axial_strain=[], 
                 shear_strain=[], axial_stress=[], torque=[], first_cycle_shift=0):
        self.runing_time = runing_time
        self.temperature = temperature
        self.axial_strain = axial_strain
        self.axial_stress = axial_stress
        self.shear_strain = shear_strain
        self.torque = torque
        self.total_runing_time = self.runing_time[-1]
        self.length = len(self.runing_time)
        self.first_cycle_shift = first_cycle_shift
        self.segment = 1
        self.film_coefficient = 0.0
        self.sink_temperature = 0.0

    def setLoadFromExperiment(self, ExperimentData, runing_time=None):
        self.runing_time = ExperimentData.runing_time
        self.temperature = ExperimentData.temperature
        self.axial_strain = ExperimentData.axial_strain
        self.axial_stress = ExperimentData.axial_stress
        self.shear_strain = ExperimentData.shear_strain
        self.torque = ExperimentData.torque
        if runing_time == None:
            self.total_runing_time = self.runing_time[-1]
        else:
            self.total_runing_time = runing_time
        self.length = len(self.runing_time)
        self.segment = 2
        
    def setLoadUniaxial(self, cycles, runing_time, temperature, axial_strain):
        for n in range(cycles):
            self.runing_time += [ n*(runing_time[-1]-runing_time[0]) + time for time in runing_time[:-1]]
        self.temperature = temperature[:-1] * cycles
        self.axial_strain = axial_strain[:-1] * cycles
        
        self.runing_time.append(cycles*(runing_time[-1]-runing_time[0]))
        self.temperature.append(temperature[-1])
        self.axial_strain.append(axial_strain[-1])
        self.length = len(self.runing_time)
        self.axial_stress = np.zeros(self.length)
        self.torque = np.zeros(self.length)
        self.total_runing_time = self.runing_time[-1]
        self.segment = 1
        
    def listToArray(self):
        self.runing_time = np.array(self.runing_time)
        self.temperature = np.array(self.temperature)
        self.axial_strain = np.array(self.axial_strain)
        self.axial_stress = np.array(self.axial_stress)
        self.shear_strain = np.array(self.shear_strain)
        self.torque = np.array(self.torque)
        
    def setLoadBiaxial(self, cycles, runing_time, temperature, axial_strain, shear_strain):
        for n in range(cycles):
            self.runing_time += [ n*(runing_time[-1]-runing_time[0]) + time for time in runing_time[:-1]]
        self.temperature += temperature[:-1] * cycles
        self.axial_strain += axial_strain[:-1] * cycles
        self.shear_strain += shear_strain[:-1] * cycles
        
        self.runing_time.append(cycles*(runing_time[-1]-runing_time[0]))
        self.temperature.append(temperature[-1])
        self.axial_strain.append(axial_strain[-1])
        self.shear_strain.append(shear_strain[-1])
        self.total_runing_time = self.runing_time[-1]

        if axial_strain[0]==0 and shear_strain[0]==0:
            self.runing_time.pop(1)
            self.temperature.pop(1)
            self.axial_strain.pop(1)
            self.shear_strain.pop(1)
        if axial_strain[0]<>0:
            self.runing_time[1] = self.first_cycle_shift
        if shear_strain[0]<>0:
            self.runing_time[1] = self.first_cycle_shift
        self.length = len(self.runing_time)
        self.axial_stress = np.zeros(self.length)
        self.torque = np.zeros(self.length)
        self.segment = 1
        
    def setConvection(self, film_coefficient, sink_temperature):
        self.film_coefficient = film_coefficient
        self.sink_temperature = sink_temperature
    
    def setThermal(self, heat_flux, film_coefficient_outer, film_coefficient_inner, 
                   emissivity, sink_temperature_inner ,sink_temperature_outer, outer_temperature):
        self.heat_flux = heat_flux
        self.film_coefficient_outer = film_coefficient_outer
        self.film_coefficient_inner = film_coefficient_inner
        self.emissivity = emissivity
        self.sink_temperature_inner = sink_temperature_inner
        self.sink_temperature_outer = sink_temperature_outer
        self.outer_temperature = outer_temperature

        
    def showLoadPath(self, directory):
        os.chdir(directory)
        show_number = 12
        save_types = ['png','pdf']
        self.listToArray()
        plt.figure(figsize=(8,6))
        ax1 = plt.gca()
        ax1.plot(self.runing_time[:show_number], self.axial_strain[:show_number], label='axial_strain')
        ax1.plot(self.runing_time[:show_number], self.shear_strain[:show_number]/np.sqrt(3), label='shear_strain_eq')
        ax1.legend(loc=1,fontsize='medium',frameon=True,numpoints=1,title='')
        ax1.set_xlabel(xylabels['runing_time'])
        ax1.set_ylabel(xylabels['axial_strain'] + ' and ' +  xylabels['shear_strain_eq'])
        ax2 = ax1.twinx() # this is the important function 
        ax2.plot(self.runing_time[:show_number], self.temperature[:show_number], label='temperature', color='red')
        ax2.set_ylabel(xylabels['temperature'])
        ax2.legend(loc=4,fontsize='medium',frameon=True,numpoints=1,title='')
        for save_type in save_types:
            save_name = '%s_%s_%s.%s' %('Loadpath','strain','temperature',save_type)
            plt.savefig(save_name, dpi=150)
        plt.close()
        
        plt.figure(figsize=(8,6))
        plt.xlabel(xylabels['runing_time'])
        plt.ylabel(xylabels['axial_strain'] + ' and ' +  xylabels['shear_strain_eq'])
        plt.plot(self.runing_time[:show_number], self.axial_strain[:show_number])
        plt.plot(self.runing_time[:show_number], np.array(self.shear_strain[:show_number])*-6.0/np.sqrt(3))
        for save_type in save_types:
            save_name = '%s_%s_%s.%s' %('Loadpath','time','strain',save_type)
            plt.savefig(save_name, dpi=150)
        plt.close()
        
        plt.figure(figsize=(8,6))
        plt.xlabel(xylabels['runing_time'])
        plt.ylabel(xylabels['temperature'])
        plt.plot(self.runing_time[:show_number], self.temperature[:show_number])
        for save_type in save_types:
            save_name = '%s_%s_%s.%s' %('Loadpath','time','temperature',save_type)
            plt.savefig(save_name, dpi=150)
        plt.close()
        
        plt.figure(figsize=(8,6))
        plt.xlabel(xylabels['axial_strain'])
        plt.ylabel(xylabels['shear_strain_eq'])
        plt.plot(self.axial_strain[:show_number], self.shear_strain[:show_number]/np.sqrt(3))
        plt.gca().set_aspect('equal')
        for save_type in save_types:
            save_name = '%s_%s_%s.%s' %('Loadpath','axial_strain','shear_strain_eq',save_type)
            plt.savefig(save_name, dpi=150)
        plt.close()

#==============================================================================
# class Job
#==============================================================================
class Job:
    """
    Job类，定义ABAQUS前处理，求解，后处理。
    """
    def __init__(self, JobName, UMAT, Step, Load, copy=True):
        self.JobName = JobName
        self.UMAT = UMAT
        self.Step = Step
        self.Load = Load
        self.copy = copy
        self.CAEName = JobName + '.cae'
        self.InputName = JobName + '.inp'
        self.AbaqusWorkDirectory = AbaqusTempDirectory + JobName + '\\'
        self.AbaqusWorkUMATDirectory = AbaqusTempDirectory + JobName + '\\UMAT\\'
        self.PythonPostProc = 'PostprocABAQUS_thermal.py'
        self.PythonPreproc = 'PreprocABAQUS_thermal.py'
        self.PythonConstants = 'Constants.py'
        self.PythonPreprocParameters = 'PreprocABAQUSParameters.py'
        self.ExpCSVFile = []
        
        self.CSVFullName = self.AbaqusWorkDirectory + JobName + '.csv'
        self.JobFullName = self.AbaqusWorkDirectory + JobName
        
    def allProc(self):
        self.createDirectory()
        self.copyFiles()
        self.creatBatchFile()
        self.createAbaqusCAE()
        self.createAbaqusInput()
#        self.createUMATFile()
        self.run()
        self.autoPostProc()
        
    def createDirectory(self):
        if not os.path.isdir(self.AbaqusWorkDirectory):
            os.makedirs(self.AbaqusWorkDirectory)
            print 'Create new directory:',self.AbaqusWorkDirectory
    
    def copyFiles(self):
        shutil.copy(PythonDirectory + self.PythonPostProc, self.AbaqusWorkDirectory)
        shutil.copy(PythonDirectory + self.PythonPreproc, self.AbaqusWorkDirectory)
        shutil.copy(PythonDirectory + self.PythonConstants, self.AbaqusWorkDirectory)

    def createAbaqusCAE(self):
        self.createDirectory()
        outfile = open(self.AbaqusWorkDirectory + self.PythonPreprocParameters,'w')
        print >>outfile, '# -*- coding: utf-8 -*-'
        print >>outfile, 'from abaqus import *'
        print >>outfile, 'from abaqusConstants import *'
        print >>outfile, 'predefined_temperature = %s' % (self.Step.predefined_temperature)
        print >>outfile, 'time_period = %s' % (self.Step.time_period)
        print >>outfile, 'initial_inc = %s' % (self.Step.initial_inc)
        print >>outfile, 'min_inc = %s' % (self.Step.min_inc)
        print >>outfile, 'max_inc = %s' % (self.Step.max_inc)
        print >>outfile, 'nonlinear = %s' % (self.Step.nonlinear)
        print >>outfile, 'JobName = %r' % (self.JobName)
        print >>outfile, 'CAEName = %r' % (self.CAEName)
        print >>outfile, 'AbaqusWorkDirectory = %r' % (self.AbaqusWorkDirectory)
        print >>outfile, 'film_coefficient = %r' % (self.Load.film_coefficient)
        print >>outfile, 'sink_temperature = %r' % (self.Load.sink_temperature)
        print >>outfile, 'heat_flux = %r' % (self.Load.heat_flux)
        print >>outfile, 'emissivity = %r' % (self.Load.emissivity)
        print >>outfile, 'sink_temperature_inner = %r' % (self.Load.sink_temperature_inner)
        print >>outfile, 'sink_temperature_outer = %r' % (self.Load.sink_temperature_outer)
        print >>outfile, 'film_coefficient_outer = %r' % (self.Load.film_coefficient_outer)
        print >>outfile, 'film_coefficient_inner = %r' % (self.Load.film_coefficient_inner)
        print >>outfile, 'outer_temperature = %r' % (self.Load.outer_temperature)    
        outfile.close()

        shutil.copy(PythonDirectory + self.PythonPreproc, self.AbaqusWorkDirectory)
        shutil.copy(PythonDirectory + self.PythonConstants, self.AbaqusWorkDirectory)
        
        os.chdir(self.AbaqusWorkDirectory)
        self.creatBatchFile()
        cmd ='AutoPreProc'
        subprocess.call(cmd,shell=True)
        
    def createAbaqusInput(self):
        self.Load.showLoadPath(self.AbaqusWorkDirectory)
        Name1 = 'amplitude_displacement.txt'
        Name2 = 'amplitude_temperature.txt'
        Name3 = 'amplitude_pressure.txt'
        Name4 = 'amplitude_rotation.txt'
        Name5 = 'amplitude_torque.txt'
        #==============================================================================
        # *Amplitude, name=DispTriangularWave
        #==============================================================================
        output_filename = self.AbaqusWorkDirectory + Name1
        outfile = open(output_filename, 'w')
        outfile.writelines('*Amplitude, name=DispTriangularWave\n')
        for i in range(0,self.Load.length,self.Load.segment):
            line = '%-20.10f,%-20.10f\n' % (self.Load.runing_time[i],self.Load.axial_strain[i])
            outfile.writelines(line)
        print 'Create ', output_filename
        outfile.close()
        #==============================================================================
        # *Amplitude, name=TempTriangularWave
        #==============================================================================
        output_filename = self.AbaqusWorkDirectory + Name2
        outfile = open(output_filename, 'w')
        outfile.writelines('*Amplitude, name=TempTriangularWave\n')
        for i in range(0,self.Load.length,self.Load.segment):
            line = '%-20.10f,%-20.10f\n' % (self.Load.runing_time[i],self.Load.temperature[i])
            outfile.writelines(line)
        print 'Create ', output_filename
        outfile.close()
        #==============================================================================
        # *Amplitude, name=PressureTriangularWave
        #==============================================================================
        output_filename = self.AbaqusWorkDirectory + Name3
        outfile = open(output_filename, 'w')
        outfile.writelines('*Amplitude, name=PressureTriangularWave\n')
        for i in range(0,self.Load.length,self.Load.segment):
            line = '%-20.10f,%-20.10f\n' % (self.Load.runing_time[i],self.Load.axial_stress[i])
            outfile.writelines(line)
        print 'Create ', output_filename
        outfile.close()
        #==============================================================================
        # *Amplitude, name=RotaTriangularWave
        #==============================================================================
        output_filename = self.AbaqusWorkDirectory + Name4
        outfile = open(output_filename, 'w')
        outfile.writelines('*Amplitude, name=RotaTriangularWave\n')
        for i in range(0,self.Load.length,self.Load.segment):
            line = '%-20.10f,%-20.10f\n' % (self.Load.runing_time[i],self.Load.shear_strain[i]/-6.0)
            outfile.writelines(line)
        print 'Create ', output_filename
        outfile.close()
        #==============================================================================
        # *Amplitude, name=TorqTriangularWave
        #==============================================================================
        output_filename = self.AbaqusWorkDirectory + Name5
        outfile = open(output_filename, 'w')
        outfile.writelines('*Amplitude, name=TorqTriangularWave\n')
        for i in range(0,self.Load.length,self.Load.segment):
            line = '%-20.10f,%-20.10f\n' % (self.Load.runing_time[i],self.Load.torque[i])
            outfile.writelines(line)
        print 'Create ', output_filename
        
        infile = open(self.AbaqusWorkDirectory+self.InputName,'r')
        list1 = infile.readlines()
        infile.close()
        
        outfile = open(self.AbaqusWorkDirectory+self.InputName, 'w')
        i=0
        while i<len(list1):
            if '*Amplitude, name=DispTriangularWave' in list1[i]:
                line = '*Include, input=%s\n' % (Name1)
                i=i+1
                print line[:-1]
            elif '*Amplitude, name=TempTriangularWave' in list1[i]:
                line = '*Include, input=%s\n' % (Name2)
                i=i+1
                print line[:-1]
            elif '*Amplitude, name=PresTriangularWave' in list1[i]:
                line = '*Include, input=%s\n' % (Name3)
                i=i+1
                print line[:-1]
            elif '*Amplitude, name=RotaTriangularWave' in list1[i]:
                line = '*Include, input=%s\n' % (Name4)
                i=i+1
                print line[:-1]
            elif '*Amplitude, name=TorqTriangularWave' in list1[i]:
                line = '*Include, input=%s\n' % (Name5)
                i=i+1
                print line[:-1]
            else:
                line = list1[i]
            i=i+1
            outfile.writelines(line)
        outfile.close()
    
    def creatBatchFile(self):
        outfile = open(self.AbaqusWorkDirectory + 'Run.bat', 'w')
        print >>outfile, (r'call set JobName=%s' % self.JobName)
        print >>outfile, (r'call set UMAT=%s' % self.UMAT.UMATDirectory + self.UMAT.UMATMainFile) # 使用原始UMAT文件夹
        print >>outfile, (r'call set UMAT=%s' % self.AbaqusWorkUMATDirectory + self.UMAT.UMATMainFile) # 使用工作目录下UMAT文件夹
        print >>outfile, r'call del *.rpy.*'
        print >>outfile, r'call del *.rpy'
        print >>outfile, r'call del %JobName%.com'
        print >>outfile, r'call del %JobName%.lck'
        print >>outfile, r'call del %JobName%.dat'
        print >>outfile, r'call del %JobName%.msg'
        print >>outfile, r'call del %JobName%.odb'
        print >>outfile, r'call del %JobName%.odb_f'
        print >>outfile, r'call del %JobName%.prt'
        print >>outfile, r'call del %JobName%.sim'
        print >>outfile, r'call del %JobName%.sta'
#        print >>outfile, r'call abaqus job=%JobName% user=%UMAT% cpus=1 int'
        print >>outfile, r'call abaqus job=%JobName% cpus=1 int'
        outfile.close()

        outfile = open(self.AbaqusWorkDirectory + 'AutoPostProc.bat', 'w')
        print >>outfile, (r'call abaqus viewer noGUI=%s' % self.PythonPostProc)
        outfile.close()
        
        outfile = open(self.AbaqusWorkDirectory + 'AutoPreProc.bat', 'w')
        print >>outfile, (r'abaqus cae noGUI=%s' % self.PythonPreproc)
        outfile.close()
        
        outfile = open(self.AbaqusWorkDirectory + 'OpenCommandline.bat', 'w')
        print >>outfile, r'call cmd'
        outfile.close()
        
    def run(self):
        self.creatBatchFile()
        os.chdir(self.AbaqusWorkDirectory)
        cmd = 'Run'
        subprocess.call(cmd,shell=True)

    def createUMATFile(self):
        fortran_file_name = self.UMAT.UMATDirectory + self.UMAT.OutputFortranFile
        text_file_name = self.AbaqusWorkDirectory + self.UMAT.OutputTextFile
        total_time = self.Load.total_runing_time
        write_umat_output_file(fortran_file_name,text_file_name,total_time)
        old = self.AbaqusWorkUMATDirectory
        if os.path.isdir(old):
            for i in range(2,100):
                new = self.AbaqusWorkUMATDirectory.replace('UMAT','UMAT'+str(i))
                if not os.path.isdir(new):
                    os.rename(old,new)
                    shutil.copytree(self.UMAT.UMATDirectory, self.AbaqusWorkUMATDirectory)
                    break
        else:
            shutil.copytree(self.UMAT.UMATDirectory, self.AbaqusWorkUMATDirectory)
            
    def autoPostProc(self):
        self.creatBatchFile()
        shutil.copy(PythonDirectory + self.PythonPostProc, self.AbaqusWorkDirectory)
        os.chdir(self.AbaqusWorkDirectory)
        cmd ='AutoPostProc'
        subprocess.call(cmd,shell=True)
        if self.copy:
            copy_suffix_files(self.AbaqusWorkDirectory,SimulationDirectory,suffixs=['csv'])
        
    def show(self):
        print 'PythonDirectory:',PythonDirectory
        print 'InputDirectory:',InputDirectory
        print 'AbaqusWorkDirectory:',self.AbaqusWorkDirectory
        print 'JobName:',self.JobName
        print 'CAEName:',self.CAEName
        print 'InputName:',self.InputName
        print 'PythonPostProc:',self.PythonPostProc
        print 'PythonPreproc:',self.PythonPreproc
        print 'ExpCSVFile:',self.ExpCSVFile