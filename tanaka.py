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

def calc_AijklBkl(A,B,dimension=3):
    C = np.zeros((dimension,dimension))
    for i in range(dimension):
        for j in range(dimension):
            for k in range(dimension):
                for l in range(dimension):
                    C[i,j] = C[i,j] + A[i,j,k,l] * B[k,l]
    return C

def calc_AijklBijkl(A,B,dimension=3):
    C = 0.0
    for i in range(dimension):
        for j in range(dimension):
            for k in range(dimension):
                for l in range(dimension):
                    C = C + A[i,j,k,l] * B[i,j,k,l]
    return C

def calc_AijBij(A,B,dimension=3):
    C = 0.0
    for i in range(dimension):
        for j in range(dimension):
            C = C + A[i,j] * B[i,j]
    return C

def calc_AijBji(A,B,dimension=3):
    C = 0.0
    for i in range(dimension):
        for j in range(dimension):
            C = C + A[i,j] * B[j,i]
    return C

def calc_AijBkl(A,B,dimension=3):
    C = np.zeros((dimension,dimension,dimension,dimension))
    for i in range(dimension):
        for j in range(dimension):
            for k in range(dimension):
                for l in range(dimension):
                    C[i,j,k,l] = A[i,j] * B[k,l]
    return C

def calc_Mean(A,dimension=3):
    Mean = 0.0
    for i in range(dimension):
        Mean += A[i,i]
    Mean *= 1.0/3.0
    return Mean

def calc_Diviatoric(A,dimension=3):
    D = A - calc_Mean(A,dimension)
    return D
    
def calc_Mises(A,dimension=3):
    D = calc_Diviatoric(A,dimension)
    Mises = np.sqrt(1.5*calc_AijBij(D,D,dimension))
    return Mises


    
def calc_Tanaka(sigma_x,tau,dimension=3,dp=0.01):
    coef = 1
    s = np.zeros((dimension,dimension))
    C = np.zeros((dimension,dimension,dimension,dimension))
    nn = np.zeros((dimension,dimension,dimension,dimension))
    phi_list = []
    for i in range(len(sigma_x)):
        s11 = sigma_x[i]
        s22 = 0.0
        s33 = 0.0
        s12 = tau[i]
        s13 = 0.0
        s23 = 0.0
        s21 = s12
        s31 = s13
        s32 = s23
    
        s[0][0] = s11
        s[0][1] = s12
        s[0][2] = s13
        s[1][0] = s21
        s[1][1] = s22
        s[1][2] = s23
        s[2][0] = s31
        s[2][1] = s32
        s[2][2] = s33
    
        n = s/np.linalg.norm(s)
        nn = calc_AijBkl(n,n,dimension)
        dC = coef * (nn - C) * dp
        C = C + dC
        """NCIJ=C_ijkl:n_kl"""
        NCIJ = calc_AijklBkl(C,n,dimension)
        """NCCN=nc_ij:nc_ji"""
        NCCN = calc_AijBji(NCIJ,NCIJ,dimension)
        """CC=c_ijkl:c_ijkl"""
        CC = calc_AijklBijkl(C,C,dimension)
        phi = np.sqrt(abs(1-NCCN/CC))
        
        phi_list.append(phi)

    return phi_list

def calc_Triaxiality(sigma_x,tau,dimension=3,dp=0.01):
    coef = 1
    s = np.zeros((dimension,dimension))
    C = np.zeros((dimension,dimension,dimension,dimension))
    nn = np.zeros((dimension,dimension,dimension,dimension))
    triaxiality_list = []
    for i in range(len(sigma_x)):
        s11 = sigma_x[i]
        s22 = 0.0
        s33 = 0.0
        s12 = tau[i]
        s13 = 0.0
        s23 = 0.0
        s21 = s12
        s31 = s13
        s32 = s23
    
        s[0][0] = s11
        s[0][1] = s12
        s[0][2] = s13
        s[1][0] = s21
        s[1][1] = s22
        s[1][2] = s23
        s[2][0] = s31
        s[2][1] = s32
        s[2][2] = s33
    
        triaxiality = calc_Mean(s,dimension)/calc_Mises(s,dimension)
        
        triaxiality_list.append(triaxiality)

    return triaxiality_list

def calc_User(sigma_x,tau,dimension=3,dp=0.01):
    coef = 1
    s = np.zeros((dimension,dimension))
    C = np.zeros((dimension,dimension,dimension,dimension))
    nn = np.zeros((dimension,dimension,dimension,dimension))
    phi_list = []
    for i in range(len(sigma_x)):
        s11 = sigma_x[i]
        s22 = 0.0
        s33 = 0.0
        s12 = tau[i]
        s13 = 0.0
        s23 = 0.0
        s21 = s12
        s31 = s13
        s32 = s23
    
        s[0][0] = s11
        s[0][1] = s12
        s[0][2] = s13
        s[1][0] = s21
        s[1][1] = s22
        s[1][2] = s23
        s[2][0] = s31
        s[2][1] = s32
        s[2][2] = s33
        
        n = s/np.linalg.norm(s)
        nn = calc_AijBkl(n,n,dimension)
        dC = coef * (nn - C) * dp
        C = C + dC
        """NCIJ=C_ijkl:n_kl"""
        NCIJ = calc_AijklBkl(C,n,dimension)
        """NCCN=nc_ij:nc_ji"""
        NCCN = calc_AijBji(NCIJ,NCIJ,dimension)
        """CC=c_ijkl:c_ijkl"""
        CC = calc_AijklBijkl(C,C,dimension)
        nnnn = calc_AijklBijkl(nn,nn,dimension)
        phi = 1-NCCN/CC
        print CC
        phi_list.append(CC)

    return phi_list

sigma_x = [-411,-281,-152,-23,106,222,304,376,437,486,521,521,514,473,430,389,374,245,116,-14,-123,-199,-263,-314,-357,-395,-411]
tau = [20,41,66,91,116,132,131,132,132,131,130,41,-47,-117,-172,-215,-229,-204,-179,-154,-123,-87,-56,-31,-9,9,20]

""" Uniaxial path """
#sigma_x = []
#tau = []
#for i in np.arange(-100,100,1):
#    sigma_x.append(i)
#    tau.append(i)
    
""" Proportional path """
#sigma_x = []
#tau = []
#for i in np.arange(1,100,1):
#    sigma_x.append(i)
#    tau.append(i)
    
""" Circle path """
sigma_x = []
tau = []
for i in np.arange(0,np.pi*100,0.01):
    sigma_x.append(np.sin(i))
    tau.append(np.cos(i+np.pi/4))

""" L path """
#sigma_x = []
#tau = []
#for i in np.arange(0,1000,0.1):
#    sigma_x.append(i)
#    tau.append(5)

phi_list = calc_Tanaka(sigma_x,tau,dimension=3,dp=0.001)
triaxiality_list = calc_Triaxiality(sigma_x,tau,dimension=3,dp=0.001)
user_list = calc_User(sigma_x,tau,dimension=3,dp=0.001)

#plt.plot(sigma_x,tau)
#plt.show()

plt.plot(sigma_x,triaxiality_list,label='Triaxiality')
plt.plot(sigma_x,phi_list,label='Tanaka')
plt.plot(sigma_x,[1,2],label='User')
lg = plt.legend(title='',loc=0)
plt.show()











#s[0][0] = '11'
#s[0][1] = '12'
#s[0][2] = '13'
#s[1][0] = '21'
#s[1][1] = '22'
#s[1][2] = '23'
#s[2][0] = '31'
#s[2][1] = '32'
#s[2][2] = '33'

#for i in range(dimension):
#    for j in range(dimension):
#            print c[:,:,i,j]