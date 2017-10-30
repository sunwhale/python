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

cc = 1
dp = 0.01
dimension = 3
c = np.zeros((dimension,dimension,dimension,dimension))
nn = np.zeros((dimension,dimension,dimension,dimension))
s = np.zeros((dimension,dimension))

#s[0][0] = '11'
#s[0][1] = '12'
#s[0][2] = '13'
#s[1][0] = '21'
#s[1][1] = '22'
#s[1][2] = '23'
#s[2][0] = '31'
#s[2][1] = '32'
#s[2][2] = '33'

s[0][0] = 1.0

#print s

count = 0

for i in range(dimension):
    for j in range(dimension):
        for k in range(dimension):
            for l in range(dimension):
                nn[i,j,k,l] = s[i][j] * s[k][l]

#for i in range(dimension):
#    for j in range(dimension):
#            print c[:,:,i,j]

dc = cc * (nn - c) * dp

print dc