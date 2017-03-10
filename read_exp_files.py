# -*- coding: utf-8 -*-
"""
Created on Fri Mar 03 14:50:05 2017

@author: Z620
"""

import numpy as np

def read_exp_files(filename, n=1):
    infile = open(filename,'r')
    list1 = infile.readlines()
    infile.close()
    list2 = [i.strip() for i in list1[::n]] #每n行取一次数据
    del list1
    list3 = [i.split(',') for i in list2]
    del list2
    for i in list3:
        for j in range(len(i)):
            if i[j]=='':
                i[j]='0'
    list4 = [[float(i[0]),float(i[1]),float(i[2]),float(i[3]),float(i[4]),float(i[5]),float(i[6]),float(i[7]),float(i[8]),float(i[9]),float(i[10]),float(i[11]),float(i[12]),float(i[13])] for i in list3[2:] ]
    del list3
    return np.array(list4)