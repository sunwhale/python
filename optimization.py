# -*- coding: utf-8 -*-
"""
Created on Thu Jan 05 11:50:47 2017

@author: j.Sun
"""

from scipy.optimize import fmin
from workbench_thermal_optimization import calc_temperature


volume_flow = 50.0
power_percent_exp_list = [0.3,0.4,0.5,0.6]

#print calc_temperature(volume_flow,power_percent_exp_list)

exp_steady = [
[67, [[0.3, 215], [0.4, 338], [0.5, 495], [0.6, 654]]], 
[50, [[0.3, 277], [0.4, 435], [0.5, 603], [0.6, 757.5]]], 
[25, [[0.3, 440], [0.4, 602.5], [0.5, 759.5], [0.6, 899]]], 
[0, [[0.3, 607], [0.4, 770], [0.5, 892], [0.6, 1015]]]
]

outfile = open('optimization.txt', 'w')

def func(x):
    volume_flow = x[0]
    sim = calc_temperature(volume_flow,power_percent_exp_list)
    exp = exp_steady[0][1]
    y = 0.0
    for i in range(4):
        y += (exp[i][1] - sim[i][1])**2
    print >>outfile, sim
    print >>outfile, x, y
    return y


if __name__ == '__main__':
    x0 = [volume_flow]    #猜一个初值 
    xopt = fmin(func, x0, maxiter=30)    #求解
    print xopt
    print >>outfile, xopt    #打印结果
    outfile.close()