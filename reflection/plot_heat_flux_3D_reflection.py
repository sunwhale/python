# -*- coding: utf-8 -*-
"""
Created on Tue Dec 18 13:49:06 2012

@author: Sun
"""

from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt


# imports specific to the plots in this example
import numpy as np
from matplotlib import cm
from mpl_toolkits.mplot3d.axes3d import get_test_data


length = 100.0
width = 80.0
txt="1.txt"
Liste1 = open(txt,"r").readlines()	
Liste2 = [i.strip() for i in Liste1]
Liste3 = [j.split('\t') for j in Liste2]
Liste4 = [[float(j) for j in i] for i in Liste3]

X = [-width/2 + i*width/127.0 for i in range(0,128)]
Y = [-length/2 + i*length/127.0 for i in range(0,128)]
X, Y = np.meshgrid(X, Y)
Z = Liste4

print Z[61][62:66]
print Z[62][62:66]
print Z[63][62:66]
print Z[64][62:66]
print Z[65][62:66]
print Z[66][62:66]

#===============================================================================
# Plot in x-y plane
#===============================================================================
fig = plt.figure()
ax = fig.gca(projection='3d')

surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0, antialiased=False)
fig.colorbar(surf, shrink=0.5, aspect=10)

#ax.plot_surface(X, Y, Z, rstride=4, cstride=4, alpha=0.2)
cset = ax.contour(X, Y, Z, zdir='z', offset=-200000, cmap=cm.coolwarm)
cset = ax.contour(X, Y, Z, zdir='x', offset=-50, cmap=cm.coolwarm)
cset = ax.contour(X, Y, Z, zdir='y', offset=50, cmap=cm.coolwarm)
ax.set_zlim3d(-200000, 400000)
ax.set_xlabel('X')
ax.set_xlim(-50, 50)
ax.set_ylabel('Y')
ax.set_ylim(-50, 50)
ax.set_zlabel('Z')
ax.set_zlim(-200000, 400000)

#plt.show()