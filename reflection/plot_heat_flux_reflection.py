# -*- coding: utf-8 -*-
"""
Created on Tue Dec 18 13:49:06 2012

@author: Sun
"""

from mpl_toolkits.mplot3d.axes3d import Axes3D
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

x = [-width/2 + i*width/127.0 for i in range(0,128)]
y = [-length/2 + i*length/127.0 for i in range(0,128)]
x, y = np.meshgrid(x, y)

#===============================================================================
# Plot in x-y plane
#===============================================================================

fig = plt.figure(figsize=plt.figaspect(0.5))
ax = fig.add_subplot(1, 2, 1, projection='3d')
surf = ax.plot_surface(x, y, Liste4, rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0, antialiased=False)
#ax.set_zlim3d(-1.01, 1.01)
fig.colorbar(surf, shrink=0.5, aspect=10)

ax = fig.add_subplot(1, 2, 2, projection='3d')
ax.plot_wireframe(x, y, Liste4, rstride=10, cstride=10)

plt.legend(loc=0)
plt.show()