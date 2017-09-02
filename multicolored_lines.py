# -*- coding: utf-8 -*-
"""
Created on Wed Mar 01 10:35:05 2017

@author: Z620
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.collections as mcoll

def multicolored_lines(x=None, y=None, z=None):
    """
    http://nbviewer.ipython.org/github/dpsanders/matplotlib-examples/blob/master/colorline.ipynb
    http://matplotlib.org/examples/pylab_examples/multicolored_line.html
    """
    if x is None:
        x = np.linspace(0, 4. * np.pi, 100)
    if y is None:
        y = np.sin(x)
    if z is None:
        z = np.linspace(0.0, 1.0, len(x))
    fig, ax = plt.subplots()
    plt.xlim(x.min(), x.max())
    plt.ylim(y.min(), y.max())
    plt.xlim(-1500, 1500)
    plt.ylim(-800, 800)
    plt.xlim(-1.0, 1)
    plt.ylim(-1.8, 1.8)
    mymap = matplotlib.colors.LinearSegmentedColormap.from_list('mycolors',['blue','green','yellow','orange','red'])
    mymap = matplotlib.colors.LinearSegmentedColormap.from_list('mycolors',['blue','red'])
    lc = colorline(x, y, z, norm=plt.Normalize(300, 650), cmap=mymap)
    cbar = plt.colorbar(lc)
    cbar.set_label('Temperature [$^{\circ}$C]')
    cbar.set_ticks(np.linspace(300,650,8))
    cbar.set_ticklabels( ('300', '350', '400', '450', '500', '550', '600', '650') )
    

def colorline(
        x, y, z=None, cmap='copper', norm=plt.Normalize(0.0, 1.0),
        linewidth=3, alpha=1.0):
    """
    http://nbviewer.ipython.org/github/dpsanders/matplotlib-examples/blob/master/colorline.ipynb
    http://matplotlib.org/examples/pylab_examples/multicolored_line.html
    Plot a colored line with coordinates x and y
    Optionally specify colors in the array z
    Optionally specify a colormap, a norm function and a line width
    """

    # Default colors equally spaced on [0,1]:
    if z is None:
        z = np.linspace(0.0, 1.0, len(x))

    # Special case if a single number:
    # to check for numerical input -- this is a hack
    if not hasattr(z, "__iter__"):
        z = np.array([z])

    z = np.asarray(z)

    segments = make_segments(x, y)
    lc = mcoll.LineCollection(segments, array=z, cmap=cmap, norm=norm,
                              linewidth=linewidth, alpha=alpha)

    ax = plt.gca()
    ax.add_collection(lc)

    return lc

def make_segments(x, y):
    """
    Create list of line segments from x and y coordinates, in the correct format
    for LineCollection: an array of the form numlines x (points per line) x 2 (x
    and y) array
    """

    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    return segments