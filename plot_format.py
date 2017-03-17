# -*- coding: utf-8 -*-
"""
Created on Thu Jan 05 11:50:47 2017

@author: j.Sun
"""

import matplotlib as mpl

def plot_format():
#==============================================================================
# http://matplotlib.org/users/customizing.html?highlight=rcparams
#==============================================================================
    mpl.rcParams['axes.linewidth'] = 1.5
    
    mpl.rcParams['xtick.major.size'] = 8
    mpl.rcParams['xtick.major.width'] = 1.5
    mpl.rcParams['xtick.minor.size'] = 4
    mpl.rcParams['xtick.minor.width'] = 1.5
    mpl.rcParams['xtick.labelsize'] = 16
    
    mpl.rcParams['ytick.major.size'] = 8
    mpl.rcParams['ytick.major.width'] = 1.5
    mpl.rcParams['ytick.minor.size'] = 4
    mpl.rcParams['ytick.minor.width'] = 1.5
    mpl.rcParams['ytick.labelsize'] = 16
    
    mpl.rcParams['figure.figsize'] = (8,6)
    mpl.rcParams['axes.labelsize'] = 20  # fontsize of the x any y labels
    mpl.rcParams['figure.subplot.left'] = 0.125
    mpl.rcParams['figure.subplot.right'] = 0.95
    mpl.rcParams['figure.subplot.bottom'] = 0.125
    mpl.rcParams['figure.subplot.top'] = 0.95
    
    mpl.rcParams['legend.fontsize'] = 16
    mpl.rcParams['legend.frameon'] = False
    mpl.rcParams['legend.numpoints'] = 1
    mpl.rcParams['legend.scatterpoints'] = 1