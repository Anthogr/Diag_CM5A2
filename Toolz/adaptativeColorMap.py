#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 13:12:25 2021

@author: Anthony Gramoullé

Evaluate if data array analyzed in this function needs an adaptative colormap 
or not.

When plotting the normalized histogram of an array, each bin containing
an amount of data is characterized by pn(<1) and sum(pn) = 1.
If the difference between max(pn) and mean(pn) is too big (I chose arbitrarily 
0.5), it means data are (very) unevenly distributed between the bins and a 
classic regular colormap might not be very useful.

So:
IF max(pn) - mean(pn) > 0.5
    We add details in colormap in bins which contain the most data
ELSE
    We keep a standard regular colormap

Inputs:
    
    - var       : Array to evaluate  
    - p1/p2     : p values to define 
                  All bins with a pn < p1 won't get added levels in colormap
                  All bins with a p1 < pn < p2 will get more levels in colormap
                  All bins with a pn > p2 will get even more levels in colormap                   
    - cmapColor : Color of the Colormap ('viridis', 'magma', etc...)
 
    
Outputs:
    
    - cmap : colormap which is either regular or irregular depending of var 
             analyzed   
    - norm : Array which contain the values of the colormap (= None if 
             colormap is regular)
"""

def adaptativeColorMap(var, p1, p2, cmapColor):
    
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    
    var = var.flatten() # 2D array -> 1D array
    var = var.compressed() # remove mask values

    n, bins = np.histogram(var, bins=10, density=True) 

    prob = n * np.diff(bins)
    
    if prob.max()-prob.mean() > 0.2:
        # Adaptative colormap (the more data available in a bin, 
        # the more increments will be generated implemented in this bin)
        lvl=[]
        
        for i in np.arange(0,len(n)):
        
            if n[i] * np.diff(bins)[i]< p1: # probability of bin n°i
                
                lvl=np.append(lvl,np.linspace(bins[i],bins[i+1],2))[0:-1]
            
            elif (p1 < n[i] * np.diff(bins)[i]) & (n[i] * np.diff(bins)[i]< p2):
                
                lvl=np.append(lvl,np.linspace(bins[i],bins[i+1],7))[0:-1] # 7 or 10
                
            else:
                
                lvl=np.append(lvl,np.linspace(bins[i],bins[i+1],20))[0:-1] # 20 or 40
        
        
        bounds = np.round(lvl,3)
        
        # diverging colormap
        # if (bounds.min() < 0) & (bounds.max() > 0):
        #     if np.abs(bounds.min()) < np.abs(bounds.max()):
        #         bounds = np.insert(bounds, 0, - bounds.max())
        #     if np.abs(bounds.min()) > np.abs(bounds.max()):
        #         bounds = np.append(bounds, np.abs(bounds.min()))
                
        cmap = plt.get_cmap(cmapColor,len(bounds)-1)
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
        
    else:
        # cmap = cmapColor
        # norm = None

        cmap = plt.get_cmap(cmapColor,100)
        bounds = np.linspace(var.min(),var.max(),100)
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
            
    return cmap, norm