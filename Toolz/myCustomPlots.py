#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 12:00:33 2021

@author: anthony
"""

def plotMapping(lon, lat, var2plot, contourLines, land_mask,
                inProjData, outProjData, cmap, norm, figTitle, figXsize, figYsize, 
                cbar_title, cbar_label_size, cbar_tick_size, title_font_size):
    
    import numpy as np
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    from Toolz import (dispGridCoastline, z_masked_overlap)
    
    fig, ax = plt.subplots(figsize=(figXsize, figYsize), subplot_kw={'projection': outProjData})

    # Display meridians and parallels
    ax.gridlines(draw_labels = True) # transform = ccrs.Geodetic(),
    
    # Choose map as global 
    # ax.set_global()
    
    # Computes lon and lat grids into projected cartesian coordinate X and Y
    # as well as maskedZ in order to be in PlateCarree standard projection
    X, Y, maskedZ = z_masked_overlap(ax, lon, lat, var2plot,       
                                     source_projection = ccrs.Geodetic())                                   
    
    # Pcolormesh
    map1 = ax.pcolormesh(lon, lat, maskedZ[0:-1,0:-1],                   
                         transform = inProjData, cmap = cmap, norm = norm, shading='nearest')
    # Contour
    if type(contourLines) is np.ndarray:
        cont1 = ax.contour(X, Y, maskedZ, contourLines,                         
                            transform = outProjData,
                            colors='k', linewidths=0.7)
        # Add labels over contour lines
        ax.clabel(cont1, fmt=' {:.1f} '.format, fontsize='x-large')  
    
    # Compute coastline
    #dispGridCoastline(lon, lat, inProjData, land_mask[0,:,:], 1.25)   

    # Display colorbar
    if cbar_title is not None:
        if len(np.unique(np.round(np.diff(norm.boundaries),2))) == 1:
            ticks = np.linspace(norm.vmin, norm.vmax, 7, endpoint=True) # control nbb of ticks for colorbar
            cbar = plt.colorbar(map1, ticks=ticks, orientation='horizontal', extend='both')
            
        elif len(np.unique(np.round(np.diff(norm.boundaries),2))) != 1:
            plt.text(0.5,-0.15,'Warning: Adaptative colormap (non-linear) !',horizontalalignment='center',
                     verticalalignment='center', transform = ax.transAxes, fontsize=16, color='r', weight='bold') 
            cbar = plt.colorbar(map1, orientation='horizontal', extend='both')
            
        cbar.ax.set_title(cbar_title,size=cbar_label_size)
        cbar.ax.tick_params(labelsize=cbar_tick_size)   
    
    # Titles/Labels and size
    plt.title(figTitle,fontsize=title_font_size)
    

def plotMappingSubPlots(nrows, ncols, lon, lat, vars2plot, contourLines, land_mask,
                inProjData, outProjData, cmap, norm, figTitle, figXsize, figYsize, 
                cbar_title, cbar_label_size, cbar_tick_size, title_font_size):
    
    import numpy as np
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    from Toolz import (dispGridCoastline, z_masked_overlap)

    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(figXsize, figYsize), subplot_kw={'projection': outProjData})
    
    for i in np.arange(0,nrows*ncols):

        plt.axes(ax[i])

        # Display meridians and parallels
        ax[i].gridlines(draw_labels = True) # transform = ccrs.Geodetic(),

        # Choose map as global 
        # ax.set_global()

        # Computes lon and lat grids into projected cartesian coordinate X and Y
        # as well as maskedZ in order to be in PlateCarree standard projection
        X, Y, maskedZ = z_masked_overlap(ax[i], lon, lat, vars2plot[i],       
                                         source_projection = ccrs.Geodetic())                                   

        # Pcolormesh
        map1 = ax[i].pcolormesh(lon, lat, maskedZ[0:-1,0:-1],                   
                             transform = inProjData, cmap = cmap, norm = norm, shading='nearest')
        # Contour
        if type(contourLines) is np.ndarray:
            cont1 = ax[i].contour(X, Y, maskedZ, contourLines,                         
                                transform = outProjData,
                                colors='k', linewidths=0.7)
            # Add labels over contour lines
            ax[i].clabel(cont1, fmt=' {:.1f} '.format, fontsize='x-large')  
    
        # Compute coastline
        #dispGridCoastline(lon, lat, inProjData, land_mask[i][0,:,:], 1.25) 

        # Titles/Labels and size
        plt.title(figTitle[i], fontsize=title_font_size)

    # Display colorbar
    cbar_ax = fig.add_axes([0.16, 0.01, 0.7, 0.04]) # list [x0, y0, width, height]
    if cbar_title is not None:
        if len(np.unique(np.round(np.diff(norm.boundaries),2))) == 1:
            ticks = np.linspace(norm.vmin, norm.vmax, 7, endpoint=True) # control nbb of ticks for colorbar
            cbar = plt.colorbar(map1, cax=cbar_ax, ticks=ticks, orientation='horizontal', extend='both')
            #cbar = fig.colorbar(map1, cax=cbar_ax, ticks=ticks, orientation='horizontal', extend='both')

        elif len(np.unique(np.round(np.diff(norm.boundaries),2))) != 1:
            plt.text(0.5,-0.15,'Warning: Adaptative colormap (non-linear) !',horizontalalignment='center',
                    verticalalignment='center', transform = ax.transAxes, fontsize=16, color='r', weight='bold') 
            cbar = plt.colorbar(map1, orientation='horizontal', extend='both')
            
        cbar.ax.set_title(cbar_title,size=cbar_label_size)
        cbar.ax.tick_params(labelsize=cbar_tick_size)   


def plotZonalAve(lat, depth, var2plot, contourLines, cmap, norm, figTitle, figXsize, figYsize, cbar_title, cbar_label_size, cbar_tick_size, title_font_size, xy_label_font_size, xy_ticks_font_size):
   
    import numpy as np
    import matplotlib.pyplot as plt 
    
    fig, ax = plt.subplots(figsize=(figXsize, figYsize))
    
    map1 = ax.pcolormesh(lat, depth, var2plot,                   
                          cmap = cmap, norm = norm, shading='auto')
    
    # Contour
    if type(contourLines) is np.ndarray:
        cont1 = ax.contour(lat, depth, var2plot, contourLines,
                            colors='k', linewidths=0.7) 
        # Add labels over contour lines
        ax.clabel(cont1, fmt=' {:.1f} '.format, fontsize='x-large')

    # Display colorbar
    if cbar_title is not None:
        if len(np.unique(np.round(np.diff(norm.boundaries),2))) == 1:
            ticks = np.linspace(norm.vmin, norm.vmax, 7, endpoint=True) # control nbb of ticks for colorbar
            cbar = plt.colorbar(map1, ticks=ticks, orientation='horizontal', extend='both')
            
        elif len(np.unique(np.round(np.diff(norm.boundaries),2))) != 1:
            plt.text(0.5,-0.15,'Warning: Adaptative colormap (non-linear) !',horizontalalignment='center',
                     verticalalignment='center', transform = ax.transAxes, fontsize=16, color='r', weight='bold') 
            cbar = plt.colorbar(map1, orientation='horizontal', extend='both')
            
        cbar.ax.set_title(cbar_title,size=cbar_label_size)
        cbar.ax.tick_params(labelsize=cbar_tick_size)   

    # Titles/Labels and size
    plt.title (figTitle , fontsize = title_font_size)
    plt.xlabel('Latitude (°N)'  , fontsize = xy_label_font_size)
    plt.ylabel('Depth (m)'      , fontsize = xy_label_font_size)
    plt.xticks(fontsize = xy_ticks_font_size)
    plt.yticks(fontsize = xy_ticks_font_size)


def plotZonalAveSubPlots(nrows, ncols, lat, depth, var2plot, contourLines, 
                        cmap, norm, figTitle, figXsize, figYsize, 
                        cbar_title, cbar_label_size, cbar_tick_size, title_font_size):
   
    import numpy as np
    import matplotlib.pyplot as plt 
    
    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(figXsize, figYsize))
    
    for i in np.arange(0,nrows*ncols):

        map1 = ax[i].pcolormesh(lat, depth, var2plot[i],                   
                            cmap = cmap, norm = norm, shading='auto')
        
        # Contour
        if type(contourLines) is np.ndarray:
            cont1 = ax[i].contour(lat, depth, var2plot[i], contourLines,
                                colors='k', linewidths=0.7) 
            if cbar_title == 'Zonal Oxygen Concentration (O2) (µmol.m$^{-3}$)':
                cont11      = ax[i].contour(lat, depth, var2plot[i], np.array([6.5, 65]), colors='w', linewidths=0.7) # Add dysoxia and anoxia contour lines
                ax[i].clabel(cont11, fmt=' {:.1f} '.format, fontsize='x-large')

            # Add labels over contour lines
            ax[i].clabel(cont1, fmt=' {:.1f} '.format, fontsize='x-large')

        # Titles
        ax[i].set_title(figTitle[i], fontsize=title_font_size)

    cbar_ax = fig.add_axes([0.16, 0.01, 0.7, 0.04]) # list [x0, y0, width, height]
    
    # Display colorbar
    if cbar_title is not None:
        if len(np.unique(np.round(np.diff(norm.boundaries),2))) == 1:
            ticks = np.linspace(norm.vmin, norm.vmax, 7, endpoint=True) # control nbb of ticks for colorbar
            cbar = plt.colorbar(map1, cax=cbar_ax, ticks=ticks, orientation='horizontal', extend='both')
            
        elif len(np.unique(np.round(np.diff(norm.boundaries),2))) != 1:
            plt.text(0.5,-0.15,'Warning: Adaptative colormap (non-linear) !',horizontalalignment='center',
                     verticalalignment='center', transform = ax[i].transAxes, fontsize=16, color='r', weight='bold') 
            cbar = plt.colorbar(map1, cax=cbar_ax, orientation='horizontal', extend='both')
            
        cbar.ax.set_title(cbar_title,size=cbar_label_size)
        cbar.ax.tick_params(labelsize=cbar_tick_size)


# def plotZonalAveSubPlots(lat, depth, var2plot, contourLines,  cmap, norm):
#    
#     import matplotlib.pyplot as plt 
#         
#     map1 = plt.pcolormesh(lat, depth, var2plot,                   
#                           cmap = cmap, norm = norm, shading='auto')
#     
#     # Contour
#     cont1 = plt.contour(lat, depth, var2plot, contourLines,
#                         colors='k', linewidths=0.7)     
#     
#     return map1, cont1


def plotMappingZoom(lon, lat, var2plot, contourLines, cbar_title, lonLim, latLim, land_mask, inProjData, outProjData, cmap, norm, figTitle, figXsize, figYsize, cbar_label_size, cbar_tick_size, title_font_size):

    import numpy as np
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    from Toolz import (dispGridCoastline, z_masked_overlap)
    
    fig, ax = plt.subplots(figsize=(figXsize, figYsize), subplot_kw={'projection': outProjData})
    
    # define limits of plot
    ax.set_extent([lonLim[0], lonLim[1], latLim[0], latLim[1]], inProjData)

    # Display meridians and parallels
    ax.gridlines(draw_labels = True) # transform = ccrs.Geodetic(),

    # Computes lon and lat grids into projected cartesian coordinate X and Y
    # as well as maskedZ in order to be in PlateCarree standard projection
    X, Y, maskedZ = z_masked_overlap(ax, lon, lat, var2plot,       
                                      source_projection = ccrs.Geodetic())                                   

    # Pcolormesh
    map1 = ax.pcolormesh(lon, lat, maskedZ[0:-1,0:-1],                   
                          transform = inProjData, cmap = cmap, norm = norm)
    
    # Contour
    if type(contourLines) is np.ndarray:
        cont1 = ax.contour(X, Y, maskedZ, contourLines,                         
                            transform = outProjData,
                            colors='k', linewidths=0.7) 
        # Add labels over contour lines
        ax.clabel(cont1,fmt=' {:.0f} '.format,fontsize='large')  
    
    # Compute coastline
    dispGridCoastline(lon, lat, inProjData, land_mask[0,:,:], 1.25)
        
    # Display colorbar
    if cbar_title is not None:
        if len(np.unique(np.round(np.diff(norm.boundaries),2))) == 1:
            ticks = np.linspace(norm.vmin, norm.vmax, 7, endpoint=True) # control nbb of ticks for colorbar
            cbar = plt.colorbar(map1, ticks=ticks, orientation='horizontal', extend='both')
            
        elif len(np.unique(np.round(np.diff(norm.boundaries),2))) != 1:
            plt.text(0.5,-0.15,'Warning: Adaptative colormap (non-linear) !',horizontalalignment='center',
                     verticalalignment='center', transform = ax.transAxes, fontsize=16, color='r', weight='bold') 
            cbar = plt.colorbar(map1, orientation='horizontal', extend='both')
            
        cbar.ax.set_title(cbar_title,size=cbar_label_size)
        cbar.ax.tick_params(labelsize=cbar_tick_size)
    
    # Titles/Labels and size
    plt.title(figTitle,fontsize=title_font_size)
    
    
def plotMappingLev(lon, lat, var2plot, cbar_title, depth, land_mask, inProjData, outProjData, cmap, norm, cbar_label_size, cbar_tick_size):
    
    import numpy as np
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    from Toolz import z_masked_overlap

    # Multiplot figs
    #----------------------------------#
    # number of plots and figure size
    nb_plots     = int(np.floor(land_mask.shape[0]/2))
    plotsPerLine = 4 # must be an integer
    nb_lines     = int(nb_plots/plotsPerLine)
    
    if nb_lines == 0:
        nb_lines += 1
    if nb_lines*plotsPerLine < nb_plots:
        nb_lines += 1
    if nb_plots >= plotsPerLine:
        nb_cols = plotsPerLine
    else:
        nb_cols = nb_plots
        
    nbOfUnusedSubplots = nb_lines*nb_cols - nb_plots
    #----------------------------------#
    
    fig, axes =  plt.subplots(figsize=(nb_cols*3.2, nb_lines*2.075),nrows=nb_lines, ncols=nb_cols)

    for k in np.arange(nb_plots):
      
        # subplot
        ax = plt.subplot(nb_lines, nb_cols,k+1, projection = outProjData)
        ax = fig.gca()
        ax.set_aspect('auto')
        
        # Display meridians et parallels
        ax.gridlines()
        
        # Choose map as global 
        ax.set_global()
        
        # Computes lon and lat grids into projected cartesian coordinate 
        # X and Y as well as Z  in order to plot data in the right 
        # projection
        X, Y, maskedZ = z_masked_overlap(ax, lon, lat, var2plot[k*2,:,:],       
                                 source_projection = ccrs.Geodetic()) 
    
        # Pcolormesh
        map1 = ax.pcolormesh(lon, lat, maskedZ[0:-1,0:-1],
                             transform = inProjData, cmap = cmap, norm = norm)
    
    
        # Compute coastline
        # dispGridCoastline(lon, lat, data_crs,                   
        #                  land_mask, 1.25)
        
        plt.title(f'{"{:.0f}".format(depth[k*2])} m')
    
    # Display colorbar
    if cbar_title is not None:
        cbar_ax = fig.add_axes([0.16, 0.01, 0.7, 0.04]) # list [x0, y0, width, height]

        if len(np.unique(np.round(np.diff(norm.boundaries),2))) == 1:
            ticks = np.linspace(norm.vmin, norm.vmax, 7, endpoint=True) # control nbb of ticks for colorbar
            cbar = fig.colorbar(map1, cax=cbar_ax, ticks=ticks, orientation='horizontal',extend='both')
            
        elif len(np.unique(np.round(np.diff(norm.boundaries),2))) != 1:
            plt.text(0.5,-0.15,'Warning: Adaptative colormap (non-linear) !',horizontalalignment='center',
                     verticalalignment='center', transform = ax.transAxes, fontsize=16, color='r', weight='bold') 
            cbar = fig.colorbar(map1, cax=cbar_ax, orientation='horizontal',extend='both')
            
        cbar.ax.set_title(cbar_title,size=cbar_label_size)
        cbar.ax.tick_params(labelsize=cbar_tick_size)
    
    # deleting unused subplots
    for i in range(-nbOfUnusedSubplots,0):
        axes[-1,i].set_axis_off()