#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
-------------------------------------------------------------------------------
Created on Mon Oct 11 14:34:22 2021

@author: anthog

This script plots several diagnostics from CMA5A2 outputs.
Depending on the options enabled:

    - Figures can be saved (format can be set as well)
    - PDF file gathering the generated figures can be created
    - Animated gif can be created for some variables (added 17/12/2021)
    - Manual value limits for plots colormap can be activated instead of the 
      automatic scaling for convenient comparison between several simulations

To run the script, jump to the 'PARAMETERS TO DEFINE' section right below after 
the 'LIBRARIES TO LOAD' section and change the variables according to your needs.

Then launch the script, you'll get through several input questions related to
the FIG/PDF/GIF settings. After this step, the script will run following your 
choices.

Note: In order to run, this script will need the 'util' and 'Toolz' 
      folders containing all the necessary functions and templates
-------------------------------------------------------------------------------
"""

#%%===========================================================================#
#                         --< LIBRARIES TO LOAD >--                           #
#=============================================================================#
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from netCDF4 import Dataset
import numpy as np
import cartopy.crs as ccrs
import warnings
from matplotlib.colors import BoundaryNorm
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import math
import xarray as xr
from matplotlib.axes import Axes
from cartopy.mpl.geoaxes import GeoAxes
GeoAxes._pcolormesh_patched = Axes.pcolormesh
import matplotlib.lines as mlines
import matplotlib.patches as mpatches

# Loading custom functions 
from Toolz import (z_masked_overlap, readncfile, dopdf, dispGridCoastline,
                   plotMapping, plotZonalAve, plotZonalAveSubPlots, plotMappingZoom, plotMappingLev,
                   adaptativeColorMap, makeGif)

#%%===========================================================================#
#                         --< PARAMETERS TO DEFINE >--                        #
#=============================================================================#
# files you want to be loaded and their path location
#---------------------------------------------------#
filePrefix   = "C30MaTotV1-3X_SE_4805_4854_1M" # [filePrefix]_grid_[X].nc with [X] = T, U, V or W
maskFile     = "C30MaTMP_mesh_mask.nc"
bathyFile    = "bathyORCA2.RupelianTotalV1.nc"
subBasinFile = "subbasins_rupelianTot.nc"
dataDirPath  = '/Users/anthony/Documents/Model/Data_CMA5A2/' 
#---------------------------------------------------#

# Manual or automatic colormap limits 
# 'y' for manual, 'n' for automatic
#---------------------------------------------------#
manual_lim_bathy   = 'n'; min_lim_bathy   = 0; max_lim_bathy   = 100; step_bathy   = 1
manual_lim_sss     = 'n'; min_lim_sss     = 0; max_lim_sss     = 100; step_sss     = 1
manual_lim_zosalin = 'n'; min_lim_zosalin = 0; max_lim_zosalin = 100; step_zosalin = 1
manual_lim_sst     = 'n'; min_lim_sst     = 0; max_lim_sst     = 100; step_sst     = 1
manual_lim_zotemp  = 'n'; min_lim_zotemp  = 0; max_lim_zotemp  = 100; step_zotemp  = 1
manual_lim_zostrf  = 'n'; min_lim_zostrf  = 0; max_lim_zostrf  = 100; step_zostrf  = 1
manual_lim_bstrf   = 'n'; min_lim_bstrf   = 0; max_lim_bstrf   = 100; step_bstrf   = 1
manual_lim_omlnh   = 'n'; min_lim_omlnh   = 0; max_lim_omlnh   = 100; step_omlnh   = 1
manual_lim_omlsh   = 'n'; min_lim_omlsh   = 0; max_lim_omlsh   = 100; step_omlsh   = 1
manual_lim_intpp   = 'n'; min_lim_intpp   = 0; max_lim_intpp   = 100; step_intpp   = 1
manual_lim_epc100  = 'n'; min_lim_epc100  = 0; max_lim_epc100  = 100; step_epc100  = 1
manual_lim_po4     = 'n'; min_lim_po4     = 0; max_lim_po4     = 100; step_po4     = 1
manual_lim_zopo4   = 'n'; min_lim_zopo4   = 0; max_lim_zopo4   = 100; step_zopo4   = 1
manual_lim_no3     = 'n'; min_lim_no3     = 0; max_lim_no3     = 100; step_no3     = 1
manual_lim_zono3   = 'n'; min_lim_zono3   = 0; max_lim_zono3   = 100; step_zono3   = 1
manual_lim_o2      = 'n'; min_lim_o2      = 0; max_lim_o2      = 100; step_o2      = 1
manual_lim_zoo2    = 'n'; min_lim_zoo2    = 0; max_lim_zoo2    = 100; step_zoo2    = 1
#---------------------------------------------------#

#%%===========================================================================#
#                             --< INPUT prompt >--                            #
#=============================================================================#
# Section that will prompt to user to choose to create figures/pdf/gif

# Figures
#---------------------------------------------------#
promptFig  = f"Create figures? (y/n)\n"
create_fig = input(promptFig)
while create_fig != 'y' and create_fig != 'n':
    print('/!\\ Answer not valid /!\\')
    create_fig = input(promptFig)

if create_fig == 'y':
    savedFigPath = os.getcwd() + '/FIG/' # get path name
    if os.system('[ ! -d "' + savedFigPath + '" ]') == 0: # folder doesn't exist, create one
        os.system('mkdir ' + savedFigPath)

    promptFig2 = f"Figures format? (png/pdf)\n"
    fig_format = input(promptFig2)
    while fig_format != 'png' and create_fig != 'pdf':
        print('/!\\ Answer not valid /!\\')
        fig_format = input(promptFig2)
#---------------------------------------------------#

# PDF
#---------------------------------------------------#
promptPdf  = f"Create pdf? (y/n)\n"
create_pdf = input(promptPdf)
while create_pdf != 'y' and create_pdf != 'n':
    print('/!\\ Answer not valid /!\\')
    create_pdf = input(promptPdf)

if create_pdf == 'y':
    savedPdfPath = os.getcwd() + '/PDF/' # get path name
    if os.system('[ ! -d "' + savedPdfPath + '" ]') == 0: # folder doesn't exist, create one
        os.system('mkdir ' + savedPdfPath)

    promptPdf2 = f"Notes to write in pdf?\n"
    sentence_to_use = input(promptPdf2)
#---------------------------------------------------#

# GIF
#---------------------------------------------------#
promptGif  = f"Create gif? (y/n)\n"
create_gif = input(promptGif)
while create_gif != 'y' and create_gif != 'n':
    print('/!\\ Answer not valid /!\\')
    create_gif = input(promptGif)

if create_gif == 'y':
    savedGifPath = os.getcwd() + '/GIF/' # get path name
    if os.system('[ ! -d "' + savedGifPath + '" ]') == 0: # folder doesn't exist, create one
        os.system('mkdir ' + savedGifPath)

    promptGif2 = f"Time in ms between each frame of the gif file? (Integer)\n"
    timeInterv = input(promptGif2)
    while timeInterv.isdigit() is False:
        print('/!\\ Answer not valid /!\\')
        timeInterv = input(promptGif2)

    timeInterv = int(timeInterv)
#---------------------------------------------------#

#%%===========================================================================#
#                             --< READING DATA >--                            #
#=============================================================================#

# Variables to load
#----------------------------------#
dict_mask = {
    'e1u':'e1u',
    'e2u':'e2u',
    'e3u_0':'e3u_0',
    'e1t':'e1t',
    'e2t':'e2t',
    'e3t_0':'e3t_0',    
    'tmask':'tmask',
    }

dict_subbasin = {
    'atlmsk':'atlmsk',
    'pacmsk':'pacmsk',
    'indmsk':'indmsk',
    }

dict_bathy = {
    'lon':'NaV_LoN',
    'lat':'nav_lat',
    'bathy':'Bathymetry',
    }

dict_T = {
    'time':'time_counter',
    'deptht':'deptht',
    'wind_speed':'windsp',
    'SSS':'sos',
    'salinity':'so',
    'SST':'tos',
    'temp':'thetao',
    'phys_seaice_frac':'siconc',
    'omlmax':'omlmax',
    }

dict_U = {
    'vozocrtx':'uo',
    }

dict_diaptr = {
    'zo_lat':'nav_lat',
    'zo_temp':'zotemglo',
    'zo_stream_function':'zomsfglo',
    }

dict_diad_T = {
    'TPP':'TPP',
    'EPC100':'EPC100',
    'DPCO2':'Dpco2',
    }

dict_ptrc_T = {
    'PO4':'PO4',
    'NO3':'NO3',
    'O2':'O2',
    'DIC':'DIC',
    }

dict_histmth = {
    'lonWind':'lon',
    'latWind':'lat',
    'u850':'u850',
    'v850':'v850',
    'z850':'z850',
    }

dict_V = {}
dict_W = {}
#----------------------------------#

# File names
#----------------------------------#
nc_file_mask     = dataDirPath + maskFile
nc_file_bathy    = dataDirPath + bathyFile
nc_file_subbasin = dataDirPath + subBasinFile
nc_file_T        = dataDirPath + filePrefix + "_grid_T.nc"
nc_file_U        = dataDirPath + filePrefix + "_grid_U.nc"
nc_file_V        = dataDirPath + filePrefix + "_grid_V.nc"
nc_file_W        = dataDirPath + filePrefix + "_grid_W.nc"
nc_file_diaptr   = dataDirPath + filePrefix + "_diaptr.nc"
nc_file_diad_T   = dataDirPath + filePrefix + "_diad_T.nc"
nc_file_ptrc_T   = dataDirPath + filePrefix + "_ptrc_T.nc"
nc_file_histmth  = dataDirPath + filePrefix + "_histmth.nc"
#----------------------------------#

# Load variables in DATA var (check if file exist before loading)
#----------------------------------#
DATA = {}

if os.path.isfile(nc_file_mask):
    DATA = readncfile(nc_file_mask, dict_mask, DATA)
    
if os.path.isfile(nc_file_bathy):
    DATA = readncfile(nc_file_bathy, dict_bathy, DATA)
    
if os.path.isfile(nc_file_subbasin):
    DATA = readncfile(nc_file_subbasin, dict_subbasin, DATA)
    
if os.path.isfile(nc_file_T):
    DATA = readncfile(nc_file_T, dict_T, DATA)
    
if os.path.isfile(nc_file_U):
    DATA = readncfile(nc_file_U, dict_U, DATA)
    
if os.path.isfile(nc_file_diaptr):
    DATA = readncfile(nc_file_diaptr, dict_diaptr, DATA)
    
if os.path.isfile(nc_file_diad_T):
    DATA = readncfile(nc_file_diad_T, dict_diad_T, DATA)
    
if os.path.isfile(nc_file_ptrc_T):
    DATA = readncfile(nc_file_ptrc_T, dict_ptrc_T, DATA)
    
if os.path.isfile(nc_file_histmth):
    DATA = readncfile(nc_file_histmth, dict_histmth, DATA)
#----------------------------------#

# Read variables in DATA and send them to work space (may be removed ?)
#----------------------------------#
for varname in DATA:
    exec(varname + " = DATA['" + varname + "']")
#----------------------------------#

#%%===========================================================================#
#                             --< CALCULATIONS >--                            #
#=============================================================================#

Mc = 12.0107 # molar mass of carbon

# Compute grids for zonal plots
#----------------------------------#
if "zo_lat" in globals() and "deptht" in globals():

    zo_lat             = np.squeeze(zo_lat)
    latGrid, depthGrid = np.meshgrid(zo_lat, deptht)
#----------------------------------#

# mask
#----------------------------------#    
if "tmask" in globals():

    tmask     = np.squeeze(tmask)
    land_mask = np.ma.masked_where(tmask==0,tmask) # puts masked value instead of 0 in land_mask
#----------------------------------#

# subbasin mask
#----------------------------------#
if "atlmsk" in globals() and "pacmsk" in globals() and "indmsk" in globals():

    atlmsk2 = np.where(atlmsk==1,1,0)
    pacmsk2 = np.where(pacmsk==1,2,0)
    indmsk2 = np.where(indmsk==1,3,0)
    
    subbasin_mask = atlmsk2+pacmsk2+indmsk2
    
    # FIX: Remove overlap between subbasins ATL and IND, when
    #      theses subbasins overlap, the plot isn't working well
    #      The overlap is arbitrarily attributed to IND ocean
    subbasin_mask = np.where(subbasin_mask==4,3,subbasin_mask)
    
    subbasin_mask = np.ma.masked_where(subbasin_mask==0,subbasin_mask)
    
    atlmsk = np.ma.masked_where(atlmsk==0,atlmsk)
    pacmsk = np.ma.masked_where(pacmsk==0,pacmsk)
    indmsk = np.ma.masked_where(indmsk==0,indmsk)
#----------------------------------#

# e2u/e3u_0
#----------------------------------#
if "e2u" in globals() and "e3u_0" in globals():

    e2u   = np.squeeze(e2u)   # e2u   : 182 x 149      CELL width (cell size along 2nd dim)
    e3u_0 = np.squeeze(e3u_0) # e3u_0 : 182 x 149 x 31 CELL height (cell size along 3rd dim)
#----------------------------------#

# Bathymetry
#----------------------------------#
if "land_mask" in globals() and "bathy" in globals():

    bathy = bathy*land_mask[0,:,:]
#----------------------------------#

# Wind speed 
#----------------------------------#
if "wind_speed" in globals():

    wind_speed = np.mean(wind_speed,axis=0) # Yearly average
#----------------------------------#

# Wind 850 hPa
#----------------------------------#
if "u850" in globals() and "v850" in globals() and "z850" in globals() and "lonWind" in globals() and "latWind" in globals():
    
    lonWindGrid, latWindGrid = np.meshgrid(lonWind,latWind)
    
    u850 = np.mean(u850,axis=0) 
    v850 = np.mean(v850,axis=0)
    z850 = np.mean(z850,axis=0)
    
    u850 = np.ma.masked_where(u850 > 1e30, u850)
    v850 = np.ma.masked_where(v850 > 1e30, v850)
    z850 = np.ma.masked_where(z850 > 1e30, z850)
#----------------------------------#

# Salinity
#----------------------------------#
if "SSS" in globals():

    SSS_mean = np.mean(SSS,axis=0) # Yearly average

if "salinity" in globals():  

    salinity_mean = np.mean(salinity,axis=0) # Yearly average
#----------------------------------#

# Temperature
#----------------------------------#
if "SST" in globals():

    SST_mean = np.mean(SST,axis=0)  # Yearly average
    
if "temp" in globals():

    temp_mean = np.mean(temp,axis=0) # Yearly average
    
if "zo_temp" in globals():
    
    zo_temp      = np.squeeze(zo_temp)
    zo_temp_mean = np.mean(zo_temp,axis=0) # Yearly average
    zo_temp_mean = np.ma.masked_where(zo_temp_mean==0, zo_temp_mean) # replace 0 by masked values
#----------------------------------#

# Zonal stream function
#----------------------------------#
if "zo_stream_function" in globals():

    zo_stream_function      = np.squeeze(zo_stream_function)
    zo_stream_function_mean = np.mean(zo_stream_function,axis=0)
#----------------------------------#

# Ocean mixed layer thickness
#----------------------------------#
if "omlmax" in globals():

    omlmaxNH = np.mean(omlmax[0:2,:,:],axis=0) # mixed layer thickness for winter in northern hemisphere (jan-mar)
    omlmaxSH = np.mean(omlmax[6:8,:,:],axis=0) # mixed layer thickness for winter in southern hemisphere (jul-sep)
#----------------------------------#

# TPP
#----------------------------------#
if "TPP" in globals():

    TPP        = TPP * Mc * 86400       # mol/m3/s -> g/m3/d    
    INTPP      = np.sum(TPP, axis=1)    # Integration over z
    INTPP_mean = np.mean(INTPP, axis=0) # Integration over z
#----------------------------------#

# EPC100
#----------------------------------#
if "EPC100" in globals():

    EPC100      = EPC100 * Mc * 86400    # mol/m2/s -> g/m2/d
    EPC100_mean = np.mean(EPC100,axis=0) # Yearly average
#----------------------------------#

# PO4
#----------------------------------#
if "PO4" in globals():

    PO4_mean = np.mean(PO4,axis=0)
#----------------------------------#

# NO3
#----------------------------------#
if "NO3" in globals():

    NO3_mean = np.mean(NO3,axis=0)
#----------------------------------#

# O2
#----------------------------------#
if "O2" in globals():

    O2_mean = np.mean(O2,axis=0)
#----------------------------------#

# COMPUTE ZONAL AVE 
#-----------------------------------------------------------------------------#
# Plots to see lat grid
# for i in np.arange(1,148): plt.plot(lat[i,:])
# for i in np.arange(96,148): plt.plot(lat[i,:])

if "zo_lat" in globals():

    varZero = np.ma.zeros(shape=(len(time), deptht.shape[0], lat.shape[0]))
    
    # Initialize variables that will contain zonal mean
    #-------------------------------#
    if "temp" in globals():    
        zo_temp2       = varZero.copy()
        zo_temp_atlmsk = varZero.copy()
        zo_temp_pacmsk = varZero.copy()
        zo_temp_indmsk = varZero.copy()

    if "salinity" in globals():
        zo_salinity        = varZero.copy()
        zo_salinity_atlmsk = varZero.copy()
        zo_salinity_pacmsk = varZero.copy()
        zo_salinity_indmsk = varZero.copy()       
        
    if "PO4" in globals():
        zo_PO4        = varZero.copy()
        zo_PO4_atlmsk = varZero.copy()
        zo_PO4_pacmsk = varZero.copy()
        zo_PO4_indmsk = varZero.copy()

    if "NO3" in globals():
        zo_NO3        = varZero.copy()
        zo_NO3_atlmsk = varZero.copy()
        zo_NO3_pacmsk = varZero.copy()
        zo_NO3_indmsk = varZero.copy()

    if "O2" in globals():
        zo_O2        = varZero.copy()
        zo_O2_atlmsk = varZero.copy()
        zo_O2_pacmsk = varZero.copy()
        zo_O2_indmsk = varZero.copy()
    #-------------------------------#
    
    # Over lat dimension, for indices from 0 to 95, lat grid is regular 
    # (lat[:,i] whatever i is constant) so we use standard mean
    #-------------------------------#
    if "temp" in globals():
        zo_temp2[:,:,0:96]       = np.ma.mean(np.squeeze(temp[:,:,0:96,:]), axis=3)
        zo_temp_atlmsk[:,:,0:96] = np.ma.mean(np.squeeze(temp[:,:,0:96,:] * atlmsk[0:96,:]), axis=3)
        zo_temp_pacmsk[:,:,0:96] = np.ma.mean(np.squeeze(temp[:,:,0:96,:] * pacmsk[0:96,:]), axis=3)
        zo_temp_indmsk[:,:,0:96] = np.ma.mean(np.squeeze(temp[:,:,0:96,:] * indmsk[0:96,:]), axis=3)
    
    if "salinity" in globals():
        zo_salinity[:,:,0:96]        = np.ma.mean(np.squeeze(salinity[:,:,0:96,:]), axis=3)
        zo_salinity_atlmsk[:,:,0:96] = np.ma.mean(np.squeeze(salinity[:,:,0:96,:] * atlmsk[0:96,:]), axis=3)
        zo_salinity_pacmsk[:,:,0:96] = np.ma.mean(np.squeeze(salinity[:,:,0:96,:] * pacmsk[0:96,:]), axis=3)
        zo_salinity_indmsk[:,:,0:96] = np.ma.mean(np.squeeze(salinity[:,:,0:96,:] * indmsk[0:96,:]), axis=3)

    if "PO4" in globals():
        zo_PO4[:,:,0:96]        = np.ma.mean(np.squeeze(PO4[:,:,0:96,:]), axis=3)
        zo_PO4_atlmsk[:,:,0:96] = np.ma.mean(np.squeeze(PO4[:,:,0:96,:] * atlmsk[0:96,:]), axis=3)
        zo_PO4_pacmsk[:,:,0:96] = np.ma.mean(np.squeeze(PO4[:,:,0:96,:] * pacmsk[0:96,:]), axis=3)
        zo_PO4_indmsk[:,:,0:96] = np.ma.mean(np.squeeze(PO4[:,:,0:96,:] * indmsk[0:96,:]), axis=3)
    
    if "NO3" in globals():
        zo_NO3[:,:,0:96]        = np.ma.mean(np.squeeze(NO3[:,:,0:96,:]), axis=3)
        zo_NO3_atlmsk[:,:,0:96] = np.ma.mean(np.squeeze(NO3[:,:,0:96,:] * atlmsk[0:96,:]), axis=3)
        zo_NO3_pacmsk[:,:,0:96] = np.ma.mean(np.squeeze(NO3[:,:,0:96,:] * pacmsk[0:96,:]), axis=3)
        zo_NO3_indmsk[:,:,0:96] = np.ma.mean(np.squeeze(NO3[:,:,0:96,:] * indmsk[0:96,:]), axis=3)
    
    if "O2" in globals():
        zo_O2[:,:,0:96]        = np.ma.mean(np.squeeze(O2[:,:,0:96,:]), axis=3)
        zo_O2_atlmsk[:,:,0:96] = np.ma.mean(np.squeeze(O2[:,:,0:96,:]  * atlmsk[0:96,:]), axis=3)
        zo_O2_pacmsk[:,:,0:96] = np.ma.mean(np.squeeze(O2[:,:,0:96,:]  * pacmsk[0:96,:]), axis=3)
        zo_O2_indmsk[:,:,0:96] = np.ma.mean(np.squeeze(O2[:,:,0:96,:]  * indmsk[0:96,:]), axis=3)
    #-------------------------------#
    
    # From 96 to 148 lat grid is not regular anymore, it oscillates so we 
    # need to define a more accurate way of meaning.
    # For this, we take variable zo_lat as a guide to compute this mean
    #-------------------------------#
    for j in np.arange(96,149): 
          
        # limits
        if j < 148:
            limit_inf = zo_lat[j] - (zo_lat[j]   - zo_lat[j-1]) /2
            limit_sup = zo_lat[j] + (zo_lat[j+1] - zo_lat[j])   /2
        else:
            limit_inf = zo_lat[j] - (zo_lat[j]   - zo_lat[j-1]) /2
            limit_sup = zo_lat[j] + 1
        
        ind = np.where(np.logical_and(lat >= limit_inf, lat <= limit_sup))
        
        if "temp" in globals():
            zo_temp2[:,:,j]        = np.ma.mean(temp[:,:,ind[0],ind[1]], axis=2)
            zo_temp_atlmsk[:,:,j]  = np.ma.mean(temp[:,:,ind[0],ind[1]] * atlmsk[ind[0],ind[1]], axis=2)
            zo_temp_pacmsk[:,:,j]  = np.ma.mean(temp[:,:,ind[0],ind[1]] * pacmsk[ind[0],ind[1]], axis=2)
            zo_temp_indmsk[:,:,j]  = np.ma.mean(temp[:,:,ind[0],ind[1]] * indmsk[ind[0],ind[1]], axis=2)
        
        if "salinity" in globals():
            zo_salinity[:,:,j]        = np.ma.mean(salinity[:,:,ind[0],ind[1]], axis=2)
            zo_salinity_atlmsk[:,:,j] = np.ma.mean(salinity[:,:,ind[0],ind[1]] * atlmsk[ind[0],ind[1]], axis=2)
            zo_salinity_pacmsk[:,:,j] = np.ma.mean(salinity[:,:,ind[0],ind[1]] * pacmsk[ind[0],ind[1]], axis=2)
            zo_salinity_indmsk[:,:,j] = np.ma.mean(salinity[:,:,ind[0],ind[1]] * indmsk[ind[0],ind[1]], axis=2)

        if "PO4" in globals():
            zo_PO4[:,:,j]        = np.ma.mean(PO4[:,:,ind[0],ind[1]], axis=2)
            zo_PO4_atlmsk[:,:,j] = np.ma.mean(PO4[:,:,ind[0],ind[1]] * atlmsk[ind[0],ind[1]], axis=2)
            zo_PO4_pacmsk[:,:,j] = np.ma.mean(PO4[:,:,ind[0],ind[1]] * pacmsk[ind[0],ind[1]], axis=2)
            zo_PO4_indmsk[:,:,j] = np.ma.mean(PO4[:,:,ind[0],ind[1]] * indmsk[ind[0],ind[1]], axis=2)
        
        if "NO3" in globals():
            zo_NO3[:,:,j]        = np.ma.mean(NO3[:,:,ind[0],ind[1]], axis=2)
            zo_NO3_atlmsk[:,:,j] = np.ma.mean(NO3[:,:,ind[0],ind[1]] * atlmsk[ind[0],ind[1]], axis=2)
            zo_NO3_pacmsk[:,:,j] = np.ma.mean(NO3[:,:,ind[0],ind[1]] * pacmsk[ind[0],ind[1]], axis=2)
            zo_NO3_indmsk[:,:,j] = np.ma.mean(NO3[:,:,ind[0],ind[1]] * indmsk[ind[0],ind[1]], axis=2)
        
        if "O2" in globals():
            zo_O2[:,:,j]        = np.ma.mean(O2[:,:,ind[0],ind[1]], axis=2)
            zo_O2_atlmsk[:,:,j] = np.ma.mean(O2[:,:,ind[0],ind[1]] * atlmsk[ind[0],ind[1]], axis=2)
            zo_O2_pacmsk[:,:,j] = np.ma.mean(O2[:,:,ind[0],ind[1]] * pacmsk[ind[0],ind[1]], axis=2)
            zo_O2_indmsk[:,:,j] = np.ma.mean(O2[:,:,ind[0],ind[1]] * indmsk[ind[0],ind[1]], axis=2)
    #-------------------------------#

    # Average over time
    #-------------------------------#
    if "temp" in globals():
        zo_temp2_mean       = np.mean(zo_temp2, axis=0)
        zo_temp_atlmsk_mean = np.mean(zo_temp_atlmsk, axis=0)
        zo_temp_pacmsk_mean = np.mean(zo_temp_pacmsk, axis=0)
        zo_temp_indmsk_mean = np.mean(zo_temp_indmsk, axis=0)

    if "salinity" in globals():
        zo_salinity_mean        = np.mean(zo_salinity, axis=0)
        zo_salinity_atlmsk_mean = np.mean(zo_salinity_atlmsk, axis=0)
        zo_salinity_pacmsk_mean = np.mean(zo_salinity_pacmsk, axis=0)
        zo_salinity_indmsk_mean = np.mean(zo_salinity_indmsk, axis=0)

    if "PO4" in globals():
        zo_PO4_mean        = np.mean(zo_PO4, axis=0)
        zo_PO4_atlmsk_mean = np.mean(zo_PO4_atlmsk, axis=0)
        zo_PO4_pacmsk_mean = np.mean(zo_PO4_pacmsk, axis=0)
        zo_PO4_indmsk_mean = np.mean(zo_PO4_indmsk, axis=0)

    if "NO3" in globals():
        zo_NO3_mean        = np.mean(zo_NO3, axis=0)
        zo_NO3_atlmsk_mean = np.mean(zo_NO3_atlmsk, axis=0)
        zo_NO3_pacmsk_mean = np.mean(zo_NO3_pacmsk, axis=0)
        zo_NO3_indmsk_mean = np.mean(zo_NO3_indmsk, axis=0)
    
    if "O2" in globals():
        zo_O2_mean        = np.mean(zo_O2,axis=0)
        zo_O2_atlmsk_mean = np.mean(zo_O2_atlmsk,axis=0)
        zo_O2_pacmsk_mean = np.mean(zo_O2_pacmsk,axis=0)
        zo_O2_indmsk_mean = np.mean(zo_O2_indmsk,axis=0)
    #-------------------------------#
#-----------------------------------------------------------------------------#  

# COMPUTE BAROTROPIC STREAM FUNCTION (bsf)
#-----------------------------------------------------------------------------#
if "vozocrtx" in globals() and "e3u_0" in globals() and "e2u" in globals():

    # puts 0 instead of masked value in vosocrtx
    vozocrtx = np.where(vozocrtx.mask==True,0,vozocrtx)   
                            
    # Zonal volume flux"/UNITS="m3/s" 
    f_m_x = vozocrtx * e3u_0 * e2u
    
    # Zonal volume flux, vertically integrated"/UNITS="m3/s"                            
    f_m_x_z = np.sum(f_m_x,axis = 1)
    
    # "Barotropic stream function"/UNITS="Sv"                                       
    bsf0 = np.empty([f_m_x_z.shape[0],f_m_x_z.shape[1],f_m_x_z.shape[2]])     
    bsf0 = np.cumsum(f_m_x_z,axis=1)*1E-6
    
    # Arbitrarily fixes a point where BSF=0. Here : eurasia, for ORCA2
    bsf_ref = np.sum(f_m_x_z[:,0:139,19],axis=1)*1E-6
                                
    # Expand vector bsf_ref (12) into a matrix bsf_ref_mat (12x149x182) so bsf 
    # calcualtion works (each element of bsf_ref is repeated over a matrix 149x182)
    bsf_ref_mat = np.expand_dims(np.expand_dims(bsf_ref, axis=1),axis=2) # expand dimensions so bsf_ref_mat is 12x1x1
    bsf_ref_mat = np.repeat(np.repeat(bsf_ref_mat,lon.shape[0],axis=1),lon.shape[1],axis=2)
    
    # "Barotropic stream function"/UNITS="Sv"
    bsf = (bsf0 - bsf_ref_mat) * land_mask[0,:,:]
    
    # Name barotropic stream function with more explicit variable name and average over time
    baro_stream_function      = bsf
    baro_stream_function_mean = np.mean(baro_stream_function,axis=0)
#-----------------------------------------------------------------------------#

#%%===========================================================================#
#                           --< PLOT PARAMETERS >--                           #
#=============================================================================#

# Common parmeters
#----------------------------------#
# Map projection for plots
projDataIn  = ccrs.PlateCarree() # Projection data are originally in
projDataOut = ccrs.Robinson(central_longitude = 0) # Final projection in plots

# figure size
figXsize = 10
figYsize = 10

# Various plot's label sizes
cbar_label_size    = 20
cbar_tick_size     = 15
title_font_size    = 20
xy_label_font_size = 20
xy_ticks_font_size = 20
plot_font_size     = 15

# Colormaps and norms
cmapColor_standard = 'viridis'
normNone           = None

# Figure title
figTitleYear  = 'Yearly average'
figTitleMonth = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December']
figTitleNone = None

# Colorbar title
cbarTitleNone = None

# Contour lines
nbContourLines = 8 # nb of contour lines plotted
#----------------------------------#

# SUBBASINS
#----------------------------------#
from palettable.cartocolors.qualitative import Safe_3 # colormap import
cmapColor_subbasin = Safe_3.mpl_colormap
# cmapColor_subbasin = mpl.cm.get_cmap('viridis')
figTitleSubbasin    = ''
contour_subbasin    = None
#----------------------------------#

# BATHYMETRY 
#----------------------------------#
cbar_title_bathy = 'Bathymetry (m)'
cmapColor_bathy  = 'YlGnBu'

#  Automatic or manual colormap limits bathymetry
if manual_lim_bathy == 'n': 
    cmap_bathy = cmapColor_bathy
    norm_bathy = None
elif manual_lim_bathy == 'y':
    cmap_bathy = mpl.cm.get_cmap(cmapColor_bathy)
    bounds     = np.arange(min_lim_bathy ,max_lim_bathy ,step_bathy)
    norm_bathy = mpl.colors.BoundaryNorm(bounds, cmap_bathy.N)

contour_bathy = None
#----------------------------------#

# WIND SPEED 850 hPa
#----------------------------------#
if "z850" in globals():
    
    figLegWind850      = 'Geopotential height 850 hPa (z850) (m)'
    cbar_title_wind850 = 'Wind 850hPa (u850, v850) (m.s$^{-1}$)'
    
    # contour_z850       = np.array([1200, 1300, 1400, 1500, 1520])
    contour_z850 = np.linspace(z850.min(),z850.max(),nbContourLines)
#----------------------------------#

# SALINITY
#----------------------------------#
cmapColor_salinity = 'magma' #gnuplot2

if "SSS" in globals(): 

    cbar_title_sss     = 'SSS (sos) (PSU)'
    
    #  Automatic or manual colormap limits SSS
    if manual_lim_sss  == 'n':        
        cmap_sss, norm_sss = adaptativeColorMap(SSS_mean, 0.015, 0.2, cmapColor_salinity)
    elif manual_lim_sss == 'y':
        cmap_sss = mpl.cm.get_cmap(cmapColor_salinity)
        bounds   = np.arange(min_lim_sss ,max_lim_sss ,step_sss)
        norm_sss = mpl.colors.BoundaryNorm(bounds, cmap_sss.N)
    
    # Contour
    if len(np.unique(np.round(np.diff(norm_sss.boundaries),2))) != 1: # norm_sss is not None:
        intervContour = np.int32(len(norm_sss.boundaries)/nbContourLines)
        contour_sss   = norm_sss.boundaries[::intervContour]
    else:
        contour_sss    = np.linspace(SSS_mean.min(),SSS_mean.max(),nbContourLines)

if "zo_salinity" in globals(): 

    cbar_title_zosalin         = 'Zonal salinity (PSU)'
    
    #  Automatic or manual colormap limits zonal salinity
    if manual_lim_zosalin  == 'n':        
        cmap_zosalin, norm_zosalin = adaptativeColorMap(zo_salinity_mean, 0.015, 0.2, cmapColor_salinity)
    elif manual_lim_zosalin == 'y':
        cmap_zosalin = mpl.cm.get_cmap(cmapColor_salinity)
        bounds       = np.arange(min_lim_zosalin ,max_lim_zosalin ,step_zosalin)
        norm_zosalin = mpl.colors.BoundaryNorm(bounds, cmap_zosalin.N)
    # Contour
    if len(np.unique(np.round(np.diff(norm_zosalin.boundaries),2))) != 1: # norm_zosalin is not None:
        intervContour = np.int32(len(norm_zosalin.boundaries)/nbContourLines)
        contour_zosalin = norm_zosalin.boundaries[::intervContour]
    else:
        contour_zosalin = np.linspace(zo_salinity_mean.min(),zo_salinity_mean.max(),nbContourLines)
#----------------------------------#

# TEMPERATURE
#----------------------------------#
cmapColor_temp           = 'inferno' #gnuplot

if "SST" in globals(): 

    cbar_title_sst     = 'SST (tos) (°C)'
    
    #  Automatic or manual colormap limits SST
    if manual_lim_sst  == 'n':        
        cmap_sst, norm_sst = adaptativeColorMap(SST_mean, 0.015, 0.2, cmapColor_temp)
    elif manual_lim_sst == 'y':
        cmap_sst = mpl.cm.get_cmap(cmapColor_temp)
        bounds   = np.arange(min_lim_sst ,max_lim_sst ,step_sst)
        norm_sst = mpl.colors.BoundaryNorm(bounds, cmap_sst.N)
    
    # Contour
    if len(np.unique(np.round(np.diff(norm_sst.boundaries),2))) != 1: # norm_sst is not None:
        intervContour = np.int32(len(norm_sst.boundaries)/nbContourLines)
        contour_sst   = norm_sst.boundaries[::intervContour]
    else:
        contour_sst = np.linspace(SST_mean.min(),SST_mean.max(),nbContourLines)

if "zo_temp" in globals(): 

    cbar_title_zotemp        = 'Zonal temperature (zotemglo) (°C)'
    
    #  Automatic or manual colormap limits zonal temperature
    if manual_lim_zotemp  == 'n':        
        cmap_zotemp, norm_zotemp = adaptativeColorMap(zo_temp_mean, 0.015, 0.2, cmapColor_temp)
    elif manual_lim_zotemp == 'y':
        cmap_zotemp = mpl.cm.get_cmap(cmapColor_temp)
        bounds      = np.arange(min_lim_zotemp ,max_lim_zotemp ,step_zotemp)
        norm_zotemp = mpl.colors.BoundaryNorm(bounds, cmap_zotemp.N)
    
    # Contour
    if len(np.unique(np.round(np.diff(norm_zotemp.boundaries),2))) != 1: # norm_zotemp is not None:
        intervContour = np.int32(len(norm_zotemp.boundaries)/nbContourLines)
        contour_zotemp   = norm_zotemp.boundaries[::intervContour]
    else:
        contour_zotemp = np.linspace(zo_temp_mean.min(),zo_temp_mean.max(),nbContourLines)
        # contour_zotemp = np.linspace(norm_zotemp.boundaries[0],norm_zotemp.boundaries[-1],nbContourLines) # norm_zotemp.boundaries[::nbContourLines]
#----------------------------------#

# STREAM FUNCTION
#----------------------------------#
if "vozocrtx" in globals():
    
    cmapColor_strf    = 'BrBG_r'
    
    cbar_title_zostrf = 'Zonal stream function (zomsfglo) (Sv = 10$^6$m$^3$.s$^{-1}$)'
    contour_zostrf    = np.array([-20, -10, 0, 10, 20])
    # Custom cmap
    bounds            = np.arange(-42,-15,4)
    bounds            = np.append(bounds,np.arange(-15,15,0.5))
    bounds            = np.append(bounds,np.arange(18,43,4))
    cmap_zostrf       = plt.get_cmap(cmapColor_strf,len(bounds)-1)
    norm_zostrf       = mpl.colors.BoundaryNorm(bounds, cmap_zostrf.N)
    
    cbar_title_bstrf  = 'Barotropic stream function (Sv = 10$^6$m$^3$.s$^{-1}$)'
    contour_bstrf     = np.array([-20, 0, 20])
    # Custom cmap
    bounds            = np.arange(-70,71,1)
    cmap_bstrf        = plt.get_cmap(cmapColor_strf,len(bounds)-1)
    norm_bstrf        = mpl.colors.BoundaryNorm(bounds, cmap_bstrf.N)
#----------------------------------#

# OCEAN MIXED LAYER
#----------------------------------#
if "omlmax" in globals():

    cmapColor_oml    = 'PuBu'
    
    fig_title_omlnh  = 'January to March average'
    cbar_title_omlnh = 'Ocean mixed layer thickness (omlmax) (m) - Northern hemisphere'
    
    #  Automatic or manual colormap limits omlNH
    if manual_lim_omlnh  == 'n':        
        cmap_omlnh, norm_omlnh = adaptativeColorMap(omlmaxNH, 0.015, 0.2, cmapColor_oml)
    elif manual_lim_omlnh == 'y':
        cmap_omlnh = mpl.cm.get_cmap(cmapColor_oml)
        bounds     = np.arange(min_lim_omlnh ,max_lim_omlnh ,step_omlnh)
        norm_omlnh = mpl.colors.BoundaryNorm(bounds, cmap_omlnh.N)
    
    if len(np.unique(np.round(np.diff(norm_omlnh.boundaries),2))) != 1: # norm_omlnh is not None:
        intervContour = np.int32(len(norm_omlnh.boundaries)/nbContourLines)
        contour_omlnh   = norm_omlnh.boundaries[::intervContour]
    else:
        contour_omlnh = np.linspace(omlmaxNH.min(),omlmaxNH.max(),nbContourLines)
        
    fig_title_omlsh  = 'July to September average'
    cbar_title_omlsh = 'Ocean mixed layer thickness (omlmax) (m) - Southern hemisphere'
   
    # Automatic or manual colormap limits omlSH
    if manual_lim_omlsh  == 'n':        
        cmap_omlsh, norm_omlsh = adaptativeColorMap(omlmaxSH, 0.015, 0.2, cmapColor_oml)
    elif manual_lim_omlsh == 'y':
        cmap_omlsh = mpl.cm.get_cmap(cmapColor_oml)
        bounds     = np.arange(min_lim_omlsh ,max_lim_omlsh ,step_omlsh)
        norm_omlsh = mpl.colors.BoundaryNorm(bounds, cmap_omlsh.N)
    
    # Contour
    if len(np.unique(np.round(np.diff(norm_omlsh.boundaries),2))) != 1: # norm_omlsh is not None:
        intervContour = np.int32(len(norm_omlsh.boundaries)/nbContourLines)
        contour_omlsh   = norm_omlsh.boundaries[::intervContour]
    else:
        contour_omlsh = np.linspace(omlmaxSH.min(),omlmaxSH.max(),nbContourLines)
#----------------------------------#

# EPC100 / TPP
#----------------------------------#
cmapColor_epctpp  = 'ocean_r'

if "TPP" in globals():

    cbar_title_intpp = 'Total Primary production of phyto depth integrated (INTPP) (g.m$^{-3}$.d$^{-1}$)'

    # Automatic or manual colormap limits INTPP
    if manual_lim_intpp == 'n':        
        # cmap_intpp = cmapColor_epctpp
        # norm_intpp = None
        cmap_intpp, norm_intpp = adaptativeColorMap(INTPP, 0.015, 0.2, cmapColor_epctpp)
    elif manual_lim_epc100 == 'y':
        cmap_intpp = mpl.cm.get_cmap(cmapColor_epctpp)
        bounds     = np.arange(min_lim_intpp ,max_lim_intpp ,step_intpp)
        norm_intpp = mpl.colors.BoundaryNorm(bounds, cmap_intpp.N)
    
    # if norm_intpp is not None:
    #     contour_intpp        = norm_intpp.boundaries[::nbContourLines]
    # else:
    #     contour_intpp        = np.linspace(INTPP.min(),INTPP.max(),nbContourLines)

    contour_intpp = None

if "EPC100" in globals():
    
    cbar_title_epc100 = 'Export of carbon particles at 100m (EPC100) (g.m$^{-2}$.d$^{-1}$)' # (mol.m$^{-2}$.s$^{-1}$)

    # Automatic or manual colormap limits EPC100
    if manual_lim_epc100 == 'n':        
        # cmap_epc100 = cmapColor_epctpp
        # norm_epc100 = None
        cmap_epc100, norm_epc100 = adaptativeColorMap(EPC100, 0.015, 0.2, cmapColor_epctpp)
    elif manual_lim_epc100 == 'y':
        cmap_epc100 = mpl.cm.get_cmap(cmapColor_epctpp)
        bounds      = np.arange(min_lim_epc100, max_lim_epc100, step_epc100)
        norm_epc100 = mpl.colors.BoundaryNorm(bounds, cmap_epc100.N)
        
    # if norm_epc100 is not None:
    #     contour_epc100        = norm_epc100.boundaries[::nbContourLines]
    # else:
    #     contour_epc100        = np.linspace(EPC100.min(),EPC100.max(),nbContourLines)
    
    contour_epc100 = None
#----------------------------------#

# PO4 / NO3 / O2
#----------------------------------#
cmapColor_po4no3o2 = 'gist_earth'

if "PO4" in globals():

    cbar_title_po4   = 'Phosphate Concentration (PO4) (µmol.m$^{-3}$)'
    cbar_title_zopo4 = 'Zonal Phosphate Concentration (PO4) (µmol.m$^{-3}$)'
    
    # Automatic or manual colormap limits PO4
    if manual_lim_po4 == 'n':
        cmap_po4, norm_po4 = adaptativeColorMap(PO4_mean, 0.015, 0.2, cmapColor_po4no3o2)
    elif manual_lim_po4 == 'y':
        cmap_po4 = mpl.cm.get_cmap(cmapColor_po4no3o2)
        bounds   = np.arange(min_lim_po4, max_lim_po4, step_po4)
        norm_po4 = mpl.colors.BoundaryNorm(bounds, cmap_po4.N)
        
    # Automatic or manual colormap limits zonal PO4
    if manual_lim_zopo4 == 'n':
        cmap_zopo4, norm_zopo4 = adaptativeColorMap(zo_PO4_mean, 0.015, 0.2, cmapColor_po4no3o2)
    elif manual_lim_zopo4 == 'y':
        cmap_zopo4 = mpl.cm.get_cmap(cmapColor_po4no3o2)
        bounds     = np.arange(min_lim_zopo4, max_lim_zopo4,step_zopo4)
        norm_zopo4 = mpl.colors.BoundaryNorm(bounds, cmap_zopo4.N)

    # Contour
    if len(np.unique(np.round(np.diff(norm_zopo4.boundaries),2))) != 1: # norm_zopo4 is not None:
        intervContour = np.int32(len(norm_zopo4.boundaries)/nbContourLines)
        contour_zopo4 = norm_zopo4.boundaries[::intervContour]
    else:
        contour_zopo4 = np.linspace(zo_PO4_mean.min(),zo_PO4_mean.max(),nbContourLines)

if "NO3" in globals():

    cbar_title_no3   = 'Nitrate Concentration (NO3) (µmol.m$^{-3}$)'
    cbar_title_zono3 = 'Zonal Nitrate Concentration (NO3) (µmol.m$^{-3}$)'
    
    # Automatic or manual colormap limits NO3
    if manual_lim_no3 == 'n':
        cmap_no3, norm_no3 = adaptativeColorMap(NO3_mean, 0.015, 0.2, cmapColor_po4no3o2)
    elif manual_lim_no3 == 'y':
        cmap_no3 = mpl.cm.get_cmap(cmapColor_po4no3o2)
        bounds   = np.arange(min_lim_no3,max_lim_no3,step_no3)
        norm_no3 = mpl.colors.BoundaryNorm(bounds, cmap_no3.N)
        
    # Automatic or manual colormap limits zonal NO3
    if manual_lim_zono3 == 'n':
        cmap_zono3, norm_zono3 = adaptativeColorMap(zo_NO3_mean, 0.015, 0.2, cmapColor_po4no3o2)
    elif manual_lim_zono3 == 'y':
        cmap_zono3 = mpl.cm.get_cmap(cmapColor_po4no3o2)
        bounds     = np.arange(min_lim_zono3,max_lim_zono3,step_zono3)
        norm_zono3 = mpl.colors.BoundaryNorm(bounds, cmap_zono3.N)
        
    # Contour
    if len(np.unique(np.round(np.diff(norm_zono3.boundaries),2))) != 1: # norm_zono3 is not None:
        intervContour = np.int32(len(norm_zono3.boundaries)/nbContourLines)
        contour_zono3 = norm_zono3.boundaries[::intervContour]
    else:
        contour_zono3 = np.linspace(zo_NO3_mean.min(),zo_NO3_mean.max(),nbContourLines)

if "O2" in globals():

    cbar_title_o2   = 'Oxygen Concentration (O2) (µmol.m$^{-3}$)'
    cbar_title_zoo2 = 'Zonal Oxygen Concentration (O2) (µmol.m$^{-3}$)'
    
    # Automatic or manual colormap limits O2
    if manual_lim_o2 == 'n':
        cmap_o2, norm_o2 = adaptativeColorMap(O2_mean, 0.015, 0.2, cmapColor_po4no3o2)
    elif manual_lim_o2 == 'y':
        cmap_o2 = mpl.cm.get_cmap(cmapColor_po4no3o2)
        bounds  = np.arange(min_lim_o2,max_lim_o2,step_o2)
        norm_o2 = mpl.colors.BoundaryNorm(bounds, cmap_o2.N)
        
    # Automatic or manual colormap limits zonal O2
    if manual_lim_zoo2 == 'n':
        cmap_zoo2, norm_zoo2 = adaptativeColorMap(zo_O2_mean, 0.015, 0.2, cmapColor_po4no3o2)
    elif manual_lim_zoo2 == 'y':
        cmap_zoo2 = mpl.cm.get_cmap(cmapColor_po4no3o2)
        bounds    = np.arange(min_lim_zoo2,max_lim_zoo2,step_zoo2)
        norm_zoo2 = mpl.colors.BoundaryNorm(bounds, cmap_zoo2.N)
    
    # Contour
    if len(np.unique(np.round(np.diff(norm_zoo2.boundaries),2))) != 1: # norm_zoo2 is not None:
        intervContour = np.int32(len(norm_zoo2.boundaries)/nbContourLines)
        contour_zoo2  = norm_zoo2.boundaries[::intervContour]
    else:
        contour_zoo2 = np.linspace(zo_O2_mean.min(),zo_O2_mean.max(),nbContourLines)
#----------------------------------#

#%%===========================================================================#
#                              --< PLOTTING >--                               #
#=============================================================================#
# Initialization variables
#----------------------------------#
filecount  = 0
savedfiles = []
#----------------------------------#

# SUBBASINS
#-----------------------------------------------------------------------------#
if "subbasin_mask" in globals():

    plotMapping(lon, lat, subbasin_mask, contour_subbasin, cbarTitleNone, land_mask, 
            projDataIn, projDataOut, cmapColor_subbasin, normNone, figTitleSubbasin,
            figXsize, figYsize, cbar_label_size, cbar_tick_size, title_font_size) # cmapColor_standard
    
    # Legend
    # cmap_subbasin = mpl.cm.get_cmap('viridis')
    
    color1 = cmapColor_subbasin(0)
    color2 = cmapColor_subbasin(0.5)
    color3 = cmapColor_subbasin(1.0)
    
    patch1 = mpatches.Patch(color=color1, label='Atlantic sub basin')
    patch2 = mpatches.Patch(color=color2, label='Pacific sub basin')
    patch3 = mpatches.Patch(color=color3, label='Indian sub basin')
    
    plt.legend(handles=[patch1, patch2, patch3],bbox_to_anchor=(0.65, 0.83), bbox_transform=plt.gcf().transFigure, prop={'size': plot_font_size})

    # Save figure
    filename = filePrefix + '_' + 'subbasin.' + fig_format
    pathFilename = savedFigPath + filename
    if create_fig == 'y':                                               
        plt.savefig(pathFilename, format=fig_format, bbox_inches='tight', facecolor='white')
        filecount += 1; lnfile = 'file' + str(filecount) + '.' + fig_format
        os.system('ln -s ' + pathFilename + ' ' + lnfile)
    savedfiles.append(filename)
#-----------------------------------------------------------------------------#

# BATHYMETRY 
#-----------------------------------------------------------------------------#
if "bathy" in globals():
    
    plotMapping(lon, lat, bathy, contour_bathy, cbar_title_bathy, land_mask, 
            projDataIn, projDataOut, cmap_bathy, norm_bathy, figTitleNone,
            figXsize, figYsize, cbar_label_size, cbar_tick_size, title_font_size)
    
    # Save figure
    filename = filePrefix + '_' + 'bathy.' + fig_format
    pathFilename = savedFigPath + filename
    if create_fig == 'y':                                               
        plt.savefig(pathFilename, format=fig_format, bbox_inches='tight', facecolor='white')
        filecount += 1; lnfile = 'file' + str(filecount) + '.' + fig_format
        os.system('ln -s ' + pathFilename + ' ' + lnfile)
    savedfiles.append(filename)
#-----------------------------------------------------------------------------#


# WIND SPEED 850 hPa
#-----------------------------------------------------------------------------#
if "u850" in globals() and "v850" in globals() and "z850" in globals() and "lonWind" in globals() and "latWind" in globals():

    latWindGrid[0,:] = 89.99 # trick so that z_masked_overlap works correctly and plt.contour can display correctly over the quiver map
    latWindGrid[-1,:] = -89.99
    
    color_array = np.sqrt((u850[::2,::2])**2 + (v850[::2,::2])**2)
    
    fig, ax = plt.subplots(figsize=(figXsize, figYsize), subplot_kw={'projection': projDataOut})
    
    # Display meridians and parallels
    ax.gridlines(draw_labels = True) # transform = ccrs.Geodetic(),
    
    
    X2, Y2, maskedZ = z_masked_overlap(ax, lonWindGrid, latWindGrid, z850, source_projection = ccrs.Geodetic())                                   
    
    # Arrow plot
    map1  = ax.quiver(lonWind[::2], latWind[::2], u850[::2,::2], v850[::2,::2], color_array, cmap = cmapColor_standard, transform = projDataIn)
    
    # Contour plot
    cont1 = ax.contour(X2, Y2, maskedZ, contour_z850, transform = projDataOut, colors='r', linewidths=1.5)
    
    # Contour labels
    ax.clabel(cont1,fmt=' {:.1f} '.format,fontsize='x-large')  
    
    # Compute coastline
    dispGridCoastline(lon, lat, projDataIn, land_mask[0,:,:], 1.25)
    
    # Colorbar
    cbar = plt.colorbar(map1, orientation='horizontal', extend='both')
    cbar.ax.set_title(cbar_title_wind850,size=cbar_label_size)
    cbar.ax.tick_params(labelsize=cbar_tick_size)
    
    # Legend does not support <cartopy.mpl.contour.GeoContourSet object at 0x7fc2739f6f10> instances. So:
    lgd_line = mlines.Line2D([], [], color='red', marker='_', markersize=1, label=figLegWind850)
    plt.legend(handles=[lgd_line],bbox_to_anchor=(0.95, 0.32), bbox_transform=plt.gcf().transFigure, prop={'size': plot_font_size})
    
    # Save figure
    filename = filePrefix + '_' + 'wind_speed_850.' + fig_format
    pathFilename = savedFigPath + filename
    if create_fig == 'y':                                               
        plt.savefig(pathFilename, format=fig_format, bbox_inches='tight', facecolor='white')
        filecount += 1; lnfile = 'file' + str(filecount) + '.' + fig_format
        os.system('ln -s ' + pathFilename + ' ' + lnfile)
    savedfiles.append(filename)
#-----------------------------------------------------------------------------#

# SSS
#-----------------------------------------------------------------------------#
if "SSS" in globals():
    
    plotMapping(lon, lat, SSS_mean, contour_sss, cbar_title_sss, land_mask, 
            projDataIn, projDataOut, cmap_sss, norm_sss, figTitleYear,
            figXsize, figYsize, cbar_label_size, cbar_tick_size, title_font_size)
    
    # Save figure
    filename = filePrefix + '_' + 'SSS.' + fig_format
    pathFilename = savedFigPath + filename
    if create_fig == 'y':                                               
        plt.savefig(pathFilename, format=fig_format, bbox_inches='tight', facecolor='white')
        filecount += 1; lnfile = 'file' + str(filecount) + '.' + fig_format
        os.system('ln -s ' + pathFilename + ' ' + lnfile)
    savedfiles.append(filename)

    # GIF
    #-------------------#
    if create_gif == 'y':
        
        for ind in np.arange(0,SSS.shape[0]):
            
            fig_title = figTitleMonth[ind]

            plotMapping(lon, lat, SSS[ind,:,:], contour_sss, cbar_title_sss, land_mask, 
                    projDataIn, projDataOut, cmap_sss, norm_sss, fig_title,
                    figXsize, figYsize, cbar_label_size, cbar_tick_size, title_font_size)

            # Save figure
            ind = str(ind); ind = ind.zfill(2) # add 0 when ind is 1 digit long (1 => 01)
            pathFilename = savedGifPath + ind + '.' + fig_format
            plt.savefig(pathFilename, bbox_inches='tight', facecolor='white')

            plt.close()

        makeGif(savedGifPath, "SSS.gif", timeInterv)
    #-------------------#
#-----------------------------------------------------------------------------#

# ZONAL AVERAGE SALINITY
#-----------------------------------------------------------------------------#
if "zo_salinity" in globals():

    plotZonalAve(latGrid, -1*depthGrid, zo_salinity_mean, contour_zosalin, cbar_title_zosalin, 
                 cmap_zosalin, norm_zosalin, figTitleYear, figXsize, figYsize, cbar_label_size, cbar_tick_size, 
                 title_font_size, xy_label_font_size, xy_ticks_font_size)
    
    filename = filePrefix + '_' + 'zoSalinity.' + fig_format
    pathFilename = savedFigPath + filename
    if create_fig == 'y':                                               
        plt.savefig(pathFilename, format=fig_format, bbox_inches='tight', facecolor='white')
        filecount += 1; lnfile = 'file' + str(filecount) + '.' + fig_format
        os.system('ln -s ' + pathFilename + ' ' + lnfile)
    savedfiles.append(filename)

        # GIF
    #-------------------#
    if create_gif == 'y':
        
        for ind in np.arange(0,zo_salinity.shape[0]):
            
            fig_title = figTitleMonth[ind]

            plotZonalAve(latGrid, -1*depthGrid, zo_salinity[ind,:,:], contour_zosalin, cbar_title_zosalin, 
                         cmap_zosalin, norm_zosalin, fig_title, figXsize, figYsize, cbar_label_size, 
                         cbar_tick_size, title_font_size, xy_label_font_size, xy_ticks_font_size)

            # Save figure
            ind = str(ind); ind = ind.zfill(2) # add 0 when ind is 1 digit long (1 => 01)
            pathFilename = savedGifPath + ind + '.' + fig_format
            plt.savefig(pathFilename, bbox_inches='tight', facecolor='white')

            plt.close()

        makeGif(savedGifPath, "zo_salinity.gif", timeInterv)
    #-------------------#
#-----------------------------------------------------------------------------#

# ZONAL AVERAGE SALINITY SUBBASINS
#-----------------------------------------------------------------------------#
if "zo_salinity_atlmsk" in globals() and "zo_salinity_pacmsk" in globals() and "zo_salinity_indmsk" in globals():
    
    fig, ax = plt.subplots(nrows=3, ncols=1, figsize=(figXsize, figYsize))
    
    plt.axes(ax[0])
    map1, cont1 = plotZonalAveSubPlots(latGrid, -1*depthGrid, zo_salinity_atlmsk_mean, 
                                       contour_zosalin, cmap_zosalin, norm_zosalin)
    
    plt.axes(ax[1])
    map2, cont2 = plotZonalAveSubPlots(latGrid, -1*depthGrid, zo_salinity_pacmsk_mean, 
                                       contour_zosalin, cmap_zosalin, norm_zosalin)
    
    plt.axes(ax[2])
    map3, cont3 = plotZonalAveSubPlots(latGrid, -1*depthGrid, zo_salinity_indmsk_mean, 
                                       contour_zosalin, cmap_zosalin, norm_zosalin)
    
    ax[0].clabel(cont1,fmt=' {:.1f} '.format,fontsize='large')   
    ax[1].clabel(cont2,fmt=' {:.1f} '.format,fontsize='large')   
    ax[2].clabel(cont3,fmt=' {:.1f} '.format,fontsize='large') 
    
    # Sub titles
    ax[0].set_title('ATLANTIC BASIN')
    ax[1].set_title('PACIFIC BASIN')
    ax[2].set_title('INDIAN BASIN')
    
    cbar_ax = fig.add_axes([0.16, 0.01, 0.7, 0.04]) # list [x0, y0, width, height]
    cbar    = fig.colorbar(map1, cax=cbar_ax, orientation='horizontal',extend='both')
    cbar.ax.set_title(cbar_title_zosalin,size=cbar_label_size)
    cbar.ax.tick_params(labelsize=cbar_tick_size)
    
    if norm_zosalin is not None and len(np.unique(np.round(np.diff(norm_zosalin.boundaries),2))) != 1:
        plt.text(0.5,-0.15,'Warning: Adaptative colormap (non-linear) !',horizontalalignment='center',
                  verticalalignment='center', transform = ax[2].transAxes, fontsize=16, color='r', weight='bold')   
    

    filename = filePrefix + '_' + 'zoSalinity_subbasins.' + fig_format
    pathFilename = savedFigPath + filename
    if create_fig == 'y':                                               
        plt.savefig(pathFilename, format=fig_format, bbox_inches='tight', facecolor='white')
        filecount += 1; lnfile = 'file' + str(filecount) + '.' + fig_format
        os.system('ln -s ' + pathFilename + ' ' + lnfile)
    savedfiles.append(filename)
#-----------------------------------------------------------------------------#

# SST
#-----------------------------------------------------------------------------#  
if "SST" in globals():
    
    plotMapping(lon, lat, SST_mean, contour_sst, cbar_title_sst, land_mask, 
            projDataIn, projDataOut, cmap_sst, norm_sst, figTitleYear,
            figXsize, figYsize, cbar_label_size, cbar_tick_size, title_font_size)
    
    # Save figure
    filename = filePrefix + '_' + 'SST.' + fig_format
    pathFilename = savedFigPath + filename
    if create_fig == 'y':                                               
        plt.savefig(pathFilename, format=fig_format, bbox_inches='tight', facecolor='white')
        filecount += 1; lnfile = 'file' + str(filecount) + '.' + fig_format
        os.system('ln -s ' + pathFilename + ' ' + lnfile)
    savedfiles.append(filename)

    # GIF
    #-------------------#
    if create_gif == 'y':
        
        for ind in np.arange(0,SST.shape[0]):
            
            fig_title = figTitleMonth[ind]

            plotMapping(lon, lat, SST[ind,:,:], contour_sst, cbar_title_sst, land_mask, 
                    projDataIn, projDataOut, cmap_sst, norm_sst, fig_title,
                    figXsize, figYsize, cbar_label_size, cbar_tick_size, title_font_size)

            # Save figure
            ind = str(ind); ind = ind.zfill(2) # add 0 when ind is 1 digit long (1 => 01)
            pathFilename = savedGifPath + ind + '.' + fig_format
            plt.savefig(pathFilename, bbox_inches='tight', facecolor='white')

            plt.close()

        makeGif(savedGifPath, "SST.gif", timeInterv)
    #-------------------#
#-----------------------------------------------------------------------------#

# ZONAL AVERAGE TEMPERATURE
#-----------------------------------------------------------------------------#
#-------------------< Test custom cmap >----------------------#
# bounds = np.arange(-2,5,1)
# bounds = np.append(bounds,np.arange(5,10,0.2))
# bounds = np.append(bounds,np.arange(10,37,2))

# cmap = plt.get_cmap(cmapColor_temp,len(bounds)-1)
# norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
#-------------------------------------------------------------#
if "zo_temp" in globals():

    plotZonalAve(latGrid, -1*depthGrid, zo_temp_mean, contour_zotemp, cbar_title_zotemp, 
                 cmap_zotemp, norm_zotemp, figTitleYear, figXsize, figYsize, cbar_label_size, cbar_tick_size, 
                 title_font_size, xy_label_font_size, xy_ticks_font_size)
    
    filename = filePrefix + '_' + 'zoTemp.' + fig_format
    pathFilename = savedFigPath + filename
    if create_fig == 'y':                                               
        plt.savefig(pathFilename, format=fig_format, bbox_inches='tight', facecolor='white')
        filecount += 1; lnfile = 'file' + str(filecount) + '.' + fig_format
        os.system('ln -s ' + pathFilename + ' ' + lnfile)
    savedfiles.append(filename)

    # GIF
    #-------------------#
    if create_gif == 'y':
        
        for ind in np.arange(0,zo_temp.shape[0]):
            
            fig_title = figTitleMonth[ind]

            plotZonalAve(latGrid, -1*depthGrid, zo_temp[ind,:,:], contour_zotemp, cbar_title_zotemp, 
                         cmap_zotemp, norm_zotemp, fig_title, figXsize, figYsize, cbar_label_size, 
                         cbar_tick_size, title_font_size, xy_label_font_size, xy_ticks_font_size)

            # Save figure
            ind = str(ind); ind = ind.zfill(2) # add 0 when ind is 1 digit long (1 => 01)
            pathFilename = savedGifPath + ind + '.' + fig_format
            plt.savefig(pathFilename, bbox_inches='tight', facecolor='white')

            plt.close()

        makeGif(savedGifPath, "zo_temp.gif", timeInterv)
    #-------------------#
#-----------------------------------------------------------------------------#

# ZONAL AVERAGE Temperature SUBBASINS
#-----------------------------------------------------------------------------#
#-------------------< Test custom cmap >----------------------#
# bounds = np.arange(-2,5,1)
# bounds = np.append(bounds,np.arange(5,10,0.2))
# bounds = np.append(bounds,np.arange(10,37,2))

# cmap = plt.get_cmap(cmapColor_temp,len(bounds)-1)
# norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
#-------------------------------------------------------------#
if "zo_temp_atlmsk" in globals() and "zo_temp_pacmsk" in globals() and "zo_temp_indmsk" in globals():

    fig, ax = plt.subplots(nrows=3, ncols=1, figsize=(figXsize, figYsize))
    
    plt.axes(ax[0])
    map1, cont1 = plotZonalAveSubPlots(latGrid, -1*depthGrid,   zo_temp_atlmsk_mean, 
                                       contour_zotemp, cmap_zotemp, norm_zotemp)
    
    plt.axes(ax[1])
    map2, cont2 = plotZonalAveSubPlots(latGrid, -1*depthGrid,   zo_temp_pacmsk_mean, 
                                       contour_zotemp, cmap_zotemp, norm_zotemp)
    
    plt.axes(ax[2])
    map3, cont3 = plotZonalAveSubPlots(latGrid, -1*depthGrid,   zo_temp_indmsk_mean, 
                                       contour_zotemp, cmap_zotemp, norm_zotemp)
    
    # Add labels over contour lines
    ax[0].clabel(cont1,fmt=' {:.1f} '.format,fontsize='large')   
    ax[1].clabel(cont2,fmt=' {:.1f} '.format,fontsize='large')   
    ax[2].clabel(cont3,fmt=' {:.1f} '.format,fontsize='large')   
    
    # Sub titles
    ax[0].set_title('ATLANTIC BASIN')
    ax[1].set_title('PACIFIC BASIN')
    ax[2].set_title('INDIAN BASIN')
    
    cbar_ax = fig.add_axes([0.16, 0.01, 0.7, 0.04]) # list [x0, y0, width, height]    
    cbar    = fig.colorbar(map1, cax=cbar_ax, orientation='horizontal',extend='both')
    cbar.ax.set_title(cbar_title_zotemp,size=cbar_label_size)
    cbar.ax.tick_params(labelsize=cbar_tick_size)
    
    if norm_zotemp is not None and len(np.unique(np.round(np.diff(norm_zotemp.boundaries),2))) != 1:
            plt.text(0.5,-0.15,'Warning: Adaptative colormap (non-linear) !',horizontalalignment='center',
                     verticalalignment='center', transform = ax[2].transAxes, fontsize=16, color='r', weight='bold')  
    
    filename = filePrefix + '_' + 'zoTemp_subbasins.' + fig_format
    pathFilename = savedFigPath + filename
    if create_fig == 'y':                                               
        plt.savefig(pathFilename, format=fig_format, bbox_inches='tight', facecolor='white')
        filecount += 1; lnfile = 'file' + str(filecount) + '.' + fig_format
        os.system('ln -s ' + pathFilename + ' ' + lnfile)
    savedfiles.append(filename)
#-----------------------------------------------------------------------------#

# ZONAL AVERAGE SST
#-----------------------------------------------------------------------------#
if "zo_temp" in globals():

    fig, ax = plt.subplots(figsize=(figXsize, figYsize))
    
    major_yticks = np.arange(0, 35+1E-9, 5)
    major_xticks = np.arange(-90, 80+1E-9, 20)
    
    ax.grid(which='both')
    plt.xlim([-90,90])
    plt.ylim([0,35])
    
    map1 = plt.step(zo_lat, zo_temp_mean[1,:], '-k', where='post')
    
    # Titles/Labels and size
    plt.title ('Zonal SST - Yearly average' , fontsize = title_font_size)
    plt.xlabel('Latitude (°N)'  , fontsize = xy_label_font_size)
    plt.ylabel('Temperature (°C)'      , fontsize = xy_label_font_size)
    plt.xticks(major_xticks , fontsize = xy_ticks_font_size)
    plt.yticks(major_yticks , fontsize = xy_ticks_font_size)
    # plt.ylim(c_lev_zosst[0],c_lev_zosst[1])
               
    filename = filePrefix + '_' + 'zoSST.' + fig_format
    pathFilename = savedFigPath + filename
    if create_fig == 'y':                                               
        plt.savefig(pathFilename, format=fig_format, bbox_inches='tight', facecolor='white')
        filecount += 1; lnfile = 'file' + str(filecount) + '.' + fig_format
        os.system('ln -s ' + pathFilename + ' ' + lnfile)
    savedfiles.append(filename)
#-----------------------------------------------------------------------------#

# ZONAL AVERAGE STREAM FUNCTION
#-----------------------------------------------------------------------------#
if "zo_stream_function" in globals():

    plotZonalAve(latGrid, -1*depthGrid, zo_stream_function_mean, contour_zostrf, cbar_title_zostrf, 
                 cmap_zostrf, norm_zostrf, figTitleYear, figXsize, figYsize, cbar_label_size, cbar_tick_size, 
                 title_font_size, xy_label_font_size, xy_ticks_font_size)
    
    filename = filePrefix + '_' + 'zoStreamFunc.' + fig_format
    pathFilename = savedFigPath + filename
    if create_fig == 'y':                                               
        plt.savefig(pathFilename, format=fig_format, bbox_inches='tight', facecolor='white')
        filecount += 1; lnfile = 'file' + str(filecount) + '.' + fig_format
        os.system('ln -s ' + pathFilename + ' ' + lnfile)
    savedfiles.append(filename)

    # GIF
    #-------------------#
    if create_gif == 'y':
        
        for ind in np.arange(0,zo_stream_function.shape[0]):
            
            fig_title = figTitleMonth[ind]

            plotZonalAve(latGrid, -1*depthGrid, zo_stream_function[ind,:,:], contour_zostrf, cbar_title_zostrf, 
                         cmap_zostrf, norm_zostrf, fig_title, figXsize, figYsize, cbar_label_size, 
                         cbar_tick_size, title_font_size, xy_label_font_size, xy_ticks_font_size)

            # Save figure
            ind = str(ind); ind = ind.zfill(2) # add 0 when ind is 1 digit long (1 => 01)
            pathFilename = savedGifPath + ind + '.' + fig_format
            plt.savefig(pathFilename, bbox_inches='tight', facecolor='white')

            plt.close()

        makeGif(savedGifPath, "zo_stream_function.gif", timeInterv)
    #-------------------# 
#-----------------------------------------------------------------------------#

# BAROTROPIC STREAM FUNCTION
#-----------------------------------------------------------------------------#
if "baro_stream_function" in globals():
    
    plotMapping(lon, lat, baro_stream_function_mean, contour_bstrf, cbar_title_bstrf, land_mask, 
            projDataIn, projDataOut, cmap_bstrf, norm_bstrf, figTitleYear,
            figXsize, figYsize, cbar_label_size, cbar_tick_size, title_font_size)
    
    # Save figure
    filename = filePrefix + '_' + 'baroStreamFunc.' + fig_format
    pathFilename = savedFigPath + filename
    if create_fig == 'y':                                               
        plt.savefig(pathFilename, format=fig_format, bbox_inches='tight', facecolor='white')
        filecount += 1; lnfile = 'file' + str(filecount) + '.' + fig_format
        os.system('ln -s ' + pathFilename + ' ' + lnfile)
    savedfiles.append(filename)

    # GIF
    #-------------------#
    if create_gif == 'y':
        
        for ind in np.arange(0,baro_stream_function.shape[0]):
            
            fig_title = figTitleMonth[ind]

            plotMapping(lon, lat, baro_stream_function[ind,:,:], contour_bstrf, cbar_title_bstrf, land_mask, 
                    projDataIn, projDataOut, cmap_bstrf, norm_bstrf, fig_title,
                    figXsize, figYsize, cbar_label_size, cbar_tick_size, title_font_size)

            # Save figure
            ind = str(ind); ind = ind.zfill(2) # add 0 when ind is 1 digit long (1 => 01)
            pathFilename = savedGifPath + ind + '.' + fig_format
            plt.savefig(pathFilename, bbox_inches='tight', facecolor='white')

            plt.close()

        makeGif(savedGifPath, "baro_stream_function.gif", timeInterv)
    #-------------------#
#-----------------------------------------------------------------------------#    

# OCEAN MIXED LAYER
#-----------------------------------------------------------------------------#
if "omlmaxNH" in globals():

    #-----------------------< Northern Hemisphere >-----------------------#    
    lon_lim_nh = [-180,180]
    lat_lim_nh = [45,90]
    
    plotMappingZoom(lon, lat, omlmaxNH, contour_omlnh, cbar_title_omlnh, 
                    lon_lim_nh, lat_lim_nh, land_mask, 
                    projDataIn, ccrs.AzimuthalEquidistant(central_longitude=0,central_latitude=90), 
                    cmap_omlnh, norm_omlnh, fig_title_omlnh, figXsize, figYsize, cbar_label_size, cbar_tick_size, title_font_size)
    
    
    # Save figure
    filename = filePrefix + '_' + 'omlmaxNH.' + fig_format
    pathFilename = savedFigPath + filename
    if create_fig == 'y':                                               
        plt.savefig(pathFilename, format=fig_format, bbox_inches='tight', facecolor='white')
        filecount += 1; lnfile = 'file' + str(filecount) + '.' + fig_format
        os.system('ln -s ' + pathFilename + ' ' + lnfile)
    savedfiles.append(filename)

    #-----------------------< Southern Hemisphere >-----------------------#    
    lon_lim_sh = [-180,180]
    lat_lim_sh = [-90,-45]
    
    plotMappingZoom(lon, lat, omlmaxSH, contour_omlsh, cbar_title_omlsh, 
                    lon_lim_sh, lat_lim_sh, land_mask, 
                    projDataIn, ccrs.AzimuthalEquidistant(central_longitude=0,central_latitude=-90), 
                    cmap_omlsh, norm_omlsh, fig_title_omlsh, figXsize, figYsize, cbar_label_size, cbar_tick_size, title_font_size)
    
    # Save figure
    filename = filePrefix + '_' + 'omlmaxSH.' + fig_format
    pathFilename = savedFigPath + filename
    if create_fig == 'y':                                               
        plt.savefig(pathFilename, format=fig_format, bbox_inches='tight', facecolor='white')
        filecount += 1; lnfile = 'file' + str(filecount) + '.' + fig_format
        os.system('ln -s ' + pathFilename + ' ' + lnfile)
    savedfiles.append(filename)
#-----------------------------------------------------------------------------#

# INTPP
#-----------------------------------------------------------------------------#
if "TPP" in globals():
    
    plotMapping(lon, lat, INTPP_mean, contour_intpp, cbar_title_intpp, land_mask, 
                projDataIn, projDataOut, cmap_intpp, norm_intpp, figTitleYear,
                figXsize, figYsize, cbar_label_size, cbar_tick_size, title_font_size)
    
    # Save figure
    filename = filePrefix + '_' + 'intpp.' + fig_format
    pathFilename = savedFigPath + filename
    if create_fig == 'y':                                               
        plt.savefig(pathFilename, format=fig_format, bbox_inches='tight', facecolor='white')
        filecount += 1; lnfile = 'file' + str(filecount) + '.' + fig_format
        os.system('ln -s ' + pathFilename + ' ' + lnfile)
    savedfiles.append(filename)

    # GIF
    #-------------------#
    if create_gif == 'y':
        
        for ind in np.arange(0,INTPP.shape[0]):
            
            fig_title = figTitleMonth[ind]

            plotMapping(lon, lat, INTPP[ind,:,:], contour_intpp, cbar_title_intpp, land_mask, 
                    projDataIn, projDataOut, cmap_intpp, norm_intpp, fig_title,
                    figXsize, figYsize, cbar_label_size, cbar_tick_size, title_font_size)

            # Save figure
            ind = str(ind); ind = ind.zfill(2) # add 0 when ind is 1 digit long (1 => 01)
            pathFilename = savedGifPath + ind + '.' + fig_format
            plt.savefig(pathFilename, bbox_inches='tight', facecolor='white')

            plt.close()

        makeGif(savedGifPath, "INTPP.gif", timeInterv)
    #-------------------#
#-----------------------------------------------------------------------------#

# EPC100
#-----------------------------------------------------------------------------#
if "EPC100" in globals():

    plotMapping(lon, lat, EPC100_mean, contour_epc100, cbar_title_epc100, land_mask, 
                projDataIn, projDataOut, cmap_epc100, norm_epc100, figTitleYear,
                figXsize, figYsize, cbar_label_size, cbar_tick_size, title_font_size)
    
    # Save figure
    filename = filePrefix + '_' + 'epc100.' + fig_format
    pathFilename = savedFigPath + filename
    if create_fig == 'y':                                               
        plt.savefig(pathFilename, format=fig_format, bbox_inches='tight', facecolor='white')
        filecount += 1; lnfile = 'file' + str(filecount) + '.' + fig_format
        os.system('ln -s ' + pathFilename + ' ' + lnfile)
    savedfiles.append(filename)

    # GIF
    #-------------------#
    if create_gif == 'y':
        
        for ind in np.arange(0,EPC100.shape[0]):
            
            fig_title = figTitleMonth[ind]

            plotMapping(lon, lat, EPC100[ind,:,:], contour_epc100, cbar_title_epc100, land_mask, 
                    projDataIn, projDataOut, cmap_epc100, norm_epc100, fig_title,
                    figXsize, figYsize, cbar_label_size, cbar_tick_size, title_font_size)

            # Save figure
            ind = str(ind); ind = ind.zfill(2) # add 0 when ind is 1 digit long (1 => 01)
            pathFilename = savedGifPath + ind + '.' + fig_format
            plt.savefig(pathFilename, bbox_inches='tight', facecolor='white')

            plt.close()

        makeGif(savedGifPath, "EPC100.gif", timeInterv)
    #-------------------#
#-----------------------------------------------------------------------------#

# PO4
#-----------------------------------------------------------------------------#
if "PO4" in globals():

    # SUBPLOTS PO4
    #-------------------------------------------------#

    plotMappingLev(lon, lat, PO4_mean, cbar_title_po4, deptht, land_mask, 
                   projDataIn, projDataOut, 
                   cmap_po4, norm_zopo4, cbar_label_size, cbar_tick_size)
    
    
    # Save figure
    filename = filePrefix + '_' + 'PO4_lev.' + fig_format
    pathFilename = savedFigPath + filename
    if create_fig == 'y':                                               
        plt.savefig(pathFilename, format=fig_format, bbox_inches='tight', facecolor='white')
        filecount += 1; lnfile = 'file' + str(filecount) + '.' + fig_format
        os.system('ln -s ' + pathFilename + ' ' + lnfile)
    savedfiles.append(filename)
    #-------------------------------------------------#

    # ZONAL AVERAGE PO4
    #-------------------------------------------------#
    plotZonalAve(latGrid, -1*depthGrid, zo_PO4_mean, contour_zopo4, cbar_title_zopo4, 
                 cmap_zopo4, norm_zopo4, figTitleYear, figXsize, figYsize, cbar_label_size, cbar_tick_size, 
                 title_font_size, xy_label_font_size, xy_ticks_font_size)
    
    filename = filePrefix + '_' + 'zoPO4.' + fig_format
    pathFilename = savedFigPath + filename
    if create_fig == 'y':                                               
        plt.savefig(pathFilename, format=fig_format, bbox_inches='tight', facecolor='white')
        filecount += 1; lnfile = 'file' + str(filecount) + '.' + fig_format
        os.system('ln -s ' + pathFilename + ' ' + lnfile)
    savedfiles.append(filename)

    # GIF
    #-------------------#
    if create_gif == 'y':
        
        for ind in np.arange(0,zo_PO4.shape[0]):
            
            fig_title = figTitleMonth[ind]

            plotZonalAve(latGrid, -1*depthGrid, zo_PO4[ind,:,:], contour_zopo4, cbar_title_zopo4, 
                         cmap_zopo4, norm_zopo4, fig_title, figXsize, figYsize, cbar_label_size, 
                         cbar_tick_size, title_font_size, xy_label_font_size, xy_ticks_font_size)

            # Save figure
            ind = str(ind); ind = ind.zfill(2) # add 0 when ind is 1 digit long (1 => 01)
            pathFilename = savedGifPath + ind + '.' + fig_format
            plt.savefig(pathFilename, bbox_inches='tight', facecolor='white')

            plt.close()

        makeGif(savedGifPath, "zo_PO4.gif", timeInterv)
    #-------------------#
    #-------------------------------------------------#
    
    # ZONAL AVERAGE PO4 SUBBASINS
    #-------------------------------------------------#
    fig, ax = plt.subplots(nrows=3, ncols=1, figsize=(figXsize, figYsize))
    
    plt.axes(ax[0])
    map1, cont1 = plotZonalAveSubPlots(latGrid, -1*depthGrid, zo_PO4_atlmsk_mean, contour_zopo4, cmap_po4, norm_zopo4)
    
    plt.axes(ax[1])
    map2, cont2 = plotZonalAveSubPlots(latGrid, -1*depthGrid, zo_PO4_pacmsk_mean, contour_zopo4, cmap_po4, norm_zopo4)
    
    plt.axes(ax[2])
    map3, cont3 = plotZonalAveSubPlots(latGrid, -1*depthGrid, zo_PO4_indmsk_mean, contour_zopo4, cmap_po4, norm_zopo4)
    
    # Add labels over contour lines
    ax[0].clabel(cont1,fmt=' {:.1f} '.format,fontsize='large')   
    ax[1].clabel(cont2,fmt=' {:.1f} '.format,fontsize='large')   
    ax[2].clabel(cont3,fmt=' {:.1f} '.format,fontsize='large')   
    
    # Sub titles
    ax[0].set_title('ATLANTIC BASIN')
    ax[1].set_title('PACIFIC BASIN')
    ax[2].set_title('INDIAN BASIN')
    
    cbar_ax = fig.add_axes([0.16, 0.01, 0.7, 0.04]) # list [x0, y0, width, height]    
    cbar    = fig.colorbar(map1, cax=cbar_ax, orientation='horizontal',extend='both')
    cbar.ax.set_title(cbar_title_zopo4,size=cbar_label_size)
    cbar.ax.tick_params(labelsize=cbar_tick_size)
    
    if norm_po4 is not None and len(np.unique(np.round(np.diff(norm_po4.boundaries),2))) != 1:
            plt.text(0.5,-0.15,'Warning: Adaptative colormap (non-linear) !',horizontalalignment='center',
                     verticalalignment='center', transform = ax[2].transAxes, fontsize=16, color='r', weight='bold')
    
    filename = filePrefix + '_' + 'zoPO4_subbasins.' + fig_format
    pathFilename = savedFigPath + filename
    if create_fig == 'y':                                               
        plt.savefig(pathFilename, format=fig_format, bbox_inches='tight', facecolor='white')
        filecount += 1; lnfile = 'file' + str(filecount) + '.' + fig_format
        os.system('ln -s ' + pathFilename + ' ' + lnfile)
    savedfiles.append(filename)
    #-------------------------------------------------#
#-----------------------------------------------------------------------------#

# NO3
#-----------------------------------------------------------------------------#
if "NO3" in globals():

    # SUBPLOTS NO3
    #-------------------------------------------------#
    plotMappingLev(lon, lat, NO3_mean, cbar_title_no3, deptht, land_mask, 
                   projDataIn, projDataOut, 
                   cmap_no3, norm_no3, cbar_label_size, cbar_tick_size)
    
    # Save figure
    filename = filePrefix + '_' + 'NO3_lev.' + fig_format
    pathFilename = savedFigPath + filename
    if create_fig == 'y':                                               
        plt.savefig(pathFilename, format=fig_format, bbox_inches='tight', facecolor='white')
        filecount += 1; lnfile = 'file' + str(filecount) + '.' + fig_format
        os.system('ln -s ' + pathFilename + ' ' + lnfile)
    savedfiles.append(filename)
    #-------------------------------------------------#
    
    # ZONAL AVERAGE NO3
    #-------------------------------------------------#
    plotZonalAve(latGrid, -1*depthGrid, zo_NO3_mean, contour_zono3, cbar_title_zono3, 
                 cmap_no3, norm_no3, figTitleYear, figXsize, figYsize, cbar_label_size, cbar_tick_size, 
                 title_font_size, xy_label_font_size, xy_ticks_font_size)
    
    filename = filePrefix + '_' + 'zoNO3.' + fig_format
    pathFilename = savedFigPath + filename
    if create_fig == 'y':                                               
        plt.savefig(pathFilename, format=fig_format, bbox_inches='tight', facecolor='white')
        filecount += 1; lnfile = 'file' + str(filecount) + '.' + fig_format
        os.system('ln -s ' + pathFilename + ' ' + lnfile)
    savedfiles.append(filename)

    # GIF
    #-------------------#
    if create_gif == 'y':
        
        for ind in np.arange(0,zo_NO3.shape[0]):
            
            fig_title = figTitleMonth[ind]

            plotZonalAve(latGrid, -1*depthGrid, zo_NO3[ind,:,:], contour_zono3, cbar_title_zono3, 
                         cmap_zono3, norm_zono3, fig_title, figXsize, figYsize, cbar_label_size, 
                         cbar_tick_size, title_font_size, xy_label_font_size, xy_ticks_font_size)

            # Save figure
            ind = str(ind); ind = ind.zfill(2) # add 0 when ind is 1 digit long (1 => 01)
            pathFilename = savedGifPath + ind + '.' + fig_format
            plt.savefig(pathFilename, bbox_inches='tight', facecolor='white')

            plt.close()

        makeGif(savedGifPath, "zo_NO3.gif", timeInterv)
    #-------------------#
    #-------------------------------------------------#
    
    # ZONAL AVERAGE NO3 SUBBASINS
    #-------------------------------------------------#
    fig, ax = plt.subplots(nrows=3, ncols=1, figsize=(figXsize, figYsize))
    
    plt.axes(ax[0])
    map1, cont1 = plotZonalAveSubPlots(latGrid, -1*depthGrid, zo_NO3_atlmsk_mean, contour_zono3, cmap_no3, norm_no3)
    
    plt.axes(ax[1])
    map2, cont2 = plotZonalAveSubPlots(latGrid, -1*depthGrid, zo_NO3_pacmsk_mean, contour_zono3, cmap_no3, norm_no3)
    
    plt.axes(ax[2])
    map3, cont3 = plotZonalAveSubPlots(latGrid, -1*depthGrid, zo_NO3_indmsk_mean, contour_zono3, cmap_no3, norm_no3)
    
    # Add labels over contour lines
    ax[0].clabel(cont1,fmt=' {:.1f} '.format,fontsize='large')   
    ax[1].clabel(cont2,fmt=' {:.1f} '.format,fontsize='large')   
    ax[2].clabel(cont3,fmt=' {:.1f} '.format,fontsize='large')   
    
    # Sub titles
    ax[0].set_title('ATLANTIC BASIN')
    ax[1].set_title('PACIFIC BASIN')
    ax[2].set_title('INDIAN BASIN')
    
    cbar_ax = fig.add_axes([0.16, 0.01, 0.7, 0.04]) # list [x0, y0, width, height]    
    cbar    = fig.colorbar(map1, cax=cbar_ax, orientation='horizontal',extend='both')
    cbar.ax.set_title(cbar_title_zono3,size=cbar_label_size)
    cbar.ax.tick_params(labelsize=cbar_tick_size)
    
    if norm_no3 is not None and len(np.unique(np.round(np.diff(norm_no3.boundaries),2))) != 1:
            plt.text(0.5,-0.15,'Warning: Adaptative colormap (non-linear) !',horizontalalignment='center',
                     verticalalignment='center', transform = ax[2].transAxes, fontsize=16, color='r', weight='bold')
    
    filename = filePrefix + '_' + 'zoNO3_subbasins.' + fig_format
    pathFilename = savedFigPath + filename
    if create_fig == 'y':                                               
        plt.savefig(pathFilename, format=fig_format, bbox_inches='tight', facecolor='white')
        filecount += 1; lnfile = 'file' + str(filecount) + '.' + fig_format
        os.system('ln -s ' + pathFilename + ' ' + lnfile)
    savedfiles.append(filename)
    #-------------------------------------------------#
#-----------------------------------------------------------------------------#

# O2
#-----------------------------------------------------------------------------#
if "O2" in globals():

    # SUBPLOTS O2
    #-------------------------------------------------#
    plotMappingLev(lon, lat, O2_mean, cbar_title_o2, deptht, land_mask, 
                   projDataIn, projDataOut, 
                   cmap_o2, norm_o2, cbar_label_size, cbar_tick_size)
    
    # Save figure
    filename = filePrefix + '_' + 'O2_lev.' + fig_format
    pathFilename = savedFigPath + filename
    if create_fig == 'y':                                               
        plt.savefig(pathFilename, format=fig_format, bbox_inches='tight', facecolor='white')
        filecount += 1; lnfile = 'file' + str(filecount) + '.' + fig_format
        os.system('ln -s ' + pathFilename + ' ' + lnfile)
    savedfiles.append(filename)
    #-------------------------------------------------#
    
    # ZONAL AVERAGE O2
    #-------------------------------------------------#
    plotZonalAve(latGrid, -1*depthGrid, zo_O2_mean, contour_zoo2, cbar_title_zoo2, 
                 cmap_o2, norm_o2, figTitleYear, figXsize, figYsize, cbar_label_size, cbar_tick_size, 
                 title_font_size, xy_label_font_size, xy_ticks_font_size)
    
    cont1 = plt.contour(latGrid, -1*depthGrid, zo_O2_mean, np.array([6.5, 65]), colors='w', linewidths=0.7) # Add dysoxia and anoxia contour lines
    plt.clabel(cont1, fmt=' {:.1f} '.format, fontsize='x-large')                                       # Add labels over contour lines
            
    filename = filePrefix + '_' + 'zoO2.' + fig_format
    pathFilename = savedFigPath + filename
    if create_fig == 'y':                                               
        plt.savefig(pathFilename, format=fig_format, bbox_inches='tight', facecolor='white')
        filecount += 1; lnfile = 'file' + str(filecount) + '.' + fig_format
        os.system('ln -s ' + pathFilename + ' ' + lnfile)
    savedfiles.append(filename)

    # GIF
    #-------------------#
    if create_gif == 'y':
        
        for ind in np.arange(0,zo_O2.shape[0]):
            
            fig_title = figTitleMonth[ind]

            plotZonalAve(latGrid, -1*depthGrid, zo_O2[ind,:,:], contour_zoo2, cbar_title_zoo2, 
                         cmap_zoo2, norm_zoo2, fig_title, figXsize, figYsize, cbar_label_size, 
                         cbar_tick_size, title_font_size, xy_label_font_size, xy_ticks_font_size)

            # Save figure
            ind = str(ind); ind = ind.zfill(2) # add 0 when ind is 1 digit long (1 => 01)
            pathFilename = savedGifPath + ind + '.' + fig_format
            plt.savefig(pathFilename, bbox_inches='tight', facecolor='white')

            plt.close()

        makeGif(savedGifPath, "zo_O2.gif", timeInterv)
    #-------------------#
    #-------------------------------------------------#
    
    # ZONAL AVERAGE O2 SUBBASINS
    #-------------------------------------------------#
    fig, ax = plt.subplots(nrows=3, ncols=1, figsize=(figXsize, figYsize))
    
    plt.axes(ax[0])
    map1, cont1 = plotZonalAveSubPlots(latGrid, -1*depthGrid, zo_O2_atlmsk_mean, contour_zoo2, cmap_o2, norm_o2)
    cont11      = plt.contour(latGrid, -1*depthGrid, zo_O2_atlmsk_mean, np.array([6.5, 65]), colors='w', linewidths=0.7) # Add dysoxia and anoxia contour lines
    
    plt.axes(ax[1])
    map2, cont2 = plotZonalAveSubPlots(latGrid, -1*depthGrid, zo_O2_pacmsk_mean, contour_zoo2, cmap_o2, norm_o2)
    cont22      = plt.contour(latGrid, -1*depthGrid, zo_O2_pacmsk_mean, np.array([6.5, 65]), colors='w', linewidths=0.7) # Add dysoxia and anoxia contour lines
    
    plt.axes(ax[2])
    map3, cont3 = plotZonalAveSubPlots(latGrid, -1*depthGrid, zo_O2_indmsk_mean, contour_zoo2, cmap_o2, norm_o2)
    cont33      = plt.contour(latGrid, -1*depthGrid, zo_O2_indmsk_mean, np.array([6.5, 65]), colors='w', linewidths=0.7) # Add dysoxia and anoxia contour lines
    
    # Add labels over contour lines
    ax[0].clabel(cont1,fmt=' {:.1f} '.format,fontsize='large')  
    ax[0].clabel(cont11, fmt=' {:.1f} '.format, fontsize='x-large')
    ax[1].clabel(cont2,fmt=' {:.1f} '.format,fontsize='large')  
    ax[1].clabel(cont22, fmt=' {:.1f} '.format, fontsize='x-large')
    ax[2].clabel(cont3,fmt=' {:.1f} '.format,fontsize='large')   
    ax[2].clabel(cont33, fmt=' {:.1f} '.format, fontsize='x-large')
    
    # Sub titles
    ax[0].set_title('ATLANTIC BASIN')
    ax[1].set_title('PACIFIC BASIN')
    ax[2].set_title('INDIAN BASIN')
    
    cbar_ax = fig.add_axes([0.16, 0.01, 0.7, 0.04]) # list [x0, y0, width, height]    
    cbar    = fig.colorbar(map1, cax=cbar_ax, orientation='horizontal',extend='both')
    cbar.ax.set_title(cbar_title_zoo2,size=cbar_label_size)
    cbar.ax.tick_params(labelsize=cbar_tick_size)
    
    if norm_o2 is not None and len(np.unique(np.round(np.diff(norm_o2.boundaries),2))) != 1:
            plt.text(0.5,-0.15,'Warning: Adaptative colormap (non-linear) !',horizontalalignment='center',
                     verticalalignment='center', transform = ax[2].transAxes, fontsize=16, color='r', weight='bold')
    
    filename = filePrefix + '_' + 'zoO2_subbasins.' + fig_format
    pathFilename = savedFigPath + filename
    if create_fig == 'y':                                               
        plt.savefig(pathFilename, format=fig_format, bbox_inches='tight', facecolor='white')
        filecount += 1; lnfile = 'file' + str(filecount) + '.' + fig_format
        os.system('ln -s ' + pathFilename + ' ' + lnfile)
    savedfiles.append(filename)
    #-------------------------------------------------#
#-----------------------------------------------------------------------------#

#%%===========================================================================#
#                            --< PDF GENERATION >--                           #
#=============================================================================#

if create_pdf == 'y':
    
    dopdf(savedPdfPath, savedfiles, 'file', filecount, filePrefix, sentence_to_use)

# %%
