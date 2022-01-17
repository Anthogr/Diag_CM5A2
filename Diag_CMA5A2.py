#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
-------------------------------------------------------------------------------
Created on Mon Oct 11 14:34:22 2021

@author: anthog

This script plots several diagnostics from CMA5A2 outputs.
Depending on the options enabled:

    - Figures format can be set (png, pdf...)
    - PDF file gathering the generated figures can be created
    - Animated gif can be created for some variables (added 17/12/2021)
    - Manual value limits for plots colormap can be activated instead of the 
      automatic scaling for convenient comparison between several simulations
    - In the case where 2 experiments are chosen, there is the ability to 
      compute the difference between them

- Launch the script by running the RUN_DIAG.sh file (command 'sh RUN_DIAG.sh')

- Enter your [ProfileName] as requested
    - if [ProfileName] is not known: a folder [ProfileName] located in userData/
      containing a [ProfileName]_profile.py will be created and program will stop.
      At this point, you will need to edit the [ProfileName]_profile.py file 
      according to your needs and relaunch the program
    - if [ProfileName] is known: program will run following the settings of the
      userData/[ProfileName]/[ProfileName]_profile.py file

    Note that the variables in [ProfileName]_profile.py file: 
        - filePrefix
        - maskFile
        - bathyFile
        - subBasinFile
        - dataDirPath
    are in list format and can contain as many elements as you want (elements 
    being the different properties of the experiments: file names, path...)
    In the specific case where you choose 2 experiments, a new prompt will appear 
    letting you choose the ability to compute differences between simu 1 & 2.

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
from shutil import copyfile
import sys
import importlib

# Loading custom functions 
from Toolz import (z_masked_overlap, readncfile, dopdf, dispGridCoastline,
                   plotMapping, plotZonalAve, plotZonalAveSubPlots, plotMappingZoom, plotMappingLev,
                   adaptativeColorMap, makeGif, myCustomInterp3)

#%%===========================================================================#
#          --< PROFILE INITIALIZATION AND LOADING INIT VARIABLES>--           #
#=============================================================================#
# Initialization part that will ask for a profile name. 
#   - If this names hasn't been used yet, the script will create a folder 
#     containing all the user related data i.e:
#       * userName_profile.py containing the initialization variables
#       * Folders FIG/PDF/GIF that will contain the script outputs
#     Those are created in userData/UserName folder
#
#   - If profile name entered already exists, the script will use the related 
#     userName_profile.py file and run the script

print('\n\n\n================================< Diag_CMA5A2 >================================')

# Initialization
#------------------------------------------------------#
# Ask for profile name
promptProfile = f"""Enter your profile name:
(This is an identifier that will allow to save your settings to run the script)\n"""
profileName   = input(promptProfile)

# Extract path name that will contain all the user related data
profilePath = os.getcwd() + '/userData/' + profileName

str_created_profile = f"""\n-------------------------------------
Profile {profileName} has been created\n
1) You should now edit
   {profilePath}/{profileName}_profile.py
2) Relaunch the script
3) When prompted for the profile name, re-enter {profileName} so the script
   will run by using the newly edited {profileName}_profile.py file
-------------------------------------
===================================================================
"""

# Check if profileName folder exists, if not create one
if os.system('[ ! -d "' + profilePath + '" ]') == 0:
    promptProfile2 = f"\nProfile {profileName} doesn't exist, do you want to create it? (y/n)\n" 
    create_profile  = input(promptProfile2)

    if create_profile == 'y':
        os.system('mkdir ' + profilePath)
        copyfile('defaultProfile.py', f'userData/{profileName}/{profileName}_profile.py') 
        sys.exit(str_created_profile)
    elif create_profile == 'n':
        sys.exit('\nProgram ended\n===============================================================================')

# And if folder and profileName_profile.py exist: load variables in profileName_profile.py
elif os.system('[ -f "' + profilePath + f'/{profileName}_profile.py' + '" ]') == 0:
    varModule = importlib.import_module(f'userData.{profileName}.{profileName}_profile')
    for varname, val in varModule.__dict__.items():
        if (varname.startswith('__') and varname.endswith('__'))==False : # get the negation of the expression between parentethis
            exec(varname + " = varModule." + varname)
    print(f'\n------------------------------------------------------------------\nVariables in {profileName}_profile.py loaded\n------------------------------------------------------------------')

else:
     sys.exit(f'/!\ Warning, {profilePath}/{profileName}_profile.py doesn''t exist')
#------------------------------------------------------#

# Check if all list variables are the same length
#------------------------------------------------------#
if not len(sentence_to_use) == len(dataDirOutput) == len(filePrefix) == len(dataDirMask) == len(maskFile) == len(dataDirBathy) == len(bathyFile) == len(subBasinFile):
    err_str = f"""\n     /!\ Warning: Lengths of list variables are not equal /!\\
------------------------------------------------------------------
Check in userData/{profileName}/{profileName}_profile.py:
- sentence_to_use - dataDirOutput - filePrefix - dataDirMask
- maskFile - dataDirBathy - bathyFile - subBasinFile
Remember you can add a str element to a list variable and leave it
blank if path or file doesn't exist
------------------------------------------------------------------
===============================================================================

"""

    sys.exit(err_str)
#------------------------------------------------------#

# Create folders for figures / PDF / GIF if specidfied 
#------------------------------------------------------#
# Figures
savedFigPath = profilePath + '/FIG/' # get path name
if os.system('[ ! -d "' + savedFigPath + '" ]') == 0: # folder doesn't exist, create one
    os.system('mkdir ' + savedFigPath)

# PDF
if create_pdf == 'y':
    savedPdfPath = profilePath + '/PDF/' # get path name
    if os.system('[ ! -d "' + savedPdfPath + '" ]') == 0: # folder doesn't exist, create one
        os.system('mkdir ' + savedPdfPath)

# GIF
if create_gif == 'y':
    savedGifPath = profilePath + '/GIF/' # get path name
    if os.system('[ ! -d "' + savedGifPath + '" ]') == 0: # folder doesn't exist, create one
        os.system('mkdir ' + savedGifPath)
#------------------------------------------------------#

# Check if nb of simulations = 2. If yes, a prompt will ask if you want to make diff between
# simulation 1 & 2
#------------------------------------------------------#
if len(filePrefix) == 2:
    promptDiff  = f"\nCreate pdf difference: [{filePrefix[0]}] - [{filePrefix[1]}]? (y/n)\n"
    create_pdf_diff = input(promptDiff)
    while create_pdf_diff != 'y' and create_pdf_diff != 'n':
        print('/!\\ Answer not valid /!\\')
        create_pdf_diff = input(promptDiff)
    if create_pdf_diff == 'y':
        promptDiff2 = f"Notes to write in pdf:\n"
        sentence_to_use.append(input(promptDiff2))
else:
    create_pdf_diff = 'n'

if create_pdf_diff == 'y':
    length_loop = 3
else:
    length_loop = len(filePrefix)
#------------------------------------------------------#

for ind_file in np.arange(0,length_loop):
    if ((ind_file == 0 or ind_file == 1) and create_pdf_diff =='y') or create_pdf_diff =='n':
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
            'bathy':'Bathymetry',
            }

        dict_T = {
            'lon':'NaV_LoN',
            'lat':'nav_lat',
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

        # Add '/' at end of path if necessary
        #----------------------------------#
        # if path not empty and doesn't end by '/': add '/' at the end
        if dataDirMask[ind_file] and dataDirMask[ind_file][-1] != '/': 
            dataDirMask[ind_file] = dataDirMask[ind_file] + '/'
        if dataDirBathy[ind_file] and dataDirBathy[ind_file][-1] != '/': 
            dataDirBathy[ind_file] = dataDirBathy[ind_file] + '/'
        if dataDirOutput[ind_file] and dataDirOutput[ind_file][-1] != '/': 
            dataDirOutput[ind_file] = dataDirOutput[ind_file] + '/'
        #----------------------------------#

        # File names
        #----------------------------------#
        nc_file_mask     = dataDirMask[ind_file]  + maskFile[ind_file]
        nc_file_bathy    = dataDirBathy[ind_file]  + bathyFile[ind_file]
        nc_file_subbasin = dataDirBathy[ind_file]  + subBasinFile[ind_file]
        nc_file_T        = dataDirOutput[ind_file] + "OCE/Analyse/SE/" + filePrefix[ind_file] + "_grid_T.nc"
        nc_file_U        = dataDirOutput[ind_file] + "OCE/Analyse/SE/" + filePrefix[ind_file] + "_grid_U.nc"
        nc_file_V        = dataDirOutput[ind_file] + "OCE/Analyse/SE/" + filePrefix[ind_file] + "_grid_V.nc"
        nc_file_W        = dataDirOutput[ind_file] + "OCE/Analyse/SE/" + filePrefix[ind_file] + "_grid_W.nc"
        nc_file_diaptr   = dataDirOutput[ind_file] + "OCE/Analyse/SE/" + filePrefix[ind_file] + "_diaptr.nc"
        nc_file_diad_T   = dataDirOutput[ind_file] + "MBG/Analyse/SE/" + filePrefix[ind_file] + "_diad_T.nc"
        nc_file_ptrc_T   = dataDirOutput[ind_file] + "MBG/Analyse/SE/" + filePrefix[ind_file] + "_ptrc_T.nc"
        nc_file_histmth  = dataDirOutput[ind_file] + "ATM/Analyse/SE/" + filePrefix[ind_file] + "_histmth.nc"
        #----------------------------------#

        # Load variables in DATA var (check if file exist before loading)
        #----------------------------------#
        DATA = {}
        notLoaded = []

        if os.path.isfile(nc_file_mask):
            DATA = readncfile(nc_file_mask, dict_mask, DATA)
        else:
            notLoaded.append(f'{nc_file_mask} doesn''t exist')

        if os.path.isfile(nc_file_bathy):
            DATA = readncfile(nc_file_bathy, dict_bathy, DATA)
        else:
            notLoaded.append(f'{nc_file_bathy} doesn''t exist')

        if os.path.isfile(nc_file_subbasin):
            DATA = readncfile(nc_file_subbasin, dict_subbasin, DATA)
        else:
            notLoaded.append(f'{nc_file_subbasin} doesn''t exist')

        if os.path.isfile(nc_file_T):
            DATA = readncfile(nc_file_T, dict_T, DATA)
        else:
            notLoaded.append(f'{nc_file_T} doesn''t exist')

        if os.path.isfile(nc_file_U):
            DATA = readncfile(nc_file_U, dict_U, DATA)
        else:
            notLoaded.append(f'{nc_file_U} doesn''t exist')

        if os.path.isfile(nc_file_diaptr):
            DATA = readncfile(nc_file_diaptr, dict_diaptr, DATA)
        else:
            notLoaded.append(f'{nc_file_diaptr} doesn''t exist')

        if os.path.isfile(nc_file_diad_T):
            DATA = readncfile(nc_file_diad_T, dict_diad_T, DATA)
        else:
            notLoaded.append(f'{nc_file_diad_T} doesn''t exist')

        if os.path.isfile(nc_file_ptrc_T):
            DATA = readncfile(nc_file_ptrc_T, dict_ptrc_T, DATA)
        else:
            notLoaded.append(f'{nc_file_ptrc_T} doesn''t exist')

        if os.path.isfile(nc_file_histmth):
            DATA = readncfile(nc_file_histmth, dict_histmth, DATA)
        else:
            notLoaded.append(f'{nc_file_histmth} doesn''t exist')

        if notLoaded:
            print('\n               /!\ Warning: File(s) not found /!\\\n------------------------------------------------------------------')
            print(*notLoaded, sep = '\n')
            print('------------------------------------------------------------------')
            if len(notLoaded) == 9:
                sys.exit('''
                  /!\ None of the files have been found /!\\
===============================================================================
''')
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
        #----------------------------------#

        # Temperature
        #----------------------------------#
        if "SST" in globals():

            SST_mean = np.mean(SST,axis=0)  # Yearly average
        
        if "zo_temp" in globals():
            
            zo_temp      = np.squeeze(zo_temp)
            zo_temp      = np.ma.masked_where(zo_temp==0, zo_temp) # replace 0 by masked values
            zo_temp_mean = np.mean(zo_temp,axis=0) # Yearly average
        #----------------------------------#

        # Zonal stream function
        #----------------------------------#
        if "zo_stream_function" in globals():

            zo_stream_function      = np.squeeze(zo_stream_function)
            zo_stream_function      = np.ma.masked_where(zo_stream_function==0, zo_stream_function) # replace 0 by masked values
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
                if "subbasin_mask" in globals():
                    zo_temp_atlmsk = varZero.copy()
                    zo_temp_pacmsk = varZero.copy()
                    zo_temp_indmsk = varZero.copy()

            if "salinity" in globals():
                zo_salinity        = varZero.copy()
                if "subbasin_mask" in globals():
                    zo_salinity_atlmsk = varZero.copy()
                    zo_salinity_pacmsk = varZero.copy()
                    zo_salinity_indmsk = varZero.copy()       
                
            if "PO4" in globals():
                zo_PO4        = varZero.copy()
                if "subbasin_mask" in globals():
                    zo_PO4_atlmsk = varZero.copy()
                    zo_PO4_pacmsk = varZero.copy()
                    zo_PO4_indmsk = varZero.copy()

            if "NO3" in globals():
                zo_NO3        = varZero.copy()
                if "subbasin_mask" in globals():
                    zo_NO3_atlmsk = varZero.copy()
                    zo_NO3_pacmsk = varZero.copy()
                    zo_NO3_indmsk = varZero.copy()

            if "O2" in globals():
                zo_O2        = varZero.copy()
                if "subbasin_mask" in globals():
                    zo_O2_atlmsk = varZero.copy()
                    zo_O2_pacmsk = varZero.copy()
                    zo_O2_indmsk = varZero.copy()
            #-------------------------------#
            
            # Lat grid is not regular, we need to to compute zonal mean
            # For this, we take variable zo_lat as a guide
            #-------------------------------#
            for j in np.arange(0,149): 
                
                # limits
                if j < 148:
                    limit_inf = zo_lat[j] - (zo_lat[j]   - zo_lat[j-1]) /2
                    limit_sup = zo_lat[j] + (zo_lat[j+1] - zo_lat[j])   /2
                else:
                    limit_inf = zo_lat[j] - (zo_lat[j]   - zo_lat[j-1]) /2
                    limit_sup = zo_lat[j] + 1
                
                ind = np.where(np.logical_and(lat >= limit_inf, lat <= limit_sup))
                
                if "salinity" in globals():
                    zo_salinity[:,:,j]        = np.ma.mean(salinity[:,:,ind[0],ind[1]], axis=2)
                    if "subbasin_mask" in globals():
                        zo_salinity_atlmsk[:,:,j] = np.ma.mean(salinity[:,:,ind[0],ind[1]] * atlmsk[ind[0],ind[1]], axis=2)
                        zo_salinity_pacmsk[:,:,j] = np.ma.mean(salinity[:,:,ind[0],ind[1]] * pacmsk[ind[0],ind[1]], axis=2)
                        zo_salinity_indmsk[:,:,j] = np.ma.mean(salinity[:,:,ind[0],ind[1]] * indmsk[ind[0],ind[1]], axis=2)

                if "temp" in globals():
                    zo_temp2[:,:,j]        = np.ma.mean(temp[:,:,ind[0],ind[1]], axis=2)
                    if "subbasin_mask" in globals():
                        zo_temp_atlmsk[:,:,j]  = np.ma.mean(temp[:,:,ind[0],ind[1]] * atlmsk[ind[0],ind[1]], axis=2)
                        zo_temp_pacmsk[:,:,j]  = np.ma.mean(temp[:,:,ind[0],ind[1]] * pacmsk[ind[0],ind[1]], axis=2)
                        zo_temp_indmsk[:,:,j]  = np.ma.mean(temp[:,:,ind[0],ind[1]] * indmsk[ind[0],ind[1]], axis=2)
                
                if "PO4" in globals():
                    zo_PO4[:,:,j]        = np.ma.mean(PO4[:,:,ind[0],ind[1]], axis=2)
                    if "subbasin_mask" in globals():
                        zo_PO4_atlmsk[:,:,j] = np.ma.mean(PO4[:,:,ind[0],ind[1]] * atlmsk[ind[0],ind[1]], axis=2)
                        zo_PO4_pacmsk[:,:,j] = np.ma.mean(PO4[:,:,ind[0],ind[1]] * pacmsk[ind[0],ind[1]], axis=2)
                        zo_PO4_indmsk[:,:,j] = np.ma.mean(PO4[:,:,ind[0],ind[1]] * indmsk[ind[0],ind[1]], axis=2)
                
                if "NO3" in globals():
                    zo_NO3[:,:,j]        = np.ma.mean(NO3[:,:,ind[0],ind[1]], axis=2)
                    if "subbasin_mask" in globals():
                        zo_NO3_atlmsk[:,:,j] = np.ma.mean(NO3[:,:,ind[0],ind[1]] * atlmsk[ind[0],ind[1]], axis=2)
                        zo_NO3_pacmsk[:,:,j] = np.ma.mean(NO3[:,:,ind[0],ind[1]] * pacmsk[ind[0],ind[1]], axis=2)
                        zo_NO3_indmsk[:,:,j] = np.ma.mean(NO3[:,:,ind[0],ind[1]] * indmsk[ind[0],ind[1]], axis=2)
                
                if "O2" in globals():
                    zo_O2[:,:,j]        = np.ma.mean(O2[:,:,ind[0],ind[1]], axis=2)
                    if "subbasin_mask" in globals():
                        zo_O2_atlmsk[:,:,j] = np.ma.mean(O2[:,:,ind[0],ind[1]] * atlmsk[ind[0],ind[1]], axis=2)
                        zo_O2_pacmsk[:,:,j] = np.ma.mean(O2[:,:,ind[0],ind[1]] * pacmsk[ind[0],ind[1]], axis=2)
                        zo_O2_indmsk[:,:,j] = np.ma.mean(O2[:,:,ind[0],ind[1]] * indmsk[ind[0],ind[1]], axis=2)
            #-------------------------------#

            # Average over time
            #-------------------------------#
            if "salinity" in globals():
                zo_salinity_mean        = np.mean(zo_salinity, axis=0)
                if "subbasin_mask" in globals():
                    zo_salinity_atlmsk_mean = np.mean(zo_salinity_atlmsk, axis=0)
                    zo_salinity_pacmsk_mean = np.mean(zo_salinity_pacmsk, axis=0)
                    zo_salinity_indmsk_mean = np.mean(zo_salinity_indmsk, axis=0)

            if "temp" in globals():
                zo_temp2_mean       = np.mean(zo_temp2, axis=0)
                if "subbasin_mask" in globals():
                    zo_temp_atlmsk_mean = np.mean(zo_temp_atlmsk, axis=0)
                    zo_temp_pacmsk_mean = np.mean(zo_temp_pacmsk, axis=0)
                    zo_temp_indmsk_mean = np.mean(zo_temp_indmsk, axis=0)

            if "PO4" in globals():
                zo_PO4_mean        = np.mean(zo_PO4, axis=0)
                if "subbasin_mask" in globals():
                    zo_PO4_atlmsk_mean = np.mean(zo_PO4_atlmsk, axis=0)
                    zo_PO4_pacmsk_mean = np.mean(zo_PO4_pacmsk, axis=0)
                    zo_PO4_indmsk_mean = np.mean(zo_PO4_indmsk, axis=0)

            if "NO3" in globals():
                zo_NO3_mean        = np.mean(zo_NO3, axis=0)
                if "subbasin_mask" in globals():
                    zo_NO3_atlmsk_mean = np.mean(zo_NO3_atlmsk, axis=0)
                    zo_NO3_pacmsk_mean = np.mean(zo_NO3_pacmsk, axis=0)
                    zo_NO3_indmsk_mean = np.mean(zo_NO3_indmsk, axis=0)
            
            if "O2" in globals():
                zo_O2_mean        = np.mean(zo_O2,axis=0)
                if "subbasin_mask" in globals():
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

    if ((ind_file == 0 or ind_file == 1) and create_pdf_diff =='y') or create_pdf_diff =='n':
        # Plot parameters for standard simu (not diff)
        #-----------------------------------------------------------------------------#
        # SUBBASINS
        #----------------------------------#
        # from palettable.cartocolors.qualitative import Safe_3 # colormap import
        # cmapColor_subbasin = Safe_3.mpl_colormap
        cmapColor_subbasin = mpl.cm.get_cmap('viridis')
        figTitleSubbasin    = ''
        contour_subbasin    = None
        #----------------------------------#

        # BATHYMETRY 
        #----------------------------------#
        if "bathy" in globals():
    
            cmapColor_bathy  = 'YlGnBu'
            cbar_title_bathy = 'Bathymetry (m)'

            #  Automatic or manual colormap limits bathymetry
            if manual_lim_bathy == 'n': 
                cmap_bathy = mpl.cm.get_cmap(cmapColor_bathy)
                bounds     = np.linspace(bathy.min(),bathy.max(),50)
                norm_bathy = mpl.colors.BoundaryNorm(bounds, cmap_bathy.N)
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

            cbar_title_sst = 'SST (tos) (°C)'
            
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
        cmapColor_strf    = 'BrBG_r'

        if "zo_stream_function" in globals():
            
            cbar_title_zostrf = 'Zonal stream function (zomsfglo) (Sv = 10$^6$m$^3$.s$^{-1}$)'
            contour_zostrf    = np.array([-20, -10, 0, 10, 20])
            # Custom cmap
            bounds            = np.arange(-42,-15,4)
            bounds            = np.append(bounds,np.arange(-15,15,0.5))
            bounds            = np.append(bounds,np.arange(18,43,4))
            cmap_zostrf       = plt.get_cmap(cmapColor_strf,len(bounds)-1)
            norm_zostrf       = mpl.colors.BoundaryNorm(bounds, cmap_zostrf.N)
            
        if "vozocrtx" in globals():

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
            
            # North hemisphere
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
            
            # South hemisphere
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
            
            # Automatic or manual colormap limits PO4
            if manual_lim_po4 == 'n':
                cmap_po4, norm_po4 = adaptativeColorMap(PO4_mean, 0.015, 0.2, cmapColor_po4no3o2)
            elif manual_lim_po4 == 'y':
                cmap_po4 = mpl.cm.get_cmap(cmapColor_po4no3o2)
                bounds   = np.arange(min_lim_po4, max_lim_po4, step_po4)
                norm_po4 = mpl.colors.BoundaryNorm(bounds, cmap_po4.N)

        if "zo_PO4" in globals():
            
            cbar_title_zopo4 = 'Zonal Phosphate Concentration (PO4) (µmol.m$^{-3}$)'
                    
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
            
            # Automatic or manual colormap limits NO3
            if manual_lim_no3 == 'n':
                cmap_no3, norm_no3 = adaptativeColorMap(NO3_mean, 0.015, 0.2, cmapColor_po4no3o2)
            elif manual_lim_no3 == 'y':
                cmap_no3 = mpl.cm.get_cmap(cmapColor_po4no3o2)
                bounds   = np.arange(min_lim_no3,max_lim_no3,step_no3)
                norm_no3 = mpl.colors.BoundaryNorm(bounds, cmap_no3.N)

        if "zo_NO3" in globals():
            
            cbar_title_zono3 = 'Zonal Nitrate Concentration (NO3) (µmol.m$^{-3}$)'

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
            
            # Automatic or manual colormap limits O2
            if manual_lim_o2 == 'n':
                cmap_o2, norm_o2 = adaptativeColorMap(O2_mean, 0.015, 0.2, cmapColor_po4no3o2)
            elif manual_lim_o2 == 'y':
                cmap_o2 = mpl.cm.get_cmap(cmapColor_po4no3o2)
                bounds  = np.arange(min_lim_o2,max_lim_o2,step_o2)
                norm_o2 = mpl.colors.BoundaryNorm(bounds, cmap_o2.N)

        if "zo_O2" in globals():
            
            cbar_title_zoo2 = 'Zonal Oxygen Concentration (O2) (µmol.m$^{-3}$)'

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
    #-----------------------------------------------------------------------------#
    
    
    # Plot parameters for diff simu
    #-----------------------------------------------------------------------------#
    elif ind_file == 2 and create_pdf_diff =='y': # plot parameters for diff

        cBarTitleDiff = f'\n [{filePrefix[0]}] - [{filePrefix[1]}]'
        
        # SALINITY
        #----------------------------------#
        cmapColor_salinity = 'RdBu_r' #gnuplot2
    
        if "SSS" in globals(): 
    
            cbar_title_sss     = 'SSS (sos) (PSU)' + cBarTitleDiff

            cmap_sss = mpl.cm.get_cmap(cmapColor_salinity)
            lim_sss  = np.round(np.max(np.abs((np.max(SSS_mean),np.min(SSS_mean)))))
            bounds = np.linspace(-lim_sss,lim_sss,100)
            norm_sss = mpl.colors.BoundaryNorm(bounds, cmap_sss.N)
            
            # Contour
            intervContour = np.int32(len(norm_sss.boundaries)/nbContourLines)
            contour_sss   = norm_sss.boundaries[::intervContour]
    
        if "zo_salinity" in globals(): 
    
            cbar_title_zosalin         = 'Zonal salinity (PSU)' + cBarTitleDiff
            
            cmap_zosalin = mpl.cm.get_cmap(cmapColor_salinity)
            lim_zosalin  = np.round(np.max(np.abs((np.max(zo_salinity_mean),np.min(zo_salinity_mean)))),2)
            bounds = np.linspace(-lim_zosalin,lim_zosalin,100)
            norm_zosalin = mpl.colors.BoundaryNorm(bounds, cmap_zosalin.N)
            
            # Contour
            intervContour = np.int32(len(norm_zosalin.boundaries)/nbContourLines)
            contour_zosalin = norm_zosalin.boundaries[::intervContour]
        #----------------------------------#
    
        # TEMPERATURE
        #----------------------------------#
        cmapColor_temp           = 'RdBu_r' #gnuplot
    
        if "SST" in globals(): 
    
            cbar_title_sst     = 'SST (tos) (°C)' + cBarTitleDiff

            cmap_sst = mpl.cm.get_cmap(cmapColor_temp)
            lim_sst  = np.round(np.max(np.abs((np.max(SST_mean),np.min(SST_mean)))),2)
            bounds = np.linspace(-lim_sst,lim_sst,100)
            norm_sst = mpl.colors.BoundaryNorm(bounds, cmap_sst.N)
        
            # Contour
            intervContour = np.int32(len(norm_sst.boundaries)/nbContourLines)
            contour_sst   = norm_sst.boundaries[::intervContour]
    
        if "zo_temp" in globals(): 
    
            cbar_title_zotemp        = 'Zonal temperature (zotemglo) (°C)' + cBarTitleDiff
                
            cmap_zotemp = mpl.cm.get_cmap(cmapColor_temp)
            lim_zotemp  = np.round(np.max(np.abs((np.max(zo_temp_mean),np.min(zo_temp_mean)))),2)
            bounds = np.linspace(-lim_zotemp,lim_zotemp,100)
            norm_zotemp = mpl.colors.BoundaryNorm(bounds, cmap_zotemp.N)
        
            # Contour
            intervContour = np.int32(len(norm_zotemp.boundaries)/nbContourLines)
            contour_zotemp   = norm_zotemp.boundaries[::intervContour]
        #----------------------------------#
    
        # STREAM FUNCTION
        #----------------------------------#
        cmapColor_strf    = 'RdBu_r'

        if "zo_stream_function" in globals():

            cbar_title_zostrf = 'Zonal stream function (zomsfglo) (Sv = 10$^6$m$^3$.s$^{-1}$)' + cBarTitleDiff
            
            cmap_zostrf = mpl.cm.get_cmap(cmapColor_strf)
            lim_zostrf  = np.round(np.max(np.abs((np.max(zo_stream_function_mean),np.min(zo_stream_function_mean)))),2)
            bounds = np.linspace(-lim_zostrf,lim_zostrf,100)
            norm_zostrf = mpl.colors.BoundaryNorm(bounds, cmap_zostrf.N)

            # Contour
            intervContour = np.int32(len(norm_zostrf.boundaries)/nbContourLines)
            contour_zostrf   = norm_zostrf.boundaries[::intervContour]

        if "vozocrtx" in globals():

            cbar_title_bstrf  = 'Barotropic stream function (Sv = 10$^6$m$^3$.s$^{-1}$)' + cBarTitleDiff

            cmap_bstrf = mpl.cm.get_cmap(cmapColor_strf)
            lim_bstrf  = np.round(np.max(np.abs((np.max(baro_stream_function_mean),np.min(baro_stream_function_mean)))),2)
            bounds = np.linspace(-lim_bstrf,lim_bstrf,100)
            norm_bstrf = mpl.colors.BoundaryNorm(bounds, cmap_bstrf.N)

            # Contour
            intervContour = np.int32(len(norm_bstrf.boundaries)/nbContourLines)
            contour_bstrf   = norm_bstrf.boundaries[::intervContour]
        #----------------------------------#
    
        # OCEAN MIXED LAYER
        #----------------------------------#
        if "omlmax" in globals():
    
            cmapColor_oml    = 'RdBu_r'
            
            # omlnh
            #---------------#
            fig_title_omlnh  = 'January to March average'
            cbar_title_omlnh = 'Ocean mixed layer thickness (omlmax) (m) - Northern hemisphere' + cBarTitleDiff

            cmap_omlnh = mpl.cm.get_cmap(cmapColor_oml)
            lim_omlnh  = np.round(np.max(np.abs((np.max(omlmaxNH),np.min(omlmaxNH)))),2)
            bounds = np.linspace(-lim_omlnh,lim_omlnh,100)
            norm_omlnh = mpl.colors.BoundaryNorm(bounds, cmap_omlnh.N)

            # Contour
            intervContour = np.int32(len(norm_omlnh.boundaries)/nbContourLines)
            contour_omlnh   = norm_omlnh.boundaries[::intervContour]
            #---------------#

            # omlsh
            #---------------#
            fig_title_omlsh  = 'July to September average'
            cbar_title_omlsh = 'Ocean mixed layer thickness (omlmax) (m) - Southern hemisphere' + cBarTitleDiff

            cmap_omlsh = mpl.cm.get_cmap(cmapColor_oml)
            lim_omlsh  = np.round(np.max(np.abs((np.max(omlmaxSH),np.min(omlmaxSH)))),2)
            bounds = np.linspace(-lim_omlsh,lim_omlsh,100)
            norm_omlsh = mpl.colors.BoundaryNorm(bounds, cmap_omlsh.N)

            # Contour
            intervContour = np.int32(len(norm_omlsh.boundaries)/nbContourLines)
            contour_omlsh   = norm_omlsh.boundaries[::intervContour]
            #---------------#
        #----------------------------------#
    
        # EPC100 / TPP
        #----------------------------------#
        cmapColor_epctpp  = 'RdBu_r'
    
        if "TPP" in globals():
    
            cbar_title_intpp = 'Total Primary production of phyto depth integrated (INTPP) (g.m$^{-3}$.d$^{-1}$)' + cBarTitleDiff

            cmap_intpp = mpl.cm.get_cmap(cmapColor_epctpp)
            lim_intpp  = np.round(np.max(np.abs((np.max(INTPP),np.min(INTPP)))),2)
            bounds = np.linspace(-lim_intpp,lim_intpp,100)
            norm_intpp = mpl.colors.BoundaryNorm(bounds, cmap_intpp.N)

            # Contour
            # intervContour = np.int32(len(norm_intpp.boundaries)/nbContourLines)
            # contour_intpp   = norm_intpp.boundaries[::intervContour]
    
            contour_intpp = None
    
        if "EPC100" in globals():
            
            cbar_title_epc100 = 'Export of carbon particles at 100m (EPC100) (g.m$^{-2}$.d$^{-1}$)' + cBarTitleDiff # (mol.m$^{-2}$.s$^{-1}$) 
            
            cmap_epc100 = mpl.cm.get_cmap(cmapColor_epctpp)
            lim_epc100  = np.round(np.max(np.abs((np.max(EPC100),np.min(EPC100)))),2)
            bounds = np.linspace(-lim_epc100,lim_epc100,100)
            norm_epc100 = mpl.colors.BoundaryNorm(bounds, cmap_epc100.N)

            # Contour
            # intervContour = np.int32(len(norm_epc100.boundaries)/nbContourLines)
            # contour_epc100   = norm_epc100.boundaries[::intervContour]
            
            contour_epc100 = None
        #----------------------------------#
    
        # PO4 / NO3 / O2
        #----------------------------------#
        cmapColor_po4no3o2 = 'RdBu_r'
        
        if "zo_PO4" in globals():
            
            cbar_title_zopo4 = 'Zonal Phosphate Concentration (PO4) (µmol.m$^{-3}$)' + cBarTitleDiff

            cmap_zopo4 = mpl.cm.get_cmap(cmapColor_po4no3o2)
            lim_zopo4  = np.round(np.max(np.abs((np.max(zo_PO4_mean),np.min(zo_PO4_mean)))),2)
            bounds = np.linspace(-lim_zopo4,lim_zopo4,100)
            norm_zopo4 = mpl.colors.BoundaryNorm(bounds, cmap_zopo4.N)

            # Contour
            intervContour = np.int32(len(norm_zopo4.boundaries)/nbContourLines)
            contour_zopo4   = norm_zopo4.boundaries[::intervContour]

                
        if "zo_NO3" in globals():
            
            cbar_title_zono3 = 'Zonal Nitrate Concentration (NO3) (µmol.m$^{-3}$)' + cBarTitleDiff
            
            cmap_zono3 = mpl.cm.get_cmap(cmapColor_po4no3o2)
            lim_zono3  = np.round(np.max(np.abs((np.max(zo_NO3_mean),np.min(zo_NO3_mean)))),2)
            bounds = np.linspace(-lim_zono3,lim_zono3,100)
            norm_zono3 = mpl.colors.BoundaryNorm(bounds, cmap_zono3.N)

            # Contour
            intervContour = np.int32(len(norm_zono3.boundaries)/nbContourLines)
            contour_zono3   = norm_zono3.boundaries[::intervContour]
                
        if "zo_O2" in globals():
            
            cbar_title_zoo2 = 'Zonal Oxygen Concentration (O2) (µmol.m$^{-3}$)' + cBarTitleDiff
            
            cmap_zoo2 = mpl.cm.get_cmap(cmapColor_po4no3o2)
            lim_zoo2  = np.round(np.max(np.abs((np.max(zo_O2_mean),np.min(zo_O2_mean)))),2)
            bounds = np.linspace(-lim_zoo2,lim_zoo2,100)
            norm_zoo2 = mpl.colors.BoundaryNorm(bounds, cmap_zoo2.N)

            # Contour
            intervContour = np.int32(len(norm_zoo2.boundaries)/nbContourLines)
            contour_zoo2   = norm_zoo2.boundaries[::intervContour]
            #----------------------------------#
    #-----------------------------------------------------------------------------#

    #%%===========================================================================#
    #                              --< PLOTTING >--                               #
    #=============================================================================#
    # Variable initialization
    #----------------------------------#
    filecount  = 0
    savedfiles = []
    #----------------------------------#
    
    # SUBBASINS
    #-----------------------------------------------------------------------------#
    if "subbasin_mask" in globals():
    
        plotMapping(lon, lat, subbasin_mask, contour_subbasin, cbarTitleNone, land_mask, 
                projDataIn, projDataOut, cmapColor_subbasin, normNone, figTitleSubbasin,
                figXsize, figYsize, cbar_label_size, cbar_tick_size, title_font_size)
        
        # Legend
        color1 = cmapColor_subbasin(0)
        color2 = cmapColor_subbasin(0.5)
        color3 = cmapColor_subbasin(1.0)
        
        patch1 = mpatches.Patch(color=color1, label='Atlantic sub basin')
        patch2 = mpatches.Patch(color=color2, label='Pacific sub basin')
        patch3 = mpatches.Patch(color=color3, label='Indian sub basin')
        
        plt.legend(handles=[patch1, patch2, patch3],bbox_to_anchor=(0.65, 0.83), bbox_transform=plt.gcf().transFigure, prop={'size': plot_font_size})
    
        # Save figure
        filename = filePrefix[ind_file] + '_' + 'subbasin.' + fig_format                                               
        pathFilename = savedFigPath + filename
        plt.savefig(pathFilename, format=fig_format, bbox_inches='tight', facecolor='white')
        if create_pdf == 'y':
            filecount += 1; lnfile = 'file' + str(filecount) + '.' + fig_format
            os.system('ln -s ' + pathFilename + ' ' + profilePath + '/' + lnfile)
            savedfiles.append(filename)
    #-----------------------------------------------------------------------------#
    
    # BATHYMETRY 
    #-----------------------------------------------------------------------------#
    if "bathy" in globals():
        
        plotMapping(lon, lat, bathy, contour_bathy, cbar_title_bathy, land_mask, 
                projDataIn, projDataOut, cmap_bathy, norm_bathy, figTitleNone,
                figXsize, figYsize, cbar_label_size, cbar_tick_size, title_font_size)
        
        # Save figure
        filename = filePrefix[ind_file] + '_' + 'bathy.' + fig_format
        pathFilename = savedFigPath + filename
        plt.savefig(pathFilename, format=fig_format, bbox_inches='tight', facecolor='white')
        if create_pdf == 'y':
            filecount += 1; lnfile = 'file' + str(filecount) + '.' + fig_format
            os.system('ln -s ' + pathFilename + ' ' + profilePath + '/' + lnfile)
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
        filename = filePrefix[ind_file] + '_' + 'wind_speed_850.' + fig_format
        pathFilename = savedFigPath + filename
        plt.savefig(pathFilename, format=fig_format, bbox_inches='tight', facecolor='white')
        if create_pdf == 'y':
            filecount += 1; lnfile = 'file' + str(filecount) + '.' + fig_format
            os.system('ln -s ' + pathFilename + ' ' + profilePath + '/' + lnfile)
            savedfiles.append(filename)
    #-----------------------------------------------------------------------------#
    
    # SSS
    #-----------------------------------------------------------------------------#
    if "SSS" in globals():
        
        plotMapping(lon, lat, SSS_mean, contour_sss, cbar_title_sss, land_mask, 
                projDataIn, projDataOut, cmap_sss, norm_sss, figTitleYear,
                figXsize, figYsize, cbar_label_size, cbar_tick_size, title_font_size)
        
        # Save figure
        filename = filePrefix[ind_file] + '_' + 'SSS.' + fig_format                                               
        pathFilename = savedFigPath + filename
        plt.savefig(pathFilename, format=fig_format, bbox_inches='tight', facecolor='white')
        if create_pdf == 'y':
            filecount += 1; lnfile = 'file' + str(filecount) + '.' + fig_format
            os.system('ln -s ' + pathFilename + ' ' + profilePath + '/' + lnfile)
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
        
        filename = filePrefix[ind_file] + '_' + 'zoSalinity.' + fig_format
        pathFilename = savedFigPath + filename
        plt.savefig(pathFilename, format=fig_format, bbox_inches='tight', facecolor='white')
        if create_pdf == 'y':                                               
            filecount += 1; lnfile = 'file' + str(filecount) + '.' + fig_format
            os.system('ln -s ' + pathFilename + ' ' + profilePath + '/' + lnfile)
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
        
    
        filename = filePrefix[ind_file] + '_' + 'zoSalinity_subbasins.' + fig_format
        pathFilename = savedFigPath + filename
        plt.savefig(pathFilename, format=fig_format, bbox_inches='tight', facecolor='white')
        if create_pdf == 'y':                                               
            filecount += 1; lnfile = 'file' + str(filecount) + '.' + fig_format
            os.system('ln -s ' + pathFilename + ' ' + profilePath + '/' + lnfile)
            savedfiles.append(filename)
    #-----------------------------------------------------------------------------#
    
    # SST
    #-----------------------------------------------------------------------------#  
    if "SST" in globals():
        
        plotMapping(lon, lat, SST_mean, contour_sst, cbar_title_sst, land_mask, 
                projDataIn, projDataOut, cmap_sst, norm_sst, figTitleYear,
                figXsize, figYsize, cbar_label_size, cbar_tick_size, title_font_size)
        
        # Save figure
        filename = filePrefix[ind_file] + '_' + 'SST.' + fig_format                                           
        pathFilename = savedFigPath + filename
        plt.savefig(pathFilename, format=fig_format, bbox_inches='tight', facecolor='white')
        if create_pdf == 'y':    
            filecount += 1; lnfile = 'file' + str(filecount) + '.' + fig_format
            os.system('ln -s ' + pathFilename + ' ' + profilePath + '/' + lnfile)
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
        
        filename = filePrefix[ind_file] + '_' + 'zoTemp.' + fig_format                                             
        pathFilename = savedFigPath + filename
        plt.savefig(pathFilename, format=fig_format, bbox_inches='tight', facecolor='white')
        if create_pdf == 'y':  
            filecount += 1; lnfile = 'file' + str(filecount) + '.' + fig_format
            os.system('ln -s ' + pathFilename + ' ' + profilePath + '/' + lnfile)
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
        
        filename = filePrefix[ind_file] + '_' + 'zoTemp_subbasins.' + fig_format
        pathFilename = savedFigPath + filename
        plt.savefig(pathFilename, format=fig_format, bbox_inches='tight', facecolor='white')
        if create_pdf == 'y':                                               
            filecount += 1; lnfile = 'file' + str(filecount) + '.' + fig_format
            os.system('ln -s ' + pathFilename + ' ' + profilePath + '/' + lnfile)
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
        
        map1 = plt.step(zo_lat, zo_temp_mean[0,:], '-k', where='post')
        
        # Titles/Labels and size
        plt.title ('Zonal SST - Yearly average' , fontsize = title_font_size)
        plt.xlabel('Latitude (°N)'  , fontsize = xy_label_font_size)
        plt.ylabel('Temperature (°C)'      , fontsize = xy_label_font_size)
        plt.xticks(major_xticks , fontsize = xy_ticks_font_size)
        plt.yticks(major_yticks , fontsize = xy_ticks_font_size)
        # plt.ylim(c_lev_zosst[0],c_lev_zosst[1])
                   
        filename = filePrefix[ind_file] + '_' + 'zoSST.' + fig_format                                              
        pathFilename = savedFigPath + filename
        plt.savefig(pathFilename, format=fig_format, bbox_inches='tight', facecolor='white')
        if create_pdf == 'y': 
            filecount += 1; lnfile = 'file' + str(filecount) + '.' + fig_format
            os.system('ln -s ' + pathFilename + ' ' + profilePath + '/' + lnfile)
            savedfiles.append(filename)
    #-----------------------------------------------------------------------------#
    
    # ZONAL AVERAGE STREAM FUNCTION
    #-----------------------------------------------------------------------------#
    if "zo_stream_function" in globals():
    
        plotZonalAve(latGrid, -1*depthGrid, zo_stream_function_mean, contour_zostrf, cbar_title_zostrf, 
                     cmap_zostrf, norm_zostrf, figTitleYear, figXsize, figYsize, cbar_label_size, cbar_tick_size, 
                     title_font_size, xy_label_font_size, xy_ticks_font_size)
        
        filename = filePrefix[ind_file] + '_' + 'zoStreamFunc.' + fig_format
        pathFilename = savedFigPath + filename
        plt.savefig(pathFilename, format=fig_format, bbox_inches='tight', facecolor='white')
        if create_pdf == 'y':                                               
            filecount += 1; lnfile = 'file' + str(filecount) + '.' + fig_format
            os.system('ln -s ' + pathFilename + ' ' + profilePath + '/' + lnfile)
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
        filename = filePrefix[ind_file] + '_' + 'baroStreamFunc.' + fig_format                                            
        pathFilename = savedFigPath + filename
        plt.savefig(pathFilename, format=fig_format, bbox_inches='tight', facecolor='white')
        if create_pdf == 'y':
            filecount += 1; lnfile = 'file' + str(filecount) + '.' + fig_format
            os.system('ln -s ' + pathFilename + ' ' + profilePath + '/' + lnfile)
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
        filename = filePrefix[ind_file] + '_' + 'omlmaxNH.' + fig_format
        pathFilename = savedFigPath + filename
        plt.savefig(pathFilename, format=fig_format, bbox_inches='tight', facecolor='white')
        if create_pdf == 'y':
            filecount += 1; lnfile = 'file' + str(filecount) + '.' + fig_format
            os.system('ln -s ' + pathFilename + ' ' + profilePath + '/' + lnfile)
            savedfiles.append(filename)
    
        #-----------------------< Southern Hemisphere >-----------------------#    
        lon_lim_sh = [-180,180]
        lat_lim_sh = [-90,-45]
        
        plotMappingZoom(lon, lat, omlmaxSH, contour_omlsh, cbar_title_omlsh, 
                        lon_lim_sh, lat_lim_sh, land_mask, 
                        projDataIn, ccrs.AzimuthalEquidistant(central_longitude=0,central_latitude=-90), 
                        cmap_omlsh, norm_omlsh, fig_title_omlsh, figXsize, figYsize, cbar_label_size, cbar_tick_size, title_font_size)
        
        # Save figure
        filename = filePrefix[ind_file] + '_' + 'omlmaxSH.' + fig_format
        pathFilename = savedFigPath + filename
        plt.savefig(pathFilename, format=fig_format, bbox_inches='tight', facecolor='white')
        if create_pdf == 'y':
            filecount += 1; lnfile = 'file' + str(filecount) + '.' + fig_format
            os.system('ln -s ' + pathFilename + ' ' + profilePath + '/' + lnfile)
            savedfiles.append(filename)
    #-----------------------------------------------------------------------------#
    
    # INTPP
    #-----------------------------------------------------------------------------#
    if "TPP" in globals():
        
        plotMapping(lon, lat, INTPP_mean, contour_intpp, cbar_title_intpp, land_mask, 
                    projDataIn, projDataOut, cmap_intpp, norm_intpp, figTitleYear,
                    figXsize, figYsize, cbar_label_size, cbar_tick_size, title_font_size)
        
        # Save figure
        filename = filePrefix[ind_file] + '_' + 'intpp.' + fig_format
        pathFilename = savedFigPath + filename
        plt.savefig(pathFilename, format=fig_format, bbox_inches='tight', facecolor='white')
        if create_pdf == 'y':
            filecount += 1; lnfile = 'file' + str(filecount) + '.' + fig_format
            os.system('ln -s ' + pathFilename + ' ' + profilePath + '/' + lnfile)
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
        filename = filePrefix[ind_file] + '_' + 'epc100.' + fig_format
        pathFilename = savedFigPath + filename
        plt.savefig(pathFilename, format=fig_format, bbox_inches='tight', facecolor='white')
        if create_pdf == 'y':
            filecount += 1; lnfile = 'file' + str(filecount) + '.' + fig_format
            os.system('ln -s ' + pathFilename + ' ' + profilePath + '/' + lnfile)
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
                    cmap_po4, norm_po4, cbar_label_size, cbar_tick_size)
        
        
        # Save figure
        filename = filePrefix[ind_file] + '_' + 'PO4_lev.' + fig_format
        pathFilename = savedFigPath + filename
        plt.savefig(pathFilename, format=fig_format, bbox_inches='tight', facecolor='white')
        if create_pdf == 'y':
            filecount += 1; lnfile = 'file' + str(filecount) + '.' + fig_format
            os.system('ln -s ' + pathFilename + ' ' + profilePath + '/' + lnfile)
            savedfiles.append(filename)
        #-------------------------------------------------#

    if "zo_PO4" in globals(): 
        # ZONAL AVERAGE PO4
        #-------------------------------------------------#
        plotZonalAve(latGrid, -1*depthGrid, zo_PO4_mean, contour_zopo4, cbar_title_zopo4, 
                    cmap_zopo4, norm_zopo4, figTitleYear, figXsize, figYsize, cbar_label_size, cbar_tick_size, 
                    title_font_size, xy_label_font_size, xy_ticks_font_size)
        
        filename = filePrefix[ind_file] + '_' + 'zoPO4.' + fig_format
        pathFilename = savedFigPath + filename
        plt.savefig(pathFilename, format=fig_format, bbox_inches='tight', facecolor='white')
        if create_pdf == 'y':
            filecount += 1; lnfile = 'file' + str(filecount) + '.' + fig_format
            os.system('ln -s ' + pathFilename + ' ' + profilePath + '/' + lnfile)
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
        if "zo_PO4_atlmsk" in globals() and "zo_PO4_pacmsk" in globals() and "zo_PO4_indmsk" in globals():

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
            
            filename = filePrefix[ind_file] + '_' + 'zoPO4_subbasins.' + fig_format
            pathFilename = savedFigPath + filename
            plt.savefig(pathFilename, format=fig_format, bbox_inches='tight', facecolor='white')
            if create_pdf == 'y':
                filecount += 1; lnfile = 'file' + str(filecount) + '.' + fig_format
                os.system('ln -s ' + pathFilename + ' ' + profilePath + '/' + lnfile)
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
        filename = filePrefix[ind_file] + '_' + 'NO3_lev.' + fig_format
        pathFilename = savedFigPath + filename
        plt.savefig(pathFilename, format=fig_format, bbox_inches='tight', facecolor='white')
        if create_pdf == 'y':
            filecount += 1; lnfile = 'file' + str(filecount) + '.' + fig_format
            os.system('ln -s ' + pathFilename + ' ' + profilePath + '/' + lnfile)
            savedfiles.append(filename)
        #-------------------------------------------------#
        
    if "zo_NO3" in globals():     
        # ZONAL AVERAGE NO3
        #-------------------------------------------------#
        plotZonalAve(latGrid, -1*depthGrid, zo_NO3_mean, contour_zono3, cbar_title_zono3, 
                    cmap_zono3, norm_zono3, figTitleYear, figXsize, figYsize, cbar_label_size, cbar_tick_size, 
                    title_font_size, xy_label_font_size, xy_ticks_font_size)
        
        filename = filePrefix[ind_file] + '_' + 'zoNO3.' + fig_format
        pathFilename = savedFigPath + filename
        plt.savefig(pathFilename, format=fig_format, bbox_inches='tight', facecolor='white')
        if create_pdf == 'y':
            filecount += 1; lnfile = 'file' + str(filecount) + '.' + fig_format
            os.system('ln -s ' + pathFilename + ' ' + profilePath + '/' + lnfile)
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
        if "zo_NO3_atlmsk" in globals() and "zo_NO3_pacmsk" in globals() and "zo_NO3_indmsk" in globals():

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
            
            filename = filePrefix[ind_file] + '_' + 'zoNO3_subbasins.' + fig_format
            pathFilename = savedFigPath + filename
            plt.savefig(pathFilename, format=fig_format, bbox_inches='tight', facecolor='white')
            if create_pdf == 'y':
                filecount += 1; lnfile = 'file' + str(filecount) + '.' + fig_format
                os.system('ln -s ' + pathFilename + ' ' + profilePath + '/' + lnfile)
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
        filename = filePrefix[ind_file] + '_' + 'O2_lev.' + fig_format
        pathFilename = savedFigPath + filename
        plt.savefig(pathFilename, format=fig_format, bbox_inches='tight', facecolor='white')
        if create_pdf == 'y':
            filecount += 1; lnfile = 'file' + str(filecount) + '.' + fig_format
            os.system('ln -s ' + pathFilename + ' ' + profilePath + '/' + lnfile)
            savedfiles.append(filename)
        #-------------------------------------------------#

    if "zo_O2" in globals():        
        # ZONAL AVERAGE O2
        #-------------------------------------------------#
        plotZonalAve(latGrid, -1*depthGrid, zo_O2_mean, contour_zoo2, cbar_title_zoo2, 
                    cmap_zoo2, norm_zoo2, figTitleYear, figXsize, figYsize, cbar_label_size, cbar_tick_size, 
                    title_font_size, xy_label_font_size, xy_ticks_font_size)
        
        cont1 = plt.contour(latGrid, -1*depthGrid, zo_O2_mean, np.array([6.5, 65]), colors='w', linewidths=0.7) # Add dysoxia and anoxia contour lines
        plt.clabel(cont1, fmt=' {:.1f} '.format, fontsize='x-large')                                       # Add labels over contour lines
                
        filename = filePrefix[ind_file] + '_' + 'zoO2.' + fig_format
        pathFilename = savedFigPath + filename
        plt.savefig(pathFilename, format=fig_format, bbox_inches='tight', facecolor='white')
        if create_pdf == 'y':
            filecount += 1; lnfile = 'file' + str(filecount) + '.' + fig_format
            os.system('ln -s ' + pathFilename + ' ' + profilePath + '/' + lnfile)
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
        if "zo_O2_atlmsk" in globals() and "zo_O2_pacmsk" in globals() and "zo_O2_indmsk" in globals():

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
            
            filename = filePrefix[ind_file] + '_' + 'zoO2_subbasins.' + fig_format
            pathFilename = savedFigPath + filename
            plt.savefig(pathFilename, format=fig_format, bbox_inches='tight', facecolor='white')
            if create_pdf == 'y':
                filecount += 1; lnfile = 'file' + str(filecount) + '.' + fig_format
                os.system('ln -s ' + pathFilename + ' ' + profilePath + '/' + lnfile)
                savedfiles.append(filename)
        #-------------------------------------------------#
    #-----------------------------------------------------------------------------#

    #%%===========================================================================#
    #                            --< PDF GENERATION >--                           #
    #=============================================================================#

    if create_pdf == 'y':

        dopdf(savedPdfPath, savedfiles, f'userData/{profileName}/file', filecount, filePrefix[ind_file], sentence_to_use[ind_file])

    #%%===========================================================================#
    #                            --< DIFF COMPUTATION >--                         #
    #=============================================================================#

    if create_pdf_diff == 'y':

        if ind_file == 0:

            # STORE VARIABLES FROM SIMULATION 1
            #----------------------------------------------------------------------#
            # Store variables from simu 1 for diff calculation that will occur 
            # on next loop (when ind_file = 1)
            # And remove other useless variables
            lon_S1       = lon
            lat_S1       = lat
            deptht_S1    = deptht
            zo_lat_S1    = zo_lat
            latGrid_S1   = latGrid
            depthGrid_S1 = depthGrid

            if "subbasin_mask" in globals():
                del subbasin_mask
                del atlmsk, pacmsk, indmsk

            if "bathy" in globals():
                del bathy

            if "tmask" in globals():
                land_mask_S1 = land_mask
                del tmask, land_mask

            if "u850" in globals() and "v850" in globals() and "z850" in globals() and "lonWind" in globals() and "latWind" in globals():
                del u850, v850, z850
                del lonWind, latWind

            if "SSS" in globals():
                SSS_S1 = SSS 
                del SSS, SSS_mean

            if "salinity" in globals():
                zo_salinity_S1 = zo_salinity
                del salinity
                del zo_salinity, zo_salinity_mean

            if "zo_salinity_atlmsk" in globals() and "zo_salinity_pacmsk" in globals() and "zo_salinity_indmsk" in globals():
                del zo_salinity_atlmsk, zo_salinity_atlmsk_mean
                del zo_salinity_pacmsk, zo_salinity_pacmsk_mean
                del zo_salinity_indmsk, zo_salinity_indmsk_mean

            if "SST" in globals():
                SST_S1 = SST
                del SST, SST_mean

            if "temp" in globals():
                zo_temp_S1 = zo_temp
                del temp
                del zo_temp, zo_temp_mean

            if "zo_temp_atlmsk" in globals() and "zo_temp_pacmsk" in globals() and "zo_temp_indmsk" in globals():
                del zo_temp_atlmsk, zo_temp_atlmsk_mean
                del zo_temp_pacmsk, zo_temp_pacmsk_mean
                del zo_temp_indmsk, zo_temp_indmsk_mean

            if "zo_temp" in globals():
                del zo_temp

            if "zo_stream_function" in globals():
                zo_stream_function_S1 = zo_stream_function
                del zo_stream_function, zo_stream_function_mean

            if "baro_stream_function" in globals():
                baro_stream_function_S1 = baro_stream_function
                del baro_stream_function, baro_stream_function_mean

            if "omlmax" in globals():
                omlmax_S1 = omlmax
                del omlmax
                del omlmaxNH, omlmaxSH

            if "TPP" in globals():
                INTPP_S1 = INTPP
                del TPP
                del INTPP, INTPP_mean

            if "EPC100" in globals():
                EPC100_S1 = EPC100
                del EPC100, EPC100_mean

            if "PO4" in globals():
                zo_PO4_S1 = zo_PO4
                del PO4, PO4_mean
                del zo_PO4, zo_PO4_mean
                if "zo_PO4_atlmsk" in globals() and "zo_PO4_pacmsk" in globals() and "zo_PO4_indmsk" in globals():
                    del zo_PO4_atlmsk, zo_PO4_atlmsk_mean
                    del zo_PO4_pacmsk, zo_PO4_pacmsk_mean
                    del zo_PO4_indmsk, zo_PO4_indmsk_mean

            if "NO3" in globals():
                zo_NO3_S1 = zo_NO3
                del NO3, NO3_mean
                del zo_NO3, zo_NO3_mean
                if "zo_NO3_atlmsk" in globals() and "zo_NO3_pacmsk" in globals() and "zo_NO3_indmsk" in globals():
                    del zo_NO3_atlmsk, zo_NO3_atlmsk_mean
                    del zo_NO3_pacmsk, zo_NO3_pacmsk_mean
                    del zo_NO3_indmsk, zo_NO3_indmsk_mean

            if "O2" in globals():
                zo_O2_S1 = zo_O2
                del O2, O2_mean
                del zo_O2, zo_O2_mean
                if "zo_O2_atlmsk" in globals() and "zo_O2_pacmsk" in globals() and "zo_O2_indmsk" in globals():
                    del zo_O2_atlmsk, zo_O2_atlmsk_mean
                    del zo_O2_pacmsk, zo_O2_pacmsk_mean
                    del zo_O2_indmsk, zo_O2_indmsk_mean
            #----------------------------------------------------------------------#

        if ind_file == 1:

            # INTERPOLATE VARIABLES & COMPUTE DIFF
            #----------------------------------------------------------------------#
            # Check if interpolation is needed <=> grid 1 different from grid 2
            comparison_lon = lon_S1 == lon      ; comparison_lat = lat_S1 == lat
            equal_arrays_lon = comparison_lon.all() ; equal_arrays_lat = comparison_lat.all()

            if ~equal_arrays_lon or ~equal_arrays_lat: # if interpolation is needed

                # Interpolate variables from simu 2 over grid from simu 1
                land_mask_S2_interp            = myCustomInterp3(lon, lat, land_mask           , lon_S1, lat_S1)

                SSS_S2_interp                  = myCustomInterp3(lon, lat, SSS                 , lon_S1, lat_S1)
                SST_S2_interp                  = myCustomInterp3(lon, lat, SST                 , lon_S1, lat_S1)
                baro_stream_function_S2_interp = myCustomInterp3(lon, lat, baro_stream_function, lon_S1, lat_S1)
                omlmax_S2_interp               = myCustomInterp3(lon, lat, omlmax              , lon_S1, lat_S1)
                INTPP_S2_interp                = myCustomInterp3(lon, lat, INTPP               , lon_S1, lat_S1)
                EPC100_S2_interp               = myCustomInterp3(lon, lat, EPC100              , lon_S1, lat_S1)
    
                zo_salinity_S2_interp        = myCustomInterp3(latGrid, depthGrid, zo_salinity       , latGrid_S1, depthGrid_S1)
                zo_temp_S2_interp            = myCustomInterp3(latGrid, depthGrid, zo_temp           , latGrid_S1, depthGrid_S1)
                zo_stream_function_S2_interp = myCustomInterp3(latGrid, depthGrid, zo_stream_function, latGrid_S1, depthGrid_S1)
                zo_PO4_S2_interp             = myCustomInterp3(latGrid, depthGrid, zo_PO4            , latGrid_S1, depthGrid_S1)
                zo_NO3_S2_interp             = myCustomInterp3(latGrid, depthGrid, zo_NO3            , latGrid_S1, depthGrid_S1)
                zo_O2_S2_interp              = myCustomInterp3(latGrid, depthGrid, zo_O2             , latGrid_S1, depthGrid_S1)

                # Make difference between simu 1 & 2
                land_mask_diff            = land_mask_S1            + land_mask_S2_interp
                SSS_diff                  = SSS_S1                  - SSS_S2_interp
                SST_diff                  = SST_S1                  - SST_S2_interp
                baro_stream_function_diff = baro_stream_function_S1 - baro_stream_function_S2_interp
                omlmax_diff               = omlmax_S1               - omlmax_S2_interp
                INTPP_diff                = INTPP_S1                - INTPP_S2_interp
                EPC100_diff               = EPC100_S1               - EPC100_S2_interp
                zo_salinity_diff          = zo_salinity_S1          - zo_salinity_S2_interp
                zo_temp_diff              = zo_temp_S1              - zo_temp_S2_interp
                zo_stream_function_diff   = zo_stream_function_S1   - zo_stream_function_S2_interp
                zo_PO4_diff               = zo_PO4_S1               - zo_PO4_S2_interp
                zo_NO3_diff               = zo_NO3_S1               - zo_NO3_S2_interp
                zo_O2_diff                = zo_O2_S1                - zo_O2_S2_interp

            else: # if interpolation is not needed

                # Make difference between simu 1 & 2 straight away
                land_mask_diff            = land_mask_S1            + land_mask
                SSS_diff                  = SSS_S1                  - SSS
                SST_diff                  = SST_S1                  - SST
                baro_stream_function_diff = baro_stream_function_S1 - baro_stream_function
                omlmax_diff               = omlmax_S1               - omlmax
                INTPP_diff                = INTPP_S1                - INTPP
                EPC100_diff               = EPC100_S1               - EPC100
                zo_salinity_diff          = zo_salinity_S1          - zo_salinity
                zo_temp_diff              = zo_temp_S1              - zo_temp
                zo_stream_function_diff   = zo_stream_function_S1   - zo_stream_function
                zo_PO4_diff               = zo_PO4_S1               - zo_PO4
                zo_NO3_diff               = zo_NO3_S1               - zo_NO3
                zo_O2_diff                = zo_O2_S1                - zo_O2
            #----------------------------------------------------------------------#

            # VARIABLE REATTRIBUTION & COMPUTE MEAN
            #----------------------------------------------------------------------#
            # Make original coord variables equal to the ones of simu 1
            # => So we can plot diff at next loop using same variable names
            lon       = lon_S1
            lat       = lat_S1
            deptht    = deptht_S1   
            zo_lat    = zo_lat_S1   
            latGrid   = latGrid_S1  
            depthGrid = depthGrid_S1

            # Make original variables equal to variables_diff
            # => So we can plot diff at next loop using same variable names
            land_mask            = land_mask_diff
            SSS                  = SSS_diff
            SST                  = SST_diff
            baro_stream_function = baro_stream_function_diff
            omlmax               = omlmax_diff
            INTPP                = INTPP_diff
            EPC100               = EPC100_diff
            zo_salinity          = zo_salinity_diff
            zo_temp              = zo_temp_diff
            zo_stream_function   = zo_stream_function_diff
            zo_PO4               = zo_PO4_diff
            zo_NO3               = zo_NO3_diff
            zo_O2                = zo_O2_diff

            # Compute mean
            SSS_mean                  = np.mean(SSS, axis = 0)
            SST_mean                  = np.mean(SST, axis = 0)
            baro_stream_function_mean = np.mean(baro_stream_function, axis = 0)
            omlmaxNH                  = np.mean(omlmax[0:2,:,:],axis=0) # mixed layer thickness for winter in northern hemisphere (jan-mar)
            omlmaxSH                  = np.mean(omlmax[6:8,:,:],axis=0) # mixed layer thickness for winter in southern hemisphere (jul-sep)
            INTPP_mean                = np.mean(INTPP, axis = 0)
            EPC100_mean               = np.mean(EPC100, axis = 0)
            zo_salinity_mean          = np.mean(zo_salinity, axis = 0)
            zo_temp_mean              = np.mean(zo_temp, axis = 0)
            # zo_temp_mean              = np.ma.masked_where(zo_temp_mean==0, zo_temp_mean) # replace 0 by masked values
            zo_stream_function_mean   = np.mean(zo_stream_function, axis = 0)
            zo_PO4_mean               = np.mean(zo_PO4, axis = 0)
            zo_NO3_mean               = np.mean(zo_NO3, axis = 0)
            zo_O2_mean                = np.mean(zo_O2, axis = 0)
            #----------------------------------------------------------------------#

            # redefine last element of filePrefix so figures and pdf 
            # created for diff will not overwrite previous ones
            filePrefix.append(f'diff_[{filePrefix[0]}]-[{filePrefix[1]}]')

            if "u850" in globals() and "v850" in globals() and "z850" in globals() and "lonWind" in globals() and "latWind" in globals():
                del u850, v850, z850
                del lonWind, latWind
                
            if "PO4" in globals():
                del PO4, PO4_mean
            
            
            if "NO3" in globals():
                del NO3, NO3_mean
            
            
            if "O2" in globals():
                del O2, O2_mean
   

# %%
