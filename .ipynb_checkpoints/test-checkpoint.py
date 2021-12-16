#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  7 11:32:29 2021

@author: anthog
"""

# %%
# conda activate oceano
import xarray as xr
import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

# MAP PROJECTIONS [https://scitools.org.uk/cartopy/docs/latest/crs/projections.html]
#projdata = ccrs.LambertCylindrical() # Andy's style
#projdata = ccrs.EqualEarth() # The trendy one
#projdata = ccrs.RotatedPole() # The funny one
projdata = ccrs.EckertIV() # Chris' style




fn = '/Users/anthog/Documents/CEREGE/DATA/CMA5A2/C20M-ICE-GA_SE_4755_4854_1M_grid_T.nc'

# ds = nc.Dataset(fn)

ds = xr.open_dataset(fn)

lon=ds.nav_lon
lat=ds.nav_lat


cf = plt.pcolormesh(lon,lat,ds.tos[5,:,:])
cf = plt.contourf(lon,lat,ds.tos[5,:,:])

ds.tos[5,:,:].plot()


fig = plt.figure()
ax = plt.subplot(projection="aitoff")
plt.title("PLOT")
plt.grid(True)

test = ax.contourf(lon,lat,ds.tos[1,:,:])

cbar = plt.colorbar(test,orientation='horizontal') 


#%%
import matplotlib.pyplot as plt

import cartopy.crs as ccrs
import iris
import iris.analysis.cartography
import iris.plot as iplt
import iris.quickplot as qplt


def main():
    # Load data
    filepath = iris.sample_data_path('orca2_votemper.nc')
    cube = iris.load_cube(filepath)

    # Choose plot projections
    projections = {}
    projections['Mollweide'] = ccrs.Mollweide()
    projections['PlateCarree'] = ccrs.PlateCarree()
    projections['NorthPolarStereo'] = ccrs.NorthPolarStereo()
    projections['Orthographic'] = ccrs.Orthographic(central_longitude=-90,
                                                    central_latitude=45)

    pcarree = projections['PlateCarree']
    # Transform cube to target projection
    new_cube, extent = iris.analysis.cartography.project(cube, pcarree,
                                                         nx=400, ny=200)

    # Plot data in each projection
    for name in sorted(projections):
        fig = plt.figure()
        fig.suptitle('ORCA2 Data Projected to {}'.format(name))
        # Set up axes and title
        ax = plt.subplot(projection=projections[name])
        # Set limits
        ax.set_global()
        # plot with Iris quickplot pcolormesh
        qplt.pcolormesh(new_cube)
        # Draw coastlines
        ax.coastlines()

        iplt.show()

if __name__ == '__main__':
    main()




#%% Aide Quentin check putNANinArray.py

import numpy as np

#-----------------------------------------------------------------------------
nc1 =  Dataset("/Users/anthog/Documents/CEREGE/DATA/CMA5A2/C20M-ICE-GA_SE_4755_4854_1M_grid_T.nc", "r", format="NETCDF4")

# grille latitude ORCA
lat_ocean = nc1.variables['nav_lat'][:][:]

# grille longitude ORCA    
lon_ocean = nc1.variables['nav_lon'][:][:]

# fonction pour lalculer la moyenne sur une saison (annuelle est comptee comme une saison)
var1 = moyenneSaison(nc1.variables[nom_var][:][:][:],saison)
#-----------------------------------------------------------------------------


import numpy as np

# Generate fake variables lon lat and sst
#------------------------------#
lon=np.array((
    [[-20,-19,-18,-17],
     [-20,-19,-18,-17],
     [-20,-19,-18,-17],
     [-20,-19,-18,-17]]))

lat=np.array((
    [[40,40,40,40],
     [41,41,41,41],
     [42,42,42,42],
     [43,43,43,43]]))

sst=np.array((
    [[2.3,3.1,5.6,2.4],
     [3.8,5.2,7.8,12.5],
     [14.6,15,8.3,2.5],
     [7.1,16.5,3.9,15.9]]))
#------------------------------#

# Find indices where for data located between -20<lon<-17 and 40<lat<43
ind=np.argwhere(((lon > -20) & (lon < -17) & (lat > 40) & (lat < 43)))

# Trick to transpose ind array so that the next line works 
ind = np.r_[ind].T

# Put nan in variable sst at indices identified above
sst[ind[0],ind[1]]=np.nan




aaa=var1[34:44,2:12].copy()

