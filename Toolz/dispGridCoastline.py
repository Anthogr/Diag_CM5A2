#!/usr/bin/env python
# zorder defined as an optional argument following https://linux.die.net/diveintopython/html/power_of_introspection/optional_arguments.html

# lon_data      = lon coordinate
# lat_data      = lat coordinate
# datacrs       = original grid projection for lon_data and lat_data (ex: ccrs.PlateCarree())
# landMask      = mask containing masked value over land cells and anything else (integer) on ocean cells
# linewidthdata = width on the line to be plotted

def dispGridCoastline(lon_data, lat_data, datacrs, landMask,linewidthdata, zorder=40):
    import numpy as np
    import matplotlib.pyplot as plt

    # land outline... the hard part
    # for each point of the cropped area, determine if it is a coastal point and plot (or not) accordingly
    land_outline_linewidth = linewidthdata
    for ilon in np.arange(0,landMask.shape[1]-1):
        for ilat in np.arange(0,landMask.shape[0]-1):
            # is there an ocean to the East or to the West?
            if (landMask.mask[ilat,ilon] != landMask.mask[ilat,ilon+1]):
                    lat1 = lat_data[ilat,ilon+1]
                    lat2 = lat_data[ilat+1,ilon+1]
                    lon1 = lon_data[ilat,ilon+1]
                    lon2 = lon_data[ilat+1,ilon+1,]
                    latpts = [lat1, lat2]; #print latpts
                    lonpts = [lon1, lon2]; #print lonpts
                    
                    # IF condition to avoid to plot horizontal line all the way from 180° up to -180° 
                    # for pieces of land spreading before and after 180° (=>split in two pieces in the plot,
                    # 1 part at the very east and the other at the very west) 
                    # (Remove the if condition to see the difference if you don't understand)
                    if (np.abs(lon2-lon1) < 50):
                        plt.plot(lonpts,latpts,'-',linewidth=land_outline_linewidth, color='k',zorder=zorder,transform=datacrs)
            # is there an ocean to the North or to the South?
            if (landMask.mask[ilat,ilon] != landMask.mask[ilat+1,ilon]):
                    lat1 = lat_data[ilat+1,ilon]
                    lat2 = lat_data[ilat+1,ilon+1,]
                    lon1 = lon_data[ilat+1,ilon]
                    lon2 = lon_data[ilat+1,ilon+1,]
                    latpts = [lat1, lat2]; #print latpts
                    lonpts = [lon1, lon2]; #print lonpts
                    
                    # IF condition to avoid to plot horizontal line all the way from 180° up to -180° 
                    # for pieces of land spreading before and after 180° (=>split in two pieces in the plot,
                    # 1 part at the very east and the other at the very west)
                    # (Remove the if condition to see the difference if you don't understand)
                    if (np.abs(lon2-lon1) < 50):
                        plt.plot(lonpts,latpts,'-',linewidth=land_outline_linewidth, color='k',zorder=zorder,transform=datacrs)
                        
    # # last point / nort pole
    # for ilon in np.arange(0,landMask.shape[1]-1):
    #     for ilat in np.arange(landMask.shape[0]-1, landMask.shape[0]):
    #         # is there an ocean to the East or to the West?
    #         if (landMask.mask[ilat,ilon] != landMask.mask[ilat,ilon+1]):
    #                 lat1 = lat_data[ilat,ilon+1]
    #                 lat2 = lat_data[ilat+1,ilon+1]
    #                 lon1 = lon_data[ilat,ilon+1]
    #                 lon2 = lon_data[ilat+1,ilon+1]
    #                 latpts = [lat1, lat2]; #print latpts
    #                 lonpts = [lon1, lon2]; #print lonpts
    #                 plt.plot(lonpts,latpts,'-',linewidth=land_outline_linewidth, color='k',zorder=zorder,transform=datacrs)
    
    # # deadling with the modulo
    # for ilat in np.arange(0,landMask.shape[0]-1):
    #     # is there an ocean to the East or to the West?
    #     if ((landMask.mask[ilat,0] == True) != (landMask.mask[ilat,-1])):
    #             lat1 = lat_data[ilat,0]
    #             lat2 = lat_data[ilat+1,0]
    #             lon1 = lon_data[ilat,0]
    #             lon2 = lon_data[ilat+1,0]
    #             latpts = [lat1, lat2]; #print latpts
    #             lonpts = [lon1, lon2]; #print lonpts
    #             plt.plot(lonpts,latpts,'-',linewidth=land_outline_linewidth, color='k',zorder=zorder,transform=datacrs)
    #     # is there an ocean to the North or to the South?
    #     if (landMask.mask[ilat,-1] != landMask.mask[ilat+1,-1]):
    #             lat1 = lat_data[ilat+1,-2]
    #             lat2 = lat_data[ilat+1,-1]
    #             lon1 = lon_data[ilat+1,-2]
    #             lon2 = lon_data[ilat+1,-1]
    #             latpts = [lat1, lat2]; #print latpts
    #             lonpts = [lon1, lon2]; #print lonpts
    #             plt.plot(lonpts,latpts,'-',linewidth=land_outline_linewidth, color='k',zorder=zorder,transform=datacrs)
