#!/usr/bin/env python
# zorder defined as an optional argument following https://linux.die.net/diveintopython/html/power_of_introspection/optional_arguments.html

def disp_grid_coastline(lon_data, lat_data, datacrs_data, topo_bathy_file,linewidthdata, zorder=40):
    # land outline... the hard part
    # for each point of the cropped area, determine if it is a coastal point and plot (or not) accordingly
    land_outline_linewidth = linewidthdata
    landseamask = np.ma.masked_where(topo_bathy_file < 0,topo_bathy_file)
    for ilon in np.arange(0,topo_bathy_file.shape[1]-1):
        for ilat in np.arange(0,topo_bathy_file.shape[0]-1):
            # is there an ocean to the East or to the West?
            if (landseamask.mask[ilat,ilon] != landseamask.mask[ilat,ilon+1]):
                    lat1 = lat_data[ilat,ilon+1]
                    lat2 = lat_data[ilat+1,ilon+1]
                    lon1 = lon_data[ilat,ilon+1]
                    lon2 = lon_data[ilat+1,ilon+1,]
                    latpts = [lat1, lat2]; #print latpts
                    lonpts = [lon1, lon2]; #print lonpts
                    plt.plot(lonpts,latpts,'-',linewidth=land_outline_linewidth, color='k',zorder=zorder,transform=datacrs_data)
            # is there an ocean to the North or to the South?
            if (landseamask.mask[ilat,ilon] != landseamask.mask[ilat+1,ilon]):
                    lat1 = lat_data[ilat+1,ilon]
                    lat2 = lat_data[ilat+1,ilon+1,]
                    lon1 = lon_data[ilat+1,ilon]
                    lon2 = lon_data[ilat+1,ilon+1,]
                    latpts = [lat1, lat2]; #print latpts
                    lonpts = [lon1, lon2]; #print lonpts
                    
                    # IF condition to avoid to plot horizontal line all the way from 180° up to -180° 
                    # for pieces of land spreading before and after 180° (=>split in two pieces in the plot,
                    # 1 part at the very east and the other at the very west)  
                    if lon2-lon1 > 0:
                        plt.plot(lonpts,latpts,'-',linewidth=land_outline_linewidth, color='k',zorder=zorder,transform=datacrs_data)
                        
    # last point / nort pole
    # for ilon in np.arange(0,topo_bathy_file.shape[1]-1):
    #     for ilat in np.arange(topo_bathy_file.shape[0]-1, topo_bathy_file.shape[0]):
    #         # is there an ocean to the East or to the West?
    #         if (landseamask.mask[ilat,ilon] != landseamask.mask[ilat,ilon+1]):
    #                 lat1 = lat_data[ilat,ilon+1]
    #                 lat2 = lat_data[ilat+1,ilon+1]
    #                 lon1 = lon_data[ilat,ilon+1]
    #                 lon2 = lon_data[ilat+1,ilon+1]
    #                 latpts = [lat1, lat2]; #print latpts
    #                 lonpts = [lon1, lon2]; #print lonpts
    #                 plt.plot(lonpts,latpts,'-',linewidth=land_outline_linewidth, color='k',zorder=zorder,transform=datacrs_data)
    # deadling with the modulo
    for ilat in np.arange(0,topo_bathy_file.shape[0]-1):
        # is there an ocean to the East or to the West?
        if ((landseamask.mask[ilat,0] == True) != (landseamask.mask[ilat,-1])):
                lat1 = lat_data[ilat,0]
                lat2 = lat_data[ilat+1,0]
                lon1 = lon_data[ilat,0]
                lon2 = lon_data[ilat+1,0]
                latpts = [lat1, lat2]; #print latpts
                lonpts = [lon1, lon2]; #print lonpts
                plt.plot(lonpts,latpts,'-',linewidth=land_outline_linewidth, color='k',zorder=zorder,transform=datacrs_data)
        # is there an ocean to the North or to the South?
        if (landseamask.mask[ilat,-1] != landseamask.mask[ilat+1,-1]):
                lat1 = lat_data[ilat+1,-2]
                lat2 = lat_data[ilat+1,-1]
                lon1 = lon_data[ilat+1,-2]
                lon2 = lon_data[ilat+1,-1]
                latpts = [lat1, lat2]; #print latpts
                lonpts = [lon1, lon2]; #print lonpts
                plt.plot(lonpts,latpts,'-',linewidth=land_outline_linewidth, color='k',zorder=zorder,transform=datacrs_data)
