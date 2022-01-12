# # -*- coding: utf-8 -*-
# #------------------------------------------------------------------------------#
# # Interpolate 3D variable (var_in) from grid_in to grid_out under format: 
# #   var_in  [L x M x N ]
# #   lon_in  [    M x N ]
# #   lat_in  [    M x N ]
# #   lon_out [    M'x N']
# #   lat_out [    M'x N']
# #
# # With Grid_in  = [ Z, lon_in , lat_in ]
# #      Grid_out = [ Z, lon_out, lat_out]
# #
# # Notice that Grid_in and Grid_out have the same 1st dimension: Z
# # (it can be time or depth for example)
# #
# # 2D interpolation will occur over 2nd and 3rd dimensions of the grid
# # and this will be done for each element of the 1st dimension of the grid 
# # to recreate 3D variable var_out
# #------------------------------------------------------------------------------#

def myCustomInterp3(lon_in,lat_in, var_in, lon_out, lat_out):

    # With scipy version on IRENE whene this was tested scipy.__version__ = 1.4.1
    # there is a strange behavior with masked array and scipy.interpolate.griddata
    # Once interpolated, depending on the input variable, masked values are either
    # converted into 0 or 1e+20 or other...
    # So before interpolation the masked array is converted into a numpy array with
    # masked values converted into NaNs. Afeter interpolation, the numpy array obtained
    # is converted back into a masked array.

    import numpy as np
    from scipy.interpolate import griddata

    coord_in  = np.array([lon_in.flatten(),lat_in.flatten()]).T
    coord_out = np.array([lon_out.flatten(),lat_out.flatten()]).T

    # var_out = np.ma.zeros(shape=(var_in.shape[0], lon_out.shape[0], lon_out.shape[1])) # Because of read above commentary
    var_out = np.zeros(shape=(var_in.shape[0], lon_out.shape[0], lon_out.shape[1]))

    # Conversion of masked array with masked values into np.array with nan is performed
    var_in = np.where(var_in.mask == True, np.nan, var_in)

    for ind in np.arange(0, var_in.shape[0]):

        var_interp = griddata(coord_in, var_in[ind,::].flatten(), coord_out, method='nearest')
        var_interp_reshaped   = var_interp.reshape(lon_out.shape)

        var_out[ind,::] = var_interp_reshaped

    # Conversion of np.array with nan into masked array with masked values
    var_out = np.ma.masked_where(np.isnan(var_out),var_out)

    return var_out
