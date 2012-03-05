'''
Created on Apr 27, 2011

@author: gaubert
'''
from netCDF4 import Dataset
import numpy

import eumetsat.common.grid_utils as grid_utils

if __name__ == '__main__':
    
    orig_file  = '/tmp/comp-tempo/20101206-EUR-L2P_GHRSST-SSTsubskin-AVHRR_METOP_A-eumetsat_sstmgr_metop02_20101206_000403-v01.7-fv01.0.nc'
    modif_file = '/tmp/comp-tempo/20101206-EUR-L2P_GHRSST-SSTsubskin-AVHRR_METOP_A-eumetsat_sstmgr_metop02_20101206_000403-v01.7-fv01.0-new.nc'
    
    orig_file  = '/tmp/comp-tempo/20101206-EUR-L2P_GHRSST-SSTsubskin-AVHRR_METOP_A-eumetsat_sstmgr_metop02_20101206_002203-v01.7-fv01.0.nc'
    modif_file = '/tmp/comp-tempo/20101206-EUR-L2P_GHRSST-SSTsubskin-AVHRR_METOP_A-eumetsat_sstmgr_metop02_20101206_002203-v01.7-fv01.0-new.nc'
    
    orig_file  = '/tmp/comp-tempo/20101206-EUR-L2P_GHRSST-SSTsubskin-AVHRR_METOP_A-eumetsat_sstmgr_metop02_20101206_004003-v01.7-fv01.0.nc'
    modif_file = '/tmp/comp-tempo/20101206-EUR-L2P_GHRSST-SSTsubskin-AVHRR_METOP_A-eumetsat_sstmgr_metop02_20101206_004003-v01.7-fv01.0-new.nc'
    
    
    orig_dset = Dataset(orig_file, 'a')
    new_dset  = Dataset(modif_file, 'a')
    
    o_lat = orig_dset.variables['lat'][:].ravel()
    o_lon = orig_dset.variables['lon'][:].ravel()
    
    n_lat = new_dset.variables['lat'][:].ravel()
    n_lon = new_dset.variables['lon'][:].ravel()
    
    distances = grid_utils.distance_array(o_lat, o_lon, n_lat, n_lon)
    
    print("Compute min, max , average\n")
    
    print("min(distance) = %f , max(distance) = %f , avg(distance) = %f \n" %( numpy.amin(distances), numpy.amax(distances), numpy.average(distances) ))
    
    
    