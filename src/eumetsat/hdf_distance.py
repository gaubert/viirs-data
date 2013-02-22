'''
Created on Aug 7, 2012

@author: gaubert
'''

import numpy
import h5py
import eumetsat.common.grid_utils as grid_utils

if __name__ == '__main__':
    
    granules_name = ['d20120224_t1821456_e1823098_b01694', 'd20120224_t1823110_e1824352_b01694', 'd20120224_t1826018_e1827260_b01694', 'd20120224_t1827272_e1828514_b01694']

    orig_filename = "/homespace/gaubert/viirs/real-data/GMODO_npp_d20120224_t1821456_e1823098_b01694_c20120306214722023793_cspp_dev.h5"
    new_filename  = "/tmp/expanded/GMODO_npp_ears_d20120224_t1821456_e1823098_b01694.h5"
    
    #new_filename  = "/tmp/expanded/GMODO_npp_ears_d20120224_t1826018_e1827260_b01694.h5"
    #orig_filename = "/homespace/gaubert/viirs/real-data/GMODO_npp_d20120224_t1826018_e1827260_b01694_c20120306215454419051_cspp_dev.h5"
    
    ofile = h5py.File(orig_filename)
    nfile = h5py.File(new_filename)
    
    o_lat = ofile['All_Data']['VIIRS-MOD-GEO_All']['Latitude'][:].ravel()
    o_lon = ofile['All_Data']['VIIRS-MOD-GEO_All']['Longitude'][:].ravel()
    
    n_lat = nfile['All_Data']['VIIRS-MOD-GEO_All']['Latitude'][:].ravel()
    n_lon = nfile['All_Data']['VIIRS-MOD-GEO_All']['Longitude'][:].ravel()
    
    
    #distances = grid_utils.distance_array(o_lat[0:3000], o_lon[0:3000], n_lat[0:3000], n_lon[0:3000])
    distances = grid_utils.distance_array(o_lat, o_lon, n_lat, n_lon)
    
    print("Compute min, max , average\n")
    
    print("min(distance) = %f , max(distance) = %f , avg(distance) = %f \n" %( numpy.amin(distances), numpy.amax(distances), numpy.average(distances) ))
    
    
    
    