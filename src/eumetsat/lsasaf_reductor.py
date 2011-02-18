'''

Created on Feb 7, 2011

@author: guillaume.aubert@eumetsat.int
'''
import sys
import h5py
import numpy
import os
import eumetsat.common.utils
import eumetsat.common.geo_utils as geo


def load_array(ds):
    a = numpy.empty(shape=ds.shape, dtype=ds.dtype)
    a[:] = ds[:]
    return a

def create_lat_land_saf_file(path):
    """ create a lat land saf file """
    
    f = h5py.File(path,'w')
    
    #create 16 bits type for Lat
    dset = f.create_dataset("LAT", (1080 , 2048), 'i4')
    #add lsasaf attributes
    dset.attrs['CLASS']          = "Data"
    dset.attrs['PRODUCT']        = "LAT" 
    dset.attrs['PRODUCT_ID']     = 255 
    dset.attrs['N_COLS']         = 2048 
    dset.attrs['N_LINES']        = 1080
    dset.attrs['NB_BYTES']       = 4
    dset.attrs['SCALING_FACTOR'] = 10000.0
    dset.attrs['OFFSET']         = 0.0
    dset.attrs['MISSING_VALUE']  = -3276800
    dset.attrs['UNITS']          = "Deg."
    dset.attrs['CAL_SLOPE']      = 0.0
    dset.attrs['CAL_OFFSET']     = 0.0
    
    f.flush()
    
    return f
    
def create_lon_land_saf_file(path):
    """ create a lon land saf file """
    
    f = h5py.File(path,'w')
    
    #create 16 bits type for Lat
    dset = f.create_dataset("LON", (1080 , 2048), 'i4')
    #add lsasaf attributes
    dset.attrs['CLASS']          = "Data"
    dset.attrs['PRODUCT']        = "LON" 
    dset.attrs['PRODUCT_ID']     = 255 
    dset.attrs['N_COLS']         = 2048 
    dset.attrs['N_LINES']        = 1080
    dset.attrs['NB_BYTES']       = 4
    dset.attrs['SCALING_FACTOR'] = 10000.0
    dset.attrs['OFFSET']         = 0.0
    dset.attrs['MISSING_VALUE']  = -3276800
    dset.attrs['UNITS']          = "Deg."
    dset.attrs['CAL_SLOPE']      = 0.0
    dset.attrs['CAL_OFFSET']     = 0.0
    
    f.flush()

    return f



if __name__ == '__main__':
    
    new_lat_lsasaf_path      = '/tmp/N_LSASAF_LAT_M02_20110202193403.hdf5'
    orig_lat_lsasaf_path     = '/homespace/gaubert/Data/LSASAF/S-LSA_-HDF5_LSASAF_EPS-AVHR_LAT_M02_20110202193403'
    
    new_lon_lsasaf_path      = '/tmp/N_LSASAF_LON_M02_20110202193403.hdf5'
    orig_lon_lsasaf_path     = '/homespace/gaubert/Data/LSASAF/S-LSA_-HDF5_LSASAF_EPS-AVHR_LON_M02_20110202193403'
    
    new_lat_f = create_lat_land_saf_file(new_lat_lsasaf_path)
    new_lon_f = create_lon_land_saf_file(new_lon_lsasaf_path)

    orig_lat_f = h5py.File(orig_lat_lsasaf_path)
    orig_lon_f = h5py.File(orig_lon_lsasaf_path)
    
    orig_lat_dset= load_array(orig_lat_f['LAT'])
    orig_lon_dset= load_array(orig_lon_f['LON'])
    
    # do the work on lats
    # apply orig scale factor
    n_ds = orig_lat_dset*0.0001
    # round value
    numpy.around(n_ds, 3, n_ds)
    #apply new scale factor
    n_ds = n_ds*1000
    #n_ds = n_ds/0.003
    #numpy.around(n_ds,0, n_ds)
    # n_ds = n_ds*0.003*10000 # to get the non scaled values
    #write values
    new_lat_f['LAT'].write_direct(n_ds)
    #close lat files
    new_lat_f.close()
    orig_lat_f.close()
    
    # do the work on lons
    # apply orig scale factor
    n_ds = orig_lon_dset*0.0001
    # round value
    numpy.around(n_ds, 3, n_ds)
    #apply new scale factor
    n_ds = n_ds*1000
    
    #n_ds = n_ds/0.006
    # round value
    #numpy.around(n_ds, 0, n_ds)
    #n_ds = n_ds*0.006*10000 # for non scaled values
    #write values
    new_lon_f['LON'].write_direct(n_ds)
    #close lat files
    new_lon_f.close()
    orig_lon_f.close()
    
    
    