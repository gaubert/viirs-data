'''
Created on Feb 9, 2011

@author: guillaume.aubert@eumetsat.int
'''

import h5py
import numpy
import eumetsat.common.geo_utils as geo

def load_array(ds):
    a = numpy.empty(shape=ds.shape, dtype=ds.dtype)
    a[:] = ds[:]
    return a


def distance_array(lat1_arr, lon1_arr, lat2_arr, lon2_arr):
    """ calculate the distance between each points """
    
    if lat1_arr.size != lat2_arr.size:
        raise Exception("Latitude arrays have different sizes")
    
    if lon1_arr.size != lon2_arr.size:
        raise Exception("Longitude arrays have different sizes")

    results = numpy.empty(lat1_arr.size)

    print("Compute distances for %d points\n" % (lat1_arr.size))
 
    for x in xrange(0, lat1_arr.size):
        results[x] = geo.distance(lat1_arr[x], lon1_arr[x], lat2_arr[x], lon2_arr[x])
        #print("Calculate distance (%f) between p1:(%f,%f) and p2:(%f,%f)\n" %(results[x], lat1_arr[x], lon1_arr[x], lat2_arr[x], lon2_arr[x]))
    
    return results




if __name__ == '__main__':
    
    new_lat_lsasaf_path      = '/tmp/N_LSASAF_LAT_M02_20110202193403.hdf5'
    orig_lat_lsasaf_path     = '/homespace/gaubert/Data/LSASAF/S-LSA_-HDF5_LSASAF_EPS-AVHR_LAT_M02_20110202193403'
    
    new_lon_lsasaf_path      = '/tmp/N_LSASAF_LON_M02_20110202193403.hdf5'
    orig_lon_lsasaf_path     = '/homespace/gaubert/Data/LSASAF/S-LSA_-HDF5_LSASAF_EPS-AVHR_LON_M02_20110202193403'
    
    orig_lat_f = h5py.File(orig_lat_lsasaf_path)
    orig_lon_f = h5py.File(orig_lon_lsasaf_path)
    
    new_lat_f = h5py.File(new_lat_lsasaf_path)
    new_lon_f = h5py.File(new_lon_lsasaf_path)
    
    # get datasets, flatten them and transform them back to floats
    lat1 = (load_array(orig_lat_f['LAT']).ravel())*0.0001
    lon1 = (load_array(orig_lon_f['LON']).ravel())*0.0001
    
    lat2 = (load_array(new_lat_f['LAT']).ravel())*0.001
    lon2 = (load_array(new_lon_f['LON']).ravel())*0.001
    
    result = distance_array(lat1, lon1, lat2, lon2)
    
    print("Compute min, max , average\n")
    
    print("min(distance) = %f , max(distance) = %f , avg(distance) = %f \n" %( numpy.amin(result), numpy.amax(result), numpy.average(result) ))
