'''
Created on May 13, 2011

@author: gaubert
'''

from netCDF4 import Dataset
import numpy




if __name__ == '__main__':
    
    dir = '/homespace/gaubert/ifremer-data'
    
    input_files = [
#                  '20110502-EUR-L2P_GHRSST-SSTsubskin-AVHRR_METOP_A-eumetsat_sstmgr_metop02_20110502_220403-v01.7-fv01.0.nc',
                   '20110426-EUR-L2P_GHRSST-SSTsubskin-AVHRR_METOP_A-eumetsat_sstmgr_metop02_20110426_111003-v01.7-fv01.0.nc',
                   '20110420-EUR-L2P_GHRSST-SSTsubskin-AVHRR_METOP_A-eumetsat_sstmgr_metop02_20110420_064903-v01.7-fv01.0.nc',
                   '20110414-EUR-L2P_GHRSST-SSTsubskin-AVHRR_METOP_A-eumetsat_sstmgr_metop02_20110414_025203-v01.7-fv01.0.nc'
                 ]
    
    for input_file in input_files: 
    
        dataset = Dataset('%s/%s' % (dir,input_file),'a')
        
        lat = dataset.variables['lat']
        lon = dataset.variables['lon']
        
        lat_data = lat[:]
        lon_data = lon[:]
        
        lat_data = numpy.around(lat_data,3)
        lon_data = numpy.around(lon_data,3)
        
        dataset.variables['lat'][:]  = lat_data 
        dataset.variables['lon'][:]  = lon_data
        
        dataset.sync()
        
        dataset.close()