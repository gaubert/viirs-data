'''
Created on Apr 27, 2011

@author: gaubert
'''
from netCDF4 import Dataset
import numpy as np
import h5py

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

def example_1():
    
    orig_file  = '/tmp/comp-tempo/20101206-EUR-L2P_GHRSST-SSTsubskin-AVHRR_METOP_A-eumetsat_sstmgr_metop02_20101206_000403-v01.7-fv01.0.nc'
    
    orig_dset = Dataset(orig_file, 'a')
    
    o_lat = orig_dset.variables['lat'][:].ravel()
    o_lon = orig_dset.variables['lon'][:].ravel()
    
    print(np.mean(o_lon))
    
    
    
    # lon_0 is the central longitude of the projection.
    # resolution = 'l' means use low resolution coastlines.
    # optional parameter 'satellite_height' may be used to
    # specify height of orbit above earth (default 35,786 km).
    m = Basemap(projection='geos',lon_0=133,resolution='l')
    m.drawcoastlines()
    m.fillcontinents(color='coral',lake_color='aqua')
    # draw parallels and meridians.
    m.drawparallels(np.arange(-90.,120.,30.))
    m.drawmeridians(np.arange(0.,420.,60.))
    m.drawmapboundary(fill_color='aqua')
    
    x, y = m(o_lat[0:100] , o_lon[0:100])
    
    #m.plot(x, y)
    
    plt.title("Full Disk Geostationary Projection")
    #plt.savefig('geos_full.png')
    plt.show()
    
    
def example_2():
    
    # retrieve data
    dir  = "/homespace/gaubert/viirs/real-data"
    #geo_file     = h5py.File("%s/%s" %(dir,"GMODO_npp_d20120224_t1821456_e1823098_b01694_c20120306214722023793_cspp_dev.h5"))
    #geo_file     = h5py.File("%s/%s" %(dir,"GMODO_npp_d20120224_t1823110_e1824352_b01694_c20120306215057805679_cspp_dev.h5"))
    #geo_file     = h5py.File("%s/%s" %(dir,"GMODO_npp_d20120224_t1826018_e1827260_b01694_c20120306215454419051_cspp_dev.h5"))
    #geo_file     = h5py.File("%s/%s" %(dir,"GMODO_npp_d20120224_t1827272_e1828514_b01694_c20120306215852012949_cspp_dev.h5"))
    geo_file      = h5py.File("/homespace/gaubert/GMTCO_npp_d20120224_t1100479_e1102121_b01689_c20120224172231282331_noaa_ops.h5")
    
    
    #lats = geo_file['All_Data']['VIIRS-MOD-GEO_All']['Latitude'][:]
    #lons = geo_file['All_Data']['VIIRS-MOD-GEO_All']['Longitude'][:]
    lats = geo_file['All_Data']['VIIRS-MOD-GEO-TC_All']['Latitude'][:]
    lons = geo_file['All_Data']['VIIRS-MOD-GEO-TC_All']['Longitude'][:]
    
    line_len = len(lats[0])
    col_len  = len(lats)
    
    print("line len %d, col len %d" % (line_len, col_len))
    print("Upper left corner point: (%f,%f)\n" % (lats[0][0], lons[0][0] ))
    print("Lower right corner point: (%f,%f)\n" % (lats[col_len-1][line_len-1], lons[col_len-1][line_len-1]))
    
    
    # draw map with markers for float locations
    #m = Basemap(projection='hammer',lon_0=180)
    lon_ref = lons[(col_len-1)/2][(line_len-1)/2]
    #lat_ref = 10
    lat_ref = lats[(col_len-1)/2][(line_len-1)/2]
    #m = Basemap(projection='ortho',lat_0=lat_ref,lon_0=lon_ref,resolution='l')
    m = Basemap(projection='nsper',lat_0=lat_ref,lon_0=lon_ref,satellite_height=2000*1000,resolution='l')
    
    #x, y = m(lons[0:10],lats[0:10])
    x,y  = m(lons,lats)
    
    
    m.drawcoastlines()
    
    m.drawmapboundary(fill_color='#99ffff')
    #m.fillcontinents(color='#cc9966',lake_color='#99ffff')
    #m.scatter(x,y,s = 1 ,color='k')
    
    m.drawgreatcircle(lons[0][0],lats[0][0],lons[0][-1],lats[0][-1],linewidth=1,color='b')
    m.drawgreatcircle(lons[0][0],lats[0][0],lons[col_len-1][0],lats[col_len-1][0],linewidth=1,color='b')
    m.drawgreatcircle(lons[col_len-1][0],lats[col_len-1][0],lons[col_len-1][line_len-1],lats[col_len-1][line_len-1],linewidth=1,color='b')
    m.drawgreatcircle(lons[0][line_len-1],lats[0][line_len-1],lons[col_len-1][line_len-1],lats[col_len-1][line_len-1],linewidth=1,color='b')
    
    plt.title('Location of VIIRS granule ')
    plt.savefig('/tmp/plot-gran3.png')



if __name__ == '__main__':
    example_2()