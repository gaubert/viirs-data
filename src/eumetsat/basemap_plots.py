'''
Created on Apr 27, 2011

@author: gaubert
'''
from netCDF4 import Dataset
import numpy 

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

def main():
    
    orig_file  = '/tmp/comp-tempo/20101206-EUR-L2P_GHRSST-SSTsubskin-AVHRR_METOP_A-eumetsat_sstmgr_metop02_20101206_000403-v01.7-fv01.0.nc'
    
    orig_dset = Dataset(orig_file, 'a')
    
    o_lat = orig_dset.variables['lat'][:].ravel()
    o_lon = orig_dset.variables['lon'][:].ravel()
    
    
    
    # lon_0 is the central longitude of the projection.
    # resolution = 'l' means use low resolution coastlines.
    # optional parameter 'satellite_height' may be used to
    # specify height of orbit above earth (default 35,786 km).
    m = Basemap(projection='geos',lon_0=-105,resolution='l')
    m.drawcoastlines()
    m.fillcontinents(color='coral',lake_color='aqua')
    # draw parallels and meridians.
    m.drawparallels(np.arange(-90.,120.,30.))
    m.drawmeridians(np.arange(0.,420.,60.))
    m.drawmapboundary(fill_color='aqua')
    
    x, y = m(o_lat[0:100] , o_lon[0:100])
    
    m.plot(x, y)
    
    plt.title("Full Disk Geostationary Projection")
    #plt.savefig('geos_full.png')
    plt.show()



if __name__ == '__main__':
    main()