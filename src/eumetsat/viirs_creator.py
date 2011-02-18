'''
Created on Jan 27, 2011

@author: guillaume.aubert@eumetsat.int
'''

import os
import random
import numpy
import sys
import h5py

import eumetsat.common.num_utils as num_utils
import eumetsat.common.filesystem_utils as fs_utils

def extract_radiance(a_out, a_in):
    """ extract the radiance from the existing file and create a radiance and radiance_factor param """
    
    print("Processing %s \n" % (os.path.basename(a_in.filename)))
    
    # get the channel name from the filename. later user the channel number to access the radiance data in each of the files
    rad_name = os.path.basename(a_in.filename).split("_")[0][2:]
    
    if rad_name in ["M10", "M11", "M12", "M13", "M14", "M15", "M16"]:
        rad_num = rad_name[1:]
    else:
        rad_num = rad_name[-1]
    
    #get the channel number
    radiance_factors = None
    #print("get radiance a_in['All_Data']['VIIRS-M%s-SDR_All']['Radiance']" %(rad_num))
    radiance = a_in['All_Data']['VIIRS-M%s-SDR_All' % (rad_num)]["Radiance"]
    if a_in['All_Data']['VIIRS-M%s-SDR_All' % (rad_num) ].get("RadianceFactors"):
        radiance_factors = a_in['All_Data']['VIIRS-M%s-SDR_All' % (rad_num)]["RadianceFactors"]
    else:
        #print("No Radiance Factor for %s\n" %(os.path.basename(a_in.filename)))
        pass
    
    # quick test convert radiance to uint16 (don't care of the value or the moment)
    if rad_num in ('3','4','5','7','13'):
        a_out["Radiance_%s" % (rad_name)] = load_array(radiance).astype('uint16')
    else:
        a_out["Radiance_%s" % (rad_name)] = radiance[:]
    
    if radiance_factors:
        a_out["RadianceFactors_%s" % (rad_name)] = radiance_factors[:]
    
    return radiance, radiance_factors

def modified_extract_radiance(a_out, a_in, keep_rad, keep_fac):
    """ extract the variance from the existing file and create a radiance and radiance_factor param """
    
    print("Processing %s \n" % (os.path.basename(a_in.filename)))
    
    # get the channel name from the filename. later user the channel number to access the radiance data in each of the files
    rad_name = os.path.basename(a_in.filename).split("_")[0][2:]
    
    if rad_name in ["M10", "M11", "M12", "M13", "M14", "M15", "M16"]:
        rad_num = rad_name[1:]
    else:
        rad_num = rad_name[-1]
    
    #get the channel number
    radiance_factors = None
    radiance = a_in['All_Data']['VIIRS-M%s-SDR_All' % (rad_num)]["Radiance"]
    if a_in['All_Data']['VIIRS-M%s-SDR_All' % (rad_num) ].get("RadianceFactors"):
        radiance_factors = a_in['All_Data']['VIIRS-M%s-SDR_All' % (rad_num)]["RadianceFactors"]
        a_out["Radiance_%s" % (rad_name)] = radiance[:]
        a_out["RadianceFactors_%s" % (rad_name)] = radiance_factors[:]
    else:
        #print("No Radiance Factor for %s. Use Default Rad and Factors\n" %(os.path.basename(a_in.filename)))
        a_out["Radiance_%s" % (rad_name)] = keep_rad[:]
        a_out["RadianceFactors_%s" % (rad_name)] = keep_fac[:]
   
    return radiance, radiance_factors

def load_array(ds):
    a = numpy.empty(shape=ds.shape, dtype=ds.dtype)
    a[:] = ds[:]
    return a


def compute_variation(geo_info,x_dim,y_dim):
    """ compute variations """
    
    #load into a numpy array
    flatten_geo_info = load_array(geo_info)

    
    flatten_geo_info = flatten_geo_info.reshape(x_dim*y_dim)
    
    variations = []
    
    icpt = 0
    for indice in xrange(1,(x_dim*y_dim), 2):
        a = round(flatten_geo_info[indice-1], 4)
        b = round(flatten_geo_info[indice], 4)
        variations.append(a-b)
        #print("variation (lat[%d]-lat[%d]) :(%f-%f)=%f\n" % ((indice-1),indice,a,b,(a-b)))
        #icpt += 1
        #if icpt == 800:
        #   return
    
    print("min variation %f - max variation %f" % (min(variations),max(variations)))

def create_tie_points_grid_indexes():
    """ Create the indexes.
         I've analysed the data and found that compared to AVHRR/3 as a simple cross track scanner, the VIIRS data are more difficult to handle mainly for the following reasons:

            * The averaging across track results in discontinuity of pixel spacing between pixels #1008 and #1009, and #2192 and #2193, respectively.
            * Along track, we have blocks of 16 scan lines, and not a regular pattern.

            Along track: Scan lines 1,16,17,32,33,48,49,64,65, ... ,721,736,737,752,753,768 
            Across track: Pixels 1,17,33,49,65,81, ... ,977,993,1008,1009,1025,1041,1057, ... ,2145,2161,2177,2192,2193,2209,2225,2241,...3153,3169,3185,3200 
    
    """
    # create scan_line_index
    scan_lines_index = []
    
    i = 0
    while i < 768:
        i += 1
        scan_lines_index.append(i)
        i += 15
        scan_lines_index.append(i)
    
    print("scan_lines = %s\n" %(scan_lines_index))
    
    pixels_index     = []
    
    i = 1
    pixels_index.append(i)
    
    while i < 3200:
        i += 16
        if i in (1009,2193):
            pixels_index.append(i-1)
        elif i == 3201:
            i = 3200
       
        pixels_index.append(i)
         
    print("pixel index = %s\n" % (pixels_index) )
    
    return { "lines" : scan_lines_index, "pixels" : pixels_index}


def create_tie_points_grid(geo_file):
    """ create the recommended tie-points grid and extract lat lon """ 
    """ Height, SolarZenithAngle, SolarAzimuthAngle, SatelliteZenithAngle, SatelliteAzimuthAngle """
    
    
    # create index of points
    indexes = create_tie_points_grid_indexes()
    
    lines   = indexes['lines']
    pixels  = indexes['pixels']
    
    out_lat    = numpy.empty(( len(indexes['lines']), len(indexes['pixels']) ))
    out_lon    = numpy.empty(( len(indexes['lines']), len(indexes['pixels']) ))
    out_sol_za = numpy.empty(( len(indexes['lines']), len(indexes['pixels']) ))
    out_sol_aa = numpy.empty(( len(indexes['lines']), len(indexes['pixels']) ))
    out_sat_za = numpy.empty(( len(indexes['lines']), len(indexes['pixels']) ))
    out_sat_aa = numpy.empty(( len(indexes['lines']), len(indexes['pixels']) ))
    
    geo_out = {}
    extract_geo_spatial_data(geo_out, geo_file)
    
    in_lat    = geo_out['Latitude']
    in_lon    = geo_out['Longitude']
    in_sol_za = geo_out['SolarZenithAngle']
    in_sol_aa = geo_out['SolarAzimuthAngle']
    in_sat_za = geo_out['SatelliteZenithAngle']
    in_sat_aa = geo_out['SatelliteAzimuthAngle']
    
    out_sat_height = geo_out['Height']
    
    l = 0
    p = 0
    for l, line_i in enumerate(lines):
        for p, pixel_i in enumerate(pixels):
            out_lat[l][p]    = in_lat[line_i-1][pixel_i-1]
            #print("out_lat[%d][%d] = in_lat[%d][%d] = %f\n" %(l, p, line_i-1, pixel_i-1, in_lat[line_i-1][pixel_i-1]))
            out_lon[l][p]    = in_lon[line_i-1][pixel_i-1]
            #print("out_lon[%d][%d] = in_lon[%d][%d] = %f\n" %(l, p, line_i, pixel_i, in_lon[line_i][pixel_i]))
            out_sol_za[l][p] = in_sol_za[line_i-1][pixel_i-1]
            
            out_sol_aa[l][p] = in_sol_aa[line_i-1][pixel_i-1]
            
            out_sat_za[l][p] = in_sat_za[line_i-1][pixel_i-1]
            
            out_sat_aa[l][p] = in_sat_aa[line_i-1][pixel_i-1]
            
    return out_lat, out_lon, out_sol_za, out_sol_aa, out_sat_za, out_sat_aa, out_sat_height
            
    
def extract_geo_spatial_data(out, geo_file):
    """ get the Lat and Lon """
    
    lat = geo_file['All_Data']['VIIRS-MOD-GEO_All']['Latitude']
    
    #compute_variation(lat,768, 3200)
    
    lon = geo_file['All_Data']['VIIRS-MOD-GEO_All']['Longitude']
    
    #compute_variation(lon,768, 3200)
    
    out['Latitude']              = lat[:]
    out['Longitude']             = lon[:]
    out['Height']                = geo_file['All_Data']['VIIRS-MOD-GEO_All']['Height'][:]
    out['SolarZenithAngle']      = geo_file['All_Data']['VIIRS-MOD-GEO_All']['SolarZenithAngle'][:]
    out['SolarAzimuthAngle']     = geo_file['All_Data']['VIIRS-MOD-GEO_All']['SolarAzimuthAngle'][:]
    out['SatelliteZenithAngle']  = geo_file['All_Data']['VIIRS-MOD-GEO_All']['SatelliteZenithAngle'][:]
    out['SatelliteAzimuthAngle'] = geo_file['All_Data']['VIIRS-MOD-GEO_All']['SatelliteAzimuthAngle'][:]
    
    
def create_aggregated_viirs_dataset():
    
    dir  = "/homespace/gaubert/viirs/Mband-SDR"
    file = "SVM01_npp_d20030125_t0847056_e0848301_b00015_c20090513182937523620_gisf_pop.h5"
    geo_file     = h5py.File("%s/%s" %(dir,"GMODO_npp_d20030125_t0847056_e0848301_b00015_c20090513182937526121_gisf_pop.h5"))

    """
    #create a fake file with integer lat lon (on a 16 bits integer lat -90:90 and lon -180:180)
    
    #create type
    i16 = numpy.dtype('<i2')
    
    #create empty array
    
    lat_arr = numpy.random.random_integers(-9000,9000, (768,3200))
    lat_arr.astype('int16')

    lon_arr = numpy.random.random_integers(-18000,18000, (768,3200))
    lon_arr.astype('int16')

    
    o_file  =  h5py.File("/tmp/prototype_file.h5" ,"w")
    
    lat_dset = o_file.create_dataset("Latitude", data=lat_arr, dtype=i16)
    lon_dset = o_file.create_dataset("Longitude", data=lon_arr, dtype=i16)
    
    lat_dset = lat_arr[:]
    lon_dset = lon_arr[:]
    
    o_file.close()
    """
    
    #create file without geo spatial info
    """
    output_file  =  h5py.File("/tmp/aggreg_sin_geospatial_file.h5" ,"w")
    
    for file in fs_utils.dirwalk(dir,"SVM*.h5"):
        input_file   =  h5py.File(file ,"r")
        
        extract_radiance(output_file, input_file)
        
    output_file.close()
    """
    
    #create file 
    output_file  =  h5py.File("/tmp/aggreg_tie_points_file.h5" ,"w")
    #create type
    f32 = numpy.dtype('<f4')
    ui16 = numpy.dtype('<u2')
    i16 = numpy.dtype('<i2')
    
    #extract tie-points grid params
    (out_lat, out_lon, out_sol_za, out_sol_aa, out_sat_za, out_sat_aa, out_height) = create_tie_points_grid(geo_file)
    output_file.create_dataset('Latitude', data = out_lat.astype('float32'), dtype = f32)
    output_file.create_dataset('Longitude', data = out_lat.astype('float32'), dtype = f32)
    
   
    
    # convert height to i16
    numpy.around(out_height, 2, out_height)
    out_height = out_height * 100
    out_height = out_height.astype('int16')
    print("out_height min %d, max %d , range %d\n" % (numpy.min(out_height), numpy.max(out_height), (numpy.max(out_height)- numpy.min(out_height)) ))
    output_file.create_dataset('Height', data = out_height.astype('int16'), dtype = i16)
    
    
    # convert solar zenith angle  => Range 0-180
    # convert solar azimuth angle => Range -180-180
    numpy.around(out_sol_za, 3, out_sol_za)
    out_sol_za = out_sol_za * 1000
    out_sol_za = out_sol_za.astype('uint16')
    print("out_sol_za min %d, max %d , range %d\n" % (numpy.min(out_sol_za), numpy.max(out_sol_za), (numpy.max(out_sol_za)- numpy.min(out_sol_za)) ))
    
    numpy.around(out_sol_aa, 2, out_sol_aa)
    out_sol_aa = out_sol_aa * 100
    out_sol_aa = out_sol_aa.astype('int16')
    print("out_sol_aa min %d, max %d , range %d\n" % (numpy.min(out_sol_aa), numpy.max(out_sol_aa), (numpy.max(out_sol_aa)- numpy.min(out_sol_aa)) ))
    
    # convert sat zenith angle  => Range 0-180
    # convert sat azimuth angle => Range -180-180
    numpy.around(out_sat_za, 3, out_sat_za)
    out_sat_za = out_sat_za * 1000
    out_sat_za = out_sat_za.astype('uint16')
    print("out_sat_za min %d, max %d , range %d\n" % (numpy.min(out_sat_za), numpy.max(out_sat_za), (numpy.max(out_sat_za)- numpy.min(out_sat_za)) ))
    
    numpy.around(out_sat_aa, 2, out_sat_aa)
    out_sat_aa = out_sat_aa * 100
    out_sat_aa = out_sat_aa.astype('int16')
    print("out_sat_aa min %d, max %d , range %d\n" % (numpy.min(out_sat_aa), numpy.max(out_sat_aa), (numpy.max(out_sat_aa)- numpy.min(out_sat_aa)) ))
    
    output_file.create_dataset('SolarZenithAngle',      data = out_sol_za, dtype = ui16)
    output_file.create_dataset('SolarAzimuthAngle',     data = out_sol_aa, dtype = i16)
    output_file.create_dataset('SatelliteZenithAngle',  data = out_sat_za, dtype = ui16)
    output_file.create_dataset('SatelliteAzimuthAngle', data = out_sat_aa, dtype = i16)
    
    #output_file.create_dataset('SolarZenithAngle',      data = out_sol_za, dtype = f32)
    #output_file.create_dataset('SolarAzimuthAngle',     data = out_sol_aa, dtype = f32)
    #output_file.create_dataset('SatelliteZenithAngle',  data = out_sat_za, dtype = f32)
    #output_file.create_dataset('SatelliteAzimuthAngle', data = out_sat_aa, dtype = f32)
    
    for file in fs_utils.dirwalk(dir,"SVM*.h5"):
        
        input_file   =  h5py.File(file ,"r")
        
        extract_radiance(output_file, input_file)    
        
    output_file.close()

def print_radiances_as_binary():
    """ print the radiances a bin values """
    dir  = "/homespace/gaubert/viirs/Mband-SDR"
    file = "SVM03_npp_d20030125_t0847056_e0848301_b00015_c20090513182937524212_gisf_pop.h5"
    
    a_in = h5py.File("%s/%s" %(dir,file))

    # get the channel name from the filename. later user the channel number to access the radiance data in each of the files
    rad_name = os.path.basename(a_in.filename).split("_")[0][2:]
    
    if rad_name in ["M10", "M11", "M12", "M13", "M14", "M15", "M16"]:
        rad_num = rad_name[1:]
    else:
        rad_num = rad_name[-1]
    
    #get the channel number
    radiance_factors = None
    radiances = load_array(a_in['All_Data']['VIIRS-M%s-SDR_All' % (rad_num)]["Radiance"])
    
    flatten_radiances = radiances.reshape(768*3200)
    
    nb_good_values = 0
    
    for i, rad in enumerate(flatten_radiances):
        if not rad <= -999.0 :
            print("Rad[%d,%f] = %s\n" %(i, rad, num_utils.float_to_bin(rad)))
            print("ReducedRad[%d,%f] = %s\n" %(i, rad, num_utils.float_to_bin(rad)[9:]))
            nb_good_values += 1
        
        if nb_good_values == 300:
            sys.exit()
        
    

if __name__ == '__main__':
    #print("bin = %s\n" % (num_utils.float_to_bin(1.0887645)))
    #sys.exit()
    
    #dir  = "/homespace/gaubert/viirs/Mband-SDR"
    #geo_file     = h5py.File("%s/%s" %(dir,"GMODO_npp_d20030125_t0847056_e0848301_b00015_c20090513182937526121_gisf_pop.h5"))
    #create_tie_points_grid(geo_file)
    
    #print_radiances_as_binary()
    create_aggregated_viirs_dataset()
    
   

            
    
    
    
    