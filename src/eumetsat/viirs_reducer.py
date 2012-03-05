'''
Created on Mar 5, 2012

@author: guillaume.aubert@eumetsat.int
'''

import os
import numpy
import sys
import h5py
import glob

class TiePointGridCreatorBase(object):
    
    
    @classmethod
    def extract_geo_spatial_reducable_data(cls, geo_file):
        """ get the Lat and Lon """
        
        out = {}
    
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
        
        return out
    
    @classmethod
    def extract_additional_geoloc_info(cls, geo_file):
        """
           extract additional geoloc info
        """
        out = {}
        out['MidTime']               = geo_file['All_Data']['VIIRS-MOD-GEO_All']['MidTime'][:]
        out['ModeGran_Geo']          = geo_file['All_Data']['VIIRS-MOD-GEO_All']['ModeGran'][:]
        out['ModeScan_Geo']          = geo_file['All_Data']['VIIRS-MOD-GEO_All']['ModeScan'][:]
        out['NumberOfScans_Geo']     = geo_file['All_Data']['VIIRS-MOD-GEO_All']['NumberOfScans'][:]
        out['PadByte1_Geo']          = geo_file['All_Data']['VIIRS-MOD-GEO_All']['PadByte1'][:]
        out['QF1_SCAN_VIIRS_SDRGEO'] = geo_file['All_Data']['VIIRS-MOD-GEO_All']['QF1_SCAN_VIIRSSDRGEO'][:]
        out['QF2_VIIRSSDRGEO']       = geo_file['All_Data']['VIIRS-MOD-GEO_All']['QF2_VIIRSSDRGEO'][:]
        out['StartTime']             = geo_file['All_Data']['VIIRS-MOD-GEO_All']['StartTime'][:]
        
        return out
        
    
    @classmethod
    def create_tie_points_grid(cls, geo_file):
        """ create the recommended tie-points grid and extract lat lon """ 
        """ Height, SolarZenithAngle, SolarAzimuthAngle, SatelliteZenithAngle, SatelliteAzimuthAngle """
    
        h5_geo_file = h5py.File(geo_file)
    
        # create index of points
        indexes = cls.create_tie_points_grid_indexes()
        
        lines   = indexes['lines']
        pixels  = indexes['pixels']
        
        out_lat    = numpy.empty(( len(indexes['lines']), len(indexes['pixels']) ))
        out_lon    = numpy.empty(( len(indexes['lines']), len(indexes['pixels']) ))
        out_sol_za = numpy.empty(( len(indexes['lines']), len(indexes['pixels']) ))
        out_sol_aa = numpy.empty(( len(indexes['lines']), len(indexes['pixels']) ))
        out_sat_za = numpy.empty(( len(indexes['lines']), len(indexes['pixels']) ))
        out_sat_aa = numpy.empty(( len(indexes['lines']), len(indexes['pixels']) ))
        
        geo_out = cls.extract_geo_spatial_reducable_data(h5_geo_file)
        
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
                
                
        res = {'lat': out_lat, 'lon': out_lon, 'sol_za' : out_sol_za, 'sol_aa' : out_sol_aa, 
                 'sat_za' : out_sat_za, 'sat_aa': out_sat_aa, 'height': out_sat_height }
        
        res.update(cls.extract_additional_geoloc_info(h5_geo_file))
        
        return res
    

class TiePointZoneGridCreator(TiePointGridCreatorBase):
    """
       Functor creating a Tie-Point Zone Grid from a geolocation file
    """
    
    @classmethod
    def create_tie_points_grid_indexes(cls):
        """
           Create tie point zone grid index
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
        
        i = 0
        pixels_index.append(i)
        
        while i < 3200:
            i += 1    
            pixels_index.append(i)
            i += 15
            pixels_index.append(i)
             
        print("len(pixel_index) = %s, pixel index = %s\n" % (len(pixels_index), pixels_index) )
        
        print("len(scan_line_index) = %s, scan_line index = %s\n" % (len(scan_lines_index), scan_lines_index) )
        
        return { "lines" : scan_lines_index, "pixels" : pixels_index}



class ReduceTiePointGridCreator(TiePointGridCreatorBase):
    """
       Functor creating a Reduced Tie-Point Grid from a geolocation file
    """
    
    @classmethod
    def create_tie_points_grid_indexes(cls):
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
            if i in (641,1009,2193,2561):
                pixels_index.append(i-1)
            elif i == 3201:
                i = 3200
           
            pixels_index.append(i)
             
        print("len(pixel_index) = %s, pixel index = %s\n" % (len(pixels_index), pixels_index) )
        
        print("len(scan_line_index) = %s, scan_line index = %s\n" % (len(scan_lines_index), scan_lines_index) )
        return { "lines" : scan_lines_index, "pixels" : pixels_index}
    
    
    
    

class VIIRSReducer(object):
    
    GEO  = 'GMODO_npp'
    M1   = 'SVM01_npp'
    M2   = 'SVM02_npp'
    M3   = 'SVM03_npp'
    M4   = 'SVM04_npp'
    M5   = 'SVM05_npp'
    M6   = 'SVM06_npp'
    M7   = 'SVM07_npp'
    M8   = 'SVM08_npp'
    M9   = 'SVM09_npp'
    M10  = 'SVM10_npp'
    M11  = 'SVM11_npp'
    M12  = 'SVM12_npp'
    M13  = 'SVM13_npp'
    M14  = 'SVM14_npp'
    M15  = 'SVM15_npp'
    M16  = 'SVM16_npp'
    
    RADIANCE_PREFIX_LIST = [M1, M2, M3, M4, M5, M6, M7, M8, M9, M10, M11, M12, M13, M14, M15, M16]
    
    PREFIX_LIST = [GEO]
    PREFIX_LIST.extend(RADIANCE_PREFIX_LIST)
    
    TYPE_LIST  = { GEO : 'geo', M1 : 'm1', M2 : 'm2', M3 : 'm3', M4 : 'm4', M5 : 'm5', M6 : 'm6', \
                   M7: 'm7', M8 : 'm8', M9 : 'm9', M10 : 'm10', M11 : 'm11', M12 : 'm12', M13 : 'm13', \
                   M14: 'm14', M15 : 'm15', M16 : 'm16'}
    
    def __init__(self):
        """
           constructor
        """
        self.list_of_files = {}
        
    def locate_files_for_granule(self, a_dir, a_id):
        """
           Find all the files related to a given granule.
           
           Params
           a_id: identifier for the granule (used to recreate the file names
        """
        
        for prefix in self.PREFIX_LIST:
            filename = '%s/%s_%s_*.h5' % (a_dir, prefix, a_id)
            files = glob.glob(filename)
            if len(files) == 0:
                raise Exception("Cannot find %s file " % (filename)) 
            else:
                self.list_of_files[self.TYPE_LIST[prefix]] = files[0]
                
        return self.list_of_files
    
    def get_tie_point_zone_grid_info(self):
        """
           return tie point zone grid info
        """
        return TiePointZoneGridCreator.create_tie_points_grid(self.list_of_files['geo'])
    
    def get_reduced_tie_point_grid_info(self):
        """
           return the reduced tie-point grid information (lat, lon, solar_azimuth_angle, solar_zenith_angle, satellite_azimuth_angle, satellite_zenith_angle, 
        """
        return ReduceTiePointGridCreator.create_tie_points_grid(self.list_of_files['geo'])
    
    def create_aggregated_viirs_dataset(self, output_filename):
        """
        """
    
        dir  = "/homespace/gaubert/viirs/Mband-SDR"
        file = "SVM01_npp_d20030125_t0847056_e0848301_b00015_c20090513182937523620_gisf_pop.h5"
        geo_file     = h5py.File("%s/%s" %(dir,"GMODO_npp_d20030125_t0847056_e0848301_b00015_c20090513182937526121_gisf_pop.h5"))
        
        #create file 
        output_file  =  h5py.File(output_filename ,"w")
        
        
        #create type
        f32 = numpy.dtype('<f4')
        ui16 = numpy.dtype('<u2')
        i16 = numpy.dtype('<i2')
        
        #return the geolocation info after the tie-point grid treatment
        geo_info = self.get_reduced_tie_point_grid_info()
        
        print("keys = %s\n" % (geo_info.keys()))
         
        output_file.create_dataset('Latitude',  data = geo_info['lat'].astype('float32'), dtype = f32)
        output_file.create_dataset('Longitude', data = geo_info['lon'].astype('float32'), dtype = f32)
        
       
        # convert solar zenith angle  => Range 0-180
        # convert solar azimuth angle => Range -180-180
        numpy.around(geo_info['sol_za'], 3, geo_info['sol_za'])
        geo_info['sol_za'] *= 1000
        geo_info['sol_za'] = geo_info['sol_za'].astype('uint16')
        out_sol_za = geo_info['sol_za']
        print("out_sol_za min %d, max %d , range %d\n" % (numpy.min(out_sol_za), numpy.max(out_sol_za), (numpy.max(out_sol_za)- numpy.min(out_sol_za)) ))
        
        numpy.around(geo_info['sol_aa'], 2, geo_info['sol_aa'])
        geo_info['sol_aa'] *= 100
        geo_info['sol_aa'] = geo_info['sol_aa'].astype('int16')
        out_sol_aa = geo_info['sol_aa']
        print("out_sol_aa min %d, max %d , range %d\n" % (numpy.min(out_sol_aa), numpy.max(out_sol_aa), (numpy.max(out_sol_aa)- numpy.min(out_sol_aa)) ))
        
        # convert sat zenith angle  => Range 0-180
        # convert sat azimuth angle => Range -180-180
        out_sat_za = geo_info['sat_za']
        numpy.around(out_sat_za, 3, out_sat_za)
        out_sat_za = out_sat_za * 1000
        out_sat_za = out_sat_za.astype('uint16')
        print("out_sat_za min %d, max %d , range %d\n" % (numpy.min(out_sat_za), numpy.max(out_sat_za), (numpy.max(out_sat_za)- numpy.min(out_sat_za)) ))
        
        out_sat_aa = geo_info['sat_aa']
        numpy.around(out_sat_aa, 2, out_sat_aa)
        out_sat_aa = out_sat_aa * 100
        out_sat_aa = out_sat_aa.astype('int16')
        print("out_sat_aa min %d, max %d , range %d\n" % (numpy.min(out_sat_aa), numpy.max(out_sat_aa), (numpy.max(out_sat_aa)- numpy.min(out_sat_aa)) ))
        
        output_file.create_dataset('SolarZenithAngle',      data = out_sol_za, dtype = ui16)
        output_file.create_dataset('SolarAzimuthAngle',     data = out_sol_aa, dtype = i16)
        output_file.create_dataset('SatelliteZenithAngle',  data = out_sat_za, dtype = ui16)
        output_file.create_dataset('SatelliteAzimuthAngle', data = out_sat_aa, dtype = i16)
        
        # extract the rest from the geolocation file
        output_file['MidTime']               = geo_info['MidTime'][:]
        output_file['StartTime']             = geo_info['StartTime'][:]
        output_file['ModeGran_Geo']          = geo_info['ModeGran_Geo'][:]
        output_file['ModeScan_Geo']          = geo_info['ModeScan_Geo'][:]
        output_file['NumberOfScans_Geo']     = geo_info['NumberOfScans_Geo'][:]
        output_file['QF1_SCAN_VIIRS_SDRGEO'] = geo_info['QF1_SCAN_VIIRS_SDRGEO'][:]
        output_file['QF2_VIIRSSDRGEO']       = geo_info['QF2_VIIRSSDRGEO'][:]
    
     
        """for rad_pref in self.RADIANCE_PREFIX_LIST:
            
            rad_fn = self.list_of_files[rad_pref]
            
            input_file   =  h5py.File(rad_fn ,"r")
            
            extract_radiance(output_file, input_file)   
        """ 
            
        output_file.close()
                

if __name__ == '__main__':
    
    reducer = VIIRSReducer()
    
    files = reducer.locate_files_for_granule('/homespace/gaubert/viirs/Mband-SDR', 'd20030125_t0847056_e0848301_b00015')
    
    info = reducer.get_reduced_tie_point_grid_info()
    
    print("nb_points for reduced grid = %d\n" % ( len(info['lat']) * len(info['lat'][0])))
    
    info = reducer.get_tie_point_zone_grid_info()
    
    print("nb_points for reduced grid = %d\n" % ( len(info['lat']) * len(info['lat'][0])))
    
    reducer.create_aggregated_viirs_dataset('/tmp/eumetsat_regional_viirs.h5')
    
    
    