# -*- coding: utf-8 -*-
'''
Created on Feb 9, 2011

@author: guillaume.aubert@eumetsat.int
'''
import math

def deg_to_rad(val):
    """ convert degree values in radians """
    return (val * math.pi)/180.00

def rad_to_deg(value):
    """ convert radians to degrees """
    (value * 180.00)/ math.pi

def distance(lat1, lon1, lat2, lon2):
    """ calculate distance between 2 points using the haversine formula 
        
        R = earth’s radius (mean radius = 6,371km)
        Δlat = lat2− lat1
        Δlong = long2− long1
        a = sin²(Δlat/2) + cos(lat1).cos(lat2).sin²(Δlong/2)
        c = 2.atan2(√a, √(1−a))
        d = R.c  
    """
    
    R     = 6371
    d_lat = deg_to_rad(lat2 - lat1)
    d_lon = deg_to_rad(lon2 - lon1)
    a     = math.sin(d_lat/2) * math.sin(d_lat/2) + \
            math.cos(deg_to_rad(lat1)) * math.cos(deg_to_rad(lat2)) * \
            math.sin(d_lon/2) * math.sin(d_lon/2)
    c     = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
    
    return R * c
            
    