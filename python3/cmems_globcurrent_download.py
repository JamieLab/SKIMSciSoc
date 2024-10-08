#!/usr/bin/env python3
"""
Created by Daniel J. Ford (d.ford@exeter.ac.uk)
Date: 08/2024

"""
import getpass
import datetime
import copernicusmarine as cmmarine
import os
import glob
import sys
sys.path.append('C:\\Users\\df391\\OneDrive - University of Exeter\\Post_Doc_ESA_Contract\\OceanICU')
sys.path.append('C:\\Users\\df391\\OneDrive - University of Exeter\\Post_Doc_ESA_Contract\\OceanICU\Data_Loading')

import data_utils as du
import cmems_glorysv12_download as cm

def load_globcurrent_monthly(loc,start_yr = 1993,end_yr = 2016):
    """
    """
    OUTPUT_DIRECTORY = loc
    yr = start_yr
    mon = 1
    while yr <= end_yr:
        date_min_v = datetime.datetime(yr,mon,1,0,0,0);
        date_max = datetime.datetime(yr,mon,27,23,59,59); date_max = date_max.strftime('%Y-%m-%d %H:%M:%S')
        OUTPUT_TEMP = os.path.join(OUTPUT_DIRECTORY,str(yr))
        du.makefolder(OUTPUT_TEMP)
        OUTPUT_FILENAME = date_min_v.strftime(f'%Y_%m_CMEMS_GLOBCURRENT.nc')
        date_min = date_min_v.strftime('%Y-%m-%dT%H:%M:%S')
        print(OUTPUT_FILENAME)
        if not du.checkfileexist(os.path.join(OUTPUT_TEMP,OUTPUT_FILENAME)):
            cmmarine.subset(
              dataset_id="cmems_obs_mob_glo_phy-cur_my_0.25deg_P1M-m",
              variables=["err_ue", "err_ugos", "err_ve", "err_vgos", "ue", "ugos", "ve", "vgos"],
              minimum_longitude=-179.875,
              maximum_longitude=179.875,
              minimum_latitude=-89.875,
              maximum_latitude=89.875,
              start_datetime=date_min,
              end_datetime=date_max,
              output_filename=OUTPUT_FILENAME,
              output_directory=OUTPUT_TEMP,
              force_download=True
            )
        mon = mon+1
        if mon == 13:
            yr = yr+1
            mon=1

def load_glorysv12_monthly(loc,start_yr = 1993,end_yr = 2016):
    """
    """
    OUTPUT_DIRECTORY = loc
    yr = start_yr
    mon = 1
    while yr <= end_yr:
        date_min_v = datetime.datetime(yr,mon,1,0,0,0);
        date_max = datetime.datetime(yr,mon,27,23,59,59); date_max = date_max.strftime('%Y-%m-%d %H:%M:%S')
        OUTPUT_TEMP = os.path.join(OUTPUT_DIRECTORY,str(yr))
        du.makefolder(OUTPUT_TEMP)
        OUTPUT_FILENAME = date_min_v.strftime(f'%Y_%m_CMEMS_GLORYSV12.nc')
        date_min = date_min_v.strftime('%Y-%m-%dT%H:%M:%S')
        print(OUTPUT_FILENAME)
        if not du.checkfileexist(os.path.join(OUTPUT_TEMP,OUTPUT_FILENAME)):
            cmmarine.subset(
              dataset_id="cmems_mod_glo_phy_my_0.083deg_P1M-m",
              variables=["uo", "vo"],
              minimum_longitude=-180,
              maximum_longitude=179.9166717529297,
              minimum_latitude=-80,
              maximum_latitude=90,
              minimum_depth=0.49402499198913574,
              maximum_depth=0.49402499198913574,
              start_datetime=date_min,
              end_datetime=date_max,
              output_filename=OUTPUT_FILENAME,
              output_directory=OUTPUT_TEMP,
              force_download=True
            )
        mon = mon+1
        if mon == 13:
            yr = yr+1
            mon=1

lon,lat = du.reg_grid(lat=0.25,lon=0.25)
# load_globcurrent_monthly('D:/SKIM_Paper_Data/SKIM/downloaded_data/CMEMS',start_yr = 1993,end_yr = 2016)
# load_glorysv12_monthly('D:/SKIM_Paper_Data/SKIM/downloaded_data/CMEMS_GLORYSV12',start_yr = 1993,end_yr = 2016)
cm.cmems_average('D:/SKIM_Paper_Data/SKIM/downloaded_data/CMEMS_GLORYSV12','D:/SKIM_Paper_Data/SKIM/downloaded_data/CMEMS_GLORYSV12_025',start_yr = 1993,end_yr=1994,log = lon,lag=lat,variable = ['uo','vo'])
