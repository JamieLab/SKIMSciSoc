#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created by Daniel J. Ford (d.ford@exeter.ac.uk)
Date: 03/2023
This script retrieves the globcurrent and wavewatch data required for the analysis from
the respective FTP servers.
Globcurrent: ftp://eftp.ifremer.fr - /data/globcurrent/v3.0/global_025_deg
Globcurrent requires signing up at: http://globcurrent.ifremer.fr/products-data/data-catalogue
Enter username and password for GlobCurrent FTP in downloadekmansurf() and downloadekmandepth()
and downloadgeostrophic() functions
WaveWatch: ftp.ifremer.fr - /ifremer/ww3/HINDCAST/GLOBAL/YYYY_CFSR/uss and /ifremer/ww3/HINDCAST/GLOBAL/YYYY_CFSR/wnd
- WaveWatch is an anouymous FTP
"""

import os
import glob
import time
import ftplib
import datetime

save_loc = 'D:/SKIM_Paper_Data/SKIM/downloaded_data/' # Location to place data structure
start_yr = 1993 #Start year of analysis
end_yr = 2016 # End year of analysis

download_ekman = False
download_geo = True
download_stokes = False
download_wind = True
download_Hs = True

def checkfileexist(file):
    #print(file)
    g = glob.glob(file)
    #print(g)
    if not g:
        return 0
    else:
        return 1

def makefolder(fold):
    if not os.path.exists(fold):
        os.makedirs(fold)

def julstr3(jul):
    if jul < 10:
        return '00'+str(jul)
    elif jul<100:
        return '0'+str(jul)
    else:
        return str(jul)

def ftpretrieve(file,ftppath,savepath,ftpserver,username='' ,password=''):
    print('Retreiving file: '+file)
    print('from FTP Server: ' + ftpserver)
    ftp = ftplib.FTP(ftpserver)
    ftp.login(username, password)
    ftp.cwd(ftppath)
    ftp.retrbinary("RETR " + file, open(os.path.join(savepath,file), 'w+b').write)
    ftp.quit()

def downloadstokes(save_loc,startyr = 1993,end_yr = 2016):
    #Function that cyles through the wavewatch archive and downloads the data.
    ftpserver='ftp.ifremer.fr'
    save_loc = save_loc+'/stokes_data_3h/'
    print('Downloading Stokes data from WaveWatch between ' + str(start_yr) + ' and ' + str(end_yr))
    yr = start_yr
    mon = 1
    while yr<end_yr+1:
        dat = datetime.datetime(yr,mon,1)
        ftp_path = dat.strftime('/ifremer/ww3/HINDCAST/GLOBAL/%Y_CFSR/uss/')
        print(ftp_path)
        file = dat.strftime('WW3-GLOB-30M_%Y%m_uss.nc')
        print(file)
        if checkfileexist(save_loc+file) == 0:
            ftpretrieve(file,ftp_path,save_loc,ftpserver)
        mon = mon+1
        if mon == 13:
            mon = 1
            yr = yr+1
    print('Stokes download complete!')

def downloadwind(save_loc,startyr = 1993,end_yr = 2016):
    #Function that cyles through the wavewatch archive and downloads the data.
    ftpserver='ftp.ifremer.fr'
    save_loc = save_loc+'/wavewatch_wnd/'
    print('Downloading Wind data from WaveWatch between ' + str(start_yr) + ' and ' + str(end_yr))
    yr = start_yr
    mon = 1
    while yr<end_yr+1:
        dat = datetime.datetime(yr,mon,1)
        ftp_path = dat.strftime('/ifremer/ww3/HINDCAST/GLOBAL/%Y_CFSR/wnd/')
        print(ftp_path)
        file = dat.strftime('WW3-GLOB-30M_%Y%m_wnd.nc')
        print(file)
        if checkfileexist(save_loc+file) == 0:
            ftpretrieve(file,ftp_path,save_loc,ftpserver)
        mon = mon+1
        if mon == 13:
            mon = 1
            yr = yr+1
    print('Wind download complete!')

def downloadHs(save_loc,startyr = 1993,end_yr = 2016):
    #Function that cyles through the wavewatch archive and downloads the data.
    ftpserver='ftp.ifremer.fr'
    save_loc = save_loc+'/wavewatch_Hs/'
    print('Downloading Significant wave height data from WaveWatch between ' + str(start_yr) + ' and ' + str(end_yr))
    yr = start_yr
    mon = 1
    while yr<end_yr+1:
        dat = datetime.datetime(yr,mon,1)
        ftp_path = dat.strftime('/ifremer/ww3/HINDCAST/GLOBAL/%Y_CFSR/hs/')
        print(ftp_path)
        file = dat.strftime('WW3-GLOB-30M_%Y%m_hs.nc')
        print(file)
        if checkfileexist(save_loc+file) == 0:
            ftpretrieve(file,ftp_path,save_loc,ftpserver)
        mon = mon+1
        if mon == 13:
            mon = 1
            yr = yr+1
    print('Significant wave height download complete!')

def downloadekmansurf(save_loc,startyr = 1993,end_yr = 2016):
    #Function that cyles through the wavewatch archive and downloads the data.
    ftpserver='eftp.ifremer.fr'
    username = ''
    password = ''
    save_loc = save_loc+'/mean_daily_ekman_currents/'
    print('Downloading Surface Ekman data from GlobCurrent between ' + str(start_yr) + ' and ' + str(end_yr))
    dat = datetime.datetime(start_yr,1,1)
    while dat.year<end_yr+1:
        jul_str = julstr3(dat.timetuple().tm_yday)

        ftp_path = '/data/globcurrent/v3.0/global_025_deg/ekman_hs/'+str(dat.year)+'/'+jul_str+'/'
        print(ftp_path)
        file = dat.strftime('%Y%m%d')+'-GLOBCURRENT-L4-CURekm_hs-ERAWS_EEM-v03.0-fv01.0.nc'
        print(file)
        if checkfileexist(save_loc+file) == 0:
            ftpretrieve(file,ftp_path,save_loc,ftpserver,username=username,password=password)
        dat = dat + datetime.timedelta(days=1)
    print('Ekman Surface download complete!')

def downloadekmandepth(save_loc,startyr = 1993,end_yr = 2016):
    #Function that cyles through the wavewatch archive and downloads the data.
    ftpserver='eftp.ifremer.fr'
    username = ''
    password = ''
    save_loc = save_loc+'/mean_daily_ekman_currents_15m/'
    print('Downloading 15m Ekman data from GlobCurrent between ' + str(start_yr) + ' and ' + str(end_yr))
    dat = datetime.datetime(start_yr,1,1)
    while dat.year<end_yr+1:
        jul_str = julstr3(dat.timetuple().tm_yday)

        ftp_path = '/data/globcurrent/v3.0/global_025_deg/ekman_15m/'+str(dat.year)+'/'+jul_str+'/'
        print(ftp_path)
        file = dat.strftime('%Y%m%d')+'-GLOBCURRENT-L4-CURekm_15m-ERAWS_EEM-v03.0-fv01.0.nc'
        print(file)
        if checkfileexist(save_loc+file) == 0:
            ftpretrieve(file,ftp_path,save_loc,ftpserver,username=username,password=password)
        dat = dat + datetime.timedelta(days=1)
    print('Ekman 15m download complete!')

def downloadgeostrophic(save_loc,startyr = 1993,end_yr = 2016):
    #Function that cyles through the wavewatch archive and downloads the data.
    ftpserver='eftp.ifremer.fr'
    username = ''
    password = ''
    save_loc = save_loc+'/mean_daily_geostrophic_currents/'
    print('Downloading Geostrophic data from GlobCurrent between ' + str(start_yr) + ' and ' + str(end_yr))
    dat = datetime.datetime(start_yr,1,1)
    while dat.year<end_yr+1:
        jul_str = julstr3(dat.timetuple().tm_yday)

        ftp_path = '/data/globcurrent/v3.0/global_025_deg/geostrophic/'+str(dat.year)+'/'+jul_str+'/'
        print(ftp_path)
        file = dat.strftime('%Y%m%d')+'000000-GLOBCURRENT-L4-CURgeo_0m-ALT_OI-v03.0-fv01.0.nc'
        print(file)
        if checkfileexist(save_loc+file) == 0:
            ftpretrieve(file,ftp_path,save_loc,ftpserver,username=username,password=password)
        dat = dat + datetime.timedelta(days=1)
    print('Geostrophic download complete!')

if download_stokes:
    downloadstokes(save_loc)
if download_ekman:
    downloadekmansurf(save_loc)
    downloadekmandepth(save_loc)
if download_wind:
    downloadwind(save_loc)
if download_geo:
    downloadgeostrophic(save_loc)
if download_Hs:
    downloadHs(save_loc)
#ftpretrieve('WW3-GLOB-30M_199305_uss.nc','/ifremer/ww3/HINDCAST/GLOBAL/1993_CFSR/uss/',save_loc,)
