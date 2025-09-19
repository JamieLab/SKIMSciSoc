#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import cPickle as pickle;
import numpy as np;
from netCDF4 import Dataset;
import python_util.parameter_sets as ps;
import python_util.skim_utilities as su;
import python_util.mask_functions as mf;
from os import path;
import datetime

ref_year = 1990
depth_l = [300,400,500,600,700,800]
# depth_l = [500]
bath_file = 'E:/SKIM/GEBCO_bathymetry_0.25x0.25deg.nc'
c = Dataset(bath_file,'r')
lat = np.array(c.variables['lat'])
lon = np.array(c.variables['lon'])
elev = np.array(c.variables['mean_depth'])
c.close()

start_yr = 1993
end_yr = 2016


for depth in depth_l:
    out_file = 'E:/SKIM/netcdf_out/'+str(depth)+'m_total_current_data_Shutler_et_al_CMEMS_'+str(start_yr)+'_'+str(end_yr)+'.nc'
    # out_file = 'E:/SKIM/netcdf_out/Jenny_MAB.nc'

    params = ps.get_global_params(cmems = True)
    # params.paramsetName = 'CMEMS_Jenny_MAB'
    # if params.paramsetName != "global": #Could actually use the shelf coordinates filename (which is also stored in params)...
    #     raise ValueError("This plotting script is intended only for the 'global' parameter set. Significant adaptation is required for use with any other datasets that may change the shelf-coordinates.");
    print(path.join("E:/SKIM", "current_data", "surface_currents_"+params.paramsetName+'_'+str(depth)+"m.p"))
    allData = pickle.load(open(path.join("E:/SKIM", "current_data", "surface_currents_"+params.paramsetName+'_'+str(depth)+"m.p"), "rb"));
    # allData = pickle.load(open(path.join("E:/SKIM", "current_data", "surface_currents_"+params.paramsetName+".p"), "rb"));

    yr = start_yr
    mon = 1
    time = []
    for j in range(len(allData)):
        time.append(datetime.datetime(yr,mon,15))
        mon = mon+1
        if mon == 13:
            yr = yr+1
            mon=1
    time_r = []
    for i in range(len(time)):
        print(time[i])
        time_r.append((time[i] - datetime.datetime(ref_year,1,15)).days)
    print(time_r)
    # vars=['totalcurrent']
    # net_var = ['total_current']
    vars = ['totalcurrent',
        'totalcurrenterr',
        'totalcurrentalong',
        'totalcurrentalongerr',
        'nEkmanAcrossShelf',
        'nEkmanAcrossShelferr',
        'nEkmanAlongShelf',
        'nEkmanAlongShelferr',
        'nGeostrophicAcrossShelf',
        'nGeostrophicAcrossShelferr',
        'nGeostrophicAlongShelf',
        'nGeostrophicAlongShelferr',
        'nStokesAcrossShelf',
        'nStokesAcrossShelferr',
        'nStokesAlongShelf'
        ,'nStokesAlongShelferr']
    net_var = ['total_current',
        'total_current_err',
        'total_along_current',
        'total_along_current_err',
        'ekman_current',
        'ekman_current_err',
        'ekman_current_along',
        'ekman_current_along_err',
        'geostrophic_current',
        'geostrophic_current_err',
        'geostrophic_current_along',
        'geostrophic_current_along_err',
        'stokes_current',
        'stokes_current_err',
        'stokes_current_along'
        ,'stokes_current_along_err']

    c = Dataset(out_file,'w',format='NETCDF4_CLASSIC')
    c.createDimension('latitude', len(lat));
    c.createDimension('longitude', len(lon));
    c.createDimension('time',len(time_r))

    latVar = c.createVariable("latitude", 'f4', ("latitude",) );
    latVar[:] = np.flip(lat)
    latVar.units = "degrees North";
    latVar.long_name = "Latitude";
    latVar.valid_min = -90;
    latVar.valid_max = 90;

    lonVar = c.createVariable("longitude", 'f4', ("longitude",) );
    lonVar[:] = lon
    lonVar.units = "degrees East";
    lonVar.long_name = "Longitude";
    lonVar.valid_min = -180;
    lonVar.valid_max = 180;

    timevar = c.createVariable("time", 'f4', ("time",) );
    timevar[:] = time_r
    timevar.units = 'Days since '+str(ref_year)+'-01-15'
    timevar.long_name = "Time";

    c.close()

    for k in range(len(vars)):
        print(vars[k])
        data_var = np.zeros((len(lon),len(lat),len(allData)))
        data_var[:] = np.nan
        for j in range(len(allData)):
            print(j)
            #print(allData[j][100])
            for i in range(len(allData[j])):
                  data_var[allData[j][i].indexX,allData[j][i].indexY,j] = eval('allData['+str(j)+']['+str(i)+'].'+vars[k])

        c = Dataset(out_file,'a',format='NETCDF4_CLASSIC')
        tot = c.createVariable(net_var[k],'f4',('longitude','latitude','time'),zlib=True)
        tot[:] = data_var
        c.close()

    ### Getting the segment length data :-)
    data_var = np.zeros((len(lon),len(lat)))
    data_var[:] = np.nan
    j=0
    for i in range(len(allData[j])):
          data_var[allData[j][i].indexX,allData[j][i].indexY] = allData[j][i].segmentDistance

    c = Dataset(out_file,'a',format='NETCDF4_CLASSIC')
    tot = c.createVariable('segment_length','f4',('longitude','latitude'),zlib=True)
    tot[:] = data_var
    tot.units = 'm'
    tot.long_name = 'Segment length across grid cell'
    c.close()



    ### Total
    # tot = c.createVariable('total_current','f4',('longitude','latitude','time'),zlib=True)
    # tot[:] = total_cur
    c = Dataset(out_file,'a',format='NETCDF4_CLASSIC')
    c['total_current'].units = 'ms-1'
    c['total_current'].long_name = 'Total surface across-shelf break current'

    c['total_current_err'].units = 'ms-1'
    c['total_current_err'].long_name = 'Total surface across-shelf break current uncertainty'
    c['total_current_err'].confidence_interval = '95 % confidence (2 sigma)'
    #

    c['total_along_current'].units = 'ms-1'
    c['total_along_current'].long_name = 'Total surface along-shelf break current'
    #

    c['total_along_current_err'].units = 'ms-1'
    c['total_along_current_err'].long_name = 'Total surface along-shelf break current uncertainty'
    c['total_along_current_err'].confidence_interval = '95 % confidence (2 sigma)'
    #
    # ### Ekman

    c['ekman_current'].units = 'ms-1'
    c['ekman_current'].long_name = 'Ekman surface across-shelf break current'
    #

    c['ekman_current_err'].units = 'ms-1'
    c['ekman_current_err'].long_name = 'Ekman surface across-shelf break current uncertainty'
    c['ekman_current_err'].confidence_interval = '95 % confidence (2 sigma)'
    #

    c['ekman_current_along'].units = 'ms-1'
    c['ekman_current_along'].long_name = 'Ekman surface along-shelf break current'
    #
    c['ekman_current_along_err'].units = 'ms-1'
    c['ekman_current_along_err'].long_name = 'Ekman surface along-shelf break current uncertainty'
    c['ekman_current_along_err'].confidence_interval = '95 % confidence (2 sigma)'
    #
    # ###Stokes

    c['stokes_current'].units = 'ms-1'
    c['stokes_current'].long_name = 'Stokes surface across-shelf break current'
    #
    c['stokes_current_err'].units = 'ms-1'
    c['stokes_current_err'].long_name = 'Stokes surface across-shelf break current uncertainty'
    c['stokes_current_err'].confidence_interval = '95 % confidence (2 sigma)'
    #
    c['stokes_current_along'].units = 'ms-1'
    c['stokes_current_along'].long_name = 'Stokes surface along-shelf break current'
    #

    c['stokes_current_along_err'].units = 'ms-1'
    c['stokes_current_along_err'].long_name = 'Stokes surface along-shelf break current uncertainty'
    c['stokes_current_along_err'].confidence_interval = '95 % confidence (2 sigma)'
    #
    # ### Geostrophic

    c['geostrophic_current'].units = 'ms-1'
    c['geostrophic_current'].long_name = 'Geostrophic surface across-shelf break current'
    #

    c['geostrophic_current_err'].units = 'ms-1'
    c['geostrophic_current_err'].long_name = 'Geostrophic surface across-shelf break current uncertainty'
    c['geostrophic_current_err'].confidence_interval = '95 % confidence (2 sigma)'
    #

    c['geostrophic_current_along'].units = 'ms-1'
    c['geostrophic_current_along'].long_name = 'Geostrophic surface along-shelf break current'
    #

    c['geostrophic_current_along_err'].units = 'ms-1'
    c['geostrophic_current_along_err'].long_name = 'Geostrophic surface along-shelf break current uncertainty'
    c['geostrophic_current_along_err'].confidence_interval = '95 % confidence (2 sigma)'
    #
    setattr(c,"date_compiled",datetime.datetime.now().strftime("%d/%m/%Y"))
    setattr(c,"complied_by",'Daniel J. Ford (d.ford@exeter.ac.uk)')
    setattr(c,"production_code","https://github.com/JamieLab/SKIMSciSoc")

    c.close()
