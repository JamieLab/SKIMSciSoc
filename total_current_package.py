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
depth = 800
bath_file = 'D:/SKIM/GEBCO_bathymetry_0.25x0.25deg.nc'
out_file = 'C:/Users/df391/OneDrive - University of Exeter/Post_Doc_Covex_Seascape/Shutler_Cross_Shelf_Transport/Data/'+str(depth)+'m_total_current_data_Shutler_et_al.nc'
c = Dataset(bath_file,'r')
lat = np.array(c.variables['lat'])
lon = np.array(c.variables['lon'])
elev = np.array(c.variables['mean_depth'])
c.close()



params = ps.get_current_params();
if params.paramsetName != "global": #Could actually use the shelf coordinates filename (which is also stored in params)...
    raise ValueError("This plotting script is intended only for the 'global' parameter set. Significant adaptation is required for use with any other datasets that may change the shelf-coordinates.");
if ('allData' in globals()) == False:
    allData = pickle.load(open(path.join("D:/SKIM", "current_data", "surface_currents_"+params.paramsetName+"_"+str(depth)+"m.p"), "rb"));

yr = 1993
mon = 1
time = []
total_cur = np.zeros((len(lon),len(lat),len(allData)))
total_cur[:] = np.nan

ekman = np.copy(total_cur)
geostrophic = np.copy(total_cur)
stokes = np.copy(total_cur)
segment = np.zeros((len(lon),len(lat)))
segment[:] = np.nan
for j in range(len(allData)):
    #print(allData[j][100])
    for i in range(len(allData[j])):
          total_cur[allData[j][i].indexX,allData[j][i].indexY,j] = allData[j][i].totalcurrent

          ekman[allData[j][i].indexX,allData[j][i].indexY,j] = allData[j][i].nEkmanAcrossShelf
          geostrophic[allData[j][i].indexX,allData[j][i].indexY,j] = allData[j][i].nGeostrophicAcrossShelf
          stokes[allData[j][i].indexX,allData[j][i].indexY,j] = allData[j][i].nStokesAcrossShelf
          if j == 0:
              segment[allData[j][i].indexX,allData[j][i].indexY] = allData[j][i].segmentDistance
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

tot = c.createVariable('total_current','f4',('longitude','latitude','time'),zlib=True)
tot[:] = total_cur
tot.units = 'ms-1'
tot.long_name = 'Total surface across-shelf break current'

tot = c.createVariable('ekman_current','f4',('longitude','latitude','time'),zlib=True)
tot[:] = ekman
tot.units = 'ms-1'
tot.long_name = 'Ekman surface across-shelf break current'

tot = c.createVariable('stokes_current','f4',('longitude','latitude','time'),zlib=True)
tot[:] = stokes
tot.units = 'ms-1'
tot.long_name = 'Stokes surface across-shelf break current'

tot = c.createVariable('geostrophic_current','f4',('longitude','latitude','time'),zlib=True)
tot[:] = geostrophic
tot.units = 'ms-1'
tot.long_name = 'Geostrophic surface across-shelf break current'

tot = c.createVariable('segment_length','f4',('longitude','latitude'),zlib=True)
tot[:] = segment
tot.units = 'm'
tot.long_name = 'Segment length across grid cell'

setattr(c,"date_compiled",datetime.datetime.now().strftime("%d/%m/%Y"))
setattr(c,"complied_by",'Daniel J. Ford (d.ford@exeter.ac.uk)')
setattr(c,"production_code","https://github.com/JamieLab/SKIMSciSoc")

c.close()
