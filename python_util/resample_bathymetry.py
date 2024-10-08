#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 10 08:47:49 2018

@author: rr
"""

from netCDF4 import Dataset;
import numpy as np;
import matplotlib.pyplot as plt;
from os import path;
import datetime

def _subgrid_mean(subx, suby, subgridWidth, subgridHeight, data):
    x0 = subx*subgridWidth;
    x1 = (subx+1)*subgridWidth;
    y0 = suby*subgridHeight;
    y1 = (suby+1)*subgridHeight;

    m = np.ma.mean(data[x0:x1, y0:y1]);
    return m;

#Calculates and returns a resampled ocean depth bathymetry.
#inputBathymetryPath - path to input data
#pixelRes - tupple containing (lon, lat) resolution for resampled grid
#outputPath - directory to save resampled bathymetry netCDF file to
def generate_resampled_bathymetry(inputBathymetryPath, pixelRes, testPlot=False, outputPath=None):
    outLonLen = int(360/pixelRes[0]);
    outLatLen = int(180/pixelRes[1]);
    lonRes = pixelRes[0];
    latRes = pixelRes[1];

    #Read in data
    ncin = Dataset(inputBathymetryPath, 'r');
    latLen = len(ncin.variables["lat"]);
    lonLen = len(ncin.variables["lon"]);
    data = ncin.variables["elevation"][:];

    #Apply mask to exclude land
    data.mask = data >= 0.0;

    #Calculate rescale values
    rescaleLat = int(latLen/outLatLen);
    #print(rescaleLat)
    rescaleLon = int(lonLen/outLonLen);

    #Empty unmasked ma array
    outData = np.ma.empty((outLatLen, outLonLen), dtype=float);
    outData.mask = np.full((outLatLen, outLonLen), False);

    #For each output grid coordinate, calculate the mean of the corresponding subgrid from the input data
    for lon in range (0, outLonLen):
        for lat in range (0, outLatLen):
            cMin = _subgrid_mean(lat, lon, rescaleLat, rescaleLon, data);
            #print lon, lat, cMin;
            if np.ma.is_masked(cMin):
                outData[lat,lon] = -1.0;
                outData.mask[lat, lon] = bool(cMin.mask);
            outData[lat,lon] = cMin;

    #Change elevation to depth
    outData = -outData;

    #Free memory
    del data;

    #Display
    if testPlot:
        plt.figure();
        plt.imshow(np.flipud(outData), cmap=plt.cm.jet);
        plt.colorbar();
        plt.title("resampled bathymetry "+str(pixelRes));
        plt.pause(1);

    if outputPath != None:
        #Create output netCDF file, write resampled bathymetry to it and close the file.
        print "Writing bathymetry to:", path.join(outputPath, "GEBCO_bathymetry_"+str(latRes)+"x"+str(lonRes)+"deg.nc"), 'w';
        ncout = Dataset(path.join(outputPath, "GEBCO_bathymetry_"+str(latRes)+"x"+str(lonRes)+"deg.nc"), 'w');
        ny, nx = outData.shape;
        print(outData.shape)
        ncout.createDimension('latitude', ny);
        ncout.createDimension('longitude', nx);
        print(np.arange(-90+(latRes/2.0), 90+(latRes/2.0), latRes).shape)
        print(np.arange(-90+(latRes/2.0), 90+(latRes/2.0), latRes))
        print(np.arange(-180+(lonRes/2.0), 180+(lonRes/2.0), lonRes).shape)
        print(np.arange(-180+(lonRes/2.0), 180+(lonRes/2.0), lonRes))
        latVar = ncout.createVariable("lat", np.dtype(float), ("latitude",) );
        latVar[:] = np.arange(-90+(latRes/2.0), 90+(latRes/2.0), latRes);
        latVar.units = "degrees North";
        latVar.long_name = "Latitude";
        latVar.valid_min = -90;
        latVar.valid_max = 90;

        lonVar = ncout.createVariable("lon", np.dtype(float), ("longitude",) );
        lonVar[:] = np.arange(-180+(lonRes/2.0), 180+(lonRes/2.0), lonRes);
        #np.arange(-179.5, 180.5);
        lonVar.units = "degrees East";
        lonVar.long_name = "Longitude";
        lonVar.valid_min = -180;
        lonVar.valid_max = 180;

        minBathyVar = ncout.createVariable("mean_depth", np.dtype(float), ("latitude", "longitude") );
        minBathyVar[:] = outData;
        minBathyVar.units = "meters (m)";
        minBathyVar.long_name = "Mean ocean depth over a " + str(pixelRes[0]) + "° grid.";
        minBathyVar.description = "description", "Mean ocean depth over a " + str(pixelRes[0]) + "° grid compiled from GEBCO bathymetry (30sec resolution) dataset (see: https://www.gebco.net/data_and_products/gridded_bathymetry_data/). Grid cells containing all land are masked/set to missing_value, where land is defined as altitude >= 0.0. Any 1°x1° cell which contained least one ocean region at the original 30sec resolution has a minimum depth calculated.";
        minBathyVar.input_file = inputBathymetryPath
        minBathyVar.valid_min = 0;
        minBathyVar.valid_max = 10000;
        minBathyVar.fill_value = -1;
        minBathyVar.missing_value = -1;

        setattr(ncout,"date_compiled",datetime.datetime.now().strftime("%d/%m/%Y"))
        setattr(ncout, "description", "Mean ocean depth over a " + str(pixelRes[0]) + "° grid compiled from GEBCO bathymetry dataset (see: https://www.gebco.net/data_and_products/gridded_bathymetry_data/). Grid cells containing all land are masked/set to missing_value, where land is defined as altitude >= 0.0.");
        setattr(ncout, "citation", "The General Bathymetric Chart of the Oceans, Gridded Bathymetry data 30sec dataset: https://www.gebco.net/data_and_products/gridded_bathymetry_data/");
        setattr(ncout, "input_file", inputBathymetryPath);
        ncout.close();

    return outData;
