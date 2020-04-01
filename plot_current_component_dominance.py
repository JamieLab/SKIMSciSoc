#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 21 14:23:19 2018

Spatial plots of the dominance of each current component.
nEkman vs nGeostrophic
nStokes vs nTotal
Equivalent calculations for skimulator data

@author: Thomas Holding
"""

import cPickle as pickle;
import numpy as np;
import matplotlib.pyplot as plt;
import python_util.parameter_sets as ps;
import python_util.skim_utilities as su;
from os import path;
from mpl_toolkits.basemap import Basemap;


def get_month_indices(months, numYears):
    indices = [];
    for year in range(numYears):
        for month in months:
            indices.append(year*12 + month);
    return indices;


#Select and return data for given set of months each year
#Returned data is reformatted into an array and stored in a QuickStruct, and contains:
#   overall mean and stddev of ekman, geostrophic and stokes current for the whole period.
#months defines the months of the year to use (e.g. 0=Jan, 11=Dec)
def get_processed_data(months, data):
    if len(data) % 12 != 0:
        raise ValueError ("get_processed_data only supports data for whole years (multiples of 12 months)");
    
    #Create the indices to use
    numYears = int(len(data)/12);
    indices = get_month_indices(months, numYears);
    
    #Create arrays to store data for ekman, geostrophic and stokes data
    length = 0;
    for i in indices:
        if len(data[i]) > length:
            length = len(data[i]);
    xcoords = np.empty((length, len(indices))); xcoords[:] = np.nan;
    ycoords = np.empty((length, len(indices))); ycoords[:] = np.nan;
    ekmanVals = np.empty((length, len(indices))); ekmanVals[:] = np.nan;
    geostrophicVals = np.empty((length, len(indices))); geostrophicVals[:] = np.nan;
    stokesVals = np.empty((length, len(indices))); stokesVals[:] = np.nan;
    stokesPassMask = np.empty((length, len(indices))); stokesPassMask[:] = np.nan; #Filters where stokes is important?
    
    #Copy data:
    for i, index in enumerate(indices):
        xcoords[:,i] = [curData.indexX for curData in data[index]];
        ycoords[:,i] = [curData.indexY for curData in data[index]];
        ekmanVals[:,i] = [curData.ekmanProportionGeoEk for curData in data[index]];
        geostrophicVals[:,i] = [curData.geostrophicProportionGeoEk for curData in data[index]];
        stokesVals[:,i] = [curData.stokesProportionOfTotal for curData in data[index]];
        stokesPassMask[:,i] = [curData.stokesMaskPass for curData in data[index]];
    
    #Calculate means and standard deviations
    results = su.QuickStruct(); #Store results in a struct
    results.xcoords = np.nanmean(xcoords, axis=1);
    results.ycoords = np.nanmean(ycoords, axis=1);
    results.ekmanProportions = np.nanmean(ekmanVals, axis=1);
    results.ekmanProportionsSD = np.nanstd(ekmanVals, axis=1);
    results.geostrophicProportions = np.nanmean(geostrophicVals, axis=1);
    results.geostrophicProportionsSD = np.nanstd(geostrophicVals, axis=1);
    results.stokesProportions = np.nanmean(stokesVals, axis=1);
    results.stokesProportionsSD = np.nanstd(stokesVals, axis=1);
    results.stokesPassMask = np.nanmean(stokesPassMask, axis=1);
    
    #Remove nans
    toKeep = (np.isnan(results.ekmanProportions)==False) | (np.isnan(results.geostrophicProportions)==False) | (np.isnan(results.stokesProportions)==False) | (np.isnan(results.stokesPassMask)==False);
    results.xcoords = results.xcoords[toKeep];
    results.ycoords = results.ycoords[toKeep];
    results.ekmanProportions = results.ekmanProportions[toKeep];
    results.ekmanProportionsSD = results.ekmanProportionsSD[toKeep];
    results.geostrophicProportions = results.geostrophicProportions[toKeep];
    results.geostrophicProportionsSD = results.geostrophicProportionsSD[toKeep];
    results.stokesProportions = results.stokesProportions[toKeep];
    results.stokesProportionsSD = results.stokesProportionsSD[toKeep];
    results.stokesPassMask = results.stokesPassMask[toKeep];
    
    return results;
    

resampleFraction = 1; #Set to 1 to disable
plot123 = True;
plot789 = True;
plotEkmanVsGeostrophic = True; #abs ekman vs across-shelf geostrophic
plotStokesTotal = True; #across-shelf stokes vs across-shelf  ???


params = ps.get_current_params();
if params.paramsetName != "global": #Could actually use the shelf coordinates filename (which is also stored in params)...
    raise ValueError("This plotting script is intended only for the 'global' parameter set. Significant adaptation is required for use with any other datasets that may change the shelf-coordinates.");
if ('allData' in globals()) == False:
    ans = raw_input("Press key to read in 'allData', ctrl+c to cancel...");
    allData = pickle.load(open(path.join("data", params.paramsetName, "current_data", "surface_currents_"+params.paramsetName+".p"), "rb"));

#Mask data to restrict analysis
#for monthIndex in range(len(allData)):
#    print "Applying mask to slice", monthIndex;
#    data[monthIndex] = params.contourMaskFunc(allData[monthIndex], params);


#global plot settings
figsize = (8.5,6.5);
ticksize = 10;

#Only keep one in ever 'rate' data points
def resample_data(data, rate):
    toKeep = range(0, len(data.ekmanProportions), rate);
    data.ekmanProportions = data.ekmanProportions[toKeep];
    data.ekmanProportionsSD = data.ekmanProportionsSD[toKeep];
    data.geostrophicProportions = data.geostrophicProportions[toKeep];
    data.geostrophicProportionsSD = data.geostrophicProportionsSD[toKeep];
    data.stokesPassMask = data.stokesPassMask[toKeep];
    data.stokesProportions = data.stokesProportions[toKeep];
    data.stokesProportionsSD = data.stokesProportionsSD[toKeep];
    data.xcoords = data.xcoords[toKeep];
    data.ycoords = data.ycoords[toKeep];
    

####################################
# Process data for Jan, Feb, March #
####################################
if plot123:
    data = get_processed_data([0, 1, 2], allData);
    if resampleFraction > 1:
        resample_data(data, resampleFraction);
    
    #Converts coordinates into lon/lat.
    lats = []; lons = [];
    for i in range(len(data.xcoords)):
        lon, lat = su.convert_index_to_lonlat(data.xcoords[i], data.ycoords[i], params.pixelRes);
        lats.append(lat); lons.append(lon);
    
    
    ###############n#############################
    # Plot Jan, Feb, Mar nekman vs ngeostrophic #
    ##################n##########################
    if plotEkmanVsGeostrophic:
        plt.figure(figsize=figsize);
        mapFig = Basemap(llcrnrlon=-180.0, llcrnrlat=-90, urcrnrlon=180.0, urcrnrlat=90.0, resolution='l', projection='cyl', lat_0 = 39.5, lon_0 = -3.25);
        mapFig.drawmapboundary(fill_color=(0.85, 0.85, 0.85));
        mapFig.fillcontinents(color=(0.5, 0.5, 0.5), zorder=1);
        mapFig.drawmeridians(np.arange(0, 360, 60), labels=[0,0,0,1], color=(0.3, 0.3, 0.3), fontsize=ticksize);
        mapFig.drawparallels(np.arange(-90, 90, 30), labels=[1,0,0,0], color=(0.3, 0.3, 0.3), fontsize=ticksize);
        mapFig.scatter(lons, lats, latlon=True, c=data.geostrophicProportions, marker='o', cmap=plt.cm.RdBu);
        
        cbar = plt.colorbar(orientation="horizontal", ticks=[0.0, 0.5, 1.0]);
        plt.clim(0, 1);
        cbar.set_label("relative strength of geostrophic versus Ekman current", fontsize=ticksize);
        cbar.ax.set_yticklabels(["0.0", "0.5", "1.0"]);
        cbar.ax.tick_params(labelsize=ticksize);
        plt.title("January, February, March, "+str(params.start_year)+"-"+str(params.end_year), fontsize=ticksize);
        #plt.savefig("plots/current_component_dominance/ekman_vs_geostrophic_123_"+str(params.start_year)+"_"+str(params.end_year)+".pdf");
        plt.savefig("plots/current_component_dominance/ekman_vs_geostrophic_123_"+str(params.start_year)+"_"+str(params.end_year)+".png");
    
    ###############n#############################
    # Plot Jan, Feb, Mar nstokes vs ntotal      #
    ##################n##########################
    if plotStokesTotal:
        passed = np.where(data.stokesPassMask != 0.0);
        
        plt.figure(figsize=figsize);
        mapFig = Basemap(llcrnrlon=-180.0, llcrnrlat=-90, urcrnrlon=180.0, urcrnrlat=90.0, resolution='l', projection='cyl', lat_0 = 39.5, lon_0 = -3.25);
        mapFig.drawmapboundary(fill_color=(0.85, 0.85, 0.85));
        mapFig.fillcontinents(color=(0.5, 0.5, 0.5), zorder=1);
        mapFig.drawmeridians(np.arange(0, 360, 60), labels=[0,0,0,1], color=(0.3, 0.3, 0.3), fontsize=ticksize);
        mapFig.drawparallels(np.arange(-90, 90, 30), labels=[1,0,0,0], color=(0.3, 0.3, 0.3), fontsize=ticksize);
        #mapFig.scatter(np.array(lons), np.array(lats), latlon=True, c=nstokesProportions, marker='o', cmap=plt.cm.RdBu);
        mapFig.scatter(np.array(lons)[passed], np.array(lats)[passed], latlon=True, c=data.stokesProportions[passed], marker='o', cmap=plt.cm.YlOrRd);
        
        print "a, max:", np.max(data.stokesProportions[passed]);
        
        #tickToUse = [0.0, 0.5, 1.0];
        ticksToUse = [0.0, 0.33, 0.66];
        cbar = plt.colorbar(orientation="horizontal", ticks=ticksToUse);
        plt.clim(ticksToUse[0], ticksToUse[-1]);
        cbar.set_label("relative strength of stokes versus total current", fontsize=ticksize);
        cbar.ax.set_yticklabels([str(v) for v in ticksToUse]);
        plt.title("January, February, March, "+str(params.start_year)+"-"+str(params.end_year), fontsize=ticksize);
        #plt.savefig("plots/current_component_dominance/stokes_vs_total_123_"+str(params.start_year)+"_"+str(params.end_year)+".pdf");
        plt.savefig("plots/current_component_dominance/stokes_vs_total_123_"+str(params.start_year)+"_"+str(params.end_year)+".png");






##################################
# Process data for Jul, Aug, Sep #
##################################
if plot789:
    data = get_processed_data([6, 7, 8], allData);
    if resampleFraction > 1:
        resample_data(data, resampleFraction);
    
    #Converts coordinates into lon/lat.
    lats = []; lons = [];
    for i in range(len(data.xcoords)):
        lon, lat = su.convert_index_to_lonlat(data.xcoords[i], data.ycoords[i], params.pixelRes);
        lats.append(lat); lons.append(lon);
    
    
    ###############n#############################
    # Plot Jan, Feb, Mar nekman vs ngeostrophic #
    ##################n##########################
    if plotEkmanVsGeostrophic:
        plt.figure(figsize=figsize);
        mapFig = Basemap(llcrnrlon=-180.0, llcrnrlat=-90, urcrnrlon=180.0, urcrnrlat=90.0, resolution='l', projection='cyl', lat_0 = 39.5, lon_0 = -3.25);
        mapFig.drawmapboundary(fill_color=(0.85, 0.85, 0.85));
        mapFig.fillcontinents(color=(0.5, 0.5, 0.5), zorder=1);
        mapFig.drawmeridians(np.arange(0, 360, 60), labels=[0,0,0,1], color=(0.3, 0.3, 0.3), fontsize=ticksize);
        mapFig.drawparallels(np.arange(-90, 90, 30), labels=[1,0,0,0], color=(0.3, 0.3, 0.3), fontsize=ticksize);
        mapFig.scatter(lons, lats, latlon=True, c=data.geostrophicProportions, marker='o', cmap=plt.cm.RdBu);
        
        cbar = plt.colorbar(orientation="horizontal", ticks=[0.0, 0.5, 1.0]);
        plt.clim(0, 1);
        cbar.set_label("relative strength of geostrophic versus Ekman current", fontsize=ticksize);
        cbar.ax.set_yticklabels(["0.0", "0.5", "1.0"]);
        plt.title("July, August, September, "+str(params.start_year)+"-"+str(params.end_year), fontsize=ticksize);
        #plt.savefig("plots/current_component_dominance/ekman_vs_geostrophic_789_"+str(params.start_year)+"_"+str(params.end_year)+".pdf");
        plt.savefig("plots/current_component_dominance/ekman_vs_geostrophic_789_"+str(params.start_year)+"_"+str(params.end_year)+".png");
    
    
    ###############n#############################
    # Plot Jan, Feb, Mar nstokes vs ntotal      #
    ##################n##########################
    if plotStokesTotal:
        passed = np.where(data.stokesPassMask != 0.0);
        
        plt.figure(figsize=figsize);
        mapFig = Basemap(llcrnrlon=-180.0, llcrnrlat=-90, urcrnrlon=180.0, urcrnrlat=90.0, resolution='l', projection='cyl', lat_0 = 39.5, lon_0 = -3.25);
        mapFig.drawmapboundary(fill_color=(0.85, 0.85, 0.85));
        mapFig.fillcontinents(color=(0.5, 0.5, 0.5), zorder=1);
        mapFig.drawmeridians(np.arange(0, 360, 60), labels=[0,0,0,1], color=(0.3, 0.3, 0.3), fontsize=ticksize);
        mapFig.drawparallels(np.arange(-90, 90, 30), labels=[1,0,0,0], color=(0.3, 0.3, 0.3), fontsize=ticksize);
        #mapFig.scatter(np.array(lons), np.array(lats), latlon=True, c=nstokesProportions, marker='o', cmap=plt.cm.RdBu);
        mapFig.scatter(np.array(lons)[passed], np.array(lats)[passed], latlon=True, c=data.stokesProportions[passed], marker='o', cmap=plt.cm.YlOrRd);
        
        print "b, max:", np.max(data.stokesProportions[passed]);
        
        #tickToUse = [0.0, 0.5, 1.0];
        ticksToUse = [0.0, 0.33, 0.66];
        cbar = plt.colorbar(orientation="horizontal", ticks=ticksToUse);
        plt.clim(ticksToUse[0], ticksToUse[-1]);
        cbar.set_label("relative strength of stokes versus total current", fontsize=ticksize);
        cbar.ax.set_yticklabels([str(v) for v in ticksToUse]);
        plt.title("July, August, September, "+str(params.start_year)+"-"+str(params.end_year), fontsize=ticksize);
        #plt.savefig("plots/current_component_dominance/stokes_vs_total_789_"+str(params.start_year)+"_"+str(params.end_year)+".pdf");
        plt.savefig("plots/current_component_dominance/stokes_vs_total_789_"+str(params.start_year)+"_"+str(params.end_year)+".png");
    
del data
