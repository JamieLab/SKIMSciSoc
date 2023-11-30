#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 21 14:23:19 2018

Spatial plots of the dominance of each current component.
nEkman vs nGeostrophic
nStokes vs nTotal
Equivalent calculations for skimulator data

@author: Thomas Holding
Modified by Daniel J. Ford (d.ford@exeter.ac.uk)
Date: 03/2023
Changes:
- Added draw_box function (L35-45)
- Added importing of Laurelle areas (L115)
- Updated data loading at (L118-123)
- Added Laurelle boxes and text to each plot (e.g: L117+)
- Updated figure plotting to produce each comparision in one figure using GridSpec
"""

import cPickle as pickle;
import numpy as np;
import matplotlib.pyplot as plt;
import python_util.parameter_sets as ps;
import python_util.skim_utilities as su;
import python_util.mask_functions as mf;
from os import path;
from mpl_toolkits.basemap import Basemap;
from matplotlib.gridspec import GridSpec


def get_month_indices(months, numYears):
    indices = [];
    for year in range(numYears):
        for month in months:
            indices.append(year*12 + month);
    return indices;

def draw_box(mapFig, lon1, lat1, lon2, lat2, color='k'):
    lons = [lon1, lon2];
    lats = [lat1, lat2];
    x, y = mapFig(lons, lats);
#    lines = [((x[0], y[0]), (x[1], y[0])),
#             ((x[0], y[0]), (x[0], y[1])),
#             ((x[1], y[1]), (x[0], y[1])),
#             ((x[1], y[1]), (x[1], y[0]))];
    linesX = [x[0], x[1], x[0], x[0], x[1], x[0], x[1], x[1]];
    linesY = [y[0], y[0], y[0], y[1], y[1], y[1], y[1], y[0]];
    mapFig.plot(linesX, linesY, color=color);


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
    totalVals = np.empty((length, len(indices))); totalVals[:] = np.nan;
    stokesPassMask = np.empty((length, len(indices))); stokesPassMask[:] = np.nan; #Filters where stokes is important?

    #Copy data:
    for i, index in enumerate(indices):
        xcoords[:,i] = [curData.indexX for curData in data[index]];
        ycoords[:,i] = [curData.indexY for curData in data[index]];
        ekmanVals[:,i] = [curData.ekmanProportionGeoEk for curData in data[index]];
        geostrophicVals[:,i] = [curData.geostrophicProportionGeoEk for curData in data[index]];
        stokesVals[:,i] = [curData.stokesProportionOfTotal for curData in data[index]];
        totalVals[:,i] = [curData.totalcurrent for curData in data[index]];
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
    results.totalcurrent = np.nanmean(totalVals,axis=1)
    results.totalcurrentSD = np.nanstd(totalVals,axis=1)

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
    results.totalcurrent = results.totalcurrent[toKeep]
    results.totalcurrentSD = results.totalcurrentSD[toKeep]

    return results;

def plot_laurelle_regions(mapFig):
    laurelle_area = mf.laurelle_area_list
    laurelle_names = mf.laurelle_names
    offsets = [[10,-15],#North Sea
        [10,-15], # English Channel
        [-50,40], # Southern Greenland
        [-15,-35], # Antartic Pennisula
        [10,-15], # Labrador Sea
        [-15,-30], # CoastOfJapan
        [-60,-30], # Cascadian Shelf
        [10,-20], # South Atlatnic Bight
        [10,-15], # Mid Atlantic Bight
        [10,-15], # Barent Sea
        [-45,-30], # Tasmanian Shelf
        [30,0], # Patagonia
        [-10,-20], # Bering Sea
        [10,25]] #Irminger Sea

    for i in range(0,len(laurelle_area)):
        draw_box(mapFig,laurelle_area[i][0][0],laurelle_area[i][0][1],laurelle_area[i][1][0],laurelle_area[i][1][1])
        x,y = mapFig((laurelle_area[i][1][0]+laurelle_area[i][0][0])/2,(laurelle_area[i][0][1]+laurelle_area[i][1][1])/2)
        x2,y2 = (offsets[i][0],offsets[i][1])
        plt.annotate(laurelle_names[i], xy=(x, y),  xycoords='data',
                xytext=(x2, y2), textcoords='offset points',
                arrowprops=dict(arrowstyle="->"), fontsize=10,fontweight='bold'
                )

resampleFraction = 1; #Set to 1 to disable
plot123 = True;
plot789 = True;
plotEkmanVsGeostrophic = True; #abs ekman vs across-shelf geostrophic
plotStokesTotal = True; #across-shelf stokes vs across-shelf  ???
plotTotal = True


params = ps.get_current_params();
if params.paramsetName != "global": #Could actually use the shelf coordinates filename (which is also stored in params)...
    raise ValueError("This plotting script is intended only for the 'global' parameter set. Significant adaptation is required for use with any other datasets that may change the shelf-coordinates.");
if ('allData' in globals()) == False:
    ans = raw_input("Press key to read in 'allData', ctrl+c to cancel...");
    allData = pickle.load(open(path.join("D:/SKIM", "current_data", "surface_currents_"+params.paramsetName+"_500m.p"), "rb"));
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



    ###############n#############################
    # Plot Jan, Feb, Mar nekman vs ngeostrophic #
    ##################n##########################
if plotEkmanVsGeostrophic:
    fig = plt.figure(figsize=(10,10))
    gs = GridSpec(2,1, figure=fig, wspace=0.33,hspace=0.10,bottom=0.05,top=0.97,left=0.05,right=1.05)

    data = get_processed_data([0, 1, 2], allData);
    if resampleFraction > 1:
        resample_data(data, resampleFraction);

    #Converts coordinates into lon/lat.
    lats = []; lons = [];
    for i in range(len(data.xcoords)):
        lon, lat = su.convert_index_to_lonlat(data.xcoords[i], data.ycoords[i], params.pixelRes);
        lats.append(lat); lons.append(lon);

    ax1 = fig.add_subplot(gs[0,0])
    mapFig = Basemap(llcrnrlon=-180.0, llcrnrlat=-90, urcrnrlon=180.0, urcrnrlat=90.0, resolution='l', projection='cyl', lat_0 = 39.5, lon_0 = -3.25);
    mapFig.drawmapboundary(fill_color=(0.85, 0.85, 0.85));
    mapFig.fillcontinents(color=(0.5, 0.5, 0.5), zorder=1);
    mapFig.drawmeridians(np.arange(0, 360, 60), labels=[0,0,0,1], color=(0.3, 0.3, 0.3), fontsize=ticksize);
    mapFig.drawparallels(np.arange(-90, 90, 30), labels=[1,0,0,0], color=(0.3, 0.3, 0.3), fontsize=ticksize);
    mapFig.scatter(lons, lats, latlon=True, c=data.geostrophicProportions, marker='o', cmap=plt.cm.BrBG);
    #lon1, lat1, lon2, lat2
    plot_laurelle_regions(mapFig)
    cbar = plt.colorbar(orientation="vertical", ticks=[0.0, 0.5, 1.0]);
    plt.clim(0, 1);
    cbar.set_label("relative strength of geostrophic versus Ekman current", fontsize=ticksize);
    cbar.ax.set_yticklabels(["0.0", "0.5", "1.0"]);
    cbar.ax.tick_params(labelsize=ticksize);
    plt.title("January, February, March, "+str(params.start_year)+"-"+str(params.end_year), fontsize=ticksize);
    #plt.savefig("plots/current_component_dominance/ekman_vs_geostrophic_123_"+str(params.start_year)+"_"+str(params.end_year)+".pdf");
    #plt.savefig("plots/current_component_dominance/ekman_vs_geostrophic_123_"+str(params.start_year)+"_"+str(params.end_year)+".png",dpi=300);

    data = get_processed_data([6, 7, 8], allData);
    if resampleFraction > 1:
        resample_data(data, resampleFraction);

    #Converts coordinates into lon/lat.
    lats = []; lons = [];
    for i in range(len(data.xcoords)):
        lon, lat = su.convert_index_to_lonlat(data.xcoords[i], data.ycoords[i], params.pixelRes);
        lats.append(lat); lons.append(lon);

    ax2 = fig.add_subplot(gs[1,0])
    mapFig = Basemap(llcrnrlon=-180.0, llcrnrlat=-90, urcrnrlon=180.0, urcrnrlat=90.0, resolution='l', projection='cyl', lat_0 = 39.5, lon_0 = -3.25);
    mapFig.drawmapboundary(fill_color=(0.85, 0.85, 0.85));
    mapFig.fillcontinents(color=(0.5, 0.5, 0.5), zorder=1);
    mapFig.drawmeridians(np.arange(0, 360, 60), labels=[0,0,0,1], color=(0.3, 0.3, 0.3), fontsize=ticksize);
    mapFig.drawparallels(np.arange(-90, 90, 30), labels=[1,0,0,0], color=(0.3, 0.3, 0.3), fontsize=ticksize);
    mapFig.scatter(lons, lats, latlon=True, c=data.geostrophicProportions, marker='o', cmap=plt.cm.BrBG);
    plot_laurelle_regions(mapFig)
    cbar = plt.colorbar(orientation="vertical", ticks=[0.0, 0.5, 1.0]);
    plt.clim(0, 1);
    cbar.set_label("relative strength of geostrophic versus Ekman current", fontsize=ticksize);
    cbar.ax.set_yticklabels(["0.0", "0.5", "1.0"]);
    plt.title("July, August, September, "+str(params.start_year)+"-"+str(params.end_year), fontsize=ticksize);
    #plt.savefig("plots/current_component_dominance/ekman_vs_geostrophic_789_"+str(params.start_year)+"_"+str(params.end_year)+".pdf");
    ax1.text(0.94,0.97,'(a)',transform=ax1.transAxes,va='top',fontsize=18,fontweight='bold')
    ax2.text(0.94,0.97,'(b)',transform=ax2.transAxes,va='top',fontsize=18,fontweight='bold')
    fig.savefig("plots/current_component_dominance/ekman_vs_geostrophic_"+str(params.start_year)+"_"+str(params.end_year)+".png",dpi=300);

    ###############n#############################
    # Plot Jan, Feb, Mar nstokes vs ntotal      #
    ##################n##########################
if plotStokesTotal:
    fig = plt.figure(figsize=(10,10))
    gs = GridSpec(2,1, figure=fig, wspace=0.33,hspace=0.10,bottom=0.05,top=0.97,left=0.05,right=1.05)

    data = get_processed_data([0, 1, 2], allData);
    if resampleFraction > 1:
        resample_data(data, resampleFraction);

    #Converts coordinates into lon/lat.
    lats = []; lons = [];
    for i in range(len(data.xcoords)):
        lon, lat = su.convert_index_to_lonlat(data.xcoords[i], data.ycoords[i], params.pixelRes);
        lats.append(lat); lons.append(lon);


    passed = np.where(data.stokesPassMask != 0.0);
    ax1 = fig.add_subplot(gs[0,0])
    #plt.figure(figsize=figsize);
    mapFig = Basemap(llcrnrlon=-180.0, llcrnrlat=-90, urcrnrlon=180.0, urcrnrlat=90.0, resolution='l', projection='cyl', lat_0 = 39.5, lon_0 = -3.25);
    mapFig.drawmapboundary(fill_color=(0.85, 0.85, 0.85));
    mapFig.fillcontinents(color=(0.5, 0.5, 0.5), zorder=1);
    mapFig.drawmeridians(np.arange(0, 360, 60), labels=[0,0,0,1], color=(0.3, 0.3, 0.3), fontsize=ticksize);
    mapFig.drawparallels(np.arange(-90, 90, 30), labels=[1,0,0,0], color=(0.3, 0.3, 0.3), fontsize=ticksize);
    #mapFig.scatter(np.array(lons), np.array(lats), latlon=True, c=nstokesProportions, marker='o', cmap=plt.cm.RdBu);
    mapFig.scatter(np.array(lons)[passed], np.array(lats)[passed], latlon=True, c=data.stokesProportions[passed], marker='o', cmap=plt.cm.YlOrRd);
    plot_laurelle_regions(mapFig)
    print "a, max:", np.max(data.stokesProportions[passed]);

    #tickToUse = [0.0, 0.5, 1.0];
    ticksToUse = [0.0, 0.33, 0.66];
    cbar = plt.colorbar(orientation="vertical", ticks=ticksToUse);
    plt.clim(ticksToUse[0], ticksToUse[-1]);
    cbar.set_label("relative strength of stokes versus total current", fontsize=ticksize);
    cbar.ax.set_yticklabels([str(v) for v in ticksToUse]);
    plt.title("January, February, March, "+str(params.start_year)+"-"+str(params.end_year), fontsize=ticksize);
    #plt.savefig("plots/current_component_dominance/stokes_vs_total_123_"+str(params.start_year)+"_"+str(params.end_year)+".pdf");
    #plt.savefig("plots/current_component_dominance/stokes_vs_total_123_"+str(params.start_year)+"_"+str(params.end_year)+".png",dpi=300);

    data = get_processed_data([6, 7, 8], allData);
    if resampleFraction > 1:
        resample_data(data, resampleFraction);

    #Converts coordinates into lon/lat.
    lats = []; lons = [];
    for i in range(len(data.xcoords)):
        lon, lat = su.convert_index_to_lonlat(data.xcoords[i], data.ycoords[i], params.pixelRes);
        lats.append(lat); lons.append(lon);
    passed = np.where(data.stokesPassMask != 0.0);

    ax2 = fig.add_subplot(gs[1,0])
    mapFig = Basemap(llcrnrlon=-180.0, llcrnrlat=-90, urcrnrlon=180.0, urcrnrlat=90.0, resolution='l', projection='cyl', lat_0 = 39.5, lon_0 = -3.25);
    mapFig.drawmapboundary(fill_color=(0.85, 0.85, 0.85));
    mapFig.fillcontinents(color=(0.5, 0.5, 0.5), zorder=1);
    mapFig.drawmeridians(np.arange(0, 360, 60), labels=[0,0,0,1], color=(0.3, 0.3, 0.3), fontsize=ticksize);
    mapFig.drawparallels(np.arange(-90, 90, 30), labels=[1,0,0,0], color=(0.3, 0.3, 0.3), fontsize=ticksize);
    #mapFig.scatter(np.array(lons), np.array(lats), latlon=True, c=nstokesProportions, marker='o', cmap=plt.cm.RdBu);
    mapFig.scatter(np.array(lons)[passed], np.array(lats)[passed], latlon=True, c=data.stokesProportions[passed], marker='o', cmap=plt.cm.YlOrRd);
    plot_laurelle_regions(mapFig)
    print "b, max:", np.max(data.stokesProportions[passed]);

    #tickToUse = [0.0, 0.5, 1.0];
    ticksToUse = [0.0, 0.33, 0.66];
    cbar = plt.colorbar(orientation="vertical", ticks=ticksToUse);
    plt.clim(ticksToUse[0], ticksToUse[-1]);
    cbar.set_label("relative strength of stokes versus total current", fontsize=ticksize);
    cbar.ax.set_yticklabels([str(v) for v in ticksToUse]);
    plt.title("July, August, September, "+str(params.start_year)+"-"+str(params.end_year), fontsize=ticksize);
    #plt.savefig("plots/current_component_dominance/stokes_vs_total_789_"+str(params.start_year)+"_"+str(params.end_year)+".pdf");
    ax1.text(0.94,0.97,'(a)',transform=ax1.transAxes,va='top',fontsize=18,fontweight='bold')
    ax2.text(0.94,0.97,'(b)',transform=ax2.transAxes,va='top',fontsize=18,fontweight='bold')
    fig.savefig("plots/current_component_dominance/stokes_vs_total_"+str(params.start_year)+"_"+str(params.end_year)+".png",dpi=300);

if plotTotal:
    fig = plt.figure(figsize=(10,10))
    gs = GridSpec(2,1, figure=fig, wspace=0.33,hspace=0.10,bottom=0.05,top=0.97,left=0.05,right=1.05)
    ####################################
    # Process data for Jan, Feb, March #
    ####################################

    data = get_processed_data([0, 1, 2], allData);
    if resampleFraction > 1:
        resample_data(data, resampleFraction);

    #Converts coordinates into lon/lat.
    lats = []; lons = [];
    for i in range(len(data.xcoords)):
        lon, lat = su.convert_index_to_lonlat(data.xcoords[i], data.ycoords[i], params.pixelRes);
        lats.append(lat); lons.append(lon);


    ax1 = fig.add_subplot(gs[0,0])
    mapFig = Basemap(llcrnrlon=-180.0, llcrnrlat=-90, urcrnrlon=180.0, urcrnrlat=90.0, resolution='l', projection='cyl', lat_0 = 39.5, lon_0 = -3.25);
    mapFig.drawmapboundary(fill_color=(0.85, 0.85, 0.85));
    mapFig.fillcontinents(color=(0.5, 0.5, 0.5), zorder=1);
    mapFig.drawmeridians(np.arange(0, 360, 60), labels=[0,0,0,1], color=(0.3, 0.3, 0.3), fontsize=ticksize);
    mapFig.drawparallels(np.arange(-90, 90, 30), labels=[1,0,0,0], color=(0.3, 0.3, 0.3), fontsize=ticksize);
    mapFig.scatter(lons, lats, latlon=True, c=data.totalcurrent, marker='o', cmap=plt.cm.RdBu);
    #lon1, lat1, lon2, lat2
    plot_laurelle_regions(mapFig)
    cbar = plt.colorbar(orientation="vertical") #ticks=[0.0, 0.5, 1.0]);
    plt.clim(-0.5, 0.5);
    cbar.set_label("Total shelf break current ($ms^{-1}$)", fontsize=ticksize);
    #cbar.ax.set_yticklabels(["0.0", "0.5", "1.0"]);
    cbar.ax.tick_params(labelsize=ticksize);
    plt.title("January, February, March, "+str(params.start_year)+"-"+str(params.end_year), fontsize=ticksize);
    #plt.savefig("plots/current_component_dominance/ekman_vs_geostrophic_123_"+str(params.start_year)+"_"+str(params.end_year)+".pdf");


    data = get_processed_data([6, 7, 8], allData);
    if resampleFraction > 1:
        resample_data(data, resampleFraction);

    #Converts coordinates into lon/lat.
    lats = []; lons = [];
    for i in range(len(data.xcoords)):
        lon, lat = su.convert_index_to_lonlat(data.xcoords[i], data.ycoords[i], params.pixelRes);
        lats.append(lat); lons.append(lon);

    ax2 = fig.add_subplot(gs[1,0])
    mapFig = Basemap(llcrnrlon=-180.0, llcrnrlat=-90, urcrnrlon=180.0, urcrnrlat=90.0, resolution='l', projection='cyl', lat_0 = 39.5, lon_0 = -3.25);
    mapFig.drawmapboundary(fill_color=(0.85, 0.85, 0.85));
    mapFig.fillcontinents(color=(0.5, 0.5, 0.5), zorder=1);
    mapFig.drawmeridians(np.arange(0, 360, 60), labels=[0,0,0,1], color=(0.3, 0.3, 0.3), fontsize=ticksize);
    mapFig.drawparallels(np.arange(-90, 90, 30), labels=[1,0,0,0], color=(0.3, 0.3, 0.3), fontsize=ticksize);
    mapFig.scatter(lons, lats, latlon=True, c=data.totalcurrent, marker='o', cmap=plt.cm.RdBu);
    #lon1, lat1, lon2, lat2
    plot_laurelle_regions(mapFig)
    cbar = plt.colorbar(orientation="vertical") #ticks=[0.0, 0.5, 1.0]);
    plt.clim(-0.5, 0.5);
    cbar.set_label("Total shelf break current ($ms^{-1}$)", fontsize=ticksize);
    #cbar.ax.set_yticklabels(["0.0", "0.5", "1.0"]);
    cbar.ax.tick_params(labelsize=ticksize);
    plt.title("July, August, September, "+str(params.start_year)+"-"+str(params.end_year), fontsize=ticksize);

    ax1.text(0.94,0.97,'(a)',transform=ax1.transAxes,va='top',fontsize=18,fontweight='bold')
    ax2.text(0.94,0.97,'(b)',transform=ax2.transAxes,va='top',fontsize=18,fontweight='bold')
    fig.savefig("plots/current_component_dominance/manuscript.png",dpi=300);


del data
