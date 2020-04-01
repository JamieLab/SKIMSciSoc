#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 21 14:23:19 2018

Spatial plots of the dominance of each current component.
Ekman vs nGeostrophic
nEkman vs nGeostrophic
nStokes vs nTotal

Note that eks is the simulator data ekman current, while 'ekman' is the wavewatch/globcurrent ekman current.

@author: rr
"""

import cPickle as pickle;
import numpy as np;
import matplotlib.pyplot as plt;
import python_util.parameter_sets as ps;
import python_util.skim_utilities as su;
from os import path;
from mpl_toolkits.basemap import Basemap;
import python_util.mask_functions as mf;
from matplotlib import ticker;
from uncertainties import unumpy as unp;
from uncertainties import ufloat;


def get_mean_month_data(months, data):
    length = 0;
    for m in months:
        if len(data[m]) > length:
            length = len(data[m]);
    
    xcoords = np.empty((length, len(months))); xcoords[:] = np.nan;
    ycoords = np.empty((length, len(months))); ycoords[:] = np.nan;

    #skimulator data
    eksOntoShelfVals = np.empty((length, len(months))); eksOntoShelfVals[:] = np.nan;
    altiOntoShelfVals = np.empty((length, len(months))); altiOntoShelfVals[:] = np.nan;
    skimOntoShelfVals = np.empty((length, len(months))); skimOntoShelfVals[:] = np.nan;
    truthOntoShelfVals = np.empty((length, len(months))); truthOntoShelfVals[:] = np.nan;
    
    for m, month in enumerate(months):
        xcoords[:,m] = [cData.indexX for cData in data[month]];
        ycoords[:,m] = [cData.indexY for cData in data[month]];
        
        eksOntoShelfVals[:,m] = [cData.nEksAcrossShelf for cData in data[month]];
        altiOntoShelfVals[:,m] = [cData.nAltiAcrossShelf for cData in data[month]];
        skimOntoShelfVals[:,m] = [cData.nSkimAcrossShelf for cData in data[month]];
        truthOntoShelfVals[:,m] = [cData.nTruthAcrossShelf for cData in data[month]];
    
    return np.nanmean(xcoords, axis=1), np.nanmean(ycoords, axis=1), eksOntoShelfVals, altiOntoShelfVals, skimOntoShelfVals, truthOntoShelfVals;

#Returns the mean and standard deviations for each bar chart (ekman component, geostrophic component, ageostrophic (ref) and ageostrophic (skim))
#Expects uarray objects from the uncertainties.unumpy package (except for mask which is a np.array)
def calc_bar_values(mask, eksOntoShelf, altiOntoShelf, totalOntoShelf):
    #Apply masks
    totalMasked = totalOntoShelf[mask,:];
    ekmanMasked = eksOntoShelf[mask,:];
    altiMasked = altiOntoShelf[mask,:];
    ageostrophicMasked = totalMasked - altiMasked - ekmanMasked;
    
    #Calculate means
    ekmanMean = np.nanmean(ekmanMasked);
    altiMean = np.nanmean(altiMasked);
    ageostrophicMean = np.nanmean(ageostrophicMasked);
    totalMean = np.nanmean(ekmanMasked);

    #Create uncertainties.ufloat objects for each mean
    uekman = ufloat(ekmanMean, np.nanstd(ekmanMasked));
    ualti = ufloat(altiMean, np.nanstd(altiMean));
    uageostrophic = ufloat(ageostrophicMean, np.nanstd(ageostrophicMean));
    utotal = ufloat(totalMean, np.nanstd(totalMasked));

    #Calculate percentage of total
    ekmanPercent = uekman / utotal * 100.0;
    altiPercent = ualti / utotal * 100.0;
    ageostrophicPercent = uageostrophic / utotal * 100.0;
    
    #Extract vals and sds seperately
    vals = [item.n for item in [ekmanPercent, altiPercent, ageostrophicPercent]];
    sds = [item.s for item in [ekmanPercent, altiPercent, ageostrophicPercent]];
    
    return vals, sds;


def add_stddev_label(bbr, xvals, barValsRef, barSDsRef):
    colourEkman="#2ca02c";
    colourGeostrophic="#1f77b4";
    colourStokes="#ff7f0e";
    colours = (colourEkman, colourGeostrophic, colourStokes);
    
    for i in range(len(xvals)):
        ypos = barValsRef[i];
        if (ypos >= 0):
            ypos += 30.0;
        else:
            ypos -= 80.0;
        bbr.text(xvals[i], ypos, "$\pm$"+str(np.around(barSDsRef[i], 2)), weight="bold", color=colours[i], fontdict={"fontsize": 8, "horizontalalignment": "center", "alpha": 0.7});


monthIndices123 = [3, 4, 5]; #Data is not necessarily in order.
monthIndices789 = [9, 10, 11];

plot123 = True;
plot789 = True;
totalMap = False;
plotbars = True;

#Map regions
norwegianShelfMask = np.array( [False]*267); #Currently off the map
northSeaMask = np.array( [False]*267); northSeaMask[0:55] = True;
celticSeaMask = np.array( [False]*267); celticSeaMask[55:150] = True;
bayOfBiscayMask = np.array( [False]*267); bayOfBiscayMask[150:215] = True;
spanishCoastMask = np.array( [False]*267); spanishCoastMask[215:267] = True;

plotSavePath = "plots/european_shelf_case_study/";

#Map extents
llcrnrlon=-18.0;
llcrnrlat=39.0;
urcrnrlon=12.0;
urcrnrlat=65.0;


params = ps.get_European_shelf_params();
if params.paramsetName != "europeanshelf":
    raise ValueError("This plotting script is intended only for the 'europeanshelf' parameter set. Significant adaptation is required for use with any other datasets that may change the shelf-coordinates.");
data = pickle.load(open(path.join("data/europeanshelf/current_data/surface_currents_"+params.paramsetName+".p"), "rb"));
for monthIndex in range(len(data)):
    print "Applying mask to slice", monthIndex;
    data[monthIndex] = mf.generic_area_mask([mf.area_EuropeanShelf], data[monthIndex], params);
    #data[monthIndex] = params.contourMaskFunc(data[monthIndex], params);

figsize = (8.5*1.3,6.5*1.3);
ticksize = 16;
ticksizeBar = 12;
barPlotYLower = -400;
barPlotYUpper = 900;


####################################
# Process data for Jan, Feb, March #
####################################
if plot123:
    xcoords, ycoords, eksOntoShelfVals, altiOntoShelfVals, skimOntoShelfVals, truthOntoShelfVals = get_mean_month_data(monthIndices123, data);
    eksOntoShelfLong = np.nanmean(eksOntoShelfVals, axis=1);
    altiOntoShelfLong = np.nanmean(altiOntoShelfVals, axis=1);
    skimOntoShelfLong = np.nanmean(skimOntoShelfVals, axis=1);
    truthOntoShelfLong = np.nanmean(truthOntoShelfVals, axis=1);

    #remove all nans
    toKeep = (unp.isnan(eksOntoShelfLong)==False) | (unp.isnan(altiOntoShelfLong)==False) | (unp.isnan(skimOntoShelfLong)==False) | (unp.isnan(truthOntoShelfLong)==False);
    xcoords = xcoords[toKeep];
    ycoords = ycoords[toKeep];
    eksOntoShelfLong = eksOntoShelfLong[toKeep];
    altiOntoShelfLong = altiOntoShelfLong[toKeep];
    skimOntoShelfLong = skimOntoShelfLong[toKeep];
    truthOntoShelfLong = truthOntoShelfLong[toKeep];
    
    #Converts coordinates into lon/lat.
    lats = []; lons = [];
    for i in range(len(xcoords)):
        lon, lat = su.convert_index_to_lonlat(xcoords[i], ycoords[i], params.pixelRes, lon0=params.originLon, lat0=params.originLat);
        lats.append(lat); lons.append(lon);
    
    
    ###############n#############################
    # Plot Jan, Feb, Mar nekman vs ngeostrophic #
    ##################n##########################
    if totalMap:
        climMax = 0.4;
        climMin = -0.15;
        plt.figure(figsize=figsize);
        mapFig = Basemap(llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat, urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat, resolution='l', projection='cyl', lat_0 = 39.5, lon_0 = -3.25);
        mapFig.drawmapboundary(fill_color=(0.85, 0.85, 0.85));
        mapFig.fillcontinents(color=(0.5, 0.5, 0.5), zorder=1);
        mapFig.drawmeridians(np.arange(0, 360, 5), labels=[0,0,0,1], color=(0.3, 0.3, 0.3), fontsize=ticksize);
        mapFig.drawparallels(np.arange(-90, 90, 5), labels=[1,0,0,0], color=(0.3, 0.3, 0.3), fontsize=ticksize);
        
#        #TEMP FINDING CODE
#        maskMin = 150;
#        maskMax = 215; #267
#        v = unp.nominal_values(truthOntoShelf);#Long);
#        v[:] = 0.0;
#        v[maskMin:maskMax]=1.0;
#        mapFig.scatter(lons, lats, latlon=True, c=v, marker='o', cmap=plt.cm.YlOrRd);
        
        mapFig.scatter(lons, lats, latlon=True, c=unp.nominal_values(truthOntoShelfLong), marker='o', cmap=plt.cm.YlOrRd);
        cbar = plt.colorbar(orientation="horizontal", ticks=[0.0, 0.5, 1.0]);
        tick_locator = ticker.MaxNLocator(nbins=5);
        cbar.locator = tick_locator;
        cbar.update_ticks();
        plt.clim(climMin, climMax);
        cbar.set_label("onto shelf current (m s$^{-1})$", fontsize=ticksize);
        plt.title("January, February, March", fontsize=ticksize);
        plt.savefig("plots/european_shelf_case_study/total_onto_shelf_123_ref_european.pdf");
        
        
        #Now with skim
        plt.figure(figsize=figsize);
        mapFig = Basemap(llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat, urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat, resolution='l', projection='cyl', lat_0 = 39.5, lon_0 = -3.25);
        mapFig.drawmapboundary(fill_color=(0.85, 0.85, 0.85));
        mapFig.fillcontinents(color=(0.5, 0.5, 0.5), zorder=1);
        mapFig.drawmeridians(np.arange(0, 360, 5), labels=[0,0,0,1], color=(0.3, 0.3, 0.3), fontsize=ticksize);
        mapFig.drawparallels(np.arange(-90, 90, 5), labels=[1,0,0,0], color=(0.3, 0.3, 0.3), fontsize=ticksize);

        mapFig.scatter(lons, lats, latlon=True, c=unp.nominal_values(skimOntoShelfLong), marker='o', cmap=plt.cm.YlOrRd);
        cbar = plt.colorbar(orientation="horizontal", ticks=[0.0, 0.5, 1.0]);
        tick_locator = ticker.MaxNLocator(nbins=5);
        cbar.locator = tick_locator;
        cbar.update_ticks();
        plt.clim(climMin, climMax);
        cbar.set_label("onto shelf current (m s$^{-1})$", fontsize=ticksize);
        plt.title("January, February, March", fontsize=ticksize);
        plt.savefig("plots/european_shelf_case_study/total_onto_shelf_123_skim_european.pdf");


    if plotbars:
        colours = [(0.3, 0.7, 0.3), (0.30, 0.30, 0.6), (0.9, 0.4, 0.0)];
        xvals = [0, 1, 2];
        tickLabels = ["", "", ""]; #["ekman", "geostrophic", "ageostrophic"];
        alpha = 0.5;
        
        #plot reference values
        fig, ((nsr, csr), (bbr, scr)) = plt.subplots(nrows=2, ncols=2, sharex=True);
        nsr.set_ylabel("percentage of total current (%)", fontsize=ticksizeBar);
        
        #Bay of buscay region
        barValsRef, barSDsRef = calc_bar_values(bayOfBiscayMask, eksOntoShelfVals, altiOntoShelfVals, truthOntoShelfVals);
        bbr.bar(xvals, barValsRef, color=colours, tick_label=tickLabels, alpha=alpha);
        #bbr.errorbar(xvals, barValsRef, yerr=barSDsRef, fmt='.', color='k', alpha=alpha);
        add_stddev_label(bbr, xvals, barValsRef, barSDsRef);
        
        #plt.xticks(fontsize=ticksizeBar, rotation=45); plt.yticks(fontsize=ticksizeBar);
        bbr.set_ylim(barPlotYLower, barPlotYUpper);
        
        #North Sea region
        barValsRef, barSDsRef = calc_bar_values(northSeaMask, eksOntoShelfVals, altiOntoShelfVals, truthOntoShelfVals);
        nsr.bar(xvals, barValsRef, color=colours, tick_label=tickLabels, alpha=alpha);
        #nsr.errorbar(xvals, barValsRef, yerr=barSDsRef, fmt='.', color='k', alpha=alpha);
        add_stddev_label(nsr, xvals, barValsRef, barSDsRef);
        #plt.xticks(fontsize=ticksizeBar, rotation=45); plt.yticks(fontsize=ticksizeBar);
        nsr.set_ylim(barPlotYLower, barPlotYUpper);
        
        #Celtic Sea region
        barValsRef, barSDsRef = calc_bar_values(celticSeaMask, eksOntoShelfVals, altiOntoShelfVals, truthOntoShelfVals);
        csr.bar(xvals, barValsRef, color=colours, tick_label=tickLabels, alpha=alpha);
        #csr.errorbar(xvals, barValsRef, yerr=barSDsRef, fmt='.', color='k', alpha=alpha);
        add_stddev_label(csr, xvals, barValsRef, barSDsRef);
        csr.set_ylim(barPlotYLower, barPlotYUpper);
        
        #Spanish Coast region
        barValsRef, barSDsRef = calc_bar_values(spanishCoastMask, eksOntoShelfVals, altiOntoShelfVals, truthOntoShelfVals);
        scr.bar(xvals, barValsRef, color=colours, tick_label=tickLabels, alpha=alpha);
        #scr.errorbar(xvals, barValsRef, yerr=barSDsRef, fmt='.', color='k', alpha=alpha);
        add_stddev_label(scr, xvals, barValsRef, barSDsRef);
        scr.set_ylim(barPlotYLower, barPlotYUpper);
        
        plt.tight_layout();
        fig.savefig(path.join(plotSavePath, "onto_shelf_current_123_ref.pdf"));
        
        
        #plot skim values
        fig, ((nsr, csr), (bbr, scr)) = plt.subplots(nrows=2, ncols=2, sharex=False);
        nsr.set_ylabel("percentage of total current (%)", fontsize=ticksizeBar);
        
        #Norwegian Shelf region
        barValsRef, barSDsRef = calc_bar_values(bayOfBiscayMask, eksOntoShelfVals, altiOntoShelfVals, truthOntoShelfVals);
        barValsSkim, barSDsSkim = calc_bar_values(bayOfBiscayMask, eksOntoShelfVals, altiOntoShelfVals, skimOntoShelfVals);
        bbr.bar(xvals, barValsRef, color=colours, tick_label=tickLabels, alpha=alpha);
        add_stddev_label(bbr, xvals, barValsRef, barSDsRef);
        #bbr.errorbar(xvals, barValsSkim, yerr=barSDsSkim, fmt='.', color='k', alpha=alpha);
        bbr.hlines(barValsSkim, np.array(xvals)-0.2, np.array(xvals)+0.2);
        bbr.set_ylim(barPlotYLower, barPlotYUpper);
        
        #North Sea region
        barValsRef, barSDsRef = calc_bar_values(northSeaMask, eksOntoShelfVals, altiOntoShelfVals, truthOntoShelfVals);
        barValsSkim, barSDsSkim = calc_bar_values(northSeaMask, eksOntoShelfVals, altiOntoShelfVals, skimOntoShelfVals);
        nsr.bar(xvals, barValsRef, color=colours, tick_label=tickLabels, alpha=alpha);
        #nsr.errorbar(xvals, barValsSkim, yerr=barSDsSkim, fmt='.', color='k', alpha=alpha);
        add_stddev_label(nsr, xvals, barValsRef, barSDsRef);
        nsr.hlines(barValsSkim, np.array(xvals)-0.2, np.array(xvals)+0.2);
        nsr.set_ylim(barPlotYLower, barPlotYUpper);
        
        #Celtic Sea region
        barValsRef, barSDsRef = calc_bar_values(celticSeaMask, eksOntoShelfVals, altiOntoShelfVals, truthOntoShelfVals);
        barValsSkim, barSDsSkim = calc_bar_values(celticSeaMask, eksOntoShelfVals, altiOntoShelfVals, skimOntoShelfVals);
        csr.bar(xvals, barValsRef, color=colours, tick_label=tickLabels, alpha=alpha);
        #csr.errorbar(xvals, barValsSkim, yerr=barSDsSkim, fmt='.', color='k', alpha=alpha);
        add_stddev_label(csr, xvals, barValsRef, barSDsRef);
        csr.hlines(barValsSkim, np.array(xvals)-0.2, np.array(xvals)+0.2);
        csr.set_ylim(barPlotYLower, barPlotYUpper);
        
        #Spanish Coast region
        barValsRef, barSDsRef = calc_bar_values(spanishCoastMask, eksOntoShelfVals, altiOntoShelfVals, truthOntoShelfVals);
        barValsSkim, barSDsSkim = calc_bar_values(spanishCoastMask, eksOntoShelfVals, altiOntoShelfVals, skimOntoShelfVals);
        scr.bar(xvals, barValsRef, color=colours, tick_label=tickLabels, alpha=alpha);
        #scr.errorbar(xvals, barValsSkim, yerr=barSDsSkim, fmt='.', color='k', alpha=alpha);
        add_stddev_label(scr, xvals, barValsRef, barSDsRef);
        scr.hlines(barValsSkim, np.array(xvals)-0.2, np.array(xvals)+0.2);
        scr.set_ylim(barPlotYLower, barPlotYUpper);
        
        plt.tight_layout();
        plt.savefig(path.join(plotSavePath, "onto_shelf_current_123_skim.pdf"));


##################################
# Process data for Jul, Aug, Sep #
##################################
if plot789:
    xcoords, ycoords, eksOntoShelfVals, altiOntoShelfVals, skimOntoShelfVals, truthOntoShelfVals = get_mean_month_data(monthIndices789, data);
    eksOntoShelfLong = np.nanmean(eksOntoShelfVals, axis=1);
    altiOntoShelfLong = np.nanmean(altiOntoShelfVals, axis=1);
    skimOntoShelfLong = np.nanmean(skimOntoShelfVals, axis=1);
    truthOntoShelfLong = np.nanmean(truthOntoShelfVals, axis=1);

    #remove all nans
    toKeep = (unp.isnan(eksOntoShelfLong)==False) | (unp.isnan(altiOntoShelfLong)==False) | (unp.isnan(skimOntoShelfLong)==False) | (unp.isnan(truthOntoShelfLong)==False);
    xcoords = xcoords[toKeep];
    ycoords = ycoords[toKeep];
    eksOntoShelfLong = eksOntoShelfLong[toKeep];
    altiOntoShelfLong = altiOntoShelfLong[toKeep];
    skimOntoShelfLong = skimOntoShelfLong[toKeep];
    truthOntoShelfLong = truthOntoShelfLong[toKeep];
    
    #Converts coordinates into lon/lat.
    lats = []; lons = [];
    for i in range(len(xcoords)):
        lon, lat = su.convert_index_to_lonlat(xcoords[i], ycoords[i], params.pixelRes, lon0=params.originLon, lat0=params.originLat);
        lats.append(lat); lons.append(lon);
    
    
    ###############n#############################
    # Plot Jan, Feb, Mar nekman vs ngeostrophic #
    ##################n##########################
    if totalMap:
        climMax = 0.4;
        climMin = -0.15;
        plt.figure(figsize=figsize);
        mapFig = Basemap(llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat, urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat, resolution='l', projection='cyl', lat_0 = 39.5, lon_0 = -3.25);
        mapFig.drawmapboundary(fill_color=(0.85, 0.85, 0.85));
        mapFig.fillcontinents(color=(0.5, 0.5, 0.5), zorder=1);
        mapFig.drawmeridians(np.arange(0, 360, 5), labels=[0,0,0,1], color=(0.3, 0.3, 0.3), fontsize=ticksize);
        mapFig.drawparallels(np.arange(-90, 90, 5), labels=[1,0,0,0], color=(0.3, 0.3, 0.3), fontsize=ticksize);
        
#        #TEMP FINDING CODE
#        maskMin = 215;
#        maskMax = 267; #267
#        v = unp.nominal_values(truthOntoShelf);
#        v[:] = 0.0;
#        v[maskMin:maskMax]=1.0;
#        mapFig.scatter(lons, lats, latlon=True, c=v, marker='o', cmap=plt.cm.YlOrRd);
        
        mapFig.scatter(lons, lats, latlon=True, c=unp.nominal_values(truthOntoShelfLong), marker='o', cmap=plt.cm.YlOrRd);
        cbar = plt.colorbar(orientation="horizontal", ticks=[0.0, 0.5, 1.0]);
        tick_locator = ticker.MaxNLocator(nbins=5);
        cbar.locator = tick_locator;
        cbar.update_ticks();
        plt.clim(climMin, climMax);
        cbar.set_label("onto shelf current (m s$^{-1})$", fontsize=ticksize);
        plt.title("July, August, September", fontsize=ticksize);
        plt.savefig("plots/european_shelf_case_study/total_onto_shelf_789_ref_european.pdf");
        
        
        #Now with skim
        plt.figure(figsize=figsize);
        mapFig = Basemap(llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat, urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat, resolution='l', projection='cyl', lat_0 = 39.5, lon_0 = -3.25);
        mapFig.drawmapboundary(fill_color=(0.85, 0.85, 0.85));
        mapFig.fillcontinents(color=(0.5, 0.5, 0.5), zorder=1);
        mapFig.drawmeridians(np.arange(0, 360, 5), labels=[0,0,0,1], color=(0.3, 0.3, 0.3), fontsize=ticksize);
        mapFig.drawparallels(np.arange(-90, 90, 5), labels=[1,0,0,0], color=(0.3, 0.3, 0.3), fontsize=ticksize);

        mapFig.scatter(lons, lats, latlon=True, c=unp.nominal_values(skimOntoShelfLong), marker='o', cmap=plt.cm.YlOrRd);
        cbar = plt.colorbar(orientation="horizontal", ticks=[0.0, 0.5, 1.0]);
        tick_locator = ticker.MaxNLocator(nbins=5);
        cbar.locator = tick_locator;
        cbar.update_ticks();
        plt.clim(climMin, climMax);
        cbar.set_label("onto shelf current (m s$^{-1})$", fontsize=ticksize);
        plt.title("July, August, September", fontsize=ticksize);
        plt.savefig("plots/european_shelf_case_study/total_onto_shelf_789_skim_european.pdf");


    if plotbars:
        colours = [(0.3, 0.7, 0.3), (0.30, 0.30, 0.6), (0.9, 0.4, 0.0)];
        xvals = [0, 1, 2];
        tickLabels = ["", "", ""]; #["ekman", "geostrophic", "ageostrophic"];
        alpha = 0.5;
        
        #plot reference values
        fig, ((nsr, csr), (bbr, scr)) = plt.subplots(nrows=2, ncols=2, sharex=True);
        nsr.set_ylabel("percentage of total current (%)", fontsize=ticksizeBar);
        
        #Bay of buscay region
        barValsRef, barSDsRef = calc_bar_values(bayOfBiscayMask, eksOntoShelfVals, altiOntoShelfVals, truthOntoShelfVals);
        bbr.bar(xvals, barValsRef, color=colours, tick_label=tickLabels, alpha=alpha);
        add_stddev_label(bbr, xvals, barValsRef, barSDsRef);
        #bbr.errorbar(xvals, barValsRef, yerr=barSDsRef, fmt='.', color='k', alpha=alpha);
        #plt.xticks(fontsize=ticksizeBar, rotation=45); plt.yticks(fontsize=ticksizeBar);
        bbr.set_ylim(barPlotYLower, barPlotYUpper);
        
        #North Sea region
        barValsRef, barSDsRef = calc_bar_values(northSeaMask, eksOntoShelfVals, altiOntoShelfVals, truthOntoShelfVals);
        nsr.bar(xvals, barValsRef, color=colours, tick_label=tickLabels, alpha=alpha);
        #nsr.errorbar(xvals, barValsRef, yerr=barSDsRef, fmt='.', color='k', alpha=alpha);
        add_stddev_label(nsr, xvals, barValsRef, barSDsRef);
        #plt.xticks(fontsize=ticksizeBar, rotation=45); plt.yticks(fontsize=ticksizeBar);
        nsr.set_ylim(barPlotYLower, barPlotYUpper);
        
        #Celtic Sea region
        barValsRef, barSDsRef = calc_bar_values(celticSeaMask, eksOntoShelfVals, altiOntoShelfVals, truthOntoShelfVals);
        csr.bar(xvals, barValsRef, color=colours, tick_label=tickLabels, alpha=alpha);
        #csr.errorbar(xvals, barValsRef, yerr=barSDsRef, fmt='.', color='k', alpha=alpha);
        add_stddev_label(csr, xvals, barValsRef, barSDsRef);
        csr.set_ylim(barPlotYLower, barPlotYUpper);
    
        #Spanish Coast region
        barValsRef, barSDsRef = calc_bar_values(spanishCoastMask, eksOntoShelfVals, altiOntoShelfVals, truthOntoShelfVals);
        scr.bar(xvals, barValsRef, color=colours, tick_label=tickLabels, alpha=alpha);
        #scr.errorbar(xvals, barValsRef, yerr=barSDsRef, fmt='.', color='k', alpha=alpha);
        add_stddev_label(scr, xvals, barValsRef, barSDsRef);
        scr.set_ylim(barPlotYLower, barPlotYUpper);
        
        plt.tight_layout();
        fig.savefig(path.join(plotSavePath, "onto_shelf_current_789_ref.pdf"));
        
        
        #plot skim values
        fig, ((nsr, csr), (bbr, scr)) = plt.subplots(nrows=2, ncols=2, sharex=False);
        nsr.set_ylabel("percentage of total current (%)", fontsize=ticksizeBar);
        
        #Norwegian Shelf region
        barValsRef, barSDsRef = calc_bar_values(bayOfBiscayMask, eksOntoShelfVals, altiOntoShelfVals, truthOntoShelfVals);
        barValsSkim, barSDsSkim = calc_bar_values(bayOfBiscayMask, eksOntoShelfVals, altiOntoShelfVals, skimOntoShelfVals);
        bbr.bar(xvals, barValsRef, color=colours, tick_label=tickLabels, alpha=alpha);
        #bbr.errorbar(xvals, barValsSkim, yerr=barSDsSkim, fmt='.', color='k', alpha=alpha);
        add_stddev_label(bbr, xvals, barValsRef, barSDsRef);
        bbr.hlines(barValsSkim, np.array(xvals)-0.2, np.array(xvals)+0.2);
        bbr.set_ylim(barPlotYLower, barPlotYUpper);
        
        #North Sea region
        barValsRef, barSDsRef = calc_bar_values(northSeaMask, eksOntoShelfVals, altiOntoShelfVals, truthOntoShelfVals);
        barValsSkim, barSDsSkim = calc_bar_values(northSeaMask, eksOntoShelfVals, altiOntoShelfVals, skimOntoShelfVals);
        nsr.bar(xvals, barValsRef, color=colours, tick_label=tickLabels, alpha=alpha);
        #nsr.errorbar(xvals, barValsSkim, yerr=barSDsSkim, fmt='.', color='k', alpha=alpha);
        add_stddev_label(nsr, xvals, barValsRef, barSDsRef);
        nsr.hlines(barValsSkim, np.array(xvals)-0.2, np.array(xvals)+0.2);
        nsr.set_ylim(barPlotYLower, barPlotYUpper);
        
        #Celtic Sea region
        barValsRef, barSDsRef = calc_bar_values(celticSeaMask, eksOntoShelfVals, altiOntoShelfVals, truthOntoShelfVals);
        barValsSkim, barSDsSkim = calc_bar_values(celticSeaMask, eksOntoShelfVals, altiOntoShelfVals, skimOntoShelfVals);
        csr.bar(xvals, barValsRef, color=colours, tick_label=tickLabels, alpha=alpha);
        #csr.errorbar(xvals, barValsSkim, yerr=barSDsSkim, fmt='.', color='k', alpha=alpha);
        add_stddev_label(csr, xvals, barValsRef, barSDsRef);
        csr.hlines(barValsSkim, np.array(xvals)-0.2, np.array(xvals)+0.2);
        csr.set_ylim(barPlotYLower, barPlotYUpper);
        
        #Spanish Coast region
        barValsRef, barSDsRef = calc_bar_values(spanishCoastMask, eksOntoShelfVals, altiOntoShelfVals, truthOntoShelfVals);
        barValsSkim, barSDsSkim = calc_bar_values(spanishCoastMask, eksOntoShelfVals, altiOntoShelfVals, skimOntoShelfVals);
        scr.bar(xvals, barValsRef, color=colours, tick_label=tickLabels, alpha=alpha);
        #scr.errorbar(xvals, barValsSkim, yerr=barSDsSkim, fmt='.', color='k', alpha=alpha);
        add_stddev_label(scr, xvals, barValsRef, barSDsRef);
        scr.hlines(barValsSkim, np.array(xvals)-0.2, np.array(xvals)+0.2);
        scr.set_ylim(barPlotYLower, barPlotYUpper);
        
        plt.tight_layout();
        plt.savefig(path.join(plotSavePath, "onto_shelf_current_789_skim.pdf"));





#
#
#bathy = Dataset("/home/rr/Files/Tasks/20180914_SKIM/data/GEBCO_bathymetry_0.25x0.25deg.nc", 'r');
#depth = np.flipud(bathy.variables["mean_depth"][:]);
#depth = su.apply_mask(depth, None, params.ilatRange, params.ilonRange);
#
#plt.figure();
#plt.imshow(depth);
#plt.scatter(coordsX, coordsY, color=colours);
#plt.title("Relative dominance: n(Geostrophic) (red) Ekman (blue)");
#
#
##n(ekman+45) vs n(geostrophic)
#nEkmanProportions = [cData.nEkmanProportion for cData in data[monthIndex]];
#nnGeostrophicProportions = [cData.nnGeostrophicProportion for cData in data[monthIndex]];
#colours = [];
#for i in range(len(coordsX)):
#    colours.append( (nnGeostrophicProportions[i], 0.0, nEkmanProportions[i]) );
#plt.figure();
#plt.imshow(depth);
#plt.scatter(coordsX, coordsY, color=colours);
#plt.title("Relative dominance: n(Geostrophic) (red) n(Ekman+45) (blue)");
#
