#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
@author: Daniel J. Ford (d.ford@exeter.ac.uk)
Created: April 2023


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
import matplotlib.transforms
font = {'weight' : 'normal',
        'size'   : 19}
matplotlib.rc('font', **font)

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

def plot_laurelle_regions(mapFig,ax):
    laurelle_area = mf.laurelle_area_list
    laurelle_names = mf.laurelle_names
    offsets = [[10,-15],#North Sea
        [10,-15], # English Channel
        [-50,50], # Southern Greenland
        [-15,-35], # Antartic Pennisula
        [10,-20], # Labrador Sea
        [-15,-30], # CoastOfJapan
        [-60,-30], # Cascadian Shelf
        [10,-25], # South Atlatnic Bight
        [10,-20], # Mid Atlantic Bight
        [10,-15], # Barent Sea
        [-45,-30], # Tasmanian Shelf
        [30,0], # Patagonia
        [-10,-20], # Bering Sea
        [10,25]] #Irminger Sea

    for i in range(0,len(laurelle_area)):
        draw_box(mapFig,laurelle_area[i][0][0],laurelle_area[i][0][1],laurelle_area[i][1][0],laurelle_area[i][1][1])
        x,y = mapFig((laurelle_area[i][1][0]+laurelle_area[i][0][0])/2,(laurelle_area[i][0][1]+laurelle_area[i][1][1])/2)
        x2,y2 = (offsets[i][0],offsets[i][1])
        ax.annotate(laurelle_names[i], xy=(x, y),  xycoords='data',
                xytext=(x2, y2), textcoords='offset points',
                arrowprops=dict(arrowstyle="->"), fontsize=14,fontweight='bold'
                )
def plot_function(ax,allData,let):
    data = get_processed_data([0, 1, 2], allData);
    if resampleFraction > 1:
        resample_data(data, resampleFraction);

    #Converts coordinates into lon/lat.
    lats = []; lons = [];
    for i in range(len(data.xcoords)):
        lon, lat = su.convert_index_to_lonlat(data.xcoords[i], data.ycoords[i], params.pixelRes);
        lats.append(lat); lons.append(lon);

    mapFig = Basemap(llcrnrlon=-180.0, llcrnrlat=-90, urcrnrlon=180.0, urcrnrlat=90.0, resolution='l', projection='cyl', lat_0 = 39.5, lon_0 = -3.25,ax=ax);
    mapFig.drawmapboundary(fill_color=(0.85, 0.85, 0.85));
    mapFig.fillcontinents(color=(0.5, 0.5, 0.5), zorder=1);
    mapFig.drawmeridians(np.arange(0, 360, 60), labels=[0,0,0,1], color=(0.3, 0.3, 0.3), fontsize=ticksize);
    mapFig.drawparallels(np.arange(-90, 90, 30), labels=[1,0,0,0], color=(0.3, 0.3, 0.3), fontsize=ticksize);
    a = mapFig.scatter(lons, lats, latlon=True, c=data.totalcurrent, marker='o', cmap=plt.cm.RdBu,vmin=-0.5,vmax=0.5);
    #lon1, lat1, lon2, lat2
    plot_laurelle_regions(mapFig,ax)
    cbar = plt.colorbar(a,orientation="vertical",ax=ax) #ticks=[0.0, 0.5, 1.0]);
    #ax.clim(-0.5, 0.5);
    cbar.set_label("Total shelf break current ($ms^{-1}$)", fontsize=ticksize);
    #cbar.ax.set_yticklabels(["0.0", "0.5", "1.0"]);
    cbar.ax.tick_params(labelsize=ticksize);
    ax.text(0.90,0.95,'('+let+')',transform=ax.transAxes,va='top',fontsize=24,fontweight='bold')
    #ax.set_title("January, February, March, "+str(params.start_year)+"-"+str(params.end_year), fontsize=ticksize);
    #plt.savefig("plots/current_component_dominance/ekman_vs_geostrophic_123_"+str(params.start_year)+"_"+str(params.end_year)+".pdf");

ticksize = 19
resampleFraction = 1
params = ps.get_current_params();
#allData = pickle.load(open(path.join("D:/SKIM", "current_data", "surface_currents_"+params.paramsetName+"_500m.p"), "rb"));

fig = plt.figure(figsize=(30,20))
gs = GridSpec(3,2, figure=fig, wspace=0.01,hspace=0.10,bottom=0.05,top=0.97,left=0.05,right=1.00)
ax1 = fig.add_subplot(gs[0,0]); ax2 = fig.add_subplot(gs[0,1]); ax3 = fig.add_subplot(gs[1,0]);ax4 = fig.add_subplot(gs[1,1]);ax5 = fig.add_subplot(gs[2,0]);ax6 = fig.add_subplot(gs[2,1]);
ax = [ax1,ax2,ax3,ax4,ax5,ax6]
l = [300,400,500,600,700,800]
let = ['a','b','c','d','e','f']
for lv in range(len(l)):
    allData = pickle.load(open(path.join("D:/SKIM", "current_data", "surface_currents_"+params.paramsetName+"_" + str(l[lv])+"m.p"), "rb"));
    plot_function(ax[lv],allData,let[lv])
fig.savefig('plots/current_component_dominance/supplementary_fig_s2.png',dpi=300)
