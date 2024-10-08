#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 13:43:32 2018

@author: rr

Modified by Daniel J. Ford (d.ford@exeter.ac.uk)
Date: 03/2023
Changes:
- Added final plotting for mansuscript figure (L436+)
- Added GridSpac to packages
- Changed output figures to .png from .pdf
- Modified data loading section (L201-207)
- Added axes as a optional arguement to plot_shelf_section, and updated figure/axes definitions (L119-123)
"""

import cPickle as pickle;
import numpy as np;
import matplotlib.pyplot as plt;
import python_util.parameter_sets as ps;
import python_util.skim_utilities as su;
import python_util.mask_functions as mf;
from matplotlib.gridspec import GridSpec
from netCDF4 import Dataset;
from os import path;
import sys

import pandas as pd;
import matplotlib.transforms
from mpl_toolkits.basemap import Basemap;

def unweighted_stats(x,y,meth):
    f = np.squeeze(np.argwhere( (np.isnan(x) == 0) & (np.isnan(y) == 0) ))
    #print(f.size)
    if f.size >= 1:
        #print(f)
        x = x[f]
        y = y[f]
        weights = np.ones((f.size))
        if f.size == 1:
            weights =weights[0]
        #print(x)
        #print(weights)
        su = np.sum(weights)
        weights = weights/su
        #print(weights)
        bi = rel_bias(y,x,weights)
        bi_m = med_rel_bias(y,x,weights)
        abi = abs_bias(y,x,weights)
        rms = rmsd(y,x,weights)
        ap = apd(y,x,weights)
        rp = rpd(y,x,weights)
        slope = np.nan
        intercept = np.nan

        struct = {'meth':meth,'rel_bias':bi,'med_rel_bias':bi_m,'abs_bias':abi,'rmsd':rms,'slope':slope,'intercept':intercept,'apd':ap,'rpd':rp,'n':f.size}
    else:
        struct = {'meth':meth,'rel_bias':np.nan,'med_rel_bias':np.nan,'abs_bias':np.nan,'rmsd':np.nan,'slope':np.nan,'intercept':np.nan,'apd':np.nan,'rpd':np.nan,'n':f.size}
    return struct
def unweight(x,y,ax,unit = '$\mu$atm'):
    """
    Function to calculate the unweighted statistics and add them to a scatter plot (well any plot really, in the bottom right corner)
    """
    stats_un = unweighted_stats(x,y,'val')
    #h2 = ax.plot(c,c*stats_un['slope']+stats_un['intercept'],'k-.',zorder=5, label = 'Unweighted')
    rmsd = '%.2f' %np.round(stats_un['rmsd'],2)
    bias = '%.2f' %np.round(stats_un['rel_bias'],2)
    sl = '%.2f' %np.round(stats_un['slope'],2)
    ip = '%.2f' %np.round(stats_un['intercept'],2)
    n = stats_un['n']
    ax.text(0.80,0.15,'RMSD = '+str(rmsd)+ unit+'\nBias = '+bias+unit+'\nN = '+str(n),transform=ax.transAxes,va='top')

def rel_bias(x,y,weights):
    bi = np.average((x)-(y), weights=weights)
    return bi

def med_rel_bias(x,y,weights):
    bi = np.median((x)-(y))
    return bi

def abs_bias(x,y,weights):
    bi = np.average(np.abs(x)-np.abs(y), weights=weights)
    return bi

def rmsd(x,y,weights):
    rms = np.sqrt(np.average((x-y)**2, weights=weights))
    return rms

def apd(x,y,weights):
    ap =  np.average( np.abs(x-y) / (np.abs(np.abs(x)+np.abs(y))/2), weights=weights) *100
    return ap

def rpd(x,y,weights):
    ap =  np.average( (x-y) / (np.abs(np.abs(x)+np.abs(y))/2), weights=weights ) *100
    return ap

def get_month_indices(months, numYears):
    indices = [];
    for year in range(numYears):
        for month in months:
            indices.append(year*12 + month);
    return indices;

def plot_shelf_region(data):
    plt.figure();
    depth = np.flipud(Dataset("data/GEBCO_bathymetry_0.25x0.25deg.nc", 'r').variables["mean_depth"][:]);
    plt.imshow(depth);
    plt.plot(data.xcoords, data.ycoords, 'y');
    midpoint = int(len(data.xcoords)/2);
    plt.plot(data.xcoords[0:midpoint], data.ycoords[0:midpoint], 'r');
    plt.plot(data.xcoords[midpoint:], data.ycoords[midpoint:], 'g');

#Select and return data for given set of months each year
#Returned data is reformatted into an array and stored in a QuickStruct, and contains:
#   overall mean and stddev of ekman, geostrophic and stokes current for the whole period.
#months defines the months of the year to use (e.g. 0=Jan, 11=Dec)
def get_processed_data_across_shelf(months, data, regionMaskBoundsList, params):
    if len(data) % 12 != 0:
        raise ValueError ("get_processed_data only supports data for whole years (multiples of 12 months)");

    #Create the indices to use
    numYears = int(len(data)/12);
    indices = get_month_indices(months, numYears);

    #Create arrays to store data for ekman, geostrophic and stokes data
    length = 0;
    for m in months:
        if len(data[m]) > length:
            length = len(data[m]);
    ekmanVals = np.empty((length, len(indices))); ekmanVals[:] = np.nan;
    geostrophicVals = np.empty((length, len(indices))); geostrophicVals[:] = np.nan;
    stokesVals = np.empty((length, len(indices))); stokesVals[:] = np.nan;
    stokesMaskPass = np.empty((length, len(indices))); stokesMaskPass[:] = np.nan;
    xcoords = np.empty((length, len(indices))); xcoords[:] = np.nan; #xcoords and ycoords are just used for checking region mask is correct
    ycoords = np.empty((length, len(indices))); ycoords[:] = np.nan;

    #Copy data:
    for i, index in enumerate(indices):
        ekmanVals[:,i] = [curData.nEkmanAcrossShelf for curData in data[index]];
        geostrophicVals[:,i] = [curData.nGeostrophicAcrossShelf for curData in data[index]];
        stokesVals[:,i] = [curData.nStokesAcrossShelf for curData in data[index]];
        stokesMaskPass[:,i] = [curData.stokesMaskPass for curData in data[index]];

        xcoords[:,i] = [curData.indexX for curData in data[index]];
        ycoords[:,i] = [curData.indexY for curData in data[index]];

    #Wherever stokesMaskPass fails, set stokes to 0
    stokesVals[stokesMaskPass==0.0] = 0.0;

    #Calculate means and standard deviations
    results = su.QuickStruct(); #Store results in a struct
    results.ekmanAcrossShelf = np.nanmean(ekmanVals, axis=1);
    results.ekmanAcrossShelfSD = np.nanstd(ekmanVals, axis=1);
    results.geostrophicAcrossShelf = np.nanmean(geostrophicVals, axis=1);
    results.geostrophicAcrossShelfSD = np.nanstd(geostrophicVals, axis=1);
    results.stokesAcrossShelf = np.nanmean(stokesVals, axis=1);
    results.stokesAcrossShelfSD = np.nanstd(stokesVals, axis=1);
    results.xcoords = np.nanmean(xcoords, axis=1);
    results.ycoords = np.nanmean(ycoords, axis=1);


    #get the regionMask and apply it
    regionMask = mf.return_area_mask(regionMaskBoundsList, data[0], params);
    results.regionMask = regionMask;
    results.ekmanAcrossShelf = results.ekmanAcrossShelf[regionMask];
    results.ekmanAcrossShelfSD = results.ekmanAcrossShelfSD[regionMask];
    results.geostrophicAcrossShelf = results.geostrophicAcrossShelf[regionMask];
    results.geostrophicAcrossShelfSD = results.geostrophicAcrossShelfSD[regionMask];
    results.stokesAcrossShelf = results.stokesAcrossShelf[regionMask];
    results.stokesAcrossShelfSD = results.stokesAcrossShelfSD[regionMask];
    results.xcoords = results.xcoords[regionMask];
    results.ycoords = results.ycoords[regionMask];

    return results;

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
    geostrophicdata = np.empty((length, len(indices))); geostrophicdata[:] = np.nan;

    #Copy data:
    for i, index in enumerate(indices):
        xcoords[:,i] = [curData.indexX for curData in data[index]];
        ycoords[:,i] = [curData.indexY for curData in data[index]];
        ekmanVals[:,i] = [curData.ekmanProportionGeoEk for curData in data[index]];
        geostrophicVals[:,i] = [curData.geostrophicProportionGeoEk for curData in data[index]];
        stokesVals[:,i] = [curData.stokesProportionOfTotal for curData in data[index]];
        totalVals[:,i] = [curData.totalcurrent for curData in data[index]];
        stokesPassMask[:,i] = [curData.stokesMaskPass for curData in data[index]];
        geostrophicdata[:,i] = [curData.nGeostrophicAcrossShelf for curData in data[index]];

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
    results.geostrophicdata =np.nanmean(geostrophicdata,axis=1)

    #Remove nans
    # toKeep = (np.isnan(results.ekmanProportions)==False) | (np.isnan(results.geostrophicProportions)==False) | (np.isnan(results.stokesProportions)==False) | (np.isnan(results.stokesPassMask)==False);
    # results.xcoords = results.xcoords[toKeep];
    # results.ycoords = results.ycoords[toKeep];
    # results.ekmanProportions = results.ekmanProportions[toKeep];
    # results.ekmanProportionsSD = results.ekmanProportionsSD[toKeep];
    # results.geostrophicProportions = results.geostrophicProportions[toKeep];
    # results.geostrophicProportionsSD = results.geostrophicProportionsSD[toKeep];
    # results.stokesProportions = results.stokesProportions[toKeep];
    # results.stokesProportionsSD = results.stokesProportionsSD[toKeep];
    # results.stokesPassMask = results.stokesPassMask[toKeep];
    # results.totalcurrent = results.totalcurrent[toKeep]
    # results.totalcurrentSD = results.totalcurrentSD[toKeep]
    # results.geostrophicdata = results.geostrophicdata[toKeep]

    return results;

def plot_shelf_section(segmentDistances, data123, data789, xlabel, yRange=None, outputFile=None, alpha=0.25, noClose=False, laruelleSections=[], laruelleSectionNames=[],ax=[],colourGeostrophic="#1f77b4",label =''):
    def plot_laruelle_data(ax):
        yRange = ax.get_ylim();
        for isection, section in enumerate(laruelleSections):
            ax.plot([section[0], section[1]], [yRange[0], yRange[0]], 'r', linewidth=2.5);

            #raise Exception("saldkfjas");
            textX = section[0] + ((section[1]-section[0]) / 2.0);
            textY = yRange[0]+0.015*(yRange[1]-yRange[0]);
            ax.text(textX, textY, laruelleSectionNames[isection], weight="bold", color='r',horizontalalignment='center');
    if ax == []:
        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(12,5));
    else:
        ax1 = ax[0]
        ax2 = ax[1]

    colourEkman="#2ca02c";
    #colourGeostrophic="#1f77b4";
    colourStokes="#ff7f0e";

    #data123
    ax1.fill_between(segmentDistances, data123.geostrophicAcrossShelf-data123.geostrophicAcrossShelfSD, data123.geostrophicAcrossShelf+data123.geostrophicAcrossShelfSD, alpha=alpha, color=colourGeostrophic);
    #ax1.fill_between(segmentDistances, data123.stokesAcrossShelf-data123.stokesAcrossShelfSD, data123.stokesAcrossShelf+data123.stokesAcrossShelfSD, alpha=alpha, color=colourStokes);
    #ax1.fill_between(segmentDistances, data123.ekmanAcrossShelf-data123.ekmanAcrossShelfSD, data123.ekmanAcrossShelf+data123.ekmanAcrossShelfSD, alpha=alpha, color=colourEkman);
    l2=ax1.plot(segmentDistances, data123.geostrophicAcrossShelf, label=label, linewidth=linewidth, color=colourGeostrophic);
    #l3=ax1.plot(segmentDistances, data123.stokesAcrossShelf, label="stokes", linewidth=linewidth, color=colourStokes);
    #l1=ax1.plot(segmentDistances, data123.ekmanAcrossShelf, label="ekman", linewidth=linewidth, color=colourEkman);
    #legendElements = [l1[0], l2[0], l3[0]];
    #ax1.legend(handles=legendElements,loc = 7);
    ax1.set_ylabel(r"Shelf break current ($m s^{-1}$)");
    ax1.set_xlabel(xlabel);
    ax1.grid()
    if yRange != None:
        ax1.set_ylim(yRange[0], yRange[1]);

    #Laruelle sections (plot horizontal lines to show the extend of each section)
    plot_laruelle_data(ax1);

    #data798
    ax2.fill_between(segmentDistances, data789.geostrophicAcrossShelf-data789.geostrophicAcrossShelfSD, data789.geostrophicAcrossShelf+data789.geostrophicAcrossShelfSD, alpha=alpha);
    #ax2.fill_between(segmentDistances, data789.stokesAcrossShelf-data789.stokesAcrossShelfSD, data789.stokesAcrossShelf+data789.stokesAcrossShelfSD, alpha=alpha);
    #ax2.fill_between(segmentDistances, data789.ekmanAcrossShelf-data789.ekmanAcrossShelfSD, data789.ekmanAcrossShelf+data789.ekmanAcrossShelfSD, alpha=alpha);
    l2=ax2.plot(segmentDistances, data789.geostrophicAcrossShelf, label=label, linewidth=linewidth);
    #l3=ax2.plot(segmentDistances, data789.stokesAcrossShelf, label="stokes", linewidth=linewidth);
    #l1=ax2.plot(segmentDistances, data789.ekmanAcrossShelf, label="ekman", linewidth=linewidth);
    #legendElements = [l1[0], l2[0], l3[0]];
    #ax2.legend(handles=legendElements,loc=7);
    ax2.set_ylabel(r"Shelf break current ($m s^{-1}$)");
    ax2.set_xlabel(xlabel);
    ax2.grid()
    if yRange != None:
        ax2.set_ylim(yRange[0], yRange[1]);

    #Laruelle sections (plot horizontal lines to show the extend of each section)
    plot_laruelle_data(ax2);

    plt.tight_layout();

    if outputFile != None:
        plt.savefig(outputFile,dpi=300);
        if noClose == False:
            plt.close;

def roll_section(data123, data789, distances, roll=170):
    def _do_roll(data, roll):
        data.xcoords = np.roll(data.xcoords, roll);
        data.ycoords = np.roll(data.ycoords, roll);
        data.ekmanAcrossShelf = np.roll(data.ekmanAcrossShelf, roll);
        data.ekmanAcrossShelfSD = np.roll(data.ekmanAcrossShelfSD, roll);
        data.geostrophicAcrossShelf = np.roll(data.geostrophicAcrossShelf, roll);
        data.geostrophicAcrossShelfSD = np.roll(data.geostrophicAcrossShelfSD, roll);
        data.stokesAcrossShelf = np.roll(data.stokesAcrossShelf, roll);
        data.stokesAcrossShelfSD = np.roll(data.stokesAcrossShelfSD, roll);
        return data;

    data123 = _do_roll(data123, roll);
    data789 = _do_roll(data789, roll);
    distances = np.roll(distances, roll);
    return data123, data789, distances;


def remove_section(data123, data789, distances, start, stop):
    def _do_remove(data, start, stop):
        data.xcoords = np.delete(data.xcoords, np.arange(start, stop));
        data.ycoords = np.delete(data.ycoords, np.arange(start, stop));
        data.ekmanAcrossShelf = np.delete(data.ekmanAcrossShelf, np.arange(start, stop));
        data.ekmanAcrossShelfSD = np.delete(data.ekmanAcrossShelfSD, np.arange(start, stop));
        data.geostrophicAcrossShelf = np.delete(data.geostrophicAcrossShelf, np.arange(start, stop));
        data.geostrophicAcrossShelfSD = np.delete(data.geostrophicAcrossShelfSD, np.arange(start, stop));
        data.stokesAcrossShelf = np.delete(data.stokesAcrossShelf, np.arange(start, stop));
        data.stokesAcrossShelfSD = np.delete(data.stokesAcrossShelfSD, np.arange(start, stop));
        return data;

    if start==None:
        start = 0;
    if stop==None:
        stop = len(data123.xcoords);

    data123 = _do_remove(data123, start, stop);
    data789 = _do_remove(data789, start, stop);
    distances = np.delete(distances, np.arange(start, stop));
    return data123, data789, distances;

def lon_lat_con(data,params):
    lats = []; lons = [];
    for i in range(len(data.xcoords)):
        lon, lat = su.convert_index_to_lonlat(data.xcoords[i], data.ycoords[i], params.pixelRes);
        lats.append(lat); lons.append(lon);
    return np.array(lats),np.array(lons)

#Global settings

# if params.paramsetName != "global": #Could actually use the shelf coordinates filename (which is also stored in params)...
#     raise ValueError("This plotting script is intended only for the 'global' parameter set. Significant adaptation is required for use with any other datasets that may change the shelf-coordinates.");
linewidth=1.5;
yRange = None;
closePlots=False;
plotTestplots=False;
ticksize = 18;
marker = 's'
m_sz = 5
v = 0.5
params = ps.get_global_params_glory(res=False)
#Read data
gridwiseData = pd.read_table(path.join('E:\SKIM', "gridwise_data", "per_grid_cell_edge_data_"+params.paramsetName+"_500m.csv"), sep=',');
gridwiseData.x = gridwiseData.x.astype(int);
gridwiseData.y = gridwiseData.y.astype(int);



allData = pickle.load(open(path.join("E:/SKIM", "current_data", "surface_currents_"+params.paramsetName+"_500m.p"), "rb"));


##### High res
high_params = ps.get_global_params_glory(res=True)
#Read data
high_gridwiseData = pd.read_table(path.join('E:\SKIM', "gridwise_data", "per_grid_cell_edge_data_"+high_params.paramsetName+"_500m.csv"), sep=',');
high_gridwiseData.x = high_gridwiseData.x.astype(int);
high_gridwiseData.y = high_gridwiseData.y.astype(int);

high_allData = pickle.load(open(path.join("E:/SKIM", "current_data", "surface_currents_"+high_params.paramsetName+"_500m.p"), "rb"));

fig = plt.figure(figsize=(45,30))
gs = GridSpec(3,2, figure=fig, wspace=0.2,hspace=0.2,bottom=0.05,top=0.97,left=0.05,right=0.93)

ax1 = fig.add_subplot(gs[0,0])
data_l = get_processed_data([0, 1, 2], allData);
lats_l,lons_l = lon_lat_con(data_l,params)

mapFig = Basemap(llcrnrlon=-180.0, llcrnrlat=-90, urcrnrlon=180.0, urcrnrlat=90.0, resolution='l', projection='cyl', lat_0 = 39.5, lon_0 = -3.25);
mapFig.drawmapboundary(fill_color=(0.85, 0.85, 0.85));
mapFig.fillcontinents(color=(0.5, 0.5, 0.5), zorder=1);
mapFig.drawmeridians(np.arange(0, 360, 60), labels=[0,0,0,1], color=(0.3, 0.3, 0.3), fontsize=ticksize);
mapFig.drawparallels(np.arange(-90, 90, 30), labels=[1,0,0,0], color=(0.3, 0.3, 0.3), fontsize=ticksize);
mapFig.scatter(lons_l, lats_l, latlon=True, c=data_l.geostrophicdata, marker=marker, cmap=plt.cm.RdBu,s = m_sz,vmin=-v,vmax=v);
ax1.set_title('January, February, March, 1993-1994',fontsize=ticksize+3,fontweight='bold')
ax1.set_ylabel('Low (1/4 deg) resolution',fontsize=ticksize+3,fontweight='bold',labelpad=100)


ax2 = fig.add_subplot(gs[1,0])
data = get_processed_data([0, 1, 2], high_allData);
lats,lons = lon_lat_con(data,high_params)
mapFig = Basemap(llcrnrlon=-180.0, llcrnrlat=-90, urcrnrlon=180.0, urcrnrlat=90.0, resolution='l', projection='cyl', lat_0 = 39.5, lon_0 = -3.25);
mapFig.drawmapboundary(fill_color=(0.85, 0.85, 0.85));
mapFig.fillcontinents(color=(0.5, 0.5, 0.5), zorder=1);
mapFig.drawmeridians(np.arange(0, 360, 60), labels=[0,0,0,1], color=(0.3, 0.3, 0.3), fontsize=ticksize);
mapFig.drawparallels(np.arange(-90, 90, 30), labels=[1,0,0,0], color=(0.3, 0.3, 0.3), fontsize=ticksize);
ap = mapFig.scatter(lons, lats, latlon=True, c=data.geostrophicdata, marker=marker, cmap=plt.cm.RdBu,s = m_sz,vmin=-v,vmax=v);
ax1.text(0.94,0.97,'(a)',transform=ax1.transAxes,va='top',fontsize=ticksize,fontweight='bold')
ax2.text(0.94,0.97,'(b)',transform=ax2.transAxes,va='top',fontsize=ticksize,fontweight='bold')
ax2.set_ylabel('High (1/12 deg) resolution',fontsize=ticksize+3,fontweight='bold',labelpad=100)
ax = fig.add_axes([0.93,0.45,0.02,0.4])
cbar = plt.colorbar(ap, cax=ax)
cbar.set_label('Across shelf current (ms$^{-1}$)',fontsize=ticksize);
plt.yticks(fontsize=ticksize)

res = 1
lon_deg = np.arange(-180+res,180,2)
lat_deg = np.arange(-90+res,89,2)

low_res = np.zeros((len(lon_deg),len(lat_deg)))
high_res = np.copy(low_res)

for i in range(len(lon_deg)):
    print(i)
    for j in range(len(lat_deg)):
        f = np.where((lons < lon_deg[i]+res) & (lons > lon_deg[i]-res) & (lats <lat_deg[j]+res) & (lats >lat_deg[j]-res) )
        high_res[i,j] = np.nanmean(data.geostrophicdata[f])

        f = np.where((lons_l < lon_deg[i]+res) & (lons_l > lon_deg[i]-res) & (lats_l <lat_deg[j]+res) & (lats_l >lat_deg[j]-res) )
        low_res[i,j] = np.nanmean(data_l.geostrophicdata[f])
dif = (high_res - low_res)
ax3 = fig.add_subplot(gs[2,0])
mapFig = Basemap(llcrnrlon=-180.0, llcrnrlat=-90, urcrnrlon=180.0, urcrnrlat=90.0, resolution='l', projection='cyl', lat_0 = 39.5, lon_0 = -3.25);
mapFig.drawmapboundary(fill_color=(0.85, 0.85, 0.85));
mapFig.fillcontinents(color=(0.5, 0.5, 0.5), zorder=1);
mapFig.drawmeridians(np.arange(0, 360, 60), labels=[0,0,0,1], color=(0.3, 0.3, 0.3), fontsize=ticksize);
mapFig.drawparallels(np.arange(-90, 90, 30), labels=[1,0,0,0], color=(0.3, 0.3, 0.3), fontsize=ticksize);
mapFig.pcolor(lon_deg,lat_deg,np.transpose(dif),latlon=True,vmin=-0.25,vmax=0.25,cmap=plt.cm.RdBu,zorder=2)
ax3.set_ylabel('Difference',fontsize=ticksize+3,fontweight='bold',labelpad=100)


ax1 = fig.add_subplot(gs[0,1])
data_l = get_processed_data([7, 8, 9], allData);
lats_l,lons_l = lon_lat_con(data_l,params)
mapFig = Basemap(llcrnrlon=-180.0, llcrnrlat=-90, urcrnrlon=180.0, urcrnrlat=90.0, resolution='l', projection='cyl', lat_0 = 39.5, lon_0 = -3.25);
mapFig.drawmapboundary(fill_color=(0.85, 0.85, 0.85));
mapFig.fillcontinents(color=(0.5, 0.5, 0.5), zorder=1);
mapFig.drawmeridians(np.arange(0, 360, 60), labels=[0,0,0,1], color=(0.3, 0.3, 0.3), fontsize=ticksize);
mapFig.drawparallels(np.arange(-90, 90, 30), labels=[1,0,0,0], color=(0.3, 0.3, 0.3), fontsize=ticksize);
mapFig.scatter(lons_l, lats_l, latlon=True, c=data_l.geostrophicdata, marker=marker, cmap=plt.cm.RdBu,s = m_sz,vmin=-v,vmax=v);
ax1.set_title('July, August, September, 1993-1994',fontsize=ticksize+3,fontweight='bold')


ax2 = fig.add_subplot(gs[1,1])
data = get_processed_data([7, 8, 9], high_allData);
lats,lons = lon_lat_con(data,high_params)
mapFig = Basemap(llcrnrlon=-180.0, llcrnrlat=-90, urcrnrlon=180.0, urcrnrlat=90.0, resolution='l', projection='cyl', lat_0 = 39.5, lon_0 = -3.25);
mapFig.drawmapboundary(fill_color=(0.85, 0.85, 0.85));
mapFig.fillcontinents(color=(0.5, 0.5, 0.5), zorder=1);
mapFig.drawmeridians(np.arange(0, 360, 60), labels=[0,0,0,1], color=(0.3, 0.3, 0.3), fontsize=ticksize);
mapFig.drawparallels(np.arange(-90, 90, 30), labels=[1,0,0,0], color=(0.3, 0.3, 0.3), fontsize=ticksize);
mapFig.scatter(lons, lats, latlon=True, c=data.geostrophicdata, marker=marker, cmap=plt.cm.RdBu,s = m_sz,vmin=-v,vmax=v);
ax1.text(0.94,0.97,'(d)',transform=ax1.transAxes,va='top',fontsize=ticksize,fontweight='bold')
ax2.text(0.94,0.97,'(e)',transform=ax2.transAxes,va='top',fontsize=ticksize,fontweight='bold')


low_res = np.zeros((len(lon_deg),len(lat_deg)))
high_res = np.copy(low_res)

for i in range(len(lon_deg)):
    print(i)
    for j in range(len(lat_deg)):
        f = np.where((lons < lon_deg[i]+res) & (lons > lon_deg[i]-res) & (lats <lat_deg[j]+res) & (lats >lat_deg[j]-res) )
        high_res[i,j] = np.nanmean(data.geostrophicdata[f])

        f = np.where((lons_l < lon_deg[i]+res) & (lons_l > lon_deg[i]-res) & (lats_l <lat_deg[j]+res) & (lats_l >lat_deg[j]-res) )
        low_res[i,j] = np.nanmean(data_l.geostrophicdata[f])
dif = (high_res - low_res)
ax4 = fig.add_subplot(gs[2,1])
mapFig = Basemap(llcrnrlon=-180.0, llcrnrlat=-90, urcrnrlon=180.0, urcrnrlat=90.0, resolution='l', projection='cyl', lat_0 = 39.5, lon_0 = -3.25);
mapFig.drawmapboundary(fill_color=(0.85, 0.85, 0.85));
mapFig.fillcontinents(color=(0.5, 0.5, 0.5), zorder=1);
mapFig.drawmeridians(np.arange(0, 360, 60), labels=[0,0,0,1], color=(0.3, 0.3, 0.3), fontsize=ticksize);
mapFig.drawparallels(np.arange(-90, 90, 30), labels=[1,0,0,0], color=(0.3, 0.3, 0.3), fontsize=ticksize);
ap = mapFig.pcolor(lon_deg,lat_deg,np.transpose(dif),latlon=True,vmin=-0.25,vmax=0.25,cmap=plt.cm.RdBu,zorder=2)
ax3.text(0.94,0.97,'(c)',transform=ax3.transAxes,va='top',fontsize=ticksize,fontweight='bold')
ax4.text(0.94,0.97,'(f)',transform=ax4.transAxes,va='top',fontsize=ticksize,fontweight='bold')

ax = fig.add_axes([0.93,0.05,0.02,0.27])
cbar = plt.colorbar(ap, cax=ax)
cbar.set_label('Difference in across-shelf flow (ms$^{-1}$)\n(high resolution - low resolution)',fontsize=ticksize);
plt.yticks(fontsize=ticksize)
fig.savefig("plots/current_along_shelf/highres_compare/manuscript_spatial_figure_global.png",dpi=300);


"""
Manuscript figure copy with high and low resoltuion data
"""
# font = {'weight' : 'normal',
#         'size'   : 18}
# matplotlib.rc('font', **font)
# fig = plt.figure(figsize=(15,15))
# gs = GridSpec(2,2, figure=fig, wspace=0.33,hspace=0.2,bottom=0.05,top=0.97,left=0.1,right=0.97)
# ax1 = fig.add_subplot(gs[0,0])
# ax2 = fig.add_subplot(gs[0,1])
# ax5 = fig.add_subplot(gs[1,0])
# ax6 = fig.add_subplot(gs[1,1])
# ax1.text(0.90,0.95,'(a)',transform=ax1.transAxes,va='top',fontsize=18,fontweight='bold')
# ax2.text(0.90,0.95,'(b)',transform=ax2.transAxes,va='top',fontsize=18,fontweight='bold')
# ax5.text(0.90,0.95,'(c)',transform=ax5.transAxes,va='top',fontsize=18,fontweight='bold')
# ax6.text(0.90,0.95,'(d)',transform=ax6.transAxes,va='top',fontsize=18,fontweight='bold')
#
# data123 = get_processed_data_across_shelf([0,1,2], allData, [mf.area_EuropeanShelf], params);
# data789 = get_processed_data_across_shelf([6,7,8], allData, [mf.area_EuropeanShelf], params);
# distances = np.array(gridwiseData["distance"])[data123.regionMask]/1000.0; #in km
#
# data123, data789, distances = remove_section(data123, data789, distances, start=370, stop=None);
# #plot_shelf_region(data123);
#
# dists1, laruNorthSeaMask = mf.calculate_km_subsection_bounds_along_shelf(data123.xcoords, data123.ycoords, distances, mf.laruelle_NorthSea, params, testPlot=plotTestplots);
# dists2, laruEnglishChannelMask = mf.calculate_km_subsection_bounds_along_shelf(data123.xcoords, data123.ycoords, distances, mf.laruelle_EnglishChannel, params, testPlot=plotTestplots);
#
# segmentDistances = np.cumsum(distances);
# dis_nw = segmentDistances[-1]
# xlabel = "distance along shelf (km, North to South)";
# #outputFile = "plots/current_along_shelf/flux_along_shelf_europeanShelf.png"
# plot_shelf_section(segmentDistances, data123, data789, xlabel=xlabel, yRange=yRange, outputFile=None, noClose=closePlots, laruelleSections=dists1+dists2, laruelleSectionNames=["North Sea", "English Channel"],ax = [ax1,ax2],label='1/4 degree');
#
#
#
#
# data123 = get_processed_data_across_shelf([0,1,2], allData, [mf.area_MidAtlanticBight], params);
# data789 = get_processed_data_across_shelf([6,7,8], allData, [mf.area_MidAtlanticBight], params);
#
# distances = np.array(gridwiseData["distance"])[data123.regionMask]/1000.0; #in km
# dists, laruMidAtlanticBightMask = mf.calculate_km_subsection_bounds_along_shelf(data123.xcoords, data123.ycoords, distances, mf.laruelle_MidAtlanticBight, params, testPlot=plotTestplots);
#
# segmentDistances = np.cumsum(distances);
# dis_sa = segmentDistances[-1]
# xlabel = "distance along shelf (km, South to North)";
# #outputFile = "plots/current_along_shelf/flux_along_shelf_midAtlanticBight.png"
# plot_shelf_section(segmentDistances, data123, data789, xlabel=xlabel, yRange=yRange, outputFile=None, noClose=closePlots, laruelleSections=dists, laruelleSectionNames=["Mid-Atlantic Bight"],ax = [ax5,ax6],label='1/4 degree');
#
#
# #High res
# data123 = get_processed_data_across_shelf([0,1,2], high_allData, [mf.area_EuropeanShelf], high_params);
# data789 = get_processed_data_across_shelf([6,7,8], high_allData, [mf.area_EuropeanShelf], high_params);
# distances = np.array(high_gridwiseData["distance"])[data123.regionMask]/1000.0; #in km
# data123, data789, distances = remove_section(data123, data789, distances, start=390*3, stop=None);
# #plot_shelf_region(data123);
#
# dists1, laruNorthSeaMask = mf.calculate_km_subsection_bounds_along_shelf(data123.xcoords, data123.ycoords, distances, mf.laruelle_NorthSea, params, testPlot=plotTestplots);
# dists2, laruEnglishChannelMask = mf.calculate_km_subsection_bounds_along_shelf(data123.xcoords, data123.ycoords, distances, mf.laruelle_EnglishChannel, params, testPlot=plotTestplots);
#
# segmentDistances = np.cumsum(distances);
# f = segmentDistances[-1] / dis_nw
# segmentDistances = segmentDistances / f
# xlabel = "distance along shelf (km, North to South)";
# #outputFile = "plots/current_along_shelf/flux_along_shelf_europeanShelf.png"
# plot_shelf_section(segmentDistances, data123, data789, xlabel=xlabel, yRange=yRange, outputFile=None, noClose=closePlots,ax = [ax1,ax2],colourGeostrophic="#ff7f0e",label='1/12 degree');
# ax1.plot([0,segmentDistances[-1]],[0,0],'k--')
# ax2.plot([0,segmentDistances[-1]],[0,0],'k--')
#
# data123 = get_processed_data_across_shelf([0,1,2], high_allData, [mf.area_MidAtlanticBight], high_params);
# data789 = get_processed_data_across_shelf([6,7,8], high_allData, [mf.area_MidAtlanticBight], high_params);
#
# distances = np.array(high_gridwiseData["distance"])[data123.regionMask]/1000.0; #in km
# dists, laruMidAtlanticBightMask = mf.calculate_km_subsection_bounds_along_shelf(data123.xcoords, data123.ycoords, distances, mf.laruelle_MidAtlanticBight, params, testPlot=plotTestplots);
#
# segmentDistances = np.cumsum(distances);
# f = segmentDistances[-1] / dis_sa
# segmentDistances = segmentDistances / f
#
# xlabel = "distance along shelf (km, South to North)";
# #outputFile = "plots/current_along_shelf/flux_along_shelf_midAtlanticBight.png"
# plot_shelf_section(segmentDistances, data123, data789, xlabel=xlabel, yRange=yRange, outputFile=None, noClose=closePlots, ax = [ax5,ax6],colourGeostrophic="#ff7f0e",label='1/12 degree');
#
# ax5.plot([0,segmentDistances[-1]],[0,0],'k--')
# ax6.plot([0,segmentDistances[-1]],[0,0],'k--')
# ax1.legend(loc = 7)
# ax2.legend(loc = 7)
# ax5.legend(loc = 7)
# ax6.legend(loc = 7)
# fig.savefig("plots/current_along_shelf/highres_compare/manuscript_spatial_figure.png",dpi=300);
# #

"""

"""

#1/4 degree data
low_data = get_processed_data(list(range(0,12)),allData)
lats_l,lons_l = lon_lat_con(data_l,params)

#1/12 degree data
high_data = get_processed_data(list(range(0,12)),high_allData)
lats,lons = lon_lat_con(data,high_params)
matches = []
for i in range(len(lats_l)):
    matches.append([])

for i in range(len(lats)):
    h = lats[i]-lats_l
    g = lons[i]-lons_l
    dis = np.sqrt(h**2 + g**2)
    f = np.where(np.min(dis) == dis)
    #print(f)
    if np.min(dis) < 1:
        matches[f[0][0]].append(i)

high_av = np.zeros((len(matches))); high_av[:] = np.nan
for i in range(len(matches)):
    if len(matches[i])>0:
        high_av[i] = np.nanmean(high_data.geostrophicdata[matches[i]])

font = {'weight' : 'normal',
        'size'   : 18}
matplotlib.rc('font', **font)
fig = plt.figure(figsize=(15,15))
gs = GridSpec(2,2, figure=fig, wspace=0.33,hspace=0.15,bottom=0.08,top=0.97,left=0.1,right=0.97)

ax2 = fig.add_subplot(gs[1,:])
# ax3 = fig.add_subplot(gs[1,1])
ax1 = fig.add_subplot(gs[0,:])
mapFig = Basemap(llcrnrlon=-180.0, llcrnrlat=-90, urcrnrlon=180.0, urcrnrlat=90.0, resolution='l', projection='cyl', lat_0 = 39.5, lon_0 = -3.25);
mapFig.drawmapboundary(fill_color=(0.85, 0.85, 0.85));
mapFig.fillcontinents(color=(0.5, 0.5, 0.5), zorder=1);
mapFig.drawmeridians(np.arange(0, 360, 60), labels=[0,0,0,1], color=(0.3, 0.3, 0.3), fontsize=ticksize);
mapFig.drawparallels(np.arange(-90, 90, 30), labels=[1,0,0,0], color=(0.3, 0.3, 0.3), fontsize=ticksize);
ap = mapFig.scatter(lons_l,lats_l,latlon=True,c = high_av-low_data.geostrophicdata,vmin=-0.25,vmax=0.25,cmap=plt.cm.RdBu,zorder=2,marker=marker, s = m_sz)
cbar = plt.colorbar(ap)
cbar.set_label('Difference in mean annual across-shelf flow (ms$^{-1}$)\n(1/12 degree - 1/4 degree)');
ax2.scatter(low_data.geostrophicdata,high_av)
ax2.set_xlabel('1/4 degree mean annual across-shelf\ncurrent 1993-1994 (ms$^{-1}$)')
ax2.set_ylabel('1/12 degree mean annual across-shelf\ncurrent 1993-1994 (ms$^{-1}$)')
c = np.array([-1.5,1.5])
ax2.plot(c,c,'k-')
ax2.set_xlim(c); ax2.set_ylim(c)
ax2.grid()
unweight(low_data.geostrophicdata,high_av,ax2,unit = 'ms$^{-1}$')
ax1.text(0.05,0.95,'(a)',transform=ax1.transAxes,va='top',fontsize=18,fontweight='bold')
ax2.text(0.05,0.95,'(b)',transform=ax2.transAxes,va='top',fontsize=18,fontweight='bold')

# #1/4 degree data
# low_data = get_processed_data([0,1,2],allData)
# lats_l,lons_l = lon_lat_con(data_l,params)
#
# #1/12 degree data
# high_data = get_processed_data([0,1,2],high_allData)
# lats,lons = lon_lat_con(data,high_params)
# matches = []
# for i in range(len(lats_l)):
#     matches.append([])
#
# for i in range(len(lats)):
#     h = lats[i]-lats_l
#     g = lons[i]-lons_l
#     dis = np.sqrt(h**2 + g**2)
#     f = np.where(np.min(dis) == dis)
#     #print(f)
#     if np.min(dis) < 1:
#         matches[f[0][0]].append(i)
#
# high_av = np.zeros((len(matches))); high_av[:] = np.nan
# for i in range(len(matches)):
#     if len(matches[i])>0:
#         high_av[i] = np.nanmean(high_data.geostrophicdata[matches[i]])
#
# ax3.scatter(low_data.geostrophicdata,high_av)
# ax3.set_xlabel('1/4 degree mean winter across-shelf\ncurrent 1993-1994 (ms$^{-1}$)')
# ax3.set_ylabel('1/12 degree mean winter across-shelf\ncurrent 1993-1994 (ms$^{-1}$)')
# c = np.array([-1.5,1.5])
# ax3.plot(c,c,'k-')
# ax3.set_xlim(c); ax3.set_ylim(c)
#
# unweight(low_data.geostrophicdata,high_av,ax3,unit = 'ms$^{-1}$')

fig.savefig("plots/current_along_shelf/highres_compare/scatter_high_lows.png",dpi=300);
# print(matches)
