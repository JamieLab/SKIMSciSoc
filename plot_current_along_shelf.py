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
import pandas as pd;
import matplotlib.transforms

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


def plot_shelf_section(segmentDistances, data123, data789, xlabel, yRange=None, outputFile=None, alpha=0.25, noClose=False, laruelleSections=[], laruelleSectionNames=[],ax=[]):
    def plot_laruelle_data(ax,laruelleSections,laruelleSectionNames):
        yRange = ax.get_ylim();
        print(laruelleSections)
        print(laruelleSectionNames)
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
    colourGeostrophic="#1f77b4";
    colourStokes="#ff7f0e";

    #data123
    ax1.fill_between(segmentDistances, data123.geostrophicAcrossShelf-data123.geostrophicAcrossShelfSD, data123.geostrophicAcrossShelf+data123.geostrophicAcrossShelfSD, alpha=alpha, color=colourGeostrophic);
    ax1.fill_between(segmentDistances, data123.stokesAcrossShelf-data123.stokesAcrossShelfSD, data123.stokesAcrossShelf+data123.stokesAcrossShelfSD, alpha=alpha, color=colourStokes);
    ax1.fill_between(segmentDistances, data123.ekmanAcrossShelf-data123.ekmanAcrossShelfSD, data123.ekmanAcrossShelf+data123.ekmanAcrossShelfSD, alpha=alpha, color=colourEkman);
    l2=ax1.plot(segmentDistances, data123.geostrophicAcrossShelf, label="geostrophic", linewidth=linewidth, color=colourGeostrophic);
    l3=ax1.plot(segmentDistances, data123.stokesAcrossShelf, label="stokes", linewidth=linewidth, color=colourStokes);
    l1=ax1.plot(segmentDistances, data123.ekmanAcrossShelf, label="ekman", linewidth=linewidth, color=colourEkman);
    legendElements = [l1[0], l2[0], l3[0]];
    ax1.legend(handles=legendElements,loc = 7);
    ax1.set_ylabel(r"Shelf break current ($m s^{-1}$)");
    ax1.set_xlabel(xlabel);
    ax1.grid()
    if yRange != None:
        ax1.set_ylim(yRange[0], yRange[1]);

    #Laruelle sections (plot horizontal lines to show the extend of each section)
    plot_laruelle_data(ax1,laruelleSections,laruelleSectionNames);

    #data798
    ax2.fill_between(segmentDistances, data789.geostrophicAcrossShelf-data789.geostrophicAcrossShelfSD, data789.geostrophicAcrossShelf+data789.geostrophicAcrossShelfSD, alpha=alpha);
    ax2.fill_between(segmentDistances, data789.stokesAcrossShelf-data789.stokesAcrossShelfSD, data789.stokesAcrossShelf+data789.stokesAcrossShelfSD, alpha=alpha);
    ax2.fill_between(segmentDistances, data789.ekmanAcrossShelf-data789.ekmanAcrossShelfSD, data789.ekmanAcrossShelf+data789.ekmanAcrossShelfSD, alpha=alpha);
    l2=ax2.plot(segmentDistances, data789.geostrophicAcrossShelf, label="geostrophic", linewidth=linewidth);
    l3=ax2.plot(segmentDistances, data789.stokesAcrossShelf, label="stokes", linewidth=linewidth);
    l1=ax2.plot(segmentDistances, data789.ekmanAcrossShelf, label="ekman", linewidth=linewidth);
    legendElements = [l1[0], l2[0], l3[0]];
    ax2.legend(handles=legendElements,loc=7);
    ax2.set_ylabel(r"Shelf break current ($m s^{-1}$)");
    ax2.set_xlabel(xlabel);
    ax2.grid()
    if yRange != None:
        ax2.set_ylim(yRange[0], yRange[1]);

    #Laruelle sections (plot horizontal lines to show the extend of each section)
    plot_laruelle_data(ax2,laruelleSections,laruelleSectionNames);

    plt.tight_layout();

    if outputFile != None:
        plt.savefig(outputFile,dpi=300);
        if noClose == False:
            plt.close;



#Control what to plot
plotAll = True;
plotEuropeanShelf = False or plotAll;
plotMidAtlanticBight = False or plotAll;
plotCoastOfJapan = False or plotAll;
plotPatagonia = False or plotAll;
plotBeringSea = False or plotAll;
plotAntarcticPeninsula = False or plotAll;
plotLabradorSea = False or plotAll;
plotTasmania = False or plotAll;
plotBarentsSea = False or plotAll;
plotSouthAtlanticBight = False or plotAll;
plotSouthGreenland = False or plotAll;
plotCascianShelf = False or plotAll;
plotIrmingerSea = False or plotAll;
plot_manuscript = False or plotAll;



#Global settings
params = ps.get_global_params(cmems=True);
# if params.paramsetName != "global": #Could actually use the shelf coordinates filename (which is also stored in params)...
#     raise ValueError("This plotting script is intended only for the 'global' parameter set. Significant adaptation is required for use with any other datasets that may change the shelf-coordinates.");
linewidth=1.5;
yRange = (-0.20, 0.35);
closePlots=False;
plotTestplots=False;


#Read data
gridwiseData = pd.read_table(path.join('E:\SKIM', "gridwise_data", "per_grid_cell_edge_data_"+params.paramsetName+"_500m.csv"), sep=',');
gridwiseData.x = gridwiseData.x.astype(int);
gridwiseData.y = gridwiseData.y.astype(int);
if ('allData' in globals()) == False:
    ans = raw_input("Press key to read in 'allData', ctrl+c to cancel...");
    allData = pickle.load(open(path.join("E:/SKIM", "current_data", "surface_currents_"+params.paramsetName+"_500m.p"), "rb"));


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


###### Plot for January, February, March (2010-2016)
yRange = None;
if plotEuropeanShelf:
    data123 = get_processed_data_across_shelf([0,1,2], allData, [mf.area_EuropeanShelf], params);
    data789 = get_processed_data_across_shelf([6,7,8], allData, [mf.area_EuropeanShelf], params);
    distances = np.array(gridwiseData["distance"])[data123.regionMask]/1000.0; #in km

    data123, data789, distances = remove_section(data123, data789, distances, start=370, stop=None);
    #plot_shelf_region(data123);

    dists1, laruNorthSeaMask = mf.calculate_km_subsection_bounds_along_shelf(data123.xcoords, data123.ycoords, distances, mf.laruelle_NorthSea, params, testPlot=plotTestplots);
    dists2, laruEnglishChannelMask = mf.calculate_km_subsection_bounds_along_shelf(data123.xcoords, data123.ycoords, distances, mf.laruelle_EnglishChannel, params, testPlot=plotTestplots);

    segmentDistances = np.cumsum(distances);
    xlabel = "distance along shelf (km, North to South)";
    outputFile = "plots/current_along_shelf/flux_along_shelf_europeanShelf.png"
    plot_shelf_section(segmentDistances, data123, data789, xlabel=xlabel, yRange=yRange, outputFile=outputFile, noClose=closePlots, laruelleSections=dists1+dists2, laruelleSectionNames=["North Sea", "English Channel"]);

if plotMidAtlanticBight:
    data123 = get_processed_data_across_shelf([0,1,2], allData, [mf.area_MidAtlanticBight], params);
    data789 = get_processed_data_across_shelf([6,7,8], allData, [mf.area_MidAtlanticBight], params);

    distances = np.array(gridwiseData["distance"])[data123.regionMask]/1000.0; #in km
    dists, laruMidAtlanticBightMask = mf.calculate_km_subsection_bounds_along_shelf(data123.xcoords, data123.ycoords, distances, mf.laruelle_MidAtlanticBight, params, testPlot=plotTestplots);

    segmentDistances = np.cumsum(distances);
    xlabel = "distance along shelf (km, South to North)";
    outputFile = "plots/current_along_shelf/flux_along_shelf_midAtlanticBight.png"
    plot_shelf_section(segmentDistances, data123, data789, xlabel=xlabel, yRange=yRange, outputFile=outputFile, noClose=closePlots, laruelleSections=dists, laruelleSectionNames=["Mid-Atlantic Bight"]);
    #plot_shelf_region(data123);

if plotCoastOfJapan:
    data123 = get_processed_data_across_shelf([0,1,2], allData, [mf.area_CoastOfJapan], params);
    data789 = get_processed_data_across_shelf([6,7,8], allData, [mf.area_CoastOfJapan], params);
    distances = np.array(gridwiseData["distance"])[data123.regionMask]/1000.0; #in km

    data123, data789, distances = remove_section(data123, data789, distances, start=0, stop=197);
    data123, data789, distances = roll_section(data123, data789, distances, roll=-110);

    #plot_shelf_region(data123);

    dists, laruCoastOfJapanMask = mf.calculate_km_subsection_bounds_along_shelf(data123.xcoords, data123.ycoords, distances, mf.laruelle_CoastOfJapan, params, testPlot=plotTestplots);

    segmentDistances = np.cumsum(distances);
    xlabel = "distance along shelf (km, Southwest to Northeast)";
    outputFile = "plots/current_along_shelf/flux_along_shelf_coastOfJapan.png"

    plot_shelf_section(segmentDistances, data123, data789, xlabel=xlabel, yRange=yRange, outputFile=outputFile, noClose=closePlots, laruelleSections=[dists[0]], laruelleSectionNames=["Coast of Japan"]);

if plotPatagonia:
    data123 = get_processed_data_across_shelf([0,1,2], allData, [mf.area_Patagonia], params);
    data789 = get_processed_data_across_shelf([6,7,8], allData, [mf.area_Patagonia], params);

    distances = np.array(gridwiseData["distance"])[data123.regionMask]/1000.0; #in km
    dists, laruPatagonianShelfMask = mf.calculate_km_subsection_bounds_along_shelf(data123.xcoords, data123.ycoords, distances, mf.laruelle_PatagonianShelf, params, testPlot=plotTestplots);

    segmentDistances = np.cumsum(distances);
    xlabel = "distance along shelf (km, West to East)";
    outputFile = "plots/current_along_shelf/flux_along_shelf_patagonia.png"
    plot_shelf_section(segmentDistances, data123, data789, xlabel=xlabel, yRange=yRange, outputFile=outputFile, noClose=closePlots, laruelleSections=dists, laruelleSectionNames=["Patagonian Shelf"]);
    #plot_shelf_region(data123);

if plotBeringSea:
    data123 = get_processed_data_across_shelf([0,1,2], allData, [mf.area_BeringSeaWest], params);
    data789 = get_processed_data_across_shelf([6,7,8], allData, [mf.area_BeringSeaWest], params);

    distances = np.array(gridwiseData["distance"])[data123.regionMask]/1000.0; #in km
    dists, laruBeringSeaMask = mf.calculate_km_subsection_bounds_along_shelf(data123.xcoords, data123.ycoords, distances, mf.laruelle_BeringSea, params, testPlot=plotTestplots);

    segmentDistances = np.cumsum(distances);
    xlabel = "distance along shelf (km, West to East)";
    outputFile = "plots/current_along_shelf/flux_along_shelf_BeringSea.png"
    plot_shelf_section(segmentDistances, data123, data789, xlabel=xlabel, yRange=yRange, outputFile=outputFile, noClose=closePlots, laruelleSections=dists, laruelleSectionNames=["Bering Sea", "Bering Sea"]);
    #plot_shelf_region(data123);

if plotAntarcticPeninsula:
    data123 = get_processed_data_across_shelf([0,1,2], allData, [mf.area_AntarcticPeninsula], params);
    data789 = get_processed_data_across_shelf([6,7,8], allData, [mf.area_AntarcticPeninsula], params);
    distances = np.array(gridwiseData["distance"])[data123.regionMask]/1000.0; #in km

    data123, data789, distances = remove_section(data123, data789, distances, start=426, stop=444); #cut out the 'loop'
    data123, data789, distances = roll_section(data123, data789, distances, roll=337);
    #plot_shelf_region(data123); plt.xlim(370,550); plt.ylim(680,560);

    dists, laruAntarcticPeninsulaMask = mf.calculate_km_subsection_bounds_along_shelf(data123.xcoords, data123.ycoords, distances, mf.laruelle_AntarcticPeninsula, params, testPlot=plotTestplots);

    segmentDistances = np.cumsum(distances);
    xlabel = "distance along shelf (km, East to West)";
    outputFile = "plots/current_along_shelf/flux_along_shelf_AntarcticPeninsula.png"
    plot_shelf_section(segmentDistances, data123, data789, xlabel=xlabel, yRange=yRange, outputFile=outputFile, noClose=closePlots, laruelleSections=dists, laruelleSectionNames=["Antarctic Peninsula", "Antarctic Peninsula"]);


if plotLabradorSea:
    data123 = get_processed_data_across_shelf([0,1,2], allData, [mf.area_LabradorSea], params);
    data789 = get_processed_data_across_shelf([6,7,8], allData, [mf.area_LabradorSea], params);
    distances = np.array(gridwiseData["distance"])[data123.regionMask]/1000.0; #in km

    data123, data798, distances = roll_section(data123, data789, distances, roll=-154);

    #plot_shelf_region(data123); plt.xlim(400, 610); plt.ylim(225, 25);

    dists, laruLabradorSeaMask = mf.calculate_km_subsection_bounds_along_shelf(data123.xcoords, data123.ycoords, distances, mf.laruelle_LabradorSea, params, testPlot=plotTestplots);

    segmentDistances = np.cumsum(distances);
    xlabel = "distance along shelf\n(km, Southwest to North to Southeast)";
    outputFile = "plots/current_along_shelf/flux_along_shelf_LabradorSea.png"
    plot_shelf_section(segmentDistances, data123, data789, xlabel=xlabel, yRange=yRange, outputFile=outputFile, noClose=closePlots, laruelleSections=dists, laruelleSectionNames=["Labrador Sea"]);


if plotTasmania:
    data123 = get_processed_data_across_shelf([0,1,2], allData, [mf.area_Tasmania], params);
    data789 = get_processed_data_across_shelf([6,7,8], allData, [mf.area_Tasmania], params);

    distances = np.array(gridwiseData["distance"])[data123.regionMask]/1000.0; #in km
    dists, laruTasmanianShelfMask = mf.calculate_km_subsection_bounds_along_shelf(data123.xcoords, data123.ycoords, distances, mf.laruelle_TasmanianShelf, params, testPlot=plotTestplots);

    segmentDistances = np.cumsum(distances);
    xlabel = "distance along shelf (km, West to East)";
    outputFile = "plots/current_along_shelf/flux_along_shelf_Tasmania.png"
    plot_shelf_section(segmentDistances, data123, data789, xlabel=xlabel, yRange=yRange, outputFile=outputFile, noClose=closePlots, laruelleSections=dists, laruelleSectionNames=["Tasmanian Shelf"]);
    #plot_shelf_region(data123);

if plotBarentsSea:
    data123 = get_processed_data_across_shelf([0,1,2], allData, [mf.area_BarentsSea], params);
    data789 = get_processed_data_across_shelf([6,7,8], allData, [mf.area_BarentsSea], params);

    distances = np.array(gridwiseData["distance"])[data123.regionMask]/1000.0; #in km
    dists, laruBarentsSeaMask = mf.calculate_km_subsection_bounds_along_shelf(data123.xcoords, data123.ycoords, distances, mf.laruelle_BarentsSea, params, testPlot=plotTestplots);

    segmentDistances = np.cumsum(distances);
    xlabel = "distance along shelf (km, East to West)";
    outputFile = "plots/current_along_shelf/flux_along_shelf_BarentsSea.png"
    plot_shelf_section(segmentDistances, data123, data789, xlabel=xlabel, yRange=yRange, outputFile=outputFile, noClose=closePlots, laruelleSections=dists, laruelleSectionNames=["Barents Sea"]);
    #plot_shelf_region(data123);

if plotSouthAtlanticBight:
    data123 = get_processed_data_across_shelf([0,1,2], allData, [mf.area_SouthAtlanticBight], params);
    data789 = get_processed_data_across_shelf([6,7,8], allData, [mf.area_SouthAtlanticBight], params);

    distances = np.array(gridwiseData["distance"])[data123.regionMask]/1000.0; #in km
    dists, laruSouthAtlanticBightMask = mf.calculate_km_subsection_bounds_along_shelf(data123.xcoords, data123.ycoords, distances, mf.laruelle_SouthAtlanticBight, params, testPlot=plotTestplots);

    segmentDistances = np.cumsum(distances);
    xlabel = "distance along shelf (km, South to North)";
    outputFile = "plots/current_along_shelf/flux_along_shelf_SouthAtlanticBight.png"
    plot_shelf_section(segmentDistances, data123, data789, xlabel=xlabel, yRange=yRange, outputFile=outputFile, noClose=closePlots, laruelleSections=dists, laruelleSectionNames=["South Atlantic Bight"]);
    #plot_shelf_region(data123);

if plotSouthGreenland:
    data123 = get_processed_data_across_shelf([0,1,2], allData, [mf.area_SouthGreenland], params);
    data789 = get_processed_data_across_shelf([6,7,8], allData, [mf.area_SouthGreenland], params);
    distances = np.array(gridwiseData["distance"])[data123.regionMask]/1000.0; #in km

    data123, data789, distances = remove_section(data123, data789, distances, start=285, stop=None);
    #plot_shelf_region(data123);

    dists, laruSouthernGreenlandMask = mf.calculate_km_subsection_bounds_along_shelf(data123.xcoords, data123.ycoords, distances, mf.laruelle_SouthernGreenland, params, testPlot=plotTestplots);

    segmentDistances = np.cumsum(distances);
    xlabel = "distance along shelf (km, West to East)";
    outputFile = "plots/current_along_shelf/flux_along_shelf_SouthGreenland.png"
    plot_shelf_section(segmentDistances, data123, data789, xlabel=xlabel, yRange=yRange, outputFile=outputFile, noClose=closePlots, laruelleSections=dists, laruelleSectionNames=["Southern Greenland"]);

if plotCascianShelf:
    data123 = get_processed_data_across_shelf([0,1,2], allData, [mf.area_CascianShelf], params);
    data789 = get_processed_data_across_shelf([6,7,8], allData, [mf.area_CascianShelf], params);
    distances = np.array(gridwiseData["distance"])[data123.regionMask]/1000.0; #in km

    dists, laruCascadianShelfMask = mf.calculate_km_subsection_bounds_along_shelf(data123.xcoords, data123.ycoords, distances, mf.laruelle_CascadianShelf, params, testPlot=plotTestplots);

    segmentDistances = np.cumsum(distances);
    xlabel = "distance along shelf (km, North to South)";
    outputFile = "plots/current_along_shelf/flux_along_shelf_CascianShelf.png"
    plot_shelf_section(segmentDistances, data123, data789, xlabel=xlabel, yRange=yRange, outputFile=outputFile, noClose=closePlots, laruelleSections=dists, laruelleSectionNames=["Cascadian Shelf"]);
    #plot_shelf_region(data123);


if plotIrmingerSea:
    data123 = get_processed_data_across_shelf([0,1,2], allData, [mf.area_IrmingerSea], params);
    data789 = get_processed_data_across_shelf([6,7,8], allData, [mf.area_IrmingerSea], params);
    distances = np.array(gridwiseData["distance"])[data123.regionMask]/1000.0; #in km

    dists, laruIrmingerSeaMask = mf.calculate_km_subsection_bounds_along_shelf(data123.xcoords, data123.ycoords, distances, mf.laruelle_IrmingerSea, params, testPlot=plotTestplots);

    segmentDistances = np.cumsum(distances);
    xlabel = "distance along shelf (km, West to East)";
    outputFile = "plots/current_along_shelf/flux_along_shelf_IrmingerSea.png"
    plot_shelf_section(segmentDistances, data123, data789, xlabel=xlabel, yRange=yRange, outputFile=outputFile, noClose=closePlots, laruelleSections=dists, laruelleSectionNames=["Irminger Sea"]);
    #plot_shelf_region(data123);

if plot_manuscript:
    font = {'weight' : 'normal',
            'size'   : 18}
    matplotlib.rc('font', **font)
    fig = plt.figure(figsize=(15,15))
    gs = GridSpec(2,2, figure=fig, wspace=0.33,hspace=0.2,bottom=0.05,top=0.97,left=0.1,right=0.97)
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[0,1])

    data123 = get_processed_data_across_shelf([0,1,2], allData, [mf.area_EuropeanShelf], params);
    data789 = get_processed_data_across_shelf([6,7,8], allData, [mf.area_EuropeanShelf], params);
    distances = np.array(gridwiseData["distance"])[data123.regionMask]/1000.0; #in km

    data123, data789, distances = remove_section(data123, data789, distances, start=370, stop=None);
    #plot_shelf_region(data123);

    dists1, laruNorthSeaMask = mf.calculate_km_subsection_bounds_along_shelf(data123.xcoords, data123.ycoords, distances, mf.laruelle_NorthSea, params, testPlot=plotTestplots);
    dists2, laruEnglishChannelMask = mf.calculate_km_subsection_bounds_along_shelf(data123.xcoords, data123.ycoords, distances, mf.laruelle_EnglishChannel, params, testPlot=plotTestplots);

    segmentDistances = np.cumsum(distances);
    xlabel = "distance along shelf (km, North to South)";
    outputFile = "plots/current_along_shelf/flux_along_shelf_europeanShelf.png"
    plot_shelf_section(segmentDistances, data123, data789, xlabel=xlabel, yRange=yRange, outputFile=None, noClose=closePlots, laruelleSections=dists1+dists2, laruelleSectionNames=["North Sea", "English Channel"],ax = [ax1,ax2]);
    ax1.text(0.90,0.95,'(a)',transform=ax1.transAxes,va='top',fontsize=18,fontweight='bold')
    ax2.text(0.90,0.95,'(b)',transform=ax2.transAxes,va='top',fontsize=18,fontweight='bold')

    # ax3 = fig.add_subplot(gs[1,0])
    # ax4 = fig.add_subplot(gs[1,1])
    #
    # data123 = get_processed_data_across_shelf([0,1,2], allData, [mf.area_Tasmania], params);
    # data789 = get_processed_data_across_shelf([6,7,8], allData, [mf.area_Tasmania], params);
    #
    # distances = np.array(gridwiseData["distance"])[data123.regionMask]/1000.0; #in km
    # dists, laruTasmanianShelfMask = mf.calculate_km_subsection_bounds_along_shelf(data123.xcoords, data123.ycoords, distances, mf.laruelle_TasmanianShelf, params, testPlot=plotTestplots);
    #
    # segmentDistances = np.cumsum(distances);
    # xlabel = "distance along shelf (km, West to East)";
    # outputFile = "plots/current_along_shelf/flux_along_shelf_Tasmania.png"
    # plot_shelf_section(segmentDistances, data123, data789, xlabel=xlabel, yRange=yRange, outputFile=None, noClose=closePlots, laruelleSections=dists, laruelleSectionNames=["Tasmanian Shelf"],ax = [ax3,ax4]);
    # ax3.text(0.90,0.95,'(c)',transform=ax3.transAxes,va='top',fontsize=18,fontweight='bold')
    # ax4.text(0.90,0.95,'(d)',transform=ax4.transAxes,va='top',fontsize=18,fontweight='bold')

    ax5 = fig.add_subplot(gs[1,0])
    ax6 = fig.add_subplot(gs[1,1])

    data123 = get_processed_data_across_shelf([0,1,2], allData, [mf.area_MidAtlanticBight], params);
    data789 = get_processed_data_across_shelf([6,7,8], allData, [mf.area_MidAtlanticBight], params);

    distances = np.array(gridwiseData["distance"])[data123.regionMask]/1000.0; #in km
    dists, laruMidAtlanticBightMask = mf.calculate_km_subsection_bounds_along_shelf(data123.xcoords, data123.ycoords, distances, mf.laruelle_MidAtlanticBight, params, testPlot=plotTestplots);

    segmentDistances = np.cumsum(distances);
    xlabel = "distance along shelf (km, South to North)";
    outputFile = "plots/current_along_shelf/flux_along_shelf_midAtlanticBight.png"
    plot_shelf_section(segmentDistances, data123, data789, xlabel=xlabel, yRange=yRange, outputFile=None, noClose=closePlots, laruelleSections=dists, laruelleSectionNames=["Mid-Atlantic Bight"],ax = [ax5,ax6]);
    ax5.text(0.90,0.95,'(c)',transform=ax5.transAxes,va='top',fontsize=18,fontweight='bold')
    ax6.text(0.90,0.95,'(d)',transform=ax6.transAxes,va='top',fontsize=18,fontweight='bold')
    fig.savefig("plots/current_along_shelf/manuscript_figure.png",dpi=300);
