#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 13:43:32 2018

@author: rr
Modified by Daniel J. Ford (d.ford@exeter.ac.uk)
Date: 03/2023
Changes:
- Added weighted mean and std as function, to remove need for DescrStatsW (cannot get to load for Python 2)
- Added absolute eks, stokes and geo currents to the output.
- Added annual and seasonal means to output table.
"""

import cPickle as pickle;
import numpy as np;
import matplotlib.pyplot as plt;
import parameter_sets as ps;
import skim_utilities as su;
from netCDF4 import Dataset;
import mask_functions as mf;
from os import path;
import pandas as pd;
import os; #mkdir


def get_month_indices(months, numYears):
    indices = [];
    for year in range(numYears):
        for month in months:
            indices.append(year*12 + month);
    return indices;


def get_temporal_index(year, imonth, startYear=1993):
    return (year-startYear)*12 + imonth;

#Select and return data for given set of months each year
#Returned data is reformatted into an array and stored in a QuickStruct, and contains:
#   overall mean and stddev of ekman, geostrophic and stokes current for the whole period.
#months defines the months of the year to use (e.g. 0=Jan, 11=Dec)
def get_processed_data_across_shelf(months, data, regionMaskBoundsList, params, filterStokes=True, gasTransferParameterisation=False):
    if len(data) % 12 != 0:
        raise ValueError ("get_processed_data only supports data for whole years (multiples of 12 months)");

    #Create the indices to use
    numYears = int(len(data)/12);
    indices = get_month_indices(months, numYears);

    #Create arrays to store data for ekman, geostrophic and stokes data
    length = 0;
    for m in months:
        #print len(data[m]);
        if len(data[m]) > length:
            length = len(data[m]);

    ekmanVals = np.empty((length, len(indices))); ekmanVals[:] = np.nan;
    geostrophicVals = np.empty((length, len(indices))); geostrophicVals[:] = np.nan;
    stokesVals = np.empty((length, len(indices))); stokesVals[:] = np.nan;
    stokesMaskPass = np.empty((length, len(indices))); stokesMaskPass[:] = np.nan;
    xcoords = np.empty((length, len(indices))); xcoords[:] = np.nan; #xcoords and ycoords are just used for checking region mask is correct
    ycoords = np.empty((length, len(indices))); ycoords[:] = np.nan;
    if gasTransferParameterisation:
        kVals = np.empty((length, len(indices))); kVals[:] = np.nan;

    #Copy data:
    for i, index in enumerate(indices):
        ekmanVals[:,i] = [curData.nEkmanAcrossShelf for curData in data[index]];
        geostrophicVals[:,i] = [curData.nGeostrophicAcrossShelf for curData in data[index]];
        stokesVals[:,i] = [curData.nStokesAcrossShelf for curData in data[index]];
        stokesMaskPass[:,i] = [curData.stokesMaskPass for curData in data[index]];
        if gasTransferParameterisation:
            kVals[:,i] = [curData.k for curData in data[index]];

        xcoords[:,i] = [curData.indexX for curData in data[index]];
        ycoords[:,i] = [curData.indexY for curData in data[index]];

    #Wherever stokesMaskPass fails, set stokes to 0
    if (filterStokes):
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

    if gasTransferParameterisation:
        results.k = np.nanmean(kVals, axis=1);
        results.kSD = np.nanstd(kVals, axis=1);

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

    results.ekmanVals = ekmanVals[regionMask,:];
    results.geostrophicVals = geostrophicVals[regionMask,:];
    results.stokesVals = stokesVals[regionMask,:];

    if gasTransferParameterisation:
        results.k = results.k[regionMask];
        results.kSD = results.kSD[regionMask];
        results.kVals = kVals[regionMask,:];

    return results;

#Select and return data for given set of months each year
#Returned data is reformatted into an array and stored in a QuickStruct, and contains:
#   overall mean and stddev of ekman, geostrophic and stokes current for the whole period.
#months defines the months of the year to use (e.g. 0=Jan, 11=Dec)
def get_processed_data_across_temporal_mask(temporalMask, data, regionMaskBoundsList, params, filterStokes=True):
    if len(data) % 12 != 0:
        raise ValueError ("get_processed_data only supports data for whole years (multiples of 12 months)");

    #Create the indices to use
    indices = np.array(temporalMask);

    #Create arrays to store data for ekman, geostrophic and stokes data
    length = 0;
    for m in indices:
        #print len(data[m]);
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
    if (filterStokes):
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

    results.ekmanVals = ekmanVals[regionMask,:];
    results.geostrophicVals = geostrophicVals[regionMask,:];
    results.stokesVals = stokesVals[regionMask,:];

    return results;

def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    # Fast and numerically precise:
    variance = np.average((values-average)**2, weights=weights)
    return (average, np.sqrt(variance))

# from statsmodels.stats.weightstats import DescrStatsW;
# a = np.array([1,2,1,2,1,2,1,2,2,2,1,3,2,3]);
# b = np.array([1,1,1,100,1,1,1,1,1,1,1,1,1,1]);



def print_statistics(data, distances, name, mask, gasTransferParameterisation=False):

    if len(data.xcoords) == 0 or len(mask) == 0:
        print name, "- No data available!\n";
        return (name, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan);

    eks = data.ekmanVals[mask,:];
    geo = data.geostrophicVals[mask,:];
    stokes = data.stokesVals[mask,:];

    eks = data.ekmanAcrossShelf[mask]
    geo = data.geostrophicAcrossShelf[mask]
    stokes = data.stokesAcrossShelf[mask]

    #distances = np.array(distances);
    if (len(distances) != len(eks)):
        raise ValueError("distances must be the same length as data values");

    #weights = np.transpose(np.array([distances]*eks.shape[1]));
    weights = np.array(distances)

    if gasTransferParameterisation:
        #ks = data.kVals[mask,:];
        ks = np.nanmean(data.kVals,axis=1)
        ks = ks[mask]
        notnan = np.where(np.isnan(ks)==False);
        kMean, kSD = weighted_avg_and_std(ks[notnan], weights=weights[notnan])
        kmean_n = len(notnan[0])
        # kStats = DescrStatsW(ks[notnan], weights=weights[notnan], ddof=0);
        # kMean = kStats.mean;
        # kSD = kStats.std;

        print name, "Gas transfer velocity stats";
        print "\tk: ", kMean, "+/-", kSD;


    totals = np.abs(eks) + np.abs(geo) + np.abs(stokes);
    eksProps = np.abs(eks) / totals;
    geoProps = np.abs(geo) / totals;
    stokesProps = np.abs(stokes) / totals;

    notnan = np.where(np.isnan(eksProps)==False);
    meanEks, stdEks = weighted_avg_and_std(eks[notnan], weights=weights[notnan])
    #eksStats = DescrStatsW(eksProps[notnan], weights=weights[notnan], ddof=0);
    # eksMean = eksStats.mean;
    # eksSD = eksStats.std;
    eksMean,eksSD = weighted_avg_and_std(eksProps[notnan], weights=weights[notnan])
    eksMean_n = len(notnan[0])

    notnan = np.where(np.isnan(geoProps)==False);
    # geoStats = DescrStatsW(geoProps[notnan], weights=weights[notnan], ddof=0);
    # geoMean = geoStats.mean;
    # geoSD = geoStats.std;
    meanGeo, stdGeo = weighted_avg_and_std(geo[notnan], weights=weights[notnan])
    geoMean, geoSD =weighted_avg_and_std(geoProps[notnan], weights=weights[notnan])
    geoMean_n = len(notnan[0])

    notnan = np.where(np.isnan(stokesProps)==False);
    # stokesStats = DescrStatsW(stokesProps[notnan], weights=weights[notnan], ddof=0);
    # stokesMean = stokesStats.mean;
    # stokesSD = stokesStats.std;
    meanStokes, stdStokes = weighted_avg_and_std(stokes[notnan], weights=weights[notnan])
    stokesMean, stokesSD =weighted_avg_and_std(stokesProps[notnan], weights=weights[notnan])
    stokesMean_n = len(notnan[0])
    print name, "mean absolute"
    print "\tEkman: ", meanEks, "+/-", stdEks;
    print "\tGeostrophic: ", meanGeo, "+/-", stdGeo;
    print "\tStokes: ", meanStokes, "+/-", stdStokes;
    print "\t", "total:", meanEks+meanGeo+meanStokes;

    print name, "mean proportions"
    print "\tEkman: ", eksMean, "+/-", eksSD;
    print "\tGeostrophic: ", geoMean, "+/-", geoSD;
    print "\tStokes: ", stokesMean, "+/-", stokesSD;
    print "\t", "total:", eksMean+geoMean+stokesMean;

    totalOntoShelf = eks+geo+stokes;
    notnan = np.where(np.isnan(totalOntoShelf)==False);
    # totalOntoShelfStats = DescrStatsW(totalOntoShelf[notnan], weights=weights[notnan], ddof=0);
    # totalMean = totalOntoShelfStats.mean;
    # totalSD = totalOntoShelfStats.std;
    totalMean, totalSD =weighted_avg_and_std(totalOntoShelf[notnan], weights=weights[notnan])
    totalMean_n = len(notnan[0])

    eksPercent = eks / totalOntoShelf * 100.0;
    geoPercent = geo / totalOntoShelf * 100.0;
    stokesPercent = stokes / totalOntoShelf * 100.0;

    notnan = np.where(np.isnan(eksPercent)==False);
    # eksPercentStats = DescrStatsW(eksPercent[notnan], weights=weights[notnan], ddof=0);
    # eksPercentMean = eksPercentStats.mean;
    # eksPercentSD = eksPercentStats.std;
    eksPercentMean, eksPercentSD =weighted_avg_and_std(eksPercent[notnan], weights=weights[notnan])
    eksPercent_n = len(notnan[0])

    notnan = np.where(np.isnan(geoPercent)==False);
    # geoPercentStats = DescrStatsW(geoPercent[notnan], weights=weights[notnan], ddof=0);
    # geoPercentMean = geoPercentStats.mean;
    # geoPercentSD = geoPercentStats.std;
    geoPercentMean, geoPercentSD =weighted_avg_and_std(geoPercent[notnan], weights=weights[notnan])
    geoPercent_n = len(notnan[0])

    notnan = np.where(np.isnan(stokesPercent)==False);
    # stokesPercentStats = DescrStatsW(stokesPercent[notnan], weights=weights[notnan], ddof=0);
    # stokesPercentMean = stokesPercentStats.mean;
    # stokesPercentSD = stokesPercentStats.std;
    stokesPercentMean, stokesPercentSD =weighted_avg_and_std(stokesPercent[notnan], weights=weights[notnan])
    stokesPercent_n = len(notnan[0])

    print name, "percentage total onto-shelf current";
    print "\tEkman: ", eksPercentMean, "+/-", eksPercentSD;
    print "\tGeostrophic: ", geoPercentMean, "+/-", geoPercentSD;
    print "\tStokes: ", stokesPercentMean, "+/-", stokesPercentSD;
    print "\t", "total onto-shelf current (m/s):", totalMean, "+/-", totalSD;
    print "";

    if gasTransferParameterisation:
        return (name, totalMean, totalSD, totalMean_n, meanEks, stdEks, eksMean_n, meanGeo,stdGeo, geoMean_n, meanStokes, stdStokes, stokesMean_n, eksMean, eksSD,eksMean_n, geoMean, geoSD, geoMean_n, stokesMean,
            stokesSD, stokesMean_n, eksPercentMean, eksPercentSD, eksMean_n, geoPercentMean, geoPercentSD, geoMean_n, stokesPercentMean, stokesPercentSD, stokesMean_n, kMean, kSD, kmean_n);
    else:
        return (name, totalMean, totalSD, totalMean_n, meanEks, stdEks, eksMean_n, meanGeo,stdGeo, geoMean_n, meanStokes, stdStokes, stokesMean_n, eksMean, eksSD,eksMean_n, geoMean, geoSD, geoMean_n, stokesMean,
            stokesSD, stokesMean_n, eksPercentMean, eksPercentSD, eksMean_n, geoPercentMean, geoPercentSD, geoMean_n, stokesPercentMean, stokesPercentSD, stokesMean_n, np.nan, np.nan,np.nan);



def calc_mean_data(params, calculateTable1=False, calculateTable2=False, calculateTable2GasTransferVelocity=True, verbose=True, outputPath=None):
    #Read data
    inputDataPath = path.join("E:/SKIM"); #params.paramsetName);
    gridwiseData = pd.read_table(path.join(inputDataPath, "gridwise_data/per_grid_cell_edge_data_"+params.paramsetName+".csv"), sep=',');
    gridwiseData.x = gridwiseData.x.astype(int);
    gridwiseData.y = gridwiseData.y.astype(int);
    allData = pickle.load(open(path.join(inputDataPath, "current_data/surface_currents_"+params.paramsetName+".p"), "rb"));

    #create path
    if path.exists(path.join(inputDataPath, "output_means")) == False:
        os.mkdir(path.join(inputDataPath, "output_means"));

    if outputPath == None:
        outputPath = path.join("data", params.paramsetName);

    if calculateTable1 == True:
        if verbose:
            print "Calculating Table 1 data for", params.paramsetName;
        t1data = [];

        ##########################
        #Painter Scottish Hebrides (2014, September)
        temporalMask = [get_temporal_index(2014, 8)];
        painterData = get_processed_data_across_temporal_mask(temporalMask, allData, [mf.verif_PainterHebredes], params);
        distances = np.array(gridwiseData["distance"])[painterData.regionMask]/1000.0; #in km
        dists, painterMask = mf.calculate_km_subsection_bounds_along_shelf(painterData.xcoords, painterData.ycoords, distances, mf.verif_PainterHebredes, params, testPlot=False);
        t1data.append(print_statistics(painterData, distances, "Painter Hebrides (Sept, 2014)", painterMask));


        ##########################
        #Yuan South Atlantic Bight (Annual analysis 2002-2008, total mean, winter(jan-mar) and summer (jul-sept))
        #total mean (across all times)
#        annualMask = range(0, len(allData));
#        yuanDataTotal = get_processed_data_across_temporal_mask(annualMask, allData, [mf.verif_YuanSouthAtlanticBight], params);
#        distances = np.array(gridwiseData["distance"])[yuanDataTotal.regionMask]/1000.0; #in km
#        dists, yuanDataTotalMask = mf.calculate_km_subsection_bounds_along_shelf(yuanDataTotal.xcoords, yuanDataTotal.ycoords, distances, mf.verif_YuanSouthAtlanticBight, params, testPlot=False);
#        t1data.append(print_statistics(yuanDataTotal, "Yuan SouthAtlanticBight (year round, 2010-2016)", yuanDataTotalMask));
#
        #winter mean (Jan-Mar)
        winterTemporalMask = [];
        for y in range(2002, 2014+1):
            winterTemporalMask.append(get_temporal_index(y, 0));
            winterTemporalMask.append(get_temporal_index(y, 1));
            winterTemporalMask.append(get_temporal_index(y, 2));
        yuanDataWinter = get_processed_data_across_temporal_mask(winterTemporalMask, allData, [mf.verif_YuanSouthAtlanticBight], params);
        distances = np.array(gridwiseData["distance"])[yuanDataWinter.regionMask]/1000.0; #in km
        dists, yuanDataWinterMask = mf.calculate_km_subsection_bounds_along_shelf(yuanDataWinter.xcoords, yuanDataWinter.ycoords, distances, mf.verif_YuanSouthAtlanticBight, params, testPlot=False);
        t1data.append(print_statistics(yuanDataWinter, distances, "Yuan SouthAtlanticBight (winter, Jan-Mar 2002-2014)", yuanDataWinterMask));

        #summer mean (Jul-Sept)
        summerTemporalMask = [];
        for y in range(2002, 2014+1):
            summerTemporalMask.append(get_temporal_index(y, 6));
            summerTemporalMask.append(get_temporal_index(y, 7));
            summerTemporalMask.append(get_temporal_index(y, 8));
        yuanDataSummer = get_processed_data_across_temporal_mask(summerTemporalMask, allData, [mf.verif_YuanSouthAtlanticBight], params);
        distances = np.array(gridwiseData["distance"])[yuanDataSummer.regionMask]/1000.0; #in km
        dists, yuanDataSummerMask = mf.calculate_km_subsection_bounds_along_shelf(yuanDataSummer.xcoords, yuanDataSummer.ycoords, distances, mf.verif_YuanSouthAtlanticBight, params, testPlot=False);
        t1data.append(print_statistics(yuanDataSummer, distances, "Yuan SouthAtlanticBight (summer, Jul-Sept 2002-2014)", yuanDataSummerMask));


        #############################
        #Fewings Mid Atlantic Bight (Summer (Jan-Mar) vs Winter(Jul-Sep))
        #Annual mean (Jan-Dec)
        temporalMask = range(get_temporal_index(2001, 5), get_temporal_index(2007, 4)+1); #June 2001 to May 2007
        data = get_processed_data_across_temporal_mask(temporalMask, allData, [mf.verif_FewingsMidAtlanticBight], params);
        distances = np.array(gridwiseData["distance"])[data.regionMask]/1000.0; #in km
        dists, dataMask = mf.calculate_km_subsection_bounds_along_shelf(data.xcoords, data.ycoords, distances, mf.verif_FewingsMidAtlanticBight, params, testPlot=False);
        t1data.append(print_statistics(data, distances, "Fewings MidAtlanticBight (annual, June 2001-May 2007)", dataMask));

        #Winter mean (October-March)
        winterTemporalMask = [];
        for y in range(2002, 2007+1):
            if y < 2007: #only goes to May 2007
                winterTemporalMask.append(get_temporal_index(y, 9)); #oct
                winterTemporalMask.append(get_temporal_index(y, 10));
                winterTemporalMask.append(get_temporal_index(y, 11));
            winterTemporalMask.append(get_temporal_index(y, 0)); #jan
            winterTemporalMask.append(get_temporal_index(y, 1));
            winterTemporalMask.append(get_temporal_index(y, 2));
        data = get_processed_data_across_temporal_mask(winterTemporalMask, allData, [mf.verif_FewingsMidAtlanticBight], params);
        distances = np.array(gridwiseData["distance"])[data.regionMask]/1000.0; #in km
        dists, dataMask = mf.calculate_km_subsection_bounds_along_shelf(data.xcoords, data.ycoords, distances, mf.verif_FewingsMidAtlanticBight, params, testPlot=False);
        t1data.append(print_statistics(data, distances, "Fewings MidAtlanticBight (winter, Oct-Mar 2002-2007)", dataMask));

        #Summer mean (April-September)
        summerTemporalMask = [];
        summerTemporalMask.append(get_temporal_index(2001, 2));
        for y in range(2002, 2007+1):
            if y > 2002: #Starts in June 2002
                summerTemporalMask.append(get_temporal_index(y, 3)); #April
                summerTemporalMask.append(get_temporal_index(y, 4));
            if y < 2007: #Ends in May 2007
                summerTemporalMask.append(get_temporal_index(y, 5)); #June
                summerTemporalMask.append(get_temporal_index(y, 6));
                summerTemporalMask.append(get_temporal_index(y, 7));
                summerTemporalMask.append(get_temporal_index(y, 8));
        data = get_processed_data_across_temporal_mask(summerTemporalMask, allData, [mf.verif_FewingsMidAtlanticBight], params);
        distances = np.array(gridwiseData["distance"])[data.regionMask]/1000.0; #in km
        dists, dataMask = mf.calculate_km_subsection_bounds_along_shelf(data.xcoords, data.ycoords, distances, mf.verif_FewingsMidAtlanticBight, params, testPlot=False);
        t1data.append(print_statistics(data, distances, "Fewings MidAtlanticBight (summer, Apr-Sept 2001-2007)", dataMask));

        ###########################
        # Woodson, californian coast (upwelling (April-Sept) and non-upwelling seasons (October-March))
        #upwelling season (April-Sept)
        upwellingTemporalMask = [];
        for y in range(2004, 2009+1):
            upwellingTemporalMask.append(get_temporal_index(y, 3)); #April
            upwellingTemporalMask.append(get_temporal_index(y, 4));
            upwellingTemporalMask.append(get_temporal_index(y, 5));
            upwellingTemporalMask.append(get_temporal_index(y, 6));
            upwellingTemporalMask.append(get_temporal_index(y, 7));
            upwellingTemporalMask.append(get_temporal_index(y, 8)); #September

        data = get_processed_data_across_temporal_mask(upwellingTemporalMask, allData, [mf.verif_WoodsonCaliforniaCoast], params);
        distances = np.array(gridwiseData["distance"])[data.regionMask]/1000.0; #in km
        dists, dataMask = mf.calculate_km_subsection_bounds_along_shelf(data.xcoords, data.ycoords, distances, mf.verif_WoodsonCaliforniaCoast, params, testPlot=False);
        t1data.append(print_statistics(data, distances, "Woodson California coast (upwelling season, Apr-Sept, 2004-2009)", dataMask));

        #Non-upwelling season (Oct-Mar)
        nonupwellingTemporalMask = [];
        for y in range(2004, 2009+1):
            nonupwellingTemporalMask.append(get_temporal_index(y, 0)); #Jan
            nonupwellingTemporalMask.append(get_temporal_index(y, 1));
            nonupwellingTemporalMask.append(get_temporal_index(y, 2)); #Mar
            nonupwellingTemporalMask.append(get_temporal_index(y, 9)); #Oct
            nonupwellingTemporalMask.append(get_temporal_index(y, 10));
            nonupwellingTemporalMask.append(get_temporal_index(y, 11)); #Dec

        data = get_processed_data_across_temporal_mask(nonupwellingTemporalMask, allData, [mf.verif_WoodsonCaliforniaCoast], params);
        distances = np.array(gridwiseData["distance"])[data.regionMask]/1000.0; #in km
        dists, dataMask = mf.calculate_km_subsection_bounds_along_shelf(data.xcoords, data.ycoords, distances, mf.verif_WoodsonCaliforniaCoast, params, testPlot=False);
        t1data.append(print_statistics(data, distances, "Woodson California coast (non-upwelling season, Oct-March, 2004-2009)", dataMask));

        ###########################
        #Waite Eastern Indian Ocean (May averaged over 2010-2016)
        mayTemporalMask = [];
        mayTemporalMask.append(get_temporal_index(2006, 4));

        data = get_processed_data_across_temporal_mask(mayTemporalMask, allData, [mf.verif_WaiteEasternIndianOcean], params);
        distances = np.array(gridwiseData["distance"])[data.regionMask]/1000.0; #in km
        dists, dataMask = mf.calculate_km_subsection_bounds_along_shelf(data.xcoords, data.ycoords, distances, mf.verif_WaiteEasternIndianOcean, params, testPlot=False);
        t1data.append(print_statistics(data, distances, "Waite Eastern Indian Ocean (May 2006)", dataMask));


        #########################
        #Wei, East China Sea, Annual, Summer () and Winter...
        #Annual
        annualMask = [];
        annualMask = range(get_temporal_index(1997, 0), get_temporal_index(2010, 0)+1);

        data = get_processed_data_across_temporal_mask(annualMask, allData, [mf.verif_WieEastChinaSea], params);
        distances = np.array(gridwiseData["distance"])[data.regionMask]/1000.0; #in km
        dists, dataMask = mf.calculate_km_subsection_bounds_along_shelf(data.xcoords, data.ycoords, distances, mf.verif_WieEastChinaSea, params, testPlot=False);
        t1data.append(print_statistics(data, distances, "Wei East China Sea (annual, Jan 1993 - Jan 2010)", dataMask));

#        #Winter (Jan-Mar)
#        winderTemporalMask = [];
#        for y in range(2010, 2016+1):
#            winderTemporalMask.append(get_temporal_index(y, 0));
#            winderTemporalMask.append(get_temporal_index(y, 1));
#            winderTemporalMask.append(get_temporal_index(y, 2));
#
#        data = get_processed_data_across_temporal_mask(winderTemporalMask, allData, [mf.verif_WieEastChinaSea], params);
#        distances = np.array(gridwiseData["distance"])[data.regionMask]/1000.0; #in km
#        dists, dataMask = mf.calculate_km_subsection_bounds_along_shelf(data.xcoords, data.ycoords, distances, mf.verif_WieEastChinaSea, params, testPlot=False);
#        t1data.append(print_statistics(data, distances, "Wei East China Sea (winter, Jan-Mar 2010-2016)", dataMask));
#
#
#        #Summer (Jul-Sept)
#        summerTemporalMask = [];
#        for y in range(2010, 2016+1):
#            summerTemporalMask.append(get_temporal_index(y, 6));
#            summerTemporalMask.append(get_temporal_index(y, 7));
#            summerTemporalMask.append(get_temporal_index(y, 8));
#
#        data = get_processed_data_across_temporal_mask(summerTemporalMask, allData, [mf.verif_WieEastChinaSea], params);
#        distances = np.array(gridwiseData["distance"])[data.regionMask]/1000.0; #in km
#        dists, dataMask = mf.calculate_km_subsection_bounds_along_shelf(data.xcoords, data.ycoords, distances, mf.verif_WieEastChinaSea, params, testPlot=False);
#        t1data.append(print_statistics(data, distances, "Wei East China Sea (Summer, Jul-Sept 2010-2016)", dataMask));



        #Write to file
        dft1 = pd.DataFrame(t1data);
        dft1.columns = ["region", "total across-shelf", "total across-shelf SD", "proportionEks", "proportionEksSD", "proportionGeo", "proportionGeoSD", "proportionStokes", "proportionStokesSD", "percentEks", "percentEksSD", "percentGeo", "percentGeoSD", "percentStokes", "percentStokesSD"];
        dft1.to_csv(path.join(outputPath, "output_means", "table1data_"+params.paramsetName+".csv"), index=False);

    if calculateTable2 == True:
        if verbose:
            print "Calculating Table 1 data for", params.paramsetName;

        t2data = [];

        #European shelf
        data_all =  get_processed_data_across_shelf([0,1,2,3,4,5,6,7,8,9,10,11], allData, [mf.area_EuropeanShelf], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        data123 = get_processed_data_across_shelf([0,1,2], allData, [mf.area_EuropeanShelf], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        data456 = get_processed_data_across_shelf([3,4,5], allData, [mf.area_EuropeanShelf], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        data789 = get_processed_data_across_shelf([6,7,8], allData, [mf.area_EuropeanShelf], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        data101112 = get_processed_data_across_shelf([9,10,11], allData, [mf.area_EuropeanShelf], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        distances = np.array(gridwiseData["distance"])[data123.regionMask]/1000.0; #in km
        dists, laruNorthSeaMask = mf.calculate_km_subsection_bounds_along_shelf(data123.xcoords, data123.ycoords, distances, mf.laruelle_NorthSea, params, testPlot=False);
        dists, laruEnglishChannelMask = mf.calculate_km_subsection_bounds_along_shelf(data123.xcoords, data123.ycoords, distances, mf.laruelle_EnglishChannel, params, testPlot=False);
        t2data.append(print_statistics(data_all, distances[laruNorthSeaMask], "North Sea (Annual)", laruNorthSeaMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        t2data.append(print_statistics(data123, distances[laruNorthSeaMask], "North Sea (Jan-Mar)", laruNorthSeaMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        t2data.append(print_statistics(data456, distances[laruNorthSeaMask], "North Sea (Apr-Jun)", laruNorthSeaMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        t2data.append(print_statistics(data789, distances[laruNorthSeaMask], "North Sea (Jul-Sep)", laruNorthSeaMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        t2data.append(print_statistics(data101112, distances[laruNorthSeaMask], "North Sea (Oct-Dec)", laruNorthSeaMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));

        t2data.append(print_statistics(data_all, distances[laruEnglishChannelMask], "English Channel (Annual)", laruEnglishChannelMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        t2data.append(print_statistics(data123, distances[laruEnglishChannelMask], "English Channel (Jan-Mar)", laruEnglishChannelMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        t2data.append(print_statistics(data456, distances[laruEnglishChannelMask], "English Channel (Apr-Jun)", laruEnglishChannelMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        t2data.append(print_statistics(data789, distances[laruEnglishChannelMask], "English Channel (Jul-Sep)", laruEnglishChannelMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        t2data.append(print_statistics(data101112, distances[laruEnglishChannelMask], "English Channel (Oct-Dec)", laruEnglishChannelMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        # t2data.append(print_statistics(data123, distances[laruEnglishChannelMask], "English Channel (Jan-Mar)", laruEnglishChannelMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        #
        # t2data.append(print_statistics(data789, distances[laruEnglishChannelMask], "English Channel (Jul-Sep)", laruEnglishChannelMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));

        #MidAtlanticBight
        data_all =  get_processed_data_across_shelf([0,1,2,3,4,5,6,7,8,9,10,11], allData, [mf.area_MidAtlanticBight], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        data123 = get_processed_data_across_shelf([0,1,2], allData, [mf.area_MidAtlanticBight], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        data456 = get_processed_data_across_shelf([3,4,5], allData, [mf.area_MidAtlanticBight], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        data789 = get_processed_data_across_shelf([6,7,8], allData, [mf.area_MidAtlanticBight], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        data101112 = get_processed_data_across_shelf([9,10,11], allData, [mf.area_MidAtlanticBight], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        # data123 = get_processed_data_across_shelf([0,1,2], allData, [mf.area_MidAtlanticBight], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        # data789 = get_processed_data_across_shelf([6,7,8], allData, [mf.area_MidAtlanticBight], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);

        distances = np.array(gridwiseData["distance"])[data123.regionMask]/1000.0; #in km
        dists, laruMidAtlanticBightMask = mf.calculate_km_subsection_bounds_along_shelf(data123.xcoords, data123.ycoords, distances, mf.laruelle_MidAtlanticBight, params, testPlot=False);
        t2data.append(print_statistics(data_all, distances[laruMidAtlanticBightMask], "Mid Atlantic Bight (Annual)", laruMidAtlanticBightMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        t2data.append(print_statistics(data123, distances[laruMidAtlanticBightMask], "Mid Atlantic Bight (Jan-Mar)", laruMidAtlanticBightMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        t2data.append(print_statistics(data456, distances[laruMidAtlanticBightMask], "Mid Atlantic Bight (Apr-Jun)", laruMidAtlanticBightMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        t2data.append(print_statistics(data789, distances[laruMidAtlanticBightMask], "Mid Atlantic Bight (Jul-Sep)", laruMidAtlanticBightMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        t2data.append(print_statistics(data101112, distances[laruMidAtlanticBightMask], "Mid Atlantic Bight (Oct-Dec)", laruMidAtlanticBightMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));

        # t2data.append(print_statistics(data123, distances[laruMidAtlanticBightMask], "Mid Atlantic Bight (Jan-Mar)", laruMidAtlanticBightMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        # t2data.append(print_statistics(data789, distances[laruMidAtlanticBightMask], "Mid Atlantic Bight (Jul-Sep)", laruMidAtlanticBightMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));

        #CoastOfJapan
        data_all =  get_processed_data_across_shelf([0,1,2,3,4,5,6,7,8,9,10,11], allData, [mf.area_CoastOfJapan], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        data123 = get_processed_data_across_shelf([0,1,2], allData, [mf.area_CoastOfJapan], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        data456 = get_processed_data_across_shelf([3,4,5], allData, [mf.area_CoastOfJapan], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        data789 = get_processed_data_across_shelf([6,7,8], allData, [mf.area_CoastOfJapan], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        data101112 = get_processed_data_across_shelf([9,10,11], allData, [mf.area_CoastOfJapan], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        # data123 = get_processed_data_across_shelf([0,1,2], allData, [mf.area_CoastOfJapan], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        # data789 = get_processed_data_across_shelf([6,7,8], allData, [mf.area_CoastOfJapan], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);

        distances = np.array(gridwiseData["distance"])[data123.regionMask]/1000.0; #in km
        def roll_data(data,roll):
            data.xcoords = np.roll(data.xcoords, roll);
            data.ycoords = np.roll(data.ycoords, roll);
            data.regionMask = np.roll(data.regionMask, roll);
            data.ekmanAcrossShelf = np.roll(data.ekmanAcrossShelf, roll);
            data.ekmanAcrossShelfSD = np.roll(data.ekmanAcrossShelfSD, roll);
            data.geostrophicAcrossShelf = np.roll(data.geostrophicAcrossShelf, roll);
            data.geostrophicAcrossShelfSD = np.roll(data.geostrophicAcrossShelfSD, roll);
            data.stokesAcrossShelf = np.roll(data.stokesAcrossShelf, roll);
            data.stokesAcrossShelfSD = np.roll(data.stokesAcrossShelfSD, roll);
            if calculateTable2GasTransferVelocity==True:
                data.k = np.roll(data.k, roll);
            return data

        #Custom roll for Japanese coast (to eliminate split)
        rollAmount = -110
        data_all= roll_data(data_all,rollAmount)
        data123 = roll_data(data123,rollAmount)
        data456 = roll_data(data456,rollAmount)
        data789 = roll_data(data789,rollAmount)
        data101112 = roll_data(data101112,rollAmount)
        distances = np.roll(distances, rollAmount);

        dists, laruCoastOfJapanMask = mf.calculate_km_subsection_bounds_along_shelf(data123.xcoords, data123.ycoords, distances, mf.laruelle_CoastOfJapan, params, testPlot=False);
        t2data.append(print_statistics(data_all, distances[laruCoastOfJapanMask], "Coast of Japan (Annual)", laruCoastOfJapanMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        t2data.append(print_statistics(data123, distances[laruCoastOfJapanMask], "Coast of Japan (Jan-Mar)", laruCoastOfJapanMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        t2data.append(print_statistics(data456, distances[laruCoastOfJapanMask], "Coast of Japan (Apr-Jun)", laruCoastOfJapanMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        t2data.append(print_statistics(data789, distances[laruCoastOfJapanMask], "Coast of Japan (Jul-Sep)", laruCoastOfJapanMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        t2data.append(print_statistics(data101112, distances[laruCoastOfJapanMask], "Coast of Japan (Oct-Dec)", laruCoastOfJapanMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        # t2data.append(print_statistics(data123, distances[laruCoastOfJapanMask], "Coast of Japan (Jan-Mar)", laruCoastOfJapanMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        # t2data.append(print_statistics(data789, distances[laruCoastOfJapanMask], "Coast of Japan (Jul-Sep)", laruCoastOfJapanMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));

        #Patagonian shelf
        data_all =  get_processed_data_across_shelf([0,1,2,3,4,5,6,7,8,9,10,11], allData, [mf.area_Patagonia], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        data123 = get_processed_data_across_shelf([0,1,2], allData, [mf.area_Patagonia], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        data456 = get_processed_data_across_shelf([3,4,5], allData, [mf.area_Patagonia], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        data789 = get_processed_data_across_shelf([6,7,8], allData, [mf.area_Patagonia], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        data101112 = get_processed_data_across_shelf([9,10,11], allData, [mf.area_Patagonia], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        # data123 = get_processed_data_across_shelf([0,1,2], allData, [mf.area_Patagonia], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        # data789 = get_processed_data_across_shelf([6,7,8], allData, [mf.area_Patagonia], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);

        distances = np.array(gridwiseData["distance"])[data123.regionMask]/1000.0; #in km
        dists, laruPatagonianShelfMask = mf.calculate_km_subsection_bounds_along_shelf(data123.xcoords, data123.ycoords, distances, mf.laruelle_PatagonianShelf, params, testPlot=False);
        t2data.append(print_statistics(data_all, distances[laruPatagonianShelfMask], "Patagonian Shelf (Annual)", laruPatagonianShelfMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        t2data.append(print_statistics(data123, distances[laruPatagonianShelfMask], "Patagonian Shelf (Jan-Mar)", laruPatagonianShelfMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        t2data.append(print_statistics(data456, distances[laruPatagonianShelfMask], "Patagonian Shelf (Apr-Jun)", laruPatagonianShelfMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        t2data.append(print_statistics(data789, distances[laruPatagonianShelfMask], "Patagonian Shelf (Jul-Sep)", laruPatagonianShelfMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        t2data.append(print_statistics(data101112, distances[laruPatagonianShelfMask], "Patagonian Shelf (Oct-Dec)", laruPatagonianShelfMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        # t2data.append(print_statistics(data123, distances[laruPatagonianShelfMask], "Patagonian Shelf (Jan-Mar)", laruPatagonianShelfMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        # t2data.append(print_statistics(data789, distances[laruPatagonianShelfMask], "Patagonian Shelf (Jul-Sep)", laruPatagonianShelfMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));


        #Bering Sea
        data_all =  get_processed_data_across_shelf([0,1,2,3,4,5,6,7,8,9,10,11], allData, [mf.area_BeringSeaWest], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        data123 = get_processed_data_across_shelf([0,1,2], allData, [mf.area_BeringSeaWest], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        data456 = get_processed_data_across_shelf([3,4,5], allData, [mf.area_BeringSeaWest], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        data789 = get_processed_data_across_shelf([6,7,8], allData, [mf.area_BeringSeaWest], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        data101112 = get_processed_data_across_shelf([9,10,11], allData, [mf.area_BeringSeaWest], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        # data123 = get_processed_data_across_shelf([0,1,2], allData, [mf.area_BeringSeaWest], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        # data789 = get_processed_data_across_shelf([6,7,8], allData, [mf.area_BeringSeaWest], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);

        distances = np.array(gridwiseData["distance"])[data123.regionMask]/1000.0; #in km
        dists, laruBeringSeaMask = mf.calculate_km_subsection_bounds_along_shelf(data123.xcoords, data123.ycoords, distances, mf.laruelle_BeringSea, params, testPlot=False);
        t2data.append(print_statistics(data_all, distances[laruBeringSeaMask], "Bering Sea (Annual)", laruBeringSeaMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        t2data.append(print_statistics(data123, distances[laruBeringSeaMask], "Bering Sea (Jan-Mar)", laruBeringSeaMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        t2data.append(print_statistics(data456, distances[laruBeringSeaMask], "Bering Sea (Apr-Jun)", laruBeringSeaMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        t2data.append(print_statistics(data789, distances[laruBeringSeaMask], "Bering Sea (Jul-Sep)", laruBeringSeaMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        t2data.append(print_statistics(data101112, distances[laruBeringSeaMask], "Bering Sea (Oct-Dec)", laruBeringSeaMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        # t2data.append(print_statistics(data123, distances[laruBeringSeaMask], "Bering Sea (Jan-Mar)", laruBeringSeaMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        # t2data.append(print_statistics(data789, distances[laruBeringSeaMask], "Bering Sea (Jul-Sep)", laruBeringSeaMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));


        #Antarctic Peninsula
        data_all =  get_processed_data_across_shelf([0,1,2,3,4,5,6,7,8,9,10,11], allData, [mf.area_AntarcticPeninsula], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        data123 = get_processed_data_across_shelf([0,1,2], allData, [mf.area_AntarcticPeninsula], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        data456 = get_processed_data_across_shelf([3,4,5], allData, [mf.area_AntarcticPeninsula], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        data789 = get_processed_data_across_shelf([6,7,8], allData, [mf.area_AntarcticPeninsula], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        data101112 = get_processed_data_across_shelf([9,10,11], allData, [mf.area_AntarcticPeninsula], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        # data123 = get_processed_data_across_shelf([0,1,2], allData, [mf.area_AntarcticPeninsula], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        # data789 = get_processed_data_across_shelf([6,7,8], allData, [mf.area_AntarcticPeninsula], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);

        distances = np.array(gridwiseData["distance"])[data123.regionMask]/1000.0; #in km
        dists, laruAntarcticPeninsulaMask = mf.calculate_km_subsection_bounds_along_shelf(data123.xcoords, data123.ycoords, distances, mf.laruelle_AntarcticPeninsula, params, testPlot=False);

        try:
            t2data.append(print_statistics(data_all, distances[laruAntarcticPeninsulaMask], "Antarctic Peninsula (Annual)", laruAntarcticPeninsulaMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
            t2data.append(print_statistics(data123, distances[laruAntarcticPeninsulaMask], "Antarctic Peninsula (Jan-Mar)", laruAntarcticPeninsulaMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
            t2data.append(print_statistics(data456, distances[laruAntarcticPeninsulaMask], "Antarctic Peninsula (Apr-Jun)", laruAntarcticPeninsulaMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
            t2data.append(print_statistics(data789, distances[laruAntarcticPeninsulaMask], "Antarctic Peninsula (Jul-Sep)", laruAntarcticPeninsulaMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
            t2data.append(print_statistics(data101112, distances[laruAntarcticPeninsulaMask], "Antarctic Peninsula (Oct-Dec)", laruAntarcticPeninsulaMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
            # t2data.append(print_statistics(data123, distances[laruAntarcticPeninsulaMask], "Antarctic Peninsula (Jan-Mar)", laruAntarcticPeninsulaMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
            # t2data.append(print_statistics(data789, distances[laruAntarcticPeninsulaMask], "Antarctic Peninsula (Jul-Sep)", laruAntarcticPeninsulaMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        except IndexError: #When distances and laruAntarcticPeninsulaMask are 0 because there were no shelf edges in the region, it fails.
                             #However, print_statistics handles this situation correctly.
            t2data.append(print_statistics(data_all, distances, "Antarctic Peninsula (Annual)", laruAntarcticPeninsulaMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
            t2data.append(print_statistics(data123, distances, "Antarctic Peninsula (Jan-Mar)", laruAntarcticPeninsulaMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
            t2data.append(print_statistics(data456, distances, "Antarctic Peninsula (Apr-Jun)", laruAntarcticPeninsulaMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
            t2data.append(print_statistics(data789, distances, "Antarctic Peninsula (Jul-Sep)", laruAntarcticPeninsulaMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
            t2data.append(print_statistics(data101112, distances, "Antarctic Peninsula (Oct-Dec)", laruAntarcticPeninsulaMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
            # t2data.append(print_statistics(data789, distances, "Antarctic Peninsula (Jul-Sep)", laruAntarcticPeninsulaMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
            # t2data.append(print_statistics(data789, distances, "Antarctic Peninsula (Jul-Sep)", laruAntarcticPeninsulaMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));


        #Labrador Sea
        data_all =  get_processed_data_across_shelf([0,1,2,3,4,5,6,7,8,9,10,11], allData, [mf.area_LabradorSea], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        data123 = get_processed_data_across_shelf([0,1,2], allData, [mf.area_LabradorSea], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        data456 = get_processed_data_across_shelf([3,4,5], allData, [mf.area_LabradorSea], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        data789 = get_processed_data_across_shelf([6,7,8], allData, [mf.area_LabradorSea], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        data101112 = get_processed_data_across_shelf([9,10,11], allData, [mf.area_LabradorSea], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        # data123 = get_processed_data_across_shelf([0,1,2], allData, [mf.area_LabradorSea], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        # data789 = get_processed_data_across_shelf([6,7,8], allData, [mf.area_LabradorSea], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        distances = np.array(gridwiseData["distance"])[data123.regionMask]/1000.0; #in km

        #Custom roll the data to eliminate gap
        rollAmount = -154
        data_all= roll_data(data_all,rollAmount)
        data123 = roll_data(data123,rollAmount)
        data456 = roll_data(data456,rollAmount)
        data789 = roll_data(data789,rollAmount)
        data101112 = roll_data(data101112,rollAmount)
        distances = np.roll(distances, rollAmount);

        dists, laruLabradorSeaMask = mf.calculate_km_subsection_bounds_along_shelf(data123.xcoords, data123.ycoords, distances, mf.laruelle_LabradorSea, params, testPlot=False);
        t2data.append(print_statistics(data_all, distances[laruLabradorSeaMask], "Labrador Sea (Annual)", laruLabradorSeaMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        t2data.append(print_statistics(data123, distances[laruLabradorSeaMask], "Labrador Sea (Jan-Mar)", laruLabradorSeaMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        t2data.append(print_statistics(data456, distances[laruLabradorSeaMask], "Labrador Sea (Apr-Jun)", laruLabradorSeaMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        t2data.append(print_statistics(data789, distances[laruLabradorSeaMask], "Labrador Sea (Jul-Sep)", laruLabradorSeaMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        t2data.append(print_statistics(data101112, distances[laruLabradorSeaMask], "Labrador Sea (Oct-Dec)", laruLabradorSeaMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        # t2data.append(print_statistics(data123, distances[laruLabradorSeaMask], "Labrador Sea (Jan-Mar)", laruLabradorSeaMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        # t2data.append(print_statistics(data789, distances[laruLabradorSeaMask], "Labrador Sea (Jul-Sep)", laruLabradorSeaMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));


        #Tasmanian Shelf
        data_all =  get_processed_data_across_shelf([0,1,2,3,4,5,6,7,8,9,10,11], allData, [mf.area_Tasmania], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        data123 = get_processed_data_across_shelf([0,1,2], allData, [mf.area_Tasmania], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        data456 = get_processed_data_across_shelf([3,4,5], allData, [mf.area_Tasmania], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        data789 = get_processed_data_across_shelf([6,7,8], allData, [mf.area_Tasmania], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        data101112 = get_processed_data_across_shelf([9,10,11], allData, [mf.area_Tasmania], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        # data123 = get_processed_data_across_shelf([0,1,2], allData, [mf.area_Tasmania], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        # data789 = get_processed_data_across_shelf([6,7,8], allData, [mf.area_Tasmania], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);

        distances = np.array(gridwiseData["distance"])[data123.regionMask]/1000.0; #in km
        dists, laruTasmanianShelfMask = mf.calculate_km_subsection_bounds_along_shelf(data123.xcoords, data123.ycoords, distances, mf.laruelle_TasmanianShelf, params, testPlot=False);
        t2data.append(print_statistics(data_all, distances[laruTasmanianShelfMask], "Tasmanian Shelf (Annual)", laruTasmanianShelfMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        t2data.append(print_statistics(data123, distances[laruTasmanianShelfMask], "Tasmanian Shelf (Jan-Mar)", laruTasmanianShelfMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        t2data.append(print_statistics(data456, distances[laruTasmanianShelfMask], "Tasmanian Shelf (Apr-Jun)", laruTasmanianShelfMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        t2data.append(print_statistics(data789, distances[laruTasmanianShelfMask], "Tasmanian Shelf (Jul-Sep)", laruTasmanianShelfMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        t2data.append(print_statistics(data101112, distances[laruTasmanianShelfMask], "Tasmanian Shelf (Oct-Dec)", laruTasmanianShelfMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        # t2data.append(print_statistics(data123, distances[laruTasmanianShelfMask], "Tasmanian Shelf (Jan-Mar)", laruTasmanianShelfMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        # t2data.append(print_statistics(data789, distances[laruTasmanianShelfMask], "Tasmanian Shelf (Jul-Sep)", laruTasmanianShelfMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));


        #Barents Sea
        data_all =  get_processed_data_across_shelf([0,1,2,3,4,5,6,7,8,9,10,11], allData, [mf.area_BarentsSea], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        data123 = get_processed_data_across_shelf([0,1,2], allData, [mf.area_BarentsSea], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        data456 = get_processed_data_across_shelf([3,4,5], allData, [mf.area_BarentsSea], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        data789 = get_processed_data_across_shelf([6,7,8], allData, [mf.area_BarentsSea], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        data101112 = get_processed_data_across_shelf([9,10,11], allData, [mf.area_BarentsSea], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        # data123 = get_processed_data_across_shelf([0,1,2], allData, [mf.area_BarentsSea], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        # data789 = get_processed_data_across_shelf([6,7,8], allData, [mf.area_BarentsSea], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);

        distances = np.array(gridwiseData["distance"])[data123.regionMask]/1000.0; #in km
        dists, laruBarentsSeaMask = mf.calculate_km_subsection_bounds_along_shelf(data123.xcoords, data123.ycoords, distances, mf.laruelle_BarentsSea, params, testPlot=False);
        t2data.append(print_statistics(data_all, distances[laruBarentsSeaMask], "Barents Sea (Annual)", laruBarentsSeaMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        t2data.append(print_statistics(data123, distances[laruBarentsSeaMask], "Barents Sea (Jan-Mar)", laruBarentsSeaMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        t2data.append(print_statistics(data456, distances[laruBarentsSeaMask], "Barents Sea (Apr-Jun)", laruBarentsSeaMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        t2data.append(print_statistics(data789, distances[laruBarentsSeaMask], "Barents Sea (Jul-Sep)", laruBarentsSeaMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        t2data.append(print_statistics(data101112, distances[laruBarentsSeaMask], "Barents Sea (Oct-Dec)", laruBarentsSeaMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        # t2data.append(print_statistics(data123, distances[laruBarentsSeaMask], "Barents Sea (Jan-Mar)", laruBarentsSeaMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        # t2data.append(print_statistics(data789, distances[laruBarentsSeaMask], "Barents Sea (Jul-Sep)", laruBarentsSeaMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));


        #South Atlantic Bight
        data_all =  get_processed_data_across_shelf([0,1,2,3,4,5,6,7,8,9,10,11], allData, [mf.area_SouthAtlanticBight], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        data123 = get_processed_data_across_shelf([0,1,2], allData, [mf.area_SouthAtlanticBight], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        data456 = get_processed_data_across_shelf([3,4,5], allData, [mf.area_SouthAtlanticBight], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        data789 = get_processed_data_across_shelf([6,7,8], allData, [mf.area_SouthAtlanticBight], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        data101112 = get_processed_data_across_shelf([9,10,11], allData, [mf.area_SouthAtlanticBight], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        # data123 = get_processed_data_across_shelf([0,1,2], allData, [mf.area_SouthAtlanticBight], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        # data789 = get_processed_data_across_shelf([6,7,8], allData, [mf.area_SouthAtlanticBight], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);

        distances = np.array(gridwiseData["distance"])[data123.regionMask]/1000.0; #in km
        dists, laruSouthAtlanticBightMask = mf.calculate_km_subsection_bounds_along_shelf(data123.xcoords, data123.ycoords, distances, mf.laruelle_SouthAtlanticBight, params, testPlot=False);
        t2data.append(print_statistics(data_all, distances[laruSouthAtlanticBightMask], "South Atlantic Bight (Annual)", laruSouthAtlanticBightMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        t2data.append(print_statistics(data123, distances[laruSouthAtlanticBightMask], "South Atlantic Bight (Jan-Mar)", laruSouthAtlanticBightMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        t2data.append(print_statistics(data456, distances[laruSouthAtlanticBightMask], "South Atlantic Bight (Apr-Jun)", laruSouthAtlanticBightMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        t2data.append(print_statistics(data789, distances[laruSouthAtlanticBightMask], "South Atlantic Bight (Jul-Sep)", laruSouthAtlanticBightMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        t2data.append(print_statistics(data101112, distances[laruSouthAtlanticBightMask], "South Atlantic Bight (Oct-Dec)", laruSouthAtlanticBightMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        # t2data.append(print_statistics(data123, distances[laruSouthAtlanticBightMask], "South Atlantic Bight (Jan-Mar)", laruSouthAtlanticBightMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        # t2data.append(print_statistics(data789, distances[laruSouthAtlanticBightMask], "South Atlantic Bight (Jul-Sep)", laruSouthAtlanticBightMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));


        #Southern Greenland
        data_all =  get_processed_data_across_shelf([0,1,2,3,4,5,6,7,8,9,10,11], allData, [mf.area_SouthGreenland], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        data123 = get_processed_data_across_shelf([0,1,2], allData, [mf.area_SouthGreenland], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        data456 = get_processed_data_across_shelf([3,4,5], allData, [mf.area_SouthGreenland], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        data789 = get_processed_data_across_shelf([6,7,8], allData, [mf.area_SouthGreenland], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        data101112 = get_processed_data_across_shelf([9,10,11], allData, [mf.area_SouthGreenland], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        # data123 = get_processed_data_across_shelf([0,1,2], allData, [mf.area_SouthGreenland], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        # data789 = get_processed_data_across_shelf([6,7,8], allData, [mf.area_SouthGreenland], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);

        distances = np.array(gridwiseData["distance"])[data123.regionMask]/1000.0; #in km
        dists, laruSouthernGreenlandMask = mf.calculate_km_subsection_bounds_along_shelf(data123.xcoords, data123.ycoords, distances, mf.laruelle_SouthernGreenland, params, testPlot=False);
        t2data.append(print_statistics(data_all, distances[laruSouthernGreenlandMask], "Southern Greenland (Annual)", laruSouthernGreenlandMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        t2data.append(print_statistics(data123, distances[laruSouthernGreenlandMask], "Southern Greenland (Jan-Mar)", laruSouthernGreenlandMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        t2data.append(print_statistics(data456, distances[laruSouthernGreenlandMask], "Southern Greenland  (Apr-Jun)", laruSouthernGreenlandMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        t2data.append(print_statistics(data789, distances[laruSouthernGreenlandMask], "Southern Greenland  (Jul-Sep)", laruSouthernGreenlandMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        t2data.append(print_statistics(data101112, distances[laruSouthernGreenlandMask], "Southern Greenland  (Oct-Dec)", laruSouthernGreenlandMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        # t2data.append(print_statistics(data123, distances[laruSouthernGreenlandMask], "Southern Greenland (Jan-Mar)", laruSouthernGreenlandMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        # t2data.append(print_statistics(data789, distances[laruSouthernGreenlandMask], "Southern Greenland (Jul-Sep)", laruSouthernGreenlandMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));


        #Cascadian Shelf
        data_all =  get_processed_data_across_shelf([0,1,2,3,4,5,6,7,8,9,10,11], allData, [mf.area_CascianShelf], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        data123 = get_processed_data_across_shelf([0,1,2], allData, [mf.area_CascianShelf], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        data456 = get_processed_data_across_shelf([3,4,5], allData, [mf.area_CascianShelf], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        data789 = get_processed_data_across_shelf([6,7,8], allData, [mf.area_CascianShelf], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        data101112 = get_processed_data_across_shelf([9,10,11], allData, [mf.area_CascianShelf], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        # data123 = get_processed_data_across_shelf([0,1,2], allData, [mf.area_CascianShelf], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        # data789 = get_processed_data_across_shelf([6,7,8], allData, [mf.area_CascianShelf], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);

        distances = np.array(gridwiseData["distance"])[data123.regionMask]/1000.0; #in km
        dists, laruCascadianShelfMask = mf.calculate_km_subsection_bounds_along_shelf(data123.xcoords, data123.ycoords, distances, mf.laruelle_CascadianShelf, params, testPlot=False);
        t2data.append(print_statistics(data_all, distances[laruCascadianShelfMask], "Cascadian Shelf (Annual)", laruCascadianShelfMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        t2data.append(print_statistics(data123, distances[laruCascadianShelfMask], "Cascadian Shelf (Jan-Mar)", laruCascadianShelfMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        t2data.append(print_statistics(data456, distances[laruCascadianShelfMask], "Cascadian Shelf  (Apr-Jun)", laruCascadianShelfMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        t2data.append(print_statistics(data789, distances[laruCascadianShelfMask], "Cascadian Shelf  (Jul-Sep)", laruCascadianShelfMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        t2data.append(print_statistics(data101112, distances[laruCascadianShelfMask], "Cascadian Shelf  (Oct-Dec)", laruCascadianShelfMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        # t2data.append(print_statistics(data123, distances[laruCascadianShelfMask], "Cascadian Shelf (Jan-Mar)", laruCascadianShelfMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        # t2data.append(print_statistics(data789, distances[laruCascadianShelfMask], "Cascadian Shelf (Jul-Sep)", laruCascadianShelfMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));


        #Irminger Sea
        data_all =  get_processed_data_across_shelf([0,1,2,3,4,5,6,7,8,9,10,11], allData, [mf.area_IrmingerSea], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        data123 = get_processed_data_across_shelf([0,1,2], allData, [mf.area_IrmingerSea], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        data456 = get_processed_data_across_shelf([3,4,5], allData, [mf.area_IrmingerSea], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        data789 = get_processed_data_across_shelf([6,7,8], allData, [mf.area_IrmingerSea], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        data101112 = get_processed_data_across_shelf([9,10,11], allData, [mf.area_IrmingerSea], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        # data123 = get_processed_data_across_shelf([0,1,2], allData, [mf.area_IrmingerSea], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);
        # data789 = get_processed_data_across_shelf([6,7,8], allData, [mf.area_IrmingerSea], params, gasTransferParameterisation=calculateTable2GasTransferVelocity);

        distances = np.array(gridwiseData["distance"])[data123.regionMask]/1000.0; #in km
        dists, laruIrmingerSeaMask = mf.calculate_km_subsection_bounds_along_shelf(data123.xcoords, data123.ycoords, distances, mf.laruelle_IrmingerSea, params, testPlot=False);
        t2data.append(print_statistics(data_all, distances[laruIrmingerSeaMask], "Irminger Sea (Annual)", laruIrmingerSeaMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        t2data.append(print_statistics(data123, distances[laruIrmingerSeaMask], "Irminger Sea (Jan-Mar)", laruIrmingerSeaMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        t2data.append(print_statistics(data456, distances[laruIrmingerSeaMask], "Irminger Sea  (Apr-Jun)", laruIrmingerSeaMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        t2data.append(print_statistics(data789, distances[laruIrmingerSeaMask], "Irminger Sea  (Jul-Sep)", laruIrmingerSeaMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        t2data.append(print_statistics(data101112, distances[laruIrmingerSeaMask], "Irminger Sea (Oct-Dec)", laruIrmingerSeaMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        # t2data.append(print_statistics(data123, distances[laruIrmingerSeaMask], "Irminger Sea (Jan-Mar)", laruIrmingerSeaMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));
        # t2data.append(print_statistics(data789, distances[laruIrmingerSeaMask], "Irminger Sea (Jul-Sep)", laruIrmingerSeaMask, gasTransferParameterisation=calculateTable2GasTransferVelocity));

        df = pd.DataFrame(t2data);
        df.columns = ["region", "total across-shelf", "total across-shelf SD","total across-shelf_n", "EksAbsolute", "EksAbsoluteStd","EksAbsolute_n", "GeoAbsolute", "GeoAbsoluteStd","GeoAbsolute_n","StokesAbsolute", "StokesAbsoluteStd","StokesAbsolute_n",
            "proportionEks", "proportionEksSD","proportionEks_n", "proportionGeo", "proportionGeoSD","proportionGeo_n", "proportionStokes", "proportionStokesSD", "proportionStokes_n", "percentEks", "percentEksSD","percetnEks_n", "percentGeo", "percentGeoSD",
            "percentGeo_n","percentStokes", "percentStokesSD","percentStokes_n", "k", "kSD","k_n"];
        df.to_csv(path.join(outputPath, "output_means", "table2data_"+params.paramsetName+'_'+str(params.start_year) +"_"+str(params.end_year)+ ".csv"), index=False);
