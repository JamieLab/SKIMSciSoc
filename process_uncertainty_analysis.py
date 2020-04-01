#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 13:43:32 2018

v1 uses celldata.nTotalMagnitude  (sum of abs(acrossShelfCurrent));
v2 uses each component and calculates a directional total (and uses stokes filtering correctly)
   also calculates the differences from baseline
v3 calculates standard deviation of the differences from baseline
v4 adds error to the total current, then calculates across shelf currents.
v5 switch for turning months 789 off.

@author: rr
"""

import cPickle as pickle;
import numpy as np;
import matplotlib.pyplot as plt;
import python_util.skim_utilities as su;
import python_util.mask_functions as mf;
from os import path;
import pandas as pd;
from statsmodels.stats.weightstats import DescrStatsW;
import python_util.parameter_sets as ps;


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
#   total current magnitude, list for each [cell, month]
#months defines the months of the year to use (e.g. 0=Jan, 11=Dec)
def get_processed_data2(months, data, regionMaskBoundsList, params, filterStokes=True):
    if len(data) % 12 != 0:
        raise ValueError ("get_processed_data only supports data for whole years (multiples of 12 months)");
    
    #Create the indices to use
    numYears = int(len(data)/12);
    indices = get_month_indices(months, numYears);
    
    #Create arrays to store data for ekman, geostrophic and stokes data
    length = len(data[0]);
    
    stokesMaskPass = np.empty((length, len(indices))); stokesMaskPass[:] = np.nan;
    ontoShelfVectorX = np.empty((length, len(indices))); ontoShelfVectorX[:] = np.nan;
    ontoShelfVectorY = np.empty((length, len(indices))); ontoShelfVectorY[:] = np.nan;
    totalCurrentX = np.empty((length, len(indices))); totalCurrentX[:] = np.nan;
    totalCurrentY = np.empty((length, len(indices))); totalCurrentY[:] = np.nan;
    xcoords = np.empty((length, len(indices))); xcoords[:] = np.nan; #xcoords and ycoords are just used for checking region mask is correct
    ycoords = np.empty((length, len(indices))); ycoords[:] = np.nan;
    
    #Copy / calculate data:
    for i, index in enumerate(indices):
        stokesMaskPass[:,i] = [curData.stokesMaskPass if np.isfinite(curData.stokesMaskPass) else np.nan for curData in data[index]];
        ontoShelfVectorX[:,i] = [curData.ontoShelfVector[0] if np.all(np.isfinite(curData.ontoShelfVector)) else np.nan for curData in data[index]];
        ontoShelfVectorY[:,i] = [curData.ontoShelfVector[1] if np.all(np.isfinite(curData.ontoShelfVector)) else np.nan for curData in data[index]];
        xcoords[:,i] = [curData.indexX if np.isfinite(curData.indexX) else np.nan for curData in data[index]];
        ycoords[:,i] = [curData.indexY if np.isfinite(curData.indexY) else np.nan for curData in data[index]];
        
        #Calculate the north and east components of the total current (taking the stokes drift mask into account)
        for j, curData in enumerate(data[index]):
            totalCurrent = curData.ekmanCurrentVector + curData.geostrophicCurrentVector;
            if stokesMaskPass[j,i] == True:
                totalCurrent += curData.stokesCurrentVector;
            if np.all(np.isfinite(totalCurrent)):
                totalCurrentX[j,i] = totalCurrent[0];
                totalCurrentY[j,i] = totalCurrent[1];
            else:
                totalCurrentX[j,i] = np.nan;
                totalCurrentY[j,i] = np.nan;

    #Calculate means and standard deviations
    results = su.QuickStruct(); #Store results in a struct
    #do not calculate means and standard deviations yet! This should be done at the noise adding step!
    
    #Construct outputs using region mask to only select data inside the region
    regionMask = mf.return_area_mask(regionMaskBoundsList, data[0], params);
    results.xcoords = xcoords[regionMask,0]; #Does not change in time, so discard time dimension
    results.ycoords = ycoords[regionMask,0]; #Does not change in time, so discard time dimension
    results.ontoShelfVectorX = ontoShelfVectorX[regionMask,0]; #Does not change in time, so discard time dimension
    results.ontoShelfVectorY = ontoShelfVectorY[regionMask,0]; #Does not change in time, so discard time dimension
    results.stokesMaskPass = stokesMaskPass[regionMask,:]; #Keep time dimension
    results.totalCurrentX = totalCurrentX[regionMask,:]; #Keep time dimension
    results.totalCurrentY = totalCurrentY[regionMask,:]; #Keep time dimension
    results.regionMask = regionMask;
    
    return results;


def copy_with_error(totalCurrentX, totalCurrentY, rmseGeoEk, rmseStokes, stokesUsed):
    totalRMSE = np.sqrt(rmseGeoEk**2 + rmseStokes**2);
    
    newTotalCurrentX = totalCurrentX.copy();
    newTotalCurrentY = totalCurrentY.copy();
    if totalRMSE != 0:
        noStokesNoiseX = np.random.normal(0.0, rmseGeoEk, totalCurrentX.shape);
        totalNoiseX = np.random.normal(0.0, totalRMSE, totalCurrentX.shape);
        noStokesNoiseY = np.random.normal(0.0, rmseGeoEk, totalCurrentY.shape);
        totalNoiseY = np.random.normal(0.0, totalRMSE, totalCurrentY.shape);
        
        #Add appropriate noise choice to the correct elements.
        whereNoStokes = np.where(stokesUsed == 0.0);
        whereTotal = np.where(stokesUsed != 0.0);
        newTotalCurrentX[whereNoStokes] += noStokesNoiseX[whereNoStokes];
        newTotalCurrentY[whereNoStokes] += noStokesNoiseY[whereNoStokes];
        newTotalCurrentX[whereTotal] += totalNoiseX[whereTotal];
        newTotalCurrentY[whereTotal] += totalNoiseY[whereTotal];
    return newTotalCurrentX, newTotalCurrentY;


def calculate_across_shelf_current(totalCurrentX, totalCurrentY, ontoShelfVectorX, ontoShelfVectorY):
    totalCurrents = np.array([totalCurrentX, totalCurrentY]);
    ontoShelfVector = np.array([ontoShelfVectorX, ontoShelfVectorY]);
    acrossShelfPointCurrent = np.dot(ontoShelfVector, totalCurrents) / np.linalg.norm(ontoShelfVector);
    return acrossShelfPointCurrent;



#Given _totalVals (cell-wise across-shelf currents), calculate the mean and SD for the region.
def calc_statistics(_totalVals, _distances, name, mask):
    if len(_totalVals) == 0 or len(mask) == 0:
        print name, "- No data available!\n";
        return (name, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan);

    totalVals = _totalVals[mask,:];
    distances = _distances[mask];

    #distances = np.array(distances);
    if (len(distances) != len(totalVals)):
        raise ValueError("distances must be the same length as data values");
    
    weights = np.transpose(np.array([distances]*totalVals.shape[1]));

    notnan = np.where(np.isnan(totalVals)==False);
    totalStats = DescrStatsW(totalVals[notnan], weights=weights[notnan], ddof=0);
    totalGrandMean = totalStats.mean;
    totalGrandSD = totalStats.std;

    #print name, "total onto-shelf current (m/s):", totalGrandMean, "+/-", totalGrandSD;
    
    #Return the mean ontoshelf current for this region (weighted by segment distance)
    #Return the SD of the same.
    return (name, totalGrandMean, totalGrandSD);



def calc_region_err_stats(data, regionMask, regionLabel, errorsGeoEk, errorsStokes, numReps, distances):
    regionMeans = []; #means for each error level
    regionSDs = [];
    repMeans = []; #2D, mean across-shelf current for each repeat and error level [err, rep]
    for iError in range(len(errorsGeoEk)):
        errMeans = []; #Collect all the mean across-shelf currents for each repeat at a single RMSE.
        errSDs = [];
        print regionLabel, iError, "of", len(errorsGeoEk);
        for repeat in range(0, numReps):
            #copy data and add errors
            totalCurrentsX, totalCurrentsY = copy_with_error(data.totalCurrentX, data.totalCurrentX, errorsGeoEk[iError], errorsStokes[iError], data.stokesMaskPass);
            
            #calculate across-shelf current from the data with errors
            acrossShelfCurrents = np.empty((totalCurrentsX.shape), dtype=float); acrossShelfCurrents[:] = np.nan;
            for index, _ in np.ndenumerate(totalCurrentsX):
                acrossShelfCurrents[index] = calculate_across_shelf_current(totalCurrentsX[index], totalCurrentsY[index], data.ontoShelfVectorX[index[0]], data.ontoShelfVectorY[index[0]]);
            
            #Calculate statistics for the region around the total across-shelf current
            #Returns mean and SD of accross shelf current for the region. Single repeat, single error level.
            name, curMean, curSD = calc_statistics(acrossShelfCurrents, distances, regionLabel, regionMask);
            errMeans.append(curMean);
            errSDs.append(curSD);
        repMeans.append(np.array(errMeans)); #list of arrays, inner array is all the means for each repeat, outer list is for each RMSE
        regionMeans.append(np.mean(errMeans)); #
        regionSDs.append(np.std(errMeans)); #Is it this or is it just the SD of the means that we want?
    
    #calculate percentage and absolute difference
    repMeans = np.array(repMeans); #dims: [err, rep]
    baseline = np.mean(repMeans[0,:]);
    absDiffs = np.abs(repMeans-baseline);
    regionDiffFromMean = np.mean(absDiffs, axis=1);
    regionDiffFromMeanSD = np.std(absDiffs, axis=1);
    regionDiffFromMeanPercentage = regionDiffFromMean/baseline * 100.0;
    
    #Return:
    #   regionMean: the mean across-shelf current across all repeats, for each RMSE.
    #   regionSDs: the SD across all repeats, for each RMSE.
    #   regionDiffFromMean: absolute difference from the baseline (e.g. different from no uncertainty/RMSE=0). Average across all repeats. I.e. how much does this error level deviate from the baseline?
    #   regionDiffFromMeanSD: SD of the absolute difference from the baseline. I.e. how much does the deviation from baseline deviate?
    #   regionDiffFromMeanPercentage: Same as regionDiffFromMean but as a percentage of the baseline.
    return regionMeans, regionSDs, regionDiffFromMean, regionDiffFromMeanSD, regionDiffFromMeanPercentage;


numRepeats = 1001;
outputPath = path.join("results", "sensitivity_to_uncertainty");
params = ps.get_global_params();
process123 = False;

proportionActualRMSEs = np.arange(0.0, 2.1, 0.1); #percentage of actual RMSE
geoEksActualRMSE = 0.135;
stokesActualRMSE = 0.03;
geoEksRMSEs = geoEksActualRMSE*proportionActualRMSEs;
stokesRMSEs = stokesActualRMSE*proportionActualRMSEs;


#Read data
inputDataPath = path.join("data", params.paramsetName);
gridwiseData = pd.read_table(path.join(inputDataPath, "gridwise_data/per_grid_cell_edge_data_"+params.paramsetName+".csv"), sep=',');
gridwiseData.x = gridwiseData.x.astype(int);
gridwiseData.y = gridwiseData.y.astype(int);
#allData = pickle.load(open(path.join(inputDataPath, "current_data/surface_currents_"+params.paramsetName+".p"), "rb"));
#allData = pickle.load(open(path.join("/home/rr/Files/Tasks/20180914_SKIM/data/global_500m/current_data/surface_currents_uncertainty_aid_global_500m.p"), "rb"));

mainDF = pd.DataFrame();
totalRMSEs = np.sqrt(np.array(geoEksRMSEs)**2 + np.array(stokesRMSEs)**2);
mainDF["RMSE"] = totalRMSEs;



############################
# North Sea
data123 = get_processed_data2([0,1,2], allData, [mf.area_EuropeanShelf], params);
data789 = get_processed_data2([6,7,8], allData, [mf.area_EuropeanShelf], params);
distances = np.array(gridwiseData["distance"])[data123.regionMask]/1000.0; #in km
dists, laruNorthSeaMask = mf.calculate_km_subsection_bounds_along_shelf(data123.xcoords, data123.ycoords, distances, mf.laruelle_NorthSea, params, testPlot=False);

regionLabel="NorthSea123";
regionMeans, regionSDs, regionDiffFromMean, regionDiffFromMeanSD, regionDiffFromMeanPercentage = calc_region_err_stats(data123, laruNorthSeaMask, regionLabel, errorsGeoEk=geoEksRMSEs, errorsStokes=stokesRMSEs, numReps=numRepeats, distances=distances);
mainDF[regionLabel+"_mean"] = regionMeans;
mainDF[regionLabel+"SD"] = regionSDs;
mainDF[regionLabel+"_diff"] = regionDiffFromMean;
mainDF[regionLabel+"_diffSD"] = regionDiffFromMeanSD;
mainDF[regionLabel+"_diffpercent"] = regionDiffFromMeanPercentage;

if process123:
    regionLabel="NorthSea789";
    regionMeans, regionSDs, regionDiffFromMean, regionDiffFromMeanSD, regionDiffFromMeanPercentage = calc_region_err_stats(data789, laruNorthSeaMask, regionLabel, errorsGeoEk=geoEksRMSEs, errorsStokes=stokesRMSEs, numReps=numRepeats, distances=distances);
    mainDF[regionLabel+"_mean"] = regionMeans;
    mainDF[regionLabel+"SD"] = regionSDs;
    mainDF[regionLabel+"_diff"] = regionDiffFromMean;
    mainDF[regionLabel+"_diffSD"] = regionDiffFromMeanSD;
    mainDF[regionLabel+"_diffpercent"] = regionDiffFromMeanPercentage;



#########################
# English Channel
dists, laruEnglishChannelMask = mf.calculate_km_subsection_bounds_along_shelf(data123.xcoords, data123.ycoords, distances, mf.laruelle_EnglishChannel, params, testPlot=False);

regionLabel="EnglishChannel123";
regionMeans, regionSDs, regionDiffFromMean, regionDiffFromMeanSD, regionDiffFromMeanPercentage = calc_region_err_stats(data123, laruEnglishChannelMask, regionLabel, errorsGeoEk=geoEksRMSEs, errorsStokes=stokesRMSEs, numReps=numRepeats, distances=distances);
mainDF[regionLabel+"_mean"] = regionMeans;
mainDF[regionLabel+"SD"] = regionSDs;
mainDF[regionLabel+"_diff"] = regionDiffFromMean;
mainDF[regionLabel+"_diffSD"] = regionDiffFromMeanSD;
mainDF[regionLabel+"_diffpercent"] = regionDiffFromMeanPercentage;

if process123:
    regionLabel="EnglishChannel789";
    regionMeans, regionSDs, regionDiffFromMean, regionDiffFromMeanSD, regionDiffFromMeanPercentage = calc_region_err_stats(data789, laruEnglishChannelMask, regionLabel, errorsGeoEk=geoEksRMSEs, errorsStokes=stokesRMSEs, numReps=numRepeats, distances=distances);
    mainDF[regionLabel+"_mean"] = regionMeans;
    mainDF[regionLabel+"SD"] = regionSDs;
    mainDF[regionLabel+"_diff"] = regionDiffFromMean;
    mainDF[regionLabel+"_diffSD"] = regionDiffFromMeanSD;
    mainDF[regionLabel+"_diffpercent"] = regionDiffFromMeanPercentage;



######################
#MidAtlanticBight
data123 = get_processed_data2([0,1,2], allData, [mf.area_MidAtlanticBight], params);
data789 = get_processed_data2([6,7,8], allData, [mf.area_MidAtlanticBight], params);
distances = np.array(gridwiseData["distance"])[data123.regionMask]/1000.0; #in km
dists, laruMidAtlanticBightMask = mf.calculate_km_subsection_bounds_along_shelf(data123.xcoords, data123.ycoords, distances, mf.laruelle_MidAtlanticBight, params, testPlot=False);

regionLabel="MidAtlanticBight123";
regionMeans, regionSDs, regionDiffFromMean, regionDiffFromMeanSD, regionDiffFromMeanPercentage = calc_region_err_stats(data123, laruMidAtlanticBightMask, regionLabel, errorsGeoEk=geoEksRMSEs, errorsStokes=stokesRMSEs, numReps=numRepeats, distances=distances);
mainDF[regionLabel+"_mean"] = regionMeans;
mainDF[regionLabel+"SD"] = regionSDs;
mainDF[regionLabel+"_diff"] = regionDiffFromMean;
mainDF[regionLabel+"_diffSD"] = regionDiffFromMeanSD;
mainDF[regionLabel+"_diffpercent"] = regionDiffFromMeanPercentage;

if process123:
    regionLabel="MidAtlanticBight789";
    regionMeans, regionSDs, regionDiffFromMean, regionDiffFromMeanSD, regionDiffFromMeanPercentage = calc_region_err_stats(data789, laruMidAtlanticBightMask, regionLabel, errorsGeoEk=geoEksRMSEs, errorsStokes=stokesRMSEs, numReps=numRepeats, distances=distances);
    mainDF[regionLabel+"_mean"] = regionMeans;
    mainDF[regionLabel+"SD"] = regionSDs;
    mainDF[regionLabel+"_diff"] = regionDiffFromMean;
    mainDF[regionLabel+"_diffSD"] = regionDiffFromMeanSD;
    mainDF[regionLabel+"_diffpercent"] = regionDiffFromMeanPercentage;


##########################
#CoastOfJapan
data123 = get_processed_data2([0,1,2], allData, [mf.area_CoastOfJapan], params);
data789 = get_processed_data2([6,7,8], allData, [mf.area_CoastOfJapan], params);
distances = np.array(gridwiseData["distance"])[data123.regionMask]/1000.0; #in km

#Custom roll for Japanese coast (to eliminate split)
data123.xcoords = np.roll(data123.xcoords, -110);
data123.ycoords = np.roll(data123.ycoords, -110);
data123.regionMask = np.roll(data123.regionMask, -110);
data123.ontoShelfVectorX = np.roll(data123.ontoShelfVectorX, -110);
data123.ontoShelfVectorY = np.roll(data123.ontoShelfVectorY, -110);
data123.stokesMaskPass = np.roll(data123.stokesMaskPass, -110);
data123.totalCurrentX = np.roll(data123.totalCurrentX, -110);
data123.totalCurrentY = np.roll(data123.totalCurrentY, -110);

if process123:
    data789.xcoords = np.roll(data789.xcoords, -110);
    data789.ycoords = np.roll(data789.ycoords, -110);
    data789.regionMask = np.roll(data789.regionMask, -110);
    data789.ontoShelfVectorX = np.roll(data789.ontoShelfVectorX, -110);
    data789.ontoShelfVectorY = np.roll(data789.ontoShelfVectorY, -110);
    data789.stokesMaskPass = np.roll(data789.stokesMaskPass, -110);
    data789.totalCurrentX = np.roll(data789.totalCurrentX, -110);
    data789.totalCurrentY = np.roll(data789.totalCurrentY, -110);

distances = np.roll(distances, -110);

dists, laruCoastOfJapanMask = mf.calculate_km_subsection_bounds_along_shelf(data123.xcoords, data123.ycoords, distances, mf.laruelle_CoastOfJapan, params, testPlot=False);

regionLabel="CoastOfJapan123";
regionMeans, regionSDs, regionDiffFromMean, regionDiffFromMeanSD, regionDiffFromMeanPercentage = calc_region_err_stats(data123, laruCoastOfJapanMask, regionLabel, errorsGeoEk=geoEksRMSEs, errorsStokes=stokesRMSEs, numReps=numRepeats, distances=distances);
mainDF[regionLabel+"_mean"] = regionMeans;
mainDF[regionLabel+"SD"] = regionSDs;
mainDF[regionLabel+"_diff"] = regionDiffFromMean;
mainDF[regionLabel+"_diffSD"] = regionDiffFromMeanSD;
mainDF[regionLabel+"_diffpercent"] = regionDiffFromMeanPercentage;

if process123:
    regionLabel="CoastOfJapan789";
    regionMeans, regionSDs, regionDiffFromMean, regionDiffFromMeanSD, regionDiffFromMeanPercentage = calc_region_err_stats(data789, laruCoastOfJapanMask, regionLabel, errorsGeoEk=geoEksRMSEs, errorsStokes=stokesRMSEs, numReps=numRepeats, distances=distances);
    mainDF[regionLabel+"_mean"] = regionMeans;
    mainDF[regionLabel+"SD"] = regionSDs;
    mainDF[regionLabel+"_diff"] = regionDiffFromMean;
    mainDF[regionLabel+"_diffSD"] = regionDiffFromMeanSD;
    mainDF[regionLabel+"_diffpercent"] = regionDiffFromMeanPercentage;


################################
#Patagonian shelf
data123 = get_processed_data2([0,1,2], allData, [mf.area_Patagonia], params);
data789 = get_processed_data2([6,7,8], allData, [mf.area_Patagonia], params);

distances = np.array(gridwiseData["distance"])[data123.regionMask]/1000.0; #in km
dists, laruPatagonianShelfMask = mf.calculate_km_subsection_bounds_along_shelf(data123.xcoords, data123.ycoords, distances, mf.laruelle_PatagonianShelf, params, testPlot=False);

regionLabel="PatagonianShelf123";
regionMeans, regionSDs, regionDiffFromMean, regionDiffFromMeanSD, regionDiffFromMeanPercentage = calc_region_err_stats(data123, laruPatagonianShelfMask, regionLabel, errorsGeoEk=geoEksRMSEs, errorsStokes=stokesRMSEs, numReps=numRepeats, distances=distances);
mainDF[regionLabel+"_mean"] = regionMeans;
mainDF[regionLabel+"SD"] = regionSDs;
mainDF[regionLabel+"_diff"] = regionDiffFromMean;
mainDF[regionLabel+"_diffSD"] = regionDiffFromMeanSD;
mainDF[regionLabel+"_diffpercent"] = regionDiffFromMeanPercentage;

if process123:
    regionLabel="PatagonianShelf789";
    regionMeans, regionSDs, regionDiffFromMean, regionDiffFromMeanSD, regionDiffFromMeanPercentage = calc_region_err_stats(data789, laruPatagonianShelfMask, regionLabel, errorsGeoEk=geoEksRMSEs, errorsStokes=stokesRMSEs, numReps=numRepeats, distances=distances);
    mainDF[regionLabel+"_mean"] = regionMeans;
    mainDF[regionLabel+"SD"] = regionSDs;
    mainDF[regionLabel+"_diff"] = regionDiffFromMean;
    mainDF[regionLabel+"_diffSD"] = regionDiffFromMeanSD;
    mainDF[regionLabel+"_diffpercent"] = regionDiffFromMeanPercentage;


##########################
#Bering Sea
data123 = get_processed_data2([0,1,2], allData, [mf.area_BeringSeaWest], params);
data789 = get_processed_data2([6,7,8], allData, [mf.area_BeringSeaWest], params);

distances = np.array(gridwiseData["distance"])[data123.regionMask]/1000.0; #in km
dists, laruBeringSeaMask = mf.calculate_km_subsection_bounds_along_shelf(data123.xcoords, data123.ycoords, distances, mf.laruelle_BeringSea, params, testPlot=False);

regionLabel="BeringSea123";
regionMeans, regionSDs, regionDiffFromMean, regionDiffFromMeanSD, regionDiffFromMeanPercentage = calc_region_err_stats(data123, laruBeringSeaMask, regionLabel, errorsGeoEk=geoEksRMSEs, errorsStokes=stokesRMSEs, numReps=numRepeats, distances=distances);
mainDF[regionLabel+"_mean"] = regionMeans;
mainDF[regionLabel+"SD"] = regionSDs;
mainDF[regionLabel+"_diff"] = regionDiffFromMean;
mainDF[regionLabel+"_diffSD"] = regionDiffFromMeanSD;
mainDF[regionLabel+"_diffpercent"] = regionDiffFromMeanPercentage;

if process123:
    regionLabel="BeringSea789";
    regionMeans, regionSDs, regionDiffFromMean, regionDiffFromMeanSD, regionDiffFromMeanPercentage = calc_region_err_stats(data789, laruBeringSeaMask, regionLabel, errorsGeoEk=geoEksRMSEs, errorsStokes=stokesRMSEs, numReps=numRepeats, distances=distances);
    mainDF[regionLabel+"_mean"] = regionMeans;
    mainDF[regionLabel+"SD"] = regionSDs;
    mainDF[regionLabel+"_diff"] = regionDiffFromMean;
    mainDF[regionLabel+"_diffSD"] = regionDiffFromMeanSD;
    mainDF[regionLabel+"_diffpercent"] = regionDiffFromMeanPercentage;


############################
#Antarctic Peninsula
data123 = get_processed_data2([0,1,2], allData, [mf.area_AntarcticPeninsula], params);
data789 = get_processed_data2([6,7,8], allData, [mf.area_AntarcticPeninsula], params);

distances = np.array(gridwiseData["distance"])[data123.regionMask]/1000.0; #in km
dists, laruAntarcticPeninsulaMask = mf.calculate_km_subsection_bounds_along_shelf(data123.xcoords, data123.ycoords, distances, mf.laruelle_AntarcticPeninsula, params, testPlot=False);

#try:
regionLabel="AntarcticPeninsula123";
regionMeans, regionSDs, regionDiffFromMean, regionDiffFromMeanSD, regionDiffFromMeanPercentage = calc_region_err_stats(data123, laruAntarcticPeninsulaMask, regionLabel, errorsGeoEk=geoEksRMSEs, errorsStokes=stokesRMSEs, numReps=numRepeats, distances=distances);
mainDF[regionLabel+"_mean"] = regionMeans;
mainDF[regionLabel+"SD"] = regionSDs;
mainDF[regionLabel+"_diff"] = regionDiffFromMean;
mainDF[regionLabel+"_diffSD"] = regionDiffFromMeanSD;
mainDF[regionLabel+"_diffpercent"] = regionDiffFromMeanPercentage;

if process123:
    regionLabel="AntarcticPeninsula789";
    regionMeans, regionSDs, regionDiffFromMean, regionDiffFromMeanSD, regionDiffFromMeanPercentage = calc_region_err_stats(data789, laruAntarcticPeninsulaMask, regionLabel, errorsGeoEk=geoEksRMSEs, errorsStokes=stokesRMSEs, numReps=numRepeats, distances=distances);
    mainDF[regionLabel+"_mean"] = regionMeans;
    mainDF[regionLabel+"SD"] = regionSDs;
    mainDF[regionLabel+"_diff"] = regionDiffFromMean;
    mainDF[regionLabel+"_diffSD"] = regionDiffFromMeanSD;
    mainDF[regionLabel+"_diffpercent"] = regionDiffFromMeanPercentage;

#except IndexError: #When distances and laruAntarcticPeninsulaMask are 0 because there were no shelf edges in the region, it fails.
#                     #However, print_statistics handles this situation correctly.
#    t2data.append(print_statistics(data789, distances, "Antarctic Peninsula (Jul-Sep)", laruAntarcticPeninsulaMask));
#    t2data.append(print_statistics(data789, distances, "Antarctic Peninsula (Jul-Sep)", laruAntarcticPeninsulaMask));

########################
#Labrador Sea
data123 = get_processed_data2([0,1,2], allData, [mf.area_LabradorSea], params);
data789 = get_processed_data2([6,7,8], allData, [mf.area_LabradorSea], params);
distances = np.array(gridwiseData["distance"])[data123.regionMask]/1000.0; #in km

#Custom roll the data to eliminate gap
rollAmount = -154;
data123.xcoords = np.roll(data123.xcoords, rollAmount);
data123.ycoords = np.roll(data123.ycoords, rollAmount);
data123.regionMask = np.roll(data123.regionMask, rollAmount);
data123.ontoShelfVectorX = np.roll(data123.ontoShelfVectorX, rollAmount);
data123.ontoShelfVectorY = np.roll(data123.ontoShelfVectorY, rollAmount);
data123.stokesMaskPass = np.roll(data123.stokesMaskPass, rollAmount);
data123.totalCurrentX = np.roll(data123.totalCurrentX, rollAmount);
data123.totalCurrentY = np.roll(data123.totalCurrentY, rollAmount);

if process123:
    data789.xcoords = np.roll(data789.xcoords, rollAmount);
    data789.ycoords = np.roll(data789.ycoords, rollAmount);
    data789.regionMask = np.roll(data789.regionMask, rollAmount);
    data789.ontoShelfVectorX = np.roll(data789.ontoShelfVectorX, rollAmount);
    data789.ontoShelfVectorY = np.roll(data789.ontoShelfVectorY, rollAmount);
    data789.stokesMaskPass = np.roll(data789.stokesMaskPass, rollAmount);
    data789.totalCurrentX = np.roll(data789.totalCurrentX, rollAmount);
    data789.totalCurrentY = np.roll(data789.totalCurrentY, rollAmount);

distances = np.roll(distances, rollAmount);

dists, laruLabradorSeaMask = mf.calculate_km_subsection_bounds_along_shelf(data123.xcoords, data123.ycoords, distances, mf.laruelle_LabradorSea, params, testPlot=False);

regionLabel="LabradorSea123";
regionMeans, regionSDs, regionDiffFromMean, regionDiffFromMeanSD, regionDiffFromMeanPercentage = calc_region_err_stats(data123, laruLabradorSeaMask, regionLabel, errorsGeoEk=geoEksRMSEs, errorsStokes=stokesRMSEs, numReps=numRepeats, distances=distances);
mainDF[regionLabel+"_mean"] = regionMeans;
mainDF[regionLabel+"SD"] = regionSDs;
mainDF[regionLabel+"_diff"] = regionDiffFromMean;
mainDF[regionLabel+"_diffSD"] = regionDiffFromMeanSD;
mainDF[regionLabel+"_diffpercent"] = regionDiffFromMeanPercentage;

if process123:
    regionLabel="LabradorSea789";
    regionMeans, regionSDs, regionDiffFromMean, regionDiffFromMeanSD, regionDiffFromMeanPercentage = calc_region_err_stats(data789, laruLabradorSeaMask, regionLabel, errorsGeoEk=geoEksRMSEs, errorsStokes=stokesRMSEs, numReps=numRepeats, distances=distances);
    mainDF[regionLabel+"_mean"] = regionMeans;
    mainDF[regionLabel+"SD"] = regionSDs;
    mainDF[regionLabel+"_diff"] = regionDiffFromMean;
    mainDF[regionLabel+"_diffSD"] = regionDiffFromMeanSD;
    mainDF[regionLabel+"_diffpercent"] = regionDiffFromMeanPercentage;


###################
#Tasmanian Shelf
data123 = get_processed_data2([0,1,2], allData, [mf.area_Tasmania], params);
data789 = get_processed_data2([6,7,8], allData, [mf.area_Tasmania], params);

distances = np.array(gridwiseData["distance"])[data123.regionMask]/1000.0; #in km
dists, laruTasmanianShelfMask = mf.calculate_km_subsection_bounds_along_shelf(data123.xcoords, data123.ycoords, distances, mf.laruelle_TasmanianShelf, params, testPlot=False);

regionLabel="TasmanianShelf123";
regionMeans, regionSDs, regionDiffFromMean, regionDiffFromMeanSD, regionDiffFromMeanPercentage = calc_region_err_stats(data123, laruTasmanianShelfMask, regionLabel, errorsGeoEk=geoEksRMSEs, errorsStokes=stokesRMSEs, numReps=numRepeats, distances=distances);
mainDF[regionLabel+"_mean"] = regionMeans;
mainDF[regionLabel+"SD"] = regionSDs;
mainDF[regionLabel+"_diff"] = regionDiffFromMean;
mainDF[regionLabel+"_diffSD"] = regionDiffFromMeanSD;
mainDF[regionLabel+"_diffpercent"] = regionDiffFromMeanPercentage;

if process123:
    regionLabel="TasmanianShelf789";
    regionMeans, regionSDs, regionDiffFromMean, regionDiffFromMeanSD, regionDiffFromMeanPercentage = calc_region_err_stats(data789, laruTasmanianShelfMask, regionLabel, errorsGeoEk=geoEksRMSEs, errorsStokes=stokesRMSEs, numReps=numRepeats, distances=distances);
    mainDF[regionLabel+"_mean"] = regionMeans;
    mainDF[regionLabel+"SD"] = regionSDs;
    mainDF[regionLabel+"_diff"] = regionDiffFromMean;
    mainDF[regionLabel+"_diffSD"] = regionDiffFromMeanSD;
    mainDF[regionLabel+"_diffpercent"] = regionDiffFromMeanPercentage;



###################
#Barents Sea
data123 = get_processed_data2([0,1,2], allData, [mf.area_BarentsSea], params);
data789 = get_processed_data2([6,7,8], allData, [mf.area_BarentsSea], params);

distances = np.array(gridwiseData["distance"])[data123.regionMask]/1000.0; #in km
dists, laruBarentsSeaMask = mf.calculate_km_subsection_bounds_along_shelf(data123.xcoords, data123.ycoords, distances, mf.laruelle_BarentsSea, params, testPlot=False);

regionLabel="BarentsSea123";
regionMeans, regionSDs, regionDiffFromMean, regionDiffFromMeanSD, regionDiffFromMeanPercentage = calc_region_err_stats(data123, laruBarentsSeaMask, regionLabel, errorsGeoEk=geoEksRMSEs, errorsStokes=stokesRMSEs, numReps=numRepeats, distances=distances);
mainDF[regionLabel+"_mean"] = regionMeans;
mainDF[regionLabel+"SD"] = regionSDs;
mainDF[regionLabel+"_diff"] = regionDiffFromMean;
mainDF[regionLabel+"_diffSD"] = regionDiffFromMeanSD;
mainDF[regionLabel+"_diffpercent"] = regionDiffFromMeanPercentage;

if process123:
    regionLabel="BarentsSea789";
    regionMeans, regionSDs, regionDiffFromMean, regionDiffFromMeanSD, regionDiffFromMeanPercentage = calc_region_err_stats(data789, laruBarentsSeaMask, regionLabel, errorsGeoEk=geoEksRMSEs, errorsStokes=stokesRMSEs, numReps=numRepeats, distances=distances);
    mainDF[regionLabel+"_mean"] = regionMeans;
    mainDF[regionLabel+"SD"] = regionSDs;
    mainDF[regionLabel+"_diff"] = regionDiffFromMean;
    mainDF[regionLabel+"_diffSD"] = regionDiffFromMeanSD;
    mainDF[regionLabel+"_diffpercent"] = regionDiffFromMeanPercentage;


###########################
#South Atlantic Bight
data123 = get_processed_data2([0,1,2], allData, [mf.area_SouthAtlanticBight], params);
data789 = get_processed_data2([6,7,8], allData, [mf.area_SouthAtlanticBight], params);

distances = np.array(gridwiseData["distance"])[data123.regionMask]/1000.0; #in km
dists, laruSouthAtlanticBightMask = mf.calculate_km_subsection_bounds_along_shelf(data123.xcoords, data123.ycoords, distances, mf.laruelle_SouthAtlanticBight, params, testPlot=False);

regionLabel="SouthAtlanticBight123";
regionMeans, regionSDs, regionDiffFromMean, regionDiffFromMeanSD, regionDiffFromMeanPercentage = calc_region_err_stats(data123, laruSouthAtlanticBightMask, regionLabel, errorsGeoEk=geoEksRMSEs, errorsStokes=stokesRMSEs, numReps=numRepeats, distances=distances);
mainDF[regionLabel+"_mean"] = regionMeans;
mainDF[regionLabel+"SD"] = regionSDs;
mainDF[regionLabel+"_diff"] = regionDiffFromMean;
mainDF[regionLabel+"_diffSD"] = regionDiffFromMeanSD;
mainDF[regionLabel+"_diffpercent"] = regionDiffFromMeanPercentage;

if process123:
    regionLabel="SouthAtlanticBight789";
    regionMeans, regionSDs, regionDiffFromMean, regionDiffFromMeanSD, regionDiffFromMeanPercentage = calc_region_err_stats(data789, laruSouthAtlanticBightMask, regionLabel, errorsGeoEk=geoEksRMSEs, errorsStokes=stokesRMSEs, numReps=numRepeats, distances=distances);
    mainDF[regionLabel+"_mean"] = regionMeans;
    mainDF[regionLabel+"SD"] = regionSDs;
    mainDF[regionLabel+"_diff"] = regionDiffFromMean;
    mainDF[regionLabel+"_diffSD"] = regionDiffFromMeanSD;
    mainDF[regionLabel+"_diffpercent"] = regionDiffFromMeanPercentage;


################################
#Southern Greenland
data123 = get_processed_data2([0,1,2], allData, [mf.area_SouthGreenland], params);
data789 = get_processed_data2([6,7,8], allData, [mf.area_SouthGreenland], params);

distances = np.array(gridwiseData["distance"])[data123.regionMask]/1000.0; #in km
dists, laruSouthernGreenlandMask = mf.calculate_km_subsection_bounds_along_shelf(data123.xcoords, data123.ycoords, distances, mf.laruelle_SouthernGreenland, params, testPlot=False);

regionLabel="SouthernGreenland123";
regionMeans, regionSDs, regionDiffFromMean, regionDiffFromMeanSD, regionDiffFromMeanPercentage = calc_region_err_stats(data123, laruSouthernGreenlandMask, regionLabel, errorsGeoEk=geoEksRMSEs, errorsStokes=stokesRMSEs, numReps=numRepeats, distances=distances);
mainDF[regionLabel+"_mean"] = regionMeans;
mainDF[regionLabel+"SD"] = regionSDs;
mainDF[regionLabel+"_diff"] = regionDiffFromMean;
mainDF[regionLabel+"_diffSD"] = regionDiffFromMeanSD;
mainDF[regionLabel+"_diffpercent"] = regionDiffFromMeanPercentage;

if process123:
    regionLabel="SouthernGreenland789";
    regionMeans, regionSDs, regionDiffFromMean, regionDiffFromMeanSD, regionDiffFromMeanPercentage = calc_region_err_stats(data789, laruSouthernGreenlandMask, regionLabel, errorsGeoEk=geoEksRMSEs, errorsStokes=stokesRMSEs, numReps=numRepeats, distances=distances);
    mainDF[regionLabel+"_mean"] = regionMeans;
    mainDF[regionLabel+"SD"] = regionSDs;
    mainDF[regionLabel+"_diff"] = regionDiffFromMean;
    mainDF[regionLabel+"_diffSD"] = regionDiffFromMeanSD;
    mainDF[regionLabel+"_diffpercent"] = regionDiffFromMeanPercentage;


#######################
#Cascadian Shelf
data123 = get_processed_data2([0,1,2], allData, [mf.area_CascianShelf], params);
data789 = get_processed_data2([6,7,8], allData, [mf.area_CascianShelf], params);

distances = np.array(gridwiseData["distance"])[data123.regionMask]/1000.0; #in km
dists, laruCascadianShelfMask = mf.calculate_km_subsection_bounds_along_shelf(data123.xcoords, data123.ycoords, distances, mf.laruelle_CascadianShelf, params, testPlot=False);

regionLabel="CascadianShelf123";
regionMeans, regionSDs, regionDiffFromMean, regionDiffFromMeanSD, regionDiffFromMeanPercentage = calc_region_err_stats(data123, laruCascadianShelfMask, regionLabel, errorsGeoEk=geoEksRMSEs, errorsStokes=stokesRMSEs, numReps=numRepeats, distances=distances);
mainDF[regionLabel+"_mean"] = regionMeans;
mainDF[regionLabel+"SD"] = regionSDs;
mainDF[regionLabel+"_diff"] = regionDiffFromMean;
mainDF[regionLabel+"_diffSD"] = regionDiffFromMeanSD;
mainDF[regionLabel+"_diffpercent"] = regionDiffFromMeanPercentage;

if process123:
    regionLabel="CascadianShelf789";
    regionMeans, regionSDs, regionDiffFromMean, regionDiffFromMeanSD, regionDiffFromMeanPercentage = calc_region_err_stats(data789, laruCascadianShelfMask, regionLabel, errorsGeoEk=geoEksRMSEs, errorsStokes=stokesRMSEs, numReps=numRepeats, distances=distances);
    mainDF[regionLabel+"_mean"] = regionMeans;
    mainDF[regionLabel+"SD"] = regionSDs;
    mainDF[regionLabel+"_diff"] = regionDiffFromMean;
    mainDF[regionLabel+"_diffSD"] = regionDiffFromMeanSD;
    mainDF[regionLabel+"_diffpercent"] = regionDiffFromMeanPercentage;


############################
#Irminger Sea
data123 = get_processed_data2([0,1,2], allData, [mf.area_IrmingerSea], params);
data789 = get_processed_data2([6,7,8], allData, [mf.area_IrmingerSea], params);

distances = np.array(gridwiseData["distance"])[data123.regionMask]/1000.0; #in km
dists, laruIrmingerSeaMask = mf.calculate_km_subsection_bounds_along_shelf(data123.xcoords, data123.ycoords, distances, mf.laruelle_IrmingerSea, params, testPlot=False);

regionLabel="IrmingerSea123";
regionMeans, regionSDs, regionDiffFromMean, regionDiffFromMeanSD, regionDiffFromMeanPercentage = calc_region_err_stats(data123, laruIrmingerSeaMask, regionLabel, errorsGeoEk=geoEksRMSEs, errorsStokes=stokesRMSEs, numReps=numRepeats, distances=distances);
mainDF[regionLabel+"_mean"] = regionMeans;
mainDF[regionLabel+"SD"] = regionSDs;
mainDF[regionLabel+"_diff"] = regionDiffFromMean;
mainDF[regionLabel+"_diffSD"] = regionDiffFromMeanSD;
mainDF[regionLabel+"_diffpercent"] = regionDiffFromMeanPercentage;

if process123:
    regionLabel="IrmingerSea789";
    regionMeans, regionSDs, regionDiffFromMean, regionDiffFromMeanSD, regionDiffFromMeanPercentage = calc_region_err_stats(data789, laruIrmingerSeaMask, regionLabel, errorsGeoEk=geoEksRMSEs, errorsStokes=stokesRMSEs, numReps=numRepeats, distances=distances);
    mainDF[regionLabel+"_mean"] = regionMeans;
    mainDF[regionLabel+"SD"] = regionSDs;
    mainDF[regionLabel+"_diff"] = regionDiffFromMean;
    mainDF[regionLabel+"_diffSD"] = regionDiffFromMeanSD;
    mainDF[regionLabel+"_diffpercent"] = regionDiffFromMeanPercentage;


#################
#Write to file
#dft1 = pd.DataFrame(t1data);
#dft1.columns = ["region", "total across-shelf", "total across-shelf SD", "proportionEks", "proportionEksSD", "proportionGeo", "proportionGeoSD", "proportionStokes", "proportionStokesSD", "percentEks", "percentEksSD", "percentGeo", "percentGeoSD", "percentStokes", "percentStokesSD"];
mainDF.to_csv(path.join(outputPath, "regional_error_analysis_"+params.paramsetName+"_v4_"+str(numRepeats)+".csv"), index=False);

import datetime;
print datetime.datetime.now();
#One repeat: 15:13 started. 15:25:48 finished. #Time: ~13 minutes. Mostly reading in of data...
#300 repeats: 15:32 started.

