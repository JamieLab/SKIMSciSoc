#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 20 11:24:33 2018

@author: rr
"""

from netCDF4 import Dataset;
import numpy as np;
import matplotlib.pyplot as plt;
import skim_utilities as su;
import parameter_sets as ps;
import pandas as pd;
import cPickle as pickle;
from string import Template;
from os import path;

def _calc_month_years(params):
    years = range(params.start_year, params.end_year+1);
    monthYears = [];
    for year in years:
        for month in range(0, 12):
            if year == params.start_year and month < params.start_month:
                continue;
            elif year == params.end_year and month >= params.end_month:
                break;
            else:
                monthYears.append( (month, year) );
    return monthYears;


def _calc_time_index(iYear, iMonth):
    return 12*iYear + iMonth;


#magnitude of across-shelf geostrophic current
def _calculate_generic_across_shelf_current(geostrophicCurrentVector, ontoShelfVector):
    pointCurrent = np.dot(ontoShelfVector, geostrophicCurrentVector) / np.linalg.norm(ontoShelfVector);
    return pointCurrent;


def _calculate_across_shelf_ekman_current(meanCurrentVector, ontoShelfVector, windspeedVector, latDegrees):
    #Calculate direction of ekman transport (90 degree from wind direction).
    
    #transportDirection = su.get_surface_ekman_direction_vector(windspeedVector, latDegrees); #90 degrees from windspeed
    transportDirection = su.get_integrated_ekman_direction(meanCurrentVector, latDegrees); #45 degrees from surface ekman current
    
    transportMagnitude = np.linalg.norm(meanCurrentVector);
    transportVector = transportDirection*transportMagnitude; #Direction dictated by wind speed and latitude, magnitude by estimated transport
    
    #On-shelf flux at a single surface point
    pointCurrent = np.dot(ontoShelfVector, transportVector) / np.linalg.norm(ontoShelfVector)
    
    #return surfaceFlux;
    return pointCurrent;

#windspeed in metres per second
#Smith 1980 C_D parameterisation, as described in Yelland1994 (introduction)
def _calculate_drag_coefficient(windspeed):
    return 0.00044 + 0.000063*windspeed;

def _calculate_wind_stress(windspeedVector, alongShelfDirection=None):
    airDensity = 1.225; #kg m^-3
    if type(alongShelfDirection) != type(None): #Windstress in the along-shelf direction
        windspeed = np.dot(windspeedVector, alongShelfDirection) / np.linalg.norm(alongShelfDirection);
    else: #Overall wind stress
        windspeed = su.calculate_magnitude(windspeedVector[0], windspeedVector[1]);
    dragCoefficient = _calculate_drag_coefficient(windspeed);
    windStress = airDensity*dragCoefficient*(windspeed**2.0);
    return windStress;


def calculate_shelf_current_data(params, inputDataPath, calculateGasTransferVelocity=False, testPlot=False, verbose=False, outputPath=None):
    #Get time span as (month, year) pairs
    monthYears = _calc_month_years(params);

    #read in gridwide data, convert x and y coordinates to proper integers.
    gridCellDataPath = path.join(inputDataPath, "gridwise_data/per_grid_cell_edge_data_"+params.paramsetName+".csv");
    gridwiseData = pd.read_table(gridCellDataPath, sep=',');
    gridwiseData.x = gridwiseData.x.astype(int);
    gridwiseData.y = gridwiseData.y.astype(int);
    
    
    #main loop to calculate for each month, year pair.
    data = [];
    dataUncertaintyAnalysis = [];
    for month, year in monthYears:
        #iYear = year - params.start_year;
        monthStr = format(month+1, "02d");
        print year, month;
        
        #read current, wind and wave height files for this (month,year)
        ekmanNc = Dataset(params.ekmanTemplate.safe_substitute(YYYY=str(year), MM=monthStr), 'r');
        geostrophicNc = Dataset(params.geostrophicTemplate.safe_substitute(YYYY=str(year), MM=monthStr), 'r');
        stokesNc = Dataset(params.stokesTemplate.safe_substitute(YYYY=str(year), MM=monthStr), 'r');
        wndNc = Dataset(params.wavewatchWndTemplate.safe_substitute(YYYY=str(year), MM=monthStr), 'r');
        hsNc = Dataset(params.wavewatchHsTemplate.safe_substitute(YYYY=str(year), MM=monthStr), 'r');
        
        #Extract the required variables from the netCDF data and apply region masks
        #data must be present in a consistent format: a 0.25 by 0.25 degree grid, using only the bounded local area.
        ekmanEast = su.apply_mask(ekmanNc.variables["ekmanU"][:, :], None, params.ilatRange, params.ilonRange);
        ekmanNorth = su.apply_mask(ekmanNc.variables["ekmanV"][:, :], None, params.ilatRange, params.ilonRange);
        geostrophicEast = su.apply_mask(geostrophicNc.variables["geostrophicU"][:, :], None, params.ilatRange, params.ilonRange);
        geostrophicNorth = su.apply_mask(geostrophicNc.variables["geostrophicV"][:, :], None, params.ilatRange, params.ilonRange);
        stokesEast = su.apply_mask(stokesNc.variables["stokesU"][:, :], None, params.ilatRange, params.ilonRange);
        stokesNorth = su.apply_mask(stokesNc.variables["stokesV"][:, :], None, params.ilatRange, params.ilonRange);
        uwnd = su.apply_mask(wndNc.variables["uwnd_mean"][:, :], None, params.ilatRange, params.ilonRange);
        vwnd = su.apply_mask(wndNc.variables["vwnd_mean"][:, :], None, params.ilatRange, params.ilonRange);
        hs = su.apply_mask(hsNc.variables["hs_mean"][:, :], None, params.ilatRange, params.ilonRange);
        
        #If calculating k we need to read in sst
        if calculateGasTransferVelocity:
            sstNc = Dataset(params.reynoldsSSTTemplate.safe_substitute(YYYY=str(year), MM=monthStr), 'r');
            sst = np.flipud(sstNc.variables["sst_mean"][0, :, :]);
            #check sst is same orientation as currents! plot both... check resolutions match too (0.25)...
            sst = su.apply_mask(sst, None, params.ilatRange, params.ilonRange);
            wndMag = np.sqrt(uwnd**2 + vwnd**2);
            wndMoment2 = wndMag**2;
            
        #Extra data required for the european shelf case study
        #Read the netCDF variables and apply region masks
        if params.paramsetName == "europeanshelf":
            skimulatorNc = Dataset(params.skimulatorTemplate.safe_substitute(YYYY=str(year), MM=monthStr), 'r');
            skimEast = su.apply_mask(skimulatorNc.variables["Ux_skim"][:, :], None, params.ilatRange, params.ilonRange);
            skimNorth = su.apply_mask(skimulatorNc.variables["Uy_skim"][:, :], None, params.ilatRange, params.ilonRange);
            eksEast = su.apply_mask(skimulatorNc.variables["u_eks"][:, :], None, params.ilatRange, params.ilonRange);
            eksNorth = su.apply_mask(skimulatorNc.variables["v_eks"][:, :], None, params.ilatRange, params.ilonRange);
            altiEast = su.apply_mask(skimulatorNc.variables["Ux_alti"][:, :], None, params.ilatRange, params.ilonRange);
            altiNorth = su.apply_mask(skimulatorNc.variables["Uy_alti"][:, :], None, params.ilatRange, params.ilonRange);
            truthEast = su.apply_mask(skimulatorNc.variables["Ux_truth"][:, :], None, params.ilatRange, params.ilonRange);
            truthNorth = su.apply_mask(skimulatorNc.variables["Uy_truth"][:, :], None, params.ilatRange, params.ilonRange);
    
        
        #####################
        #Begin calculating shelf current data.
        monthData = [];
        monthDataUncertaintyAnalysis = []; #attach data for uncertainty analysis
        #ekmanBoundaryTotal = 0.0;
        #geostrophicBoundaryTotal = 0.0;
        #stokesBoundaryTotal = 0.0;
        #filterInfo = [];
        #tTotal = 0.0;
        #Each dataframe row corresponds to a grid location. Use the grid data to calculate grid specific current data
        # and store in 'cellData'
        for i, cell in enumerate(gridwiseData.itertuples()):
            cellData = su.QuickStruct();
            cellDataUncertaintyAnalysis = su.QuickStruct(); #Extra data for the uncertainty analysis.
            #cellData.row = cell; #More info for debugging
            cellData.indexX = cell.x;
            cellData.indexY = cell.y;
            cellDataUncertaintyAnalysis.indexX = cell.x;
            cellDataUncertaintyAnalysis.indexY = cell.y;
            
            #print cell;
            indexX = cell.x;
            indexY = cell.y;
            latDegrees = cell.lat;
            ontoShelfVector = np.array([cell.onshelfX, cell.onshelfY]);
            segmentDistance = cell.distance;
            sigWaveHeight = hs[indexY, indexX];
            
            ekmanCurrentVector = np.array([ekmanEast[indexY, indexX], ekmanNorth[indexY, indexX]])
            geostrophicCurrentVector = np.array([geostrophicEast[indexY, indexX], geostrophicNorth[indexY, indexX]]);
            stokesCurrentVector = np.array([stokesEast[indexY, indexX], stokesNorth[indexY, indexX]]);
            windspeedVector = np.array([uwnd[indexY, indexX], vwnd[indexY, indexX]]);
            
            cellData.segmentDistance = segmentDistance;
            
            
            ####################### Calculating values used for plotting #######################
            #Only calculate if there are actually values at this position
            if np.all(np.isfinite(ekmanCurrentVector)) and np.all(np.isfinite(ontoShelfVector)) and np.all(np.isfinite(windspeedVector)) \
              and np.all(np.isfinite(geostrophicCurrentVector)) and np.all(np.isfinite(stokesCurrentVector)):
                #Calculate the across-shelf current for each component
                nEkmanAcrossShelf = _calculate_across_shelf_ekman_current(ekmanCurrentVector, ontoShelfVector, windspeedVector, latDegrees);
                cellData.nEkmanAcrossShelf = nEkmanAcrossShelf;
                
                nGeostrophicAcrossShelf = _calculate_generic_across_shelf_current(geostrophicCurrentVector, ontoShelfVector);
                cellData.nGeostrophicAcrossShelf = nGeostrophicAcrossShelf;
                
                nStokesAcrossShelf = _calculate_generic_across_shelf_current(stokesCurrentVector, ontoShelfVector);
                cellData.nStokesAcrossShelf = nStokesAcrossShelf;
                
                #Store the abs (magnitude) of each across-shelf component
                nEkmanMagnitude = np.abs(nEkmanAcrossShelf); #absolute value of the magnitude across-shelf ekman current
                cellData.nEkmanMagnitude = nEkmanMagnitude;
                
                nGeostrophicMagnitude = np.abs(nGeostrophicAcrossShelf); #absolute value of the magnitude of across-shelf geostrophic current
                cellData.nGeostrophicMagnitude = nGeostrophicMagnitude;
                
                nStokesMagnitude = np.abs(nStokesAcrossShelf);
                cellData.nStokesMagnitude = nStokesMagnitude;
                
                #Calculate the Ekman and geostrophic proportions of ekman+geostrophic
                cellData.ekmanProportionGeoEk = nEkmanMagnitude / (nEkmanMagnitude+nGeostrophicMagnitude);
                cellData.geostrophicProportionGeoEk = nGeostrophicMagnitude / (nEkmanMagnitude+nGeostrophicMagnitude);
                
                #Calculate total magnitude (ekman+geostrophic+stokes) and stokes proportion of total
                nTotalMagnitude = nEkmanMagnitude+nGeostrophicMagnitude+nStokesMagnitude;
                cellData.nTotalMagnitude = nTotalMagnitude;
                cellData.stokesProportionOfTotal = nStokesMagnitude / nTotalMagnitude;
                
                #Calculate whether stokes is likely to be important (i.e. where the mixed layer depth is shallow so won't be swamped at depth by other components)
                windStress = _calculate_wind_stress(windspeedVector);
                cellData.windStress = windStress;
                cellData.sigWaveHeight = sigWaveHeight;
                cellData.stokesMaskPass = (np.abs(windStress) < params.stokesMaskWindStressThreshold) and (sigWaveHeight > params.stokesMaskSigWaveHeightThreshold);
                
                cellDataUncertaintyAnalysis.indexX = cell.x;
                cellDataUncertaintyAnalysis.indexY = cell.y;
                cellDataUncertaintyAnalysis.ontoShelfVector = ontoShelfVector;
                cellDataUncertaintyAnalysis.segmentDistance = segmentDistance;
                cellDataUncertaintyAnalysis.geostrophicCurrentVector = geostrophicCurrentVector;
                cellDataUncertaintyAnalysis.stokesCurrentVector = stokesCurrentVector;
                cellDataUncertaintyAnalysis.windspeedVector = windspeedVector;
                cellDataUncertaintyAnalysis.stokesMaskPass = cellData.stokesMaskPass;
                
                #modify ekman current direction based on corialis
                ekmanDirection = su.get_integrated_ekman_direction(ekmanCurrentVector, latDegrees); #45 degrees from surface ekman current
                ekmanMagnitude = np.linalg.norm(ekmanCurrentVector);
                adjustedEkmanVector = ekmanDirection*ekmanMagnitude; #Direction dictated by wind speed and latitude, magnitude by estimated transport
                cellDataUncertaintyAnalysis.ekmanCurrentVector = adjustedEkmanVector;
                
            else:
                cellData.nEkmanAcrossShelf = np.nan;
                cellData.nGeostrophicAcrossShelf = np.nan;
                cellData.nStokesAcrossShelf = np.nan;
                cellData.nEkmanMagnitude = np.nan;
                cellData.nGeostrophicMagnitude = np.nan;
                cellData.nStokesMagnitude = np.nan;
                cellData.ekmanProportionGeoEk = np.nan;
                cellData.geostrophicProportionGeoEk = np.nan;
                cellData.nTotalMagnitude = np.nan;
                cellData.stokesProportionOfTotal = np.nan;
                cellData.windStress = np.nan;
                cellData.sigWaveHeight = np.nan;
                cellData.stokesMaskPass = np.nan;
                
                cellDataUncertaintyAnalysis.indexX = np.nan;
                cellDataUncertaintyAnalysis.indexY = np.nan;
                cellDataUncertaintyAnalysis.ontoShelfVector = np.nan;
                cellDataUncertaintyAnalysis.segmentDistance = np.nan;
                cellDataUncertaintyAnalysis.ekmanCurrentVector = np.nan;
                cellDataUncertaintyAnalysis.geostrophicCurrentVector = np.nan;
                cellDataUncertaintyAnalysis.stokesCurrentVector = np.nan;
                cellDataUncertaintyAnalysis.windspeedVector = np.nan;
                cellDataUncertaintyAnalysis.stokesMaskPass = np.nan;
            
            
            #European case study requires some values calculated from the skimulator and reference data
            if params.paramsetName == "europeanshelf":
                #Create current vectors from the skimulator dataset
                eksVector = np.array([eksEast[indexY, indexX], eksNorth[indexY, indexX]], dtype=float);
                altiVector = np.array([altiEast[indexY, indexX], altiNorth[indexY, indexX]], dtype=float); #skimulater dataset geostrophic
                skimVector = np.array([skimEast[indexY, indexX], skimNorth[indexY, indexX]], dtype=float);
                truthVector = np.array([truthEast[indexY, indexX], truthNorth[indexY, indexX]], dtype=float); #skimulator dataset total?
    
                #Only calculate if there are actually values at this position
                if np.all(np.isfinite(eksVector)) and np.all(np.isfinite(altiVector)) and np.all(np.isfinite(skimVector)) \
                  and np.all(np.isfinite(truthVector)) and np.all(np.isfinite(ontoShelfVector)) and np.all(np.isfinite(windspeedVector)):
                    #repeat calculations with skimulator data...
                    nEksAcrossShelf = _calculate_across_shelf_ekman_current(eksVector, ontoShelfVector, windspeedVector, latDegrees);
                    cellData.nEksAcrossShelf = nEksAcrossShelf;
                    
                    nAltiAcrossShelf = _calculate_generic_across_shelf_current(altiVector, ontoShelfVector); #geostrophic
                    cellData.nAltiAcrossShelf = nAltiAcrossShelf;
                    
                    nSkimAcrossShelf = _calculate_generic_across_shelf_current(skimVector, ontoShelfVector); #Where is ekman???
                    cellData.nSkimAcrossShelf = nSkimAcrossShelf;
                    
                    nTruthAcrossShelf = _calculate_generic_across_shelf_current(truthVector, ontoShelfVector); #total reference?
                    cellData.nTruthAcrossShelf = nTruthAcrossShelf;
                else:
                    cellData.nEksAcrossShelf = np.nan;
                    cellData.nAltiAcrossShelf = np.nan;
                    cellData.nSkimAcrossShelf = np.nan;
                    cellData.nTruthAcrossShelf = np.nan;

            #Calculate Nightingale gas transfer velocity
            if calculateGasTransferVelocity:
                schmidt = 2116.8 + (-136.25 * sst[indexY, indexX]) + (4.7353 * sst[indexY, indexX]**2) + (-0.092307 * sst[indexY, indexX]**3) + (0.0007555 * sst[indexY, indexX]**4);
                k = 0.222 * wndMoment2[indexY, indexX] + (0.333 * wndMag[indexY, indexX]) * np.sqrt(600.0/schmidt);
                cellData.k = k;
                cellData.schmidt = schmidt;
    
            #Only store if either the skimulator or surface current data was used to calculate shelf current data.
            monthData.append(cellData);
            monthDataUncertaintyAnalysis.append(cellDataUncertaintyAnalysis); #Extra data for the uncertainty analysis.
        
        data.append(monthData);
        dataUncertaintyAnalysis.append(monthDataUncertaintyAnalysis);
        print "Data len:", len(data);
    
    if outputPath != None:
        outputFile = path.join(inputDataPath, "current_data", "surface_currents_"+params.paramsetName+".p");
        print "Pickling output to", outputFile;
        pickle.dump(data, open(outputFile, 'wb'));
        
        outputFile = path.join(inputDataPath, "current_data", "surface_currents_uncertainty_aid_"+params.paramsetName+".p");
        print "Pickling output data to aid uncertainty analysis to", outputFile;
        pickle.dump(dataUncertaintyAnalysis, open(outputFile, 'wb'));
    
    return data;