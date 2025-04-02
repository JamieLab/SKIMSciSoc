#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 15:06:22 2019

@author: Tom Holding
"""

import resample_bathymetry;
import shelf_coord_creator;
import generate_cell_data;
import calculate_shelf_currents;
import calculate_shelf_currents_errprop;
from os import path;
import os; #for mkdir



##############################################
# Calculates across-shelf current components
#   Params - parameter set to use from parameter_sets
# Execution switches to control which parts of the analysis will run:
#   generateBathymetry
#   generateShelfCoordinates
#   generateCellData
#   calculateShelfCurrents
#   withGasTransferCalc will calculate k along with shelf currents
def shelf_current_analysis(params, generateBathymetry=False, generateShelfCoordinates=False, generateCellData=False, calculateShelfCurrents=False, withGasTransferCalc=True, outputPathOverride=None, errorPropagation=None,cmems='False',shape=False,shape_file_shallow='',shape_file_deep=''):
    #Global settings
    if outputPathOverride != None:
        runDataPath = outputPathOverride;
    else:
        runDataPath = path.join("data", params.paramsetName);

    #Create required folders
    if path.exists(runDataPath) == False:
        os.makedirs(runDataPath);
    if path.exists(path.join(runDataPath, "shelf_coordinates")) == False:
        os.makedirs(path.join(runDataPath, "shelf_coordinates"));
    if path.exists(path.join(runDataPath, "gridwise_data")) == False:
        os.makedirs(path.join(runDataPath, "gridwise_data"));
    if path.exists(path.join(runDataPath, "current_data")) == False:
        os.makedirs(path.join(runDataPath, "current_data"));



    ###################
    # Make bathymetry #
    ###################
    if generateBathymetry:
        print "Resampling bathymetry.";
        gebcoDataPath = "D:/Data/Bathymetry/GEBCO_2014_2D.nc";
        #gebcoDataPath = "/home/rr/data/GEBCO_bathymetry_30sec/GEBCO_2014_2D.nc";
        resample_bathymetry.generate_resampled_bathymetry(gebcoDataPath, params.pixelRes, testPlot=True, outputPath=runDataPath);


    #########################################
    # Generate shelf-edge coordinates file. #
    #########################################
    if generateShelfCoordinates:
        print "Generating shelf path coordinate data.";
        #Create contour plots and select paths to use as the shelf edge coordinates.
        shelf_coord_creator.get_shelf_edge_data(params, testPlots=True, outputPath=runDataPath,shape=shape,shape_file_shallow=shape_file_shallow,shape_file_deep=shape_file_deep);


    ######################
    # Generate cell data #
    ######################
    if generateCellData:
        generate_cell_data.generate_cell_data(params, runDataPath, testPlots=True, verbose=True, outputPath=runDataPath);


    ###################################
    # Calculate across-shelf currents #
    ###################################
    if calculateShelfCurrents:
        if errorPropagation == None:
            calculate_shelf_currents.calculate_shelf_current_data(params, runDataPath, calculateGasTransferVelocity=withGasTransferCalc, testPlot=True, verbose=True, outputPath=runDataPath,cmems = cmems);
        else:
            calculate_shelf_currents_errprop.calculate_shelf_current_data(params, runDataPath, calculateGasTransferVelocity=withGasTransferCalc, testPlot=True, verbose=True, outputPath=runDataPath, errorProp=errorPropagation);
