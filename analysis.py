#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 13:06:53 2019

@author: Tom Holding
Modified by Daniel J. Ford (d.ford@exeter.ac.uk)
Date: 03/2023
Changes:
- Updated output directory
"""

import copy;
import python_util.parameter_sets as ps;
import python_util.analysis_driver;
import python_util.skim_utilities as su;

import python_util.calc_means;

def transform_params(_params, shallowContourDepth):#, errorLevel, repeat):
    params = copy.deepcopy(_params);
    params.paramsetName = params.paramsetName+"_"+str(shallowContourDepth)+"m";
    #params.paramsetNameErrorOutput = params.paramsetName+"_"+format(errorLevel, ".3f")+"_rep"+str(repeat);
    params.contourPathFile = params.contourPathFile.split(".")[0]+"_"+str(shallowContourDepth)+"m.p";
    params.contourPathFileDeep = params.contourPathFileDeep.split(".")[0]+"_"+str(shallowContourDepth)+"m.p";
    params.shallowDepth = float(shallowContourDepth);
    params.deepDepth = params.shallowDepth+100.0;
    return params;


######################################
### Test run: 500m shallow, 600m deep.
### Extracts data for methods plots
# params = ps.get_example_method_plot_params();
# python_util.analysis_driver.shelf_current_analysis(params, generateBathymetry=True, generateShelfCoordinates=True, generateCellData=True, calculateShelfCurrents=True,outputPathOverride='D:/SKIM');
# params = ps.get_example_method_plot_params_64();
# python_util.analysis_driver.shelf_current_analysis(params, generateBathymetry=True, generateShelfCoordinates=True, generateCellData=True, calculateShelfCurrents=True,outputPathOverride='D:/SKIM');

###Unused?
# params = ps.get_global_params();
# python_util.analysis_driver.shelf_current_analysis(params, generateBathymetry=True, generateShelfCoordinates=True, generateCellData=True, calculateShelfCurrents=True,outputPathOverride='D:/SKIM');
# #python_util.calc_means.calc_mean_data(params, calculateTable1=False, calculateTable2=True, calculateTable2GasTransferVelocity=True, verbose=True);
# python_util.calc_means.calc_mean_data(params, calculateTable1=False, calculateTable2=True, calculateTable2GasTransferVelocity=True, verbose=True);



######################################
### European shelf case study run:
#params = ps.get_European_shelf_params();
#python_util.analysis_driver.shelf_current_analysis(params, generateBathymetry=False, generateShelfCoordinates=True, generateCellData=True, calculateShelfCurrents=True);



# #####################################
# # Sensitivity to shelf contour depth.
masterParams = ps.get_global_params(cmems=True);
# masterParams = ps.get_global_params_glory(res=False)
# contourDepths = [500];#
contourDepths = [300, 400, 500, 600, 700, 800];
#contourDepths = [600, 700, 800];
# errorLevels = [0.05, 0.10, 0.135, 0.15, 0.20];
# numRepeats = 3;
#
### Uncomment as neccessary for multiple repeats and uncertainty analysis

#for repeat in range(numRepeats):
for depth in contourDepths:
    #if (repeat == 0) and (depth == 300):
    #        continue;
    #for errorLevel in errorLevels:
        params = transform_params(masterParams, depth);#, errorLevel, repeat);
        print params.paramsetName;
        python_util.analysis_driver.shelf_current_analysis(params, generateBathymetry=False, generateShelfCoordinates=True, generateCellData=True, calculateShelfCurrents=True, withGasTransferCalc=True,outputPathOverride='E:/SKIM',cmems='True');
        # python_util.calc_means.calc_mean_data(params, calculateTable1=True, calculateTable2=False, calculateTable2GasTransferVelocity=True, verbose=True,outputPath='E:/SKIM');
        #python_util.calc_means.calc_mean_data(params, calculateTable1=False, calculateTable2=True, calculateTable2GasTransferVelocity=True, verbose=True);

        ##errorProp = [(errorLevel, errorLevel, 0.0)]; #tupples for (Ekman, geostrophic, stokes)
        ##python_util.analysis_driver.shelf_current_analysis(params, generateBathymetry=False, generateShelfCoordinates=False, generateCellData=False, calculateShelfCurrents=True, withGasTransferCalc=True, errorPropagation=errorProp);
        ##python_util.calc_means.calc_mean_data(params, calculateTable1=False, calculateTable2=True, calculateTable2GasTransferVelocity=True, verbose=True);

# params = transform_params(masterParams, 300);#, errorLevel, repeat);
# params.paramsetName = 'CMEMS_Jenny_MAB'
# # params.minContourPathSizeShallow = 3
# # params.minContourPathSizeDeep = 3
# # params.numLineApproximationFunction = su.fixed_length_lines;
# params.start_year = 1993;
# params.end_year = 2016;
# print params.paramsetName;
# python_util.analysis_driver.shelf_current_analysis(params, generateBathymetry=False, generateShelfCoordinates=True, generateCellData=True, calculateShelfCurrents=True, withGasTransferCalc=True,outputPathOverride='E:/SKIM',cmems='True',shape=False,shape_file_shallow = 'E:/SKIM_Paper_Data/SKIM/MAB_300m_Contour_v3/MAB_300m_Contour_v3.shp',shape_file_deep = 'E:/SKIM_Paper_Data/SKIM/MAB_400m_Contour_v3/MAB_400m_Contour_v3.shp');
