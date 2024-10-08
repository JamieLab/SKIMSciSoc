#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 20 11:23:51 2018

@author: rr
"""

import matplotlib.pyplot as plt;
from matplotlib.gridspec import GridSpec
from netCDF4 import Dataset;
import numpy as np;
import skim_utilities as su;
import parameter_sets as ps;
import cPickle as pickle;
from os import path;
import geopy.distance;
from random import random;
import pandas as pd;


def _flatten_shelf_coordinate_data(shelfCoordData, lineXCoordsList, lineYCoordsList, lineParamsList, pointLineIndexList):
    newData = su.QuickStruct();

    #Simple attributes just copied over
    newData.pathIndicesUsed = shelfCoordData.pathIndicesUsed;
    newData.contourDepth = shelfCoordData.contourDepth;

    #lists of lists concatinated
    newData.xs = [item for sublist in shelfCoordData.xs for item in sublist];
    newData.ys = [item for sublist in shelfCoordData.ys for item in sublist];
    newData.lineXCoords = [item for sublist in lineXCoordsList for item in sublist];
    newData.lineYCoords = [item for sublist in lineYCoordsList for item in sublist];
    newData.lineParamsList = [item for sublist in lineParamsList for item in sublist];

    #list of arrays stacked
    newData.coordinates = np.vstack(shelfCoordData.coordinatesLists);

    #pointLineIndex must keep track of line numbers to generate a cumulative sum
    newData.pointLineIndex = [];
    offset = 0;
    for sublist in pointLineIndexList:
        for val in sublist:
            newData.pointLineIndex.append(val+offset);
        offset = newData.pointLineIndex[-1]+1;

    return newData;

def _unflatten(pointLineIndexList, flattenedData):
    unflattened = [];
    startPoint=0;
    for sublist in pointLineIndexList:
        endPoint = startPoint + sublist[-1]+1;
        unflattened.append(flattenedData[startPoint : endPoint]);

        startPoint = endPoint;

    return unflattened;


#Calculates cell specific data (normals, coordinates, distances etc.)
#   inputDataPath - path to the parameter set's input data path (e.g. data/<paramsetName>)
#params = ps.get_current_params();
#inputDataPath = path.join("data", params.paramsetName);
#testPlots=True;
#verbose = False;
#outputPath = inputDataPath;
def generate_cell_data(params, inputDataPath, testPlots=False, verbose=False, outputPath=None):
#    params = ps.get_current_params();
#    ilatRange = params.ilatRange; #(50,225);
#    ilonRange = params.ilonRange; #(600, 850)
#    threshold=params.shallowDepth;
#    deepThreshold = params.deepDepth;
#    pixelRes = params.pixelRes; #spatial resolution (in degrees) of the current data
#    numLines = params.numberOfBoundryApproximationLines; #Number of lines to use to approximate the shelf edge

    #Load bathymetry - used to calculate gradients of normals
    bathy = Dataset(inputDataPath+"/GEBCO_bathymetry_"+str(params.pixelRes[1])+"x"+str(params.pixelRes[1])+"deg.nc", 'r');
    lat = np.flipud(bathy.variables["lat"][:]);
    lon = bathy.variables["lon"][:];
    depth = np.flipud(bathy.variables["mean_depth"][:]);
    depth = su.apply_mask(depth, None, params.ilatRange, params.ilonRange);

    #Read shelf-edge data
    shelfCoordData = pickle.load(open(path.join(inputDataPath, "shelf_coordinates", params.contourPathFile), "rb"));
    shelfCoordDataDeep = pickle.load(open(path.join(inputDataPath, "shelf_coordinates", params.contourPathFileDeep), "rb"));
    numPaths = len(shelfCoordData.coordinatesLists);
    numLines = shelfCoordData.numLines;

    #These will be used to create the cell data dataframe
    columns = ["x", "y", "lon", "lat", "distance", "distanceProp", "onshelfX", "onshelfY", "alongshelfX", "alongshelfY"];
    dtypes = [np.int16, np.int16, np.float32, np.float32, np.float32, np.float32, np.float32, np.float32, np.float32, np.float32];
    df = su.create_empty_df(columns, dtypes);

    #Calculate straight line segments for each shelf path
    lineXCoordsList = []; lineYCoordsList = []; lineParamsList = []; pointLineIndexList = [];
    for pathIndex in range(numPaths):
        if verbose:
            print "Generating straight line approximations for path #"+str(pathIndex+1)+" of", numPaths;

        shelfEdgeCoords = np.array(shelfCoordData.coordinatesLists[pathIndex], dtype=float);

        #approximate the shelf edge by generating a series of straight line segments along the shelf edge coordinates.
        lineXCoords, lineYCoords, lineParams, pointLineIndex = su.split_into_n_lines(shelfEdgeCoords[:,0], shelfEdgeCoords[:,1], numLines[pathIndex]);
        lineXCoordsList.append(lineXCoords);
        lineYCoordsList.append(lineYCoords);
        lineParamsList.append(lineParams);
        pointLineIndexList.append(pointLineIndex);

    #Remove so these don't accidentally get used.
    del (lineXCoords, lineYCoords, lineParams, pointLineIndex, shelfEdgeCoords);



    #combined all paths into one for each shallow and deep shelf boundary.
    #This is required because deep and shallow path data do not necessarily line up. Otherwise paths must be manually selected.
    print('Flatten')
    flatShallowCoordData = _flatten_shelf_coordinate_data(shelfCoordData, lineXCoordsList, lineYCoordsList, lineParamsList, pointLineIndexList);
    flatDeepCoordData = _flatten_shelf_coordinate_data(shelfCoordDataDeep, [], [], [], []); #No straight line approximation data for deep contour


    #Now calculate the shelf normals for the flattened data
    #Calculate the onto-shelf direction vectors for each line segment of the shelf edge approximation
    onShelfDirectionVectors, lineCentrePoints = su.calculate_on_shelf_direction(flatShallowCoordData.lineXCoords, flatShallowCoordData.lineYCoords, flatShallowCoordData.coordinates, flatDeepCoordData.coordinates);


    #Test plot to show shelf segments, the straight-line approximations of these segments, as well as their normal direction.
    if testPlots == True:
        normalsFig = plt.figure(figsize=(14,7))
        gs = GridSpec(1,1, figure=normalsFig, wspace=0.33,hspace=0.10,bottom=0.07,top=0.97,left=0.05,right=1)
        normalsAx = normalsFig.add_subplot(gs[0,0])
        a = normalsAx.imshow(depth,vmin=0,vmax=7000);

        for i in range(1, len(flatShallowCoordData.lineXCoords)):
            minCol = 0.25; scal=1.0-minCol; colour = (minCol+scal*random(), minCol+scal*random(), minCol+scal*random());
            segmentMask = np.array(flatShallowCoordData.pointLineIndex)==i;
            normalsAx.plot(flatShallowCoordData.coordinates[segmentMask,0], flatShallowCoordData.coordinates[segmentMask,1], color=(colour[0]*0.7, colour[1]*0.7, colour[2]*0.7), alpha=0.8, linewidth=2.0);
            normalsAx.plot(flatShallowCoordData.lineXCoords[i], flatShallowCoordData.lineYCoords[i], color = colour);
            normalsAx.arrow(lineCentrePoints[i][0], lineCentrePoints[i][1], onShelfDirectionVectors[i][0], onShelfDirectionVectors[i][1], head_width=6.0, head_length=6.0, fc=colour, ec=colour);
            #normalsAx.plot([lineCentrePoints[i][0], lineCentrePoints[i][0]+onShelfDirectionVectors[i][0]], [lineCentrePoints[i][1], lineCentrePoints[i][1]+onShelfDirectionVectors[i][1]], color=colour);
        #plt.title("Shelf segments, straight line approximations and onto-shelf directions");
        cbar = plt.colorbar(a)
        cbar.set_label("Depth (m)");
        normalsAx.set_xlabel('Longitude ($^{o}$E)')
        normalsAx.set_ylabel('Latitude ($^{o}$N)')
        xtickslocs = [0,180,360,540,720,900,1080,1260,1439]
        print(xtickslocs)
        labels = [str(np.around(lon[int(i)])) for i in xtickslocs]
        normalsAx.xaxis.set_ticks(xtickslocs)
        normalsAx.set_xticklabels(labels)
        ytickslocs = [0,180,360,540,719]
        labels = [str(np.around(lat[int(i)])) for i in ytickslocs]
        normalsAx.yaxis.set_ticks(ytickslocs)
        normalsAx.set_yticklabels(labels)
        normalsFig.savefig('plots/shelf_normals.png',dpi=300)
        plt.pause(1);


    #Unflatten onShelfDirectionVectors and lineCentrePoints to each shelf segment/coordinate pair
    #This is required so that adjacent coordinates belonging to different paths are not interpretted as continous shelf sections.
    #This also allows them to be used with the path-structured shelfCoordData data
    print('Unflattening')
    onShelfDirectionVectors = _unflatten(pointLineIndexList, onShelfDirectionVectors);
    lineCentrePoints = _unflatten(pointLineIndexList, lineCentrePoints);







    #Main calculations, loop through each path
    for pathIndex in range(numPaths):
        shelfEdgeCoords = np.array(shelfCoordData.coordinatesLists[pathIndex], dtype=float);

        #calculate where each of the shelf edge lines cross pixels.
        gridIntercepts, abnormalIntercepts = su.get_grid_cell_intercepts_with_edge_lines(shelfEdgeCoords, params);
        print(len(gridIntercepts))
        print(len(abnormalIntercepts))
        print(abnormalIntercepts)
        #Sanity check.
        # if len(abnormalIntercepts) != 0:
        #     raise RuntimeError("There were more than zero abnormal grid intercepts.");

        #Calculate distance proportions between each grid intercept point
        for intercept in gridIntercepts:
            #lon lat of diagonal distance
            diagonalDistance = geopy.distance.distance((intercept.b1lat, intercept.b1lon), (intercept.b2lat, intercept.b2lon)).meters;

            dist = geopy.distance.distance((intercept.i1lat, intercept.i1lon), (intercept.i2lat, intercept.i2lon)).meters;
            #dist = su.calculate_distance(intercept.ix1, intercept.iy1, intercept.ix2, intercept.iy2);
            distProportion = dist/diagonalDistance;
            intercept.dist = dist;
            intercept.distProportion = distProportion;

            #check by plotting grid intercepts
    #        fig, ax = plt.subplots(1);
    #        for intercept in gridIntercepts:
    #            su.plot_debug_intersect(intercept, ax);

        #bring together the shelf edge line mapping to approximation lines, onshelf direction vectors and grid intercept data
        #to construct a dataframe containing all this information in a per-pixel manner
        for iCell in range(0, len(gridIntercepts)):
            x = gridIntercepts[iCell].x;
            y = gridIntercepts[iCell].y;
            lon = gridIntercepts[iCell].gridCentreLon;
            lat = gridIntercepts[iCell].gridCentreLat;
            distance = gridIntercepts[iCell].dist;
            distanceProp = gridIntercepts[iCell].distProportion;


            correspondingApproxLine = pointLineIndexList[pathIndex][gridIntercepts[iCell].edgeShelfLineIndex]; #shelf edge line used to index into pointLineIndex to get the approx line number that corresponds to the shelf edge line.
            normalisingCoefficient = np.linalg.norm(onShelfDirectionVectors[pathIndex][correspondingApproxLine]);
            onshelfX = onShelfDirectionVectors[pathIndex][correspondingApproxLine][0]/normalisingCoefficient;
            onshelfY = onShelfDirectionVectors[pathIndex][correspondingApproxLine][1]/normalisingCoefficient;

            alongShelfVect = np.array([gridIntercepts[iCell].alongDirX, gridIntercepts[iCell].alongDirY]);
            alongShelfVect = alongShelfVect / np.linalg.norm(alongShelfVect);

            #Add row to dataframe
            row = (x, y, lon, lat, distance, distanceProp, onshelfX, onshelfY, alongShelfVect[0], alongShelfVect[1]);
            ###columns = ["x", "y", "lon", "lat", "distance", "distanceProp", "onshelfX", "onshelfY", "alongshelfX", "alongshelfY"];

            ser = pd.Series(row, index=columns);
            df=df.append(ser, ignore_index=True);
        print "Grid cells processed: ", len(df);






    if outputPath != None:
        print "Writing gridwise data."
        df.to_csv(path.join(outputPath, "gridwise_data", "per_grid_cell_edge_data_"+params.paramsetName+".csv"));

        #write other useful plotting info to file. This is used to produce many of the plots
        # (e.g. grid intercepts for the method demo plot)
        extraInfo = {};
        extraInfo["lineXCoords"] = lineXCoordsList;
        extraInfo["lineYCoords"] = lineYCoordsList;
        extraInfo["lineParams"] = lineParamsList;
        extraInfo["pointLineIndex"] = pointLineIndexList;
        extraInfo["onShelfDirectionVectors"] = onShelfDirectionVectors;
        extraInfo["lineCentrePoints"] = lineCentrePoints;
        extraInfo["gridIntercepts"] = gridIntercepts;

        pickle.dump(extraInfo, open(path.join(outputPath, "gridwise_data", "extra_info_"+params.paramsetName+".p"), "wb"));
