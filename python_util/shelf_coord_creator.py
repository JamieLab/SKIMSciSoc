#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 20 11:23:51 2018

@author: rr
"""

import skim_utilities as su;
import numpy as np;
from netCDF4 import Dataset;
import matplotlib.pyplot as plt;
from os import path;
import cPickle as pickle;

#Generates contours for a particular depth on a bathymetry.
#Extracts lists of lists for each path in the contour containing coordinate data.
#Saves a pickled QuickStruct contraining this data and also returns the QuickStruct object (for convenience).
#bathymetry - 2D matrix containing bathymetry
#contourDepth - depth at which to draw contours
#minLengthThreshold - paths with fewer coordinate pairs than this will be ignored
#numLinesFunc - function that calculates the number of lines that should be used to approximate this part of the path
#lon, lat - the longitude and latitude values that define the grid.
def _extract_contour_paths(bathymetry, lon, lat, contourDepth, minLengthThreshold, numLinesFunc, verbose=False):
    #Generate contours
    plt.ioff();
    cs = plt.contour(range(0, len(lon)), range(0, len(lat)), bathymetry, [contourDepth], colors=[(1,0.5,0.5)]);
    plt.close();
    plt.ion();

    #Filter paths by length, store x and y data for each path that passes.
    #Also calculate the number of straight line approximations needed for each path.
    xs = []; ys = [];
    pathIndicesUsed = [];
    numLines = [];
    for pathIndex in range(len(cs.collections[0].get_paths())):
        #Extract current path and vertices
        p = cs.collections[0].get_paths()[pathIndex];
        v = p.vertices;
        if len(v) >= minLengthThreshold: #If the threshold is passed, add this path
            xs.append(list(v[:,0]));
            ys.append(list(v[:,1]));
            pathIndicesUsed.append(pathIndex);
            numLines.append(numLinesFunc(xs[-1], ys[-1]))

    #Zip x and y coordinates into a single matrix (useful format for some calculations)
    coordinatesLists = [];
    for i in range(len(xs)):
        coordinatesLists.append (zip(xs[i],ys[i]));

    #Store data as a struct and return
    shelfCoordinateData = su.QuickStruct();
    shelfCoordinateData.xs = xs;
    shelfCoordinateData.ys = ys;
    shelfCoordinateData.coordinatesLists = coordinatesLists;
    shelfCoordinateData.pathIndicesUsed = pathIndicesUsed;
    shelfCoordinateData.contourDepth = contourDepth;
    shelfCoordinateData.numLines = numLines;

    if verbose:
        print "Total number of lines:", np.sum(numLines);

    return shelfCoordinateData;

#Plots
def _test_plot_coordinate_data(coordData, depth, title):
    plt.figure();
    plt.imshow(depth);
    for i in range(len(coordData.xs)):
        plt.scatter(coordData.xs[i], coordData.ys[i], s=0.5, color=(1,0,0));
    plt.title(title);



#Stores the shallow and deep shelf edge coordinates in a file.
#Stores the number of straigh line approximations to use in a file
def get_shelf_edge_data(params, testPlots=False, outputPath=None):
    ilatRange = params.ilatRange; #(50,225);
    print(ilatRange)
    ilonRange = params.ilonRange; #(600, 850)

    #Load depth dataset
    bathy = Dataset(outputPath+"/GEBCO_bathymetry_"+str(params.pixelRes[0])+"x"+str(params.pixelRes[1])+"deg.nc", 'r');
    depth = np.flipud(bathy.variables["mean_depth"][:]);
    depth = su.apply_mask(depth, None, ilatRange, ilonRange);
    lat = bathy.variables["lat"][ilatRange[0]:ilatRange[1]];
    lon = bathy.variables["lon"][ilonRange[0]:ilonRange[1]];
    print(depth.shape)
    print(lat.shape)
    print(lon.shape)

    #Get shallow shelf coordinates
    shallowShelfCoordinateData = _extract_contour_paths(depth, lon, lat, params.shallowDepth, params.minContourPathSizeShallow, params.numLineApproximationFunction, verbose=True);
    deepShelfCoordinateData = _extract_contour_paths(depth, lon, lat, params.deepDepth, params.minContourPathSizeDeep, params.numLineApproximationFunction);

    if testPlots:
        _test_plot_coordinate_data(shallowShelfCoordinateData, depth, "Shallow shelf used");
        _test_plot_coordinate_data(deepShelfCoordinateData, depth, "Deep shelf used");
        #plt.pause(1);

    #Save coordinate data
    if outputPath != None:
        print "Writing shallow shelf coordinate data to: ", path.join(outputPath, "shelf_coordinates", params.contourPathFile);
        pickle.dump(shallowShelfCoordinateData, open(path.join(outputPath, "shelf_coordinates", params.contourPathFile), "wb"));
        print "Writing deep shelf coordinate data to: ", path.join(outputPath, "shelf_coordinates", params.contourPathFileDeep);
        pickle.dump(deepShelfCoordinateData, open(path.join(outputPath, "shelf_coordinates", params.contourPathFileDeep), "wb"));
