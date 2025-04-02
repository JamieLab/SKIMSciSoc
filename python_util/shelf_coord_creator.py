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
import shapefile

#Generates contours for a particular depth on a bathymetry.
#Extracts lists of lists for each path in the contour containing coordinate data.
#Saves a pickled QuickStruct contraining this data and also returns the QuickStruct object (for convenience).
#bathymetry - 2D matrix containing bathymetry
#contourDepth - depth at which to draw contours
#minLengthThreshold - paths with fewer coordinate pairs than this will be ignored
#numLinesFunc - function that calculates the number of lines that should be used to approximate this part of the path
#lon, lat - the longitude and latitude values that define the grid.
def _extract_contour_paths(bathymetry, lon, lat, contourDepth, minLengthThreshold, numLinesFunc, verbose=False,shape=False,shape_file=''):
    #Generate contours
    plt.ioff();
    if shape:
        lat_temp = np.flipud(lat)
        # lon = lon.data
        sf = shapefile.Reader(shape_file)
        shapes = sf.shapes()
        points = shapes[0].points
        parts = shapes[0].parts
        # Need to convert points into lon lat indexes
        x_out = []
        y_out = []
        for j in range(len(parts)):
            print(j)
            x = []
            y = []
            #print(lon)
            if j < len(parts)-1:
                points_t = points[parts[j]:parts[j+1]]
            else:
                points_t = points[parts[j]:]
            for i in range(len(points_t)):
                xt = points_t[i][0]
                f = np.where(np.abs(xt-lon) == np.min(np.abs(xt-lon)))
                print(f)
                print(lon[f[0]])
                print(xt)
                print(xt-lon[f[0]])

                xv = f[0] + ((xt-lon[f[0]])/0.25)
                print(xv)
                x.append(xv.data[0])

                yt = points_t[i][1]
                f = np.where(np.abs(yt-lat_temp) == np.min(np.abs(yt-lat_temp)))
                yv = f[0] - ((yt-lat_temp[f[0]])/0.25)
                print(yv)
                y.append(yv.data[0])
            x_out.append(x)
            y_out.append(y)
        #print(x)
        #print(y)
        #plt.contour(range(0, len(lon)), range(0, len(lat)), bathymetry, [contourDepth], colors=[(1,0.5,0.5)]);
        for i in range(len(x_out)):
            plt.plot(x_out[i],y_out[i])
            plt.scatter(x_out[i],y_out[i])
        plt.show()

    else:
        cs = plt.contour(range(0, len(lon)), range(0, len(lat)), bathymetry, [contourDepth], colors=[(1,0.5,0.5)]);
    plt.close();
    plt.ion();

    #Filter paths by length, store x and y data for each path that passes.
    #Also calculate the number of straight line approximations needed for each path.
    xs = []; ys = [];
    pathIndicesUsed = [];
    numLines = [];
    if shape:
        # print(len(x_out))
        for i in range(len(x_out)):
            v = (x_out[i],y_out[i])
            print(v)
            # print(minLengthThreshold)
            # print(len(v))
            if len(v[0])>=minLengthThreshold:
                xs.append(v[0][:])
                ys.append(v[1][:])
                pathIndicesUsed = i
                numLines.append(numLinesFunc(xs[-1], ys[-1]))
    else:
        for pathIndex in range(len(cs.collections[0].get_paths())):
            #Extract current path and vertices
            p = cs.collections[0].get_paths()[pathIndex];
            v = p.vertices;
            if len(v) >= minLengthThreshold: #If the threshold is passed, add this path
                xs.append(v[:,0]);
                ys.append(v[:,1]);
                pathIndicesUsed.append(pathIndex);
                numLines.append(numLinesFunc(xs[-1], ys[-1]))
    print(xs)
    print(ys)
    #Zip x and y coordinates into a single matrix (useful format for some calculations)
    coordinatesLists = [];
    for i in range(1):#range(len(xs)):
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
def _test_plot_coordinate_data(coordData, depth, title,coordDatadeep=False):
    plt.figure();
    plt.imshow(depth);
    for i in range(len(coordData.xs)):
        plt.scatter(coordData.xs[i], coordData.ys[i], s=0.5, color=(1,0,0));
    if coordDatadeep:
        for i in range(len(coordDatadeep.xs)):
            plt.scatter(coordDatadeep.xs[i], coordDatadeep.ys[i], s=0.5, color=(0,0,1));
    plt.title(title);



#Stores the shallow and deep shelf edge coordinates in a file.
#Stores the number of straigh line approximations to use in a file
def get_shelf_edge_data(params, testPlots=False, outputPath=None,shape=False,shape_file_shallow='',shape_file_deep=''):
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
    shallowShelfCoordinateData = _extract_contour_paths(depth, lon, lat, params.shallowDepth, params.minContourPathSizeShallow, params.numLineApproximationFunction, verbose=True,shape=shape,shape_file=shape_file_shallow);
    deepShelfCoordinateData = _extract_contour_paths(depth, lon, lat, params.deepDepth, params.minContourPathSizeDeep, params.numLineApproximationFunction,shape=shape,shape_file=shape_file_deep);

    if testPlots:
        _test_plot_coordinate_data(shallowShelfCoordinateData, depth, "Shallow shelf used",coordDatadeep=deepShelfCoordinateData);
        _test_plot_coordinate_data(deepShelfCoordinateData, depth, "Deep shelf used");
        plt.pause(1);

    #Save coordinate data
    if outputPath != None:
        print "Writing shallow shelf coordinate data to: ", path.join(outputPath, "shelf_coordinates", params.contourPathFile);
        pickle.dump(shallowShelfCoordinateData, open(path.join(outputPath, "shelf_coordinates", params.contourPathFile), "wb"));
        print "Writing deep shelf coordinate data to: ", path.join(outputPath, "shelf_coordinates", params.contourPathFileDeep);
        pickle.dump(deepShelfCoordinateData, open(path.join(outputPath, "shelf_coordinates", params.contourPathFileDeep), "wb"));
