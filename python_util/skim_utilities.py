#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 24 13:28:59 2018

Misc utilities

@author: tom holding
Modified by Daniel J. Ford (d.ford@exeter.ac.uk)
Date: 03/2023
Changes:
- Updated surface current -> ekamn transport angle from 45 to 59.25 described in Rio and Picot (2014; https://doi.org/10.1002/2014GL061773)
"""

import numpy as np;
import matplotlib.pyplot as plt;
import inspect;
import matplotlib.patches as patches; #Rectangle
import math;


class QuickStruct:
    def __init__(self):
        pass;
    def to_multiline_string(self):
        s = "";
        attrNames = [name for name in dir(self) if not name.startswith('__')];
        for name in attrNames:
            val = getattr(self, name);
            if inspect.ismethod(val) == False:
                s += name + ":\t" + str(val) + "\n";
        return s;
    def to_string(self):
        s = "";
        attrNames = [name for name in dir(self) if not name.startswith('__')];
        for name in attrNames:
            val = getattr(self, name);
            if inspect.ismethod(val) == False:
                s += name + ":" + str(val) + ", ";
        return s;

    def __str__(self):
        return self.to_string();

    def value_string(self):
        s = "";
        attrNames = [name for name in dir(self) if not name.startswith('__')];
        for name in attrNames:
            val = getattr(self, name);
            if inspect.ismethod(val) == False:
                s+= str(val) + ", ";
        return s;

def _other_num_line_approximations(xCoords, yCoords):
    numCoordPairs = len(xCoords);
    if numCoordPairs != len(yCoords):
        raise ValueError("X and Y coordinate lists are not the same length.");
    return int(np.round(numCoordPairs/math.log(numCoordPairs, 1.35)));

def num_line_approximations(xCoords, yCoords):
    numCoordPairs = len(xCoords);
    if numCoordPairs != len(yCoords):
        raise ValueError("X and Y coordinate lists are not the same length.");
    return int(np.round(numCoordPairs*0.05))+1;

def fixed_length_lines(xCoords,yCoords):
    return int(30)

#Magnitude from x/y component of 2D vectors
def calculate_magnitude(xComponent, yComponent):
    return np.sqrt(xComponent**2 + yComponent**2);

#returns vector direction (clockwise from north) in degrees
def calculate_direction(xComponent, yComponent):
    return np.rad2deg(np.arctan2(xComponent, yComponent)); #clockwise from north

def calculate_distance(x1, y1, x2, y2):
    return np.sqrt( (x2-x1)**2 + (y2-y1)**2 );

#Calculate line parameters from the vector; dx, dy, and a point the line crosses; x1, y1
def calculate_line_parameters(dx, dy, x1, y1):
    m = dy/dx;
    c = y1 - m*x1;
    return m, c;

#Calculates shortest distance of point (xp, yp) from line (x1, y1) to (x2, y2)
def distance_from_line(xp, yp, x1, y1, x2, y2):
    from numpy.linalg import norm;

    point = np.array([xp, yp]);
    linePoint1 = np.array([x1, y1]);
    linePoint2 = np.array([x2, y2]);
    distance = np.cross(linePoint2-linePoint1, point-linePoint1)/norm(linePoint2-linePoint1);
    return np.abs(distance);

#Returns the two normal vectors of a 2D line (left hand, right hand)
def calculate_normal_vectors(x1, y1, x2, y2):
    dx = x2-x1;
    dy = y2-y1;
    return [(-dy, dx), (dy, -dx)];


def rotate_vector(vector, angleDeg):
    angleRad = np.deg2rad(angleDeg);
    newX = np.cos(angleRad)*vector[0] - np.sin(angleRad)*vector[1];
    newY = np.sin(angleRad)*vector[0] + np.cos(angleRad)*vector[1];
    return np.array([newX, newY]);

#returns the intersect point of two straight lines
def intersect_point(params1, params2):
    xi = (params2[1]-params1[1]) / (params1[0]-params2[0]);
    yi = params1[0]*xi + params1[1];
    return xi, yi;

#returns the coordinates at which point a line crosses a verticle line
def coords_at_const_x(lineParams, constX):
    return float(constX), lineParams[0]*float(constX)+lineParams[1];

#returns the coordinates at which point a line crosses a verticle line
def coords_at_const_y(lineParams, constY):
    return (float(constY)-lineParams[1])/lineParams[0], float(constY);

#returns whether a line intersects a box and the two sets of coordinates where it does as a tuple (intersected, p1, p2)
#box defined by lower left corner (bx1, by1) and upper right corner (bx2, by2).
#line defined by it's end points (lx1, ly1) and (lx2, ly2).
#assumes 'line' isn't horizontal or verticle
def line_intersects_box(bx1, by1, bx2, by2, lx1, ly1, lx2, ly2):
    def to_quick_struct(interceptCoords, label):
        out = QuickStruct();
        out.ix=interceptCoords[0];
        out.iy=interceptCoords[1];
        out.label = label;
        return out;

    if lx1>lx2:
        lx1, lx2 = lx2, lx1;
        ly1, ly2 = ly2, ly1;
    ldx = lx2-lx1;
    ldy = ly2-ly1;
    lineParams = calculate_line_parameters(ldx, ldy, lx1, ly1);

    #There must be between 0 and 2 intercepts. Find them
    intercepts = [];
    leftIntercept = coords_at_const_x(lineParams, bx1);
    rightIntercept = coords_at_const_x(lineParams, bx2);
    bottomIntercept = coords_at_const_y(lineParams, by1);
    topIntercept = coords_at_const_y(lineParams, by2);

    #Note that horizontal and verticle intercepts must use different open/closed bounds to permit perfect diagonals
    if (leftIntercept[1] > by1 and leftIntercept[1] <= by2):
        intercepts.append(to_quick_struct(leftIntercept, "W"));
    if (rightIntercept[1] > by1 and rightIntercept[1] <= by2):
        intercepts.append(to_quick_struct(rightIntercept, "E"));
    if (bottomIntercept[0] >= bx1 and bottomIntercept[0] < bx2):
        intercepts.append(to_quick_struct(bottomIntercept, "S"));
    if (topIntercept[0] >= bx1 and topIntercept[0] < bx2):
        intercepts.append(to_quick_struct(topIntercept, "N"));

    return intercepts;

#returns a QuickStruct containing information about whether a bounded line intersects an unbounded line:
#   intersectInfo.isIntersected (did an intersect occur within the bounds of the bounded line?)
#   intersectInfo.xi and .yi the coordinates of the intersect (regardless of whether it was in the bounded area or not)
def is_intersected(boundedX1, boundedY1, boundedX2, boundedY2, intersecterParams):
    dx = boundedX2 - boundedX1;
    dy = boundedY2 - boundedY1;
    mBounded, cBounded = calculate_line_parameters(dx, dy, boundedX1, boundedY1);

    #Calculate intersect point
    xi, yi = intersect_point((intersecterParams[0], intersecterParams[1]), (mBounded, cBounded));
    intersectInfo = QuickStruct();
    intersectInfo.xi = xi;
    intersectInfo.yi = yi;

    #Check to see if intersect point is within bounded line range
    if (xi >= np.min((boundedX1, boundedX2))) and (xi <= np.max((boundedX1, boundedX2))):
        intersectInfo.intersected = True;
    else:
        intersectInfo.intersected = False;

    return intersectInfo;

#Returns a list of QuickStruct objects which contain information about the intersects between an intersector and a set of lines.
#lineCoords defines the set of bounded lines to be intersected and takes the form of a 2D (n by 2) matrix of x and y coordinates.
#intersector defines the unbounded intersecting line in the form of a tuple containing the parameters for a straight line (m, c)
#distanceBaseline defines the point from which intersect distances will be calculated from.
#label gives an optional label to each intersect.
def calculate_intersect_info(lineCoords, intersector, distanceBaseline, label=None):
    infoList = [];
    for i in range(1, len(lineCoords)):
        intersectInfo = is_intersected(lineCoords[i-1,0], lineCoords[i-1,1], lineCoords[i,0], lineCoords[i,1], intersector);

        #If the lines intersect then calculate the distance from the intersect point to the centre of the approximate edge line
        if intersectInfo.intersected == True:
            #Add the intersected line (useful for debugging)
            intersectInfo.intersectedLineX = (lineCoords[i-1,0], lineCoords[i,0]);
            intersectInfo.intersectedLineY = (lineCoords[i-1,1], lineCoords[i,1]);

            #Calculate distance from xi, yi to centrePointx, centrePointY
            intersectInfo.dist = calculate_distance(intersectInfo.xi, distanceBaseline[0], intersectInfo.yi, distanceBaseline[1]);
            intersectInfo.offsetVector = np.array([intersectInfo.xi - distanceBaseline[0], intersectInfo.yi - distanceBaseline[1]]);

            #Add the intersect info to a list.
            intersectInfo.type = label;
            infoList.append(intersectInfo);
    return infoList;

#returns a list of direction vectors which point to the on-shelf direction from a straight line segment which approximates the shelf edge
#lineXCoords and lineYCoords are lists which define the (X1, X2) and (Y1, Y2) coordinates of the shelf edge approximation lines, respectively.
#shelfEdgeCoords and shelfEdgeDeepCoords are x,y coordinates of the shelf edge contour (2d n by 2 matrix) for shallow and deep contours respectively.
def calculate_on_shelf_direction(lineXCoords, lineYCoords, shelfEdgeCoords, shelfEdgeDeepCoords):
    #for each line find the inwardFlux direction (onto shelf)
    onShelfDirectionVectors = []; #Normal vectors in the 'onto shelf' direction
    lineCentrePoints = []; #Centre points of approx lines
    for lineNum in range(len(lineXCoords)):
        print "Line number", lineNum, "of", len(lineXCoords);
        #Get normals for line
        lineNormals = calculate_normal_vectors(lineXCoords[lineNum][0], lineYCoords[lineNum][0], lineXCoords[lineNum][1], lineYCoords[lineNum][1]);

        #find centre of approximation line
        approxLineCentrePointX = (lineXCoords[lineNum][1] - lineXCoords[lineNum][0])/2.0 + lineXCoords[lineNum][0];
        approxLineCentrePointY = (lineYCoords[lineNum][1] - lineYCoords[lineNum][0])/2.0 + lineYCoords[lineNum][0];

        #calculate perpendicular 'normal line' that crosses the centre of the line
        mNormal, cNormal = calculate_line_parameters(lineNormals[0][0], lineNormals[0][1], approxLineCentrePointX, approxLineCentrePointY);

        #find any shelf edge lines that crosses the normal line and calculate distance 'centre' point (the intersect of normal and approximation lines).
        intersectInfoList1 = calculate_intersect_info(shelfEdgeCoords, (mNormal, cNormal), (approxLineCentrePointX, approxLineCentrePointY), "shallow");
        #Repeat for deep shelf edge coordinate lines
        intersectInfoList2 = calculate_intersect_info(shelfEdgeDeepCoords, (mNormal, cNormal), (approxLineCentrePointX, approxLineCentrePointY), "deep");

        #Combine lists into one
        intersectInfoList = intersectInfoList1 + intersectInfoList2;

        #Assign positive or negative direction to distances from centre of the approximation line (positive = left hand, negative = right hand)
        for iInfo, intersectInfo in enumerate(intersectInfoList):
            offsetBearing = calculate_direction(intersectInfo.offsetVector[0], intersectInfo.offsetVector[1]); #from north (clockwise)
            normalBearing = calculate_direction(lineNormals[0][0], lineNormals[0][1]);
            reverseNormalBearing = calculate_direction(lineNormals[1][0], lineNormals[1][1]);

            #negate distance of offset bearing matches reverseNormalBearing.
            if np.isclose(offsetBearing, reverseNormalBearing, rtol=0, atol=1e-06):
                intersectInfoList[iInfo].dist = -intersectInfoList[iInfo].dist;
            elif np.isclose(offsetBearing, normalBearing, rtol=0, atol=1e-06) == False:
                raise ValueError("Could not match intersect offset bearing with either normal bearing!");

        #Using intersectInfo distances, identify which normal to use for the current approximation line
        #Find closest deep and closest shallow intercepts
        closeShallow = None;
        for intersectInfo in intersectInfoList:
            if (closeShallow == None) or (np.abs(intersectInfo.dist) < np.abs(closeShallow.dist)):
                closeShallow = intersectInfo;
        closeDeep = None;
        for intersectInfo in intersectInfoList:
            if (closeDeep == None) or (np.abs(intersectInfo.dist) < np.abs(closeDeep.dist)):
                closeDeep = intersectInfo;

        #Determine on shelf direction
        #if closeShallow == None or closeDeep == None:
        #    onShelfDirectionVectors.append(None);
       # else:
        if (closeDeep.dist < closeShallow.dist):
            onShelfDirectionVectors.append(lineNormals[0]);
        else:
            onShelfDirectionVectors.append(lineNormals[1]);
        lineCentrePoints.append( (approxLineCentrePointX, approxLineCentrePointY) );

    return onShelfDirectionVectors, lineCentrePoints;

#returns a list of direction vectors which point to the on-shelf direction from a straight line segment which approximates the shelf edge
#lineXCoords and lineYCoords are lists which define the (X1, X2) and (Y1, Y2) coordinates of the shelf edge approximation lines, respectively.
#shelfEdgeCoords and shelfEdgeDeepCoords are x,y coordinates of the shelf edge contour (2d n by 2 matrix) for shallow and deep contours respectively.
#pointLineIndex is the index of the straight line approximation that each point along the shelf belongs to
#bathymetry is a bathymetry on the same grid scale as the coordinates
def calculate_on_shelf_direction2(lineXCoords, lineYCoords, shelfEdgeCoords, pointLineIndex, bathymetry):
    #for each shelf coordinate find the onto-shelf direction
    onShelfDirectionVectors = []; #Normal vectors in the 'onto shelf' direction
    lineCentrePoints = []; #Centre points of approx lines

    for lineNum in range(len(lineXCoords)):


        #Calculate both direction vectors corresponding to the straigh line approximation's normal
        lineNormals = calculate_normal_vectors(lineXCoords[lineNum][0], lineYCoords[lineNum][0], lineXCoords[lineNum][1], lineYCoords[lineNum][1]);



        #Get the shelf point which belong to the current straight line approximation
        curEdgeCoords = shelfEdgeCoords[poinLineIndex==lineNum,:];

        #for each shelf coordinate, calculate which normal direction is 'onto shelf' (i.e. positive bathymetry gradient)
        count0=0; count1=0;
        for icoord in range(len(curEdgeCoords)):
            pass;
            #bathymetry at shelf coordinate

            #offset in direction 0

            #offset in direction 1

            #increment counts

        #Select the normal which is onto-shelf for most shelf points that belong to this line

        #calculate centre points of approximation line




    for ishelf in range(len(pointLineIndex)):
        lineNum = pointLineIndex[ishelf];
        approxLineCentrePointX = (lineXCoords[lineNum][1] - lineXCoords[lineNum][0])/2.0 + lineXCoords[lineNum][0];
        approxLineCentrePointY = (lineYCoords[lineNum][1] - lineYCoords[lineNum][0])/2.0 + lineYCoords[lineNum][0];

        #Calculate normals to the approximation line.
        lineNormals = calculate_normal_vectors(lineXCoords[lineNum][0], lineYCoords[lineNum][0], lineXCoords[lineNum][1], lineYCoords[lineNum][1]);
        #Convert to line parameters. Note that we can use either normal to get the gradient and intercept of this line
        mNormal, cNormal = calculate_line_parameters(lineNormals[0][0], lineNormals[0][1], approxLineCentrePointX, approxLineCentrePointY);

        #Using the gradient of the normal line, sample either side of the shelf coordinate
        shelfX, shelfY = shelfEdgeCoords[ishelf, :];


        #Calculate centre of shelf segment


        #Calculate position one cell from each shelf coord in the direction of each normal vector

        #assess which one has positive slow

        #select positive slow

        #calculate and add lineCentrePoints (centre of approximation line)

        #....


        #Get normals for line
        lineNormals = calculate_normal_vectors(lineXCoords[lineNum][0], lineYCoords[lineNum][0], lineXCoords[lineNum][1], lineYCoords[lineNum][1]);

        #find centre of approximation line
        approxLineCentrePointX = (lineXCoords[lineNum][1] - lineXCoords[lineNum][0])/2.0 + lineXCoords[lineNum][0];
        approxLineCentrePointY = (lineYCoords[lineNum][1] - lineYCoords[lineNum][0])/2.0 + lineYCoords[lineNum][0];

        #calculate perpendicular 'normal line' that crosses the centre of the current straight line approximation
        #note that we can use either normal to get the gradient and intercept of this line
        mNormal, cNormal = calculate_line_parameters(lineNormals[0][0], lineNormals[0][1], approxLineCentrePointX, approxLineCentrePointY);

        #find any shelf edge lines that crosses the normal line and calculate distance 'centre' point (the intersect of normal and approximation lines).
        intersectInfoList1 = calculate_intersect_info(shelfEdgeCoords, (mNormal, cNormal), (approxLineCentrePointX, approxLineCentrePointY), "shallow");
        #Repeat for deep shelf edge coordinate lines
        intersectInfoList2 = calculate_intersect_info(shelfEdgeDeepCoords, (mNormal, cNormal), (approxLineCentrePointX, approxLineCentrePointY), "deep");

        #Combine lists into one
        intersectInfoList = intersectInfoList1 + intersectInfoList2;

        #Assign positive or negative direction to distances from centre of the approximation line (positive = left hand, negative = right hand)
        for iInfo, intersectInfo in enumerate(intersectInfoList):
            offsetBearing = calculate_direction(intersectInfo.offsetVector[0], intersectInfo.offsetVector[1]); #from north (clockwise)
            normalBearing = calculate_direction(lineNormals[0][0], lineNormals[0][1]);
            reverseNormalBearing = calculate_direction(lineNormals[1][0], lineNormals[1][1]);

            #negate distance of offset bearing matches reverseNormalBearing.
            if np.isclose(offsetBearing, reverseNormalBearing, rtol=0, atol=1e-06):
                intersectInfoList[iInfo].dist = -intersectInfoList[iInfo].dist;
            elif np.isclose(offsetBearing, normalBearing, rtol=0, atol=1e-06) == False:
                raise ValueError("Could not match intersect offset bearing with either normal bearing!");

        #Using intersectInfo distances, identify which normal to use for the current approximation line
        #Find closest deep and closest shallow intercepts
        closeShallow = None;
        for intersectInfo in intersectInfoList:
            if (closeShallow == None) or (np.abs(intersectInfo.dist) < np.abs(closeShallow.dist)):
                closeShallow = intersectInfo;
        closeDeep = None;
        for intersectInfo in intersectInfoList:
            if (closeDeep == None) or (np.abs(intersectInfo.dist) < np.abs(closeDeep.dist)):
                closeDeep = intersectInfo;

        #Determine on shelf direction
        #if closeShallow == None or closeDeep == None:
        #    onShelfDirectionVectors.append(None);
       # else:
        if (closeDeep.dist < closeShallow.dist):
            onShelfDirectionVectors.append(lineNormals[0]);
        else:
            onShelfDirectionVectors.append(lineNormals[1]);
        lineCentrePoints.append( (approxLineCentrePointX, approxLineCentrePointY) );

    return onShelfDirectionVectors, lineCentrePoints;


def apply_mask(data, mask=None, ilatRange=None, ilonRange=None):
    d = data.copy();
    if mask is not None:
        d.mask = d.mask | (mask!=1);

    if ilatRange != None and ilonRange != None:
        d = d[ilatRange[0]:ilatRange[1], ilonRange[0]:ilonRange[1]];
    return d;

#Adds a striaght line to the active plot.
#TODO: Delete, no longer needed
def draw_line(xRange, params):
    yRange = [(params[0]*xRange[0] + params[1]), (params[0]*xRange[1] + params[1])];
    plt.plot(xRange, yRange, 'k');

#splits 'coordinates' into n contiguous sections and fits a line to each section
#Each line is constrained by the need to cross through the end of the previous line.
#returns a tuple containing:
#   list of x ranges that the parameters apply to as a typle: (start, end)
#   tuple containing fitted parameters (m, c) for y=mx+c
def fit_n_lines(x, y, n):
    from scipy.optimize import curve_fit;
    if len(x) != len(y):
        raise ValueError("x and y must have the same length.");

    #Split data into n subsets
    xList = np.array_split(x, n);
    yList = np.array_split(y, n);
    xRange = [(xList[0][0], xList[0][-1])];
    for i in range(1, len(xList)):
        xRange.append( (xList[i-1][-1], xList[i][-1]) );


    #Fit each subset of coordinates
    def f(x, m, c): #Straight line, y=mx+c
        return m*x + c;
    def f2(x, m): #Straight line clamped at known point, y=mx+c where c = z - m*x_0
        return m*x + (f2.z - m*f2.x0);
    parameters = [];
    for i in range(0, n):
        if i==0:
            m, c = curve_fit(f, xList[i], yList[i])[0];
            parameters.append((m, c));
            f2.x0 = xRange[i][1];
            f2.z = m*f2.x0+c;
        else:
            m = curve_fit(f2, xList[i], yList[i])[0];
            c = f2.z - m*f2.x0;
            parameters.append((m, c));
            f2.x0 = xRange[i][1];
            f2.z = m*f2.x0+c;
    return xRange, parameters;


#splits 'coordinates' into n contiguous sections by distance along the line,
#returns a tuple containing:
#   list of x ranges that the parameters apply to as a typle: (start, end)
#   tuple containing fitted parameters of the line (m, c) for y=mx+c
def split_into_n_lines(x, y, n):
    if len(x) != len(y):
        raise ValueError("x and y must have the same length.");

    def calc_dist(x1, y1, x2, y2):
        return np.sqrt( (x1-x2)**2 + (y1-y2)**2 );

    def interpolate_line(proportion, x1, y1, x2, y2):
        proportion = float(proportion);
        xi = proportion*x1 + (1.0-proportion)*x2;
        yi = proportion*y1 + (1.0-proportion)*y2;
        return xi, yi;

    #Calculate total line length
    length = 0.0;
    for i in range(1, len(x)):
        length += calc_dist(x[i-1], y[i-1], x[i], y[i]);

    #Calculate length of each line
    lineLengthThreshold = length/n;
    if lineLengthThreshold < 1.5: # > the grid cell diagonal (~1.4)
        raise ValueError("lineLength is too small (<2.0 grid cells). Halted because results would be unreliable.");

    parameters = [];
    xRange = [];
    yRange = [];
    #iterate through coordinates summing distances until we exceed lineLength. Calculate coordinates which correspond to the length
    #Calculate parameters for straight line between start and stop coordinates
    growingLength = 0;
    curX1 = x[0];
    curY1 = y[0];
    pointLineIndex = [0]; #index corresponding to the line that each point 'belongs' to
    lineNum=0;
    for i in range(1, len(x)):
        #Store coordinates for convenience
        x1 = x[i-1]; y1 = y[i-1];
        x2 = x[i]; y2 = y[i];

        newLength = calc_dist(x1, y1, x2, y2);


        if growingLength+newLength < lineLengthThreshold and i != (len(x)-1): #If there new line can be added without exceeding the length threshold then do so
                                                                            #AND it isn't the last line (last line always added)
            growingLength += newLength;
            pointLineIndex.append(lineNum);
            continue;

        else: #if adding growingLength would supassed the length threshold
            #calculate coordinates along the current line that will give the correct length
            requiredDistance = lineLengthThreshold-float(growingLength); #How far along the current line do we need to go in order to reach the length threshold?
            #Perform linear interpolation to get the correct point on the line
            proportion = requiredDistance/float(newLength);
            xi, yi = interpolate_line(proportion, x1, y1, x2, y2);

            #The interpolated points are the end point of the current growing line
            curX2 = xi;
            curY2 = yi;

            #Now we know start end end of the line, add the ranges and parameters
            xRange.append( (curX1, curX2) );
            yRange.append( (curY1, curY2) );
            m = float(curY2-curY1)/float(curX2-curX1);
            c = float(curY2 - m*curX2);
            parameters.append( (m, c) );

            #Calculate values for the next growing line
            remainingDistance = float(newLength) - requiredDistance;
            growingLength = remainingDistance;
            curX1 = curX2; #end point of last line is start of the current line.
            curY1 = curY2;

            if (proportion < 0.5) or (i == len(x)-1): #check which line the point should belong to. #Note last line always belongs to lineNum not lineNum+1
                pointLineIndex.append(lineNum);
            else:
                pointLineIndex.append(lineNum+1);
            lineNum += 1;

    return xRange, yRange, parameters, np.array(pointLineIndex);


def coord_values(start, stop, resolution):
    return np.linspace(start, stop-resolution, int((stop-start)/resolution))

def coord_values_from_idx(start, stop, resolution):
    mask = np.arange(start, stop);
    fullCoords = coord_values(0, 360, resolution);
    return fullCoords[mask];

#returns the indices corresponding to a given longitude and latitude.
#resolution is the grid resolution (lonres, latres)
#ilat0 and ilon0 are optional offsets (for non-global grids) specifying the
#   index corresponding to the first valid point on global grid.
def lat_lon_to_coords(lat, lon, resolution, ilat0=0, ilon0=0):
    if (lat % resolution[0] != 0 or lon % resolution[1] != 0):
        raise ValueError("provided lon/lat are not perfect multiples of the supplied resolution.");

    ilat = (lat/resolution[0]) - ilat0;
    ilon = (lon/resolution[1]) - ilon0;
    return (ilat, ilon);

#Find where shelf edge lines intersect grid cells.
#Returns a list of QuickStruct objects which contains information on each intersect found:
#the grid cell coordinates in which an intersect was found
#the grid cell bounding box
#edge line which intersects the grid cell
#the number of other intersects which occured for this grid cell
#coordinates of the intersect
#a label for the intersect indicating which edge of the bounding box was intersected (e.g. N, S, E, W)
def get_grid_cell_intercepts_with_edge_lines(shelfEdgeCoords, params):
    pixelRes = params.pixelRes;
    fullIntercepts = [];
    abnormalIntercepts = [];
    for iLine in range(1, len(shelfEdgeCoords)):
        lx1 = shelfEdgeCoords[iLine-1][0];
        ly1 = shelfEdgeCoords[iLine-1][1];
        lx2 = shelfEdgeCoords[iLine][0];
        ly2 = shelfEdgeCoords[iLine][1];

        #lowerLon = np.floor(min(lx1, lx2)/pixelRes[1]) * pixelRes[1];
        #upperLon = np.ceil(max(lx1, lx2)/pixelRes[1]) * pixelRes[1];
        #lowerLat = np.floor(min(ly1, ly2)/pixelRes[0]) * pixelRes[0];
        #upperLat = np.ceil(max(ly1, ly2)/pixelRes[0]) * pixelRes[0];

        lowerX = int(np.floor(min(lx1, lx2)));
        upperX = int(np.ceil(max(lx1, lx2)));
        lowerY = int(np.floor(min(ly1, ly2)));
        upperY = int(np.ceil(max(ly1, ly2)));

        #for cellLat in coord_values(lowerLat, upperLat, pixelRes[0]):
        #    for cellLon in coord_values(lowerLon, upperLon, pixelRes[1]):
        for y in range(lowerY, upperY):
            for x in range(lowerX, upperX):
                bx1 = x;
                bx2 = x+1;
                by1 = y;
                by2 = y+1;

                intercepts = line_intersects_box(bx1, by1, bx2, by2, lx1, ly1, lx2, ly2);
                if len(intercepts) == 2:
                    struct = QuickStruct();
                    struct.ix1 = intercepts[0].ix;
                    struct.ix2 = intercepts[1].ix;
                    struct.iy1 = intercepts[0].iy;
                    struct.iy2 = intercepts[1].iy;

                    struct.alongDirX = struct.ix2-struct.ix1;
                    struct.alongDirY = struct.iy2-struct.iy1;

                    struct.label1 = intercepts[0].label;
                    struct.label2 = intercepts[1].label;
                    struct.bx1 = bx1;
                    struct.bx2 = bx2;
                    struct.by1 = by1;
                    struct.by2 = by2;
                    struct.lx1 = lx1;
                    struct.lx2 = lx2;
                    struct.ly1 = ly1;
                    struct.ly2 = ly2;

                    struct.x = x; #Cell coordinates
                    struct.y = y;
                    #Calculates lon/lat of centre of cell
                    struct.gridCentreLon, struct.gridCentreLat = convert_index_to_lonlat(x, y, pixelRes, lon0=params.originLon, lat0=params.originLat, centre=True);
                    struct.b1lon, struct.b1lat = convert_index_to_lonlat(x, y, pixelRes, lon0=params.originLon, lat0=params.originLat, centre=False);
                    struct.b2lon = struct.b1lon + pixelRes[1];
                    struct.b2lat = struct.b1lat + pixelRes[0];
                    struct.i1lon, struct.i1lat = convert_index_to_lonlat(struct.ix1, struct.iy1, pixelRes, lon0=params.originLon, lat0=params.originLat);
                    struct.i2lon, struct.i2lat = convert_index_to_lonlat(struct.ix2, struct.iy2, pixelRes, lon0=params.originLon, lat0=params.originLat);


                    struct.numIntersectsInBox = len(intercepts); #used later for differentiating between full intercepts, line ending in box, and glancing passes.
                    struct.edgeShelfLineIndex = iLine;
                    fullIntercepts.append(struct);

                else: #greater than or less than two intercets
                    for intercept in intercepts:
                        intercept.bx1 = bx1;
                        intercept.bx2 = bx2;
                        intercept.by1 = by1;
                        intercept.by2 = by2;
                        intercept.lx1 = lx1;
                        intercept.lx2 = lx2;
                        intercept.ly1 = ly1;
                        intercept.ly2 = ly2;
                        lon, lat = convert_index_to_lonlat(x, y, pixelRes, lon0=params.originLon, lat0=params.originLat);
                        struct.gridCentreLat = lat+(pixelRes[0]/2.0);
                        struct.gridCentreLon = lon+(pixelRes[1]/2.0);
                        struct.edgeShelfLineIndex = iLine;
                        intercept.numIntersectsInBox = len(intercepts); #used later for differentiating between full intercepts, line ending in box, and glancing passes.
                        abnormalIntercepts.append(intercept);
    return fullIntercepts, abnormalIntercepts;

#creates a new set of lines from an original set of lines by expanding in the direction indicated.
def extract_parallel_transect(originalLineParams, originalLineXCoords, originalLineYCoords, originalLineCentrePoints, directionList, distance, reverseDirection=True):
    #create list of parallel lines to the shelf edge contour
    distance = 2.0;
    reverseDirection = True;

    parallelLineParams = [];
    parallelLineCentrePoints = [];
    for i in range(len(originalLineParams)):
        direction = directionList[i] / np.linalg.norm(directionList[i]);
        if reverseDirection:
            direction = -direction;
        newX, newY = originalLineCentrePoints[i]+(direction*distance);

        m = originalLineParams[i][0];
        c = newY - m*newX;
        parallelLineParams.append((m, c));
        parallelLineCentrePoints.append( (newX, newY) );

    #Find intercepts between each line and calculate the X and Y start and end points for each line.
    parallelEdgeLinesX = [];
    parallelEdgeLinesY = [];
    for i, parallelLineParam in enumerate(parallelLineParams):
        if i==0: #start of the first line is the start of the first edge line + the offset
            direction = directionList[0] / np.linalg.norm(directionList[0]);
            if reverseDirection:
                direction = -direction;
            x0 = originalLineXCoords[0][0]+(direction[0]*distance);
            y0 = originalLineYCoords[0][0]+(direction[1]*distance);
        else: #start of the current line is the end of the last parallel line
            x0 = parallelEdgeLinesX[i-1][1];
            y0 = parallelEdgeLinesY[i-1][1];

        if i==len(parallelLineParams)-1: #end of the last line is the end of the last edge line + the offset
            direction = directionList[i] / np.linalg.norm(directionList[i]);
            if reverseDirection:
                direction = -direction;
            x1 = originalLineXCoords[i][1]+(direction[0]*distance);
            y1 = originalLineYCoords[i][1]+(direction[1]*distance);
        else: #end of the current line is the intersect with the next line
            x1, y1 = intersect_point(parallelLineParams[i], parallelLineParams[i+1]);

        parallelEdgeLinesX.append( (x0, x1));
        parallelEdgeLinesY.append( (y0, y1));

    return parallelEdgeLinesX, parallelEdgeLinesY, parallelLineCentrePoints, parallelLineParams;

#combines paired X with paird Y coordinates into a single numpy array
def combind_XY_coords(xCoords, yCoords):
    out = np.zeros( (len(xCoords)+1, 2), dtype=float);
    for i in range(len(xCoords)):
        out[i,0] = xCoords[i][0];
        out[i,1] = yCoords[i][0];
    out[i+1, 0] = xCoords[i][1];
    out[i+1, 1] = yCoords[i][1];

    return out;


def plot_debug_intersect(intercept, ax, lineColour='k'):
    ax.plot([intercept.lx1, intercept.lx2], [intercept.ly1, intercept.ly2], lineColour);
    box = patches.Rectangle( (intercept.bx1, intercept.by1), intercept.bx2-intercept.bx1, intercept.by2-intercept.by1, edgecolor='r', facecolor=None);
    ax.add_patch(box);

#Create an empty dataframe with specific dtypes
#From https://stackoverflow.com/questions/36462257/create-empty-dataframe-in-pandas-specifying-column-types
def create_empty_df(columns, dtypes, index=None):
    import pandas as pd;
    assert len(columns)==len(dtypes)
    df = pd.DataFrame(index=index)
    for c,d in zip(columns, dtypes):
        df[c] = pd.Series(dtype=d)
    return df

##new and tested.
##y, x: the y and x index
##pixelRes: tupple containing (latitudeResolution, longitudeResolution)
##lon0, lat0: the longitude and latitude at index [0, 0]
##Returns the lonlat of the bottom left grid cell corner unless centre==True when centre is used
def convert_index_to_lonlat(x, y, pixelRes, lon0=-180, lat0=89.75, centre=False):
    #print(x)
    #print(pixelRes)
    lon = lon0 + (x*pixelRes[1])
    lat = lat0 - (y*pixelRes[0]);
    if centre == False:
        return lon, lat;
    else:
        return lon+(pixelRes[1]/2.0), lat+(pixelRes[0]/2.0);


#new and tested.
#lon, lat: the lon and lat index
#pixelRes: tupple containing (latitudeResolution, longitudeResolution)
#lon0, lat0: the longitude and latitude at index [0, 0]
#lonlat refers to the bottom left grid cell corner unless centre==True when centre is expected
def convert_lonlat_to_index(_lon, _lat, pixelRes, lon0=-180, lat0=89.75, centre=False):
    lon = _lon;
    lat = _lat;
    if centre:
        lon -= (pixelRes[1]/2.0);
        lat -= (pixelRes[0]/2.0);

    x = (lon-lon0) / pixelRes[1];
    y = -(lat-lat0) / pixelRes[0];
    return int(x), int(y);


#Returns the normalised direction vector for the direction of ekman current at the surface (90 degrees right/left from the windspeed).
def get_surface_ekman_direction_vector(windspeedVector, latDegrees, flippedY=True):
    if latDegrees == 0.0: #If at equator, no change... This is a bit of a simplification.
        return windspeedVector/np.linalg.norm(windspeedVector);

    normals = calculate_normal_vectors(0.0, 0.0, windspeedVector[0], windspeedVector[1]);
    if flippedY:
        if latDegrees > 0.0: #Northern hemisphere
            transportDirection = np.array(normals[0]);
        else: #Southern hemisphere, TODO: equator is special case.
            transportDirection = np.array(normals[1]);
    else:
        if latDegrees > 0.0: #Northern hemisphere
            transportDirection = np.array(normals[1]);
        else: #Southern hemisphere, TODO: equator is special case.
            transportDirection = np.array(normals[0]);
    return transportDirection / np.linalg.norm(transportDirection);

#Returns the normalised direction vector for the direction of ekman current at the surface (90 degrees right/left from the windspeed).
def get_integrated_ekman_direction(ekmanVector, latDegrees, flippedY=True):
    if latDegrees == 0.0: #If at equator, no change... This is a bit of a simplification.
        return ekmanVector/np.linalg.norm(ekmanVector);

    rotation = 59.25; #degrees. Net Ekman transport is at 45 degrees to that of the surface flow
    if flippedY:
        rotation = -rotation;

    if latDegrees > 0.0: #if northern hemisphere
        rotation = -rotation; #move to the right, otherwise keep rotation to the left

    transportDirection = rotate_vector(ekmanVector, rotation);
    return transportDirection/np.linalg.norm(transportDirection);
