#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 14:28:00 2018

@author: rr
"""

import skim_utilities as su;
import numpy as np;

#Our areas
area_EuropeanShelf = ((-15.0, 40.0), (20.0, 70.0)); #European shelf
area_LabradorSea = ((-70.0, 40.0), (-40.0, 70.0)); #Labrador Sea
area_MidAtlanticBight = ((-78.0, 32.0), (-67.0, 43.0)); #Mid-Atlantic Bight
area_CoastOfJapan = ((127.0, 29.0), (149.0, 44.0)); #Coast of Japan
area_Patagonia = ((-87.0, -58.0), (-59.0, -41.0)); #Patagonia
area_Tasmania = ((141.0, -46.0), (153.0, -38.0)); #Tasmania
area_BarentsSea = ((0.0, 60.0), (75.0, 90.0)); #Barents Sea
area_SouthAtlanticBight = ((-82.0, 23.5), (-74.0, 36.5)); #South Atlantic Bight
area_SouthGreenland = ((-57.0, 59.0), (-18.5, 67.5)); #South Greenland
area_AntarcticPeninsula = ((-80.0, -75.0), (-50.0, -61.0)) #Antarctic Peninsula
area_BeringSeaWest = ((-180.0, 45.0), (-150.0, 61.0)); #Berine Sea (West)
area_BeringSeaEast = ((178.0, 45.0), (180.0, 61.0)); #Bering Sea (East)
area_CascianShelf = ((-140.0, 40.0), (-120.0, 54.0)); #Cascian Shelf
#area_IrmingerSea = ((-57.0, 59.0), (-18.5, 67.5)); #Irminger Sea
area_IrmingerSea = ((-47.0, 58.0), (-17.5, 67.0)); #Irminger Sea

#Laruelle 2018 areas (approx.) (table 2)
#Laruelle, G. G., Cai, W. J., Hu, X., Gruber, N., Mackenzie, F. T., & Regnier, P. (2018). Continental shelves as a variable but increasing global sink for atmospheric carbon dioxide. Nature communications, 9(1), 454.
laruelle_NorthSea = ((-5.0, 56), (7.0, 64));
laruelle_EnglishChannel = ((-12.0, 46.5), (1.0, 51.0));
laruelle_SouthernGreenland = ((-53.0, 59.0), (-44.5, 62));
laruelle_AntarcticPeninsula = ((-69.0, -67.0), (-61.5, -62.5));
laruelle_LabradorSea = ((-67.0, 41.5), (-46.5, 48.25));
laruelle_CoastOfJapan = ((128.0, 27.5), (146.0, 43.0));
laruelle_CascadianShelf = ((-130.0, 42.0), (-119.75, 51.5));
laruelle_SouthAtlanticBight = ((-83.4, 25.0), (-77.0, 34.0));
laruelle_MidAtlanticBight = ((-77.0, 33.5), (-68.25, 41.25));
laruelle_BarentsSea = ((14.0, 67.5), (19.0, 71.0));
laruelle_TasmanianShelf = ((142.0, -45.25), (151.25, -39.0));
laruelle_PatagonianShelf = ((-69.0, -57.25), (-65.0, -53.0));
laruelle_BeringSea = ((-170.0, 50.0), (-160.0, 56.5));
#laruelle_WestIceland = ((-26.25, 61.25), (-20.5, 65.5)); #Do they mean this by Irminger Sea??
#laruelle_IrmingerSea = ((-34.0, 62.75), (-19.25, 65.5));
laruelle_IrmingerSea = ((-32.0, 61.25), (-25.0, 65.0));


#Verification of results against literature (table 1)
verif_PainterHebredes = ((-10.25, 55.25), (-5.75, 59.5));
verif_YuanSouthAtlanticBight = ((-82.5, 27.0), (-76.0, 34.0));
verif_FewingsMidAtlanticBight = ((-72.0, 39.0), (-65.5, 46.0)); #This had to be expanded considerably to reach our definition of the shelf edge
verif_WoodsonCaliforniaCoast = ((-122.5, 36.25), (-121.0, 37.25));
verif_WaiteEasternIndianOcean = ((112.0, -34.0), (117.0, -30.0));
verif_WieEastChinaSea = ((120.0, 21.0), (133.0, 34.0));


def point_in_area(llx, lly, urx, ury, px, py):
    if (px >= llx) & (py < lly) & (px < urx) & (py >= ury):
        return True;
    else:
        return False;
    
def empty_mask_function(data, params):
    return data;

def generic_area_mask(allowedAreas, data, params):
    keep = np.empty(len(data), dtype=bool);
    for pointIndex in range(0, len(data)):
        px = data[pointIndex].indexX;
        py = data[pointIndex].indexY;
        
        #Check each area to see if the point falls into one of them
        passed = False;
        for area in allowedAreas:
            #convert lon/lat to coordinates
            llx, lly = su.convert_lonlat_to_index(area[0][0], area[0][1], params.pixelRes, lon0=params.originLon, lat0=params.originLat);
            urx, ury = su.convert_lonlat_to_index(area[1][0], area[1][1], params.pixelRes, lon0=params.originLon, lat0=params.originLat);
            if point_in_area(llx, lly, urx, ury, px, py) == True:
                passed=True;
                break;
        keep[pointIndex] = passed;
        
    #Only keep points that fell into one of the above areas
    keep = np.array(keep);
    data = list(np.array(data)[keep]);
    
    return data;

def return_area_mask(allowedAreas, data, params):
    keep = np.empty(len(data), dtype=bool);
    for pointIndex in range(0, len(data)):
        px = data[pointIndex].indexX;
        py = data[pointIndex].indexY;
        
        #Check each area to see if the point falls into one of them
        passed = False;
        for area in allowedAreas:
            #convert lon/lat to coordinates
            llx, lly = su.convert_lonlat_to_index(area[0][0], area[0][1], params.pixelRes, lon0=params.originLon, lat0=params.originLat);
            urx, ury = su.convert_lonlat_to_index(area[1][0], area[1][1], params.pixelRes, lon0=params.originLon, lat0=params.originLat);
            if point_in_area(llx, lly, urx, ury, px, py) == True:
                passed=True;
                break;
        keep[pointIndex] = passed;
        
    keep = np.array(keep);
    return keep;

def global_shelf_mask_func(data, params):
    allowedAreas = []; #((llLon, llLat), (urLon, urLat))
    allowedAreas.append( area_EuropeanShelf ); #European shelf
    allowedAreas.append( area_LabradorSea ); #Labrador Sea
    allowedAreas.append( area_MidAtlanticBight ); #Mid-Atlantic Bight
    allowedAreas.append( area_CoastOfJapan ); #Coast of Japan
    allowedAreas.append( area_Patagonia ); #Patagonia
    allowedAreas.append( area_Tasmania ); #Tasmania
    allowedAreas.append( area_BarentsSea ); #Barents Sea
    allowedAreas.append( area_SouthAtlanticBight ); #South Atlantic Bight
    allowedAreas.append( area_SouthGreenland ); #South Greenland
    allowedAreas.append( area_AntarcticPeninsula ); #Antarctic Peninsula
    allowedAreas.append( area_BeringSeaWest ); #Bering Sea (West)
    allowedAreas.append( area_BeringSeaEast ); #Bering Sea (East)
    allowedAreas.append( area_IrmingerSea ); #Irminger Sea
    allowedAreas.append( area_CascianShelf ); #Cascian Shelf
    
#    allowedAreas.append( ((141.0, -45.0), (151.0, -39.0)) ); #Tasmania
#    allowedAreas.append( ((0.0, 60.0), (75.0, 90.0)) ); #Barents Sea
#    allowedAreas.append( ((-82.0, 26.0), (-74.0, 36.5)) ); #South Atlantic Bight
#    allowedAreas.append( ((-57.0, 59.0), (-18.5, 67.5)) ); #South Greenland
#    allowedAreas.append( ((-80.0, -75.0), (-50.0, -61.0)) ); #Antarctic Peninsula
#    allowedAreas.append( ((-180.0, 45.0), (-150.0, 61.0)) ); #Bering Sea (West)
#    allowedAreas.append( ((178.0, 45.0), (180.0, 61.0)) ); #Bering Sea (East)
#    allowedAreas.append( ((-57.0, 59.0), (-18.5, 67.5)) ); #Irminger Sea
#    allowedAreas.append( ((-140.0, 40.0), (-120.0, 54.0)) ); #Cascian Shelf
    

    keep = np.empty(len(data), dtype=bool);
    for pointIndex in range(0, len(data)):
        px = data[pointIndex].indexX;
        py = data[pointIndex].indexY;
        
        #Check each area to see if the point falls into one of them
        passed = False;
        for area in allowedAreas:
            #convert lon/lat to coordinates
            llx, lly = su.convert_lonlat_to_index(area[0][0], area[0][1], params.pixelRes);
            urx, ury = su.convert_lonlat_to_index(area[1][0], area[1][1], params.pixelRes);
            if point_in_area(llx, lly, urx, ury, px, py) == True:
                passed=True;
                break;
        keep[pointIndex] = passed;
        
    #Only keep points that fell into one of the above areas
    keep = np.array(keep);
    data = list(np.array(data)[keep]);
    
    return data;
    
#overrideMaskGapCheck will inhibit checks that the mask represents a contiguous section of the shelf edge
#allowMulti will return a list of tuples instead, allowing regions with gaps
def calculate_km_subsection_bounds_along_shelf(wholeShelfX, wholeShelfY, wholeShelfDistances, subsectionBounds, params, testPlot=False):
    #Calculate subset mask
    subsetMask = [False]*len(wholeShelfX);
    for pointIndex in range(0, len(wholeShelfX)):
        px = wholeShelfX[pointIndex];
        py = wholeShelfY[pointIndex];
        
        #convert lon/lat to coordinates
        llx, lly = su.convert_lonlat_to_index(subsectionBounds[0][0], subsectionBounds[0][1], params.pixelRes, lon0=params.originLon, lat0=params.originLat);
        urx, ury = su.convert_lonlat_to_index(subsectionBounds[1][0], subsectionBounds[1][1], params.pixelRes, lon0=params.originLon, lat0=params.originLat);
        
        #Check each area to see if the point falls into one of them and set the mask
        if point_in_area(llx, lly, urx, ury, px, py) == True:
            subsetMask[pointIndex] = True;

    #find the first and last instance of 'True' for each contiguous span of 'True'
    ifirsts = [];
    ilasts = [];
    prev=False;
    for i in range(0, len(subsetMask)):
        if subsetMask[i] == True and prev==False:
            ifirsts.append(i);
        if subsetMask[i] == False and prev==True:
            ilasts.append(i-1);
        prev=subsetMask[i];
    if prev==True: #last one won't trigger an end of a contiguous block
        ilasts.append(i);
    if len(ilasts) != len(ifirsts): #sanity check
        raise RuntimeError("Unequal starts and ends of contiguous blocks in subset filter.");
    
    #convert mask to array
    subsetMask = np.array(subsetMask);
    cumulativeDists = np.cumsum(wholeShelfDistances);
    
    #calculate intercept distances for each contiguous span
    distances = [];
    for i in range(0, len(ifirsts)):
        dist1 = cumulativeDists[ifirsts[i]];
        dist2 = cumulativeDists[ilasts[i]];
        distances.append( (dist1, dist2) );
    
    #plot to perform manual sanity check
    if testPlot:
        import matplotlib.pyplot as plt;
        from netCDF4 import Dataset;
        
        depth = np.flipud(Dataset("data/GEBCO_bathymetry_0.25x0.25deg.nc", 'r').variables["mean_depth"][:]);
        
        plt.figure();
        plt.imshow(depth);
        for i in range(len(distances)):
            plt.plot(wholeShelfX[ifirsts[i]:ilasts[i]+1], wholeShelfY[ifirsts[i]:ilasts[i]+1], 'y', linewidth=5, alpha=1.0); #plot just subsection with thicker line
            #plt.plot(wholeShelfX[iFirst:iLast+1], wholeShelfY[iFirst:iLast+1], 'y', linewidth=5, alpha=1.0); #plot just subsection with thicker line
        plt.plot(wholeShelfX, wholeShelfY, 'r');
        plt.xlim(wholeShelfX.min(), wholeShelfX.max());
        plt.ylim(wholeShelfY.max(), wholeShelfY.min());
    
    #return the intercept distances along the shelf boundary
    return distances, subsetMask;



#def global_shelf_mask_func(data, params):
#    allowedAreas = []; #((llLon, llLat), (urLon, urLat))
#    allowedAreas.append( ((-15.0, 40.0), (20.0, 70.0)) ); #European shelf
#    allowedAreas.append( ((-70.0, 50.0), (-40.0, 70.0)) ); #Labrador Sea
#    allowedAreas.append( ((-76.0, 34.0), (-68.0, 42.0)) ); #Mid-Atlantic Bight
#    allowedAreas.append( ((128.0, 30.0), (149.0, 44.0)) ); #Coast of Japan
#    allowedAreas.append( ((141.0, -45.0), (151.0, -39.0)) ); #Tasmania
#    allowedAreas.append( ((0.0, 60.0), (75.0, 90.0)) ); #Barents Sea
#    allowedAreas.append( ((-82.0, 26.0), (-74.0, 36.5)) ); #South Atlantic Bight
#    allowedAreas.append( ((-57.0, 59.0), (-18.5, 67.5)) ); #South Greenland
#    allowedAreas.append( ((-80.0, -75.0), (-50.0, -61.0)) ); #Antarctic Peninsula
#    allowedAreas.append( ((-180.0, 45.0), (-150.0, 61.0)) ); #Bering Sea (West)
#    allowedAreas.append( ((178.0, 45.0), (180.0, 61.0)) ); #Bering Sea (East)
#    allowedAreas.append( ((-57.0, 59.0), (-18.5, 67.5)) ); #Irminger Sea
#    allowedAreas.append( ((-140.0, 40.0), (-120.0, 54.0)) ); #Cascian Shelf
#    allowedAreas.append( ((-87.0, -58.0), (-59.0, -41.0)) ); #Patagonia
#
#    for pathIndex in range(len(data.xs)):
#        keep = np.empty(len(data.xs[pathIndex]), dtype=bool);
#        for pointIndex in range(0, len(data.xs[pathIndex])):
#            px = data.xs[pathIndex][pointIndex];
#            py = data.ys[pathIndex][pointIndex];
#            
#            #Check each area to see if the point falls into one of them
#            passed = False;
#            for area in allowedAreas:
#                #convert lon/lat to coordinates
#                llx, lly = su.convert_lonlat_to_index(area[0][0], area[0][1], params.pixelRes);
#                urx, ury = su.convert_lonlat_to_index(area[1][0], area[1][1], params.pixelRes);
#                if point_in_area(llx, lly, urx, ury, px, py) == True:
#                    passed=True;
#                    break;
#            keep[pointIndex] = passed;
#        
#        #Only keep points that fell into one of the above areas
#        keep = np.array(keep);
#        data.xs[pathIndex] = list(np.array(data.xs[pathIndex])[keep]);
#        data.ys[pathIndex] = list(np.array(data.ys[pathIndex])[keep]);
#        data.coordinatesLists[pathIndex] = list(np.array(data.coordinatesLists[pathIndex])[keep]);
#    
#    return data;
            


