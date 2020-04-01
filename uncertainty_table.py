#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May  8 15:14:16 2019
V2: corrected percentage calculation.
V3: used difference from baseline (RMSD=0) instead of mean values (which would be the same due to convergence to the mean)
    percentages outputted as absolutes...
@author: rr
"""

import numpy as np;
import pandas as pd;

regionNames = ["NorthSea", "EnglishChannel", "MidAtlanticBight", "CoastOfJapan", "PatagonianShelf", "BeringSea", "AntarcticPeninsula", "LabradorSea", "TasmanianShelf", "BarentsSea", "SouthAtlanticBight", "SouthernGreenland", "CascadianShelf", "IrmingerSea"];
regionNames123 = [regionName+"123" for regionName in regionNames];
regionLabels = ["North Sea", "English Channel", "Mid-Atlantic Bight", "Coast of Japan", "Patagonian Shelf", "Bering Sea", "Antarctic Peninsula", "Labrador Sea", "Tasmanian Shelf", "Barents Sea", "South Atlantic Bight", "Southern Greenland", "Cascadian Shelf", "Irminger Sea"];

#df = pd.read_table("/home/rr/Files/Tasks/20180914_SKIM/plots/uncertainty_analysis/regional_error_analysis_global_small.csv", sep=",");
df = pd.read_table("/home/rr/Files/Tasks/20180914_SKIM/results/sensitivity_to_uncertainty/regional_error_analysis_global_v4_1000_reps.csv", sep=",");
df.columns

#Extract the base 0 error across-shelf current for each region.
baseValues = [df[regionName+"123_mean"][0] for regionName in regionNames];


#Calculate percentage intervals...
for i, regionName in enumerate(regionNames):
    percentageInterval = df[regionName+"123_diffSD"] / baseValues[i] * 100.0;
    df[regionName+"123_diffpercentSD"] = percentageInterval;


newdf = pd.DataFrame();
newdf["Error"] = df["RMSE"]


for i, regionName in enumerate(regionNames):
    means = df[regionName+"123_mean"].values;
    sds = df[regionName+"123SD"].values;
    strs = [format(means[v], ".2E") + " +/- " + format(sds[v], ".2E") for v in range(len(means))];
    newdf[regionLabels[i]+" diff"] = strs;


for i, regionName in enumerate(regionNames):
    means = np.abs(df[regionName+"123_diffpercent"].values);
    sds = np.abs(df[regionName+"123_diffpercentSD"].values);
    strs = [format(means[v], ".2f") + " +/- " + format(sds[v], ".2f") for v in range(len(means))];
    newdf[regionLabels[i]+" perc"] = strs;

newdft = newdf.transpose();
#newdft["Region"] = regionLabels
newdft.shape


newdft.to_csv("/home/rr/Files/Tasks/20180914_SKIM/plots/uncertainty_analysis/regional_error_analysis_global_small_table_v4_3.csv", index=True);
