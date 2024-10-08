#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  3 10:46:32 2018

@author: rr
Modified by Daniel J. Ford (d.ford@exeter.ac.uk)
Date: 03/2023
Changes:
- Updated file paths for local system (L31-37)
- Changed Reynolds sst to SST-CCI (L39)
"""

import skim_utilities as su;
import mask_functions as mf;
from string import Template;

#baseline parameters, used as the basis for all other parameter sets.
def get_baseline_params():
    params = su.QuickStruct();
    params.numLineApproximationFunction = su.num_line_approximations; #Function to determine the number of line segments to use
    params.minContourPathSizeShallow = 150; #Contour paths which contain less than this number of coordinate pairs will be ignored
    params.minContourPathSizeDeep = 100; #Contour paths which contain less than this number of coordinate pairs will be ignored.
                                        #To be on the safe side this should be less restrictive than the shallow threshold

    #Thresholds used to filter for regions when stokes has an appreciable effect.
    params.stokesMaskWindStressThreshold = 0.03; #N m^-2 (wind must be below this threshold)
    params.stokesMaskSigWaveHeightThreshold = 2.0; #metres (wave height must be above this threshold)


    #file paths
    params.ekmanTemplate = Template("D:/SKIM_Paper_Data/SKIM/processed_data/monthly_means/ekman_monthly_mean/${YYYY}${MM}_ekman_surface_monthly_mean.nc");
    params.geostrophicTemplate = Template("D:/SKIM_Paper_Data/SKIM/processed_data/monthly_means/geostrophic_monthly_mean/${YYYY}${MM}_geostrophic_monthly_mean.nc");
    params.stokesTemplate = Template("D:/SKIM_Paper_Data/SKIM/processed_data/monthly_means/stokes_monthly_mean/${YYYY}${MM}_stokes_monthly_mean.nc");
    params.wavewatchWndTemplate = Template("D:/SKIM_Paper_Data/SKIM/processed_data/monthly_means/wavewatch_wnd_monthly_means/${YYYY}${MM}_wind_monthly_mean.nc");
    params.wavewatchHsTemplate = Template("D:/SKIM_Paper_Data/SKIM/processed_data/monthly_means/wavewatch_Hs_monthly_means/${YYYY}${MM}_Hs_monthly_mean.nc");
    params.skimulatorTemplate = Template("D:/SKIM_Paper_Data/SKIM/processed_data/monthly_means/simulated_skim/${YYYY}${MM}_monthly_mean_scisoc_atne.nc");
    #params.reynoldsSSTTemplate = Template("D:/SKIM_Paper_Data/SKIM/reynolds_sst/reynolds_avhrr_only_monthly_0.25_calculated_tmh/${YYYY}/${YYYY}${MM}01_OCF-SST-GLO-1M-100-REYNOLDS_0.25deg_TMH.nc");
    params.reynoldsSSTTemplate = Template("D:/Data/SST-CCI/MONTHLY_025_DEG/${YYYY}/${YYYY}${MM}_ESA_CCI_MONTHLY_SST_025_deg.nc");
    #params.ekmanTemplate = Template("/media/verwirrt/Backup_THolding/THoldingBackup/SKIM/processed_data/monthly_means/ekman_monthly_mean/${YYYY}${MM}_ekman_surface_monthly_mean.nc");
    #params.geostrophicTemplate = Template("/media/verwirrt/Backup_THolding/THoldingBackup/SKIM/processed_data/monthly_means/geostrophic_monthly_mean/${YYYY}${MM}_geostrophic_monthly_mean.nc");
    #params.stokesTemplate = Template("/media/verwirrt/Backup_THolding/THoldingBackup/SKIM/processed_data/monthly_means/stokes_monthly_mean/${YYYY}${MM}_stokes_monthly_mean.nc");
    #params.wavewatchWndTemplate = Template("/media/verwirrt/Backup_THolding/THoldingBackup/SKIM/processed_data/monthly_means/wavewatch_wnd_monthly_means/wnd_monthly_means_${YYYY}${MM}.nc");
    #params.wavewatchHsTemplate = Template("/media/verwirrt/Backup_THolding/THoldingBackup/SKIM/processed_data/monthly_means/wavewatch_Hs_monthly_means/Hs_monthly_means_${YYYY}${MM}.nc");
    #params.skimulatorTemplate = Template("/media/verwirrt/Backup_THolding/THoldingBackup/SKIM/processed_data/monthly_means/simulated_skim/$(YYYY)$(MM)_monthly_mean_scisoc_atne.nc");

    return params;

def get_baseline_params_cmems():
    params = su.QuickStruct();
    params.numLineApproximationFunction = su.num_line_approximations; #Function to determine the number of line segments to use
    params.minContourPathSizeShallow = 150; #Contour paths which contain less than this number of coordinate pairs will be ignored
    params.minContourPathSizeDeep = 100; #Contour paths which contain less than this number of coordinate pairs will be ignored.
                                        #To be on the safe side this should be less restrictive than the shallow threshold

    #Thresholds used to filter for regions when stokes has an appreciable effect.
    params.stokesMaskWindStressThreshold = 0.03; #N m^-2 (wind must be below this threshold)
    params.stokesMaskSigWaveHeightThreshold = 2.0; #metres (wave height must be above this threshold)


    #file paths
    params.ekmanTemplate = Template("E:/SKIM_Paper_Data/SKIM/downloaded_data/CMEMS/${YYYY}/${YYYY}_${MM}_CMEMS_GLOBCURRENT.nc");
    params.geostrophicTemplate = Template("E:/SKIM_Paper_Data/SKIM/downloaded_data/CMEMS/${YYYY}/${YYYY}_${MM}_CMEMS_GLOBCURRENT.nc");
    params.stokesTemplate = Template("E:/SKIM_Paper_Data/SKIM/processed_data/monthly_means/stokes_monthly_mean/${YYYY}${MM}_stokes_monthly_mean.nc");
    params.wavewatchWndTemplate = Template("E:/SKIM_Paper_Data/SKIM/processed_data/monthly_means/wavewatch_wnd_monthly_means/${YYYY}${MM}_wind_monthly_mean.nc");
    params.wavewatchHsTemplate = Template("E:/SKIM_Paper_Data/SKIM/processed_data/monthly_means/wavewatch_Hs_monthly_means/${YYYY}${MM}_Hs_monthly_mean.nc");
    params.skimulatorTemplate = Template("E:/SKIM_Paper_Data/SKIM/processed_data/monthly_means/simulated_skim/${YYYY}${MM}_monthly_mean_scisoc_atne.nc");
    #params.reynoldsSSTTemplate = Template("D:/SKIM_Paper_Data/SKIM/reynolds_sst/reynolds_avhrr_only_monthly_0.25_calculated_tmh/${YYYY}/${YYYY}${MM}01_OCF-SST-GLO-1M-100-REYNOLDS_0.25deg_TMH.nc");
    params.reynoldsSSTTemplate = Template("D:/Data/SST-CCI/MONTHLY_025_DEG/${YYYY}/${YYYY}${MM}_ESA_CCI_MONTHLY_SST_025_deg.nc");

    return params;


#returns QuickStruct which contains the parameters used for the original development of the algorithms
def __old__get_European_shelf_params(): #Whole European shelf
    params = get_baseline_params();
    params.pixelRes = (0.25, 0.25);
    params.timeRes = "monthly";
    params.originLon = -20.0; #at (0,0)
    params.maxLon = 22.5;
    params.originLat = 72.0; #at (0,0)
    params.maxLat = 35.5;
    params.contourPathFile = "european_shelf_contour_paths.p";
    params.contourPathFileDeep = "european_shelf_contour_paths_deep.p";
    params.paramsetName = "old_europeanshelf";

    originLon, originLat = su.convert_lonlat_to_index(params.originLon, params.originLat, params.pixelRes);
    maxLon, maxLat = su.convert_lonlat_to_index(params.maxLon, params.maxLat, params.pixelRes);
    params.ilonRange = (originLon, maxLon); #(600, 850)
    params.ilatRange = (originLat, maxLat); #(50, 225)
    params.shallowDepth = 500.0;
    params.deepDepth = params.shallowDepth+300.0;
    params.contourPathsShallow = [0];
    params.contourPathsDeep = [0];
    params.numberOfBoundryApproximationLines = [64]; #how many lines should be used to approximate the shelf edge boundry

    params.start_year = 2014;
    params.end_year = 2015;
    params.start_month = 0;
    params.end_month = 12;

    params.contourMaskFunc = mf.empty_mask_function;

    return params;

def get_example_method_plot_params(): #Whole European shelf
    params = get_baseline_params();
    params.pixelRes = (0.25, 0.25);
    params.timeRes = "monthly";
    params.originLon = -20.0; #at (0,0)
    params.maxLon = 0.0;
    params.originLat = 59.0; #at (0,0)
    params.maxLat = 45.0;
    params.contourPathFile = "method_example_contour_paths.p";
    params.contourPathFileDeep = "method_example_contour_paths_deep.p";
    params.paramsetName = "method_example";

    originLon, originLat = su.convert_lonlat_to_index(params.originLon, params.originLat, params.pixelRes);
    maxLon, maxLat = su.convert_lonlat_to_index(params.maxLon, params.maxLat, params.pixelRes);
    params.ilonRange = (originLon, maxLon); #(600, 850)
    params.ilatRange = (originLat, maxLat); #(50, 225)
    params.shallowDepth = 500.0;
    params.deepDepth = params.shallowDepth+100.0;
    params.contourPathsShallow = [1];
    params.contourPathsDeep = [2];
    params.numberOfBoundryApproximationLines = [8]; #how many lines should be used to approximate the shelf edge boundry

    params.start_year = 2014;
    params.end_year = 2015;
    params.start_month = 0;
    params.end_month = 12;

    params.contourMaskFunc = mf.empty_mask_function;

    return params;

def get_example_method_plot_params_64(): #Whole European shelf
    params = get_baseline_params();
    params.pixelRes = (0.25, 0.25);
    params.timeRes = "monthly";
    params.originLon = -20.0; #at (0,0)
    params.maxLon = 0.0;
    params.originLat = 59.0; #at (0,0)
    params.maxLat = 45.0;
    params.contourPathFile = "method_example_64_contour_paths.p";
    params.contourPathFileDeep = "method_example_64_contour_paths_deep.p";
    params.paramsetName = "method_example_64";

    originLon, originLat = su.convert_lonlat_to_index(params.originLon, params.originLat, params.pixelRes);
    maxLon, maxLat = su.convert_lonlat_to_index(params.maxLon, params.maxLat, params.pixelRes);
    params.ilonRange = (originLon, maxLon); #(600, 850)
    params.ilatRange = (originLat, maxLat); #(50, 225)
    params.shallowDepth = 500.0;
    params.deepDepth = params.shallowDepth+100.0;
    params.contourPathsShallow = [1];
    params.contourPathsDeep = [2];
    #params.numberOfBoundryApproximationLines = [64]; #how many lines should be used to approximate the shelf edge boundry
    params.numLineApproximationFunction = lambda x, y : 64;
    params.start_year = 2014;
    params.end_year = 2015;
    params.start_month = 0;
    params.end_month = 12;

    params.contourMaskFunc = mf.empty_mask_function;

    return params;

#returns QuickStruct which contains the parameters used for the original development of the algorithms
def get_European_shelf_params(): #Whole European shelf
    params = get_baseline_params();
    params.pixelRes = (0.25, 0.25);
    params.timeRes = "monthly";
    params.originLon = -25.0; #at (0,0)
    params.maxLon = 0.0;
    params.originLat = 62.0; #at (0,0)
    params.maxLat = 40.0;
    params.contourPathFile = "european_shelf_contour_paths.p";
    params.contourPathFileDeep = "european_shelf_contour_paths_deep.p";
    params.paramsetName = "europeanshelf";

    originLon, originLat = su.convert_lonlat_to_index(params.originLon, params.originLat, params.pixelRes);
    maxLon, maxLat = su.convert_lonlat_to_index(params.maxLon, params.maxLat, params.pixelRes);
    params.ilonRange = (originLon, maxLon); #(600, 850)
    params.ilatRange = (originLat, maxLat); #(50, 225)
    params.shallowDepth = 500.0;
    params.deepDepth = params.shallowDepth+100.0;
    params.contourPathsShallow = [1];
    params.contourPathsDeep = [1];
    #params.numberOfBoundryApproximationLines = [64]; #how many lines should be used to approximate the shelf edge boundry
    params.numLineApproximationFunction = lambda x, y : 64; #Function to determine the number of line segments to use

    #params.start_year = 2014;
    #params.end_year = 2015;
    #params.start_month = 0;
    #params.end_month = 12;
    params.start_year = 2011;
    params.end_year = 2012;
    params.start_month = 9;
    params.end_month = 9;

    params.contourMaskFunc = mf.empty_mask_function;
    params.yeariMonths = [(params.start_year, imonth) for imonth in range(params.start_month, 12)] + [(params.end_year, imonth) for imonth in range(0, params.end_month)];

    return params;


#Global analysis using multiple years. This is the main analysis.
def get_global_params(cmems = False):
    if cmems:
        params = get_baseline_params_cmems();
    else:
        params = get_baseline_params();
    params.pixelRes = (0.25, 0.25);
    params.timeRes = "monthly";
    params.originLon = -180.0; #at (0,0)
    params.maxLon = 179.75;
    params.originLat = 89.75; #at (0,0)
    params.maxLat = -90.0;

    originLon, originLat = su.convert_lonlat_to_index(params.originLon, params.originLat, params.pixelRes);
    maxLon, maxLat = su.convert_lonlat_to_index(params.maxLon, params.maxLat, params.pixelRes);
    params.ilonRange = (originLon, maxLon); #(600, 850)
    params.ilatRange = (originLat, maxLat); #(50, 225)
    params.shallowDepth = 500.0;
    params.deepDepth = params.shallowDepth+100.0;
    params.contourPathsShallow = [13, 15, 39, 68, 75, 101, 107, 110, 113, 115, 136, 175, 257, 177, 180];
    params.contourPathsDeep = [2, 3, 21, 41, 45, 67, 66, 74, 77, 79, 97, 129, 168, 131, 133];
    params.numberOfBoundryApproximationLines = [32, 16, 32, 16, 4, 8, 8, 32, 16, 8, 32, 32, 8, 8, 8]; #how many lines should be used to approximate the shelf edge boundry

    params.contourMaskFunc = mf.global_shelf_mask_func;

    params.start_year = 1993;
    params.end_year = 2016;
    params.start_month = 0;
    params.end_month = 12;
    if cmems:
        params.contourPathFile = "cmems_global_shelf_contour_paths.p";
        params.contourPathFileDeep = "cmems_global_shelf_contour_paths_deep.p";
        params.paramsetName = "cmems_global";
    else:
        params.contourPathFile = "global_shelf_contour_paths.p";
        params.contourPathFileDeep = "global_shelf_contour_paths_deep.p";
        params.paramsetName = "global";

    return params

#Global analysis using multiple years. This is the main analysis.
def get_global_params_glory(res = False):
    params = su.QuickStruct();
    params.numLineApproximationFunction = su.num_line_approximations; #Function to determine the number of line segments to use
    params.minContourPathSizeShallow = 150; #Contour paths which contain less than this number of coordinate pairs will be ignored
    params.minContourPathSizeDeep = 100; #Contour paths which contain less than this number of coordinate pairs will be ignored.
                                        #To be on the safe side this should be less restrictive than the shallow threshold

    #Thresholds used to filter for regions when stokes has an appreciable effect.
    params.stokesMaskWindStressThreshold = 0.03; #N m^-2 (wind must be below this threshold)
    params.stokesMaskSigWaveHeightThreshold = 2.0; #metres (wave height must be above this threshold)


    #file paths
    params.ekmanTemplate = Template("D:/SKIM_Paper_Data/SKIM/downloaded_data/CMEMS/${YYYY}/${YYYY}_${MM}_CMEMS_GLOBCURRENT.nc");

    params.stokesTemplate = Template("D:/SKIM_Paper_Data/SKIM/processed_data/monthly_means/stokes_monthly_mean/${YYYY}${MM}_stokes_monthly_mean.nc");
    params.wavewatchWndTemplate = Template("D:/SKIM_Paper_Data/SKIM/processed_data/monthly_means/wavewatch_wnd_monthly_means/${YYYY}${MM}_wind_monthly_mean.nc");
    params.wavewatchHsTemplate = Template("D:/SKIM_Paper_Data/SKIM/processed_data/monthly_means/wavewatch_Hs_monthly_means/${YYYY}${MM}_Hs_monthly_mean.nc");
    params.skimulatorTemplate = Template("D:/SKIM_Paper_Data/SKIM/processed_data/monthly_means/simulated_skim/${YYYY}${MM}_monthly_mean_scisoc_atne.nc");
    #params.reynoldsSSTTemplate = Template("D:/SKIM_Paper_Data/SKIM/reynolds_sst/reynolds_avhrr_only_monthly_0.25_calculated_tmh/${YYYY}/${YYYY}${MM}01_OCF-SST-GLO-1M-100-REYNOLDS_0.25deg_TMH.nc");
    params.reynoldsSSTTemplate = Template("D:/Data/SST-CCI/MONTHLY_025_DEG/${YYYY}/${YYYY}${MM}_ESA_CCI_MONTHLY_SST_025_deg.nc");
    if res:

        params.pixelRes = (1.0/12, 1.0/12);
        params.timeRes = "monthly";
        params.originLon = -180.0; #at (0,0)
        params.maxLon = 180;
        params.originLat = 90; #at (0,0)
        params.maxLat = -80.0;
        params.geostrophicTemplate = Template("D:/SKIM_Paper_Data/SKIM/downloaded_data/CMEMS_GLORYSV12/${YYYY}/${YYYY}_${MM}_CMEMS_GLORYSV12.nc");
    else:
        params.pixelRes = (0.25, 0.25);
        params.timeRes = "monthly";
        params.originLon = -180.0; #at (0,0)
        params.maxLon = 179.75;
        params.originLat = 89.75; #at (0,0)
        params.maxLat = -90.0;
        params.geostrophicTemplate = Template("D:/SKIM_Paper_Data/SKIM/downloaded_data/CMEMS_GLORYSV12_025/${YYYY}/${YYYY}_${MM}_CMEMS_GLORYSV12_0.25_deg.nc");

    originLon, originLat = su.convert_lonlat_to_index(params.originLon, params.originLat, params.pixelRes,lon0 = -180,lat0=90.0-params.pixelRes[0]);
    maxLon, maxLat = su.convert_lonlat_to_index(params.maxLon, params.maxLat, params.pixelRes,lon0 = -180,lat0=90.0-params.pixelRes[0]);
    params.ilonRange = (originLon, maxLon); #(600, 850)
    params.ilatRange = (originLat, maxLat); #(50, 225)
    params.shallowDepth = 500.0;
    params.deepDepth = params.shallowDepth+100.0;
    params.contourPathsShallow = [13, 15, 39, 68, 75, 101, 107, 110, 113, 115, 136, 175, 257, 177, 180];
    params.contourPathsDeep = [2, 3, 21, 41, 45, 67, 66, 74, 77, 79, 97, 129, 168, 131, 133];
    params.numberOfBoundryApproximationLines = [32, 16, 32, 16, 4, 8, 8, 32, 16, 8, 32, 32, 8, 8, 8]; #how many lines should be used to approximate the shelf edge boundry

    params.contourMaskFunc = mf.global_shelf_mask_func;

    params.start_year = 1993;
    params.end_year = 1994;
    params.start_month = 0;
    params.end_month = 12;
    if res:
        params.contourPathFile = "cmems_glory_high_global_shelf_contour_paths.p";
        params.contourPathFileDeep = "cmems_glory_high_global_shelf_contour_paths_deep.p";
        params.paramsetName = "cmems_glory_high_global";
    else:
        params.contourPathFile = "cmems_glory_low_global_shelf_contour_paths.p";
        params.contourPathFileDeep = "cmems_glory_low_global_shelf_contour_paths_deep.p";
        params.paramsetName = "cmems_glory_low_global";

    return params
#Global analysis using multiple years. This is the main analysis.
def get_global_test_params():
    params = get_baseline_params();
    params.pixelRes = (0.25, 0.25);
    params.timeRes = "monthly";
    params.originLon = -180.0; #at (0,0)
    params.maxLon = 179.75;
    params.originLat = 89.75; #at (0,0)
    params.maxLat = -90.0;

    originLon, originLat = su.convert_lonlat_to_index(params.originLon, params.originLat, params.pixelRes);
    maxLon, maxLat = su.convert_lonlat_to_index(params.maxLon, params.maxLat, params.pixelRes);
    params.ilonRange = (originLon, maxLon); #(600, 850)
    params.ilatRange = (originLat, maxLat); #(50, 225)
    params.shallowDepth = 500.0;
    params.deepDepth = params.shallowDepth+100.0;
#    params.contourPathsShallow = [13, 15, 39, 68, 75, 101, 107, 110, 113, 115, 136, 175, 257, 177, 180];
#    params.contourPathsDeep = [2, 3, 21, 41, 45, 67, 66, 74, 77, 79, 97, 129, 168, 131, 133];
#    params.numberOfBoundryApproximationLines = [32, 16, 32, 16, 4, 8, 8, 32, 16, 8, 32, 32, 8, 8, 8]; #how many lines should be used to approximate the shelf edge boundry

    params.contourMaskFunc = mf.global_shelf_mask_func;

    params.start_year = 1993;
    params.end_year = 1993;
    params.start_month = 0;
    params.end_month = 12;
    params.contourPathFile = "global_test_shelf_contour_paths.p";
    params.contourPathFileDeep = "global_test_shelf_contour_paths_deep.p";
    params.paramsetName = "global_test";

    return params

##Global analysis using multiple years with a shallower 300m shelf contour (as opposed) and 800m deep
#def get_global_multiyear_300m_shelf_params():
#    params = get_global_multiyear_params();
#
#    params.shallowDepth = 300.0;
#    params.deepDepth = params.shallowDepth+300.0;
#    params.contourPathFile = "global_300m_shelf_contour_paths.p";
#    params.contourPathFileDeep = "global_300m_shelf_contour_paths_deep.p";
#    params.contourPathsShallow = [0, 1, 2, 10, 11, 15, 22, 24, 25, 26, 31, 35, 36, 40, 42, 44, 52, 53, 69, 79, 80, 88, 89, 98, 103, 112, 122, 127, 130, 132, 134, 140, 141, 147, 151, 154, 156, 165, 171, 176, 181, 182, 183, 187, 190, 195, 197, 207, 218, 224, 226, 230, 234, 235, 236, 237, 238, 239, 244, 248, 251, 253, 254, 256, 257, 259, 262, 265, 266, 269, 276, 277, 280, 281, 292, 295, 306, 311, 312, 326, 328, 329, 330, 332, 334, 339, 340, 341, 347, 350, 354, 356, 359, 365, 366, 376, 379, 407, 413, 415, 417, 426, 427, 430, 433, 440, 454, 455, 465, 473, 479, 487, 497, 503, 514, 532, 537, 547, 549, 556, 568, 570, 580, 585, 586, 592, 593, 600, 604, 605, 606, 608, 614, 631, 633, 638, 643, 649, 651, 663, 665, 672, 676, 677, 678, 680, 684, 696, 702, 704, 705, 731, 774, 779, 811, 814, 819, 822, 827, 832, 845, 851, 855, 860, 865, 892, 898, 904, 911, 918, 949, 952, 955, 961, 967, 968, 975, 976, 978, 980, 985, 995, 998];
#    params.contourPathsDeep =  [0, 1, 2, 3, 4, 21, 25, 26, 27, 31, 32, 33, 36, 37, 38, 40, 41, 44, 45, 46, 49, 50, 51, 52, 56, 57, 59, 60, 61, 63, 64, 65, 66, 67, 68, 70, 71, 72, 74, 75, 77, 78, 79, 83, 87, 89, 90, 95, 96, 97, 98, 100, 109, 111, 113, 115, 116, 117, 118, 119, 120, 121, 122, 124, 125, 127, 128, 129, 131, 132, 133, 134, 137, 138, 139, 140, 141, 142, 144, 149, 150, 151, 152, 154, 155, 157, 158, 161, 162, 165, 166, 168, 169, 172, 173, 174, 176, 178, 185, 188, 190, 191, 198, 201, 202, 203, 210, 217, 219, 222, 223, 224, 225, 226, 230, 232, 235, 240, 241, 242, 243, 249, 253, 255, 257, 259, 269, 272, 273, 278, 289, 290, 298, 301, 303, 304, 305, 315, 316, 318, 324, 327, 328, 330, 331, 343, 345, 346, 355, 363, 365, 370, 373, 380, 383, 394, 400, 406, 410, 416, 421, 423, 429, 433, 435, 436, 440, 441, 442, 443, 444, 445, 448, 451, 452, 454, 458, 463, 464, 478, 484, 488, 491, 496, 502, 506, 508, 511, 513, 515, 524, 529, 530, 533, 542, 546, 549, 550, 555, 557, 559, 563, 566, 569, 573, 575, 576, 584, 585, 588, 590, 600, 602, 603, 605, 608, 612, 613, 615, 616, 619, 621, 625];
#    params.numberOfBoundryApproximationLines = [3, 5, 9, 5, 2, 2, 3, 8, 2, 12, 9, 4, 3, 95, 20, 7, 3, 3, 3, 24, 16, 3, 5, 2, 23, 2, 3, 2, 4, 2, 2, 4, 2, 3, 4, 7, 13, 2, 2, 26, 5, 5, 6, 6, 9, 8, 11, 7, 2, 2, 2, 3, 6, 2, 2, 2, 3, 3, 3, 2, 3, 2, 45, 15, 2, 2, 12, 9, 8, 2, 3, 4, 4, 3, 3, 31, 2, 3, 7, 2, 5, 5, 10, 14, 2, 3, 2, 2, 2, 5, 24, 9, 6, 9, 5, 3, 3, 2, 3, 4, 2, 2, 3, 2, 2, 2, 2, 2, 3, 2, 2, 3, 3, 4, 2, 2, 2, 2, 3, 2, 4, 2, 4, 2, 5, 3, 3, 9, 3, 6, 2, 2, 5, 11, 2, 3, 2, 3, 4, 3, 3, 2, 2, 2, 2, 3, 3, 2, 3, 2, 4, 2, 3, 2, 2, 3, 4, 2, 2, 2, 3, 3, 2, 2, 2, 3, 2, 2, 2, 2, 2, 2, 4, 2, 2, 2, 4, 4, 3, 3, 2, 2, 2]; #how many lines should be used to approximate the shelf edge boundry
#
#    params.paramsetName = "global_300m";
#
#    return params


#Set the parameters to be used by
def get_current_params():
    #return get_example_method_plot_params();

    #return get_European_shelf_params(); #European shelf
    return get_global_params();
