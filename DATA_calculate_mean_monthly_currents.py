#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 20 08:44:21 2018

@author: Tom Holding
Modified by Daniel J. Ford (d.ford@exeter.ac.uk)
Date: 03/2023
- Updated templates
- Added process_wind_month() function to average wind speed from WaveWatch
- Modified data loading components for process_stokes_month() and process_wind_month() due to netCDF files containing all data for a month.
"""

from netCDF4 import Dataset;
import numpy as np;
import matplotlib.pyplot as plt;
import calendar;
from string import Template;
import python_util.skim_utilities as su;

processEkman = False;
processEkman15m = False;
processGeostrophic = False;
processStokes = True;
processwind = False;
processHs = False

processMonthly = True;
processDaily = False;

ekmanInputTemplate = Template("D:/SKIM_Paper_Data/SKIM/downloaded_data/mean_daily_ekman_currents/${YYYY}${MM}${DD}-GLOBCURRENT-L4-CURekm_hs-ERAWS_EEM-v03.0-fv01.0.nc");
ekman15mInputTemplate = Template("D:/SKIM_Paper_Data/SKIM/downloaded_data/mean_daily_ekman_currents_15m/${YYYY}${MM}${DD}-GLOBCURRENT-L4-CURekm_15m-ERAWS_EEM-v03.0-fv01.0.nc");
geostrophicInputTemplate = Template("D:/SKIM_Paper_Data/SKIM/downloaded_data/mean_daily_geostrophic_currents/${YYYY}${MM}${DD}000000-GLOBCURRENT-L4-CURgeo_0m-ALT_OI-v03.0-fv01.0.nc");
stokesInputTemplate = Template("D:/SKIM_Paper_Data/SKIM/downloaded_data/stokes_data_3h/WW3-GLOB-30M_${YYYY}${MM}_uss.nc");
windInputTemplate = Template("D:/SKIM_Paper_Data/SKIM/downloaded_data/wavewatch_wnd/WW3-GLOB-30M_${YYYY}${MM}_wnd.nc");
hsInputTemplate = Template("D:/SKIM_Paper_Data/SKIM/downloaded_data/wavewatch_Hs/WW3-GLOB-30M_${YYYY}${MM}_hs.nc");

ekmanMonthlyOutputTemplate = Template("D:/SKIM_Paper_Data/SKIM/processed_data/monthly_means/ekman_monthly_mean/${YYYY}${MM}_ekman_surface_monthly_mean.nc");
ekman15mMonthlyOutputTemplate = Template("D:/SKIM_Paper_Data/SKIM/processed_data/monthly_means/ekman_15m_monthly_mean/${YYYY}${MM}_ekman_15m_monthly_mean.nc");
geostrophicMonthlyOutputTemplate = Template("D:/SKIM_Paper_Data/SKIM/processed_data/monthly_means/geostrophic_monthly_mean/${YYYY}${MM}_geostrophic_monthly_mean.nc");
stokesMonthlyOutputTemplate = Template("D:/SKIM_Paper_Data/SKIM/processed_data/monthly_means/stokes_monthly_mean/${YYYY}${MM}_stokes_monthly_mean.nc");
windMonthlyOutputTemplate = Template("D:/SKIM_Paper_Data/SKIM/processed_data/monthly_means/wavewatch_wnd_monthly_means/${YYYY}${MM}_wind_monthly_mean.nc");
hsMonthlyOutputTemplate = Template("D:/SKIM_Paper_Data/SKIM/processed_data/monthly_means/wavewatch_Hs_monthly_means/${YYYY}${MM}_Hs_monthly_mean.nc");
# Not used?
# mixedLayerDepthPath = "/home/rr/data/mixed_layer_depth_ifremer/mld_DT02_c1m_reg2.0.nc";
# shelfMaskPath = "/home/rr/Files/Tasks/20180914_SKIM/data/shelf_mask_1000.0_0.0.nc";
# mixedLayerDepth = Dataset(mixedLayerDepthPath, 'r').variables["mld"][:];
shelfMask = None; #np.flipud(Dataset(shelfMaskPath, 'r').variables["shelf_mask"][:]);


def apply_mask(data, mask=None, ilatRange=None, ilonRange=None):
    d = data.copy();
    if mask is not None:
        d.mask = d.mask | (mask!=1);

    if ilatRange != None and ilonRange != None:
        d = d[ilatRange[0]:ilatRange[1], ilonRange[0]:ilonRange[1]];
    return d;

def calc_time_index(iYear, iMonth):
    return 12*iYear + iMonth;

def create_netCDF(filename, variableNames, latRes, lonRes):
    ncout = Dataset(filename, 'w');
    ncout.createDimension('lat', ilatRange[1]-ilatRange[0]);
    ncout.createDimension('lon', ilonRange[1]-ilonRange[0]);

    latVar = ncout.createVariable("lat", np.dtype(float), ("lat",) );
    latVar[:] = np.arange(-90+(latRes/2.0), 90+(latRes/2.0), latRes)[::-1];
    latVar.units = "degrees North";
    latVar.long_name = "Latitude";
    latVar.valid_min = -90;
    latVar.valid_max = 90;

    lonVar = ncout.createVariable("lon", np.dtype(float), ("lon",) );
    lonVar[:] = np.arange(-180+(lonRes/2.0), 180+(lonRes/2.0), lonRes)#[::-1];
    lonVar.units = "degrees East";
    lonVar.long_name = "Longitude";
    lonVar.valid_min = -180;
    lonVar.valid_max = 180;

    for varName in variableNames:
        var = ncout.createVariable(varName, np.dtype(float), ("lat", "lon"));
        var.valid_min = -10000.0;
        var.valid_max = 10000.0;

    return ncout;

def process_ekman_month(year, month,time_cor = 2):
    monthStr = format(month+1, "02d");
    resolution = 0.25;

    ekmanU = np.ma.empty( shape=(calendar.monthrange(year,month+1)[1], int(180.0/resolution), int(360.0/resolution)) );
    ekmanU_err = np.ma.empty( shape=(calendar.monthrange(year,month+1)[1], int(180.0/resolution), int(360.0/resolution)) );
    ekmanV = np.ma.empty( shape=(calendar.monthrange(year,month+1)[1], int(180.0/resolution), int(360.0/resolution)) );
    ekmanV_err = np.ma.empty( shape=(calendar.monthrange(year,month+1)[1], int(180.0/resolution), int(360.0/resolution)) );
    ekmanCurrent = np.ma.empty( shape=(calendar.monthrange(year,month+1)[1], int(180.0/resolution), int(360.0/resolution)) );
    ekmanCurrent_err = np.ma.empty( shape=(calendar.monthrange(year,month+1)[1], int(180.0/resolution), int(360.0/resolution)) );
    #for each day in the month, load the input netCDF file and sum the results
    for day in range(0, calendar.monthrange(year,month+1)[1]):
        dayStr = format(day+1, "02d");
        inputData = Dataset(ekmanInputTemplate.safe_substitute(YYYY=str(year), MM=monthStr, DD=dayStr));
        u = np.flipud(inputData.variables["eastward_ekman_current_velocity"][0,:,:]);
        u_err =  np.abs(np.flipud(inputData.variables["eastward_ekman_current_velocity_error"][0,:,:])); # Need np.abs as some uncertainties are negative...
        v = np.flipud(inputData.variables["northward_ekman_current_velocity"][0,:,:]);
        v_err = np.abs(np.flipud(inputData.variables["northward_ekman_current_velocity_error"][0,:,:])); # Need np.abs as some uncertainties are negative...
        ekmanU[day,:,:] = u;
        ekmanU_err[day,:,:] = u_err;
        ekmanV[day,:,:] = v;
        ekmanV_err[day,:,:] = v_err;
        ecur = np.sqrt(u**2 + v**2);
        ekmanCurrent[day,:,:] = ecur
        comb_err = ((np.sqrt((((u_err/u)*2)*u)**2 + (((v_err/v)*2)*v)**2)/ecur)*0.5)*ecur # Convert to percentage error, multiply by two to propagate the squared function, convert back to absolute error for the addition where we assume they are independent so combine in quadrature
        ekmanCurrent_err[day,:,:] = comb_err
        inputData.close()
    #calculate mean, do any final processing
    ekmanUmean = np.ma.mean(ekmanU, axis=0);
    ekmanUmeanMasked = apply_mask(ekmanUmean, shelfMask, ilatRange, ilonRange);

    ekmanUmean_err = np.ma.mean(ekmanU_err,axis=0) / np.sqrt(calendar.monthrange(year,month+1)[1]/time_cor)
    ekmanUmeanMasked_err = apply_mask(ekmanUmean_err, shelfMask, ilatRange, ilonRange);

    ekmanVmean = np.ma.mean(ekmanV, axis=0);
    ekmanVmeanMasked = apply_mask(ekmanVmean, shelfMask, ilatRange, ilonRange);

    ekmanVmean_err = np.ma.mean(ekmanV_err,axis=0) / np.sqrt(calendar.monthrange(year,month+1)[1]/time_cor)
    ekmanVmeanMasked_err = apply_mask(ekmanVmean_err, shelfMask, ilatRange, ilonRange);

    ekmanCurrentMean = np.ma.mean(ekmanCurrent, axis=0);
    ekmanCurrentMeanMasked = apply_mask(ekmanCurrentMean, shelfMask, ilatRange, ilonRange);
    ### Code below does a more complex propagation (taking into account correlated periods, and uncorrelated periods) which gets the same as the mean with the sqrt...
    # ekmancurrent_temp = np.ma.empty( shape=(int(np.ceil(calendar.monthrange(year,month+1)[1]/time_cor)), int(180.0/resolution), int(360.0/resolution))) ;
    # time_val = np.arange(0,calendar.monthrange(year,month+1)[1],time_cor)
    # for i in range(0,len(time_val)-1):
    #     ekmancurrent_temp[i,:,:] = np.ma.mean(ekmanCurrent_err[time_val[i]:time_val[i+1],:,:],axis=0)
    # ekmancurrent_temp[-1,:,:]= np.ma.mean(ekmanCurrent_err[time_val[i+1]:time_val[-1],:,:],axis=0)
    # ekmancurrent_temp = np.sqrt(np.ma.sum(ekmancurrent_temp**2,axis=0))/(calendar.monthrange(year,month+1)[1]/time_cor)
    # ekmanCurrentMean_err = ekmancurrent_temp
    ekmanCurrentMean_err = np.ma.mean(ekmanCurrent_err,axis=0) / np.sqrt(calendar.monthrange(year,month+1)[1]/time_cor)
    ekmanCurrentMeanMasked_err = apply_mask(ekmanCurrentMean_err, shelfMask, ilatRange, ilonRange);

    #write to netCDF
    ekmanMonthlyOutputPath = ekmanMonthlyOutputTemplate.safe_substitute(YYYY=str(year), MM=monthStr);
    ncout = create_netCDF(ekmanMonthlyOutputPath, ["ekmanU",'ekmanUerr',"ekmanV",'ekmanVerr',"ekmanCurrent",'ekmanCurrenterr'], latRes=resolution, lonRes=resolution);

    ncout.variables["ekmanU"][:] = ekmanUmeanMasked;
    ncout.variables["ekmanU"].long_name = "Mean monthly Ekman surface current velocity (eastward component)";
    ncout.variables["ekmanU"].units = "m s-1";

    ncout.variables["ekmanUerr"][:] = ekmanUmeanMasked_err;
    ncout.variables["ekmanUerr"].long_name = "Mean monthly Ekman surface current velocity (eastward component) uncertainty";
    ncout.variables["ekmanUerr"].units = "m s-1";
    ncout.variables["ekmanUerr"].time_correlation = 'Assumed uncertainties are correlated for ' +str(time_cor) + ' days';

    ncout.variables["ekmanV"][:] = ekmanVmeanMasked;
    ncout.variables["ekmanV"].long_name = "Mean monthly Ekman surface current velocity (northward component)";
    ncout.variables["ekmanV"].units = "m s-1";

    ncout.variables["ekmanVerr"][:] = ekmanVmeanMasked_err;
    ncout.variables["ekmanVerr"].long_name = "Mean monthly Ekman surface current velocity (northward component) uncertainty";
    ncout.variables["ekmanVerr"].units = "m s-1";
    ncout.variables["ekmanVerr"].time_correlation = 'Assumed uncertainties are correlated for ' +str(time_cor) + ' days';

    ncout.variables["ekmanCurrent"][:] = ekmanCurrentMeanMasked;
    ncout.variables["ekmanCurrent"].long_name = "Mean monthly Ekman surface current velocity";
    ncout.variables["ekmanCurrent"].units = "m s-1";

    ncout.variables["ekmanCurrenterr"][:] = ekmanCurrentMeanMasked_err;
    ncout.variables["ekmanCurrenterr"].long_name = "Mean monthly Ekman surface current velocity uncertainty";
    ncout.variables["ekmanCurrenterr"].units = "m s-1";
    ncout.variables["ekmanCurrenterr"].time_correlation = 'Assumed uncertainties are correlated for ' +str(time_cor) + ' days';
    ncout.close();


def process_ekman_15m_month(year, month):
    monthStr = format(month+1, "02d");
    resolution = 0.25;

    ekmanU = np.ma.empty( shape=(calendar.monthrange(year,month+1)[1], int(180.0/resolution), int(360.0/resolution)) );
    ekmanV = np.ma.empty( shape=(calendar.monthrange(year,month+1)[1], int(180.0/resolution), int(360.0/resolution)) );
    ekmanCurrent = np.ma.empty( shape=(calendar.monthrange(year,month+1)[1], int(180.0/resolution), int(360.0/resolution)) );
    #for each day in the month, load the input netCDF file and sum the results
    for day in range(0, calendar.monthrange(year,month+1)[1]):
        dayStr = format(day+1, "02d");
        inputData = Dataset(ekman15mInputTemplate.safe_substitute(YYYY=str(year), MM=monthStr, DD=dayStr));
        u = np.flipud(inputData.variables["eastward_ekman_current_velocity"][0,:,:]);
        v = np.flipud(inputData.variables["northward_ekman_current_velocity"][0,:,:]);
        ekmanU[day,:,:] = u;
        ekmanV[day,:,:] = v;
        ekmanCurrent[day,:,:] = np.sqrt(u**2 + v**2);

    #calculate mean, do any final processing
    ekmanUmean = np.ma.mean(ekmanU, axis=0);
    ekmanUmeanMasked = apply_mask(ekmanUmean, shelfMask, ilatRange, ilonRange);
    ekmanVmean = np.ma.mean(ekmanV, axis=0);
    ekmanVmeanMasked = apply_mask(ekmanVmean, shelfMask, ilatRange, ilonRange);
    ekmanCurrentMean = np.ma.mean(ekmanCurrent, axis=0);
    ekmanCurrentMeanMasked = apply_mask(ekmanCurrentMean, shelfMask, ilatRange, ilonRange);

    #write to netCDF
    ekman15mMonthlyOutputPath = ekman15mMonthlyOutputTemplate.safe_substitute(YYYY=str(year), MM=monthStr);
    ncout = create_netCDF(ekman15mMonthlyOutputPath, ["ekmanU", "ekmanV", "ekmanCurrent"], latRes=resolution, lonRes=resolution);

    ncout.variables["ekmanU"][:] = ekmanUmeanMasked;
    ncout.variables["ekmanU"].long_name = "Mean monthly Ekman 15m current velocity (eastward component)";
    ncout.variables["ekmanU"].units = "m s-1";

    ncout.variables["ekmanV"][:] = ekmanVmeanMasked;
    ncout.variables["ekmanV"].long_name = "Mean monthly Ekman 15m current velocity (westward component)";
    ncout.variables["ekmanV"].units = "m s-1";

    ncout.variables["ekmanCurrent"][:] = ekmanCurrentMeanMasked;
    ncout.variables["ekmanCurrent"].long_name = "Mean monthly Ekman 15m current velocity";
    ncout.variables["ekmanCurrent"].units = "m s-1";
    ncout.close();


def process_geostrophic_month(year, month,time_cor=2):
    monthStr = format(month+1, "02d");
    resolution = 0.25;

    geostrophicU = np.ma.empty( shape=(calendar.monthrange(year,month+1)[1], int(180.0/resolution), int(360.0/resolution)) );
    geostrophicU_err = np.ma.empty( shape=(calendar.monthrange(year,month+1)[1], int(180.0/resolution), int(360.0/resolution)) );
    geostrophicV = np.ma.empty( shape=(calendar.monthrange(year,month+1)[1], int(180.0/resolution), int(360.0/resolution)) );
    geostrophicV_err = np.ma.empty( shape=(calendar.monthrange(year,month+1)[1], int(180.0/resolution), int(360.0/resolution)) );
    geostrophicCurrent = np.ma.empty( shape=(calendar.monthrange(year,month+1)[1], int(180.0/resolution), int(360.0/resolution)) );
    geostrophicCurrent_err = np.ma.empty( shape=(calendar.monthrange(year,month+1)[1], int(180.0/resolution), int(360.0/resolution)) );
    #for each day in the month, load the input netCDF file and sum the results
    for day in range(0, calendar.monthrange(year,month+1)[1]):
        dayStr = format(day+1, "02d");
        inputData = Dataset(geostrophicInputTemplate.safe_substitute(YYYY=str(year), MM=monthStr, DD=dayStr));
        u = np.flipud(inputData.variables["eastward_geostrophic_current_velocity"][0,:,:]);
        u_err = np.abs(np.flipud(inputData.variables["eastward_geostrophic_current_velocity_error"][0,:,:]));
        v = np.flipud(inputData.variables["northward_geostrophic_current_velocity"][0,:,:]);
        v_err = np.abs(np.flipud(inputData.variables["northward_geostrophic_current_velocity_error"][0,:,:]));
        inputData.close()
        geostrophicU[day,:,:] = u;
        geostrophicU_err[day,:,:] = u_err;
        geostrophicV[day,:,:] = v;
        geostrophicV_err[day,:,:] = v_err;
        ecur = np.sqrt(u**2 + v**2);
        geostrophicCurrent[day,:,:] = ecur
        comb_err = ((np.sqrt((((u_err/u)*2)*u)**2 + (((v_err/v)*2)*v)**2)/ecur)*0.5)*ecur # Convert to percentage error, multiply by two to propagate the squared function, convert back to absolute error for the addition where we assume they are independent so combine in quadrature
        geostrophicCurrent_err[day,:,:] = comb_err

    #calculate mean, do any final processing

    geostrophicUmean = np.ma.mean(geostrophicU, axis=0);
    geostrophicUmeanMasked = apply_mask(geostrophicUmean, shelfMask, ilatRange, ilonRange);

    geostrophicUmean_err = np.ma.mean(geostrophicU_err,axis=0) / np.sqrt(calendar.monthrange(year,month+1)[1]/time_cor)
    geostrophicUmeanMasked_err = apply_mask(geostrophicUmean_err, shelfMask, ilatRange, ilonRange);

    geostrophicVmean = np.ma.mean(geostrophicV, axis=0);
    geostrophicVmeanMasked = apply_mask(geostrophicVmean, shelfMask, ilatRange, ilonRange);

    geostrophicVmean_err = np.ma.mean(geostrophicV_err,axis=0) / np.sqrt(calendar.monthrange(year,month+1)[1]/time_cor)
    geostrophicVmeanMasked_err = apply_mask(geostrophicVmean_err, shelfMask, ilatRange, ilonRange);

    geostrophicCurrentMean = np.ma.mean(geostrophicCurrent, axis=0);
    geostrophicCurrentMeanMasked = apply_mask(geostrophicCurrentMean, shelfMask, ilatRange, ilonRange);

    geostrophicCurrentMean_err = np.ma.mean(geostrophicCurrent_err,axis=0) / np.sqrt(calendar.monthrange(year,month+1)[1]/time_cor)
    geostrophicCurrentMeanMasked_err = apply_mask(geostrophicCurrentMean_err, shelfMask, ilatRange, ilonRange);

    #write to netCDF
    geostrophicMonthlyOutputPath = geostrophicMonthlyOutputTemplate.safe_substitute(YYYY=str(year), MM=monthStr);
    ncout = create_netCDF(geostrophicMonthlyOutputPath, ["geostrophicU","geostrophicUerr", "geostrophicV","geostrophicVerr", "geostrophicCurrent","geostrophicCurrenterr"], latRes=resolution, lonRes=resolution);

    ncout.variables["geostrophicU"][:] = geostrophicUmeanMasked;
    ncout.variables["geostrophicU"].long_name = "Mean monthly geostrophic surface current velocity (eastward component)";
    ncout.variables["geostrophicU"].units = "m s-1";

    ncout.variables["geostrophicUerr"][:] = geostrophicUmeanMasked_err;
    ncout.variables["geostrophicUerr"].long_name = "Mean monthly geostrophic surface current velocity (eastward component) uncertainty";
    ncout.variables["geostrophicUerr"].units = "m s-1";
    ncout.variables["geostrophicUerr"].time_correlation = 'Assumed uncertainties are correlated for ' +str(time_cor) + ' days';

    ncout.variables["geostrophicV"][:] = geostrophicVmeanMasked;
    ncout.variables["geostrophicV"].long_name = "Mean monthly geostrophic surface current velocity (northward component)";
    ncout.variables["geostrophicV"].units = "m s-1";

    ncout.variables["geostrophicVerr"][:] = geostrophicVmeanMasked_err;
    ncout.variables["geostrophicVerr"].long_name = "Mean monthly geostrophic surface current velocity (northward component) uncertainty";
    ncout.variables["geostrophicVerr"].units = "m s-1";
    ncout.variables["geostrophicVerr"].time_correlation = 'Assumed uncertainties are correlated for ' +str(time_cor) + ' days';

    ncout.variables["geostrophicCurrent"][:] = geostrophicCurrentMeanMasked;
    ncout.variables["geostrophicCurrent"].long_name = "Mean monthly geostrophic surface current velocity";
    ncout.variables["geostrophicCurrent"].units = "m s-1";

    ncout.variables["geostrophicCurrenterr"][:] = geostrophicCurrentMeanMasked_err;
    ncout.variables["geostrophicCurrenterr"].long_name = "Mean monthly geostrophic surface current velocity uncertainty";
    ncout.variables["geostrophicCurrenterr"].units = "m s-1";
    ncout.variables["geostrophicCurrenterr"].time_correlation = 'Assumed uncertainties are correlated for ' +str(time_cor) + ' days';
    ncout.close();


def process_stokes_month(year, month,percent_err = 0.2,time_cor=2):
    monthStr = format(month+1, "02d");
    resolution = 0.5;

    stokesU = np.ma.empty( shape=(calendar.monthrange(year,month+1)[1]*8, int(180.0/resolution), int(360.0/resolution)) );
    stokesU[:,:,:] = np.nan
    stokesV = np.ma.empty( shape=(calendar.monthrange(year,month+1)[1]*8, int(180.0/resolution), int(360.0/resolution)) );
    stokesV[:,:,:] = np.nan
    stokesCurrent = np.ma.empty( shape=(calendar.monthrange(year,month+1)[1]*8, int(180.0/resolution), int(360.0/resolution)) );
    stokesCurrent[:,:,:] = np.nan
    inputData = Dataset(stokesInputTemplate.safe_substitute(YYYY=str(year), MM=monthStr));
    u = inputData.variables["uuss"][:,:,:];

    v = inputData.variables["vuss"][:,:,:];

    stokesU[:,24:341,:] = u;
    u_err = np.abs(stokesU) * percent_err
    stokesV[:,24:341,:] = v;
    v_err = np.abs(stokesV) * percent_err
    ecur = np.sqrt(stokesU**2 + stokesV**2);
    comb_err = ((np.sqrt((((u_err/stokesU)*2)*stokesU)**2 + (((v_err/stokesV)*2)*stokesV)**2)/ecur)*0.5)*ecur # Convert to percentage error, multiply by two to propagate the squared function, convert back to absolute error for the addition where we assume they are independent so combine in quadrature

    stokesCurrent[:,:,:] = np.sqrt(stokesU**2 + stokesV**2);

    #calculate mean, do any final processing
    stokesUmean = np.ma.mean(stokesU, axis=0);
    stokesUmean = np.flipud(np.repeat(np.repeat(stokesUmean, 2, axis=0), 2, axis=1));
    stokesUmeanMasked = apply_mask(stokesUmean, shelfMask, ilatRange, ilonRange);

    stokesUmeanerr = np.ma.mean(u_err, axis=0) / np.sqrt(calendar.monthrange(year,month+1)[1]/time_cor);
    stokesUmeanerr = np.flipud(np.repeat(np.repeat(stokesUmeanerr, 2, axis=0), 2, axis=1));
    stokesUmeanMaskederr = apply_mask(stokesUmeanerr, shelfMask, ilatRange, ilonRange);

    stokesVmean = np.ma.mean(stokesV, axis=0);
    stokesVmean = np.flipud(np.repeat(np.repeat(stokesVmean, 2, axis=0), 2, axis=1));
    stokesVmeanMasked = apply_mask(stokesVmean, shelfMask, ilatRange, ilonRange);

    stokesVmeanerr = np.ma.mean(v_err, axis=0) / np.sqrt(calendar.monthrange(year,month+1)[1]/time_cor);
    stokesVmeanerr = np.flipud(np.repeat(np.repeat(stokesVmeanerr, 2, axis=0), 2, axis=1));
    stokesVmeanMaskederr = apply_mask(stokesVmeanerr, shelfMask, ilatRange, ilonRange);

    stokesCurrentMean = np.ma.mean(stokesCurrent, axis=0);
    stokesCurrentMean = np.flipud(np.repeat(np.repeat(stokesCurrentMean, 2, axis=0), 2, axis=1));
    stokesCurrentMeanMasked = apply_mask(stokesCurrentMean, shelfMask, ilatRange, ilonRange);

    stokesCurrentmeanerr = np.ma.mean(comb_err, axis=0) / np.sqrt(calendar.monthrange(year,month+1)[1]/time_cor);
    stokesCurrentmeanerr = np.flipud(np.repeat(np.repeat(stokesCurrentmeanerr, 2, axis=0), 2, axis=1));
    stokesCurrentmeanMaskederr = apply_mask(stokesCurrentmeanerr, shelfMask, ilatRange, ilonRange);

    #write to netCDF
    stokesMonthlyOutputPath = stokesMonthlyOutputTemplate.safe_substitute(YYYY=str(year), MM=monthStr);
    ncout = create_netCDF(stokesMonthlyOutputPath, ["stokesU",'stokesUerr', "stokesV",'stokesVerr', "stokesCurrent",'stokesCurrenterr'], latRes=resolution/2.0, lonRes=resolution/2.0);

    ncout.variables["stokesU"][:] = stokesUmeanMasked;
    ncout.variables["stokesU"].long_name = "Mean monthly stokes surface current velocity (eastward component)";
    ncout.variables["stokesU"].units = "m s-1";

    ncout.variables["stokesUerr"][:] = stokesUmeanMaskederr;
    ncout.variables["stokesUerr"].long_name = "Mean monthly stokes surface current velocity (eastward component) uncertainty";
    ncout.variables["stokesUerr"].units = "m s-1";
    ncout.variables["stokesUerr"].time_correlation = 'Assumed uncertainties are correlated for ' +str(time_cor) + ' days';
    ncout.variables["stokesUerr"].fixed_uncertainty = 'Data produced assuming a 20 % uncertainty on the 3 hourly observations from: https://doi.org/10.1016/j.ocemod.2012.12.001';

    ncout.variables["stokesV"][:] = stokesVmeanMasked;
    ncout.variables["stokesV"].long_name = "Mean monthly stokes surface current velocity (northward component)";
    ncout.variables["stokesV"].units = "m s-1";

    ncout.variables["stokesVerr"][:] = stokesVmeanMaskederr;
    ncout.variables["stokesVerr"].long_name = "Mean monthly stokes surface current velocity (northward component) uncertainty";
    ncout.variables["stokesVerr"].units = "m s-1";
    ncout.variables["stokesVerr"].time_correlation = 'Assumed uncertainties are correlated for ' +str(time_cor) + ' days';
    ncout.variables["stokesVerr"].fixed_uncertainty = 'Data produced assuming a 20 % uncertainty on the 3 hourly observations from: https://doi.org/10.1016/j.ocemod.2012.12.001';

    ncout.variables["stokesCurrent"][:] = stokesCurrentMeanMasked;
    ncout.variables["stokesCurrent"].long_name = "Mean monthly stokes surface current velocity";
    ncout.variables["stokesCurrent"].units = "m s-1";

    ncout.variables["stokesCurrenterr"][:] = stokesCurrentmeanMaskederr;
    ncout.variables["stokesCurrenterr"].long_name = "Mean monthly stokes surface current velocity uncertainty";
    ncout.variables["stokesCurrenterr"].units = "m s-1";
    ncout.variables["stokesCurrenterr"].time_correlation = 'Assumed uncertainties are correlated for ' +str(time_cor) + ' days';
    ncout.variables["stokesCurrenterr"].fixed_uncertainty = 'Data produced assuming a 20 % uncertainty on the 3 hourly observations from: https://doi.org/10.1016/j.ocemod.2012.12.001';
    ncout.close();

def process_wind_month(year, month):
    monthStr = format(month+1, "02d");
    resolution = 0.5;

    stokesU = np.ma.empty( shape=(calendar.monthrange(year,month+1)[1]*8, int(180.0/resolution), int(360.0/resolution)) );
    stokesV = np.ma.empty( shape=(calendar.monthrange(year,month+1)[1]*8, int(180.0/resolution), int(360.0/resolution)) );
    stokesCurrent = np.ma.empty( shape=(calendar.monthrange(year,month+1)[1]*8, int(180.0/resolution), int(360.0/resolution)) );
    #for each day in the month, load the input netCDF file and sum the results
    # for day in range(0, calendar.monthrange(year,month+1)[1]):
    #     dayStr = format(day+1, "02d");
    inputData = Dataset(windInputTemplate.safe_substitute(YYYY=str(year), MM=monthStr));
    u = inputData.variables["uwnd"][:,:,:];
    v = inputData.variables["vwnd"][:,:,:];
    stokesU[:,24:341,:] = u;
    stokesV[:,24:341,:] = v;
    stokesCurrent[:,24:341,:] = np.sqrt(u**2 + v**2);

    #calculate mean, do any final processing
    stokesUmean = np.ma.mean(stokesU, axis=0);
    stokesUmean = np.flipud(np.repeat(np.repeat(stokesUmean, 2, axis=0), 2, axis=1));
    stokesUmeanMasked = apply_mask(stokesUmean, shelfMask, ilatRange, ilonRange);
    stokesVmean = np.ma.mean(stokesV, axis=0);
    stokesVmean = np.flipud(np.repeat(np.repeat(stokesVmean, 2, axis=0), 2, axis=1));
    stokesVmeanMasked = apply_mask(stokesVmean, shelfMask, ilatRange, ilonRange);
    stokesCurrentMean = np.ma.mean(stokesCurrent, axis=0);
    stokesCurrentMean = np.flipud(np.repeat(np.repeat(stokesCurrentMean, 2, axis=0), 2, axis=1));
    stokesCurrentMeanMasked = apply_mask(stokesCurrentMean, shelfMask, ilatRange, ilonRange);

    #write to netCDF
    stokesMonthlyOutputPath = windMonthlyOutputTemplate.safe_substitute(YYYY=str(year), MM=monthStr);
    ncout = create_netCDF(stokesMonthlyOutputPath, ["uwnd_mean", "vwind_mean", "windspeed_mean"], latRes=resolution/2.0, lonRes=resolution/2.0);

    ncout.variables["uwnd_mean"][:] = stokesUmeanMasked;
    ncout.variables["uwnd_mean"].long_name = "Mean monthly wind velocity (eastward component)";
    ncout.variables["uwnd_mean"].units = "m s-1";

    ncout.variables["vwind_mean"][:] = stokesVmeanMasked;
    ncout.variables["vwind_mean"].long_name = "Mean monthly wind velocity (westward component)";
    ncout.variables["vwind_mean"].units = "m s-1";

    ncout.variables["windspeed_mean"][:] = stokesCurrentMeanMasked;
    ncout.variables["windspeed_mean"].long_name = "Mean monthly wind velocity";
    ncout.variables["windspeed_mean"].units = "m s-1";
    ncout.close();

def process_Hs_month(year, month):
    monthStr = format(month+1, "02d");
    resolution = 0.5;

    stokesU = np.ma.empty( shape=(calendar.monthrange(year,month+1)[1]*8, int(180.0/resolution), int(360.0/resolution)) );

    #for each day in the month, load the input netCDF file and sum the results
    # for day in range(0, calendar.monthrange(year,month+1)[1]):
    #     dayStr = format(day+1, "02d");
    inputData = Dataset(hsInputTemplate.safe_substitute(YYYY=str(year), MM=monthStr));
    u= inputData.variables["hs"][:,:,:];
    stokesU[:,24:341,:] = u;


    #calculate mean, do any final processing
    stokesUmean = np.ma.mean(stokesU, axis=0);
    stokesUmean = np.flipud(np.repeat(np.repeat(stokesUmean, 2, axis=0), 2, axis=1));
    stokesUmeanMasked = apply_mask(stokesUmean, shelfMask, ilatRange, ilonRange);

    #write to netCDF
    stokesMonthlyOutputPath = hsMonthlyOutputTemplate.safe_substitute(YYYY=str(year), MM=monthStr);
    ncout = create_netCDF(stokesMonthlyOutputPath, ["hs_mean"], latRes=resolution/2.0, lonRes=resolution/2.0);

    ncout.variables["hs_mean"][:] = stokesUmeanMasked;
    ncout.variables["hs_mean"].long_name = "Mean monthly significant wave height";
    ncout.variables["hs_mean"].units = "m";

    ncout.close();


#years = [2002];
years = range(1993, 2017);
months = range(0, 12);
ilatRange = (0, 180*4);
ilonRange = (0, 360*4);

for iYear, year in enumerate(years):
    for month in months:
        if processMonthly:
            print "Processing monthly:", year, month;
            if processEkman:
                process_ekman_month(year, month);
            if processEkman15m:
                process_ekman_15m_month(year, month);
            if processGeostrophic:
                process_geostrophic_month(year, month);
            if processStokes:
                process_stokes_month(year, month);
            if processwind:
                process_wind_month(year,month)
            if processHs:
                process_Hs_month(year,month)
        if processDaily:
            for day in range(0, calendar.monthrange(year,month+1)[1]):
                print year, month, day;
                pass;

#        ncout = create_ncdf("data/monthly_mean_currents_"+str(year)+monthStr+".nc", 12*len(years), 1.0, 0.25, 0.25, ilatRange, ilonRange);
#
#        ekmanMonthMean = np.ma.empty( (calendar.monthrange(year,month+1)[1], 720, 1440) );
#        ekmanDirMonthMean = np.ma.empty( (calendar.monthrange(year,month+1)[1], 720, 1440) );
#        ekmanUMonthMean = np.ma.empty( (calendar.monthrange(year,month+1)[1], 720, 1440) );
#        ekmanVMonthMean = np.ma.empty( (calendar.monthrange(year,month+1)[1], 720, 1440) );
#        ekman15MonthMean = np.ma.empty( (calendar.monthrange(year,month+1)[1], 720, 1440) );
#        ekman15DirMonthMean = np.ma.empty( (calendar.monthrange(year,month+1)[1], 720, 1440) );
#        ekman15UMonthMean = np.ma.empty( (calendar.monthrange(year,month+1)[1], 720, 1440) );
#        ekman15VMonthMean = np.ma.empty( (calendar.monthrange(year,month+1)[1], 720, 1440) );
#        geostrophicMonthMean = np.ma.empty( (calendar.monthrange(year,month+1)[1], 720, 1440) );
#        geostrophicDirMonthMean = np.ma.empty( (calendar.monthrange(year,month+1)[1], 720, 1440) );
#        geostrophicUMonthMean = np.ma.empty( (calendar.monthrange(year,month+1)[1], 720, 1440) );
#        geostrophicVMonthMean = np.ma.empty( (calendar.monthrange(year,month+1)[1], 720, 1440) );
#        stokesMonthMean = np.ma.empty( (calendar.monthrange(year,month+1)[1], 720, 1440) );
#        stokesDirMonthMean = np.ma.empty( (calendar.monthrange(year,month+1)[1], 720, 1440) );
#        stokesUMonthMean = np.ma.empty( (calendar.monthrange(year,month+1)[1], 720, 1440) );
#        stokesVMonthMean = np.ma.empty( (calendar.monthrange(year,month+1)[1], 720, 1440) );
#        totalMonthMean = np.ma.empty( (calendar.monthrange(year,month+1)[1], 720, 1440) );
#        totalDirMonthMean = np.ma.empty( (calendar.monthrange(year,month+1)[1], 720, 1440) );
#        totalUMonthMean = np.ma.empty( (calendar.monthrange(year,month+1)[1], 720, 1440) );
#        totalVMonthMean = np.ma.empty( (calendar.monthrange(year,month+1)[1], 720, 1440) );
#        #for each day in the month, accumulate the data into a single 3d array
#        for day in range(0, calendar.monthrange(year,month+1)[1]):
#            print year, month, day;
#            dayStr = format(day+1, "02d");
#            ekmanNc = Dataset(ekmanTemplate.safe_substitute(YYYY=str(year), MM=monthStr, DD=dayStr), 'r');
#            ekman15Nc = Dataset(ekman15Template.safe_substitute(YYYY=str(year), MM=monthStr, DD=dayStr), 'r');
#            geostrophicNc = Dataset(geostrophicTemplate.safe_substitute(YYYY=str(year), MM=monthStr, DD=dayStr), 'r');
#            stokesNc = Dataset(stokesTemplate.safe_substitute(YYYY=str(year), MM=monthStr, DD=dayStr), 'r');
#
#            ekmanMonthMean[day, :, :] = np.flipud(su.calculate_magnitude(ekmanNc.variables["eastward_ekman_current_velocity"][0,:,:], ekmanNc.variables["northward_ekman_current_velocity"][0,:,:]));
#            ekmanDirMonthMean[day, :, :] = np.flipud(su.calculate_direction(ekmanNc.variables["eastward_ekman_current_velocity"][0,:,:], ekmanNc.variables["northward_ekman_current_velocity"][0,:,:]));
#            ekmanUMonthMean[day, :, :] = np.flipud(ekmanNc.variables["eastward_ekman_current_velocity"][0,:,:]);
#            ekmanVMonthMean[day, :, :] = np.flipud(ekmanNc.variables["northward_ekman_current_velocity"][0,:,:]);
#
#            ekman15MonthMean[day, :, :] = np.flipud(su.calculate_magnitude(ekman15Nc.variables["eastward_ekman_current_velocity"][0,:,:], ekman15Nc.variables["northward_ekman_current_velocity"][0,:,:]));
#            ekman15DirMonthMean[day, :, :] = np.flipud(su.calculate_direction(ekman15Nc.variables["eastward_ekman_current_velocity"][0,:,:], ekman15Nc.variables["northward_ekman_current_velocity"][0,:,:]));
#            ekman15UMonthMean[day, :, :] = np.flipud(ekman15Nc.variables["eastward_ekman_current_velocity"][0,:,:]);
#            ekman15VMonthMean[day, :, :] = np.flipud(ekman15Nc.variables["northward_ekman_current_velocity"][0,:,:]);
#
#            geostrophicMonthMean[day, :, :] = np.flipud(su.calculate_magnitude(geostrophicNc.variables["eastward_geostrophic_current_velocity"][0,:,:], geostrophicNc.variables["northward_geostrophic_current_velocity"][0,:,:]));
#            geostrophicDirMonthMean[day, :, :] = np.flipud(su.calculate_direction(geostrophicNc.variables["eastward_geostrophic_current_velocity"][0,:,:], geostrophicNc.variables["northward_geostrophic_current_velocity"][0,:,:]));
#            geostrophicUMonthMean[day, :, :] = np.flipud(geostrophicNc.variables["eastward_geostrophic_current_velocity"][0,:,:]);
#            geostrophicVMonthMean[day, :, :] = np.flipud(geostrophicNc.variables["northward_geostrophic_current_velocity"][0,:,:]);
#
#            stokesEastStretched =  np.repeat(np.repeat(stokesNc.variables["uuss"][:,:,:], 2, axis=1), 2, axis=2);
#            stokesNorthStretched =  np.repeat(np.repeat(stokesNc.variables["vuss"][:,:,:], 2, axis=1), 2, axis=2);
#            stokesMonthMean[day, :, :] = np.flipud(su.calculate_magnitude(stokesEastStretched, stokesNorthStretched));
#            stokesDirMonthMean[day, :, :] = np.flipud(su.calculate_direction(stokesEastStretched, stokesNorthStretched));
#            stokesUMonthMean[day, :, :] = np.flipud(stokesEastStretched[0,:,:]);
#            stokesVMonthMean[day, :, :] = np.flipud(stokesNorthStretched[0,:,:]);
#
#            totalEastwardCurrentVelocity = ekmanNc.variables["eastward_ekman_current_velocity"][0,:,:] + geostrophicNc.variables["eastward_geostrophic_current_velocity"][0,:,:] + stokesEastStretched[0,:,:]
#            totalNorthwardCurrentVelocity = ekmanNc.variables["northward_ekman_current_velocity"][0,:,:] + geostrophicNc.variables["northward_geostrophic_current_velocity"][0,:,:] + stokesNorthStretched[0,:,:]
#            totalMonthMean[day, :, :] = np.flipud(su.calculate_magnitude(totalEastwardCurrentVelocity, totalNorthwardCurrentVelocity));
#            totalDirMonthMean[day, :, :] = np.flipud(su.calculate_direction(totalEastwardCurrentVelocity, totalNorthwardCurrentVelocity));
#            totalUMonthMean[day, :, :] = np.flipud(totalEastwardCurrentVelocity);
#            totalVMonthMean[day, :, :] = np.flipud(totalNorthwardCurrentVelocity);
#
#
#
#        #Calculate means and write to in netcdf file
#        ekmanMonthMean = np.ma.mean(ekmanMonthMean, axis=0);
#        ekmanMonthMeanMasked = su.apply_mask(ekmanMonthMean, shelfMask, ilatRange, ilonRange);
#        ncout.variables["ekman_velocity_mean"][calc_time_index(iYear, month),:,:] = ekmanMonthMeanMasked;
#        ekmanDirMonthMean = np.ma.mean(ekmanDirMonthMean, axis=0);
#        ekmanDirMonthMeanMasked = su.apply_mask(ekmanDirMonthMean, shelfMask, ilatRange, ilonRange);
#        ncout.variables["ekman_direction_mean"][calc_time_index(iYear, month),:,:] = ekmanDirMonthMeanMasked;
#        ekmanUMonthMean = np.ma.mean(ekmanUMonthMean, axis=0);
#        ekmanUMonthMeanMasked = su.apply_mask(ekmanUMonthMean, shelfMask, ilatRange, ilonRange);
#        ncout.variables["ekman_eastward_velocity_mean"][calc_time_index(iYear, month),:,:] = ekmanUMonthMeanMasked;
#        ekmanVMonthMean = np.ma.mean(ekmanVMonthMean, axis=0);
#        ekmanVMonthMeanMasked = su.apply_mask(ekmanVMonthMean, shelfMask, ilatRange, ilonRange);
#        ncout.variables["ekman_northward_velocity_mean"][calc_time_index(iYear, month),:,:] = ekmanVMonthMeanMasked;
#
#        ekman15MonthMean = np.ma.mean(ekman15MonthMean, axis=0);
#        ekman15MonthMeanMasked = su.apply_mask(ekman15MonthMean, shelfMask, ilatRange, ilonRange);
#        ncout.variables["ekman15_velocity_mean"][calc_time_index(iYear, month),:,:] = ekman15MonthMeanMasked;
#        ekman15DirMonthMean = np.ma.mean(ekman15DirMonthMean, axis=0);
#        ekman15DirMonthMeanMasked = su.apply_mask(ekman15DirMonthMean, shelfMask, ilatRange, ilonRange);
#        ncout.variables["ekman15_direction_mean"][calc_time_index(iYear, month),:,:] = ekman15DirMonthMeanMasked;
#        ekman15UMonthMean = np.ma.mean(ekman15UMonthMean, axis=0);
#        ekman15UMonthMeanMasked = su.apply_mask(ekman15UMonthMean, shelfMask, ilatRange, ilonRange);
#        ncout.variables["ekman15_eastward_velocity_mean"][calc_time_index(iYear, month),:,:] = ekman15UMonthMeanMasked;
#        ekman15VMonthMean = np.ma.mean(ekman15VMonthMean, axis=0);
#        ekman15VMonthMeanMasked = su.apply_mask(ekman15VMonthMean, shelfMask, ilatRange, ilonRange);
#        ncout.variables["ekman15_northward_velocity_mean"][calc_time_index(iYear, month),:,:] = ekman15VMonthMeanMasked;
#
#        geostrophicMonthMean = np.ma.mean(geostrophicMonthMean, axis=0);
#        geostrophicMonthMeanMasked = su.apply_mask(geostrophicMonthMean, shelfMask, ilatRange, ilonRange);
#        ncout.variables["geostrophic_velocity_mean"][calc_time_index(iYear, month),:,:] = geostrophicMonthMeanMasked;
#        geostrophicDirMonthMean = np.ma.mean(geostrophicDirMonthMean, axis=0);
#        geostrophicDirMonthMeanMasked = su.apply_mask(geostrophicDirMonthMean, shelfMask, ilatRange, ilonRange);
#        ncout.variables["geostrophic_direction_mean"][calc_time_index(iYear, month),:,:] = geostrophicDirMonthMeanMasked;
#        geostrophicUMonthMean = np.ma.mean(geostrophicUMonthMean, axis=0);
#        geostrophicUMonthMeanMasked = su.apply_mask(geostrophicUMonthMean, shelfMask, ilatRange, ilonRange);
#        ncout.variables["geostrophic_eastward_velocity_mean"][calc_time_index(iYear, month),:,:] = geostrophicUMonthMeanMasked;
#        geostrophicVMonthMean = np.ma.mean(geostrophicVMonthMean, axis=0);
#        geostrophicVMonthMeanMasked = su.apply_mask(geostrophicVMonthMean, shelfMask, ilatRange, ilonRange);
#        ncout.variables["geostrophic_northward_velocity_mean"][calc_time_index(iYear, month),:,:] = geostrophicVMonthMeanMasked;
#
#        stokesMonthMean = np.ma.mean(stokesMonthMean, axis=0);
#        stokesMonthMeanMasked = su.apply_mask(stokesMonthMean, shelfMask, ilatRange, ilonRange);
#        ncout.variables["stokes_velocity_mean"][calc_time_index(iYear, month),:,:] = np.flipud(stokesMonthMeanMasked);
#        stokesDirMonthMean = np.ma.mean(stokesDirMonthMean, axis=0);
#        stokesDirMonthMeanMasked = su.apply_mask(stokesDirMonthMean, shelfMask, ilatRange, ilonRange);
#        ncout.variables["stokes_direction_mean"][calc_time_index(iYear, month),:,:] = np.flipud(stokesDirMonthMeanMasked);
#        stokesUMonthMean = np.ma.mean(stokesUMonthMean, axis=0);
#        stokesUMonthMeanMasked = su.apply_mask(stokesUMonthMean, shelfMask, ilatRange, ilonRange);
#        ncout.variables["stokes_eastward_velocity_mean"][calc_time_index(iYear, month),:,:] = stokesUMonthMeanMasked;
#        stokesVMonthMean = np.ma.mean(stokesVMonthMean, axis=0);
#        stokesVMonthMeanMasked = su.apply_mask(stokesVMonthMean, shelfMask, ilatRange, ilonRange);
#        ncout.variables["stokes_northward_velocity_mean"][calc_time_index(iYear, month),:,:] = stokesVMonthMeanMasked;
#
#        totalMonthMean = np.ma.mean(totalMonthMean, axis=0);
#        totalMonthMeanMasked = su.apply_mask(totalMonthMean, shelfMask, ilatRange, ilonRange);
#        ncout.variables["total_velocity_mean"][calc_time_index(iYear, month),:,:] = totalMonthMeanMasked;
#        totalDirMonthMean = np.ma.mean(totalDirMonthMean, axis=0);
#        totalDirMonthMeanMasked = su.apply_mask(totalDirMonthMean, shelfMask, ilatRange, ilonRange);
#        ncout.variables["total_direction_mean"][calc_time_index(iYear, month),:,:] = totalDirMonthMeanMasked;
#        totalUMonthMean = np.ma.mean(totalUMonthMean, axis=0);
#        totalUMonthMeanMasked = su.apply_mask(totalUMonthMean, shelfMask, ilatRange, ilonRange);
#        ncout.variables["total_eastward_velocity_mean"][calc_time_index(iYear, month),:,:] = totalUMonthMeanMasked;
#        totalVMonthMean = np.ma.mean(totalVMonthMean, axis=0);
#        totalVMonthMeanMasked = su.apply_mask(totalVMonthMean, shelfMask, ilatRange, ilonRange);
#        ncout.variables["total_northward_velocity_mean"][calc_time_index(iYear, month),:,:] = totalVMonthMeanMasked;
#
#
#ncout.close();
#
#plt.figure(); plt.imshow(ekmanMonthMeanMasked); plt.colorbar(); plt.title("Monthly mean Ekman current magnitude");
#plt.figure(); plt.imshow(ekmanDirMonthMeanMasked); plt.colorbar(); plt.title("Monthly mean Ekman current direction");
#if shelfMask != None:
#    plt.figure(); plt.imshow(shelfMask); plt.colorbar(); plt.title("Shelf mask");


#plt.figure(); plt.imshow(totalMonthMean);
#totalMonthMeanMasked = su.apply_mask(totalMonthMean, shelfMask, ilatRange, ilonRange);
#plt.figure(); plt.imshow(totalMonthMeanMasked);
#plt.imshow(shelfMask);
