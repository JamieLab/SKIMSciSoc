#!/usr/bin/env python3

import pandas as pd
import numpy as np
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import matplotlib.transforms
import scipy.stats
from scipy.optimize import curve_fit

font = {'weight' : 'normal',
        'size'   : 24}
matplotlib.rc('font', **font)
let = list(map(chr, range(97, 123)))

def func(x, a, b,c):
    return a*(x**2) + b*x + c

def bootstrap_linear(x,y):
    r2 = []
    for i in range(0,100):
        a = list(np.random.choice(np.arange(0,len(x)), len(x)-2, replace=False))
        stat = scipy.stats.linregress(x[a],y[a])
        r2.append(stat.rvalue**2)
    r2 = np.array(r2)
    print('Linear bootstrap')
    print(np.mean(r2))
    print(np.std(r2))
    return np.std(r2)

def bootstrap_quad(x,y):
    r2 = []
    for i in range(0,100):
        a = list(np.random.choice(np.arange(0,len(x)), len(x)-2, replace=False))
        popt, pcov = curve_fit(func, x[a],y[a])
        coef,pvalue = scipy.stats.pearsonr(y[a],func(x[a],*popt))
        r2.append(coef**2)
    r2 = np.array(r2)
    print('Quad bootstrap')
    print(np.mean(r2))
    print(np.std(r2))
    return np.std(r2)

def plot_scatter(xdata = None, dataerr = None,ydata=None,ydata_err=None,xlabel = None,ax = None,lett = 'a',size=18):
    ax.errorbar(xdata,ydata,yerr=ydata_err,xerr=dataerr,linestyle='',color='k')
    c = ax.scatter(xdata,ydata,zorder=3)
    ax.set_ylabel('d$\Delta$pCO$_2$/dt ($\mu$atm yr$^{-1}$)')
    ax.set_xlabel(xlabel)
    #plt.colorbar(c,ax=axs[0])
    ax.grid()
    a = np.arange(np.amin(xdata),np.amax(xdata),np.abs(np.amin(xdata)-np.amax(xdata))/100 )
    stat = scipy.stats.linregress(xdata,ydata)
    lin_boot = bootstrap_linear(xdata,ydata)

    popt, pcov = curve_fit(func, xdata,ydata)
    coef,pvalue = scipy.stats.pearsonr(ydata,func(xdata,*popt))
    quad_boot = bootstrap_quad(xdata,ydata)
    print('Linear = ' + str(stat.rvalue**2))
    print('Quad = ' + str(coef**2))
    if stat.rvalue**2+lin_boot > coef**2:
        ax.plot(a,a*stat.slope + stat.intercept,'r--',linewidth=2)
        r = stat.rvalue**2
        pval = stat.pvalue
        boot = lin_boot
    else:
        ax.plot(a,func(a,*popt),'b--',linewidth=2)
        r = coef**2
        pval = pvalue
        boot = quad_boot
    if pval < 0.01:
        ax.text(0.05,0.95,'r$^{2}$ = ' + str(np.round(coef**2,2)) +' $\pm$ ' + str(np.round(quad_boot,2)) +'\np-value < 0.01' + '\nn = ' + str(len(data['EksAbsolute'])),transform=ax.transAxes,va='top',fontsize=22)
    else:
        ax.text(0.05,0.95,'r$^{2}$ = ' + str(np.round(r,2)) +' $\pm$ ' + str(np.round(boot,2)) +'\np-value = ' + str(np.round(pval,2)) + '\nn = ' + str(len(xdata)),transform=ax.transAxes,va='top',fontsize=22)
    if pval < 0.05:
        ax.set_facecolor('lightgray')
    if lett is str():
        ax.text(0.90,0.95,'('+lett+')',transform=ax.transAxes,va='top',fontsize=size,fontweight='bold')
        return lett
    else:
        if lett[0] == 0:
            letter = let[lett[1]]
        else:
            letter = let[lett[0]-1] + let[lett[1]]
        ax.text(0.85,0.95,'('+letter+')',transform=ax.transAxes,va='top',fontsize=size,fontweight='bold')
        lett[1] = lett[1]+1
        if lett[1]>25:
            lett[0] = lett[0]+1
            lett[1] = 0
        return lett

def plot_scatter_grid(data,ax,text,t = [0,0]):

    ax[0].annotate(text, xy=(0, 0.5), xytext=(-ax[0].yaxis.labelpad - 20, 0),
                xycoords=ax[0].yaxis.label, textcoords='offset points',
                size='large', ha='right', va='center',weight='bold')
    ax1 = ax[0]
    t = plot_scatter(xdata=data['EksAbsolute'],dataerr=data['EksAbsoluteStd']/np.sqrt(data['EksAbsolute_n']),ydata =data['dpCO2 trend wide shelf mean'],ydata_err=data['dpCO2 trend wide shelf std']/np.sqrt(data['dpCO2 trend wide shelf_n']), xlabel='Mean shelf \n break Ekman current (ms$^{-1}$)',ax=ax1,lett=t)#
    ax2 = ax[1]
    t = plot_scatter(xdata=np.abs(data['EksAbsolute']),dataerr=data['EksAbsoluteStd']/np.sqrt(data['EksAbsolute_n']),ydata =data['dpCO2 trend wide shelf mean'],ydata_err=data['dpCO2 trend wide shelf std']/np.sqrt(data['dpCO2 trend wide shelf_n']),xlabel='Absolute mean\n shelf break Ekman current ($ms^{-1}$)',ax=ax2,lett=t)#
    ax1 = ax[2]
    t = plot_scatter(xdata=data['total across-shelf'],dataerr=data['total across-shelf SD']/np.sqrt(data['total across-shelf_n']),ydata =data['dpCO2 trend wide shelf mean'],ydata_err=data['dpCO2 trend wide shelf std']/np.sqrt(data['dpCO2 trend wide shelf_n']),xlabel='Mean shelf\n break total current (ms$^{-1}$)',ax=ax1,lett=t)#
    ax2 = ax[3]
    t = plot_scatter(xdata=np.abs(data['total across-shelf']),dataerr=data['total across-shelf SD']/np.sqrt(data['total across-shelf_n']),ydata =data['dpCO2 trend wide shelf mean'],ydata_err=data['dpCO2 trend wide shelf std']/np.sqrt(data['dpCO2 trend wide shelf_n']),xlabel='Absolute mean shelf\n break total current (ms$^{-1}$)',ax=ax2,lett=t)#
    ax1 = ax[4]
    t = plot_scatter(xdata=data['GeoAbsolute'],dataerr=data['GeoAbsoluteStd']/np.sqrt(data['GeoAbsolute_n']),ydata =data['dpCO2 trend wide shelf mean'],ydata_err=data['dpCO2 trend wide shelf std']/np.sqrt(data['dpCO2 trend wide shelf_n']),xlabel='Mean shelf\n break Geostropic current (ms$^{-1}$)',ax=ax1,lett=t)#
    ax1 = ax[5]
    t = plot_scatter(xdata=np.abs(data['GeoAbsolute']),dataerr=data['GeoAbsoluteStd']/np.sqrt(data['GeoAbsolute_n']),ydata =data['dpCO2 trend wide shelf mean'],ydata_err=data['dpCO2 trend wide shelf std']/np.sqrt(data['dpCO2 trend wide shelf_n']),xlabel='Absolute mean shelf\n break Geostrophic current (ms$^{-1}$)',ax=ax1,lett=t)#
    ax1 = ax[6]
    t = plot_scatter(xdata=data['k'],dataerr=data['kSD']/data['k_n'],ydata =data['dpCO2 trend wide shelf mean'],ydata_err=data['dpCO2 trend wide shelf std']/np.sqrt(data['dpCO2 trend wide shelf_n']),xlabel='Mean air-sea CO$_2$ gas\ntransfer velocity (cm $hr^{-1}$)',ax=ax1,lett=t)
    ax1 = ax[7]
    t = plot_scatter(xdata=data['proportionEks']*100,dataerr=data['proportionEksSD']*100/np.sqrt(data['proportionEks_n']),ydata =data['dpCO2 trend wide shelf mean'],ydata_err=data['dpCO2 trend wide shelf std']/np.sqrt(data['dpCO2 trend wide shelf_n']),xlabel='Ekman proportion of\ntotal shelf break current (%)',ax=ax1,lett=t)
    ax1 = ax[8]
    t = plot_scatter(xdata=data['proportionGeo']*100,dataerr=data['proportionGeoSD']*100/np.sqrt(data['proportionGeo_n']),ydata =data['dpCO2 trend wide shelf mean'],ydata_err=data['dpCO2 trend wide shelf std']/np.sqrt(data['dpCO2 trend wide shelf_n']),xlabel='Geostrophic proportion of\ntotal shelf break current (%)',ax=ax1,lett=t)
    ax1 = ax[9]
    t = plot_scatter(xdata=data['proportionStokes']*100,dataerr=data['proportionStokesSD']*100/np.sqrt(data['proportionStokes_n']),ydata =data['dpCO2 trend wide shelf mean'],ydata_err=data['dpCO2 trend wide shelf std']/np.sqrt(data['dpCO2 trend wide shelf_n']),xlabel='Stokes proportion of\ntotal shelf break current (%)',ax=ax1,lett=t)

    ax1 = ax[10]
    t = plot_scatter(xdata=data['sst'],dataerr=data['sstSD']/np.sqrt(data['sst_n']),ydata =data['dpCO2 trend wide shelf mean'],ydata_err=data['dpCO2 trend wide shelf std']/np.sqrt(data['dpCO2 trend wide shelf_n']),xlabel='Sea Surface Temperature ($^o$C)',ax=ax1,lett=t)
    # ax1 = ax[11]
    # t = plot_scatter(xdata=data['eks_trend'],dataerr=data['eks_trend_std']/np.sqrt(data['eks_trend_n']),ydata =data['dpCO2 trend wide shelf mean'],ydata_err=data['dpCO2 trend wide shelf std']/np.sqrt(data['dpCO2 trend wide shelf_n']),xlabel='Mean shelf \n break Ekman current trend (ms$^{-1}$ yr$^{-1}$)',ax=ax1,lett=t)
    # ax1 = ax[12]
    # t = plot_scatter(xdata=data['geo_trend'],dataerr=data['geo_trend_std']/np.sqrt(data['eks_trend_n']),ydata =data['dpCO2 trend wide shelf mean'],ydata_err=data['dpCO2 trend wide shelf std']/np.sqrt(data['dpCO2 trend wide shelf_n']),xlabel='Mean shelf \n break Geostrophic current trend (ms$^{-1}$ yr$^{-1}$)',ax=ax1,lett=t)
    #
    # ax1 = ax[13]
    # t = plot_scatter(xdata=data['stokes_trend'],dataerr=data['stokes_trend_std']/np.sqrt(data['stokes_trend_n']),ydata =data['dpCO2 trend wide shelf mean'],ydata_err=data['dpCO2 trend wide shelf std']/np.sqrt(data['dpCO2 trend wide shelf_n']),xlabel='Mean shelf \n break Stokes current trend (ms$^{-1}$ yr$^{-1}$)',ax=ax1,lett=t)

    return t

def plot_scatter_grid_trend(data,ax,text,t = [0,0],size=22):
    ax[0].annotate(text, xy=(0, 0.5), xytext=(-ax[0].yaxis.labelpad - 20, 0),
                xycoords=ax[0].yaxis.label, textcoords='offset points', ha='right', va='center',weight='bold')
    ax1 = ax[0]
    t = plot_scatter(xdata=data['eks_trend'],dataerr=data['eks_trend_std']/np.sqrt(data['eks_trend_n']),ydata =data['dpCO2 trend wide shelf mean'],ydata_err=data['dpCO2 trend wide shelf std']/np.sqrt(data['dpCO2 trend wide shelf_n']),xlabel='Mean shelf \n break Ekman current trend (ms$^{-1}$ yr$^{-1}$)',ax=ax1,lett=t,size=size)
    ax1 = ax[1]
    t = plot_scatter(xdata=data['geo_trend'],dataerr=data['geo_trend_std']/np.sqrt(data['eks_trend_n']),ydata =data['dpCO2 trend wide shelf mean'],ydata_err=data['dpCO2 trend wide shelf std']/np.sqrt(data['dpCO2 trend wide shelf_n']),xlabel='Mean shelf \n break Geostrophic current trend (ms$^{-1}$ yr$^{-1}$)',ax=ax1,lett=t,size=size)

    ax1 = ax[2]
    t = plot_scatter(xdata=data['stokes_trend'],dataerr=data['stokes_trend_std']/np.sqrt(data['stokes_trend_n']),ydata =data['dpCO2 trend wide shelf mean'],ydata_err=data['dpCO2 trend wide shelf std']/np.sqrt(data['dpCO2 trend wide shelf_n']),xlabel='Mean shelf \n break Stokes current trend (ms$^{-1}$ yr$^{-1}$)',ax=ax1,lett=t,size=size)

    return t

def load_data(file):
    print(file)
    data = pd.read_table(file,sep=',')
    return data


in_fold = 'C:/Users/df391/OneDrive - University of Exeter/Post_Doc_Covex_Seascape/Shutler_Cross_Shelf_Transport/csv/'
in_files = ['table2data_cmems_global_500m_annual_1995_2015.csv','table2data_cmems_global_500m_winter_1995_2015.csv','table2data_cmems_global_500m_spring_1995_2015.csv','table2data_cmems_global_500m_summer_1995_2015.csv','table2data_cmems_global_500m_autumn_1995_2015.csv']
season = ['Annual','Winter\nBoreal: Jan-Mar\nAustral: Jul-Sep','Spring\nBoreal: Apr-Jun\nAustral: Oct-Dec','Summer\nBoreal: Jul-Sep\nAustral: Jan-Mar','Autumn\nBoreal: Oct-Dec\nAustral: Apr-Jun']
plot_location = 'plots/data_means/'
"""
"""
# # #print(data.columns)
# # """
# # Linear Regression scatter plots for information - pretty much every combination of the data in table 1.
# # """
# # fig = plt.figure(figsize=(70,35))
# # gs = GridSpec(5,10, figure=fig, wspace=0.33,hspace=0.25,bottom=0.05,top=0.95,left=0.08,right=0.975)
# # for j in range(0,len(in_files)):
# #     ax = [fig.add_subplot(gs[j,i]) for i in range(0,10)]
# #     data = load_data(in_fold+in_files[j])
# #     plot_scatter_grid(data,ax,season[j])
# # # # print(scipy.stats.linregress(np.abs(data['EksAbsolute']),data['dpCO2 trend wide shelf mean']))
# # # # print(scipy.stats.linregress(np.abs(data['EksAbsolute']),data['dpCO2 trend wide shelf mean']).rvalue**2)
# # # # print(scipy.stats.linregress(np.abs(data['total across-shelf']),data['dpCO2 trend wide shelf mean']))
# # # # print(scipy.stats.linregress(np.abs(data['total across-shelf']),data['dpCO2 trend wide shelf mean']).rvalue**2)
# # # # print(scipy.stats.linregress(np.abs(data['EksAbsolute']/(data['total across-shelf']-data['EksAbsolute'])),data['dpCO2 trend wide shelf mean']))
# # # # print(scipy.stats.linregress(np.abs(data['EksAbsolute']/(data['total across-shelf']-data['EksAbsolute'])),data['dpCO2 trend wide shelf mean']).rvalue**2)
# # fig.savefig('plots/scatter_plots.png',format='png',dpi=300)
# #
# font = {'weight' : 'normal',
#         'size'   : 18}
# matplotlib.rc('font', **font)
# def opacity_bar(han):
#     [bar.set_alpha(0.5) for bar in han[2]]
#     [cap.set_alpha(0.5) for cap in han[1]]
#
# # """
# # Figure for manuscript!
# # """
# ylim = np.array([-1.75,4])
#
# data = load_data(in_fold+'table2data_cmems_global_500m_winter_1995_2015.csv')
# nam = []
# for i in data['region']:
#     nam.append(i.split('(')[0])
# print(nam)
# cm = plt.cm.get_cmap('tab20',14)
# print(cm)
# fig2 = plt.figure(figsize=(19,7))
# gs = GridSpec(1,2, figure=fig2, wspace=0.20,hspace=0.2,bottom=0.15,top=0.95,left=0.05,right=0.95)
# ax1 = fig2.add_subplot(gs[0,0])
# c = ax1.errorbar(data['EksAbsolute'],data['dpCO2 trend wide shelf mean'],yerr=data['dpCO2 trend wide shelf std']/np.sqrt(data['dpCO2 trend wide shelf_n']),xerr=data['EksAbsoluteStd']/np.sqrt(data['EksAbsolute_n']),linestyle='',color='b')
# opacity_bar(c)
# ax1.scatter(data['EksAbsolute'],data['dpCO2 trend wide shelf mean'],zorder=3,c='b')
# ax1.set_ylabel('d$\Delta$pCO$_2$/dt ($\mu$atm yr$^{-1}$)')
# ax1.set_xlabel('Mean wintertime\n shelf break Ekman current ($ms^{-1}$)')
# offset_vals = [
#     [-0.02,0.25],# North Sea
#     [-0.01,-0.5], # English Channel
#     [0.01,0.25], # Mid Atlantic Bight
#     [0.005,0.25], # Coast of Japan
#     [0.005,-0.25], # Patagonia
#     [0.01,0], # Bering
#     [0.01,0.25],# Antartica
#     [-0.02,-.5], # Labrador Sea
#     [0.01,0], # Tasmania
#     [0.01,0.25], # Barents
#     [0.01,0.5], # South Atlantic Bight
#     [-0.02,0.75], #Greenland
#     [0.01,0.5], # Cascadian
#     [0.01,0.3] # Irminger
# ]
# for i in range(len(nam)):
#     ax1.annotate(nam[i],xy= (data['EksAbsolute'][i],data['dpCO2 trend wide shelf mean'][i]), xycoords='data',xytext = (data['EksAbsolute'][i]+offset_vals[i][0],data['dpCO2 trend wide shelf mean'][i]+offset_vals[i][1]),
#        arrowprops=dict(arrowstyle="->"), fontsize=14)
#
# ax1.grid()
# ax1.text(0.90,0.95,'(a)',transform=ax1.transAxes,va='top',fontsize=18,fontweight='bold')
# a = np.arange(np.min(data['EksAbsolute']),np.max(data['EksAbsolute'])+0.0025,0.1/100)
# popt, pcov = curve_fit(func, data['EksAbsolute'],data['dpCO2 trend wide shelf mean'])
# coef,pvalue = scipy.stats.pearsonr(data['dpCO2 trend wide shelf mean'],func(data['EksAbsolute'],*popt))
# quad_boot = bootstrap_quad(data['EksAbsolute'],data['dpCO2 trend wide shelf mean'])
# #stat = scipy.stats.linregress(np.abs(data['EksAbsolute']),data['dpCO2 trend wide shelf mean'])
# #ax1.plot(a,a*stat.slope + stat.intercept,'k--')
# ax1.plot(a,func(a,*popt),'k--',linewidth=2)
# if pvalue < 0.01:
#     ax1.text(0.05,0.95,'r$^{2}$ = ' + str(np.round(coef**2,2)) +' $\pm$ ' + str(np.round(quad_boot,2)) +'\np-value < 0.01' + '\nn = ' + str(len(data['EksAbsolute'])),transform=ax1.transAxes,va='top',fontsize=18)
# else:
#     ax1.text(0.05,0.95,'r$^{2}$ = ' + str(np.round(coef**2,2)) +' $\pm$ ' + str(np.round(quad_boot,2)) +'\np-value = ' + str(np.round(pvalue,2)) + '\nn = ' + str(len(data['EksAbsolute'])),transform=ax1.transAxes,va='top',fontsize=18)
# ax1.set_ylim(ylim)
# ax1.set_xlim([-0.045,0.06])
#
# data = load_data(in_fold+'table2data_cmems_global_500m_autumn_1995_2015.csv')
# ax2 = fig2.add_subplot(gs[0,1])
# c = ax2.errorbar(data['EksAbsolute'],data['dpCO2 trend wide shelf mean'],yerr=data['dpCO2 trend wide shelf std']/np.sqrt(data['dpCO2 trend wide shelf_n']),xerr=data['EksAbsoluteStd']/np.sqrt(data['EksAbsolute_n']),linestyle='',color='b')
# opacity_bar(c)
# cbar_scat = ax2.scatter(data['EksAbsolute'],data['dpCO2 trend wide shelf mean'],zorder=3,c='b')
# ax2.set_ylabel('d$\Delta$pCO$_2$/dt ($\mu$atm yr$^{-1}$)')
# ax2.set_xlabel('Mean autumntime\n shelf break Ekman current ($ms^{-1}$)')
# ax2.grid()
#
# offset_vals = [
#     [0.01,0],# North Sea
#     [-0.01,-0.75], # English Channel
#     [0,1], # Mid Atlantic Bight
#     [-0.005,-0.75], # Coast of Japan
#     [0.01,0], # Patagonia
#     [0.01,0], # Bering
#     [0.005,0.5],# Antartica
#     [-0.015,-.5], # Labrador Sea
#     [0.005,-1], # Tasmania
#     [0.01,0.25], # Barents
#     [-0.01,0.5], # South Atlantic Bight
#     [-0.02,0.75], #Greenland
#     [0.01,0.5], # Cascadian
#     [-0.01,0.3] # Irminger
# ]
# for i in range(len(nam)):
#     ax2.annotate(nam[i],xy= (data['EksAbsolute'][i],data['dpCO2 trend wide shelf mean'][i]), xycoords='data',xytext = (data['EksAbsolute'][i]+offset_vals[i][0],data['dpCO2 trend wide shelf mean'][i]+offset_vals[i][1]),
#        arrowprops=dict(arrowstyle="->"), fontsize=14)
# a = np.arange(np.min(data['EksAbsolute'])-0.0025,np.max(data['EksAbsolute'])+0.0025,0.1/100)
# popt, pcov = curve_fit(func, data['EksAbsolute'],data['dpCO2 trend wide shelf mean'])
# coef,pvalue = scipy.stats.pearsonr(data['dpCO2 trend wide shelf mean'],func(data['EksAbsolute'],*popt))
# quad_boot = bootstrap_quad(data['EksAbsolute'],data['dpCO2 trend wide shelf mean'])
# # stat = scipy.stats.linregress(np.abs(data['EksAbsolute']),data['dpCO2 trend wide shelf mean'])
# # ax2.plot(a,a*stat.slope + stat.intercept,'k--')
# ax2.plot(a,func(a,*popt),'k--',linewidth=2)
# if pvalue < 0.01:
#     ax2.text(0.05,0.95,'r$^{2}$ = ' + str(np.round(coef**2,2)) +' $\pm$ ' + str(np.round(quad_boot,2)) +'\np-value < 0.01' + '\nn = ' + str(len(data['EksAbsolute'])),transform=ax2.transAxes,va='top',fontsize=18)
# else:
#     ax2.text(0.05,0.95,'r$^{2}$ = ' + str(np.round(coef**2,2)) +' $\pm$ ' + str(np.round(quad_boot,2)) +'\np-value = ' + str(np.round(pvalue,2)) + '\nn = ' + str(len(data['EksAbsolute'])),transform=ax2.transAxes,va='top',fontsize=18)
# ax2.text(0.90,0.95,'(b)',transform=ax2.transAxes,va='top',fontsize=18,fontweight='bold')
# ax2.set_ylim(ylim)
# ax2.set_xlim([-0.04,0.055])
# # print(data['region'])
#
# # cax = fig2.add_axes([0.85, .15, .01, 0.8])
# # cbar = fig2.colorbar(cbar_scat, orientation='vertical',cax=cax)
# # cbar.set_ticks(np.arange(0,len(data))+0.5)
# # cbar.set_ticklabels(nam)
#
# fig2.savefig('plots/scatter_plots_manu.png',format='png',dpi=300)
# plt.close(fig2)

# """
# GAM approach
# """
#
# font = {'weight' : 'normal',
#         'size'   : 14}
# def bootstrap_gam(x, y,ens=100):
#     r2 = []
#     for i in range(0,100):
#         a = list(np.random.choice(np.arange(0,len(x)), len(x)-2, replace=False))
#         y2 = y[a]
#         if x.shape[1] == 1:
#             x2 = x[a]
#         else:
#             x2 = x[a,:]
#         gam = GAM(s(0)+s(1),fit_intercept=False,lam=[0.6]*x2.shape[1]).gridsearch(x2, y2,n_splines=np.arange(4,10,1))
#         r,pvalue = scipy.stats.pearsonr(gam.predict(x2),y2)
#         r2.append(r**2)
#         #r2.append(gam.statistics_['pseudo_r2']['explained_deviance'])
#     return np.std(r2)
#
# from pygam import GAM, s, f,l
# matplotlib.rc('font', **font)
# out = []
# for j in ['annual','winter','spring','summer','autumn']:
#     data = load_data(in_fold+'table2data_cmems_global_500m_'+j+'_1995_2015.csv')
#     #from sklearn.linear_model import LinearRegression
#
#     #
#     x = np.array(data[['EksAbsolute','k']])
#     print(x)
#     y = np.array(data['dpCO2 trend wide shelf mean'])
#     lam = np.logspace(-2, 1, 10)
#     lams = [lam] * 2
#     # popt, pcov = curve_fit(func, data['EksAbsolute'],data['dpCO2 trend wide shelf mean'])
#     # gam = GAM(s(0) + s(1),fit_intercept=False).fit(x,y)
#     print(x.shape)
#     gam = GAM(s(0)+s(1),fit_intercept=False,lam=[0.6]*x.shape[1]).gridsearch(x, y,n_splines=np.arange(4,10,1))
#     st = bootstrap_gam(x,y)
#     #st = 0
#     gam.summary()
#     #print(gam.predict(x))
#     r2,pvalue = scipy.stats.pearsonr(gam.predict(x),y)
#
#     #print(gam.statistics_)
#     print(j+'='+str(gam.statistics_['pseudo_r2']['explained_deviance']))
#     fig, axs = plt.subplots(1,x.shape[1],figsize=(14,7));
#     if x.shape[1] == 1:
#         axs = [axs]
#     fig_label = [j+'time Ekman\nshelf break current','Mean wintertime gas transfer']
#     label = ['EksAbsolute','k']
#     for i, ax in enumerate(axs):
#         XX = gam.generate_X_grid(term=i)
#         ax.plot(XX[:, i], gam.partial_dependence(term=i, X=XX),label='GAM relationship')
#         ax.plot(XX[:, i], gam.partial_dependence(term=i, X=XX, width=.95)[1], c='r', ls='--')
#         #ax.scatter(x[label[i]],y)
#         ax.set_title(fig_label[i])
#         a = ax.get_xlim()
#         #a = np.arange(a[0],a[1],(a[1]-a[0])/100)
#         #popt, pcov = curve_fit(func, x[label[i]],y)
#         #ax.plot(a,func(a,*popt),'k--',linewidth=2,label='Quadratic relationship')
#         ax.legend()
#         ax.grid()
#     out.append([j,gam.statistics_['pseudo_r2']['explained_deviance'],st,r2**2,pvalue])
#
#     fig.savefig('plots/'+j+'_GAM_testing.png',dpi=300)
# print(out)
# plt.show()


"""
Supplementary figure with scatter plots
"""
font = {'weight' : 'normal',
        'size'   : 18}
matplotlib.rc('font', **font)
fig = plt.figure(figsize=(77,35))

# plot_scatter_grid(df_means,ax,'')
# plt.show()
gs = GridSpec(5,11, figure=fig, wspace=0.33,hspace=0.2,bottom=0.05,top=0.95,left=0.08,right=0.975)
t = [0,0]
for j in range(0,len(in_files)):

    ax = [fig.add_subplot(gs[j,i]) for i in range(0,11)]
    data = load_data(in_fold+in_files[j])
    t = plot_scatter_grid(data,ax,season[j],t)
fig.savefig('plots/scatter_plots.png',format='png',dpi=300)
fig.savefig('plots/scatter_plots.pdf',format='pdf',dpi=300)

# font = {'weight' : 'normal',
#         'size'   : 22}
# matplotlib.rc('font', **font)
# fig = plt.figure(figsize=(28,35))
#
# # plot_scatter_grid(df_means,ax,'')
# # plt.show()
# gs = GridSpec(5,3, figure=fig, wspace=0.33,hspace=0.2,bottom=0.05,top=0.95,left=0.17,right=0.975)
# t = [0,0]
# for j in range(0,len(in_files)):
#
#     ax = [fig.add_subplot(gs[j,i]) for i in range(0,3)]
#     data = load_data(in_fold+in_files[j])
#     t = plot_scatter_grid_trend(data,ax,season[j],t)
# fig.savefig('plots/scatter_plots_trend.png',format='png',dpi=300)
