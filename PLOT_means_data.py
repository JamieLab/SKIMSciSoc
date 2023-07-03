#!/usr/bin/env python3

import pandas as pd
import numpy as np
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import matplotlib.transforms
import scipy.stats

font = {'weight' : 'normal',
        'size'   : 24}
matplotlib.rc('font', **font)

def plot_scatter(xdata = None, dataerr = None,ydata=None,ydata_err=None,xlabel = None,ax = None,let = 'a'):
    ax.errorbar(xdata,ydata,yerr=ydata_err,xerr=dataerr,linestyle='',color='k')
    c = ax.scatter(xdata,ydata,zorder=3)
    ax.set_ylabel('d$\Delta$pCO$_2$/dt ($\mu$atm yr$^{-1}$)')
    ax.set_xlabel(xlabel)
    #plt.colorbar(c,ax=axs[0])
    ax.grid()
    a = np.array([np.amin(xdata),np.amax(xdata)])
    stat = scipy.stats.linregress(xdata,ydata)
    ax.plot(a,a*stat.slope + stat.intercept,'r--',linewidth=2)
    ax.text(0.05,0.95,'R$^{2}$ = ' + str(np.round(stat.rvalue**2,2)) +'\np-value = ' + str(np.round(stat.pvalue,2)) + '\nn = ' + str(len(xdata)),transform=ax.transAxes,va='top',fontsize=22)
    #ax.text(0.90,0.95,'('+let+')',transform=ax.transAxes,va='top',fontsize=18,fontweight='bold')

def plot_scatter_grid(data,ax,text):

    ax[0].annotate(text, xy=(0, 0.5), xytext=(-ax[0].yaxis.labelpad - 20, 0),
                xycoords=ax[0].yaxis.label, textcoords='offset points',
                size='large', ha='right', va='center',weight='bold')
    ax1 = ax[0]
    plot_scatter(xdata=data['EksAbsolute'],dataerr=data['EksAbsoluteStd'],ydata =data['dpCO2 trend wide shelf mean'],ydata_err=data['dpCO2 trend wide shelf std'], xlabel='Mean shelf \n break Ekman current (ms$^{-1}$)',ax=ax1,let='a')#
    ax2 = ax[1]
    plot_scatter(xdata=np.abs(data['EksAbsolute']),dataerr=data['EksAbsoluteStd'],ydata =data['dpCO2 trend wide shelf mean'],ydata_err=data['dpCO2 trend wide shelf std'],xlabel='Absolute mean\n shelf break Ekman current ($ms^{-1}$)',ax=ax2,let='b')#
    ax1 = ax[2]
    plot_scatter(xdata=data['total across-shelf'],dataerr=data['total across-shelf SD'],ydata =data['dpCO2 trend wide shelf mean'],ydata_err=data['dpCO2 trend wide shelf std'],xlabel='Mean shelf\n break total current (ms$^{-1}$)',ax=ax1,let='c')#
    ax2 = ax[3]
    plot_scatter(xdata=np.abs(data['total across-shelf']),dataerr=data['total across-shelf SD'],ydata =data['dpCO2 trend wide shelf mean'],ydata_err=data['dpCO2 trend wide shelf std'],xlabel='Absolute mean shelf\n break total current (ms$^{-1}$)',ax=ax2,let='d')#
    ax1 = ax[4]
    plot_scatter(xdata=data['GeoAbsolute'],dataerr=data['GeoAbsoluteStd'],ydata =data['dpCO2 trend wide shelf mean'],ydata_err=data['dpCO2 trend wide shelf std'],xlabel='Mean shelf\n break geostropic current (ms$^{-1}$)',ax=ax1,let='e')#
    ax1 = ax[5]
    plot_scatter(xdata=np.abs(data['GeoAbsolute']),dataerr=data['GeoAbsoluteStd'],ydata =data['dpCO2 trend wide shelf mean'],ydata_err=data['dpCO2 trend wide shelf std'],xlabel='Absolute mean shelf\n break geostrophic current (ms$^{-1}$)',ax=ax1,let='f')#
    ax1 = ax[6]
    plot_scatter(xdata=data['k_ms-1_10-5'],dataerr=data['k_ms-1_10-5_SD'],ydata =data['dpCO2 trend wide shelf mean'],ydata_err=data['dpCO2 trend wide shelf std'],xlabel='Mean gas transfer (10$^{-5}$ ms$^{-1}$)',ax=ax1,let='h')
    ax1 = ax[7]
    plot_scatter(xdata=data['proportionEks']*100,dataerr=data['proportionEksSD']*100,ydata =data['dpCO2 trend wide shelf mean'],ydata_err=data['dpCO2 trend wide shelf std'],xlabel='Ekman proportion of\ntotal shelf break current (%)',ax=ax1,let='i')
    ax1 = ax[8]
    plot_scatter(xdata=data['proportionGeo']*100,dataerr=data['proportionGeoSD']*100,ydata =data['dpCO2 trend wide shelf mean'],ydata_err=data['dpCO2 trend wide shelf std'],xlabel='Geostrophic proportion of\ntotal shelf break current (%)',ax=ax1,let='j')
    ax1 = ax[9]
    plot_scatter(xdata=data['proportionStokes']*100,dataerr=data['proportionStokesSD']*100,ydata =data['dpCO2 trend wide shelf mean'],ydata_err=data['dpCO2 trend wide shelf std'],xlabel='Stokes proportion of\ntotal shelf break current (%)',ax=ax1,let='k')

def load_data(file):
    data = pd.read_excel(file)
    return data
in_fold = 'C:/Users/df391/OneDrive - University of Exeter/Post_Doc_Covex_Seascape/Shutler_Cross_Shelf_Transport/'
in_files = ['table2data_global_500m_annual.xlsx','table2data_global_500m_winter.xlsx','table2data_global_500m_spring.xlsx','table2data_global_500m_summer.xlsx','table2data_global_500m_autumn.xlsx']
season = ['Annual','Winter\nBoreal: Jan-Mar\nAustral: Jul-Sep','Spring\nBoreal: Apr-Jun\nAustral: Oct-Dec','Summer\nBoreal: Jul-Sep\nAustral: Jan-Mar','Autumn\nBoreal: Oct-Dec\nAustral: Apr-Jun']
plot_location = 'plots/data_means/'
#
# #print(data.columns)
# """
# Linear Regression scatter plots for information - pretty much every combination of the data in table 1.
# """
fig = plt.figure(figsize=(70,35))
gs = GridSpec(5,10, figure=fig, wspace=0.33,hspace=0.25,bottom=0.05,top=0.95,left=0.08,right=0.975)
for j in range(0,len(in_files)):
    ax = [fig.add_subplot(gs[j,i]) for i in range(0,10)]
    data = load_data(in_fold+in_files[j])
    plot_scatter_grid(data,ax,season[j])
# # print(scipy.stats.linregress(np.abs(data['EksAbsolute']),data['dpCO2 trend wide shelf mean']))
# # print(scipy.stats.linregress(np.abs(data['EksAbsolute']),data['dpCO2 trend wide shelf mean']).rvalue**2)
# # print(scipy.stats.linregress(np.abs(data['total across-shelf']),data['dpCO2 trend wide shelf mean']))
# # print(scipy.stats.linregress(np.abs(data['total across-shelf']),data['dpCO2 trend wide shelf mean']).rvalue**2)
# # print(scipy.stats.linregress(np.abs(data['EksAbsolute']/(data['total across-shelf']-data['EksAbsolute'])),data['dpCO2 trend wide shelf mean']))
# # print(scipy.stats.linregress(np.abs(data['EksAbsolute']/(data['total across-shelf']-data['EksAbsolute'])),data['dpCO2 trend wide shelf mean']).rvalue**2)
fig.savefig('plots/scatter_plots.png',format='png',dpi=300)
#
font = {'weight' : 'normal',
        'size'   : 18}
matplotlib.rc('font', **font)
def opacity_bar(han):
    [bar.set_alpha(0.5) for bar in han[2]]
    [cap.set_alpha(0.5) for cap in han[1]]

# """
# Figure for manuscript!
# """
data = pd.read_excel(in_fold+'table2data_global_500m_winter.xlsx')
cm = plt.cm.get_cmap('RdYlBu')
fig2 = plt.figure(figsize=(15,7))
gs = GridSpec(1,2, figure=fig2, wspace=0.33,hspace=0.2,bottom=0.15,top=0.95,left=0.1,right=0.95)
ax1 = fig2.add_subplot(gs[0,0])
c = ax1.errorbar(np.abs(data['EksAbsolute']),data['dpCO2 trend wide shelf mean'],yerr=data['dpCO2 trend wide shelf std'],xerr=data['EksAbsoluteStd'],linestyle='',color='b')
opacity_bar(c)
ax1.scatter(np.abs(data['EksAbsolute']),data['dpCO2 trend wide shelf mean'],zorder=3,c=np.sign(data['EksAbsolute']),vmin=-0.01,vmax=0.01,cmap=cm)
ax1.set_ylabel('d$\Delta$pCO$_2$/dt ($\mu$atm yr$^{-1}$)')
ax1.set_xlabel('Absolute mean wintertime\n shelf break Ekman current ($ms^{-1}$)')
ax1.annotate('North Sea',xy= (0.03418,1.81), xycoords='data',xytext = (0.05,1.5),
    arrowprops=dict(arrowstyle="->"), fontsize=12)
ax1.annotate('English Channel',xy= (0.01471,-0.03), xycoords='data',xytext = (0.03,-0.8),
    arrowprops=dict(arrowstyle="->"), fontsize=12)
ax1.grid()
ax1.text(0.90,0.95,'(a)',transform=ax1.transAxes,va='top',fontsize=18,fontweight='bold')
a = np.array([0,0.04])
stat = scipy.stats.linregress(np.abs(data['EksAbsolute']),data['dpCO2 trend wide shelf mean'])
ax1.plot(a,a*stat.slope + stat.intercept,'k--')
ax1.text(0.05,0.95,'R$^{2}$ = ' + str(np.round(stat.rvalue**2,2)) +'\np-value = ' + str(np.round(stat.pvalue,2)) + '\nn = ' + str(len(data['EksAbsolute'])),transform=ax1.transAxes,va='top',fontsize=18)


data = pd.read_excel(in_fold+'table2data_global_500m_autumn.xlsx')
ax2 = fig2.add_subplot(gs[0,1])
c = ax2.errorbar(np.abs(data['EksAbsolute']),data['dpCO2 trend wide shelf mean'],yerr=data['dpCO2 trend wide shelf std'],xerr=data['EksAbsoluteStd'],linestyle='',color='b')
opacity_bar(c)
ax2.scatter(np.abs(data['EksAbsolute']),data['dpCO2 trend wide shelf mean'],zorder=3,c = np.sign(data['EksAbsolute']),vmin=-0.01,vmax=0.01,cmap=cm)
ax2.set_ylabel('d$\Delta$pCO$_2$/dt ($\mu$atm yr$^{-1}$)')
ax2.set_xlabel('Absolute mean autumntime\n shelf break Ekman current ($ms^{-1}$)')
ax2.grid()
a = np.array([0,0.04])
ax2.annotate('North Sea',xy= (0.026763115,1.81), xycoords='data',xytext = (0.05,1.5),
    arrowprops=dict(arrowstyle="->"), fontsize=12)
ax2.annotate('English Channel',xy= (0.013474338,-0.03), xycoords='data',xytext = (0.03,-0.8),
    arrowprops=dict(arrowstyle="->"), fontsize=12)
stat = scipy.stats.linregress(np.abs(data['EksAbsolute']),data['dpCO2 trend wide shelf mean'])
ax2.plot(a,a*stat.slope + stat.intercept,'k--')
ax2.text(0.05,0.95,'R$^{2}$ = ' + str(np.round(stat.rvalue**2,2)) +'\np-value = ' + str(np.round(stat.pvalue,2)) + '\nn = ' + str(len(data['EksAbsolute'])),transform=ax2.transAxes,va='top',fontsize=18)
ax2.text(0.90,0.95,'(b)',transform=ax2.transAxes,va='top',fontsize=18,fontweight='bold')
fig2.savefig('plots/scatter_plots_manu.png',format='png',dpi=300)
#
# # from sklearn.linear_model import LinearRegression
# # import statsmodels.api as sm
# # #
# # inp = sm.add_constant(np.transpose(np.array([np.abs(data['EksAbsolute']),np.abs(data['GeoAbsolute']),np.abs(data['total across-shelf']),data['k_ms-1_10-5']])))
# # out = sm.GLM(np.array(data['dpCO2 trend wide shelf mean']),inp,family=sm.families.Poisson())
# # result = out.fit()
# # print(result.summary())

# data = pd.read_excel(in_fold+in_files[1])
# data2 = pd.read_excel(in_fold+in_files[4])
# df_concat = pd.concat((data,data2))
# by_row_index = df_concat.groupby(df_concat.index)
# df_means = by_row_index.mean()
# print(df_means)
# fig = plt.figure(figsize=(70,35))
#
# # plot_scatter_grid(df_means,ax,'')
# # plt.show()
# gs = GridSpec(5,10, figure=fig, wspace=0.33,hspace=0.2,bottom=0.05,top=0.95,left=0.08,right=0.975)
# for j in range(0,len(in_files)):
#
#     ax = [fig.add_subplot(gs[0,i]) for i in range(0,10)]
#     data = load_data(in_fold+in_files[j])
#     plot_scatter_grid(data,ax,season[j])
# fig.savefig('plots/scatter_plots.png',format='png',dpi=300)
