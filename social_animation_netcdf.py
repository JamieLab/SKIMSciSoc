#!/usr/bin/env python2
# -*- coding: utf-8 -*-

# import cPickle as pickle;
import numpy as np;
import matplotlib.pyplot as plt;
# import python_util.parameter_sets as ps;
# import python_util.skim_utilities as su;
# import python_util.mask_functions as mf;
from os import path;
from mpl_toolkits.basemap import Basemap;
from matplotlib.gridspec import GridSpec
from matplotlib.animation import FuncAnimation, PillowWriter,FFMpegWriter,ImageMagickWriter
import matplotlib as mpl
from netCDF4 import Dataset
mpl.rcParams['animation.ffmpeg_path'] = r'C:\\Users\\df391\\OneDrive - University of Exeter\\Python\ffmpeg\\bin\\ffmpeg.exe'
mpl.rcParams['animation.convert_path'] = r'C:\\Program Files\\ImageMagick-7.1.1-Q16-HDRI\\magick.exe'

def animate(j):
    fig.clear()
    #print(i)
    gs = GridSpec(1,1, figure=fig, wspace=0.2,hspace=0.2,bottom=0.1,top=0.95,left=0.1,right=0.98)
    ax = fig.add_subplot(gs[0,0])
    mapFig = Basemap(llcrnrlon=-180.0, llcrnrlat=-90, urcrnrlon=180.0, urcrnrlat=90.0, resolution='l', projection='cyl', lat_0 = 39.5, lon_0 = -3.25,ax=ax);
    mapFig.drawmapboundary(fill_color=(0.85, 0.85, 0.85));
    mapFig.fillcontinents(color=(0.5, 0.5, 0.5), zorder=1);
    mapFig.drawmeridians(np.arange(0, 360, 60), labels=[0,0,0,1], color=(0.3, 0.3, 0.3), fontsize=ticksize);
    mapFig.drawparallels(np.arange(-90, 90, 30), labels=[1,0,0,0], color=(0.3, 0.3, 0.3), fontsize=ticksize);
    tc = total_cur[:,:,j]
    f = np.where(np.isnan(tc) == 0)
    a = mapFig.scatter(long[f], latg[f], c = tc[f],marker='o',latlon=True,cmap=plt.cm.RdBu,vmin=-0.5,vmax=0.5);
    cbar = plt.colorbar(a,orientation="vertical") #ticks=[0.0, 0.5, 1.0]);
    cbar.set_label("Total cross-shelf break current ($ms^{-1}$)", fontsize=ticksize);
    print(j)
    global year
    global mon

    ax.set_title('Year: ' + str(int(year)) + ' - Month: ' +str(int(mon)))
    mon = mon + 1
    if mon == 13:
        year = year+1
        mon = 1

file = 'F:/ConvexSeascapeSurvey/SurfaceFlows/Shutler_et_al_500m_total_current_data_CMEMS_1993_2016_v2-0.nc'

c = Dataset(file,'r')
lat = np.array(c['latitude'])
lon = np.array(c['longitude'])
latg,long = np.meshgrid(lat,lon)
print(long.shape)
total_cur = np.array(c['total_current'])
print(total_cur.shape)
c.close()

ticksize = 10;
year = 1993
mon =1
fig = plt.figure(figsize=(10,5))
#
ani = FuncAnimation(fig, animate, interval=40, blit=False, repeat=False, frames=total_cur.shape[2])
# #plt.show()
ani.save('plots/animated_netcdf.mp4',dpi=300, writer=FFMpegWriter(fps=6))
#ani.save('plots/animated.gif',dpi=300,writer=ImageMagickWriter(fps=12))
