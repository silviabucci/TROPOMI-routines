#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 12:12:03 2020

@author: Silvia Bucci
"""
from cartopy import feature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.axes_grid1 import make_axes_locatable


def color_scale_rainbow():

    #define a new colorscale
    cdict = {'red': ((0., 1, 1),
                     (0.05, 1, 1),
                     (0.11, 0, 0),
                     (0.66, 1, 1),
                     (0.89, 1, 1),
                     (1, 0.5, 0.5)),
             'green': ((0., 1, 1),
                       (0.05, 1, 1),
                       (0.11, 0, 0),
                       (0.375, 1, 1),
                       (0.64, 1, 1),
                       (0.91, 0, 0),
                       (1, 0, 0)),
             'blue': ((0., 1, 1),
                      (0.05, 1, 1),
                      (0.11, 1, 1),
                      (0.34, 1, 1),
                      (0.65, 0, 0),
                      (1, 0, 0))}

    my_cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    return my_cmap

def cartopy_only_map(txt='',unit='',save=False,figname='fig',dir0='./',clim=[0,1],fig=plt.figure(),ax=plt.axis,boundaries=[-180,180,-90,90],xticks=None,yticks=None,cmap='jet',proj=None):
    ax.set_global()
    print(boundaries)
    minlon=boundaries[0]
    maxlon=boundaries[1]
    minlat=boundaries[2]
    maxlat=boundaries[3]
    ax.set_title(txt,size=8)
    ax.add_feature(feature.NaturalEarthFeature(
            category='cultural',
            name='admin_1_states_provinces_lines',
            scale='50m',
            facecolor='none'))    
    ax.set_extent([minlon, maxlon, minlat, maxlat],crs=proj)
    ax.coastlines()
    if np.sum(xticks)==None:
        xticks=np.arange(minlon,maxlon,(maxlon-minlon)/10.)
    if np.sum(yticks)==None:
        yticks=np.arange(minlat,maxlat,(maxlat-minlat)/5.)
    ax.set_xticks(xticks, crs=proj)
    ax.set_yticks(yticks, crs=proj)        
    lon_formatter = LongitudeFormatter(zero_direction_label=False)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    ax.grid()
    return fig,ax

        

def cartopy_map_new(X,Y,rr,txt='',unit='',save=False,figname='fig',dir0='./',clim=[0,1],fig=plt.figure(),ax=plt.axis,boundaries=[-180,180,-90,90],xticks=None,yticks=None,cmap='jet',proj=None):
    minlon=boundaries[0]
    maxlon=boundaries[1]
    minlat=boundaries[2]
    maxlat=boundaries[3]
    wx=np.where(np.logical_and(X>=minlon,X<maxlon))[0]
    wy=np.where(np.logical_and(Y>=minlat,Y<maxlat))[0]
    x=X[wx]
    y=Y[wy]
    var=rr[wx,:]
    var=var[:,wy]
    cs=ax.pcolormesh(x,y,var.T,cmap=cmap,vmin=clim[0],vmax=clim[1],transform=proj)
    ax.set_title(txt,size=8)
    ax.add_feature(feature.NaturalEarthFeature(
            category='cultural',
            name='admin_1_states_provinces_lines',
            scale='50m',
            facecolor='none'))
    ax.set_extent([minlon, maxlon, minlat, maxlat],crs=proj)
    ax.coastlines() 
    gl = ax.gridlines(crs=proj, draw_labels=True,
                  linewidth=2, color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    if np.sum(np.array(xticks))!=None:
        gl.xlocator = mticker.FixedLocator(xticks)
    if np.sum(np.array(yticks))!=None:
        gl.ylocator = mticker.FixedLocator(yticks)
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'color': 'k','size':8}
    gl.ylabel_style = {'color': 'k','size':8}
    
    divider = make_axes_locatable(ax)
    ax_cb = divider.append_axes('right',size="5%", pad=0.1, axes_class=plt.Axes)
    fig.add_axes(ax_cb)
    cbar=plt.colorbar(cs, cax=ax_cb,orientation='vertical')
    cbar.ax.tick_params(labelsize=7)
    cbar.ax.set_ylabel(unit, position=(30,1), fontsize=7, rotation=0)

    return
