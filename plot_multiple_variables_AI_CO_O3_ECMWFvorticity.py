#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 23:38:30 2020

@author: Silvia Bucci
"""
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
import pickle,gzip
from extra_functions import color_scale_rainbow
import cartopy
import cartopy.io.shapereader as shpreader
from cartopy import feature
import cartopy.crs as ccrs
mymap=color_scale_rainbow()

pathAI='/home/sbucci/routines/s5p_AI/'
pathO3='/home/sbucci/routines/s5p_O3/'
pathCO='/home/sbucci/routines/s5p_CO/'
path_V='/home/sbucci/routines/satie_copy/'

v=True


def cartopy_map(X,Y,rr,txt='',unit='',save=False,figname='fig',dir0='./',clim=[0,1],fig=plt.figure(),ax=plt.axis,cmap='jet',xax=1,yax=1):

    maxlon=-60
    minlon=-180
    minlat=-65
    maxlat=-35
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
    import matplotlib.ticker as mticker
    gl.xlocator = mticker.FixedLocator([-165,-135,-105,-75])
    gl.ylocator = mticker.FixedLocator([-60,-50,-40])
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'color': 'k','size':8}
    gl.ylabel_style = {'color': 'k','size':8}
    if xax==0:
        gl.xlabels_bottom = False
    if yax==0:
        gl.ylabels_left = False
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax)
    ax_cb = divider.append_axes('right',size="5%", pad=0.1, axes_class=plt.Axes)
    fig.add_axes(ax_cb)
    cbar=plt.colorbar(cs, cax=ax_cb,orientation='vertical')
    cbar.ax.tick_params(labelsize=7)
    cbar.ax.set_ylabel(unit, position=(30,1), fontsize=7, rotation=0)
    return


longrid=np.arange(-180,180,0.25)
latgrid=np.arange(-70,0,0.25)
source_range=np.array([[longrid.min(),longrid.max()],[latgrid.min(),latgrid.max()]])
dateseq=[datetime(2020,1,17),datetime(2020,2,3)]        

#Load Vorticity values from ECMWF saved from B. Legras routines
#It requires the ECMWF_N, mki2d and constants modules from Legras
filev=pickle.load(gzip.open(path_V+'OPZ-extract-1.pkl','rb'))
filev.update(pickle.load(gzip.open(path_V+'OPZ-extract-2.pkl','rb')))

for idate in dateseq:
    if idate==datetime(2020,1,17):
        livelb=49 #chose the vertical levels to extract from Vorticity field
    if idate==datetime(2020,2,3):
        livelb=39
    else:
        livelb=49
    #read TROPOMI values saved from read_sentinel5p_AI_CO_O3.py
    H=pickle.load(gzip.open(pathO3+'s5p_O3_'+idate.strftime('%Y-%m-%d.pkl'),'rb'))
    o3=H['O3']
    H=pickle.load(gzip.open(pathCO+'s5p_CO_'+idate.strftime('%Y-%m-%d.pkl'),'rb'))
    co=H['CO']
    H=pickle.load(gzip.open(pathAI+'s5p_AI_'+idate.strftime('%Y-%m-%d.pkl'),'rb'))
    ai=H['AI']    
    lon=H['lon']
    lat=H['lat']
    #Compute time index to read from the vorticity file
    il=idate-datetime(2020,1,4,6)
    istep=il.days*2+il.seconds/3600./12.-0.5+1
    #Read the variables from the vorticity file
    latv=filev[int(istep)].attr['lats']
    lonv=filev[int(istep)].attr['lons']
    levv=filev[int(istep)].attr['levs']
    V=filev[istep].var['VO']
    nlow2=np.where(lonv<180)[0]            
    nhigh2=np.where(lonv>=180)[0]
    v0=V[:,:,nlow2]
    v1=V[:,:,nhigh2]
    vn=np.append(v0,v1,axis=2)
    lonv2=np.append(lonv[nlow2],lonv[nhigh2]-360)            
        
    if len(nlow2)>0: plt.contour(lonv[nlow2],latv,V[livelb,:,nlow2].T)
    #Set the min and max values for the variables plot
    vminCO=0.0
    vmaxCO=0.04  
    vminO3=0.11
    vmaxO3=0.16
    vminAI=0
    vmaxAI=2
    vminv=-2e-5
    vmaxv=1e-4
    conversionCO=6.022141e19*1e4*1e-22 #convert CO in units 10^22 mol m^{-2}
    conversionO3=2241.15 #convert O3 in DU
    conversionV=1e4 #convert Vorticity in 10^{-5} s^{-1}
    for livel in [livelb]:
        minb={} 
        minb[0,0]=vminCO*conversionCO
        minb[0,1]=vminO3*conversionO3
        minb[1,1]=vminAI
        minb[1,0]=vminv*conversionV
        maxb={}
        maxb[0,0]=vmaxCO*conversionCO
        maxb[0,1]=vmaxO3*conversionO3
        maxb[1,1]=vmaxAI
        maxb[1,0]=vmaxv*conversionV   
        varb={}
        varb[0,0]=co*conversionCO
        varb[0,1]=o3*conversionO3
        varb[1,1]=ai
        varb[1,0]=V[livel,:,:].T*conversionV
        titb={}
        titb[0,0]='Total CO column'
        titb[0,1]='Total O3 column'
        titb[1,1]='Aerosol Absorbing Index'
        titb[1,0]=u'Vorticity'
        lonb={}
        lonb[0,0]=lon
        lonb[0,1]=lon
        lonb[1,1]=lon
        lonb[1,0]=lonv-360.
        latb={}
        latb[0,0]=lat
        latb[0,1]=lat
        latb[1,1]=lat
        latb[1,0]=latv
        unitb={}
        unitb[0,0]='(10$^{22}$) mol m$^{-2}$'
        unitb[0,1]='DU'
        unitb[1,1]='AI 380-340 nm'        
        unitb[1,0]='(10$^{-5}$ s$^{-1}$) '+'{:2.1f}'.format(filev[int(istep)].attr['zscale'][livel])+' km'
        cmapb={}
        cmapb[0,0]=mymap
        cmapb[0,1]=mymap
        cmapb[1,1]=mymap
        cmapb[1,0]='seismic'    
        xaxb={}
        xaxb[0,0]=0
        xaxb[0,1]=0
        xaxb[1,1]=1
        xaxb[1,0]=1
        yaxb={}
        yaxb[0,0]=1
        yaxb[0,1]=0
        yaxb[1,1]=0
        yaxb[1,0]=1
        #DO the 4 panels plots
        ncol=2 ; nlin=2
        proj = ccrs.PlateCarree()
        fig,axs=plt.subplots(ncol,nlin,figsize=[17,12],subplot_kw=dict(projection=proj),sharex=True,sharey=True)
        fig.subplots_adjust(bottom=0.3, top=0.5, left=0.05, right=0.45,hspace=0.1,wspace=0.15)
        fig.suptitle(idate.strftime('%Y-%m-%d'),y=0.53,x=0.25,size=10)
        ic=0
        for i in np.arange(ncol):
            for j in np.arange(nlin):   
                print(i,j)
                ic=ic+1
                var=varb[i,j]                    
                minv=minb[i,j]                    
                maxv=maxb[i,j]                    
                unit=unitb[i,j]
                tit=titb[i,j]+' '+unit
                cmap=cmapb[i,j]
                xax=xaxb[i,j]
                yax=yaxb[i,j]
                
                cartopy_map(lonb[i,j],latb[i,j],varb[i,j],txt=tit,clim=[minb[i,j],maxb[i,j]],fig=fig,\
                            ax=axs[i,j],cmap=cmap,xax=xax,yax=yax)

    fig.savefig('./s5p_AI_CO_O3_and_ECMWF_Vorticity_'+idate.strftime('%Y-%m-%d')+'{:2.1f}'.format(filev[26].attr['zscale'][livel])+'.png',bbox_inches='tight',dpi=300)            
