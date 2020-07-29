#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 23:38:30 2020

@author: Silvia Bucci 
"""
import os
import pandas as pd
from datetime import datetime,timedelta
import numpy as np
import matplotlib.pyplot as plt
import pickle,gzip
import cartopy.crs as ccrs
from netCDF4 import Dataset
from extra_functions import color_scale_rainbow,cartopy_map_new
mymap=color_scale_rainbow()

pathAI='/dsk3/data/s5p/AI/'
pathO3='/dsk3/data/s5p/O3/'
pathCO='/dsk3/data/s5p/CO/'

#choose here wich variable to read
plotAI=True
plotO3=False
plotCO=False
#Choose here is you want to save the AI features position by click
click=False
track = []    

#Geographical domain
longrid=np.arange(-180,180,0.25)
latgrid=np.arange(-85,0,0.25)
source_range=np.array([[longrid.min(),longrid.max()],[latgrid.min(),latgrid.max()]])
#Temporal domain
dateseq=pd.date_range(pd.datetime(2019,12,31,0,0),pd.datetime(2020,2,29,0,0),freq='24H')
dateseq=dateseq.to_pydatetime()



#
def on_click(event):
    print('latitude=%f longitude=%f'%(event.xdata,event.ydata))

def on_key_ai(event):
    print(dir(event))
    xx=[np.abs(i-event.xdata) for i in xedges]
    yy=[np.abs(i-event.ydata) for i in yedges]
    a=np.where(np.array(xx)==np.min(xx))[0]
    b=np.where(np.array(yy)==np.min(yy))[0]
    print('you pressed ',event.key,event.xdata,event.ydata,xedges[a],yedges[b],ai_m[a,b][0],(ai_tim[a,b][0])/nstep)
    track.append([event.key,idate,event.xdata,event.ydata,ai_m[a,b][0],(ai_tim[a,b][0])/nstep])    

#%%
#READ THE VARIABLES!    
#Initialize variables
co0_tot=0    
co_tot=0    

o3_tot=0    
o30_tot=0   

ai0_tot=0    
ai_tot=0    
 

for idate in dateseq:
    nstep=0
    nstepco=0
    nstepo3=0
    co00=0
    o300=0
    ai00=0
    co_t=0
    o3_t=0
    ai_t=0
    co_tim=0
    o3_tim=0
    ai_tim=0
    if plotCO:
        path0=pathCO
        for filename in os.listdir(path0):        
            if 'CO_____'+idate.strftime('%Y%m%d') in filename:
                file=Dataset(os.path.join(path0,filename),'r')
                print(filename)
                utc=file.groups['PRODUCT'].variables['time_utc'][:][0,:]
                utc=[datetime.strptime(i[0:19], '%Y-%m-%d'+'T'+'%H:%M:%S') for i in utc]
                wt=np.where(np.logical_and(np.array(utc)>=idate,np.array(utc)<idate+timedelta(seconds=86400)))[0]
                if len(wt)>0:
                    nstepco=nstepco+1
                    print(nstepco)
                    lons = file.groups['PRODUCT'].variables['longitude'][:][0,:,:]
                    lonsf=lons.flatten()
                    lats = file.groups['PRODUCT'].variables['latitude'][:][0,:,:]
                    co = file.groups['PRODUCT'].variables['carbonmonoxide_total_column'][0,:,:]
                    co_units = file.groups['PRODUCT'].variables['carbonmonoxide_total_column'].units
                    lon_0 = lons.mean()
                    lat_0 = lats.mean()
                    co[co>1e5]=0
                    coo=co.filled(-1e18)
                    wo=[coo.flatten()>0]
                    co0,xedges,yedges=np.histogram2d(lonsf.flatten()[wo],lats.flatten()[wo],bins=(longrid,latgrid))
                    co_h,xedges,yedges=np.histogram2d(lonsf.flatten()[wo],lats.flatten()[wo],bins=(longrid,latgrid),weights=coo.flatten()[wo])
                    co_h=np.ma.array(co_h,mask=co_h<-1e10)
                    co00=co00+co0
                    co_t=co_t+co_h                    
                    coco=(co0/co0)*utc[int(len(utc)/2)].hour
                    coco[np.isnan(coco)]=0
                    co_tim=co_tim+coco
    
        co_tot=co_tot+co_t
        co0_tot=co0_tot+co00
        if np.sum(co_t)!=0:
            co_m=co_t/co00
            co_m[co00==0]=0
            vmin=0.0
            vmax=0.05           
            fig,ax=plt.subplots(figsize=[18,10])
            cartopy_map_new(xedges[:-1],yedges[:-1],co_m,cmap=mymap,clim=[vmin,vmax],fig=fig,ax=ax,boundaries=[source_range[0,0],source_range[0,1],source_range[1,0],source_range[1,1]])        
                    # Add Title
            plt.title('CO total column '+idate.strftime('%Y-%m-%d'))
            plt.savefig('/home/sbucci/Documents/Australia/sentinel/s5p_CO_'+idate.strftime('%Y-%m-%d'),bbox_inches='tight')            
            plt.pause(7)
            plt.close('all')
            plt.clf()
            file.close()
            HH={}
            HH['CO']=co_m
            HH['lon']=xedges[:-1]
            HH['lat']=yedges[:-1]
            pickle.dump(HH,gzip.open('s5p_CO_'+idate.strftime('%Y-%m-%d')+'.pkl','wb'))
    else:
        print('NO CO HERE')

    
    if plotAI:
        path0=pathAI
        for filename in os.listdir(path0):        
            if 'AER_AI_'+idate.strftime('%Y%m%d') in filename:
                file=Dataset(os.path.join(path0,filename),'r')
                print(filename)
                utc=file.groups['PRODUCT'].variables['time_utc'][:][0,:]
                utc=[datetime.strptime(i[0:19], '%Y-%m-%d'+'T'+'%H:%M:%S') for i in utc]
                wt=np.where(np.logical_and(np.array(utc)>=idate,np.array(utc)<idate+timedelta(seconds=86400)))[0]
                if len(wt)>0:
                    nstep=nstep+1
                    print(nstep)
                    lons = file.groups['PRODUCT'].variables['longitude'][:][0,:,:]
                    lonsf=lons.flatten()
                    lats = file.groups['PRODUCT'].variables['latitude'][:][0,:,:]
                    ai = file.groups['PRODUCT'].variables['aerosol_index_340_380'][0,:,:]
                    ai_units = file.groups['PRODUCT'].variables['aerosol_index_340_380'].units
                    lon_0 = lons.mean()
                    lat_0 = lats.mean()
                    ai[ai>1e5]=0
                    aio=ai.filled(-1e18)
                    wo=[aio.flatten()>0]
                    ai0,xedges,yedges=np.histogram2d(lonsf.flatten()[wo],lats.flatten()[wo],bins=(longrid,latgrid))
                    ai_h,xedges,yedges=np.histogram2d(lonsf.flatten()[wo],lats.flatten()[wo],bins=(longrid,latgrid),weights=aio.flatten()[wo])
                    ai_h=np.ma.array(ai_h,mask=ai_h<-1e10)
                    ai00=ai00+ai0
                    ai_t=ai_t+ai_h                    
                    aiai=(ai0/ai0)*utc[int(len(utc)/2)].hour
                    aiai[np.isnan(aiai)]=0
                    ai_tim=ai_tim+aiai
    
        ai_tot=ai_tot+ai_t
        ai0_tot=ai0_tot+ai00
        if np.sum(ai_t)!=0:
            ai_m=ai_t/ai00
            ai_m[ai00==0]=0
            vmin=0
            vmax=2
            proj = ccrs.PlateCarree()
            fig,ax2=plt.subplots(1,1,figsize=[18,10],subplot_kw=dict(projection=proj))
            cartopy_map_new(xedges[:-1],yedges[:-1],ai_m,txt='AI in atmosphere '+idate.strftime('%Y-%m-%d'),cmap=mymap,clim=[vmin,vmax],fig=fig,ax=ax2,\
                            boundaries=[source_range[0,0],source_range[0,1],source_range[1,0],source_range[1,1]],proj=proj)            

            if click:
                cid1 = fig.canvas.mpl_connect('button_press_event', on_click)
                cid2 = fig.canvas.mpl_connect('key_press_event', on_key_ai)        
            plt.savefig('/home/sbucci/Documents/Australia/sentinel/s5p_AI_'+idate.strftime('%Y-%m-%d'),bbox_inches='tight')            
            plt.pause(7)
            plt.close('all')
            plt.clf()
            file.close()
            HH={}
            HH['AI']=ai_m
            HH['lon']=xedges[:-1]
            HH['lat']=yedges[:-1]
            pickle.dump(HH,gzip.open('s5p_AI_'+idate.strftime('%Y-%m-%d')+'.pkl','wb'))
    else:
        print('NO AI HERE')

    if plotO3:
        print('OZONE')
        path0=pathO3
        for filename in os.listdir(pathO3):                             
            if np.logical_or('O3_____'+idate.strftime('%Y%m%d') in filename,'O3_____'+(idate+timedelta(seconds=86400)).strftime('%Y%m%d') in filename):
                file=Dataset(os.path.join(path0,filename),'r')
                print(filename)
                utc=file.groups['PRODUCT'].variables['time_utc'][:][0,:]
                utc=datetime(2010,1,1,0,0,0)+timedelta(seconds=np.float(file.groups['PRODUCT'].variables['time'][:][0]))
                wt=np.where(np.logical_and(np.array(utc)>=idate,np.array(utc)<idate+timedelta(seconds=86400)))[0]                
                if np.logical_and(np.array(utc)>=idate,np.array(utc)<idate+timedelta(seconds=90000)):
                    nstepo3=nstepo3+1
                    if nstepo3<16:
                        print(nstepo3)
                        lonso3 = file.groups['PRODUCT'].variables['longitude'][:][0,:,:]
                        lonso3f=lonso3.flatten()
                        lato3s = file.groups['PRODUCT'].variables['latitude'][:][0,:,:]
                        o3 = file.groups['PRODUCT'].variables['ozone_total_vertical_column'][0,:,:]
                        qa_o3=file.groups['PRODUCT']['qa_value'][0,:,:]                        
                        co_units = file.groups['PRODUCT'].variables['ozone_total_vertical_column'].units
                        lon_0 = lonso3.mean()
                        lat_0 = lato3s.mean()
                        o3[o3>1e25]=0
                        o3o=o3.filled(-1e18)
                        o3[qa_o3<0.7]=0
                        wo=[o3o.flatten()>0]
                        o30,xedges,yedges=np.histogram2d(lonso3f.flatten()[wo],lato3s.flatten()[wo],bins=(longrid,latgrid))
                        o3_h,xedges,yedges=np.histogram2d(lonso3f.flatten()[wo],lato3s.flatten()[wo],bins=(longrid,latgrid),weights=o3o.flatten()[wo])
                        o3_h=np.ma.array(o3_h,mask=o3_h<=0)
                        o300=o300+o30
                        o3_t=o3_t+o3_h.filled(0)
                        coco=(o30/o30)*utc.hour
                        coco[np.isnan(coco)]=0
                        o3_tim=o3_tim+coco
        o3_tot=o3_tot+o3_t
        o30_tot=o30_tot+o300
        if np.sum(o3_t)!=0:
            vmin=0.11
            vmax=0.16
            o3_m=o3_t/o300
            o3_m[o300==0]=0            
            source_range=np.array([[120,300],[-70,-20]])            
            nlow=xedges[:-1]<120
            nhigh=xedges[:-1]>=120       
            xedges[:-1]=np.append(xedges[:-1][nhigh],xedges[:-1][nlow]+360)
            o3_ms=np.append(o3_m[nhigh,:],o3_m[nlow,:],axis=0)
            proj = ccrs.PlateCarree()
            fig,ax2=plt.subplots(1,1,figsize=[18,10],subplot_kw=dict(projection=proj))
            cartopy_map_new(xedges[:-1],yedges[:-1],o3_ms,cmap=mymap,clim=[vmin,vmax],fig=fig,ax=ax2,\
                            boundaries=[source_range[0,0],source_range[0,1],source_range[1,0],source_range[1,1]],proj=proj)
            plt.title('O3 in atmosphere '+idate.strftime('%Y-%m-%d'))
            plt.close('all')
            plt.clf()
            file.close()
            HH={}
            HH['O3']=o3_m
            HH['lon']=xedges[:-1]
            HH['lat']=yedges[:-1]
            pickle.dump(HH,gzip.open('s5p_O3_'+idate.strftime('%Y-%m-%d')+'.pkl','wb'))            

#%%       
#Save here the composite of AI along the whole period            
if plotAI:
    H={}
    H['AI']=ai_tot
    H['AI0']=ai0_tot
    H['AIlon']=xedges[-1:]
    H['AIlat']=yedges[-1:]
    pickle.dump(H,gzip.open('AI_composite_tot.pkl','wb'))            
#Save here the track of the features selected by click
if plotAI:
    if click:
        tr=np.array(track)
        pickle.dump(tr,gzip.open('AI.pkl','wb'))

