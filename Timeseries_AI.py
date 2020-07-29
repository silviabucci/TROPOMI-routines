#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 23:38:30 2020

@author: Silvia Bucci
"""
import pandas as pd
from datetime import datetime,timedelta
import os
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import pickle,gzip
from mpl_toolkits.axes_grid1 import make_axes_locatable
from extra_functions import color_scale_rainbow,cartopy_only_map
from matplotlib.colors import LogNorm
mymap=color_scale_rainbow()
pathAI='/dsk3/data/s5p/AI/'

#Define the domain and its resolution
longrid=np.arange(-180,180,0.25)
latgrid=np.arange(-85,0,0.25)
source_range=np.array([[longrid.min(),longrid.max()],[latgrid.min(),latgrid.max()]])
dateseq=pd.date_range(pd.datetime(2019,12,15,0,0),pd.datetime(2020,2,29,0,0),freq='24H')
dateseq=dateseq.to_pydatetime()

#Load the area of the domain (saved by save_area_km2.py) 
area=pickle.load(gzip.open('area_for_vortex.pkl','rb'))
area0=area[:-1,:-1].T/1e6
area_value=[]
ai_value=[]
date_value=[]
#Read saved data from read_sentinel5p_AI_CO_O3
#%%
for idate in dateseq:    
    HH=pickle.load(gzip.open('../s5p_AI/s5p_AI_'+idate.strftime('%Y-%m-%d')+'.pkl','rb'))
    ai_m=HH['AI']        
    xedges=HH['lon']    
    yedges=HH['lat']
    #compute area when AI>3
    war=ai_m>3
    date_value.append(idate)

    wlat=np.where(yedges<-70)[0]            
    area0[:,wlat]=0         
    area_value.append(np.sum(area0[war]))

#%%
#Compute timeseries of 95% of AI values near the Australian coast
ai_value=[]
ai_date=[]
for idate in dateseq:
    path0='/dsk3/data/s5p/AI/'
    for filename in os.listdir(path0):        
        if 'AER_AI_'+idate.strftime('%Y%m%d') in filename:
            file=Dataset(os.path.join(path0,filename),'r')
            print(filename)
            utc=file.groups['PRODUCT'].variables['time_utc'][:][0,:]
            utc=[datetime.strptime(i[0:19], '%Y-%m-%d'+'T'+'%H:%M:%S') for i in utc]
            wt=np.where(np.logical_and(np.array(utc)>=idate,np.array(utc)<idate+timedelta(seconds=86400)))[0]
            if len(wt)>0:
                lons = file.groups['PRODUCT'].variables['longitude'][:][0,:,:]
                lonsf=lons.flatten()
                lats = file.groups['PRODUCT'].variables['latitude'][:][0,:,:]
                var = file.groups['PRODUCT'].variables['aerosol_index_354_388'][0,:,:]
                var_units = file.groups['PRODUCT'].variables['aerosol_index_354_388'].units
                lon_0 = lons.mean()
                lat_0 = lats.mean()
                var[var>1e5]=0
                var0=var.filled(-1e18)
                wo=[var0.flatten()>0]
                lon0=150
                lon1=155
                lat0=-40
                lat1=-20
                mm=np.where(np.logical_and(np.logical_and(np.array(lons)>=lon0,np.array(lons)<=lon1),np.logical_and(np.array(lats)>=lat0,np.array(lats)<=lat1)))[0]
                try:
                    var[var<0]=0
                    var=np.ma.array(var,mask=var==0)
                    print(mm[0])
                    ai_value.append(np.percentile(var[mm],95))
                    ai_date.append(idate)
                except:
                    print('')    
#%% 

line=np.array(ai_value)
dateline=np.array(ai_date)
dateline2=[]
line2=[]
dateline=np.array(dateline)
line=np.array(line)
a,b=np.unique(dateline,return_index=True)
for i in range(len(a)):
    if np.sum(dateline==a[i])<2:        
        line2.append(line[dateline==a[i]][0])
    else:
        line2.append(np.max(line[dateline==a[i]]))
dateline2=np.array(a)
line2=np.array(line2)

#%%
#Plot timeserie
fs=15
fig, ax = plt.subplots(figsize=[12,9],ncols=1,nrows=2,sharex=True); ax[0].plot(date_value,np.array(area_value),'ob-') ;ax[0].set_ylabel('million km2',fontsize=fs) ; ax[0].set_title('Surface of plumes with AI > 3',fontsize=fs)
fig.subplots_adjust(wspace=0.15)
ax[0].set_xlim(datetime(2019,12,15),datetime(2020,1,31))
ax[0].set_ylim(0,8)
ax[1].set_xlim(datetime(2019,12,15),datetime(2020,1,31))
ax[1].set_ylim(0,12)
ax[1].plot(dateline2,line2,'ob-',color='k')  ;ax[1].set_ylabel('AI',fontsize=fs) ; ax[1].set_title('95% of the AI values near the Australian coast (150E-155E 40S-20S)',fontsize=fs)

fig.autofmt_xdate()
dayss=[datetime(2020,1,2),datetime(2020,1,4),datetime(2020,1,7),datetime(2020,1,10),datetime(2020,1,13)]
wt=[np.where(np.array(date_value)==da)[0][0] for da in dayss]
for ii in range(0,len(wt),1):
    ax[0].text(date_value[wt[ii]],area_value[wt[ii]]+0.5,date_value[wt[ii]].strftime('%d-%b'),fontsize=fs)                     

dayss=[datetime(2019,12,31),datetime(2020,1,2),datetime(2020,1,5)]
wt=[np.where(np.array(dateline2)==da)[0][0] for da in dayss]        
for ii in range(0,len(wt),1):
    ax[1].text(dateline2[wt[ii]],line2[wt[ii]]+0.6,dateline2[wt[ii]].strftime('%d-%b'),fontsize=fs)                             

import matplotlib 
matplotlib.rc('xtick', labelsize=fs) 
matplotlib.rc('ytick', labelsize=fs)         
#plt.savefig('/home/sbucci/routines/AI_surface_and_percentile.png',bbox_inches='tight',dpi=300)        
#%%
#Read composite of AI and centroid of the blobs from read_sentinel5p_AI_CO_O3
import cartopy.crs as ccrs

AI=pickle.load(gzip.open('AI_composite_tot.pkl','rb'))
#read vortex 1 and 2 positions saved from read_Sentinel5p_AI_CO_O3
a=pickle.load(gzip.open('AI_1v.pkl','rb'))
lona=np.array(a)[:,2] ;lata=np.array(a)[:,3]
b=pickle.load(gzip.open('AI_2v.pkl','rb'))
lonb=np.array(b)[:,2] ;latb=np.array(b)[:,3]

AIm=AI['AI']/AI['AI0']

source_range=np.array([[130,130+360],[-80,-10]])
proj = ccrs.PlateCarree(central_longitude=(source_range[0,1]-source_range[0,0])/2+source_range[0,0])
fig,ax2=plt.subplots(1,1,figsize=[10,10],subplot_kw=dict(projection=proj))
xticks=np.arange(-180,180,60)

#Plot composite AI    
fig,ax=cartopy_only_map(cmap=mymap,fig=fig,ax=ax2,\
                boundaries=[source_range[0,0],source_range[0,1],source_range[1,0],source_range[1,1]],proj=proj,xticks=xticks)
cs=ax.pcolormesh(xedges+(180-source_range[0,0]),yedges,AIm[:,:].T,cmap=mymap,vmin=0.1,vmax=5,norm=LogNorm())
plt.title('Mean AI and bubbles centroid 2020-01-06 to 2020-02-29')
ms=4
lona=np.array(lona);lata=np.array(lata)
lonb=np.array(lonb);latb=np.array(latb)

lonb2=np.append(lonb[0:2]-360.,lonb[2:])
lona2=np.append(lona[:-11],lona[-11:]-360.)

ax.plot(lona2+(180-source_range[0,0]),lata,'ob-',c='r',markersize=ms)    
ax.plot(lonb2+(180-source_range[0,0]),latb,'ob-',c='dimgray',markersize=ms)    

wb=lonb2>130
wp=np.min(np.where(wb==1)[0])-1
wb[wp]=True
ax.plot(lonb2[wb]-360+(180-source_range[0,0]),latb[wb],'ob-',c='dimgray',markersize=ms)    

wa=lona2<-180-(180-source_range[0,0])
wp=np.min(np.where(wa==1)[0])-1
wa[wp]=True
ax.plot(lona2[wa]+360+(180-source_range[0,0]),lata[wa],'ob-',c='r',markersize=ms)        

ax.plot(lonb2[0]+(180-source_range[0,0]),latb[0],'ob-',c='dimgray',label='Bubble 2, from 4th January',markersize=ms)    
ax.plot(lona2[0]+(180-source_range[0,0]),lata[0],'+b-',c='r',label='Bubble 1, from 1st January',markersize=ms)
ax.plot(lonb2[0]+(180-source_range[0,0]),latb[0],'s',c='dimgray',markersize=ms*1.6)    
ax.plot(lona2[0]+(180-source_range[0,0]),lata[0],'s',c='r',markersize=ms*1.6)
ax.plot(lonb2[0]+(180-source_range[0,0]),latb[0],'s',c='k',markersize=ms*1.2)    
ax.plot(lona2[0]+(180-source_range[0,0]),lata[0],'s',c='k',markersize=ms*1.2)
ax.plot(lonb2+(180-source_range[0,0]),latb,'o',c='k',markersize=ms)
ax.plot(lona2+(180-source_range[0,0]),lata,'+',c='k',markersize=ms*1.2)
ax.plot(lonb2[wb]-360+(180-source_range[0,0]),latb[wb],'o',c='k',markersize=ms)    
ax.plot(lona2[wa]+360+(180-source_range[0,0]),lata[wa],'+',c='k',markersize=ms*1.2)        


ax.legend(bbox_to_anchor=(0.05, 0.90, 0.96,0.8))
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05,axes_class=plt.Axes)
cbar0=plt.colorbar(cs, cax=cax,label='AI')
#plt.savefig('/home/sbucci/routines/plumes_composite.png',bbox_inches='tight',dpi=300)
plt.show()




