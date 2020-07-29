#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 12:27:27 2019

@author: silvia bucci
"""
import numpy 
import pickle,gzip
from mpl_toolkits.basemap import Basemap
xedges=numpy.arange(-180,180,0.25)
yedges=numpy.arange(-85,0,0.25)
source_range=numpy.array([[xedges.min(),xedges.max()],[yedges.min(),yedges.max()]])

area=numpy.zeros([len(yedges),len(xedges)])
for i in range(len(xedges)-1):
    for j in range(len(yedges)-1):
        coordinates=numpy.array([
        [xedges[i],yedges[j+1]], 
        [xedges[i], yedges[j]], 
        [xedges[i+1], yedges[j]], 
        [xedges[i+1], yedges[j+1]],
        [xedges[i], yedges[j+1]]])               
        lats=coordinates[:,1]
        lons=coordinates[:,0]
        lat1=numpy.min(lats)
        lat2=numpy.max(lats)
        lon1=numpy.min(lons)
        lon2=numpy.max(lons)
        bmap=Basemap(projection='cea',llcrnrlat=lat1,llcrnrlon=lon1,urcrnrlat=lat2,urcrnrlon=lon2)
        xs,ys=bmap(lons,lats)
        ar=numpy.abs(0.5*numpy.sum(ys[:-1]*numpy.diff(xs)-xs[:-1]*numpy.diff(ys)))                
        ar=ar/1e6
        area[j,i]=ar
        print(i,j,ar) 
pickle.dump(area,gzip.open('area_for_vortex.pkl','wb'))
