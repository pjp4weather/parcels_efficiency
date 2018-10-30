#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 29 15:10:48 2018

@author: paul
"""

import netCDF4
import numpy as np 
from timeit import default_timer as timer
import matplotlib.pyplot as plt


particles = 1000000
time_steps_arr = np.arange(100,1000,100)
timer_arr = np.zeros(len(time_steps_arr))

appending = False
known_shape = False

#%%
for i in range(len(time_steps_arr)):
    print time_steps_arr[i]    
    a = np.random.uniform(-1.,1.,(particles,time_steps_arr[i]))
    
    dataset = netCDF4.Dataset("test_output.nc", "w", format="NETCDF4")
    
    if known_shape:
        dataset.createDimension("obs", particles)
        dataset.createDimension("time", time_steps_arr[i])
    else:    
        dataset.createDimension("obs", None)
        dataset.createDimension("time", None)
    
    var1 = dataset.createVariable("Var1","f8",('obs','time'))

    start = timer()
    if appending:
        for j in range(time_steps_arr[i]):
            var1[:,j] =  a[:,j]
    else:
        var1[:,:] =  a[:,:]
    end = timer()
    
    timer_arr[i] = end - start
    
    dataset.close()
    del a
    del var1
#%%

plt.title("Writing time")
plt.plot(time_steps_arr,timer_arr,marker="o",label="Appending:"+str(appending)+", shape_known: "+str(known_shape))
plt.xlabel("time steps")
plt.ylabel("time [s]")
plt.legend()