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


particles = 100000
time_steps_arr = np.arange(100,1000,100)
timer_arr = np.zeros(len(time_steps_arr))

appending = False
known_shape = False


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

np.save('A_'+str(appending)+"_S_"+str(known_shape)+"_N_"+str(particles),np.array([time_steps_arr,timer_arr]))

#%%
plt.close("all")
plt.figure(figsize=(12,6))
particles = np.array([1000,10000,100000])
c = ['g','b','r']
ls = ["-","-","--"]

for j in range(3):
    for i in range(1,3):
        if i==0:
            appending = True
            known_shape = False
        if i==1:
            appending = False
            known_shape = False
        if i ==2:
            known_shape = True
            
        steps, time = np.load('A_'+str(appending)+"_S_"+str(known_shape)+"_N_"+str(particles[j])+".npy")
        #plt.loglog(steps,time,marker="o",label="Appending:"+str(appending)+", shape_known: "+str(known_shape))
        
        plt.loglog(steps,time,label="Appending: False, limited: "+str(known_shape) + ", particles: " + str(particles[j]),c=c[j],ls=ls[i])
    
particles = 1000
appending = True
known_shape = False
steps, time = np.load('A_'+str(appending)+"_S_"+str(known_shape)+"_N_"+str(particles)+".npy")
plt.loglog(steps,time,c="k",lw=2.,label="Appending: "+str(appending)+", limited: "+str(known_shape)+ ", particles: " + str(particles))

plt.title("Writing time",fontsize=14)
plt.xlabel("time steps",fontsize=14)
plt.ylabel("time [s]",fontsize=14)
plt.yticks(fontsize=14)
plt.xticks(fontsize=14)
plt.ylim(0.001,100)
plt.xlim(100,1000)
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize=12)
plt.subplots_adjust(left=0.1, right=0.5, top=0.9, bottom=0.1)
