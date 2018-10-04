#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  2 17:07:22 2018

@author: paul
"""

import numpy as np
import os
import netCDF4
from timeit import default_timer as timer
import multiprocessing as mp
import glob

#%%
def convert(particleset,name,multiProcess):
    """
    Main method that converts pickle files to one netcdf file
    """
    start_time = timer()
    
    # list of id directories
    id_list = os.listdir("out")
    n_id = len(id_list)
    
    # list of files
    file_list = os.listdir("out/0")
    n_time = len(file_list)
    
    # init netcdf file
    ids,time,lat,lon,z,dataset = init_netcdf(name,particleset,n_id,n_time)
    
    # read data and write to netcdf file
    if multiProcess:
        ids[:,:], time[:,:],lat[:,:],lon[:,:],z[:,:] = read_mp()
    else:
         ids[:,:], time[:,:],lat[:,:],lon[:,:],z[:,:] = read_loop(file_list,n_id,n_time)
            
    end_time = timer()
    convert_time = end_time-start_time
    
    dataset.close()
    
    return convert_time


def init_netcdf(name,particleset,n_id,n_time):
    """
    initialize a netcdf file
    """
    dataset = netCDF4.Dataset(name,"w", format="NETCDF4")
    dataset.createDimension("traj", n_id) 
    dataset.createDimension("obs", n_time)
    coords = ("traj", "obs")
    #.parcels_version = parcels_version
    dataset.parcels_mesh = particleset.fieldset.gridset.grids[0].mesh
    
    # Create ID variable according to CF conventions
    ids = dataset.createVariable("trajectory", "i4", coords)#, chunksizes=self.chunksizes)
    ids.long_name = "Unique identifier for each particle"
    ids.cf_role = "trajectory_id"
    
    time = dataset.createVariable("time", "f8", coords, fill_value=np.nan)#, chunksizes=self.chunksizes)
    time.long_name = ""
    time.standard_name = "time"
    time.axis = "T"
    
    lat = dataset.createVariable("lat", "f4", coords, fill_value=np.nan)#, chunksizes=self.chunksizes)
    lat.long_name = ""
    lat.standard_name = "latitude"
    lat.units = "degrees_north"
    lat.axis = "Y"
    
    lon = dataset.createVariable("lon", "f4", coords, fill_value=np.nan)#, chunksizes=self.chunksizes)
    lon.long_name = ""
    lon.standard_name = "longitude"
    lon.units = "degrees_east"
    lon.axis = "X"
    
    z = dataset.createVariable("z", "f4", coords, fill_value=np.nan)#, chunksizes=self.chunksizes)
    z.long_name = ""
    z.standard_name = "depth"
    z.units = "m"
    z.positive = "down"
    
    return ids,time,lat,lon,z, dataset


def load(file_name):
    """
    mathod to load one pickle file
    """
    arr = np.load(file_name)
    return arr[0],arr[1],arr[2],arr[3],arr[4]  


def sort_list(path):
    """
    method to sort a list of all pathes to output by 
        1. ID
        2. time stamp
    """
    splitted = path.split("/")
    return int(splitted[1]), float(splitted[2][:-4])


def read_loop(file_list,n_id,n_time):
    """
    read pickles using a loop over all files and return one array 
    for each variable
    """
    #generate sorted list
    file_list = [x[:-4] for x in file_list]
    file_list = map(float,file_list)
    file_list.sort()
    
    # merge allocate arrays
    id_m, time_m,lon_m,lat_m,z_m =\
            map(lambda n,m: np.zeros((n,m)), \
                [n_id,n_id,n_id,n_id,n_id],[n_time,n_time,n_time,n_time,n_time,])
    
    # loop over all files
    for i in range(n_id):
        for j in range(n_time):
            arr = np.load("out/"+str(i)+"/"+str(file_list[j])+".npy")
            id_m[i,j] = arr[0]
            time_m[i,j] = arr[1]
            lat_m[i,j] = arr[2]
            lon_m[i,j] = arr[3]
            z_m[i,j] = arr[4]
            
    return id_m, time_m, lat_m, lon_m, z_m

def read_mp():
    """
    read pickles using multiprocessing
    """
    # init a pool that runs on all cores of the machine
    p = mp.Pool()
    
    # get a sorted list of all pickle file pathes
    file_list_mp = glob.glob("out/*/*.npy")
    file_list_mp.sort(key=sort_list)
    
    # read files using multiprocessing
    a = p.map(load, file_list_mp)
    a = np.array(a)
    return a[:,0],a[:,1],a[:,2],a[:,3],a[:,4]


#%%  
    
if __name__ =="__main__":
    p = mp.Pool()
    
    file1 = "/home/paul/Dropbox/Studium/UU/SOAC/Project/code/out/0/0.0.npy"
    file2 = "/home/paul/Dropbox/Studium/UU/SOAC/Project/code/out/0/3600.0.npy"
    file_list = glob.glob("out/*/*.npy")
    file_list.sort(key=sort_list)
    
    start = timer()
    a = p.map(load, file_list)
    a = np.array(a)
    
    end = timer()
    print end -start
    p.close()
    