#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  2 17:07:22 2018

@author: paul
"""

import numpy as np
import os
from timeit import default_timer as timer
import multiprocessing as mp
import glob


def get_netcdf(pfile,n_id,n_time):
    """
    initialize a netcdf file
    """
    
    ids = pfile.id
    time = pfile.time
    lat = pfile.lat
    lon = pfile.lon
    z = pfile.z
    return ids,time,lat,lon,z, pfile

# =============================================================================
# convert from pickles saved per time step (WITH all ids)
# =============================================================================

def convert_from_time(pfile,multiProcess):
    """
    Main method that converts pickle files to one netcdf file
    """
    
    def sort_list(path):
        """
        method to sort a list of all pathes to output by 
            1. time stamp
        """
        splitted = path.split("/")
        return float(splitted[1][:-4])    
    
    def read(file_list,n_id,n_time):
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
        for j in range(n_time):
            arr = np.load("out/"+str(file_list[j])+".npy")
            id_m[:,j] = arr[0,:]
            time_m[:,j] = arr[1,:]
            lat_m[:,j] = arr[2,:]
            lon_m[:,j] = arr[3,:]
            z_m[:,j] = arr[4,:]
                
        return id_m, time_m, lat_m, lon_m, z_m
    
    start_time = timer()
    
    # list of files
    time_list = os.listdir("out")
    n_time = len(time_list)
    

    n_id = len(np.load("out/"+time_list[0])[0,:])
    
    # init netcdf file
    ids,time,lat,lon,z,dataset = get_netcdf(pfile,n_id,n_time)
    
    # read data and write to netcdf file
    ids[:,:], time[:,:],lat[:,:],lon[:,:],z[:,:] = read(time_list,n_id,n_time)
            
    end_time = timer()
    convert_time = end_time-start_time
    
    return convert_time
 
# =============================================================================
# convert from pickles saved per time step and id
# =============================================================================
def load(file_str):
    """
    mathod to load one pickle file
    """
    #print file_str
    arr = np.load(file_str)
    return arr[0],arr[1],arr[2],arr[3],arr[4] 


def convert_from_id(pfile,multiProcess):
    """
    Main method that converts pickle files to one netcdf file
    """
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
        a = a.reshape(n_id,n_time,5)
        
        return a[:,:,0],a[:,:,1],a[:,:,2],a[:,:,3],a[:,:,4]
    
    start_time = timer()
    
    # list of id directories
    id_list = os.listdir("out")
    n_id = len(id_list)
    
    # list of files
    file_list = os.listdir("out/0")
    n_time = len(file_list)
    
    # init netcdf file
    ids,time,lat,lon,z,dataset = get_netcdf(pfile,n_id,n_time)
    
    # read data and write to netcdf file
    if multiProcess:
        ids[:,:], time[:,:],lat[:,:],lon[:,:],z[:,:] = read_mp()
    
    else:
         ids[:,:], time[:,:],lat[:,:],lon[:,:],z[:,:] = read_loop(file_list,n_id,n_time)
            
    end_time = timer()
    convert_time = end_time-start_time
    
    return convert_time