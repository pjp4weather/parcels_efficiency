#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
MODULE picke2nc

This module contains methods that convert outputs from the pickle format to
netcdf format (ParticleFile). Pickle files can be saved per time step and id 
(id_tstep_pickle) or per time step with ALL ids in one pickle (tstep_pickle)

            
@author: paul
"""

import numpy as np
import os
from timeit import default_timer as timer
import multiprocessing as mp
import glob

def get_pfile(pfile):
    """
    Get the ParticleFile.
    
    Parameters:
    -----------
    pfile: ParticleFile instance
    
    Returns:
    -------
    Unwrapped (?) variables (id, time,lat, lon, z)
    """
    
    ids = pfile.id
    time = pfile.time
    lat = pfile.lat
    lon = pfile.lon
    z = pfile.z
    
    return ids,time,lat,lon,z

# =============================================================================
# convert from pickles saved per time step (WITH all ids)
# =============================================================================

def convert_tstep_pickle(pfile,multiProcess):
    """
    Writes outputs from pickle files from time step pickles 
    (all ids in one file per time step) to ParticleFile instance
    
    Parameters:
    -----------
    pfile: ParticleFile instance
    
    multiProcess: boolean 
        Use multiprocessing for reading pickles. Currently, not possible.
    
    Returns
    -------
    conversion time : float
    """
    
    def sort_list(path):
        """
        method to sort a list of all pathes to output by their name that is the
        time stamp.
        """
        splitted = path.split("/")
        return float(splitted[1][:-4])    
    
    def read(file_list,id_fill):
        """
        read pickles using a loop over all files and return one array 
        for each variable. 
        
        Note: id and time arrays are filled in advance and 
        not read for each time step from the array because they are known in 
        advance.
        
        Parameters
        -----------
        file_list : list of strings
            List that  contains all file names in the output directory
        
        Returns
        -------
        Merged arrays for id, time, lat, lon, z
        
        """
        #generate sorted list
        file_list = [x[:-4] for x in file_list]
        file_list = map(float,file_list)
        file_list.sort()
        
        # infer array size of dimension id from the highest id in NPY file from 
        # last time step
        data_dict  = np.load("out/"+str(file_list[-1])+".npy").item()
        
        n_id = int(max(data_dict["ids"])+1)
        n_time = len(file_list)
        
        # dictionary for merging
        merge_dict = {}
        for var in data_dict.keys():
            merge_dict[var] = np.zeros((n_id,n_time))
            
            if var!="ids":
                merge_dict[var][:] = np.nan
            else:
                merge_dict[var][:] = id_fill
        
        # initiated indeces for time axis
        time_index = np.zeros(n_id,dtype=int)
        
        # loop over all files
        for i in range(n_time):
            data_dict = np.load("out/"+str(file_list[i])+".npy").item()
            
            # don't convert to netdcf if all values are nan for a time step
            if np.isnan(data_dict["ids"]).all():
                for key in merge_dict.keys():
                    merge_dict[key] = merge_dict[key][:,:-1]
                continue
            
            # get ids that going to be filled
            id_ind =  np.array(data_dict["ids"],dtype=int)
            
            # get the corresponding time indeces
            t_ind = time_index[id_ind]
            
            # write into merge array
            for key in merge_dict.keys():
                merge_dict[key][id_ind,t_ind] = data_dict[key]
           
            # new time index for ids that where filled with values
            time_index[id_ind]  = time_index[id_ind]  + 1
        
        # remove rows that are completely filled with nan values
        out_dict = {}
        for var in merge_dict.keys():
            out_dict[var] = merge_dict[var][~np.isnan(merge_dict["lat"]).all(axis=1)]
        
        return out_dict
        
    start_time = timer()
    
    # list of files
    time_list = os.listdir("out")
    
    
    data_dict = read(time_list,pfile.id._FillValue)

    for var in data_dict.keys():
        if var !="ids":
            getattr(pfile, var)[:,:] = data_dict[var]
        else:
            getattr(pfile, "id")[:,:] = data_dict[var]
           
    end_time = timer()
    convert_time = end_time-start_time
    
    if os.path.exists("out"):
        os.system("rm -rf out")
    
    return convert_time
 

# =============================================================================
# convert from pickles saved per time step and id
# =============================================================================
def load(file_str):
    """
    mehod to load one pickle file id_tspep_pickle
    """
    #print file_str
    arr = np.load(file_str)
    return arr[0],arr[1],arr[2],arr[3],arr[4] 


def convert_id_tspep_pickle(pfile,multiProcess):
    """
    Writes outputs from pickle files from time step pickles 
    (all ids in one file per time step) to ParticleFile instance
    
    Parameters:
    -----------
    pfile: ParticleFile instance
    
    multiProcess: boolean 
        Use multiprocessing for reading pickles. 
    
    Returns
    -------
    conversion time : float
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
        
        Parameters
        -----------
        file_list : list of strings
            List that  contains all file names in the output directory
        
        n_id : int
            Number of ids
        
        n_time : int
            Number if time steps
            
        Returns
        -------
        Merged arrays for id, time, lat, lon, z
        
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
        read pickles using a multiprocessing
        
            
        Returns
        -------
        Merged arrays for id, time, lat, lon, z
        
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
    ids,time,lat,lon,z,dataset = get_pfile(pfile)
    
    # read data and write to netcdf file
    if multiProcess:
        ids[:,:], time[:,:],lat[:,:],lon[:,:],z[:,:] = read_mp()
    
    else:
         ids[:,:], time[:,:],lat[:,:],lon[:,:],z[:,:] = read_loop(file_list,n_id,n_time)
            
    end_time = timer()
    convert_time = end_time-start_time
    
    return convert_time