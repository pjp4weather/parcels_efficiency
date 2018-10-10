#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  4 16:06:43 2018

@author: paul
"""

from parcels import FieldSet, ParticleSet, JITParticle, AdvectionRK4, plotTrajectoriesFile
import numpy as np
from datetime import timedelta
import matplotlib.pyplot as plt
import os
from pickle2nc import convert_tstep_pickle, convert_id_tspep_pickle
from timeit import default_timer as timer
import xarray as xr

class tests():
    """
    class that contains all tests that should be performed with a new writing 
    routine
    """
    
    def __init__(self, write_routine,multi_process=False, particle_number = 500):
        
        self.write_routine = write_routine
        self.multi_process = multi_process
        
        if self.write_routine == "write_pickle_per_tstep":
            self.convert = convert_tstep_pickle
            
        elif self.write_routine =="write_pickle_per_id_tstep":
            self.convert = convert_id_tspep_pickle
            self.multi_process = multi_process
            
        if write_routine != "write_pickle_per_id_tstep" and multi_process:
            print "multi processing for conversion just available for 'write_1pickle'. set to 'False'"
            self.multi_process = multi_process
        
        self.test_particle_number = particle_number
    
    def timing(self):
        """
        Method to compute the time needed for the integration, writing of data
        and conversion to netcdf
        """
        np.random.seed(0)
        start = timer()

        os.system("rm -rf out")
        os.mkdir("out")
        
        integration_time_days = 6
        iters = integration_time_days * 24
                
        write_time_arr = np.zeros(iters)
        
        fieldset = FieldSet.from_parcels("/home/paul/parcels_examples/MovingEddies_data/moving_eddies")
        pset = ParticleSet.from_list(fieldset=fieldset,   # the fields on which the particles are advected
                                     pclass=JITParticle,  # the type of particles (JITParticle or ScipyParticle)
                                     lon = np.random.uniform(low=3.3e5, high=3.4e5, size=self.test_particle_number),  # a vector of release longitudes 
                                     lat=np.random.uniform(low=1.5e5, high=2.8e5, size=self.test_particle_number)  )  # a vector of release latitudes
        output_name = "output.nc"
        pfile = pset.ParticleFile(name=output_name, outputdt=timedelta(hours=1))
        
        # evaluate to the chosen writing routine
        try:
            write = eval("pfile."+self.write_routine)
        except:
            raise AttributeError("'ParticleFile' does not has writing routine '"+self.write_routine +"'")
        
        for i in range(iters):
            
            pset.execute(AdvectionRK4,                 # the kernel (which defines how particles move)
                         runtime=timedelta(hours=1),    # the total length of the run
                         dt=timedelta(minutes=5))      # the timestep of the kernel
            
            start_writing = timer()
            write(pset, pset[0].time)
            end_writing = timer()
            
            write_time_arr[i] = end_writing - start_writing
        
        self.time_convert = 0
        
        if self.write_routine !="write":
            self.time_convert = self.convert(pfile,self.multi_process)
        
        plotTrajectoriesFile(output_name)
        
        end = timer()
        
        # get process time
        self.time_total = end - start
        self.time_writing = np.sum(write_time_arr)
        self.time_integration = self.time_total - self.time_writing - self.time_convert
        
    def plot_timing(self):
        """
        plot the results of the timing test
        """
        #plt.close("all")
        
        fig,ax = plt.subplots(1,1)
        
        ind = np.arange(0,4)
        tt, ti,tw, tc  = ax.bar(ind,[self.time_total,self.time_integration,self.time_writing,self.time_convert])
        
        tt.set_facecolor('k')
        ti.set_facecolor('tomato')
        tw.set_facecolor('lime')
        tc.set_facecolor('steelblue')
        ax.set_xticks(ind)
        ax.set_xticklabels(['Total','Integration', 'Writing', 'Converting'])
        ax.set_ylabel("s")
        plt.ylim(0,100)
        ax.set_title("Number of particles: " + str(self.test_particle_number) +", " +self.write_routine +", mp: " +str(self.multi_process))
       
        
    def equality(self):
        """
        Check if new writing routine writes exactly the same output as the 
        default one
        """
        np.random.seed(0)

        os.system("rm -rf out")
        os.mkdir("out")
        
        integration_time_days = 1
        iters = integration_time_days * 12
        
        fieldset = FieldSet.from_parcels("/home/paul/parcels_examples/MovingEddies_data/moving_eddies")
        pset = ParticleSet.from_list(fieldset=fieldset,   # the fields on which the particles are advected
                                     pclass=JITParticle,  # the type of particles (JITParticle or ScipyParticle)
                                     lon = np.random.uniform(low=3.3e5, high=3.4e5, size=self.test_particle_number),  # a vector of release longitudes 
                                     lat=np.random.uniform(low=1.5e5, high=2.8e5, size=self.test_particle_number)  )  # a vector of release latitudes
       
        # output file from for default writing routine
        output_name = "output.nc"
        pfile = pset.ParticleFile(name=output_name, outputdt=timedelta(hours=1))
        
        # output file from new writing routine
        output_name_compare = "output_compare.nc"
        pfile_compare = pset.ParticleFile(name=output_name_compare, outputdt=timedelta(hours=1))
        
        # evaluate to the chosen writing routine
        try:
            write = eval("pfile_compare."+self.write_routine)
        except:
            raise AttributeError("'ParticleFile' does not has writing routine '"+self.write_routine +"'")
        
        for i in range(iters):
            pset.execute(AdvectionRK4,                 # the kernel (which defines how particles move)
                         runtime=timedelta(hours=1),    # the total length of the run
                         dt=timedelta(minutes=5))      # the timestep of the kernel
            
            # default writing routine
            pfile.write(pset, pset[0].time)
            
            # new writing routine
            write(pset, pset[0].time)

        self.time_convert = self.convert(pfile_compare,self.multi_process)
        
        self._comapre_outputs(output_name,output_name_compare)

    
    def _comapre_outputs(self,file1,file2):
        """
        Check if the two input files are identical
        
        """
        read_file1 = xr.open_dataset(file1)
        read_file2 = xr.open_dataset(file2)
        
        compare = (read_file1==read_file2).all()
        
        print compare
            
#%%        
if __name__ == "__main__":
    """
    choose  write_routine from:
        write, 
        write_pickle_per_tstep, 
        write_pickle_per_id_tstep
    """
    write_routine = "write_pickle_per_tstep"
    
    t = tests(write_routine,multi_process=False,particle_number=500)
    
    # Test timing
    t.timing()
    t.plot_timing()
    
    # Test if writing routine is identical
    t.equality()
