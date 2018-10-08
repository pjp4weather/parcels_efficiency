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
from pickle2nc import convert_from_time, convert_from_id
from timeit import default_timer as timer



class tests():
    """
    class that contains all tests that should be performed with a new writing 
    routine
    """
    
    def __init__(self, write_routine,multi_process=False, particle_number = 500):
        
        self.write_routine = write_routine
    
        if self.write_routine == "write_pickle":
            self.convert = convert_from_time
            
        elif self.write_routine =="write_1pickle":
            self.convert = convert_from_id
            self.multi_process = multi_process
            
        if write_routine != "write_1pickle" and multi_process:
            raise AttributeError("multi processing for conversion just available for 'write_1pickle'")
        
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
        
#%%        
if __name__ == "__main__":
    
    write_routine = "write"
    
    t = tests(write_routine,multi_process=False,particle_number=500)
    
    t.timing()
    t.plot_timing()
