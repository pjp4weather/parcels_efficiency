#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  4 16:06:43 2018

@author: paul
"""
from parcels import FieldSet, ParticleSet, JITParticle, AdvectionRK4, plotTrajectoriesFile, Variable
import numpy as np
from datetime import timedelta
import matplotlib.pyplot as plt
from timeit import default_timer as timer
import xarray as xr

plt.close("all")

def DeleteParticle(particle, fieldset, time, dt):
    print("Particle lost !! (%g %g %g %g)" % (particle.lon, particle.lat, particle.depth, particle.time))
    particle.delete()
    
def Ageing(particle, fieldset, time, dt):
    particle.age += dt 

class PlasticParticle(JITParticle):
    age = Variable('age', dtype=np.float32, initial=0.)

class tests():
    """
    class that contains all tests that should be performed with a new writing 
    routine
    """
    
    def __init__(self, particle_number = 500):        
        self.test_particle_number = particle_number                    
    def timing(self,integration_time_days=6):
        """
        Method to compute the time needed for the integration, writing of data
        and conversion to netcdf
        """
        np.random.seed(0)
        start = timer()
        
        iters = integration_time_days * 24
                
        write_time_arr = np.zeros(iters)
        
        fieldset = FieldSet.from_parcels("/home/paul/parcels_examples/MovingEddies_data/moving_eddies")
        pset = ParticleSet.from_list(fieldset=fieldset,   # the fields on which the particles are advected
                                     pclass=PlasticParticle,  # the type of particles (JITParticle or ScipyParticle)
                                     lon = np.random.uniform(low=3.3e5, high=3.4e5, size=self.test_particle_number),  # a vector of release longitudes 
                                     lat=np.random.uniform(low=1.5e5, high=2.8e5, size=self.test_particle_number))    # a vector of release latitudes
                                     
        output_name = "output.nc"
        pfile = pset.ParticleFile(name=output_name, outputdt=timedelta(hours=1))
        
        kernel = AdvectionRK4 + pset.Kernel(Ageing)
        
        for i in range(iters):
            
            pset.execute(kernel,                 # the kernel (which defines how particles move)
                         runtime=timedelta(hours=1),    # the total length of the run
                         dt=timedelta(minutes=5),     # the timestep of the kernel
                         )  
                        
            
            start_writing = timer()
            pfile.write(pset, pset[0].time)
            end_writing = timer()
            
            write_time_arr[i] = end_writing - start_writing
        
        end = timer()
        # get process time
        self.time_total = end - start
        self.time_writing = np.sum(write_time_arr)
        self.time_integration = self.time_total - self.time_writing
        
    def plot_timing(self):
        """
        plot the results of the timing test
        """
        #plt.close("all")
        
        fig,ax = plt.subplots(1,1)
        fs = 14
        ind = np.arange(0,3)
        tt, ti,tw  = ax.bar(ind, [self.time_total,self.time_integration,self.time_writing])
        
        tt.set_facecolor('k')
        ti.set_facecolor('tomato')
        tw.set_facecolor('lime')
        ax.set_xticks(ind)
        ax.set_xticklabels(['Total','Integration', 'Writing'],fontsize=fs)
        ax.set_ylabel("time [s]",fontsize=fs)
        ax.set_title("Number of particles: " + str(self.test_particle_number))

#%%        
if __name__ == "__main__":
    
    tt = tests(particle_number=100)
    tt.timing(integration_time_days=6)
    tt.plot_timing()
    plotTrajectoriesFile("output.nc")
    dataset = xr.open_dataset("output.nc")
    print dataset