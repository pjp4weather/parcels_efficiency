#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  4 16:06:43 2018

@author: paul
"""

from parcels import FieldSet, ParticleSet, Variable, JITParticle, AdvectionRK4, plotTrajectoriesFile
import numpy as np
import math
from datetime import timedelta
from operator import attrgetter
import matplotlib.pyplot as plt
import os
from pickle_to_netcdf import convert
from timeit import default_timer as timer

start = timer()

os.system("rm -rf out")
os.mkdir("out")

write_pickle = False
multi_process = True




fieldset = FieldSet.from_parcels("/home/paul/parcels_examples/MovingEddies_data/moving_eddies")

#pset = ParticleSet.from_list(fieldset=fieldset,   # the fields on which the particles are advected
#                             pclass=JITParticle,  # the type of particles (JITParticle or ScipyParticle)
#                             lon=[3.3e5,  3.3e5], # a vector of release longitudes 
#                             lat=[1e5, 2.8e5])    # a vector of release latitudes

particle_number = 100
np.random.seed(0)
pset = ParticleSet.from_list(fieldset=fieldset,   # the fields on which the particles are advected
                             pclass=JITParticle,  # the type of particles (JITParticle or ScipyParticle)
                             lon = np.random.uniform(low=3.3e5, high=3.4e5, size=particle_number),  # a vector of release longitudes 
                             lat=np.random.uniform(low=1.5e5, high=2.8e5, size=particle_number)  )  # a vector of release latitudes

pfile = pset.ParticleFile(name="EddyParticles.nc", outputdt=timedelta(hours=1))

iters = 24*6
write_time_arr = np.zeros(iters)

for i in range(24*6):
    pset.execute(AdvectionRK4,                 # the kernel (which defines how particles move)
                 runtime=timedelta(hours=1),    # the total length of the run
                 dt=timedelta(minutes=5))      # the timestep of the kernel
    
    start_writing = timer()
    pfile.write(pset, pset[0].time,pickle = write_pickle)
    end_writing = timer()
    
    write_time_arr[i] = end_writing - start_writing

if write_pickle:
    time_convert = convert(pset,"EddyParticles2.nc",multi_process)
    plotTrajectoriesFile("EddyParticles2.nc")

else:
    time_convert = 0.
    plotTrajectoriesFile("EddyParticles.nc")

end = timer()


# get process time
time_total = end - start
time_writing = np.sum(write_time_arr)
time_integration = time_total - time_writing - time_convert

    

# =============================================================================
# Plot
# =============================================================================
plt.close("all")

fig,ax = plt.subplots(1,1)

ind = np.arange(0,4)
tt, ti,tw, tc  = ax.bar(ind,[time_total,time_integration,time_writing,time_convert])

tt.set_facecolor('k')
ti.set_facecolor('tomato')
tw.set_facecolor('lime')
tc.set_facecolor('steelblue')
ax.set_xticks(ind)
ax.set_xticklabels(['Total','Integration', 'Writing', 'Converting'])
ax.set_ylabel("s")
plt.ylim(0,18)
ax.set_title("Number of particles: " + str(particle_number) +", pickle: " + str(write_pickle) +", mp: " +str(multi_process))