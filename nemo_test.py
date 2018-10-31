from parcels import FieldSet, ParticleFile, ParticleSet, JITParticle, Variable
from parcels import ErrorCode
import numpy as np
from glob import glob
from parcels import timer
import time as timelib
from datetime import timedelta as delta


print 'Running file: %s' % __file__

def get_nemo_fieldset():
    data_dir = '/data2/imau/oceanparcels/hydrodynamic_data/NEMO-MEDUSA/ORCA0083-N006/'
    ufiles = sorted(glob(data_dir+'means/ORCA0083-N06_200?????d05U.nc'))
    vfiles = sorted(glob(data_dir+'means/ORCA0083-N06_200?????d05V.nc'))
    mesh_mask = data_dir + 'domain/coordinates.nc'

    filenames = {'U': ufiles,
                 'V': vfiles,
                 'mesh_mask': mesh_mask}

    variables = {'U': 'uo',
                 'V': 'vo'}
    dimensions = {'lon': 'glamf', 'lat': 'gphif', 'time': 'time_counter'}
    fieldset = FieldSet.from_nemo(filenames, variables, dimensions)
    fieldset.U.vmax = 5
    fieldset.V.vmax = 5
    return fieldset

timer.root = timer.Timer('Main')
timer.fieldset = timer.Timer('FieldSet', parent=timer.root)
fieldset = get_nemo_fieldset()
#set_diff(fieldset)
#set_cmems(fieldset)
#set_cmems_unbeaching(fieldset)
#set_cfsr(fieldset, .01)
timer.fieldset.stop()


class PlasticParticle(JITParticle):
    age = Variable('age', dtype=np.float32, initial=0.)

def DeleteParticle(particle, fieldset, time, dt):
    print("Particle lost !! (%g %g %g %g)" % (particle.lon, particle.lat, particle.depth, particle.time))
    particle.delete()

def AdvectionRK4(particle, fieldset, time, dt):
    (u1, v1) = fieldset.UV[time, particle.lon, particle.lat, particle.depth]
    lon1, lat1 = (particle.lon + u1*.5*dt, particle.lat + v1*.5*dt)

    (u2, v2) = fieldset.UV[time + .5 * dt, lon1, lat1, particle.depth]
    lon2, lat2 = (particle.lon + u2*.5*dt, particle.lat + v2*.5*dt)

    (u3, v3) = fieldset.UV[time + .5 * dt, lon2, lat2, particle.depth]
    lon3, lat3 = (particle.lon + u3*dt, particle.lat + v3*dt)

    (u4, v4) = fieldset.UV[time + dt, lon3, lat3, particle.depth]
    particle.lon += (u1 + 2*u2 + 2*u3 + u4) / 6. * dt
    particle.lat += (v1 + 2*v2 + 2*v3 + v4) / 6. * dt

def Ageing(particle, fieldset, time, dt):
    particle.age += dt 


# Release particles
vec = np.linspace(0,1,11)
xsi, eta = np.meshgrid(vec, vec)

# Rotterdam
lonCorners = [2.96824026, 3.22713804, 3.26175451, 3.002671]
latCorners = [51.60693741, 51.58454132, 51.73711395, 51.759758] 
lon_r = (1-xsi)*(1-eta) * lonCorners[0] + xsi*(1-eta) * lonCorners[1] + \
        xsi*eta * lonCorners[2] + (1-xsi)*eta * lonCorners[3]
lat_r = (1-xsi)*(1-eta) * latCorners[0] + xsi*(1-eta) * latCorners[1] + \
        xsi*eta * latCorners[2] + (1-xsi)*eta * latCorners[3]

lonCorners = [1.37941658, 1.63887346, 1.67183721, 1.41217935]
latCorners = [51.58309555, 51.56196213, 51.71636581, 51.73773575]
lon_t = (1-xsi)*(1-eta) * lonCorners[0] + xsi*(1-eta) * lonCorners[1] + \
        xsi*eta * lonCorners[2] + (1-xsi)*eta * lonCorners[3]
lat_t = (1-xsi)*(1-eta) * latCorners[0] + xsi*(1-eta) * latCorners[1] + \
        xsi*eta * latCorners[2] + (1-xsi)*eta * latCorners[3]

lons = np.concatenate((lon_r.flatten(), lon_t.flatten()))
lats = np.concatenate((lat_r.flatten(), lat_t.flatten()))
times = np.arange(np.datetime64('2000-01-05'), np.datetime64('2001-01-05'))

pset = ParticleSet.from_list(fieldset, PlasticParticle,
                             lon=np.tile(lons, [len(times)]),
                             lat=np.tile(lats, [len(times)]),
                             time=np.repeat(times, len(lons)))

kernel = AdvectionRK4 #+ pset.Kernel(Ageing)

new_write = True

timer.particlefile = timer.Timer('ParticleFile', parent=timer.root)
outfile = './'+__file__[:-3]
pfile = ParticleFile(outfile, pset)

if new_write:
    pfile.write_pickle_per_tstep(pset, pset[0].time)
else:
    pfile.write(pset, pset[0].time)

timer.particlefile.stop()
tic = timelib.time()

ndays = 100
timer.run = timer.Timer('Execution', parent=timer.root, start=False)

for d in range(ndays/2):
    day = 2 * d
    print('running %d / %d [time %g s]: %d particles ' % (day+1, ndays, timelib.time()-tic, len(pset)))
    timer.run.start()
    pset.execute(kernel, runtime=delta(days=2), dt=900, verbose_progress=False, recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle})
    timer.run.stop()
    timer.particlefile.start()
    
    if new_write:
        pfile.write_pickle_per_tstep(pset, pset[0].time)
    else:
        pfile.write(pset, pset[0].time)
    
    timer.particlefile.stop()
    
timer.root.stop()
timer.root.print_tree()

