#simple script to convert gadget-type files to tipsy
#this script only converts gas-particles for the use of pure hydrodynamics
#USAGE: python gadgetConvertor.py filename

import pynbody
import sys
import pytipsy
import numpy as np

#file loading
f = pynbody.load(sys.argv[1])
c = sys.argv[1] + '_conv'

#number of particles in each family
ng = len(f.gas)
nd = len(f.dark)
ns = len(f.star)

#get size of system
xmax = f.gas['x'].max()
ymax = f.gas['y'].max()
zmax = f.gas['z'].max()

#get header information
header = {}
header['time'] = 0 #set to zero as IC
header['n'] = ng + nd + ns
header['ndim'] = 3 #set to dimensionality of your system
header['ngas'] = ng
header['ndark'] = nd
header['nstar'] = ns

#gas particle information
catg = {'mass':np.zeros(ng), 'pos':np.zeros((ng,3)), 'x':np.zeros(ng),'y':np.zeros(ng),'z':np.zeros(ng),
        'vel':np.zeros((ng,3)), 'vx':np.zeros(ng), 'vy':np.zeros(ng), 'vz':np.zeros(ng), 'dens':np.zeros(ng), 
        'tempg':np.zeros(ng), 'h':np.zeros(ng), 'zmetal':np.zeros(ng), 'phi':np.zeros(ng)}

for i in range(ng):
    catg['mass'][i] = f.gas['mass'][i]
    catg['pos'][i] = f.gas['pos'][i]
    catg['x'][i] = f.gas['pos'][i][0] - 0.5*xmax #Move to set the middle of the system at (0,0,0)
    catg['y'][i] = f.gas['pos'][i][1] - 0.5*ymax
    catg['z'][i] = f.gas['pos'][i][2] - 0.5*zmax
    catg['vel'][i] = f.gas['vel'][i]
    catg['vx'][i] = f.gas['vel'][i][0]
    catg['vy'][i] = f.gas['vel'][i][1]
    catg['vz'][i] = f.gas['vel'][i][2]
    catg['dens'][i] = f.gas['rho'][i]
    catg['h'][i] = f.gas['smooth'][i]
    catg['tempg'][i] = f.gas['temp'][i]
    catg['zmetal'][i] = 0
    catg['phi'][i] = 0

pytipsy.wtipsy(c, header, catg, {}, {}, STANDARD=True)
