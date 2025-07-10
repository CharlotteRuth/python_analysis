#9/23/18
#Charlotte Christensen
#Identify all gas that was ever in the halo/disk of a specific galaxy and its progenitor

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.colors as colors
import matplotlib.cm as cm
import matplotlib as mpl
import pynbody
import os, glob, pickle
import os.path
import socket

if (socket.gethostname() == "quirm.math.grinnell.edu"):
    prefix = '/home/christenc/Data/Sims/'
    outprefix = '/home/christenc/Figures/marvel/'
    dataprefix = '/home/christenc/Code/Datafiles/'        
else:
    prefix = '/home/christensen/Storage2/UW/MolecH/Cosmo/'
    outprefix = '/home/christensen/Storage2/UW/MolecH/Cosmo/'
    dataprefix = '/home/christensen/Code/Datafiles/'

stemname = 'cptmarvel.cosmo25cmb'
basename = 'cptmarvel.cosmo25cmb.4096g5HbwK1BH'
tfile = prefix + 'cptmarvel.cosmo25cmb/cptmarvel.cosmo25cmb.4096g5HbwK1BH/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096'

grp = '1'
laststep = 4096
currentstep = laststep

#Read in the merger file.
#Then create a dictionary where the keys are the timesteps and the values are the halo ids associated with those steps
timestepsDict = {}
timestepsDict[currentstep] = []
timestepsDict[currentstep].append(grp)

f = open(prefix + stemname + '/' + basename + '/' + basename + '.grp' + grp + '.test.halo_step.out', 'r')
while 1:
    line = f.readline()
    if not line:
        break
    tokens = line.split()
    if tokens[0] == 'step':
        currentstep = int(tokens[1])
        timestepsDict[currentstep] = []
    else:
        timestepsDict[currentstep].append(tokens[0])    
f.close()

#Make an empty dictionary that will hold the gas and star particles ever part of the halo or its progenitor
#Iord
#Accreted as star/gas
#Timestep first identified as being in halo
#Mass when first identified as being in halo
#Question: How will I distinguish between satellite or not? Do I want to include satellites? What about gas that is idenitfied as a satellite but is actually a disk clump?
allgas = {}

timesteps = timestepsDict.keys()
timesteps.sort()
for timestep in timesteps:
    print('{:0>6}'.format(timestep),timestepsDict[timestep])
    tfile = prefix + stemname + '/' + basename + '/' + basename + '.' + '{:0>6}'.format(timestep) + '/' + basename + '.' + '{:0>6}'.format(timestep)
    if os.path.isfile(tfile):
        s = pynbody.load(tfile)
        h = s.halos()
        for halo in timestepsDict[timestep]:
            gas = h[int(halo)].gas
            stars = h[int(halo)].star
            print(len(gas))
            gas['accrtime'] = h[int(halo)].properties['time'].in_units('Gyr')
            if len(allgas) == 0:
                #allgas = gas
                newgas = gas
            else :
                newgas_iord = numpy.isin(gas['iord'],allgas['iord'],invert=True)
                newgas = gas[newgas]
                #allgas = [allgas,gas[newgas]]
            for i in range(len(newgas)):
                     #print(gas[i]['iord'])
                     #print(int(newgas[i]['iord']))
                     allgas[int(newgas[i]['iord'])] = {}
                     #print(float(gas[i]['mass'].in_units('Msol')))
                     allgas[int(newgas[i]['iord'])]['accrmass'] = float(gas[i]['mass'].in_units('Msol'))
                     #print(float(gas[i]['accrtime']))
                     allgas[int(newgas[i]['iord'])]['accrtime'] = float(gas[i]['accrtime'])
                     allgas[int(newgas[i]['iord'])]['accrtype'] = 'gas'
    else:
        print('File does not exist')

        #Check for unique gas particles. Add those particles to the list along with their mass and the current time
        #Check for star particles whose iord OR progenitor iord never appear on the list. Add this to the list along with their mass and current time. Flag that they are accreted as a star, I guess.
