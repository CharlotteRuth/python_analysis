#7/12/18
#Charlotte Christensen
#Plot the fraction of baryons that currently are or ever were in the halo and disk

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.colors as colors
import matplotlib.cm as cm
import matplotlib as mpl
import pynbody
import os, glob, pickle
import os.path
from astropy.io import fits

#entropy_hdu = fits.open('/home/christenc/Data/Sims/cptmarvel.cosmo25cmb/cptmarvel.cosmo25cmb.4096g5HbwK1BH/grp1.allgas.entropy.fits')
#header = entropy_hdu[1].header
#entropy = entropy_hdu[1].data
#entropy_hdu.close()

#halo_mass_z0 = 
#halo_mass_ever

f_bar = 0.16510

basename = prefix + 'cptmarvel.cosmo25cmb/cptmarvel.cosmo25cmb.4096g5HbwK1BH'
tfile = prefix + 'cptmarvel.cosmo25cmb/cptmarvel.cosmo25cmb.4096g5HbwK1BH/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096'
objs_cm = []
f=open(tfile + '.data', 'rb')
while 1:
    try:
        objs_cm.append(pickle.load(f))
    except EOFError:
        break
        
f.close()
objs_pd_cm = pd.DataFrame(objs_cm)

halo_ever = np.empty(len(objs_pd_cm['haloid']))
disk_ever = np.empty(len(objs_pd_cm['haloid']))
smass_accr = np.empty(len(objs_pd_cm['haloid']))

ind = 0
for i in objs_pd_cm['haloid']:
    census_data_file = open(basename + '/grp' + i + '.eject_quant_census.txt')
    temp = census_data_file.readline()
    census_data_file.close()
    halo_ever[ind]  = float(temp.split()[0])
    disk_ever[ind]  = float(temp.split()[1])
    smass_accr[ind]  = float(temp.split()[2])
    ind = + 1

objs_pd_cm['mstar']
objs_pd_cm['mgas']
objs_pd_cm['mISM']
objs_pd_cm['mvir']*f_bar

print("Fraction of cosmic baryons accreted: ",(halo_ever[0] + smass_accr[0])/objs_pd_cm['mvir'][0]/f_bar)
print("Fraction of halo baryons retained: ",np.array(objs_pd_cm['mstar'][0] +  objs_pd_cm['mgas'][0])/(halo_ever[0] + smass_accr[0]))
print("Fraction of halo baryons to disk: ",(disk_ever[0] + smass_accr[0])/(halo_ever[0] + smass_accr[0]))
print("Fraction of disk baryons retained: ",np.array(objs_pd_cm['mstar'][0] +  objs_pd_cm['mISM'][0])/(disk_ever[0] + smass_accr[0]))
