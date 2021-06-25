
import numpy as np
import pynbody
import pynbody.plot as pp
import pynbody.plot.sph
import pynbody.filt as filt
import pynbody.units as units
import pynbody.analysis.profile as profile
import sys, os, glob, pickle
import socket
import pandas as pd

def bulk_processing(tfile,halo_nums):

    ZSOLAR = 0.0130215
    XSOLO = 0.84E-2 #What pynbody uses
    XSOLH = 0.706
    
    s = pynbody.load(tfile)
    s.physical_units()
    h = s.halos()

    (s.star[filt.LowPass('OxMassFrac',1e-7)])['OxMassFrac'] = 1e-7
    (s.star[filt.LowPass('FeMassFrac',1e-8)])['FeMassFrac'] = 1e-8
    (s.star[filt.LowPass('metals',1e-7*ZSOLAR)])['metals'] = 1e-7*ZSOLAR
    maxHI = np.amax(s.gas['HI'])
        
    f=open(tfile+'ex.data','wb')

    
    for halo_num in halo_nums:
        print(halo_num)
        halo_ahf = h[int(halo_num)] #h.load_copy(int(halo_num))

        halo_ahf = halo_ahf[pynbody.filt.Sphere(5*np.std(np.sqrt(halo_ahf['x'].in_units('kpc')**2 + halo_ahf['y'].in_units('kpc')**2 + halo_ahf['z'].in_units('kpc')**2)), (np.mean(halo_ahf['x'].in_units('kpc')),np.mean(halo_ahf['y'].in_units('kpc')),np.mean(halo_ahf['z'].in_units('kpc'))))] #Filter out outlying particles. This was happening with Sonia. A massive dark matter particle was being included in the halos
        pynbody.analysis.halo.center(halo_ahf,vel = False)
        rvir = pynbody.array.SimArray(np.max(np.sqrt(halo_ahf['x'].in_units('kpc')**2 + halo_ahf['y'].in_units('kpc')**2 + halo_ahf['z'].in_units('kpc')**2)),'kpc')
        halo = halo_ahf
        #halo = s[pynbody.filt.Sphere(rvir, (0,0,0))] #This code can be used to select all material within a virial radius, whether AHF considers it part of the halo or not. It presents problems for satellites as then the host halo material is included.
        
        stars = halo.star
        stars.physical_units()

        #For calculating HI and halo mass
        currenttime = halo.properties['time'].in_units('Gyr')
            
        #For stellar half-mass raidus
        #profile_stellar = pynbody.analysis.profile.Profile(stars,ndin = 2,min = 0, max = np.ceil(rvir), nbins = int(np.ceil(rvir)/0.01))
        #index = np.argmin(np.abs(profile_stellar['mass_enc']/max(profile_stellar['mass_enc']) - 0.5))
        #r_half = profile_stellar['rbins'].in_units('kpc')[index]

        pickle.dump({'haloid': halo_num,
             'z': s.properties['z'],
             'time': currenttime,
             'mvir':  np.sum(halo['mass'].in_units('Msol')),
             'rvir':  rvir,
             'mgas':  np.sum(halo.gas['mass'].in_units('Msol')),
             'mstar': np.sum(halo.star['mass'].in_units('Msol')),
             },f, pickle.HIGHEST_PROTOCOL)

    f.close()

    
def bulk_processing_read(tfile):

    objs = []
    f=open(tfile + '_ex.data', 'rb')
    while 1:
        try:
            objs.append(pickle.load(f))
        except EOFError:
            break
        
    f.close()

    objs_pd = pd.DataFrame(objs)
    
if __name__ == '__main__':    
    if (socket.gethostname() == "quirm.math.grinnell.edu"):
        prefix = '/home/christenc/Data/Sims/'
    else:
        prefix = '/home/christensen/Storage2/UW/MolecH/Cosmo/'

#Sandra
#   tfile_name = 'h148.cosmo50PLK.3072g3HbwK1BH.004096'
#    tfile = prefix + 'h148.cosmo50PLK.3072g/h148.cosmo50PLK.3072g3HbwK1BH/h148.cosmo50PLK.3072g3HbwK1BH.004096/h148.cosmo50PLK.3072g3HbwK1BH.004096'
#    halo_nums = ['1','2','3','5','6','9','10','11','12','14','18','23','26','28','31','34','36','42','57','64','77','94','125','160','252','264','271','304']
#    bulk_processing(tfile,halo_nums)

#Sonia    
    tfile_name = 'h242.cosmo50PLK.3072gst5HbwK1BH.004096'
    tfile = prefix + 'h242.cosmo50PLK.3072g/h242.cosmo50PLK.3072gst5HbwK1BH/h242.cosmo50PLK.3072gst5HbwK1BH.004096/h242.cosmo50PLK.3072gst5HbwK1BH.004096'      
    halo_nums = ['9','11','24','29','30','33','39','40','45','75','76']
    #M200 AHF halo_nums = ['1','8','10','21','26','30','34','38','42','69','70','401','421']
    bulk_processing(tfile,halo_nums)

    
