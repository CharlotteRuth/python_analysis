#Charlotte Christensen

#1/29/18
#Generate the properties of dwarf galaxies in the Marvel runs and save them. Based on https://pynbody.github.io/pynbody/tutorials/processing.html

#Run with
#%run /home/christensen/Code/python/python_analysis/bulk_processing_marvel.py
#or, on quirm,
#%run /home/christenc/Code/python/python_analysis/bulk_processing_marvel.py
#ipython --pylab


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
        
    f=open(tfile+'.data','wb')

    for halo_num in halo_nums:
        print(halo_num)
        halo_ahf = h[int(halo_num)] #h.load_copy(int(halo_num))
        
        pynbody.analysis.halo.center(halo_ahf,vel= False)           
        rvir = pynbody.array.SimArray(np.sqrt(np.max(halo_ahf['x'].in_units('kpc')**2 + halo_ahf['y'].in_units('kpc')**2 + halo_ahf['z'].in_units('kpc')**2)),'kpc')
        halo = s[pynbody.filt.Sphere(rvir, (0,0,0))]
        
        stars = halo_ahf.star
        stars.physical_units()

        #For calculating HI and halo mass
        currenttime = halo.properties['time'].in_units('Gyr')
        ISM = halo.gas[pynbody.filt.LowPass('temp', '1e4 K')]
        cool = halo.gas[pynbody.filt.BandPass('temp', '1e4 K', '1e5 K')]
        warm = halo.gas[pynbody.filt.BandPass('temp', '1e5 K', '1e6 K')]
        hot = halo.gas[pynbody.filt.HighPass('temp', '1e6 K')]
        if len(halo.g) != 0 :
            mHI = np.sum(halo.gas[filt.LowPass('temp',2e4)].gas['HI']*halo.gas[filt.LowPass('temp',2e4)].gas['mass'].in_units('Msol'))
            mwarm_coolon = np.sum((warm[pynbody.filt.HighPass('coolontime', str(currenttime)+' Gyr')])['mass'].in_units('Msol')),
            mCool_coolon = np.sum((cool[pynbody.filt.HighPass('coolontime', str(currenttime)+' Gyr')])['mass'].in_units('Msol')),
            mHot_coolon =  np.sum((hot[pynbody.filt.HighPass('coolontime', str(currenttime)+' Gyr')])['mass'].in_units('Msol')),
        else :
            mHI = 0
        
        #For calculating the metallicity

        minmetal = np.amin(s.stars['metals'][np.where(s.stars['metals'] > 0)])
        s.stars['metals'][np.where(s.stars['metals'] == 0)] = minmetal

        #For stellar half-mass raidus
        profile_stellar = pynbody.analysis.profile.Profile(stars,ndin = 2,min = 0, max = np.ceil(rvir), nbins = int(np.ceil(rvir)/0.01))
        index = np.argmin(np.abs(profile_stellar['mass_enc']/max(profile_stellar['mass_enc']) - 0.5))
        r_half = profile_stellar['rbins'].in_units('kpc')[index]

        #For calculating the gas-phase oxygen metallicity
        H2form = 1 #These currently need to be copied over by hand from the paramfile. It would be better to have them read in
        cstar = 0.1
        tempcut = 1e3
        denscut = 0.1 #amu/cc
        deltat = 1e6*3.15569e7 #seconds
        grav =  6.67408e-11 #m3 kg-1 s-2
        halo.g['tdyn'] = 1.0/np.sqrt(4.0*np.pi*grav*halo.g['rho'].in_units('kg m^-3')) #in seconds
        sfgas = halo.g[filt.HighPass('rho',str(denscut)+' m_p cm^-3')]
        if len(sfgas) != 0 :
            sfgas = sfgas[filt.LowPass('temp',str(tempcut)+' K')]
            if len(sfgas) != 0 :
                if H2form:
                    sfgas['sfeff'] = 1.0 - np.exp(-1.0*cstar*deltat/sfgas['tdyn']*2.0*sfgas['H2']/(2.0*sfgas['H2'] + sfgas['HI']))
                else :
                    sfgas['sfeff'] = 1.0 - np.exp(-1.0*cstar*deltat/sfgas['tdyn'])
                oxh_sfr = np.log10(np.sum((10**sfgas['oxh'])*(XSOLO / XSOLH)/15.9*sfgas['sfeff']*sfgas['mass'])/np.sum(sfgas['sfeff']*sfgas['mass']))
            else :
                oxh_sfr = 0
        else :
            oxh_sfr = 0

        coldcut = 1e4
        coldgas = halo.g[filt.LowPass('temp',str(coldcut)+' K')]
        if len(coldgas) != 0:
            oxh_cold = np.log10(np.sum((10**coldgas['oxh'])*(XSOLO / XSOLH)/15.9*coldgas['mass'])/np.sum(coldgas['mass']))
        else :
            oxh_cold = 0
            
        pickle.dump({'haloid': halo_num,
             'z': s.properties['z'],
             'time': currenttime,
             'mvir':  np.sum(halo['mass'].in_units('Msol')),
             'rvir':  rvir,
             'mgas':  np.sum(halo.gas['mass'].in_units('Msol')),
             'mISM':  np.sum(ISM['mass'].in_units('Msol')),
             'mwarm': np.sum(warm['mass'].in_units('Msol')),
             'mCool': np.sum(cool['mass'].in_units('Msol')),
             'mHot': np.sum(hot['mass'].in_units('Msol')),
             'mZISM':  np.sum(ISM['mass'].in_units('Msol')*ISM['metals']),
             'mZwarm': np.sum(warm['mass'].in_units('Msol')*warm['metals']),
             'mZCool': np.sum(cool['mass'].in_units('Msol')*cool['metals']),
             'mZHot': np.sum(hot['mass'].in_units('Msol')*hot['metals']),
             'mZstar': np.sum(stars['mass'].in_units('Msol')*stars['metals']),
             'mwarm_coolon': mwarm_coolon,
             'mCool_coolon': mCool_coolon,
             'mHot_coolon': mHot_coolon,
             'mHI':   mHI,
             'mstar': np.sum(halo.star['mass'].in_units('Msol')),
             'mstarform': np.sum(halo.star['massform'].in_units('Msol')),
             'vmax':  h[int(halo_num)].properties['Vmax'],
             'FeH_mean': np.mean(stars['feh']), #Note that this is the mean of a log
             'FeH_std': np.std(stars['feh']),  #Note that this is a standard deviation of a log
             'Zstellar_mean': np.mean(np.log10(stars['metals']/ZSOLAR)),
             'M_V': -2.5*np.log10(np.sum(10.**(-0.4*halo.star['v_mag']))), # Padova Simple stellar populations, Maiz-Apellaniz 2006 + Bessell 1990
             'M_B': -2.5*np.log10(np.sum(10.**(-0.4*halo.star['b_mag']))), # Padova Simple stellar populations, Maiz-Apellaniz 2006 + Bessell 1990
             'r_half': r_half, #Stellar half-mass radius
             'SFR': float(np.sum(stars[filt.LowPass('age','2e8 yr')].star['massform'].in_units('Msol')))/2e8 , # star formation rate averaged over last 200 Myr
             'sigma_v': np.sqrt(float(np.std(halo.star['vx'].in_units('km s^-1')))**2 + float(np.std(halo.star['vy'].in_units('km s^-1')))**2 + float(np.std(halo.star['vz'].in_units('km s^-1')))**2)/np.sqrt(3), #1-d velocity dispersion calculated as 3-d velocity dispersion/sqrt(3)
             'oxh_sfr': oxh_sfr,
             'oxh_cold': oxh_cold
             },f, pickle.HIGHEST_PROTOCOL)

    f.close()

    
def bulk_processing_read(tfile):

    objs = []
    f=open(tfile + '.data', 'rb')
    while 1:
        try:
            objs.append(pickle.load(f))
        except EOFError:
            break
        
    f.close()

    objs_pd = pd.DataFrame(objs)
    
if __name__ == '__main__':    
    if (socket.gethostname() == "quirm"):
        prefix = '/home/christenc/Data/Sims/'
    else:
        prefix = '/home/christensen/Storage2/UW/MolecH/Cosmo/'
 
    tfile_name = 'h148.cosmo50PLK.3072g3HbwK1BH.004096'
    tfile = prefix + 'h148.cosmo50PLK.3072g/h148.cosmo50PLK.3072g3HbwK1BH/h148.cosmo50PLK.3072g3HbwK1BH.004096/h148.cosmo50PLK.3072g3HbwK1BH.004096'
    halo_nums = ['1','2','3','5','6','9','10','11','12','14','18','23','26','28','31','34','36','42','57','64','77','94','125','160','252','264','271','304']
    bulk_processing(tfile,halo_nums)
    
    tfile_name = 'h329.cosmo50PLK.3072gst5HbwK1BH.004096'
    tfile = prefix + 'h329.cosmo50PLK.3072g/h329.cosmo50PLK.3072gst5HbwK1BH/h329.cosmo50PLK.3072gst5HbwK1BH.004096/h329.cosmo50PLK.3072gst5HbwK1BH.004096'
    halo_nums = ['1','9','32','126','129']
    bulk_processing(tfile,halo_nums)
    
    tfile_name = 'h229.cosmo50PLK.3072gst5HbwK1BH.004096'
    tfile = prefix + 'h229.cosmo50PLK.3072g/h229.cosmo50PLK.3072gst5HbwK1BH/h229.cosmo50PLK.3072gst5HbwK1BH.004096/h229.cosmo50PLK.3072gst5HbwK1BH.004096'
    halo_nums = ['1','2','4','7','17','21','22','27','51','52','70','104','203']
    bulk_processing(tfile,halo_nums)
    
    tfile_name = 'h242.cosmo50PLK.3072gst5HbwK1BH.004096'
    tfile = prefix + 'h242.cosmo50PLK.3072g/h242.cosmo50PLK.3072gst5HbwK1BH/h242.cosmo50PLK.3072gst5HbwK1BH.004096/h242.cosmo50PLK.3072gst5HbwK1BH.004096'      
    halo_nums = ['1','9','11','24','29','30','33','39','40','45','75','76']
    bulk_processing(tfile,halo_nums)
        
    tfile = prefix + 'cptmarvel.cosmo25cmb/cptmarvel.cosmo25cmb.4096g5HbwK1BH/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_name = 'cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096'
    halo_nums = ['1','2','4','5','6','7','10','11','13','14','27']
    bulk_processing(tfile,halo_nums) #yield: 0.02788242
    
    halo_nums = ['1','3','7','8','10','11','12','16','17','18','30','34','36']
    tfile = prefix + 'rogue.cosmo25cmb/rogue.cosmo25cmb.4096g5HbwK1BH/rogue.cosmo25cmb.4096g5HbwK1BH.004096/rogue.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_name = 'rogue.cosmo25cmb.4096g5HbwK1BH.004096'
    bulk_processing(tfile,halo_nums) #yield: 0.02706356

    tfile = prefix + 'elektra.cosmo25cmb/elektra.cosmo25cmb.4096g5HbwK1BH/elektra.cosmo25cmb.4096g5HbwK1BH.004096/elektra.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_name = 'elektra.cosmo25cmb.4096g5HbwK1BH.004096'
    halo_nums = ['1','2','3','4','5','9','10','11','12','17','18','37']
    bulk_processing(tfile,halo_nums) #yield: 0.02773383

    tfile = prefix + 'storm.cosmo25cmb/storm.cosmo25cmb.4096g5HbwK1BH/storm.cosmo25cmb.4096g5HbwK1BH.004096/storm.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_name = 'storm.cosmo25cmb.4096g5HbwK1BH'
    halo_nums = ['1','2','3','4','5','6','7','8','10','11','13','14','16','17','24','34','35','49','109','125','192','208']   
    bulk_processing(tfile,halo_nums) #yield: 0.02758432

    tfile = prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/h986.cosmo50cmb.3072g14HBWK.00512/h986.cosmo50cmb.3072g14HBWK.00512'
    tfile_name = 'h986.cosmo50cmb.3072g14HBWK'
    halo_nums = ['1','2','3','8','15','16']
    #bulk_processing(tfile,halo_nums)
