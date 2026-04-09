#Charlotte Christensen

#1/29/18
#Generate the properties of dwarf galaxies in the Marvel runs and save them. Based on https://pynbody.github.io/pynbody/tutorials/processing.html

#Run with
#%run /home/christensen/Code/python/python_analysis/bulk_processing_marvel.py
#or, on quirm,
#%run /home/christenc/Code/python/python_analysis/bulk_processing_marvel.py
#or, on dodo:
#%run /home/christenc/Code/python_analysis/bulk_processing_marvel.py 
#ipython --pylab


import numpy as np
import pynbody
import pynbody.filt as filt
import pynbody.units as units
import pynbody.analysis.profile as profile
import sys, os, glob, pickle
import socket
import pandas as pd

def bulk_processing(tfile,halo_nums,add=False):

    ZSOLAR = 0.0130215
    XSOLO = 0.84E-2 #What pynbody uses
    XSOLH = 0.706
    
    s = pynbody.load(tfile)
    s.physical_units()
    h = s.halos()

    print("File Opened")
    #(s.star[filt.LowPass('OxMassFrac',1e-7)])['OxMassFrac'] = 1e-7
    #(s.star[filt.LowPass('FeMassFrac',1e-8)])['FeMassFrac'] = 1e-8
    #(s.star[filt.LowPass('metals',1e-7*ZSOLAR)])['metals'] = 1e-7*ZSOLAR
    #maxHI = np.amax(s.gas['HI'])
    #print("Stars filtered")
    maxHI = 0.7433429956436157
    
    print(add)

    if add:
        if os.path.isfile(tfile + '.data'):
            data = bulk_processing_read(tfile)
        else:
            data = []
    else:
        data = []

    f=open(tfile+'.data','wb')
    
    for halo_num in halo_nums:
        if len(data)!=0:
            print(data.keys())
            print(data['haloid'] == halo_num)
        if len(data)==0 or (sum(data['haloid'] == halo_num) ==0):
            print("Halo Number: ",halo_num)
            #halo_ahf = h[int(halo_num)] #h.load_copy(int(halo_num))
            halo_ahf = h.load_copy(int(halo_num))

            #halo_ahf = halo_ahf[pynbody.filt.Sphere(5*np.std(np.sqrt(halo_ahf['x'].in_units('kpc')**2 + halo_ahf['y'].in_units('kpc')**2 + halo_ahf['z'].in_units('kpc')**2)), 
            #(np.mean(halo_ahf['x'].in_units('kpc')),np.mean(halo_ahf['y'].in_units('kpc')),np.mean(halo_ahf['z'].in_units('kpc'))))] 
            #Filter out outlying particles. I thought this was happening with Sonia but I was wrong. A massive dark matter particle was being included in the halos
            pynbody.analysis.halo.center(halo_ahf,vel = False)
            rvir = pynbody.array.SimArray(np.max(np.sqrt(halo_ahf['x'].in_units('kpc')**2 + halo_ahf['y'].in_units('kpc')**2 + halo_ahf['z'].in_units('kpc')**2)),'kpc')
            halo = halo_ahf
            #halo = s[pynbody.filt.Sphere(rvir, (0,0,0))] 
            #This code can be used to select all material within a virial radius, whether AHF considers it part of the halo or not. It presents problems for satellites as then the host halo material is included.
        
            stars = halo.star
            stars.physical_units()
            (stars[filt.LowPass('OxMassFrac',1e-7)])['OxMassFrac'] = 1e-7
            (stars[filt.LowPass('FeMassFrac',1e-8)])['FeMassFrac'] = 1e-8
            (stars[filt.LowPass('metals',1e-7*ZSOLAR)])['metals'] = 1e-7*ZSOLAR

            #For calculating HI and halo mass
            currenttime = halo.properties['time'].in_units('Gyr')
            if len(halo.g) != 0 :
                ISM = halo.gas[pynbody.filt.LowPass('temp', '1e4 K')]
                cool = halo.gas[pynbody.filt.BandPass('temp', '1e4 K', '1e5 K')]
                warm = halo.gas[pynbody.filt.BandPass('temp', '1e5 K', '1e6 K')]
                hot = halo.gas[pynbody.filt.HighPass('temp', '1e6 K')]
                mISM = np.sum(ISM['mass'].in_units('Msol'))
                mwarm =  np.sum(warm['mass'].in_units('Msol'))
                mCool =  np.sum(cool['mass'].in_units('Msol'))
                mHot =  np.sum(hot['mass'].in_units('Msol'))
                mZISM =   np.sum(ISM['mass'].in_units('Msol')*ISM['metals'])
                mZwarm = np.sum(warm['mass'].in_units('Msol')*warm['metals'])
                mZCool = np.sum(cool['mass'].in_units('Msol')*cool['metals'])
                mZHot = np.sum(hot['mass'].in_units('Msol')*hot['metals'])
                mHI = np.sum(halo.gas[filt.LowPass('temp',2e4)].gas['HI']*halo.gas[filt.LowPass('temp',2e4)].gas['mass'].in_units('Msol'))
#            mwarm_coolon = np.sum((warm[pynbody.filt.HighPass('coolontime', str(currenttime)+' Gyr')])['mass'].in_units('Msol'))
#            mCool_coolon = np.sum((cool[pynbody.filt.HighPass('coolontime', str(currenttime)+' Gyr')])['mass'].in_units('Msol'))
#            mHot_coolon =  np.sum((hot[pynbody.filt.HighPass('coolontime', str(currenttime)+' Gyr')])['mass'].in_units('Msol'))

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
            else :
                mISM = 0
                mwarm = 0
                mcool = 0
                mHot = 0
                mZISM = 0
                mZwarm = 0
                mZCool = 0
                mZHot = 0
                mHI = 0
                mwarm_coolon = 0
                mCool_coolon = 0
                mHot_coolon =  0
                oxh_sfr = 0
                oxh_cold = 0
            
            #For stellar half-mass raidus
            profile_stellar = pynbody.analysis.profile.Profile(stars,ndin = 2,min = 0, max = np.ceil(rvir), nbins = int(np.ceil(rvir)/0.01))
            index = np.argmin(np.abs(profile_stellar['mass_enc']/max(profile_stellar['mass_enc']) - 0.5))
            r_half = profile_stellar['rbins'].in_units('kpc')[index]

            #sfh, bins = pynbody.plot.stars.sfh(halo, filename=None, massform=True, clear=False, legend=False, subplot=False, trange=False, bins=128)
            sfh, bin_edges = np.histogram(stars['tform'], weights = stars['massform'], bins = 128, range = (0,14.2))
            bins = (bin_edges[:-1] - bin_edges[1:])/2

            print('Mvir')
            mvir =  np.sum(halo['mass'].in_units('Msol'))
            print('Mgas')
            mgas = np.sum(halo.gas['mass'].in_units('Msol'))
            print('mZstar')
            mZstar = np.sum(stars['mass'].in_units('Msol')*stars['metals'])
            print('mstar')
            mstar = np.sum(halo.star['mass'].in_units('Msol'))
            print('mstarform')
            mstarform = np.sum(halo.star['massform'].in_units('Msol'))
            print('Vmax')
            vmax = h[int(halo_num)].properties['Vmax']
            print('Fe')
            FeH_mean = np.mean(stars['feh'])
            FeH_std = np.std(stars['feh'])
            print('ZstellarMean')
            Zstellar_mean =  np.mean(np.log10(stars['metals']/ZSOLAR))
            print('M_V')
            M_V =  -2.5*np.log10(np.sum(10.**(-0.4*halo.star['v_mag'])))
            print('M_B')
            M_B =  -2.5*np.log10(np.sum(10.**(-0.4*halo.star['b_mag'])))
            print('SFR')
            sfr  = float(np.sum(stars[filt.LowPass('age','2e8 yr')].star['massform'].in_units('Msol')))/2e8
            print('sigma_v')
            sigma_v =  np.sqrt(float(np.std(halo.star['vx'].in_units('km s^-1')))**2 + float(np.std(halo.star['vy'].in_units('km \
s^-1')))**2 + float(np.std(halo.star['vz'].in_units('km s^-1')))**2)/np.sqrt(3)

            print(np.mean(np.log10(stars['metals']/ZSOLAR)))
            pickle.dump({'haloid': halo_num,
                         'z': s.properties['z'],
                         'time': currenttime,
                         'mvir': mvir,
                         'rvir':  rvir,
                         'mgas':  mgas,
                         'mISM':  mISM,
                         'mwarm': mwarm,
                         'mCool': mCool,
                         'mHot': mHot,
                         'mZISM':  mZISM,
                         'mZwarm': mZwarm,
                         'mZCool': mZCool,
                         'mZHot': mZHot,
                         'mZstar': mZstar,
                         'mwarm_coolon': 0, #mwarm_coolon,
                         'mCool_coolon': 0, #mCool_coolon,
                         'mHot_coolon': 0, #mHot_coolon,
                         'mHI':   mHI,
                         'mstar': mstar,
                         'mstarform': mstarform,
                         'vmax':  vmax,
                         'FeH_mean': FeH_mean, #Note that this is the mean of a log
                         'FeH_std': FeH_std,  #Note that this is a standard deviation of a log
                         'Zstellar_mean': Zstellar_mean,
                         'M_V': M_V, # Padova Simple stellar populations, Maiz-Apellaniz 2006 + Bessell 1990
                         'M_B': M_B, # Padova Simple stellar populations, Maiz-Apellaniz 2006 + Bessell 1990
                         'r_half': r_half, #Stellar half-mass radius
                         'sfh':sfh,
                         'sfhbins':bins,
                         'SFR': sfr , # star formation rate averaged over last 200 Myr
                         'sigma_v': sigma_v, #1-d velocity dispersion calculated as 3-d velocity dispersion/sqrt(3)
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
    return objs_pd
    
if __name__ == '__main__':    
    if (socket.gethostname() == "quirm.math.grinnell.edu"):
        prefix = '/home/christenc/Data/Sims/'
    elif (socket.gethostname() == "dodo.physics.rutgers.edu"):
        prefix = '/home/christenc/REPOSITORY/e12Gals/'
    else:
        prefix = '/home/christensen/Storage2/UW/MolecH/Cosmo/'

#Sandra 
    tfile_name = 'h148.cosmo50PLK.3072g3HbwK1BH.004096'
    tfile = prefix + 'h148.cosmo50PLK.3072g/h148.cosmo50PLK.3072g3HbwK1BH/h148.cosmo50PLK.3072g3HbwK1BH.004096/h148.cosmo50PLK.3072g3HbwK1BH.004096'
    halo_nums = ['1','2','3','5','6','9','10','11','12','14','18','23','26','28','31','34','36','42','57','64','77','94','125','160','252','264','271','304']
    #bulk_processing(tfile,halo_nums)

#Sandra, Hi res
    tfile_name = 'h148.cosmo50PLK.6144g3HbwK1BH'
    M200_AHF_tfile = prefix + 'h148.cosmo50PLK.6144g3HbwK1BH/h148.cosmo50PLK.6144g3HbwK1BH.004096/ahf_200/h148.cosmo50PLK.6144g3HbwK1BH.004096'
    M200_AHF_halo_nums = ['2',    '3',    '4',    '5',    '6',    '7',    '9',   '11', '12',   '13',   '16', '17',   '19',   '23',   '24',   '25',   '26',   '29',   '30',   '32',   '34',   '35', '36',   '39',   '44',   '46',   '49',   '50',   '51',   '59',   '64',   '71',   '72','73',   '85',   '86',  '133',  '170',  '181',  '195',  '200',  '250',  '262',  '289', '319',  '348',  '392',  '426',  '454',  '532',  '644',  '741', '1007', '1109', '1367', '1394', '1918', '3252', '5768']
    #bulk_processing(tfile,halo_nums)
    bulk_processing(M200_AHF_tfile,M200_AHF_halo_nums)

#Ruth    
    tfile_name = 'h229.cosmo50PLK.3072gst5HbwK1BH.004096'
    tfile = prefix + 'h229.cosmo50PLK.3072g/h229.cosmo50PLK.3072gst5HbwK1BH/h229.cosmo50PLK.3072gst5HbwK1BH.004096/h229.cosmo50PLK.3072gst5HbwK1BH.004096'
    halo_nums = ['1','2','4','7','17','21','22','27','51','52','70','104','203']
    #M200_AHF_tfile = prefix + 'h229.cosmo50PLK.3072g/h229.cosmo50PLK.3072gst5HbwK1BH/h229.cosmo50PLK.3072gst5HbwK1BH.004096/ahf_200/h229.cosmo50PLK.3072gst5HbwK1BH.004096'    
    #M200_AHF_halo_nums = ['1','2','3','6','14','15','18','20','22','33','47','48','49','52','57','62','64','89','92','193','528','886','1302','5313']
    #bulk_processing(M200_AHF_tfile,M200_AHF_halo_nums)

#Sonia    
    tfile_name = 'h242.cosmo50PLK.3072gst5HbwK1BH.004096'
    tfile = prefix + 'h242.cosmo50PLK.3072g/h242.cosmo50PLK.3072gst5HbwK1BH/h242.cosmo50PLK.3072gst5HbwK1BH.004096/h242.cosmo50PLK.3072gst5HbwK1BH.004096'      
    halo_nums = ['1','9','11','24','29','30','33','39','40','45','75','76']
    #M200_AHF_tfile = prefix + 'h242.cosmo50PLK.3072g/h242.cosmo50PLK.3072gst5HbwK1BH/h242.cosmo50PLK.3072gst5HbwK1BH.004096/ahf_200/h242.cosmo50PLK.3072gst5HbwK1BH.004096' 
    #M200_AHF_halo_nums = ['1','8','10','21','26','30','34','38','42','44','63','69','70','81','138','192','401','421','1633','5495','8816']
    #bulk_processing(M200_AHF_tfile,M200_AHF_halo_nums)

#Elena
    tfile_name = 'h329.cosmo50PLK.3072gst5HbwK1BH.004096'
    tfile = prefix + 'h329.cosmo50PLK.3072g/h329.cosmo50PLK.3072gst5HbwK1BH/h329.cosmo50PLK.3072gst5HbwK1BH.004096/h329.cosmo50PLK.3072gst5HbwK1BH.004096'
    halo_nums = ['1','9','32','126','129']
    #bulk_processing(tfile,halo_nums)
    
#Elena, Hi res
    tfile_name = 'h329.cosmo50PLK.6144g5HbwK1BH'
    M200_AHF_tfile = prefix + 'h329.cosmo50PLK.6144g5HbwK1BH/h329.cosmo50PLK.6144g5HbwK1BH.004096/ahf_200/h329.cosmo50PLK.6144g5HbwK1BH.004096'
    M200_AHF_halo_nums = ['1',  '7',  '12',  '15',  '18',  '28',  '41',  '64', '307', '406', '543', '887']
    #bulk_processing(M200_AHF_tfile,M200_AHF_halo_nums)

    tfile = prefix + 'cptmarvel.cosmo25cmb/cptmarvel.cosmo25cmb.4096g5HbwK1BH/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_name = 'cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096'
    halo_nums = ['1','2','4','5','6','7','10','11','13','14','27']
    #bulk_processing(tfile,halo_nums) #yield: 0.02788242
    
    halo_nums = ['1','3','7','8','10','11','12','16','17','18','30','34','36']
    tfile = prefix + 'rogue.cosmo25cmb/rogue.cosmo25cmb.4096g5HbwK1BH/rogue.cosmo25cmb.4096g5HbwK1BH.004096/rogue.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_name = 'rogue.cosmo25cmb.4096g5HbwK1BH.004096'
    #bulk_processing(tfile,halo_nums) #yield: 0.02706356

    tfile = prefix + 'elektra.cosmo25cmb/elektra.cosmo25cmb.4096g5HbwK1BH/elektra.cosmo25cmb.4096g5HbwK1BH.004096/elektra.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_name = 'elektra.cosmo25cmb.4096g5HbwK1BH.004096'
    halo_nums = ['1','2','3','4','5','9','10','11','12','17','18','37']
    #bulk_processing(tfile,halo_nums) #yield: 0.02773383

    tfile = prefix + 'storm.cosmo25cmb/storm.cosmo25cmb.4096g5HbwK1BH/storm.cosmo25cmb.4096g5HbwK1BH.004096/storm.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_name = 'storm.cosmo25cmb.4096g5HbwK1BH'
    halo_nums = ['1','2','3','4','5','6','7','8','10','11','13','14','16','17','24','34','35','49','109','125','192','208']   
    #bulk_processing(tfile,halo_nums) #yield: 0.02758432

    tfile = prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/h986.cosmo50cmb.3072g14HBWK.00512/h986.cosmo50cmb.3072g14HBWK.00512'
    tfile_name = 'h986.cosmo50cmb.3072g14HBWK'
    halo_nums = ['1','2','3','8','15','16']
    #bulk_processing(tfile,halo_nums)
