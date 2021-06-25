#Calculate the tidal factor to repliate the work on Bennet+ 2019

#Run with
#%run /home/christensen/Code/python/python_analysis/tidal_index.py
#or, on quirm,
#%run /home/christenc/Code/python/python_analysis/tidal_index.py
#ipython --pylab

import numpy as np
import pynbody
import socket
import matplotlib.pyplot as plt
import sys, os, glob, pickle
import pandas as pd

def tidal_index(tfile,halo_num,pointcolor = [1,0,0]):

    use_smass = 1
    s = pynbody.load(tfile)
    s.physical_units()
    h = s.halos(dummy = True)

    mainhalo = h[halo_num]
    mass = []
    dist = []
    for halo in h:
        if use_smass:
            mass.append(halo.properties['M_star'])
        else:
            mass.append(halo.properties['mass'])
        dist.append(np.sqrt((halo.properties['Xc'] - mainhalo.properties['Xc'])**2 + (halo.properties['Yc'] - mainhalo.properties['Yc'])**2 + (halo.properties['Zc'] - mainhalo.properties['Zc'])**2))

    dist = np.array(dist)/mainhalo.properties['h']/1000 #Convert to Mpc
    mass = np.array(mass)
    dist = dist[1:]
    mass = mass[1:]
    tidal_factor_arr = mass/dist**3

    #cmap, norm = mpl.colors.from_levels_and_colors([0, 2, 5, 6], ['red', 'green', 'blue'])
    #colors = np.asarray([(0, 0, 1, (log10(m)-6)/8) for m in mass])
    #lphas = np.linspace(0.1, 2, len(mass)) #(log10(mass)-6)/8
    if use_smass:
        minlogmass = 3
        maxlogmass = 10
    else:
        minlogmass = 7
        maxlogmass = 13        
    alphas = (np.log10(mass)-minlogmass)/(maxlogmass - minlogmass)
    alphas = np.where(alphas > 0, alphas, alphas*0)
    alphas = np.where(alphas < 1, alphas, alphas*0 + 1)
    rgba_colors = np.zeros((len(dist),4))
    # for red the first column needs to be one
    rgba_colors[:,0] = pointcolor[0]
    rgba_colors[:,1] = pointcolor[1]
    rgba_colors[:,2] = pointcolor[2]
    # the fourth column needs to be your alphas
    rgba_colors[:, 3] = alphas

    pointstyle = plt.scatter(dist,np.log10(tidal_factor_arr) - 10.96,color = rgba_colors,edgecolor = pointcolor)
    plt.xscale('log')
    #plt.yscale('log')
    plt.xlabel(r'Distance [Mpc]')
    plt.ylabel(r'Tidal Index')
    #plt.axis([1e-2,1,1e11,1e13])
    plt.show()
    
    sorted_ind = tidal_factor_arr.argsort()
    tidal_factor  = np.log10(np.sum(tidal_factor_arr[sorted_ind[-5:]])) + -10.96 #as defined in Karachentsev+ 2013
    print(tfile + ': ' + str(tidal_factor))
    print('Mass: ',mass[sorted_ind[-5:]])
    print('Distance: ',dist[sorted_ind[-5:]])
    print('Tidal Index: ',np.log10(tidal_factor_arr[sorted_ind[-5:]]) + -10.96)
    return pointstyle

if __name__ == '__main__':    
    if (socket.gethostname() == "quirm.math.grinnell.edu"):
        prefix = '/home/christenc/Data/Sims/'
    else:
        prefix = '/home/christensen/Storage2/UW/MolecH/Cosmo/'
    plt.figure(1)
    plt.clf()
    
#Sandra 
    tfile_name = 'h148.cosmo50PLK.3072g3HbwK1BH.004096'
    tfile = prefix + 'h148.cosmo50PLK.3072g/h148.cosmo50PLK.3072g3HbwK1BH/h148.cosmo50PLK.3072g3HbwK1BH.004096/h148.cosmo50PLK.3072g3HbwK1BH.004096'
    pointcolor = [1,0,0] #red
    sandra = tidal_index(tfile,1,pointcolor = pointcolor)
    objs_dat = []
    f=open(tfile + '.MAP.data', 'rb')
    while 1:
        try:
            objs_dat.append(pickle.load(f))
        except EOFError:
            break
        
    f.close()
    objs_pd_sandra = pd.DataFrame(objs_dat)
    
#Ruth    
    tfile_name = 'h229.cosmo50PLK.3072gst5HbwK1BH.004096'
    tfile = prefix + 'h229.cosmo50PLK.3072g/h229.cosmo50PLK.3072gst5HbwK1BH/h229.cosmo50PLK.3072gst5HbwK1BH.004096/h229.cosmo50PLK.3072gst5HbwK1BH.004096'
    pointcolor = [0,1,0] #green
    ruth = tidal_index(tfile,1,pointcolor = pointcolor)
    objs_dat = []
    f=open(tfile + '.MAP.data', 'rb')
    while 1:
        try:
            objs_dat.append(pickle.load(f))
        except EOFError:
            break
        
    f.close()
    objs_pd_ruth = pd.DataFrame(objs_dat)

#Sonia    
    tfile_name = 'h242.cosmo50PLK.3072gst5HbwK1BH.004096'
    tfile = prefix + 'h242.cosmo50PLK.3072g/h242.cosmo50PLK.3072gst5HbwK1BH/snapshots/h242.cosmo50PLK.3072gst5HbwK1BH.004096'      
    pointcolor = [0,0,1] #blue
    sonia = tidal_index(tfile,1,pointcolor = pointcolor)
    objs_dat = []
    f=open(tfile + '.MAP.data', 'rb')
    while 1:
        try:
            objs_dat.append(pickle.load(f))
        except EOFError:
            break
        
    f.close()
    objs_pd_sonia = pd.DataFrame(objs_dat)
    
#Elena
    tfile_name = 'h329.cosmo50PLK.3072gst5HbwK1BH.004096'
    tfile = prefix + 'h329.cosmo50PLK.3072g/h329.cosmo50PLK.3072gst5HbwK1BH/snapshots/h329.cosmo50PLK.3072gst5HbwK1BH.004096'
    pointcolor = [0.5,0,0.5] #purple
    elena = tidal_index(tfile,1,pointcolor = pointcolor)
    objs_dat = []
    f=open(tfile + '.MAP.data', 'rb')
    while 1:
        try:
            objs_dat.append(pickle.load(f))
        except EOFError:
            break
        
    f.close()
    objs_pd_elena = pd.DataFrame(objs_dat)

    
    plt.legend([sandra,ruth,sonia,elena],['Sandra','Ruth','Sonia','Elena'])
    plt.axis([4e-2,1,-1,2.5])
    plt.savefig('/home/christenc/Figures/marvel/marvel.smass_tidalFactor.png')

    tidal_factor_sandra = 0.5770758479496578
    tidal_factor_ruth = -0.507694004121344
    tidal_factor_sonia = 1.2806522407890597
    tidal_factor_elena = -0.10882281891372614

    #********** Calculate the number of satellites with M_K < -12 and the fraction that are star forming
    rvir_sandra = objs_pd_sandra[objs_pd_sandra['haloid'] == 1]['Rvir']
    objs_pd_sandra_sat = objs_pd_sandra[objs_pd_sandra['h1dist'] < float(rvir_sandra)]
    objs_pd_sandra_sat = objs_pd_sandra_sat[objs_pd_sandra_sat['h1dist'] > 0]
    Kmag_sandra = -2.5*np.log10(objs_pd_sandra_sat['M_star'].tolist()) + 3.27 #Vega 2mass K band mag for sun
    N_Mv_ltm12_sandra = len(where(objs_pd_sandra_sat['V_mag'] < -12)[0])
    sSFR_sandra = 10**(np.log10(objs_pd_sandra_sat['sSFR'].tolist()))
    frac_ltm12_sf_sandra = float(len(where(sSFR_sandra[where(objs_pd_sandra_sat['V_mag'] < -12)] > 1e-11)[0]))/N_Mv_ltm12_sandra

    rvir_ruth = objs_pd_ruth[objs_pd_ruth['haloid'] == 1]['Rvir']
    objs_pd_ruth_sat = objs_pd_ruth[objs_pd_ruth['h1dist'] < float(rvir_ruth)]
    objs_pd_ruth_sat = objs_pd_ruth_sat[objs_pd_ruth_sat['h1dist'] > 0]
    Kmag_ruth = -2.5*np.log10(objs_pd_ruth_sat['M_star'].tolist()) + 3.27 #Vega 2mass K band mag for sun
    N_Mv_ltm12_ruth = len(where(objs_pd_ruth_sat['V_mag'] < -12)[0])
    sSFR_ruth = 10**(np.log10(objs_pd_ruth_sat['sSFR'].tolist()))
    frac_ltm12_sf_ruth = float(len(where(sSFR_ruth[where(objs_pd_ruth_sat['V_mag'] < -12)] > 1e-11)[0]))/N_Mv_ltm12_ruth

    rvir_sonia = objs_pd_sonia[objs_pd_sonia['haloid'] == 1]['Rvir']
    objs_pd_sonia_sat = objs_pd_sonia[objs_pd_sonia['h1dist'] < float(rvir_sonia)]
    objs_pd_sonia_sat = objs_pd_sonia_sat[objs_pd_sonia_sat['h1dist'] > 0]
    Kmag_sonia = -2.5*np.log10(objs_pd_sonia_sat['M_star'].tolist()) + 3.27 #Vega 2mass K band mag for sun
    N_Mv_ltm12_sonia = len(where(objs_pd_sonia_sat['V_mag'] < -12)[0])
    sSFR_sonia = 10**(np.log10(objs_pd_sonia_sat['sSFR'].tolist()))
    frac_ltm12_sf_sonia = float(len(where(sSFR_sonia[where(objs_pd_sonia_sat['V_mag'] < -12)] > 1e-11)[0]))/N_Mv_ltm12_sonia
    
    rvir_elena = objs_pd_elena[objs_pd_elena['haloid'] == 1]['Rvir']
    objs_pd_elena_sat = objs_pd_elena[objs_pd_elena['h1dist'] < float(rvir_elena)]
    objs_pd_elena_sat = objs_pd_elena_sat[objs_pd_elena_sat['h1dist'] > 0]
    Kmag_elena = -2.5*np.log10(objs_pd_elena_sat['M_star'].tolist()) + 3.27 #Vega 2mass K band mag for sun
    N_Mv_ltm12_elena = len(where(objs_pd_elena_sat['V_mag'] < -12)[0])
    sSFR_elena = 10**(np.log10(objs_pd_elena_sat['sSFR'].tolist()))
    frac_ltm12_sf_elena = float(len(where(sSFR_elena[where(objs_pd_elena_sat['V_mag'] < -12)] > 1e-11)[0]))/N_Mv_ltm12_elena

    frac_ltm12_sf = [frac_ltm12_sf_sandra, frac_ltm12_sf_ruth, frac_ltm12_sf_sonia, frac_ltm12_sf_elena]
    tidal_factor = [tidal_factor_sandra,tidal_factor_ruth,tidal_factor_sonia,tidal_factor_elena]
    N_Mv_ltm8 = [30, 12, 14, 6] #Estimated from Hollis plot for Mv

    plt.clf()
    plt.figure(1,figsize=(plt_width,plt_width)) 
    plt.scatter(tidal_factor,N_Mv_ltm8)
    plt.xlabel(r'Tidal Index')
    plt.ylabel(r'# of satellites: M_V < -8')
    plt.axis([-0.6,3.5,0,35])
    plt.savefig('/home/christenc/Figures/marvel/marvel.numsat_tidalFactor.png')

    plt.clf()
    plt.figure(1,figsize=(plt_width,plt_width)) 
    plt.scatter(tidal_factor,frac_ltm12_sf)
    plt.xlabel(r'Tidal Index')
    plt.ylabel(r'SF fraction: M_V < -8')
    plt.axis([-0.6,3.5,0,1.1])
    plt.savefig('/home/christenc/Figures/marvel/marvel.fracSFsat_tidalFactor.png')    
