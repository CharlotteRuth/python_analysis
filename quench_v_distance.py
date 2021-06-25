#Calculate the relationship between quenching time, stellar mass, and distance to host, replicating Weisz+ 2015, Figure 1

#Run with
#%run /home/christensen/Code/python/python_analysis/quench_v_distance.py
#or, on quirm,
#%run /home/christenc/Code/python/python_analysis/quench_v_distance.py
#ipython --pylab

import numpy as np
import pynbody
import socket
import matplotlib.pyplot as plt
import pandas as pd
import sys, os, glob, pickle
from scipy.interpolate import interp1d
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
#import distance_to_nearest_host


#Calculates the distance to the nearest massive galaxy

def distance_to_nearest_host(data):
    distances = []
    hostrvirs = []
    for i in range(len(data)):
        s = data['sim'].tolist()[i]
       
        if s=='h148' or s=='h229' or s=='h242' or s=='h329': # if sat simulation, find distance to halo 1
            h1dist = data['h1dist'].tolist()[i]*0.6776942783267969
            distances.append(h1dist)
       
            h1rvir = data['Rvir'][(data.sim==s) & (data.haloid==1)].tolist()[0]*0.6776942783267969
            hostrvirs.append(h1rvir)
           
        else: # if field simulation, find distance to nearest massive DM halo (currently > 0.5e12.5 Msol)
            if s=='cptmarvel':
                path = '/home/akinshol/Data/Sims/cptmarvel.cosmo25cmb.4096g5HbwK1BH/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096'

            if s=='elektra':
                path = '/home/akinshol/Data/Sims/elektra.cosmo25cmb.4096g5HbwK1BH/elektra.cosmo25cmb.4096g5HbwK1BH.004096.dir/elektra.cosmo25cmb.4096g5HbwK1BH.004096'

            if s=='rogue':
                path = '/home/akinshol/Data/Sims/rogue.cosmo25cmb.4096g5HbwK1BH/rogue.cosmo25cmb.4096g5HbwK1BH.004096.dir/rogue.cosmo25cmb.4096g5HbwK1BH.004096'

            if s=='storm':
                path = '/home/akinshol/Data/Sims/storm.cosmo25cmb.4096g5HbwK1BH/storm.cosmo25cmb.4096g5HbwK1BH.004096/storm.cosmo25cmb.4096g5HbwK1BH.004096'
           
            coords = []
            with open(path+'.coords','rb') as f:
                while True:
                    try:
                        coords.append(pickle.load(f)) #,encoding='latin1'))
                    except EOFError:
                        break
            coords = pd.DataFrame(coords)
           
            threshold = 5*10**(11) # this threshold can be adjusted,
            # i tried to pick something similar to the virial masses of the host in the JL simulations
           
            coords = coords[coords.mass > threshold]
           
            halocoords = np.array([data['Xc'].tolist()[i],data['Yc'].tolist()[i],data['Zc'].tolist()[i]])

            x = np.array(coords['Xc'])
            y = np.array(coords['Yc'])
            z = np.array(coords['Zc'])
            Rvir = np.array(coords['Rv'])

           
            c = np.array([x,y,z])
            c = np.transpose(c)
            dist = np.sqrt(np.sum((halocoords-c)**2, axis=1))*0.6776942783267969
            distances.append(np.min(dist))
            hostrvirs.append(Rvir[np.argmin(dist)]*0.6776942783267969)
           
    return np.array(distances),np.array(hostrvirs)



if __name__ == '__main__':    
    if (socket.gethostname() == "quirm.math.grinnell.edu"):
        prefix = '/home/christenc/Data/Sims/'
        outprefix = '/home/christenc/Figures/marvel/'
    else:
        prefix = '/home/christensen/Storage2/UW/MolecH/Cosmo/'
        outprefix = '/home/christensen/Storage2/UW/MolecH/Cosmo/'
    plt.figure(1)
    plt.clf()

    presentation = True #False
    if presentation:
        outbase = outprefix + 'marvel_pres_'
        plt.style.use(['default','/home/christenc/.config/matplotlib/presentation.mplstyle'])
        plt_width = 8 #inches
        aspect_ratio = 3.0/4.0
        legendsize = 16
        dpi = 100
        markersize = 40
    else:
        outbase = outprefix + 'marvel'
        plt.style.use(['default','/home/christenc/.config/matplotlib/article.mplstyle'])
        plt_width = 3.5 #inches
        aspect_ratio = 3.0/4.0
        legendsize = 5
        dpi = 300
        markersize = 12
    
#Sandra 
    tfile_name = 'h148.cosmo50PLK.3072g3HbwK1BH.004096'
    tfile = prefix + 'h148.cosmo50PLK.3072g/h148.cosmo50PLK.3072g3HbwK1BH/h148.cosmo50PLK.3072g3HbwK1BH.004096/h148.cosmo50PLK.3072g3HbwK1BH.004096'
    pointcolor = [1,0,0] #red
    objs_dat = []
    f=open(tfile + '.MAP.data', 'rb')
    while 1:
        try:
            objs_dat.append(pickle.load(f))
        except EOFError:
            break        
    f.close()
    objs_pd_sandra = pd.DataFrame(objs_dat)
    objs_pd_sandra['sim'] = ['h148'] * len(objs_pd_sandra)
    
#Ruth    
    tfile_name = 'h229.cosmo50PLK.3072gst5HbwK1BH.004096'
    tfile = prefix + 'h229.cosmo50PLK.3072g/h229.cosmo50PLK.3072gst5HbwK1BH/h229.cosmo50PLK.3072gst5HbwK1BH.004096/h229.cosmo50PLK.3072gst5HbwK1BH.004096'
    pointcolor = [0,1,0] #green
    objs_dat = []
    f=open(tfile + '.MAP.data', 'rb')
    while 1:
        try:
            objs_dat.append(pickle.load(f))
        except EOFError:
            break        
    f.close()
    objs_pd_ruth = pd.DataFrame(objs_dat)
    objs_pd_ruth['sim'] = ['h229'] * len(objs_pd_ruth)

#Sonia    
    tfile_name = 'h242.cosmo50PLK.3072gst5HbwK1BH.004096'
    tfile = prefix + 'h242.cosmo50PLK.3072g/h242.cosmo50PLK.3072gst5HbwK1BH/snapshots/h242.cosmo50PLK.3072gst5HbwK1BH.004096'      
    pointcolor = [0,0,1] #blue
    objs_dat = []
    f=open(tfile + '.MAP.data', 'rb')
    while 1:
        try:
            objs_dat.append(pickle.load(f))
        except EOFError:
            break        
    f.close()
    objs_pd_sonia = pd.DataFrame(objs_dat)
    objs_pd_sonia['sim'] = ['h242'] * len(objs_pd_sonia)
    
#Elena
    tfile_name = 'h329.cosmo50PLK.3072gst5HbwK1BH.004096'
    tfile = prefix + 'h329.cosmo50PLK.3072g/h329.cosmo50PLK.3072gst5HbwK1BH/snapshots/h329.cosmo50PLK.3072gst5HbwK1BH.004096'
    pointcolor = [0.5,0,0.5] #purple
    objs_dat = []
    f=open(tfile + '.MAP.data', 'rb')
    while 1:
        try:
            objs_dat.append(pickle.load(f))
        except EOFError:
            break        
    f.close()
    objs_pd_elena = pd.DataFrame(objs_dat)
    objs_pd_elena['sim'] = ['h329'] * len(objs_pd_elena)

    #Cpt Marvel
    tfile = prefix + 'cptmarvel.cosmo25cmb/cptmarvel.cosmo25cmb.4096g5HbwK1BH/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096'
    objs_cm = []
    f=open(tfile + '.MAP.data', 'rb')
    while 1:
        try:
            objs_cm.append(pickle.load(f))
        except EOFError:
            break        
    f.close()
    objs_pd_cm = pd.DataFrame(objs_cm)
    objs_pd_cm['sim'] = ['cptmarvel'] * len(objs_pd_cm)

#Elektra
    tfile = prefix + 'elektra.cosmo25cmb/elektra.cosmo25cmb.4096g5HbwK1BH/elektra.cosmo25cmb.4096g5HbwK1BH.004096/elektra.cosmo25cmb.4096g5HbwK1BH.004096'
    objs_e = []
    f=open(tfile + '.MAP.data', 'rb')
    while 1:
        try:
            objs_e.append(pickle.load(f))
        except EOFError:
            break        
    f.close()
    objs_pd_e = pd.DataFrame(objs_e)
    objs_pd_e['sim'] = ['elektra'] * len(objs_pd_e)

#Rogue
    tfile = prefix + 'rogue.cosmo25cmb/rogue.cosmo25cmb.4096g5HbwK1BH/rogue.cosmo25cmb.4096g5HbwK1BH.004096/rogue.cosmo25cmb.4096g5HbwK1BH.004096'
    objs_r = []
    f=open(tfile + '.MAP.data', 'rb')
    while 1:
        try:
            objs_r.append(pickle.load(f))
        except EOFError:
            break        
    f.close()
    objs_pd_r = pd.DataFrame(objs_r)
    objs_pd_r['sim'] = ['rogue'] * len(objs_pd_r)

#Storm
    tfile = prefix + 'storm.cosmo25cmb/storm.cosmo25cmb.4096g5HbwK1BH/storm.cosmo25cmb.4096g5HbwK1BH.004096/storm.cosmo25cmb.4096g5HbwK1BH.004096'
    objs_s = []
    f=open(tfile + '.MAP.data', 'rb')
    while 1:
        try:
            objs_s.append(pickle.load(f))
        except EOFError:
            break        
    f.close()
    objs_pd_s = pd.DataFrame(objs_s)
    objs_pd_s['sim'] = ['storm'] * len(objs_pd_s)

    objs_pd = [objs_pd_sandra,objs_pd_ruth,objs_pd_sonia,objs_pd_elena,objs_pd_s,objs_pd_r,objs_pd_e,objs_pd_cm]
    #objs_pd = [objs_pd_sandra,objs_pd_ruth,objs_pd_sonia,objs_pd_elena]
    objs_pd = pd.concat(objs_pd)
    tau90 = np.empty(len(objs_pd))
    distances, rvir = distance_to_nearest_host(objs_pd)

    
    plt.clf()
    plt.figure(1)
    ind = 0
    for index, row in objs_pd.iterrows():
        xarr = row['sfhbins'][1:] - (row['sfhbins'][1] - row['sfhbins'][0])
        yarr = np.cumsum(row['sfh'])/max(np.cumsum(row['sfh']))
        plt.clf()
        plt.plot(xarr,yarr)
        plt.plot([0,14],[0.9,0.9])
        #if np.isnan(interp(0.9)):
            #print(row)
            #print(row['sfhbins'],row['sfh'])
        if (yarr[0] >= 0.9):
            tau90[ind] = xarr[0]
        else:
            interp = interp1d(yarr, xarr) #, kind='cubic')
            print(ind,interp(0.9))
            tau90[ind] = float(interp(0.9))
        ind = ind + 1

    cmx = plt.get_cmap("viridis") 
    cNorm  = colors.Normalize(vmin=0, vmax = 14)
    #scalarMap = cm.ScalarMappable(norm=cNorm, cmap=cmx)
    plt.clf()
    fig1 = plt.figure(1,figsize = (plt_width,plt_width*aspect_ratio))
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax1 = fig1.add_subplot(gs[0])
    ax1sub = fig1.add_subplot(gs[1])    
    q = ax1.scatter(distances[objs_pd['SFR'] < 1e-11],objs_pd['M_star'][objs_pd['SFR'] < 1e-11],s = (objs_pd['Rvir'][objs_pd['SFR'] < 1e-11]*2).tolist(), c = tau90[objs_pd['SFR'] < 1e-11], cmap = cmx, norm = cNorm,edgecolor = 'k',marker = 'D')
    sf = ax1.scatter(distances[objs_pd['SFR'] >= 1e-11],objs_pd['M_star'][objs_pd['SFR'] >= 1e-11],s = (objs_pd['Rvir'][objs_pd['SFR'] >= 1e-11]*2).tolist(), c = tau90[objs_pd['SFR'] >= 1e-11], cmap = cmx, norm = cNorm,edgecolor = 'k')
    #plt.scatter(objs_pd_e['h1dist'],objs_pd_e['M_star'])
    ax1.legend([q,sf],['Quenched','Star forming'],scatterpoints = 1,facecolor = 'white',loc = 3,framealpha = 0,frameon = False)    
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.axis([10, 5e3, 1e2, 1e11])
    ax1.set_ylabel(r'M$_*$/M$_\odot$')
    ax1.set_xlabel(r'Distance to massive galaxy (kpc)')
    cb = mpl.colorbar.ColorbarBase(ax1sub, cmap=cmx, norm=cNorm)
    cb.set_label(r"$\tau_{90}$ (Gyr)")

    plt.savefig(outbase + '.distance_smass_t90.png')    
