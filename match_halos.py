# Matches halos that were identified with different runs of AHF


#Run with
#%run /home/christensen/Code/python/python_analysis/match_halos.py
#or, on quirm,
#%run /home/christenc/Code/python/python_analysis/match_halos.py
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
import matplotlib as mpl

# Read in a pickled data file for python 3
def pickle_read(file):

    objs = []
    f=open(file, 'rb')
    while 1:
        try:
            u = pickle._Unpickler(f)
            u.encoding = 'latin1'
            p = u.load()
            objs.append(p)
            #objs.append(pickle.load(f))
        except EOFError:
            break
        
    f.close()

    return pd.DataFrame(objs)


if __name__ == '__main__':    
    if (socket.gethostname() == "quirm.math.grinnell.edu"):
        prefix = '/home/christenc/Data/Sims/'
        outprefix = '/home/christenc/Figures/marvel/'
    else:
        prefix = '/home/christensen/Storage2/UW/MolecH/Cosmo/'
        outprefix = '/home/christensen/Storage2/UW/MolecH/Cosmo/'
    plt.figure(1)
    plt.clf()
    
    dataprefix = '/home/christenc/Code/Datafiles/'        
    f = open(dataprefix+'mstar_vs_mhalo_4Charlotte.txt', 'r')
    fdmdata = []
    for line in f:
        line = line.strip()
        columns = line.split()
        if len(columns) == 12 and columns[0] != 'Volume':
            source = {}
            source['simname'] = columns[0]
            source['halogrp_z0'] = int(columns[1])
            source['halogrp_Mpeak'] = int(columns[2])
            source['Mpeak_snap'] = float(columns[3])
            source['Mpeak'] = float(columns[4])
            source['Mhalo_z0'] = float(columns[5])
            source['Mstar_z0'] = float(columns[6])
            source['Mstar_z0_photo'] = float(columns[7])
            source['Mstar_Mpeak'] = float(columns[8])
            source['Mstar_Mpeak_z0'] = float(columns[9])
            source['Vmag'] = float(columns[10])
            source['type'] = columns[11]            
            fdmdata.append(source)
    f.close()
    fdmdata = pd.DataFrame(fdmdata)
    
#Ruth    
    tfile_name = 'h229.cosmo50PLK.3072gst5HbwK1BH.004096'
    tfile = prefix + 'h229.cosmo50PLK.3072g/h229.cosmo50PLK.3072gst5HbwK1BH/h229.cosmo50PLK.3072gst5HbwK1BH.004096/h229.cosmo50PLK.3072gst5HbwK1BH.004096'
    pointcolor = [0,1,0] #green
    objs_pd_ruth = pickle_read(tfile + '.MAP.data')
    objs_pd_ruth['sim'] = ['h229'] * len(objs_pd_ruth)

#Sonia    
    tfile_name = 'h242.cosmo50PLK.3072gst5HbwK1BH.004096'
    tfile = prefix + 'h242.cosmo50PLK.3072g/h242.cosmo50PLK.3072gst5HbwK1BH/h242.cosmo50PLK.3072gst5HbwK1BH.004096/h242.cosmo50PLK.3072gst5HbwK1BH.004096'      
    pointcolor = [0,0,1] #blue
    objs_pd_sonia = pickle_read(tfile + '.MAP.data')
    objs_pd_sonia['sim'] = ['h242'] * len(objs_pd_sonia)

#Cpt Marvel
    tfile = prefix + 'cptmarvel.cosmo25cmb/cptmarvel.cosmo25cmb.4096g5HbwK1BH/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096'
    objs_pd_cm = pickle_read(tfile + '.MAP.data')
    objs_pd_cm['sim'] = ['cptmarvel'] * len(objs_pd_cm)

#Elektra
    tfile = prefix + 'elektra.cosmo25cmb/elektra.cosmo25cmb.4096g5HbwK1BH/elektra.cosmo25cmb.4096g5HbwK1BH.004096/elektra.cosmo25cmb.4096g5HbwK1BH.004096'
    objs_pd_e = pickle_read(tfile + '.MAP.data')
    objs_pd_e['sim'] = ['elektra'] * len(objs_pd_e)

#Rogue
    tfile = prefix + 'rogue.cosmo25cmb/rogue.cosmo25cmb.4096g5HbwK1BH/rogue.cosmo25cmb.4096g5HbwK1BH.004096/rogue.cosmo25cmb.4096g5HbwK1BH.004096'
    objs_pd_r = pickle_read(tfile + '.MAP.data')
    objs_pd_r['sim'] = ['rogue'] * len(objs_pd_r)
    
#Storm
    tfile = prefix + 'storm.cosmo25cmb/storm.cosmo25cmb.4096g5HbwK1BH/storm.cosmo25cmb.4096g5HbwK1BH.004096/storm.cosmo25cmb.4096g5HbwK1BH.004096'
    objs_pd_s = pickle_read(tfile + '.MAP.data')
    objs_pd_s['sim'] = ['storm'] * len(objs_pd_s)
    
    objs_pd = [objs_pd_ruth]
    #objs_pd = [objs_pd_sandra,objs_pd_ruth,objs_pd_sonia,objs_pd_elena]
    objs_pd = pd.concat(objs_pd)
    smass_tol = 0.5 #fractional tolerance for the stellar mass
    vmass_tol = 0.9

    match_id = {'m200_haloid': 0}
    objs_pd = objs_pd.join(pd.DataFrame(columns=match_id))

    for sim in pd.unique(objs_pd['sim']):
        objs_pd_sim = objs_pd[objs_pd['sim'] == sim]
        for halo in np.sort(objs_pd_sim['haloid']):
            possible_match = (float(objs_pd_sim[objs_pd_sim['haloid'] == halo]['M_star']) < fdmdata[fdmdata['simname']==sim]['Mstar_z0']*(1 + smass_tol)) & (float(objs_pd_sim[objs_pd_sim['haloid'] == halo]['M_star']) > fdmdata[fdmdata['simname']==sim]['Mstar_z0']*(1 - smass_tol))& (float(objs_pd_sim[objs_pd_sim['haloid'] == halo]['mass']) < fdmdata[fdmdata['simname']==sim]['Mhalo_z0']*(1 + vmass_tol)) & (float(objs_pd_sim[objs_pd_sim['haloid'] == halo]['mass']) > fdmdata[fdmdata['simname']==sim]['Mhalo_z0']*(1 - vmass_tol))
            if sum(possible_match) == 0:
                print(sim, halo, 'XXX', float(objs_pd[(objs_pd['sim'] == sim) & (objs_pd['haloid'] == halo)]['M_star']), float(objs_pd[(objs_pd['sim'] == sim) & (objs_pd['haloid'] == halo)]['mass']), 'No Match')
            else:
                #print(fdmdata[fdmdata['simname']==sim][possible_match]['halogrp_z0'])
                arg_best_match = np.argmin(np.abs(fdmdata[fdmdata['simname']==sim][possible_match]['Mstar_z0'] - float(objs_pd_sim[objs_pd_sim['haloid'] == halo]['M_star'])))
                objs_pd.loc[(objs_pd['sim'] == sim) & (objs_pd['haloid'] == halo),'m200_haloid'] = fdmdata.loc[arg_best_match]['halogrp_z0']
                print(sim, halo, fdmdata.loc[arg_best_match]['halogrp_z0'], float(objs_pd[(objs_pd['sim'] == sim) & (objs_pd['haloid'] == halo)]['M_star']), fdmdata.loc[arg_best_match]['Mstar_z0'], float(objs_pd[(objs_pd['sim'] == sim) & (objs_pd['haloid'] == halo)]['mass']), fdmdata.loc[arg_best_match]['Mhalo_z0'])
                fdmdata.loc[arg_best_match,'Mstar_z0'] = 0 #set stellar mass to zero so it can't be used again for another halo

#h229 halo 33 -- stellar mass in data files is different from stat file so no match is found
#242 halo fdm: 26, 30, 42
#cptmarvel -- All correct
#electra -- All correct
#rogue -- All correct!
#storm -- All correct!
