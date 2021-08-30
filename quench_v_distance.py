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
import matplotlib as mpl
#import pickle_read
#import distance_to_nearest_host

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

# Matches the halos in the objs_pd file with the halo information from Ferah (m200)
# see "match_halos.py" for development
def match_halos(objs_pd, fdmdata):
    smass_tol = 0.5 #fractional tolerance for the stellar mass
    vmass_tol = 0.9

    match_id = {'m200_haloid': 0}
    objs_pd = objs_pd.join(pd.DataFrame(columns=match_id))
    
    for sim in pd.unique(objs_pd['sim']):
        objs_pd_sim = objs_pd[objs_pd['sim'] == sim].copy()
        for halo in np.sort(objs_pd_sim['haloid']):
            possible_match = (float(objs_pd_sim[objs_pd_sim['haloid'] == halo]['M_star']) < fdmdata[fdmdata['simname']==sim]['Mstar_z0']*(1 + smass_tol)) & (float(objs_pd_sim[objs_pd_sim['haloid'] == halo]['M_star']) > fdmdata[fdmdata['simname']==sim]['Mstar_z0']*(1 - smass_tol))& (float(objs_pd_sim[objs_pd_sim['haloid'] == halo]['mass']) < fdmdata[fdmdata['simname']==sim]['Mhalo_z0']*(1 + vmass_tol)) & (float(objs_pd_sim[objs_pd_sim['haloid'] == halo]['mass']) > fdmdata[fdmdata['simname']==sim]['Mhalo_z0']*(1 - vmass_tol))
            if sum(possible_match) == 0:
                print(sim, halo, 'XXX', float(objs_pd[(objs_pd['sim'] == sim) & (objs_pd['haloid'] == halo)]['M_star']), float(objs_pd[(objs_pd['sim'] == sim) & (objs_pd['haloid'] == halo)]['mass']), 'No Match')
            else:
                #print(fdmdata[fdmdata['simname']==sim][possible_match]['halogrp_z0'])
                arg_best_match = np.argmin(np.abs(fdmdata[fdmdata['simname']==sim][possible_match]['Mstar_z0'] - float(objs_pd_sim[objs_pd_sim['haloid'] == halo]['M_star'])))
                objs_pd.loc[(objs_pd['sim'] == sim) & (objs_pd['haloid'] == halo),'m200_haloid'] = fdmdata.loc[arg_best_match]['halogrp_z0']
                #print(sim, halo, fdmdata.loc[arg_best_match]['halogrp_z0'], float(objs_pd[(objs_pd['sim'] == sim) & (objs_pd['haloid'] == halo)]['M_star']), fdmdata.loc[arg_best_match]['Mstar_z0'], float(objs_pd[(objs_pd['sim'] == sim) & (objs_pd['haloid'] == halo)]['mass']), fdmdata.loc[arg_best_match]['Mhalo_z0'])
                fdmdata.loc[arg_best_match,'Mstar_z0'] = 0 #set stellar mass to zero so it

    return objs_pd
                
#Calculates the distance to the nearest massive galaxy
def distance_to_nearest_host(tfiles,data):
    distances = []
    hostrvirs = []
    min_massiveHalo = 10**11.5
    sprev = ''
    tfile_it = -1    
    for i in range(len(data)):
        s = data['sim'].tolist()[i]

        print(s,data['haloid'].tolist()[i])       
        if s=='h148' or s=='h229' or s=='h242' or s=='h329': # if sat simulation, find distance to halo 1
            if s != sprev:
                tfile_it = tfile_it + 1
                sprev = s
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

            if s != sprev:
                sim = pynbody.load(tfiles[tfiles_it])
                tfile_it = tfile_it + 1
                sprev = s                
                h_dummy = sim.halos(dummy = True)
                loc = []
                rvir = []
                for AHFhalo in h_dummy:
                    properties = AHFhalo.properties            
                    if (properties['mass'] > min_massiveHalo):
#                        print('Halo id: ',properties['halo_id'])
                        loc.append(np.array([properties['Xc']/properties['h'], properties['Yc']/properties['h'], properties['Zc']/properties['h']]))
                        rvir.append(properties['Rvir'])
                
                loc = np.array(loc)
                rvir = np.array(rvir)

            properties = h_dummy[int(data['haloid'].tolist()[i])].properties
            distances.append(min(((properties['Xc'] - loc[:,0])**2 + (properties['Yc'] - loc[:,1])**2 + (properties['Zc'] - loc[:,2])**2)**(0.5)))
            minind = np.where((((properties['Xc'] - loc[:,0])**2 + (properties['Yc'] - loc[:,1])**2 + (properties['Zc'] - loc[:,2])**2)**(0.5)) == distances[-1])
            hostrvirs.append(rvir[minind]*0.6776942783267969)

    return np.array(distances),np.array(hostrvirs)



if __name__ == '__main__':    
    if (socket.gethostname() == "quirm.math.grinnell.edu"):
        prefix = '/home/christenc/Data/Sims/'
        outprefix = '/home/christenc/Figures/marvel/marvelJL'
    else:
        prefix = '/home/christensen/Storage2/UW/MolecH/Cosmo/'
        outprefix = '/home/christensen/Storage2/UW/MolecH/Cosmo/marvelJL'
    plt.figure(1)
    plt.clf()

    presentation = False
    if presentation:
        outbase = outprefix + '_pres_'
        plt.style.use(['default','/home/christenc/.config/matplotlib/presentation.mplstyle'])
        plt_width = 8 #inches
        aspect_ratio = 3.0/4.0
        legendsize = 16
        dpi = 200
        markersize = 100
        ms_scale = 1
        lw = mpl.rcParams['lines.linewidth'] - 1
        edgewidth = 2
    else:
        outbase = outprefix #+ 'marvel'
        plt.style.use(['default','/home/christenc/.config/matplotlib/article.mplstyle'])
        plt_width = 3.5 #inches
        aspect_ratio = 3.0/4.0
        legendsize = 5
        dpi = 300
        markersize = 25
        ms_scale = 0.25
        lw = mpl.rcParams['lines.linewidth']
        edgewidth = 0.5
        
    if (socket.gethostname() == "ozma.grinnell.edu"):
        dataprefix = '/home/christensen/Code/Datafiles/' 
    else:
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

    tfile_base_cm = 'cptmarvel.cosmo25cmb.4096g5HbwK1BH'
    tfile_cm = prefix + 'cptmarvel.cosmo25cmb/cptmarvel.cosmo25cmb.4096g5HbwK1BH/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096' #'cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096'

    tfile_r = prefix + 'rogue.cosmo25cmb/rogue.cosmo25cmb.4096g5HbwK1BH/rogue.cosmo25cmb.4096g5HbwK1BH.004096/rogue.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_base_r = 'rogue.cosmo25cmb.4096g5HbwK1BH'

    tfile_e = prefix + 'elektra.cosmo25cmb/elektra.cosmo25cmb.4096g5HbwK1BH/elektra.cosmo25cmb.4096g5HbwK1BH.004096/elektra.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_base_e = 'elektra.cosmo25cmb.4096g5HbwK1BH'

    tfile_s = prefix + 'storm.cosmo25cmb/storm.cosmo25cmb.4096g5HbwK1BH/storm.cosmo25cmb.4096g5HbwK1BH.004096/storm.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_base_s = 'storm.cosmo25cmb.4096g5HbwK1BH'
    
    tfile_base_1 = 'h148.cosmo50PLK.3072g3HbwK1BH'
    tfile_1 = prefix + 'h148.cosmo50PLK.3072g/h148.cosmo50PLK.3072g3HbwK1BH/snapshots_200bkgdens/h148.cosmo50PLK.3072g3HbwK1BH.004096' #
    tfile_1 = prefix + 'h148.cosmo50PLK.3072g/h148.cosmo50PLK.3072g3HbwK1BH/h148.cosmo50PLK.3072g3HbwK1BH.004096/h148.cosmo50PLK.3072g3HbwK1BH.004096'
    
    tfile_base_1hr = 'h148.cosmo50PLK.6144g3HbwK1BH/'
    tfile_1hr = prefix + 'h148.cosmo50PLK.6144g/h148.cosmo50PLK.6144g3HbwK1BH/h148.cosmo50PLK.6144g3HbwK1BH.004096/ahf_200/h148.cosmo50PLK.6144g3HbwK1BH.004096'
    
    tfile_base_2 = 'h229.cosmo50PLK.3072gst5HbwK1BH'
    tfile_2 = prefix + 'h229.cosmo50PLK.3072g/h229.cosmo50PLK.3072gst5HbwK1BH/snapshots_200bkgdens/h229.cosmo50PLK.3072gst5HbwK1BH.004096' #
    tfile_2 = prefix + 'h229.cosmo50PLK.3072g/h229.cosmo50PLK.3072gst5HbwK1BH/h229.cosmo50PLK.3072gst5HbwK1BH.004096/h229.cosmo50PLK.3072gst5HbwK1BH.004096'    
    
    tfile_base_3 = 'h242.cosmo50PLK.3072gst5HbwK1BH'
    tfile_3 = prefix + 'h242.cosmo50PLK.3072g/h242.cosmo50PLK.3072gst5HbwK1BH/snapshots_200bkgdens/h242.cosmo50PLK.3072gst5HbwK1BH.004096' #
    tfile_3 = prefix + 'h242.cosmo50PLK.3072g/h242.cosmo50PLK.3072gst5HbwK1BH/h242.cosmo50PLK.3072gst5HbwK1BH.004096/h242.cosmo50PLK.3072gst5HbwK1BH.004096'

    tfile_base_4 = 'h329.cosmo50PLK.3072gst5HbwK1BH'
    tfile_4 = prefix + 'h329.cosmo50PLK.3072g/h329.cosmo50PLK.3072gst5HbwK1BH/snapshots_200bkgdens/h329.cosmo50PLK.3072gst5HbwK1BH.004096' #
    tfile_4 = prefix + 'h329.cosmo50PLK.3072g/h329.cosmo50PLK.3072gst5HbwK1BH/h329.cosmo50PLK.3072gst5HbwK1BH.004096/h329.cosmo50PLK.3072gst5HbwK1BH.004096'
    
    tfile_base_4hr = 'h329.cosmo50PLK.6144g5HbwK1BH'
    tfile_4hr = prefix + 'h329.cosmo50PLK.6144g/h329.cosmo50PLK.6144g5HbwK1BH/h329.cosmo50PLK.6144g5HbwK1BH.004096/ahf_200/h329.cosmo50PLK.6144g5HbwK1BH.004096' #    
    
    tfiles = [tfile_cm, tfile_e, tfile_r, tfile_s, tfile_1, tfile_1hr, tfile_2, tfile_3, tfile_4, tfile_4hr]
    tfile_base = [tfile_base_cm, tfile_base_e, tfile_base_r, tfile_base_s, tfile_base_1, tfile_base_1hr, tfile_base_2, tfile_base_3, tfile_base_4, tfile_base_4hr]

    objs_pd = None 
    for tfile, base in zip(tfiles, tfile_base):
        objs_dat = []
        #print(tfile)
        f=open(tfile + '.MAP.data', 'rb')
        while 1:
            try:
                objs_dat.append(pickle.load(f))
            except EOFError:
                break        
        f.close()
        if len(objs_dat) == 1:
            temp = pd.DataFrame(objs_dat[0])
        else:
            temp = pd.DataFrame(objs_dat)
        simname = base.split('.')[0]
        if (base.split('.')[2])[0] == '6':
            simname = simname+'_6144'
            temp['M_star'] = temp['mstar']
            temp['mass'] = temp['mvir']
            
        temp['sim'] = [simname]*len(temp)
        if not 'massiveDist' in temp:
            temp = distance_to_nearest_host(temp,[tfile])
            temp.to_pickle(tfile + '.MAP.data')

        if objs_pd is None: 
            objs_pd = temp
        else:
            objs_pd = objs_pd.append(temp, ignore_index = True)     
    
    fdmdata_mod = fdmdata.copy()
    objs_pd = match_halos(objs_pd, fdmdata_mod)
    
    ind = 0
    tau90 = np.empty(len(objs_pd))            
    for index, row in objs_pd.iterrows():
        if len(row['sfhbins']) != len(row['sfh']):
            xarr = row['sfhbins'][1:] - (row['sfhbins'][1] - row['sfhbins'][0])
        else:
            xarr = row['sfhbins'][:]
        yarr = np.cumsum(row['sfh'])/max(np.cumsum(row['sfh']))
        if (yarr[0] >= 0.9):
            tau90[ind] = xarr[0]
        else:
            interp = interp1d(yarr, xarr) #, kind='cubic')
            if np.isnan(interp(0.9)):
                tau90[ind] = 0
            else:
                tau90[ind] = float(interp(0.9))
        ind = ind + 1            

    objs_pd['tau90'] = tau90
    
    halo_label = {'halo_label': ""}
    objs_pd = objs_pd.join(pd.DataFrame(columns=halo_label))
    objs_pd.loc[~objs_pd['m200_haloid'].isnull(),'halo_label'] = objs_pd[~objs_pd['m200_haloid'].isnull()]['sim']+objs_pd[~objs_pd['m200_haloid'].isnull()]['m200_haloid'].astype(str)
    objs_pd = objs_pd.set_index('halo_label')

    fdmdata = fdmdata.join(pd.DataFrame(columns=halo_label))
    fdmdata['halo_label'] = fdmdata['simname']+fdmdata['halogrp_z0'].astype(str)
    fdmdata = fdmdata.set_index('halo_label') 

    objs_pd_comb = pd.concat([objs_pd,fdmdata], join="inner", axis=1)

    '''
    plt.clf()
    plt.figure(1)
    ind = 0
    for index, row in objs_pd.iterrows():
        if (type(row['sfhbins']) == float):
            tau90[ind] = 0
        else:
            if len(row['sfhbins']) != len(row['sfh']):
                xarr = row['sfhbins'][1:] - (row['sfhbins'][1] - row['sfhbins'][0])
            else:
                xarr = row['sfhbins'][:]
            yarr = np.cumsum(row['sfh'])/max(np.cumsum(row['sfh']))
            plt.clf()
            plt.plot(xarr,yarr)
            plt.plot([0,14],[0.9,0.9])
            if (yarr[0] >= 0.9):
                tau90[ind] = xarr[0]
            else:
                interp = interp1d(yarr, xarr) #, kind='cubic')
                print(ind,interp(0.9))
                tau90[ind] = float(interp(0.9))
        ind = ind + 1
    '''


    #scalarMap = cm.ScalarMappable(norm=cNorm, cmap=cmx)
    plt.clf()
    fig1 = plt.figure(1,figsize = (plt_width,plt_width*aspect_ratio))
    fig1.clear()
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax1 = fig1.add_subplot(gs[0])
    ax1sub = fig1.add_subplot(gs[1])

    cmx = plt.get_cmap("viridis") 
    cNorm  = colors.Normalize(vmin=0, vmax = 14)    
    q = ax1.scatter(objs_pd[objs_pd['SFR'] < 1e-11]['massiveDist'],objs_pd['M_star'][objs_pd['SFR'] < 1e-11],s = (objs_pd['Rvir'][objs_pd['SFR'] < 1e-11]*2*ms_scale).tolist(), c = tau90[objs_pd['SFR'] < 1e-11], cmap = cmx, norm = cNorm,edgecolor = 'k',marker = 'D', linewidths = edgewidth)
    sf = ax1.scatter(objs_pd[objs_pd['SFR'] >= 1e-11]['massiveDist'],objs_pd['M_star'][objs_pd['SFR'] >= 1e-11],s = (objs_pd['Rvir'][objs_pd['SFR'] >= 1e-11]*2*ms_scale).tolist(), c = tau90[objs_pd['SFR'] >= 1e-11], cmap = cmx, norm = cNorm,edgecolor = 'k', linewidths = edgewidth)
    #plt.scatter(objs_pd_e['h1dist'],objs_pd_e['M_star'])
    lgnd = ax1.legend([q,sf],['Quenched','Star forming'],scatterpoints = 1,facecolor = 'white',loc = 3,framealpha = 0,frameon = False)
    lgnd.legendHandles[0]._sizes = [markersize]
    lgnd.legendHandles[1]._sizes = [markersize]
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.axis([10, 1e4, 1e2, 1e10])
    ax1.set_ylabel(r'M$_*$/M$_\odot$')
    ax1.set_xlabel(r'Distance to massive galaxy (kpc)')
    cb = mpl.colorbar.ColorbarBase(ax1sub, cmap=cmx, norm=cNorm)
    cb.set_label(r"$\tau_{90}$ (Gyr)")
    fig1.set_size_inches(plt_width,plt_width*aspect_ratio)
    fig1.tight_layout()
    fig1.show()
    fig1.savefig(outbase + '_distance_smass_t90.png',dpi = dpi)    

'''
    fig1= plt.figure() 
    ax1 = fig1.add_subplot(1,1,1)
    ax1.plot(objs_pd_comb['distances'],objs_pd_comb['Mhalo_z0']/objs_pd_comb['Mpeak'],'x')
    ax1.plot(objs_pd_comb['Mhalo_z0'],objs_pd_comb['Mstar_z0'],'o')
    ax1.plot(objs_pd_comb['Mpeak'],objs_pd_comb['Mstar_z0'],'o')
    ax1.set_yscale('log')
    ax1.set_xscale('log')
'''
    
