#Charlotte Christensen

#8/13/19
#Plot the fraction of stars formed before 4 and 8 Gyr in satellites and field galaxies to compare with Digby+ 2019


#%run /home/christenc/Code/python/python_analysis/sftimes_env
import matplotlib as mpl
#mpl.use('tkagg') #Also can try mpl.use('Agg') #for using mpl over ssh    
import matplotlib.pylab as plt
import matplotlib.mlab as mlab
import matplotlib.colors as colors
import matplotlib.cm as cm
import pynbody
import math
import numpy as np
import socket
import matplotlib.gridspec as gridspec
import pandas as pd
import sys, os, glob, pickle
from scipy.interpolate import interp1d
sys.path.append(os.path.abspath("/home/christenc/Code/python/python_analysis/"))
from modules.user_tools import task

# Matches the halos in the objs_pd file with the halo information from Ferah (m200)
# see "match_halos.py" for development
def match_halos(objs_pd, fdmdata):
    smass_tol = 0.5 #fractional tolerance for the stellar mass
    vmass_tol = 0.9

    match_id = {'m200_haloid': 0}
    if not 'm200_haloid' in objs_pd.keys():
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
                index_best_match = (fdmdata[fdmdata['simname']==sim][possible_match]).index[arg_best_match]
                objs_pd.loc[(objs_pd['sim'] == sim) & (objs_pd['haloid'] == halo),'m200_haloid'] = fdmdata.loc[index_best_match]['halogrp_z0']
                #print(sim, halo, fdmdata.loc[index_best_match]['halogrp_z0'], float(objs_pd[(objs_pd['sim'] == sim) & (objs_pd['haloid'] == halo)]['M_star']), fdmdata.loc[index_best_match]['Mstar_z0'], float(objs_pd[(objs_pd['sim'] == sim) & (objs_pd['haloid'] == halo)]['mass']), fdmdata.loc[index_best_match]['Mhalo_z0'])
                fdmdata.loc[index_best_match,'Mstar_z0'] = 0 #set stellar mass to zero so it isn't used again

    return objs_pd

def sftimes_env(tfiles,outfile_base,tfile_base,*halo_nums):
    presentation = False
    if presentation:
        outfile_base = outfile_base + '_pres'
        plt.style.use(['default','/home/christenc/.config/matplotlib/presentation.mplstyle'])
        plt_width = 8 #inches
        aspect_ratio = 3.0/4.0
        legendsize = 16
        dpi = 200
        markersize = 100
        ms_scale = 1
        lw = mpl.rcParams['lines.linewidth'] - 1
        edgewidth = 1
    else:
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

    use_m200 = 1
    if use_m200:
        ext = '.m200.'
    else:
        ext = '.MAP.'

    # Data from Digby+ 2019
    digbyfield = pd.read_csv(dataprefix + 'DigbyField.csv')
    digbysat = pd.read_csv(dataprefix + 'DigbySatellite.csv')
        
    f = open(dataprefix+'mstar_vs_mhalo_extended_corrected5.txt', 'r')
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

    objs_pd = None 
    for tfile, base in zip(tfiles, tfile_base):
        objs_dat = []
        print(tfile)
        
        f=open(tfile + ext + 'data', 'rb')
        while 1:
            try:
                objs_dat.append(pickle.load(f))
            except EOFError:
                break        
        f.close()
        
        #objs_dat = pd.read_csv(tfile + '.MAP.data.csv')
        if len(objs_dat) == 1:
            temp = pd.DataFrame(objs_dat[0])
        else:
            temp = pd.DataFrame(objs_dat)
        simname = base.split('.')[0]
        if (base.split('.')[2])[0] == '6':
            simname = simname+'_6144'
        if not ('M_star' in temp.keys()):
            temp['M_star'] = temp['mstar']
        if not ('mass' in temp.keys()):            
            temp['mass'] = temp['mvir']
            
        temp['sim'] = [simname]*len(temp)

        if objs_pd is None: 
            objs_pd = temp
        else:
            objs_pd = objs_pd.append(temp, ignore_index = True)
        

    ind = 0
    f1Gyr = np.empty(len(objs_pd))
    f4Gyr = np.empty(len(objs_pd))
    f8Gyr = np.empty(len(objs_pd))
    
    for index, row in objs_pd.iterrows():
#        print(row['sim'],row['haloid'])
        sfh = row['sfh']
        sfhbins = row['sfhbins']
        
        if len(sfhbins) != len(sfh):
            xarr = sfhbins[1:] - (sfhbins[1] - sfhbins[0])
        else:
            xarr = sfhbins[:]
        yarr = np.cumsum(sfh)/max(np.cumsum(sfh))
        interp = interp1d(xarr, yarr)
        if min(xarr) > 1:
            f1Gyr[ind] = 0
        else:
            if max(xarr) < 1:
                f1Gyr[ind] = 1.0
            else:
                f1Gyr[ind] = float(interp(1))
                
        if min(xarr) > 4:
            f4Gyr[ind] = 0
        else:                
            if max(xarr) < 4:
                f4Gyr[ind] = 1.0
            else:
                f4Gyr[ind] = float(interp(4))
                
        if min(xarr) > 8:
            f8Gyr[ind] = 0
        else:
            if max(xarr) < 8:
                f8Gyr[ind] = 1.0
            else:                       
                f8Gyr[ind] = float(interp(8))
#        print(ind,f1Gyr[ind],f4Gyr[ind],f8Gyr[ind])
            
        ind = ind + 1

    objs_pd['f1Gyr'] = f1Gyr
    objs_pd['f4Gyr'] = f4Gyr
    objs_pd['f8Gyr'] = f8Gyr

    #Match halos between my and Ferah's data
    fdmdata_mod = fdmdata.copy()
    objs_pd = match_halos(objs_pd, fdmdata_mod)
    #remove entries (rows) from objs_pd that have no match in Ferah's data
    index_rm = objs_pd[(objs_pd['m200_haloid']).isnull()].index
    objs_pd = objs_pd.drop(index_rm)
    #objs_pd = distance_to_nearest_host(objs_pd, tfiles)

    halo_label = {'halo_label': ""}
    objs_pd = objs_pd.join(pd.DataFrame(columns=halo_label))
    objs_pd.loc[~objs_pd['m200_haloid'].isnull(),'halo_label'] = objs_pd[~objs_pd['m200_haloid'].isnull()]['sim']+objs_pd[~objs_pd['m200_haloid'].isnull()]['m200_haloid'].astype(str)
    objs_pd = objs_pd.set_index('halo_label')

    fdmdata = fdmdata.join(pd.DataFrame(columns=halo_label))
    fdmdata['halo_label'] = fdmdata['simname']+fdmdata['halogrp_z0'].astype(str)
    fdmdata = fdmdata.set_index('halo_label') 

    objs_pd_comb = pd.concat([objs_pd,fdmdata], join="inner", axis=1)    
    
    plt.close('all')    
    fig1 = plt.figure(figsize=(plt_width*2,plt_width*aspect_ratio*2))
    gs = fig1.add_gridspec(2,2, hspace=0, wspace=0) #, vspace=0)
    axs1 = gs.subplots(sharex = True, sharey = True) #, constrained_layout=True)
    axs1 = axs1.flatten()
    ax1a = axs1[0]
    ax1b = axs1[1]
    ax1c = axs1[2]
    ax1d = axs1[3]
    ax1a.scatter(objs_pd_comb['mstar'][objs_pd_comb['type']=='Central'],objs_pd_comb['f4Gyr'][objs_pd_comb['type']=='Central'], c = 'b', s = markersize, alpha = 0.6)
    ax1a.scatter(objs_pd_comb['mstar'][objs_pd_comb['type']=='Backsplash'],objs_pd_comb['f4Gyr'][objs_pd_comb['type']=='Backsplash'], c = 'b', s = markersize, alpha = 0.6)
    ax1a.errorbar(digbyfield['Mstar'][digbyfield['oMSTO']=='y'],digbyfield['f4'][digbyfield['oMSTO']=='y'], yerr = np.array([digbyfield['f4_lerr'][digbyfield['oMSTO']=='y'].tolist(),digbyfield['f4_herr'][digbyfield['oMSTO']=='y'].tolist()]), c = 'r', markersize = markersize*0.2, alpha = 0.8, marker = 'o', ls = 'None')
    ax1a.errorbar(digbyfield['Mstar'][digbyfield['oMSTO']=='n'],digbyfield['f4'][digbyfield['oMSTO']=='n'], yerr = np.array([digbyfield['f4_lerr'][digbyfield['oMSTO']=='n'].tolist(),digbyfield['f4_herr'][digbyfield['oMSTO']=='n'].tolist()]), c = 'r', markersize = markersize*0.1, alpha = 0.1, marker = 'o', ls = 'None')    
    ax1a.set_xscale('log')
    ax1a.set_ylabel(r'$f_{4Gyr}$')
    ax1a.text(1e8,0.9,'Field')
    #ax1b.scatter(objs_pd_comb['mstar'][objs_pd_comb['type']=='Backsplash'],objs_pd_comb['f4Gyr'][objs_pd_comb['type']=='Backsplash'], c = 'g', s = markersize, alpha = 0.6)    
    ax1b.scatter(objs_pd_comb['mstar'][objs_pd_comb['type']=='Satellite'],objs_pd_comb['f4Gyr'][objs_pd_comb['type']=='Satellite'], c = 'g', s = markersize, alpha = 0.6)
    ax1b.errorbar(digbysat['Mstar'][digbysat['oMSTO']=='y'],digbysat['f4'][digbysat['oMSTO']=='y'], yerr = np.array([digbysat['f4_lerr'][digbysat['oMSTO']=='y'].tolist(),digbysat['f4_herr'][digbysat['oMSTO']=='y'].tolist()]), c = 'r', markersize = markersize*0.2, alpha = 0.8, marker = 'o', ls = 'None')
    ax1b.errorbar(digbysat['Mstar'][digbysat['oMSTO']=='n'],digbysat['f4'][digbysat['oMSTO']=='n'], yerr = np.array([digbysat['f4_lerr'][digbysat['oMSTO']=='n'].tolist(),digbysat['f4_herr'][digbysat['oMSTO']=='n'].tolist()]), c = 'r', markersize = markersize*0.1, alpha = 0.1, marker = 'o', ls = 'None')    
    ax1a.set_xscale('log')    
    ax1b.set_xscale('log')
    ax1b.text(1e8,0.9,'Satellites')
    ax1c.scatter(objs_pd_comb['mstar'][objs_pd_comb['type']=='Central'],objs_pd_comb['f8Gyr'][objs_pd_comb['type']=='Central'], c = 'b', s = markersize, alpha = 0.6)
    ax1c.scatter(objs_pd_comb['mstar'][objs_pd_comb['type']=='Backsplash'],objs_pd_comb['f8Gyr'][objs_pd_comb['type']=='Backsplash'], c = 'b', s = markersize, alpha = 0.6)
    ax1c.errorbar(digbyfield['Mstar'][digbyfield['oMSTO']=='y'],digbyfield['f8'][digbyfield['oMSTO']=='y'], yerr = np.array([digbyfield['f8_lerr'][digbyfield['oMSTO']=='y'].tolist(),digbyfield['f8_herr'][digbyfield['oMSTO']=='y'].tolist()]), c = 'r', markersize = markersize*0.2, alpha = 0.8, marker = 'o', ls = 'None')
    ax1c.errorbar(digbyfield['Mstar'][digbyfield['oMSTO']=='n'],digbyfield['f8'][digbyfield['oMSTO']=='n'], yerr = np.array([digbyfield['f8_lerr'][digbyfield['oMSTO']=='n'].tolist(),digbyfield['f8_herr'][digbyfield['oMSTO']=='n'].tolist()]), c = 'r', markersize = markersize*0.1, alpha = 0.1, marker = 'o', ls = 'None')     
    ax1c.set_xscale('log')
    ax1c.set_ylabel(r'$f_{8Gyr}$')
    ax1c.set_xlabel(r'$M_*/M_\odot$')
    #ax1d.scatter(objs_pd_comb['mstar'][objs_pd_comb['type']=='Backsplash'],objs_pd_comb['f8Gyr'][objs_pd_comb['type']=='Backsplash'], c = 'g', s = markersize, alpha = 0.6)    
    ax1d.scatter(objs_pd_comb['mstar'][objs_pd_comb['type']=='Satellite'],objs_pd_comb['f8Gyr'][objs_pd_comb['type']=='Satellite'], c = 'g', s = markersize, alpha = 0.6)
    ax1d.errorbar(digbysat['Mstar'][digbysat['oMSTO']=='y'],digbysat['f8'][digbysat['oMSTO']=='y'], yerr = np.array([digbysat['f8_lerr'][digbysat['oMSTO']=='y'].tolist(),digbysat['f8_herr'][digbysat['oMSTO']=='y'].tolist()]), c = 'r', markersize = markersize*0.2, alpha = 0.8, marker = 'o', ls = 'None')
    ax1d.errorbar(digbysat['Mstar'][digbysat['oMSTO']=='n'],digbysat['f8'][digbysat['oMSTO']=='n'], yerr = np.array([digbysat['f8_lerr'][digbysat['oMSTO']=='n'].tolist(),digbysat['f8_herr'][digbysat['oMSTO']=='n'].tolist()]), c = 'r', markersize = markersize*0.1, alpha = 0.1, marker = 'o', ls = 'None')     
    ax1d.set_xlabel(r'$M_*/M_\odot$')    
    ax1d.set_xscale('log')
    ax1d.axis([5e3, 5e9, -0.05, 1.05])
#    fig1.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.8, wspace=0, hspace=0)
    fig1.tight_layout()
    fig1.show()
    fig1.savefig(outfile_base + '_f4_and_f8.png',dpi = dpi)

    plt.close('all')
    fig1 = plt.figure(1,figsize=(plt_width*2,plt_width*aspect_ratio))
    gs = fig1.add_gridspec(1,2,wspace=0)
    axs1 = gs.subplots(sharey = True) #, constrained_layout=True) 
    axs1 = axs1.flatten()
    ax1a = axs1[0]
    ax1b = axs1[1]    
    ax1a.scatter(objs_pd_comb['mstar'][objs_pd_comb['type']=='Central'],objs_pd_comb['f8Gyr'][objs_pd_comb['type']=='Central'] - objs_pd_comb['f4Gyr'][objs_pd_comb['type']=='Central'], c = 'b', s = markersize, alpha = 0.6)
    ax1a.scatter(objs_pd_comb['mstar'][objs_pd_comb['type']=='Backsplash'],objs_pd_comb['f8Gyr'][objs_pd_comb['type']=='Backsplash'] - objs_pd_comb['f4Gyr'][objs_pd_comb['type']=='Backsplash'], c = 'b', s = markersize, alpha = 0.6)
    ax1a.errorbar(digbyfield['Mstar'][digbyfield['oMSTO']=='y'],digbyfield['f8'][digbyfield['oMSTO']=='y']-digbyfield['f4'][digbyfield['oMSTO']=='y'], yerr = np.array([digbyfield['f4_lerr'][digbyfield['oMSTO']=='y'].tolist(),digbyfield['f4_herr'][digbyfield['oMSTO']=='y'].tolist()]), c = 'r', markersize = markersize*0.2, alpha = 0.8, marker = 'o', ls = 'None')
    ax1a.errorbar(digbyfield['Mstar'][digbyfield['oMSTO']=='n'],digbyfield['f8'][digbyfield['oMSTO']=='n']-digbyfield['f4'][digbyfield['oMSTO']=='n'], yerr = np.array([digbyfield['f4_lerr'][digbyfield['oMSTO']=='n'].tolist(),digbyfield['f4_herr'][digbyfield['oMSTO']=='n'].tolist()]), c = 'r', markersize = markersize*0.1, alpha = 0.1, marker = 'o', ls = 'None')     
    ax1a.set_xscale('log')
    ax1a.set_ylabel(r'($f_{8Gyr} - f_{4Gyr}$)')
    ax1a.set_xlabel(r'$M_*/M_\odot$')    
    ax1a.text(1e8,0.9,'Field')
    ax1a.axis([5e3, 5e9, -0.05, 1.05])    
    #ax1b.scatter(objs_pd_comb['mstar'][objs_pd_comb['type']=='Backsplash'],objs_pd_comb['f8Gyr'][objs_pd_comb['type']=='Backsplash'] - objs_pd_comb['f4Gyr'][objs_pd_comb['type']=='Backsplash'], c = 'g', s = markersize, alpha = 0.6)    
    ax1b.scatter(objs_pd_comb['mstar'][objs_pd_comb['type']=='Satellite'],objs_pd_comb['f8Gyr'][objs_pd_comb['type']=='Satellite'] - objs_pd_comb['f4Gyr'][objs_pd_comb['type']=='Satellite'], c = 'g', s = markersize, alpha = 0.6)
    ax1b.errorbar(digbysat['Mstar'][digbysat['oMSTO']=='y'],digbysat['f8'][digbysat['oMSTO']=='y']-digbysat['f4'][digbysat['oMSTO']=='y'], yerr = np.array([digbysat['f4_lerr'][digbysat['oMSTO']=='y'].tolist(),digbysat['f4_herr'][digbysat['oMSTO']=='y'].tolist()]), c = 'r', markersize = markersize*0.2, alpha = 0.8, marker = 'o', ls = 'None')
    ax1b.errorbar(digbysat['Mstar'][digbysat['oMSTO']=='n'],digbysat['f8'][digbysat['oMSTO']=='n']-digbysat['f4'][digbysat['oMSTO']=='n'], yerr = np.array([digbysat['f4_lerr'][digbysat['oMSTO']=='n'].tolist(),digbysat['f4_herr'][digbysat['oMSTO']=='n'].tolist()]), c = 'r', markersize = markersize*0.1, alpha = 0.1, marker = 'o', ls = 'None')       
    ax1b.set_xscale('log')
    ax1b.text(1e8,0.9,'Satellites')
    ax1b.set_xlabel(r'$M_*/M_\odot$')
    ax1b.axis([5e3, 5e9, -0.05, 1.05])
    fig1.tight_layout()
    fig1.show()
    fig1.savefig(outfile_base + '_f4_f8.png',dpi = dpi)

    plt.close('all')
    fig1 = plt.figure(1,figsize=(plt_width*2,plt_width*aspect_ratio))
    gs = fig1.add_gridspec(1,2,wspace=0)
    axs1 = gs.subplots(sharey = True) #, constrained_layout=True) 
    axs1 = axs1.flatten()
    ax1a = axs1[0]
    ax1b = axs1[1]    
    ax1a.scatter(objs_pd_comb['mstar'][objs_pd_comb['type']=='Central'],objs_pd_comb['f1Gyr'][objs_pd_comb['type']=='Central'], c = 'b', s = markersize, alpha = 0.6)
    ax1a.scatter(objs_pd_comb['mstar'][objs_pd_comb['type']=='Backsplash'],objs_pd_comb['f1Gyr'][objs_pd_comb['type']=='Backsplash'], c = 'b', s = markersize, alpha = 0.6)
    ax1a.scatter(digbyfield['Mstar'][digbyfield['oMSTO']=='y'],np.array((digbyfield['f1'][digbyfield['oMSTO']=='y']).tolist()).astype('float'), c = 'r', s = markersize, alpha = 0.8, marker = 'o')
    #ax1a.errorbar(digbyfield['Mstar'][digbyfield['oMSTO']=='y'],np.array((digbyfield['f4'][digbyfield['oMSTO']=='y']).tolist()).astype('float'), yerr = [digbyfield['f1_lerr'][digbyfield['oMSTO']=='y'].tolist(),digbyfield['f1_herr'][digbyfield['oMSTO']=='y'].tolist()], c = 'r', markersize = markersize*0.2, alpha = 0.8, marker = 'o', ls = 'None')
    ax1a.set_xscale('log')
    ax1a.set_ylabel(r'$f_{1Gyr}$')
    ax1a.set_xlabel(r'$M_*/M_\odot$')
    ax1a.axis([5e3, 5e9, -0.05, 1.05])    
    ax1a.text(1e8,0.9,'Field')
    #ax1b.scatter(objs_pd_comb['mstar'][objs_pd_comb['type']=='Backsplash'],objs_pd_comb['f1Gyr'][objs_pd_comb['type']=='Backsplash'], c = 'g', s = markersize, alpha = 0.6)    
    ax1b.scatter(objs_pd_comb['mstar'][objs_pd_comb['type']=='Satellite'],objs_pd_comb['f1Gyr'][objs_pd_comb['type']=='Satellite'], c = 'g', s = markersize, alpha = 0.6)
    ax1b.scatter(digbysat['Mstar'][digbysat['oMSTO']=='y'],np.array((digbysat['f1'][digbysat['oMSTO']=='y']).tolist()).astype('float'), c = 'r', s = markersize, alpha = 0.8, marker = 'o', ls = 'None')    
    #ax1b.errorbar(digbysat['Mstar'][digbysat['oMSTO']=='y'],np.array((digbysat['f4'][digbysat['oMSTO']=='y']).tolist()).astype('float'), yerr = np.array([digbysat['f1_lerr'][digbysat['oMSTO']=='y'].tolist(),digbysat['f1_herr'][digbysat['oMSTO']=='y'].tolist()]), c = 'r', markersize = markersize*0.2, alpha = 0.8, marker = 'o', ls = 'None')
    ax1b.set_xscale('log')
    ax1b.text(1e8,0.9,'Satellites')
    ax1b.set_xlabel(r'$M_*/M_\odot$')
    ax1b.axis([5e3, 5e9, -0.05, 1.05])
    fig1.tight_layout()
    fig1.show()
    fig1.savefig(outfile_base + '_f1.png',dpi = dpi)    
    
if __name__ == '__main__':
    if (socket.gethostname() == "quirm.math.grinnell.edu"):
        prefix = '/home/christenc/Data/Sims/'
        prefix_outfile = '/home/christenc/Figures/marvel/'
    else:
        prefix = '/home/christensen/Storage2/UW/MolecH/Cosmo/'
        prefix_outfile = '/home/christensen/Plots/marvel/'

    tfile_base_cm = 'cptmarvel.cosmo25cmb.4096g5HbwK1BH'
    tfile_cm = prefix + 'cptmarvel.cosmo25cmb/cptmarvel.cosmo25cmb.4096g5HbwK1BH/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096/ahf_200/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096' #'cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096'

    tfile_r = prefix + 'rogue.cosmo25cmb/rogue.cosmo25cmb.4096g5HbwK1BH/rogue.cosmo25cmb.4096g5HbwK1BH.004096/ahf_200/rogue.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_base_r = 'rogue.cosmo25cmb.4096g5HbwK1BH'

    tfile_e = prefix + 'elektra.cosmo25cmb/elektra.cosmo25cmb.4096g5HbwK1BH/elektra.cosmo25cmb.4096g5HbwK1BH.004096/ahf_200/elektra.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_base_e = 'elektra.cosmo25cmb.4096g5HbwK1BH'

    tfile_s = prefix + 'storm.cosmo25cmb/storm.cosmo25cmb.4096g5HbwK1BH/storm.cosmo25cmb.4096g5HbwK1BH.004096/ahf_200/storm.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_base_s = 'storm.cosmo25cmb.4096g5HbwK1BH'
    
    tfile_base_1 = 'h148.cosmo50PLK.3072g3HbwK1BH'
    tfile_1 = prefix + 'h148.cosmo50PLK.3072g/h148.cosmo50PLK.3072g3HbwK1BH/snapshots_200bkgdens/h148.cosmo50PLK.3072g3HbwK1BH.004096' #
    tfile_1 = prefix + 'h148.cosmo50PLK.3072g/h148.cosmo50PLK.3072g3HbwK1BH/h148.cosmo50PLK.3072g3HbwK1BH.004096/ahf_200/h148.cosmo50PLK.3072g3HbwK1BH.004096'
    
    tfile_base_1hr = 'h148.cosmo50PLK.6144g3HbwK1BH/'
    tfile_1hr = prefix + 'h148.cosmo50PLK.6144g/h148.cosmo50PLK.6144g3HbwK1BH/h148.cosmo50PLK.6144g3HbwK1BH.004096/ahf_200/h148.cosmo50PLK.6144g3HbwK1BH.004096'
    
    tfile_base_2 = 'h229.cosmo50PLK.3072gst5HbwK1BH'
    tfile_2 = prefix + 'h229.cosmo50PLK.3072g/h229.cosmo50PLK.3072gst5HbwK1BH/snapshots_200bkgdens/h229.cosmo50PLK.3072gst5HbwK1BH.004096' #
    tfile_2 = prefix + 'h229.cosmo50PLK.3072g/h229.cosmo50PLK.3072gst5HbwK1BH/h229.cosmo50PLK.3072gst5HbwK1BH.004096/ahf_200/h229.cosmo50PLK.3072gst5HbwK1BH.004096'    
    
    tfile_base_3 = 'h242.cosmo50PLK.3072gst5HbwK1BH'
    tfile_3 = prefix + 'h242.cosmo50PLK.3072g/h242.cosmo50PLK.3072gst5HbwK1BH/snapshots_200bkgdens/h242.cosmo50PLK.3072gst5HbwK1BH.004096' #
    tfile_3 = prefix + 'h242.cosmo50PLK.3072g/h242.cosmo50PLK.3072gst5HbwK1BH/h242.cosmo50PLK.3072gst5HbwK1BH.004096/ahf_200/h242.cosmo50PLK.3072gst5HbwK1BH.004096'

    tfile_base_4 = 'h329.cosmo50PLK.3072gst5HbwK1BH'
    tfile_4 = prefix + 'h329.cosmo50PLK.3072g/h329.cosmo50PLK.3072gst5HbwK1BH/snapshots_200bkgdens/h329.cosmo50PLK.3072gst5HbwK1BH.004096' #
    tfile_4 = prefix + 'h329.cosmo50PLK.3072g/h329.cosmo50PLK.3072gst5HbwK1BH/h329.cosmo50PLK.3072gst5HbwK1BH.004096/ahf_200/h329.cosmo50PLK.3072gst5HbwK1BH.004096'
    
    tfile_base_4hr = 'h329.cosmo50PLK.6144g5HbwK1BH'
    tfile_4hr = prefix + 'h329.cosmo50PLK.6144g/h329.cosmo50PLK.6144g5HbwK1BH/h329.cosmo50PLK.6144g5HbwK1BH.004096/ahf_200/h329.cosmo50PLK.6144g5HbwK1BH.004096' #    
    
    outfile_base = prefix_outfile + 'marvelJL'
    #SMHMv__distance([tfile_1,tfile_2,tfile_3,tfile_4],outfile_base,[tfile_base_1,tfile_base_2,tfile_base_3,tfile_base_4])
    tfiles = [tfile_cm, tfile_e, tfile_r, tfile_s, tfile_1, tfile_1hr, tfile_2, tfile_3, tfile_4, tfile_4hr]
    tfile_base = [tfile_base_cm, tfile_base_e, tfile_base_r, tfile_base_s, tfile_base_1, tfile_base_1hr, tfile_base_2, tfile_base_3, tfile_base_4, tfile_base_4hr]
    sftimes_env([tfile_cm, tfile_e, tfile_r, tfile_s, tfile_1, tfile_1hr, tfile_2, tfile_3, tfile_4, tfile_4hr],outfile_base,[tfile_base_cm, tfile_base_e, tfile_base_r, tfile_base_s, tfile_base_1, tfile_base_1hr, tfile_base_2, tfile_base_3, tfile_base_4, tfile_base_4hr])

