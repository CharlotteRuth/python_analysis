#Charlotte Christensen

#8/13/19
#Plot the SMHM relation for the marvel and Justice League runs, coloring points according to distance/tau_90

#SMHM for environment paper

#%run /home/christenc/Code/python/python_analysis/SMHM_v_distance
import matplotlib as mpl
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
def distance_to_nearest_host(data,tfiles):
    massiveDist = []
    min_massiveHalo = 10**11.5
    sprev = ''
    tfile_it = 0
    for i in range(len(data)):
        s = data['sim'].tolist()[i]
        print(s,data['haloid'].tolist()[i])       
        if s=='h148' or s=='h229' or s=='h242' or s=='h329': # if sat simulation, find distance to halo 1
            h1dist = data['h1dist'].tolist()[i]*0.6776942783267969
            massiveDist.append(h1dist)
       
            h1rvir = data['Rvir'][(data.sim==s) & (data.haloid==1)].tolist()[0]*0.6776942783267969
            #hostrvirs.append(h1rvir)
           
        else: # if field simulation, find distance to nearest massive DM halo (currently > 0.5e12.5 Msol)
            if s != sprev:
                sim = pynbody.load(tfiles[tfile_it])
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
                #for index, halo in data.iterrows():

            properties = h_dummy[int(data['haloid'].tolist()[i])].properties
            massiveDist.append(min(((properties['Xc'] - loc[:,0])**2 + (properties['Yc'] - loc[:,1])**2 + (properties['Zc'] - loc[:,2])**2)**(0.5)))
            minind = np.where((((properties['Xc'] - loc[:,0])**2 + (properties['Yc'] - loc[:,1])**2 + (properties['Zc'] - loc[:,2])**2)**(0.5)) == massiveDist[-1])
            #hostrvirs.append(rvir[minind]*0.6776942783267969)
            if s != sprev:
                tfile_it = tfile_it + 1
            sprev = s

    data['massiveDist'] = massiveDist
    return data

def SMHM_v_distance(tfiles,outfile_base,tfile_base,*halo_nums):
    presentation = True
    if presentation:
        outfile_base = outfile_base + '_pres'
        plt.style.use(['default','/home/christenc/.config/matplotlib/presentation.mplstyle'])
        plt_width = 16 #8 #inches
        aspect_ratio = 3.0/4.0
        legendsize = 16
        dpi = 100
        lw = mpl.rcParams['lines.linewidth'] - 1        
    else:
        plt.style.use(['default','/home/christenc/.config/matplotlib/article.mplstyle'])
        plt_width = 3.5 #inches
        aspect_ratio = 3.0/4.0
        legendsize = 5
        dpi = 300
        lw = mpl.rcParams['lines.linewidth']    

    if (socket.gethostname() == "ozma.grinnell.edu"):
        dataprefix = '/home/christensen/Code/Datafiles/' 
    else:
        dataprefix = '/home/christenc/Code/Datafiles/'
    f = open(dataprefix+'Read2017.txt','r')
    readdata = []
    for line in f:
        line = line.strip()
        columns = line.split()
        if len(columns) == 22 and columns[0] != 'Galaxy':
            source = {}
            source['Galaxy'] = columns[0]
            source['vmax'] = float(columns[1])
            source['vmax_err'] = float(columns[2])
            source['i'] = float(columns[3])
            source['i_err_n'] = float(columns[4])
            source['i_err_p'] = float(columns[5])
            source['D'] = float(columns[6])
            source['D_err'] = float(columns[7])
            source['M_*'] = float(columns[8])
            source['M_*_err'] = float(columns[9])
            source['Mgas'] = float(columns[10])
            source['R*'] = float(columns[11])
            source['R_gas'] = float(columns[12])
            source['Rmin'] = float(columns[13])
            source['Rmax'] = float(columns[14])
            source['M200'] = float(columns[15])
            source['M200_err_n'] = float(columns[16])
            source['M200_err_p'] = float(columns[17])
            source['c'] = float(columns[18])
            source['c_err_n'] = float(columns[19])
            source['c_err_p'] = float(columns[20])
            source['Xi'] = float(columns[21])
            readdata.append(source)
    f.close()
    readdata = pd.DataFrame(readdata)

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
    
    objs_pd = None 
    for tfile, base in zip(tfiles, tfile_base):
        objs_dat = []
        print(tfile)
        f=open(tfile + '.MAP.data', 'rb')
        while 1:
            try:
                objs_dat.append(pickle.load(f))
            except EOFError:
                break        
        f.close()
        temp = pd.DataFrame(objs_dat)
        temp['sim'] = [base.split('.')[0]]*len(objs_dat)

        if objs_pd is None: 
            objs_pd = temp
        else:
            objs_pd = objs_pd.append(temp, ignore_index = True)       

    fdmdata_mod = fdmdata.copy()
    objs_pd = match_halos(objs_pd, fdmdata_mod)    
            
    objs_pd = distance_to_nearest_host(objs_pd, tfiles)
    
    """
        h_dummy = s.halos(dummy = True)
        loc = []        
        for AHFhalo in h_dummy:
            properties = AHFhalo.properties            
            if (properties['mass'] > min_massiveHalo):
                loc.append(np.array([properties['Xc']/properties['h'], properties['Yc']/properties['h'], properties['Zc']/properties['h']]))
        loc = np.array(loc)
        massiveDist = [] #Distance to nearest massive (M_vir > min_massiveHalo) galaxy
        for halo in objs_dat:
            massiveDist.append(min(((halo['Xc'] - loc[:,0])**2 + (halo['Yc'] - loc[:,1])**2 + (halo['Zc'] - loc[:,2])**2)**(0.5)))
        temp = pd.DataFrame(objs_dat)
        temp['massiveDist'] = massiveDist
        if objs_pd is None: 
            objs_pd = temp
        else:
            objs_pd = objs_pd.append(temp, ignore_index = True)
    """

    ind = 0
    tau90 = np.empty(len(objs_pd))            
    for index, row in objs_pd.iterrows():
        xarr = row['sfhbins'][1:] - (row['sfhbins'][1] - row['sfhbins'][0])
        yarr = np.cumsum(row['sfh'])/max(np.cumsum(row['sfh']))
        if (yarr[0] >= 0.9):
            tau90[ind] = xarr[0]
        else:
            interp = interp1d(yarr, xarr) #, kind='cubic')
            if np.isnan(interp(0.9)):
                print(row)
                print(row['sfhbins'],row['sfh'])
                tau90[ind] = 0
            else:
                print(ind,interp(0.9))
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
        
    #SMHM colored by tau_90
    #fig1.clear()
    fig1 = plt.figure(1,figsize=(plt_width,plt_width*aspect_ratio))
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax1 = fig1.add_subplot(gs[0])
    ax1sub = fig1.add_subplot(gs[1])

    cmx = plt.get_cmap("viridis") 
    cNorm  = colors.Normalize(vmin=0, vmax = 14)    
    sat_v_tau90 = ax1.scatter(objs_pd_comb['mass'][objs_pd_comb['type']=='Satellite'],objs_pd_comb['M_star'][objs_pd_comb['type']=='Satellite'],c = objs_pd_comb['tau90'][objs_pd_comb['type']=='Satellite'], cmap = cmx, norm = cNorm,edgecolor = 'k')
    cen_v_tau90 = ax1.scatter(objs_pd_comb['mass'][objs_pd_comb['type']=='Central'],objs_pd_comb['M_star'][objs_pd_comb['type']=='Central'],c = objs_pd_comb['tau90'][objs_pd_comb['type']=='Central'], cmap = cmx, norm = cNorm,edgecolor = 'k',marker = 'D')
    ax1.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_ylabel(r'M$_*$/M$_\odot$')
    ax1.set_xlabel(r'M$_{\mathrm{vir}}$/M$_\odot$')
    ax1.axis([2e6, 2e11, 2e2, 5e9])
    ax1.legend([cen_v_tau90,sat_v_tau90],['Central','Satellite'],scatterpoints = 1,facecolor = 'white',loc = 0,framealpha = 0,frameon = False) 
    cb = mpl.colorbar.ColorbarBase(ax1sub, cmap=cmx, norm=cNorm)
    cb.set_label(r"$\tau_{90}$ (Gyr)")
    plt.savefig(outfile_base + '_SMHM_t90.png')

    #SMHM colored by tau_90
    fig1.clear()
    fig1a = plt.figure(10,figsize=(plt_width,plt_width*aspect_ratio))
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax1 = fig1a.add_subplot(gs[0])
    ax1sub = fig1a.add_subplot(gs[1])

    cmx = plt.get_cmap("viridis") 
    cNorm  = colors.Normalize(vmin=0, vmax = 14)    
    sat_v_tau90 = ax1.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],objs_pd_comb['Mstar_Mpeak'][objs_pd_comb['type']=='Satellite'],c = objs_pd_comb['tau90'][objs_pd_comb['type']=='Satellite'], cmap = cmx, norm = cNorm,edgecolor = 'k')
    cen_v_tau90 = ax1.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],objs_pd_comb['Mstar_Mpeak'][objs_pd_comb['type']=='Central'],c = objs_pd_comb['tau90'][objs_pd_comb['type']=='Central'], cmap = cmx, norm = cNorm,edgecolor = 'k',marker = 'D')
    ax1.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_ylabel(r'M$_*$/M$_\odot$')
    ax1.set_xlabel(r'M$_{\mathrm{vir, peak}}$/M$_\odot$')
    # ax1.axis([2e6, 2e11, 2e2, 5e9])
    ax1.axis([1e8, 2e11, 2e2, 5e9])
    ax1.legend([cen_v_tau90,sat_v_tau90],['Central','Satellite'],scatterpoints = 1,facecolor = 'white',loc = 0,framealpha = 0,frameon = False) 
    cb = mpl.colorbar.ColorbarBase(ax1sub, cmap=cmx, norm=cNorm)
    cb.set_label(r"$\tau_{90}$ (Gyr)")
    plt.savefig(outfile_base + '_SMHM_t90_Mpeak.png')
    
    #SMHM colored by distance to massive galaxy
    #fig2.clear()
    fig2 = plt.figure(2,figsize=(plt_width,plt_width*aspect_ratio))
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax2 = fig2.add_subplot(gs[0])
    ax2sub = fig2.add_subplot(gs[1])

    cNorm  = colors.Normalize(vmin=1.5, vmax = 4)    
    ax2.scatter(objs_pd['mass'],objs_pd['M_star'],c = np.log10(objs_pd['massiveDist']), cmap = cmx, norm = cNorm,edgecolor = 'k')
    ax2.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_ylabel(r'M$_*$/M$_\odot$')
    ax2.set_xlabel(r'M$_{\mathrm{vir}}$/M$_\odot$')
    ax2.axis([2e6, 2e11, 2e2, 5e9])
    cb = mpl.colorbar.ColorbarBase(ax2sub, cmap=cmx, norm=cNorm)
    cb.set_label(r"Log(Distance to Massive Galaxy/1 kpc)")
    plt.savefig(outfile_base + '_SMHM_rMassGal.png')

    #SMHM colored by distance to massive galaxy
    fig2.clear()
    fig2 = plt.figure(2,figsize=(plt_width,plt_width*aspect_ratio))
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax2 = fig2.add_subplot(gs[0])
    ax2sub = fig2.add_subplot(gs[1])

    cNorm  = colors.Normalize(vmin=1.5, vmax = 4)    
    ax2.scatter(objs_pd_comb['Mpeak'],objs_pd_comb['Mstar_Mpeak'],c = np.log10(objs_pd_comb['massiveDist']), cmap = cmx, norm = cNorm,edgecolor = 'k')
    ax2.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_ylabel(r'M$_*$/M$_\odot$')
    ax2.set_xlabel(r'M$_{\mathrm{vir, peak}}$/M$_\odot$')
    #ax2.axis([2e6, 2e11, 2e2, 5e9])
    ax2.axis([1e8, 2e11, 2e2, 5e9])
    cb = mpl.colorbar.ColorbarBase(ax2sub, cmap=cmx, norm=cNorm)
    cb.set_label(r"Log(Distance to Massive Galaxy/1 kpc)")
    plt.savefig(outfile_base + '_SMHM_rMassGal_Mpeak.png')
    
    #SMHM colored by HI mass
    #fig3.clear()
    fig3 = plt.figure(3,figsize=(plt_width,plt_width*aspect_ratio))
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax3 = fig3.add_subplot(gs[0])
    ax3sub = fig3.add_subplot(gs[1])
    
    cNorm  = colors.Normalize(vmin=2, vmax = 10)
    ax3.scatter(objs_pd['mass'],objs_pd['M_star'] ,edgecolor = 'k',facecolor = 'none')    
    ax3.scatter(objs_pd['mass'],objs_pd['M_star'],c = np.log10(np.array(objs_pd['mHI'].tolist())), cmap = cmx, norm = cNorm,edgecolor = 'k')
    ax3.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax3.set_xscale('log')
    ax3.set_yscale('log')
    ax3.set_ylabel(r'M$_*$/M$_\odot$')
    ax3.set_xlabel(r'M$_{\mathrm{vir}}$/M$_\odot$')
    ax3.axis([2e6, 2e11, 2e2, 5e9])
    cb = mpl.colorbar.ColorbarBase(ax3sub, cmap=cmx, norm=cNorm)
    cb.set_label(r"M$_{\mathrm{HI}}$ [M$_\odot$]")
    plt.savefig(outfile_base + '_SMHM_mHI.png')

    fig3.clear()
    fig3 = plt.figure(3,figsize=(plt_width,plt_width*aspect_ratio))
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax3 = fig3.add_subplot(gs[0])
    ax3sub = fig3.add_subplot(gs[1])
    
    cNorm  = colors.Normalize(vmin=2, vmax = 10)
    ax3.scatter(objs_pd_comb['Mpeak'],objs_pd_comb['Mstar_Mpeak'] ,edgecolor = 'k',facecolor = 'none')    
    ax3.scatter(objs_pd_comb['Mpeak'],objs_pd_comb['Mstar_Mpeak'],c = np.log10(np.array(objs_pd_comb['mHI'].tolist())), cmap = cmx, norm = cNorm,edgecolor = 'k')
    ax3.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax3.set_xscale('log')
    ax3.set_yscale('log')
    ax3.set_ylabel(r'M$_*$/M$_\odot$')
    ax3.set_xlabel(r'M$_{\mathrm{vir, peak}}$/M$_\odot$')
    ax3.axis([2e6, 2e11, 2e2, 5e9])
    cb = mpl.colorbar.ColorbarBase(ax3sub, cmap=cmx, norm=cNorm)
    cb.set_label(r"M$_{\mathrm{HI}}$ [M$_\odot$]")
    plt.savefig(outfile_base + '_SMHM_mHI_Mpeak.png')
   

    #SMHM colored by HI fraction
    #fig4.clear()
    fig4 = plt.figure(4,figsize=(plt_width,plt_width*aspect_ratio))
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax4 = fig4.add_subplot(gs[0])
    ax4sub = fig4.add_subplot(gs[1])
    
    cNorm  = colors.Normalize(vmin=0, vmax = 1)
    ax4.scatter(objs_pd['mass'],objs_pd['M_star'],c = objs_pd['mHI']/(objs_pd['M_star'] + objs_pd['mHI']), cmap = cmx, norm = cNorm,edgecolor = 'k')
    ax4.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax4.set_xscale('log')
    ax4.set_yscale('log')
    ax4.set_ylabel(r'M$_*$/M$_\odot$')
    ax4.set_xlabel(r'M$_{\mathrm{vir}}$/M$_\odot$')
    ax4.axis([2e6, 2e11, 2e2, 5e9])
    cb = mpl.colorbar.ColorbarBase(ax4sub, cmap=cmx, norm=cNorm)
    cb.set_label(r"M$_{\mathrm{HI}}$/(M$_*$ + M$_{\mathrm{HI}}$)")
    plt.savefig(outfile_base + '_SMHM_fHI.png')

    fig4.clear()
    fig4 = plt.figure(4,figsize=(plt_width,plt_width*aspect_ratio))
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax4 = fig4.add_subplot(gs[0])
    ax4sub = fig4.add_subplot(gs[1])
    
    cNorm  = colors.Normalize(vmin=0, vmax = 1)
    ax4.scatter(objs_pd_comb['Mpeak'],objs_pd_comb['Mstar_Mpeak'],c = objs_pd_comb['mHI']/(objs_pd_comb['M_star'] + objs_pd_comb['mHI']), cmap = cmx, norm = cNorm,edgecolor = 'k')
    ax4.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax4.set_xscale('log')
    ax4.set_yscale('log')
    ax4.set_ylabel(r'M$_*$/M$_\odot$')
    ax4.set_xlabel(r'M$_{\mathrm{vir, peak}}$/M$_\odot$')
    ax4.axis([1e8, 2e11, 2e2, 5e9])
    cb = mpl.colorbar.ColorbarBase(ax4sub, cmap=cmx, norm=cNorm)
    cb.set_label(r"M$_{\mathrm{HI}}$/(M$_*$ + M$_{\mathrm{HI}}$)")
    plt.savefig(outfile_base + '_SMHM_fHI_Mpeak.png')   
    
    #Baryonic mass in disk vs. halo mass, colored by distance to massive galaxy
    #fig5.clear()
    fig5 = plt.figure(5,figsize=(plt_width,plt_width*aspect_ratio))
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax5 = fig5.add_subplot(gs[0])
    ax5sub = fig5.add_subplot(gs[1])

    cNorm  = colors.Normalize(vmin=1.5, vmax = 4)    
    ax5.scatter(objs_pd['mass'],(objs_pd['M_star'] + objs_pd['mHI']),c = np.log10(objs_pd['massiveDist']), cmap = cmx, norm = cNorm,edgecolor = 'k')
    ax5.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax5.set_xscale('log')
    ax5.set_yscale('log')
    ax5.set_ylabel(r'M$_*$ + M$_{\mathrm{HI}}$ [M$_\odot$]')
    ax5.set_xlabel(r'M$_{\mathrm{vir}}$ [M$_\odot$]')
    ax5.axis([2e6, 2e11, 2e2, 5e9])
    cb = mpl.colorbar.ColorbarBase(ax5sub, cmap=cmx, norm=cNorm)
    cb.set_label(r"Log(Distance to Massive Galaxy/1 kpc)")
    plt.savefig(outfile_base + '_BMHM_rMassGal.png')

    fig5.clear()
    fig5 = plt.figure(5,figsize=(plt_width,plt_width*aspect_ratio))
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax5 = fig5.add_subplot(gs[0])
    ax5sub = fig5.add_subplot(gs[1])

    cNorm  = colors.Normalize(vmin=1.5, vmax = 4)    
    ax5.scatter(objs_pd_comb['Mpeak'],(objs_pd_comb['M_star'] + objs_pd_comb['mHI']),c = np.log10(objs_pd_comb['massiveDist']), cmap = cmx, norm = cNorm,edgecolor = 'k')
    ax5.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax5.set_xscale('log')
    ax5.set_yscale('log')
    ax5.set_ylabel(r'M$_*$ + M$_{\mathrm{HI}}$ [M$_\odot$]')
    ax5.set_xlabel(r'M$_{\mathrm{vir, peak}}$ [M$_\odot$]')
    ax5.axis([1e8, 2e11, 2e2, 5e9])
    cb = mpl.colorbar.ColorbarBase(ax5sub, cmap=cmx, norm=cNorm)
    cb.set_label(r"Log(Distance to Massive Galaxy/1 kpc)")
    plt.savefig(outfile_base + '_BMHM_rMassGal_Mpeak.png')
    
    #Baryon fraction vs distance
    #fig6.clear()
    fig6 = plt.figure(6,figsize=(plt_width,plt_width*aspect_ratio))
    gs = gridspec.GridSpec(1,1)
    ax6 = fig6.add_subplot(gs[0])

    ax6.scatter(objs_pd['massiveDist'],(objs_pd['M_star'] + objs_pd['mHI'])/objs_pd['mass'],facecolor = 'none',edgecolor = 'k')
    fbary = ax6.scatter(objs_pd['massiveDist'],(objs_pd['M_star'] + objs_pd['mHI'])/objs_pd['mass'],edgecolor = 'k',facecolor = 'k',alpha = 0.3)
    ax6.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    #ax6.scatter(objs_pd['massiveDist'],objs_pd['mHI']/objs_pd['mass'],facecolor = 'none',edgecolor = 'k')
    fstar = ax6.scatter(objs_pd['massiveDist'],objs_pd['M_star']/objs_pd['mass'],facecolor = 'none',edgecolor = 'red')
    ax6.set_xscale('log')
    ax6.set_yscale('log')
    ax6.set_ylabel(r'M/M$_{\mathrm{vir}}$')
    ax6.set_xlabel(r'Log(Distance to Massive Galaxy/1 kpc)')
    ax6.axis([20, 7e3, 1e-7, 0.2])
    ax6.legend([fbary,fstar],[r'(M$_*$ + M$_{\mathrm{HI}})$/M$_{\mathrm{vir}}$',r'M$_*$/M$_{\mathrm{vir}}$'],loc = 3)
    plt.savefig(outfile_base + '_fbary_rMassGal.png')
    
    
if __name__ == '__main__':
    if (socket.gethostname() == "quirm.math.grinnell.edu"):
        prefix = '/home/christenc/Data/Sims/'
        prefix_outfile = '/home/christenc/Figures/marvel/'
    else:
        prefix = '/home/christensen/Storage2/UW/MolecH/Cosmo/'
        prefix_outfile = '/home/christensen/Plots/marvel/'

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
    
    tfile_base_2 = 'h229.cosmo50PLK.3072gst5HbwK1BH'
    tfile_2 = prefix + 'h229.cosmo50PLK.3072g/h229.cosmo50PLK.3072gst5HbwK1BH/snapshots_200bkgdens/h229.cosmo50PLK.3072gst5HbwK1BH.004096' #
    tfile_2 = prefix + 'h229.cosmo50PLK.3072g/h229.cosmo50PLK.3072gst5HbwK1BH/h229.cosmo50PLK.3072gst5HbwK1BH.004096/h229.cosmo50PLK.3072gst5HbwK1BH.004096'    
    
    tfile_base_3 = 'h242.cosmo50PLK.3072gst5HbwK1BH'
    tfile_3 = prefix + 'h242.cosmo50PLK.3072g/h242.cosmo50PLK.3072gst5HbwK1BH/snapshots_200bkgdens/h242.cosmo50PLK.3072gst5HbwK1BH.004096' #
    tfile_3 = prefix + 'h242.cosmo50PLK.3072g/h242.cosmo50PLK.3072gst5HbwK1BH/h242.cosmo50PLK.3072gst5HbwK1BH.004096/h242.cosmo50PLK.3072gst5HbwK1BH.004096'

    tfile_base_4 = 'h329.cosmo50PLK.3072gst5HbwK1BH'
    tfile_4 = prefix + 'h329.cosmo50PLK.3072g/h329.cosmo50PLK.3072gst5HbwK1BH/snapshots_200bkgdens/h329.cosmo50PLK.3072gst5HbwK1BH.004096' #
    tfile_4 = prefix + 'h329.cosmo50PLK.3072g/h329.cosmo50PLK.3072gst5HbwK1BH/h329.cosmo50PLK.3072gst5HbwK1BH.004096/h329.cosmo50PLK.3072gst5HbwK1BH.004096'
    
    outfile_base = prefix_outfile + 'marvelJL'
    #SMHM_v_distance([tfile_1,tfile_2,tfile_3,tfile_4],outfile_base,[tfile_base_1,tfile_base_2,tfile_base_3,tfile_base_4])
    SMHM_v_distance([tfile_cm, tfile_e, tfile_r, tfile_s, tfile_1,tfile_2,tfile_3,tfile_4],outfile_base,[tfile_base_cm, tfile_base_e, tfile_base_r, tfile_base_s, tfile_base_1,tfile_base_2,tfile_base_3,tfile_base_4])
