#Charlotte Christensen

#8/13/19
#Plot the SMHM relation for the marvel and Justice League runs, coloring points according to distance/tau_90

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

    min_massiveHalo = 10**11.5

    if (socket.gethostname() == "quirm.math.grinnell.edu"):
        dataprefix = '/home/christenc/Code/Datafiles/'
    else:
        dataprefix = '/home/christensen/Code/Datafiles/'
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
    
    
    objs_pd = None 
    for tfile in tfiles:
        objs_dat = []
        print(tfile)
        f=open(tfile + '.data', 'rb')
        while 1:
            try:
                objs_dat.append(pickle.load(f))
            except EOFError:
                break        
        f.close()
        s = pynbody.load(tfile)
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

    #SMHM colored by tau_90
    #fig1.clear()
    fig1 = plt.figure(1,figsize=(plt_width,plt_width*aspect_ratio))
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax1 = fig1.add_subplot(gs[0])
    ax1sub = fig1.add_subplot(gs[1])

    cmx = plt.get_cmap("viridis") 
    cNorm  = colors.Normalize(vmin=0, vmax = 14)    
    sat_v_tau90 = ax1.scatter(objs_pd['mass'][~np.isnan(objs_pd['hostDist'])],objs_pd['M_star'][~np.isnan(objs_pd['hostDist'])],c = tau90[~np.isnan(objs_pd['hostDist'])], cmap = cmx, norm = cNorm,edgecolor = 'k')
    cen_v_tau90 = ax1.scatter(objs_pd['mass'][np.isnan(objs_pd['hostDist'])],objs_pd['M_star'][np.isnan(objs_pd['hostDist'])],c = tau90[np.isnan(objs_pd['hostDist'])], cmap = cmx, norm = cNorm,edgecolor = 'k',marker = 'D')
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
    #tfile_1 = prefix + 'h148.cosmo50PLK.3072g/h148.cosmo50PLK.3072g3HbwK1BH/h148.cosmo50PLK.3072g3HbwK1BH.004096/h148.cosmo50PLK.3072g3HbwK1BH.004096'    
    
    tfile_base_2 = 'h229.cosmo50PLK.3072gst5HbwK1BH'
    tfile_2 = prefix + 'h229.cosmo50PLK.3072g/h229.cosmo50PLK.3072gst5HbwK1BH/snapshots_200bkgdens/h229.cosmo50PLK.3072gst5HbwK1BH.004096' #
    #tfile_2 = prefix + 'h229.cosmo50PLK.3072g/h229.cosmo50PLK.3072gst5HbwK1BH/h229.cosmo50PLK.3072gst5HbwK1BH.004096/h229.cosmo50PLK.3072gst5HbwK1BH.004096'    
    
    tfile_base_3 = 'h242.cosmo50PLK.3072gst5HbwK1BH'
    tfile_3 = prefix + 'h242.cosmo50PLK.3072g/h242.cosmo50PLK.3072gst5HbwK1BH/snapshots_200bkgdens/h242.cosmo50PLK.3072gst5HbwK1BH.004096' #
    #tfile_3 = prefix + 'h242.cosmo50PLK.3072g/h242.cosmo50PLK.3072gst5HbwK1BH/h242.cosmo50PLK.3072gst5HbwK1BH.004096/h242.cosmo50PLK.3072gst5HbwK1BH.004096'

    tfile_base_4 = 'h329.cosmo50PLK.3072gst5HbwK1BH'
    tfile_4 = prefix + 'h329.cosmo50PLK.3072g/h329.cosmo50PLK.3072gst5HbwK1BH/snapshots_200bkgdens/h329.cosmo50PLK.3072gst5HbwK1BH.004096' #
    #tfile_4 = prefix + 'h329.cosmo50PLK.3072g/h329.cosmo50PLK.3072gst5HbwK1BH/h329.cosmo50PLK.3072gst5HbwK1BH.004096/h329.cosmo50PLK.3072gst5HbwK1BH.004096'
    
    outfile_base = prefix_outfile + 'marvelJL'
    SMHM_v_distance([tfile_1,tfile_2,tfile_3,tfile_4],outfile_base,[tfile_base_1,tfile_base_2,tfile_base_3,tfile_base_4])
    #SMHM_v_distance([tfile_cm, tfile_e, tfile_r, tfile_s, tfile_1,tfile_2,tfile_3,tfile_4],outfile_base,[tfile_base_cm, tfile_base_e, tfile_base_r, tfile_base_s, tfile_base_1,tfile_base_2,tfile_base_3,tfile_base_4])
