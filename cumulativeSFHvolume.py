#Charlotte Christensen

#10/12/16
#Calculate and plot the cumulative star formation histories of every halo above a given mass range in a simulation volume

#%run /home/christenc/Code/python/python_analysis/cumulativeSFHvolume

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

def cumulativeSFHvolume(tfiles,outfile_base,tfile_base,*halo_nums):
    presentation = True
    if presentation:
        outfile_base = outfile_base + '_pres'
        plt.style.use(['default','/home/christenc/.config/matplotlib/presentation.mplstyle'])
        plt_width = 8 #inches
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
        
    #min_mass = 1e9 #Minimum halo mass for analysis in solar masses
    min_nstar =  100 #100 #Minimum number of stars for
    #min_vmass = 1e8
    min_noncontamFrac = 0.9 #Mass fraction in low mass DM particles

    maxtime = 13.7
        
    mass_color_scale = 0
    use_virial_mass = 0 #1
    if mass_color_scale:
        if use_virial_mass:
            min_vmass = 7e7 #1e8
            max_vmass = 1e11 #3e12
        else: #use stellar mass instead
            min_vmass = min_nstar*5e3 #1e8
            max_vmass = 2.5e9 #3e12          
        outfile_base = outfile_base + '_masscolor'
        cmx = plt.get_cmap("viridis_r") 
        #values = range(0,20)
        #cNorm  = colors.Normalize(vmin=0, vmax = values[-1]+1)
        #scalarMap = cm.ScalarMappable(norm=cNorm, cmap=cmx)

        cNorm  = colors.Normalize(vmin=np.log10(min_vmass), vmax = np.log10(max_vmass))
        scalarMap = cm.ScalarMappable(norm=cNorm, cmap=cmx)
    else:
        cmx = plt.get_cmap("tab20")
        values = range(0,20)
        cNorm  = colors.Normalize(vmin=0, vmax = values[-1]+1)
        scalarMap = cm.ScalarMappable(norm=cNorm, cmap=cmx)

        
    #dtime = .1
    #s = pynbody.load(tfile)
    #s.physical_units()
    #h = s.halos()
    #h_dummy = s.halos(dummy = True)
    #properties = h_dummy[1].properties

    objs_pd = None 
    for tfile in tfiles:
        objs_dat = []
        f=open(tfile + '.MAP.data', 'rb')
        while 1:
            try:
                objs_dat.append(pickle.load(f))
            except EOFError:
                break       
        f.close()
        sim = ''
        if tfile == '/home/christenc/Data/Sims/h329.cosmo50PLK.3072g/h329.cosmo50PLK.3072gst5HbwK1BH/snapshots_200bkgdens/h329.cosmo50PLK.3072gst5HbwK1BH.004096/h329.cosmo50PLK.3072gst5HbwK1BH.004096':
            sim = 'Elena'
        if tfile == '/home/christenc/Data/Sims/h242.cosmo50PLK.3072g/h242.cosmo50PLK.3072gst5HbwK1BH/snapshots_200bkgdens/h242.cosmo50PLK.3072gst5HbwK1BH.004096':
            sim = 'Sonia'
        if tfile == '/home/christenc/Data/Sims/h229.cosmo50PLK.3072g/h229.cosmo50PLK.3072gst5HbwK1BH/snapshots_200bkgdens/h229.cosmo50PLK.3072gst5HbwK1BH.004096':
            sim = 'Ruth'
        if tfile == '/home/christenc/Data/Sims/h148.cosmo50PLK.3072g/h148.cosmo50PLK.3072g3HbwK1BH/snapshots_200bkgdens/h148.cosmo50PLK.3072g3HbwK1BH.004096':
            sim = 'Sandra'
            
        if objs_pd is None: 
            objs_pd = pd.DataFrame(objs_dat)
            objs_pd['sim'] = [sim]*len(objs_dat)
        else:
            temp = pd.DataFrame(objs_dat)
            temp['sim'] = [sim]*len(objs_dat)
            objs_pd = objs_pd.append(temp, ignore_index = True)
        print(objs_pd.shape)
           
    fig1 = plt.figure(1,figsize=(plt_width,plt_width*aspect_ratio))
    fig2 = plt.figure(2,figsize=(plt_width,plt_width*aspect_ratio))
    if mass_color_scale:
        gs =  gridspec.GridSpec(1,2,width_ratios=[15,1])
        ax1 = fig1.add_subplot(gs[0])
        ax1sub = fig1.add_subplot(gs[1])
        ax2 = fig2.add_subplot(gs[0])
        ax2sub = fig2.add_subplot(gs[1])
    else:
        ax1 = fig1.add_subplot(111)
        ax2 = fig2.add_subplot(111)

    i = 0    
    labels = []
    legendlines = []
    #if not halo_nums:
    #    halo_nums_all = range(1,len(h))
    #else:
    #    halo_nums_all = halo_nums
    
    #Loop through all the halos, breaking when the virial mass is below min_mass
    #for ihalo in halo_nums_all: #len(h)-1):
    for index, halo in objs_pd.iterrows(): #len(h)-1):    
        #properties = h_dummy[int(ihalo)].properties

        #if not halo_nums:
        if (halo['n_star'] < min_nstar):
            continue
        if (halo['fMhires'] < min_noncontamFrac):
            continue
        if (halo['n_star'] > 10*(halo['n_particles'] - halo['n_star'] - halo['n_gas'])):
            continue
        #else:
            #valid_halo = str(properties['halo_id']) in halo_nums

        #if int(objs_pd['HIgasfrac'][objs_pd['haloid'] == int(ihalo)] >= 0.25):
        if (halo['sSFR'] > 1e-11): 
            quenched = 0
            lw_plot = lw*2
        else:
            quenched = 1
            lw_plot = lw

        linestyle = '-'
        #if (halo['hostVirialR'] > 0):
        #    linestyle = '--'
        #else:
        #    linestyle = '-'
                
        labels.append(str(halo['haloid']))
        #halo = h.load_copy(int(ihalo))
        #print(properties['halo_id'],np.log10(properties['mass']),properties['fMhires'])
        #sfh, bin_edges = np.histogram(halo.star['tform'].in_units('1e9 yr'),range = [0,13.7],weights = halo.star['massform'].in_units('Msol'),bins = int(round(maxtime/dtime)))
        #sfhcum = np.cumsum(sfh/dtime/1e9)
        dtime = halo['sfhbins'][1] - halo['sfhbins'][0] #Gyr
        sfh = halo['sfh']
        sfhcum = np.cumsum(halo['sfh']*dtime)
            
        if mass_color_scale:
            if use_virial_mass:
                colorVal = scalarMap.to_rgba(np.log10(halo['mass']))
            else:
                colorVal = scalarMap.to_rgba(np.log10(halo['M_star']))            
            legendline = ax1.plot((halo['sfhbins'][:-1] + halo['sfhbins'][1:])/2,sfh/dtime/1e9,color = colorVal,linestyle = linestyle,linewidth = lw_plot)
            ax2.plot((halo['sfhbins'][:-1] + halo['sfhbins'][1:])/2,sfhcum/max(sfhcum),color = colorVal,linestyle = linestyle,linewidth = lw_plot)
        else:
            colorVal = scalarMap.to_rgba(i % 20)
            legendline = ax1.plot((halo['sfhbins'][:-1] + halo['sfhbins'][1:])/2,sfh/dtime/1e9,color = colorVal,linestyle = linestyle,linewidth = lw_plot) #,color = plt.cm.tab20(i)) #, cmap = 'tab20')
            ax2.plot((halo['sfhbins'][:-1] + halo['sfhbins'][1:])/2,sfhcum/max(sfhcum),color = colorVal,linestyle = linestyle,linewidth = lw_plot) #,c = i) #, cmap = 'tab20')

        legendlines.append(legendline)
        fig1.show()
        fig2.show()
        i += 1
            
    ax1.set_xlabel('Time [Gyr]')
    ax1.set_ylabel(r'M$_\odot$ yr$^{-1}$')
    #ax1.axis([0.1, 20, 0, 1])
    #ax1.set_xscale('log')
    #ax1.set_title(tfile_base)
    if mass_color_scale:
        cb = mpl.colorbar.ColorbarBase(ax1sub, cmap=cmx, norm=cNorm)
        cb.set_label('Log Virial Mass [M$_\odot$]')
    else:
        ax1.legend(labels,loc = 1,fontsize = 'xx-small')
    plt.tight_layout()
    fig1.show()
    fig1.savefig(outfile_base + '_sfh.png',dpi = dpi)
    
    ax2.set_xlabel('Time [Gyr]')
    ax2.set_ylabel(r'$M_*( t_{form} < t)/M_*$')
    ax2.axis([0, 13.7, 0, 1])
    #ax2.set_title(tfile_base)
    if mass_color_scale:
        cb = mpl.colorbar.ColorbarBase(ax2sub, cmap=cmx, norm=cNorm)
        cb.set_label('Log Virial Mass [M$_\odot$]')
    else:
        ax2.legend(labels,loc = 4,fontsize = 8) #'xx-small')
    #ax2.set_xscale('log')
    plt.tight_layout()
    fig2.show()
    fig2.savefig(outfile_base + '_sfhcum.png',dpi = dpi)
    fig1.clear()
    fig2.clear()

if __name__ == '__main__':
    if (socket.gethostname() == "quirm.math.grinnell.edu"):
        prefix = '/home/christenc/Data/Sims/'
        prefix_outfile = '/home/christenc/Figures/marvel/'
    else:
        prefix = '/home/christensen/Storage2/UW/MolecH/Cosmo/'
        prefix_outfile = '/home/christensen/Plots/marvel/'

    halo_nums = ['1','2','4','5','6','7','10','11','13','14','27']
    tfile_base = 'cptmarvel.cosmo25cmb.4096g5HbwK1BH'
    tfile = prefix + 'cptmarvel.cosmo25cmb/cptmarvel.cosmo25cmb.4096g5HbwK1BH/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096' #'cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096'
    outfile_base = prefix_outfile + tfile_base
    #cumulativeSFHenviron([tfile],outfile_base,tfile_base) #,halo_nums)

    halo_nums = ['1','3','7','8','10','11','12','16','17','18','30','32','34','36','61','123']
    tfile = prefix + 'rogue.cosmo25cmb/rogue.cosmo25cmb.4096g5HbwK1BH/rogue.cosmo25cmb.4096g5HbwK1BH.004096/rogue.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_base = 'rogue.cosmo25cmb.4096g5HbwK1BH'
    outfile_base = prefix_outfile + tfile_base
    #cumulativeSFHenviron([tfile],outfile_base,tfile_base) #,halo_nums)

    halo_nums = ['1','2','3','4','5','8','9','10','11','12','17','18','37','75']
    tfile = prefix + 'elektra.cosmo25cmb/elektra.cosmo25cmb.4096g5HbwK1BH/elektra.cosmo25cmb.4096g5HbwK1BH.004096/elektra.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_base = 'elektra.cosmo25cmb.4096g5HbwK1BH'
    outfile_base = prefix_outfile + tfile_base
    #cumulativeSFHenviron([tfile],outfile_base,tfile_base) #,halo_nums)

    halo_nums = ['1','2','3','4','5','6','7','8','10','11','13','14','15','16','17','23','24','28','34','35','43','49','50','60','109','124','125','192','208'] 
    tfile = prefix + 'storm.cosmo25cmb/storm.cosmo25cmb.4096g5HbwK1BH/storm.cosmo25cmb.4096g5HbwK1BH.004096/storm.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_base = 'storm.cosmo25cmb.4096g5HbwK1BH'
    outfile_base = prefix_outfile + tfile_base
    #cumulativeSFHenviron([tfile],outfile_base,tfile_base) #,halo_nums)
    
    tfile_base = 'h148.cosmo50PLK.3072g3HbwK1BH'
    tfile = prefix + 'h148.cosmo50PLK.3072g/h148.cosmo50PLK.3072g3HbwK1BH/snapshots_200bkgdens/h148.cosmo50PLK.3072g3HbwK1BH.004096'    
    outfile_base = prefix_outfile + tfile_base
    cumulativeSFHvolume([tfile],outfile_base,tfile_base)
    
    tfile_base = 'h229.cosmo50PLK.3072gst5HbwK1BH'
    tfile = prefix + 'h229.cosmo50PLK.3072g/h229.cosmo50PLK.3072gst5HbwK1BH/snapshots_200bkgdens/h229.cosmo50PLK.3072gst5HbwK1BH.004096'    
    outfile_base = prefix_outfile + tfile_base
    cumulativeSFHvolume([tfile],outfile_base,tfile_base)
    
    tfile_base = 'h242.cosmo50PLK.3072gst5HbwK1BH'
    tfile = prefix + 'h242.cosmo50PLK.3072g/h242.cosmo50PLK.3072gst5HbwK1BH/snapshots_200bkgdens/h242.cosmo50PLK.3072gst5HbwK1BH.004096'
    outfile_base = prefix_outfile + tfile_base
    cumulativeSFHvolume([tfile],outfile_base,tfile_base)

    tfile_base = 'h329.cosmo50PLK.3072gst5HbwK1BH'
    tfile = prefix + 'h329.cosmo50PLK.3072g/h329.cosmo50PLK.3072gst5HbwK1BH/snapshots_200bkgdens/h329.cosmo50PLK.3072gst5HbwK1BH.004096'
    outfile_base = prefix_outfile + tfile_base
    cumulativeSFHvolume([tfile],outfile_base,tfile_base)
    

