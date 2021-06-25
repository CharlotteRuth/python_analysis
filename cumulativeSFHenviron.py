#Charlotte Christensen

# 8/2/19
#For the given simulations, calculate the star formation histories at different distances from a massive galaxy

#%run /home/christenc/Code/python/python_analysis/cumulativeSFHenviron
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

def cumulativeSFHenviron(tfiles,outfile_base,tfile_base,*halo_nums):
    presentation = True
    if presentation:
        outfile_base = outfile_base + '_pres'
        plt.style.use(['default','/home/christenc/.config/matplotlib/presentation.mplstyle'])
        plt_width = 16 #8 #inches
        aspect_ratio = 1.0/4.0
        legendsize = 16
        dpi = 100
        lw = mpl.rcParams['lines.linewidth'] - 1        
    else:
        plt.style.use(['default','/home/christenc/.config/matplotlib/article.mplstyle'])
        plt_width = 3.5 #inches
        aspect_ratio = 1.0/4.0
        legendsize = 5
        dpi = 300
        lw = mpl.rcParams['lines.linewidth']
        
    #min_mass = 1e9 #Minimum halo mass for analysis in solar masses
    min_nstar =  100 #100 #Minimum number of stars for
    min_noncontamFrac = 0.9
    min_massiveHalo = 10**11.5

    maxtime = 13.7
    
    mass_color_scale = 0
    use_virial_mass = 1 #1
    if mass_color_scale:
        if use_virial_mass:
            min_vmass = 7e7 #1e8
            max_vmass = 1e11 #3e12
            outfile_base = outfile_base + '_vmass'
        else: #use stellar mass instead
            min_vmass = min_nstar*5e3 #1e8
            max_vmass = 2.5e9 #3e12
            outfile_base = outfile_base + '_smass'
        outfile_base = outfile_base + '_masscolor'
        cmx = plt.get_cmap("viridis_r") 
            #values = range(0,20)
            #cNorm  = colors.Normalize(vmin=0, vmax = values[-1]+1)
            #scalarMap = cm.ScalarMappable(norm=cNorm, cmap=cmx)
            
        cNorm  = colors.Normalize(vmin=np.log10(min_vmass), vmax = math.ceil(np.log10(max_vmass)))
        scalarMap = cm.ScalarMappable(norm=cNorm, cmap=cmx)     
    else:
        cmx = plt.get_cmap("tab20")
        values = range(0,20)
        cNorm  = colors.Normalize(vmin=0, vmax = values[-1]+1)
        scalarMap = cm.ScalarMappable(norm=cNorm, cmap=cmx)

    objs_pd = None 
    for tfile in tfiles:
        objs_dat = []
        print(tfile)
        f=open(tfile + '.MAP.data', 'rb')
        while 1:
            try:
                objs_dat.append(pickle.load(f))
            except EOFError:
                break        
        f.close()
        s = pynbody.load(tfile)
        h_dummy = s.halos(dummy = True)
        loc = []
        fMhires_ahf = []
        nstar_ahf = []
        haloid_ahf = []
        for AHFhalo in h_dummy:
            properties = AHFhalo.properties
            fMhires_ahf.append(AHFhalo.properties['fMhires'])
            nstar_ahf.append(AHFhalo.properties['n_star'])
            haloid_ahf.append(AHFhalo.properties['halo_id'])
            if (properties['mass'] > min_massiveHalo):
                loc.append(np.array([properties['Xc']/properties['h'], properties['Yc']/properties['h'], properties['Zc']/properties['h']]))
        valid_halos_ahf = (np.array(haloid_ahf))[(np.array(fMhires_ahf) >= min_noncontamFrac) & (np.array(nstar_ahf) >= min_nstar)]
        loc = np.array(loc)
        massiveDist = [] #Distance to nearest massive (M_vir > min_massiveHalo) galaxy
        for halo in objs_dat:
            massiveDist.append(min(((halo['Xc'] - loc[:,0])**2 + (halo['Yc'] - loc[:,1])**2 + (halo['Zc'] - loc[:,2])**2)**(0.5)))
            
        sim = ''
        if tfile == '/home/christenc/Data/Sims/cptmarvel.cosmo25cmb/cptmarvel.cosmo25cmb.4096g5HbwK1BH/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096':
            sim = 'Cpt. Marvel'
        if tfile == '/home/christenc/Data/Sims/rogue.cosmo25cmb/rogue.cosmo25cmb.4096g5HbwK1BH/rogue.cosmo25cmb.4096g5HbwK1BH.004096/rogue.cosmo25cmb.4096g5HbwK1BH.004096':
            sim = 'Rogue'
        if tfile == '/home/christenc/Data/Sims/elektra.cosmo25cmb/elektra.cosmo25cmb.4096g5HbwK1BH/elektra.cosmo25cmb.4096g5HbwK1BH.004096/elektra.cosmo25cmb.4096g5HbwK1BH.004096':
            sim = 'Elektra'
        if tfile == '/home/christenc/Data/Sims/storm.cosmo25cmb/storm.cosmo25cmb.4096g5HbwK1BH/storm.cosmo25cmb.4096g5HbwK1BH.004096/storm.cosmo25cmb.4096g5HbwK1BH.004096':
            sim = 'Storm'
        if tfile == '/home/christenc/Data/Sims/h329.cosmo50PLK.3072g/h329.cosmo50PLK.3072gst5HbwK1BH/snapshots_200bkgdens/h329.cosmo50PLK.3072gst5HbwK1BH.004096':
            sim = 'Elena'
        if tfile == '/home/christenc/Data/Sims/h242.cosmo50PLK.3072g/h242.cosmo50PLK.3072gst5HbwK1BH/snapshots_200bkgdens/h242.cosmo50PLK.3072gst5HbwK1BH.004096':
            sim = 'Sonia'
        if tfile == '/home/christenc/Data/Sims/h229.cosmo50PLK.3072g/h229.cosmo50PLK.3072gst5HbwK1BH/snapshots_200bkgdens/h229.cosmo50PLK.3072gst5HbwK1BH.004096':
            sim = 'Ruth'
        if tfile == '/home/christenc/Data/Sims/h148.cosmo50PLK.3072g/h148.cosmo50PLK.3072g3HbwK1BH/snapshots_200bkgdens/h148.cosmo50PLK.3072g3HbwK1BH.004096':
            sim = 'Sandra'

        temp = pd.DataFrame(objs_dat)
        temp['sim'] = [sim]*len(objs_dat)
        temp['massiveDist'] = massiveDist
        valid_halos = np.array(temp[(temp['n_star'] >= min_nstar) &(temp['fMhires']>=min_noncontamFrac )]['haloid'])
        intersec = np.intersect1d(valid_halos,valid_halos_ahf)
        if len(valid_halos) is not len(intersec):
            print(tfile)
            print('Number of valid halos in *.data: ',len(valid_halos))
            print('Extra halos in *.dat')
            print(np.setdiff1d(valid_halos,intersec))
            print('Number of valid halos in AHF: ',len(valid_halos_ahf))
            print('Extra halos in AHF')
            print(np.setdiff1d(valid_halos_ahf,intersec))
        if objs_pd is None: 
            objs_pd = temp
        else:
            objs_pd = objs_pd.append(temp, ignore_index = True)
        
    fig1 = plt.figure(1,figsize=(plt_width,plt_width*aspect_ratio))
    if mass_color_scale:
        gs =  gridspec.GridSpec(1,2,width_ratios=[45,1])
    else:
        gs =  gridspec.GridSpec(1,1)
    
    gs00 = gridspec.GridSpecFromSubplotSpec(1, 4, subplot_spec=gs[0],wspace=0.0)
    ax1 = fig1.add_subplot(gs00[0])
    ax2 = fig1.add_subplot(gs00[1],sharey=ax1)
    ax3 = fig1.add_subplot(gs00[2],sharey=ax1)
    ax4 = fig1.add_subplot(gs00[3],sharey=ax1)
    if mass_color_scale:
        ax1sub = fig1.add_subplot(gs[1])
    axes = [ax1,ax2,ax3,ax4]
        
    legendlines = []

    maxradii_arr = [150,300,750,3000]
    minradii_arr = [0.1,maxradii_arr[0],maxradii_arr[1],maxradii_arr[2]]
    titles = ['r < 150 kpc','150 kpc < r < 300 kpc','300 kpc < r < 750 kpc','750 kpc < r']
    for minradii, maxradii, ax, title in zip(minradii_arr, maxradii_arr, axes, titles):
        i = 0
        labels = []
        print(maxradii, len(objs_pd[(objs_pd['h1dist'] < maxradii) & (objs_pd['h1dist'] >= minradii)]))
        print(max(objs_pd[(objs_pd['h1dist'] < maxradii) & (objs_pd['h1dist'] >= minradii)]['h1dist']))
        print(min(objs_pd[(objs_pd['h1dist'] < maxradii) & (objs_pd['h1dist'] >= minradii)]['h1dist']))
        for index, halo in objs_pd[(objs_pd['h1dist'] < maxradii) & (objs_pd['h1dist'] >= minradii)].iterrows(): #len(h)-1):
            if (halo['n_star'] < min_nstar):
                continue
            if (halo['fMhires'] < min_noncontamFrac):
                continue
            if (halo['n_star'] > 10*(halo['n_particles'] - halo['n_star'] - halo['n_gas'])):
                continue            

            if (halo['sSFR'] > 1e-11):
                quenched = 0
                lw_plot = lw #*2
                linestyle = '-'
            else:
                quenched = 1
                lw_plot = lw
                linestyle = '--'
            
#            if (halo['hostVirialR'] > 0):
#                linestyle = '--'
#            else:
#                linestyle = '-'
                
            labels.append(str(halo['haloid']) + ' ' + halo['sim'])
            dtime = halo['sfhbins'][1] - halo['sfhbins'][0] #Gyr
            sfhcum = np.cumsum(halo['sfh']*dtime)
            
            if mass_color_scale:
                if use_virial_mass:
                    colorVal = scalarMap.to_rgba(np.log10(halo['mass']))
                else:
                    colorVal = scalarMap.to_rgba(np.log10(halo['M_star']))
                legendline = ax.plot((halo['sfhbins'][:-1] + halo['sfhbins'][1:])/2,sfhcum/max(sfhcum),color = colorVal,linestyle = linestyle,linewidth = lw_plot)
            else:
                colorVal = scalarMap.to_rgba((i % 20) + 1)
                legendline = ax.plot((halo['sfhbins'][:-1] + halo['sfhbins'][1:])/2,sfhcum/max(sfhcum),color = colorVal,linestyle = linestyle,linewidth = lw_plot) #,c = i) #, cmap = 'tab20')

            ax.label_outer()
            legendlines.append(legendline)
            fig1.show()
            ax.set_ylim([0,1])
            ax.set_title(title,fontsize= 'x-small')
            i = i + 1
        if not mass_color_scale:
            ax.legend(labels,loc = 1,fontsize = 8) #'xx-small')
            
    ax2.set_xlabel('Time [Gyr]')
    ax1.set_ylabel(r'$M_*( t_{form} < t)/M_*$')
    #ax1.axis([0.1, 20, 0, 1])
    #ax1.set_xscale('log')
    #ax1.set_title(tfile_base)
    if mass_color_scale:
        cb = mpl.colorbar.ColorbarBase(ax1sub, cmap=cmx, norm=cNorm)
        if use_virial_mass:
            cb.set_label('Log Virial Mass [M$_\odot$]')
        else:
            cb.set_label('Log Stellar Mass [M$_\odot$]')

    plt.tight_layout()
    fig1.show()
    fig1.savefig(outfile_base + '_environ_sfhcum.png',dpi = dpi)
    fig1.clear()

    
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
    cumulativeSFHenviron([tfile_1,tfile_2,tfile_3,tfile_4],outfile_base,[tfile_base_1,tfile_base_2,tfile_base_3,tfile_base_4])
    #cumulativeSFHenviron([tfile_cm, tfile_e, tfile_r, tfile_s, tfile_1,tfile_2,tfile_3,tfile_4],outfile_base,[tfile_base_cm, tfile_base_e, tfile_base_r, tfile_base_s, tfile_base_1,tfile_base_2,tfile_base_3,tfile_base_4])

        

