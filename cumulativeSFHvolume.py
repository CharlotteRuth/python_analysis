#Charlotte Christensen

# 10/12/16
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

def cumulativeSFHvolume(tfile,outfile_base,tfile_base,*halo_nums):
    presentation = True
    if presentation:
        outfile_base = outfile_base + '_pres'
        plt.style.use(['default','/home/christenc/.config/matplotlib/presentation.mplstyle'])
        plt_width = 8 #inches
        aspect_ratio = 3.0/4.0
        legendsize = 16
        dpi = 100
    else:
        plt.style.use(['default','/home/christenc/.config/matplotlib/article.mplstyle'])
        plt_width = 3.5 #inches
        aspect_ratio = 3.0/4.0
        legendsize = 5
        dpi = 300
    
    #min_mass = 1e9 #Minimum halo mass for analysis in solar masses
    min_nstar =  100 #Minimum number of stars for
    min_vmass = 1e8
    min_noncontamFrac = 0.9 #Mass fraction in low mass DM particles

    if halo_nums:
        halo_nums = halo_nums[0]

    maxtime = 13.7
    dtime = .1
    s = pynbody.load(tfile)
    s.physical_units()
    h = s.halos()
    h_dummy = s.halos(dummy = True)

    mass_color_scale = 1    
    if mass_color_scale:
        outfile_base = outfile_base + '_masscolor'
        cmx = plt.get_cmap("viridis_r") 
        #values = range(0,20)
        #cNorm  = colors.Normalize(vmin=0, vmax = values[-1]+1)
        #scalarMap = cm.ScalarMappable(norm=cNorm, cmap=cmx)

        cNorm  = colors.Normalize(vmin=np.log10(min_vmass), vmax = math.ceil(np.log10(h[1].properties['mass'])))
        scalarMap = cm.ScalarMappable(norm=cNorm, cmap=cmx)
    else:
        cmx = plt.get_cmap("tab20")
        values = range(0,20)
        cNorm  = colors.Normalize(vmin=0, vmax = values[-1]+1)
        scalarMap = cm.ScalarMappable(norm=cNorm, cmap=cmx)

        
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
    
    labels = []

    i = 0
    #Loop through all the halos, breaking when the virial mass is below min_mass
    for ihalo in range(1,len(h)): #len(h)-1):
        if not halo_nums:
            valid_halo = h_dummy[ihalo].properties['n_star'] >= min_nstar and h[ihalo].properties['fMhires'] >= min_noncontamFrac and h[ihalo].properties['mass'] >= min_vmass
        else:
            valid_halo = str(h_dummy[ihalo].properties['halo_id']) in halo_nums

            
        if valid_halo:
            labels.append(str(ihalo))
            halo = h.load_copy(ihalo)
            print(h[ihalo].properties['halo_id'],np.log10(h[ihalo].properties['mass']),h[ihalo].properties['fMhires'])
            sfh, bin_edges = np.histogram(halo.star['tform'].in_units('1e9 yr'),range = [0,13.7],weights = halo.star['massform'].in_units('Msol'),bins = int(round(maxtime/dtime)))
            sfhcum = np.cumsum(sfh/dtime/1e9)
            
            if mass_color_scale:
                colorVal = scalarMap.to_rgba(np.log10(h[ihalo].properties['mass']))
                ax1.plot((bin_edges[:-1] + bin_edges[1:])/2,sfh/dtime/1e9,color = colorVal)
                ax2.plot((bin_edges[:-1] + bin_edges[1:])/2,sfhcum/max(sfhcum),color = colorVal)
            else:
                colorVal = scalarMap.to_rgba(i % 20)
                ax1.plot((bin_edges[:-1] + bin_edges[1:])/2,sfh/dtime/1e9,color = colorVal) #,color = plt.cm.tab20(i)) #, cmap = 'tab20')
                ax2.plot((bin_edges[:-1] + bin_edges[1:])/2,sfhcum/max(sfhcum),color = colorVal) #,c = i) #, cmap = 'tab20')
                
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
        ax1.legend(labels,loc = 1)
    plt.tight_layout()
    fig1.show()
    fig1.savefig(outfile_base + '_sfh.png',dpi = dpi)

    fig1.clear()

    
    ax2.set_xlabel('Time [Gyr]')
    ax2.set_ylabel(r'$M_*( t_{form} < t)/M_*$')
    ax2.axis([0, 13.7, 0, 1])
    #ax2.set_title(tfile_base)
    if mass_color_scale:
        cb = mpl.colorbar.ColorbarBase(ax2sub, cmap=cmx, norm=cNorm)
        cb.set_label('Log Virial Mass [M$_\odot$]')
    else:
        ax2.legend(labels,loc = 4)
    #ax2.set_xscale('log')
    plt.tight_layout()
    fig2.show()
    fig2.savefig(outfile_base + '_sfhcum.png',dpi = dpi)

    fig2.clear()
        

if __name__ == '__main__':
    if (socket.gethostname() == "quirm"):
        prefix = '/home/christenc/Data/Sims/'
        prefix_outfile = '/home/christenc/Figures/marvel/'
    else:
        prefix = '/home/christensen/Storage2/UW/MolecH/Cosmo/'
        prefix_outfile = '/home/christensen/Plots/marvel/'

    halo_nums = ['1','2','4','5','6','7','10','11','13','14','27']
    tfile_base = 'cptmarvel.cosmo25cmb.4096g5HbwK1BH'
    tfile = prefix + 'cptmarvel.cosmo25cmb/cptmarvel.cosmo25cmb.4096g5HbwK1BH/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096' #'cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096'
    outfile_base = prefix_outfile + tfile_base
    cumulativeSFHvolume(tfile,outfile_base,tfile_base,halo_nums)

    halo_nums = ['1','3','7','8','10','11','12','16','17','18','30','32','34','36','61','123']
    tfile = prefix + 'rogue.cosmo25cmb/rogue.cosmo25cmb.4096g5HbwK1BH/rogue.cosmo25cmb.4096g5HbwK1BH.004096/rogue.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_base = 'rogue.cosmo25cmb.4096g5HbwK1BH'
    outfile_base = prefix_outfile + tfile_base
    cumulativeSFHvolume(tfile,outfile_base,tfile_base,halo_nums)

    halo_nums = ['1','2','3','4','5','8','9','10','11','12','17','18','37','75']
    tfile = prefix + 'elektra.cosmo25cmb/elektra.cosmo25cmb.4096g5HbwK1BH/elektra.cosmo25cmb.4096g5HbwK1BH.004096/elektra.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_base = 'elektra.cosmo25cmb.4096g5HbwK1BH'
    outfile_base = prefix_outfile + tfile_base
    cumulativeSFHvolume(tfile,outfile_base,tfile_base,halo_nums)

    halo_nums = ['1','2','3','4','5','6','7','8','10','11','13','14','15','16','17','23','24','28','34','35','43','49','50','60','109','124','125','192','208'] 
    tfile = prefix + 'storm.cosmo25cmb/storm.cosmo25cmb.4096g5HbwK1BH/storm.cosmo25cmb.4096g5HbwK1BH.004096/storm.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_base = 'storm.cosmo25cmb.4096g5HbwK1BH'
    outfile_base = prefix_outfile + tfile_base
    cumulativeSFHvolume(tfile,outfile_base,tfile_base,halo_nums)
    
    tfile_base = 'h148.cosmo50PLK.3072g3HbwK1BH'
    tfile = prefix + 'h148.cosmo50PLK.3072g/h148.cosmo50PLK.3072g3HbwK1BH/h148.cosmo50PLK.3072g3HbwK1BH.004096/h148.cosmo50PLK.3072g3HbwK1BH.004096'    
    outfile_base = prefix_outfile + tfile_base
    #cumulativeSFHvolume(tfile,outfile_base,tfile_base)
     
    tfile_base = 'h329.cosmo50PLK.3072gst5HbwK1BH.004096'
    tfile = prefix + 'h329.cosmo50PLK.3072g/h329.cosmo50PLK.3072gst5HbwK1BH/h329.cosmo50PLK.3072gst5HbwK1BH.004096/h329.cosmo50PLK.3072gst5HbwK1BH.004096'
    outfile_base = prefix_outfile + tfile_base
    #cumulativeSFHvolume(tfile,outfile_base,tfile_base)
    
    tfile_base = 'h229.cosmo50PLK.3072gst5HbwK1BH.004096'
    tfile = prefix + 'h229.cosmo50PLK.3072g/h229.cosmo50PLK.3072gst5HbwK1BH/h229.cosmo50PLK.3072gst5HbwK1BH.004096/h229.cosmo50PLK.3072gst5HbwK1BH.004096'    
    outfile_base = prefix_outfile + tfile_base
    #cumulativeSFHvolume(tfile,outfile_base,tfile_base)
    
    tfile_base = 'h242.cosmo50PLK.3072gst5HbwK1BH.004096'
    tfile = prefix + 'h242.cosmo50PLK.3072g/h242.cosmo50PLK.3072gst5HbwK1BH/h242.cosmo50PLK.3072gst5HbwK1BH.004096/h242.cosmo50PLK.3072gst5HbwK1BH.004096'
    outfile_base = prefix_outfile + tfile_base
    #cumulativeSFHvolume(tfile,outfile_base,tfile_base)

