#Plot the spatial distribution of the halos in a simulation

#Run with
#%run /home/christensen/Code/python/python_analysis/plt_halos.py
#or, on quirm,
#%run /home/christenc/Code/python/python_analysis/plt_halos.py
#ipython --pylab

import matplotlib.pylab as plt
import matplotlib.mlab as mlab
import matplotlib.colors as colors
import matplotlib.cm as cm
import pynbody
import math
import numpy as np
import socket

def plt_halos(tfile,outfilebase,tfile_base,*halo_nums):

    #min_mass = 1e9 #Minimum halo mass for analysis in solar masses
    min_nstar =  100 #Minimum number of stars for
    min_vmass = 1e8
    min_noncontamFrac = 0.9 #Mass fraction in low mass DM particles

    if halo_nums:
        halo_nums = halo_nums[0]
        
    s = pynbody.load(tfile)
    s.physical_units()
    h = s.halos()
    h_dummy = s.halos(dummy = True)

    mass_color_scale = 0

    if mass_color_scale:
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
    
    fig1 = plt.figure(1)
    fig2 = plt.figure(2)
    ax1 = fig1.add_subplot(111)
    ax2 = fig2.add_subplot(111)
    
    labels = []

    i = 0
    #Loop through all the halos, breaking when the virial mass is below min_mass
    xmin = h_dummy[1].properties['Xc'] - h_dummy[1].properties['Rvir']
    xmax = h_dummy[1].properties['Xc'] + h_dummy[1].properties['Rvir']
    ymin = h_dummy[1].properties['Yc'] - h_dummy[1].properties['Rvir']
    ymax = h_dummy[1].properties['Yc'] + h_dummy[1].properties['Rvir']
    zmin = h_dummy[1].properties['Zc'] - h_dummy[1].properties['Rvir']
    zmax = h_dummy[1].properties['Zc'] + h_dummy[1].properties['Rvir']    
    for ihalo in range(1,len(h)): #len(h)-1):
        if not halo_nums:
            valid_halo = h_dummy[ihalo].properties['n_star'] >= min_nstar and h[ihalo].properties['fMhires'] >= min_noncontamFrac and h[ihalo].properties['mass'] >= min_vmass
        else:
            valid_halo = str(h_dummy[ihalo].properties['halo_id']) in halo_nums

            
        if valid_halo:
            labels.append(str(ihalo))
            print(h_dummy[ihalo].properties['halo_id'])
            
            if mass_color_scale:
                colorVal = scalarMap.to_rgba(np.log10(h_dummy[ihalo].properties['mass']))
            else:
                colorVal = scalarMap.to_rgba(i)

            if h_dummy[ihalo].properties['Xc'] - h_dummy[ihalo].properties['Rvir'] < xmin:
                xmin = h_dummy[ihalo].properties['Xc'] - h_dummy[ihalo].properties['Rvir']
            if h_dummy[ihalo].properties['Xc'] + h_dummy[ihalo].properties['Rvir'] > xmax:
                xmax = h_dummy[ihalo].properties['Xc'] + h_dummy[ihalo].properties['Rvir']

            if h_dummy[ihalo].properties['Yc'] - h_dummy[ihalo].properties['Rvir'] < ymin:
                ymin = h_dummy[ihalo].properties['Yc'] - h_dummy[ihalo].properties['Rvir']
            if h_dummy[ihalo].properties['Yc'] + h_dummy[ihalo].properties['Rvir'] > ymax:
                ymax = h_dummy[ihalo].properties['Yc'] + h_dummy[ihalo].properties['Rvir']

            if h_dummy[ihalo].properties['Zc'] - h_dummy[ihalo].properties['Rvir'] < zmin:
                zmin = h_dummy[ihalo].properties['Zc'] - h_dummy[ihalo].properties['Rvir']
            if h_dummy[ihalo].properties['Zc'] + h_dummy[ihalo].properties['Rvir'] > zmax:
                zmax = h_dummy[ihalo].properties['Zc'] + h_dummy[ihalo].properties['Rvir']
                
            circlexy = plt.Circle((h_dummy[ihalo].properties['Xc'],h_dummy[ihalo].properties['Yc']),h_dummy[ihalo].properties['Rvir'],color = colorVal,fill=False)
            circlexz = plt.Circle((h_dummy[ihalo].properties['Xc'],h_dummy[ihalo].properties['Zc']),h_dummy[ihalo].properties['Rvir'],color = colorVal,fill=False)
            ax1.add_artist(circlexy)
            ax1.set_xlim(xmin,xmax)
            ax1.set_ylim(ymin,ymax)
            ax2.add_artist(circlexz)
            ax2.set_xlim(xmin,xmax)
            ax2.set_ylim(zmin,zmax)            
            fig1.show()
            fig2.show()
            i += 1
            
    ax1.set_xlabel('X')
    ax1.set_ylabel('Y')
    ax1.set_title(tfile_base)
    ax1.legend(labels,loc = 1)
    fig1.show()
    fig1.savefig(outfilebase + '.halos_xy.png')

    fig1.clear()

    
    ax2.set_xlabel('X')
    ax2.set_ylabel('Z')
    ax2.set_title(tfile_base)
    ax2.legend(labels,loc = 4)
    fig2.show()
    fig2.savefig(outfilebase + '.halos_xz.png')

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
    #plt_halos(tfile,outfile_base,tfile_base,halo_nums)

    halo_nums = ['1','3','7','8','10','11','12','16','17','18','30','32','34','36','61','123']
    tfile = prefix + 'rogue.cosmo25cmb/rogue.cosmo25cmb.4096g5HbwK1BH/rogue.cosmo25cmb.4096g5HbwK1BH.004096/rogue.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_base = 'rogue.cosmo25cmb.4096g5HbwK1BH'
    outfile_base = prefix_outfile + tfile_base
    #plt_halos(tfile,outfile_base,tfile_base,halo_nums)

    halo_nums = ['1','2','3','4','5','8','9','10','11','12','17','18','37','75']
    tfile = prefix + 'elektra.cosmo25cmb/elektra.cosmo25cmb.4096g5HbwK1BH/elektra.cosmo25cmb.4096g5HbwK1BH.004096/elektra.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_base = 'elektra.cosmo25cmb.4096g5HbwK1BH'
    outfile_base = prefix_outfile + tfile_base
    #plt_halos(tfile,outfile_base,tfile_base,halo_nums)

    halo_nums = ['1','2','3','4','5','6','7','8','10','11','13','14','15','16','17','23','24','28','34','35','43','49','50','60','109','124','125','192','208'] 
    tfile = prefix + 'storm.cosmo25cmb/storm.cosmo25cmb.4096g5HbwK1BH/storm.cosmo25cmb.4096g5HbwK1BH.004096/storm.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_base = 'storm.cosmo25cmb.4096g5HbwK1BH'
    outfile_base = prefix_outfile + tfile_base
    #plt_halos(tfile,outfile_base,tfile_base,halo_nums)
    
    tfile_base = 'h148.cosmo50PLK.3072g3HbwK1BH'
    tfile = prefix + 'h148.cosmo50PLK.3072g/h148.cosmo50PLK.3072g3HbwK1BH/h148.cosmo50PLK.3072g3HbwK1BH.004096/h148.cosmo50PLK.3072g3HbwK1BH.004096'    
    outfile_base = prefix_outfile + tfile_base
    halo_nums = ['1','2','3','5','6','9','10','11','12','14','18','23','26','28','31','34','36','42','57','64','77','94','125','160','252','264','271','304']
    plt_halos(tfile,outfile_base,tfile_base,halo_nums)
     
    tfile_base = 'h329.cosmo50PLK.3072gst5HbwK1BH.004096'
    tfile = prefix + 'h329.cosmo50PLK.3072g/h329.cosmo50PLK.3072gst5HbwK1BH/h329.cosmo50PLK.3072gst5HbwK1BH.004096/h329.cosmo50PLK.3072gst5HbwK1BH.004096'    
    outfile_base = prefix_outfile + tfile_base
    halo_nums = ['1','9','32','126','129']
    plt_halos(tfile,outfile_base,tfile_base,halo_nums)
    
    tfile_base = 'h229.cosmo50PLK.3072gst5HbwK1BH.004096'
    tfile = prefix + 'h229.cosmo50PLK.3072g/h229.cosmo50PLK.3072gst5HbwK1BH/h229.cosmo50PLK.3072gst5HbwK1BH.004096/h229.cosmo50PLK.3072gst5HbwK1BH.004096'    
    outfile_base = prefix_outfile + tfile_base
    halo_nums = ['1','2','4','7','17','21','22','27','51','52','70','104','203']
    plt_halos(tfile,outfile_base,tfile_base,halo_nums)
    
    tfile_base = 'h242.cosmo50PLK.3072gst5HbwK1BH.004096'
    tfile = prefix + 'h242.cosmo50PLK.3072g/h242.cosmo50PLK.3072gst5HbwK1BH/h242.cosmo50PLK.3072gst5HbwK1BH.004096/h242.cosmo50PLK.3072gst5HbwK1BH.004096'    
    outfile_base = prefix_outfile + tfile_base
    halo_nums = ['1','9','11','24','29','30','33','39','40','45','75','76']
    plt_halos(tfile,outfile_base,tfile_base,halo_nums)

