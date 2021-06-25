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
import matplotlib as mpl
import pynbody
import math
import numpy as np
import socket
import pandas as pd
import sys, os, glob, pickle

def plt_halos(tfile,outfile_base,tfile_base,central_halo = '1',*halo_nums):
    plt_background = True
    presentation = True #False
    if presentation:
        outfile_base = outfile_base + 'marvel_pres'
        plt.style.use(['default','/home/christenc/.config/matplotlib/presentation.mplstyle'])
        plt_width = 8 #inches
        aspect_ratio = 1.0
        legendsize = 5 #16
        dpi = 100
        markersize = 40
        lw = mpl.rcParams['lines.linewidth'] - 2
    else:
        plt.style.use(['default','/home/christenc/.config/matplotlib/article.mplstyle'])
        plt_width = 3.5 #inches
        aspect_ratio = 1.0
        legendsize = 5
        dpi = 300
        markersize = 12
        lw = mpl.rcParams['lines.linewidth'] - 1
    
    #min_mass = 1e9 #Minimum halo mass for analysis in solar masses
    min_nstar =  1 #100 #Minimum number of stars for plotting
    min_vmass = 1e8
    min_noncontamFrac = 0.9 #Mass fraction in low mass DM particles

    #if halo_nums:
    #    halo_nums = halo_nums[0]

    objs_dat = []
    f=open(tfile + '.MAP.data', 'rb')
    while 1:
        try:
            objs_dat.append(pickle.load(f))
        except EOFError:
            break
        
    f.close()
    objs_pd = pd.DataFrame(objs_dat)
    
    s = pynbody.load(tfile)
    s.physical_units()
    h = s.halos()
    h_dummy = s.halos(dummy = True)

    i = 0
    #Loop through all the halos, breaking when the virial mass is below min_mass or if the number of stars is too small
    properties = h_dummy[1].properties
    xmin_store = xmin = -1*properties['Rvir']/properties['h'] #properties['Xc'] - properties['Rvir']
    xmax_store = xmax = properties['Rvir']/properties['h'] #properties['Xc'] + properties['Rvir']
    ymin_store = ymin = -1*properties['Rvir']/properties['h'] #properties['Yc'] - properties['Rvir']
    ymax_store = ymax = properties['Rvir']/properties['h'] #properties['Yc'] + properties['Rvir']
    zmin_store = zmin = -1*properties['Rvir']/properties['h'] #properties['Zc'] - properties['Rvir']
    zmax_store = zmax = properties['Rvir']/properties['h'] #properties['Zc'] + properties['Rvir']
    xcen = properties['Xc']/properties['h']
    ycen = properties['Yc']/properties['h']
    zcen = properties['Zc']/properties['h']

    if not halo_nums:
        halo_nums_all = range(1,len(h))
    else:
        halo_nums_all = halo_nums

    #Loop through to find the extent of all the halos    
    for ihalo in halo_nums_all: #len(h)-1):
        #print(ihalo)
        properties = h_dummy[int(ihalo)].properties
        if not halo_nums:
            valid_halo = properties['n_star'] >= min_nstar and properties['fMhires'] >= min_noncontamFrac and properties['mass'] >= min_vmass
        else:
            valid_halo = 1 #str(properties['halo_id']) in halo_nums[i]
            
        if valid_halo:
            print('valid halo: '+str(ihalo))
                    
            print('Xc: ' + str(properties['Xc']/properties['h'] - xcen) + ', Yc: ' +  str(properties['Yc']/properties['h']- ycen) + ', Zc: ' + str(properties['Zc']/properties['h'] - zcen))

            if properties['Xc']/properties['h'] - properties['Rvir']/properties['h'] - xcen < xmin:
                xmin = properties['Xc']/properties['h'] - properties['Rvir']/properties['h'] - xcen
            if properties['Xc']/properties['h'] + properties['Rvir']/properties['h'] - xcen > xmax:
                xmax = properties['Xc']/properties['h'] + properties['Rvir']/properties['h'] - xcen

            if properties['Yc']/properties['h'] - properties['Rvir']/properties['h'] - ycen < ymin:
                ymin = properties['Yc']/properties['h'] - properties['Rvir']/properties['h'] - ycen
            if properties['Yc']/properties['h'] + properties['Rvir']/properties['h'] - ycen > ymax:
                ymax = properties['Yc']/properties['h'] + properties['Rvir']/properties['h'] - ycen

            if properties['Zc']/properties['h'] - properties['Rvir']/properties['h'] - zcen < zmin:
                zmin = properties['Zc']/properties['h'] - properties['Rvir']/properties['h'] - zcen
            if properties['Zc']/properties['h'] + properties['Rvir']/properties['h'] - zcen > zmax:
                zmax = properties['Zc'] /properties['h']+ properties['Rvir']/properties['h'] - zcen

            xmax_store = max([xmax_store,xmax])
            xmin_store = min([xmin_store,xmin])
            ymax_store = max([ymax_store,ymax])
            ymin_store = min([ymin_store,ymin])
            zmax_store = max([zmax_store,zmax])
            zmin_store = min([zmin_store,zmin])
            i += 1

    extent = max([xmax_store - xmin_store,ymax_store - ymin_store,zmax_store - zmin_store])*1.05
    xcenbox = (xmax_store + xmin_store)/2
    ycenbox = (ymax_store + ymin_store)/2
    zcenbox = (zmax_store + zmin_store)/2

    fig1 = plt.figure(1,figsize=(plt_width,plt_width))
    fig2 = plt.figure(2,figsize=(plt_width,plt_width))
    ax1 = fig1.add_subplot(111)
    ax2 = fig2.add_subplot(111)        
    if plt_background:
        pynbody.analysis.halo.center(h[1],mode = 'com')
        cube = pynbody.filt.Cuboid(str(xcenbox - extent/2) +' kpc', y1=str(ycenbox - extent/2) +' kpc', z1=str(zcenbox - extent/2)+' kpc', x2=str(xcenbox + extent/2)+' kpc', y2=str(ycenbox + extent/2)+' kpc', z2=str(zcenbox + extent/2)+' kpc')
        pynbody.plot.image(s[cube].d, width = 2*(xcenbox + extent/2),x1 = xcenbox - extent/2, y2 = ycenbox + extent/2, y1 = ycenbox - extent/2,units = 'Msol kpc^-2',subplot = ax1,cmap = cm.gray,show_cbar =False,vmin = 100) #, threaded = False)
        s.rotate_x(-90)
        cubexz = pynbody.filt.Cuboid(str(xcenbox - extent/2) +' kpc', z1=str(ycenbox - extent/2) +' kpc', y1=str(zcenbox - extent/2)+' kpc', x2=str(xcenbox + extent/2)+' kpc', z2=str(ycenbox + extent/2)+' kpc', y2=str(zcenbox + extent/2)+' kpc')
        pynbody.plot.image(s[cubexz].d, width = 2*(xcenbox + extent/2),x1 = xcenbox - extent/2, y2 = zcenbox + extent/2, y1 = zcenbox - extent/2,units = 'Msol kpc^-2',subplot = ax2,cmap = cm.gray,show_cbar =False,vmin = 100) #, threaded = False)

    ax1.set_xlim(xcenbox - extent/2,xcenbox + extent/2)
    ax1.set_ylim(ycenbox - extent/2,ycenbox + extent/2)
    ax2.set_xlim(xcenbox - extent/2,xcenbox + extent/2)
    ax2.set_ylim(zcenbox - extent/2,zcenbox + extent/2)

    i = 0
    labels = []
    legendlines = []

    mass_color_scale = 1
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
        
    for ihalo in halo_nums_all: #len(h)-1):
        #print(ihalo)
        properties = h_dummy[int(ihalo)].properties
        if not halo_nums:
            valid_halo = properties['n_star'] >= min_nstar and properties['fMhires'] >= min_noncontamFrac and properties['mass'] >= min_vmass
        else:
            valid_halo = 1 #str(properties['halo_id']) in halo_nums[i]
            
        if valid_halo:
            print('valid halo: '+str(ihalo))
            labels.append(str(ihalo))

            quenched = 1
            lw_plot = lw
            linestyle = '--'
            try:
                if int(objs_pd['sSFR'][objs_pd['haloid'] == int(ihalo)] >= 1e-11):
                    quenched = 0
                    #lw_plot = lw*2
                    linestyle = '-'
                if int(objs_pd['HIgasfrac'][objs_pd['haloid'] == int(ihalo)] >= 0.25):
                    quenched = 0
                    #lw_plot = lw*2
                    linestyle = '-'
            except:
                print('halo: '+ihalo+' is not found')
            
            if mass_color_scale:
                colorVal = scalarMap.to_rgba(np.log10(properties['mass']))
            else:
                colorVal = scalarMap.to_rgba(i % 20)
                
            if tfile_base == 'h242.cosmo50PLK.3072gst5HbwK1BH.004096' or tfile_base == 'h148.cosmo50PLK.3072g3HbwK1BH.004096' or tfile_base == 'elektra.cosmo25cmb.4096g5HbwK1BH':
                if properties['hostHalo'] == -1:
                    #linestyle = '-'
                    lw_plot = lw*2
                else:
                    #linestyle = '--'
                    lw_plot = lw
            else:
                if properties['hostHalo'] == 0:
                    #linestyle = '-'
                    lw_plot = lw*2
                else:
                    #linestyle = '--'
                    lw_plot = lw
            circlexy = plt.Circle((properties['Xc']/properties['h'] - xcen,properties['Yc']/properties['h'] - ycen),properties['Rvir']/properties['h'],color = colorVal,fill=False,linestyle = linestyle,lw = lw_plot)
            circlexz = plt.Circle((properties['Xc']/properties['h'] - xcen,properties['Zc']/properties['h'] - zcen),properties['Rvir']/properties['h'],color = colorVal,fill=False,linestyle = linestyle,lw = lw_plot)
            legendlines.append(circlexy)
            ax1.add_patch(mpl.patches.Circle((properties['Xc']/properties['h'] - xcen,properties['Yc']/properties['h'] - ycen),properties['Rvir']/properties['h'],color = colorVal,fill=False,linestyle = linestyle,lw = lw_plot))
            ax2.add_patch(mpl.patches.Circle((properties['Xc']/properties['h'] - xcen,properties['Zc']/properties['h'] - zcen),properties['Rvir']/properties['h'],color = colorVal,fill=False,linestyle = linestyle,lw = lw_plot))
           
            fig1.show()
            fig2.show()
            i += 1
        
    ax1.set_xlabel('X [kpc]')
    ax1.set_ylabel('Y [kpc]')
    ax1.set_title(tfile_base)
    if not mass_color_scale:
        ax1.legend(legendlines,labels,loc = 2,fontsize = 'xx-small') #loc = 2 for Elena
    fig1.show()
    if not mass_color_scale:
        fig1.savefig(outfile_base + '.halos_xy.png')
    else:
        fig1.savefig(outfile_base + '.halos_xy.mscale.png')
    fig1.clear()

    ax2.set_xlabel('X [kpc]')
    ax2.set_ylabel('Z [kpc]')
    ax2.set_title(tfile_base)
    if not mass_color_scale:    
        ax2.legend(legendlines,labels,loc = 2,fontsize = 'xx-small') #loc = 2 for Elena
    fig2.show()
    if not mass_color_scale:
        fig2.savefig(outfile_base + '.halos_xz.png')
    else:
        fig2.savefig(outfile_base + '.halos_xz.mscale.png')
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
    #plt_halos(tfile,outfile_base,tfile_base) #,halo_nums)

    halo_nums = ['1','3','7','8','10','11','12','16','17','18','30','32','34','36','61','123']
    tfile = prefix + 'rogue.cosmo25cmb/rogue.cosmo25cmb.4096g5HbwK1BH/rogue.cosmo25cmb.4096g5HbwK1BH.004096/rogue.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_base = 'rogue.cosmo25cmb.4096g5HbwK1BH'
    outfile_base = prefix_outfile + tfile_base
    #plt_halos(tfile,outfile_base,tfile_base) #,halo_nums)

    halo_nums = ['1','2','3','4','5','8','9','10','11','12','17','18','37','75']
    tfile = prefix + 'elektra.cosmo25cmb/elektra.cosmo25cmb.4096g5HbwK1BH/elektra.cosmo25cmb.4096g5HbwK1BH.004096/elektra.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_base = 'elektra.cosmo25cmb.4096g5HbwK1BH'
    outfile_base = prefix_outfile + tfile_base
    #plt_halos(tfile,outfile_base,tfile_base) #,halo_nums)

    halo_nums = ['1','2','3','4','5','6','7','8','10','11','13','14','15','16','17','23','24','28','34','35','43','49','50','60','109','124','125','192','208'] 
    tfile = prefix + 'storm.cosmo25cmb/storm.cosmo25cmb.4096g5HbwK1BH/storm.cosmo25cmb.4096g5HbwK1BH.004096/storm.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_base = 'storm.cosmo25cmb.4096g5HbwK1BH'
    outfile_base = prefix_outfile + tfile_base
    #plt_halos(tfile,outfile_base,tfile_base) #,halo_nums)

    #Sandra
    tfile_base = 'h148.cosmo50PLK.3072g3HbwK1BH.004096'
    tfile = prefix + 'h148.cosmo50PLK.3072g/h148.cosmo50PLK.3072g3HbwK1BH/snapshots_200bkgdens/h148.cosmo50PLK.3072g3HbwK1BH.004096'
    #h148.cosmo50PLK.3072g3HbwK1BH.004096/h148.cosmo50PLK.3072g3HbwK1BH.004096'    
    outfile_base = prefix_outfile + tfile_base
    halo_nums = ['1','2','3','5','6','9','10','11','12','14','18','23','26','28','31','34','36','42','57','64','77','94','125','160','252','264','271','304']
    #plt_halos(tfile,outfile_base,tfile_base) #,halo_nums)

    #Ruth
    tfile_base = 'h229.cosmo50PLK.3072gst5HbwK1BH.004096'
    tfile = prefix + 'h229.cosmo50PLK.3072g/h229.cosmo50PLK.3072gst5HbwK1BH/snapshots_200bkgdens/h229.cosmo50PLK.3072gst5HbwK1BH.004096'
    #h229.cosmo50PLK.3072gst5HbwK1BH.004096/h229.cosmo50PLK.3072gst5HbwK1BH.004096'    
    outfile_base = prefix_outfile + tfile_base
    halo_nums = ['1','2','4','7','17','21','22','27','51','52','70','104','203']
    #plt_halos(tfile,outfile_base,tfile_base) #,halo_nums)

    #Sonia
    tfile_base = 'h242.cosmo50PLK.3072gst5HbwK1BH.004096'
    tfile = prefix + 'h242.cosmo50PLK.3072g/h242.cosmo50PLK.3072gst5HbwK1BH/h242.cosmo50PLK.3072gst5HbwK1BH.004096/h242.cosmo50PLK.3072gst5HbwK1BH.004096'
    #snapshots_200bkgdens/h242.cosmo50PLK.3072gst5HbwK1BH.004096'
    #h242.cosmo50PLK.3072gst5HbwK1BH.004096/h242.cosmo50PLK.3072gst5HbwK1BH.004096'    
    outfile_base = prefix_outfile + tfile_base
    halo_nums = ['1','9','11','24','29','30','33','39','40','45','75','76']
    plt_halos(tfile,outfile_base,tfile_base) #,halo_nums)

    #Elena
    tfile_base = 'h329.cosmo50PLK.3072gst5HbwK1BH.004096'
    tfile = prefix + 'h329.cosmo50PLK.3072g/h329.cosmo50PLK.3072gst5HbwK1BH/snapshots_200bkgdens/h329.cosmo50PLK.3072gst5HbwK1BH.004096'
    #h329.cosmo50PLK.3072gst5HbwK1BH.004096/h329.cosmo50PLK.3072gst5HbwK1BH.004096'    
    outfile_base = prefix_outfile + tfile_base
    halo_nums = ['1','9','32','126','129']
    halo_nums = ['1','9','31','32','40','63','99','126','129','135','170']    
    #plt_halos(tfile,outfile_base,tfile_base) #,halo_nums)
