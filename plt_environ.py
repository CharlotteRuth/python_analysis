#Plot the dark matter distribution surrounding the main halo

#Run with
#%run /home/christensen/Code/python/python_analysis/plt_environ.py
#or, on quirm,
#%run /home/christenc/Code/python/python_analysis/plt_environ.py
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


def plt_environ(tfile,outfile_base,tfile_base,central_halo = '1'):
    plt_background = True
    presentation = False #True #False
    if presentation:
        outfile_base = outfile_base + 'marvel_pres'
        plt.style.use(['default','/home/christenc/.config/matplotlib/presentation.mplstyle'])
        plt_width = 8 #inches
        aspect_ratio = 1.0
        legendsize = 5 #16
        dpi = 100
        markersize = 40
        lw = mpl.rcParams['lines.linewidth'] - 1
    else:
        plt.style.use(['default','/home/christenc/.config/matplotlib/article.mplstyle'])
        plt_width = 3.5 #inches
        aspect_ratio = 1.0
        legendsize = 5
        dpi = 300
        markersize = 12
        lw = mpl.rcParams['lines.linewidth']

    s = pynbody.load(tfile)
    s.physical_units()
    h = s.halos()
    h_dummy = s.halos(dummy = True)
        
    extent = 2000
    xcenbox = 0
    ycenbox = 0
    zcenbox = 0
    fig1 = plt.figure(1,figsize=(plt_width,plt_width))
    fig2 = plt.figure(2,figsize=(plt_width,plt_width))
    ax1 = fig1.add_subplot(111)
    ax2 = fig2.add_subplot(111)
    
    pynbody.analysis.halo.center(h[1],mode = 'com')
    cube = pynbody.filt.Cuboid(str(xcenbox - extent/2) +' kpc', y1=str(ycenbox - extent/2) +' kpc', z1=str(zcenbox - extent/2)+' kpc', x2=str(xcenbox + extent/2)+' kpc', y2=str(ycenbox + extent/2)+' kpc', z2=str(zcenbox + extent/2)+' kpc')
    pynbody.plot.image(s[cube].d, width = 2*(xcenbox + extent/2),x1 = xcenbox - extent/2, y2 = ycenbox + extent/2, y1 = ycenbox - extent/2,units = 'Msol kpc^-2',subplot = ax1,cmap = cm.gray,show_cbar =False,vmin = 100)
    s.rotate_x(-90)
    cubexz = pynbody.filt.Cuboid(str(xcenbox - extent/2) +' kpc', z1=str(ycenbox - extent/2) +' kpc', y1=str(zcenbox - extent/2)+' kpc', x2=str(xcenbox + extent/2)+' kpc', z2=str(ycenbox + extent/2)+' kpc', y2=str(zcenbox + extent/2)+' kpc')
    pynbody.plot.image(s[cubexz].d, width = 2*(xcenbox + extent/2),x1 = xcenbox - extent/2, y2 = zcenbox + extent/2, y1 = zcenbox - extent/2,units = 'Msol kpc^-2',subplot = ax2,cmap = cm.gray,show_cbar =False,vmin = 100)

    ax1.set_xlim(xcenbox - extent/2,xcenbox + extent/2)
    ax1.set_ylim(ycenbox - extent/2,ycenbox + extent/2)
    ax2.set_xlim(xcenbox - extent/2,xcenbox + extent/2)
    ax2.set_ylim(zcenbox - extent/2,zcenbox + extent/2)

    ax1.set_xlabel('X [kpc]')
    ax1.set_ylabel('Y [kpc]')
    ax1.set_title(tfile_base)
    fig1.show()
    fig1.savefig(outfile_base + '.environ_xy.png')
    fig1.clear()

    ax2.set_xlabel('X [kpc]')
    ax2.set_ylabel('Z [kpc]')
    ax2.set_title(tfile_base)
    fig2.show()
    fig2.savefig(outfile_base + '.environ_xz.png')
    fig2.clear()
        
if __name__ == '__main__':
    if (socket.gethostname() == "quirm.math.grinnell.edu"):
        prefix = '/home/christenc/Data/Sims/'
        prefix_outfile = '/home/christenc/Figures/marvel/'
    else:
        prefix = '/home/christensen/Storage2/UW/MolecH/Cosmo/'
        prefix_outfile = '/home/christensen/Plots/marvel/'

    #Sandra
    tfile_base = 'h148.cosmo50PLK.3072g3HbwK1BH.004096'
    tfile = prefix + 'h148.cosmo50PLK.3072g/h148.cosmo50PLK.3072g3HbwK1BH/h148.cosmo50PLK.3072g3HbwK1BH.004096/h148.cosmo50PLK.3072g3HbwK1BH.004096'    
    outfile_base = prefix_outfile + tfile_base
    plt_environ(tfile,outfile_base,tfile_base) #,halo_nums)

    #Ruth
    tfile_base = 'h229.cosmo50PLK.3072gst5HbwK1BH.004096'
    tfile = prefix + 'h229.cosmo50PLK.3072g/h229.cosmo50PLK.3072gst5HbwK1BH/h229.cosmo50PLK.3072gst5HbwK1BH.004096/h229.cosmo50PLK.3072gst5HbwK1BH.004096'    
    outfile_base = prefix_outfile + tfile_base
    plt_environ(tfile,outfile_base,tfile_base) #,halo_nums)

    #Sonia
    tfile_base = 'h242.cosmo50PLK.3072gst5HbwK1BH.004096'
    tfile = prefix + 'h242.cosmo50PLK.3072g/h242.cosmo50PLK.3072gst5HbwK1BH/h242.cosmo50PLK.3072gst5HbwK1BH.004096/h242.cosmo50PLK.3072gst5HbwK1BH.004096'    
    outfile_base = prefix_outfile + tfile_base
    plt_environ(tfile,outfile_base,tfile_base) #,halo_nums)

    #Elena
    tfile_base = 'h329.cosmo50PLK.3072gst5HbwK1BH.004096'
    tfile = prefix + 'h329.cosmo50PLK.3072g/h329.cosmo50PLK.3072gst5HbwK1BH/h329.cosmo50PLK.3072gst5HbwK1BH.004096/h329.cosmo50PLK.3072gst5HbwK1BH.004096'    
    outfile_base = prefix_outfile + tfile_base
    plt_environ(tfile,outfile_base,tfile_base) #,halo_nums)
