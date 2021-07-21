# Charlotte Christensen
# 6/25/2021

# Read in Ferah's data about halos for the SMHM relation
import numpy as np
import socket
import sys, os, glob, pickle
import pandas as pd
import matplotlib.pylab as plt

if __name__ == '__main__':    
    if (socket.gethostname() == "quirm.math.grinnell.edu"):
        prefix = '/home/christenc/Data/Sims/'
        outprefix = '/home/christenc/Figures/marvel/'
        dataprefix = '/home/christenc/Code/Datafiles/'        
    else:
        prefix = '/home/christensen/Storage2/UW/MolecH/Cosmo/'
        outprefix = '/home/christensen/Storage2/UW/MolecH/Cosmo/'
        dataprefix = '/home/christensen/Code/Datafiles/'

    presentation = True
    if presentation:
        outbase = outprefix + 'marvel_pres_'
        plt.style.use(['default','/home/christenc/.config/matplotlib/presentation.mplstyle'])
        plt_width = 8 #inches
        aspect_ratio = 3.0/4.0
        legendsize = 12
        labelsize = 18
        dpi = 100
    else:
        outbase = outprefix + 'marvel'
        plt.style.use(['default','/home/christenc/.config/matplotlib/article.mplstyle'])
        plt_width = 3.5 #inches
        aspect_ratio = 3.0/4.0
        legendsize = 5
        labelsize = 5
        dpi = 300  

    dataprefix = '/home/christenc/Code/Datafiles/'
        
    f = open(dataprefix+'mstar_vs_mhalo_4Charlotte.txt', 'r')
    fdmdata = []
    for line in f:
        line = line.strip()
        columns = line.split()
        if len(columns) == 12 and columns[0] != 'Volume':
            source = {}
            source['simname'] = columns[0]
            source['halogrp_z0'] = columns[1]
            source['halogrp_Mpeak'] = columns[2]
            source['Mpeak_snap'] = columns[3]
            source['Mpeak'] = columns[4]
            source['Mhalo_z0'] = columns[5]
            source['Mstar_z0'] = columns[6]
            source['Mstar_z0_photo'] = columns[7]
            source['Mstar_Mpeak'] = columns[8]
            source['Mstar_Mpeak_z0'] = columns[9]
            source['Vmag'] = columns[10]
            source['type'] = columns[11]            
            fdmdata.append(source)
    f.close()
    fdmdata = pd.DataFrame(fdmdata)

    
