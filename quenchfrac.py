import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors as colors
import socket
import pandas as pd
import sys, os, glob, pickle

if (socket.gethostname() == "quirm.math.grinnell.edu"):
    prefix = '/home/christenc/Code/Justice_League_Code/Data/'
    outprefix = '/home/christenc/Figures/marvel/'
    dataprefix = '/home/christenc/Code/Datafiles/'
else:
    prefix = '/home/christensen/Storage2/UW/MolecH/Cosmo/'
    outprefix = '/home/christensen/Storage2/UW/MolecH/Cosmo/'
    dataprefix = '/home/christensen/Code/Datafiles/'

presentation = True
if presentation:
    outbase = outprefix + 'jl_pres_'
    plt.style.use(['default','/home/christenc/.config/matplotlib/presentation.mplstyle'])
    plt_width = 8 #inches
    aspect_ratio = 3.0/4.0
    legendsize = 16
    dpi = 100
    markersize = 40
else:
    outbase = outprefix + 'jl'
    plt.style.use(['default','/home/christenc/.config/matplotlib/article.mplstyle'])
    plt_width = 3.5 #inches
    aspect_ratio = 3.0/4.0
    legendsize = 5
    dpi = 300
    markersize = 12

objs_sandra = []
f=open(prefix + 'z0_data/h148.data', 'rb')
#f=open(prefix + 'timesteps_data/h148.data', 'rb')
while 1:
    try:
        objs_sandra.append(pickle.load(f))
    except EOFError:
        break
f.close()

objs_elena = []
f=open(prefix + 'z0_data/h329.data', 'rb')
#f=open(prefix + 'timesteps_data/h329.data', 'rb')
while 1:
    try:
        objs_elena.append(pickle.load(f))
    except EOFError:
        break
f.close()

