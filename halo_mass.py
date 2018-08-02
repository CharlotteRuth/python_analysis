import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors as colors
import socket
import pandas as pd
import sys, os, glob, pickle

if (socket.gethostname() == "quirm"):
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
    legendsize = 16
    dpi = 100
else:
    outbase = outprefix + 'marvel'
    plt.style.use(['default','/home/christenc/.config/matplotlib/article.mplstyle'])
    plt_width = 3.5 #inches
    aspect_ratio = 3.0/4.0
    legendsize = 5
    dpi = 300
    
"""
SMALL_SIZE = 8
MEDIUM_SIZE = 10
BIGGER_SIZE = 12
plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
"""
f_bar = 0.16510

#Markers and colors for the simulations
markers = ['o','o','o','o','o','o','o','o']
colormap = 'tab20'
colors_cm = 'b'
colors_r = 'b'
colors_e = 'b'
colors_s = 'b'
colors_sandra = 'b'
colors_ruth = 'b'
colors_sonia = 'b'
colors_elena = 'b'
colors_cm_f = 'none'
colors_r_f = 'none'
colors_e_f = 'none'
colors_s_f = 'none'
colors_sandra_f = 'none'
colors_ruth_f = 'none'
colors_sonia_f = 'none'
colors_elena_f = 'none'
cNorm  = colors.Normalize(vmin=0, vmax=1)

"""
markers = ['o','o','o','o','o','o','o','o']    
colormap = 'tab20' 
colors_cm = 'b'
colors_r = 'r'
colors_e = 'orange'
colors_s = 'purple'
colors_sandra = 'g'
colors_ruth = 'm'
colors_sonia = 'lime'
colors_elena = 'c'
colors_cm_f = 'b'
colors_r_f = 'r'
colors_e_f = 'orange'
colors_s_f = 'purple'
colors_sandra_f = 'g'
colors_ruth_f = 'm'
colors_sonia_f = 'lime'
colors_elena_f = 'c'
"""

"""
markers = ['s','d','^','p','+','x',(3,2,60),(5,2,72)]
colormap = 'tab20'
values = range(0,20)
cNorm  = colors.Normalize(vmin=0, vmax = values[-1]+1)
colors_cm = np.linspace(0,len(objs_pd_cm)-1,num = len(objs_pd_cm)) % 20
colors_r = np.linspace(0,len(objs_pd_r)-1,num = len(objs_pd_r)) % 20
colors_e = np.linspace(0,len(objs_pd_e)-1,num = len(objs_pd_e)) % 20
colors_s = np.linspace(0,len(objs_pd_s)-1,num = len(objs_pd_s)) % 20
colors_sandra = np.linspace(0,len(objs_pd_sandra)-1,num = len(objs_pd_sandra)) % 20
colors_ruth = np.linspace(0,len(objs_pd_ruth)-1,num = len(objs_pd_ruth)) % 20
colors_sonia = np.linspace(0,len(objs_pd_sonia)-1,num = len(objs_pd_sonia)) % 20
colors_elena = np.linspace(0,len(objs_pd_elena)-1,num = len(objs_pd_elena)) % 20
colors_cm_f = np.linspace(0,len(objs_pd_cm)-1,num = len(objs_pd_cm)) % 20
colors_r_f = np.linspace(0,len(objs_pd_r)-1,num = len(objs_pd_r)) % 20
colors_e_f = np.linspace(0,len(objs_pd_e)-1,num = len(objs_pd_e)) % 20
colors_s_f = np.linspace(0,len(objs_pd_s)-1,num = len(objs_pd_s)) % 20
colors_sandra_f = np.linspace(0,len(objs_pd_sandra)-1,num = len(objs_pd_sandra)) % 20
colors_ruth_f = np.linspace(0,len(objs_pd_ruth)-1,num = len(objs_pd_ruth)) % 20
colors_sonia_f = np.linspace(0,len(objs_pd_sonia)-1,num = len(objs_pd_sonia)) % 20
colors_elena_f = np.linspace(0,len(objs_pd_elena)-1,num = len(objs_pd_elena)) % 20
"""


#Cpt Marvel
tfile = prefix + 'cptmarvel.cosmo25cmb/cptmarvel.cosmo25cmb.4096g5HbwK1BH/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096'
objs_cm = []
f=open(tfile + '.data', 'rb')
while 1:
    try:
        objs_cm.append(pickle.load(f))
    except EOFError:
        break
        
f.close()
objs_pd_cm = pd.DataFrame(objs_cm)

#Elektra
tfile = prefix + 'elektra.cosmo25cmb/elektra.cosmo25cmb.4096g5HbwK1BH/elektra.cosmo25cmb.4096g5HbwK1BH.004096/elektra.cosmo25cmb.4096g5HbwK1BH.004096'
objs_e = []
f=open(tfile + '.data', 'rb')
while 1:
    try:
        objs_e.append(pickle.load(f))
    except EOFError:
        break
        
f.close()
objs_pd_e = pd.DataFrame(objs_e)

#Rogue
tfile = prefix + 'rogue.cosmo25cmb/rogue.cosmo25cmb.4096g5HbwK1BH/rogue.cosmo25cmb.4096g5HbwK1BH.004096/rogue.cosmo25cmb.4096g5HbwK1BH.004096'
objs_r = []
f=open(tfile + '.data', 'rb')
while 1:
    try:
        objs_r.append(pickle.load(f))
    except EOFError:
        break
        
f.close()
objs_pd_r = pd.DataFrame(objs_r)

#Storm
tfile = prefix + 'storm.cosmo25cmb/storm.cosmo25cmb.4096g5HbwK1BH/storm.cosmo25cmb.4096g5HbwK1BH.004096/storm.cosmo25cmb.4096g5HbwK1BH.004096'
objs_s = []
f=open(tfile + '.data', 'rb')
while 1:
    try:
        objs_s.append(pickle.load(f))
    except EOFError:
        break
        
f.close()
objs_pd_s = pd.DataFrame(objs_s)
#markers = ['s','d','^','p']

plt.figure(1,figsize=(plt_width,plt_width*aspect_ratio))
cptmarvel_plt = plt.scatter(np.array(objs_pd_cm['mstar']),np.array(objs_pd_cm['mHI'])/np.array(objs_pd_cm['mgas']),marker = markers[0], label = 'Cpt. Marvel' ,zorder = 2,c = colors_cm, cmap = colormap, norm=cNorm)
rogue_plt = plt.scatter(np.array(objs_pd_r['mstar']),np.array(objs_pd_r['mHI'])/np.array(objs_pd_r['mgas']),marker = markers[1], label = 'Rogue', zorder = 3,c = colors_r, cmap = colormap, norm=cNorm)
elektra_plt = plt.scatter(np.array(objs_pd_e['mstar']),np.array(objs_pd_e['mHI'])/np.array(objs_pd_e['mgas']),marker = markers[2], label = 'Elektra',zorder = 4,c = colors_e, cmap = colormap, norm=cNorm)
storm_plt = plt.scatter(np.array(objs_pd_s['mstar']),np.array(objs_pd_s['mHI'])/np.array(objs_pd_s['mgas']),marker = markers[3], label = 'Storm', zorder = 5,c = colors_s, cmap = colormap, norm=cNorm)
#sandra_plt = plt.scatter(np.array(objs_pd_sandra['mstar']),np.array(objs_pd_sandra['mHI'])/np.array(objs_pd_sandra['mgas']),marker = markers[4], label = 'Sandra', zorder = 5,c = colors_sandra, cmap = colormap, norm=cNorm)
#ruth_plt = plt.scatter(np.array(objs_pd_ruth['mstar']),np.array(objs_pd_ruth['mHI'])/np.array(objs_pd_ruth['mgas']),marker = markers[5], label = 'Ruth', zorder = 5,c = colors_ruth, cmap = colormap, norm=cNorm)
#elena_plt = plt.scatter(np.array(objs_pd_elena['mstar']),np.array(objs_pd_elena['mHI'])/np.array(objs_pd_elena['mgas']),marker = markers[6], label = 'Elena', zorder = 5,c = colors_elena, cmap = colormap, norm=cNorm)
plt.xlabel(r'M$_*$/M$_\odot$')
plt.ylabel(r'M$_{HI}$/M$_{gas}$')
plt.xscale('log')
#plt.axis([1e6, 3e9, 0, 8.5])
plt.axis([1e3, 1e10, 0, 0.3])
plt.tight_layout()
plt.show()
plt.savefig(outbase + '.HI_mgas_vMstar.png',dpi = dpi)
plt.clf()

plt.figure(1,figsize=(plt_width,plt_width*aspect_ratio))
cptmarvel_plt = plt.scatter(np.array(objs_pd_cm['SFR'])/np.array(objs_pd_cm['mstar']),np.array(objs_pd_cm['mHI'])/np.array(objs_pd_cm['mgas']),marker = markers[0], label = 'Cpt. Marvel' ,zorder = 2,c = colors_cm, cmap = colormap, norm=cNorm)
rogue_plt = plt.scatter(np.array(objs_pd_r['SFR'])/np.array(objs_pd_r['mstar']),np.array(objs_pd_r['mHI'])/np.array(objs_pd_r['mgas']),marker = markers[1], label = 'Rogue', zorder = 3,c = colors_r, cmap = colormap, norm=cNorm)
elektra_plt = plt.scatter(np.array(objs_pd_e['SFR'])/np.array(objs_pd_e['mstar']),np.array(objs_pd_e['mHI'])/np.array(objs_pd_e['mgas']),marker = markers[2], label = 'Elektra',zorder = 4,c = colors_e, cmap = colormap, norm=cNorm)
storm_plt = plt.scatter(np.array(objs_pd_s['SFR'])/np.array(objs_pd_s['mstar']),np.array(objs_pd_s['mHI'])/np.array(objs_pd_s['mgas']),marker = markers[3], label = 'Storm', zorder = 5,c = colors_s, cmap = colormap, norm=cNorm)
#sandra_plt = plt.scatter(np.array(objs_pd_sandra['mstar']),np.array(objs_pd_sandra['mHI'])/np.array(objs_pd_sandra['mgas']),marker = markers[4], label = 'Sandra', zorder = 5,c = colors_sandra, cmap = colormap, norm=cNorm)
#ruth_plt = plt.scatter(np.array(objs_pd_ruth['mstar']),np.array(objs_pd_ruth['mHI'])/np.array(objs_pd_ruth['mgas']),marker = markers[5], label = 'Ruth', zorder = 5,c = colors_ruth, cmap = colormap, norm=cNorm)
#elena_plt = plt.scatter(np.array(objs_pd_elena['mstar']),np.array(objs_pd_elena['mHI'])/np.array(objs_pd_elena['mgas']),marker = markers[6], label = 'Elena', zorder = 5,c = colors_elena, cmap = colormap, norm=cNorm)
plt.xlabel(r'SSFR [yr$^{-1}$]')
plt.ylabel(r'M$_{HI}$/M$_{gas}$')
plt.xscale('log')
#plt.axis([1e6, 3e9, 0, 8.5])
plt.axis([6e-12, 1e-9, 0, 0.3])
plt.tight_layout()
plt.show()
plt.savefig(outbase + '.HI_mgas_vSSFR.png',dpi = dpi)
plt.clf()

plt.figure(1,figsize=(plt_width,plt_width*aspect_ratio))
cptmarvel_plt = plt.scatter(np.array(objs_pd_cm['mvir']),np.array(objs_pd_cm['mHI'])/np.array(objs_pd_cm['mgas']),marker = markers[0], label = 'Cpt. Marvel' ,zorder = 2,c = colors_cm, cmap = colormap, norm=cNorm)
rogue_plt = plt.scatter(np.array(objs_pd_r['mvir']),np.array(objs_pd_r['mHI'])/np.array(objs_pd_r['mgas']),marker = markers[1], label = 'Rogue', zorder = 3,c = colors_r, cmap = colormap, norm=cNorm)
elektra_plt = plt.scatter(np.array(objs_pd_e['mvir']),np.array(objs_pd_e['mHI'])/np.array(objs_pd_e['mgas']),marker = markers[2], label = 'Elektra',zorder = 4,c = colors_e, cmap = colormap, norm=cNorm)
storm_plt = plt.scatter(np.array(objs_pd_s['mvir']),np.array(objs_pd_s['mHI'])/np.array(objs_pd_s['mgas']),marker = markers[3], label = 'Storm', zorder = 5,c = colors_s, cmap = colormap, norm=cNorm)
#sandra_plt = plt.scatter(np.array(objs_pd_sandra['mvir']),np.array(objs_pd_sandra['mHI'])/np.array(objs_pd_sandra['mgas']),marker = markers[4], label = 'Sandra', zorder = 5,c = colors_sandra, cmap = colormap, norm=cNorm)
#ruth_plt = plt.scatter(np.array(objs_pd_ruth['mvir']),np.array(objs_pd_ruth['mHI'])/np.array(objs_pd_ruth['mgas']),marker = markers[5], label = 'Ruth', zorder = 5,c = colors_ruth, cmap = colormap, norm=cNorm)
#elena_plt = plt.scatter(np.array(objs_pd_elena['mvir']),np.array(objs_pd_elena['mHI'])/np.array(objs_pd_elena['mgas']),marker = markers[6], label = 'Elena', zorder = 5,c = colors_elena, cmap = colormap, norm=cNorm)
plt.xlabel(r'M$_{vir}$/M$_\odot$')
plt.ylabel(r'M$_{HI}$/M$_{gas}$')
plt.xscale('log')
#plt.axis([1e6, 3e9, 0, 8.5])
plt.axis([8e8, 1e11, 0, 0.3])
plt.tight_layout()
plt.show()
plt.savefig(outbase + '.HI_mgas_vMvir.png',dpi = dpi)
plt.clf()

plt.figure(1,figsize=(plt_width,plt_width*aspect_ratio))
cptmarvel_plt = plt.scatter(np.array(objs_pd_cm['mstar']),np.array(objs_pd_cm['mHI'])/np.array(objs_pd_cm['mvir']),marker = markers[0], label = 'Cpt. Marvel' ,zorder = 2,c = colors_cm, cmap = colormap, norm=cNorm)
rogue_plt = plt.scatter(np.array(objs_pd_r['mstar']),np.array(objs_pd_r['mHI'])/np.array(objs_pd_r['mvir']),marker = markers[1], label = 'Rogue', zorder = 3,c = colors_r, cmap = colormap, norm=cNorm)
elektra_plt = plt.scatter(np.array(objs_pd_e['mstar']),np.array(objs_pd_e['mHI'])/np.array(objs_pd_e['mvir']),marker = markers[2], label = 'Elektra',zorder = 4,c = colors_e, cmap = colormap, norm=cNorm)
storm_plt = plt.scatter(np.array(objs_pd_s['mstar']),np.array(objs_pd_s['mHI'])/np.array(objs_pd_s['mvir']),marker = markers[3], label = 'Storm', zorder = 5,c = colors_s, cmap = colormap, norm=cNorm)
#sandra_plt = plt.scatter(np.array(objs_pd_sandra['mstar']),np.array(objs_pd_sandra['mHI'])/np.array(objs_pd_sandra['mvir']),marker = markers[4], label = 'Sandra', zorder = 5,c = colors_sandra, cmap = colormap, norm=cNorm)
#ruth_plt = plt.scatter(np.array(objs_pd_ruth['mstar']),np.array(objs_pd_ruth['mHI'])/np.array(objs_pd_ruth['mvir']),marker = markers[5], label = 'Ruth', zorder = 5,c = colors_ruth, cmap = colormap, norm=cNorm)
#elena_plt = plt.scatter(np.array(objs_pd_elena['mstar']),np.array(objs_pd_elena['mHI'])/np.array(objs_pd_elena['mvir']),marker = markers[6], label = 'Elena', zorder = 5,c = colors_elena, cmap = colormap, norm=cNorm)
plt.xlabel(r'M$_*$/M$_\odot$')
plt.ylabel(r'M$_{HI}$/M$_{vir}$')
plt.xscale('log')
#plt.axis([1e6, 3e9, 0, 8.5])
plt.axis([1e3, 1e10, 0, 0.025])
plt.tight_layout()
plt.show()
plt.savefig(outbase + '.HI_mvir_vMstar.png',dpi = dpi)
plt.clf()

plt.figure(1,figsize=(plt_width,plt_width*aspect_ratio))
cptmarvel_plt = plt.scatter(np.array(objs_pd_cm['SFR'])/np.array(objs_pd_cm['mstar']),np.array(objs_pd_cm['mHI'])/np.array(objs_pd_cm['mvir']),marker = markers[0], label = 'Cpt. Marvel' ,zorder = 2,c = colors_cm, cmap = colormap, norm=cNorm)
rogue_plt = plt.scatter(np.array(objs_pd_r['SFR'])/np.array(objs_pd_r['mstar']),np.array(objs_pd_r['mHI'])/np.array(objs_pd_r['mvir']),marker = markers[1], label = 'Rogue', zorder = 3,c = colors_r, cmap = colormap, norm=cNorm)
elektra_plt = plt.scatter(np.array(objs_pd_e['SFR'])/np.array(objs_pd_e['mstar']),np.array(objs_pd_e['mHI'])/np.array(objs_pd_e['mvir']),marker = markers[2], label = 'Elektra',zorder = 4,c = colors_e, cmap = colormap, norm=cNorm)
storm_plt = plt.scatter(np.array(objs_pd_s['SFR'])/np.array(objs_pd_s['mstar']),np.array(objs_pd_s['mHI'])/np.array(objs_pd_s['mvir']),marker = markers[3], label = 'Storm', zorder = 5,c = colors_s, cmap = colormap, norm=cNorm)
#sandra_plt = plt.scatter(np.array(objs_pd_sandra['mstar']),np.array(objs_pd_sandra['mHI'])/np.array(objs_pd_sandra['mvir']),marker = markers[4], label = 'Sandra', zorder = 5,c = colors_sandra, cmap = colormap, norm=cNorm)
#ruth_plt = plt.scatter(np.array(objs_pd_ruth['mstar']),np.array(objs_pd_ruth['mHI'])/np.array(objs_pd_ruth['mvir']),marker = markers[5], label = 'Ruth', zorder = 5,c = colors_ruth, cmap = colormap, norm=cNorm)
#elena_plt = plt.scatter(np.array(objs_pd_elena['mstar']),np.array(objs_pd_elena['mHI'])/np.array(objs_pd_elena['mvir']),marker = markers[6], label = 'Elena', zorder = 5,c = colors_elena, cmap = colormap, norm=cNorm)
plt.xlabel(r'SSFR [yr$^{-1}$]')
plt.ylabel(r'M$_{HI}$/M$_{vir}$')
plt.xscale('log')
#plt.axis([1e6, 3e9, 0, 8.5])
plt.axis([6e-12, 1e-9, 0, 0.025])
plt.tight_layout()
plt.show()
plt.savefig(outbase + '.HI_mvir_vSSFR.png',dpi = dpi)
plt.clf()

plt.figure(1,figsize=(plt_width,plt_width*aspect_ratio))
cptmarvel_plt = plt.scatter(np.array(objs_pd_cm['mvir']),np.array(objs_pd_cm['mHI'])/np.array(objs_pd_cm['mvir']),marker = markers[0], label = 'Cpt. Marvel' ,zorder = 2,c = colors_cm, cmap = colormap, norm=cNorm)
rogue_plt = plt.scatter(np.array(objs_pd_r['mvir']),np.array(objs_pd_r['mHI'])/np.array(objs_pd_r['mvir']),marker = markers[1], label = 'Rogue', zorder = 3,c = colors_r, cmap = colormap, norm=cNorm)
elektra_plt = plt.scatter(np.array(objs_pd_e['mvir']),np.array(objs_pd_e['mHI'])/np.array(objs_pd_e['mvir']),marker = markers[2], label = 'Elektra',zorder = 4,c = colors_e, cmap = colormap, norm=cNorm)
storm_plt = plt.scatter(np.array(objs_pd_s['mvir']),np.array(objs_pd_s['mHI'])/np.array(objs_pd_s['mvir']),marker = markers[3], label = 'Storm', zorder = 5,c = colors_s, cmap = colormap, norm=cNorm)
#sandra_plt = plt.scatter(np.array(objs_pd_sandra['mvir']),np.array(objs_pd_sandra['mHI'])/np.array(objs_pd_sandra['mvir']),marker = markers[4], label = 'Sandra', zorder = 5,c = colors_sandra, cmap = colormap, norm=cNorm)
#ruth_plt = plt.scatter(np.array(objs_pd_ruth['mvir']),np.array(objs_pd_ruth['mHI'])/np.array(objs_pd_ruth['mvir']),marker = markers[5], label = 'Ruth', zorder = 5,c = colors_ruth, cmap = colormap, norm=cNorm)
#elena_plt = plt.scatter(np.array(objs_pd_elena['mvir']),np.array(objs_pd_elena['mHI'])/np.array(objs_pd_elena['mvir']),marker = markers[6], label = 'Elena', zorder = 5,c = colors_elena, cmap = colormap, norm=cNorm)
plt.xlabel(r'M$_{vir}$/M$_\odot$')
plt.ylabel(r'M$_{HI}$/M$_{vir}$')
plt.xscale('log')
#plt.axis([1e6, 3e9, 0, 8.5])
plt.axis([8e8, 1e11, 0, 0.025])
plt.tight_layout()
plt.show()
plt.savefig(outbase + '.HI_mvir_vMvir.png',dpi = dpi)
plt.clf()

plt.figure(1,figsize=(plt_width,plt_width*aspect_ratio))
cptmarvel_plt = plt.scatter(np.array(objs_pd_cm['mstar']),np.array(objs_pd_cm['mgas'])/np.array(objs_pd_cm['mvir']),marker = markers[0], label = 'Cpt. Marvel' ,zorder = 2,c = colors_cm, cmap = colormap, norm=cNorm)
rogue_plt = plt.scatter(np.array(objs_pd_r['mstar']),np.array(objs_pd_r['mgas'])/np.array(objs_pd_r['mvir']),marker = markers[1], label = 'Rogue', zorder = 3,c = colors_r, cmap = colormap, norm=cNorm)
elektra_plt = plt.scatter(np.array(objs_pd_e['mstar']),np.array(objs_pd_e['mgas'])/np.array(objs_pd_e['mvir']),marker = markers[2], label = 'Elektra',zorder = 4,c = colors_e, cmap = colormap, norm=cNorm)
storm_plt = plt.scatter(np.array(objs_pd_s['mstar']),np.array(objs_pd_s['mgas'])/np.array(objs_pd_s['mvir']),marker = markers[3], label = 'Storm', zorder = 5,c = colors_s, cmap = colormap, norm=cNorm)
plt.xlabel(r'M$_*$/M$_\odot$')
plt.ylabel(r'M$_{gas}$/M$_{vir}$')
plt.xscale('log')
plt.axis([1e3, 1e10, 0, 0.1])
plt.tight_layout()
plt.show()
plt.savefig(outbase + '_mgas_mvir_vMstar.png',dpi = dpi)
plt.clf()

plt.figure(1,figsize=(plt_width,plt_width*aspect_ratio))
cptmarvel_plt = plt.scatter(np.array(objs_pd_cm['SFR'])/np.array(objs_pd_cm['mstar']),np.array(objs_pd_cm['mgas'])/np.array(objs_pd_cm['mvir']),marker = markers[0], label = 'Cpt. Marvel' ,zorder = 2,c = colors_cm, cmap = colormap, norm=cNorm)
rogue_plt = plt.scatter(np.array(objs_pd_r['SFR'])/np.array(objs_pd_r['mstar']),np.array(objs_pd_r['mgas'])/np.array(objs_pd_r['mvir']),marker = markers[1], label = 'Rogue', zorder = 3,c = colors_r, cmap = colormap, norm=cNorm)
elektra_plt = plt.scatter(np.array(objs_pd_e['SFR'])/np.array(objs_pd_e['mstar']),np.array(objs_pd_e['mgas'])/np.array(objs_pd_e['mvir']),marker = markers[2], label = 'Elektra',zorder = 4,c = colors_e, cmap = colormap, norm=cNorm)
storm_plt = plt.scatter(np.array(objs_pd_s['SFR'])/np.array(objs_pd_s['mstar']),np.array(objs_pd_s['mgas'])/np.array(objs_pd_s['mvir']),marker = markers[3], label = 'Storm', zorder = 5,c = colors_s, cmap = colormap, norm=cNorm)
plt.xlabel(r'SSFR [yr$^{-1}$]')
plt.ylabel(r'M$_{gas}$/M$_{vir}$')
plt.xscale('log')
plt.axis([6e-12, 1e-9, 0, 0.1])
plt.tight_layout()
plt.show()
plt.savefig(outbase + '_mgas_mvir_vSSFR.png',dpi = dpi)
plt.clf()

plt.figure(1,figsize=(plt_width,plt_width*aspect_ratio))
cptmarvel_plt = plt.scatter(np.array(objs_pd_cm['mvir']),np.array(objs_pd_cm['mgas'])/np.array(objs_pd_cm['mvir']),marker = markers[0], label = 'Cpt. Marvel' ,zorder = 2,c = colors_cm, cmap = colormap, norm=cNorm)
rogue_plt = plt.scatter(np.array(objs_pd_r['mvir']),np.array(objs_pd_r['mgas'])/np.array(objs_pd_r['mvir']),marker = markers[1], label = 'Rogue', zorder = 3,c = colors_r, cmap = colormap, norm=cNorm)
elektra_plt = plt.scatter(np.array(objs_pd_e['mvir']),np.array(objs_pd_e['mgas'])/np.array(objs_pd_e['mvir']),marker = markers[2], label = 'Elektra',zorder = 4,c = colors_e, cmap = colormap, norm=cNorm)
storm_plt = plt.scatter(np.array(objs_pd_s['mvir']),np.array(objs_pd_s['mgas'])/np.array(objs_pd_s['mvir']),marker = markers[3], label = 'Storm', zorder = 5,c = colors_s, cmap = colormap, norm=cNorm)
plt.xlabel(r'M$_{vir}$/M$_\odot$')
plt.ylabel(r'M$_{gas}$/M$_{vir}$')
plt.xscale('log')
plt.axis([8e8, 1e11, 0, 0.1])
plt.tight_layout()
plt.show()
plt.savefig(outbase + '_mgas_mvir_vMvir.png',dpi = dpi)
plt.clf()


plt.figure(1,figsize=(plt_width,plt_width*aspect_ratio))
cptmarvel_plt = plt.scatter(objs_pd_cm['mstar'],(objs_pd_cm['mgas'] - objs_pd_cm['mHI'])/objs_pd_cm['mvir'],marker = markers[0], label = 'Cpt. Marvel' ,zorder = 2,c = colors_cm, cmap = colormap, norm=cNorm)
rogue_plt = plt.scatter(objs_pd_r['mstar'],(objs_pd_r['mgas'] - objs_pd_r['mHI'])/objs_pd_r['mvir'],marker = markers[1], label = 'Rogue' ,zorder = 2,c = colors_r, cmap = colormap, norm=cNorm)
elektra_plt = plt.scatter(objs_pd_e['mstar'],(objs_pd_e['mgas'] - objs_pd_e['mHI'])/objs_pd_e['mvir'],marker = markers[2], label = 'Elektra' ,zorder = 2,c = colors_e, cmap = colormap, norm=cNorm)
storm_plt = plt.scatter(objs_pd_s['mstar'],(objs_pd_s['mgas'] - objs_pd_s['mHI'])/objs_pd_s['mvir'],marker = markers[3], label = 'Storm' ,zorder = 2,c = colors_s, cmap = colormap, norm=cNorm)
plt.xlabel(r'M$_*$/M$_\odot$')
plt.ylabel(r'M$_{halo}$/M$_{vir}$')
plt.xscale('log')
plt.axis([1e3, 1e10, 0, 0.1])
plt.tight_layout()
plt.show()
plt.savefig(outbase + '_mhgas_mvir_vMstar.png',dpi = dpi)
plt.clf()

plt.figure(1,figsize=(plt_width,plt_width*aspect_ratio))
cptmarvel_plt = plt.scatter(objs_pd_cm['SFR']/objs_pd_cm['mstar'],(objs_pd_cm['mgas'] - objs_pd_cm['mHI'])/objs_pd_cm['mvir'],marker = markers[0], label = 'Cpt. Marvel' ,zorder = 2,c = colors_cm, cmap = colormap, norm=cNorm)
rogue_plt = plt.scatter(objs_pd_r['SFR']/objs_pd_r['mstar'],(objs_pd_r['mgas'] - objs_pd_r['mHI'])/objs_pd_r['mvir'],marker = markers[1], label = 'Rogue' ,zorder = 2,c = colors_r, cmap = colormap, norm=cNorm)
elektra_plt = plt.scatter(objs_pd_e['SFR']/objs_pd_e['mstar'],(objs_pd_e['mgas'] - objs_pd_e['mHI'])/objs_pd_e['mvir'],marker = markers[2], label = 'Elektra' ,zorder = 2,c = colors_e, cmap = colormap, norm=cNorm)
storm_plt = plt.scatter(objs_pd_s['SFR']/objs_pd_s['mstar'],(objs_pd_s['mgas'] - objs_pd_s['mHI'])/objs_pd_s['mvir'],marker = markers[3], label = 'Storm' ,zorder = 2,c = colors_s, cmap = colormap, norm=cNorm)
plt.xlabel(r'SSFR [yr$^{-1}$]')
plt.ylabel(r'M$_{halo}$/M$_{vir}$')
plt.xscale('log')
plt.axis([6e-12, 1e-9, 0, 0.1])
plt.tight_layout()
plt.show()
plt.savefig(outbase + '_mhgas_mvir_vSSFR.png',dpi = dpi)
plt.clf()

plt.figure(1,figsize=(plt_width,plt_width*aspect_ratio))
cptmarvel_plt = plt.scatter(objs_pd_cm['mvir'],(objs_pd_cm['mgas'] - objs_pd_cm['mHI'])/objs_pd_cm['mvir'],marker = markers[0], label = 'Cpt. Marvel' ,zorder = 2,c = colors_cm, cmap = colormap, norm=cNorm)
rogue_plt = plt.scatter(objs_pd_r['mvir'],(objs_pd_r['mgas'] - objs_pd_r['mHI'])/objs_pd_r['mvir'],marker = markers[1], label = 'Rogue' ,zorder = 2,c = colors_r, cmap = colormap, norm=cNorm)
elektra_plt = plt.scatter(objs_pd_e['mvir'],(objs_pd_e['mgas'] - objs_pd_e['mHI'])/objs_pd_e['mvir'],marker = markers[2], label = 'Elektra' ,zorder = 2,c = colors_e, cmap = colormap, norm=cNorm)
storm_plt = plt.scatter(objs_pd_s['mvir'],(objs_pd_s['mgas'] - objs_pd_s['mHI'])/objs_pd_s['mvir'],marker = markers[3], label = 'Storm' ,zorder = 2,c = colors_s, cmap = colormap, norm=cNorm)
plt.xlabel(r'M$_{vir}$/M$_\odot$')
plt.ylabel(r'M$_{gas}$/M$_{vir}$')
plt.xscale('log')
plt.axis([8e8, 1e11, 0, 0.1])
plt.tight_layout()
plt.show()
plt.savefig(outbase + '_mhgas_mvir_vMvir.png',dpi = dpi)
plt.clf()

######################################33

cm_mstar_log = np.empty(len(objs_pd_cm['mstar']))
for i in np.arange(0,len(objs_pd_cm['mstar'])):
    cm_mstar_log[i] = np.log10(objs_pd_cm['mstar'][i])
e_mstar_log = np.empty(len(objs_pd_e['mstar']))
for i in np.arange(0,len(objs_pd_e['mstar'])):
    e_mstar_log[i] = np.log10(objs_pd_e['mstar'][i])
r_mstar_log = np.empty(len(objs_pd_r['mstar']))
for i in np.arange(0,len(objs_pd_r['mstar'])):
    r_mstar_log[i] = np.log10(objs_pd_r['mstar'][i])
s_mstar_log = np.empty(len(objs_pd_s['mstar']))
for i in np.arange(0,len(objs_pd_s['mstar'])):
    s_mstar_log[i] = np.log10(objs_pd_s['mstar'][i])

cm_mvir_log = np.empty(len(objs_pd_cm['mvir']))
for i in np.arange(0,len(objs_pd_cm['mvir'])):
    cm_mvir_log[i] = np.log10(objs_pd_cm['mvir'][i])
e_mvir_log = np.empty(len(objs_pd_e['mvir']))
for i in np.arange(0,len(objs_pd_e['mvir'])):
    e_mvir_log[i] = np.log10(objs_pd_e['mvir'][i])
r_mvir_log = np.empty(len(objs_pd_r['mvir']))
for i in np.arange(0,len(objs_pd_r['mvir'])):
    r_mvir_log[i] = np.log10(objs_pd_r['mvir'][i])
s_mvir_log = np.empty(len(objs_pd_s['mvir']))
for i in np.arange(0,len(objs_pd_s['mvir'])):
    s_mvir_log[i] = np.log10(objs_pd_s['mvir'][i])

cm_SSFR = objs_pd_cm['SFR']/objs_pd_cm['mstar']
e_SSFR = objs_pd_e['SFR']/objs_pd_e['mstar']
r_SSFR = objs_pd_r['SFR']/objs_pd_r['mstar']
s_SSFR = objs_pd_s['SFR']/objs_pd_s['mstar']

#--------------------------
plt.figure(1,figsize=(plt_width,plt_width*aspect_ratio))
hot_plt = plt.bar(cm_mstar_log,(objs_pd_cm['mgas'] - (objs_pd_cm['mwarm'] + objs_pd_cm['mCool'] + objs_pd_cm['mISM']))/objs_pd_cm['mvir']/f_bar,0.1,(objs_pd_cm['mwarm'] + objs_pd_cm['mCool'] + objs_pd_cm['mISM'] + objs_pd_cm['mstar'])/objs_pd_cm['mvir']/f_bar,color = 'goldenrod',alpha = 0.5)
plt.bar(e_mstar_log,(objs_pd_e['mgas'] - (objs_pd_e['mwarm'] + objs_pd_e['mCool'] + objs_pd_e['mISM']))/objs_pd_e['mvir']/f_bar,0.1,(objs_pd_e['mwarm'] + objs_pd_e['mCool'] + objs_pd_e['mISM'] + objs_pd_e['mstar'])/objs_pd_e['mvir']/f_bar,color = 'goldenrod',alpha = 0.5)
plt.bar(r_mstar_log,(objs_pd_r['mgas'] - (objs_pd_r['mwarm'] + objs_pd_r['mCool'] + objs_pd_r['mISM']))/objs_pd_r['mvir']/f_bar,0.1,(objs_pd_r['mwarm'] + objs_pd_r['mCool'] + objs_pd_r['mISM'] + objs_pd_r['mstar'])/objs_pd_r['mvir']/f_bar,color = 'goldenrod',alpha = 0.5)
plt.bar(s_mstar_log,(objs_pd_s['mgas'] - (objs_pd_s['mwarm'] + objs_pd_s['mCool'] + objs_pd_s['mISM']))/objs_pd_s['mvir']/f_bar,0.1,(objs_pd_s['mwarm'] + objs_pd_s['mCool'] + objs_pd_s['mISM']+ objs_pd_s['mstar'])/objs_pd_s['mvir']/f_bar,color = 'goldenrod',alpha = 0.5)

warm_plt =plt.bar(cm_mstar_log,objs_pd_cm['mwarm']/objs_pd_cm['mvir']/f_bar,0.1,(objs_pd_cm['mCool'] + objs_pd_cm['mISM'] + objs_pd_cm['mstar'])/objs_pd_cm['mvir']/f_bar,color= 'g',alpha = 0.5)
plt.bar(e_mstar_log,objs_pd_e['mwarm']/objs_pd_e['mvir']/f_bar,0.1,(objs_pd_e['mCool'] + objs_pd_e['mISM'] + objs_pd_e['mstar'])/objs_pd_e['mvir']/f_bar,color= 'g',alpha = 0.5)
plt.bar(r_mstar_log,objs_pd_r['mwarm']/objs_pd_r['mvir']/f_bar,0.1,(objs_pd_r['mCool'] + objs_pd_r['mISM'] + objs_pd_r['mstar'])/objs_pd_r['mvir']/f_bar,color= 'g',alpha = 0.5)
plt.bar(s_mstar_log,objs_pd_s['mwarm']/objs_pd_s['mvir']/f_bar,0.1,(objs_pd_s['mCool'] + objs_pd_s['mISM'] + objs_pd_s['mstar'])/objs_pd_s['mvir']/f_bar,color= 'g',alpha = 0.5)

cool_plt =plt.bar(cm_mstar_log,objs_pd_cm['mCool']/objs_pd_cm['mvir']/f_bar,0.1,(objs_pd_cm['mISM'] + objs_pd_cm['mstar'])/objs_pd_cm['mvir']/f_bar,color= 'darkviolet',alpha = 0.5)
plt.bar(e_mstar_log,objs_pd_e['mCool']/objs_pd_e['mvir']/f_bar,0.1,(objs_pd_e['mISM'] + objs_pd_e['mstar'])/objs_pd_e['mvir']/f_bar,color= 'darkviolet',alpha = 0.5)
plt.bar(r_mstar_log,objs_pd_r['mCool']/objs_pd_r['mvir']/f_bar,0.1,(objs_pd_r['mISM'] + objs_pd_r['mstar'])/objs_pd_r['mvir']/f_bar,color= 'darkviolet',alpha = 0.5)
plt.bar(s_mstar_log,objs_pd_s['mCool']/objs_pd_s['mvir']/f_bar,0.1,(objs_pd_s['mISM'] + objs_pd_s['mstar'])/objs_pd_s['mvir']/f_bar,color= 'darkviolet',alpha = 0.5)

ISM_plt =plt.bar(cm_mstar_log,objs_pd_cm['mISM']/objs_pd_cm['mvir']/f_bar,0.1,objs_pd_cm['mstar']/objs_pd_cm['mvir']/f_bar,color= 'b',alpha = 0.5)
plt.bar(e_mstar_log,objs_pd_e['mISM']/objs_pd_e['mvir']/f_bar,0.1,objs_pd_e['mstar']/objs_pd_e['mvir']/f_bar,color= 'b',alpha = 0.5)
plt.bar(r_mstar_log,objs_pd_r['mISM']/objs_pd_r['mvir']/f_bar,0.1,objs_pd_r['mstar']/objs_pd_r['mvir']/f_bar,color= 'b',alpha = 0.5)
plt.bar(s_mstar_log,objs_pd_s['mISM']/objs_pd_s['mvir']/f_bar,0.1,objs_pd_s['mstar']/objs_pd_s['mvir']/f_bar,color= 'b',alpha = 0.5)

stars_plt =plt.bar(cm_mstar_log,objs_pd_cm['mstar']/objs_pd_cm['mvir']/f_bar,width=0.1,color= 'r',alpha = 0.5)
plt.bar(e_mstar_log,objs_pd_e['mstar']/objs_pd_e['mvir']/f_bar,width=0.1,color= 'r',alpha = 0.5)
plt.bar(r_mstar_log,objs_pd_r['mstar']/objs_pd_r['mvir']/f_bar,width=0.1,color= 'r',alpha = 0.5)
plt.bar(s_mstar_log,objs_pd_s['mstar']/objs_pd_s['mvir']/f_bar,width=0.1,color= 'r',alpha = 0.5)
plt.xlabel(r'Log(M$_*$/M$_\odot)$')
plt.ylabel(r'$M/(f_{bary} \times M_{vir})$')
plt.legend([hot_plt,warm_plt,cool_plt,ISM_plt,stars_plt],['$10^6 K < T$','$10^5 K < T < 10^6 K$','$10^4 K < T < 10^5 K$','$T < 10^4 K$','Stars'],loc = 2,framealpha=0.5,fontsize = legendsize)
plt.tight_layout()
plt.show()
plt.savefig(outbase + '_gascensus_vMstar.png',dpi = dpi)
plt.close()

#--------------------------------------
plt.figure(1,figsize=(plt_width,plt_width*aspect_ratio))
hot_plt = plt.bar(cm_mvir_log,(objs_pd_cm['mgas'] - (objs_pd_cm['mwarm'] + objs_pd_cm['mCool'] + objs_pd_cm['mISM']))/objs_pd_cm['mvir']/f_bar,0.05,(objs_pd_cm['mwarm'] + objs_pd_cm['mCool'] + objs_pd_cm['mISM'] + objs_pd_cm['mstar'])/objs_pd_cm['mvir']/f_bar,color = 'goldenrod',alpha = 0.5)
plt.bar(e_mvir_log,(objs_pd_e['mgas'] - (objs_pd_e['mwarm'] + objs_pd_e['mCool'] + objs_pd_e['mISM']))/objs_pd_e['mvir']/f_bar,0.05,(objs_pd_e['mwarm'] + objs_pd_e['mCool'] + objs_pd_e['mISM'] + objs_pd_e['mstar'])/objs_pd_e['mvir']/f_bar,color = 'goldenrod',alpha = 0.5)
plt.bar(r_mvir_log,(objs_pd_r['mgas'] - (objs_pd_r['mwarm'] + objs_pd_r['mCool'] + objs_pd_r['mISM']))/objs_pd_r['mvir']/f_bar,0.05,(objs_pd_r['mwarm'] + objs_pd_r['mCool'] + objs_pd_r['mISM'] + objs_pd_r['mstar'])/objs_pd_r['mvir']/f_bar,color = 'goldenrod',alpha = 0.5)
plt.bar(s_mvir_log,(objs_pd_s['mgas'] - (objs_pd_s['mwarm'] + objs_pd_s['mCool'] + objs_pd_s['mISM']))/objs_pd_s['mvir']/f_bar,0.05,(objs_pd_s['mwarm'] + objs_pd_s['mCool'] + objs_pd_s['mISM']+ objs_pd_s['mstar'])/objs_pd_s['mvir']/f_bar,color = 'goldenrod',alpha = 0.5)

warm_plt = plt.bar(cm_mvir_log,objs_pd_cm['mwarm']/objs_pd_cm['mvir']/f_bar,0.05,(objs_pd_cm['mCool'] + objs_pd_cm['mISM'] + objs_pd_cm['mstar'])/objs_pd_cm['mvir']/f_bar,color= 'g',alpha = 0.5)
plt.bar(e_mvir_log,objs_pd_e['mwarm']/objs_pd_e['mvir']/f_bar,0.05,(objs_pd_e['mCool'] + objs_pd_e['mISM'] + objs_pd_e['mstar'])/objs_pd_e['mvir']/f_bar,color= 'g',alpha = 0.5)
plt.bar(r_mvir_log,objs_pd_r['mwarm']/objs_pd_r['mvir']/f_bar,0.05,(objs_pd_r['mCool'] + objs_pd_r['mISM'] + objs_pd_r['mstar'])/objs_pd_r['mvir']/f_bar,color= 'g',alpha = 0.5)
plt.bar(s_mvir_log,objs_pd_s['mwarm']/objs_pd_s['mvir']/f_bar,0.05,(objs_pd_s['mCool'] + objs_pd_s['mISM'] + objs_pd_s['mstar'])/objs_pd_s['mvir']/f_bar,color= 'g',alpha = 0.5)

cool_plt = plt.bar(cm_mvir_log,objs_pd_cm['mCool']/objs_pd_cm['mvir']/f_bar,0.05,(objs_pd_cm['mISM'] + objs_pd_cm['mstar'])/objs_pd_cm['mvir']/f_bar,color= 'darkviolet',alpha = 0.5)
plt.bar(e_mvir_log,objs_pd_e['mCool']/objs_pd_e['mvir']/f_bar,0.05,(objs_pd_e['mISM'] + objs_pd_e['mstar'])/objs_pd_e['mvir']/f_bar,color= 'darkviolet',alpha = 0.5)
plt.bar(r_mvir_log,objs_pd_r['mCool']/objs_pd_r['mvir']/f_bar,0.05,(objs_pd_r['mISM'] + objs_pd_r['mstar'])/objs_pd_r['mvir']/f_bar,color= 'darkviolet',alpha = 0.5)
plt.bar(s_mvir_log,objs_pd_s['mCool']/objs_pd_s['mvir']/f_bar,0.05,(objs_pd_s['mISM'] + objs_pd_s['mstar'])/objs_pd_s['mvir']/f_bar,color= 'darkviolet',alpha = 0.5)

ISM_plt = plt.bar(cm_mvir_log,objs_pd_cm['mISM']/objs_pd_cm['mvir']/f_bar,0.05,objs_pd_cm['mstar']/objs_pd_cm['mvir']/f_bar,color= 'b',alpha = 0.5)
plt.bar(e_mvir_log,objs_pd_e['mISM']/objs_pd_e['mvir']/f_bar,0.05,objs_pd_e['mstar']/objs_pd_e['mvir']/f_bar,color= 'b',alpha = 0.5)
plt.bar(r_mvir_log,objs_pd_r['mISM']/objs_pd_r['mvir']/f_bar,0.05,objs_pd_r['mstar']/objs_pd_r['mvir']/f_bar,color= 'b',alpha = 0.5)
plt.bar(s_mvir_log,objs_pd_s['mISM']/objs_pd_s['mvir']/f_bar,0.05,objs_pd_s['mstar']/objs_pd_s['mvir']/f_bar,color= 'b',alpha = 0.5)

stars_plt = plt.bar(cm_mvir_log,objs_pd_cm['mstar']/objs_pd_cm['mvir']/f_bar,width=0.05,color= 'r',alpha = 0.5)
plt.bar(e_mvir_log,objs_pd_e['mstar']/objs_pd_e['mvir']/f_bar,width=0.05,color= 'r',alpha = 0.5)
plt.bar(r_mvir_log,objs_pd_r['mstar']/objs_pd_r['mvir']/f_bar,width=0.05,color= 'r',alpha = 0.5)
plt.bar(s_mvir_log,objs_pd_s['mstar']/objs_pd_s['mvir']/f_bar,width=0.05,color= 'r',alpha = 0.5)
plt.xlabel(r'Log(M$_{vir}$/M$_\odot)$')
plt.ylabel(r'$M/(f_{bary} \times M_{vir})$')
#plt.yscale('log')
#plt.axis([8, 11.2, 1e-2, 1])
plt.legend([hot_plt,warm_plt,cool_plt,ISM_plt,stars_plt],['$10^6 K < T$','$10^5 K < T < 10^6 K$','$10^4 K < T < 10^5 K$','$T < 10^4 K$','Stars'],loc = 2,framealpha=0.5,fontsize = legendsize)
plt.tight_layout()
plt.show()
plt.savefig(outbase + '_gascensus_vMvir.png',dpi = dpi)
plt.close()

#-----------------------------------
plt.figure(1,figsize=(plt_width,plt_width*aspect_ratio))
hot_plt = plt.bar(cm_SSFR,(objs_pd_cm['mgas'] - (objs_pd_cm['mwarm'] + objs_pd_cm['mCool'] + objs_pd_cm['mISM']))/objs_pd_cm['mvir']/f_bar,1e-11,(objs_pd_cm['mwarm'] + objs_pd_cm['mCool'] + objs_pd_cm['mISM'] + objs_pd_cm['mstar'])/objs_pd_cm['mvir']/f_bar,color = 'goldenrod',alpha = 0.5)
plt.bar(e_SSFR,(objs_pd_e['mgas'] - (objs_pd_e['mwarm'] + objs_pd_e['mCool'] + objs_pd_e['mISM']))/objs_pd_e['mvir']/f_bar,1e-11,(objs_pd_e['mwarm'] + objs_pd_e['mCool'] + objs_pd_e['mISM'] + objs_pd_e['mstar'])/objs_pd_e['mvir']/f_bar,color = 'goldenrod',alpha = 0.5)
plt.bar(r_SSFR,(objs_pd_r['mgas'] - (objs_pd_r['mwarm'] + objs_pd_r['mCool'] + objs_pd_r['mISM']))/objs_pd_r['mvir']/f_bar,1e-11,(objs_pd_r['mwarm'] + objs_pd_r['mCool'] + objs_pd_r['mISM'] + objs_pd_r['mstar'])/objs_pd_r['mvir']/f_bar,color = 'goldenrod',alpha = 0.5)
plt.bar(s_SSFR,(objs_pd_s['mgas'] - (objs_pd_s['mwarm'] + objs_pd_s['mCool'] + objs_pd_s['mISM']))/objs_pd_s['mvir']/f_bar,1e-11,(objs_pd_s['mwarm'] + objs_pd_s['mCool'] + objs_pd_s['mISM']+ objs_pd_s['mstar'])/objs_pd_s['mvir']/f_bar,color = 'goldenrod',alpha = 0.5)

warm_plt = plt.bar(cm_SSFR,objs_pd_cm['mwarm']/objs_pd_cm['mvir']/f_bar,1e-11,(objs_pd_cm['mCool'] + objs_pd_cm['mISM'] + objs_pd_cm['mstar'])/objs_pd_cm['mvir']/f_bar,color= 'g',alpha = 0.5)
plt.bar(e_SSFR,objs_pd_e['mwarm']/objs_pd_e['mvir']/f_bar,1e-11,(objs_pd_e['mCool'] + objs_pd_e['mISM'] + objs_pd_e['mstar'])/objs_pd_e['mvir']/f_bar,color= 'g',alpha = 0.5)
plt.bar(r_SSFR,objs_pd_r['mwarm']/objs_pd_r['mvir']/f_bar,1e-11,(objs_pd_r['mCool'] + objs_pd_r['mISM'] + objs_pd_r['mstar'])/objs_pd_r['mvir']/f_bar,color= 'g',alpha = 0.5)
plt.bar(s_SSFR,objs_pd_s['mwarm']/objs_pd_s['mvir']/f_bar,1e-11,(objs_pd_s['mCool'] + objs_pd_s['mISM'] + objs_pd_s['mstar'])/objs_pd_s['mvir']/f_bar,color= 'g',alpha = 0.5)

cool_plt = plt.bar(cm_SSFR,objs_pd_cm['mCool']/objs_pd_cm['mvir']/f_bar,1e-11,(objs_pd_cm['mISM'] + objs_pd_cm['mstar'])/objs_pd_cm['mvir']/f_bar,color= 'darkviolet',alpha = 0.5)
plt.bar(e_SSFR,objs_pd_e['mCool']/objs_pd_e['mvir']/f_bar,1e-11,(objs_pd_e['mISM'] + objs_pd_e['mstar'])/objs_pd_e['mvir']/f_bar,color= 'darkviolet',alpha = 0.5)
plt.bar(r_SSFR,objs_pd_r['mCool']/objs_pd_r['mvir']/f_bar,1e-11,(objs_pd_r['mISM'] + objs_pd_r['mstar'])/objs_pd_r['mvir']/f_bar,color= 'darkviolet',alpha = 0.5)
plt.bar(s_SSFR,objs_pd_s['mCool']/objs_pd_s['mvir']/f_bar,1e-11,(objs_pd_s['mISM'] + objs_pd_s['mstar'])/objs_pd_s['mvir']/f_bar,color= 'darkviolet',alpha = 0.5)

ISM_plt = plt.bar(cm_SSFR,objs_pd_cm['mISM']/objs_pd_cm['mvir']/f_bar,1e-11,objs_pd_cm['mstar']/objs_pd_cm['mvir']/f_bar,color= 'b',alpha = 0.5)
plt.bar(e_SSFR,objs_pd_e['mISM']/objs_pd_e['mvir']/f_bar,1e-11,objs_pd_e['mstar']/objs_pd_e['mvir']/f_bar,color= 'b',alpha = 0.5)
plt.bar(r_SSFR,objs_pd_r['mISM']/objs_pd_r['mvir']/f_bar,1e-11,objs_pd_r['mstar']/objs_pd_r['mvir']/f_bar,color= 'b',alpha = 0.5)
plt.bar(s_SSFR,objs_pd_s['mISM']/objs_pd_s['mvir']/f_bar,1e-11,objs_pd_s['mstar']/objs_pd_s['mvir']/f_bar,color= 'b',alpha = 0.5)

stars_plt = plt.bar(cm_SSFR,objs_pd_cm['mstar']/objs_pd_cm['mvir']/f_bar,width=1e-11,color= 'r',alpha = 0.5)
plt.bar(e_SSFR,objs_pd_e['mstar']/objs_pd_e['mvir']/f_bar,width=1e-11,color= 'r',alpha = 0.5)
plt.bar(r_SSFR,objs_pd_r['mstar']/objs_pd_r['mvir']/f_bar,width=1e-11,color= 'r',alpha = 0.5)
plt.bar(s_SSFR,objs_pd_s['mstar']/objs_pd_s['mvir']/f_bar,width=1e-11,color= 'r',alpha = 0.5)
plt.xlabel(r'sSFR')
plt.ylabel(r'$M/(f_{bary}\times M_{vir})$')
plt.legend([hot_plt,warm_plt,cool_plt,ISM_plt,stars_plt],['$10^6 K < T$','$10^5 K < T < 10^6 K$','$10^4 K < T < 10^5 K$','$T < 10^4 K$','Stars'],loc = 1,framealpha=0.5,fontsize = legendsize)
plt.axis([0, 1e-9, 0, 1])
plt.tight_layout()
plt.show()
plt.savefig(outbase + '_gascensus_vSSFR.png',dpi = dpi)
plt.close()

# ==============================================
ZSOLAR = 0.0130215
Zyield = 0.02788242 #0.015
objs_pd_cm_mZstar = objs_pd_cm['mstar']*ZSOLAR*10**objs_pd_cm['Zstellar_mean']
objs_pd_e_mZstar = objs_pd_e['mstar']*ZSOLAR*10**objs_pd_e['Zstellar_mean']
objs_pd_r_mZstar = objs_pd_r['mstar']*ZSOLAR*10**objs_pd_r['Zstellar_mean']
objs_pd_s_mZstar = objs_pd_s['mstar']*ZSOLAR*10**objs_pd_s['Zstellar_mean']

plt.figure(1,figsize=(plt_width,plt_width*aspect_ratio))
hot_plt = plt.bar(cm_mstar_log,objs_pd_cm['mZHot']/objs_pd_cm['mstarform']/Zyield,0.1,(objs_pd_cm['mZwarm'] + objs_pd_cm['mZCool'] + objs_pd_cm['mZISM'] + objs_pd_cm_mZstar)/objs_pd_cm['mstarform']/Zyield,color = 'goldenrod',alpha = 0.5)
plt.bar(e_mstar_log,objs_pd_e['mZHot']/objs_pd_e['mstarform']/Zyield,0.1,(objs_pd_e['mZwarm'] + objs_pd_e['mZCool'] + objs_pd_e['mZISM'] + objs_pd_e_mZstar)/objs_pd_e['mstarform']/Zyield,color = 'goldenrod',alpha = 0.5)
plt.bar(r_mstar_log,objs_pd_r['mZHot']/objs_pd_r['mstarform']/Zyield,0.1,(objs_pd_r['mZwarm'] + objs_pd_r['mZCool'] + objs_pd_r['mZISM'] + objs_pd_r_mZstar)/objs_pd_r['mstarform']/Zyield,color = 'goldenrod',alpha = 0.5)
plt.bar(s_mstar_log,objs_pd_s['mZHot']/objs_pd_s['mstarform']/Zyield,0.1,(objs_pd_s['mZwarm'] + objs_pd_s['mZCool'] + objs_pd_s['mZISM']+ objs_pd_s_mZstar)/objs_pd_s['mstarform']/Zyield,color = 'goldenrod',alpha = 0.5)

warm_plt = plt.bar(cm_mstar_log,objs_pd_cm['mZwarm']/objs_pd_cm['mstarform']/Zyield,0.1,(objs_pd_cm['mZCool'] + objs_pd_cm['mZISM'] + objs_pd_cm_mZstar)/objs_pd_cm['mstarform']/Zyield,color= 'g',alpha = 0.5)
plt.bar(e_mstar_log,objs_pd_e['mZwarm']/objs_pd_e['mstarform']/Zyield,0.1,(objs_pd_e['mZCool'] + objs_pd_e['mZISM'] + objs_pd_e_mZstar)/objs_pd_e['mstarform']/Zyield,color= 'g',alpha = 0.5)
plt.bar(r_mstar_log,objs_pd_r['mZwarm']/objs_pd_r['mstarform']/Zyield,0.1,(objs_pd_r['mZCool'] + objs_pd_r['mZISM'] + objs_pd_r_mZstar)/objs_pd_r['mstarform']/Zyield,color= 'g',alpha = 0.5)
plt.bar(s_mstar_log,objs_pd_s['mZwarm']/objs_pd_s['mstarform']/Zyield,0.1,(objs_pd_s['mZCool'] + objs_pd_s['mZISM'] + objs_pd_s_mZstar)/objs_pd_s['mstarform']/Zyield,color= 'g',alpha = 0.5)

cool_plt = plt.bar(cm_mstar_log,objs_pd_cm['mZCool']/objs_pd_cm['mstarform']/Zyield,0.1,(objs_pd_cm['mZISM'] + objs_pd_cm_mZstar)/objs_pd_cm['mstarform']/Zyield,color= 'darkviolet',alpha = 0.5)
plt.bar(e_mstar_log,objs_pd_e['mZCool']/objs_pd_e['mstarform']/Zyield,0.1,(objs_pd_e['mZISM'] + objs_pd_e_mZstar)/objs_pd_e['mstarform']/Zyield,color= 'darkviolet',alpha = 0.5)
plt.bar(r_mstar_log,objs_pd_r['mZCool']/objs_pd_r['mstarform']/Zyield,0.1,(objs_pd_r['mZISM'] + objs_pd_r_mZstar)/objs_pd_r['mstarform']/Zyield,color= 'darkviolet',alpha = 0.5)
plt.bar(s_mstar_log,objs_pd_s['mZCool']/objs_pd_s['mstarform']/Zyield,0.1,(objs_pd_s['mZISM'] + objs_pd_s_mZstar)/objs_pd_s['mstarform']/Zyield,color= 'darkviolet',alpha = 0.5)

ISM_plt = plt.bar(cm_mstar_log,objs_pd_cm['mZISM']/objs_pd_cm['mstarform']/Zyield,0.1,objs_pd_cm_mZstar/objs_pd_cm['mstarform']/Zyield,color= 'b',alpha = 0.5)
plt.bar(e_mstar_log,objs_pd_e['mZISM']/objs_pd_e['mstarform']/Zyield,0.1,objs_pd_e_mZstar/objs_pd_e['mstarform']/Zyield,color= 'b',alpha = 0.5)
plt.bar(r_mstar_log,objs_pd_r['mZISM']/objs_pd_r['mstarform']/Zyield,0.1,objs_pd_r_mZstar/objs_pd_r['mstarform']/Zyield,color= 'b',alpha = 0.5)
plt.bar(s_mstar_log,objs_pd_s['mZISM']/objs_pd_s['mstarform']/Zyield,0.1,objs_pd_s_mZstar/objs_pd_s['mstarform']/Zyield,color= 'b',alpha = 0.5)

stars_plt = plt.bar(cm_mstar_log,objs_pd_cm_mZstar/objs_pd_cm['mstarform']/Zyield,width=0.1,color= 'r',alpha = 0.5)
plt.bar(e_mstar_log,objs_pd_e_mZstar/objs_pd_e['mstarform']/Zyield,width=0.1,color= 'r',alpha = 0.5)
plt.bar(r_mstar_log,objs_pd_r_mZstar/objs_pd_r['mstarform']/Zyield,width=0.1,color= 'r',alpha = 0.5)
plt.bar(s_mstar_log,objs_pd_s_mZstar/objs_pd_s['mstarform']/Zyield,width=0.1,color= 'r',alpha = 0.5)
plt.xlabel(r'Log(M$_*/M_\odot)$')
plt.ylabel(r'$M_Z/(y \times M_{*})$')
plt.plot([2.5,9.5],[1,1],color = 'k',linestyle = '--')
plt.axis([2.5, 9.5, 2e-3, 20])
#plt.axis([3.5, 9.5, 0, 1])
plt.yscale('log')
plt.tight_layout()
plt.legend([hot_plt,warm_plt,cool_plt,ISM_plt,stars_plt],['$10^6 K < T$','$10^5 K < T < 10^6 K$','$10^4 K < T < 10^5 K$','$T < 10^4 K$','Stars'],loc = 1,ncol = 2,framealpha=0.5,fontsize = legendsize)
plt.tight_layout()
plt.show()
plt.savefig(outbase + '_Zcensus_vMstar.png',dpi = dpi)
plt.close()


plt.figure(1,figsize=(plt_width,plt_width*aspect_ratio))
hot_plt = plt.bar(cm_mvir_log,objs_pd_cm['mZHot']/objs_pd_cm['mstarform']/Zyield,0.05,(objs_pd_cm['mZwarm'] + objs_pd_cm['mZCool'] + objs_pd_cm['mZISM'] + objs_pd_cm_mZstar)/objs_pd_cm['mstarform']/Zyield,color = 'goldenrod',alpha = 0.5)
plt.bar(e_mvir_log,objs_pd_e['mZHot']/objs_pd_e['mstarform']/Zyield,0.05,(objs_pd_e['mZwarm'] + objs_pd_e['mZCool'] + objs_pd_e['mZISM'] + objs_pd_e_mZstar)/objs_pd_e['mstarform']/Zyield,color = 'goldenrod',alpha = 0.5)
plt.bar(r_mvir_log,objs_pd_r['mZHot']/objs_pd_r['mstarform']/Zyield,0.05,(objs_pd_r['mZwarm'] + objs_pd_r['mZCool'] + objs_pd_r['mZISM'] + objs_pd_r_mZstar)/objs_pd_r['mstarform']/Zyield,color = 'goldenrod',alpha = 0.5)
plt.bar(s_mvir_log,objs_pd_s['mZHot']/objs_pd_s['mstarform']/Zyield,0.05,(objs_pd_s['mZwarm'] + objs_pd_s['mZCool'] + objs_pd_s['mZISM']+ objs_pd_s_mZstar)/objs_pd_s['mstarform']/Zyield,color = 'goldenrod',alpha = 0.5)

warm_plt = plt.bar(cm_mvir_log,objs_pd_cm['mZwarm']/objs_pd_cm['mstarform']/Zyield,0.05,(objs_pd_cm['mZCool'] + objs_pd_cm['mZISM'] + objs_pd_cm_mZstar)/objs_pd_cm['mstarform']/Zyield,color= 'g',alpha = 0.5)
plt.bar(e_mvir_log,objs_pd_e['mZwarm']/objs_pd_e['mstarform']/Zyield,0.05,(objs_pd_e['mZCool'] + objs_pd_e['mZISM'] + objs_pd_e_mZstar)/objs_pd_e['mstarform']/Zyield,color= 'g',alpha = 0.5)
plt.bar(r_mvir_log,objs_pd_r['mZwarm']/objs_pd_r['mstarform']/Zyield,0.05,(objs_pd_r['mZCool'] + objs_pd_r['mZISM'] + objs_pd_r_mZstar)/objs_pd_r['mstarform']/Zyield,color= 'g',alpha = 0.5)
plt.bar(s_mvir_log,objs_pd_s['mZwarm']/objs_pd_s['mstarform']/Zyield,0.05,(objs_pd_s['mZCool'] + objs_pd_s['mZISM'] + objs_pd_s_mZstar)/objs_pd_s['mstarform']/Zyield,color= 'g',alpha = 0.5)

cool_plt = plt.bar(cm_mvir_log,objs_pd_cm['mZCool']/objs_pd_cm['mstarform']/Zyield,0.05,(objs_pd_cm['mZISM'] + objs_pd_cm_mZstar)/objs_pd_cm['mstarform']/Zyield,color= 'darkviolet',alpha = 0.5)
plt.bar(e_mvir_log,objs_pd_e['mZCool']/objs_pd_e['mstarform']/Zyield,0.05,(objs_pd_e['mZISM'] + objs_pd_e_mZstar)/objs_pd_e['mstarform']/Zyield,color= 'darkviolet',alpha = 0.5)
plt.bar(r_mvir_log,objs_pd_r['mZCool']/objs_pd_r['mstarform']/Zyield,0.05,(objs_pd_r['mZISM'] + objs_pd_r_mZstar)/objs_pd_r['mstarform']/Zyield,color= 'darkviolet',alpha = 0.5)
plt.bar(s_mvir_log,objs_pd_s['mZCool']/objs_pd_s['mstarform']/Zyield,0.05,(objs_pd_s['mZISM'] + objs_pd_s_mZstar)/objs_pd_s['mstarform']/Zyield,color= 'darkviolet',alpha = 0.5)

ISM_plt = plt.bar(cm_mvir_log,objs_pd_cm['mZISM']/objs_pd_cm['mstarform']/Zyield,0.05,objs_pd_cm_mZstar/objs_pd_cm['mstarform']/Zyield,color= 'b',alpha = 0.5)
plt.bar(e_mvir_log,objs_pd_e['mZISM']/objs_pd_e['mstarform']/Zyield,0.05,objs_pd_e_mZstar/objs_pd_e['mstarform']/Zyield,color= 'b',alpha = 0.5)
plt.bar(r_mvir_log,objs_pd_r['mZISM']/objs_pd_r['mstarform']/Zyield,0.05,objs_pd_r_mZstar/objs_pd_r['mstarform']/Zyield,color= 'b',alpha = 0.5)
plt.bar(s_mvir_log,objs_pd_s['mZISM']/objs_pd_s['mstarform']/Zyield,0.05,objs_pd_s_mZstar/objs_pd_s['mstarform']/Zyield,color= 'b',alpha = 0.5)

stars_plt = plt.bar(cm_mvir_log,objs_pd_cm_mZstar/objs_pd_cm['mstarform']/Zyield,width=0.05,color= 'r',alpha = 0.5)
plt.bar(e_mvir_log,objs_pd_e_mZstar/objs_pd_e['mstarform']/Zyield,width=0.05,color= 'r',alpha = 0.5)
plt.bar(r_mvir_log,objs_pd_r_mZstar/objs_pd_r['mstarform']/Zyield,width=0.05,color= 'r',alpha = 0.5)
plt.bar(s_mvir_log,objs_pd_s_mZstar/objs_pd_s['mstarform']/Zyield,width=0.05,color= 'r',alpha = 0.5)
plt.yscale('log')
plt.xlabel(r'Log($M_{vir}/M_\odot)$')
plt.ylabel(r'$M/(f_{bary} \times M_{vir})$')
plt.plot([8,11.2],[1,1],color = 'k',linestyle = '--')
plt.axis([8.1, 11.2, 2e-3, 20])
plt.tight_layout()
plt.legend([hot_plt,warm_plt,cool_plt,ISM_plt,stars_plt],['$10^6 K < T$','$10^5 K < T < 10^6 K$','$10^4 K < T < 10^5 K$','$T < 10^4 K$','Stars'],loc = 1,ncol = 2,fontsize = legendsize)
plt.show()
plt.savefig(outbase + '_Zcensus_vMvir.png',dpi = dpi)
plt.close()
