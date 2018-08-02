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

min_nstar =  100 
dInitStarMass   = 1.82699e-13*2.310e15
min_mstar = min_nstar*dInitStarMass

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

#Sandra
tfile = prefix + 'h148.cosmo50PLK.3072g/h148.cosmo50PLK.3072g3HbwK1BH/h148.cosmo50PLK.3072g3HbwK1BH.004096/h148.cosmo50PLK.3072g3HbwK1BH.004096'
objs_sandra = []
f=open(tfile + '.data', 'rb')
while 1:
    try:
        objs_sandra.append(pickle.load(f))
    except EOFError:
        break
        
f.close()
objs_pd_sandra = pd.DataFrame(objs_sandra)


#Ruth
tfile = prefix + 'h229.cosmo50PLK.3072g/h229.cosmo50PLK.3072gst5HbwK1BH/h229.cosmo50PLK.3072gst5HbwK1BH.004096/h229.cosmo50PLK.3072gst5HbwK1BH.004096'
objs_ruth = []
f=open(tfile + '.data', 'rb')
while 1:
    try:
        objs_ruth.append(pickle.load(f))
    except EOFError:
        break
        
f.close()
objs_pd_ruth = pd.DataFrame(objs_ruth)


#Sonia
tfile = prefix + 'h242.cosmo50PLK.3072g/h242.cosmo50PLK.3072gst5HbwK1BH/h242.cosmo50PLK.3072gst5HbwK1BH.004096/h242.cosmo50PLK.3072gst5HbwK1BH.004096'
objs_sonia = []
f=open(tfile + '.data', 'rb')
while 1:
    try:
        objs_sonia.append(pickle.load(f))
    except EOFError:
        break
        
f.close()
objs_pd_sonia = pd.DataFrame(objs_sonia)

#Elena
tfile = prefix + 'h329.cosmo50PLK.3072g/h329.cosmo50PLK.3072gst5HbwK1BH/h329.cosmo50PLK.3072gst5HbwK1BH.004096/h329.cosmo50PLK.3072gst5HbwK1BH.004096'
objs_elena = []
f=open(tfile + '.data', 'rb')
while 1:
    try:
        objs_elena.append(pickle.load(f))
    except EOFError:
        break
        
f.close()
objs_pd_elena = pd.DataFrame(objs_elena)


#Markers and colors for the simulations
values = range(0,20)
cNorm  = colors.Normalize(vmin=0, vmax = values[-1]+1)
"""
markers = ['o','o','o','o','o','o','o','o']
colormap = 'tab20'
colors_cm = 'c'
colors_r = 'c'
colors_e = 'c'
colors_s = 'c'
colors_sandra = 'c'
colors_ruth = 'c'
colors_sonia = 'c'
colors_elena = 'c'
colors_cm_f = 'none'
colors_r_f = 'none'
colors_e_f = 'none'
colors_s_f = 'none'
colors_sandra_f = 'none'
colors_ruth_f = 'none'
colors_sonia_f = 'none'
colors_elena_f = 'none'
"""

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


#Recent SFR vs Average historical SFR
#McQuinn+ 2015, Figures 4 and 5
f = open(dataprefix+'McQuinn2015.txt', 'r')
mcquinndata = []
for line in f:
    line = line.strip()
    columns = line.split()
    if len(columns) == 22 and columns[0] != 'Galaxy':
        source = {}
        source['gal'] = columns[0]
        source['SFRrecent'] = float(columns[1])
        source['SFRrecent_perr'] = float(columns[2])
        source['SFRrecent_merr'] = float(columns[3])
        source['SFRlifetime'] = float(columns[4])
        source['SFRlifetime_perr'] = float(columns[5])
        source['SFRlifetime_merr'] = float(columns[6])
        source['recent_by_lifetime'] = float(columns[7])
        source['recent_by_lifetime_perr'] = float(columns[8])
        source['recent_by_lifetime_merr'] = float(columns[9])
        source['mstar'] = float(columns[10])
        source['mstar_perr'] = float(columns[11])
        source['mstar_merr'] = float(columns[12])
        source['mHI'] = float(columns[13])
        source['mHI_perr'] = float(columns[14])
        source['mHI_merr'] = float(columns[15])
        source['mHI_by_mstar'] = float(columns[16])
        source['mHI_by_mstar_perr'] = float(columns[17])
        source['mHI_by_mstar_merr'] = float(columns[18])
        source['tau'] = float(columns[19])
        source['tau_perr'] = float(columns[20])
        source['tau_merr'] = float(columns[21])
        mcquinndata.append(source)
f.close()
mcquinndata = pd.DataFrame(mcquinndata)

"""
plt.figure(1)
#Fix upper limits in McQuinn data (marked as a zero in the minus error term) before publishing
mcquinn_plt = plt.errorbar(1e-3*mcquinndata['SFRlifetime'],1e-3*mcquinndata['SFRrecent'],xerr=[1e-3*mcquinndata['SFRlifetime_merr'],1e-3*mcquinndata['SFRlifetime_perr']],yerr=[1e-3*mcquinndata['SFRrecent_merr'],1e-3*mcquinndata['SFRrecent_perr']],color = "grey",fmt="o",zorder = 1)
cptmarvel_plt = plt.scatter(np.array(objs_pd_cm['mstar'])/np.array(objs_pd_cm['time'])/1e9,np.array(objs_pd_cm['SFR']),marker = markers[0], label = 'Cpt. Marvel',zorder = 2,c = colors_cm, cmap = colormap,norm = cNorm)
rogue_plt = plt.scatter(np.array(objs_pd_r['mstar'])/np.array(objs_pd_r['time'])/1e9,np.array(objs_pd_r['SFR']),marker = markers[1],zorder = 3,c = colors_r, cmap = colormap, norm=cNorm)
elektra_plt = plt.scatter(np.array(objs_pd_e['mstar'])/np.array(objs_pd_e['time'])/1e9,np.array(objs_pd_e['SFR']),marker = markers[2],zorder = 4,c = colors_e, cmap = colormap, norm=cNorm)
storm_plt = plt.scatter(np.array(objs_pd_s['mstar'])/np.array(objs_pd_s['time'])/1e9,np.array(objs_pd_s['SFR']),marker = markers[3],zorder = 5,c = colors_s, cmap = colormap, norm=cNorm)
sandra_plt = plt.scatter(np.array(objs_pd_sandra['mstar'])/np.array(objs_pd_sandra['time'])/1e9,np.array(objs_pd_sandra['SFR']),marker = markers[4],zorder = 5,c = colors_sandra, cmap = colormap, norm=cNorm)
ruth_plt = plt.scatter(np.array(objs_pd_ruth['mstar'])/np.array(objs_pd_ruth['time'])/1e9,np.array(objs_pd_ruth['SFR']),marker = markers[5],zorder = 5,c = colors_ruth, cmap = colormap, norm=cNorm)
sonia_plt = plt.scatter(np.array(objs_pd_sonia['mstar'])/np.array(objs_pd_sonia['time'])/1e9,np.array(objs_pd_sonia['SFR']),marker = markers[6],zorder = 5,c = colors_sonia, cmap = colormap, norm=cNorm)
elena_plt = plt.scatter(np.array(objs_pd_elena['mstar'])/np.array(objs_pd_elena['time'])/1e9,np.array(objs_pd_elena['SFR']),marker = markers[7],zorder = 5,c = colors_elena, cmap = colormap, norm=cNorm)
plt.plot([1e-5,100],[1e-5,100],color = 'k')
plt.xlabel(r'Lifetime SFR [M$_\odot$ yr$^{-1}$]')
plt.ylabel(r'Recent SFR [M$_\odot$ yr$^{-1}$]')
plt.xscale('log')
plt.yscale('log')
#plt.axis([1e-5, 3e-1, 1e-5, 3e-1])
plt.axis([1e-5, 40, 1e-5, 40])
#plt.legend([mcquinn_plt, cptmarvel_plt, rogue_plt, elektra_plt, storm_plt],['McQuinn+ 2015','Cpt. Marvel','Rogue','Elektra','Storm'],loc=2)
plt.legend([mcquinn_plt, cptmarvel_plt, rogue_plt, elektra_plt, storm_plt, sandra_plt, ruth_plt, sonia_plt, elena_plt],['McQuinn+ 2015','Cpt. Marvel','Rogue','Elektra','Storm','Sandra','Ruth','Sonia','Elena'],loc=2)
plt.show()
plt.savefig(outbase + '.sfr.png')
plt.close()

#Star forming main sequence
#McQuinn+ 2015, Figure 7 and 8
plt.figure(1,figsize=(plt_width,plt_width))
#mcquinn_plt = plt.errorbar(mcquinndata['mstar']*1e6,1e-3*mcquinndata['SFRrecent']/mcquinndata['mstar']/1e6,xerr=[mcquinndata['mstar_merr']*1e6,mcquinndata['mstar_perr']*1e6],yerr=[1e-3*mcquinndata['SFRrecent_merr']/mcquinndata['mstar']/1e6,1e-3*mcquinndata['SFRrecent_perr']/mcquinndata['mstar']/1e6],color = "grey",fmt="o",zorder = 1)
cptmarvel_plt = plt.scatter(np.array(objs_pd_cm['mstar']),np.array(objs_pd_cm['SFR'])/np.array(objs_pd_cm['mstar']),marker = markers[0], label = 'Cpt. Marvel' ,zorder = 2,edgecolors = 'k', cmap = colormap, norm=cNorm,facecolors = colors_cm_f)
rogue_plt = plt.scatter(np.array(objs_pd_r['mstar']),np.array(objs_pd_r['SFR'])/np.array(objs_pd_r['mstar']),marker = markers[1],zorder = 3,edgecolors = 'k', cmap = colormap, norm=cNorm, facecolors = colors_r_f)
elektra_plt = plt.scatter(np.array(objs_pd_e['mstar']),np.array(objs_pd_e['SFR'])/np.array(objs_pd_e['mstar']),marker = markers[2],zorder = 4,edgecolors = 'k', cmap = colormap, norm=cNorm, facecolors = colors_e_f)
storm_plt = plt.scatter(np.array(objs_pd_s['mstar']),np.array(objs_pd_s['SFR'])/np.array(objs_pd_s['mstar']),marker = markers[3],zorder = 5,edgecolors = 'k', cmap = colormap, norm=cNorm, facecolors = colors_s_f)
#sandra_plt = plt.scatter(np.array(objs_pd_sandra['mstar']),np.array(objs_pd_sandra['SFR'])/np.array(objs_pd_sandra['mstar']),marker = markers[4],zorder = 5,edgecolors = colors_sandra, cmap = colormap, norm=cNorm, facecolors = colors_sandra_f)
#ruth_plt = plt.scatter(np.array(objs_pd_ruth['mstar']),np.array(objs_pd_ruth['SFR'])/np.array(objs_pd_ruth['mstar']),marker = markers[5],zorder = 5,edgecolors = colors_ruth, cmap = colormap, norm=cNorm, facecolors = colors_ruth_f)
#sonia_plt = plt.scatter(np.array(objs_pd_sonia['mstar']),np.array(objs_pd_sonia['SFR'])/np.array(objs_pd_sonia['mstar']),marker = markers[6],zorder = 5,edgecolors = colors_sonia, cmap = colormap, norm=cNorm, facecolors = colors_sonia_f)
#elena_plt = plt.scatter(np.array(objs_pd_elena['mstar']),np.array(objs_pd_elena['SFR'])/np.array(objs_pd_elena['mstar']),marker = markers[7],zorder = 5,edgecolors = colors_elena, cmap = colormap, norm=cNorm, facecolors = colors_elena_f)
plt.xlabel(r'M$_*$/M$_\odot$')
plt.ylabel(r'SSFR [yr$^{-1}$]')
plt.xscale('log')
plt.yscale('log')
#plt.axis([1e4, 3e9, 5e-12, 1e-9])
#plt.axis([1e4, 1e12, 1e-12, 1e-9])
plt.axis([316227, 3162277660, 1e-13, 1e-8])
#plt.legend([mcquinn_plt, cptmarvel_plt, rogue_plt, elektra_plt, storm_plt],['McQuinn+ 2015','Cpt. Marvel','Rogue','Elektra','Storm'],loc=3)
#plt.legend([mcquinn_plt, cptmarvel_plt],['McQuinn+ 2015','Simulations'],loc=4)
#plt.legend([mcquinn_plt, cptmarvel_plt, rogue_plt, elektra_plt, storm_plt, sandra_plt, ruth_plt, sonia_plt, elena_plt],['McQuinn+ 2015','Cpt. Marvel','Rogue','Elektra','Storm','Sandra','Ruth','Sonia','Elena'],loc=3)
#plt.legend([mcquinn_plt,cptmarvel_plt],['McQuinn et al. 2015','Simulations'],loc = 4,fontsize = legendsize)
plt.tight_layout()
plt.show()
plt.savefig(outbase + '.sfms_ssfr.png')
plt.close()

#Star forming main sequence
#McQuinn+ 2015, Figure 7 and 8
plt.figure(1)
mcquinn_plt = plt.errorbar(mcquinndata['mstar']*1e6,1e-3*mcquinndata['SFRrecent'],xerr=[mcquinndata['mstar_merr']*1e6,mcquinndata['mstar_perr']*1e6],yerr=[1e-3*mcquinndata['SFRrecent_merr'],1e-3*mcquinndata['SFRrecent_perr']],color = "grey",fmt="o",zorder = 1)
cptmarvel_plt = plt.scatter(np.array(objs_pd_cm['mstar']),np.array(objs_pd_cm['SFR']),marker = markers[0], label = 'Cpt. Marvel' ,zorder = 2,edgecolors = colors_cm, cmap = colormap, norm=cNorm,facecolors = colors_cm_f)
rogue_plt = plt.scatter(np.array(objs_pd_r['mstar']),np.array(objs_pd_r['SFR']),marker = markers[1],zorder = 3,edgecolors = colors_r, cmap = colormap, norm=cNorm, facecolors = colors_r_f)
elektra_plt = plt.scatter(np.array(objs_pd_e['mstar']),np.array(objs_pd_e['SFR']),marker = markers[2],zorder = 4,edgecolors = colors_e, cmap = colormap, norm=cNorm, facecolors = colors_e_f)
storm_plt = plt.scatter(np.array(objs_pd_s['mstar']),np.array(objs_pd_s['SFR']),marker = markers[3],zorder = 5,edgecolors = colors_s, cmap = colormap, norm=cNorm, facecolors = colors_s_f)
#sandra_plt = plt.scatter(np.array(objs_pd_sandra['mstar']),np.array(objs_pd_sandra['SFR'])/np.array(objs_pd_sandra['mstar']),marker = markers[4],zorder = 5,edgecolors = colors_sandra, cmap = colormap, norm=cNorm, facecolors = colors_sandra_f)
#ruth_plt = plt.scatter(np.array(objs_pd_ruth['mstar']),np.array(objs_pd_ruth['SFR'])/np.array(objs_pd_ruth['mstar']),marker = markers[5],zorder = 5,edgecolors = colors_ruth, cmap = colormap, norm=cNorm, facecolors = colors_ruth_f)
#sonia_plt = plt.scatter(np.array(objs_pd_sonia['mstar']),np.array(objs_pd_sonia['SFR'])/np.array(objs_pd_sonia['mstar']),marker = markers[6],zorder = 5,edgecolors = colors_sonia, cmap = colormap, norm=cNorm, facecolors = colors_sonia_f)
#elena_plt = plt.scatter(np.array(objs_pd_elena['mstar']),np.array(objs_pd_elena['SFR'])/np.array(objs_pd_elena['mstar']),marker = markers[7],zorder = 5,edgecolors = colors_elena, cmap = colormap, norm=cNorm, facecolors = colors_elena_f)
plt.xlabel(r'M$_*$/M$_\odot$')
plt.ylabel(r'SFR [M$_\odot$ yr$^{-1}$]')
plt.xscale('log')
plt.yscale('log')
#plt.legend([mcquinn_plt, cptmarvel_plt, rogue_plt, elektra_plt, storm_plt],['McQuinn+ 2015','Cpt. Marvel','Rogue','Elektra','Storm'],loc=3)
#plt.legend([mcquinn_plt, cptmarvel_plt],['McQuinn+ 2015','Simulations'],loc=4)
#plt.legend([mcquinn_plt, cptmarvel_plt, rogue_plt, elektra_plt, storm_plt, sandra_plt, ruth_plt, sonia_plt, elena_plt],['McQuinn+ 2015','Cpt. Marvel','Rogue','Elektra','Storm','Sandra','Ruth','Sonia','Elena'],loc=3)
plt.legend([mcquinn_plt,cptmarvel_plt],['McQuinn et al. 2015','Simulations'],loc = 2,fontsize = legendsize)
plt.tight_layout()
plt.show()
plt.savefig(outbase + '.sfms.png')
plt.close()

#HI fraction vs stellar mass
#McQuinn+ 2015, Figures 9 and 10 for SHIELD galaxies
plt.figure(1)
mcquinn_plt = plt.errorbar(mcquinndata['mstar']*1e6,mcquinndata['mHI_by_mstar'],xerr=[mcquinndata['mstar_merr']*1e6,mcquinndata['mstar_perr']*1e6],yerr=[mcquinndata['mHI_by_mstar_merr'],mcquinndata['mHI_by_mstar_perr']],color = "grey",fmt="o",zorder = 1)
cptmarvel_plt = plt.scatter(np.array(objs_pd_cm['mstar']),np.array(objs_pd_cm['mHI'])/np.array(objs_pd_cm['mstar']),marker = markers[0], label = 'Cpt. Marvel' ,zorder = 2,c = colors_cm, cmap = colormap, norm=cNorm)
rogue_plt = plt.scatter(np.array(objs_pd_r['mstar']),np.array(objs_pd_r['mHI'])/np.array(objs_pd_r['mstar']),marker = markers[1], label = 'Rogue', zorder = 3,c = colors_r, cmap = colormap, norm=cNorm)
elektra_plt = plt.scatter(np.array(objs_pd_e['mstar']),np.array(objs_pd_e['mHI'])/np.array(objs_pd_e['mstar']),marker = markers[2], label = 'Elektra',zorder = 4,c = colors_e, cmap = colormap, norm=cNorm)
storm_plt = plt.scatter(np.array(objs_pd_s['mstar']),np.array(objs_pd_s['mHI'])/np.array(objs_pd_s['mstar']),marker = markers[3], label = 'Storm', zorder = 5,c = colors_s, cmap = colormap, norm=cNorm)
sandra_plt = plt.scatter(np.array(objs_pd_sandra['mstar']),np.array(objs_pd_sandra['mHI'])/np.array(objs_pd_sandra['mstar']),marker = markers[4], label = 'Sandra', zorder = 5,c = colors_sandra, cmap = colormap, norm=cNorm)
ruth_plt = plt.scatter(np.array(objs_pd_ruth['mstar']),np.array(objs_pd_ruth['mHI'])/np.array(objs_pd_ruth['mstar']),marker = markers[5], label = 'Ruth', zorder = 5,c = colors_ruth, cmap = colormap, norm=cNorm)
sonia_plt = plt.scatter(np.array(objs_pd_sonia['mstar']),np.array(objs_pd_sonia['mHI'])/np.array(objs_pd_sonia['mstar']),marker = markers[6], label = 'Sonia', zorder = 5,c = colors_sonia, cmap = colormap, norm=cNorm)
elena_plt = plt.scatter(np.array(objs_pd_elena['mstar']),np.array(objs_pd_elena['mHI'])/np.array(objs_pd_elena['mstar']),marker = markers[7], label = 'Elena', zorder = 5,c = colors_elena, cmap = colormap, norm=cNorm)
plt.xlabel(r'M$_*$/M$_\odot$')
plt.ylabel(r'M$_{HI}$/M$_*$')
plt.xscale('log')
#plt.axis([1e6, 3e9, 0, 8.5])
plt.axis([1e6, 1e12, 0, 8.5])
#plt.legend([mcquinn_plt, cptmarvel_plt, rogue_plt, elektra_plt, storm_plt],['McQuinn+ 2015','Cpt. Marvel','Rogue','Elektra','Storm'],loc=1)
plt.legend([mcquinn_plt, cptmarvel_plt, rogue_plt, elektra_plt, storm_plt, sandra_plt, ruth_plt, sonia_plt, elena_plt],['McQuinn+ 2015','Cpt. Marvel','Rogue','Elektra','Storm','Sandra','Ruth','Sonia','Elena'],loc=1)
plt.show()
plt.savefig(outbase + '.gasfrac_mstar.png')
plt.clf()

plt.figure(1)
mcquinn_plt = plt.errorbar(mcquinndata['mstar']*1e6,mcquinndata['mHI_by_mstar'],xerr=[mcquinndata['mstar_merr']*1e6,mcquinndata['mstar_perr']*1e6],yerr=[mcquinndata['mHI_by_mstar_merr'],mcquinndata['mHI_by_mstar_perr']],color = "grey",fmt="o",zorder = 1)
cptmarvel_plt = plt.scatter(np.array(objs_pd_cm['mstar']),np.array(objs_pd_cm['mHI'])/np.array(objs_pd_cm['mstar']),marker = markers[0], label = 'Cpt. Marvel' ,zorder = 2,c = colors_cm, cmap = colormap, norm=cNorm)
rogue_plt = plt.scatter(np.array(objs_pd_r['mstar']),np.array(objs_pd_r['mHI'])/np.array(objs_pd_r['mstar']),marker = markers[1], label = 'Rogue', zorder = 3,c = colors_r, cmap = colormap, norm=cNorm)
elektra_plt = plt.scatter(np.array(objs_pd_e['mstar']),np.array(objs_pd_e['mHI'])/np.array(objs_pd_e['mstar']),marker = markers[2], label = 'Elektra',zorder = 4,c = colors_e, cmap = colormap, norm=cNorm)
storm_plt = plt.scatter(np.array(objs_pd_s['mstar']),np.array(objs_pd_s['mHI'])/np.array(objs_pd_s['mstar']),marker = markers[3], label = 'Storm', zorder = 5,c = colors_s, cmap = colormap, norm=cNorm)
sandra_plt = plt.scatter(np.array(objs_pd_sandra['mstar']),np.array(objs_pd_sandra['mHI'])/np.array(objs_pd_sandra['mstar']),marker = markers[4], label = 'Storm', zorder = 5,c = colors_sandra, cmap = colormap, norm=cNorm)
ruth_plt = plt.scatter(np.array(objs_pd_ruth['mstar']),np.array(objs_pd_ruth['mHI'])/np.array(objs_pd_ruth['mstar']),marker = markers[5], label = 'Ruth', zorder = 5,c = colors_ruth, cmap = colormap, norm=cNorm)
sonia_plt = plt.scatter(np.array(objs_pd_sonia['mstar']),np.array(objs_pd_sonia['mHI'])/np.array(objs_pd_sonia['mstar']),marker = markers[6], label = 'Sonia', zorder = 5,c = colors_sonia, cmap = colormap, norm=cNorm)
elena_plt = plt.scatter(np.array(objs_pd_elena['mstar']),np.array(objs_pd_elena['mHI'])/np.array(objs_pd_elena['mstar']),marker = markers[7], label = 'Elena', zorder = 5,c = colors_elena, cmap = colormap, norm=cNorm)
plt.xlabel(r'M$_*$/M$_\odot$')
plt.ylabel(r'M$_{HI}$/M$_*$')
plt.xscale('log')
#plt.axis([1e3, 3e9, 0, 200])
plt.axis([1e3, 1e12, 0, 200])
#plt.legend([mcquinn_plt, cptmarvel_plt, rogue_plt, elektra_plt, storm_plt],['McQuinn+ 2015','Cpt. Marvel','Rogue','Elektra','Storm'],loc=1)
plt.legend([mcquinn_plt, cptmarvel_plt, rogue_plt, elektra_plt, storm_plt, sandra_plt, ruth_plt, sonia_plt, elena_plt],['McQuinn+ 2015','Cpt. Marvel','Rogue','Elektra','Storm','Sandra','Ruth','Sonia','Elena'],loc=1)
plt.show()
plt.savefig(outbase + '.gasfrac_mstar_all.png')
plt.clf()

#Specific SFR vs gas ratio
#McQuinn+ 2015, Figure 11 and 12
plt.figure(1)
mcquinn_plt = plt.errorbar(mcquinndata['mHI_by_mstar'],1e-9*mcquinndata['SFRrecent']/mcquinndata['mstar'],xerr=[mcquinndata['mHI_by_mstar_merr'],mcquinndata['mHI_by_mstar_perr']],yerr=[1e-9*mcquinndata['SFRrecent_merr']/mcquinndata['mstar'],1e-9*mcquinndata['SFRrecent_perr']/mcquinndata['mstar_perr']],color = "grey",fmt="o",zorder = 1)
cptmarvel_plt = plt.scatter(np.array(objs_pd_cm['mHI'])/np.array(objs_pd_cm['mstar']),np.array(objs_pd_cm['SFR'])/np.array(objs_pd_cm['mstar']),marker =  markers[0], label = 'Cpt. Marvel' ,zorder = 2,c = colors_cm, cmap = colormap, norm=cNorm)
rogue_plt = plt.scatter(np.array(objs_pd_r['mHI'])/np.array(objs_pd_r['mstar']),np.array(objs_pd_r['SFR'])/np.array(objs_pd_r['mstar']),marker = markers[1], label = 'Rogue', zorder = 3,c = colors_r, cmap = colormap, norm=cNorm)
elektra_plt = plt.scatter(np.array(objs_pd_e['mHI'])/np.array(objs_pd_e['mstar']),np.array(objs_pd_e['SFR'])/np.array(objs_pd_e['mstar']),marker = markers[2], label = 'Elektra',zorder = 4,c = colors_e, cmap = colormap, norm=cNorm)
storm_plt = plt.scatter(np.array(objs_pd_s['mHI'])/np.array(objs_pd_s['mstar']),np.array(objs_pd_s['SFR'])/np.array(objs_pd_s['mstar']),marker = markers[3], label = 'Storm', zorder = 5,c = colors_s, cmap = colormap, norm=cNorm)
sandra_plt = plt.scatter(np.array(objs_pd_sandra['mHI'])/np.array(objs_pd_sandra['mstar']),np.array(objs_pd_sandra['SFR'])/np.array(objs_pd_sandra['mstar']),marker = markers[4], label = 'Sandra', zorder = 5,c = colors_sandra, cmap = colormap, norm=cNorm)
ruth_plt = plt.scatter(np.array(objs_pd_ruth['mHI'])/np.array(objs_pd_ruth['mstar']),np.array(objs_pd_ruth['SFR'])/np.array(objs_pd_ruth['mstar']),marker = markers[5], label = 'Ruth', zorder = 5,c = colors_ruth, cmap = colormap, norm=cNorm)
sonia_plt = plt.scatter(np.array(objs_pd_sonia['mHI'])/np.array(objs_pd_sonia['mstar']),np.array(objs_pd_sonia['SFR'])/np.array(objs_pd_sonia['mstar']),marker = markers[6], label = 'Sonia', zorder = 5,c = colors_sonia, cmap = colormap, norm=cNorm)
elena_plt = plt.scatter(np.array(objs_pd_elena['mHI'])/np.array(objs_pd_elena['mstar']),np.array(objs_pd_elena['SFR'])/np.array(objs_pd_elena['mstar']),marker = markers[7], label = 'Elena', zorder = 5,c = colors_elena, cmap = colormap, norm=cNorm)
plt.xlabel(r'M$_{HI}$/M$_*$')
plt.ylabel(r'SSFR [yr$^{-1}$]')
plt.yscale('log')
plt.axis([0, 8.5, 3e-12, 3e-9])
#plt.legend([mcquinn_plt, cptmarvel_plt, rogue_plt, elektra_plt, storm_plt],['McQuinn+ 2015','Cpt. Marvel','Rogue','Elektra','Storm'],loc=4)
plt.legend([mcquinn_plt, cptmarvel_plt, rogue_plt, elektra_plt, storm_plt, sandra_plt, ruth_plt, sonia_plt, elena_plt],['McQuinn+ 2015','Cpt. Marvel','Rogue','Elektra','Storm','Sandra','Ruth','Sonia','Elena'],loc=4)
plt.show()
plt.savefig(outbase + '.ssfr_gasfrac.png')
plt.clf()

plt.figure(1)
mcquinn_plt = plt.errorbar(mcquinndata['mHI_by_mstar'],1e-9*mcquinndata['SFRrecent']/mcquinndata['mstar'],xerr=[mcquinndata['mHI_by_mstar_merr'],mcquinndata['mHI_by_mstar_perr']],yerr=[1e-9*mcquinndata['SFRrecent_merr']/mcquinndata['mstar_merr'],1e-9*mcquinndata['SFRrecent_perr']/mcquinndata['mstar_perr']],color = "grey",fmt="o",zorder = 1)
cptmarvel_plt = plt.scatter(np.array(objs_pd_cm['mHI'])/np.array(objs_pd_cm['mstar']),np.array(objs_pd_cm['SFR'])/np.array(objs_pd_cm['mstar']),marker = markers[0], label = 'Cpt. Marvel' ,zorder = 2,c = colors_cm, cmap = colormap, norm=cNorm)
rogue_plt = plt.scatter(np.array(objs_pd_r['mHI'])/np.array(objs_pd_r['mstar']),np.array(objs_pd_r['SFR'])/np.array(objs_pd_r['mstar']),marker = markers[1], label = 'Rogue', zorder = 3,c = colors_r, cmap = colormap, norm=cNorm)
elektra_plt = plt.scatter(np.array(objs_pd_e['mHI'])/np.array(objs_pd_e['mstar']),np.array(objs_pd_e['SFR'])/np.array(objs_pd_e['mstar']),marker = markers[2], label = 'Elektra',zorder = 4,c = colors_e, cmap = colormap, norm=cNorm)
storm_plt = plt.scatter(np.array(objs_pd_s['mHI'])/np.array(objs_pd_s['mstar']),np.array(objs_pd_s['SFR'])/np.array(objs_pd_s['mstar']),marker = markers[3], label = 'Storm', zorder = 5,c = colors_s, cmap = colormap, norm=cNorm)
sandra_plt = plt.scatter(np.array(objs_pd_sandra['mHI'])/np.array(objs_pd_sandra['mstar']),np.array(objs_pd_sandra['SFR'])/np.array(objs_pd_sandra['mstar']),marker = markers[4], label = 'Sandra', zorder = 5,c = colors_sandra, cmap = colormap, norm=cNorm)
ruth_plt = plt.scatter(np.array(objs_pd_ruth['mHI'])/np.array(objs_pd_ruth['mstar']),np.array(objs_pd_ruth['SFR'])/np.array(objs_pd_ruth['mstar']),marker = markers[5], label = 'Ruth', zorder = 5,c = colors_ruth, cmap = colormap, norm=cNorm)
sonia_plt = plt.scatter(np.array(objs_pd_sonia['mHI'])/np.array(objs_pd_sonia['mstar']),np.array(objs_pd_sonia['SFR'])/np.array(objs_pd_sonia['mstar']),marker = markers[6], label = 'Sonia', zorder = 5,c = colors_sonia, cmap = colormap, norm=cNorm)
elena_plt = plt.scatter(np.array(objs_pd_elena['mHI'])/np.array(objs_pd_elena['mstar']),np.array(objs_pd_elena['SFR'])/np.array(objs_pd_elena['mstar']),marker = markers[7], label = 'Elena', zorder = 5,c = colors_elena, cmap = colormap, norm=cNorm)
plt.xlabel(r'M$_{HI}$/M$_*$')
plt.ylabel(r'SSFR [yr$^{-1}$]')
plt.yscale('log')
plt.axis([0, 150, 3e-12, 3e-9])
#plt.legend([mcquinn_plt, cptmarvel_plt, rogue_plt, elektra_plt, storm_plt],['McQuinn+ 2015','Cpt. Marvel','Rogue','Elektra','Storm'],loc=4)
plt.legend([mcquinn_plt, cptmarvel_plt, rogue_plt, elektra_plt, storm_plt, sandra_plt, ruth_plt, sonia_plt, elena_plt],['McQuinn+ 2015','Cpt. Marvel','Rogue','Elektra','Storm','Sandra','Ruth','Sonia','Elena'],loc=4)
plt.show()
plt.savefig(outbase + '.ssfr_gasfrac_all.png')
plt.clf()
"""

#Half stellar mass radius vs stellar mass
#
plt.figure(1,figsize=(plt_width,plt_width))
cptmarvel_plt = plt.scatter(np.array(objs_pd_cm['mstar']),np.array(objs_pd_cm['r_half']),marker = markers[0], label = 'Cpt. Marvel' ,zorder = 2,c = colors_cm, cmap = colormap, norm=cNorm)
rogue_plt = plt.scatter(np.array(objs_pd_r['mstar']),np.array(objs_pd_r['r_half']),marker = markers[1], label = 'Rogue', zorder = 3,c = colors_r, cmap = colormap, norm=cNorm)
elektra_plt = plt.scatter(np.array(objs_pd_e['mstar']),np.array(objs_pd_e['r_half']),marker = markers[2], label = 'Elektra',zorder = 4,c = colors_e, cmap = colormap, norm=cNorm)
storm_plt = plt.scatter(np.array(objs_pd_s['mstar']),np.array(objs_pd_s['r_half']),marker = markers[3], label = 'Storm', zorder = 5,c = colors_s, cmap = colormap, norm=cNorm)
sandra_plt = plt.scatter(np.array(objs_pd_sandra['mstar']),np.array(objs_pd_sandra['r_half']),marker = markers[4], label = 'Sandra', zorder = 5,c = colors_sandra, cmap = colormap, norm=cNorm)
ruth_plt = plt.scatter(np.array(objs_pd_ruth['mstar']),np.array(objs_pd_ruth['r_half']),marker = markers[5], label = 'Ruth', zorder = 5,c = colors_ruth, cmap = colormap, norm=cNorm)
sonia_plt = plt.scatter(np.array(objs_pd_sonia['mstar']),np.array(objs_pd_sonia['r_half']),marker = markers[6], label = 'Sonia', zorder = 5,c = colors_sonia, cmap = colormap, norm=cNorm)
elena_plt = plt.scatter(np.array(objs_pd_elena['mstar']),np.array(objs_pd_elena['r_half']),marker = markers[7], label = 'Elena', zorder = 5,c = colors_elena, cmap = colormap, norm=cNorm)
plt.xlabel(r'M$_*$/M$_\odot$')
plt.ylabel(r'r$_{[1/2]}$ [kpc]')
plt.xscale('log')
plt.yscale('log')
#plt.axis([2e3, 3e9, 3e-2, 6])
plt.axis([2e3, 1e12, 3e-2, 6])
#plt.legend(['Cpt. Marvel','Rogue','Elektra','Storm'],loc=4)
#plt.legend([cptmarvel_plt, rogue_plt, elektra_plt, storm_plt, sandra_plt, ruth_plt, sonia_plt, elena_plt],['Cpt. Marvel','Rogue','Elektra','Storm','Sandra','Ruth','Sonia','Elena'],loc=4)
plt.tight_layout() 
plt.show()
plt.savefig(outbase + '.size_mstar.png',dpi = dpi)
plt.clf()

#1-D stellar velocity dispersion vs stellar mass
#
plt.figure(1,figsize=(plt_width,plt_width))
plt.scatter(np.array(objs_pd_cm['mstar']),np.array(objs_pd_cm['sigma_v']),marker = markers[0], label = 'Cpt. Marvel',zorder = 2,c = colors_cm, cmap = colormap, norm=cNorm)
plt.scatter(np.array(objs_pd_r['mstar']),np.array(objs_pd_r['sigma_v']),marker = markers[1], label = 'Rogue', zorder = 3,c = colors_r, cmap = colormap, norm=cNorm)
plt.scatter(np.array(objs_pd_e['mstar']),np.array(objs_pd_e['sigma_v']),marker = markers[2], label = 'Elektra',zorder = 4,c = colors_e, cmap = colormap, norm=cNorm)
plt.scatter(np.array(objs_pd_s['mstar']),np.array(objs_pd_s['sigma_v']),marker = markers[3], label = 'Storm', zorder = 5,c = colors_s, cmap = colormap, norm=cNorm)
plt.scatter(np.array(objs_pd_sandra['mstar']),np.array(objs_pd_sandra['sigma_v']),marker = markers[4], label = 'Sandra', zorder = 5,c = colors_sandra, cmap = colormap, norm=cNorm)
plt.scatter(np.array(objs_pd_ruth['mstar']),np.array(objs_pd_ruth['sigma_v']),marker = markers[5], label = 'Ruth', zorder = 5,c = colors_ruth, cmap = colormap, norm=cNorm)
plt.scatter(np.array(objs_pd_sonia['mstar']),np.array(objs_pd_sonia['sigma_v']),marker = markers[6], label = 'Sonia', zorder = 5,c = colors_sonia, cmap = colormap, norm=cNorm)
plt.scatter(np.array(objs_pd_elena['mstar']),np.array(objs_pd_elena['sigma_v']),marker = markers[7], label = 'Elena', zorder = 5,c = colors_elena, cmap = colormap, norm=cNorm)
plt.xlabel(r'M$_*$/M$_\odot$')
plt.ylabel(r'$\sigma_*$ [km/s]')
plt.xscale('log')
plt.yscale('log')
#plt.legend(['Cpt. Marvel','Rogue','Elektra','Storm'],loc=2)
plt.legend([cptmarvel_plt, rogue_plt, elektra_plt, storm_plt, sandra_plt, ruth_plt, sonia_plt, elena_plt],['Cpt. Marvel','Rogue','Elektra','Storm','Sandra','Ruth','Sonia','Elena'],loc=2)
plt.show()
plt.savefig(outbase + '.sigma_mstar.png',dpi = dpi)
plt.clf()

#Stellar Luminosities vs Metallicities
#Kirby+13 , Figure 8
magVsolar = 4.83
f = open(dataprefix+'Kirby2013.txt', 'r')
kirbydata = []
for line in f:
    line = line.strip()
    columns = line.split()
    if len(columns) == 17 and columns[0] != 'Galaxy':
        source = {}
        source['gal'] = columns[0]
        source['nstars'] = columns[1]
        source['LogL_V'] = float(columns[2])
        source['LogL_V_err'] = float(columns[3])
        source['LogM_star'] = columns[4]
        source['LogM_star_err'] = columns[5]
        source['LogFe'] = float(columns[6])
        source['LogFe_err'] = float(columns[7])
        source['sigma'] = columns[8]
        source['sigma2'] = columns[9]
        source['median'] = columns[10]
        source['mad'] = columns[11]
        source['iqr'] = columns[12]
        source['skewness'] = columns[13]
        source['skewness_err'] = columns[14]
        source['Kurtosis'] = columns[15]
        source['Kurtosis_err'] = columns[16]
        kirbydata.append(source)
f.close()
kirbydata = pd.DataFrame(kirbydata)
objs_pd_cm['LogL_V'] = ((magVsolar - objs_pd_cm['M_V'])/2.5)
objs_pd_r['LogL_V'] = ((magVsolar - objs_pd_r['M_V'])/2.5)
objs_pd_e['LogL_V'] = ((magVsolar - objs_pd_e['M_V'])/2.5)
objs_pd_s['LogL_V'] = ((magVsolar - objs_pd_s['M_V'])/2.5)
objs_pd_sandra['LogL_V'] = ((magVsolar - objs_pd_sandra['M_V'])/2.5)
objs_pd_ruth['LogL_V'] = ((magVsolar - objs_pd_ruth['M_V'])/2.5)
objs_pd_sonia['LogL_V'] = ((magVsolar - objs_pd_sonia['M_V'])/2.5)
objs_pd_elena['LogL_V'] = ((magVsolar - objs_pd_elena['M_V'])/2.5)
                               
plt.figure(1)
#kirbydata.plot(x = 'LogL_V',y = 'LogFe',yerr = 'LogFe_err',xerr = 'LogL_V_err',fmt="o",color = "grey")
kirby_plt = plt.errorbar((kirbydata['LogL_V']).tolist(),kirbydata['LogFe'].tolist(),yerr = kirbydata['LogFe_err'].tolist(),xerr = kirbydata['LogL_V_err'].tolist(),color = "grey",fmt="o",zorder = 1)
plt.scatter(np.array(objs_pd_cm['LogL_V'].values),np.array(objs_pd_cm['FeH_mean']),marker = markers[0], label = 'Cpt. Marvel',zorder = 2,edgecolors = colors_cm, cmap = colormap, norm=cNorm,facecolors = colors_cm_f)
plt.scatter(np.array(objs_pd_r['LogL_V']),np.array(objs_pd_r['FeH_mean']),marker = markers[1], label = 'Rogue', zorder = 3,edgecolors = colors_r, cmap = colormap, norm=cNorm,facecolors = colors_r_f)
plt.scatter(np.array(objs_pd_e['LogL_V']),np.array(objs_pd_e['FeH_mean']),marker = markers[2], label = 'Elektra',zorder = 4,edgecolors = colors_e, cmap = colormap, norm=cNorm,facecolors = colors_e_f)
plt.scatter(np.array(objs_pd_s['LogL_V']),np.array(objs_pd_s['FeH_mean']),marker = markers[3], label = 'Storm', zorder = 5,edgecolors = colors_s, cmap = colormap, norm=cNorm,facecolors = colors_s_f)
plt.scatter(np.array(objs_pd_sandra['LogL_V']),np.array(objs_pd_sandra['FeH_mean']),marker = markers[4], label = 'Sandra', zorder = 5,edgecolors = colors_sandra, cmap = colormap, norm=cNorm,facecolors = colors_sandra_f)
plt.scatter(np.array(objs_pd_ruth['LogL_V']),np.array(objs_pd_ruth['FeH_mean']),marker = markers[5], label = 'Ruth', zorder = 5,edgecolors = colors_ruth, cmap = colormap, norm=cNorm,facecolors = colors_ruth_f)
plt.scatter(np.array(objs_pd_sonia['LogL_V']),np.array(objs_pd_sonia['FeH_mean']),marker = markers[6], label = 'Sonia', zorder = 5,edgecolors = colors_sonia, cmap = colormap, norm=cNorm,facecolors = colors_sonia_f)
plt.scatter(np.array(objs_pd_elena['LogL_V']),np.array(objs_pd_elena['FeH_mean']),marker = markers[7], label = 'Elena', zorder = 5,edgecolors = colors_elena, cmap = colormap, norm=cNorm,facecolors = colors_elena_f)
plt.xlabel(r'L$_V$/L$_{V,solar}$')
plt.ylabel(r'<[Fe/H]>')
#plt.legend([kirby_plt, cptmarvel_plt, rogue_plt, elektra_plt, storm_plt],['Kirby+ 2013','Cpt. Marvel','Rogue','Elektra','Storm'],loc=2)
plt.legend([ cptmarvel_plt, rogue_plt, elektra_plt, storm_plt, sandra_plt, ruth_plt, sonia_plt, elena_plt],['Cpt. Marvel','Rogue','Elektra','Storm','Sandra','Ruth','Sonia','Elena'],loc=2)
#plt.axis([3.5, 9.5, -3.5, -0.5])
plt.axis([3.5, 12, -3.5, 0.1])
plt.show()
plt.savefig(outbase + '.FeH_Lv.png')
plt.clf()

plt.figure(1)
#kirbydata.plot(x = 'LogL_V',y = 'LogFe',yerr = 'LogFe_err',xerr = 'LogL_V_err',fmt="o",color = "grey")
kirby_plt = plt.errorbar((kirbydata['LogL_V']).tolist(),kirbydata['LogFe'].tolist(),yerr = kirbydata['LogFe_err'].tolist(),xerr = kirbydata['LogL_V_err'].tolist(),color = "grey",fmt="o",zorder = 1)
cptmarvel_plt = plt.scatter(np.array(objs_pd_cm['LogL_V'].values),np.array(objs_pd_cm['Zstellar_mean']),marker = markers[0], label = 'Cpt. Marvel',zorder = 2,edgecolors = colors_cm, cmap = colormap, norm=cNorm,facecolors = colors_cm_f)
rogue_plt = plt.scatter(np.array(objs_pd_r['LogL_V'].values),np.array(objs_pd_r['Zstellar_mean']),marker = markers[1], label = 'Rogue', zorder = 3,edgecolors = colors_r, cmap = colormap, norm=cNorm,facecolors = colors_r_f)
elektra_plt = plt.scatter(np.array(objs_pd_e['LogL_V'].values),np.array(objs_pd_e['Zstellar_mean']),marker = markers[2], label = 'Elektra',zorder = 4,edgecolors = colors_e, cmap = colormap, norm=cNorm,facecolors = colors_e_f)
storm_plt = plt.scatter(np.array(objs_pd_s['LogL_V'].values),np.array(objs_pd_s['Zstellar_mean']),marker = markers[3], label = 'Storm', zorder = 5,edgecolors = colors_s, cmap = colormap, norm=cNorm,facecolors = colors_s_f)
#sandra_plt = plt.scatter(np.array(objs_pd_sandra['LogL_V'].values),np.array(objs_pd_sandra['Zstellar_mean']),marker = markers[4], label = 'Sandra', zorder = 5,edgecolors = colors_sandra, cmap = colormap, norm=cNorm,facecolors = colors_sandra_f)
#plt.scatter(np.array(objs_pd_ruth['LogL_V']),np.array(objs_pd_ruth['Zstellar_mean']),marker = markers[5], label = 'Ruth', zorder = 5,edgecolors = colors_ruth, cmap = colormap, norm=cNorm,facecolors = colors_ruth_f)
#plt.scatter(np.array(objs_pd_sonia['LogL_V']),np.array(objs_pd_sonia['Zstellar_mean']),marker = markers[6], label = 'Sonia', zorder = 5,edgecolors = colors_sonia, cmap = colormap, norm=cNorm,facecolors = colors_sonia_f)
#plt.scatter(np.array(objs_pd_elena['LogL_V']),np.array(objs_pd_elena['Zstellar_mean']),marker = markers[7], label = 'Elena', zorder = 5,edgecolors = colors_elena, cmap = colormap, norm=cNorm,facecolors = colors_elena_f)
plt.xlabel(r'L$_V$/L$_{V,solar}$')
plt.ylabel(r'<Z/Z$_\odot$>')
#plt.legend([kirby_plt, cptmarvel_plt, rogue_plt, elektra_plt, storm_plt],['Kirby+ 2013','Cpt. Marvel','Rogue','Elektra','Storm'],loc=2)
plt.legend([kirby_plt,cptmarvel_plt],['Kirby+ 2013','Simulations'],loc=4)
#plt.legend([kirby_plt, cptmarvel_plt, rogue_plt, elektra_plt, storm_plt, sandra_plt, ruth_plt, sonia_plt, elena_plt],['McQuinn+ 2015','Cpt. Marvel','Rogue','Elektra','Storm','Sandra','Ruth','Sonia','Elena'],loc=2)
plt.axis([4.3, 9.5, -3.0, -0.5])
#plt.axis([3.5, 11.5, -3.0, 0.5])
plt.tight_layout()
plt.show()
plt.savefig(outbase + '.Z_Lv.png')
plt.clf()

plt.figure(1)
plt.scatter(np.array(objs_pd_cm['mstar']),np.array(objs_pd_cm['FeH_mean']),marker = markers[0], label = 'Cpt. Marvel',zorder = 2,c = colors_cm, cmap = colormap, norm=cNorm)
plt.scatter(np.array(objs_pd_r['mstar']),np.array(objs_pd_r['FeH_mean']),marker = markers[1], label = 'Rogue', zorder = 3,c = colors_r, cmap = colormap, norm=cNorm)
plt.scatter(np.array(objs_pd_e['mstar']),np.array(objs_pd_e['FeH_mean']),marker = markers[2], label = 'Elektra',zorder = 4,c = colors_e, cmap = colormap, norm=cNorm)
plt.scatter(np.array(objs_pd_s['mstar']),np.array(objs_pd_s['FeH_mean']),marker = markers[3], label = 'Storm', zorder = 5,c = colors_s, cmap = colormap, norm=cNorm)
plt.scatter(np.array(objs_pd_sandra['mstar']),np.array(objs_pd_sandra['FeH_mean']),marker = markers[4], label = 'Sandra', zorder = 5,c = colors_sandra, cmap = colormap, norm=cNorm)
plt.scatter(np.array(objs_pd_ruth['mstar']),np.array(objs_pd_ruth['FeH_mean']),marker = markers[5], label = 'Ruth', zorder = 5,c = colors_ruth, cmap = colormap, norm=cNorm)
plt.scatter(np.array(objs_pd_sonia['mstar']),np.array(objs_pd_sonia['FeH_mean']),marker = markers[6], label = 'Sonia', zorder = 5,c = colors_sonia, cmap = colormap, norm=cNorm)
plt.scatter(np.array(objs_pd_elena['mstar']),np.array(objs_pd_elena['FeH_mean']),marker = markers[7], label = 'Elena', zorder = 5,c = colors_elena, cmap = colormap, norm=cNorm)
plt.xlabel(r'M$_*$/M$_\odot$')
plt.ylabel(r'<[Fe/H]>')
plt.xscale('log')
#plt.legend(['Cpt. Marvel','Rogue','Elektra','Storm'],loc=2)
plt.legend([cptmarvel_plt, rogue_plt, elektra_plt, storm_plt, sandra_plt, ruth_plt, sonia_plt, elena_plt],['Cpt. Marvel','Rogue','Elektra','Storm','Sandra','Ruth','Sonia','Elena'],loc=2)
#plt.axis([2e3, 3e9, -3.5, -0.5])
plt.axis([2e3, 1e12, -3.5, 0.1])
plt.show()
plt.savefig(outbase + '.FeH_mstar.png')
plt.clf()

#Stellar Luminosities vs Spread in Metallicities
#Kirby+ , Figure 4
plt.figure(1)
plt.scatter(np.array(objs_pd_cm['LogL_V']),np.array(objs_pd_cm['FeH_std']),marker = markers[0], label = 'Cpt. Marvel',zorder = 2,c = colors_cm, cmap = colormap, norm=cNorm)
plt.scatter(np.array(objs_pd_r['LogL_V']),np.array(objs_pd_r['FeH_std']),marker = markers[1], label = 'Rogue', zorder = 3,c = colors_r, cmap = colormap, norm=cNorm)
plt.scatter(np.array(objs_pd_e['LogL_V']),np.array(objs_pd_e['FeH_std']),marker = markers[2], label = 'Elektra',zorder = 4,c = colors_e, cmap = colormap, norm=cNorm)
plt.scatter(np.array(objs_pd_s['LogL_V']),np.array(objs_pd_s['FeH_std']),marker = markers[3], label = 'Storm', zorder = 5,c = colors_s, cmap = colormap, norm=cNorm)
plt.scatter(np.array(objs_pd_sandra['LogL_V']),np.array(objs_pd_sandra['FeH_std']),marker = markers[4], label = 'Sandra', zorder = 5,c = colors_sandra, cmap = colormap, norm=cNorm)
plt.scatter(np.array(objs_pd_ruth['LogL_V']),np.array(objs_pd_ruth['FeH_std']),marker = markers[5], label = 'Ruth', zorder = 5,c = colors_ruth, cmap = colormap, norm=cNorm)
plt.scatter(np.array(objs_pd_sonia['LogL_V']),np.array(objs_pd_sonia['FeH_std']),marker = markers[6], label = 'Sonia', zorder = 5,c = colors_sonia, cmap = colormap, norm=cNorm)
plt.scatter(np.array(objs_pd_elena['LogL_V']),np.array(objs_pd_elena['FeH_std']),marker = markers[7], label = 'Elena', zorder = 5,c = colors_elena, cmap = colormap, norm=cNorm)
plt.xlabel(r'L$_V$/L$_{V,solar}$')
plt.ylabel(r'$\sigma_{[Fe/H]}$')
#plt.legend(['Cpt. Marvel','Rogue','Elektra','Storm'],loc=1)
plt.legend([ cptmarvel_plt, rogue_plt, elektra_plt, storm_plt, sandra_plt, ruth_plt, sonia_plt, elena_plt],['Cpt. Marvel','Rogue','Elektra','Storm','Sandra','Ruth','Sonia','Elena'],loc=1)
plt.show()
plt.savefig(outbase + '.sigmaFeH.png')
plt.clf()

#Gas-phase mass metallicity relationship
XSOLO = 0.84E-2 #What pynbody uses
XSOLH = 0.706
f = open(dataprefix+'Lee06.txt', 'r')
leedata = []
for line in f:
    line = line.strip()
    columns = line.split()
    if len(columns) == 15 and columns[0] != 'Galaxy':
        source = {}
        source['gal'] = columns[0]
        source['source'] = columns[1]
        source['aorkey'] = columns[2]
        source['F'] = float(columns[3])
        source['F_err'] = float(columns[4])
        source['mag'] = float(columns[5])
        source['magerr'] = float(columns[6])
        source['disref'] = columns[7]
        source['mag_ab'] = float(columns[8])
        source['mag_ab_err'] = float(columns[9])
        source['O_H'] = float(columns[10])
        source['O_H_err'] = float(columns[11])
        source['O_Href'] = columns[12]
        source['B'] = float(columns[13])
        source['logMstar'] = float(columns[14])
        leedata.append(source)
f.close()
leedata = pd.DataFrame(leedata)

f = open(dataprefix+'Tremonti04.txt', 'r')
tremontidata = []
for line in f:
    line = line.strip()
    columns = line.split()
    if len(columns) == 6 and columns[0] != 'LogMstar':
        source = {}
        source['logMstar'] = float(columns[0])
        source['P2_5'] = float(columns[1])
        source['P16'] = float(columns[2])
        source['P50'] = float(columns[3])
        source['P84'] = float(columns[4])
        source['P97_5'] = float(columns[5])
        tremontidata.append(source)
f.close()
tremontidata = pd.DataFrame(tremontidata)
a = 319.7570 
b = -107.13160 
c = 12.208670 
d = -0.4606539 
        
plt.figure(1,figsize=(plt_width,plt_width*aspect_ratio))
lee_plt = plt.errorbar(10**leedata['logMstar'],leedata['O_H'],yerr=leedata['O_H_err'],color = "grey",fmt="o",zorder = 1)
plt.plot(10**tremontidata['logMstar'],a + b*tremontidata['P2_5'] + c*tremontidata['P2_5']**2 + d*tremontidata['P2_5']**3,color = "grey",linestyle = ":")
plt.plot(10**tremontidata['logMstar'],a + b*tremontidata['P97_5'] + c*tremontidata['P97_5']**2 + d*tremontidata['P97_5']**3,color = "grey",linestyle = ":")
plt.plot(10**tremontidata['logMstar'],a + b*tremontidata['P16'] + c*tremontidata['P16']**2 + d*tremontidata['P16']**3,color = "grey",linestyle = "--")
plt.plot(10**tremontidata['logMstar'],a + b*tremontidata['P84'] + c*tremontidata['P84']**2 + d*tremontidata['P84']**3,color = "grey",linestyle = "--")
tremonti_plt, = plt.plot(10**tremontidata['logMstar'],a + b*tremontidata['P50'] + c*tremontidata['P50']**2 + d*tremontidata['P50']**3,color = "grey",linestyle = "-")
cptmarvel_plt = plt.scatter(objs_pd_cm['mstar'],objs_pd_cm['oxh_sfr'] + 12,marker = markers[0], label = 'Cpt. Marvel',zorder = 2,edgecolors = colors_cm, cmap = colormap, norm=cNorm,facecolors = colors_cm_f,linewidth = 2)
plt.scatter(objs_pd_r['mstar'],objs_pd_r['oxh_sfr'] + 12,marker = markers[1], label = 'Rogue', zorder = 3,edgecolors = colors_r, cmap = colormap, norm=cNorm,facecolors = colors_r_f,linewidth = 2)
plt.scatter(objs_pd_e['mstar'],objs_pd_e['oxh_sfr'] + 12,marker = markers[2], label = 'Elektra',zorder = 4,edgecolors = colors_e, cmap = colormap, norm=cNorm,facecolors = colors_e_f,linewidth = 2)
plt.scatter(objs_pd_s['mstar'],objs_pd_s['oxh_sfr'] + 12,marker = markers[3], label = 'Storm', zorder = 5,edgecolors = colors_s, cmap = colormap, norm=cNorm,facecolors = colors_s_f,linewidth = 2)
#plt.scatter(np.array(objs_pd_sandra['mstar'])[np.where(objs_pd_sandra['oxh_sfr'] != 0)],np.array(objs_pd_sandra['oxh_sfr'] + 12)[np.where(objs_pd_sandra['oxh_sfr'] != 0)],marker = markers[4], label = 'Sandra', zorder = 5,c = colors_s[np.where(objs_pd_sandra['oxh_sfr'] != 0)], cmap = colormap, norm=cNorm)
#plt.scatter(np.array(objs_pd_ruth['mstar'])[np.where(objs_pd_ruth['oxh_sfr'] != 0)],np.array(objs_pd_ruth['oxh_sfr'] + 12)[np.where(objs_pd_ruth['oxh_sfr'] != 0)],marker = markers[5], label = 'Ruth', zorder = 5,c = colors_s[np.where(objs_pd_ruth['oxh_sfr'] != 0)], cmap = colormap, norm=cNorm)
#plt.scatter(np.array(objs_pd_sonia['mstar'])[np.where(objs_pd_sonia['oxh_sfr'] != 0)],np.array(objs_pd_sonia['oxh_sfr'] + 12)[np.where(objs_pd_sonia['oxh_sfr'] != 0)],marker = markers[6], label = 'Sonia', zorder = 5,c = colors_s[np.where(objs_pd_sonia['oxh_sfr'] != 0)], cmap = colormap, norm=cNorm)
#plt.scatter(np.array(objs_pd_elena['mstar'])[np.where(objs_pd_elena['oxh_sfr'] != 0)],np.array(objs_pd_elena['oxh_sfr'] + 12)[np.where(objs_pd_elena['oxh_sfr'] != 0)],marker = markers[7], label = 'Elena', zorder = 5,c = colors_s[np.where(objs_pd_elena['oxh_sfr'] != 0)], cmap = colormap, norm=cNorm)
#plt.legend([mcquinn_plt, cptmarvel_plt, rogue_plt, elektra_plt, storm_plt, sandra_plt, ruth_plt, sonia_plt, elena_plt],['McQuinn+ 2015','Cpt. Marvel','Rogue','Elektra','Storm','Sandra','Ruth','Sonia','Elena'],loc=2)
plt.legend([tremonti_plt,lee_plt,cptmarvel_plt],['Tremonti et al. 2004','Lee et al. 2006','Simulations'],loc = 4)
plt.xlabel(r'M$_*$/M$_\odot$')
plt.ylabel(r'12 + log(O/H)')
plt.xscale('log')
#plt.axis([1e3, 1e12, 6.0, 9.5])
plt.axis([5e3, 5e9, 6.0, 9.0])
plt.tight_layout()
plt.show()
plt.savefig(outbase + '.mzr_sfe.png',dpi = dpi)
plt.clf()

plt.figure(1)
plt.scatter(objs_pd_cm[objs_pd_cm['oxh_cold'] != 0]['mstar'],objs_pd_cm[objs_pd_cm['oxh_cold'] != 0]['oxh_cold'] + 12,marker = markers[0], label = 'Cpt. Marvel',zorder = 2,c = colors_cm[np.where(objs_pd_cm['oxh_cold'] != 0)[0]], cmap = colormap, norm=cNorm)
plt.scatter(np.array(objs_pd_r['mstar'])[np.where(objs_pd_r['oxh_cold'] != 0)],np.array(objs_pd_r['oxh_cold'] + 12)[np.where(objs_pd_r['oxh_cold'] != 0)],marker = markers[1], label = 'Rogue', zorder = 3,c = colors_r[np.where(objs_pd_r['oxh_cold'] != 0)], cmap = colormap, norm=cNorm)
plt.scatter(np.array(objs_pd_e['mstar'])[np.where(objs_pd_e['oxh_cold'] != 0)],np.array(objs_pd_e['oxh_cold'] + 12)[np.where(objs_pd_e['oxh_cold'] != 0)],marker = markers[2], label = 'Elektra',zorder = 4,c = colors_e[np.where(objs_pd_e['oxh_cold'] != 0)], cmap = colormap, norm=cNorm)
plt.scatter(np.array(objs_pd_s['mstar'])[np.where(objs_pd_s['oxh_cold'] != 0)],np.array(objs_pd_s['oxh_cold'] + 12)[np.where(objs_pd_s['oxh_cold'] != 0)],marker = markers[3], label = 'Storm', zorder = 5,c = colors_s[np.where(objs_pd_s['oxh_cold'] != 0)], cmap = colormap, norm=cNorm)
#plt.scatter(np.array(objs_pd_sandra['mstar'])[np.where(objs_pd_sandra['oxh_cold'] != 0)],np.array(objs_pd_sandra['oxh_cold'] + 12)[np.where(objs_pd_sandra['oxh_cold'] != 0)],marker = markers[4], label = 'Sandra', zorder = 5,c = colors_s[np.where(objs_pd_sandra['oxh_cold'] != 0)], cmap = colormap, norm=cNorm)
#plt.scatter(np.array(objs_pd_ruth['mstar'])[np.where(objs_pd_ruth['oxh_cold'] != 0)],np.array(objs_pd_ruth['oxh_cold'] + 12)[np.where(objs_pd_ruth['oxh_cold'] != 0)],marker = markers[5], label = 'Ruth', zorder = 5,c = colors_s[np.where(objs_pd_ruth['oxh_cold'] != 0)], cmap = colormap, norm=cNorm)
#plt.scatter(np.array(objs_pd_sonia['mstar'])[np.where(objs_pd_sonia['oxh_cold'] != 0)],np.array(objs_pd_sonia['oxh_cold'] + 12)[np.where(objs_pd_sonia['oxh_cold'] != 0)],marker = markers[6], label = 'Sonia', zorder = 5,c = colors_s[np.where(objs_pd_sonia['oxh_cold'] != 0)], cmap = colormap, norm=cNorm)
#plt.scatter(np.array(objs_pd_elena['mstar'])[np.where(objs_pd_elena['oxh_cold'] != 0)],np.array(objs_pd_elena['oxh_cold'] + 12)[np.where(objs_pd_elena['oxh_cold'] != 0)],marker = markers[7], label = 'Elena', zorder = 5,c = colors_s[np.where(objs_pd_elena['oxh_cold'] != 0)], cmap = colormap, norm=cNorm)
plt.xlabel(r'M$_*$/M$_\odot$')
plt.ylabel(r'12 + log(O/H)')
plt.xscale('log')
plt.legend([tremonti_plt,lee_plt,cptmarvel_plt],['Tremonti et al. 2004','Lee et al. 2006','Simulations'],loc = 4)
#plt.legend(['Cpt. Marvel','Rogue','Elektra','Storm'],loc=2)
#plt.legend([mcquinn_plt, cptmarvel_plt, rogue_plt, elektra_plt, storm_plt, sandra_plt, ruth_plt, sonia_plt, elena_plt],['McQuinn+ 2015','Cpt. Marvel','Rogue','Elektra','Storm','Sandra','Ruth','Sonia','Elena'],loc=2)
#plt.axis([7e5, 3e9, 6.5, 8.2])
plt.axis([3e5, 1e12, 7, 9.5])
plt.show()
plt.savefig(outbase + '.mzr_cold.png')
plt.clf()

