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

outbase = outprefix + 'marvel'

f_bar = 0.16510
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

plt.figure(1)
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
plt.show()
plt.savefig(outbase + '.HI_mgas_vMstar.png')
plt.clf()

plt.figure(1)
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
plt.show()
plt.savefig(outbase + '.HI_mgas_vSSFR.png')
plt.clf()

plt.figure(1)
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
plt.show()
plt.savefig(outbase + '.HI_mgas_vMvir.png')
plt.clf()

plt.figure(1)
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
plt.show()
plt.savefig(outbase + '.HI_mvir_vMstar.png')
plt.clf()

plt.figure(1)
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
plt.show()
plt.savefig(outbase + '.HI_mvir_vSSFR.png')
plt.clf()

plt.figure(1)
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
plt.show()
plt.savefig(outbase + '.HI_mvir_vMvir.png')
plt.clf()

plt.figure(1)
cptmarvel_plt = plt.scatter(np.array(objs_pd_cm['mstar']),np.array(objs_pd_cm['mgas'])/np.array(objs_pd_cm['mvir']),marker = markers[0], label = 'Cpt. Marvel' ,zorder = 2,c = colors_cm, cmap = colormap, norm=cNorm)
rogue_plt = plt.scatter(np.array(objs_pd_r['mstar']),np.array(objs_pd_r['mgas'])/np.array(objs_pd_r['mvir']),marker = markers[1], label = 'Rogue', zorder = 3,c = colors_r, cmap = colormap, norm=cNorm)
elektra_plt = plt.scatter(np.array(objs_pd_e['mstar']),np.array(objs_pd_e['mgas'])/np.array(objs_pd_e['mvir']),marker = markers[2], label = 'Elektra',zorder = 4,c = colors_e, cmap = colormap, norm=cNorm)
storm_plt = plt.scatter(np.array(objs_pd_s['mstar']),np.array(objs_pd_s['mgas'])/np.array(objs_pd_s['mvir']),marker = markers[3], label = 'Storm', zorder = 5,c = colors_s, cmap = colormap, norm=cNorm)
plt.xlabel(r'M$_*$/M$_\odot$')
plt.ylabel(r'M$_{gas}$/M$_{vir}$')
plt.xscale('log')
plt.axis([1e3, 1e10, 0, 0.1])
plt.show()
plt.savefig(outbase + '.mgas_mvir_vMstar.png')
plt.clf()

plt.figure(1)
cptmarvel_plt = plt.scatter(np.array(objs_pd_cm['SFR'])/np.array(objs_pd_cm['mstar']),np.array(objs_pd_cm['mgas'])/np.array(objs_pd_cm['mvir']),marker = markers[0], label = 'Cpt. Marvel' ,zorder = 2,c = colors_cm, cmap = colormap, norm=cNorm)
rogue_plt = plt.scatter(np.array(objs_pd_r['SFR'])/np.array(objs_pd_r['mstar']),np.array(objs_pd_r['mgas'])/np.array(objs_pd_r['mvir']),marker = markers[1], label = 'Rogue', zorder = 3,c = colors_r, cmap = colormap, norm=cNorm)
elektra_plt = plt.scatter(np.array(objs_pd_e['SFR'])/np.array(objs_pd_e['mstar']),np.array(objs_pd_e['mgas'])/np.array(objs_pd_e['mvir']),marker = markers[2], label = 'Elektra',zorder = 4,c = colors_e, cmap = colormap, norm=cNorm)
storm_plt = plt.scatter(np.array(objs_pd_s['SFR'])/np.array(objs_pd_s['mstar']),np.array(objs_pd_s['mgas'])/np.array(objs_pd_s['mvir']),marker = markers[3], label = 'Storm', zorder = 5,c = colors_s, cmap = colormap, norm=cNorm)
plt.xlabel(r'SSFR [yr$^{-1}$]')
plt.ylabel(r'M$_{gas}$/M$_{vir}$')
plt.xscale('log')
plt.axis([6e-12, 1e-9, 0, 0.1])
plt.show()
plt.savefig(outbase + '.mgas_mvir_vSSFR.png')
plt.clf()

plt.figure(1)
cptmarvel_plt = plt.scatter(np.array(objs_pd_cm['mvir']),np.array(objs_pd_cm['mgas'])/np.array(objs_pd_cm['mvir']),marker = markers[0], label = 'Cpt. Marvel' ,zorder = 2,c = colors_cm, cmap = colormap, norm=cNorm)
rogue_plt = plt.scatter(np.array(objs_pd_r['mvir']),np.array(objs_pd_r['mgas'])/np.array(objs_pd_r['mvir']),marker = markers[1], label = 'Rogue', zorder = 3,c = colors_r, cmap = colormap, norm=cNorm)
elektra_plt = plt.scatter(np.array(objs_pd_e['mvir']),np.array(objs_pd_e['mgas'])/np.array(objs_pd_e['mvir']),marker = markers[2], label = 'Elektra',zorder = 4,c = colors_e, cmap = colormap, norm=cNorm)
storm_plt = plt.scatter(np.array(objs_pd_s['mvir']),np.array(objs_pd_s['mgas'])/np.array(objs_pd_s['mvir']),marker = markers[3], label = 'Storm', zorder = 5,c = colors_s, cmap = colormap, norm=cNorm)
plt.xlabel(r'M$_{vir}$/M$_\odot$')
plt.ylabel(r'M$_{gas}$/M$_{vir}$')
plt.xscale('log')
plt.axis([8e8, 1e11, 0, 0.1])
plt.show()
plt.savefig(outbase + '.mgas_mvir_vMvir.png')
plt.clf()


plt.figure(1)
cptmarvel_plt = plt.scatter(objs_pd_cm['mstar'],(objs_pd_cm['mgas'] - objs_pd_cm['mHI'])/objs_pd_cm['mvir'],marker = markers[0], label = 'Cpt. Marvel' ,zorder = 2,c = colors_cm, cmap = colormap, norm=cNorm)
rogue_plt = plt.scatter(objs_pd_r['mstar'],(objs_pd_r['mgas'] - objs_pd_r['mHI'])/objs_pd_r['mvir'],marker = markers[1], label = 'Rogue' ,zorder = 2,c = colors_r, cmap = colormap, norm=cNorm)
elektra_plt = plt.scatter(objs_pd_e['mstar'],(objs_pd_e['mgas'] - objs_pd_e['mHI'])/objs_pd_e['mvir'],marker = markers[2], label = 'Elektra' ,zorder = 2,c = colors_e, cmap = colormap, norm=cNorm)
storm_plt = plt.scatter(objs_pd_s['mstar'],(objs_pd_s['mgas'] - objs_pd_s['mHI'])/objs_pd_s['mvir'],marker = markers[3], label = 'Storm' ,zorder = 2,c = colors_s, cmap = colormap, norm=cNorm)
plt.xlabel(r'M$_*$/M$_\odot$')
plt.ylabel(r'M$_{halo}$/M$_{vir}$')
plt.xscale('log')
plt.axis([1e3, 1e10, 0, 0.1])
plt.show()
plt.savefig(outbase + '.mhgas_mvir_vMstar.png')
plt.clf()

plt.figure(1)
cptmarvel_plt = plt.scatter(objs_pd_cm['SFR']/objs_pd_cm['mstar'],(objs_pd_cm['mgas'] - objs_pd_cm['mHI'])/objs_pd_cm['mvir'],marker = markers[0], label = 'Cpt. Marvel' ,zorder = 2,c = colors_cm, cmap = colormap, norm=cNorm)
rogue_plt = plt.scatter(objs_pd_r['SFR']/objs_pd_r['mstar'],(objs_pd_r['mgas'] - objs_pd_r['mHI'])/objs_pd_r['mvir'],marker = markers[1], label = 'Rogue' ,zorder = 2,c = colors_r, cmap = colormap, norm=cNorm)
elektra_plt = plt.scatter(objs_pd_e['SFR']/objs_pd_e['mstar'],(objs_pd_e['mgas'] - objs_pd_e['mHI'])/objs_pd_e['mvir'],marker = markers[2], label = 'Elektra' ,zorder = 2,c = colors_e, cmap = colormap, norm=cNorm)
storm_plt = plt.scatter(objs_pd_s['SFR']/objs_pd_s['mstar'],(objs_pd_s['mgas'] - objs_pd_s['mHI'])/objs_pd_s['mvir'],marker = markers[3], label = 'Storm' ,zorder = 2,c = colors_s, cmap = colormap, norm=cNorm)
plt.xlabel(r'SSFR [yr$^{-1}$]')
plt.ylabel(r'M$_{halo}$/M$_{vir}$')
plt.xscale('log')
plt.axis([6e-12, 1e-9, 0, 0.1])
plt.show()
plt.savefig(outbase + '.mhgas_mvir_vSSFR.png')
plt.clf()

plt.figure(1)
cptmarvel_plt = plt.scatter(objs_pd_cm['mvir'],(objs_pd_cm['mgas'] - objs_pd_cm['mHI'])/objs_pd_cm['mvir'],marker = markers[0], label = 'Cpt. Marvel' ,zorder = 2,c = colors_cm, cmap = colormap, norm=cNorm)
rogue_plt = plt.scatter(objs_pd_r['mvir'],(objs_pd_r['mgas'] - objs_pd_r['mHI'])/objs_pd_r['mvir'],marker = markers[1], label = 'Rogue' ,zorder = 2,c = colors_r, cmap = colormap, norm=cNorm)
elektra_plt = plt.scatter(objs_pd_e['mvir'],(objs_pd_e['mgas'] - objs_pd_e['mHI'])/objs_pd_e['mvir'],marker = markers[2], label = 'Elektra' ,zorder = 2,c = colors_e, cmap = colormap, norm=cNorm)
storm_plt = plt.scatter(objs_pd_s['mvir'],(objs_pd_s['mgas'] - objs_pd_s['mHI'])/objs_pd_s['mvir'],marker = markers[3], label = 'Storm' ,zorder = 2,c = colors_s, cmap = colormap, norm=cNorm)
plt.xlabel(r'M$_{vir}$/M$_\odot$')
plt.ylabel(r'M$_{gas}$/M$_{vir}$')
plt.xscale('log')
plt.axis([8e8, 1e11, 0, 0.1])
plt.show()
plt.savefig(outbase + '.mhgas_mvir_vMvir.png')
plt.clf()

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

plt.figure(1)
hot_plt = plt.bar(cm_mstar_log,(objs_pd_cm['mgas'])/objs_pd_cm['mvir']/f_bar,width=0.1,color = 'k')
plt.bar(e_mstar_log,(objs_pd_e['mgas'])/objs_pd_e['mvir']/f_bar,width=0.1,color = 'k')
plt.bar(r_mstar_log,(objs_pd_r['mgas'])/objs_pd_r['mvir']/f_bar,width=0.1,color = 'k')
plt.bar(s_mstar_log,(objs_pd_s['mgas'])/objs_pd_s['mvir']/f_bar,width=0.1,color = 'k')
warm_plt = plt.bar(cm_mstar_log,(objs_pd_cm['mwarm'] + objs_pd_cm['mCool'] + objs_pd_cm['mISM'] + objs_pd_cm['mstar'])/objs_pd_cm['mvir']/f_bar,width=0.1,color= 'r')
plt.bar(e_mstar_log,(objs_pd_e['mwarm'] + objs_pd_e['mCool'] + objs_pd_e['mISM'] + objs_pd_e['mstar'])/objs_pd_e['mvir']/f_bar,width=0.1,color= 'r')
plt.bar(r_mstar_log,(objs_pd_r['mwarm'] + objs_pd_r['mCool'] + objs_pd_r['mISM'] + objs_pd_r['mstar'])/objs_pd_r['mvir']/f_bar,width=0.1,color= 'r')
plt.bar(s_mstar_log,(objs_pd_s['mwarm'] + objs_pd_s['mCool'] + objs_pd_s['mISM']+ objs_pd_s['mstar'])/objs_pd_s['mvir']/f_bar,width=0.1,color= 'r')
cool_plt = plt.bar(cm_mstar_log,(objs_pd_cm['mCool'] + objs_pd_cm['mISM'] + objs_pd_cm['mstar'])/objs_pd_cm['mvir']/f_bar,width=0.1,color= 'g')
plt.bar(e_mstar_log,(objs_pd_e['mCool'] + objs_pd_e['mISM'] + objs_pd_e['mstar'])/objs_pd_e['mvir']/f_bar,width=0.1,color= 'g')
plt.bar(r_mstar_log,(objs_pd_r['mCool'] + objs_pd_r['mISM'] + objs_pd_r['mstar'])/objs_pd_r['mvir']/f_bar,width=0.1,color= 'g')
plt.bar(s_mstar_log,(objs_pd_s['mCool'] + objs_pd_s['mISM'] + objs_pd_s['mstar'])/objs_pd_s['mvir']/f_bar,width=0.1,color= 'g')
ISM_plt = plt.bar(cm_mstar_log,(objs_pd_cm['mISM'] + objs_pd_cm['mstar'])/objs_pd_cm['mvir']/f_bar,width=0.1,color= 'b')
plt.bar(e_mstar_log,(objs_pd_e['mISM'] + objs_pd_e['mstar'])/objs_pd_e['mvir']/f_bar,width=0.1,color= 'b')
plt.bar(r_mstar_log,(objs_pd_r['mISM'] + objs_pd_r['mstar'])/objs_pd_r['mvir']/f_bar,width=0.1,color= 'b')
plt.bar(s_mstar_log,(objs_pd_s['mISM'] + objs_pd_s['mstar'])/objs_pd_s['mvir']/f_bar,width=0.1,color= 'b')
plt.xlabel(r'Log(M$_*$/M$_\odot)$')
plt.ylabel(r'$M/M_{vir}/f_{bary}$')
plt.legend([hot_plt,warm_plt,cool_plt,ISM_plt],['$10^6 < $T','$10^5 <$ T $< 10^6$','$10^4 <$ T $< 10^5$','T $< 10^4$'],loc = 2)
plt.show()
plt.savefig(outbase + '.gascensus_vMstar.png')
plt.close()

plt.figure(1)
hot_plt = plt.bar(cm_mvir_log,(objs_pd_cm['mgas'])/objs_pd_cm['mvir']/f_bar,width=0.05,color = 'k')
plt.bar(e_mvir_log,(objs_pd_e['mgas'])/objs_pd_e['mvir']/f_bar,width=0.05,color = 'k')
plt.bar(r_mvir_log,(objs_pd_r['mgas'])/objs_pd_r['mvir']/f_bar,width=0.05,color = 'k')
plt.bar(s_mvir_log,(objs_pd_s['mgas'])/objs_pd_s['mvir']/f_bar,width=0.05,color = 'k')
warm_plt = plt.bar(cm_mvir_log,(objs_pd_cm['mwarm'] + objs_pd_cm['mCool'] + objs_pd_cm['mISM'] + objs_pd_cm['mstar'])/objs_pd_cm['mvir']/f_bar,width=0.05,color= 'r')
plt.bar(e_mvir_log,(objs_pd_e['mwarm'] + objs_pd_e['mCool'] + objs_pd_e['mISM'] + objs_pd_e['mstar'])/objs_pd_e['mvir']/f_bar,width=0.05,color= 'r')
plt.bar(r_mvir_log,(objs_pd_r['mwarm'] + objs_pd_r['mCool'] + objs_pd_r['mISM'] + objs_pd_r['mstar'])/objs_pd_r['mvir']/f_bar,width=0.05,color= 'r')
plt.bar(s_mvir_log,(objs_pd_s['mwarm'] + objs_pd_s['mCool'] + objs_pd_s['mISM']+ objs_pd_s['mstar'])/objs_pd_s['mvir']/f_bar,width=0.05,color= 'r')
cool_plt = plt.bar(cm_mvir_log,(objs_pd_cm['mCool'] + objs_pd_cm['mISM'] + objs_pd_cm['mstar'])/objs_pd_cm['mvir']/f_bar,width=0.05,color= 'g')
plt.bar(e_mvir_log,(objs_pd_e['mCool'] + objs_pd_e['mISM'] + objs_pd_e['mstar'])/objs_pd_e['mvir']/f_bar,width=0.05,color= 'g')
plt.bar(r_mvir_log,(objs_pd_r['mCool'] + objs_pd_r['mISM'] + objs_pd_r['mstar'])/objs_pd_r['mvir']/f_bar,width=0.05,color= 'g')
plt.bar(s_mvir_log,(objs_pd_s['mCool'] + objs_pd_s['mISM'] + objs_pd_s['mstar'])/objs_pd_s['mvir']/f_bar,width=0.05,color= 'g')
ISM_plt = plt.bar(cm_mvir_log,(objs_pd_cm['mISM'] + objs_pd_cm['mstar'])/objs_pd_cm['mvir']/f_bar,width=0.05,color= 'b')
plt.bar(e_mvir_log,(objs_pd_e['mISM'] + objs_pd_e['mstar'])/objs_pd_e['mvir']/f_bar,width=0.05,color= 'b')
plt.bar(r_mvir_log,(objs_pd_r['mISM'] + objs_pd_r['mstar'])/objs_pd_r['mvir']/f_bar,width=0.05,color= 'b')
plt.bar(s_mvir_log,(objs_pd_s['mISM'] + objs_pd_s['mstar'])/objs_pd_s['mvir']/f_bar,width=0.05,color= 'b')
plt.xlabel(r'Log(M$_{vir}$/M$_\odot)$')
plt.ylabel(r'$M/M_{vir}/f_{bary}$')
plt.legend([hot_plt,warm_plt,cool_plt,ISM_plt],['$10^6 < $T','$10^5 <$ T $< 10^6$','$10^4 <$ T $< 10^5$','T $< 10^4$'],loc = 2)
plt.show()
plt.savefig(outbase + '.gascensus_vMvir.png')
plt.close()

plt.figure(1)
hot_plt = plt.bar(cm_SSFR,(objs_pd_cm['mgas'])/objs_pd_cm['mvir']/f_bar,width=1e-11,color = 'k')
plt.bar(e_SSFR,(objs_pd_e['mgas'])/objs_pd_e['mvir']/f_bar,width=1e-11,color = 'k')
plt.bar(r_SSFR,(objs_pd_r['mgas'])/objs_pd_r['mvir']/f_bar,width=1e-11,color = 'k')
plt.bar(s_SSFR,(objs_pd_s['mgas'])/objs_pd_s['mvir']/f_bar,width=1e-11,color = 'k')
warm_plt = plt.bar(cm_SSFR,(objs_pd_cm['mwarm'] + objs_pd_cm['mCool'] + objs_pd_cm['mISM'] + objs_pd_cm['mstar'])/objs_pd_cm['mvir']/f_bar,width=1e-11,color= 'r')
plt.bar(e_SSFR,(objs_pd_e['mwarm'] + objs_pd_e['mCool'] + objs_pd_e['mISM'] + objs_pd_e['mstar'])/objs_pd_e['mvir']/f_bar,width=1e-11,color= 'r')
plt.bar(r_SSFR,(objs_pd_r['mwarm'] + objs_pd_r['mCool'] + objs_pd_r['mISM'] + objs_pd_r['mstar'])/objs_pd_r['mvir']/f_bar,width=1e-11,color= 'r')
plt.bar(s_SSFR,(objs_pd_s['mwarm'] + objs_pd_s['mCool'] + objs_pd_s['mISM']+ objs_pd_s['mstar'])/objs_pd_s['mvir']/f_bar,width=1e-11,color= 'r')
cool_plt = plt.bar(cm_SSFR,(objs_pd_cm['mCool'] + objs_pd_cm['mISM'] + objs_pd_cm['mstar'])/objs_pd_cm['mvir']/f_bar,width=1e-11,color= 'g')
plt.bar(e_SSFR,(objs_pd_e['mCool'] + objs_pd_e['mISM'] + objs_pd_e['mstar'])/objs_pd_e['mvir']/f_bar,width=1e-11,color= 'g')
plt.bar(r_SSFR,(objs_pd_r['mCool'] + objs_pd_r['mISM'] + objs_pd_r['mstar'])/objs_pd_r['mvir']/f_bar,width=1e-11,color= 'g')
plt.bar(s_SSFR,(objs_pd_s['mCool'] + objs_pd_s['mISM'] + objs_pd_s['mstar'])/objs_pd_s['mvir']/f_bar,width=1e-11,color= 'g')
ISM_plt = plt.bar(cm_SSFR,(objs_pd_cm['mISM'] + objs_pd_cm['mstar'])/objs_pd_cm['mvir']/f_bar,width=1e-11,color= 'b')
plt.bar(e_SSFR,(objs_pd_e['mISM'] + objs_pd_e['mstar'])/objs_pd_e['mvir']/f_bar,width=1e-11,color= 'b')
plt.bar(r_SSFR,(objs_pd_r['mISM'] + objs_pd_r['mstar'])/objs_pd_r['mvir']/f_bar,width=1e-11,color= 'b')
plt.bar(s_SSFR,(objs_pd_s['mISM'] + objs_pd_s['mstar'])/objs_pd_s['mvir']/f_bar,width=1e-11,color= 'b')
plt.xlabel(r'sSFR$')
plt.ylabel(r'$M/M_{vir}/f_{bary}$')
plt.legend([hot_plt,warm_plt,cool_plt,ISM_plt],['$10^6 < $T','$10^5 <$ T $< 10^6$','$10^4 <$ T $< 10^5$','T $< 10^4$'],loc = 2)
plt.show()
plt.savefig(outbase + '.gascensus_vSSFR.png')
plt.close()
