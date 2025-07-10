# Charlotte Christensen
# 11/11/19
# Examine the star formation history and halo assembly history for a particular halo

import matplotlib as mpl
import matplotlib.pylab as plt
import matplotlib.mlab as mlab
import matplotlib.colors as colors
import matplotlib.cm as cm
#mpl.use('Qt5Agg')
import pynbody
import math
import numpy as np
import socket
import matplotlib.gridspec as gridspec
import pandas as pd
import sys, os, glob, pickle
sys.path.append(os.path.abspath("/home/christenc/Code/python/python_analysis/"))
from modules.user_tools import task

presentation = False
if presentation:
    outfile_base = outfile_base + '_pres'
    plt.style.use(['default','/home/christenc/.config/matplotlib/presentation.mplstyle'])
    plt_width = 16 #8 #inches
    aspect_ratio = 3.0/4.0
    legendsize = 16
    dpi = 100
    lw = mpl.rcParams['lines.linewidth'] - 1        
else:
    plt.style.use(['default','/home/christenc/.config/matplotlib/article.mplstyle'])
    plt_width = 8 #inches
    aspect_ratio = 3.0/4.0
    legendsize = 5
    dpi = 300
    lw = mpl.rcParams['lines.linewidth']

if (socket.gethostname() == "ozma.grinnell.edu"):
    dataprefix = '/home/christensen/Code/Datafiles/'
    prefix = '/home/christensen/Storage2/UW/MolecH/Cosmo/'
    prefix_outfile = '/home/christensen/Plots/marvel/'    
else:
    dataprefix = '/home/christenc/Code/Datafiles/'
    prefix = '/home/christenc/Data/Sims/'
    prefix_outfile = '/home/christenc/Figures/marvel/'    

use_m200 = 1
if use_m200:
    ext = '.m200.'
else:
    ext = '.MAP.'
    
assem_history = task(
    np.load,
    start_text="Loading assembly history data",
    end_text="Loaded assembly history data",
    fail_text="Failed to load assembly history data",
    exit_on_fail=True
)(dataprefix + "/assembly_histories.npy", encoding="latin1", allow_pickle=True).item()
assem_history = pd.DataFrame.from_dict(assem_history)

tfile_cm = prefix + 'cptmarvel.cosmo25cmb/cptmarvel.cosmo25cmb.4096g5HbwK1BH/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096/ahf_200/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096'
tfile_r = prefix + 'rogue.cosmo25cmb/rogue.cosmo25cmb.4096g5HbwK1BH/rogue.cosmo25cmb.4096g5HbwK1BH.004096/ahf_200/rogue.cosmo25cmb.4096g5HbwK1BH.004096'
tfile_e = prefix + 'elektra.cosmo25cmb/elektra.cosmo25cmb.4096g5HbwK1BH/elektra.cosmo25cmb.4096g5HbwK1BH.004096/ahf_200/elektra.cosmo25cmb.4096g5HbwK1BH.004096'
tfile_s = prefix + 'storm.cosmo25cmb/storm.cosmo25cmb.4096g5HbwK1BH/storm.cosmo25cmb.4096g5HbwK1BH.004096/ahf_200/storm.cosmo25cmb.4096g5HbwK1BH.004096'
tfile_1 = prefix + 'h148.cosmo50PLK.3072g/h148.cosmo50PLK.3072g3HbwK1BH/h148.cosmo50PLK.3072g3HbwK1BH.004096/ahf_200/h148.cosmo50PLK.3072g3HbwK1BH.004096'
tfile_2 = prefix + 'h229.cosmo50PLK.3072g/h229.cosmo50PLK.3072gst5HbwK1BH/h229.cosmo50PLK.3072gst5HbwK1BH.004096/ahf_200/h229.cosmo50PLK.3072gst5HbwK1BH.004096' 
tfile_3 = prefix + 'h242.cosmo50PLK.3072g/h242.cosmo50PLK.3072gst5HbwK1BH/h242.cosmo50PLK.3072gst5HbwK1BH.004096/ahf_200/h242.cosmo50PLK.3072gst5HbwK1BH.004096'
tfile_4 = prefix + 'h329.cosmo50PLK.3072g/h329.cosmo50PLK.3072gst5HbwK1BH/h329.cosmo50PLK.3072gst5HbwK1BH.004096/ahf_200/h329.cosmo50PLK.3072gst5HbwK1BH.004096'
tfiles = [tfile_cm, tfile_e, tfile_r, tfile_s, tfile_1, tfile_2, tfile_3, tfile_4]

objs_pd = None 
for tfile in tfiles:
    objs_dat = []
    print(tfile)
    f=open(tfile + ext + 'data', 'rb')
    while 1:
        try:
            objs_dat.append(pickle.load(f))
        except EOFError:
            break        
    f.close()
        
    if len(objs_dat) == 1:
        temp = pd.DataFrame(objs_dat[0])
    else:
        temp = pd.DataFrame(objs_dat)
    base = tfile.split('/')[-1]
    simname = base.split('.')[0]
    if (base.split('.')[2])[0] == '6':
        simname = simname+'_6144'
    if not ('M_star' in temp.keys()):
        temp['M_star'] = temp['mstar']
    if not ('mass' in temp.keys()):            
        temp['mass'] = temp['mvir']            
    temp['sim'] = [simname]*len(temp)

    if objs_pd is None: 
        objs_pd = temp
    else:
        objs_pd = objs_pd.append(temp, ignore_index = True)        

# -------------------------------------------------
    
simname = 'storm'
haloid = 11
tfile = prefix + 'storm.cosmo25cmb/storm.cosmo25cmb.4096g5HbwK1BH/storm.cosmo25cmb.4096g5HbwK1BH.004096/ahf_200/storm.cosmo25cmb.4096g5HbwK1BH.004096'

simname = 'h148'
haloid = 55
haloid = 41
haloid = 15
tfile = prefix + 'h148.cosmo50PLK.3072g/h148.cosmo50PLK.3072g3HbwK1BH/h148.cosmo50PLK.3072g3HbwK1BH.004096/ahf_200/h148.cosmo50PLK.3072g3HbwK1BH.004096'

#simname = 'h229'
#tfile = prefix + 'h229.cosmo50PLK.3072g/h229.cosmo50PLK.3072gst5HbwK1BH/h229.cosmo50PLK.3072gst5HbwK1BH.004096/ahf_200/h229.cosmo50PLK.3072gst5HbwK1BH.004096'    
#haloid = 528

sim = pynbody.load(tfile)
h  = sim.halos()
simhalo = h[haloid]
cond = ((assem_history['Volume'] == simname) & (assem_history['halo grp @ z=0'] == haloid))
cond2 = ((objs_pd['haloid'] == str(haloid)) & (objs_pd['sim'] == simname))

# Plot the assembly history
halo = objs_pd[cond2]
time = np.array((halo['sfhbins'].tolist())[0])
dtime = time[1] - time[0]
sfhcum =  np.cumsum(np.array((halo['sfh']).tolist()[0])*dtime)
if len(time) != len(sfhcum):
    time = time[0:-1] + dtime/2

sfh2, time2 = np.histogram(simhalo.star['tform'].in_units('Gyr'), bins = len(sfhcum), weights = simhalo.star['massform'].in_units('Msol'))
sfhcum2 = np.cumsum(sfh2)
time2 = time2[0:-1] - (time2[1] - time2[0])/2

barymass = np.array((assem_history['Gas Mass History'][cond].tolist())[0]) + np.array((assem_history['Stellar Mass History'][cond].tolist())[0])
plt.clf()
fig1 = plt.figure(1,figsize=(plt_width,plt_width*aspect_ratio))
fig1.set_size_inches(plt_width,plt_width*aspect_ratio)
fig1.clear()
ax1 = fig1.add_subplot()
ax1.plot((assem_history['Time'][cond].tolist())[0],(assem_history['Assembly History'][cond].tolist())[0]/max((assem_history['Assembly History'][cond].tolist())[0]),color = 'b')
ax1.plot(time,sfhcum/max(sfhcum),color = 'g', linestyle = '-')
ax1.plot(time2,sfhcum2/max(sfhcum2),color = 'g', linestyle = '--')
ax1.plot((assem_history['Time'][cond].tolist())[0],(assem_history['Stellar Mass History'][cond].tolist())[0]/max(sfhcum2),color = 'g', linestyle = '-.')
ax1.plot((assem_history['Time'][cond].tolist())[0],(assem_history['Gas Mass History'][cond].tolist())[0]/max((assem_history['Gas Mass History'][cond].tolist())[0]),color = 'k')
ax1.plot(np.array((assem_history['Time'][cond].tolist())[0]), barymass/max(barymass), color = 'r')
