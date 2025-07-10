
import os
import pandas as pd
import pynbody
import tangos
import matplotlib.colors as colors
from matplotlib import path
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap
import sys, os, glob, pickle, socket
from scipy.interpolate import interp1d

def snap_to_time_code(snap):
    return 13.8*np.int(snap)/4096

sys.path.append(os.path.abspath("/home/christenc/Code/python/python_analysis/"))
from modules.user_tools import task, save_table_to_txt

if (socket.gethostname() == "ozma.grinnell.edu"):
    dataprefix = '/home/christensen/Code/Datafiles/' 
else:
    dataprefix = '/home/christenc/Code/Datafiles/'

if (socket.gethostname() == "quirm.math.grinnell.edu"):
    prefix = '/home/christenc/Data/Sims/'
    prefix_outfile = '/home/christenc/Figures/marvel/'
else:
    prefix = '/home/christensen/Storage2/UW/MolecH/Cosmo/'
    prefix_outfile = '/home/christensen/Plots/marvel/'

hubble = 0.6776942783267969

presentation = False
if presentation:
    outfile_base = outfile_base + '_pres'
    plt.style.use(['default','/home/christenc/.config/matplotlib/presentation.mplstyle'])
    plt_width = 8 #inches
    aspect_ratio = 3.0/4.0
    legendsize = 16
    dpi = 200
    markersize = 100
    ms_scale = 1
    lw = mpl.rcParams['lines.linewidth'] - 1
    edgewidth = 1
else:
    plt.style.use(['default','/home/christenc/.config/matplotlib/article.mplstyle'])
    plt_width = 3.5 #inches
    aspect_ratio = 3.0/4.0
    legendsize = 5
    dpi = 300
    markersize = 25
    ms_scale = 0.25
    lw = mpl.rcParams['lines.linewidth']
    edgewidth = 0.5

tangos.core.init_db('/home/christenc/Storage/tangos_db/JL_r200.db')

filepath = '/home/christenc/Data/Sims/h148.cosmo50PLK.3072g/h148.cosmo50PLK.3072g3HbwK1BH/'
filename = filepath.split('/')[-2]
simname = 'Sandra'
simshort = 'h148'

"""
filepath = '/home/christenc/Data/Sims/h229.cosmo50PLK.3072g/h229.cosmo50PLK.3072gst5HbwK1BH/'
filename = filepath.split('/')[-2]
simname = 'Ruth'
simshort = 'h229'
"""
"""
filepath = '/home/christenc/Data/Sims/h242.cosmo50PLK.3072g/h242.cosmo50PLK.3072gst5HbwK1BH/'
filename = filepath.split('/')[-2]
simname = 'Sonia'
simshort = 'h242'
"""
"""
filepath = '/home/christenc/Data/Sims/h329.cosmo50PLK.3072g/h329.cosmo50PLK.3072gst5HbwK1BH/'
filename = filepath.split('/')[-2]
simname = 'Elena'
simshort = 'h329'
"""

#tangos.core.init_db('/home/christenc/Storage/tangos_db/Marvel_r200.db')
"""
filepath = '/home/christenc/Data/Sims/storm.cosmo25cmb/storm.cosmo25cmb.4096g5HbwK1BH/'
filename = filepath.split('/')[-2]
simname = 'Storm'
simshort = 'storm'
"""
"""
filepath = '/home/christenc/Data/Sims/elektra.cosmo25cmb/elektra.cosmo25cmb.4096g5HbwK1BH/'
filename = filepath.split('/')[-2]
simname = 'Elektra'
simshort = 'elektra'

filepath = '/home/christenc/Data/Sims/rogue.cosmo25cmb/rogue.cosmo25cmb.4096g5HbwK1BH/'
filename = filepath.split('/')[-2]
simname = 'Rogue'
simshort = 'rogue'
"""
"""
filepath = '/home/christenc/Data/Sims/cptmarvel.cosmo25cmb/cptmarvel.cosmo25cmb.4096g5HbwK1BH/'
filename = filepath.split('/')[-2]
simname = 'Cptmarvel'
simshort = 'cptmarvel'
"""

stat_data = np.load("/home/christenc/Code/Datafiles/consolidated_M200_" + simshort + ".npy", encoding="latin1", allow_pickle=True)

stat_data = task(
    np.load,
    start_text="Loading M200 stat data",
    end_text="Loaded M200 stat data",
    fail_text="Failed to load M200 stat data",
    exit_on_fail=True
)("/home/christenc/Code/Datafiles/consolidated_M200_" + simshort + ".npy", encoding="latin1", allow_pickle=True).item()

timesteps = tangos.get_simulation("snapshots_200crit_" + simshort).timesteps
snap_nums = [re.findall(r'.00+[\d]+', str(snap))[0][3:] for snap in timesteps]

path = '/home/christenc/Data/Sims/'
trace_file_paths = {
        "cptmarvel": path + 'cptmarvel.cosmo25cmb/cptmarvel.cosmo25cmb.4096g5HbwK1BH/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.M200.trace_back.hdf5'
        ,"rogue": path + 'rogue.cosmo25cmb/rogue.cosmo25cmb.4096g5HbwK1BH/rogue.cosmo25cmb.4096g5HbwK1BH.004096.M200.trace_back.hdf5'
        ,"elektra": path + 'elektra.cosmo25cmb/elektra.cosmo25cmb.4096g5HbwK1BH/elektra.cosmo25cmb.4096g5HbwK1BH.004096.M200.trace_back.hdf5'
        ,"storm": path + 'storm.cosmo25cmb/storm.cosmo25cmb.4096g5HbwK1BH/storm.cosmo25cmb.4096g5HbwK1BH.004096.M200.trace_back.hdf5'
        ,"h148_6144": path + 'h148.cosmo50PLK.6144g/h148.cosmo50PLK.6144g3HbwK1BH/h148.cosmo50PLK.6144g3HbwK1BH.004096.M200.trace_back.hdf5'
        ,"h329_6144": path + 'h329.cosmo50PLK.6144g/h329.cosmo50PLK.6144g5HbwK1BH/h329.cosmo50PLK.6144g5HbwK1BH.004096.M200.trace_back.hdf5'
        ,"h148": path + 'h148.cosmo50PLK.3072g/h148.cosmo50PLK.3072g3HbwK1BH/h148.cosmo50PLK.3072g3HbwK1BH.004096.M200.trace_back.hdf5'        
        ,"h229": path + 'h229.cosmo50PLK.3072g/h229.cosmo50PLK.3072gst5HbwK1BH/h229.cosmo50PLK.3072gst5HbwK1BH.004096.M200.trace_back.hdf5'
        ,"h242": path + 'h242.cosmo50PLK.3072g/h242.cosmo50PLK.3072gst5HbwK1BH/h242.cosmo50PLK.3072gst5HbwK1BH.004096.M200.trace_back.hdf5'
        ,"h329": path + 'h329.cosmo50PLK.3072g/h329.cosmo50PLK.3072gst5HbwK1BH/h329.cosmo50PLK.3072gst5HbwK1BH.004096.M200.trace_back.hdf5'
}

sim_trace = pd.read_hdf( trace_file_paths[simshort] )

objs_dat = []
f=open(filepath + filename + '.00' + snap_nums[-1] + '/ahf_200/' + filename + '.00' + snap_nums[-1] + '.m200.dist.data', 'rb')
#f=open(filepath + filename + '.' + snap_nums[-1] + '/' + filename + '.00' + snap_nums[-1] + '.m200.dist.data', 'rb')
while 1:
    try:
        objs_dat.append(pickle.load(f))
    except EOFError:
        break        
f.close()
if len(objs_dat) == 1:
    objs_pd_200 = pd.DataFrame(objs_dat[0])
else:
    objs_pd_200 = pd.DataFrame(objs_dat)
    
tangos_halo = tangos.get_halo("snapshots_200crit_" + simshort + '/' + filename + '.004096/halo_1')
#time, scalefactor, id_main, mvir_main, rvir_main = tangos_halo.calculate_for_progenitors("t()", "a()", "halo_number()", "Mvir", "Rvir")
#timestep = tangos.get_simulation("snapshots_200crit_h329").timesteps[-1]
#grpz0 = timestep.calculate_all("Grp") #Redshift zero data for all galaxies in the tangos database

# Create a function to go from time to scalefactor
#t_to_a = interp1d(np.flip(np.append(time,0)), np.flip(np.append(scalefactor,0)), kind='cubic')

alltimes = np.array([snap_to_time_code(snap) for snap in list(objs_pd_200.loc[1]['mpb_grp'].keys())])
"""
alltimes = []
for timestep in timesteps:
    ts = np.array(timestep.calculate_all("t()")[0])
    if len(ts) !=0:
        alltimes.append(ts[0])
    else: alltimes.append(0)
alltimes = np.array(alltimes)
"""

# If all the data is consistent then the following should be equal:
# np.array(snap_nums)  -- minus the two starting zeros
# stat_data.keys()
# sim_trace.keys() -- minus the '004096'

# Calculate the distance when the gravitational force was at the peak and the time that it happened
"""
for i in arange(0,len(objs_pd_200)):  
    if np.max(np.array(list(objs_pd_200.iloc[i]['mpb_gforce'].values()))) != 0:
        #print(np.max(np.array(list(objs_pd_200.iloc[i]['mpb_gforce'].values()))))
        ind = np.argmax(np.array(list(objs_pd_200.iloc[i]['mpb_gforce'].values())))
        objs_pd_200.iloc[i,objs_pd_200.columns.get_loc('min_dist')] = np.sqrt(np.array(list(objs_pd_200.iloc[i]['mpb_gforce_mass'].values()))/np.array(list(objs_pd_200.iloc[i]['mpb_gforce'].values())))[ind]
        objs_pd_200.iloc[i,objs_pd_200.columns.get_loc('min_dist_t')] = alltimes[ind]
        print(objs_pd_200.iloc[i]['min_dist_t'],objs_pd_200.iloc[i]['min_dist'])
    else:
        objs_pd_200.iloc[i,objs_pd_200.columns.get_loc('min_dist_t')] = objs_pd_200.iloc[i]['time']
        objs_pd_200.iloc[i,objs_pd_200.columns.get_loc('min_dist')] = objs_pd_200.iloc[i]['massiveDist']
"""
#objs_pd_200.to_pickle(filepath + filename + '.00' + snap_nums[-1] + '/ahf_200/' + filename + '.00' + snap_nums[-1] + '.m200.dist.data')

plt.clf()
fig1 = plt.figure(1)
ax1 = fig1.add_subplot(1,1,1)
for index, row in objs_pd_200.iloc[1:].iterrows():
    ts_link = np.array(list(row['mpb_tindex_mass'].values())) != 0
    ax1.plot(alltimes[ts_link],((np.array(list(row['mpb_tindex_mass'].values()))/np.array(list(row['mpb_tindex'].values())))**(1/3.0))[ts_link],alpha = 0.4)
row = objs_pd_200.iloc[0]
ax1.plot(alltimes[ts_link],((np.array(list(row['mpb_tindex_mass'].values()))/np.array(list(row['mpb_tindex'].values())))**(1/3.0))[ts_link],alpha = 0.4)
#ax1.plot(time,rvir_main,color = 'k')    
#ax.axis([0, 14, 2e5, 3e9])
ax1.set_yscale('log')
ax1.set_xlabel('Time [Gyr]')
ax1.set_ylabel('Physical Distance to Most Impactful Companion')

# Distance to main halo
fig2 = plt.figure(2)
ax2 = fig2.add_subplot(1,1,1)
for index, row in objs_pd_200.iloc[1:].iterrows():
    ts_link = np.array(list(row['mpb_massive_dist'].values())) != 0
    ax2.plot(alltimes[ts_link],(np.array(list(row['mpb_massive_dist'].values())))[ts_link],alpha = 0.4)
row = objs_pd_200.iloc[0]
ax2.plot(alltimes[ts_link],(np.array(list(row['mpb_massive_dist'].values())))[ts_link],alpha = 0.4)
ax2.plot(time,rvir_main,color = 'k')    
#ax.axis([0, 14, 2e5, 3e9])
ax2.set_yscale('log')
ax2.set_xlabel('Time [Gyr]')
ax2.set_ylabel('Physical Distance to Closest Massive (> 1e11 Msun) Companion')

# Distance to main halo
fig3 = plt.figure(3)
ax3 = fig3.add_subplot(1,1,1)
for index, row in objs_pd_200.iloc[1:].iterrows():
    ts_link = np.array(list(row['mpb_massive_dist'].values())) != 0
    ax3.plot(alltimes[ts_link],(np.array(list(row['mpb_massive_dist'].values())) -np.sqrt(np.array(list(row['mpb_tindex_mass'].values()))/np.array(list(row['mpb_tindex'].values()))))[ts_link] ,alpha = 0.4)
#ax.axis([0, 14, 2e5, 3e9])
ax3.plot(time,rvir_main,color = 'k')   
ax3.set_xlabel('Time [Gyr]')
ax3.set_ylabel('Distance between massive halo and most impactful halo')

# Tidal Index
fig4 = plt.figure(4)
ax4 = fig4.add_subplot(1,1,1)
for index, row in objs_pd_200.iloc[1:].iterrows():
    ts_link = np.array(list(row['mpb_tindex_mass'].values())) != 0
    ax4.plot(alltimes[ts_link],(np.array(list(row['mpb_tindex'].values())))[ts_link] ,alpha = 0.4)
row = objs_pd_200.iloc[0]
ax4.plot(alltimes[ts_link],(np.array(list(row['mpb_tindex'].values())))[ts_link] ,alpha = 0.4)    
plt.gca().set_prop_cycle(None)
for index, row in objs_pd_200.iloc[1:].iterrows():
    print(row['min_dist_t'],row['min_dist'])
    ax4.plot([row['max_tindex_t'],row['max_tindex_t']],[row['max_tindex'],row['max_tindex']],marker = 'o')  
row = objs_pd_200.iloc[0]
ax4.plot([row['min_dist_t'],row['min_dist_t']],[row['max_tindex'],row['max_tindex']],marker = 'o') 
ax4.set_xlabel('Time [Gyr]')
ax4.set_ylabel('Tidal Index')
ax4.set_yscale('log')

fig5 = plt.figure(5)
ax5 = fig5.add_subplot(1,1,1)
for index, row in objs_pd_200.iloc[1:].iterrows():
    ts_link = np.array(list(row['mpb_tindex_mass'].values())) != 0
    ax5.plot(alltimes[ts_link],((np.array(list(row['mpb_tindex_mass'].values()))/np.array(list(row['mpb_tindex'].values())))**(1/3.0))[ts_link],alpha = 0.4)
plt.gca().set_prop_cycle(None)
for index, row in objs_pd_200.iloc[1:].iterrows():
    ts_link = np.array(list(row['mpb_tindex_mass'].values())) != 0
    ax5.plot(alltimes[ts_link],(np.array(list(row['mpb_massive_dist'].values())))[ts_link],alpha = 0.4,linestyle = '--')
plt.gca().set_prop_cycle(None)
for index, row in objs_pd_200.iloc[1:].iterrows():
#    print(row['min_dist_t'],row['min_dist'])
#    ax5.plot([row['min_dist_t'],row['min_dist_t']],[row['min_dist'],row['min_dist']],marker = 'o')
    print(row['max_tindex_t'],row['max_tindex'])
    ax5.plot([row['max_tindex_t'],row['max_tindex_t']],[row['max_tindex_dist'],row['max_tindex_dist']],marker = 'o') 
    
#row = objs_pd_200.iloc[0]
#ax5.plot(alltimes[ts_link],np.sqrt(np.array(list(row['mpb_tindex_mass'].values()))/np.array(list(row['mpb_tindex'].values())))[ts_link],alpha = 0.4)
#ax5.plot(alltimes[ts_link],(np.array(list(row['mpb_massive_dist'].values())))[ts_link],alpha = 0.4,linestyle = '--')
ax5.plot(time,rvir_main,color = 'k')    
#ax.axis([0, 14, 2e5, 3e9])
ax5.set_yscale('log')
ax5.set_xlabel('Time [Gyr]')
ax5.set_ylabel('Physical Distance to Main Halo Progenitor/Most Impactful Halo')


for index, row in objs_pd_200.iloc[1:].iterrows():
    plt.clf()
    fig1 = plt.figure(1)
    ax1a = fig1.add_subplot(1,1,1)
    ax1a.set_xlabel('Time [Gyr]')
    ax1a.set_ylabel('Physical Distance to Most Impactful Companion')
    ts_link = np.array(list(row['mpb_tindex_mass'].values())) != 0
    ax1a.plot(alltimes[ts_link],((np.array(list(row['mpb_tindex_mass'].values()))/np.array(list(row['mpb_tindex'].values())))**(1/3.0))[ts_link],alpha = 0.4)
    ax1b = ax1a.twinx()
    ax1b.plot(alltimes[ts_link],np.array(list(row['mpb_mvir'].values()))[ts_link]/max(np.array(list(row['mpb_mvir'].values()))[ts_link]),linestyle = '--')
    ax1b.axis([0, 14, 0, 1])
    wait = input('continue')
    ax1a.clear()
    ax1b.clear()
    fig1.clear()
    close(fig1)
row = objs_pd_200.iloc[0]
ax1.plot(alltimes[ts_link],((np.array(list(row['mpb_tindex_mass'].values()))/np.array(list(row['mpb_tindex'].values())))**(1/3.0))[ts_link],alpha = 0.4)
#ax1.plot(time,rvir_main,color = 'k')    
#ax.axis([0, 14, 2e5, 3e9])

for index, row in objs_pd_200.iloc[1:].iterrows():
    plt.clf()
    fig1 = plt.figure(1)
    ax1a = fig1.add_subplot(1,1,1)
    ax1a.set_xlabel('Time [Gyr]')
    ax1a.set_ylabel('Physical Distance to Massive Companion')
    ax1a.set_title(row['sim']+' '+row['haloid'])
    ts_link = np.array(list(row['mpb_massive_dist'].values())) != 0
    ax1a.plot(alltimes[ts_link],np.array(list(row['mpb_massive_dist'].values()))[ts_link],alpha = 0.4)
    ax1b = ax1a.twinx()
    ax1b.plot(alltimes[ts_link],np.array(list(row['mpb_mvir'].values()))[ts_link]/max(np.array(list(row['mpb_mvir'].values()))[ts_link]),linestyle = '--')
    ax1b.axis([0, 14, 0, 1])
    wait = input('continue')
    ax1a.clear()
    ax1b.clear()
    fig1.clear()
    close(fig1)

for index, row in objs_pd_200.iloc[1:].iterrows():
    plt.clf()
    fig1 = plt.figure(1)
    ax1a = fig1.add_subplot(1,1,1)
    ax1a.set_xlabel('Time [Gyr]')
    ax1a.set_ylabel('Tidal Index')
    ax1a.set_title(row['sim']+' '+row['haloid'])
    ts_link = np.array(list(row['mpb_massive_dist'].values())) != 0
    ax1a.plot(alltimes[ts_link],np.array(list(row['mpb_tindex'].values()))[ts_link],alpha = 0.4)
    ax1b = ax1a.twinx()
    ax1b.plot(alltimes[ts_link],np.array(list(row['mpb_mvir'].values()))[ts_link]/max(np.array(list(row['mpb_mvir'].values()))[ts_link]),linestyle = '--')
    ax1b.axis([0, 14, 0, 1])

    fig2 = plt.figure(2)
    ax2a = fig2.add_subplot(1,1,1)
    ax2a.set_xlabel('Time [Gyr]')
    ax2a.set_ylabel('Physical Distance to Massive Companion')
    ax2a.set_title(row['sim']+' '+row['haloid'])
    ts_link = np.array(list(row['mpb_massive_dist'].values())) != 0
    ax2a.plot(alltimes[ts_link],np.array(list(row['mpb_massive_dist'].values()))[ts_link],alpha = 0.4)
    ax2a.plot(alltimes[ts_link],((np.array(list(row['mpb_tindex_mass'].values()))/np.array(list(row['mpb_tindex'].values())))**(1/3.0))[ts_link],alpha = 0.4)
    ax2b = ax2a.twinx()
    ax2b.plot(alltimes[ts_link],np.array(list(row['mpb_mvir'].values()))[ts_link]/max(np.array(list(row['mpb_mvir'].values()))[ts_link]),linestyle = '--')
    ax2b.axis([0, 14, 0, 1])
    
    wait = input('continue')
    ax1a.clear()
    ax1b.clear()
    fig1.clear()
    close(fig1)
    ax2a.clear()
    ax2b.clear()
    fig2.clear()
    close(fig2)
