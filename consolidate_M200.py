# %% Initialization Cell

import time

print("Loading libraries...") #end and flush removed
t1 = time.time()

import sys, shutil, os, gc, copy, re
import numpy as np
import pandas as pd

# Now import custom packages
sys.path.append(os.path.abspath("."))
#import modules.cubehelix as cubehelix
from modules.user_tools import task, print_dict
#import modules.halo_trace as ht

t2 = time.time()
print("Libraries loaded in "+str( round(t2-t1, 3) )+"s.")

#clear any unused variables in RAM
gc.collect()

print("Working directory:", os.getcwd())

# %%

def time_to_snap_code(time):
    snap_code = str(np.int( np.round( time * 4096/13.8 ) ))
    return (4-len(snap_code))*"0" + snap_code

def snap_code_to_time(snap_code):
    return np.int(snap_code)*(13.8/4096.)

def load_txt(fname):
    ''' Load the file using std open'''
    f = open(fname,'r')

    data = []
    for line in f.readlines():
        data.append(line.replace('\n','').split(' '))

    f.close()

    return data

# %%

DCJL_sims = ["h148_6144","h329_6144","h148","h229","h242","h329"]
Marvel_sims = ["cptmarvel","rogue","elektra","storm"]
sims_to_include = DCJL_sims + Marvel_sims #["h329_6144"]

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

M200_dir = "/home/christen/myhome1/M200_amiga_stat"
M200_dir = "/home/christenc/Code/Datafiles/M200_amiga_stat"
folders = os.listdir( M200_dir )
folders = np.sort(folders)[::-1]
stat_data = {}

M200_dir_sf = "/home/christen/myhome1/M200_amiga_stat/004096/"
M200_dir_sf = "/home/christenc/Code/Datafiles/M200_amiga_stat/004096/"
M200_stat_files = os.listdir( M200_dir_sf )
stat_filepaths = {}
for file in M200_stat_files:
    if "_old" in file: continue
    for sim in sims_to_include:
        if sim.split("_")[0] in file:
            if len(sim.split("_")) > 1:
                if sim.split("_")[1] in file:
                    stat_filepaths[sim] = M200_dir_sf+file
            else:
                stat_filepaths[sim] = M200_dir_sf+file

# print_dict(stat_filepaths)
# sys.exit()

h = 0.73
for snapshot in folders:
    files = os.listdir( M200_dir + "/" + snapshot )
    for file in files:
        if "_old" in file: continue
        sim = file.split(".")[0]
        if "6144" in file:
            sim = sim+"_6144"
        if sim not in sims_to_include: continue
        if sim not in stat_data.keys():
            stat_data[sim] = {}
        if snapshot not in stat_data[sim].keys():
            stat_data[sim][snapshot] = {}

        stat_file_path = M200_dir + "/" + snapshot + "/" + file
        data = task(
            np.genfromtxt,
            start_text=f"Loading stat file for {sim} {snapshot} from {stat_file_path}",
            end_text=f"Loaded stat file for {sim} {snapshot} from {stat_file_path}",
            fail_text=f"Failed to load stat file for {sim} {snapshot} from {stat_file_path}",
            exit_on_fail=True
        )(stat_file_path, dtype="float", usecols=(0,1,3,5,6,7,8,13,14,15,20), skip_header=0)
        grps, Ntots, Nstars, Mvirs, Rvirs, Mgass, Mstars, xcs, ycs, zcs, sats = np.transpose(data)

        for i in range(len(grps)):
            grp = grps[i]
            Ntot = Ntots[i]
            Nstar = Nstars[i]
            Mvir = Mvirs[i]
            Mgas = Mgass[i]
            Mstar = Mstars[i]
            sat = sats[i]
            xc = xcs[i]/h
            yc = ycs[i]/h
            zc = zcs[i]/h
            Rvir = Rvirs[i]

            # Get Nearest(Rvir)
            # Identify either the MW or the largest halo within 5 Mpc
            if sim in DCJL_sims:
                index_of_MW = np.argwhere( grps == 1 )
                xc_mw = xcs[index_of_MW]/h
                yc_mw = ycs[index_of_MW]/h
                zc_mw = zcs[index_of_MW]/h
                Rvir_of_MW = Rvirs[index_of_MW]

                d_mw = np.sqrt( (xc-xc_mw)**2. + (yc-yc_mw)**2. + (zc-zc_mw)**2. ) / Rvir_of_MW

            elif sim in Marvel_sims:
                d_mw = 0


            stat_data[sim][snapshot][grp] = [ Mvir, Mstar, sat, d_mw, Rvir, Ntot, Nstar, Mgas ]

for sim in stat_data.keys():
    task(
        np.save,
        start_text=f"Saving {sim} data",
        end_text=f"Saved {sim} data",
        fail_text=f"Failed to save {sim} data",
        exit_on_fail=True,
        no_spinner=False
    )(f"consolidated_M200_{sim}.npy", stat_data[sim], allow_pickle=True)

task(
    np.save,
    start_text=f"Saving data",
    end_text=f"Saved data",
    fail_text=f"Failed to save data",
    exit_on_fail=True,
    no_spinner=False
)(f"consolidated_M200.npy", stat_data, allow_pickle=True)
