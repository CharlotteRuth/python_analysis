# Charlotte Christensen
# February 16, 2022
# Calculate the closest massive neighbor to a galaxy

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

def snap_to_time_code(snap):
    return 13.8*np.int(snap)/4096

sys.path.append(os.path.abspath("/home/christenc/Code/python/python_analysis/"))
from modules.user_tools import task, save_table_to_txt

M200_dir = "/home/christenc/Code/Datafiles/M200_amiga_stat/"

use_halo_trace = 1 # Use the data from Elaad's halo_trace.py, not tangos linking

read_stat_data = 0 # Read the halo information from the compilation of the amiga.stat files. Unfortunately, this doesn't work when looking for distances to halos outside the resolved region.

read_AHF = 0 # Read the properties from AHF_halos

read_tangosDB = 1 # Read the properties from the tangos DB

mass_threshold = 1e11 # The minimum mass in solar masses for what is considered a ``massive galaxy"

tangos.core.init_db('/home/christenc/Storage/tangos_db/JL_r200.db')
tangos.core.init_db('/home/christenc/Storage/tangos_db/JL_r200_N100min.db')

"""
filepath = '/home/christenc/Data/Sims/h148.cosmo50PLK.3072g/h148.cosmo50PLK.3072g3HbwK1BH/'
filename = filepath.split('/')[-2]
simname = 'Sandra'
simshort = 'h148'
"""
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

#tangos.core.init_db('/home/christenc/Storage/tangos_db/Marvel_r200.db') # Doesn't have information on unresolved halos
tangos.core.init_db('/home/christenc/Storage/tangos_db/Marvel_r200_N100min.db')

filepath = '/home/christenc/Data/Sims/storm.cosmo25cmb/storm.cosmo25cmb.4096g5HbwK1BH/'
filename = filepath.split('/')[-2]
simname = 'Storm'
simshort = 'storm'

"""
filepath = '/home/christenc/Data/Sims/elektra.cosmo25cmb/elektra.cosmo25cmb.4096g5HbwK1BH/'
filename = filepath.split('/')[-2]
simname = 'Elektra'
simshort = 'elektra'

filepath = '/home/christenc/Data/Sims/rogue.cosmo25cmb/rogue.cosmo25cmb.4096g5HbwK1BH/'
filename = filepath.split('/')[-2]
simname = 'Rogue'
simshort = 'rogue'


filepath = '/home/christenc/Data/Sims/cptmarvel.cosmo25cmb/cptmarvel.cosmo25cmb.4096g5HbwK1BH/'
filename = filepath.split('/')[-2]
simname = 'Cptmarvel'
simshort = 'cptmarvel'
"""


tangos.core.init_db('/home/christenc/Storage/tangos_db/JLmint_r200_N100.db')
"""
filepath = '/home/christenc/Data/Sims/h148.cosmo50PLK.6144g/h148.cosmo50PLK.6144g3HbwK1BH/'
filename = filepath.split('/')[-2]
simname = 'Sandra_mint'
simshort = 'h148mint'
"""
"""
filepath = '/home/christenc/Data/Sims/h329.cosmo50PLK.6144g/h329.cosmo50PLK.6144g5HbwK1BH/'
filename = filepath.split('/')[-2]
simname = 'Elena_mint'
simshort = 'h329mint'
"""

timesteps = tangos.get_simulation("snapshots_200crit_" + simshort).timesteps
snap_nums = [re.findall(r'.00+[\d]+', str(snap))[0][1:] for snap in timesteps]

stat_data = task(
    np.load,
    start_text="Loading M200 stat data",
    end_text="Loaded M200 stat data",
    fail_text="Failed to load M200 stat data",
    exit_on_fail=True
)("/home/christenc/Code/Datafiles/consolidated_M200_" + simshort + ".npy", encoding="latin1", allow_pickle=True).item()

path = '/home/christenc/Data/Sims/'
trace_file_paths = {
        "cptmarvel": path + 'cptmarvel.cosmo25cmb/cptmarvel.cosmo25cmb.4096g5HbwK1BH/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.M200.trace_back.hdf5'
        ,"rogue": path + 'rogue.cosmo25cmb/rogue.cosmo25cmb.4096g5HbwK1BH/rogue.cosmo25cmb.4096g5HbwK1BH.004096.M200.trace_back.hdf5'
        ,"elektra": path + 'elektra.cosmo25cmb/elektra.cosmo25cmb.4096g5HbwK1BH/elektra.cosmo25cmb.4096g5HbwK1BH.004096.M200.trace_back.hdf5'
        ,"storm": path + 'storm.cosmo25cmb/storm.cosmo25cmb.4096g5HbwK1BH/storm.cosmo25cmb.4096g5HbwK1BH.004096.M200.trace_back.hdf5'
        ,"h148mint": path + 'h148.cosmo50PLK.6144g/h148.cosmo50PLK.6144g3HbwK1BH/h148.cosmo50PLK.6144g3HbwK1BH.004096.M200.trace_back.hdf5'
        ,"h329mint": path + 'h329.cosmo50PLK.6144g/h329.cosmo50PLK.6144g5HbwK1BH/h329.cosmo50PLK.6144g5HbwK1BH.004096.M200.trace_back.hdf5'
        ,"h148": path + 'h148.cosmo50PLK.3072g/h148.cosmo50PLK.3072g3HbwK1BH/h148.cosmo50PLK.3072g3HbwK1BH.004096.M200.trace_back.hdf5'        
        ,"h229": path + 'h229.cosmo50PLK.3072g/h229.cosmo50PLK.3072gst5HbwK1BH/h229.cosmo50PLK.3072gst5HbwK1BH.004096.M200.trace_back.hdf5'
        ,"h242": path + 'h242.cosmo50PLK.3072g/h242.cosmo50PLK.3072gst5HbwK1BH/h242.cosmo50PLK.3072gst5HbwK1BH.004096.M200.trace_back.hdf5'
        ,"h329": path + 'h329.cosmo50PLK.3072g/h329.cosmo50PLK.3072gst5HbwK1BH/h329.cosmo50PLK.3072gst5HbwK1BH.004096.M200.trace_back.hdf5'
}

sim_trace = pd.read_hdf( trace_file_paths[simshort] )

objs_dat = []
f=open(filepath + filename + '.' + snap_nums[-1] + '/ahf_200/' + filename + '.' + snap_nums[-1] + '.m200.data', 'rb')
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
#objs_pd_200.set_index('haloid') #doesn't appear to be necessary but would allow access by grp

# If all the data is consistent then the following should be equal:
# np.array(snap_nums)  -- minus the two starting zeros
# stat_data.keys()
# sim_trace.keys() -- minus the '004096'

# Also, compare sim_trace.index to objs_pd_200['haloid'] to ensure that all of the latter are contained in the former

sim = pynbody.load(filepath + filename + '.' + snap_nums[-1] + '/ahf_200/' + filename + '.' + snap_nums[-1])
h = sim.properties['h']
tangos_halo = tangos.get_halo("snapshots_200crit_" + simshort + '/' + filename + '.004096/halo_1')
time, scalefactor, id_main, mvir_main, rvir_main = tangos_halo.calculate_for_progenitors("t()", "a()", "halo_number()", "Mvir", "Rvir")
timestep = tangos.get_simulation("snapshots_200crit_" + simshort).timesteps[-1]
grp_key = "halo_number()" #-- 6/7/22 had to change from "Grp" because now using AHF data, not amiga.stat
grpz0 = (timestep.calculate_all(grp_key))[0] #Redshift zero data for all galaxies in the tangos database 

# Check if the Mvir and Rvir data is in h^-1 units (as from AHF) or absolute (as from amiga) and adjust accordingly
mvir_all = (timestep.calculate_all("Mvir"))[0]
ind_match = int(np.argwhere(grpz0 == int(objs_pd_200['haloid'][0]))[0]) # index where the halo id matches
if (objs_pd_200['mvir'][0] < 1.01*mvir_all[ind_match]) and (objs_pd_200['mvir'][0] < 0.99*mvir_all[ind_match]):
    print("Using data from .amiga.stat")
    amiga_dat = True
else:
    print("Using data from AHF:", objs_pd_200['mvir'][0], mvir_all[ind_match]/h)
    amiga_dat = False
    mvir_main = mvir_main/h
    rvir_main = rvir_main/h
    

#alltimes = np.flip((tangos.get_halo("snapshots_200crit_" + simshort + '/' + '%.004096/halo_0')).calculate_for_progenitors("t()")[0])

alltimes = []
for timestep in timesteps:
    #print(np.array(timestep.calculate_all("t()")[0]))
    ts = np.array(timestep.calculate_all("t()")[0])
    if len(ts) !=0:
        alltimes.append(ts[0])
    else: alltimes.append(0)

#alltimes = []
#for timestep in timesteps: alltimes.append(timestep.calculate_all("t()")[0][0]) #h138, h329, h242
#for timestep in timesteps: alltimes.append(timestep.calculate_all("t()")[0])    
alltimes = np.array(alltimes)

if use_halo_trace:
    snap_nums = np.union1d(snap_nums, sim_trace.keys())
    snap_steps = np.append(np.flip(sim_trace.keys()),'004096')
    alltimes = np.array([snap_to_time_code(snap) for snap in snap_nums])
else:
    snap_steps = snap_nums

# For each halo in the list, calculate the main progenitor at all time steps
"""
mpb.clear()
mpb_mvir.clear()
mpb_dist.clear()
"""
#mpb = []
mpb_grp = []
mpb_mvir = []
mpb_main_dist = []
#mpb_close_id = []
mpb_close_grp = []
mpb_tindex = []
mpb_tindex_grp = []
mpb_tindex_mass = []
mpb_massive_dist = []
for index, row in objs_pd_200.iterrows():
    print(row['haloid'])

    #ids = np.empty(len(snap_nums))*0
    grp = np.empty(len(snap_nums))*0
    mvir = np.empty(len(snap_nums))*0
    if use_halo_trace:
        it = 0
        for snap_step in snap_nums[:-1]:
            try:
                grp[it] = sim_trace.loc[int(row['haloid'])][snap_step]
                if grp[it] == -1 or np.isnan(grp[it]):
                    grp[it] = 0
                else:
                # Structure of stat_data[snap][haloid] is [ Mvir, Mstar, sat, d_mw, Rvir, Ntot, Nstar, Mgas ]
                    mvir[it] = stat_data[snap_step][grp[it]][0]
            except:
                pass
            #print(snap_step,grp[it],mvir[it])  
            it = it + 1
            
        grp[-1] = int(row['haloid']) #last value is the z=0 value
        mvir[-1] = row['mvir']
    else:
        #haloind = (np.where(grpz0[0] == int(row['haloid'])))[0]

        haloind = np.where(timesteps[-1].calculate_all(grp_key)[0] == int(row['haloid']))
        if len(haloind[0]) == 1: # If the halo is in the tangos database
            tangos_halo = timesteps[-1][int(haloind[0][0]) + 1]
            """
            # otherwise, if doing haloind = (np.where(grpz0[0] == int(row['haloid'])))[0]
            haloind = haloind[0] + 1
            tangos_halo = tangos.get_halo("snapshots_200crit_" + simshort + '/' + '%.004096/halo_' + str(haloind))
            """
            
            if int(tangos_halo.calculate(grp_key)) != int(row['haloid']):
                print("Error for halo: ", tangos_halo.calculate(grp_key))
        
            if tangos_halo == None:
                #ids[-1] = int(haloind)
                grp[-1] = int(row['haloid'])        
                mvir[-1] = row['mvir']
            else:
                time_short, grp_short, mvir_short, z_short  = tangos_halo.calculate_for_progenitors("t()", grp_key, "Mvir", "z()")
                sort_times = numpy.searchsorted(alltimes,time_short)
                if not amiga_dat:
                    mvir_short = mvir_short/h
                #ids[sort_times] = ids_short
                grp[sort_times] = grp_short        
                mvir[sort_times] = mvir_short

    #mpb_temp = {snap_nums[i]: int(ids[i]) for i in range(len(snap_nums))}
    #mpb.append(mpb_temp.copy())
    mpb_grp_temp = {snap_nums[i]: int(grp[i]) for i in range(len(snap_nums))}
    mpb_grp.append(mpb_grp_temp.copy())    
    mpb_mvir_temp = {snap_nums[i]: mvir[i] for i in range(len(snap_nums))}
    mpb_mvir.append(mpb_mvir_temp.copy())
    mpb_main_dist_temp = {snap_nums[i]: 0 for i in range(len(snap_nums))}
    mpb_main_dist.append(mpb_main_dist_temp.copy())
    #mpb_close_id_temp = {snap_nums[i]: 0 for i in range(len(snap_nums))}
    #mpb_close_id.append(mpb_close_id_temp.copy())
    mpb_close_grp_temp = {snap_nums[i]: 0 for i in range(len(snap_nums))}
    mpb_close_grp.append(mpb_close_grp_temp.copy())    
    mpb_tindex_temp = {snap_nums[i]: 0 for i in range(len(snap_nums))}
    mpb_tindex.append(mpb_tindex_temp.copy())
    mpb_tindex_grp_temp = {snap_nums[i]: 0 for i in range(len(snap_nums))}
    mpb_tindex_grp.append(mpb_tindex_grp_temp.copy())
    mpb_tindex_mass_temp = {snap_nums[i]: 0 for i in range(len(snap_nums))}
    mpb_tindex_mass.append(mpb_tindex_mass_temp.copy())
    mpb_massive_dist_temp = {snap_nums[i]: 0 for i in range(len(snap_nums))}
    mpb_massive_dist.append(mpb_massive_dist_temp.copy())    
    #print(mpb)



mindist = {'min_dist': ""}
mindistt = {'min_dist_t': ""}
maxtindex = {'max_tindex': ""}
maxtindext = {'max_tindex_t': ""}
maxtindextd = {'max_tindex_dist': ""}
objs_pd_200 = objs_pd_200.join(pd.DataFrame(columns=mindist))
objs_pd_200 = objs_pd_200.join(pd.DataFrame(columns=mindistt))
objs_pd_200 = objs_pd_200.join(pd.DataFrame(columns=maxtindex))
objs_pd_200 = objs_pd_200.join(pd.DataFrame(columns=maxtindext))
objs_pd_200 = objs_pd_200.join(pd.DataFrame(columns=maxtindextd))
#objs_pd_200['mpb_id'] = mpb
objs_pd_200['mpb_grp'] = mpb_grp
objs_pd_200['mpb_mvir'] = mpb_mvir
objs_pd_200['mpb_main_dist'] = mpb_main_dist
#mpb_close_id = mpb_main_dist.copy()
#mpb_tindex = mpb_main_dist.copy()
#mpb_tindex_id = mpb_main_dist.copy()
objs_pd_200['mpb_close_grp'] = mpb_close_grp
objs_pd_200['mpb_tindex'] = mpb_tindex
objs_pd_200['mpb_tindex_grp'] = mpb_tindex_grp
objs_pd_200['mpb_tindex_mass'] = mpb_tindex_mass
objs_pd_200['mpb_massive_dist'] = mpb_massive_dist

# If the tangos halo tracing and Elaads halo_trace are consistent than the following should be the same when use_halo_trace = 0
i = 3
haloid = objs_pd_200['haloid'][i]
print(objs_pd_200[objs_pd_200['haloid']==haloid]['mpb_grp'][i].values()) 
print(sim_trace.loc[int(haloid)])

i = 0
haloid = objs_pd_200['haloid'][i]
ind = np.array(list(objs_pd_200[objs_pd_200['haloid']==haloid]['mpb_mvir'][i].values())) > mass_threshold
if np.sum(ind) > 0 :
    print("Time when main progenitor greater than mass threshold: ",
          snap_to_time_code((np.array(list(objs_pd_200[objs_pd_200['haloid']==haloid]['mpb_mvir'][i].keys())))[ind][0]))
else:
    print("Main progenitor never more massive than threshold")


# Loop through the time steps, calculating the distance to the nearest halo
# timestep = snap_nums[0]

#sim = pynbody.load(filepath + 'snapshots_200crit_' + simshort + '/' + filename + '.' + snap_steps[-1])
#h_dummy = sim.halos(dummy = True)
#h = h_dummy[1].properties['h']

a_arr = []
for snap_step in snap_steps:
    print(filepath + 'snapshots_200crit_' + simshort + '/' + filename + '.' + snap_step)
    halomasses = np.array([row['mpb_mvir'][snap_step] for index, row in objs_pd_200.iterrows()])
    #if max(halomasses) == 0: continue
    minmass = 1e4 #np.min(halomasses[np.nonzero(halomasses)])

    if read_tangosDB: # Read data from the tangos DB
        try:
            timestep = tangos.get_timestep("snapshots_200crit_" + simshort + '/' + filename + '.' + snap_step)
            z = float(re.findall(r'z=.* t',str(timestep))[0][2:-2])
            a_s = 1/(1 + z)
            Xc, Yc, Zc, Rvir, mass, halo_id = timestep.calculate_all('Xc','Yc','Zc','Rvir','Mvir',grp_key)
            if amiga_dat:
                h = 1 # Because the hubble constant is already accounted for in Mvir, Rvir, Xc etc.
            Xc = Xc/h; Yc = Yc/h; Zc = Zc/h; Rvir = Rvir/h; mass = mass/h
            
            loc = np.array([Xc[mass > minmass], Yc[mass > minmass], Zc[mass > minmass]]).T
            rvir = Rvir[mass > minmass]
            mvir = mass[mass > minmass]
            a = a_s + np.zeros(np.sum([mass > minmass]))
            ids = halo_id[mass > minmass]
            n_halos = len(Xc)
        except:
            objs_pd_200.iloc[i]['mpb_close_grp'][snap_step] = 0
            objs_pd_200.iloc[i]['mpb_tindex'][snap_step] = 0
            objs_pd_200.iloc[i]['mpb_tindex_grp'][snap_step] = 0
            objs_pd_200.iloc[i]['mpb_tindex_mass'][snap_step] = 0
            objs_pd_200.iloc[i]['mpb_massive_dist'][snap_step] = 0
            continue
    else:
        loc = []
        rvir = []
        mvir = []
        a = []
        ids = []
        if read_AHF: #!read_stat_data:
            sim = pynbody.load(filepath + 'snapshots_200crit_' + simshort + '/' + filename + '.' + snap_step)
            h_dummy = sim.halos(dummy = True)
        else:
            stat_file_path = M200_dir + snap_step + '/' + filename + '.' + snap_step + '.M200.amiga.stat'
            data = task(
                np.genfromtxt,
                start_text=f"Loading stat file from {stat_file_path}",
                end_text=f"Loaded stat file from {stat_file_path}",
                fail_text=f"Failed to load stat file from {stat_file_path}",
                exit_on_fail=True
            )(stat_file_path, dtype="float", usecols=(0,1,3,5,6,8,13,14,15,20), skip_header=0)
            grps, Ntots, Nstars, Mvirs, Rvirs, Mstars, xcs, ycs, zcs, sats = np.transpose(data)
            AHFhalo_properties = [{'halo_id': grps[i], 'mass': Mvirs[i], 'Rvir': Rvirs[i], 'Xc': xcs[i], 'Yc': ycs[i], 'Zc': zcs[i]} for i in range(len(grps))]
            h = 1 # Because the hubble constant is already accounted for in Mvir, Rvir, Xc etc.

        # Ignore the least massive halos
        for AHFhalo in h_dummy:
            properties = AHFhalo.properties
            #    for properties in AHFhalo_properties:
            if (properties['mass']/h > minmass):
                loc.append(np.array([properties['Xc']/h, properties['Yc']/h, properties['Zc']/h]))
                rvir.append(properties['Rvir']/h)
                mvir.append(properties['mass']/(h))
                a.append(properties['a'])
                ids.append(properties['halo_id'])

    loc = np.array(loc)
    rvir = np.array(rvir)            
    mvir = np.array(mvir)
    a = np.array(a)
    ids = np.array(ids)
    a_arr.append(a[0])
    
    #row = objs_pd_200.iloc[0]
    #for index, row in objs_pd_200.iterrows():
    mpb_grp = objs_pd_200.iloc[0]['mpb_grp'][snap_step]
    if mpb_grp != 0:
        if read_tangosDB: # Read data from the tangos DB
            properties_main = {
                "Xc": timestep.halos[mpb_grp - 1].calculate("Xc"),
                "Yc": timestep.halos[mpb_grp - 1].calculate("Yc"),
                "Zc": timestep.halos[mpb_grp - 1].calculate("Zc"),
                "mass": timestep.halos[mpb_grp - 1].calculate("Mvir") 
            }
            #if not amiga_dat:
            #    Xc = Xc/h; Yc = Yc/h; Zc = Zc/h; Rvir = Rvir/h; mass = mass/h
        else:
            properties_main =  h_dummy[mpb_grp].properties
    else: properties_main = None
    print("ID, ID (z=0), (mass): Massive ID (mass, dist), T. Index ID (mass, dist, index), Main dist")
    for i in arange(0,len(objs_pd_200)):
        #haloid = row['haloid'] #replace .iloc[i] with [objs_pd_200['haloid'] == haloid]
        mpb_grp = objs_pd_200.iloc[i]['mpb_grp'][snap_step]
        if mpb_grp != 0:
            if read_tangosDB: # Read data from the tangos DB
                if mpb_grp > n_halos:
                    continue # This halo is too small to be tracked by tangos
                properties = {
                    "Xc": timestep.halos[mpb_grp - 1].calculate("Xc"),
                    "Yc": timestep.halos[mpb_grp - 1].calculate("Yc"),
                    "Zc": timestep.halos[mpb_grp - 1].calculate("Zc"),
                    "mass": timestep.halos[mpb_grp - 1].calculate("Mvir") 
                }
            else:            
                properties = h_dummy[mpb_grp].properties
            #print(i, objs_pd_200.iloc[i]['haloid'], mpb_grp, properties['halo_id'])

            if properties_main != None:
                objs_pd_200.iloc[i]['mpb_main_dist'][snap_step] = ((properties['Xc']/h - properties_main['Xc']/h)**2 + (properties['Yc']/h - properties_main['Yc']/h)**2 + (properties['Zc']/h - properties_main['Zc']/h)**2)**(0.5)*a[0]
                #objs_pd_200.iloc[i]['mpb_main_dist'][snap_step] = ((properties['Xc']/h - properties_main['Xc']/properties_main['h'])**2 + (properties['Yc']/h - properties_main['Yc']/properties_main['h'])**2 + (properties['Zc']/h - properties_main['Zc']/properties_main['h'])**2)**(0.5)*a[0]
                

            moremass = (mvir > minmass/h)            
            if np.sum(moremass) > 0:
                dists = ((properties['Xc']/h - loc[moremass,0])**2 + (properties['Yc']/h - loc[moremass,1])**2 + (properties['Zc']/h - loc[moremass,2])**2)**(0.5)*a[0]
                
                tindex = mvir[moremass]/dists**3 # Tidal Index
                #tindex = mvir/np.sqrt((properties['Xc']/h - loc[:,0])**2 + (properties['Yc']/h - loc[:,1])**2 + (properties['Zc']/h - loc[:,2])**2)/a #Grav potential
                tindex_max = np.nanmax(tindex[np.isfinite(tindex)])
                objs_pd_200.iloc[i]['mpb_tindex'][snap_step] = tindex_max
                objs_pd_200.iloc[i]['mpb_tindex_grp'][snap_step] = ids[np.where(tindex == tindex_max)[0][0]]
                objs_pd_200.iloc[i]['mpb_tindex_mass'][snap_step] = mvir[np.where(tindex == tindex_max)[0][0]]
                
                #print('Dist: ',ids[minind],objs_pd_200.iloc[i]['mpb_close_id'][snap_step],objs_pd_200.iloc[i]['mpb_massive_dist'][snap_step],np.min(dists[np.nonzero(dists)]))                
                #print('Grav: ',ids[np.where(tindex == tindex_max)[0][0]],objs_pd_200.iloc[i]['mpb_tindex_id'][snap_step],objs_pd_200.iloc[i]['mpb_tindex'][snap_step],np.nanmax(tindex[np.isfinite(tindex)]))

            moremass = (mvir > mass_threshold) #/h)            
            if np.sum(moremass) > 0: 
                dists = ((properties['Xc']/h - loc[moremass,0])**2 + (properties['Yc']/h - loc[moremass,1])**2 + (properties['Zc']/h - loc[moremass,2])**2)**(0.5)
                objs_pd_200.iloc[i]['mpb_massive_dist'][snap_step] = np.min(dists[np.nonzero(dists)])
                minind = np.where((((properties['Xc']/h - loc[:,0])**2 + (properties['Yc']/h - loc[:,1])**2 + (properties['Zc']/h - loc[:,2])**2)**(0.5)) == objs_pd_200.iloc[i]['mpb_massive_dist'][snap_step])[0][0]
                objs_pd_200.iloc[i]['mpb_close_grp'][snap_step] = ids[minind]
                objs_pd_200.iloc[i]['mpb_massive_dist'][snap_step] = objs_pd_200.iloc[i]['mpb_massive_dist'][snap_step]*a[0]
                
                if ids[minind] != ids[np.where(tindex == tindex_max)[0][0]]:
                    print("{0}, {1} ({2:4.2e}): {3:3d} ({4:4.2e}, {5:4.2e}) vs {6:3d} ({7:4.2e}, {8:4.2e}, {9:4.2e}), {10:4.2e}. ".format(mpb_grp,objs_pd_200.iloc[i]['haloid'],properties['mass'],int(ids[minind]),mvir[minind],objs_pd_200.iloc[i]['mpb_massive_dist'][snap_step],int(ids[np.where(tindex == tindex_max)[0][0]]),mvir[np.where(tindex == tindex_max)[0][0]],(mvir[np.where(tindex == tindex_max)[0][0]]/objs_pd_200.iloc[i]['mpb_tindex'][snap_step])**(1/3.0),objs_pd_200.iloc[i]['mpb_tindex'][snap_step],objs_pd_200.iloc[i]['mpb_main_dist'][snap_step]))
                else:
                    print("{0}, {1} ({2:4.2e}): {3:3d} ({4:4.2e}, {5:4.2e}) == {6:3d} ({7:4.2e}, {8:4.2e}, {9:4.2e}), {10:4.2e}. ".format(mpb_grp,objs_pd_200.iloc[i]['haloid'],properties['mass'],int(ids[minind]),mvir[minind],objs_pd_200.iloc[i]['mpb_massive_dist'][snap_step],int(ids[np.where(tindex == tindex_max)[0][0]]),mvir[np.where(tindex == tindex_max)[0][0]],sqrt(mvir[np.where(tindex == tindex_max)[0][0]]/objs_pd_200.iloc[i]['mpb_tindex'][snap_step])**(1/3.0),objs_pd_200.iloc[i]['mpb_tindex'][snap_step],objs_pd_200.iloc[i]['mpb_main_dist'][snap_step]))
                    #print(snap_step,objs_pd_200.iloc[i]['haloid'],ids[minind],ids[np.where(tindex == tindex_max)[0][0]])
            else:
                print("{0}, {1} ({2:4.2e}):  -- (--, --) vs. {3:3d} ({4:4.2e}, {5:4.2e}, {6:4.2e}) vs, {7:4.2e}. ".format(mpb_grp,objs_pd_200.iloc[i]['haloid'],properties['mass'],int(ids[np.where(tindex == tindex_max)[0][0]]),mvir[np.where(tindex == tindex_max)[0][0]],sqrt(mvir[np.where(tindex == tindex_max)[0][0]]/objs_pd_200.iloc[i]['mpb_tindex'][snap_step])**(1/3.0),objs_pd_200.iloc[i]['mpb_tindex'][snap_step],objs_pd_200.iloc[i]['mpb_main_dist'][snap_step]))

"""
            else:
                objs_pd_200.iloc[i]['mpb_dist'][snap_step] = 0
                objs_pd_200.iloc[i]['mpb_close_grp'][snap_step] = 0
                objs_pd_200.iloc[i]['mpb_tindex'][snap_step] = 0
                objs_pd_200.iloc[i]['mpb_tindex_grp'][snap_step] = 0
                objs_pd_200.iloc[i]['mpb_tindex_mass'][snap_step] = 0
#                objs_pd_200.iloc[i]['mpb_massive_dist'][snap_step] = 0
"""

# Calculate the distance when the gravitational force was at the peak and the time that it happened
for i in arange(0,len(objs_pd_200)):
    if np.max(np.array(list(objs_pd_200.iloc[i]['mpb_massive_dist'].values()))) != 0:
        dist_arr = np.array(list(objs_pd_200.iloc[i]['mpb_massive_dist'].values()))
        snap_arr = np.array(list(objs_pd_200.iloc[i]['mpb_massive_dist'].keys()))

        cond = (dist_arr > 0) & (alltimes > 4)
        ind = np.argmin(dist_arr[cond])
        objs_pd_200.iloc[i,objs_pd_200.columns.get_loc('min_dist')] = np.min(dist_arr[cond])
        objs_pd_200.iloc[i,objs_pd_200.columns.get_loc('min_dist_t')] = snap_to_time_code((snap_arr[cond])[ind])
        print("{0:3d}: {1:3.2f}, {2:3.2f}".format(i,objs_pd_200.iloc[i]['min_dist_t'],objs_pd_200.iloc[i]['min_dist']))
    else:
        objs_pd_200.iloc[i]['min_dist_t'] = 0 #objs_pd_200.iloc[i]['time']
        objs_pd_200.iloc[i]['min_dist'] = 0 #objs_pd_200.iloc[i]['massiveDist']    
    
    if np.max(np.array(list(objs_pd_200.iloc[i]['mpb_tindex'].values()))) != 0:
        ind = np.argmax(np.array(list(objs_pd_200.iloc[i]['mpb_tindex'].values()))[alltimes > 4])
        objs_pd_200.iloc[i,objs_pd_200.columns.get_loc('max_tindex')] = (np.array(list(objs_pd_200.iloc[i]['mpb_tindex'].values()))[alltimes > 4])[ind]
        objs_pd_200.iloc[i,objs_pd_200.columns.get_loc('max_tindex_dist')] = ((np.array(list(objs_pd_200.iloc[i]['mpb_tindex_mass'].values()))[alltimes > 4]/np.array(list(objs_pd_200.iloc[i]['mpb_tindex'].values()))[alltimes > 4])[ind])**(1/3.0)
        objs_pd_200.iloc[i,objs_pd_200.columns.get_loc('max_tindex_t')] = snap_to_time_code((np.array(list(objs_pd_200.iloc[i]['mpb_tindex'].keys()))[alltimes > 4])[ind])
        print("{0:3d}: {1:3.2f}, {2:3.2f}, {3:3.2f}".format(i,objs_pd_200.iloc[i]['max_tindex_t'],objs_pd_200.iloc[i]['max_tindex_dist'],objs_pd_200.iloc[i]['max_tindex'])) #np.log10(np.max(np.array(list(objs_pd_200.iloc[i]['mpb_tindex'].values()))[alltimes > 4]))))
    else:
        objs_pd_200.iloc[i]['max_tindex_t'] = 0 # objs_pd_200.iloc[i]['time']
        objs_pd_200.iloc[i]['max_tindex'] = 0 # objs_pd_200.iloc[i]['massiveDist']
        objs_pd_200.iloc[i]['max_tindex_dist'] = 0
    
#objs_pd_200.to_pickle(filepath + filename + '.' + snap_nums[-1] + '/ahf_200/' + filename + '.' + snap_nums[-1] + '.m200.dist.data')
    
for index, row in objs_pd_200.iterrows(): print(row['mpb_tindex_grp']['004096'])

# Distance to the halo that provides a maximum gravitational force
# d = sqrt(Mvir/G_force)
plt.clf()
fig1 = plt.figure(1)
ax1 = fig1.add_subplot(1,1,1)
for index, row in objs_pd_200.iloc[1:].iterrows():
    ts_link = np.array(list(row['mpb_tindex_mass'].values())) != 0
    ax1.plot(alltimes[ts_link],((np.array(list(row['mpb_tindex_mass'].values()))/np.array(list(row['mpb_tindex'].values())))**(1/3.0))[ts_link],alpha = 0.4)
row = objs_pd_200.iloc[0]
ax1.plot(alltimes[ts_link],((np.array(list(row['mpb_tindex_mass'].values()))/np.array(list(row['mpb_tindex'].values())))**(1/3.0))[ts_link],alpha = 0.4)
plt.gca().set_prop_cycle(None)
for index, row in objs_pd_200.iloc[1:].iterrows():
    ax1.plot([row['max_tindex_t'],row['max_tindex_t']],[row['max_tindex_dist'],row['max_tindex_dist']],marker = 'o')  
row = objs_pd_200.iloc[0]
ax1.plot([row['max_tindex_t'],row['max_tindex_t']],[row['max_tindex_dist'],row['max_tindex_dist']],marker = 'o') 
ax1.plot(time,rvir_main,color = 'k')    
#ax.axis([0, 14, 2e5, 3e9])
ax1.set_yscale('log')
ax1.set_xlabel('Time [Gyr]')
ax1.set_ylabel('Physical Distance to Most Impactful Companion')

# Distance to massive halo
fig2 = plt.figure(2)
ax2 = fig2.add_subplot(1,1,1)
for index, row in objs_pd_200.iloc[1:].iterrows():
    ts_link = np.array(list(row['mpb_massive_dist'].values())) != 0
    ax2.plot(alltimes[ts_link],(np.array(list(row['mpb_massive_dist'].values())))[ts_link],alpha = 0.4)
row = objs_pd_200.iloc[0]
ax2.plot(alltimes[ts_link],(np.array(list(row['mpb_massive_dist'].values())))[ts_link],alpha = 0.4)
plt.gca().set_prop_cycle(None)
for index, row in objs_pd_200.iloc[1:].iterrows():
    ax2.plot([row['min_dist_t'],row['min_dist_t']],[row['min_dist'],row['min_dist']],marker = 'o')  
#row = objs_pd_200.iloc[22]
#ax2.plot([row['min_dist_t'],row['min_dist_t']],[row['min_dist'],row['min_dist']],marker = 'o') 
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
ts_link = np.array(list((((objs_pd_200.loc[objs_pd_200['haloid']=='8','mpb_tindex_mass']).values)[0]).values()))!=0
ax4.plot(alltimes[ts_link],(np.array(list((objs_pd_200.loc[objs_pd_200['haloid']=='8','mpb_tindex'].values)[0].values())))[ts_link] ,alpha = 0.4, color = 'k')
ax4.plot([row['max_tindex_t'],row['max_tindex_t']],[row['max_tindex'],row['max_tindex']],marker = 'o') 
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
