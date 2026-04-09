#!/usr/bin/env python

import pandas as pd
import tangos
import tangos.examples.mergers as mergers

import numpy as np
import h5py
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
h = 0.6776942783267969

from tangos_halo_module.halo_properties import track_halo_property, get_timesteps, ID_to_sim_halo_snap, infall_final_n_particles, infall_final_coordinates, apocentric_distance, disruption_time, accretion_time, orbit_interpolation, infall_velocity, quenching_time, max_sSFR_time, max_mass_time 
from tangos_halo_module.path import get_file_path, get_halo_snap_num, read_file
from tangos_halo_module.halos import ID_to_tangos_halo, get_survivors, get_main_progenitor_branch, get_zombies, get_host, get_survivor_IDs, get_zombie_IDs, blockPrint, enablePrint, tangos_to_pynbody_halo

ids = np.loadtxt('ids100.txt')
ids=[int(i) for i in ids]
print(ids)

# Initialize dataframe
column_names = ['ID', 'IDs', 'SFR', 'sSFR', 'stellar_mass', 'new_stellar_mass']
df = pd.DataFrame(columns = column_names)
df.to_csv('sSFR_data100.csv', index=False)

# Load routine for getting tracks
def sSFR_track(ID=0, simulation=0, status=0, tangos_halo=0, halo_id=0, snap_num=0, resolution=1000):
    
#     from tangos_halo_module.halos import tangos_to_pynbody_halo, get_main_progenitor_branch

    if ID != 0:
        print('Started ', ID)
        simulation, status, halo_id, snap_num = ID_to_sim_halo_snap(ID=ID)
        print(simulation, status, halo_id, snap_num)
        tangos_halo = ID_to_tangos_halo(ID=ID, resolution=resolution)
    
    # Satellite MPB ids
    MPB, MPB_ids = get_main_progenitor_branch(tangos_halo=tangos_halo, simulation=simulation, resolution=resolution)
    print(MPB_ids)
    print(MPB)
    # Store and recall
    # Host MPB and ids
    host = get_host(simulation=simulation, resolution=resolution)
    host_MPB, host_ids = get_main_progenitor_branch(tangos_halo=host, simulation=simulation, 
                                                    resolution=resolution)

    # Only consider halo while separate from host
    MPB_ids[MPB_ids==host_ids] = 0
    print(MPB_ids)

    for i in range(len(MPB_ids)):
        if MPB_ids[i]!=0:
            if len(str(MPB[i]))>6:
                MPB_ids[i]=get_halo_snap_num(tangos_halo=MPB[i])[2]
            else: 
                MPB_ids[i]=0
    print('MPB new IDs: ', MPB_ids)
    IDs_track = np.zeros(0)
    SFR_track = np.zeros(0)
    sSFR_track = np.zeros(0)
    mass_track = np.zeros(0)
    new_mass_track = np.zeros(0)

    for index, idx in enumerate(MPB_ids):
        print(idx)
        if idx==0:
            IDs_track = np.concatenate((IDs_track, [-1]), axis=None)
            SFR_track = np.concatenate((SFR_track, [-1]), axis=None)
            sSFR_track = np.concatenate((sSFR_track, [-1]), axis=None)
            mass_track = np.concatenate((mass_track, [-1]), axis=None)
            new_mass_track = np.concatenate((new_mass_track, [-1]), axis=None)
        elif len(str(MPB[index]))<6:
            IDs_track = np.concatenate((IDs_track, [idx]), axis=None)
            SFR_track = np.concatenate((SFR_track, [-1]), axis=None)
            sSFR_track = np.concatenate((sSFR_track, [-1]), axis=None)
            mass_track = np.concatenate((mass_track, [-1]), axis=None)
            new_mass_track = np.concatenate((new_mass_track, [-1]), axis=None)
        else:
            tang_halo = ID_to_tangos_halo(ID=idx, resolution=resolution)
            print(tang_halo)
            pyn_halo = tangos_to_pynbody_halo(tangos_halo=tang_halo, simulation=simulation, resolution=resolution)
            mass = np.asarray(pyn_halo.s['mass'].in_units('Msol'))
            age = pyn_halo.properties['time'].in_units('Gyr') - pyn_halo.s['tform'].in_units('Gyr') #SimArray in Gyr
            new_mass = np.sum(mass[age <= 0.1]) # Formed within last 100 Myr
            total_mass = np.sum(mass)
            SFR = new_mass/1e8
            sSFR = SFR/total_mass
            IDs_track = np.concatenate((IDs_track, [idx]), axis=None)
            SFR_track = np.concatenate((SFR_track, [SFR]), axis=None)
            sSFR_track = np.concatenate((sSFR_track, [sSFR]), axis=None)
            mass_track = np.concatenate((mass_track, [total_mass]), axis=None)
            new_mass_track = np.concatenate((new_mass_track, [new_mass]), axis=None)
    print(SFR_track, sSFR_track, mass_track, new_mass_track)
    dat = {'ID': ID,
           'IDs': [IDs_track],
           'SFR': [SFR_track], 
           'sSFR': [sSFR_track], 
           'stellar_mass': [mass_track], 
           'new_stellar_mass': [new_mass_track],}
    df = pd.DataFrame(data=dat)
    return df
    
df = pd.read_csv('sSFR_data100.csv')
completed_ids = np.asarray(df['ID'])
print('Completed Count before run starts: ', len(completed_ids))

# Set counters and run on loop
count=0
done=0
while done<124:
    try:
        print('Watch me go for the ', count, 'th time...')
        for idx in ids:
            if idx in completed_ids:
                print('Already Exists ', idx)
                done=len(completed_ids)
                print('Completed Count --> ', done)
                pass
            else:
                print('Started ', idx)
                sat_df = sSFR_track(ID=idx, resolution=100)
                df = pd.concat([df, sat_df])
                df.to_csv('sSFR_data100.csv', index=False)
                completed_ids = np.concatenate((completed_ids, [idx]), axis=None)
                print('Finished ', idx)
                done=len(completed_ids)
                print('Completed Count --> ', done)
    except: 
        count+=1
        continue
        
print('Finally, we are done!')