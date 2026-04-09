#!/usr/bin/env python

import pandas as pd
import tangos
import tangos.examples.mergers as mergers
import pynbody

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

idx=[str(int(i)) for i in ids]
idx.sort()
h148_ids = idx[:77]
h329_ids = idx[77:]
print(h148_ids)
print(h329_ids)

for sim in ['h329', 'h148']:
    print('started simulation ', sim)
    
    df = pd.read_csv(str(sim)+'MPB_ids.csv')
    IDs = np.asarray(df.columns)[1:]
    print('my IDs are ', IDs)
    ts = get_timesteps(simulation=sim, resolution=100)[3]
    print('the timesteps are ', ts)
    
    
    
    #Read in new dfs
    df_SFR = pd.read_csv(str(sim)+'MPB_SFR.csv')
    df_sSFR = pd.read_csv(str(sim)+'MPB_sSFR.csv')
    df_mass = pd.read_csv(str(sim)+'MPB_mass.csv')
    df_new_mass = pd.read_csv(str(sim)+'MPB_new_mass.csv')
    
    #Iterate through each timestep
    for index in range(0, len(ts)):
        #index is #row
        snapnum = ts[index]
        print('starting with timestep ', snapnum)
        
        halo_ids = np.asarray(df.loc[index][1:])
        print('for timestep ', snapnum, ' my halo IDs are ', halo_ids)
    
        if sim=='h329':
            snapshot = pynbody.load('/home/dp1144/tangos_simulations/'+ str(sim) +'/snapshots/'+ str(sim) +'.cosmo50PLK.6144g5HbwK1BH.00'+ str(snapnum))
        else:
            snapshot = pynbody.load('/home/dp1144/tangos_simulations/'+ str(sim) +'/snapshots/'+ str(sim) +'.cosmo50PLK.6144g3HbwK1BH.00'+ str(snapnum))
        snapshot.physical_units()
        all_halos = snapshot.halos()
        print(snapshot)
        print(all_halos)
        
        count=0 #This is #column-1
        
        #Iterate through each halo
        for idx in halo_ids:
            print('started [', index, ', ', count, ']')
            ID = IDs[count]
            print('my ID is ', ID)
            print('')
            if idx==0:
                #set everything to -1 and move on
                df_SFR.iat[index, count+1] = -1
                df_sSFR.iat[index, count+1] = -1
                df_mass.iat[index, count+1] = -1
                df_new_mass.iat[index, count+1] = -1
            else:
                idx=int(idx)
                simulation, status, halo_id, snap_num = ID_to_sim_halo_snap(ID=idx)
                halo=all_halos[int(halo_id)]
                print(halo)
                #Do the calculation
                mass = np.asarray(halo.s['mass'].in_units('Msol'))
                age = halo.properties['time'].in_units('Gyr') - halo.s['tform'].in_units('Gyr') #SimArray in Gyr
                new_mass = np.sum(mass[age <= 0.1]) # Formed within last 100 Myr
                total_mass = np.sum(mass)
                SFR = new_mass/1e8
                if SFR==0:
                    sSFR=0
                else:
                    sSFR = SFR/total_mass
                
                print('SFR: ', SFR)
                print('sSFR: ', sSFR)
                print('mass: ', total_mass)
                print('new mass: ', new_mass)
                
                #add data in dfs
                df_SFR.iat[index, count+1] = SFR
                df_sSFR.iat[index, count+1] = sSFR
                df_mass.iat[index, count+1] = total_mass
                df_new_mass.iat[index, count+1] = new_mass
            print('finished [', index, ', ', count+1, ']')     

            df_SFR.to_csv(str(sim)+'MPB_SFR.csv', index=False)
            df_sSFR.to_csv(str(sim)+'MPB_sSFR.csv', index=False)
            df_mass.to_csv(str(sim)+'MPB_mass.csv', index=False)
            df_new_mass.to_csv(str(sim)+'MPB_new_mass.csv', index=False)
            
            print('saved [', index, ', ', count, ']')
            
            #Move on to next id == next column
            count+=1
        
        print('finished timestep ', snapnum)
    print('finished simulation ', sim)
print('FINALLY WE ARE DONE!!!')

    
