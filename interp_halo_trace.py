import time

print("Loading libraries...") #end and flush removed
t1 = time.time()

import sys, os, gc, copy
from scipy.optimize import curve_fit
import numpy as np
import pandas as pd
#import pymp # Open MP
import matplotlib as mpl
#mpl.use('Agg') #for using mpl over ssh
#mpl.use('Qt5Agg')
import matplotlib.pyplot as plt

sys.path.append(os.path.abspath("."))
from modules.user_tools import task, save_table_to_txt #, create_gif
#import modules.cubehelix as cubehelix

t2 = time.time()
print("Libraries loaded in "+str( round(t2-t1, 3) )+"s.")

#clear any unused variables in RAM
gc.collect()

print("Working directory:", os.getcwd())

def time_to_snap_code(time):
    snap_code = str(np.int( np.round( time * 4096/13.8 ) ))
    return (4-len(snap_code))*"0" + snap_code

def snap_code_to_time(snap_code):
    return np.int(snap_code)*(13.8/4096.)

create_assembly_histories = True

if create_assembly_histories:
    stat_data = task(
        np.load,
        start_text="Loading M200 stat data",
        end_text="Loaded M200 stat data",
        fail_text="Failed to load M200 stat data",
        exit_on_fail=True
    )("/home/Code/Datafiles/christenc/consolidated_M200.npy", encoding="latin1", allow_pickle=True).item()

    DCJL_sims = ["h148_6144","h329_6144","h148","h229","h242","h329"]
    Marvel_sims = ["cptmarvel","rogue","elektra","storm"]
    sims_to_include = Marvel_sims + DCJL_sims
    #sims_to_include = ['h148'] #['storm']
    #sims_to_include = ["h148_6144","h329_6144"]

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

    times = []
    volumes = []
    halo_grps = []
    assembly_histories = []
    stellar_mass_histories = []
    gas_mass_histories = []
    Rvir_histories = []
    masses = []
    tracing_fractions = []
    times_per_sim = { sim: [] for sim in sims_to_include }

    for sim in sims_to_include:
        print(f"Getting Mstar and Mpeak for {sim}...")
        sim_trace = pd.read_hdf( trace_file_paths[sim] )

        # satellites_sim = satellites[sim]
        # v_band_mags_sim = v_band_mags[sim]

        halo_grps_sim = list(range(1000)) #list(v_band_mags[sim].keys())

        available_snapshots = list(np.sort(list(stat_data[sim].keys()))[::-1])
        traced_snapshots = list(sim_trace.columns)

        snapshots = [ snapshot for snapshot in traced_snapshots if snapshot in available_snapshots ]
        all_snapshots = [ "4096" ] + [ snapshot[2:] for snapshot in snapshots ]

        # times = [] #np.vectorize(snap_code_to_time)(snapshots)


        for halo_grp in halo_grps_sim:
            try:
                halo_grp_trace = [halo_grp] + list(sim_trace.loc[halo_grp, snapshots])
            except:
                continue

            # print(all_snapshots)
            # print(snapshots)
            # print(halo_grp_trace)
            # sys.exit()

            # Now get the grp from the earliest snapshot

            '''
            final_grp = halo_grp_trace[-1]
            if final_grp > 0 and sim != "elektra":
                if final_grp in inverted_cat.keys():
                    halo_grp_trace.append( int(np.min(inverted_cat[final_grp])) )
                    snaps.append( final_snapshots[sim] )
            '''
            # I couldn't figure out inverted_cat -- CRC

            # Now build Mstar and Mvir histories from stat_data

            time = []
            Mvir = []
            Mstar = []
            Mgas = []
            Rvir = []
            tracing_length = 0
            for i in range(len(all_snapshots)):
                grp = halo_grp_trace[i]
                if grp == -1: continue
                snapshot = all_snapshots[i]
                snap = stat_data[sim][ '00' + snapshot ]
                if grp not in snap.keys():
                    print(f"Halo {grp} is not in stat_data[{sim}][{snapshot}]")
                    continue
                # print(stat_data[sim][ snapshot ])
                Mvir.append( snap[ grp ][0] if grp >= 1 else 0 )
                Mstar.append( snap[ grp ][1] if grp >= 1 else 0 )
                Rvir.append( snap[ grp ][4] if grp >= 1 else 0 )
                Mgas.append( snap[ grp ][7] if grp >= 1 else 0 )
                time.append( snap_code_to_time(snapshot) )
                # if grp > 1:
                #     d_mw.append( float(stat_data[sim][ snapshot ][ grp ][3]) )
                if grp >= 1:
                    tracing_length += 1
            tracing_fraction = float(tracing_length)/len(all_snapshots)
            # time = np.vectorize(snapshot_to_time)(all_snapshots)

            # print(sim, halo_grp, d_mw)

            # Mmax_index = np.argmax(Mvir)
            # snapshot_Mpeak = all_snapshots[Mmax_index]
            # halo_grp_Mpeak = halo_grp_trace[Mmax_index]
            # Mmax = Mvir[Mmax_index]
            # z0_mass = Mvir[0]
            # stellar_mass_z0 = 0.4*Mstar[0]

            times.append( time[::-1] )
            volumes.append( sim )
            halo_grps.append( halo_grp )
            # times.append( time )
            masses.append( Mvir[-1] )
            assembly_histories.append( Mvir[::-1] )
            stellar_mass_histories.append( Mstar[::-1] )
            gas_mass_histories.append( Mgas[::-1] )
            Rvir_histories.append( Rvir[::-1] )
            tracing_fractions.append( tracing_fraction )

    # column_titles = [ "Volume", "halo grp @ z=0", "M_halo @ z=0", "M_peak", "M_star @ z=0", "M_star @ M_peak", "V band mag" ]
    #
    # [ "Volume", "halo grp @ z=0", "halo grp @ M_peak", "snapshot @ M_peak", "M_peak", "M_halo @ z=0", "M_star @ z=0", "M_star,photometric @ z=0", "M_star @ M_peak", "M_star,photometric @ M_peak", "V band mag", "Satellite type" ]
    #
    table = {
        "Time": times,
        "Volume": volumes,
        "halo grp @ z=0": halo_grps,
        # "Age @ peak": peak_ages,
        "Assembly History": assembly_histories,
        "Stellar Mass History": stellar_mass_histories,
        "Gas Mass History": gas_mass_histories,
        "Rvir History": Rvir_histories,
        "M_halo @ z=0": masses,
        "Tracing Fraction": tracing_fractions
    }

    task(
        np.save,
        start_text=f"Saving data",
        end_text=f"Saved data",
        fail_text=f"Failed to save data",
        exit_on_fail=True,
        no_spinner=False
    )("/home/christenc/Code/Datafiles/assembly_histories.npy", table, allow_pickle=True)

    print("Done!")
