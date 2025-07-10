# Charlotte Christensen
# February 16th, 2022
# This file compares the merger tree used by Deb (created by Tangos link) to that produced by Elaad's halo_trace.py

import pandas as pd

simshort = "h329_6144"
simshort = "h148_6144"
simshort = "h148"

path = '/home/christenc/Data/Sims/'
trace_file_paths = {
        "cptmarvel": path + 'cptmarvel.cosmo25cmb/cptmarvel.cosmo25cmb.4096g5HbwK1BH/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.M200.trace_back.hdf5'
        ,"rogue": path + 'rogue.cosmo25cmb/rogue.cosmo25cmb.4096g5HbwK1BH/rogue.cosmo25cmb.4096g5HbwK1BH.004096.M200.trace_back.hdf5'
        ,"elektra": path + 'elektra.cosmo25cmb/elektra.cosmo25cmb.4096g5HbwK1BH/elektra.cosmo25cmb.4096g5HbwK1BH.004096.M200.trace_back.hdf5'
        ,"storm": path + 'storm.cosmo25cmb/storm.cosmo25cmb.4096g5HbwK1BH/storm.cosmo25cmb.4096g5HbwK1BH.004096.M200.trace_back.hdf5'
        ,"h148_6144": path + 'h148.cosmo50PLK.6144g/h148.cosmo50PLK.6144g3HbwK1BH/h148.cosmo50PLK.6144g3HbwK1BH.004096.trace_back.hdf5'
        ,"h329_6144": path + 'h329.cosmo50PLK.6144g/h329.cosmo50PLK.6144g5HbwK1BH/h329.cosmo50PLK.6144g5HbwK1BH.004096.trace_back.hdf5'
        ,"h148": path + 'h148.cosmo50PLK.3072g/h148.cosmo50PLK.3072g3HbwK1BH/h148.cosmo50PLK.3072g3HbwK1BH.004096.M200.trace_back.hdf5'        
        ,"h229": path + 'h229.cosmo50PLK.3072g/h229.cosmo50PLK.3072gst5HbwK1BH/h229.cosmo50PLK.3072gst5HbwK1BH.004096.M200.trace_back.hdf5'
        ,"h242": path + 'h242.cosmo50PLK.3072g/h242.cosmo50PLK.3072gst5HbwK1BH/h242.cosmo50PLK.3072gst5HbwK1BH.004096.M200.trace_back.hdf5'
        ,"h329": path + 'h329.cosmo50PLK.3072g/h329.cosmo50PLK.3072gst5HbwK1BH/h329.cosmo50PLK.3072gst5HbwK1BH.004096.M200.trace_back.hdf5'
}

sim_trace = pd.read_hdf( trace_file_paths[simshort] )

path_data = "/home/christenc/Code/Datafiles/"
all_sat = pd.read_csv(path_data + simshort + "MPB_ids.csv")

all_sat = all_sat.set_index('snapshot')
all_sat = all_sat.transpose()

survive_sat = all_sat[all_sat[4096] != 0]

# All traced halos
# sim_trace.index.tolist()

for index, row in survive_sat.iterrows():
    print("====",str(row[4096])[7:],"====")
    try:
        halo_trace = sim_trace.loc[int(str(row[4096])[7:])]
        for key, value in halo_trace.items():
            if str(int(value)) == str(row[int(key)])[7:]:
                print(key,int(value),str(row[int(key)])[7:], "Match")
            else:
                print(key,int(value),str(row[int(key)])[7:])
    except:
        print("Missing ", str(row[4096])[7:])
