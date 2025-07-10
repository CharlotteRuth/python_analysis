# Charlotte Christensen
# 3/2/22
# This program identifies the iords of all gas particles that are in the progenitors at the time of accretion onto the main halo


#mpl.use('tkagg') #Also can try mpl.use('Agg') #for using mpl over ssh    
import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
import pynbody
import numpy as np
import pandas as pd
import socket, sys, os, glob, pickle
sys.path.append(os.path.abspath("/home/christenc/Code/python/python_analysis/"))

datapath = "/home/christenc/Code/students/DebPathak/mergers/"
simpath = "/home/christenc/Data/Sims/"

halo_data = pd.read_csv(datapath + "Data100.csv")

steps = [str(row['infall_ID'])[3:7] for index, row in halo_data.iterrows()]
haloids = [str(row['infall_ID'])[7:] for index, row in halo_data.iterrows()]
halo_data['infall_step'] = steps
halo_data['infall_haloid'] = haloids

filenames = {'h148': simpath + "h148.cosmo50PLK.3072g/h148.cosmo50PLK.3072g3HbwK1BH/snapshots_200crit_h148/h148.cosmo50PLK.3072g3HbwK1BH",
                 'h229': simpath + "h229.cosmo50PLK.3072g/h229.cosmo50PLK.3072gst5HbwK1BH/snapshots_200crit_h229/h229.cosmo50PLK.3072gst5HbwK1BH",
                 'h242': simpath + "h242.cosmo50PLK.3072g/h242.cosmo50PLK.3072gst5HbwK1BH/snapshots_200crit_h242/h242.cosmo50PLK.3072gst5HbwK1BH",
                 'h329': simpath + "h329.cosmo50PLK.3072g/h329.cosmo50PLK.3072gst5HbwK1BH/snapshots_200crit_h329/h329.cosmo50PLK.3072gst5HbwK1BH"
                 }

gas_iord = []
gas_temp = []
gas_rho = []
gas_halo = []
gas_mass = []
gas_OxMassFrac = []
gas_FeMassFrac = []
gas_uniqID = []


# simshort = "h148"
# halos = pd.read_csv(datapath + "Halo_Files/MPB_ids/" + simshort + "MPB_ids.csv")

for simname in halo_data['Simulation'].unique():
    print("Simulation: ", simname)
    for step in sort(halo_data[halo_data['Simulation'] == simname]['infall_step'].unique()):
        print("Step: ", step, sum((halo_data["Simulation"] == simname) & (halo_data["infall_step"] == step)))
        filename = glob.glob(filenames[simname] + ".00" + step)
        sim = pynbody.load(filename[0])
        sim.physical_units()
        h = sim.halos()
        for index, row in halo_data[(halo_data["Simulation"] == simname) & (halo_data["infall_step"] == step)].iterrows():
            haloid = row['infall_haloid']
            uniqID = row['ID']
            print("Halo: ", haloid)
            halo = h.load_copy(int(haloid))
            if len(halo.gas) != 0:
                gas_iord = np.append(gas_iord, halo.gas['iord'])
                gas_temp = np.append(gas_temp, halo.gas['temp'])
                gas_rho = np.append(gas_rho, halo.gas['rho'])
                gas_mass = np.append(gas_mass, halo.gas['mass'])
                gas_OxMassFrac = np.append(gas_OxMassFrac, halo.gas['OxMassFrac'])
                gas_FeMassFrac = np.append(gas_FeMassFrac, halo.gas['FeMassFrac'])                
                gas_halo = np.append(gas_halo, np.full(len(halo.gas), haloid))
                gas_uniqID = np.append(gas_uniqID, np.full(len(halo.gas), uniqID))

tstep_disrupt = np.array([int(str(uniqID)[3:7]) for uniqID in gas_uniqID])
                
gas_part =  pd.DataFrame({'iord' : gas_iord, 'mass' : gas_mass, 'temp' : gas_temp, 'rho' : gas_rho, 'OxMassFrac' : gas_OxMassFrac, 'FeMassFrac' : gas_FeMassFrac, 'infall_halo' : gas_halo, 'uniqID' : gas_uniqID, 'Simulation' : gas_sims, 'tstep_disrupt' : tstep_disrupt})       
gas_part.to_csv("/home/christenc/Code/Datafiles/gas_accr_NearMintJL.csv", index = False)
