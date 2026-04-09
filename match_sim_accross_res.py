import pynbody
import numpy as np
import sys
sys.path.append('/home/christenc/Code')
import scipy.spatial
import pickle
import os

simkey = "h329"
base_lr = simkey + ".cosmo50PLK.3072gst"
base_hr = simkey + ".cosmo50PLK.6144g"
gas_key = "5"
source_halos = [1308, 859, 449, 174, 154, 127, 117, 115, 92, 90, 53, 48, 38, 37, 30, 29, 26, 18, 15, 14, 13, 10, 7]

simkey = "h148"
base_lr = simkey + ".cosmo50PLK.3072g"
base_hr = simkey + ".cosmo50PLK.6144g"
gas_key = "3"
source_halos = [1291, 1084, 1067, 925, 801, 699, 647, 637, 603, 557, 447, 438, 358, 341, 336, 327, 282, 251, 249, 230, 188, 148, 131, 128, 122, 114, 109, 94, 86, 75, 72, 65, 59, 55, 52, 51, 43, 41, 38, 37, 35, 34, 33, 31, 29, 28, 27, 23, 20, 15, 13, 12, 11, 10, 7, 6, 4, 3, 2]
source_halos = [282, 65]

file_path_lr = "/data/REPOSITORY/e12Gals/" + base_lr + gas_key + "HbwK1BH/"
file_path_hr = "/data/REPOSITORY/e12Gals/" + base_hr + gas_key + "HbwK1BH/"

#sim_lr_file = file_path_lr + base_lr + "5HbwK1BH.000071/" + base_lr + "5HbwK1BH.000071"
#sim_hr_file = file_path_hr + base_hr + "5HbwK1BH.000071/" + base_hr + "5HbwK1BH.000071"

sim_hr_file = file_path_hr + base_hr + ".tbin"
sim_hr = pynbody.load(sim_hr_file)
hr_mass, hr_mass_ct = np.unique(sim_hr.dark['mass'], return_counts = True)
sim_hr['iord'] = np.arange(len(sim_hr))
pot_match = sim_hr.dark[sim_hr.dark['mass'] == hr_mass[0]]


#-----------------------------------------
#A = np.array([sim_lr.dark[unmatched]['x'].tolist(),sim_lr.dark[unmatched]['y'].tolist(),sim_lr.dark[unmatched]['z'].tolist()]).T

if not os.path.exists(file_path_hr + "tree.pkl"):
     print("Creating tree")
     B = np.array([pot_match['x'].tolist(),pot_match['y'].tolist(),pot_match['z'].tolist()]).T
     tree = scipy.spatial.KDTree(B)
     with open(file_path_hr + "tree.pkl", "wb") as f:
          pickle.dump(tree, f)
else:
     print("Reading tree")
     with open(file_path_hr + "tree.pkl", "rb") as f:
          tree = pickle.load(f)


#-----------------------------------------

sim_lr_file = file_path_lr + base_lr + "st.tbin"
sim_lr = pynbody.load(sim_lr_file)
lr_mass, lr_mass_ct = np.unique(sim_lr.dark['mass'], return_counts = True)
sim_lr['iord'] = np.arange(len(sim_lr))
"""
sim_lr['matchiord'] = sim_lr['iord']*0
for mass in lr_mass[2:]: #Doesn't start at zero because the lr is souped
     (sim_lr.dark[sim_lr.dark['mass'] == mass])['matchiord'] = (sim_hr.dark[sim_hr.dark['mass'] == mass])['iord']
"""
#unmatched = pynbody.filt.LowPass('matchiord', 1)

sim_lr_file_z0 = file_path_lr + base_lr + gas_key + "HbwK1BH.004096/ahf_200/" + base_lr + gas_key + "HbwK1BH.004096"
sim_hr_file_z0 = file_path_hr + base_hr + gas_key + "HbwK1BH.004096/ahf_200/" + base_hr + gas_key + "HbwK1BH.004096"

sim_lr_z0 = pynbody.load(sim_lr_file_z0)
sim_hr_z0 = pynbody.load(sim_hr_file_z0)
pynbody.config['halo-class-priority'] = [pynbody.halo.ahf.AHFCatalogue] # Read the AHF halo cataloge
h_lr_dummy = sim_lr_z0.halos(dummy = True)
h_lr = sim_lr_z0.halos()

ct = 0
for halo in h_lr_dummy:
    ct = ct + 1
    try:
        if halo.properties['n_star'] > 1:
             #            dm = h_lr.load_copy(ct).dark
             print(ct, halo.properties['#ID'], halo.properties['npart'], halo.properties['n_star']) #, len(dm))    
    except:
        continue

match_halos = []

#match_halos = [543, ]

# For each halo of interest (e.g., those other than halo 1 with stars) in the lr sim
    # Determine the iords of the dm particles in that halo
    # Use those iords to select those particles from the lr tbin file
    # Loop through those particles
         # Match particles to the iords from the high res sim using the KD tree
         # Append to list
    # Use matched iords to select particles in the z = 0 high res simulation

for halo_id in source_halos:
     halo = h_lr[halo_id]
     dm = halo.dark
     halo_parts = sim_lr[dm['iord'].tolist()]
     A = np.array([halo_parts['x'].tolist(),halo_parts['y'].tolist(),halo_parts['z'].tolist()]).T
     print(len(A))

     assoc=np.array([])
     i = 0
     for I1,point in enumerate(A):
          if (i%1000)==0:
               print("{:.2f} ".format(i/len(A)), end='')
          i = i+1
          temp, I2 = tree.query(point,k=1)
          # note that B[I2] should be the same as tree.data[I2,:] and np.array([float(pot_match[I2]['x']), float(pot_match[I2]['y']), float(pot_match[I2]['z'])])
          if len(assoc) == 0:
               assoc = [int(halo_parts[I1]['iord']), int(pot_match[I2]['iord'])]
          else:
               assoc = np.vstack([assoc, [int(halo_parts[I1]['iord']), int(pot_match[I2]['iord'])]])

     len_gas = len(sim_hr.gas)
     len_gas_z0 = len(sim_hr_z0.gas)
     hr_part_z0 = sim_hr_z0[len_gas_z0 + np.sort(assoc[:,1]) - len_gas]
     hr_part = sim_hr[np.sort(assoc[:,1])]

     uniq, uniq_cts = np.unique(hr_part_z0['amiga.grp'], return_counts = True)
     print("\n", uniq, uniq_cts)
     if np.sum(np.where(uniq > 1)) > 0:
          match_halo = uniq[(np.where(uniq > 1)[0])[np.argmax(uniq_cts[np.where(uniq > 1)])]]
     else:
          match_halo = -1

     match_halos.append(match_halo)

np.savetxt(simkey + "_matched_halos.txt", np.vstack((np.array(source_halo),np.array(match_halo))).T,fmt="%d")

# Check how well the selection works
"""
     result_lr_IC, bins_IC = np.histogram(halo_parts.dark['x'].in_units('kpc'))
     result_hr_IC, bins_IC = np.histogram(hr_part.dark['x'].in_units('kpc'), bins=bins_IC)
     result_lr_z0, bins_z0_lr = np.histogram(halo.dark['x'].in_units('kpc'))
     result_hr_z0, bins_z0_hr = np.histogram(hr_part_z0.dark['x'].in_units('kpc'))
     print(bins_IC, result_lr_IC, result_hr_IC)
     print(bins_z0_lr, bins_z0_hr, result_lr_z0, result_hr_z0)
"""


