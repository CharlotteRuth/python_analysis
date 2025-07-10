# Retrieve metallicity information for Deb

import tangos
import matplotlib.pyplot as plt
import pandas
import numpy as np

halolist = pandas.read_csv("~/infall_halo_IDs.csv")
#export TANGOS_DB_CONNECTION=/home/christenc/Storage/tangos_db/JL_r200_N100min.db
halolist = pandas.read_csv("/home/christenc/Code/Datafiles/HaloFiles_Deb/Data100_new.csv")
halolist = pandas.read_csv("/home/christenc/Code/Datafiles/HaloFiles_Deb/Data100_new_edited.csv")
#export TANGOS_DB_CONNECTION=/home/christenc/Storage/tangos_db/JLmint_r200_N100.db
#halolist = pandas.read_csv("/home/christenc/Code/Datafiles/HaloFiles_Deb/Datamint_new_offline_wFeHfix_edited.csv")


meanOxH = np.zeros(len(halolist['ID']))
meanOxFe = np.zeros(len(halolist['ID']))
meanFeH = np.zeros(len(halolist['ID']))
massstar = np.zeros(len(halolist['ID']))

tangos.all_simulations() # list all simulations

#key='h329'
#simname = "snapshots_200crit_"+key+"mint" # Stored on emu; Define the simulation you want to look at
#current_sim = key

for index, data in halolist.iterrows():
        simname = "snapshots_200crit_" + data['Simulation']#+"mint"
        #halo = tangos.get_halo(simname + "/%" + str(data['infall_snap_num'])+"/halo_" + str(data['infall_halo_ID']))
        #simname = "snapshots_200crit_" + data['Simulation']
        snap = str(data['ID'])[3:7]
        halo_id = str(data['ID'])[7:]
        halo = tangos.get_halo(simname + "/%" + snap + "/halo_" + halo_id)
        meanOxH[index] = halo['meanOxH']
        meanOxFe[index] = halo['meanOxFe']
        meanFeH[index] = halo['meanFeH']
        massstar[index] = halo['M_star']
        
halolist['meanOxH'] = meanOxH
halolist['meanFeH'] = meanFeH
halolist['meanOxFe'] = meanOxFe
halolist['mass_star2']=massstar

print(massstar)
print(meanFeH)
#halolist.to_csv("~/infall_halo_IDs_update.csv")
#halolist.to_csv("/home/christenc/Code/Datafiles/HaloFiles_Deb/Datamint_new_edited_Z.csv")
halolist.to_csv("/home/christenc/Code/Datafiles/HaloFiles_Deb/Data100_new_edited_Z.csv")

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.plot(halolist[halolist['Status'] == 'Survivor']['meanFeH'],halolist[halolist['Status'] == 'Survivor']['meanOxFe'],'ko')
ax1.plot(halolist[halolist['Status'] == 'Zombie']['meanFeH'],halolist[halolist['Status'] == 'Zombie']['meanOxFe'],'ro')
ax1.set_xlabel("[Fe/H]")
ax1.set_ylabel("[O/Fe]")
ax1.set_xlim([-6,0.1])
ax1.set_ylim([0,0.75])
plt.savefig("/home/christenc/alpha_fe.png")
plt.show()
