import pynbody
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
import sys
sys.path.append('/home/christenc/Code')

simkey = "h329"
base_lr = simkey + ".cosmo50PLK.3072gst"
base_hr = simkey + ".cosmo50PLK.6144g"
gas_key = "5"

simkey = "h148"
base_lr = simkey + ".cosmo50PLK.3072g"
base_hr = simkey + ".cosmo50PLK.6144g"
gas_key = "3"

file_path_lr = "/data/REPOSITORY/e12Gals/" + base_lr + gas_key + "HbwK1BH/"
file_path_hr = "/data/REPOSITORY/e12Gals/" + base_hr + gas_key + "HbwK1BH/"

sim_lr_file_z0 = file_path_lr + base_lr + gas_key + "HbwK1BH.004096/ahf_200/" + base_lr + gas_key + "HbwK1BH.004096"
sim_hr_file_z0 = file_path_hr + base_hr + gas_key + "HbwK1BH.004096/ahf_200/" + base_hr + gas_key + "HbwK1BH.004096"

sim_lr_z0 = pynbody.load(sim_lr_file_z0)
sim_hr_z0 = pynbody.load(sim_hr_file_z0)
pynbody.config['halo-class-priority'] = [pynbody.halo.ahf.AHFCatalogue] 

h_lr_dummy = sim_lr_z0.halos(dummy = True)
h_lr = sim_lr_z0.halos()
h_hr_dummy = sim_hr_z0.halos(dummy = True)
h_hr = sim_hr_z0.halos()

match_arr = np.loadtxt(file_path_lr + simkey + "_mint_nmint_halos.txt", dtype = int)

lr_halos = match_arr[match_arr[:,1] > 0, 0] # Remove halos for which there is no match in the high res run
hr_halos = match_arr[match_arr[:,1] > 0, 1]

lr_halos_df = pd.DataFrame({'HaloID':lr_halos})
hr_halos_df = pd.DataFrame({'HaloID':hr_halos})

#hr_halos_df['HaloID'] = lr_halos
#lr_halos_df['HaloID'] = hr_halos
hr_halos_df['Mstar'] = np.empty(len(hr_halos))
lr_halos_df['Mstar'] = np.empty(len(lr_halos))
hr_halos_df['Mvir'] = np.empty(len(hr_halos))
lr_halos_df['Mvir'] = np.empty(len(lr_halos))
hr_halos_df['fMhires'] = np.empty(len(hr_halos))
lr_halos_df['fMhires'] = np.empty(len(lr_halos))
#hr_halos_df['SFH'] = np.empty(len(hr_halos))
#lr_halos_df['SFH'] = np.empty(len(lr_halos))

#lr_halos_df = lr_halos_df.set_index(0)
#hr_halos_df = hr_halos_df.set_index(0)

lr_SFH = []
for index, row in lr_halos_df.iterrows():
    lr_halos_df.loc[index,'fMhires'] = h_lr[row['HaloID']].properties['fMhires']    
    lr_halos_df.loc[index,'Mvir'] = h_lr[row['HaloID']].properties['mass']
    lr_halos_df.loc[index,'Mstar'] = h_lr[row['HaloID']].properties['M_star']
    halo = h_lr.load_copy(row['HaloID'])
    if len(halo.star) > 0:
        sfhist = plt.hist(halo.star['tform'].in_units('Gyr'), weights = halo.star['massform'].in_units('Msol'), range = (0,13.8), bins = 138)
        lr_SFH.append(sfhist[0])
    else:
        lr_SFH.append(np.empty(138)*0)
    #lr_halos_df.loc[index,'SFH'] = sfhist

lr_halos_df['sfh'] = lr_SFH
pickle.dump(lr_halos_df, open( simkey + "_matched_nmint.pkl", "wb" ) )
hr_SFH = []
for index, row in hr_halos_df.iterrows():
    hr_halos_df.loc[index,'fMhires'] = h_hr[row['HaloID']].properties['fMhires']
    hr_halos_df.loc[index,'Mvir'] = h_hr[row['HaloID']].properties['mass']
    hr_halos_df.loc[index,'Mstar'] = h_hr[row['HaloID']].properties['M_star']
    halo = h_hr.load_copy(row['HaloID'])
    if len(halo.star) > 0:
        sfhist = plt.hist(halo.star[halo.star['mass'] > 0]['tform'].in_units('Gyr'), weights = halo.star[halo.star['mass'] > 0]['massform'].in_units('Msol'), range = (0,13.8), bins = 138)
        hr_SFH.append(sfhist[0])
    else:
        hr_SFH.append(np.empty(138)*0)
    #hr_halos_df.loc[index,'SFH'] = sfhist

hr_halos_df['sfh'] = hr_SFH
pickle.dump(hr_halos_df, open( simkey + "_matched_mint.pkl", "wb" ) )
fig1 = plt.figure(1)
ax = plt.subplot(111)
ax.plot(lr_halos_df[lr_halos_df['fMhires'] > 0.9]['Mvir'], np.array(hr_halos_df[lr_halos_df['fMhires'] > 0.9]['Mvir'])/np.array(lr_halos_df[lr_halos_df['fMhires'] > 0.9]['Mvir']), "o")
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlabel("NM Halo Mass")
ax.set_ylabel("(M Halo Mass)/(NM Malo Mass)")
fig1.show()
fig1.savefig(simkey + "_Mint_NMint_hmass.png")

fig2 = plt.figure(2)
ax = plt.subplot(111)
ax.plot(lr_halos_df[lr_halos_df['fMhires'] > 0.9]['Mvir'], np.array(hr_halos_df[lr_halos_df['fMhires'] > 0.9]['Mstar'])/np.array(lr_halos_df[lr_halos_df['fMhires'] > 0.9]['Mstar']), "bo")
ax.plot(lr_halos_df[(lr_halos_df['fMhires'] > 0.9) & (lr_halos_df['Mstar'] == 0)]['Mvir'],np.empty(np.sum([(lr_halos_df['fMhires'] > 0.9) & (lr_halos_df['Mstar'] == 0)])) + 9 , "k^")
ax.plot(lr_halos_df[(lr_halos_df['fMhires'] > 0.9) & (hr_halos_df['Mstar'] == 0)]['Mvir'],np.empty(np.sum([(lr_halos_df['fMhires'] > 0.9) & (hr_halos_df['Mstar'] == 0)])) + 0.4 , "kv")
ax.plot([1e7,1e11],[1,1],'k')
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlabel("NM Halo Mass")
ax.set_ylabel("(M Mstar)/(NM Mstar)")
ax.set_xlim(1e7,1e11)
fig2.show()
fig2.savefig(simkey + "_Mint_NMint_smass.png")

fig3 = plt.figure(3, figsize=(8,10.5))
if np.sum(lr_halos_df['fMhires'] > 0.9) > 36:
   gs = gridspec.GridSpec(int(np.ceil(np.sum(lr_halos_df['fMhires'] > 0.9)/3)),3)
elif np.sum(lr_halos_df['fMhires'] > 0.9) > 18:
   gs = gridspec.GridSpec(int(np.ceil(np.sum(lr_halos_df['fMhires'] > 0.9)/2)),2)
else: 
    gs = gridspec.GridSpec(np.sum(lr_halos_df['fMhires'] > 0.9),1)
ct = 0
for index, row in lr_halos_df.iterrows():
    if row['fMhires'] < 0.9:
        continue
    if (row['Mstar'] == 0) and (hr_halos_df.loc[index, 'Mstar'] == 0):
        continue
    ax = fig3.add_subplot(gs[ct])
    sfh_nmint = row['sfh']/1e8*100
    sfh_nmint[sfh_nmint < 0] = 0
    sfh_mint = hr_halos_df.loc[index, 'sfh']/1e8*100
    sfh_mint[sfh_mint < 0] = 0
    row['sfh'][row['sfh'] < 0] = 0
    if hr_halos_df.loc[index, 'Mstar'] > 0:
        ax.plot(np.linspace(0, 13.8, num = 138, endpoint=False) + 0.05, sfh_mint, color = 'darkorange')
    else:
        sfh_mint = sfh_mint*0
    if row['Mstar'] > 0:
        ax.plot(np.linspace(0, 13.8, num = 138, endpoint=False) + 0.05, sfh_nmint,'--',color = 'cornflowerblue')
    else:
        sfh_nmint = sfh_nmint*0
    ax.text(8, 0.65*max(np.nanmax(sfh_mint),np.nanmax(sfh_nmint)), "{}, {}".format(row['HaloID'], hr_halos_df.loc[index, 'HaloID']))
    print(row['HaloID'],(np.nanmax(sfh_mint),np.nanmax(sfh_nmint)))
    ct = ct + 1
fig3.add_subplot(111, frame_on=False)
plt.tick_params(labelcolor="none", bottom=False, left=False)
plt.xlabel("Time [Gyr]")
plt.ylabel("SFR [100*Msol/yr]")
plt.title("Mint (orange) vs Near Mint (blue)")
#plt.tight_layout()
fig3.show()
fig3.savefig(simkey + "_Mint_NMint_SFH.png")
    
fig4 = plt.figure(4)
ax = plt.subplot(111)
ax.plot(lr_halos_df[lr_halos_df['fMhires'] > 0.9]['Mvir'], (np.array(hr_halos_df[lr_halos_df['fMhires'] > 0.9]['Mstar'])-np.array(lr_halos_df[lr_halos_df['fMhires'] > 0.9]['Mstar']))/np.array(lr_halos_df[lr_halos_df['fMhires'] > 0.9]['Mstar']), "bo")
ax.plot(lr_halos_df[(lr_halos_df['fMhires'] > 0.9) & (lr_halos_df['Mstar'] == 0)]['Mvir'],np.empty(np.sum([(lr_halos_df['fMhires'] > 0.9) & (lr_halos_df['Mstar'] == 0)])) + 4 , "k^")
#ax.plot(lr_halos_df[(lr_halos_df['fMhires'] > 0.9) & (hr_halos_df['Mstar'] == 0)]['Mvir'],np.empty(np.sum([(lr_halos_df['fMhires'] > 0.9) & (hr_halos_df['Mstar'] == 0)])) + 0.4 , "kv")
ax.plot([1e7,1e11],[0,0],'k')
#ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlabel("NM Halo Mass")
ax.set_ylabel("(M Mstar - NM Mstar)/(NM Mstar)")
ax.set_ylim(-1.5,7)
fig4.show()
fig4.savefig(simkey + "_Mint_NMint_frac_smass.png")

fig5 = plt.figure(5)
ax = plt.subplot(111)
ax.plot(lr_halos_df[lr_halos_df['fMhires'] > 0.9]['Mvir'], (np.array(hr_halos_df[lr_halos_df['fMhires'] > 0.9]['Mvir'])-np.array(lr_halos_df[lr_halos_df['fMhires'] > 0.9]['Mvir']))/np.array(lr_halos_df[lr_halos_df['fMhires'] > 0.9]['Mvir']), "bo")
#ax.plot(lr_halos_df[(lr_halos_df['fMhires'] > 0.9) & (lr_halos_df['Mstar'] == 0)]['Mvir'],np.empty(np.sum([(lr_halos_df['fMhires'] > 0.9) & (lr_halos_df['Mstar'] == 0)])) + 9 , "k^")
#ax.plot(lr_halos_df[(lr_halos_df['fMhires'] > 0.9) & (hr_halos_df['Mstar'] == 0)]['Mvir'],np.empty(np.sum([(lr_halos_df['fMhires'] > 0.9) & (hr_halos_df['Mstar'] == 0)])) + 0.4 , "kv")
ax.plot([1e7,1e11],[0,0],'k')
#ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlabel("NM Halo Mass")
ax.set_ylabel("(M Mhalo - NM Mhalo)/(NM Mhalo)")
ax.set_ylim(-1.5,10)
fig5.show()
fig5.savefig(simkey + "_Mint_NMint_frac_hmass.png")

fig6 = plt.figure(6)
ax = plt.subplot(111)
ax.plot(lr_halos_df[lr_halos_df['fMhires'] > 0.9]['Mvir'], (lr_halos_df[lr_halos_df['fMhires'] > 0.9]['Mstar']), "o")
ax.plot(hr_halos_df[lr_halos_df['fMhires'] > 0.9]['Mvir'], (hr_halos_df[lr_halos_df['fMhires'] > 0.9]['Mstar']), "o")
#ax.plot(lr_halos_df[(lr_halos_df['fMhires'] > 0.9) & (lr_halos_df['Mstar'] == 0)]['Mvir'],np.empty(np.sum([(lr_halos_df['fMhires'] > 0.9) & (lr_halos_df['Mstar'] == 0)])) + 9 , "k^")
#ax.plot(lr_halos_df[(lr_halos_df['fMhires'] > 0.9) & (hr_halos_df['Mstar'] == 0)]['Mvir'],np.empty(np.sum([(lr_halos_df['fMhires'] > 0.9) & (hr_halos_df['Mstar'] == 0)])) + 0.4 , "kv")
#ax.plot([1e7,1e11],[0,0],'k')
#ax.set_yscale('log')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel("Halo Mass [Msol]")
ax.set_ylabel("Stellar Mass [Msol]")
fig6.show()
fig6.savefig(simkey + "_Mint_NMint_matchedSMHM.png")

