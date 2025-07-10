import os
import pandas as pd
import pynbody
import tangos
import numpy as np
import sys, os, glob, pickle

import matplotlib.colors as colors
from matplotlib import path
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap
import sys, os, glob, pickle, socket


filepath = '/home/christenc/Data/Sims/h148.cosmo50PLK.3072g/h148.cosmo50PLK.3072g3HbwK1BH/'
filename = filepath.split('/')[-2]
simname = 'Sandra'
simshort = 'h148'

filepath = '/home/christenc/Data/Sims/h229.cosmo50PLK.3072g/h229.cosmo50PLK.3072gst5HbwK1BH/'
filename = filepath.split('/')[-2]
simname = 'Ruth'
simshort = 'h229'

filepath = '/home/christenc/Data/Sims/h242.cosmo50PLK.3072g/h242.cosmo50PLK.3072gst5HbwK1BH/'
filename = filepath.split('/')[-2]
simname = 'Sonia'
simshort = 'h242'

filepath = '/home/christenc/Data/Sims/h329.cosmo50PLK.3072g/h329.cosmo50PLK.3072gst5HbwK1BH/'
filename = filepath.split('/')[-2]
simname = 'Elena'
simshort = 'h329'

#tangos.init_db(filepath + simshort + '.db')
tangos.core.init_db('/home/christenc/Storage/tangos_db/JL_r200.db')

filepath = '/home/christenc/Data/Sims/storm.cosmo25cmb/storm.cosmo25cmb.4096g5HbwK1BH/'
filename = filepath.split('/')[-2]
simname = 'Storm'
simshort = 'storm'

#tangos.init_db(filepath + simshort + '.db')
#tangos.core.init_db('/home/christenc/Storage/tangos_db/Marvel_r200.db')
#tangos.init_db('/home/christenc/Storage/tangos_db/JL_r200.db')
tangos.init_db('/data2/REPOSITORY/tangos.db')

#simname="h329.cosmo50PLK.3072gst"
#simname="h242.cosmo50PLK.3072gst"
#simname="h229.cosmo50PLK.3072gst"
simname="h148.cosmo50PLK.3072gst"
timestep = tangos.get_simulation(simname).timesteps[-1]
#halo_ids = timestep.calculate_all("Grp")
halo_ids = timestep.calculate_all("halo_number()")

mvir_max = np.empty(len(halo_ids[0]))
t_mvir_max = np.empty(len(halo_ids[0]))
halogrp_mvir_max = np.empty(len(halo_ids[0]))
mvirs = np.empty(len(halo_ids[0]))

hubble = 0.6776942783267969

i = 0
for halo in halo_ids[0]:
    tangos_halo = tangos.get_halo(simname+'/%.004096/halo_' + str(halo))
    if not tangos_halo is None:
        time, ids, mvir = tangos_halo.calculate_for_progenitors("t()", "halo_number()", "Mhalo")
        mvir = mvir/hubble
        #time, ids, mvir = tangos_halo.calculate_for_progenitors("t()", "halo_number()", "Mvir")
        if len(mvir) > 0:
            print(halo,time[np.argmax(mvir)],mvir[0]/max(mvir),max(mvir),mvir[0])
            mvir_max[i] = max(mvir)
            t_mvir_max[i] = time[np.argmax(mvir)]
            halogrp_mvir_max[i] = ids[np.argmax(mvir)]
            mvirs[i] = mvir[0]
    i = i + 1

halo_data = np.array([halo_ids[0],halogrp_mvir_max,mvir_max,t_mvir_max,mvirs])
df = pd.DataFrame(data=halo_data.T,index=halo_ids[0],columns=["halogrp_z0","halogrp_Mpeak","Mpeak","t_Mpeak","Mhalo_z0"])
df.to_csv(simname+'_Mpeak.csv',index=False)

timesteps = tangos.get_simulation("snapshots_200crit_" + simshort).timesteps
snap_nums = [re.findall(r'.00+[\d]+', str(snap))[0][3:] for snap in timesteps]
"""
snap_nums_h329 = ('0139', '0188', '0192', '0225', '0275',
       '0288', '0347', '0384', '0456', '0480', '0576', '0637', '0672',
       '0768', '0776', '0864', '0960', '0974', '1056', '1106', '1152',
       '1248', '1269', '1344', '1440', '1475', '1536',
       '1740', '1824', '1920', '2016', '2088', '2112', '2208', '2304',
       '2400', '2496', '2554', '2592', '2688', '2784', '2880', '2976',
       '3072', '3168', '3195', '3264', '3360', '3456', '3552', '3606',
       '3648', '3744', '3936', '4032', '4096')
snap_nums = snap_nums_h329
"""
#time, ids = tangos_halo.calculate_for_progenitors("t()", "halo_number()")
#snap_nums = snap_nums[len(snap_nums) - len(time):]

#xc, yc, zc =

objs_dat = []
f=open(filepath + filename + '.00' + snap_nums[-1] + '/ahf_200/' + filename + '.00' + snap_nums[-1] + '.m200.data', 'rb')
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
    
# This isn't the crit 200 data
objs_dat = []
f=open(filepath + filename + '.00' + snap_nums[-1] + '/' + filename + '.00' + snap_nums[-1] + '.MAP.data', 'rb')
while 1:
    try:
        objs_dat.append(pickle.load(f))
    except EOFError:
        break        
f.close()
if len(objs_dat) == 1:
    objs_pd = pd.DataFrame(objs_dat[0])
else:
    objs_pd = pd.DataFrame(objs_dat)    

tangos_halo = tangos.get_halo("snapshots_200crit_" + simshort + '/' + filename + '.004096/halo_1')
time, scalefactor, id_main, mvir_main, rvir_main = tangos_halo.calculate_for_progenitors("t()", "a()", "halo_number()", "Mvir", "Rvir")
timestep = tangos.get_simulation("snapshots_200crit_" + simshort).timesteps[-1]
grpz0 = timestep.calculate_all("Grp") #Redshift zero data for all galaxies in the tangos database

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

# For each halo in the list, calculate the main progenitor at all time steps
"""
mpb.clear()
mpb_mvir.clear()
mpb_dist.clear()
"""
mpb = []
mpb_grp = []
mpb_mvir = []
mpb_dist = []
mpb_close_id = []
mpb_close_grp = []
mpb_gforce = []
mpb_gforce_id = []
mpb_gforce_mass = []
mpb_main_dist = []
for index, row in objs_pd_200.iterrows():
    print(row['haloid'])
    haloind = (np.where(grpz0[0] == int(row['haloid'])))[0]
    if len(haloind) == 1:
        haloind = haloind[0] + 1
        tangos_halo = tangos.get_halo("snapshots_200crit_" + simshort + '/' + '%.004096/halo_' + str(haloind))
        if tangos_halo == None:
            ids = np.empty(len(alltimes))*0
            ids[-1] = int(haloind)
            grp = np.empty(len(alltimes))*0
            grp[-1] = int(row['haloid'])        
            mvir = np.empty(len(alltimes))*0
            mvir[-1] = row['mass']
        else:
            time_short, ids_short, grp_short, mvir_short, z_short  = tangos_halo.calculate_for_progenitors("t()", "halo_number()", "Grp", "Mvir", "z()")
            sort_times = numpy.searchsorted(alltimes,time_short)
            ids = np.empty(len(alltimes))*0
            ids[sort_times] = ids_short
            grp = np.empty(len(alltimes))*0
            grp[sort_times] = grp_short        
            mvir = np.empty(len(alltimes))*0
            mvir[sort_times] = mvir_short
    else: # halo not in tangos database
        ids = np.empty(len(alltimes))*0
        grp = np.empty(len(alltimes))*0
        mvir = np.empty(len(alltimes))*0
    mpb_temp = {snap_nums[i]: int(ids[i]) for i in range(len(snap_nums))}
    mpb.append(mpb_temp.copy())
    mpb_grp_temp = {snap_nums[i]: int(grp[i]) for i in range(len(snap_nums))}
    mpb_grp.append(mpb_grp_temp.copy())    
    mpb_mvir_temp = {snap_nums[i]: mvir[i] for i in range(len(snap_nums))}
    mpb_mvir.append(mpb_mvir_temp.copy())
    mpb_dist_temp = {snap_nums[i]: 0 for i in range(len(snap_nums))}
    mpb_dist.append(mpb_dist_temp.copy())
    mpb_close_id_temp = {snap_nums[i]: 0 for i in range(len(snap_nums))}
    mpb_close_id.append(mpb_close_id_temp.copy())
    mpb_close_grp_temp = {snap_nums[i]: 0 for i in range(len(snap_nums))}
    mpb_close_grp.append(mpb_close_grp_temp.copy())    
    mpb_gforce_temp = {snap_nums[i]: 0 for i in range(len(snap_nums))}
    mpb_gforce.append(mpb_gforce_temp.copy())
    mpb_gforce_id_temp = {snap_nums[i]: 0 for i in range(len(snap_nums))}
    mpb_gforce_id.append(mpb_gforce_id_temp.copy())
    mpb_gforce_mass_temp = {snap_nums[i]: 0 for i in range(len(snap_nums))}
    mpb_gforce_mass.append(mpb_gforce_mass_temp.copy())
    mpb_main_dist_temp = {snap_nums[i]: 0 for i in range(len(snap_nums))}
    mpb_main_dist.append(mpb_main_dist_temp.copy())    
    #print(mpb)

mindist = {'min_dist': ""}
mindistt = {'min_dist_t': ""}
objs_pd_200 = objs_pd_200.join(pd.DataFrame(columns=mindist))
objs_pd_200 = objs_pd_200.join(pd.DataFrame(columns=mindistt))

objs_pd_200['mpb_id'] = mpb
objs_pd_200['mpb_grp'] = mpb_grp
objs_pd_200['mpb_mvir'] = mpb_mvir
objs_pd_200['mpb_dist'] = mpb_dist
#mpb_close_id = mpb_dist.copy()
#mpb_gforce = mpb_dist.copy()
#mpb_gforce_id = mpb_dist.copy()
objs_pd_200['mpb_close_id'] = mpb_close_id
objs_pd_200['mpb_gforce'] = mpb_gforce
objs_pd_200['mpb_gforce_id'] = mpb_gforce_id
objs_pd_200['mpb_gforce_mass'] = mpb_gforce_mass
objs_pd_200['mpb_main_dist'] = mpb_main_dist

# Loop through the time steps, calculating the distance to the nearest halo
# timestep = snap_nums[0]
for timestep in snap_nums:
    print(filepath + 'snapshots_200crit_' + simshort + '/' + filename + '.00' + timestep)
    halomasses = np.array([row['mpb_mvir'][timestep] for index, row in objs_pd_200.iterrows()])
    if max(halomasses) == 0: continue
    minmass = np.min(halomasses[np.nonzero(halomasses)])
    
    sim = pynbody.load(filepath + 'snapshots_200crit_' + simshort + '/' + filename + '.00' + timestep)
    h_dummy = sim.halos(dummy = True)
    loc = []
    rvir = []
    mvir = []
    a = []
    ids = []
    for AHFhalo in h_dummy:
        properties = AHFhalo.properties    
        if (properties['mass']/properties['h'] > minmass):
            loc.append(np.array([properties['Xc']/properties['h'], properties['Yc']/properties['h'], properties['Zc']/properties['h']]))
            rvir.append(properties['Rvir']/properties['h'])
            mvir.append(properties['mass']/(properties['h']))
            a.append(properties['a'])
            ids.append(properties['halo_id'])
    loc = np.array(loc)
    rvir = np.array(rvir)            
    mvir = np.array(mvir)
    a = np.array(a)
    ids = np.array(ids)
    
    #row = objs_pd_200.iloc[0]
    #for index, row in objs_pd_200.iterrows():
    mpb_id = objs_pd_200.iloc[0]['mpb_grp'][timestep]
    if mpb_id != 0:
        properties_main =  h_dummy[mpb_id].properties
    else: properties_main = None
    for i in arange(0,len(objs_pd_200)):
        #haloid = row['haloid'] #replace .iloc[i] with [objs_pd_200['haloid'] == haloid]
        mpb_id = objs_pd_200.iloc[i]['mpb_grp'][timestep]
        if mpb_id != 0:
            properties = h_dummy[mpb_id].properties
            #print(i, objs_pd_200.iloc[i]['haloid'], mpb_id, properties['halo_id'])

            if properties_main != None:
                objs_pd_200.iloc[i]['mpb_main_dist'][timestep] = ((properties['Xc']/properties['h'] - properties_main['Xc']/properties_main['h'])**2 + (properties['Yc']/properties['h'] - properties_main['Yc']/properties_main['h'])**2 + (properties['Zc']/properties['h'] - properties_main['Zc']/properties_main['h'])**2)**(0.5)*a[0]
            
            moremass = (mvir > properties['mass']/properties['h'])            
            if np.sum(moremass) > 0:
                dists = ((properties['Xc']/properties['h'] - loc[moremass,0])**2 + (properties['Yc']/properties['h'] - loc[moremass,1])**2 + (properties['Zc']/properties['h'] - loc[moremass,2])**2)**(0.5)
                objs_pd_200.iloc[i]['mpb_dist'][timestep] = np.min(dists[np.nonzero(dists)])
                minind = np.where((((properties['Xc']/properties['h'] - loc[:,0])**2 + (properties['Yc']/properties['h'] - loc[:,1])**2 + (properties['Zc']/properties['h'] - loc[:,2])**2)**(0.5)) == objs_pd_200.iloc[i]['mpb_dist'][timestep])[0][0]
                objs_pd_200.iloc[i]['mpb_close_id'][timestep] = ids[minind]
                objs_pd_200.iloc[i]['mpb_dist'][timestep] = objs_pd_200.iloc[i]['mpb_dist'][timestep]*a[0]
                
                gforce = mvir/((properties['Xc']/properties['h'] - loc[:,0])**2 + (properties['Yc']/properties['h'] - loc[:,1])**2 + (properties['Zc']/properties['h'] - loc[:,2])**2)/a[0]**2 #Grav force
                #gforce = mvir/np.sqrt((properties['Xc']/properties['h'] - loc[:,0])**2 + (properties['Yc']/properties['h'] - loc[:,1])**2 + (properties['Zc']/properties['h'] - loc[:,2])**2)/a #Grav potential
                gforce_max = np.nanmax(gforce[moremass])
                objs_pd_200.iloc[i]['mpb_gforce'][timestep] = gforce_max
                objs_pd_200.iloc[i]['mpb_gforce_id'][timestep] = ids[np.where(gforce == gforce_max)[0][0]]
                objs_pd_200.iloc[i]['mpb_gforce_mass'][timestep] = mvir[np.where(gforce == gforce_max)[0][0]]
                
                #print('Dist: ',ids[minind],objs_pd_200.iloc[i]['mpb_close_id'][timestep],objs_pd_200.iloc[i]['mpb_dist'][timestep],np.min(dists[np.nonzero(dists)]))                
                #print('Grav: ',ids[np.where(gforce == gforce_max)[0][0]],objs_pd_200.iloc[i]['mpb_gforce_id'][timestep],objs_pd_200.iloc[i]['mpb_gforce'][timestep],np.nanmax(gforce[np.isfinite(gforce)]))
                if ids[minind] != ids[np.where(gforce == gforce_max)[0][0]]:
                    print("{0}, {1} ({2:4.2e}): {3:3d} ({4:4.2e}, {5:4.2e}) vs {6:3d} ({7:4.2e}, {8:4.2e}), {9:4.2e}. ".format(mpb_id,objs_pd_200.iloc[i]['haloid'],properties['mass'],ids[minind],mvir[minind],objs_pd_200.iloc[i]['mpb_dist'][timestep],ids[np.where(gforce == gforce_max)[0][0]],mvir[np.where(gforce == gforce_max)[0][0]],sqrt(mvir[np.where(gforce == gforce_max)[0][0]]/objs_pd_200.iloc[i]['mpb_gforce'][timestep]/a[0]**2),objs_pd_200.iloc[i]['mpb_main_dist'][timestep]))
                else:
                    print("{0}, {1} ({2:4.2e}): {3:3d} ({4:4.2e}, {5:4.2e}) == {6:3d} ({7:4.2e}, {8:4.2e}), {9:4.2e}. ".format(mpb_id,objs_pd_200.iloc[i]['haloid'],properties['mass'],ids[minind],mvir[minind],objs_pd_200.iloc[i]['mpb_dist'][timestep],ids[np.where(gforce == gforce_max)[0][0]],mvir[np.where(gforce == gforce_max)[0][0]],sqrt(mvir[np.where(gforce == gforce_max)[0][0]]/objs_pd_200.iloc[i]['mpb_gforce'][timestep]/a[0]**2),objs_pd_200.iloc[i]['mpb_main_dist'][timestep]))
                    #print(timestep,objs_pd_200.iloc[i]['haloid'],ids[minind],ids[np.where(gforce == gforce_max)[0][0]])
            else:
                objs_pd_200.iloc[i]['mpb_dist'][timestep] = 0
                objs_pd_200.iloc[i]['mpb_close_id'][timestep] = 0
                objs_pd_200.iloc[i]['mpb_gforce'][timestep] = 0
                objs_pd_200.iloc[i]['mpb_gforce_id'][timestep] = 0
                objs_pd_200.iloc[i]['mpb_gforce_mass'][timestep] = 0
#                objs_pd_200.iloc[i]['mpb_main_dist'][timestep] = 0

# Calculate the distance when the gravitational force was at the peak and the time that it happened
for i in arange(0,len(objs_pd_200)):  
    if np.max(np.array(list(objs_pd_200.iloc[i]['mpb_gforce'].values()))) != 0:
        #print(np.max(np.array(list(objs_pd_200.iloc[i]['mpb_gforce'].values()))))
        ind = np.argmax(np.array(list(objs_pd_200.iloc[i]['mpb_gforce'].values())))
        objs_pd_200.iloc[i,objs_pd_200.columns.get_loc('min_dist')] = np.sqrt(np.array(list(objs_pd_200.iloc[i]['mpb_gforce_mass'].values()))/np.array(list(objs_pd_200.iloc[i]['mpb_gforce'].values())))[ind]
        objs_pd_200.iloc[i,objs_pd_200.columns.get_loc('min_dist_t')] = alltimes[ind]
        print(objs_pd_200.iloc[i]['min_dist_t'],objs_pd_200.iloc[i]['min_dist'])
    else:
        objs_pd_200.iloc[i]['min_dist_t'] = objs_pd_200.iloc[i]['time']
        objs_pd_200.iloc[i]['min_dist'] = objs_pd_200.iloc[i]['massiveDist']        
    
objs_pd_200.to_pickle(filepath + filename + '.00' + snap_nums[-1] + '/ahf_200/' + filename + '.00' + snap_nums[-1] + '.m200.dist.data')
    
for index, row in objs_pd_200.iterrows(): print(row['mpb_gforce_id']['4096'])

scalefactor = 1 #I have to fix this later
# Distance to the halo that provides a maximum gravitational force
plt.clf()
fig1 = plt.figure(1)
ax1 = fig1.add_subplot(1,1,1)
for index, row in objs_pd_200.iloc[:].iterrows():
    ts_link = np.array(list(row['mpb_gforce_mass'].values())) != 0
    ax1.plot(alltimes[ts_link],(np.sqrt(np.array(list(row['mpb_gforce_mass'].values()))/np.array(list(row['mpb_gforce'].values())))/np.flip(scalefactor))[ts_link],alpha = 0.4)
row = objs_pd_200.iloc[0]
#ax1.plot(alltimes[ts_link],(np.sqrt(np.array(list(row['mpb_gforce_mass'].values()))/np.array(list(row['mpb_gforce'].values())))/np.flip(scalefactor))[ts_link],alpha = 0.4)
ax1.plot(time,rvir_main*scalefactor,color = 'k')    
#ax.axis([0, 14, 2e5, 3e9])
ax1.set_yscale('log')
ax1.set_xlabel('Time [Gyr]')
ax1.set_ylabel('Distance to Most Impactful Companion')

# Distance to main halo
fig2 = plt.figure(2)
ax2 = fig2.add_subplot(1,1,1)
for index, row in objs_pd_200.iloc[:].iterrows():
    ts_link = np.array(list(row['mpb_gforce_mass'].values())) != 0
    ax2.plot(alltimes[ts_link],(np.array(list(row['mpb_main_dist'].values()))/np.flip(scalefactor))[ts_link],alpha = 0.4)
row = objs_pd_200.iloc[0]
#ax2.plot(alltimes[ts_link],(np.array(list(row['mpb_main_dist'].values()))/np.flip(scalefactor))[ts_link],alpha = 0.4)
ax2.plot(time,rvir_main*scalefactor,color = 'k')    
#ax.axis([0, 14, 2e5, 3e9])
ax2.set_yscale('log')
ax2.set_xlabel('Time [Gyr]')
ax2.set_ylabel('Distance to Main Halo Progenitor')

# Distance to main halo
fig3 = plt.figure(3)
ax3 = fig3.add_subplot(1,1,1)
for index, row in objs_pd_200.iloc[:].iterrows():
    ts_link = np.array(list(row['mpb_gforce_mass'].values())) != 0
    ax3.plot(alltimes[ts_link],(np.array(list(row['mpb_main_dist'].values()))/np.flip(scalefactor) -np.sqrt(np.array(list(row['mpb_gforce_mass'].values()))/np.array(list(row['mpb_gforce'].values())))/np.flip(scalefactor))[ts_link] ,alpha = 0.4)
#ax.axis([0, 14, 2e5, 3e9])
ax3.plot(time,rvir_main*scalefactor,color = 'k')   
ax3.set_xlabel('Time [Gyr]')
ax3.set_ylabel('Distance between main progenitor and most impactful halo')

fig4 = plt.figure(4)
ax4 = fig4.add_subplot(1,1,1)
for index, row in objs_pd_200.iloc[:].iterrows():
    ts_link = np.array(list(row['mpb_gforce_mass'].values())) != 0
    ax4.plot(alltimes[ts_link],(np.array(list(row['mpb_gforce'].values())))[ts_link] ,alpha = 0.4)
#ax.axis([0, 14, 2e5, 3e9])
ax4.set_xlabel('Time [Gyr]')
ax4.set_ylabel('Gravitational Force')
ax4.set_yscale('log')
#ax4.axis([0, 14, 2e5, 3e9])
