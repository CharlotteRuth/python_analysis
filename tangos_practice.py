import os
import pandas as pd
import pynbody
import tangos
import matplotlib.colors as colors
from matplotlib import path
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap
import sys, os, glob, pickle

filepath = '/home/christenc/Data/Sims/h148.cosmo50PLK.3072g/h148.cosmo50PLK.3072g3HbwK1BH/'
filename = filepath.split('/')[-2]
simname = 'Sandra'
simshort = 'h148'

filepath = '/home/christenc/Data/Sims/h229.cosmo50PLK.3072g/h229.cosmo50PLK.3072gst5HbwK1BH/'
filename = filepath.split('/')[-2]
simname = 'Ruth'
simshort = 'h229'



filepath = '/home/christenc/Data/Sims/h329.cosmo50PLK.3072g/h329.cosmo50PLK.3072gst5HbwK1BH/'
filename = filepath.split('/')[-2]
simname = 'Elena'
simshort = 'h329'


#tangos.init_db(filepath + simshort + '.db')
tangos.init_db('/home/christenc/Storage/tangos_db/JL_r200.db')

timesteps = tangos.get_simulation("snapshots_200crit_" + simshort).timesteps
snap_nums = [re.findall(r'.00+[\d]+', str(snap))[0][3:] for snap in timesteps]
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
time, a, id_main, mvir_main, rvir_main = tangos_halo.calculate_for_progenitors("t()", "a()", "halo_number()", "Mvir", "Rvir")
                                                                        
# For each halo in the list, calculate the main progenitor at all time steps
"""
mpb.clear()
mpb_mvir.clear()
mpb_dist.clear()
"""
mpb = []
mpb_mvir = []
mpb_dist = []
mpb_close_id = []
mpb_gforce = []
mpb_gforce_id = []
for index, row in objs_pd.iterrows():
    print(row['haloid'])
    tangos_halo = tangos.get_halo("snapshots_200crit_" + simshort + '/' + '.004096/halo_' + str(row['haloid']))
    if tangos_halo == None:
        ids = np.empty(len(snap_nums))*0
        ids[-1] = int(row['haloid'])
        mvir = np.empty(len(snap_nums))*0
        mvir[-1] = row['mvir']
    else:
        time_short, ids_short, mvir_short = tangos_halo.calculate_for_progenitors("t()", "halo_number()", "Mvir")
        ids = np.empty(len(snap_nums))*0
        ids[-1*len(ids_short):] = np.flip(ids_short)
        mvir = np.empty(len(snap_nums))*0
        mvir[-1*len(ids_short):] = np.flip(mvir_short)        
    mpb_temp = {snap_nums[i]: int(ids[i]) for i in range(len(snap_nums))}
    mpb.append(mpb_temp.copy())
    mpb_mvir_temp = {snap_nums[i]: mvir[i] for i in range(len(snap_nums))}
    mpb_mvir.append(mpb_mvir_temp.copy())
    mpb_dist_temp = {snap_nums[i]: 0 for i in range(len(snap_nums))}
    mpb_dist.append(mpb_dist_temp.copy())
    mpb_close_id_temp = {snap_nums[i]: 0 for i in range(len(snap_nums))}
    mpb_close_id.append(mpb_close_id_temp.copy())
    mpb_gforce_temp = {snap_nums[i]: 0 for i in range(len(snap_nums))}
    mpb_gforce.append(mpb_gforce_temp.copy())
    mpb_gforce_id_temp = {snap_nums[i]: 0 for i in range(len(snap_nums))}
    mpb_gforce_id.append(mpb_gforce_id_temp.copy())     
    #print(mpb)

objs_pd['mpb_id'] = mpb
objs_pd['mpb_mvir'] = mpb_mvir
objs_pd['mpb_dist'] = mpb_dist
#mpb_close_id = mpb_dist.copy()
#mpb_gforce = mpb_dist.copy()
#mpb_gforce_id = mpb_dist.copy()
objs_pd['mpb_close_id'] = mpb_close_id
objs_pd['mpb_gforce'] = mpb_gforce
objs_pd['mpb_gforce_id'] = mpb_gforce_id
mindist = {'min_dist': ""}
mindistt = {'min_dist_t': ""}
objs_pd = objs_pd.join(pd.DataFrame(columns=mindist))
objs_pd = objs_pd.join(pd.DataFrame(columns=mindistt))
# Loop through the time steps, calculating the distance to the nearest halo
# timestep = snap_nums[0]
for timestep in snap_nums:
    print(filepath + 'snapshots/' + filename + '.00' + timestep)
    halomasses = np.array([row['mpb_mvir'][timestep] for index, row in objs_pd.iterrows()])
    minmass = np.min(halomasses[np.nonzero(halomasses)])
    
    sim = pynbody.load(filepath + 'snapshots/' + filename + '.00' + timestep)
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
            mvir.append(properties['mass']/(properties['mass']))
            a.append(properties['a'])
            ids.append(properties['halo_id'])
    loc = np.array(loc)
    rvir = np.array(rvir)            
    mvir = np.array(mvir)
    a = np.array(a)
    ids = np.array(ids)
    
    #row = objs_pd.iloc[0]
    #for index, row in objs_pd.iterrows():
    for i in arange(0,len(objs_pd)):
        mpb_id = objs_pd.iloc[i]['mpb_id'][timestep]
        if mpb_id != 0:
            properties = h_dummy[mpb_id].properties
            if np.sum(mvir > properties['mass']/properties['h']) > 0:
                dists = ((properties['Xc']/properties['h'] - loc[mvir > properties['mass']/properties['h'],0])**2 + (properties['Yc']/properties['h'] - loc[mvir > properties['mass']/properties['h'],1])**2 + (properties['Zc']/properties['h'] - loc[mvir > properties['mass']/properties['h'],2])**2)**(0.5)
                objs_pd.iloc[i]['mpb_dist'][timestep] = np.min(dists[np.nonzero(dists)])
                minind = np.where((((properties['Xc']/properties['h'] - loc[:,0])**2 + (properties['Yc']/properties['h'] - loc[:,1])**2 + (properties['Zc']/properties['h'] - loc[:,2])**2)**(0.5)) == objs_pd.iloc[i]['mpb_dist'][timestep])[0][0]
                objs_pd.iloc[i]['mpb_close_id'][timestep] = ids[minind]
                
                gforce = mvir/((properties['Xc']/properties['h'] - loc[:,0])**2 + (properties['Yc']/properties['h'] - loc[:,1])**2 + (properties['Zc']/properties['h'] - loc[:,2])**2)/a**2
                gforce_max = np.nanmax(gforce[np.isfinite(gforce)])
                objs_pd.iloc[i]['mpb_gforce'][timestep] = gforce_max
                objs_pd.iloc[i]['mpb_gforce_id'][timestep] = ids[np.where(gforce == gforce_max)[0][0]]
                
                #print('Dist: ',ids[minind],objs_pd.iloc[i]['mpb_close_id'][timestep],objs_pd.iloc[i]['mpb_dist'][timestep],np.min(dists[np.nonzero(dists)]))                
                #print('Grav: ',ids[np.where(gforce == gforce_max)[0][0]],objs_pd.iloc[i]['mpb_gforce_id'][timestep],objs_pd.iloc[i]['mpb_gforce'][timestep],np.nanmax(gforce[np.isfinite(gforce)]))
                if ids[minind] != ids[np.where(gforce == gforce_max)[0][0]]:
                    print("{0}, {1} ({8:4.2e}): {2:3d} ({3:4.2e}, {4:4.2e}) vs {5:3d} ({6:4.2e}, {7:4.2e}). ".format(timestep,objs_pd.iloc[i]['haloid'],ids[minind],mvir[minind],objs_pd.iloc[i]['mpb_dist'][timestep],ids[np.where(gforce == gforce_max)[0][0]],mvir[np.where(gforce == gforce_max)[0][0]],sqrt(mvir[np.where(gforce == gforce_max)[0][0]]/objs_pd.iloc[i]['mpb_gforce'][timestep]/a[0]**2),properties['mass']))
                    #print(timestep,objs_pd.iloc[i]['haloid'],ids[minind],ids[np.where(gforce == gforce_max)[0][0]])
            else:
                objs_pd.iloc[i]['mpb_dist'][timestep] = 0
                objs_pd.iloc[i]['mpb_close_id'][timestep] = 0
                objs_pd.iloc[i]['mpb_gforce'][timestep] = 0
                objs_pd.iloc[i]['mpb_gforce_id'][timestep] = 0

for index, row in objs_pd.iterrows(): print(row['mpb_gforce_id']['4096'])

dist_z0 = []
for index, row in objs_pd.iterrows(): dist_z0.append(row['mpb_dist']['4096'])
dist_z0 = np.array(dist_z0)
plt.clf()
fig = plt.figure(1)
ax = fig.add_subplot(1,1,1)
ax.plot(objs_pd.massiveDist,dist_z0,linestyle ="", marker = "o",alpha = 0.2)
    
plt.clf()
fig = plt.figure(1)
ax = fig.add_subplot(1,1,1)
for index, row in objs_pd.iloc[:].iterrows():
    ax.plot(np.flip(time),np.array(list(row['mpb_dist'].values())),alpha = 0.2)
ax.plot(time,rvir_main,color = 'k')        
#ax.plot(time,mvir_main/(rvir_main**2*a**2),color = 'k')    
#ax.axis([0, 14, 2e5, 3e9])
ax.set_yscale('log')

plt.clf()
fig = plt.figure(1)
ax = fig.add_subplot(1,1,1)
for index, row in objs_pd.iloc[:20].iterrows():
    ax.plot(np.flip(time),np.array(list(row['mpb_dist'].values())) - np.flip(rvir_main/a),alpha = 0.8)
ax.plot((0,14),(0,0),color = 'k')        
ax.axis([0, 14, -500, 500])
#ax.plot(time,mvir_main/(rvir_main**2*a**2),color = 'k')    
#ax.axis([0, 14, 2e5, 3e9])
#ax.set_yscale('log')
    

plt.clf()
fig = plt.figure(1)
ax = fig.add_subplot(1,1,1)
for index, row in objs_pd.iloc[:].iterrows():
    ax.plot(np.flip(time),np.array(list(row['mpb_gforce'].values())),alpha = 0.2)
ax.plot(time,mvir_main/(rvir_main**2*a**2),color = 'k')    
ax.axis([0, 14, 2e5, 3e9])
ax.set_yscale('log')

plt.clf()
fig = plt.figure(1)
ax = fig.add_subplot(1,1,1)
for index, row in objs_pd.iloc[:].iterrows():
    ax.plot(np.flip(time),np.array(list(row['mpb_gforce_id'].values())),alpha = 0.2)
ax.set_yscale('log')
#ax.plot(time,id_main,color = 'k')
#row = objs_pd.iloc[0]


plt.clf()
fig = plt.figure(1)
ax = fig.add_subplot(1,1,1)
for index, row in objs_pd.iterrows():
    ax.plot(np.flip(time),np.abs(np.array(list(row['mpb_close_id'].values()))-np.array(list(row['mpb_gforce_id'].values()))),linestyle ="", marker = "o",alpha = 0.2)
#ax.set_yscale('log')
#ax.set_xscale('log')

