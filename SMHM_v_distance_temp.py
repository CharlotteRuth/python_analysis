
# Charlotte Christensen

# 8/13/19
# Plot the SMHM relation for the marvel and Justice League runs, coloring points according to distance/tau_90

# SMHM for environment paper

#%run /home/christenc/Code/python/python_analysis/SMHM_v_distance
import matplotlib as mpl
#mpl.use('tkagg') #Also can try mpl.use('Agg') #for using mpl over ssh    
#mpl.use('Agg')
import matplotlib.pylab as plt
import matplotlib.mlab as mlab
import matplotlib.colors as colors
import matplotlib.cm as cm
import pynbody
import math
import numpy as np
import socket
import matplotlib.gridspec as gridspec
import pandas as pd
import sys, os, glob, pickle
from scipy.interpolate import interp1d
sys.path.append(os.path.abspath("/home/christenc/Code/python/python_analysis/"))
from modules.user_tools import task
import Simpy
import tangos
import importlib
import os
# EX: importlib.reload(SMHM_v_distance) use after import SMHM_v_distance

def discrete_cmap(N, base_cmap=None):
    """Create an N-bin discrete colormap from the specified input map"""

    # Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:

    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    return base.from_list(cmap_name, color_list, N)

#A function to return the percentiles of data binned along the xaxis
def bin_plt(xdata, ydata, xmin=None, xmax = None, nbins = 10, perc = np.array([10,25,50,75,90])):
    if xmin==None:
        xmin = xdata[np.isfinite(xdata)].min()
    if xmax==None:
        xmax = xdata[np.isfinite(xdata)].max()
    dwidth = (xmax - xmin)/nbins/2
    xaxis = np.linspace(xmin, xmax, num = nbins, endpoint = False) + dwidth
    data = np.zeros((nbins,len(perc))) - 1
    for i in range(nbins):
        if np.sum((xdata >= xaxis[i] - dwidth) & (xdata < xaxis[i] + dwidth)) > 2:
            data[i,:] = np.percentile(ydata[(xdata >= xaxis[i] - dwidth) & (xdata < xaxis[i] + dwidth)],perc)
    return xaxis, data
    
# Matches the halos in the objs_pd file with the halo information from Ferah (m200)
# see "match_halos.py" for development
def match_halos(objs_pd, fdmdata):
    smass_tol = 0.5 #fractional tolerance for the stellar mass
    vmass_tol = 0.9

    match_id = {'m200_haloid': 0}
    if not 'm200_haloid' in objs_pd.keys():
        objs_pd = objs_pd.join(pd.DataFrame(columns=match_id))
    
    for sim in pd.unique(objs_pd['sim']):
        objs_pd_sim = objs_pd[objs_pd['sim'] == sim].copy()
        for halo in np.sort(objs_pd_sim['haloid']):
            possible_match = (float(objs_pd_sim[objs_pd_sim['haloid'] == halo]['M_star']) < fdmdata[fdmdata['simname']==sim]['Mstar_z0']*(1 + smass_tol)) & \
                (float(objs_pd_sim[objs_pd_sim['haloid'] == halo]['M_star']) > fdmdata[fdmdata['simname']==sim]['Mstar_z0']*(1 - smass_tol)) & \
                 (float(objs_pd_sim[objs_pd_sim['haloid'] == halo]['mass']) < fdmdata[fdmdata['simname']==sim]['Mhalo_z0']*(1 + vmass_tol)) & \
                  (float(objs_pd_sim[objs_pd_sim['haloid'] == halo]['mass']) > fdmdata[fdmdata['simname']==sim]['Mhalo_z0']*(1 - vmass_tol))
            if sum(possible_match) == 0:
                pass
#                #print(sim, halo, 'XXX', float(objs_pd[(objs_pd['sim'] == sim) & (objs_pd['haloid'] == halo)]['M_star']), float(objs_pd[(objs_pd['sim'] == sim) & (objs_pd['haloid'] == halo)]['mass']), 'No Match')
            else:
                #print(fdmdata[fdmdata['simname']==sim][possible_match]['halogrp_z0'])
                arg_best_match = np.argmin(np.abs(fdmdata[fdmdata['simname']==sim][possible_match]['Mstar_z0'] - float(objs_pd_sim[objs_pd_sim['haloid'] == halo]['M_star'])))
                index_best_match = (fdmdata[fdmdata['simname']==sim][possible_match]).index[arg_best_match]
                objs_pd.loc[(objs_pd['sim'] == sim) & (objs_pd['haloid'] == halo),'m200_haloid'] = fdmdata.loc[index_best_match]['halogrp_z0']
                print(sim, halo, fdmdata.loc[index_best_match]['halogrp_z0'], float(objs_pd[(objs_pd['sim'] == sim) & (objs_pd['haloid'] == halo)]['M_star']), fdmdata.loc[index_best_match]['Mstar_z0'], \
                          float(objs_pd[(objs_pd['sim'] == sim) & (objs_pd['haloid'] == halo)]['mass']), fdmdata.loc[index_best_match]['Mhalo_z0'])
                fdmdata.loc[index_best_match,'Mstar_z0'] = 0 #set stellar mass to zero so it isn't used again

    return objs_pd

#Calculates the distance to the nearest massive galaxy
def distance_to_nearest_host(data,tfiles):
    massiveDist = []
    min_massiveHalo = 10**11.5
    sprev = ''
    tfile_it = -1
    for i in range(len(data)):
        s = data['sim'].tolist()[i]

        #print(s,data['haloid'].tolist()[i])       
        if 0: #s=='h148' or s=='h229' or s=='h242' or s=='h329': # or s=='h148_6144' or s=='h329_6144': # if sat simulation, find distance to halo 1
            if s != sprev:
                tfile_it = tfile_it + 1
                sprev = s
            h1dist = data['h1dist'].tolist()[i]*0.6776942783267969
            massiveDist.append(h1dist)       
            h1rvir = data['Rvir'][(data.sim==s) & (data.haloid==1)].tolist()[0]*0.6776942783267969           
        else: # if field simulation, find distance to nearest massive DM halo (currently > 0.5e12.5 Msol)
            if s != sprev:
                sim = pynbody.load(tfiles[tfile_it])
                tfile_it = tfile_it + 1
                sprev = s
                h_dummy = sim.halos(dummy = True)
                loc = []
                rvir = []
                for AHFhalo in h_dummy:
                    properties = AHFhalo.properties            
                    if (properties['mass'] > min_massiveHalo):
#                        print('Halo id: ',properties['halo_id'])
                        loc.append(np.array([properties['Xc']/properties['h'], properties['Yc']/properties['h'], properties['Zc']/properties['h']]))
                        rvir.append(properties['Rvir'])
                
                loc = np.array(loc)
                rvir = np.array(rvir)
                #for index, halo in data.iterrows():

            properties = h_dummy[int(data['haloid'].tolist()[i])].properties
            massiveDist.append(min(((properties['Xc']/properties['h'] - loc[:,0])**2 + (properties['Yc']/properties['h'] - loc[:,1])**2 + (properties['Zc']/properties['h'] - loc[:,2])**2)**(0.5)))
            minind = np.where((((properties['Xc'] - loc[:,0])**2 + (properties['Yc'] - loc[:,1])**2 + (properties['Zc'] - loc[:,2])**2)**(0.5)) == massiveDist[-1])
            #hostrvirs.append(rvir[minind]*0.6776942783267969)

    data['massiveDist'] = massiveDist
    return data

def SMHM_v_distance_data(tfiles,outfile_base,tfile_base,*halo_nums):
    if (socket.gethostname() == "ozma.grinnell.edu"):
        dataprefix = '/home/christensen/Code/Datafiles/' 
    else:
        dataprefix = '/home/christenc/Code/Datafiles/'

    use_m200 = 1
    if use_m200:
        ext = '.m200.dist.'
    else:
        ext = '.MAP.'

    # --------------------- Ferah Data giving M_peak ---------------------------
    #f = open(dataprefix+'mstar_vs_mhalo_4Charlotte.txt', 'r')
    

    show_vmax = 0
    if show_vmax:    
        f = open(dataprefix+'Mstar_vs_Mhalo_withVmaxpeak.txt', 'r')
        fdmdata2 = []
        for line in f:
            line = line.strip()
            columns = line.split()
            if len(columns) == 13 and columns[0] != 'Volume':
                source = {}
                source['simname'] = columns[0]
                source['halogrp_z0'] = int(columns[1])
                source['halogrp_Mpeak'] = int(columns[2])
                source['Mpeak_snap'] = float(columns[3])
                source['Mpeak'] = float(columns[4])
                source['Mhalo_z0'] = float(columns[5])
                source['Mstar_z0'] = float(columns[6])
                source['Mstar_z0_photo'] = float(columns[7])
                source['Mstar_Mpeak'] = float(columns[8])
                source['Mstar_Mpeak_z0'] = float(columns[9])
                source['Vmag'] = float(columns[10])
                source['type'] = columns[11]
                source['Vmax'] = float(columns[12])
                fdmdata2.append(source)
        f.close()
        fdmdata = pd.DataFrame(fdmdata2)
    else:
        f = open(dataprefix+'mstar_vs_mhalo_extended_corrected5.txt', 'r')
        fdmdata = []
        for line in f:
            line = line.strip()
            columns = line.split()
            if len(columns) == 12 and columns[0] != 'Volume':
                source = {}
                source['simname'] = columns[0]
                source['halogrp_z0'] = int(columns[1])
                source['halogrp_Mpeak'] = int(columns[2])
                source['Mpeak_snap'] = float(columns[3])
                source['Mpeak'] = float(columns[4])
                source['Mhalo_z0'] = float(columns[5])
                source['Mstar_z0'] = float(columns[6])
                source['Mstar_z0_photo'] = float(columns[7])
                source['Mstar_Mpeak'] = float(columns[8])
                source['Mstar_Mpeak_z0'] = float(columns[9])
                source['Vmag'] = float(columns[10])
                source['type'] = columns[11]            
                fdmdata.append(source)
        f.close()
        fdmdata = pd.DataFrame(fdmdata)

    # ----------------------
    #read dataprefix+'/assembly_histories.npy'
    assem_history = task(
        np.load,
        start_text="Loading assembly history data",
        end_text="Loaded assembly history data",
        fail_text="Failed to load assembly history data",
        exit_on_fail=True
    )(dataprefix + "/assembly_histories.npy", encoding="latin1", allow_pickle=True).item()
    assem_history = pd.DataFrame.from_dict(assem_history)
    
    #times, volumes, halo_grps, stellar_mass_histories, Rvir_histories, masses, tracing_fractions = ( assem_history[key] for key in ["Time","Volume","halo grp @ z=0","Stellar Mass History","Rvir History","M_halo @ z=0","Tracing Fraction"] )
     
    '''
    assem_history = []
    f=open(dataprefix+'/assembly_histories.npy', 'rb')
    while 1:
        try:
            assem_history.append(pickle.load(f))
        except EOFError:
            break        
    f.close() 
    '''  

    #read dataprefix+'/reduced_time_series_data.npy'
    """
    time_series = []
    f=open(dataprefix+'/reduced_time_series_data.npy', 'rb')
    while 1:
        try:
            time_series.append(pickle.load(f))
        except EOFError:
            break        
    f.close() 
    """
    
    # Does not include data for Mint runs
    data = task(
        np.load,
        start_text="Loading data",
        end_text="Loaded data",
        fail_text="Failed to load data",
        exit_on_fail=True
    )(dataprefix + "reduced_time_series_data.npy", encoding="latin1", allow_pickle=True).item()   

    # Includes concentrations from einasto profiles; Does not have data for Mint runs
    # Ex: concentrations['rogue'][1]['concentration']
    concentrations = np.load(dataprefix + 'concentrations_and_einasto.npy', allow_pickle = True).item() 

    # 
    objs_pd = None 
    for tfile, base in zip(tfiles, tfile_base):
        objs_dat = []
        
        f=open(tfile + ext + 'data', 'rb') # For distance data, use calc_closest_dist.py
        print(f)
        while 1:
            try:
                objs_dat.append(pickle.load(f))
            except EOFError:
                break        
        f.close()
        
        #objs_dat = pd.read_csv(tfile + '.MAP.data.csv')
        if len(objs_dat) == 1:
            temp = pd.DataFrame(objs_dat[0])
        else:
            temp = pd.DataFrame(objs_dat)
        simname = base.split('.')[0]
        if (base.split('.')[2])[0] == '6':
            simname = simname+'_6144'
        if not ('M_star' in temp.keys()):
            temp['M_star'] = temp['mstar']
        if not ('mass' in temp.keys()):            
            temp['mass'] = temp['mvir']
            
        temp['sim'] = [simname]*len(temp)
        if not 'massiveDist' in temp:
            temp = distance_to_nearest_host(temp,[tfile])
            temp.to_pickle(tfile + ext + 'data')
            temp.to_csv(tfile + ext + 'data.csv', index=False)

        if objs_pd is None: 
            objs_pd = temp
        else:
            objs_pd = objs_pd.append(temp, ignore_index = True)

    #Match halos between my and Ferah's data
    fdmdata_mod = fdmdata.copy()
    objs_pd = match_halos(objs_pd, fdmdata_mod)
    #remove entries (rows) from objs_pd that have no match in Ferah's data
    index_rm = objs_pd[(objs_pd['m200_haloid']).isnull()].index
    objs_pd = objs_pd.drop(index_rm)
    #objs_pd = distance_to_nearest_host(objs_pd, tfiles)
    
    """
        h_dummy = s.halos(dummy = True)
        loc = []        
        for AHFhalo in h_dummy:
            properties = AHFhalo.properties            
            if (properties['mass'] > min_massiveHalo):
                loc.append(np.array([properties['Xc']/properties['h'], properties['Yc']/properties['h'], properties['Zc']/properties['h']]))
        loc = np.array(loc)
        massiveDist = [] #Distance to nearest massive (M_vir > min_massiveHalo) galaxy
        for halo in objs_dat:
            massiveDist.append(min(((halo['Xc'] - loc[:,0])**2 + (halo['Yc'] - loc[:,1])**2 + (halo['Zc'] - loc[:,2])**2)**(0.5)))
        temp = pd.DataFrame(objs_dat)
        temp['massiveDist'] = massiveDist
        if objs_pd is None: 
            objs_pd = temp
        else:
            objs_pd = objs_pd.append(temp, ignore_index = True)
    """

    tau90_vir = np.zeros(len(objs_pd))
    tau50_vir = np.zeros(len(objs_pd))
    t_max = np.zeros(len(objs_pd)) # Time when mvir was at a maximum
    
    smass_z0_5 = np.zeros(len(objs_pd))
    vmass_z0_5 = np.zeros(len(objs_pd))
    smass_z1 = np.zeros(len(objs_pd))
    vmass_z1 = np.zeros(len(objs_pd))
    smass_z2 = np.zeros(len(objs_pd))
    vmass_z2 = np.zeros(len(objs_pd))    
    smass_z4 = np.zeros(len(objs_pd))
    vmass_z4 = np.zeros(len(objs_pd))
    
    time_z0_5 = Simpy.cosmology.getTime(0.5) # Find stellar mass to halo mass at these different
    time_z1 = Simpy.cosmology.getTime(1)    # redshifts, only plotting it if the galaxy has not yet
    time_z2 = Simpy.cosmology.getTime(2)    # reached peak halo mass
    time_z4 = Simpy.cosmology.getTime(4)
    
    ind = 0
    for index, row in objs_pd.iterrows():
        cond = ((assem_history['Volume'] == row['sim']) & (assem_history['halo grp @ z=0'] == int(row['haloid'])))
        if np.sum(cond) == 0:
            print(row['sim'],row['haloid'])
            tau90_vir[ind] = -1
            tau50_vir[ind] = -1
            t_max[ind] = -1
            ind = ind + 1
            continue
        xarr = np.array(((assem_history['Time'][cond]).tolist())[0])
        assem_history_halo = np.array(((assem_history['Assembly History'][cond]).tolist())[0])
        star_history_halo = np.array(((assem_history['Stellar Mass History'][cond]).tolist())[0])
        yarr = (assem_history_halo)/np.max(assem_history_halo)

        #fig1 = plt.figure()
        #fig1.clear()
        #ax1 = fig1.add_subplot()
        #ax1.plot(xarr, yarr)
        #ax1.set_title(row['sim']+", "+row['haloid'])
        #ax1.axvline(x = (fdmdata.loc[row['sim']+row['haloid']]['Mpeak_snap']*13.8/4096), color = 'k')
        #fig1.show()
        #wait = input('continue')
        #plt.clf()

        end = np.argmax(yarr)
        t_max[ind] = xarr[end]
        if len(xarr) != 0:
            xarr = xarr[0:end + 1]
            yarr = yarr[0:end + 1]
            star_history_halo = star_history_halo[0:end + 1]
        if (yarr[0] >= 0.9):
            tau90_vir[ind] = xarr[0]
        else:
            interp = interp1d(yarr, xarr) #, kind='cubic')
            if np.isnan(interp(0.9)):
                tau90_vir[ind] = 0
                print('Problem: ',row['sim'],row['haloid'])
            else:
                tau90_vir[ind] = float(interp(0.9))
        if (yarr[0] >= 0.5):
            tau50_vir[ind] = xarr[0]
        else:
            interp = interp1d(yarr, xarr) #, kind='cubic')
            if np.isnan(interp(0.5)):
                tau50_vir[ind] = 0
                print('Problem: ',row['sim'],row['haloid'])
            else:
                tau50_vir[ind] = float(interp(0.5))
        if tau90_vir[ind] >= 14: #13.8:
            print(row['sim'],row['haloid'],tau90_vir[ind])
        if tau50_vir[ind] >= 14: #13.8:
            print(row['sim'],row['haloid'],tau50_vir[ind])

        if len(xarr) > 1:
            interp = interp1d(xarr, yarr*np.max(assem_history_halo))
            interp_star = interp1d(xarr, star_history_halo)
        else:
            continue
        if (np.min(xarr) < time_z4) and (np.max(xarr) > time_z4):
           vmass_z4[ind] = interp(time_z4)
           smass_z4[ind] = interp_star(time_z4)
        if (np.min(xarr) < time_z2) and (np.max(xarr) > time_z2):
           vmass_z2[ind] = interp(time_z2)
           smass_z2[ind] = interp_star(time_z2)
        if (np.min(xarr) < time_z1) and (np.max(xarr) > time_z1):
           vmass_z1[ind] = interp(time_z1)
           smass_z1[ind] = interp_star(time_z1)           
        if (np.min(xarr) < time_z0_5) and (np.max(xarr) > time_z0_5):
           vmass_z0_5[ind] = interp(time_z0_5)
           smass_z0_5[ind] = interp_star(time_z0_5)
        ind = ind + 1
    #        #    if (np.min(xarr) < time_z0) and (np.max(xarr) > time_z0):
    #        #       vmass_z0[ind] = interp(time_z0)
        
    objs_pd['tau90_vir'] = tau90_vir
    objs_pd['tau50_vir'] = tau50_vir
    objs_pd['vmass_z4'] = vmass_z4
    objs_pd['vmass_z2'] = vmass_z2
    objs_pd['vmass_z1'] = vmass_z1
    objs_pd['vmass_z0_5'] = vmass_z0_5
    #objs_pd['vmass_z0'] = vmass0
    objs_pd['smass_z4'] = smass_z4
    objs_pd['smass_z2'] = smass_z2
    objs_pd['smass_z1'] = smass_z1
    objs_pd['smass_z0_5'] = smass_z0_5       

    tau90 = np.zeros(len(objs_pd))
    tau50 = np.zeros(len(objs_pd))
    concent = np.zeros(len(objs_pd))
    ind = 0    
    for index, row in objs_pd.iterrows():
#        print(row['sim'],row['haloid'])
#        row['sfh'] = row['sfh'].replace('  ',' ')
#         row['sfh'] = row['sfh'].replace('   ',' ')
#         row['sfh'] = row['sfh'].replace('    ',' ')
#         row['sfh'] = row['sfh'].replace('     ',' ')
#        row['sfh'] = row['sfh'].replace('      ',' ')
#         row['sfh'] = row['sfh'].replace('       ',' ')
#         row['sfh'] = row['sfh'].replace('        ',' ')
#         row['sfh'] = row['sfh'].replace('         ',' ')
#         row['sfh'] = row['sfh'].replace('          ',' ')
#         row['sfh'] = row['sfh'].replace('           ',' ')
#         row['sfh'] = row['sfh'].replace('            ',' ')        
#         row['sfh'] = row['sfh'].replace('             ',' ')
#         sfh_str = ((row['sfh'])[2:-1].replace('\n','')).split(' ')
#         sfh = np.array([float(s) for s in sfh_str if s != ''])
#         #sfh = np.array([float(x) for x in sfh_str])
#         row['sfhbins'] = row['sfhbins'].replace('  ',' ')
#         row['sfhbins'] = row['sfhbins'].replace('   ',' ')
#         row['sfhbins'] = row['sfhbins'].replace('    ',' ')
#         row['sfhbins'] = row['sfhbins'].replace('     ',' ')
#         sfhbins_str = ((row['sfhbins'])[2:-1].replace('\n','')).split(' ')
#         sfhbins = np.array([float(s) for s in sfhbins_str if s != ''])
#         #sfhbins = np.array([float(x) for x in sfhbins_str])
        sfh = row['sfh']
        sfhbins = row['sfhbins']
        
        if len(sfhbins) != len(sfh):
            xarr = sfhbins[1:] - (sfhbins[1] - sfhbins[0])
        else:
            xarr = sfhbins[:]
        #print(min(xarr),max(xarr),len(xarr))
        yarr = np.cumsum(sfh)/max(np.cumsum(sfh))
        if (yarr[0] >= 0.9):
            tau90[ind] = xarr[0]
        else:
            interp = interp1d(yarr, xarr) #, kind='cubic')
            if np.isnan(interp(0.9)):
                tau90[ind] = 0
            else:
                tau90[ind] = float(interp(0.9))
        if (yarr[0] >= 0.5):
            tau50[ind] = xarr[0]
        else:
            interp = interp1d(yarr, xarr) #, kind='cubic')
            if np.isnan(interp(0.5)):
                tau50[ind] = 0
            else:
                tau50[ind] = float(interp(0.5))

        #if False:
        if (row['sim'] == 'cptmarvel') or (row['sim'] == 'elektra') or (row['sim'] == 'storm') or (row['sim'] == 'rogue'):
            try: # concentrations[row['sim']][int(row['haloid'])]:
                if 'concentration' in concentrations[row['sim']][int(row['haloid'])].keys():
                    concent[ind] = concentrations[row['sim']][int(row['haloid'])]['concentration']
                    print(row['sim'], int(row['haloid']), concentrations[row['sim']][int(row['haloid'])]['concentration']) #concent[ind])
                else:
                        concent[ind] = concentrations[row['sim']][int(row['haloid'])]['scale radius']['r_s']/row['rvir'] # This rvir should be r_200
                        print(row['sim'], int(row['haloid']), concentrations[row['sim']][int(row['haloid'])]['scale radius']['r_s']/row['rvir'], concent[ind])
            except:
                concent[ind] = -1
        else:
            concent[ind] = -1
            if t_max[ind] == -1:
                ind = ind + 1
                continue
            step_max = '{:0>6}'.format(round(t_max[ind]/Simpy.cosmology.getTime(0)*4096))
            
            if row['sim'] == 'h329':
                tfile = '/home/christenc/Data/Sims/h329.cosmo50PLK.3072g/h329.cosmo50PLK.3072gst5HbwK1BH/snapshots_200crit_h329/h329.cosmo50PLK.3072gst5HbwK1BH.' + step_max
            elif row['sim'] == 'h242':
                tfile = '/home/christenc/Data/Sims/h242.cosmo50PLK.3072g/h242.cosmo50PLK.3072gst5HbwK1BH/snapshots_200crit_h242/h242.cosmo50PLK.3072gst5HbwK1BH.' + step_max
            elif row['sim'] == 'h229':
                tfile = '/home/christenc/Data/Sims/h229.cosmo50PLK.3072g/h229.cosmo50PLK.3072gst5HbwK1BH/snapshots_200crit_h229/h229.cosmo50PLK.3072gst5HbwK1BH.' + step_max
            elif row['sim'] == 'h148':
                tfile = '/home/christenc/Data/Sims/h148.cosmo50PLK.3072g/h148.cosmo50PLK.3072g3HbwK1BH/snapshots_200crit_h148/h148.cosmo50PLK.3072g3HbwK1BH.' + step_max
                """
            elif row['sim'] == 'cptmarvel':
                tfile = '/home/christenc/Data/Sims/cptmarvel.cosmo25cmb/cptmarvel.cosmo25cmb.4096g5HbwK1BH/snapshots_200crit_cptmarvel/cptmarvel.cosmo25cmb.4096g5HbwK1BH.' + step_max
            elif row['sim'] == 'elektra':
                tfile = '/home/christenc/Data/Sims/elektra.cosmo25cmb/elektra.cosmo25cmb.4096g5HbwK1BH/snapshots_200crit_elektra/elektra.cosmo25cmb.4096g5HbwK1BH.' + step_max
            elif row['sim'] == 'storm':
                tfile = '/home/christenc/Data/Sims/storm.cosmo25cmb/storm.cosmo25cmb.4096g5HbwK1BH/snapshots_200crit_storm/storm.cosmo25cmb.4096g5HbwK1BH.' + step_max
            elif row['sim'] == 'rogue':
                tfile = '/home/christenc/Data/Sims/rogue.cosmo25cmb/rogue.cosmo25cmb.4096g5HbwK1BH/snapshots_200crit_rogue/rogue.cosmo25cmb.4096g5HbwK1BH.' + step_max
              """  
            if os.path.exists(tfile) :
                s = pynbody.load(tfile)
                h_props = s.halos(dummy = True)
                concent[ind] = h_props[int(row['haloid'])].properties['cNFW']
            else:
                concent[ind] = -1
            print(row['sim'], int(row['haloid']), concent[ind], 'NFW')
        ind = ind + 1

    objs_pd['tau90'] = tau90
    objs_pd['tau50'] = tau50
    concent_pretangos = concent.copy()
    vmax = np.empty(len(concent))
    vmax_max = np.empty(len(concent))
    vmax_mmax = np.empty(len(concent))
    bmax = np.empty(len(concent))
    bz0 = np.empty(len(concent))
    b_peak = np.empty(len(concent))
    
    ind = 0
    for index, row in objs_pd.iterrows():
        if (row['sim'] == 'cptmarvel') or (row['sim'] == 'elektra') or (row['sim'] == 'rogue') or (row['sim'] == 'storm'):
            #tangos_db = '/home/christenc/Storage/tangos_db/Marvel_r200.db'
            tangos_db = '/home/christenc/Storage/tangos_db/Marvel_r200_N100min.db'
            sim_key = row['sim']
        elif (row['sim'] == 'h148') or (row['sim'] == 'h229') or (row['sim'] == 'h242') or (row['sim'] == 'h329'):
            tangos_db = '/home/christenc/Storage/tangos_db/JL_r200_N100min.db'
            sim_key = row['sim']
        elif (row['sim'] == 'h148_6144') or (row['sim'] == 'h329_6144'):
            tangos_db = '/home/christenc/Storage/tangos_db/JLmint_r200_N100.db'
            sim_key = row['sim'][0:4]+'mint'
        else:
            print("Database not available")
        tangos.init_db(tangos_db)
        
        halo = tangos.get_halo("snapshots_200crit_" + sim_key + "/%4096/halo_" + row['haloid'])
        if 'cNFW' in halo.keys():
            concent_prog, mvir_prog, t_prog = halo.calculate_for_progenitors("cNFW", "Mvir", "t()")
            concent[ind] = concent_prog[np.argmax(mvir_prog)]
            print('cNFW: ', concent_prog[np.argmax(mvir_prog)],concent_pretangos[ind])
        if 'Vmax' in halo.keys():
            vmax_prog, t_prog = halo.calculate_for_progenitors("Vmax", "t()")
            #vmax[ind] = vmax_prog[np.argmax(mvir_prog)]
            vmax[ind] = vmax_prog[0]
            vmax_max[ind] = np.max(vmax_prog)
            vmax_mmax[ind] = vmax_prog[np.argmax(mvir_prog)]
            print('Vmax: ', vmax_prog[np.argmax(mvir_prog)], np.max(vmax_prog), vmax_prog[0])
        smass_prog, gmass_prog, mvir_prog = halo.calculate_for_progenitors("M_star", "M_gas", "Mvir")
        bmax[ind] = np.max(smass_prog + gmass_prog)
        b_peak[ind] = smass_prog[np.argmax(mvir_prog)] + gmass_prog[np.argmax(mvir_prog)]
        bz0[ind] = smass_prog[0] + gmass_prog[0]
        ind = ind + 1

    objs_pd['concentration'] = concent
    objs_pd['vmax'] = vmax
    objs_pd['vmax_max'] = vmax_max
    objs_pd['vmax_mmax'] = vmax_mmax
    objs_pd['bmass_max'] = bmax
    objs_pd['bmass_z0'] = bz0
    objs_pd['bmass_peak'] = b_peak
        
    lum = {'Lum': ""}
    objs_pd = objs_pd.join(pd.DataFrame(columns=lum))
    objs_pd['Lum'] = 10**((objs_pd['M_V'] - 4.81)/-2.5)

    #b_v = {'B-V': ""}
    #objs_pd_comb = objs_pd_comb.join(pd.DataFrame(columns=b_v))
    objs_pd['B-V'] = objs_pd['M_B'] - objs_pd['M_V']
    
    halo_label = {'halo_label': ""}
    objs_pd = objs_pd.join(pd.DataFrame(columns=halo_label))
    objs_pd.loc[~objs_pd['m200_haloid'].isnull(),'halo_label'] = objs_pd[~objs_pd['m200_haloid'].isnull()]['sim']+objs_pd[~objs_pd['m200_haloid'].isnull()]['m200_haloid'].astype(str)
    objs_pd = objs_pd.set_index('halo_label')

    fdmdata = fdmdata.join(pd.DataFrame(columns=halo_label))
    fdmdata['halo_label'] = fdmdata['simname']+fdmdata['halogrp_z0'].astype(str)
    fdmdata = fdmdata.set_index('halo_label') 

    objs_pd_comb = pd.concat([objs_pd,fdmdata], join="inner", axis=1)

    #Removing Mint for now
    #objs_pd_comb = objs_pd_comb[(objs_pd_comb['sim'] != 'h148_6144') & (objs_pd_comb['sim'] != 'h329_6144')]

    #Calculate the metals in the gas
    objs_pd_comb['mZgas'] = objs_pd_comb['mZISM'] + objs_pd_comb['mZCool'] + objs_pd_comb['mZwarm'] + objs_pd_comb['mZHot']

    objs_pd_comb_bk = objs_pd_comb
    #objs_pd_comb = objs_pd_comb[(objs_pd_comb['sim'] == 'h329_6144') | (objs_pd_comb['sim'] == 'h148_6144') | (objs_pd_comb['sim'] == 'cptmarvel') | (objs_pd_comb['sim'] == 'elektra') | (objs_pd_comb['sim'] == 'rogue') | (objs_pd_comb['sim'] == 'storm')]
    """
    objs_pd_comb = objs_pd_comb_bk
    objs_pd_comb = objs_pd_comb[(objs_pd_comb['sim'] == 'h329_6144') | (objs_pd_comb['sim'] == 'h148_6144') | (objs_pd_comb['sim'] == 'h242') | (objs_pd_comb['sim'] == 'h229') | \
                                    (objs_pd_comb['sim'] == 'cptmarvel') | (objs_pd_comb['sim'] == 'elektra') | (objs_pd_comb['sim'] == 'rogue') | (objs_pd_comb['sim'] == 'storm')]

    objs_pd_comb = objs_pd_comb_bk
    objs_pd_comb = objs_pd_comb[(objs_pd_comb['sim'] == 'h329_6144') | (objs_pd_comb['sim'] == 'h148_6144') | (objs_pd_comb['sim'] == 'h329') | (objs_pd_comb['sim'] == 'h148')]
    """
    objs_pd_comb_wr = objs_pd_comb[['sim','type','haloid','mass','Mpeak','mHI','Mpeak_snap','Mstar_Mpeak','Mstar_z0','max_tindex','min_dist','massiveDist','concentration','tau90_vir','tau90','vmax']].copy()
    objs_pd_comb_wr.to_csv(dataprefix+'SMHM_vmax_env.csv')

    return objs_pd_comb


def SMHM_v_distance_plts(tfiles,outfile_base,tfile_base,objs_pd_comb,*halo_nums):
    presentation = False
    if presentation:
        outfile_base = outfile_base + '_pres'
        plt.style.use(['default','/home/christenc/.config/matplotlib/presentation.mplstyle'])
        plt_width = 8 #inches
        aspect_ratio = 3.0/4.0
        legendsize = 16
        dpi = 200
        markersize = 100
        ms_scale = 1
        lw = mpl.rcParams['lines.linewidth'] - 1
        edgewidth = 1
    else:
        plt.style.use(['default','/home/christenc/.config/matplotlib/article.mplstyle'])
        plt_width = 3.5 #inches
        aspect_ratio = 3.0/4.0
        legendsize = 5
        dpi = 300
        markersize = 25
        ms_scale = 0.25
        lw = mpl.rcParams['lines.linewidth']
        edgewidth = 0.5
    if (socket.gethostname() == "ozma.grinnell.edu"):
        dataprefix = '/home/christensen/Code/Datafiles/' 
    else:
        dataprefix = '/home/christenc/Code/Datafiles/'

    #Set pointsizes so that lower res simulations are smaller
    objs_pd_comb['p_size'] = np.ones(len(objs_pd_comb))*markersize
    mask = objs_pd_comb[(objs_pd_comb['simname']=='h148') | (objs_pd_comb['simname']=='h229') | (objs_pd_comb['simname']=='h242') | (objs_pd_comb['simname']=='h329')].index
    objs_pd_comb.loc[mask,'p_size'] = markersize*0.5        

# -------------------- READ OBS DATA ------------------------------
     # Data from Digby+ 2019
    digbyfield = pd.read_csv(dataprefix + 'DigbyField.csv')
    digbysat = pd.read_csv(dataprefix + 'DigbySatellite.csv')       

    # From McConnachie data
    MW_M31_dist = [58,184,75,40,110,218,110,104,133,180,162,174,45,1350,1066,64,40,218,161,107,681,2030,1862,149,252,520,2266,2436,1945,2187,803,155,422,179,50,23,142,187,2078,1301,2860,1930,452,\
                       474,415,86,28,1435,1430,61,882,1367,2583,2288,2387,43,836]
    MW_M31_MHI_LV = [0,0,0,0,0,0,0,0,0,0,0,0,0,0.5588855,3.6373396,0,0,0,0,0,0,0.1292644,0.2144484,0.0085202,0.595621,1.3706083,0.3962263,0.9891736,0.624481,4.7199784,0.205775,0,2.1045441,0,\
                         0.3095693,0,0,0.0016497,0.8856712,3.2901259,0.2195973,0.583369,1.2762723,0.8449909,0.1581908,0.0978189,0,2.6453962,1.7043947,1.0156822,0,0.1163385,1.8810399,2.3728621,55.782561,0,2.2978932,]
    Calc_V_Mag = [-11.66,-12.37,-9.97,-8.12,-9.14,-12.61,-7.63,-6.9,-6.4,-6.7,-8.43,-9.4,-8.7,-10.48,-10.32,-3.71,-5.5,-8.59,-4.92,-9.11,-11.19,-11.24,-11.52,-13.44,-15,-14.38,-13.98,-15.84,-15.55,-9.5,-14.51,\
                      -5.84,-8,-5.25,-18.12,-16.45,-14.65,-14.75,-18.46,-15.53,-14.08,-18.56,-15.21,-12.3,-9.89,-11.07,-1.5,-13.85,-13.88,-16.83,-9.54,-12.47,-13.16,-12.39,-13.16,-2.7,-13.75]

    MW_M31_dist = np.array(MW_M31_dist)
    MW_M31_MHI_LV = np.array(MW_M31_MHI_LV)
    Calc_V_Mag = np.array(Calc_V_Mag)

    mcconnachie_data = pd.read_csv(dataprefix + 'McConnachie2012.csv', dtype ={'Min(D_MW, D_M31)': np.float64, 'M_V': np.float64, 'M_H I': np.float64}) #, skiprows = 1 )
    
    f = open(dataprefix+'Read2017.txt','r')
    readdata = []
    for line in f:
        line = line.strip()
        columns = line.split()
        if len(columns) == 22 and columns[0] != 'Galaxy':
            source = {}
            source['Galaxy'] = columns[0]
            source['vmax'] = float(columns[1])
            source['vmax_err'] = float(columns[2])
            source['i'] = float(columns[3])
            source['i_err_n'] = float(columns[4])
            source['i_err_p'] = float(columns[5])
            source['D'] = float(columns[6])
            source['D_err'] = float(columns[7])
            source['M_*'] = float(columns[8])
            source['M_*_err'] = float(columns[9])
            source['Mgas'] = float(columns[10])
            source['R*'] = float(columns[11])
            source['R_gas'] = float(columns[12])
            source['Rmin'] = float(columns[13])
            source['Rmax'] = float(columns[14])
            source['M200'] = float(columns[15])
            source['M200_err_n'] = float(columns[16])
            source['M200_err_p'] = float(columns[17])
            source['c'] = float(columns[18])
            source['c_err_n'] = float(columns[19])
            source['c_err_p'] = float(columns[20])
            source['Xi'] = float(columns[21])
            readdata.append(source)
    f.close()
    readdata = pd.DataFrame(readdata)

    # Read data from Justin Read's abundance matching, provided by Ferah
    f = open(dataprefix+'mstar-mhalo-field.txt', 'r')
    read_abunmatch = []
    for line in f:
        line = line.strip()
        columns = line.split()
        if len(columns) == 3:
            source = {}
            source['M200'] = float(columns[0])
            source['Mstar_low'] = float(columns[1])
            source['Mstar_high'] = float(columns[2])
            read_abunmatch.append(source)
    f.close()
    read_abunmatch = pd.DataFrame(read_abunmatch)

    # Shaded regions in fig 6 of Nadler et al. 2020                                                                                                                  
    upper_fid_Mhalo, upper_fid_Mstar = np.transpose(np.genfromtxt(dataprefix+"webplotdigitizer/fiducial upper bound.csv", dtype="float", usecols=(0,1), skip_header=0))
    lower_fid_Mhalo, lower_fid_Mstar = np.transpose(np.genfromtxt(dataprefix+"webplotdigitizer/fiducial lower bound.csv", dtype="float", usecols=(0,1), skip_header=0))

    upper_fid_inner_Mhalo, upper_fid_inner_Mstar = np.transpose(np.genfromtxt(dataprefix+"webplotdigitizer/fid inner upper bound.csv", dtype="float", usecols=(0,1), skip_header=0))
    lower_fid_inner_Mhalo, lower_fid_inner_Mstar = np.transpose(np.genfromtxt(dataprefix+"webplotdigitizer/fid inner lower bound.csv", dtype="float", usecols=(0,1), skip_header=0))

    # Jethwa et al. 2018 (HO + scatter)                                                                                                                           
#    Jethwa_upper_inner_Mhalo, Jethwa_upper_inner_Mstar = np.transpose(np.genfromtxt(dataprefix+"webplotdigitizer/Jethwa inner upper.csv", dtype="float", usecols=(0,1), skip_header=0))
#    Jethwa_lower_inner_Mhalo, Jethwa_lower_inner_Mstar = np.transpose(np.genfromtxt(dataprefix+"webplotdigitizer/Jethwa inner lower.csv", dtype="float", usecols=(0,1), skip_header=0))

#    Jethwa_upper_outer_Mhalo, Jethwa_upper_outer_Mstar = np.transpose(np.genfromtxt(dataprefix+"webplotdigitizer/Jethwa outer upper.csv", dtype="float", usecols=(0,1), skip_header=0))
#    Jethwa_lower_outer_Mhalo, Jethwa_lower_outer_Mstar = np.transpose(np.genfromtxt(dataprefix+"webplotdigitizer/Jethwa outer lower.csv", dtype="float", usecols=(0,1), skip_header=0))

    jethwa_inner = np.loadtxt(dataprefix+"jethwa_inner.csv",skiprows=1,delimiter=',')
    Jethwa_upper_inner_Mhalo = jethwa_inner[:,0]
    Jethwa_lower_inner_Mhalo = jethwa_inner[:,0]    
    Jethwa_upper_inner_Mstar = 10**jethwa_inner[:,1]
    Jethwa_lower_inner_Mstar = 10**jethwa_inner[:,2]

    jethwa_outer = np.loadtxt(dataprefix+"jethwa_outer.csv",skiprows=1,delimiter=',')
    Jethwa_upper_outer_Mhalo = jethwa_outer[:,0]
    Jethwa_lower_outer_Mhalo = jethwa_outer[:,0]
    Jethwa_upper_outer_Mstar = 10**jethwa_outer[:,1]
    Jethwa_lower_outer_Mstar = 10**jethwa_outer[:,2]    

    Read_upper_dashed_Mhalo, Read_upper_dashed_Mstar = np.transpose(np.genfromtxt(dataprefix+"webplotdigitizer/Read upper dashed.csv", dtype="float", usecols=(0,1), skip_header=0))
    Read_lower_dashed_Mhalo, Read_lower_dashed_Mstar = np.transpose(np.genfromtxt(dataprefix+"webplotdigitizer/Read lower dashed.csv", dtype="float", usecols=(0,1), skip_header=0))
    Read_upper_solid_Mhalo, Read_upper_solid_Mstar = np.transpose(np.genfromtxt(dataprefix+"webplotdigitizer/Read upper solid.csv", dtype="float", usecols=(0,1), skip_header=0))
    Read_lower_solid_Mhalo, Read_lower_solid_Mstar = np.transpose(np.genfromtxt(dataprefix+"webplotdigitizer/Read lower solid.csv", dtype="float", usecols=(0,1), skip_header=0))
    # ---------------------- END READ OBS DATA ----------------------------------------

    # ---------------------- ROMULUS DATA ---------------------------------------------
    objs_rom = pd.read_csv('/home/christenc/Storage/Cosmo/cosmo25/romulus25.m200.data.csv')
    nbins = 10 # 5, for 11.4
    mvir_bins_edges = np.linspace(9.4,13.4,nbins + 1)
    mvir_bins = (mvir_bins_edges[1:] + mvir_bins_edges[:-1])/2 
    mvirz0_bins_edges = np.linspace(9.4,13.4,nbins + 1)
    mvirz0_bins = (mvirz0_bins_edges[1:] + mvirz0_bins_edges[:-1])/2 
    dist_bins = [0,150,300,750,25000]
    mstar_bins0 = np.empty([5,nbins])
    mstar_bins1 = np.empty([5,nbins])
    mstar_bins2 = np.empty([5,nbins])
    mstar_bins3 = np.empty([5,nbins])
    mstarz0_bins0 = np.empty([5,nbins])
    mstarz0_bins1 = np.empty([5,nbins])
    mstarz0_bins2 = np.empty([5,nbins])
    mstarz0_bins3 = np.empty([5,nbins])    
    mbar_bins0 = np.empty([5,nbins])
    mbar_bins1 = np.empty([5,nbins])
    mbar_bins2 = np.empty([5,nbins])
    mbar_bins3 = np.empty([5,nbins])
    fHI_bins0 = np.empty([5,nbins])
    fHI_bins1 = np.empty([5,nbins])
    fHI_bins2 = np.empty([5,nbins])
    fHI_bins3 = np.empty([5,nbins]) 
    dist_binsc = np.empty(len(dist_bins) - 1)
    for i in range(len(dist_bins) - 1):
        dist_binsc[i] = np.median(objs_rom[(objs_rom['massiveDist'] > dist_bins[i]) & (objs_rom['massiveDist'] < dist_bins[i + 1])]['massiveDist'])
    for i in range(nbins):
        if len(objs_rom[(np.log10(objs_rom['Mpeak'])>mvir_bins_edges[i]) & (np.log10(objs_rom['Mpeak'])<mvir_bins_edges[i + 1]) & (objs_rom['massiveDist'] > dist_bins[0]) & (objs_rom['massiveDist'] < dist_bins[1])]['Mstar_Mpeak']) > 0:
            cond = (np.log10(objs_rom['Mpeak'])>mvir_bins_edges[i]) & (np.log10(objs_rom['Mpeak'])<mvir_bins_edges[i + 1]) & (objs_rom['massiveDist'] > dist_bins[0]) & (objs_rom['massiveDist'] < dist_bins[1])
            if (np.sum(cond)) > 3:
                mstar_bins0[:,i] = np.percentile(objs_rom[cond]['Mstar_Mpeak'],[10,25,50,75,90])
                mbar_bins0[:,i] = np.percentile(objs_rom[cond]['mstar'] + objs_rom[cond]['mHI'],[10,25,50,75,90])
                fHI_bins0[:,i] = np.percentile(objs_rom[cond]['mHI']/(objs_rom[cond]['mstar'] + objs_rom[cond]['mHI']),[10,25,50,75,90])
        if len(objs_rom[(np.log10(objs_rom['Mpeak'])>mvir_bins_edges[i]) & (np.log10(objs_rom['Mpeak'])<mvir_bins_edges[i + 1]) & (objs_rom['massiveDist'] > dist_bins[1]) & (objs_rom['massiveDist'] < dist_bins[2])]['Mstar_Mpeak']) > 0:
            cond = (np.log10(objs_rom['Mpeak'])>mvir_bins_edges[i]) & (np.log10(objs_rom['Mpeak'])<mvir_bins_edges[i + 1]) & (objs_rom['massiveDist'] > dist_bins[1]) & (objs_rom['massiveDist'] < dist_bins[2])
            if (np.sum(cond)) > 3:
                mstar_bins1[:,i] = np.percentile(objs_rom[cond]['Mstar_Mpeak'],[10,25,50,75,90])
                mbar_bins1[:,i] = np.percentile(objs_rom[cond]['mstar'] + objs_rom[cond]['mHI'],[10,25,50,75,90])
                fHI_bins1[:,i] = np.percentile(objs_rom[cond]['mHI']/(objs_rom[cond]['mstar'] + objs_rom[cond]['mHI']),[10,25,50,75,90])
        if len(objs_rom[(np.log10(objs_rom['Mpeak'])>mvir_bins_edges[i]) & (np.log10(objs_rom['Mpeak'])<mvir_bins_edges[i + 1]) & (objs_rom['massiveDist'] > dist_bins[2]) & (objs_rom['massiveDist'] < dist_bins[3])]['Mstar_Mpeak']) > 0:
            cond = (np.log10(objs_rom['Mpeak'])>mvir_bins_edges[i]) & (np.log10(objs_rom['Mpeak'])<mvir_bins_edges[i + 1]) & (objs_rom['massiveDist'] > dist_bins[2]) & (objs_rom['massiveDist'] < dist_bins[3])
            if (np.sum(cond)) > 3:
                mstar_bins2[:,i] = np.percentile(objs_rom[cond]['Mstar_Mpeak'],[10,25,50,75,90])
                mbar_bins2[:,i] = np.percentile(objs_rom[cond]['mstar'] + objs_rom[cond]['mHI'],[10,25,50,75,90])
                fHI_bins2[:,i] = np.percentile(objs_rom[cond]['mHI']/(objs_rom[cond]['mstar'] + objs_rom[cond]['mHI']),[10,25,50,75,90])                                                
        if len(objs_rom[(np.log10(objs_rom['Mpeak'])>mvir_bins_edges[i]) & (np.log10(objs_rom['Mpeak'])<mvir_bins_edges[i + 1]) & (objs_rom['massiveDist'] > dist_bins[3]) & (objs_rom['massiveDist'] < dist_bins[4])]['Mstar_Mpeak']) > 0:
            cond = (np.log10(objs_rom['Mpeak'])>mvir_bins_edges[i]) & (np.log10(objs_rom['Mpeak'])<mvir_bins_edges[i + 1]) & (objs_rom['massiveDist'] > dist_bins[3]) & (objs_rom['massiveDist'] < dist_bins[4])
            if (np.sum(cond)) > 3:
                mstar_bins3[:,i] = np.percentile(objs_rom[cond]['Mstar_Mpeak'],[10,25,50,75,90])
                mbar_bins3[:,i] = np.percentile(objs_rom[cond]['mstar'] + objs_rom[cond]['mHI'],[10,25,50,75,90])
                fHI_bins3[:,i] = np.percentile(objs_rom[cond]['mHI']/(objs_rom[cond]['mstar'] + objs_rom[cond]['mHI']),[10,25,50,75,90])
    for i in range(nbins):
        if len(objs_rom[(np.log10(objs_rom['mvir'])>mvirz0_bins_edges[i]) & (np.log10(objs_rom['mvir'])<mvirz0_bins_edges[i + 1]) & (objs_rom['massiveDist'] > dist_bins[0]) & (objs_rom['massiveDist'] < dist_bins[1])]['mstar']) > 0:
            cond = (np.log10(objs_rom['mvir'])>mvirz0_bins_edges[i]) & (np.log10(objs_rom['mvir'])<mvirz0_bins_edges[i + 1]) & (objs_rom['massiveDist'] > dist_bins[0]) & (objs_rom['massiveDist'] < dist_bins[1])
            if (np.sum(cond)) > 3:
                mstarz0_bins0[:,i] = np.percentile(objs_rom[cond]['mstar'],[10,25,50,75,90])
        if len(objs_rom[(np.log10(objs_rom['mvir'])>mvirz0_bins_edges[i]) & (np.log10(objs_rom['mvir'])<mvirz0_bins_edges[i + 1]) & (objs_rom['massiveDist'] > dist_bins[1]) & (objs_rom['massiveDist'] < dist_bins[2])]['mstar']) > 0:
            cond = (np.log10(objs_rom['mvir'])>mvirz0_bins_edges[i]) & (np.log10(objs_rom['mvir'])<mvirz0_bins_edges[i + 1]) & (objs_rom['massiveDist'] > dist_bins[1]) & (objs_rom['massiveDist'] < dist_bins[2])
            if (np.sum(cond)) > 3:
                mstarz0_bins1[:,i] = np.percentile(objs_rom[cond]['mstar'],[10,25,50,75,90])
        if len(objs_rom[(np.log10(objs_rom['mvir'])>mvirz0_bins_edges[i]) & (np.log10(objs_rom['mvir'])<mvirz0_bins_edges[i + 1]) & (objs_rom['massiveDist'] > dist_bins[2]) & (objs_rom['massiveDist'] < dist_bins[3])]['mstar']) > 0:
            cond = (np.log10(objs_rom['mvir'])>mvirz0_bins_edges[i]) & (np.log10(objs_rom['mvir'])<mvirz0_bins_edges[i + 1]) & (objs_rom['massiveDist'] > dist_bins[2]) & (objs_rom['massiveDist'] < dist_bins[3])
            if (np.sum(cond)) > 3:
                mstarz0_bins2[:,i] = np.percentile(objs_rom[cond]['mstar'],[10,25,50,75,90])
        if len(objs_rom[(np.log10(objs_rom['mvir'])>mvirz0_bins_edges[i]) & (np.log10(objs_rom['mvir'])<mvirz0_bins_edges[i + 1]) & (objs_rom['massiveDist'] > dist_bins[3]) & (objs_rom['massiveDist'] < dist_bins[4])]['mstar']) > 0:
            cond = (np.log10(objs_rom['mvir'])>mvirz0_bins_edges[i]) & (np.log10(objs_rom['mvir'])<mvirz0_bins_edges[i + 1]) & (objs_rom['massiveDist'] > dist_bins[3]) & (objs_rom['massiveDist'] < dist_bins[4])
            if (np.sum(cond)) > 3:
                mstarz0_bins3[:,i] = np.percentile(objs_rom[cond]['mstar'],[10,25,50,75,90])
            
    # --------------------- END READ ROMULUS ------------------------------------        


    
    use_peak_mstar = False
    if use_peak_mstar:
        mstar_key = 'Mstar_Mpeak'
        mstar_axes_label = r"M$_{\mathrm{*, peak}}$/M$_\odot$"
        mstar_ratio_axes_label = r"M$_{\mathrm{*, peak}}$/M$_{\mathrm{vir, peak}}$"        
    else:
        mstar_key = 'Mstar_z0'
        mstar_axes_label = r"M$_{*, z = 0}$/M$_\odot$"
        mstar_ratio_axes_label = r"M$_{*, z = 0}$/M$_{\mathrm{vir, peak}}$"

    time_z0_5 = Simpy.cosmology.getTime(0.5) # Find stellar mass to halo mass at these different
    time_z1 = Simpy.cosmology.getTime(1)    # redshifts, only plotting it if the galaxy has not yet
    time_z2 = Simpy.cosmology.getTime(2)    # reached peak halo mass
    time_z4 = Simpy.cosmology.getTime(4)
        
    #SMHM colored by tau_90
    #ig1.clear()
    plt.clf()
    fig1 = plt.figure(1,figsize=(plt_width,plt_width*aspect_ratio))
    fig1.set_size_inches(plt_width,plt_width*aspect_ratio)
    fig1.clear()
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax1 = fig1.add_subplot(gs[0])
    ax1sub = fig1.add_subplot(gs[1])
    cmx = plt.get_cmap("viridis") 
    cNorm  = colors.Normalize(vmin=0, vmax = 14)
    #Read
    ax1.fill_between(read_abunmatch['M200'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_low'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_high'][read_abunmatch['M200']> 1e9],color = 'grey', alpha = 0.5)
    #Nadler outer
    ax1.fill_between(upper_fid_Mhalo, # Mhalo values 
        np.interp( upper_fid_Mhalo, lower_fid_Mhalo, lower_fid_Mstar ), # Mstar lower bound 
        upper_fid_Mstar, # Mstar upper bound  
        facecolor="#33669A",
        linewidth=0,
        zorder=0,
        alpha=0.3
	)
    #Nadler inner
    ax1.fill_between(
        upper_fid_inner_Mhalo, # Mhalo values     
        np.interp( upper_fid_inner_Mhalo, lower_fid_inner_Mhalo, lower_fid_inner_Mstar ), # Mstar lower bound    
        upper_fid_inner_Mstar, # Mstar upper bound 
        facecolor="#33669A",
        linewidth=0,
        zorder=0,
        alpha=0.3
    )
    #Jethwa outer
    ax1.fill_between(
        Jethwa_upper_outer_Mhalo, # Mhalo values  
        np.interp( Jethwa_upper_outer_Mhalo, Jethwa_lower_outer_Mhalo, Jethwa_lower_outer_Mstar ), # Mstar lower bound   
        Jethwa_upper_outer_Mstar, # Mstar upper bound  
        facecolor="#66CC66",
        linewidth=0,
        zorder=0,
        alpha=0.3
    )
    #Jethwa inner
    ax1.fill_between(
        Jethwa_upper_inner_Mhalo, # Mhalo values   
        np.interp( Jethwa_upper_inner_Mhalo, Jethwa_lower_inner_Mhalo, Jethwa_lower_inner_Mstar ), # Mstar lower bound  
        Jethwa_upper_inner_Mstar, # Mstar upper bound 
        facecolor="#66CC66",
        linewidth=0,
        zorder=0,
        alpha=0.3
    )
    read, = ax1.plot( [], [],  color="grey", linestyle="solid", linewidth=4) #, label="Read+ 2017" )
    nadler, = ax1.plot( [], [],  color="#33669A", linestyle="solid", linewidth=8) #, label="Nadler+ 2020" )
    jethwa, = ax1.plot( [], [],  color="#66CC66", linestyle="solid", linewidth=8) #, label="Jethwa+ 2018" )    
    #ax1.plot( Read_upper_dashed_Mhalo, Read_upper_dashed_Mstar,  linestyle='--', color="#FF9966", linewidth=2 )
    #ax1.plot( Read_lower_dashed_Mhalo, Read_lower_dashed_Mstar,  linestyle='--', color="#FF9966", linewidth=2 )
    #ax1.plot( Read_upper_solid_Mhalo, Read_upper_solid_Mstar,  linestyle='solid', color="#FF9966", linewidth=4 )
    #ax1.plot( Read_lower_solid_Mhalo, Read_lower_solid_Mstar,  linestyle='solid', color="#FF9966", linewidth=4 )
    cen_plt = ax1.scatter(objs_pd_comb['mass'][objs_pd_comb['type']=='Central'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Central'],c = objs_pd_comb['tau90'][objs_pd_comb['type']=='Central'], cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    cen_plt = ax1.scatter(objs_pd_comb['mass'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb[mstar_key][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = objs_pd_comb['tau90'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], linewidths = edgewidth)
    sat_plt = ax1.scatter(objs_pd_comb['mass'][objs_pd_comb['type']=='Satellite'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Satellite'],c = objs_pd_comb['tau90'][objs_pd_comb['type']=='Satellite'], cmap = cmx, norm = cNorm, edgecolor = 'k',marker = "*", s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    #ax1.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_ylabel(mstar_axes_label)
    ax1.set_xlabel(r'M$_{\mathrm{vir}}$/M$_\odot$')
    ax1.axis([2e6, 2e11, 2e2, 5e9])
    legend1 = ax1.legend([cen_plt,sat_plt],['Central','Satellite'],scatterpoints = 1,facecolor = 'white',loc = 4,framealpha = 0,frameon = False)
    legend2 = ax1.legend([read,nadler,jethwa],['Read+ 2017','Nadler+ 2020','Jethwa+ 2018'],facecolor = 'white',loc = 0,framealpha = 0,frameon = False)
    ax1.add_artist(legend1)
    ax1.add_artist(legend2)
    cb = mpl.colorbar.ColorbarBase(ax1sub, cmap=cmx, norm=cNorm)
    cb.set_label(r"$\tau_{90}$ (Gyr)")
    fig1.tight_layout()
    fig1.show()
    fig1.savefig(outfile_base + '_SMHM_t90.png',dpi = dpi)

    #SMHM colored by tau_90
    #plt.clf()
    fig1 = plt.figure(1,figsize=(plt_width,plt_width*aspect_ratio))
    fig1.set_size_inches(plt_width,plt_width*aspect_ratio)
    fig1.clear()
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax1 = fig1.add_subplot(gs[0])
    ax1sub = fig1.add_subplot(gs[1])
    cmx = plt.get_cmap("viridis") 
    cNorm  = colors.Normalize(vmin=0, vmax = 14)    
    ax1.fill_between(read_abunmatch['M200'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_low'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_high'][read_abunmatch['M200']> 1e9],color = 'grey',linewidth=0,alpha=0.6) #Read
    ax1.fill_between(upper_fid_Mhalo, np.interp( upper_fid_Mhalo, lower_fid_Mhalo, lower_fid_Mstar ), upper_fid_Mstar, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3)#Nadler outer
    ax1.fill_between(upper_fid_inner_Mhalo, np.interp( upper_fid_inner_Mhalo, lower_fid_inner_Mhalo, lower_fid_inner_Mstar ),upper_fid_inner_Mstar, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3) #Nadler inner
    ax1.fill_between(Jethwa_upper_outer_Mhalo, np.interp( Jethwa_upper_outer_Mhalo, Jethwa_lower_outer_Mhalo, Jethwa_lower_outer_Mstar ),Jethwa_upper_outer_Mstar, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3)#Jethwa outer
    ax1.fill_between(Jethwa_upper_inner_Mhalo,np.interp( Jethwa_upper_inner_Mhalo, Jethwa_lower_inner_Mhalo, Jethwa_lower_inner_Mstar ), Jethwa_upper_inner_Mstar, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3) #Jethwa inner
    read, = ax1.plot( [], [],  color="grey", linestyle="solid", linewidth=4) #, label="Read+ 2017" )
    nadler, = ax1.plot( [], [],  color="#33669A", linestyle="solid", linewidth=8) #, label="Nadler+ 2020" )
    jethwa, = ax1.plot( [], [],  color="#66CC66", linestyle="solid", linewidth=8) #, label="Jethwa+ 2018" ) 
    cen_plt = ax1.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Central'],c = objs_pd_comb['tau90'][objs_pd_comb['type']=='Central'], cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    sat_plt = ax1.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Satellite'],c = objs_pd_comb['tau90'][objs_pd_comb['type']=='Satellite'], cmap = cmx, norm = cNorm,edgecolor = 'k',marker="*", s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    #ax1.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_ylabel(mstar_axes_label)
    ax1.set_xlabel(r'M$_{\mathrm{vir, peak}}$/M$_\odot$')
    ax1.axis([1e8, 2e11, 2e2, 5e9])
    legend1 = ax1.legend([cen_plt,sat_plt],['Central','Satellite'],scatterpoints = 1,facecolor = 'white',loc = 0,framealpha = 0,frameon = False)
    legend2 = ax1.legend([read,nadler,jethwa],['Read+ 2017','Nadler+ 2020','Jethwa+ 2018'],facecolor = 'white',loc = 4,framealpha = 0,frameon = False)
    ax1.add_artist(legend1)
    ax1.add_artist(legend2)    
    cb = mpl.colorbar.ColorbarBase(ax1sub, cmap=cmx, norm=cNorm)
    cb.set_label(r"$\tau_{90}$ (Gyr)")
    fig1.tight_layout()
    fig1.show()    
    fig1.savefig(outfile_base + '_SMHM_t90_Mpeak.png',dpi = dpi)

    fig1 = plt.figure(1,figsize=(plt_width,plt_width*aspect_ratio))
    fig1.set_size_inches(plt_width,plt_width*aspect_ratio)
    fig1.clear()
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax1 = fig1.add_subplot(gs[0])
    ax1sub = fig1.add_subplot(gs[1])
    cmx = plt.get_cmap("viridis") 
    cNorm  = colors.Normalize(vmin=0, vmax = 14)    
    ax1.fill_between(read_abunmatch['M200'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_low'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_high'][read_abunmatch['M200']> 1e9],color = 'grey',linewidth=0,alpha=0.6) #Read
    ax1.fill_between(upper_fid_Mhalo, np.interp( upper_fid_Mhalo, lower_fid_Mhalo, lower_fid_Mstar ), upper_fid_Mstar, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3)#Nadler outer
    ax1.fill_between(upper_fid_inner_Mhalo, np.interp( upper_fid_inner_Mhalo, lower_fid_inner_Mhalo, lower_fid_inner_Mstar ),upper_fid_inner_Mstar, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3) #Nadler inner
    ax1.fill_between(Jethwa_upper_outer_Mhalo, np.interp( Jethwa_upper_outer_Mhalo, Jethwa_lower_outer_Mhalo, Jethwa_lower_outer_Mstar ),Jethwa_upper_outer_Mstar, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3)#Jethwa outer
    ax1.fill_between(Jethwa_upper_inner_Mhalo,np.interp( Jethwa_upper_inner_Mhalo, Jethwa_lower_inner_Mhalo, Jethwa_lower_inner_Mstar ), Jethwa_upper_inner_Mstar, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3) #Jethwa inner
    read, = ax1.plot( [], [],  color="grey", linestyle="solid", linewidth=4) #, label="Read+ 2017" )
    nadler, = ax1.plot( [], [],  color="#33669A", linestyle="solid", linewidth=8) #, label="Nadler+ 2020" )
    jethwa, = ax1.plot( [], [],  color="#66CC66", linestyle="solid", linewidth=8) #, label="Jethwa+ 2018" ) 
    cen_plt = ax1.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Central'],c = objs_pd_comb['tau90_vir'][objs_pd_comb['type']=='Central'], cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax1.scatter(objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb[mstar_key][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = objs_pd_comb['tau90_vir'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], linewidths = edgewidth)
    sat_plt = ax1.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Satellite'],c = objs_pd_comb['tau90_vir'][objs_pd_comb['type']=='Satellite'], cmap = cmx, norm = cNorm,edgecolor = 'k',marker="*", s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    #ax1.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_ylabel(mstar_axes_label)
    ax1.set_xlabel(r'M$_{\mathrm{vir, peak}}$/M$_\odot$')
    ax1.axis([1e8, 2e11, 2e2, 5e9])
    legend1 = ax1.legend([cen_plt,sat_plt],['Central','Satellite'],scatterpoints = 1,facecolor = 'white',loc = 0,framealpha = 0,frameon = False)
    legend2 = ax1.legend([read,nadler,jethwa],['Read+ 2017','Nadler+ 2020','Jethwa+ 2018'],facecolor = 'white',loc = 4,framealpha = 0,frameon = False)
    ax1.add_artist(legend1)
    ax1.add_artist(legend2)    
    cb = mpl.colorbar.ColorbarBase(ax1sub, cmap=cmx, norm=cNorm)
    cb.set_label(r"$\tau_{90, vir}$ (Gyr)")
    fig1.tight_layout()
    fig1.show()    
    fig1.savefig(outfile_base + '_SMHM_t90mvir_Mpeak.png',dpi = dpi)

    #SMHM colored by tau_90
    #plt.clf()
    fig1 = plt.figure(1,figsize=(plt_width*2,plt_width*aspect_ratio))
    fig1.set_size_inches(plt_width*2,plt_width*aspect_ratio)
    fig1.clear()
    gs = fig1.add_gridspec(1,2,wspace=0)
    axs1 = gs.subplots(sharey = True) #, constrained_layout=True) 
    axs1 = axs1.flatten()
    ax1a = axs1[0]
    ax1b = axs1[1]       
    cmx = plt.get_cmap("cool_r") 
    cNorm  = colors.Normalize(vmin=0, vmax = 14)    
    ax1a.fill_between(read_abunmatch['M200'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_low'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_high'][read_abunmatch['M200']> 1e9],color = 'grey',linewidth=0,alpha=0.6) #Read
    ax1a.fill_between(upper_fid_Mhalo, np.interp( upper_fid_Mhalo, lower_fid_Mhalo, lower_fid_Mstar ), upper_fid_Mstar, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3)#Nadler outer
    ax1a.fill_between(upper_fid_inner_Mhalo, np.interp( upper_fid_inner_Mhalo, lower_fid_inner_Mhalo, lower_fid_inner_Mstar ),upper_fid_inner_Mstar, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3) #Nadler inner
    ax1a.fill_between(Jethwa_upper_outer_Mhalo, np.interp( Jethwa_upper_outer_Mhalo, Jethwa_lower_outer_Mhalo, Jethwa_lower_outer_Mstar ),Jethwa_upper_outer_Mstar, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3)#Jethwa outer
    ax1a.fill_between(Jethwa_upper_inner_Mhalo,np.interp( Jethwa_upper_inner_Mhalo, Jethwa_lower_inner_Mhalo, Jethwa_lower_inner_Mstar ), Jethwa_upper_inner_Mstar, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3) #Jethwa inner
    read, = ax1a.plot( [], [],  color="grey", linestyle="solid", linewidth=4) #, label="Read+ 2017" )
    nadler, = ax1a.plot( [], [],  color="#33669A", linestyle="solid", linewidth=8) #, label="Nadler+ 2020" )
    jethwa, = ax1a.plot( [], [],  color="#66CC66", linestyle="solid", linewidth=8) #, label="Jethwa+ 2018" ) 
    cen_plt = ax1a.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Central'],c = objs_pd_comb['tau90'][objs_pd_comb['type']=='Central'], cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    sat_plt = ax1a.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Satellite'],c = objs_pd_comb['tau90'][objs_pd_comb['type']=='Satellite'], cmap = cmx, norm = cNorm,edgecolor = 'k',marker="*", s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    #ax1.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax1a.set_xscale('log')
    ax1a.set_yscale('log')
    ax1a.set_ylabel(mstar_axes_label)
    ax1a.set_xlabel(r'M$_{\mathrm{vir, peak}}$/M$_\odot$')
    ax1a.axis([1e8, 2e11, 2e2, 5e9])
    legend1 = ax1a.legend([read,nadler,jethwa],['Read+ 2017','Nadler+ 2020','Jethwa+ 2018'],facecolor = 'white',loc = 4,framealpha = 0,frameon = False)
    ax1a.add_artist(legend1)
    ax1a.text(2e8,7e8,r'Stellar $\tau$')

    ax1b.fill_between(read_abunmatch['M200'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_low'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_high'][read_abunmatch['M200']> 1e9],color = 'grey',linewidth=0,alpha=0.6) #Read
    ax1b.fill_between(upper_fid_Mhalo, np.interp( upper_fid_Mhalo, lower_fid_Mhalo, lower_fid_Mstar ), upper_fid_Mstar, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3)#Nadler outer
    ax1b.fill_between(upper_fid_inner_Mhalo, np.interp( upper_fid_inner_Mhalo, lower_fid_inner_Mhalo, lower_fid_inner_Mstar ),upper_fid_inner_Mstar, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3) #Nadler inner
    ax1b.fill_between(Jethwa_upper_outer_Mhalo, np.interp( Jethwa_upper_outer_Mhalo, Jethwa_lower_outer_Mhalo, Jethwa_lower_outer_Mstar ),Jethwa_upper_outer_Mstar, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3)#Jethwa outer
    ax1b.fill_between(Jethwa_upper_inner_Mhalo,np.interp( Jethwa_upper_inner_Mhalo, Jethwa_lower_inner_Mhalo, Jethwa_lower_inner_Mstar ), Jethwa_upper_inner_Mstar, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3) #Jethwa inner
    read, = ax1b.plot( [], [],  color="grey", linestyle="solid", linewidth=4) #, label="Read+ 2017" )
    nadler, = ax1b.plot( [], [],  color="#33669A", linestyle="solid", linewidth=8) #, label="Nadler+ 2020" )
    jethwa, = ax1b.plot( [], [],  color="#66CC66", linestyle="solid", linewidth=8) #, label="Jethwa+ 2018" ) 
    cen_plt = ax1b.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Central'],c = objs_pd_comb['tau90_vir'][objs_pd_comb['type']=='Central'], cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax1b.scatter(objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb[mstar_key][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = objs_pd_comb['tau90_vir'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], linewidths = edgewidth)
    sat_plt = ax1b.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Satellite'],c = objs_pd_comb['tau90_vir'][objs_pd_comb['type']=='Satellite'], cmap = cmx, norm = cNorm,edgecolor = 'k',marker="*", s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    #ax1.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax1b.set_xscale('log')
    ax1b.set_yscale('log')
    ax1b.set_xlabel(r'M$_{\mathrm{vir, peak}}$/M$_\odot$')
    ax1b.axis([1e8, 2e11, 2e2, 5e9])
    legend1 = ax1b.legend([cen_plt,sat_plt],['Central','Satellite/Backsplash'],scatterpoints = 1,facecolor = 'white',loc = 4,framealpha = 0,frameon = False)
    ax1b.add_artist(legend1)
    ax1b.text(2e8,7e8,r'Virial Mass $\tau$')
    fig1.tight_layout()
    sm = plt.cm.ScalarMappable(cmap=cmx, norm=cNorm)
    cbar = plt.colorbar(sm, location="right",ax=[ax1a,ax1b],aspect = 20)  
    cbar.set_label(r"$\tau_{90}$ (Gyr)")
    fig1.show()
    
    fig1 = plt.figure(1,figsize=(plt_width,plt_width*aspect_ratio))
    fig1.set_size_inches(plt_width,plt_width*aspect_ratio)
    fig1.clear()
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax1 = fig1.add_subplot(gs[0])
    ax1sub = fig1.add_subplot(gs[1])
    cmx = plt.get_cmap("viridis") 
    cNorm  = colors.Normalize(vmin=0, vmax = 14)
    #Read
    ax1.fill_between(read_abunmatch['M200'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_low'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_high'][read_abunmatch['M200']> 1e9],color = 'grey', alpha = 0.5)
    #Nadler outer
    ax1.fill_between(upper_fid_Mhalo, # Mhalo values 
        np.interp( upper_fid_Mhalo, lower_fid_Mhalo, lower_fid_Mstar ), # Mstar lower bound 
        upper_fid_Mstar, # Mstar upper bound  
        facecolor="#33669A",
        linewidth=0,
        zorder=0,
        alpha=0.3
	)
    #Nadler inner
    ax1.fill_between(
        upper_fid_inner_Mhalo, # Mhalo values     
        np.interp( upper_fid_inner_Mhalo, lower_fid_inner_Mhalo, lower_fid_inner_Mstar ), # Mstar lower bound    
        upper_fid_inner_Mstar, # Mstar upper bound 
        facecolor="#33669A",
        linewidth=0,
        zorder=0,
        alpha=0.3
    )
    #Jethwa outer
    ax1.fill_between(
        Jethwa_upper_outer_Mhalo, # Mhalo values  
        np.interp( Jethwa_upper_outer_Mhalo, Jethwa_lower_outer_Mhalo, Jethwa_lower_outer_Mstar ), # Mstar lower bound   
        Jethwa_upper_outer_Mstar, # Mstar upper bound  
        facecolor="#66CC66",
        linewidth=0,
        zorder=0,
        alpha=0.3
    )
    #Jethwa inner
    ax1.fill_between(
        Jethwa_upper_inner_Mhalo, # Mhalo values   
        np.interp( Jethwa_upper_inner_Mhalo, Jethwa_lower_inner_Mhalo, Jethwa_lower_inner_Mstar ), # Mstar lower bound  
        Jethwa_upper_inner_Mstar, # Mstar upper bound 
        facecolor="#66CC66",
        linewidth=0,
        zorder=0,
        alpha=0.3
    )
    read, = ax1.plot( [], [],  color="grey", linestyle="solid", linewidth=4) #, label="Read+ 2017" )
    nadler, = ax1.plot( [], [],  color="#33669A", linestyle="solid", linewidth=8) #, label="Nadler+ 2020" )
    jethwa, = ax1.plot( [], [],  color="#66CC66", linestyle="solid", linewidth=8) #, label="Jethwa+ 2018" )    
    #ax1.plot( Read_upper_dashed_Mhalo, Read_upper_dashed_Mstar,  linestyle='--', color="#FF9966", linewidth=2 )
    #ax1.plot( Read_lower_dashed_Mhalo, Read_lower_dashed_Mstar,  linestyle='--', color="#FF9966", linewidth=2 )
    #ax1.plot( Read_upper_solid_Mhalo, Read_upper_solid_Mstar,  linestyle='solid', color="#FF9966", linewidth=4 )
    #ax1.plot( Read_lower_solid_Mhalo, Read_lower_solid_Mstar,  linestyle='solid', color="#FF9966", linewidth=4 )
    cen_plt = ax1.scatter(objs_pd_comb['mass'][objs_pd_comb['type']=='Central'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Central'],c = objs_pd_comb['tau90'][objs_pd_comb['type']=='Central'], cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    cen_plt = ax1.scatter(objs_pd_comb['mass'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb[mstar_key][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = objs_pd_comb['tau90'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], linewidths = edgewidth)
    sat_plt = ax1.scatter(objs_pd_comb['mass'][objs_pd_comb['type']=='Satellite'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Satellite'],c = objs_pd_comb['tau90'][objs_pd_comb['type']=='Satellite'], cmap = cmx, norm = cNorm, edgecolor = 'k',marker = "*", s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    #ax1.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_ylabel(mstar_axes_label)
    ax1.set_xlabel(r'M$_{\mathrm{vir}}$/M$_\odot$')
    ax1.axis([2e6, 2e11, 2e2, 5e9])
    legend1 = ax1.legend([cen_plt,sat_plt],['Central','Satellite'],scatterpoints = 1,facecolor = 'white',loc = 4,framealpha = 0,frameon = False)
    legend2 = ax1.legend([read,nadler,jethwa],['Read+ 2017','Nadler+ 2020','Jethwa+ 2018'],facecolor = 'white',loc = 0,framealpha = 0,frameon = False)
    ax1.add_artist(legend1)
    ax1.add_artist(legend2)
    cb = mpl.colorbar.ColorbarBase(ax1sub, cmap=cmx, norm=cNorm)
    cb.set_label(r"$\tau_{90}$ (Gyr)")
    fig1.tight_layout()
    fig1.show()
    fig1.savefig(outfile_base + '_SMHM_t90.png',dpi = dpi)

    #SMHM colored by tau_90
    #plt.clf()
    fig1 = plt.figure(1,figsize=(plt_width,plt_width*aspect_ratio))
    fig1.set_size_inches(plt_width,plt_width*aspect_ratio)
    fig1.clear()
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax1 = fig1.add_subplot(gs[0])
    ax1sub = fig1.add_subplot(gs[1])
    cmx = plt.get_cmap("viridis") 
    cNorm  = colors.Normalize(vmin=0, vmax = 14)    
    ax1.fill_between(read_abunmatch['M200'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_low'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_high'][read_abunmatch['M200']> 1e9],color = 'grey',linewidth=0,alpha=0.6) #Read
    ax1.fill_between(upper_fid_Mhalo, np.interp( upper_fid_Mhalo, lower_fid_Mhalo, lower_fid_Mstar ), upper_fid_Mstar, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3)#Nadler outer
    ax1.fill_between(upper_fid_inner_Mhalo, np.interp( upper_fid_inner_Mhalo, lower_fid_inner_Mhalo, lower_fid_inner_Mstar ),upper_fid_inner_Mstar, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3) #Nadler inner
    ax1.fill_between(Jethwa_upper_outer_Mhalo, np.interp( Jethwa_upper_outer_Mhalo, Jethwa_lower_outer_Mhalo, Jethwa_lower_outer_Mstar ),Jethwa_upper_outer_Mstar, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3)#Jethwa outer
    ax1.fill_between(Jethwa_upper_inner_Mhalo,np.interp( Jethwa_upper_inner_Mhalo, Jethwa_lower_inner_Mhalo, Jethwa_lower_inner_Mstar ), Jethwa_upper_inner_Mstar, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3) #Jethwa inner
    read, = ax1.plot( [], [],  color="grey", linestyle="solid", linewidth=4) #, label="Read+ 2017" )
    nadler, = ax1.plot( [], [],  color="#33669A", linestyle="solid", linewidth=8) #, label="Nadler+ 2020" )
    jethwa, = ax1.plot( [], [],  color="#66CC66", linestyle="solid", linewidth=8) #, label="Jethwa+ 2018" ) 
    cen_plt = ax1.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Central'],c = objs_pd_comb['tau90'][objs_pd_comb['type']=='Central'], cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    sat_plt = ax1.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Satellite'],c = objs_pd_comb['tau90'][objs_pd_comb['type']=='Satellite'], cmap = cmx, norm = cNorm,edgecolor = 'k',marker="*", s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    #ax1.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_ylabel(mstar_axes_label)
    ax1.set_xlabel(r'M$_{\mathrm{vir, peak}}$/M$_\odot$')
    ax1.axis([1e8, 2e11, 2e2, 5e9])
    legend1 = ax1.legend([cen_plt,sat_plt],['Central','Satellite'],scatterpoints = 1,facecolor = 'white',loc = 0,framealpha = 0,frameon = False)
    legend2 = ax1.legend([read,nadler,jethwa],['Read+ 2017','Nadler+ 2020','Jethwa+ 2018'],facecolor = 'white',loc = 4,framealpha = 0,frameon = False)
    ax1.add_artist(legend1)
    ax1.add_artist(legend2)    
    cb = mpl.colorbar.ColorbarBase(ax1sub, cmap=cmx, norm=cNorm)
    cb.set_label(r"$\tau_{90}$ (Gyr)")
    fig1.tight_layout()
    fig1.show()    
    fig1.savefig(outfile_base + '_SMHM_t90_Mpeak.png',dpi = dpi)

    fig1 = plt.figure(1,figsize=(plt_width,plt_width*aspect_ratio))
    fig1.set_size_inches(plt_width,plt_width*aspect_ratio)
    fig1.clear()
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax1 = fig1.add_subplot(gs[0])
    ax1sub = fig1.add_subplot(gs[1])
    cmx = plt.get_cmap("viridis") 
    cNorm  = colors.Normalize(vmin=0, vmax = 14)    
    ax1.fill_between(read_abunmatch['M200'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_low'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_high'][read_abunmatch['M200']> 1e9],color = 'grey',linewidth=0,alpha=0.6) #Read
    ax1.fill_between(upper_fid_Mhalo, np.interp( upper_fid_Mhalo, lower_fid_Mhalo, lower_fid_Mstar ), upper_fid_Mstar, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3)#Nadler outer
    ax1.fill_between(upper_fid_inner_Mhalo, np.interp( upper_fid_inner_Mhalo, lower_fid_inner_Mhalo, lower_fid_inner_Mstar ),upper_fid_inner_Mstar, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3) #Nadler inner
    ax1.fill_between(Jethwa_upper_outer_Mhalo, np.interp( Jethwa_upper_outer_Mhalo, Jethwa_lower_outer_Mhalo, Jethwa_lower_outer_Mstar ),Jethwa_upper_outer_Mstar, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3)#Jethwa outer
    ax1.fill_between(Jethwa_upper_inner_Mhalo,np.interp( Jethwa_upper_inner_Mhalo, Jethwa_lower_inner_Mhalo, Jethwa_lower_inner_Mstar ), Jethwa_upper_inner_Mstar, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3) #Jethwa inner
    read, = ax1.plot( [], [],  color="grey", linestyle="solid", linewidth=4) #, label="Read+ 2017" )
    nadler, = ax1.plot( [], [],  color="#33669A", linestyle="solid", linewidth=8) #, label="Nadler+ 2020" )
    jethwa, = ax1.plot( [], [],  color="#66CC66", linestyle="solid", linewidth=8) #, label="Jethwa+ 2018" ) 
    cen_plt = ax1.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Central'],c = objs_pd_comb['tau90_vir'][objs_pd_comb['type']=='Central'], cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax1.scatter(objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb[mstar_key][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = objs_pd_comb['tau90_vir'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], linewidths = edgewidth)
    sat_plt = ax1.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Satellite'],c = objs_pd_comb['tau90_vir'][objs_pd_comb['type']=='Satellite'], cmap = cmx, norm = cNorm,edgecolor = 'k',marker="*", s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    #ax1.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_ylabel(mstar_axes_label)
    ax1.set_xlabel(r'M$_{\mathrm{vir, peak}}$/M$_\odot$')
    ax1.axis([1e8, 2e11, 2e2, 5e9])
    legend1 = ax1.legend([cen_plt,sat_plt],['Central','Satellite'],scatterpoints = 1,facecolor = 'white',loc = 0,framealpha = 0,frameon = False)
    legend2 = ax1.legend([read,nadler,jethwa],['Read+ 2017','Nadler+ 2020','Jethwa+ 2018'],facecolor = 'white',loc = 4,framealpha = 0,frameon = False)
    ax1.add_artist(legend1)
    ax1.add_artist(legend2)    
    cb = mpl.colorbar.ColorbarBase(ax1sub, cmap=cmx, norm=cNorm)
    cb.set_label(r"$\tau_{90, vir}$ (Gyr)")
    fig1.tight_layout()
    fig1.show()    
    fig1.savefig(outfile_base + '_SMHM_t90mvir_Mpeak.png',dpi = dpi)

    #SMHM colored by tau_90
    #plt.clf()
    plt.close('all')
    fig1 = plt.figure(1,figsize=(plt_width*2,plt_width*aspect_ratio))
    gs = fig1.add_gridspec(1,2,wspace=0)
    axs1 = gs.subplots(sharey = True) #, constrained_layout=True) 
    axs1 = axs1.flatten()
    ax1a = axs1[0]
    ax1b = axs1[1]       
    cmx = plt.get_cmap("cool_r") 
    cNorm  = colors.Normalize(vmin=0, vmax = 14)    
    ax1a.fill_between(read_abunmatch['M200'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_low'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_high'][read_abunmatch['M200']> 1e9],color = 'grey',linewidth=0,alpha=0.6) #Read
    ax1a.fill_between(upper_fid_Mhalo, np.interp( upper_fid_Mhalo, lower_fid_Mhalo, lower_fid_Mstar ), upper_fid_Mstar, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3)#Nadler outer
    ax1a.fill_between(upper_fid_inner_Mhalo, np.interp( upper_fid_inner_Mhalo, lower_fid_inner_Mhalo, lower_fid_inner_Mstar ),upper_fid_inner_Mstar, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3) #Nadler inner
    ax1a.fill_between(Jethwa_upper_outer_Mhalo, np.interp( Jethwa_upper_outer_Mhalo, Jethwa_lower_outer_Mhalo, Jethwa_lower_outer_Mstar ),Jethwa_upper_outer_Mstar, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3)#Jethwa outer
    ax1a.fill_between(Jethwa_upper_inner_Mhalo,np.interp( Jethwa_upper_inner_Mhalo, Jethwa_lower_inner_Mhalo, Jethwa_lower_inner_Mstar ), Jethwa_upper_inner_Mstar, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3) #Jethwa inner
    read, = ax1a.plot( [], [],  color="grey", linestyle="solid", linewidth=4) #, label="Read+ 2017" )
    nadler, = ax1a.plot( [], [],  color="#33669A", linestyle="solid", linewidth=8) #, label="Nadler+ 2020" )
    jethwa, = ax1a.plot( [], [],  color="#66CC66", linestyle="solid", linewidth=8) #, label="Jethwa+ 2018" ) 
    cen_plt = ax1a.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Central'],c = objs_pd_comb['tau90'][objs_pd_comb['type']=='Central'], cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    sat_plt = ax1a.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Satellite'],c = objs_pd_comb['tau90'][objs_pd_comb['type']=='Satellite'], cmap = cmx, norm = cNorm,edgecolor = 'k',marker="*", s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    #ax1.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax1a.set_xscale('log')
    ax1a.set_yscale('log')
    ax1a.set_ylabel(mstar_axes_label)
    ax1a.set_xlabel(r'M$_{\mathrm{vir, peak}}$/M$_\odot$')
    ax1a.axis([1e8, 2e11, 2e2, 5e9])
    legend1 = ax1a.legend([read,nadler,jethwa],['Read+ 2017','Nadler+ 2020','Jethwa+ 2018'],facecolor = 'white',loc = 4,framealpha = 0,frameon = False)
    ax1a.add_artist(legend1)
    ax1a.text(2e8,7e8,r'Stellar $\tau$')

    ax1b.fill_between(read_abunmatch['M200'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_low'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_high'][read_abunmatch['M200']> 1e9],color = 'grey',linewidth=0,alpha=0.6) #Read
    ax1b.fill_between(upper_fid_Mhalo, np.interp( upper_fid_Mhalo, lower_fid_Mhalo, lower_fid_Mstar ), upper_fid_Mstar, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3)#Nadler outer
    ax1b.fill_between(upper_fid_inner_Mhalo, np.interp( upper_fid_inner_Mhalo, lower_fid_inner_Mhalo, lower_fid_inner_Mstar ),upper_fid_inner_Mstar, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3) #Nadler inner
    ax1b.fill_between(Jethwa_upper_outer_Mhalo, np.interp( Jethwa_upper_outer_Mhalo, Jethwa_lower_outer_Mhalo, Jethwa_lower_outer_Mstar ),Jethwa_upper_outer_Mstar, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3)#Jethwa outer
    ax1b.fill_between(Jethwa_upper_inner_Mhalo,np.interp( Jethwa_upper_inner_Mhalo, Jethwa_lower_inner_Mhalo, Jethwa_lower_inner_Mstar ), Jethwa_upper_inner_Mstar, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3) #Jethwa inner
    read, = ax1b.plot( [], [],  color="grey", linestyle="solid", linewidth=4) #, label="Read+ 2017" )
    nadler, = ax1b.plot( [], [],  color="#33669A", linestyle="solid", linewidth=8) #, label="Nadler+ 2020" )
    jethwa, = ax1b.plot( [], [],  color="#66CC66", linestyle="solid", linewidth=8) #, label="Jethwa+ 2018" ) 
    cen_plt = ax1b.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Central'],c = objs_pd_comb['tau90_vir'][objs_pd_comb['type']=='Central'], cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax1b.scatter(objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb[mstar_key][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = objs_pd_comb['tau90_vir'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], linewidths = edgewidth)
    sat_plt = ax1b.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Satellite'],c = objs_pd_comb['tau90_vir'][objs_pd_comb['type']=='Satellite'], cmap = cmx, norm = cNorm,edgecolor = 'k',marker="*", s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    #ax1.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax1b.set_xscale('log')
    ax1b.set_yscale('log')
    ax1b.set_xlabel(r'M$_{\mathrm{vir, peak}}$/M$_\odot$')
    ax1b.axis([1e8, 2e11, 2e2, 5e9])
    legend1 = ax1b.legend([cen_plt,sat_plt],['Central','Satellite/Backsplash'],scatterpoints = 1,facecolor = 'white',loc = 4,framealpha = 0,frameon = False)
    ax1b.add_artist(legend1)
    ax1b.text(2e8,7e8,r'Virial Mass $\tau$')
    fig1.tight_layout()
    sm = plt.cm.ScalarMappable(cmap=cmx, norm=cNorm)
    cbar = plt.colorbar(sm, location="right",ax=[ax1a,ax1b],aspect = 20)  
    cbar.set_label(r"$\tau_{90}$ (Gyr)")
    fig1.show()
    fig1.set_size_inches(plt_width*2,plt_width*aspect_ratio)
    fig1.show()
    fig1.savefig(outfile_base + '_SMHM_t90_t90mvir_Mpeak.png',dpi = dpi)
    
    fig1 = plt.figure(1,figsize=(plt_width,plt_width*aspect_ratio))
    fig1.set_size_inches(plt_width,plt_width*aspect_ratio)
    fig1.clear()
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax1 = fig1.add_subplot(gs[0])
    ax1sub = fig1.add_subplot(gs[1])
    cmx = plt.get_cmap("viridis") 
    cNorm  = colors.Normalize(vmin=0, vmax = 10)    
    ax1.fill_between(read_abunmatch['M200'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_low'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_high'][read_abunmatch['M200']> 1e9],color = 'grey',linewidth=0,alpha=0.6) #Read
    ax1.fill_between(upper_fid_Mhalo, np.interp( upper_fid_Mhalo, lower_fid_Mhalo, lower_fid_Mstar ), upper_fid_Mstar, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3)#Nadler outer
    ax1.fill_between(upper_fid_inner_Mhalo, np.interp( upper_fid_inner_Mhalo, lower_fid_inner_Mhalo, lower_fid_inner_Mstar ),upper_fid_inner_Mstar, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3) #Nadler inner
    ax1.fill_between(Jethwa_upper_outer_Mhalo, np.interp( Jethwa_upper_outer_Mhalo, Jethwa_lower_outer_Mhalo, Jethwa_lower_outer_Mstar ),Jethwa_upper_outer_Mstar, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3)#Jethwa outer
    ax1.fill_between(Jethwa_upper_inner_Mhalo,np.interp( Jethwa_upper_inner_Mhalo, Jethwa_lower_inner_Mhalo, Jethwa_lower_inner_Mstar ), Jethwa_upper_inner_Mstar, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3) #Jethwa inner
    read, = ax1.plot( [], [],  color="grey", linestyle="solid", linewidth=4) #, label="Read+ 2017" )
    nadler, = ax1.plot( [], [],  color="#33669A", linestyle="solid", linewidth=8) #, label="Nadler+ 2020" )
    jethwa, = ax1.plot( [], [],  color="#66CC66", linestyle="solid", linewidth=8) #, label="Jethwa+ 2018" ) 
    cen_plt = ax1.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Central'],c = objs_pd_comb['tau50_vir'][objs_pd_comb['type']=='Central'], cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax1.scatter(objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb[mstar_key][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = objs_pd_comb['tau50_vir'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], linewidths = edgewidth)
    sat_plt = ax1.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Satellite'],c = objs_pd_comb['tau50_vir'][objs_pd_comb['type']=='Satellite'], cmap = cmx, norm = cNorm,edgecolor = 'k',marker="*", s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    #ax1.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_ylabel(mstar_axes_label)
    ax1.set_xlabel(r'M$_{\mathrm{vir, peak}}$/M$_\odot$')
    ax1.axis([1e8, 2e11, 2e2, 5e9])
    legend1 = ax1.legend([cen_plt,sat_plt],['Central','Satellite'],scatterpoints = 1,facecolor = 'white',loc = 0,framealpha = 0,frameon = False)
    legend2 = ax1.legend([read,nadler,jethwa],['Read+ 2017','Nadler+ 2020','Jethwa+ 2018'],facecolor = 'white',loc = 4,framealpha = 0,frameon = False)
    ax1.add_artist(legend1)
    ax1.add_artist(legend2)    
    cb = mpl.colorbar.ColorbarBase(ax1sub, cmap=cmx, norm=cNorm)
    cb.set_label(r"$\tau_{50, vir}$ (Gyr)")
    fig1.tight_layout()
    fig1.show()    
    fig1.savefig(outfile_base + '_SMHM_t50mvir_Mpeak.png',dpi = dpi)        
    
    #SMHM colored by distance to massive galaxy
    #plt.clf()
    fig2 = plt.figure(2,figsize=(plt_width,plt_width*aspect_ratio))
    fig2.set_size_inches(plt_width,plt_width*aspect_ratio)
    fig2.clear()
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax2 = fig2.add_subplot(gs[0])
    ax2sub = fig2.add_subplot(gs[1])
    cNorm  = colors.LogNorm(vmin=10**1.5, vmax = 10**4)
    ax2.fill_between(read_abunmatch['M200'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_low'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_high'][read_abunmatch['M200']> 1e9],color = 'grey',linewidth=0,alpha=0.6)
    ax2.fill_between(upper_fid_Mhalo, np.interp( upper_fid_Mhalo, lower_fid_Mhalo, lower_fid_Mstar ), upper_fid_Mstar, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3)#Nadler outer
    ax2.fill_between(upper_fid_inner_Mhalo, np.interp( upper_fid_inner_Mhalo, lower_fid_inner_Mhalo, lower_fid_inner_Mstar ),upper_fid_inner_Mstar, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3) #Nadler inner
    ax2.fill_between(Jethwa_upper_outer_Mhalo, np.interp( Jethwa_upper_outer_Mhalo, Jethwa_lower_outer_Mhalo, Jethwa_lower_outer_Mstar ),Jethwa_upper_outer_Mstar, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3)#Jethwa outer
    ax2.fill_between(Jethwa_upper_inner_Mhalo,np.interp( Jethwa_upper_inner_Mhalo, Jethwa_lower_inner_Mhalo, Jethwa_lower_inner_Mstar ), Jethwa_upper_inner_Mstar, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3) #Jethwa inner
    read, = ax1.plot( [], [],  color="grey", linestyle="solid", linewidth=4) #, label="Read+ 2017" )
    nadler, = ax1.plot( [], [],  color="#33669A", linestyle="solid", linewidth=8) #, label="Nadler+ 2020" )
    jethwa, = ax1.plot( [], [],  color="#66CC66", linestyle="solid", linewidth=8) #, label="Jethwa+ 2018" )     
    ax2.scatter(objs_pd_comb['mass'][objs_pd_comb['type']=='Central'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Central'],c = (objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax2.scatter(objs_pd_comb['mass'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb[mstar_key][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = (objs_pd_comb['massiveDist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], linewidths = edgewidth)    
    ax2.scatter(objs_pd_comb['mass'][objs_pd_comb['type']=='Satellite'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Satellite'],c = (objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k',marker = '*', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    #ax2.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_ylabel(mstar_axes_label)
    ax2.set_xlabel(r'M$_{\mathrm{vir}}$/M$_\odot$')
    ax2.axis([2e6, 2e11, 2e2, 5e9])
    legend1 = ax2.legend([cen_plt,sat_plt],['Central','Satellite'],scatterpoints = 1,facecolor = 'white',loc = 4,framealpha = 0,frameon = False)
    legend2 = ax2.legend([read,nadler,jethwa],['Read+ 2017','Nadler+ 2020','Jethwa+ 2018'],facecolor = 'white',loc = 0,framealpha = 0,frameon = False)
    ax2.add_artist(legend1)
    ax2.add_artist(legend2)    
    cb = mpl.colorbar.ColorbarBase(ax2sub, cmap=cmx, norm=cNorm)
    cb.set_label(r"D$_{\mathrm{massive}}$ [kpc]")
    fig2.tight_layout()
    fig2.show()    
    fig2.savefig(outfile_base + '_SMHM_rMassGal.png',dpi = dpi)

    #SMHM colored by distance to massive galaxy
    #plt.clf()
    #fig2.clear()
    plt.close('all')
    fig2 = plt.figure(2,figsize=(plt_width,plt_width*aspect_ratio))
    fig2.set_size_inches(plt_width,plt_width*aspect_ratio)
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax2 = fig2.add_subplot(gs[0])
    ax2sub = fig2.add_subplot(gs[1])
    cNorm  = colors.LogNorm(vmin=10**1.5, vmax = 10**4)
    ax2.fill_between(read_abunmatch['M200'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_low'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_high'][read_abunmatch['M200']> 1e9],color = 'grey',linewidth=0,alpha=0.6)
    ax2.fill_between(upper_fid_Mhalo, np.interp( upper_fid_Mhalo, lower_fid_Mhalo, lower_fid_Mstar ), upper_fid_Mstar, facecolor='y',linewidth=0,zorder=0,alpha=0.3)#Nadler outer  facecolor="#33669A"
    ax2.fill_between(upper_fid_inner_Mhalo, np.interp( upper_fid_inner_Mhalo, lower_fid_inner_Mhalo, lower_fid_inner_Mstar ),upper_fid_inner_Mstar, facecolor='y',linewidth=0,zorder=0,alpha=0.3) #Nadler inner  facecolor="#33669A"
    ax2.fill_between(Jethwa_upper_outer_Mhalo, np.interp( Jethwa_upper_outer_Mhalo, Jethwa_lower_outer_Mhalo, Jethwa_lower_outer_Mstar ),Jethwa_upper_outer_Mstar, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3)#Jethwa outer
    ax2.fill_between(Jethwa_upper_inner_Mhalo,np.interp( Jethwa_upper_inner_Mhalo, Jethwa_lower_inner_Mhalo, Jethwa_lower_inner_Mstar ), Jethwa_upper_inner_Mstar, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3) #Jethwa inner
    read, = ax2.plot( [], [],  color="grey", linestyle="solid", linewidth=4) #, label="Read+ 2017" )
    nadler, = ax2.plot( [], [],  color="#33669A", linestyle="solid", linewidth=8) #, label="Nadler+ 2020" )
    jethwa, = ax2.plot( [], [],  color="#66CC66", linestyle="solid", linewidth=8) #, label="Jethwa+ 2018" )    
    ax2.scatter(objs_rom['Mpeak'],objs_rom['Mstar_Mpeak'],c = (objs_rom['massiveDist']),cmap = cmx, norm = cNorm,alpha = 0.4, s = markersize*0.5)
    ax2.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],objs_pd_comb['Mstar_Mpeak'][objs_pd_comb['type']=='Central'],c = (objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax2.scatter(objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb['Mstar_Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = (objs_pd_comb['massiveDist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], linewidths = edgewidth)     
    ax2.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],objs_pd_comb['Mstar_Mpeak'][objs_pd_comb['type']=='Satellite'],c = (objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k', marker = "*", s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    #ax2.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_ylabel(r'M$_{\mathrm{*, peak}}$/M$_\odot$')
    ax2.set_xlabel(r'M$_{\mathrm{vir, peak}}$/M$_\odot$')
    ax2.axis([1e8, 2e11, 2e2, 5e9]) #ax2.axis([1e8, 2e13, 2e2, 5e9])
    legend1 = ax2.legend([cen_plt,sat_plt],['Central','Satellite'],scatterpoints = 1,facecolor = 'white',loc = 0,framealpha = 0,frameon = False)
    legend2 = ax2.legend([read,nadler,jethwa],['Read+ 2017','Nadler+ 2020','Jethwa+ 2018'],facecolor = 'white',loc = 4,framealpha = 0,frameon = False)
    ax2.add_artist(legend1)
    ax2.add_artist(legend2)      
    cb = mpl.colorbar.ColorbarBase(ax2sub, cmap=cmx, norm=cNorm)
    cb.set_label(r"D$_{\mathrm{massive}}$ [kpc]")
    fig2.tight_layout()
    fig2.show() 
    fig2.savefig(outfile_base + '_SMHM_rMassGal_Mpeak.png',dpi = dpi)

    plt.close('all')
    fig2 = plt.figure(2,figsize=(plt_width,plt_width*aspect_ratio))
    fig2.set_size_inches(plt_width,plt_width*aspect_ratio)
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax2 = fig2.add_subplot(gs[0])
    ax2sub = fig2.add_subplot(gs[1])
    cmx = plt.get_cmap("viridis")
    cNorm  = colors.LogNorm(vmin=10**1.5, vmax = 10**4)
    sm = plt.cm.ScalarMappable(cmap=cmx, norm=cNorm)
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_ylabel(r'M$_{*, peak}$/M$_\odot$')
    ax2.set_xlabel(r'M$_{\mathrm{vir, peak}}$/M$_\odot$')
    #ax2.axis([2e9, 2e11, 3e3, 1e10])
    ax2.axis([1e8, 2e11, 2e2, 5e9])
    ax2.fill_between(10**(mvir_bins[mstar_bins0[2,:]>1e5]),mstar_bins0[0,mstar_bins0[2,:]>1e5],mstar_bins0[4,mstar_bins0[2,:]>1e5],color = sm.to_rgba(math.log10(dist_binsc[0])),alpha = 0.3,linewidth=0)
    ax2.fill_between(10**(mvir_bins[mstar_bins0[2,:]>1e5]),mstar_bins0[1,mstar_bins0[2,:]>1e5],mstar_bins0[3,mstar_bins0[2,:]>1e5],color = sm.to_rgba(math.log10(dist_binsc[0])),alpha = 0.3,linewidth=0)    
    ax2.plot(10**(mvir_bins[mstar_bins0[2,:]>1e5]),mstar_bins0[2,mstar_bins0[2,:]>1e5],c = sm.to_rgba(math.log10(dist_binsc[0])))

    ax2.fill_between(10**(mvir_bins[mstar_bins1[2,:]>1e5]),mstar_bins1[0,mstar_bins1[2,:]>1e5],mstar_bins1[4,mstar_bins1[2,:]>1e5],color = sm.to_rgba(math.log10(dist_binsc[1])),alpha = 0.3,linewidth=0)
    ax2.fill_between(10**(mvir_bins[mstar_bins1[2,:]>1e5]),mstar_bins1[1,mstar_bins1[2,:]>1e5],mstar_bins1[3,mstar_bins1[2,:]>1e5],color = sm.to_rgba(math.log10(dist_binsc[1])),alpha = 0.3,linewidth=0)    
    ax2.plot(10**(mvir_bins[mstar_bins1[2,:]>1e5]),mstar_bins1[2,mstar_bins1[2,:]>1e5],c = sm.to_rgba(math.log10(dist_binsc[1])))

    ax2.fill_between(10**(mvir_bins[mstar_bins2[2,:]>1e5]),mstar_bins2[0,mstar_bins2[2,:]>1e5],mstar_bins2[4,mstar_bins2[2,:]>1e5],color = sm.to_rgba(math.log10(dist_binsc[2])),alpha = 0.3,linewidth=0)
    ax2.fill_between(10**(mvir_bins[mstar_bins2[2,:]>1e5]),mstar_bins2[1,mstar_bins2[2,:]>1e5],mstar_bins2[3,mstar_bins2[2,:]>1e5],color = sm.to_rgba(math.log10(dist_binsc[2])),alpha = 0.3,linewidth=0)    
    ax2.plot(10**(mvir_bins[mstar_bins2[2,:]>1e5]),mstar_bins2[2,mstar_bins2[2,:]>1e5],c = sm.to_rgba(math.log10(dist_binsc[2])))    
    
    ax2.fill_between(10**(mvir_bins[mstar_bins3[2,:]>1e5]),mstar_bins3[0,mstar_bins3[2,:]>1e5],mstar_bins3[4,mstar_bins3[2,:]>1e5],color = sm.to_rgba(math.log10(dist_binsc[3])),alpha = 0.3,linewidth=0)
    ax2.fill_between(10**(mvir_bins[mstar_bins3[2,:]>1e5]),mstar_bins3[1,mstar_bins3[2,:]>1e5],mstar_bins3[3,mstar_bins3[2,:]>1e5],color = sm.to_rgba(math.log10(dist_binsc[3])),alpha = 0.3,linewidth=0)    
    ax2.plot(10**(mvir_bins[mstar_bins3[2,:]>1e5]),mstar_bins3[2,mstar_bins3[2,:]>1e5],c = sm.to_rgba(math.log10(dist_binsc[3])))
    ax2.scatter(objs_pd_comb['Mpeak'],objs_pd_comb[mstar_key],c = (objs_pd_comb['massiveDist']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'], linewidths = edgewidth)
    cb = mpl.colorbar.ColorbarBase(ax2sub, cmap=cmx, norm=cNorm)
    cb.set_label(r"D$_{\mathrm{massive}}$ [kpc]")
    fig2.tight_layout()    
    fig2.show()
    fig2.savefig(outfile_base + '_SMHM_rMassGal_Mpeak_rom.png',dpi = dpi)

    plt.close('all')
    fig2 = plt.figure(2,figsize=(plt_width*2,plt_width*aspect_ratio))
    gs = fig2.add_gridspec(1,2,wspace=0)
    axs2 = gs.subplots(sharey = True) #sharey = True) #, constrained_layout=True) 
    axs2 = axs2.flatten()
    ax2a = axs2[0]
    ax2b = axs2[1]
    cNorm  = colors.LogNorm(vmin=10**1.5, vmax = 10**4)
    cNormLog  = colors.Normalize(vmin=1.5, vmax = 4)
    sm = plt.cm.ScalarMappable(cmap=cmx, norm=cNormLog)
    ax2a.set_xscale('log')
    ax2a.set_yscale('log')
    ax2a.set_ylabel(r'M$_{*, z = 0}$/M$_\odot$')
    ax2a.set_xlabel(r'M$_{\mathrm{vir, z = 0}}$/M$_\odot$')
    #ax2a.axis([2e9, 2e11, 3e3, 1e10])
    ax2a.axis([1e6, 2e11, 1e3, 5e9])
    ax2a.fill_between(10**(mvirz0_bins[mstarz0_bins0[2,:] > 1e5]),mstarz0_bins0[0,mstarz0_bins0[2,:] > 1e5],mstarz0_bins0[4,mstarz0_bins0[2,:] > 1e5],color = sm.to_rgba(math.log10(dist_binsc[0])),alpha = 0.3,linewidth=0)
    ax2a.fill_between(10**(mvirz0_bins[mstarz0_bins0[2,:] > 1e5]),mstarz0_bins0[1,mstarz0_bins0[2,:] > 1e5],mstarz0_bins0[3,mstarz0_bins0[2,:] > 1e5],color = sm.to_rgba(math.log10(dist_binsc[0])),alpha = 0.3,linewidth=0)    
    ax2a.plot(10**(mvirz0_bins[mstarz0_bins0[2,:] > 1e5]),mstarz0_bins0[2,mstarz0_bins0[2,:] > 1e5],c = sm.to_rgba(math.log10(dist_binsc[0])))

    ax2a.fill_between(10**(mvirz0_bins[mstarz0_bins1[2,:] > 1e5]),mstarz0_bins1[0,mstarz0_bins1[2,:] > 1e5],mstarz0_bins1[4,mstarz0_bins1[2,:] > 1e5],color = sm.to_rgba(math.log10(dist_binsc[1])),alpha = 0.3,linewidth=0)
    ax2a.fill_between(10**(mvirz0_bins[mstarz0_bins1[2,:] > 1e5]),mstarz0_bins1[1,mstarz0_bins1[2,:] > 1e5],mstarz0_bins1[3,mstarz0_bins1[2,:] > 1e5],color = sm.to_rgba(math.log10(dist_binsc[1])),alpha = 0.3,linewidth=0)    
    ax2a.plot(10**(mvirz0_bins[mstarz0_bins1[2,:] > 1e5]),mstarz0_bins1[2,mstarz0_bins1[2,:] > 1e5],c = sm.to_rgba(math.log10(dist_binsc[1])))

    ax2a.fill_between(10**(mvirz0_bins[mstarz0_bins2[2,:] > 1e5]),mstarz0_bins2[0,mstarz0_bins2[2,:] > 1e5],mstarz0_bins2[4,mstarz0_bins2[2,:] > 1e5],color = sm.to_rgba(math.log10(dist_binsc[2])),alpha = 0.3,linewidth=0)
    ax2a.fill_between(10**(mvirz0_bins[mstarz0_bins2[2,:] > 1e5]),mstarz0_bins2[1,mstarz0_bins2[2,:] > 1e5],mstarz0_bins2[3,mstarz0_bins2[2,:] > 1e5],color = sm.to_rgba(math.log10(dist_binsc[2])),alpha = 0.3,linewidth=0)    
    ax2a.plot(10**(mvirz0_bins[mstarz0_bins2[2,:] > 1e5]),mstarz0_bins2[2,mstarz0_bins2[2,:] > 1e5],c = sm.to_rgba(math.log10(dist_binsc[2])))    
    
    ax2a.fill_between(10**(mvirz0_bins[mstarz0_bins3[2,:] > 1e5]),mstarz0_bins3[0,mstarz0_bins3[2,:] > 1e5],mstarz0_bins3[4,mstarz0_bins3[2,:] > 1e5],color = sm.to_rgba(math.log10(dist_binsc[3])),alpha = 0.3,linewidth=0)
    ax2a.fill_between(10**(mvirz0_bins[mstarz0_bins3[2,:] > 1e5]),mstarz0_bins3[1,mstarz0_bins3[2,:] > 1e5],mstarz0_bins3[3,mstarz0_bins3[2,:] > 1e5],color = sm.to_rgba(math.log10(dist_binsc[3])),alpha = 0.3,linewidth=0)    
    ax2a.plot(10**(mvirz0_bins[mstarz0_bins3[2,:] > 1e5]),mstarz0_bins3[2,mstarz0_bins3[2,:] > 1e5],c = sm.to_rgba(math.log10(dist_binsc[3])))
    cen_plt = ax2a.scatter(objs_pd_comb['Mhalo_z0'][objs_pd_comb['type']=='Central'],objs_pd_comb['Mstar_z0'][objs_pd_comb['type']=='Central'],c = (objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm, edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax2a.scatter(objs_pd_comb['Mhalo_z0'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb['Mstar_z0'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = (objs_pd_comb['massiveDist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]), cmap = cmx, norm = cNorm, edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]*1.5, marker = '*', linewidths = edgewidth)    
    sat_plt = ax2a.scatter(objs_pd_comb['Mhalo_z0'][objs_pd_comb['type']=='Satellite'],objs_pd_comb['Mstar_z0'][objs_pd_comb['type']=='Satellite'],c = (objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm, edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, marker = '*', linewidths = edgewidth)
    
    ax2b.set_xscale('log')
    ax2b.set_yscale('log')
    ax2b.set_ylabel(r'M$_{\mathrm{*, peak}}$/M$_\odot$')
    ax2b.set_xlabel(r'M$_{\mathrm{vir, peak}}$/M$_\odot$')
    #ax2b.axis([2e9, 2e11, 3e3, 1e10])
    ax2b.axis([1e8, 2e11, 1e3, 5e9])
    ax2b.fill_between(10**(mvir_bins[mstar_bins0[2,:] > 1e5]),mstar_bins0[0,mstar_bins0[2,:] > 1e5],mstar_bins0[4,mstar_bins0[2,:] > 1e5],color = sm.to_rgba(math.log10(dist_binsc[0])),alpha = 0.3,linewidth=0)
    ax2b.fill_between(10**(mvir_bins[mstar_bins0[2,:] > 1e5]),mstar_bins0[1,mstar_bins0[2,:] > 1e5],mstar_bins0[3,mstar_bins0[2,:] > 1e5],color = sm.to_rgba(math.log10(dist_binsc[0])),alpha = 0.3,linewidth=0)    
    ax2b.plot(10**(mvir_bins[mstar_bins0[2,:] > 1e5]),mstar_bins0[2,mstar_bins0[2,:] > 1e5],c = sm.to_rgba(math.log10(dist_binsc[0])))

    ax2b.fill_between(10**(mvir_bins[mstar_bins1[2,:] > 1e5]),mstar_bins1[0,mstar_bins1[2,:] > 1e5],mstar_bins1[4,mstar_bins1[2,:] > 1e5],color = sm.to_rgba(math.log10(dist_binsc[1])),alpha = 0.3,linewidth=0)
    ax2b.fill_between(10**(mvir_bins[mstar_bins1[2,:] > 1e5]),mstar_bins1[1,mstar_bins1[2,:] > 1e5],mstar_bins1[3,mstar_bins1[2,:] > 1e5],color = sm.to_rgba(math.log10(dist_binsc[1])),alpha = 0.3,linewidth=0)    
    ax2b.plot(10**(mvir_bins[mstar_bins1[2,:] > 1e5]),mstar_bins1[2,mstar_bins1[2,:] > 1e5],c = sm.to_rgba(math.log10(dist_binsc[1])))

    ax2b.fill_between(10**(mvir_bins[mstar_bins2[2,:] > 1e5]),mstar_bins2[0,mstar_bins2[2,:] > 1e5],mstar_bins2[4,mstar_bins2[2,:] > 1e5],color = sm.to_rgba(math.log10(dist_binsc[2])),alpha = 0.3,linewidth=0)
    ax2b.fill_between(10**(mvir_bins[mstar_bins2[2,:] > 1e5]),mstar_bins2[1,mstar_bins2[2,:] > 1e5],mstar_bins2[3,mstar_bins2[2,:] > 1e5],color = sm.to_rgba(math.log10(dist_binsc[2])),alpha = 0.3,linewidth=0)    
    ax2b.plot(10**(mvir_bins[mstar_bins2[2,:] > 1e5]),mstar_bins2[2,mstar_bins2[2,:] > 1e5],c = sm.to_rgba(math.log10(dist_binsc[2])))    
    
    ax2b.fill_between(10**(mvir_bins[mstar_bins3[2,:] > 1e5]),mstar_bins3[0,mstar_bins3[2,:] > 1e5],mstar_bins3[4,mstar_bins3[2,:] > 1e5],color = sm.to_rgba(math.log10(dist_binsc[3])),alpha = 0.3,linewidth=0)
    ax2b.fill_between(10**(mvir_bins[mstar_bins3[2,:] > 1e5]),mstar_bins3[1,mstar_bins3[2,:] > 1e5],mstar_bins3[3,mstar_bins3[2,:] > 1e5],color = sm.to_rgba(math.log10(dist_binsc[3])),alpha = 0.3,linewidth=0)    
    ax2b.plot(10**(mvir_bins[mstar_bins3[2,:] > 1e5]),mstar_bins3[2,mstar_bins3[2,:] > 1e5],c = sm.to_rgba(math.log10(dist_binsc[3])))
    ax2b.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Central'],c = (objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm, edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax2b.scatter(objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb[mstar_key][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = (objs_pd_comb['massiveDist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]), cmap = cmx, norm = cNorm, edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]*1.5, marker = '*', linewidths = edgewidth)    
    ax2b.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Satellite'],c = (objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm, edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, marker = '*', linewidths = edgewidth)
    fig2.tight_layout()
    legend2 = ax2a.legend([cen_plt,sat_plt],['Central','Satellite/Backsplash'],scatterpoints = 1,facecolor = 'white',loc = 0,framealpha = 0,frameon = False)
    
    sm = plt.cm.ScalarMappable(cmap=cmx, norm=cNorm)
    cbar = plt.colorbar(sm, location="right",ax=[ax2a,ax2b],aspect = 20) #pad=0.04)
    #cbar = mpl.colorbar.ColorbarBase(ax2sub, cmap=cmx, norm=cNorm, location="right",ax=[ax2a,ax2b]) #pad=0.04)
    cbar.set_label(r'D$_{\mathrm{massive}}$ [kpc]') #, rotation=0)
    fig2.show()
    fig2.set_size_inches(plt_width*2,plt_width*aspect_ratio)
    fig2.show()
    fig2.savefig(outfile_base + '_SMHM2_rMassGal_MpeakMz0_rom.png',dpi = dpi)
    
    plt.close('all')
    #fig2.clear()
    fig2 = plt.figure(2,figsize=(plt_width,plt_width*aspect_ratio))
    fig2.set_size_inches(plt_width,plt_width*aspect_ratio)
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax2 = fig2.add_subplot(gs[0])
    ax2sub = fig2.add_subplot(gs[1])
    cNorm  = colors.LogNorm(vmin=10**1.5, vmax = 10**4)
    ax2.fill_between(read_abunmatch['M200'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_low'][read_abunmatch['M200']> 1e9]/read_abunmatch['M200'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_high'][read_abunmatch['M200']> 1e9]/read_abunmatch['M200'][read_abunmatch['M200']> 1e9],color = 'grey',linewidth=0,alpha = 0.6)
    ax2.fill_between(upper_fid_Mhalo, np.interp( upper_fid_Mhalo, lower_fid_Mhalo, lower_fid_Mstar )/upper_fid_Mhalo, upper_fid_Mstar/upper_fid_Mhalo, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3)#Nadler outer
    ax2.fill_between(upper_fid_inner_Mhalo, np.interp( upper_fid_inner_Mhalo, lower_fid_inner_Mhalo, lower_fid_inner_Mstar)/upper_fid_inner_Mhalo,upper_fid_inner_Mstar/upper_fid_inner_Mhalo, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3) #Nadler inner
    ax2.fill_between(Jethwa_upper_outer_Mhalo, np.interp( Jethwa_upper_outer_Mhalo, Jethwa_lower_outer_Mhalo, Jethwa_lower_outer_Mstar )/Jethwa_upper_outer_Mhalo,Jethwa_upper_outer_Mstar/Jethwa_upper_outer_Mhalo, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3)#Jethwa outer
    ax2.fill_between(Jethwa_upper_inner_Mhalo,np.interp( Jethwa_upper_inner_Mhalo, Jethwa_lower_inner_Mhalo, Jethwa_lower_inner_Mstar )/Jethwa_upper_inner_Mhalo, Jethwa_upper_inner_Mstar/Jethwa_upper_inner_Mhalo, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3) #Jethwa inner
    read, = ax2.plot( [], [],  color="grey", linestyle="solid", linewidth=4) #, label="Read+ 2017" )
    nadler, = ax2.plot( [], [],  color="#33669A", linestyle="solid", linewidth=8) #, label="Nadler+ 2020" )
    jethwa, = ax2.plot( [], [],  color="#66CC66", linestyle="solid", linewidth=8) #, label="Jethwa+ 2018" )        
    ax2.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Central']/objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],c = (objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax2.scatter(objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb[mstar_key][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]/objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = (objs_pd_comb['massiveDist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], linewidths = edgewidth)
    ax2.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Satellite']/objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],c = (objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k', marker = "*", s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    #ax2.scatter(readdata['M200']*1e10,readdata['M_*']*1e7/(readdata['M200']*1e10),marker = '+',c = 'k')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_ylabel(mstar_ratio_axes_label)
    ax2.set_xlabel(r'M$_{\mathrm{vir, peak}}$/M$_\odot$')
    ax2.axis([1e8, 2e11, 1e-6, 0.05])
    legend1 = ax2.legend([cen_plt,sat_plt],['Central','Satellite'],scatterpoints = 1,facecolor = 'white',loc = 4,framealpha = 0,frameon = False)
    legend2 = ax2.legend([read,nadler,jethwa],['Read+ 2017','Nadler+ 2020','Jethwa+ 2018'],facecolor = 'white',loc = 0,framealpha = 0,frameon = False)
    ax2.add_artist(legend1)
    ax2.add_artist(legend2)      
    cb = mpl.colorbar.ColorbarBase(ax2sub, cmap=cmx, norm=cNorm)
    cb.set_label(r"D$_{\mathrm{massive}}$ [kpc]")
    fig2.tight_layout()
    fig2.show() 
    fig2.savefig(outfile_base + '_SMHMr_rMassGal_Mpeak.png',dpi = dpi)

    plt.close('all')
    #fig2.clear()
    #fig2 = plt.figure(2, figsize=(plt_width*2,plt_width*aspect_ratio))
    #gs = gridspec.GridSpec(1,3,width_ratios=[15,15,1])
    #ax2a = fig2.add_subplot(gs[0])
    #ax2b = fig2.add_subplot(gs[1])
    #ax2sub = fig2.add_subplot(gs[2])
    fig2 = plt.figure(figsize=(plt_width,plt_width*aspect_ratio*2))
    gs = fig2.add_gridspec(2, hspace=0)
    axs2 = gs.subplots(sharex = True) #, constrained_layout=True)
    axs2 = axs2.flatten()
    ax2a = axs2[0]
    ax2b = axs2[1]
    cNorm  = colors.LogNorm(vmin=10**1.5, vmax = 10**4)
    ax2a.fill_between(read_abunmatch['M200'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_low'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_high'][read_abunmatch['M200']> 1e9],color = 'grey',linewidth=0,alpha = 0.6)
    ax2a.fill_between(upper_fid_Mhalo, np.interp( upper_fid_Mhalo, lower_fid_Mhalo, lower_fid_Mstar ), upper_fid_Mstar, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3)#Nadler outer
    ax2a.fill_between(upper_fid_inner_Mhalo, np.interp( upper_fid_inner_Mhalo, lower_fid_inner_Mhalo, lower_fid_inner_Mstar ),upper_fid_inner_Mstar, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3) #Nadler inner
    ax2a.fill_between(Jethwa_upper_outer_Mhalo, np.interp( Jethwa_upper_outer_Mhalo, Jethwa_lower_outer_Mhalo, Jethwa_lower_outer_Mstar ),Jethwa_upper_outer_Mstar, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3)#Jethwa outer
    ax2a.fill_between(Jethwa_upper_inner_Mhalo,np.interp( Jethwa_upper_inner_Mhalo, Jethwa_lower_inner_Mhalo, Jethwa_lower_inner_Mstar ), Jethwa_upper_inner_Mstar, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3) #Jethwa inner
    read, = ax2a.plot( [], [],  color="grey", linestyle="solid", linewidth=4) #, label="Read+ 2017" )
    nadler, = ax2a.plot( [], [],  color="#33669A", linestyle="solid", linewidth=8) #, label="Nadler+ 2020" )
    jethwa, = ax2a.plot( [], [],  color="#66CC66", linestyle="solid", linewidth=8) #, label="Jethwa+ 2018" )     
    cen_plt = ax2a.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Central'],c = (objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax2a.scatter(objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb[mstar_key][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = (objs_pd_comb['massiveDist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], linewidths = edgewidth)    
    sat_plt =  ax2a.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Satellite'],c = (objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k', marker = "*", s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    #ax2.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax2a.set_xscale('log')
    ax2a.set_yscale('log')
    ax2a.set_ylabel(mstar_axes_label)
    #ax2a.set_xlabel(r'M$_{\mathrm{vir, peak}}$/M$_\odot$')
    ax2a.axis([1e8, 2e11, 2e2, 5e9])
    #legend1 = ax2a.legend([cen_plt,sat_plt],['Central','Satellite'],scatterpoints = 1,facecolor = 'white',loc = 0,framealpha = 0,frameon = False)
    legend2 = ax2a.legend([read,nadler,jethwa],['Read+ 2017','Nadler+ 2020','Jethwa+ 2018'],facecolor = 'white',loc = 4,framealpha = 0,frameon = False)
    #ax2.add_artist(legend1)
    ax2a.add_artist(legend2)
    ax2a.label_outer()

    cNorm  = colors.LogNorm(vmin=10**1.5, vmax = 10**4)
    ax2b.fill_between(read_abunmatch['M200'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_low'][read_abunmatch['M200']> 1e9]/read_abunmatch['M200'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_high'][read_abunmatch['M200']> 1e9]/read_abunmatch['M200'][read_abunmatch['M200']> 1e9],color = 'grey',linewidth=0,alpha = 0.6)
    ax2b.fill_between(upper_fid_Mhalo, np.interp( upper_fid_Mhalo, lower_fid_Mhalo, lower_fid_Mstar )/upper_fid_Mhalo, upper_fid_Mstar/upper_fid_Mhalo, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3)#Nadler outer
    ax2b.fill_between(upper_fid_inner_Mhalo, np.interp( upper_fid_inner_Mhalo, lower_fid_inner_Mhalo, lower_fid_inner_Mstar)/upper_fid_inner_Mhalo,upper_fid_inner_Mstar/upper_fid_inner_Mhalo, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3) #Nadler inner
    ax2b.fill_between(Jethwa_upper_outer_Mhalo, np.interp( Jethwa_upper_outer_Mhalo, Jethwa_lower_outer_Mhalo, Jethwa_lower_outer_Mstar )/Jethwa_upper_outer_Mhalo,Jethwa_upper_outer_Mstar/Jethwa_upper_outer_Mhalo, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3)#Jethwa outer
    ax2b.fill_between(Jethwa_upper_inner_Mhalo,np.interp( Jethwa_upper_inner_Mhalo, Jethwa_lower_inner_Mhalo, Jethwa_lower_inner_Mstar )/Jethwa_upper_inner_Mhalo, Jethwa_upper_inner_Mstar/Jethwa_upper_inner_Mhalo, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3) #Jethwa inner
    read, = ax2b.plot( [], [],  color="grey", linestyle="solid", linewidth=4) #, label="Read+ 2017" )
    nadler, = ax2b.plot( [], [],  color="#33669A", linestyle="solid", linewidth=8) #, label="Nadler+ 2020" )
    jethwa, = ax2b.plot( [], [],  color="#66CC66", linestyle="solid", linewidth=8) #, label="Jethwa+ 2018" )        
    ax2b.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Central']/objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],c = (objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax2b.scatter(objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb[mstar_key][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]/objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = (objs_pd_comb['massiveDist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], linewidths = edgewidth)    
    ax2b.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Satellite']/objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],c = (objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k', marker = "*", s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    #ax2.scatter(readdata['M200']*1e10,readdata['M_*']*1e7/(readdata['M200']*1e10),marker = '+',c = 'k')
    ax2b.set_xscale('log')
    ax2b.set_yscale('log')
    ax2b.set_ylabel(mstar_ratio_axes_label)
    ax2b.set_xlabel(r'M$_{\mathrm{vir, peak}}$/M$_\odot$')
    ax2b.axis([1e8, 2e11, 1e-6, 0.05])
    legend1 = ax2b.legend([cen_plt,sat_plt],['Central','Satellite'],scatterpoints = 1,facecolor = 'white',loc = 4,framealpha = 0,frameon = False)
    #legend2 = ax2b.legend([read,nadler,jethwa],['Read+ 2017','Nadler+ 2020','Jethwa+ 2018'],facecolor = 'white',loc = 0,framealpha = 0,frameon = False)
    ax2b.add_artist(legend1)
    #ax2b.add_artist(legend2)
    ax2b.label_outer()
    fig2.tight_layout()
    
    sm = plt.cm.ScalarMappable(cmap=cmx, norm=cNorm)
    cbar = plt.colorbar(sm, location="right",ax=[ax2a,ax2b],aspect = 40) #pad=0.04)
    #cbar = mpl.colorbar.ColorbarBase(ax2sub, cmap=cmx, norm=cNorm, location="right",ax=[ax2a,ax2b]) #pad=0.04)
    cbar.set_label(r'D$_{\mathrm{massive}}$ [kpc]') #, rotation=0)
    #fig2.tight_layout()

    #sm = plt.cm.ScalarMappable(cmap=cmx, norm=cNorm)
    #sm.set_array([])
    #cbar = plt.colorbar(sm, location="top", ax=axs2.ravel().tolist(), fraction = 0.1) #pad=0.04)
    #cbar.set_label(r'D$_{\mathrm{massive}}$ [kpc]') #, rotation=0)
    #cbar.ax.get_xaxis().labelpad = 15
    fig2.show()
    fig2.savefig(outfile_base + '_SMHM2_rMassGal_Mpeak.png',dpi = dpi)
    
    #SMHM colored by HI mass
    #plt.clf()
    plt.close('all')
    fig3 = plt.figure(3,figsize=(plt_width,plt_width*aspect_ratio))
    fig3.clear()
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax3 = fig3.add_subplot(gs[0])
    ax3sub = fig3.add_subplot(gs[1])    
    cNorm  = colors.Normalize(vmin=2.0, vmax = 10.0)
    ax3.fill_between(read_abunmatch['M200'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_low'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_high'][read_abunmatch['M200']> 1e9],color = 'grey',linewidth=0,alpha = 0.6)
    ax3.fill_between(upper_fid_Mhalo, np.interp( upper_fid_Mhalo, lower_fid_Mhalo, lower_fid_Mstar ), upper_fid_Mstar, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3)#Nadler outer
    ax3.fill_between(upper_fid_inner_Mhalo, np.interp( upper_fid_inner_Mhalo, lower_fid_inner_Mhalo, lower_fid_inner_Mstar ),upper_fid_inner_Mstar, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3) #Nadler inner
    ax3.fill_between(Jethwa_upper_outer_Mhalo, np.interp( Jethwa_upper_outer_Mhalo, Jethwa_lower_outer_Mhalo, Jethwa_lower_outer_Mstar ),Jethwa_upper_outer_Mstar, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3)#Jethwa outer
    ax3.fill_between(Jethwa_upper_inner_Mhalo,np.interp( Jethwa_upper_inner_Mhalo, Jethwa_lower_inner_Mhalo, Jethwa_lower_inner_Mstar ), Jethwa_upper_inner_Mstar, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3) #Jethwa inner
    read, = ax3.plot( [], [],  color="grey", linestyle="solid", linewidth=4) #, label="Read+ 2017" )
    nadler, = ax3.plot( [], [],  color="#33669A", linestyle="solid", linewidth=8) #, label="Nadler+ 2020" )
    jethwa, = ax3.plot( [], [],  color="#66CC66", linestyle="solid", linewidth=8) #, label="Jethwa+ 2018" )
    ax3.scatter(objs_pd_comb['mass'][objs_pd_comb['type']=='Central'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Central'] ,edgecolor = 'k',facecolor = 'none', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)    
    ax3.scatter(objs_pd_comb['mass'][objs_pd_comb['type']=='Central'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Central'],c = (np.array(objs_pd_comb['mHI'][objs_pd_comb['type']=='Central'].tolist())), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax3.scatter(objs_pd_comb['mass'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb[mstar_key][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')] ,edgecolor = 'k',facecolor = 'none', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], linewidths = edgewidth)    
    ax3.scatter(objs_pd_comb['mass'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb[mstar_key][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = (np.array(objs_pd_comb['mHI'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')].tolist())), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], linewidths = edgewidth)      
    ax3.scatter(objs_pd_comb['mass'][objs_pd_comb['type']=='Satellite'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Satellite'] ,edgecolor = 'k',facecolor = 'none',marker = '*', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)    
    ax3.scatter(objs_pd_comb['mass'][objs_pd_comb['type']=='Satellite'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Satellite'],c = (np.array(objs_pd_comb['mHI'][objs_pd_comb['type']=='Satellite'].tolist())), cmap = cmx, norm = cNorm,edgecolor = 'k',marker = '*', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)  
    #ax3.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax3.set_xscale('log')
    ax3.set_yscale('log')
    ax3.set_ylabel(r'M$_*$/M$_\odot$')
    ax3.set_xlabel(r'M$_{\mathrm{vir}}$/M$_\odot$')
    ax3.axis([2e6, 2e11, 2e2, 5e9])
    legend1 = ax3.legend([cen_plt,sat_plt],['Central','Satellite'],scatterpoints = 1,facecolor = 'white',loc = 4,framealpha = 0,frameon = False)
    legend2 = ax3.legend([read,nadler,jethwa],['Read+ 2017','Nadler+ 2020','Jethwa+ 2018'],facecolor = 'white',loc = 0,framealpha = 0,frameon = False)
    ax3.add_artist(legend1)
    ax3.add_artist(legend2)
    fig3.show()
    cb = mpl.colorbar.ColorbarBase(ax3sub, cmap=cmx, norm=cNorm)
    cb.set_label(r"M$_{\mathrm{HI}}$ [M$_\odot$]")
    fig3.set_size_inches(plt_width,plt_width*aspect_ratio)
    fig3.tight_layout()
    fig3.show()
    fig3.savefig(outfile_base + '_SMHM_mHI.png',dpi = dpi)
    
    #plt.clf()
    plt.close('all')
    fig3 = plt.figure(3,figsize=(plt_width,plt_width*aspect_ratio))
    fig3.clear()
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax3 = fig3.add_subplot(gs[0])
    ax3sub = fig3.add_subplot(gs[1])    
    cNorm  = colors.Normalize(vmin=2, vmax = 10)
    ax3.fill_between(read_abunmatch['M200'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_low'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_high'][read_abunmatch['M200']> 1e9],color = 'grey',linewidth=0,alpha = 0.6)
    ax3.fill_between(upper_fid_Mhalo, np.interp( upper_fid_Mhalo, lower_fid_Mhalo, lower_fid_Mstar ), upper_fid_Mstar, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3)#Nadler outer
    ax3.fill_between(upper_fid_inner_Mhalo, np.interp( upper_fid_inner_Mhalo, lower_fid_inner_Mhalo, lower_fid_inner_Mstar ),upper_fid_inner_Mstar, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3) #Nadler inner
    ax3.fill_between(Jethwa_upper_outer_Mhalo, np.interp( Jethwa_upper_outer_Mhalo, Jethwa_lower_outer_Mhalo, Jethwa_lower_outer_Mstar ),Jethwa_upper_outer_Mstar, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3)#Jethwa outer
    ax3.fill_between(Jethwa_upper_inner_Mhalo,np.interp( Jethwa_upper_inner_Mhalo, Jethwa_lower_inner_Mhalo, Jethwa_lower_inner_Mstar ), Jethwa_upper_inner_Mstar, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3) #Jethwa inner
    read, = ax3.plot( [], [],  color="grey", linestyle="solid", linewidth=4) #, label="Read+ 2017" )
    nadler, = ax3.plot( [], [],  color="#33669A", linestyle="solid", linewidth=8) #, label="Nadler+ 2020" )
    jethwa, = ax3.plot( [], [],  color="#66CC66", linestyle="solid", linewidth=8) #, label="Jethwa+ 2018" )
    ax3.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Central'] ,edgecolor = 'k',facecolor = 'none', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)    
    ax3.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Central'],c = (np.array(objs_pd_comb['mHI'][objs_pd_comb['type']=='Central'].tolist())), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax3.scatter(objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb[mstar_key][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')] ,edgecolor = 'k',facecolor = 'none', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], linewidths = edgewidth)    
    ax3.scatter(objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb[mstar_key][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = (np.array(objs_pd_comb['mHI'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')].tolist())), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], linewidths = edgewidth)
    ax3.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Satellite'] ,edgecolor = 'k',facecolor = 'none', marker = '*', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)    
    ax3.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Satellite'],c = (np.array(objs_pd_comb['mHI'][objs_pd_comb['type']=='Satellite'].tolist())), cmap = cmx, norm = cNorm,edgecolor = 'k', marker = '*', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    #ax3.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax3.set_xscale('log')
    ax3.set_yscale('log')
    ax3.set_ylabel(r'M$_*$/M$_\odot$')
    ax3.set_xlabel(r'M$_{\mathrm{vir, peak}}$/M$_\odot$')
    ax3.axis([1e8, 2e11, 2e2, 5e9])
    legend1 = ax3.legend([cen_plt,sat_plt],['Central','Satellite'],scatterpoints = 1,facecolor = 'white',loc = 0,framealpha = 0,frameon = False)
    legend2 = ax3.legend([read,nadler,jethwa],['Read+ 2017','Nadler+ 2020','Jethwa+ 2018'],facecolor = 'white',loc = 4,framealpha = 0,frameon = False)
    ax3.add_artist(legend1)
    ax3.add_artist(legend2)
    fig3.show()
    cb = mpl.colorbar.ColorbarBase(ax3sub, cmap=cmx, norm=cNorm)
    cb.set_label(r"M$_{\mathrm{HI}}$ [M$_\odot$]")
    fig3.set_size_inches(plt_width,plt_width*aspect_ratio)
    fig3.tight_layout()
    fig3.show()
    fig3.savefig(outfile_base + '_SMHM_mHI_Mpeak.png',dpi = dpi)
   
    #SMHM colored by HI fraction ###############################
    #plt.clf()
    plt.close('all')
    fig4 = plt.figure(4,figsize=(plt_width,plt_width*aspect_ratio))
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax4 = fig4.add_subplot(gs[0])
    ax4sub = fig4.add_subplot(gs[1])    
    cNorm  = colors.Normalize(vmin=0, vmax = 1)
    ax4.fill_between(read_abunmatch['M200'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_low'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_high'][read_abunmatch['M200']> 1e9],color = 'grey',linewidth=0,alpha = 0.6)
    ax4.fill_between(upper_fid_Mhalo, np.interp( upper_fid_Mhalo, lower_fid_Mhalo, lower_fid_Mstar ), upper_fid_Mstar, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3)#Nadler outer
    ax4.fill_between(upper_fid_inner_Mhalo, np.interp( upper_fid_inner_Mhalo, lower_fid_inner_Mhalo, lower_fid_inner_Mstar ),upper_fid_inner_Mstar, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3) #Nadler inner
    ax4.fill_between(Jethwa_upper_outer_Mhalo, np.interp( Jethwa_upper_outer_Mhalo, Jethwa_lower_outer_Mhalo, Jethwa_lower_outer_Mstar ),Jethwa_upper_outer_Mstar, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3)#Jethwa outer
    ax4.fill_between(Jethwa_upper_inner_Mhalo,np.interp( Jethwa_upper_inner_Mhalo, Jethwa_lower_inner_Mhalo, Jethwa_lower_inner_Mstar ), Jethwa_upper_inner_Mstar, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3) #Jethwa inner
    read, = ax4.plot( [], [],  color="grey", linestyle="solid", linewidth=4) #, label="Read+ 2017" )
    nadler, = ax4.plot( [], [],  color="#33669A", linestyle="solid", linewidth=8) #, label="Nadler+ 2020" )
    jethwa, = ax4.plot( [], [],  color="#66CC66", linestyle="solid", linewidth=8) #, label="Jethwa+ 2018" )
    #Adding 1 to avoid divide by zero error
    ax4.scatter(objs_pd_comb['mass'][objs_pd_comb['type']=='Central'],objs_pd_comb['mstar'][objs_pd_comb['type']=='Central'],c = objs_pd_comb['mHI'][objs_pd_comb['type']=='Central']/(objs_pd_comb['mstar'][objs_pd_comb['type']=='Central'] + objs_pd_comb['mHI'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax4.scatter(objs_pd_comb['mass'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb['mstar'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = objs_pd_comb['mHI'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]/(objs_pd_comb['mstar'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')] + objs_pd_comb['mHI'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], linewidths = edgewidth)
    ax4.scatter(objs_pd_comb['mass'][objs_pd_comb['type']=='Satellite'],objs_pd_comb['mstar'][objs_pd_comb['type']=='Satellite'],c = objs_pd_comb['mHI'][objs_pd_comb['type']=='Satellite']/(objs_pd_comb['mstar'][objs_pd_comb['type']=='Satellite'] + objs_pd_comb['mHI'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k', marker = '*', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    #ax4.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax4.set_xscale('log')
    ax4.set_yscale('log')
    ax4.set_ylabel(r'M$_*$/M$_\odot$')
    ax4.set_xlabel(r'M$_{\mathrm{vir}}$/M$_\odot$')
    ax4.axis([2e6, 2e11, 2e2, 5e9])
    legend1 = ax4.legend([cen_plt,sat_plt],['Central','Satellite'],scatterpoints = 1,facecolor = 'white',loc = 4,framealpha = 0,frameon = False)
    legend2 = ax4.legend([read,nadler,jethwa],['Read+ 2017','Nadler+ 2020','Jethwa+ 2018'],facecolor = 'white',loc = 0,framealpha = 0,frameon = False)
    ax4.add_artist(legend1)
    ax4.add_artist(legend2)
    fig4.show()
    cb = mpl.colorbar.ColorbarBase(ax4sub, cmap=cmx, norm=cNorm)
    cb.set_label(r"M$_{\mathrm{HI}}$/(M$_*$ + M$_{\mathrm{HI}}$)")
    fig4.set_size_inches(plt_width,plt_width*aspect_ratio)
    fig4.tight_layout()
    fig4.show()
    fig4.savefig(outfile_base + '_SMHM_fHI.png',dpi = dpi)

    plt.close('all')
    fig4 = plt.figure(4,figsize=(plt_width,plt_width*aspect_ratio))
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax4 = fig4.add_subplot(gs[0])
    ax4sub = fig4.add_subplot(gs[1])    
    cNorm  = colors.Normalize(vmin=0, vmax = 1)
    ax4.fill_between(read_abunmatch['M200'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_low'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_high'][read_abunmatch['M200']> 1e9],color = 'grey',linewidth=0,alpha = 0.6)
    ax4.fill_between(upper_fid_Mhalo, np.interp( upper_fid_Mhalo, lower_fid_Mhalo, lower_fid_Mstar ), upper_fid_Mstar, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3)#Nadler outer
    ax4.fill_between(upper_fid_inner_Mhalo, np.interp( upper_fid_inner_Mhalo, lower_fid_inner_Mhalo, lower_fid_inner_Mstar ),upper_fid_inner_Mstar, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3) #Nadler inner
    ax4.fill_between(Jethwa_upper_outer_Mhalo, np.interp( Jethwa_upper_outer_Mhalo, Jethwa_lower_outer_Mhalo, Jethwa_lower_outer_Mstar ),Jethwa_upper_outer_Mstar, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3)#Jethwa outer
    ax4.fill_between(Jethwa_upper_inner_Mhalo,np.interp( Jethwa_upper_inner_Mhalo, Jethwa_lower_inner_Mhalo, Jethwa_lower_inner_Mstar ), Jethwa_upper_inner_Mstar, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3) #Jethwa inner
    read, = ax4.plot( [], [],  color="grey", linestyle="solid", linewidth=4) #, label="Read+ 2017" )
    nadler, = ax4.plot( [], [],  color="#33669A", linestyle="solid", linewidth=8) #, label="Nadler+ 2020" )
    jethwa, = ax4.plot( [], [],  color="#66CC66", linestyle="solid", linewidth=8) #, label="Jethwa+ 2018" ) 
    ax4.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],objs_pd_comb['mstar'][objs_pd_comb['type']=='Central'],c = objs_pd_comb['mHI'][objs_pd_comb['type']=='Central']/(objs_pd_comb['mstar'][objs_pd_comb['type']=='Central'] + objs_pd_comb['mHI'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax4.scatter(objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb['mstar'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = objs_pd_comb['mHI'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]/(objs_pd_comb['mstar'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')] + objs_pd_comb['mHI'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], linewidths = edgewidth)    
    ax4.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],objs_pd_comb['mstar'][objs_pd_comb['type']=='Satellite'],c = objs_pd_comb['mHI'][objs_pd_comb['type']=='Satellite']/(objs_pd_comb['mstar'][objs_pd_comb['type']=='Satellite'] + objs_pd_comb['mHI'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k',marker = '*', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    #ax4.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax4.set_xscale('log')
    ax4.set_yscale('log')
    ax4.set_ylabel(r'M$_*$/M$_\odot$')
    ax4.set_xlabel(r'M$_{\mathrm{vir, peak}}$/M$_\odot$')
    ax4.axis([1e8, 2e11, 2e2, 5e9])
    legend1 = ax4.legend([cen_plt,sat_plt],['Central','Satellite'],scatterpoints = 1,facecolor = 'white',loc = 0,framealpha = 0,frameon = False)
    legend2 = ax4.legend([read,nadler,jethwa],['Read+ 2017','Nadler+ 2020','Jethwa+ 2018'],facecolor = 'white',loc = 4,framealpha = 0,frameon = False)
    ax4.add_artist(legend1)
    ax4.add_artist(legend2)
    fig4.show()
    cb = mpl.colorbar.ColorbarBase(ax4sub, cmap=cmx, norm=cNorm)
    cb.set_label(r"M$_{\mathrm{HI}}$/(M$_*$ + M$_{\mathrm{HI}}$)")
    fig4.set_size_inches(plt_width,plt_width*aspect_ratio)
    fig4.tight_layout()
    fig4.show()
    fig4.savefig(outfile_base + '_SMHM_fHI_Mpeak.png',dpi = dpi)   
    
    #Baryonic mass in disk vs. halo mass, colored by distance to massive galaxy ###############################
    plt.close('all')
    fig5 = plt.figure(5,figsize=(plt_width,plt_width*aspect_ratio))
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax5 = fig5.add_subplot(gs[0])
    #ax5sub = fig5.add_subplot(gs[1])
    cNorm  = colors.LogNorm(vmin=10**1.5, vmax = 10**4)
    ax5.scatter(objs_pd_comb['mass'][objs_pd_comb['type']=='Central'],(objs_pd_comb['Mstar_z0'][objs_pd_comb['type']=='Central'] + objs_pd_comb['mHI'][objs_pd_comb['type']=='Central']),c = (objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central']*1, linewidths = edgewidth)
    ax5.scatter(objs_pd_comb['mass'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],(objs_pd_comb['mstar'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')] + objs_pd_comb['mHI'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]),c = (objs_pd_comb['massiveDist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]*1, linewidths = edgewidth)
    ax5.scatter(objs_pd_comb['mass'][objs_pd_comb['type']=='Satellite'],(objs_pd_comb['Mstar_z0'][objs_pd_comb['type']=='Satellite'] + objs_pd_comb['mHI'][objs_pd_comb['type']=='Satellite']),c = (objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k', marker = '*', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    #ax5.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax5.set_xscale('log')
    ax5.set_yscale('log')
    ax5.set_ylabel(r'M$_*$ + M$_{\mathrm{HI}}$ [M$_\odot$]')
    ax5.set_xlabel(r'M$_{\mathrm{vir}}$ [M$_\odot$]')
    ax5.axis([2e6, 2e11, 1e3, 1e10])
    fig5.tight_layout()
    ax5.legend([cen_plt,sat_plt],['Central','Satellite'],scatterpoints = 1,facecolor = 'white',loc = 0,framealpha = 0,frameon = False)
    sm = plt.cm.ScalarMappable(cmap=cmx, norm=cNorm)
    cbar = plt.colorbar(sm, ax=[ax5],aspect = 20) #pad=0.04)
    cbar.set_label(r'D$_{\mathrm{massive}}$ [kpc]') #, rotation=0)
    #cb = mpl.colorbar.ColorbarBase(ax5sub, cmap=cmx, norm=cNorm)
    #cb.set_label(r"D$_{\mathrm{massive}}$ [kpc]")
    fig5.show()
    fig5.set_size_inches(plt_width,plt_width*aspect_ratio)
    fig5.show()
    fig5.savefig(outfile_base + '_BMHM_rMassGal.png',dpi = dpi)

    plt.close('all')
    fig5 = plt.figure(5,figsize=(plt_width,plt_width*aspect_ratio))
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax5 = fig5.add_subplot(gs[0])
    #ax5sub = fig5.add_subplot(gs[1])
    cNorm  = colors.LogNorm(vmin=10**1.5, vmax = 10**4)    
    ax5.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],(objs_pd_comb['Mstar_z0'][objs_pd_comb['type']=='Central'] + objs_pd_comb['mHI'][objs_pd_comb['type']=='Central']),c = (objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax5.scatter(objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],(objs_pd_comb['mstar'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')] + objs_pd_comb['mHI'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]),c = (objs_pd_comb['massiveDist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]*1.5, marker="*", linewidths = edgewidth)    
    ax5.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],(objs_pd_comb['Mstar_z0'][objs_pd_comb['type']=='Satellite'] + objs_pd_comb['mHI'][objs_pd_comb['type']=='Satellite']),c = (objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k', marker="*", s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    #ax5.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax5.set_xscale('log')
    ax5.set_yscale('log')
    ax5.set_ylabel(r'M$_*$ + M$_{\mathrm{HI}}$ [M$_\odot$]')
    ax5.set_xlabel(r'M$_{\mathrm{vir, peak}}$ [M$_\odot$]')
    ax5.axis([1e8, 2e11, 1e3, 1e10])
    fig5.tight_layout()
    ax5.legend([cen_plt,sat_plt],['Central','Satellite/Backsplash'],scatterpoints = 1,facecolor = 'white',loc = 0,framealpha = 0,frameon = False)
    sm = plt.cm.ScalarMappable(cmap=cmx, norm=cNorm)
    cbar = plt.colorbar(sm, ax=[ax5],aspect = 20) #pad=0.04)
    cbar.set_label(r'D$_{\mathrm{massive}}$ [kpc]') #, rotation=0)    
    #cb = mpl.colorbar.ColorbarBase(ax5sub, cmap=cmx, norm=cNorm)
    #cb.set_label(r"D$_{\mathrm{massive}}$ [kpc]")
    fig5.set_size_inches(plt_width,plt_width*aspect_ratio)
    fig5.show()
    fig5.savefig(outfile_base + '_BMHM_rMassGal_Mpeak.png',dpi = dpi)

    plt.close('all')
    f_bary = 0.15637
    fig5 = plt.figure(5,figsize=(plt_width,plt_width*aspect_ratio))
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax5 = fig5.add_subplot(gs[0])
    cmx = plt.get_cmap("viridis") 
    #ax5sub = fig5.add_subplot(gs[1])
    cNorm  = colors.LogNorm(vmin=10, vmax = 2e3)    
    cen_plt = ax5.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],(objs_pd_comb['Mstar_z0'][objs_pd_comb['type']=='Central'] + objs_pd_comb['mHI'][objs_pd_comb['type']=='Central']),c = (objs_pd_comb['min_dist'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    sat_plt = ax5.scatter(objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],(objs_pd_comb['mstar'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')] + objs_pd_comb['mHI'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]),c = (objs_pd_comb['min_dist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]*1.5, marker="*", linewidths = edgewidth)    
    ax5.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],(objs_pd_comb['Mstar_z0'][objs_pd_comb['type']=='Satellite'] + objs_pd_comb['mHI'][objs_pd_comb['type']=='Satellite']),c = (objs_pd_comb['min_dist'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k', marker="*", s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    log_mvir_arr = np.linspace(7,12,num = 500)
    mvir_arr = 10**log_mvir_arr    
    fbary_ln = ax5.plot(mvir_arr, f_bary*mvir_arr, 'k--')
    fbary_1_ln = ax5.plot(mvir_arr, 0.1*f_bary*mvir_arr, 'k')
    fbary_01_ln = ax5.plot(mvir_arr, 0.01*f_bary*mvir_arr, 'k:')
    fbary_001_ln = ax5.plot(mvir_arr, 0.001*f_bary*mvir_arr, 'k-.')
    #ax5.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax5.set_xscale('log')
    ax5.set_yscale('log')
    ax5.set_ylabel(r'M$_*$ + M$_{\mathrm{HI}}$ [M$_\odot$]')
    ax5.set_xlabel(r'M$_{\mathrm{vir, peak}}$ [M$_\odot$]')
    ax5.axis([1e8, 2e11, 1e3, 1e10])
    fig5.tight_layout()
    legend1 = ax5.legend([sat_plt,cen_plt],['Satellite/Backsplash','Central'],scatterpoints = 1,facecolor = 'white',loc = 'upper left',framealpha = 0,frameon = False)
    ax5.legend(['f$_{\mathrm{bary}}$', '10% f$_{\mathrm{bary}}$', '1% f$_{\mathrm{bary}}$', '0.1% f$_{\mathrm{bary}}$'],facecolor = 'white',loc = 0,framealpha = 0,frameon = False)
    ax5.add_artist(legend1)
    #ax5.add_artist(legend2)
    sm = plt.cm.ScalarMappable(cmap=cmx, norm=cNorm)
    cbar = plt.colorbar(sm, ax=[ax5],aspect = 20) #pad=0.04)
    cbar.set_label(r'Min(D$_{\mathrm{massive}}$) [kpc]') #, rotation=0)    
    fig5.set_size_inches(plt_width,plt_width*aspect_ratio)
    fig5.show()
    fig5.savefig(outfile_base + '_BMHM_rClosestGal_Mpeak.png',dpi = dpi)

    # Baryonic Tully Fisher, including 
    plt.close('all')
    f_bary = 0.15637
    fig5 = plt.figure(5,figsize=(plt_width,plt_width*aspect_ratio))
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax5 = fig5.add_subplot(gs[0])
    cmx = plt.get_cmap("viridis") 
    #ax5sub = fig5.add_subplot(gs[1])
    cNorm  = colors.LogNorm(vmin=10, vmax = 2e3)    
    cen_plt = ax5.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],(objs_pd_comb['bmass_peak'][objs_pd_comb['type']=='Central']),c = (objs_pd_comb['min_dist'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    sat_plt = ax5.scatter(objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],(objs_pd_comb['bmass_peak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]),c = (objs_pd_comb['min_dist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]*1.5, marker="*", linewidths = edgewidth)    
    ax5.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],(objs_pd_comb['bmass_peak'][objs_pd_comb['type']=='Satellite']),c = (objs_pd_comb['min_dist'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k', marker="*", s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    log_mvir_arr = np.linspace(7,12,num = 500)
    mvir_arr = 10**log_mvir_arr    
    fbary_ln = ax5.plot(mvir_arr, f_bary*mvir_arr, 'k--')
    fbary_1_ln = ax5.plot(mvir_arr, 0.1*f_bary*mvir_arr, 'k')
    fbary_01_ln = ax5.plot(mvir_arr, 0.01*f_bary*mvir_arr, 'k:')
    fbary_001_ln = ax5.plot(mvir_arr, 0.001*f_bary*mvir_arr, 'k-.')
    #ax5.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax5.set_xscale('log')
    ax5.set_yscale('log')
    ax5.set_ylabel(r'M$_*$ + M$_{\mathrm{gas}}$ [M$_\odot$]')
    ax5.set_xlabel(r'M$_{\mathrm{vir, peak}}$ [M$_\odot$]')
    ax5.axis([1e8, 2e11, 1e3, 1e10])
    fig5.tight_layout()
    legend1 = ax5.legend([sat_plt,cen_plt],['Satellite/Backsplash','Central'],scatterpoints = 1,facecolor = 'white',loc = 'upper left',framealpha = 0,frameon = False)
    ax5.legend(['f$_{\mathrm{bary}}$', '10% f$_{\mathrm{bary}}$', '1% f$_{\mathrm{bary}}$', '0.1% f$_{\mathrm{bary}}$'],facecolor = 'white',loc = 0,framealpha = 0,frameon = False)
    ax5.add_artist(legend1)
    #ax5.add_artist(legend2)
    sm = plt.cm.ScalarMappable(cmap=cmx, norm=cNorm)
    cbar = plt.colorbar(sm, ax=[ax5],aspect = 20) #pad=0.04)
    cbar.set_label(r'Min(D$_{\mathrm{massive}}$) [kpc]') #, rotation=0)    
    fig5.set_size_inches(plt_width,plt_width*aspect_ratio)
    fig5.show()
    fig5.savefig(outfile_base + '_B2MHM_rClosestGal_Mpeak.png',dpi = dpi)
    
    plt.close('all')
    fig5 = plt.figure(5,figsize=(plt_width,plt_width*aspect_ratio))
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax5 = fig5.add_subplot(gs[0])
    #ax5sub = fig5.add_subplot(gs[1])
    cNorm  = colors.LogNorm(vmin=1e2, vmax = 1e8)
    ax5.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],(objs_pd_comb['Mstar_z0'][objs_pd_comb['type']=='Central'] + objs_pd_comb['mHI'][objs_pd_comb['type']=='Central']),c = (objs_pd_comb['max_tindex'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax5.scatter(objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],(objs_pd_comb['mstar'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')] + objs_pd_comb['mHI'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]),c = (objs_pd_comb['max_tindex'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]*1.5, marker="*", linewidths = edgewidth)    
    ax5.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],(objs_pd_comb['Mstar_z0'][objs_pd_comb['type']=='Satellite'] + objs_pd_comb['mHI'][objs_pd_comb['type']=='Satellite']),c = (objs_pd_comb['max_tindex'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k', marker="*", s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    #ax5.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax5.set_xscale('log')
    ax5.set_yscale('log')
    ax5.set_ylabel(r'M$_*$ + M$_{\mathrm{HI}}$ [M$_\odot$]')
    ax5.set_xlabel(r'M$_{\mathrm{vir, peak}}$ [M$_\odot$]')
    ax5.axis([1e8, 2e11, 1e3, 1e10])
    fig5.tight_layout()
    ax5.legend([cen_plt,sat_plt],['Central','Satellite/Backsplash'],scatterpoints = 1,facecolor = 'white',loc = 0,framealpha = 0,frameon = False)
    sm = plt.cm.ScalarMappable(cmap=cmx, norm=cNorm)
    cbar = plt.colorbar(sm, ax=[ax5],aspect = 20) #pad=0.04)
    cbar.set_label(r'Tidal Index') #, rotation=0)    
    #cb = mpl.colorbar.ColorbarBase(ax5sub, cmap=cmx, norm=cNorm)
    #cb.set_label(r"D$_{\mathrm{massive}}$ [kpc]")
    fig5.set_size_inches(plt_width,plt_width*aspect_ratio)
    fig5.show()
    fig5.savefig(outfile_base + '_BMHM_rtindex_Mpeak.png',dpi = dpi)    
    
    plt.close('all')
    fig5 = plt.figure(5,figsize=(plt_width,plt_width*aspect_ratio))
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax5 = fig5.add_subplot(gs[0])
    cNorm  = colors.LogNorm(vmin=10**1.5, vmax = 10**4)
    sm = plt.cm.ScalarMappable(cmap=cmx, norm=cNorm)
    ax5.fill_between(10**mvir_bins,mbar_bins0[0,:],mbar_bins0[4,:],color = sm.to_rgba(dist_binsc[0]),alpha = 0.3,linewidth=0)
    ax5.fill_between(10**mvir_bins,mbar_bins0[1,:],mbar_bins0[3,:],color = sm.to_rgba(dist_binsc[0]),alpha = 0.3,linewidth=0)    
    ax5.plot(10**mvir_bins,mbar_bins0[2,:],c = sm.to_rgba(dist_binsc[0]))

    ax5.fill_between(10**mvir_bins,mbar_bins1[0,:],mbar_bins1[4,:],color = sm.to_rgba(dist_binsc[1]),alpha = 0.3,linewidth=0)
    ax5.fill_between(10**mvir_bins,mbar_bins1[1,:],mbar_bins1[3,:],color = sm.to_rgba(dist_binsc[1]),alpha = 0.3,linewidth=0)    
    ax5.plot(10**mvir_bins,mbar_bins1[2,:],c = sm.to_rgba(dist_binsc[1]))

    ax5.fill_between(10**mvir_bins,mbar_bins2[0,:],mbar_bins2[4,:],color = sm.to_rgba(dist_binsc[2]),alpha = 0.3,linewidth=0)
    ax5.fill_between(10**mvir_bins,mbar_bins2[1,:],mbar_bins2[3,:],color = sm.to_rgba(dist_binsc[2]),alpha = 0.3,linewidth=0)    
    ax5.plot(10**mvir_bins,mbar_bins2[2,:],c = sm.to_rgba(dist_binsc[2]))    
    
    ax5.fill_between(10**mvir_bins,mbar_bins3[0,:],mbar_bins3[4,:],color = sm.to_rgba(dist_binsc[3]),alpha = 0.3,linewidth=0)
    ax5.fill_between(10**mvir_bins,mbar_bins3[1,:],mbar_bins3[3,:],color = sm.to_rgba(dist_binsc[3]),alpha = 0.3,linewidth=0)
    ax5.plot(10**mvir_bins,mbar_bins3[2,:],c = sm.to_rgba(dist_binsc[3]))
    ax5.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],(objs_pd_comb['mstar'][objs_pd_comb['type']=='Central'] + objs_pd_comb['mHI'][objs_pd_comb['type']=='Central']),c = (objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax5.scatter(objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],(objs_pd_comb['mstar'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')] + objs_pd_comb['mHI'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]),c = (objs_pd_comb['massiveDist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], linewidths = edgewidth)    
    ax5.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],(objs_pd_comb['mstar'][objs_pd_comb['type']=='Satellite'] + objs_pd_comb['mHI'][objs_pd_comb['type']=='Satellite']),c = (objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k', marker="*", s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    ax5.set_xscale('log')
    ax5.set_yscale('log')
    ax5.set_ylabel(r'M$_*$ + M$_{\mathrm{HI}}$ [M$_\odot$]')
    ax5.set_xlabel(r'M$_{\mathrm{vir, peak}}$ [M$_\odot$]')
    ax5.axis([1e8, 2e11, 1e3, 1e10]) #ax5.axis([5e9, 2e13, 1e5, 1e12])
    fig5.tight_layout()
    ax5.legend([cen_plt,sat_plt],['Central','Satellite'],scatterpoints = 1,facecolor = 'white',loc = 0,framealpha = 0,frameon = False)
    sm = plt.cm.ScalarMappable(cmap=cmx, norm=cNorm)
    cbar = plt.colorbar(sm, ax=[ax5],aspect = 20) #pad=0.04)
    cbar.set_label(r'D$_{\mathrm{massive}}$ [kpc]') #, rotation=0)    
    fig5.set_size_inches(plt_width,plt_width*aspect_ratio)
    fig5.show()
    fig5.savefig(outfile_base + '_BMHM_rMassGal_Mpeak_rom.png',dpi = dpi)
    
    # Baryon fraction vs distance ###############################
    plt.close('all')
    fig6 = plt.figure(6,figsize=(plt_width,plt_width*aspect_ratio))
    fig6.clear()
    gs = gridspec.GridSpec(1,1)
    ax6 = fig6.add_subplot(gs[0])
    #Baryonic mass fraction of central galaxies
    ax6.scatter(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central'],(objs_pd_comb['M_star'][objs_pd_comb['type']=='Central'] + objs_pd_comb['mHI'][objs_pd_comb['type']=='Central'])/objs_pd_comb['mass'][objs_pd_comb['type']=='Central'],edgecolor = 'k', facecolor = 'none',          s = (objs_pd_comb['rvir'][objs_pd_comb['type']=='Central']*2*ms_scale).tolist(), linewidths = edgewidth)
    ax6.scatter(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central'],(objs_pd_comb['M_star'][objs_pd_comb['type']=='Central'] + objs_pd_comb['mHI'][objs_pd_comb['type']=='Central'])/objs_pd_comb['mass'][objs_pd_comb['type']=='Central'],edgecolor = 'k', facecolor = 'k',alpha = 0.3, s = (objs_pd_comb['rvir'][objs_pd_comb['type']=='Central']*2*ms_scale).tolist(), linewidths = edgewidth)
    ax6.scatter(objs_pd_comb['massiveDist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],(objs_pd_comb['M_star'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')] + objs_pd_comb['mHI'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')])/objs_pd_comb['mass'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],edgecolor = 'k', facecolor = 'none',          s = (objs_pd_comb['rvir'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]*2*ms_scale).tolist(), linewidths = edgewidth)
    ax6.scatter(objs_pd_comb['massiveDist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],(objs_pd_comb['M_star'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')] + objs_pd_comb['mHI'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')])/objs_pd_comb['mass'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],edgecolor = 'k', facecolor = 'k',alpha = 0.3, s = (objs_pd_comb['rvir'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]*2*ms_scale).tolist(), linewidths = edgewidth)
    #Baryonic mass fraction of satellites
    ax6.scatter(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite'],        (objs_pd_comb['M_star'][objs_pd_comb['type']=='Satellite'] + objs_pd_comb['mHI'][objs_pd_comb['type']=='Satellite'])/objs_pd_comb['mass'][objs_pd_comb['type']=='Satellite'],facecolor = 'none',edgecolor = 'k',             marker = '*', s = (objs_pd_comb['rvir'][objs_pd_comb['type']=='Satellite']*2*ms_scale*1.5).tolist(), linewidths = edgewidth)
    fbary = ax6.scatter(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite'],(objs_pd_comb['M_star'][objs_pd_comb['type']=='Satellite'] + objs_pd_comb['mHI'][objs_pd_comb['type']=='Satellite'])/objs_pd_comb['mass'][objs_pd_comb['type']=='Satellite'],facecolor = 'k',   edgecolor = 'k',alpha = 0.3, marker = '*', s = (objs_pd_comb['rvir'][objs_pd_comb['type']=='Satellite']*2*ms_scale*1.5).tolist(), linewidths = edgewidth)
    #Stellar mass fraction of central galaxies
    ax6.scatter(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central'],objs_pd_comb['M_star'][objs_pd_comb['type']=='Central']/objs_pd_comb['mass'][objs_pd_comb['type']=='Central'],    edgecolor = 'red',facecolor = 'none',               s = (objs_pd_comb['rvir'][objs_pd_comb['type']=='Central']*2*ms_scale).tolist(),     linewidths = edgewidth)
    ax6.scatter(objs_pd_comb['massiveDist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb['M_star'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]/objs_pd_comb['mass'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],    edgecolor = 'red',facecolor = 'none',               s = (objs_pd_comb['rvir'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]*2*ms_scale).tolist(),     linewidths = edgewidth)    
    #Stellar mass fraction of satellites
    fstar = ax6.scatter(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite'],objs_pd_comb['M_star'][objs_pd_comb['type']=='Satellite']/objs_pd_comb['mass'][objs_pd_comb['type']=='Satellite'],edgecolor = 'red',facecolor = 'none', marker = '*', s = (objs_pd_comb['rvir'][objs_pd_comb['type']=='Satellite']*2*ms_scale*1.5).tolist(), linewidths = edgewidth)
    #ax6.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    #ax6.scatter(objs_pd['massiveDist'],objs_pd['mHI']/objs_pd['mass'],facecolor = 'none',edgecolor = 'k')
    ax6.set_xscale('log')
    ax6.set_yscale('log')
    ax6.set_ylabel(r'M/M$_{\mathrm{vir}}$')
    ax6.set_xlabel(r"Distance to massive galaxy [kpc]")
    ax6.axis([17, 8e3, 5e-6, 0.2])
    ax6.legend([fstar,fbary],[r'M$_*$/M$_{\mathrm{vir}}$',r'(M$_*$ + M$_{\mathrm{HI}})$/M$_{\mathrm{vir}}$'],loc = 3)
    fig6.show()
    fig6.set_size_inches(plt_width,plt_width*aspect_ratio)
    fig6.tight_layout()
    fig6.show()
    fig6.savefig(outfile_base + '_fbary_rMassGal.png',dpi = dpi)

    plt.close('all')
    fig6 = plt.figure(6,figsize=(plt_width,plt_width*aspect_ratio))
    gs = gridspec.GridSpec(1,1)
    ax6 = fig6.add_subplot(gs[0])
    dmin = 1.4
    dmax = 4
    nbins = 6
    xaxis, data = bin_plt(np.log10(objs_pd_comb['massiveDist']), objs_pd_comb['M_star']/objs_pd_comb['Mpeak'], xmin = dmin, xmax = dmax, nbins = nbins)
    ax6.fill_between(10**xaxis, data[:,0], data[:,4], color='r', alpha=0.1, linewidth=0)
    ax6.fill_between(10**xaxis, data[:,1], data[:,3], color='r', alpha=0.1, linewidth=0)
    ax6.plot(10**xaxis, data[:,2], c='r')
    xaxis, data = bin_plt(np.log10(objs_pd_comb['massiveDist']), (objs_pd_comb['M_star'] + objs_pd_comb['mHI'])/objs_pd_comb['Mpeak'], xmax = dmax, xmin = dmin, nbins = nbins)
    ax6.fill_between(10**xaxis, data[:,0], data[:,4], color='k', alpha=0.1, linewidth=0)
    ax6.fill_between(10**xaxis, data[:,1], data[:,3], color='k', alpha=0.1, linewidth=0)
    ax6.plot(10**xaxis, data[:,2], c='k')           
    #Baryonic mass fraction of central galaxies
    ax6.scatter(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central'],(objs_pd_comb['M_star'][objs_pd_comb['type']=='Central'] + objs_pd_comb['mHI'][objs_pd_comb['type']=='Central'])/objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],edgecolor = 'k', facecolor = 'none',          s = ((objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'])**0.33/10*ms_scale).tolist(), linewidths = edgewidth)
    ax6.scatter(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central'],(objs_pd_comb['M_star'][objs_pd_comb['type']=='Central'] + objs_pd_comb['mHI'][objs_pd_comb['type']=='Central'])/objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],edgecolor = 'k', facecolor = 'k',alpha = 0.3, s = ((objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'])**0.33/10*ms_scale).tolist(), linewidths = edgewidth)
    ax6.scatter(objs_pd_comb['massiveDist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],(objs_pd_comb['M_star'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')] + objs_pd_comb['mHI'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')])/objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],edgecolor = 'k', facecolor = 'none', marker = '*' ,s = ((objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')])**0.33/10*ms_scale*1.5).tolist(), linewidths = edgewidth)
    ax6.scatter(objs_pd_comb['massiveDist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],(objs_pd_comb['M_star'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')] + objs_pd_comb['mHI'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')])/objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],edgecolor = 'k', facecolor = 'k',marker = '*',alpha = 0.3, s = ((objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')])**0.33/10*ms_scale*1.5).tolist(), linewidths = edgewidth)
    #Baryonic mass fraction of satellites
    ax6.scatter(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite'],        (objs_pd_comb['M_star'][objs_pd_comb['type']=='Satellite'] + objs_pd_comb['mHI'][objs_pd_comb['type']=='Satellite'])/objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],facecolor = 'none',edgecolor = 'k',             marker = '*', s = ((objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'])**0.33/10*ms_scale*1.5).tolist(), linewidths = edgewidth)
    fbary = ax6.scatter(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite'],(objs_pd_comb['M_star'][objs_pd_comb['type']=='Satellite'] + objs_pd_comb['mHI'][objs_pd_comb['type']=='Satellite'])/objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],facecolor = 'k',   edgecolor = 'k',alpha = 0.3, marker = '*', s = ((objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'])**0.33/10*ms_scale*1.5).tolist(), linewidths = edgewidth)
    #Stellar mass fraction of central galaxies
    ax6.scatter(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central'],          objs_pd_comb['M_star'][objs_pd_comb['type']=='Central']/objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],    edgecolor = 'red',facecolor = 'none', s = ((objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'])**0.33/10*ms_scale).tolist(),  linewidths = edgewidth)
    
    ax6.scatter(objs_pd_comb['massiveDist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],          objs_pd_comb['M_star'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]/objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],    edgecolor = 'red',facecolor = 'none',  marker = '*', s = ((objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')])**0.33/10*ms_scale).tolist(),     linewidths = edgewidth)
    #Stellar mass fraction of satellites
    fstar = ax6.scatter(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite'],objs_pd_comb['M_star'][objs_pd_comb['type']=='Satellite']/objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],edgecolor = 'red',facecolor = 'none', marker = '*', s = ((objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'])**0.33/10*ms_scale*1.5).tolist(), linewidths = edgewidth)
    #ax6.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax6.set_xscale('log')
    ax6.set_yscale('log')
    ax6.set_ylabel(r'M/M$_{\mathrm{vir, peak}}$')
    ax6.set_xlabel(r"Distance to massive galaxy at $z = 0$ [kpc]")
    ax6.axis([17, 8e3, 1e-6, 0.08])
    ax6.legend([fstar,fbary],[r'M$_*$/M$_{\mathrm{vir, peak}}$',r'(M$_*$ + M$_{\mathrm{HI}})$/M$_{\mathrm{vir, peak}}$'],loc = 3, framealpha = 0)
    fig6.tight_layout()
    fig6.show()
    fig6.set_size_inches(plt_width,plt_width*aspect_ratio)
    fig6.show()
    fig6.savefig(outfile_base + '_fbary_rMassGal_Mpeak.png',dpi = dpi)

    plt.close('all')
    fig6 = plt.figure(6,figsize=(plt_width,plt_width*aspect_ratio))
    gs = gridspec.GridSpec(1,1)
    ax6 = fig6.add_subplot(gs[0])
    dmin = 1.3
    dmax = 3.3
    nbins = 6
    xaxis, data = bin_plt(np.log10(np.array(objs_pd_comb['min_dist'].tolist())), np.array((objs_pd_comb['M_star']/objs_pd_comb['Mpeak']).tolist()), xmin = dmin, xmax = dmax, nbins = nbins)
    ax6.fill_between(10**xaxis, data[:,0], data[:,4], color='r', alpha=0.1, linewidth=0)
    ax6.fill_between(10**xaxis, data[:,1], data[:,3], color='r', alpha=0.1, linewidth=0)
    ax6.plot(10**xaxis, data[:,2], c='r')
    xaxis, data = bin_plt(np.log10(np.array(objs_pd_comb['min_dist'].tolist())), (objs_pd_comb['M_star'] + objs_pd_comb['mHI'])/objs_pd_comb['Mpeak'], xmax = dmax, xmin = dmin, nbins = nbins)
    ax6.fill_between(10**xaxis, data[:,0], data[:,4], color='k', alpha=0.1, linewidth=0)
    ax6.fill_between(10**xaxis, data[:,1], data[:,3], color='k', alpha=0.1, linewidth=0)
    ax6.plot(10**xaxis, data[:,2], c='k')   
    #Baryonic mass fraction of central galaxies
    ax6.scatter(objs_pd_comb['min_dist'][objs_pd_comb['type']=='Central'],(objs_pd_comb['M_star'][objs_pd_comb['type']=='Central'] + objs_pd_comb['mHI'][objs_pd_comb['type']=='Central'])/objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],edgecolor = 'k', facecolor = 'none',          s = ((objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'])**0.33/10*ms_scale).tolist(), linewidths = edgewidth)
    ax6.scatter(objs_pd_comb['min_dist'][objs_pd_comb['type']=='Central'],(objs_pd_comb['M_star'][objs_pd_comb['type']=='Central'] + objs_pd_comb['mHI'][objs_pd_comb['type']=='Central'])/objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],edgecolor = 'k', facecolor = 'k',alpha = 0.3, s = ((objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'])**0.33/10*ms_scale).tolist(), linewidths = edgewidth)
    ax6.scatter(objs_pd_comb['min_dist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],(objs_pd_comb['M_star'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')] + objs_pd_comb['mHI'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')])/objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],edgecolor = 'k', facecolor = 'none', marker = '*' ,s = ((objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')])**0.33/10*ms_scale*1.5).tolist(), linewidths = edgewidth)
    ax6.scatter(objs_pd_comb['min_dist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],(objs_pd_comb['M_star'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')] + objs_pd_comb['mHI'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')])/objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],edgecolor = 'k', facecolor = 'k',marker = '*',alpha = 0.3, s = ((objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')])**0.33/10*ms_scale*1.5).tolist(), linewidths = edgewidth)
    #Baryonic mass fraction of satellites
    ax6.scatter(objs_pd_comb['min_dist'][objs_pd_comb['type']=='Satellite'],        (objs_pd_comb['M_star'][objs_pd_comb['type']=='Satellite'] + objs_pd_comb['mHI'][objs_pd_comb['type']=='Satellite'])/objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],facecolor = 'none',edgecolor = 'k',             marker = '*', s = ((objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'])**0.33/10*ms_scale*1.5).tolist(), linewidths = edgewidth)
    fbary = ax6.scatter(objs_pd_comb['min_dist'][objs_pd_comb['type']=='Satellite'],(objs_pd_comb['M_star'][objs_pd_comb['type']=='Satellite'] + objs_pd_comb['mHI'][objs_pd_comb['type']=='Satellite'])/objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],facecolor = 'k',   edgecolor = 'k',alpha = 0.3, marker = '*', s = ((objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'])**0.33/10*ms_scale*1.5).tolist(), linewidths = edgewidth)
    #Stellar mass fraction of central galaxies
    ax6.scatter(objs_pd_comb['min_dist'][objs_pd_comb['type']=='Central'],          objs_pd_comb['M_star'][objs_pd_comb['type']=='Central']/objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],    edgecolor = 'red',facecolor = 'none', s = ((objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'])**0.33/10*ms_scale).tolist(),  linewidths = edgewidth)
    
    ax6.scatter(objs_pd_comb['min_dist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],          objs_pd_comb['M_star'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]/objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],    edgecolor = 'red',facecolor = 'none',  marker = '*', s = ((objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')])**0.33/10*ms_scale).tolist(),     linewidths = edgewidth)
    #Stellar mass fraction of satellites
    fstar = ax6.scatter(objs_pd_comb['min_dist'][objs_pd_comb['type']=='Satellite'],objs_pd_comb['M_star'][objs_pd_comb['type']=='Satellite']/objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],edgecolor = 'red',facecolor = 'none', marker = '*', s = ((objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'])**0.33/10*ms_scale*1.5).tolist(), linewidths = edgewidth)
    #ax6.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax6.set_xscale('log')
    ax6.set_yscale('log')
    ax6.set_ylabel(r'M/M$_{\mathrm{vir, peak}}$')
    ax6.set_xlabel(r"Minimum distance to massive galaxy [kpc]")    
    # ax6.axis([10, 4e3, 1e-6, 0.2])
    ax6.axis([10, 1.5e3, 1e-6, 0.2])
    ax6.legend([fstar,fbary],[r'M$_*$/M$_{\mathrm{vir, peak}}$',r'(M$_*$ + M$_{\mathrm{HI}})$/M$_{\mathrm{vir, peak}}$'],loc = 3)
    fig6.tight_layout()
    fig6.show()
    fig6.set_size_inches(plt_width,plt_width*aspect_ratio)
    fig6.show()
    fig6.savefig(outfile_base + '_fbary_rClosestDist_Mpeak.png',dpi = dpi)
    
    #Color vs virial or stellar mass
    plt.close('all')
    fig7 = plt.figure(7,figsize=(plt_width,plt_width*aspect_ratio))
    fig7.clear()
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax7 = fig7.add_subplot(gs[0])
    ax7sub = fig7.add_subplot(gs[1]) 
    cNorm  = colors.LogNorm(vmin=10**1.5, vmax = 10**4)    
    ax7.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],(objs_pd_comb['B-V'][objs_pd_comb['type']=='Central']),c = (objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax7.scatter(objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],(objs_pd_comb['B-V'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]),c = (objs_pd_comb['massiveDist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]*1.5, marker = '*', linewidths = edgewidth)
    ax7.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],(objs_pd_comb['B-V'][objs_pd_comb['type']=='Satellite']),c = (objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k', marker = '*', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    ax7.set_xscale('log')
    ax7.set_ylabel(r'B-V')
    ax7.set_xlabel(r'M$_{\mathrm{vir, peak}}$ [M$_\odot$]')
    ax7.axis([1e8, 1e11, 0, 0.8])
    ax7.legend([cen_plt,sat_plt],['Central','Satellite/Backsplash'],scatterpoints = 1,facecolor = 'white',loc = 3,framealpha = 0,frameon = False)
    cb = mpl.colorbar.ColorbarBase(ax7sub, cmap=cmx, norm=cNorm)
    cb.set_label(r"D$_{\mathrm{massive}}$ [kpc]")
    fig7.set_size_inches(plt_width,plt_width*aspect_ratio)
    fig7.tight_layout()
    fig7.show()
    fig7.savefig(outfile_base + '_BV_Mpeak.png',dpi = dpi)

    plt.close('all')
    fig7 = plt.figure(7,figsize=(plt_width,plt_width*aspect_ratio))
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax7 = fig7.add_subplot(gs[0])
    ax7sub = fig7.add_subplot(gs[1])
    cNorm  = colors.LogNorm(vmin=10**1.5, vmax = 10**4)
    ax7.scatter(objs_pd_comb['mass'][objs_pd_comb['type']=='Central'],(objs_pd_comb['B-V'][objs_pd_comb['type']=='Central']),c = (objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax7.scatter(objs_pd_comb['mass'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],(objs_pd_comb['B-V'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]),c = (objs_pd_comb['massiveDist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]*1.5, marker = '*', linewidths = edgewidth)
    ax7.scatter(objs_pd_comb['mass'][objs_pd_comb['type']=='Satellite'],(objs_pd_comb['B-V'][objs_pd_comb['type']=='Satellite']),c = (objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k', marker = '*', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    ax7.set_xscale('log')
    ax7.set_ylabel(r'B-V')
    ax7.set_xlabel(r'M$_{\mathrm{vir, z = 0}}$ [M$_\odot$]')
    ax7.axis([1e8, 1e11, 0, 0.8])
    ax7.legend([cen_plt,sat_plt],['Central','Satellite/Backsplash'],scatterpoints = 1,facecolor = 'white',loc = 3,framealpha = 0,frameon = False)
    cb = mpl.colorbar.ColorbarBase(ax7sub, cmap=cmx, norm=cNorm)
    cb.set_label(r"D$_{\mathrm{massive}}$ [kpc]")
    fig7.set_size_inches(plt_width,plt_width*aspect_ratio)
    fig7.tight_layout()
    fig7.show()
    fig7.savefig(outfile_base + '_BV_Mvir.png',dpi = dpi)

    plt.close('all')
    fig7 = plt.figure(7,figsize=(plt_width,plt_width*aspect_ratio))
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax7 = fig7.add_subplot(gs[0])
    ax7sub = fig7.add_subplot(gs[1]) 
    cNorm  = colors.LogNorm(vmin=10**1.5, vmax = 10**4)
    ax7.scatter(objs_pd_comb[mstar_key][objs_pd_comb['type']=='Central'],(objs_pd_comb['B-V'][objs_pd_comb['type']=='Central']),c = (objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax7.scatter(objs_pd_comb[mstar_key][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],(objs_pd_comb['B-V'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]),c = (objs_pd_comb['massiveDist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]*1.5, marker = '*', linewidths = edgewidth)
    ax7.scatter(objs_pd_comb[mstar_key][objs_pd_comb['type']=='Satellite'],(objs_pd_comb['B-V'][objs_pd_comb['type']=='Satellite']),c = (objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k', marker = '*', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    ax7.set_xscale('log')
    ax7.set_ylabel(r'B-V')
    ax7.set_xlabel(r'M$_{\mathrm{star}}$ [M$_\odot$]')
    ax7.axis([1e3, 5e9, 0.05, 0.8])
    ax7.legend([cen_plt,sat_plt],['Central','Satellite/Backsplash'],scatterpoints = 1,facecolor = 'white',loc = 4,framealpha = 0,frameon = False)
    cb = mpl.colorbar.ColorbarBase(ax7sub, cmap=cmx, norm=cNorm)
    cb.set_label(r"D$_{\mathrm{massive}}$ [kpc]")
    fig7.set_size_inches(plt_width,plt_width*aspect_ratio)
    fig7.tight_layout()
    fig7.show()
    fig7.savefig(outfile_base + '_BV_Mstar.png',dpi = dpi)

    plt.close('all')
    fig7 = plt.figure(7,figsize=(plt_width,plt_width*aspect_ratio))
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax7 = fig7.add_subplot(gs[0])
    ax7sub = fig7.add_subplot(gs[1]) 
    cNorm  = colors.LogNorm(vmin=10**1.5, vmax = 10**4)
    ax7.scatter(objs_pd_comb['mstar'][objs_pd_comb['type']=='Central'],(objs_pd_comb['B-V'][objs_pd_comb['type']=='Central']),c = (objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax7.scatter(objs_pd_comb['mstar'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],(objs_pd_comb['B-V'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]),c = (objs_pd_comb['massiveDist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]*1.5, marker = '*', linewidths = edgewidth)
    ax7.scatter(objs_pd_comb['mstar'][objs_pd_comb['type']=='Satellite'],(objs_pd_comb['B-V'][objs_pd_comb['type']=='Satellite']),c = (objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k', marker = '*', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    #ax5.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax7.set_xscale('log')
    ax7.set_ylabel(r'B-V')
    ax7.set_xlabel(r'M$_{\mathrm{star, peak}}$ [M$_\odot$]')
    ax7.axis([1e3, 5e9, 0.05, 0.8])
    ax7.legend([cen_plt,sat_plt],['Central','Satellite/Backsplash'],scatterpoints = 1,facecolor = 'white',loc = 4,framealpha = 0,frameon = False)
    cb = mpl.colorbar.ColorbarBase(ax7sub, cmap=cmx, norm=cNorm)
    cb.set_label(r"D$_{\mathrm{massive}}$ [kpc]")
    fig7.set_size_inches(plt_width,plt_width*aspect_ratio)
    fig7.tight_layout()
    fig7.show()
    fig7.savefig(outfile_base + '_BV_Mstarpeak.png',dpi = dpi)

    plt.close('all')
    fig7 = plt.figure(7,figsize=(plt_width*2,plt_width*aspect_ratio))
    gs = fig7.add_gridspec(1,2,wspace=0)
    axs7 = gs.subplots(sharey = True) #, constrained_layout=True) 
    axs7 = axs7.flatten()
    ax7a = axs7[0]
    ax7b = axs7[1]
    cNorm  = colors.LogNorm(vmin=10**1.5, vmax = 10**4)
    ax7a.scatter(objs_pd_comb[mstar_key][objs_pd_comb['type']=='Central'],(objs_pd_comb['B-V'][objs_pd_comb['type']=='Central']),c = (objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax7a.scatter(objs_pd_comb[mstar_key][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],(objs_pd_comb['B-V'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]),c = (objs_pd_comb['massiveDist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]*1.5, marker = '*', linewidths = edgewidth)
    ax7a.scatter(objs_pd_comb[mstar_key][objs_pd_comb['type']=='Satellite'],(objs_pd_comb['B-V'][objs_pd_comb['type']=='Satellite']),c = (objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k', marker = '*', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    #ax5.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax7a.set_ylabel(r'B-V')
    ax7a.set_xscale('log')
    ax7a.set_xlabel(r'M$_{\mathrm{star}}$ [M$_\odot$]')
    ax7a.axis([1e3, 5e9, 0.05, 0.85])    
    ax7b.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],(objs_pd_comb['B-V'][objs_pd_comb['type']=='Central']),c = (objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax7b.scatter(objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],(objs_pd_comb['B-V'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]),c = (objs_pd_comb['massiveDist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]*1.5, marker = '*', linewidths = edgewidth)
    ax7b.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],(objs_pd_comb['B-V'][objs_pd_comb['type']=='Satellite']),c = (objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k', marker = '*', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    ax7b.set_xscale('log')
    ax7b.set_xlabel(r'M$_{\mathrm{vir, peak}}$ [M$_\odot$]')
    ax7b.axis([1e8, 1e11, 0.05, 0.85])
    ax7b.legend([cen_plt,sat_plt],['Central','Satellite/Backsplash'],scatterpoints = 1,facecolor = 'white',loc = 3,framealpha = 0,frameon = False)
    fig7.tight_layout()
    sm = plt.cm.ScalarMappable(cmap=cmx, norm=cNorm)
    cbar = plt.colorbar(sm, location="right",ax=[ax7a,ax7b],aspect = 20) #pad=0.04)
    #cbar = mpl.colorbar.ColorbarBase(ax2sub, cmap=cmx, norm=cNorm, location="right",ax=[ax2a,ax2b]) #pad=0.04)
    cbar.set_label(r'D$_{\mathrm{massive}}$ [kpc]') #, rotation=0)
    fig7.show()
    fig7.set_size_inches(plt_width*2,plt_width*aspect_ratio)
    fig7.show()
    fig7.savefig(outfile_base + '_BV_Mvir_Mstar.png',dpi = dpi)

# Observational comparisons        
    plt.close('all')
    fig8 = plt.figure(8,figsize=(plt_width,plt_width*aspect_ratio*2))
    gs = fig8.add_gridspec(2,hspace=0)
    axs8 = gs.subplots(sharex = True) #, constrained_layout=True) 
    axs8 = axs8.flatten()
    ax8b = axs8[0]
    ax8a = axs8[1]
    #cmxr = plt.get_cmap("viridis_r")
    cmxr = plt.get_cmap("winter")
    cNorm  = colors.Normalize(vmin=-20, vmax = -2)
    sm = plt.cm.ScalarMappable(cmap=cmxr, norm=cNorm)
    cen_plt = ax8a.scatter(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central'],(objs_pd_comb['B-V'][objs_pd_comb['type']=='Central']),c = objs_pd_comb['M_V'][objs_pd_comb['type']=='Central'], cmap = cmxr, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    sat_plt = ax8a.scatter(objs_pd_comb['massiveDist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],(objs_pd_comb['B-V'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]),c = objs_pd_comb['M_V'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], cmap = cmxr, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]*1.5, marker = "*", linewidths = edgewidth)
    ax8a.scatter(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite'],(objs_pd_comb['B-V'][objs_pd_comb['type']=='Satellite']),c = objs_pd_comb['M_V'][objs_pd_comb['type']=='Satellite'], cmap = cmxr, norm = cNorm,edgecolor = 'k', marker = '*', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    #ax5.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax8a.set_ylabel(r'B-V')
    ax8a.set_xscale('log')
    ax8a.axis([17, 7e3, 0.09, 0.81])
    ax8a.legend([cen_plt,sat_plt],['Central','Satellite/Backsplash'],scatterpoints = 1,facecolor = 'white',loc = 3,framealpha = 0,frameon = False)
    #ax8a.set_xlabel(r'Distance to massive galaxy (kpc)')

    #ax8b.scatter(MW_M31_dist, MW_M31_MHI_LV, marker='x', c=Calc_V_Mag, vmin=-20, vmax = -8, cmap=cmxr, s = markersize)
    M_V_sun = 4.81
    mHI_LV = mcconnachie_data['M_H I']*1e6/(10**(0.4*(M_V_sun-mcconnachie_data['M_V'])))
    ax8b.scatter(mcconnachie_data['Min(D_MW, D_M31)'], mHI_LV, marker='x', c=mcconnachie_data['M_V'], vmin=-20, vmax = -8, cmap=cmxr, s = markersize)
    mc_nogas = np.where(mcconnachie_data['M_H I'] == 0)[0]
    ax8b.scatter(np.array((mcconnachie_data['Min(D_MW, D_M31)'].to_list()))[mc_nogas],np.zeros(len(mc_nogas)) + 1e-4, marker='x', c=np.array((mcconnachie_data['M_V'].to_list()))[mc_nogas], vmin=-20, vmax = -8, cmap=cmxr, s = markersize)
    colorVal_mc = np.empty(len(mc_nogas))
    i = 0
    for m_v in np.array((mcconnachie_data['M_V'].to_list()))[mc_nogas]:
        colorVal_mc[i] = float(m_v)
        i = i + 1
    colorVal_mc = sm.to_rgba(colorVal_mc)
    for (dist, c) in zip(np.array((mcconnachie_data['Min(D_MW, D_M31)'].to_list()))[mc_nogas], colorVal_mc):
        ax8b.errorbar(dist, 1e-4,yerr =  6e-5, ecolor = c, c=c, uplims = True, marker = '.',zorder = 1, capsize = 3, barsabove = False)
    #ax8b.scatter(MW_M31_dist, MW_M31_MHI_LV, marker=',', c='k') #, s = markersize*0.5)
    
    nogas = (objs_pd_comb['mHI'] == 0)# & np.isnan(objs_pd_comb['Lum'])) #~np.isfinite(objs_pd_comb[objs_pd_comb['mHI'] != 0]['Lum']/objs_pd_comb[objs_pd_comb['mHI'] != 0]['mHI'])
    ax8b.scatter(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central'],(objs_pd_comb['mHI'][objs_pd_comb['type']=='Central']/objs_pd_comb['Lum'][objs_pd_comb['type']=='Central']),c = objs_pd_comb['M_V'][objs_pd_comb['type']=='Central'], cmap = cmxr, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    colorVal = np.empty(len(objs_pd_comb['M_V']))
    i = 0
    for m_v in objs_pd_comb['M_V']:
        colorVal[i] = float(m_v)
        i = i + 1
    colorVal = sm.to_rgba(colorVal)
    for (dist, c, ms) in zip(objs_pd_comb['massiveDist'][(objs_pd_comb['type']=='Central') & nogas], colorVal[(objs_pd_comb['type']=='Central') & nogas], objs_pd_comb['p_size'][(objs_pd_comb['type']=='Central') & nogas]):
        ax8b.errorbar(dist, 1e-4,yerr =  8e-5,ecolor = c, c=c, marker = "o",uplims = True,zorder = 1, capsize = 3, barsabove = False)
        #ax8b.scatter(dist, 1e-4,c=c, marker = "o",edgecolor = 'k')
    cen_plt = ax8b.scatter(objs_pd_comb['massiveDist'][(objs_pd_comb['type']=='Central') & nogas],(objs_pd_comb['mHI'][(objs_pd_comb['type']=='Central') & nogas])*0 + 1e-4,c = objs_pd_comb['M_V'][(objs_pd_comb['type']=='Central') & nogas],cmap = cmxr, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Central') & nogas], linewidths = edgewidth,zorder = 2)
    
    sat_plt = ax8b.scatter(objs_pd_comb['massiveDist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],(objs_pd_comb['mHI'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]/objs_pd_comb['Lum'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]),c = objs_pd_comb['M_V'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], cmap = cmxr, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]*1.5, marker = "*", linewidths = edgewidth)
    for (dist, c) in zip(objs_pd_comb['massiveDist'][((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')) & nogas], colorVal[((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')) & nogas]):
        ax8b.errorbar(dist, 1e-4,yerr =  8e-5,ecolor = c, c=c, marker = "*", uplims = True,zorder = 1, capsize = 3, barsabove = False)
    ax8b.scatter(objs_pd_comb['massiveDist'][((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')) & nogas],(objs_pd_comb['mHI'][((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')) & nogas])*0 + 1e-4,c = objs_pd_comb['M_V'][((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')) & nogas], cmap = cmxr, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')) & nogas], linewidths = edgewidth,zorder = 2)
    
    ax8b.scatter(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite'],(objs_pd_comb['mHI'][objs_pd_comb['type']=='Satellite']/objs_pd_comb['Lum'][objs_pd_comb['type']=='Satellite']),c = objs_pd_comb['M_V'][objs_pd_comb['type']=='Satellite'], cmap = cmxr, norm = cNorm,edgecolor = 'k', marker = '*', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    ax8b.scatter(objs_pd_comb['massiveDist'][(objs_pd_comb['type']=='Satellite') & nogas],(objs_pd_comb['mHI'][(objs_pd_comb['type']=='Satellite') & nogas])*0 + 1e-4,c = objs_pd_comb['M_V'][(objs_pd_comb['type']=='Satellite') & nogas], cmap = cmxr, norm = cNorm,edgecolor = 'k', marker = '*', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Satellite') & nogas]*1.5, linewidths = edgewidth)
    for (dist, c) in zip(objs_pd_comb['massiveDist'][(objs_pd_comb['type']=='Satellite') & nogas], colorVal[(objs_pd_comb['type']=='Satellite') & nogas]):
        ax8b.errorbar(dist, 1e-4,yerr =  8e-5,ecolor = c, c=c, uplims = True, marker = '*',zorder = 1, capsize = 3)        
    ax8b.scatter(objs_pd_comb['massiveDist'][(objs_pd_comb['type']=='Satellite') & nogas],(objs_pd_comb['mHI'][(objs_pd_comb['type']=='Satellite') & nogas])*0 + 1e-4,c = objs_pd_comb['M_V'][(objs_pd_comb['type']=='Satellite') & nogas], cmap = cmxr, norm = cNorm,edgecolor = 'k', marker = '*', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Satellite') & nogas]*1.5, linewidths = edgewidth,zorder = 2)
    
    ax8b.set_xscale('log')
    ax8b.set_yscale('log')
    ax8b.set_xlabel(r'D$_{\mathrm{massive}}$ [kpc]')
    ax8b.set_ylabel(r'M$_{HI}$/L$_V$ [M$_{\odot}$/L$_{\odot}$]')
    ax8b.axis([17, 7e3, 1e-5, 50])
    fig8.tight_layout()
    cbar = plt.colorbar(sm, location="right",ax=[ax8a,ax8b],aspect = 40) #pad=0.04)
    #cbar = mpl.colorbar.ColorbarBase(ax2sub, cmap=cmxr, norm=cNorm, location="right",ax=[ax2a,ax2b]) #pad=0.04)
    cbar.set_label(r'V Magnitude') #, rotation=0)
    cbar.ax.invert_yaxis()
    fig8.show()
    fig8.set_size_inches(plt_width,plt_width*aspect_ratio*2)
    fig8.show()
    fig8.savefig(outfile_base + '_BV_HI_dist.png',dpi = dpi)

# Metallicity comparison for Arora 2021
    plt.close('all')
    fig9 = plt.figure(9,figsize=(plt_width,plt_width*aspect_ratio*2))
    gs = fig9.add_gridspec(2,hspace=0)
    axs9 = gs.subplots(sharex = True) #, constrained_layout=True)     
    axs9 = axs9.flatten()
    ax9a = axs9[0]
    ax9b = axs9[1]
    cNorm  = colors.LogNorm(vmin=10**1.5, vmax = 10**4)
    ax9a.scatter(objs_pd_comb[mstar_key][objs_pd_comb['type']=='Central'],(objs_pd_comb['mZgas'][objs_pd_comb['type']=='Central']),c = (objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax9a.scatter(objs_pd_comb[mstar_key][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],(objs_pd_comb['mZgas'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]),c = (objs_pd_comb['massiveDist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], linewidths = edgewidth)
    ax9a.scatter(objs_pd_comb[mstar_key][objs_pd_comb['type']=='Satellite'],(objs_pd_comb['mZgas'][objs_pd_comb['type']=='Satellite']),c = (objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k', marker = '*', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    #ax5.scatter(readdata['M200']*1e10,readdata['M_*']*1e9,marker = '+',c = 'k')
    ax9a.set_ylabel(r'M$_{\mathrm{Z}}$ [M$_\odot$]')
    ax9a.set_xscale('log')
    ax9a.set_yscale('log')
    ax9a.axis([1e3, 5e9, 1, 1e8])    
    ax9b.scatter(objs_pd_comb[mstar_key][objs_pd_comb['type']=='Central'],(objs_pd_comb['mZgas'][objs_pd_comb['type']=='Central']/objs_pd_comb['mgas'][objs_pd_comb['type']=='Central']),c = (objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax9b.scatter(objs_pd_comb[mstar_key][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],(objs_pd_comb['mZgas'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')])/(objs_pd_comb['mgas'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]),c = (objs_pd_comb['massiveDist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], linewidths = edgewidth)
    ax9b.scatter(objs_pd_comb[mstar_key][objs_pd_comb['type']=='Satellite'],(objs_pd_comb['mZgas'][objs_pd_comb['type']=='Satellite'])/(objs_pd_comb['mgas'][objs_pd_comb['type']=='Satellite']),c = (objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k', marker = '*', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    ax9b.set_xscale('log')
    ax9b.set_yscale('log')
    ax9b.set_ylabel(r'Z$_{\mathrm{gas}}$')
    ax9b.set_xlabel(r'M$_{\mathrm{star}}$ [M$_\odot$]')
    #ax9b.axis([1e8, 1e11, 0.05, 0.85])
    fig9.tight_layout()
    sm = plt.cm.ScalarMappable(cmap=cmx, norm=cNorm)
    cbar = plt.colorbar(sm, location="right",ax=[ax9a,ax9b],aspect = 20) #pad=0.04)
    #cbar = mpl.colorbar.ColorbarBase(ax2sub, cmap=cmx, norm=cNorm, location="right",ax=[ax2a,ax2b]) #pad=0.04)
    cbar.set_label(r'D$_{\mathrm{massive}}$ [kpc]') #, rotation=0)
    fig9.show()
    fig9.set_size_inches(plt_width,plt_width*aspect_ratio*2)
    fig9.show()
    fig9.savefig(outfile_base + '_Mstar_Z.png',dpi = dpi)

    #Gas fraction
    plt.close('all')
    fig10 = plt.figure(10,figsize=(plt_width,plt_width*aspect_ratio))
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax10 = fig10.add_subplot(gs[0])
    ax10sub = fig10.add_subplot(gs[1]) 
    cNorm  = colors.LogNorm(vmin=10**1.5, vmax = 10**4)
    ax10.scatter(objs_rom['mstar'],objs_rom['mHI']/(objs_rom['mstar'] + objs_rom['mHI']),c = (objs_rom['massiveDist']),cmap = cmx, norm = cNorm,alpha = 0.4, s = markersize*0.5)
    ax10.scatter(objs_pd_comb['mstar'][objs_pd_comb['type']=='Central'],objs_pd_comb['mHI'][objs_pd_comb['type']=='Central']/(objs_pd_comb['mstar'][objs_pd_comb['type']=='Central'] + objs_pd_comb['mHI'][objs_pd_comb['type']=='Central']),c = (objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax10.scatter(objs_pd_comb['mstar'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb['mHI'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]/(objs_pd_comb['mstar'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')] + objs_pd_comb['mHI'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]),c = (objs_pd_comb['massiveDist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], linewidths = edgewidth)
    ax10.scatter(objs_pd_comb['mstar'][objs_pd_comb['type']=='Satellite'],objs_pd_comb['mHI'][objs_pd_comb['type']=='Satellite']/(objs_pd_comb['mstar'][objs_pd_comb['type']=='Satellite'] + objs_pd_comb['mHI'][objs_pd_comb['type']=='Satellite']),c = (objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k',marker = "*", s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    ax10.set_xscale('log')
    ax10.set_xlabel(mstar_axes_label)
    ax10.set_ylabel(r'$f_{\mathrm{HI}}$')
    ax10.axis([1e3, 6e9, -0.05, 1.05])    
    cb = mpl.colorbar.ColorbarBase(ax10sub, cmap=cmx, norm=cNorm)
    cb.set_label(r"D$_{\mathrm{massive}}$ [kpc]")
    fig10.tight_layout()    
    fig10.show()
    fig10.savefig(outfile_base + '_smass_fHI.png',dpi = dpi)
    """
    plt.close('all')
    fig10 = plt.figure(10,figsize=(plt_width,plt_width*aspect_ratio))
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax10 = fig10.add_subplot(gs[0])
    ax10sub = fig10.add_subplot(gs[1]) 
    cNorm  = colors.LogNorm(vmin=10**1.5, vmax = 10**4)
    ax10.fill_between(10**((mvir_bins[0:-1] + mvir_bins[1:])/2),fHI_bins0[0,:],fHI_bins0[4,:],color = sm.to_rgba(math.log10(dist_binsc[0])),alpha = 0.1,linewidth=0)
    ax10.fill_between(10**((mvir_bins[0:-1] + mvir_bins[1:])/2),fHI_bins0[1,:],fHI_bins0[3,:],color = sm.to_rgba(math.log10(dist_binsc[0])),alpha = 0.1,linewidth=0)    
    ax10.plot(10**((mvir_bins[0:-1] + mvir_bins[1:])/2),fHI_bins0[2,:],c = sm.to_rgba(math.log10(dist_binsc[0])))

    ax10.fill_between(10**((mvir_bins[0:-1] + mvir_bins[1:])/2),fHI_bins1[0,:],fHI_bins1[4,:],color = sm.to_rgba(math.log10(dist_binsc[1])),alpha = 0.1,linewidth=0)
    ax10.fill_between(10**((mvir_bins[0:-1] + mvir_bins[1:])/2),fHI_bins1[1,:],fHI_bins1[3,:],color = sm.to_rgba(math.log10(dist_binsc[1])),alpha = 0.1,linewidth=0)    
    ax10.plot(10**((mvir_bins[0:-1] + mvir_bins[1:])/2),fHI_bins1[2,:],c = sm.to_rgba(math.log10(dist_binsc[1])))

    ax10.fill_between(10**((mvir_bins[0:-1] + mvir_bins[1:])/2),fHI_bins2[0,:],fHI_bins2[4,:],color = sm.to_rgba(math.log10(dist_binsc[2])),alpha = 0.1,linewidth=0)
    ax10.fill_between(10**((mvir_bins[0:-1] + mvir_bins[1:])/2),fHI_bins2[1,:],fHI_bins2[3,:],color = sm.to_rgba(math.log10(dist_binsc[2])),alpha = 0.1,linewidth=0)    
    ax10.plot(10**((mvir_bins[0:-1] + mvir_bins[1:])/2),fHI_bins2[2,:],c = sm.to_rgba(math.log10(dist_binsc[2])))    
    
    ax10.fill_between(10**((mvir_bins[0:-1] + mvir_bins[1:])/2),fHI_bins3[0,:],fHI_bins3[4,:],color = sm.to_rgba(math.log10(dist_binsc[3])),alpha = 0.1,linewidth=0)
    ax10.fill_between(10**((mvir_bins[0:-1] + mvir_bins[1:])/2),fHI_bins3[1,:],fHI_bins3[3,:],color = sm.to_rgba(math.log10(dist_binsc[3])),alpha = 0.1,linewidth=0)    
    ax10.plot(10**((mvir_bins[0:-1] + mvir_bins[1:])/2),fHI_bins3[2,:],c = sm.to_rgba(math.log10(dist_binsc[3])))    
    #ax10.scatter(objs_rom['mvir'],objs_rom['mHI']/(objs_rom['mstar'] + objs_rom['mHI']),c = (objs_rom['massiveDist']),cmap = cmx, norm = cNorm,alpha = 0.4, s = markersize*0.5)
    ax10.scatter(objs_pd_comb['mvir'][objs_pd_comb['type']=='Central'],objs_pd_comb['mHI'][objs_pd_comb['type']=='Central']/(objs_pd_comb['mstar'][objs_pd_comb['type']=='Central'] + objs_pd_comb['mHI'][objs_pd_comb['type']=='Central']),c = (objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax10.scatter(objs_pd_comb['mvir'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb['mHI'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]/(objs_pd_comb['mstar'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')] + objs_pd_comb['mHI'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]),c = (objs_pd_comb['massiveDist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], linewidths = edgewidth)
    ax10.scatter(objs_pd_comb['mvir'][objs_pd_comb['type']=='Satellite'],objs_pd_comb['mHI'][objs_pd_comb['type']=='Satellite']/(objs_pd_comb['mstar'][objs_pd_comb['type']=='Satellite'] + objs_pd_comb['mHI'][objs_pd_comb['type']=='Satellite']),c = (objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k',marker = "*", s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    ax10.set_xscale('log')
    ax10.set_xlabel(r'M$_{vir, z = 0}$/M$_\odot$')
    ax10.set_ylabel(r'$f_{\mathrm{HI}}$')
    ax10.axis([1e9, 2e11, -0.05, 1.05])    
    cb = mpl.colorbar.ColorbarBase(ax10sub, cmap=cmx, norm=cNorm)
    cb.set_label(r"D$_{\mathrm{massive}}$ [kpc]")
    fig10.tight_layout()    
    fig10.show()
    fig10.savefig(outfile_base + '_vmass_fHI_rom.png',dpi = dpi)
    """
    plt.close('all')
    fig11 = plt.figure(2,figsize=(plt_width,plt_width*aspect_ratio))
    fig11.set_size_inches(plt_width,plt_width*aspect_ratio)
    fig11.clear()
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax11 = fig11.add_subplot(gs[0])
    ax11sub = fig11.add_subplot(gs[1])
    cNorm  = colors.LogNorm(vmin=10**1.5, vmax = 10**4)
    ax11.scatter(objs_pd_comb['tau90_vir'][objs_pd_comb['type']=='Central'],objs_pd_comb['tau90'][objs_pd_comb['type']=='Central'],c = (objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax11.scatter(objs_pd_comb['tau90_vir'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb['tau90'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = (objs_pd_comb['massiveDist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], linewidths = edgewidth)    
    ax11.scatter(objs_pd_comb['tau90_vir'][objs_pd_comb['type']=='Satellite'],objs_pd_comb['tau90'][objs_pd_comb['type']=='Satellite'],c = (objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite'], linewidths = edgewidth)
    
    ax11.scatter(objs_pd_comb['tau50_vir'][objs_pd_comb['type']=='Central'],objs_pd_comb['tau50'][objs_pd_comb['type']=='Central'],c = (objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth, marker = '+')
    ax11.scatter(objs_pd_comb['tau50_vir'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb['tau50'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = (objs_pd_comb['massiveDist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], linewidths = edgewidth, marker = '+')    
    ax11.scatter(objs_pd_comb['tau50_vir'][objs_pd_comb['type']=='Satellite'],objs_pd_comb['tau50'][objs_pd_comb['type']=='Satellite'],c = (objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k',s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite'], linewidths = edgewidth, marker = '+')    
    #ax11.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax11.set_xlabel(r'$\tau_{\mathrm{vir, t}}$')
    ax11.set_ylabel(r'$\tau_{t}$')
    ax11.axis([-0.5, 14, -0.5, 14])
    #legend1 = ax11.legend([cen_plt,sat_plt],['Central','Satellite'],scatterpoints = 1,facecolor = 'white',loc = 4,framealpha = 0,frameon = False)
    #ax11.add_artist(legend1)
    cb = mpl.colorbar.ColorbarBase(ax11sub, cmap=cmx, norm=cNorm)
    cb.set_label(r"D$_{\mathrm{massive}}$ [kpc]")
    fig11.tight_layout()
    fig11.show()    
    fig11.savefig(outfile_base + '_tau90_tau50.png',dpi = dpi)

    plt.close('all') #339924903
    fig11 = plt.figure(2,figsize=(plt_width,plt_width*aspect_ratio))
    fig11.set_size_inches(plt_width,plt_width*aspect_ratio)
    fig11.clear()
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax11 = fig11.add_subplot(gs[0])
    ax11sub = fig11.add_subplot(gs[1])
    #cNorm  = colors.LogNorm(vmin=10**1.5, vmax = 10**4)
    cNorm  = colors.LogNorm(vmin=1e4, vmax = 1e8)
    tau90_plt = ax11.scatter(objs_pd_comb['tau90_vir'][objs_pd_comb['type']=='Central'],objs_pd_comb['concentration'][objs_pd_comb['type']=='Central'],c = (objs_pd_comb['max_tindex'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax11.scatter(objs_pd_comb['tau90_vir'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb['concentration'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = (objs_pd_comb['max_tindex'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], linewidths = edgewidth)    
    ax11.scatter(objs_pd_comb['tau90_vir'][objs_pd_comb['type']=='Satellite'],objs_pd_comb['concentration'][objs_pd_comb['type']=='Satellite'],c = (objs_pd_comb['max_tindex'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite'], linewidths = edgewidth)
    """
    tau50_plt = ax11.scatter(objs_pd_comb['tau50_vir'][objs_pd_comb['type']=='Central'],objs_pd_comb['concentration'][objs_pd_comb['type']=='Central'],c = (objs_pd_comb['max_tindex'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth, marker = '+')
    ax11.scatter(objs_pd_comb['tau50_vir'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb['concentration'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = (objs_pd_comb['max_tindex'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], linewidths = edgewidth, marker = '+')    
    ax11.scatter(objs_pd_comb['tau50_vir'][objs_pd_comb['type']=='Satellite'],objs_pd_comb['concentration'][objs_pd_comb['type']=='Satellite'],c = (objs_pd_comb['max_tindex'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k',s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite'], linewidths = edgewidth, marker = '+')    
    """
    ax11.set_xlabel(r'$\tau_{\mathrm{vir, t}}$')
    ax11.set_ylabel('Concentration')
    ax11.axis([-0.5, 14, -0.5, 20])
    #legend1 = ax11.legend([tau50_plt, tau90_plt],[r'50%',r'90%'],frameon = False) 
    #legend1 = ax11.legend([cen_plt,sat_plt],['Central','Satellite'],scatterpoints = 1,facecolor = 'white',loc = 4,framealpha = 0,frameon = False)
    #ax11.add_artist(legend1)
    cb = mpl.colorbar.ColorbarBase(ax11sub, cmap=cmx, norm=cNorm)
    #cb.set_label(r"D$_{\mathrm{massive}}$ [kpc]")
    cb.set_label(r"Tidal Index")
    fig11.tight_layout()
    fig11.show()    
    fig11.savefig(outfile_base + '_tau90_tau50_concent.png',dpi = dpi)

    plt.close('all')
    fig11 = plt.figure(2,figsize=(plt_width,plt_width*aspect_ratio))
    fig11.set_size_inches(plt_width,plt_width*aspect_ratio)
    fig11.clear()
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax11 = fig11.add_subplot(gs[0])
    ax11sub = fig11.add_subplot(gs[1])
    #cNorm  = colors.LogNorm(vmin=10**1.5, vmax = 10**4)
    cNorm  = colors.LogNorm(vmin=1e4, vmax = 1e8)
    tau50_plt = ax11.scatter(objs_pd_comb['tau50'][objs_pd_comb['type']=='Central'],objs_pd_comb['concentration'][objs_pd_comb['type']=='Central'],c = (objs_pd_comb['max_tindex'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax11.scatter(objs_pd_comb['tau50'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb['concentration'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = (objs_pd_comb['max_tindex'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], linewidths = edgewidth)    
    ax11.scatter(objs_pd_comb['tau50'][objs_pd_comb['type']=='Satellite'],objs_pd_comb['concentration'][objs_pd_comb['type']=='Satellite'],c = (objs_pd_comb['max_tindex'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite'], linewidths = edgewidth)
    
    tau90_plt = ax11.scatter(objs_pd_comb['tau90'][objs_pd_comb['type']=='Central'],objs_pd_comb['concentration'][objs_pd_comb['type']=='Central'],c = (objs_pd_comb['max_tindex'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth, marker = '+')
    ax11.scatter(objs_pd_comb['tau90'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb['concentration'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = (objs_pd_comb['max_tindex'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], linewidths = edgewidth, marker = '+')    
    ax11.scatter(objs_pd_comb['tau90'][objs_pd_comb['type']=='Satellite'],objs_pd_comb['concentration'][objs_pd_comb['type']=='Satellite'],c = (objs_pd_comb['max_tindex'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k',s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite'], linewidths = edgewidth, marker = '+')    
    ax11.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax11.set_xlabel(r'$\tau_{\mathrm{vir, t}}$')
    ax11.set_ylabel('Concentration')
    #ax11.set_xlabel(r'$\tau_{t}$')
    #ax11.axis([-0.5, 14, -0.5, 30])
    legend1 = ax11.legend([tau50_plt, tau90_plt],[r'50%',r'90%'],frameon = False) 
    #legend1 = ax11.legend([cen_plt,sat_plt],['Central','Satellite'],scatterpoints = 1,facecolor = 'white',loc = 4,framealpha = 0,frameon = False)
    #ax11.add_artist(legend1)
    cb = mpl.colorbar.ColorbarBase(ax11sub, cmap=cmx, norm=cNorm)
    #cb.set_label(r"D$_{\mathrm{massive}}$ [kpc]")
    cb.set_label(r"Tidal Index")
    fig11.tight_layout()
    fig11.show()    
    fig11.savefig(outfile_base + '_star_tau90_tau50_concent.png',dpi = dpi)

    plt.close('all')
    fig11 = plt.figure(2,figsize=(plt_width,plt_width*aspect_ratio))
    fig11.set_size_inches(plt_width,plt_width*aspect_ratio)
    fig11.clear()
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax11 = fig11.add_subplot(gs[0])
    ax11sub = fig11.add_subplot(gs[1])
    cNorm  = colors.Normalize(vmin=8, vmax = 11.5)
    ax11.scatter(np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central']),objs_pd_comb['tau90'][objs_pd_comb['type']=='Central']/objs_pd_comb['tau90_vir'][objs_pd_comb['type']=='Central'],c = (objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax11.scatter(np.log10(objs_pd_comb['massiveDist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]),objs_pd_comb['tau90'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]/objs_pd_comb['tau90_vir'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = (objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], linewidths = edgewidth)    
    ax11.scatter(np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite']),objs_pd_comb['tau90'][objs_pd_comb['type']=='Satellite']/objs_pd_comb['tau90_vir'][objs_pd_comb['type']=='Satellite'],c = (objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite'], linewidths = edgewidth)
    
    ax11.scatter(np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central']),objs_pd_comb['tau50'][objs_pd_comb['type']=='Central']/objs_pd_comb['tau50_vir'][objs_pd_comb['type']=='Central'],c = (objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth, marker = '+')
    ax11.scatter(np.log10(objs_pd_comb['massiveDist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]),objs_pd_comb['tau50'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]/objs_pd_comb['tau50_vir'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = (objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], linewidths = edgewidth, marker = '+')    
    ax11.scatter(np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite']),objs_pd_comb['tau50'][objs_pd_comb['type']=='Satellite']/objs_pd_comb['tau50_vir'][objs_pd_comb['type']=='Satellite'],c = (objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k',s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite'], linewidths = edgewidth, marker = '+')    
    #ax11.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax11.set_xlabel(r'D$_{\mathrm{massive}}$ [kpc]')
    ax11.set_ylabel(r'$\tau_{t}/\tau_{\mathrm{vir, t}}$')
    ax11.axis([1.4, 4, 0, 4])
    #legend1 = ax11.legend([cen_plt,sat_plt],['Central','Satellite'],scatterpoints = 1,facecolor = 'white',loc = 4,framealpha = 0,frameon = False)
    #ax11.add_artist(legend1)
    cb = mpl.colorbar.ColorbarBase(ax11sub, cmap=cmx, norm=cNorm)
    cb.set_label(r"Log(M$_{\mathrm{peak}}$/M$_\odot$)")
    fig11.tight_layout()
    fig11.show()    
    fig11.savefig(outfile_base + '_tau90_tau50_dist.png',dpi = dpi)

    plt.close('all')
    fig11 = plt.figure(2,figsize=(plt_width,plt_width*aspect_ratio))
    fig11.set_size_inches(plt_width,plt_width*aspect_ratio)
    fig11.clear()
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax11 = fig11.add_subplot(gs[0])
    ax11sub = fig11.add_subplot(gs[1])
    cNorm  = colors.Normalize(vmin=8, vmax = 11.5)
    ax11.scatter(np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central']),objs_pd_comb['tau90'][objs_pd_comb['type']=='Central'],c = (np.array(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'].tolist())), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth, alpha = 0.4)
    ax11.scatter(np.log10(objs_pd_comb['massiveDist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]),objs_pd_comb['tau90'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = (np.array(objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')].tolist())), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], linewidths = edgewidth, alpha = 0.4)    
    ax11.scatter(np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite']),objs_pd_comb['tau90'][objs_pd_comb['type']=='Satellite'],c = (np.array(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'].tolist())), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite'], linewidths = edgewidth, alpha = 0.4)
    
    ax11.scatter(np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central']),objs_pd_comb['tau90_vir'][objs_pd_comb['type']=='Central'],c = (np.array(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'].tolist())), cmap = cmx, norm = cNorm,s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth*2, marker = '+')
    ax11.scatter(np.log10(objs_pd_comb['massiveDist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]),objs_pd_comb['tau90_vir'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = (np.array(objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')].tolist())), cmap = cmx, norm = cNorm,s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], linewidths = edgewidth*2, marker = '+')    
    ax11.scatter(np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite']),objs_pd_comb['tau90_vir'][objs_pd_comb['type']=='Satellite'],c = (np.array(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'].tolist())), cmap = cmx, norm = cNorm,s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite'], linewidths = edgewidth*2, marker = '+')
    '''
    ax11.scatter(np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central']),objs_pd_comb['tau50'][objs_pd_comb['type']=='Central'],c = 'None',edgecolor = 'r', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax11.scatter(np.log10(objs_pd_comb['massiveDist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]),objs_pd_comb['tau50'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = 'None',edgecolor = 'r', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], linewidths = edgewidth)    
    ax11.scatter(np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite']),objs_pd_comb['tau50'][objs_pd_comb['type']=='Satellite'],c = 'None',edgecolor = 'r', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite'], linewidths = edgewidth)
    
    ax11.scatter(np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central']),objs_pd_comb['tau50_vir'][objs_pd_comb['type']=='Central'],c = 'r',s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth, marker = '+')
    ax11.scatter(np.log10(objs_pd_comb['massiveDist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]),objs_pd_comb['tau50_vir'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = 'r',s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], linewidths = edgewidth, marker = '+')    
    ax11.scatter(np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite']),objs_pd_comb['tau50_vir'][objs_pd_comb['type']=='Satellite'],c = 'r',s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite'], linewidths = edgewidth, marker = '+')    
    '''
    ax11.set_xlabel(r'D$_{\mathrm{massive}}$ [kpc]')
    ax11.set_ylabel(r'$\tau_{90}$')
    ax11.axis([1.4, 4, -0.5, 14])
    #legend1 = ax11.legend([cen_plt,sat_plt],['Central','Satellite'],scatterpoints = 1,facecolor = 'white',loc = 4,framealpha = 0,frameon = False)
    #ax11.add_artist(legend1)
    cb = mpl.colorbar.ColorbarBase(ax11sub, cmap=cmx, norm=cNorm)
    cb.set_label(r"Log(M$_{\mathrm{*, z = 0}}$/M$_\odot$)")
    fig11.tight_layout()
    fig11.show()    
    fig11.savefig(outfile_base + '_tau90_dist.png',dpi = dpi)

    plt.close('all')
    fig11 = plt.figure(2,figsize=(plt_width,plt_width*aspect_ratio))
    fig11.set_size_inches(plt_width,plt_width*aspect_ratio)
    fig11.clear()
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax11 = fig11.add_subplot(gs[0])
    ax11sub = fig11.add_subplot(gs[1])
    cNorm  = colors.Normalize(vmin=8, vmax = 11.5)
    ax11.scatter(np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central']),objs_pd_comb['tau50'][objs_pd_comb['type']=='Central'],c = (np.array(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'].tolist())), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth, alpha = 0.4)
    ax11.scatter(np.log10(objs_pd_comb['massiveDist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]),objs_pd_comb['tau50'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = (np.array(objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')].tolist())), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], linewidths = edgewidth, alpha = 0.4)    
    ax11.scatter(np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite']),objs_pd_comb['tau50'][objs_pd_comb['type']=='Satellite'],c = (np.array(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'].tolist())), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite'], linewidths = edgewidth, alpha = 0.4)
    
    ax11.scatter(np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central']),objs_pd_comb['tau50_vir'][objs_pd_comb['type']=='Central'],c = (np.array(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'].tolist())), cmap = cmx, norm = cNorm,s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth*2, marker = '+')
    ax11.scatter(np.log10(objs_pd_comb['massiveDist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]),objs_pd_comb['tau50_vir'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = (np.array(objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')].tolist())), cmap = cmx, norm = cNorm,s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], linewidths = edgewidth*2, marker = '+')    
    ax11.scatter(np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite']),objs_pd_comb['tau50_vir'][objs_pd_comb['type']=='Satellite'],c = (np.array(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'].tolist())), cmap = cmx, norm = cNorm,s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite'], linewidths = edgewidth*2, marker = '+')

    ax11.set_xlabel(r'D$_{\mathrm{massive}}$ [kpc]')
    ax11.set_ylabel(r'$\tau_{50}$')
    ax11.axis([1.4, 4, -0.5, 14])
    #legend1 = ax11.legend([cen_v_tau50,sat_v_tau50],['Central','Satellite'],scatterpoints = 1,facecolor = 'white',loc = 4,framealpha = 0,frameon = False)
    #ax11.add_artist(legend1)
    cb = mpl.colorbar.ColorbarBase(ax11sub, cmap=cmx, norm=cNorm)
    cb.set_label(r"Log(M$_{\mathrm{*, z = 0}}$/M$_\odot$)")
    fig11.tight_layout()
    fig11.show()    
    fig11.savefig(outfile_base + '_tau50_dist.png',dpi = dpi)

    plt.close('all')    
    fig11 = plt.figure(figsize=(plt_width*2,plt_width*aspect_ratio*2))
    gs = fig11.add_gridspec(2,2, hspace=0, wspace=0) #, vspace=0)
    axs11 = gs.subplots(sharex = True, sharey = True) #, constrained_layout=True)
    axs11 = axs11.flatten()
    ax11a = axs11[0]
    ax11b = axs11[1]
    ax11c = axs11[2]
    ax11d = axs11[3]
    cmx = plt.get_cmap("winter_r")
    dcmap = discrete_cmap(256, cmx)
    cNorm  = colors.LogNorm(vmin=1e8, vmax = 1e11)
    dmin = 1.4
    dmax = 4
    nbins = 6
    quench_cond = np.array(list(zip(*objs_pd_comb['SFR']))[0])/objs_pd_comb['M_star'] < 1e-11
    sf_cond =  np.array(list(zip(*objs_pd_comb['SFR']))[0])/objs_pd_comb['M_star'] >= 1e-11    
    ax11a.scatter(objs_pd_comb['massiveDist'][sf_cond * (objs_pd_comb['type']=='Central')],objs_pd_comb['tau50'][sf_cond * (objs_pd_comb['type']=='Central')],c = (objs_pd_comb['Mpeak'][sf_cond * (objs_pd_comb['type']=='Central')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][sf_cond * (objs_pd_comb['type']=='Central')]*1.5, linewidths = edgewidth)
    ax11a.scatter(objs_pd_comb['massiveDist'][sf_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))],objs_pd_comb['tau50'][sf_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))],c = (objs_pd_comb['Mpeak'][sf_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][sf_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))]*1.5*1.5, marker = '*', linewidths = edgewidth)    
    ax11a.scatter(objs_pd_comb['massiveDist'][sf_cond * (objs_pd_comb['type']=='Satellite')],objs_pd_comb['tau50'][sf_cond * (objs_pd_comb['type']=='Satellite')],c = (objs_pd_comb['Mpeak'][sf_cond * (objs_pd_comb['type']=='Satellite')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][sf_cond * (objs_pd_comb['type']=='Satellite')]*1.5*1.5, marker = '*', linewidths = edgewidth)
    
    ax11a.scatter(objs_pd_comb['massiveDist'][quench_cond * (objs_pd_comb['type']=='Central')],objs_pd_comb['tau50'][quench_cond * (objs_pd_comb['type']=='Central')],edgecolor =  dcmap(cNorm(objs_pd_comb['Mpeak'][quench_cond * (objs_pd_comb['type']=='Central')])), cmap = cmx, norm = cNorm, facecolor = 'none', s = objs_pd_comb['p_size'][quench_cond * (objs_pd_comb['type']=='Central')]*1.5, linewidths = edgewidth*2)
    ax11a.scatter(objs_pd_comb['massiveDist'][quench_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))],objs_pd_comb['tau50'][quench_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))],edgecolor = dcmap(cNorm(objs_pd_comb['Mpeak'][quench_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))])), cmap = cmx, norm = cNorm,facecolor = 'none', s = objs_pd_comb['p_size'][quench_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))]*1.5*1.5, marker = '*', linewidths = edgewidth*2)    
    ax11a.scatter(objs_pd_comb['massiveDist'][quench_cond * (objs_pd_comb['type']=='Satellite')],objs_pd_comb['tau50'][quench_cond * (objs_pd_comb['type']=='Satellite')],edgecolor =  dcmap(cNorm(objs_pd_comb['Mpeak'][quench_cond * (objs_pd_comb['type']=='Satellite')])), cmap = cmx, norm = cNorm,facecolor = 'none', s = objs_pd_comb['p_size'][quench_cond * (objs_pd_comb['type']=='Satellite')]*1.5*1.5, marker = '*', linewidths = edgewidth*2)
    
    xaxis, data = bin_plt(np.log10(objs_pd_comb['massiveDist'][(objs_pd_comb['Mpeak'] > 3.16e9) & (objs_pd_comb['tau90_vir'] > 0)]), objs_pd_comb['tau50'][ (objs_pd_comb['Mpeak'] > 3.16e9) & (objs_pd_comb['tau90_vir'] > 0)], xmax = dmax, nbins = nbins)
    xaxis2, data2 = bin_plt(np.log10(objs_pd_comb['massiveDist'][(objs_pd_comb['type'] == 'Central') & (objs_pd_comb['Mpeak'] > 3.16e9) & (objs_pd_comb['tau90_vir'] > 0)]), objs_pd_comb['tau50'][(objs_pd_comb['type'] == 'Central') & (objs_pd_comb['Mpeak'] > 3.16e9) & (objs_pd_comb['tau90_vir'] > 0)], xmax = dmax, nbins = nbins)
    ax11a.fill_between(10**xaxis, data[:,0], data[:,4], color='k', alpha=0.1, linewidth=0)
    ax11a.fill_between(10**xaxis, data[:,1], data[:,3], color='k', alpha=0.1, linewidth=0)
    ax11a.plot(10**xaxis, data[:,2], c='k')
    #ax11a.plot(10**xaxis2, data2[:,2], "--k")
    ax11a.set_xscale('log')
    ax11a.set_ylabel(r'$\tau_{50}$')
    ax11a.text(30,12,'Stars')
    ax11a.axis([17, 1e4, -0.5, 14])

    ax11b.scatter(objs_pd_comb['massiveDist'][sf_cond * (objs_pd_comb['type']=='Central')],objs_pd_comb['tau50_vir'][sf_cond * (objs_pd_comb['type']=='Central')],c = (objs_pd_comb['Mpeak'][sf_cond * (objs_pd_comb['type']=='Central')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][sf_cond * (objs_pd_comb['type']=='Central')]*1.5, linewidths = edgewidth)
    ax11b.scatter(objs_pd_comb['massiveDist'][sf_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))],objs_pd_comb['tau50_vir'][sf_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))],c = (objs_pd_comb['Mpeak'][sf_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][sf_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))]*1.5*1.5, marker = '*', linewidths = edgewidth)    
    ax11b.scatter(objs_pd_comb['massiveDist'][sf_cond * (objs_pd_comb['type']=='Satellite')],objs_pd_comb['tau50_vir'][sf_cond * (objs_pd_comb['type']=='Satellite')],c = (objs_pd_comb['Mpeak'][sf_cond * (objs_pd_comb['type']=='Satellite')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][sf_cond * (objs_pd_comb['type']=='Satellite')]*1.5*1.5, marker = '*', linewidths = edgewidth)
    
    ax11b.scatter(objs_pd_comb['massiveDist'][quench_cond * (objs_pd_comb['type']=='Central')],objs_pd_comb['tau50_vir'][quench_cond * (objs_pd_comb['type']=='Central')],edgecolor =  dcmap(cNorm(objs_pd_comb['Mpeak'][quench_cond * (objs_pd_comb['type']=='Central')])), cmap = cmx, norm = cNorm, facecolor = 'none', s = objs_pd_comb['p_size'][quench_cond * (objs_pd_comb['type']=='Central')]*1.5, linewidths = edgewidth*2)
    ax11b.scatter(objs_pd_comb['massiveDist'][quench_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))],objs_pd_comb['tau50_vir'][quench_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))],edgecolor = dcmap(cNorm(objs_pd_comb['Mpeak'][quench_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))])), cmap = cmx, norm = cNorm,facecolor = 'none', s = objs_pd_comb['p_size'][quench_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))]*1.5*1.5, marker = '*', linewidths = edgewidth*2)    
    ax11b.scatter(objs_pd_comb['massiveDist'][quench_cond * (objs_pd_comb['type']=='Satellite')],objs_pd_comb['tau50_vir'][quench_cond * (objs_pd_comb['type']=='Satellite')],edgecolor =  dcmap(cNorm(objs_pd_comb['Mpeak'][quench_cond * (objs_pd_comb['type']=='Satellite')])), cmap = cmx, norm = cNorm,facecolor = 'none', s = objs_pd_comb['p_size'][quench_cond * (objs_pd_comb['type']=='Satellite')]*1.5*1.5, marker = '*', linewidths = edgewidth*2)    

    xaxis, data = bin_plt(np.log10(objs_pd_comb['massiveDist'][(objs_pd_comb['Mpeak'] > 3.16e9) & (objs_pd_comb['tau90_vir'] > 0)]), objs_pd_comb['tau50_vir'][(objs_pd_comb['Mpeak'] > 3.16e9) & (objs_pd_comb['tau90_vir'] > 0)], xmax = dmax, nbins = nbins)
    xaxis2, data2 = bin_plt(np.log10(objs_pd_comb['massiveDist'][(objs_pd_comb['type'] == 'Central') & (objs_pd_comb['Mpeak'] > 3.16e9) & (objs_pd_comb['tau90_vir'] > 0)]), objs_pd_comb['tau50_vir'][(objs_pd_comb['type'] == 'Central') & (objs_pd_comb['Mpeak'] > 3.16e9) & (objs_pd_comb['tau90_vir'] > 0)], xmax = dmax, nbins = nbins)    
    ax11b.fill_between(10**xaxis, data[:,0], data[:,4], color='k', alpha=0.1, linewidth=0)
    ax11b.fill_between(10**xaxis, data[:,1], data[:,3], color='k', alpha=0.1, linewidth=0)
    ax11b.plot(10**xaxis, data[:,2], c='k')
    #ax11b.plot(10**xaxis2, data2[:,2], "--k")
    ax11b.set_xscale('log')
    ax11b.text(30,12,'Total Mass')
    
    ax11c.scatter(objs_pd_comb['massiveDist'][sf_cond * (objs_pd_comb['type']=='Central')],objs_pd_comb['tau90'][sf_cond * (objs_pd_comb['type']=='Central')],c = (objs_pd_comb['Mpeak'][sf_cond * (objs_pd_comb['type']=='Central')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][sf_cond * (objs_pd_comb['type']=='Central')]*1.5, linewidths = edgewidth)
    ax11c.scatter(objs_pd_comb['massiveDist'][sf_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))],objs_pd_comb['tau90'][sf_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))],c = (objs_pd_comb['Mpeak'][sf_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][sf_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))]*1.5*1.5, marker = '*', linewidths = edgewidth)    
    ax11c.scatter(objs_pd_comb['massiveDist'][sf_cond * (objs_pd_comb['type']=='Satellite')],objs_pd_comb['tau90'][sf_cond * (objs_pd_comb['type']=='Satellite')],c = (objs_pd_comb['Mpeak'][sf_cond * (objs_pd_comb['type']=='Satellite')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][sf_cond * (objs_pd_comb['type']=='Satellite')]*1.5*1.5, marker = '*', linewidths = edgewidth)
    
    ax11c.scatter(objs_pd_comb['massiveDist'][quench_cond * (objs_pd_comb['type']=='Central')],objs_pd_comb['tau90'][quench_cond * (objs_pd_comb['type']=='Central')],edgecolor =  dcmap(cNorm(objs_pd_comb['Mpeak'][quench_cond * (objs_pd_comb['type']=='Central')])), cmap = cmx, norm = cNorm, facecolor = 'none', s = objs_pd_comb['p_size'][quench_cond * (objs_pd_comb['type']=='Central')]*1.5, linewidths = edgewidth*2)
    ax11c.scatter(objs_pd_comb['massiveDist'][quench_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))],objs_pd_comb['tau90'][quench_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))],edgecolor = dcmap(cNorm(objs_pd_comb['Mpeak'][quench_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))])), cmap = cmx, norm = cNorm,facecolor = 'none', s = objs_pd_comb['p_size'][quench_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))]*1.5*1.5, marker = '*', linewidths = edgewidth*2)    
    ax11c.scatter(objs_pd_comb['massiveDist'][quench_cond * (objs_pd_comb['type']=='Satellite')],objs_pd_comb['tau90'][quench_cond * (objs_pd_comb['type']=='Satellite')],edgecolor =  dcmap(cNorm(objs_pd_comb['Mpeak'][quench_cond * (objs_pd_comb['type']=='Satellite')])), cmap = cmx, norm = cNorm,facecolor = 'none', s = objs_pd_comb['p_size'][quench_cond * (objs_pd_comb['type']=='Satellite')]*1.5*1.5, marker = '*', linewidths = edgewidth*2)
    
    xaxis, data = bin_plt(np.log10(objs_pd_comb['massiveDist'][(objs_pd_comb['Mpeak'] > 3.16e9) & (objs_pd_comb['tau90_vir'] > 0)]), objs_pd_comb['tau90'][(objs_pd_comb['Mpeak'] > 3.16e9) & (objs_pd_comb['tau90_vir'] > 0)], xmax = dmax, nbins = nbins)
    xaxis2, data2 = bin_plt(np.log10(objs_pd_comb['massiveDist'][(objs_pd_comb['type'] == 'Central') & (objs_pd_comb['Mpeak'] > 3.16e9) & (objs_pd_comb['tau90_vir'] > 0)]), objs_pd_comb['tau90'][(objs_pd_comb['type'] == 'Central') & (objs_pd_comb['Mpeak'] > 3.16e9) & (objs_pd_comb['tau90_vir'] > 0)], xmax = dmax, nbins = nbins)
    ax11c.fill_between(10**xaxis, data[:,0], data[:,4], color='k', alpha=0.1, linewidth=0)
    ax11c.fill_between(10**xaxis, data[:,1], data[:,3], color='k', alpha=0.1, linewidth=0)
    ax11c.plot(10**xaxis, data[:,2], c='k')
    #ax11c.plot(10**xaxis2, data2[:,2], "--k")
    ax11c.set_xscale('log')
    ax11c.set_ylabel(r'$\tau_{90}$')
    ax11c.set_xlabel(r'D$_{\mathrm{massive}}$ [kpc]')

    cen_plt = ax11d.scatter(objs_pd_comb['massiveDist'][sf_cond * (objs_pd_comb['type']=='Central')],objs_pd_comb['tau90_vir'][sf_cond * (objs_pd_comb['type']=='Central')],c = (objs_pd_comb['Mpeak'][sf_cond * (objs_pd_comb['type']=='Central')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][sf_cond * (objs_pd_comb['type']=='Central')]*1.5, linewidths = edgewidth)
    sf_plt = ax11d.scatter(objs_pd_comb['massiveDist'][sf_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))],objs_pd_comb['tau90_vir'][sf_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))],c = (objs_pd_comb['Mpeak'][sf_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][sf_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))]*1.5*1.5, marker = '*', linewidths = edgewidth)    
    sat_plt = ax11d.scatter(objs_pd_comb['massiveDist'][sf_cond * (objs_pd_comb['type']=='Satellite')],objs_pd_comb['tau90_vir'][sf_cond * (objs_pd_comb['type']=='Satellite')],c = (objs_pd_comb['Mpeak'][sf_cond * (objs_pd_comb['type']=='Satellite')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][sf_cond * (objs_pd_comb['type']=='Satellite')]*1.5*1.5, marker = '*', linewidths = edgewidth)
    
    ax11d.scatter(objs_pd_comb['massiveDist'][quench_cond * (objs_pd_comb['type']=='Central')],objs_pd_comb['tau90_vir'][quench_cond * (objs_pd_comb['type']=='Central')],edgecolor =  dcmap(cNorm(objs_pd_comb['Mpeak'][quench_cond * (objs_pd_comb['type']=='Central')])), cmap = cmx, norm = cNorm, facecolor = 'none', s = objs_pd_comb['p_size'][quench_cond * (objs_pd_comb['type']=='Central')]*1.5, linewidths = edgewidth*2)
    quench_plt = ax11d.scatter(objs_pd_comb['massiveDist'][quench_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))],objs_pd_comb['tau90_vir'][quench_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))],edgecolor = dcmap(cNorm(objs_pd_comb['Mpeak'][quench_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))])), cmap = cmx, norm = cNorm,facecolor = 'none', s = objs_pd_comb['p_size'][quench_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))]*1.5*1.5, marker = '*', linewidths = edgewidth*2)    
    ax11d.scatter(objs_pd_comb['massiveDist'][quench_cond * (objs_pd_comb['type']=='Satellite')],objs_pd_comb['tau90_vir'][quench_cond * (objs_pd_comb['type']=='Satellite')],edgecolor =  dcmap(cNorm(objs_pd_comb['Mpeak'][quench_cond * (objs_pd_comb['type']=='Satellite')])), cmap = cmx, norm = cNorm,facecolor = 'none', s = objs_pd_comb['p_size'][quench_cond * (objs_pd_comb['type']=='Satellite')]*1.5*1.5, marker = '*', linewidths = edgewidth*2)    
    
    xaxis, data = bin_plt(np.log10(objs_pd_comb['massiveDist'][(objs_pd_comb['Mpeak'] > 3.16e9) & (objs_pd_comb['tau90_vir'] > 0)]), objs_pd_comb['tau90_vir'][(objs_pd_comb['Mpeak'] > 3.16e9) & (objs_pd_comb['tau90_vir'] > 0)],  xmax = dmax, nbins = nbins)
    xaxis2, data2 = bin_plt(np.log10(objs_pd_comb['massiveDist'][(objs_pd_comb['type'] == 'Central') & (objs_pd_comb['Mpeak'] > 3.16e9) & (objs_pd_comb['tau90_vir'] > 0)]), objs_pd_comb['tau90_vir'][(objs_pd_comb['type'] == 'Central') &  (objs_pd_comb['Mpeak'] > 3.16e9) & (objs_pd_comb['tau90_vir'] > 0)],  xmax = dmax, nbins = nbins)
    ax11d.fill_between(10**xaxis, data[:,0], data[:,4], color='k', alpha=0.1, linewidth=0)
    ax11d.fill_between(10**xaxis, data[:,1], data[:,3], color='k', alpha=0.1, linewidth=0)
    ax11d.plot(10**xaxis, data[:,2], c='k')
    #ax11d.plot(10**xaxis2, data2[:,2], "--k")
    ax11d.set_xscale('log')
    ax11d.axis([17, 1e4, -0.5, 14])
    ax11d.set_xlabel(r'D$_{\mathrm{massive}}$ [kpc]')
    ax11b.legend([cen_plt,sat_plt],['Central','Satellite/Backsplash'],scatterpoints = 1,facecolor = 'white',framealpha = 0,frameon = False)
    ax11d.legend([sf_plt,quench_plt],['SF','Quenched'],scatterpoints = 1,facecolor = 'white',framealpha = 0,frameon = False)
    fig11.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.8, wspace=0, hspace=0)
    cb_ax = fig11.add_axes([0.83, 0.1, 0.02, 0.8])
    sm = plt.cm.ScalarMappable(cmap=cmx, norm=cNorm)
    cbar = fig11.colorbar(sm,cax = cb_ax) #pad=0.04)
    cbar.set_label(r'Log(M$_{\mathrm{vir, peak}}$/M$_\odot$)') #, rotation=0)
    fig11.show()
    fig11.savefig(outfile_base + '_rMassGal_taus.png',dpi = dpi)

    plt.close('all')    
    fig11 = plt.figure(figsize=(plt_width*2,plt_width*aspect_ratio*2))
    gs = fig11.add_gridspec(2,2, hspace=0, wspace=0) #, vspace=0)
    axs11 = gs.subplots(sharex = True, sharey = True) #, constrained_layout=True)
    axs11 = axs11.flatten()
    ax11a = axs11[0]
    ax11b = axs11[1]
    ax11c = axs11[2]
    ax11d = axs11[3]
    dcmap = discrete_cmap(256, cmx)
    cNorm  = colors.LogNorm(vmin=1e8, vmax = 1e11)
    dmin = 1
    dmax = 3.2
    nbins = 6
    quench_cond = np.array(list(zip(*objs_pd_comb['SFR']))[0])/objs_pd_comb['M_star'] < 1e-11
    sf_cond =  np.array(list(zip(*objs_pd_comb['SFR']))[0])/objs_pd_comb['M_star'] >= 1e-11    
    ax11a.scatter(objs_pd_comb['min_dist'][sf_cond * (objs_pd_comb['type']=='Central')],objs_pd_comb['tau50'][sf_cond * (objs_pd_comb['type']=='Central')],c = (objs_pd_comb['Mpeak'][sf_cond * (objs_pd_comb['type']=='Central')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][sf_cond * (objs_pd_comb['type']=='Central')]*1.5, linewidths = edgewidth)
    ax11a.scatter(objs_pd_comb['min_dist'][sf_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))],objs_pd_comb['tau50'][sf_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))],c = (objs_pd_comb['Mpeak'][sf_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][sf_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))]*1.5*1.5, marker = '*', linewidths = edgewidth)    
    ax11a.scatter(objs_pd_comb['min_dist'][sf_cond * (objs_pd_comb['type']=='Satellite')],objs_pd_comb['tau50'][sf_cond * (objs_pd_comb['type']=='Satellite')],c = (objs_pd_comb['Mpeak'][sf_cond * (objs_pd_comb['type']=='Satellite')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][sf_cond * (objs_pd_comb['type']=='Satellite')]*1.5*1.5, marker = '*', linewidths = edgewidth)
    
    ax11a.scatter(objs_pd_comb['min_dist'][quench_cond * (objs_pd_comb['type']=='Central')],objs_pd_comb['tau50'][quench_cond * (objs_pd_comb['type']=='Central')],edgecolor =  dcmap(cNorm(objs_pd_comb['Mpeak'][quench_cond * (objs_pd_comb['type']=='Central')])), cmap = cmx, norm = cNorm, facecolor = 'none', s = objs_pd_comb['p_size'][quench_cond * (objs_pd_comb['type']=='Central')]*1.5, linewidths = edgewidth*2)
    ax11a.scatter(objs_pd_comb['min_dist'][quench_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))],objs_pd_comb['tau50'][quench_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))],edgecolor = dcmap(cNorm(objs_pd_comb['Mpeak'][quench_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))])), cmap = cmx, norm = cNorm,facecolor = 'none', s = objs_pd_comb['p_size'][quench_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))]*1.5*1.5, marker = '*', linewidths = edgewidth*2)    
    ax11a.scatter(objs_pd_comb['min_dist'][quench_cond * (objs_pd_comb['type']=='Satellite')],objs_pd_comb['tau50'][quench_cond * (objs_pd_comb['type']=='Satellite')],edgecolor =  dcmap(cNorm(objs_pd_comb['Mpeak'][quench_cond * (objs_pd_comb['type']=='Satellite')])), cmap = cmx, norm = cNorm,facecolor = 'none', s = objs_pd_comb['p_size'][quench_cond * (objs_pd_comb['type']=='Satellite')]*1.5*1.5, marker = '*', linewidths = edgewidth*2)
    
    xaxis, data = bin_plt(np.log10((objs_pd_comb['min_dist'][(objs_pd_comb['Mpeak'] > 3.16e9) & (objs_pd_comb['tau90_vir'] > 0)]).tolist()), objs_pd_comb['tau50'][(objs_pd_comb['Mpeak'] > 3.16e9) & (objs_pd_comb['tau90_vir'] > 0)], xmax = dmax, nbins = nbins)
    xaxis2, data2 = bin_plt(np.log10((objs_pd_comb['min_dist'][(objs_pd_comb['type'] == 'Central') & (objs_pd_comb['Mpeak'] > 3.16e9) & (objs_pd_comb['tau90_vir'] > 0)]).tolist()), objs_pd_comb['tau50'][(objs_pd_comb['type'] == 'Central') & (objs_pd_comb['Mpeak'] > 3.16e9) & (objs_pd_comb['tau90_vir'] > 0)], xmax = dmax, nbins = nbins)
    ax11a.fill_between(10**xaxis, data[:,0], data[:,4], color='k', alpha=0.1, linewidth=0)
    ax11a.fill_between(10**xaxis, data[:,1], data[:,3], color='k', alpha=0.1, linewidth=0)
    ax11a.plot(10**xaxis, data[:,2], c='k')
    #ax11a.plot(10**xaxis2, data2[:,2], "--k")
    ax11a.set_xscale('log')
    ax11a.set_ylabel(r'$\tau_{50}$')
    ax11a.text(30,12,'Stars')
    ax11a.axis([10, 4e3, -0.5, 14])

    ax11b.scatter(objs_pd_comb['min_dist'][sf_cond * (objs_pd_comb['type']=='Central')],objs_pd_comb['tau50_vir'][sf_cond * (objs_pd_comb['type']=='Central')],c = (objs_pd_comb['Mpeak'][sf_cond * (objs_pd_comb['type']=='Central')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][sf_cond * (objs_pd_comb['type']=='Central')]*1.5, linewidths = edgewidth)
    ax11b.scatter(objs_pd_comb['min_dist'][sf_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))],objs_pd_comb['tau50_vir'][sf_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))],c = (objs_pd_comb['Mpeak'][sf_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][sf_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))]*1.5*1.5, marker = '*', linewidths = edgewidth)    
    ax11b.scatter(objs_pd_comb['min_dist'][sf_cond * (objs_pd_comb['type']=='Satellite')],objs_pd_comb['tau50_vir'][sf_cond * (objs_pd_comb['type']=='Satellite')],c = (objs_pd_comb['Mpeak'][sf_cond * (objs_pd_comb['type']=='Satellite')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][sf_cond * (objs_pd_comb['type']=='Satellite')]*1.5*1.5, marker = '*', linewidths = edgewidth)
    
    ax11b.scatter(objs_pd_comb['min_dist'][quench_cond * (objs_pd_comb['type']=='Central')],objs_pd_comb['tau50_vir'][quench_cond * (objs_pd_comb['type']=='Central')],edgecolor =  dcmap(cNorm(objs_pd_comb['Mpeak'][quench_cond * (objs_pd_comb['type']=='Central')])), cmap = cmx, norm = cNorm, facecolor = 'none', s = objs_pd_comb['p_size'][quench_cond * (objs_pd_comb['type']=='Central')]*1.5, linewidths = edgewidth*2)
    ax11b.scatter(objs_pd_comb['min_dist'][quench_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))],objs_pd_comb['tau50_vir'][quench_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))],edgecolor = dcmap(cNorm(objs_pd_comb['Mpeak'][quench_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))])), cmap = cmx, norm = cNorm,facecolor = 'none', s = objs_pd_comb['p_size'][quench_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))]*1.5*1.5, marker = '*', linewidths = edgewidth*2)    
    ax11b.scatter(objs_pd_comb['min_dist'][quench_cond * (objs_pd_comb['type']=='Satellite')],objs_pd_comb['tau50_vir'][quench_cond * (objs_pd_comb['type']=='Satellite')],edgecolor =  dcmap(cNorm(objs_pd_comb['Mpeak'][quench_cond * (objs_pd_comb['type']=='Satellite')])), cmap = cmx, norm = cNorm,facecolor = 'none', s = objs_pd_comb['p_size'][quench_cond * (objs_pd_comb['type']=='Satellite')]*1.5*1.5, marker = '*', linewidths = edgewidth*2)    

    xaxis, data = bin_plt(np.log10((objs_pd_comb['min_dist'][(objs_pd_comb['Mpeak'] > 3.16e9) & (objs_pd_comb['tau90_vir'] > 0)]).tolist()), objs_pd_comb['tau50_vir'][(objs_pd_comb['Mpeak'] > 3.16e9) & (objs_pd_comb['tau90_vir'] > 0)], xmax = dmax, nbins = nbins)
    xaxis2, data2 = bin_plt(np.log10((objs_pd_comb['min_dist'][(objs_pd_comb['type'] == 'Central') & (objs_pd_comb['Mpeak'] > 3.16e9) & (objs_pd_comb['tau90_vir'] > 0)]).tolist()), objs_pd_comb['tau50_vir'][(objs_pd_comb['type'] == 'Central') & (objs_pd_comb['Mpeak'] > 3.16e9) & (objs_pd_comb['tau90_vir'] > 0)], xmax = dmax, nbins = nbins)
    ax11b.fill_between(10**xaxis, data[:,0], data[:,4], color='k', alpha=0.1, linewidth=0)
    ax11b.fill_between(10**xaxis, data[:,1], data[:,3], color='k', alpha=0.1, linewidth=0)
    ax11b.plot(10**xaxis, data[:,2], c='k')
    #ax11b.plot(10**xaxis2, data2[:,2], "--k")
    ax11b.set_xscale('log')
    ax11b.text(30,12,'Total Mass')
    
    ax11c.scatter(objs_pd_comb['min_dist'][sf_cond * (objs_pd_comb['type']=='Central')],objs_pd_comb['tau90'][sf_cond * (objs_pd_comb['type']=='Central')],c = (objs_pd_comb['Mpeak'][sf_cond * (objs_pd_comb['type']=='Central')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][sf_cond * (objs_pd_comb['type']=='Central')]*1.5, linewidths = edgewidth)
    ax11c.scatter(objs_pd_comb['min_dist'][sf_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))],objs_pd_comb['tau90'][sf_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))],c = (objs_pd_comb['Mpeak'][sf_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][sf_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))]*1.5*1.5, marker = '*', linewidths = edgewidth)    
    ax11c.scatter(objs_pd_comb['min_dist'][sf_cond * (objs_pd_comb['type']=='Satellite')],objs_pd_comb['tau90'][sf_cond * (objs_pd_comb['type']=='Satellite')],c = (objs_pd_comb['Mpeak'][sf_cond * (objs_pd_comb['type']=='Satellite')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][sf_cond * (objs_pd_comb['type']=='Satellite')]*1.5*1.5, marker = '*', linewidths = edgewidth)
    
    ax11c.scatter(objs_pd_comb['min_dist'][quench_cond * (objs_pd_comb['type']=='Central')],objs_pd_comb['tau90'][quench_cond * (objs_pd_comb['type']=='Central')],edgecolor =  dcmap(cNorm(objs_pd_comb['Mpeak'][quench_cond * (objs_pd_comb['type']=='Central')])), cmap = cmx, norm = cNorm, facecolor = 'none', s = objs_pd_comb['p_size'][quench_cond * (objs_pd_comb['type']=='Central')]*1.5, linewidths = edgewidth*2)
    ax11c.scatter(objs_pd_comb['min_dist'][quench_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))],objs_pd_comb['tau90'][quench_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))],edgecolor = dcmap(cNorm(objs_pd_comb['Mpeak'][quench_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))])), cmap = cmx, norm = cNorm,facecolor = 'none', s = objs_pd_comb['p_size'][quench_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))]*1.5*1.5, marker = '*', linewidths = edgewidth*2)    
    ax11c.scatter(objs_pd_comb['min_dist'][quench_cond * (objs_pd_comb['type']=='Satellite')],objs_pd_comb['tau90'][quench_cond * (objs_pd_comb['type']=='Satellite')],edgecolor =  dcmap(cNorm(objs_pd_comb['Mpeak'][quench_cond * (objs_pd_comb['type']=='Satellite')])), cmap = cmx, norm = cNorm,facecolor = 'none', s = objs_pd_comb['p_size'][quench_cond * (objs_pd_comb['type']=='Satellite')]*1.5*1.5, marker = '*', linewidths = edgewidth*2)
    
    xaxis, data = bin_plt(np.log10((objs_pd_comb['min_dist'][(objs_pd_comb['Mpeak'] > 3.16e9) & (objs_pd_comb['tau90_vir'] > 0)]).tolist()), objs_pd_comb['tau90'][(objs_pd_comb['Mpeak'] > 3.16e9) & (objs_pd_comb['tau90_vir'] > 0)], xmax = dmax, nbins = nbins)
    xaxis2, data2 = bin_plt(np.log10((objs_pd_comb['min_dist'][(objs_pd_comb['type'] == 'Central') & (objs_pd_comb['Mpeak'] > 3.16e9) & (objs_pd_comb['tau90_vir'] > 0)]).tolist()), objs_pd_comb['tau90'][(objs_pd_comb['type'] == 'Central') & (objs_pd_comb['Mpeak'] > 3.16e9) & (objs_pd_comb['tau90_vir'] > 0)], xmax = dmax, nbins = nbins)    
    ax11c.fill_between(10**xaxis, data[:,0], data[:,4], color='k', alpha=0.1, linewidth=0)
    ax11c.fill_between(10**xaxis, data[:,1], data[:,3], color='k', alpha=0.1, linewidth=0)
    ax11c.plot(10**xaxis, data[:,2], c='k')
    #ax11c.plot(10**xaxis2, data2[:,2], "--k")
    ax11c.set_xscale('log')
    ax11c.set_ylabel(r'$\tau_{90}$')
    ax11c.set_xlabel(r'Min(D$_{\mathrm{massive}}$) [kpc]')

    cen_plt = ax11d.scatter(objs_pd_comb['min_dist'][sf_cond * (objs_pd_comb['type']=='Central')],objs_pd_comb['tau90_vir'][sf_cond * (objs_pd_comb['type']=='Central')],c = (objs_pd_comb['Mpeak'][sf_cond * (objs_pd_comb['type']=='Central')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][sf_cond * (objs_pd_comb['type']=='Central')]*1.5, linewidths = edgewidth)
    sf_plt = ax11d.scatter(objs_pd_comb['min_dist'][sf_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))],objs_pd_comb['tau90_vir'][sf_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))],c = (objs_pd_comb['Mpeak'][sf_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][sf_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))]*1.5*1.5, marker = '*', linewidths = edgewidth)    
    sat_plt = ax11d.scatter(objs_pd_comb['min_dist'][sf_cond * (objs_pd_comb['type']=='Satellite')],objs_pd_comb['tau90_vir'][sf_cond * (objs_pd_comb['type']=='Satellite')],c = (objs_pd_comb['Mpeak'][sf_cond * (objs_pd_comb['type']=='Satellite')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][sf_cond * (objs_pd_comb['type']=='Satellite')]*1.5*1.5, marker = '*', linewidths = edgewidth)
    
    ax11d.scatter(objs_pd_comb['min_dist'][quench_cond * (objs_pd_comb['type']=='Central')],objs_pd_comb['tau90_vir'][quench_cond * (objs_pd_comb['type']=='Central')],edgecolor =  dcmap(cNorm(objs_pd_comb['Mpeak'][quench_cond * (objs_pd_comb['type']=='Central')])), cmap = cmx, norm = cNorm, facecolor = 'none', s = objs_pd_comb['p_size'][quench_cond * (objs_pd_comb['type']=='Central')]*1.5, linewidths = edgewidth*2)
    quench_plt = ax11d.scatter(objs_pd_comb['min_dist'][quench_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))],objs_pd_comb['tau90_vir'][quench_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))],edgecolor = dcmap(cNorm(objs_pd_comb['Mpeak'][quench_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))])), cmap = cmx, norm = cNorm,facecolor = 'none', s = objs_pd_comb['p_size'][quench_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))]*1.5*1.5, marker = '*', linewidths = edgewidth*2)    
    ax11d.scatter(objs_pd_comb['min_dist'][quench_cond * (objs_pd_comb['type']=='Satellite')],objs_pd_comb['tau90_vir'][quench_cond * (objs_pd_comb['type']=='Satellite')],edgecolor =  dcmap(cNorm(objs_pd_comb['Mpeak'][quench_cond * (objs_pd_comb['type']=='Satellite')])), cmap = cmx, norm = cNorm,facecolor = 'none', s = objs_pd_comb['p_size'][quench_cond * (objs_pd_comb['type']=='Satellite')]*1.5*1.5, marker = '*', linewidths = edgewidth*2)    
    
    xaxis, data = bin_plt(np.log10((objs_pd_comb['min_dist'][(objs_pd_comb['Mpeak'] > 3.16e9) & (objs_pd_comb['tau90_vir'] > 0)]).tolist()), objs_pd_comb['tau90_vir'][(objs_pd_comb['Mpeak'] > 3.16e9) & (objs_pd_comb['tau90_vir'] > 0)],  xmax = dmax, nbins = nbins)
    xaxis2, data2 = bin_plt(np.log10((objs_pd_comb['min_dist'][(objs_pd_comb['type'] == 'Central') & (objs_pd_comb['Mpeak'] > 3.16e9) & (objs_pd_comb['tau90_vir'] > 0)]).tolist()), objs_pd_comb['tau90_vir'][(objs_pd_comb['type'] == 'Central') & (objs_pd_comb['Mpeak'] > 3.16e9) & (objs_pd_comb['tau90_vir'] > 0)],  xmax = dmax, nbins = nbins)
    ax11d.fill_between(10**xaxis, data[:,0], data[:,4], color='k', alpha=0.1, linewidth=0)
    ax11d.fill_between(10**xaxis, data[:,1], data[:,3], color='k', alpha=0.1, linewidth=0)
    ax11d.plot(10**xaxis, data[:,2], c='k')
    #ax11d.plot(10**xaxis2, data2[:,2], "--k")
    ax11d.set_xscale('log')
    ax11d.axis([10,4e3,-0.5,14])
    ax11d.set_xlabel(r'Min(D$_{\mathrm{massive}}$) [kpc]')
    ax11b.legend([cen_plt,sat_plt],['Central','Satellite/Backsplash'],scatterpoints = 1,facecolor = 'white',framealpha = 0,frameon = False)
    ax11d.legend([sf_plt,quench_plt],['SF','Quenched'],scatterpoints = 1,facecolor = 'white',framealpha = 0,frameon = False)
    fig11.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.8, wspace=0, hspace=0)
    cb_ax = fig11.add_axes([0.83, 0.1, 0.02, 0.8])
    sm = plt.cm.ScalarMappable(cmap=cmx, norm=cNorm)
    cbar = fig11.colorbar(sm,cax = cb_ax) #pad=0.04)
    cbar.set_label(r'Log(M$_{\mathrm{vir, peak}}$/M$_\odot$)') #, rotation=0)
    fig11.show()
    fig11.savefig(outfile_base + '_rClosest_taus.png',dpi = dpi)

    #############333
    plt.close('all')    
    fig11 = plt.figure(figsize=(plt_width*2,plt_width*aspect_ratio*2))
    gs = fig11.add_gridspec(2,2, hspace=0, wspace=0) #, vspace=0)
    axs11 = gs.subplots(sharex = True, sharey = True) #, constrained_layout=True)
    axs11 = axs11.flatten()
    ax11a = axs11[0]
    ax11b = axs11[1]
    ax11c = axs11[2]
    ax11d = axs11[3]
    cmx = plt.get_cmap("winter_r")
    dcmap = discrete_cmap(256, cmx)
    cNorm  = colors.LogNorm(vmin=1e8, vmax = 1e11)
    dmin = 1
    dmax = 3.2
    nbins = 6
    quench_cond = np.array(list(zip(*objs_pd_comb['SFR']))[0])/objs_pd_comb['M_star'] < 1e-11
    sf_cond =  np.array(list(zip(*objs_pd_comb['SFR']))[0])/objs_pd_comb['M_star'] >= 1e-11    
    ax11a.scatter(objs_pd_comb['min_dist'][quench_cond * (objs_pd_comb['type']=='Central')],objs_pd_comb['tau50'][quench_cond * (objs_pd_comb['type']=='Central')],edgecolor =  dcmap(cNorm(objs_pd_comb['Mpeak'][quench_cond * (objs_pd_comb['type']=='Central')])), cmap = cmx, norm = cNorm, facecolor = 'none', s = objs_pd_comb['p_size'][quench_cond * (objs_pd_comb['type']=='Central')]*1.5, linewidths = edgewidth*2)
    ax11a.scatter(objs_pd_comb['min_dist'][quench_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))],objs_pd_comb['tau50'][quench_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))],edgecolor = dcmap(cNorm(objs_pd_comb['Mpeak'][quench_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))])), cmap = cmx, norm = cNorm,facecolor = 'none', s = objs_pd_comb['p_size'][quench_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))]*1.5*1.5, marker = '*', linewidths = edgewidth*2)    
    ax11a.scatter(objs_pd_comb['min_dist'][quench_cond * (objs_pd_comb['type']=='Satellite')],objs_pd_comb['tau50'][quench_cond * (objs_pd_comb['type']=='Satellite')],edgecolor =  dcmap(cNorm(objs_pd_comb['Mpeak'][quench_cond * (objs_pd_comb['type']=='Satellite')])), cmap = cmx, norm = cNorm,facecolor = 'none', s = objs_pd_comb['p_size'][quench_cond * (objs_pd_comb['type']=='Satellite')]*1.5*1.5, marker = '*', linewidths = edgewidth*2)
    
    xaxis, data = bin_plt(np.log10((objs_pd_comb['min_dist'][(objs_pd_comb['Mpeak'] > 3.16e9) & (objs_pd_comb['tau90_vir'] > 0)]).tolist()), objs_pd_comb['tau50'][(objs_pd_comb['Mpeak'] > 3.16e9) & (objs_pd_comb['tau90_vir'] > 0)], xmax = dmax, nbins = nbins)
    xaxis2, data2 = bin_plt(np.log10((objs_pd_comb['min_dist'][(objs_pd_comb['type'] == 'Central') & (objs_pd_comb['Mpeak'] > 3.16e9) & (objs_pd_comb['tau90_vir'] > 0)]).tolist()), objs_pd_comb['tau50'][(objs_pd_comb['type'] == 'Central') & (objs_pd_comb['Mpeak'] > 3.16e9) & (objs_pd_comb['tau90_vir'] > 0)], xmax = dmax, nbins = nbins)
    ax11a.fill_between(10**xaxis, data[:,0], data[:,4], color='k', alpha=0.1, linewidth=0)
    ax11a.fill_between(10**xaxis, data[:,1], data[:,3], color='k', alpha=0.1, linewidth=0)
    ax11a.plot(10**xaxis, data[:,2], c='k')
    #ax11a.plot(10**xaxis2, data2[:,2], "--k")
    ax11a.set_xscale('log')
    ax11a.set_ylabel(r'$\tau_{50}$')
    ax11a.text(30,12,'Stars')
    ax11a.axis([10, 4e3, -0.5, 14])

    ax11b.scatter(objs_pd_comb['min_dist'][sf_cond * (objs_pd_comb['type']=='Central')],objs_pd_comb['tau50_vir'][sf_cond * (objs_pd_comb['type']=='Central')],c = (objs_pd_comb['Mpeak'][sf_cond * (objs_pd_comb['type']=='Central')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][sf_cond * (objs_pd_comb['type']=='Central')]*1.5, linewidths = edgewidth)
    ax11b.scatter(objs_pd_comb['min_dist'][sf_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))],objs_pd_comb['tau50_vir'][sf_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))],c = (objs_pd_comb['Mpeak'][sf_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][sf_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))]*1.5*1.5, marker = '*', linewidths = edgewidth)    
    ax11b.scatter(objs_pd_comb['min_dist'][sf_cond * (objs_pd_comb['type']=='Satellite')],objs_pd_comb['tau50_vir'][sf_cond * (objs_pd_comb['type']=='Satellite')],c = (objs_pd_comb['Mpeak'][sf_cond * (objs_pd_comb['type']=='Satellite')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][sf_cond * (objs_pd_comb['type']=='Satellite')]*1.5*1.5, marker = '*', linewidths = edgewidth)
    
    ax11b.scatter(objs_pd_comb['min_dist'][quench_cond * (objs_pd_comb['type']=='Central')],objs_pd_comb['tau50_vir'][quench_cond * (objs_pd_comb['type']=='Central')],edgecolor =  dcmap(cNorm(objs_pd_comb['Mpeak'][quench_cond * (objs_pd_comb['type']=='Central')])), cmap = cmx, norm = cNorm, facecolor = 'none', s = objs_pd_comb['p_size'][quench_cond * (objs_pd_comb['type']=='Central')]*1.5, linewidths = edgewidth*2)
    ax11b.scatter(objs_pd_comb['min_dist'][quench_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))],objs_pd_comb['tau50_vir'][quench_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))],edgecolor = dcmap(cNorm(objs_pd_comb['Mpeak'][quench_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))])), cmap = cmx, norm = cNorm,facecolor = 'none', s = objs_pd_comb['p_size'][quench_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))]*1.5*1.5, marker = '*', linewidths = edgewidth*2)    
    ax11b.scatter(objs_pd_comb['min_dist'][quench_cond * (objs_pd_comb['type']=='Satellite')],objs_pd_comb['tau50_vir'][quench_cond * (objs_pd_comb['type']=='Satellite')],edgecolor =  dcmap(cNorm(objs_pd_comb['Mpeak'][quench_cond * (objs_pd_comb['type']=='Satellite')])), cmap = cmx, norm = cNorm,facecolor = 'none', s = objs_pd_comb['p_size'][quench_cond * (objs_pd_comb['type']=='Satellite')]*1.5*1.5, marker = '*', linewidths = edgewidth*2)    

    xaxis, data = bin_plt(np.log10((objs_pd_comb['min_dist'][(objs_pd_comb['Mpeak'] > 3.16e9) & (objs_pd_comb['tau90_vir'] > 0)]).tolist()), objs_pd_comb['tau50_vir'][(objs_pd_comb['Mpeak'] > 3.16e9) & (objs_pd_comb['tau90_vir'] > 0)], xmax = dmax, nbins = nbins)
    xaxis2, data2 = bin_plt(np.log10((objs_pd_comb['min_dist'][(objs_pd_comb['type'] == 'Central') & (objs_pd_comb['Mpeak'] > 3.16e9) & (objs_pd_comb['tau90_vir'] > 0)]).tolist()), objs_pd_comb['tau50_vir'][(objs_pd_comb['type'] == 'Central') & (objs_pd_comb['Mpeak'] > 3.16e9) & (objs_pd_comb['tau90_vir'] > 0)], xmax = dmax, nbins = nbins)
    ax11b.fill_between(10**xaxis, data[:,0], data[:,4], color='k', alpha=0.1, linewidth=0)
    ax11b.fill_between(10**xaxis, data[:,1], data[:,3], color='k', alpha=0.1, linewidth=0)
    ax11b.plot(10**xaxis, data[:,2], c='k')
    #ax11b.plot(10**xaxis2, data2[:,2], "--k")
    ax11b.set_xscale('log')
    ax11b.text(30,12,'Total Mass')
    
    ax11c.scatter(objs_pd_comb['min_dist'][sf_cond * (objs_pd_comb['type']=='Central')],objs_pd_comb['tau90'][sf_cond * (objs_pd_comb['type']=='Central')],c = (objs_pd_comb['Mpeak'][sf_cond * (objs_pd_comb['type']=='Central')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][sf_cond * (objs_pd_comb['type']=='Central')]*1.5, linewidths = edgewidth)
    ax11c.scatter(objs_pd_comb['min_dist'][sf_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))],objs_pd_comb['tau90'][sf_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))],c = (objs_pd_comb['Mpeak'][sf_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][sf_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))]*1.5*1.5, marker = '*', linewidths = edgewidth)    
    ax11c.scatter(objs_pd_comb['min_dist'][sf_cond * (objs_pd_comb['type']=='Satellite')*(objs_pd_comb['Mpeak_snap'] == 13.8)],objs_pd_comb['tau90'][sf_cond * (objs_pd_comb['type']=='Satellite')*(objs_pd_comb['Mpeak_snap'] == 13.8)],c = (objs_pd_comb['Mpeak'][sf_cond * (objs_pd_comb['type']=='Satellite')*(objs_pd_comb['Mpeak_snap'] == 13.8)]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][sf_cond * (objs_pd_comb['type']=='Satellite')*(objs_pd_comb['Mpeak_snap'] == 13.8)]*1.5*1.5, marker = '*', linewidths = edgewidth)
    
    ax11c.scatter(objs_pd_comb['min_dist'][quench_cond * (objs_pd_comb['type']=='Central')*(objs_pd_comb['Mpeak_snap'] == 13.8)],objs_pd_comb['tau90'][quench_cond * (objs_pd_comb['type']=='Central')*(objs_pd_comb['Mpeak_snap'] == 13.8)],edgecolor =  dcmap(cNorm(objs_pd_comb['Mpeak'][quench_cond * (objs_pd_comb['type']=='Central')*(objs_pd_comb['Mpeak_snap'] == 13.8)])), cmap = cmx, norm = cNorm, facecolor = 'none', s = objs_pd_comb['p_size'][quench_cond * (objs_pd_comb['type']=='Central')*(objs_pd_comb['Mpeak_snap'] == 13.8)]*1.5, linewidths = edgewidth*2)
    ax11c.scatter(objs_pd_comb['min_dist'][quench_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))*(objs_pd_comb['Mpeak_snap'] == 13.8)],objs_pd_comb['tau90'][quench_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))*(objs_pd_comb['Mpeak_snap'] == 13.8)],edgecolor = dcmap(cNorm(objs_pd_comb['Mpeak'][quench_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))*(objs_pd_comb['Mpeak_snap'] == 13.8)])), cmap = cmx, norm = cNorm,facecolor = 'none', s = objs_pd_comb['p_size'][quench_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))*(objs_pd_comb['Mpeak_snap'] == 13.8)]*1.5*1.5, marker = '*', linewidths = edgewidth*2)    
    ax11c.scatter(objs_pd_comb['min_dist'][quench_cond * (objs_pd_comb['type']=='Satellite')*(objs_pd_comb['Mpeak_snap'] == 13.8)],objs_pd_comb['tau90'][quench_cond * (objs_pd_comb['type']=='Satellite')*(objs_pd_comb['Mpeak_snap'] == 13.8)],edgecolor =  dcmap(cNorm(objs_pd_comb['Mpeak'][quench_cond * (objs_pd_comb['type']=='Satellite')*(objs_pd_comb['Mpeak_snap'] == 13.8)])), cmap = cmx, norm = cNorm,facecolor = 'none', s = objs_pd_comb['p_size'][quench_cond * (objs_pd_comb['type']=='Satellite')*(objs_pd_comb['Mpeak_snap'] == 13.8)]*1.5*1.5, marker = '*', linewidths = edgewidth*2)
    
    xaxis, data = bin_plt(np.log10((objs_pd_comb['min_dist'][(objs_pd_comb['Mpeak'] > 3.16e9) & (objs_pd_comb['tau90_vir'] > 0)]).tolist()), objs_pd_comb['tau90'][(objs_pd_comb['Mpeak'] > 3.16e9) & (objs_pd_comb['tau90_vir'] > 0)], xmax = dmax, nbins = nbins)
    xaxis2, data2 = bin_plt(np.log10((objs_pd_comb['min_dist'][(objs_pd_comb['type'] == 'Central') & (objs_pd_comb['Mpeak'] > 3.16e9) & (objs_pd_comb['tau90_vir'] > 0)]).tolist()), objs_pd_comb['tau90'][(objs_pd_comb['type'] == 'Central') & (objs_pd_comb['Mpeak'] > 3.16e9) & (objs_pd_comb['tau90_vir'] > 0)], xmax = dmax, nbins = nbins)    
    ax11c.fill_between(10**xaxis, data[:,0], data[:,4], color='k', alpha=0.1, linewidth=0)
    ax11c.fill_between(10**xaxis, data[:,1], data[:,3], color='k', alpha=0.1, linewidth=0)
    ax11c.plot(10**xaxis, data[:,2], c='k')
    #ax11c.plot(10**xaxis2, data2[:,2], "--k")
    ax11c.set_xscale('log')
    ax11c.set_ylabel(r'$\tau_{90}$')
    ax11c.set_xlabel(r'Min(D$_{\mathrm{massive}}$) [kpc]')

    cen_plt = ax11d.scatter(objs_pd_comb['min_dist'][sf_cond * (objs_pd_comb['type']=='Central')],objs_pd_comb['tau90_vir'][sf_cond * (objs_pd_comb['type']=='Central')],c = (objs_pd_comb['Mpeak'][sf_cond * (objs_pd_comb['type']=='Central')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][sf_cond * (objs_pd_comb['type']=='Central')]*1.5, linewidths = edgewidth)
    sf_plt = ax11d.scatter(objs_pd_comb['min_dist'][sf_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))],objs_pd_comb['tau90_vir'][sf_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))],c = (objs_pd_comb['Mpeak'][sf_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][sf_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))]*1.5*1.5, marker = '*', linewidths = edgewidth)    
    sat_plt = ax11d.scatter(objs_pd_comb['min_dist'][sf_cond * (objs_pd_comb['type']=='Satellite')],objs_pd_comb['tau90_vir'][sf_cond * (objs_pd_comb['type']=='Satellite')],c = (objs_pd_comb['Mpeak'][sf_cond * (objs_pd_comb['type']=='Satellite')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][sf_cond * (objs_pd_comb['type']=='Satellite')]*1.5*1.5, marker = '*', linewidths = edgewidth)
    
    ax11d.scatter(objs_pd_comb['min_dist'][quench_cond * (objs_pd_comb['type']=='Central')],objs_pd_comb['tau90_vir'][quench_cond * (objs_pd_comb['type']=='Central')],edgecolor =  dcmap(cNorm(objs_pd_comb['Mpeak'][quench_cond * (objs_pd_comb['type']=='Central')])), cmap = cmx, norm = cNorm, facecolor = 'none', s = objs_pd_comb['p_size'][quench_cond * (objs_pd_comb['type']=='Central')]*1.5, linewidths = edgewidth*2)
    quench_plt = ax11d.scatter(objs_pd_comb['min_dist'][quench_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))],objs_pd_comb['tau90_vir'][quench_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))],edgecolor = dcmap(cNorm(objs_pd_comb['Mpeak'][quench_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))])), cmap = cmx, norm = cNorm,facecolor = 'none', s = objs_pd_comb['p_size'][quench_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))]*1.5*1.5, marker = '*', linewidths = edgewidth*2)    
    ax11d.scatter(objs_pd_comb['min_dist'][quench_cond * (objs_pd_comb['type']=='Satellite')],objs_pd_comb['tau90_vir'][quench_cond * (objs_pd_comb['type']=='Satellite')],edgecolor =  dcmap(cNorm(objs_pd_comb['Mpeak'][quench_cond * (objs_pd_comb['type']=='Satellite')])), cmap = cmx, norm = cNorm,facecolor = 'none', s = objs_pd_comb['p_size'][quench_cond * (objs_pd_comb['type']=='Satellite')]*1.5*1.5, marker = '*', linewidths = edgewidth*2)    
    
    xaxis, data = bin_plt(np.log10((objs_pd_comb['min_dist'][(objs_pd_comb['Mpeak'] > 3.16e9) & (objs_pd_comb['tau90_vir'] > 0)]).tolist()), objs_pd_comb['tau90_vir'][(objs_pd_comb['Mpeak'] > 3.16e9) & (objs_pd_comb['tau90_vir'] > 0)],  xmax = dmax, nbins = nbins)
    xaxis2, data2 = bin_plt(np.log10((objs_pd_comb['min_dist'][(objs_pd_comb['type'] == 'Central') & (objs_pd_comb['Mpeak'] > 3.16e9) & (objs_pd_comb['tau90_vir'] > 0)]).tolist()), objs_pd_comb['tau90_vir'][(objs_pd_comb['type'] == 'Central') & (objs_pd_comb['Mpeak'] > 3.16e9) & (objs_pd_comb['tau90_vir'] > 0)],  xmax = dmax, nbins = nbins)
    ax11d.fill_between(10**xaxis, data[:,0], data[:,4], color='k', alpha=0.1, linewidth=0)
    ax11d.fill_between(10**xaxis, data[:,1], data[:,3], color='k', alpha=0.1, linewidth=0)
    ax11d.plot(10**xaxis, data[:,2], c='k')
    #ax11d.plot(10**xaxis2, data2[:,2], "--k")
    ax11d.set_xscale('log')
    ax11d.axis([10,4e3,-0.5,14])
    ax11d.set_xlabel(r'Min(D$_{\mathrm{massive}}$) [kpc]')
    ax11b.legend([cen_plt,sat_plt],['Central','Satellite/Backsplash'],scatterpoints = 1,facecolor = 'white',framealpha = 0,frameon = False)
    ax11d.legend([sf_plt,quench_plt],['SF','Quenched'],scatterpoints = 1,facecolor = 'white',framealpha = 0,frameon = False)
    fig11.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.8, wspace=0, hspace=0)
    cb_ax = fig11.add_axes([0.83, 0.1, 0.02, 0.8])
    sm = plt.cm.ScalarMappable(cmap=cmx, norm=cNorm)
    cbar = fig11.colorbar(sm,cax = cb_ax) #pad=0.04)
    cbar.set_label(r'Log(M$_{\mathrm{vir, peak}}$/M$_\odot$)') #, rotation=0)
    fig11.show()
    fig11.savefig(outfile_base + '_rClosest_taus_nonpeak.png',dpi = dpi)
    
    plt.close('all')
    fig11 = plt.figure(11,figsize=(plt_width*2,plt_width*aspect_ratio))
    fig11.set_size_inches(plt_width*2,plt_width*aspect_ratio)
    fig11.clear()
    gs = fig11.add_gridspec(1,2,wspace=0)
    axs11 = gs.subplots(sharey = True) #, constrained_layout=True) 
    axs11 = axs11.flatten()
    ax11a = axs11[0]
    ax11b = axs11[1]
    cmx = plt.get_cmap("viridis")
    cNorm  = colors.LogNorm(vmin=10, vmax = 2e3)

    xaxis, data = bin_plt((np.log10(objs_pd_comb['Mpeak'])[np.isfinite(objs_pd_comb['tau50']/objs_pd_comb['tau50_vir'])]), (objs_pd_comb['tau50']/objs_pd_comb['tau50_vir'])[np.isfinite(objs_pd_comb['tau50']/objs_pd_comb['tau50_vir'])], xmin = np.log10(1e8), xmax = np.log10(2e11), nbins = nbins)    
    ax11a.plot([1e8,1e12],[1,1], color = 'k', linestyle = 'dashed')
    ax11a.fill_between(10**xaxis, data[:,0], data[:,4], color='k', alpha=0.1, linewidth=0, zorder=1)
    ax11a.fill_between(10**xaxis, data[:,1], data[:,3], color='k', alpha=0.1, linewidth=0, zorder=1)
    tau50 = ax11a.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],objs_pd_comb['tau50'][objs_pd_comb['type']=='Central']/objs_pd_comb['tau50_vir'][objs_pd_comb['type']=='Central'],c = (objs_pd_comb['min_dist'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth, zorder=2)
    ax11a.scatter(objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb['tau50'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]/objs_pd_comb['tau50_vir'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = (objs_pd_comb['min_dist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], linewidths = edgewidth, zorder=2)    
    ax11a.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],objs_pd_comb['tau50'][objs_pd_comb['type']=='Satellite']/objs_pd_comb['tau50_vir'][objs_pd_comb['type']=='Satellite'],c =(objs_pd_comb['min_dist'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k',s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth, zorder=2, marker = '*')
    ax11a.plot(10**xaxis, data[:,2], c='k', zorder=2)    
    ax11a.text(1e10,1e-1,r'$\tau$ = 50')    
    ax11a.set_ylabel(r'$\tau_{t}/\tau_{\mathrm{vir, t}}$')
    ax11a.set_xlabel(r'M$_{\mathrm{vir, peak}}$ [M$_\odot$]')
    ax11a.set_xscale('log')
    ax11a.set_yscale('log')
    #ax11a.axis([1e8, 2e11,-0.5, 4])
    ax11a.axis([1e8, 2e11, 3e-2, 4])
    
    xaxis, data = bin_plt((np.log10(objs_pd_comb['Mpeak']))[np.isfinite(objs_pd_comb['tau90']/objs_pd_comb['tau90_vir'])], (objs_pd_comb['tau90']/objs_pd_comb['tau90_vir'])[np.isfinite(objs_pd_comb['tau90']/objs_pd_comb['tau90_vir'])], xmin = np.log10(1e8), xmax = np.log10(2e11), nbins = nbins)
    ax11b.plot([1e8,1e12],[1,1], color = 'k', linestyle = 'dashed')
    ax11b.fill_between(10**xaxis, data[:,0], data[:,4], color='k', alpha=0.1, linewidth=0, zorder=1)
    ax11b.fill_between(10**xaxis, data[:,1], data[:,3], color='k', alpha=0.1, linewidth=0, zorder=1)
    cen_plt = ax11b.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],objs_pd_comb['tau90'][objs_pd_comb['type']=='Central']/objs_pd_comb['tau90_vir'][objs_pd_comb['type']=='Central'],c = (objs_pd_comb['min_dist'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth, zorder=2)
    ax11b.scatter(objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb['tau90'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]/objs_pd_comb['tau90_vir'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = (objs_pd_comb['min_dist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], linewidths = edgewidth, zorder=2)    
    sat_plt = ax11b.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],objs_pd_comb['tau90'][objs_pd_comb['type']=='Satellite']/objs_pd_comb['tau90_vir'][objs_pd_comb['type']=='Satellite'],c = (objs_pd_comb['min_dist'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth, zorder=2, marker = '*')
    ax11b.plot(10**xaxis, data[:,2], c='k', zorder=2)
    ax11b.text(1e10,1e-1,r'$\tau$ = 90')
    #ax11b.set_ylabel(r'$\tau_{t}/\tau_{\mathrm{vir, t}}$')
    ax11b.set_xlabel(r'M$_{\mathrm{vir, peak}}$ [M$_\odot$]')
    ax11b.set_xscale('log')
    ax11b.set_yscale('log')
    #ax11n.axis([1e8, 2e11,-0.5, 4])
    ax11b.axis([1e8, 2e11, 3e-2, 4])
    ax11b.legend([cen_plt,sat_plt],['Central','Satellite/Backsplash'],scatterpoints = 1,facecolor = 'white',loc = 0,framealpha = 0,frameon = False)
    fig11.tight_layout()

    #legend1 = ax11.legend([tau50,tau90],['t = 50','t = 90'],scatterpoints = 1,facecolor = 'white',loc = 4,framealpha = 0,frameon = False)
    #ax11.add_artist(legend1)
    #cb = mpl.colorbar.ColorbarBase(ax11sub, cmap=cmx, norm=cNorm)
    sm = plt.cm.ScalarMappable(cmap=cmx, norm=cNorm)
    cbar = plt.colorbar(sm, location="right",ax=[ax11a,ax11b],aspect = 20)    
    cbar.set_label(r"Min(D$_{\mathrm{massive}}$) [kpc]")
    fig11.show()    
    fig11.savefig(outfile_base + '_rClosest_tau90_tau50_Mpeak.png',dpi = dpi)

    plt.close('all')
    fig11 = plt.figure(11,figsize=(plt_width*2,plt_width*aspect_ratio))
    fig11.set_size_inches(plt_width*2,plt_width*aspect_ratio)
    fig11.clear()
    gs = fig11.add_gridspec(1,2,wspace=0)
    axs11 = gs.subplots(sharey = True) #, constrained_layout=True) 
    axs11 = axs11.flatten()
    ax11a = axs11[0]
    ax11b = axs11[1]    
    cNorm  = colors.LogNorm(vmin=10**1.5, vmax = 10**4)

    xaxis, data = bin_plt((np.log10(objs_pd_comb['Mpeak'])[np.isfinite(objs_pd_comb['tau50']/objs_pd_comb['tau50_vir'])]), (objs_pd_comb['tau50']/objs_pd_comb['tau50_vir'])[np.isfinite(objs_pd_comb['tau50']/objs_pd_comb['tau50_vir'])], xmin = np.log10(1e8), xmax = np.log10(2e11), nbins = nbins)    
    ax11a.plot([1e8,1e12],[1,1], color = 'k', linestyle = 'dashed')
    ax11a.fill_between(10**xaxis, data[:,0], data[:,4], color='k', alpha=0.1, linewidth=0, zorder=1)
    ax11a.fill_between(10**xaxis, data[:,1], data[:,3], color='k', alpha=0.1, linewidth=0, zorder=1)
    tau50 = ax11a.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],objs_pd_comb['tau50'][objs_pd_comb['type']=='Central']/objs_pd_comb['tau50_vir'][objs_pd_comb['type']=='Central'],c = (objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth, zorder=2)
    ax11a.scatter(objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb['tau50'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]/objs_pd_comb['tau50_vir'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = (objs_pd_comb['massiveDist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], linewidths = edgewidth, zorder=2)    
    ax11a.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],objs_pd_comb['tau50'][objs_pd_comb['type']=='Satellite']/objs_pd_comb['tau50_vir'][objs_pd_comb['type']=='Satellite'],c =(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k',s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth, zorder=2, marker = '*')
    ax11a.plot(10**xaxis, data[:,2], c='k', zorder=2)    
    ax11a.text(1e10,1e-1,r'$\tau$ = 50')    
    ax11a.set_ylabel(r'$\tau_{t}/\tau_{\mathrm{vir, t}}$')
    ax11a.set_xlabel(r'M$_{\mathrm{vir, peak}}$ [M$_\odot$]')
    ax11a.set_xscale('log')
    ax11a.set_yscale('log')
    #ax11a.axis([1e8, 2e11,-0.5, 4])
    ax11a.axis([1e8, 2e11, 3e-2, 4])
    
    xaxis, data = bin_plt((np.log10(objs_pd_comb['Mpeak']))[np.isfinite(objs_pd_comb['tau90']/objs_pd_comb['tau90_vir'])], (objs_pd_comb['tau90']/objs_pd_comb['tau90_vir'])[np.isfinite(objs_pd_comb['tau90']/objs_pd_comb['tau90_vir'])], xmin = np.log10(1e8), xmax = np.log10(2e11), nbins = nbins)
    ax11b.plot([1e8,1e12],[1,1], color = 'k', linestyle = 'dashed')
    ax11b.fill_between(10**xaxis, data[:,0], data[:,4], color='k', alpha=0.1, linewidth=0, zorder=1)
    ax11b.fill_between(10**xaxis, data[:,1], data[:,3], color='k', alpha=0.1, linewidth=0, zorder=1)
    cen_plt = ax11b.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],objs_pd_comb['tau90'][objs_pd_comb['type']=='Central']/objs_pd_comb['tau90_vir'][objs_pd_comb['type']=='Central'],c = (objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth, zorder=2)
    ax11b.scatter(objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb['tau90'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]/objs_pd_comb['tau90_vir'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = (objs_pd_comb['massiveDist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], linewidths = edgewidth, zorder=2)    
    sat_plt = ax11b.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],objs_pd_comb['tau90'][objs_pd_comb['type']=='Satellite']/objs_pd_comb['tau90_vir'][objs_pd_comb['type']=='Satellite'],c = (objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth, zorder=2, marker = '*')
    ax11b.plot(10**xaxis, data[:,2], c='k', zorder=2)
    ax11b.text(1e10,1e-1,r'$\tau$ = 90')
    #ax11b.set_ylabel(r'$\tau_{t}/\tau_{\mathrm{vir, t}}$')
    ax11b.set_xlabel(r'M$_{\mathrm{vir, peak}}$ [M$_\odot$]')
    ax11b.set_xscale('log')
    ax11b.set_yscale('log')
    #ax11n.axis([1e8, 2e11,-0.5, 4])
    ax11b.axis([1e8, 2e11, 3e-2, 4])
    ax11b.legend([cen_plt,sat_plt],['Central','Satellite/Backsplash'],scatterpoints = 1,facecolor = 'white',loc = 0,framealpha = 0,frameon = False)
    fig11.tight_layout()

    #legend1 = ax11.legend([tau50,tau90],['t = 50','t = 90'],scatterpoints = 1,facecolor = 'white',loc = 4,framealpha = 0,frameon = False)
    #ax11.add_artist(legend1)
    #cb = mpl.colorbar.ColorbarBase(ax11sub, cmap=cmx, norm=cNorm)
    sm = plt.cm.ScalarMappable(cmap=cmx, norm=cNorm)
    cbar = plt.colorbar(sm, location="right",ax=[ax11a,ax11b],aspect = 20)    
    cbar.set_label(r"D$_{\mathrm{massive}}$ [kpc]")
    fig11.show()    
    fig11.savefig(outfile_base + '_rMassGal_tau90_tau50_Mpeak.png',dpi = dpi)
    
    plt.close('all')
    fig11 = plt.figure(2,figsize=(plt_width,plt_width*aspect_ratio))
    fig11.set_size_inches(plt_width,plt_width*aspect_ratio)
    fig11.clear()
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax11 = fig11.add_subplot(gs[0])
    ax11sub = fig11.add_subplot(gs[1])
    cNorm  = colors.LogNorm(vmin=10**1.5, vmax = 10**4)
    tau90 = ax11.scatter(objs_pd_comb['mstar'][objs_pd_comb['type']=='Central'],objs_pd_comb['tau90'][objs_pd_comb['type']=='Central'],c = (objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax11.scatter(objs_pd_comb['mstar'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb['tau90'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = (objs_pd_comb['massiveDist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], linewidths = edgewidth)    
    ax11.scatter(objs_pd_comb['mstar'][objs_pd_comb['type']=='Satellite'],objs_pd_comb['tau90'][objs_pd_comb['type']=='Satellite'],c = (objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite'], linewidths = edgewidth, marker = '*')
    ax11.scatter(digbyfield['Mstar'][digbyfield['oMSTO']=='y'],np.array((digbyfield['t90'][digbyfield['oMSTO']=='y']).tolist()).astype('float'), c = 'k', s = markersize, facecolor = 'None',marker = 'o')
    ax11.scatter(digbysat['Mstar'][digbysat['oMSTO']=='y'],np.array((digbysat['t90'][digbysat['oMSTO']=='y']).tolist()).astype('float'), c = 'k', s = markersize, facecolor = 'None',marker = '*', ls = 'None') 
    #ax11.errorbar(digbysat['Mstar'][digbysat['oMSTO']=='y'],digbysat['t90'][digbysat['oMSTO']=='y'], yerr = np.array([digbysat['t90_lerr'][digbysat['oMSTO']=='y'].tolist(),digbysat['t90_herr'][digbysat['oMSTO']=='y'].tolist()]), c = 'k', markersize = markersize*0.2, alpha = 0.8, marker = '*', ls = 'None')
    #ax11.errorbar(digbyfield['Mstar'][digbyfield['oMSTO']=='y'],digbyfield['t90'][digbyfield['oMSTO']=='y'], yerr = np.array([digbyfield['t90_lerr'][digbyfield['oMSTO']=='y'].tolist(),digbyfield['t90_herr'][digbyfield['oMSTO']=='y'].tolist()]), c = 'k', markersize = markersize*0.2, alpha = 0.8, marker = 'o', ls = 'None')
    
    cb = mpl.colorbar.ColorbarBase(ax11sub, cmap=cmx, norm=cNorm)
    cb.set_label(r"D$_{\mathrm{massive}}$ [kpc]")
    ax11.set_ylabel(r'$\tau_{90}$')
    ax11.set_xlabel(r'M$_{\mathrm{vir, peak}}$ [M$_\odot$]')
    ax11.set_xscale('log')    
    fig11.tight_layout()
    fig11.show()
    fig11.savefig(outfile_base + '_tau90_Mpeak_obs.png',dpi = dpi)


        #SMHM colored by distance to massive galaxy
    #plt.clf()
    cmx = plt.get_cmap("viridis") 
    fig12 = plt.figure(12,figsize=(plt_width,plt_width*aspect_ratio))
    fig12.set_size_inches(plt_width,plt_width*aspect_ratio)
    fig12.clear()
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax12 = fig12.add_subplot(gs[0])
    ax12sub = fig12.add_subplot(gs[1])
    cNorm  = colors.LogNorm(vmin=10**1.5, vmax = 10**4)
    ax12.fill_between(read_abunmatch['M200'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_low'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_high'][read_abunmatch['M200']> 1e9],color = 'grey',linewidth=0,alpha=0.6)
    ax12.fill_between(upper_fid_Mhalo, np.interp( upper_fid_Mhalo, lower_fid_Mhalo, lower_fid_Mstar ), upper_fid_Mstar, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3)#Nadler outer
    ax12.fill_between(upper_fid_inner_Mhalo, np.interp( upper_fid_inner_Mhalo, lower_fid_inner_Mhalo, lower_fid_inner_Mstar ),upper_fid_inner_Mstar, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3) #Nadler inner
    ax12.fill_between(Jethwa_upper_outer_Mhalo, np.interp( Jethwa_upper_outer_Mhalo, Jethwa_lower_outer_Mhalo, Jethwa_lower_outer_Mstar ),Jethwa_upper_outer_Mstar, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3)#Jethwa outer
    ax12.fill_between(Jethwa_upper_inner_Mhalo,np.interp( Jethwa_upper_inner_Mhalo, Jethwa_lower_inner_Mhalo, Jethwa_lower_inner_Mstar ), Jethwa_upper_inner_Mstar, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3) #Jethwa inner
    read, = ax12.plot( [], [],  color="grey", linestyle="solid", linewidth=4) #, label="Read+ 2017" )
    nadler, = ax12.plot( [], [],  color="#33669A", linestyle="solid", linewidth=8) #, label="Nadler+ 2020" )
    jethwa, = ax12.plot( [], [],  color="#66CC66", linestyle="solid", linewidth=8) #, label="Jethwa+ 2018" )     
    cen_plt = ax12.scatter(objs_pd_comb['mass'][objs_pd_comb['type']=='Central'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Central'],c = objs_pd_comb['min_dist'][objs_pd_comb['type']=='Central'], cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax12.scatter(objs_pd_comb['mass'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb[mstar_key][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = objs_pd_comb['min_dist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]*1.5,marker = '*', linewidths = edgewidth)    
    sat_plt = ax12.scatter(objs_pd_comb['mass'][objs_pd_comb['type']=='Satellite'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Satellite'],c = objs_pd_comb['min_dist'][objs_pd_comb['type']=='Satellite'], cmap = cmx, norm = cNorm,edgecolor = 'k',marker = '*', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    #ax12.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax12.set_xscale('log')
    ax12.set_yscale('log') 
    ax12.set_ylabel(r'M$_*$/M$_\odot$')
    ax12.set_xlabel(r'M$_{\mathrm{vir}}$/M$_\odot$')
    ax12.axis([2e6, 2e11, 2e2, 5e9])
    legend1 = ax12.legend([cen_plt,sat_plt],['Central','Satellite'],scatterpoints = 1,facecolor = 'white',loc = 4,framealpha = 0,frameon = False)
    legend2 = ax12.legend([read,nadler,jethwa],['Read+ 2017','Nadler+ 2020','Jethwa+ 2018'],facecolor = 'white',loc = 0,framealpha = 0,frameon = False)
    ax12.add_artist(legend1)
    ax12.add_artist(legend2)    
    cb = mpl.colorbar.ColorbarBase(ax12sub, cmap=cmx, norm=cNorm)
    cb.set_label(r"Min(D$_{\mathrm{massive}}$) [kpc]")
    fig12.tight_layout()
    fig12.show()    
    fig12.savefig(outfile_base + '_SMHM_rClosestDist.png',dpi = dpi)

    #SMHM colored by distance to massive galaxy
    #plt.clf()
    #fig12.clear()
    plt.close('all')
    cNorm  = colors.LogNorm(vmin=10, vmax = 2e3)
    fig12 = plt.figure(12,figsize=(plt_width,plt_width*aspect_ratio))
    fig12.set_size_inches(plt_width,plt_width*aspect_ratio)
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax12 = fig12.add_subplot(gs[0])
    ax12sub = fig12.add_subplot(gs[1])
    ax12.fill_between(read_abunmatch['M200'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_low'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_high'][read_abunmatch['M200']> 1e9],color = 'grey',linewidth=0,alpha=0.6)
    ax12.fill_between(upper_fid_Mhalo, np.interp( upper_fid_Mhalo, lower_fid_Mhalo, lower_fid_Mstar ), upper_fid_Mstar, facecolor='y',linewidth=0,zorder=0,alpha=0.3)#Nadler outer  facecolor="#33669A"
    ax12.fill_between(upper_fid_inner_Mhalo, np.interp( upper_fid_inner_Mhalo, lower_fid_inner_Mhalo, lower_fid_inner_Mstar ),upper_fid_inner_Mstar, facecolor='y',linewidth=0,zorder=0,alpha=0.3) #Nadler inner  facecolor="#33669A"
    ax12.fill_between(Jethwa_upper_outer_Mhalo, np.interp( Jethwa_upper_outer_Mhalo, Jethwa_lower_outer_Mhalo, Jethwa_lower_outer_Mstar ),Jethwa_upper_outer_Mstar, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3)#Jethwa outer
    ax12.fill_between(Jethwa_upper_inner_Mhalo,np.interp( Jethwa_upper_inner_Mhalo, Jethwa_lower_inner_Mhalo, Jethwa_lower_inner_Mstar ), Jethwa_upper_inner_Mstar, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3) #Jethwa inner
    read, = ax12.plot( [], [],  color="grey", linestyle="solid", linewidth=4) #, label="Read+ 2017" )
    nadler, = ax12.plot( [], [],  color="#33669A", linestyle="solid", linewidth=8) #, label="Nadler+ 2020" )
    jethwa, = ax12.plot( [], [],  color="#66CC66", linestyle="solid", linewidth=8) #, label="Jethwa+ 2018" )    
    #ax12.scatter(objs_rom['Mpeak'],objs_rom['mstar'],c = objs_rom['massiveDist'],cmap = cmx, norm = cNorm,alpha = 0.4, s = markersize*0.5)
    ax12.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Central'],c = objs_pd_comb['min_dist'][objs_pd_comb['type']=='Central'], cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax12.scatter(objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb[mstar_key][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = objs_pd_comb['min_dist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]*1.5, marker = "*", linewidths = edgewidth)     
    ax12.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Satellite'],c = objs_pd_comb['min_dist'][objs_pd_comb['type']=='Satellite'], cmap = cmx, norm = cNorm,edgecolor = 'k', marker = "*", s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    #ax12.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax12.set_xscale('log')
    ax12.set_yscale('log')
    ax12.set_ylabel(mstar_axes_label)
    ax12.set_xlabel(r'M$_{\mathrm{vir, peak}}$/M$_\odot$')
    ax12.axis([1e8, 2e11, 2e2, 5e9]) #ax12.axis([1e8, 2e13, 2e2, 5e9])
    legend1 = ax12.legend([cen_plt,sat_plt],['Central','Satellite'],scatterpoints = 1,facecolor = 'white',loc = 0,framealpha = 0,frameon = False)
    legend2 = ax12.legend([read,nadler,jethwa],['Read+ 2017','Nadler+ 2020','Jethwa+ 2018'],facecolor = 'white',loc = 4,framealpha = 0,frameon = False)
    ax12.add_artist(legend1)
    ax12.add_artist(legend2)      
    cb = mpl.colorbar.ColorbarBase(ax12sub, cmap=cmx, norm=cNorm)
    cb.set_label(r"Min(D$_{\mathrm{massive}}$) [kpc]")
    fig12.tight_layout()
    fig12.show() 
    fig12.savefig(outfile_base + '_SMHM_rClosestDist_Mpeak.png',dpi = dpi)

    """
    plt.close('all')
    fig12 = plt.figure(2,figsize=(plt_width,plt_width*aspect_ratio))
    fig12.set_size_inches(plt_width,plt_width*aspect_ratio)
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax12 = fig12.add_subplot(gs[0])
    ax12sub = fig12.add_subplot(gs[1])
    sm = plt.cm.ScalarMappable(cmap=cmx, norm=cNorm)
    ax12.set_xscale('log')
    ax12.set_yscale('log')
    ax12.set_ylabel(r'M$_{*, peak}$/M$_\odot$')
    ax12.set_xlabel(r'M$_{\mathrm{vir, peak}}$/M$_\odot$')
    #ax12.axis([2e9, 2e11, 3e3, 1e10])
    ax12.axis([1e8, 2e11, 2e2, 5e9])
    ax12.fill_between(10**((mvir_bins[0:-1] + mvir_bins[1:])/2),mstar_bins0[0,:],mstar_bins0[4,:],color = sm.to_rgba((dist_binsc[0])),alpha = 0.3,linewidth=0)
    ax12.fill_between(10**((mvir_bins[0:-1] + mvir_bins[1:])/2),mstar_bins0[1,:],mstar_bins0[3,:],color = sm.to_rgba((dist_binsc[0])),alpha = 0.3,linewidth=0)    
    ax12.plot(10**((mvir_bins[0:-1] + mvir_bins[1:])/2),mstar_bins0[2,:],c = sm.to_rgba((dist_binsc[0])))

    ax12.fill_between(10**((mvir_bins[0:-1] + mvir_bins[1:])/2),mstar_bins1[0,:],mstar_bins1[4,:],color = sm.to_rgba((dist_binsc[1])),alpha = 0.3,linewidth=0)
    ax12.fill_between(10**((mvir_bins[0:-1] + mvir_bins[1:])/2),mstar_bins1[1,:],mstar_bins1[3,:],color = sm.to_rgba((dist_binsc[1])),alpha = 0.3,linewidth=0)    
    ax12.plot(10**((mvir_bins[0:-1] + mvir_bins[1:])/2),mstar_bins1[2,:],c = sm.to_rgba((dist_binsc[1])))

    ax12.fill_between(10**((mvir_bins[0:-1] + mvir_bins[1:])/2),mstar_bins2[0,:],mstar_bins2[4,:],color = sm.to_rgba((dist_binsc[2])),alpha = 0.3,linewidth=0)
    ax12.fill_between(10**((mvir_bins[0:-1] + mvir_bins[1:])/2),mstar_bins2[1,:],mstar_bins2[3,:],color = sm.to_rgba((dist_binsc[2])),alpha = 0.3,linewidth=0)    
    ax12.plot(10**((mvir_bins[0:-1] + mvir_bins[1:])/2),mstar_bins2[2,:],c = sm.to_rgba((dist_binsc[2])))    
    
    ax12.fill_between(10**((mvir_bins[0:-1] + mvir_bins[1:])/2),mstar_bins3[0,:],mstar_bins3[4,:],color = sm.to_rgba((dist_binsc[3])),alpha = 0.3,linewidth=0)
    ax12.fill_between(10**((mvir_bins[0:-1] + mvir_bins[1:])/2),mstar_bins3[1,:],mstar_bins3[3,:],color = sm.to_rgba((dist_binsc[3])),alpha = 0.3,linewidth=0)    
    ax12.plot(10**((mvir_bins[0:-1] + mvir_bins[1:])/2),mstar_bins3[2,:],c = sm.to_rgba((dist_binsc[3])))
    ax12.scatter(objs_pd_comb['Mpeak'],objs_pd_comb[mstar_key],c = objs_pd_comb['min_dist'], cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'], linewidths = edgewidth)
    cb = mpl.colorbar.ColorbarBase(ax12sub, cmap=cmx, norm=cNorm)
    cb.set_label(r"Min(D$_{\mathrm{massive}}$) [kpc]")
    fig12.tight_layout()    
    fig12.show()
    fig12.savefig(outfile_base + '_SMHM_rClosestDist_Mpeak_rom.png',dpi = dpi)
    """
    plt.close('all')
    #fig12.clear()
    fig12 = plt.figure(12,figsize=(plt_width,plt_width*aspect_ratio))
    fig12.set_size_inches(plt_width,plt_width*aspect_ratio)
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax12 = fig12.add_subplot(gs[0])
    ax12sub = fig12.add_subplot(gs[1])
    ax12.fill_between(read_abunmatch['M200'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_low'][read_abunmatch['M200']> 1e9]/read_abunmatch['M200'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_high'][read_abunmatch['M200']> 1e9]/read_abunmatch['M200'][read_abunmatch['M200']> 1e9],color = 'grey',linewidth=0,alpha = 0.6)
    ax12.fill_between(upper_fid_Mhalo, np.interp( upper_fid_Mhalo, lower_fid_Mhalo, lower_fid_Mstar )/upper_fid_Mhalo, upper_fid_Mstar/upper_fid_Mhalo, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3)#Nadler outer
    ax12.fill_between(upper_fid_inner_Mhalo, np.interp( upper_fid_inner_Mhalo, lower_fid_inner_Mhalo, lower_fid_inner_Mstar)/upper_fid_inner_Mhalo,upper_fid_inner_Mstar/upper_fid_inner_Mhalo, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3) #Nadler inner
    ax12.fill_between(Jethwa_upper_outer_Mhalo, np.interp( Jethwa_upper_outer_Mhalo, Jethwa_lower_outer_Mhalo, Jethwa_lower_outer_Mstar )/Jethwa_upper_outer_Mhalo,Jethwa_upper_outer_Mstar/Jethwa_upper_outer_Mhalo, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3)#Jethwa outer
    ax12.fill_between(Jethwa_upper_inner_Mhalo,np.interp( Jethwa_upper_inner_Mhalo, Jethwa_lower_inner_Mhalo, Jethwa_lower_inner_Mstar )/Jethwa_upper_inner_Mhalo, Jethwa_upper_inner_Mstar/Jethwa_upper_inner_Mhalo, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3) #Jethwa inner
    read, = ax12.plot( [], [],  color="grey", linestyle="solid", linewidth=4) #, label="Read+ 2017" )
    nadler, = ax12.plot( [], [],  color="#33669A", linestyle="solid", linewidth=8) #, label="Nadler+ 2020" )
    jethwa, = ax12.plot( [], [],  color="#66CC66", linestyle="solid", linewidth=8) #, label="Jethwa+ 2018" )        
    ax12.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Central']/objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],c = (objs_pd_comb['min_dist'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax12.scatter(objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb[mstar_key][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]/objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = (objs_pd_comb['min_dist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]*1.5, marker = "*", linewidths = edgewidth)
    ax12.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Satellite']/objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],c = (objs_pd_comb['min_dist'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k', marker = "*", s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    #ax12.scatter(readdata['M200']*1e10,readdata['M_*']*1e7/(readdata['M200']*1e10),marker = '+',c = 'k')
    ax12.set_xscale('log')
    ax12.set_yscale('log')
    ax12.set_ylabel(mstar_ratio_axes_label)
    ax12.set_xlabel(r'M$_{\mathrm{vir, peak}}$/M$_\odot$')
    ax12.axis([1e8, 2e11, 1e-6, 0.05])
    legend1 = ax12.legend([cen_plt,sat_plt],['Central','Satellite/Backsplash'],scatterpoints = 1,facecolor = 'white',loc = 4,framealpha = 0,frameon = False)
    legend2 = ax12.legend([read,nadler,jethwa],['Read+ 2017','Nadler+ 2020','Jethwa+ 2018'],facecolor = 'white',loc = 0,framealpha = 0,frameon = False)
    ax12.add_artist(legend1)
    ax12.add_artist(legend2)      
    cb = mpl.colorbar.ColorbarBase(ax12sub, cmap=cmx, norm=cNorm)
    cb.set_label(r"Min(D$_{\mathrm{massive}}$) [kpc]")
    fig12.tight_layout()
    fig12.show() 
    fig12.savefig(outfile_base + '_SMHMr_rClosestDist_Mpeak.png',dpi = dpi)

    plt.close('all')
    #fig12.clear()
    #fig2 = plt.figure(2, figsize=(plt_width*2,plt_width*aspect_ratio))
    #gs = gridspec.GridSpec(1,3,width_ratios=[15,15,1])
    #ax2a = fig12.add_subplot(gs[0])
    #ax2b = fig12.add_subplot(gs[1])
    #ax2sub = fig12.add_subplot(gs[2])
    cNorm  = colors.LogNorm(vmin=10, vmax = 2e3)
    fig12 = plt.figure(figsize=(plt_width,plt_width*aspect_ratio*2))
    gs = fig12.add_gridspec(2, hspace=0)
    axs12 = gs.subplots(sharex = True) #, constrained_layout=True)
    axs12 = axs12.flatten()
    ax12a = axs12[0]
    ax12b = axs12[1]
    ax12a.fill_between(read_abunmatch['M200'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_low'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_high'][read_abunmatch['M200']> 1e9],color = 'grey',linewidth=0,alpha = 0.6)
    ax12a.fill_between(upper_fid_Mhalo, np.interp( upper_fid_Mhalo, lower_fid_Mhalo, lower_fid_Mstar ), upper_fid_Mstar, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3)#Nadler outer
    ax12a.fill_between(upper_fid_inner_Mhalo, np.interp( upper_fid_inner_Mhalo, lower_fid_inner_Mhalo, lower_fid_inner_Mstar ),upper_fid_inner_Mstar, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3) #Nadler inner
    ax12a.fill_between(Jethwa_upper_outer_Mhalo, np.interp( Jethwa_upper_outer_Mhalo, Jethwa_lower_outer_Mhalo, Jethwa_lower_outer_Mstar ),Jethwa_upper_outer_Mstar, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3)#Jethwa outer
    ax12a.fill_between(Jethwa_upper_inner_Mhalo,np.interp( Jethwa_upper_inner_Mhalo, Jethwa_lower_inner_Mhalo, Jethwa_lower_inner_Mstar ), Jethwa_upper_inner_Mstar, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3) #Jethwa inner
    read, = ax12b.plot( [], [],  color="grey", linestyle="solid", linewidth=4) #, label="Read+ 2017" )
    nadler, = ax12b.plot( [], [],  color="#33669A", linestyle="solid", linewidth=8) #, label="Nadler+ 2020" )
    jethwa, = ax12b.plot( [], [],  color="#66CC66", linestyle="solid", linewidth=8) #, label="Jethwa+ 2018" )     
    cen_plt = ax12a.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Central'],c = (objs_pd_comb['min_dist'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax12a.scatter(objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb[mstar_key][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = (objs_pd_comb['min_dist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]*1.5, marker = "*", linewidths = edgewidth)    
    sat_plt =  ax12a.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Satellite'],c = (objs_pd_comb['min_dist'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k', marker = "*", s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    #ax12.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax12a.set_xscale('log')
    ax12a.set_yscale('log')
    ax12a.set_ylabel(mstar_axes_label)
    #ax12a.set_xlabel(r'M$_{\mathrm{vir, peak}}$/M$_\odot$')
    ax12a.axis([1e8, 2e11, 2e2, 5e9])
    #legend1 = ax12a.legend([cen_plt,sat_plt],['Central','Satellite'],scatterpoints = 1,facecolor = 'white',loc = 0,framealpha = 0,frameon = False)
    #ax12.add_artist(legend1)
    ax12a.label_outer()

    ax12b.fill_between(read_abunmatch['M200'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_low'][read_abunmatch['M200']> 1e9]/read_abunmatch['M200'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_high'][read_abunmatch['M200']> 1e9]/read_abunmatch['M200'][read_abunmatch['M200']> 1e9],color = 'grey',linewidth=0,alpha = 0.6)
    ax12b.fill_between(upper_fid_Mhalo, np.interp( upper_fid_Mhalo, lower_fid_Mhalo, lower_fid_Mstar )/upper_fid_Mhalo, upper_fid_Mstar/upper_fid_Mhalo, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3)#Nadler outer
    ax12b.fill_between(upper_fid_inner_Mhalo, np.interp( upper_fid_inner_Mhalo, lower_fid_inner_Mhalo, lower_fid_inner_Mstar)/upper_fid_inner_Mhalo,upper_fid_inner_Mstar/upper_fid_inner_Mhalo, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3) #Nadler inner
    ax12b.fill_between(Jethwa_upper_outer_Mhalo, np.interp( Jethwa_upper_outer_Mhalo, Jethwa_lower_outer_Mhalo, Jethwa_lower_outer_Mstar )/Jethwa_upper_outer_Mhalo,Jethwa_upper_outer_Mstar/Jethwa_upper_outer_Mhalo, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3)#Jethwa outer
    ax12b.fill_between(Jethwa_upper_inner_Mhalo,np.interp( Jethwa_upper_inner_Mhalo, Jethwa_lower_inner_Mhalo, Jethwa_lower_inner_Mstar )/Jethwa_upper_inner_Mhalo, Jethwa_upper_inner_Mstar/Jethwa_upper_inner_Mhalo, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3) #Jethwa inner
    read, = ax12b.plot( [], [],  color="grey", linestyle="solid", linewidth=4) #, label="Read+ 2017" )
    nadler, = ax12b.plot( [], [],  color="#33669A", linestyle="solid", linewidth=8) #, label="Nadler+ 2020" )
    jethwa, = ax12b.plot( [], [],  color="#66CC66", linestyle="solid", linewidth=8) #, label="Jethwa+ 2018" )        
    ax12b.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Central']/objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],c = (objs_pd_comb['min_dist'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax12b.scatter(objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb[mstar_key][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]/objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = (objs_pd_comb['min_dist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]*1.5, marker = "*", linewidths = edgewidth)    
    ax12b.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Satellite']/objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],c = (objs_pd_comb['min_dist'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k', marker = "*", s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    #ax12.scatter(readdata['M200']*1e10,readdata['M_*']*1e7/(readdata['M200']*1e10),marker = '+',c = 'k')
    ax12b.set_xscale('log')
    ax12b.set_yscale('log')
    ax12b.set_ylabel(mstar_ratio_axes_label)
    ax12b.set_xlabel(r'M$_{\mathrm{vir, peak}}$/M$_\odot$')
    ax12b.axis([1e8, 2e11, 1e-6, 0.05])
    legend1 = ax12a.legend([cen_plt,sat_plt],['Central','Satellite/Backsplash'],scatterpoints = 1,facecolor = 'white',loc = 0,framealpha = 0,frameon = False)
    legend2 = ax12b.legend([read,nadler,jethwa],['Read+ 2017','Nadler+ 2020','Jethwa+ 2018'],facecolor = 'white',loc = 4,framealpha = 0,frameon = False)
    #legend2 = ax12b.legend([read,nadler,jethwa],['Read+ 2017','Nadler+ 2020','Jethwa+ 2018'],facecolor = 'white',loc = 0,framealpha = 0,frameon = False)
    ax12a.add_artist(legend1)
    ax12b.add_artist(legend2)
    ax12b.label_outer()
    fig12.tight_layout()
    
    sm = plt.cm.ScalarMappable(cmap=cmx, norm=cNorm)
    cbar = plt.colorbar(sm, location="right",ax=[ax12a,ax12b],aspect = 40) #pad=0.04)
    #cbar = mpl.colorbar.ColorbarBase(ax12sub, cmap=cmx, norm=cNorm, location="right",ax=[ax12a,ax12b]) #pad=0.04)
    cbar.set_label(r'Min(D$_{\mathrm{massive}}$ [kpc])') #, rotation=0)
    #fig12.tight_layout()
    #sm = plt.cm.ScalarMappable(cmap=cmx, norm=cNorm)
    #sm.set_array([])
    #cbar = plt.colorbar(sm, location="top", ax=axs12.ravel().tolist(), fraction = 0.1) #pad=0.04)
    #cbar.set_label(r'Min(D$_{\mathrm{massive}}$) [kpc]') #, rotation=0)
    #cbar.ax.get_xaxis().labelpad = 15
    fig12.show()
    fig12.savefig(outfile_base + '_SMHM2_rClosestDist_Mpeak.png',dpi = dpi)

    # SMHM colored by maximum tidal index
    cmx_r = plt.get_cmap("viridis_r") 
    fig13 = plt.figure(13,figsize=(plt_width,plt_width*aspect_ratio))
    fig13.set_size_inches(plt_width,plt_width*aspect_ratio)
    fig13.clear()
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax13 = fig13.add_subplot(gs[0])
    ax13sub = fig13.add_subplot(gs[1])
    cNorm  = colors.LogNorm(vmin=1e4, vmax = 1e8)
    ax13.fill_between(read_abunmatch['M200'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_low'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_high'][read_abunmatch['M200']> 1e9],color = 'grey',linewidth=0,alpha=0.6)
    ax13.fill_between(upper_fid_Mhalo, np.interp( upper_fid_Mhalo, lower_fid_Mhalo, lower_fid_Mstar ), upper_fid_Mstar, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3)#Nadler outer
    ax13.fill_between(upper_fid_inner_Mhalo, np.interp( upper_fid_inner_Mhalo, lower_fid_inner_Mhalo, lower_fid_inner_Mstar ),upper_fid_inner_Mstar, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3) #Nadler inner
    ax13.fill_between(Jethwa_upper_outer_Mhalo, np.interp( Jethwa_upper_outer_Mhalo, Jethwa_lower_outer_Mhalo, Jethwa_lower_outer_Mstar ),Jethwa_upper_outer_Mstar, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3)#Jethwa outer
    ax13.fill_between(Jethwa_upper_inner_Mhalo,np.interp( Jethwa_upper_inner_Mhalo, Jethwa_lower_inner_Mhalo, Jethwa_lower_inner_Mstar ), Jethwa_upper_inner_Mstar, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3) #Jethwa inner
    read, = ax13.plot( [], [],  color="grey", linestyle="solid", linewidth=4) #, label="Read+ 2017" )
    nadler, = ax13.plot( [], [],  color="#33669A", linestyle="solid", linewidth=8) #, label="Nadler+ 2020" )
    jethwa, = ax13.plot( [], [],  color="#66CC66", linestyle="solid", linewidth=8) #, label="Jethwa+ 2018" )     
    cen_plt = ax13.scatter(objs_pd_comb['mass'][objs_pd_comb['type']=='Central'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Central'],c = objs_pd_comb['max_tindex'][objs_pd_comb['type']=='Central'], cmap = cmx_r, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax13.scatter(objs_pd_comb['mass'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb[mstar_key][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = objs_pd_comb['max_tindex'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], cmap = cmx_r, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]*1.5, marker = "*", linewidths = edgewidth)    
    sat_plt = ax13.scatter(objs_pd_comb['mass'][objs_pd_comb['type']=='Satellite'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Satellite'],c = objs_pd_comb['max_tindex'][objs_pd_comb['type']=='Satellite'], cmap = cmx_r, norm = cNorm,edgecolor = 'k',marker = '*', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    #ax13.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax13.set_xscale('log')
    ax13.set_yscale('log') 
    ax13.set_ylabel(r'M$_*$/M$_\odot$')
    ax13.set_xlabel(r'M$_{\mathrm{vir}}$/M$_\odot$')
    ax13.axis([2e6, 2e11, 2e2, 5e9])
    legend1 = ax13.legend([cen_plt,sat_plt],['Central','Satellite'],scatterpoints = 1,facecolor = 'white',loc = 4,framealpha = 0,frameon = False)
    legend2 = ax13.legend([read,nadler,jethwa],['Read+ 2017','Nadler+ 2020','Jethwa+ 2018'],facecolor = 'white',loc = 0,framealpha = 0,frameon = False)
    ax13.add_artist(legend1)
    ax13.add_artist(legend2)    
    cb = mpl.colorbar.ColorbarBase(ax13sub, cmap=cmx_r, norm=cNorm)
    cb.set_label(r"Tidal Index")
    fig13.tight_layout()
    fig13.show()    
    fig13.savefig(outfile_base + '_SMHM_rtindex.png',dpi = dpi)

    #SMHM colored by maximum tidal index
    #plt.clf()
    #fig13.clear()
    plt.close('all')
    fig13 = plt.figure(13,figsize=(plt_width,plt_width*aspect_ratio))
    fig13.set_size_inches(plt_width,plt_width*aspect_ratio)
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax13 = fig13.add_subplot(gs[0])
    ax13sub = fig13.add_subplot(gs[1])
    cNorm  = colors.LogNorm(vmin=1e4, vmax = 1e8)
    ax13.fill_between(read_abunmatch['M200'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_low'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_high'][read_abunmatch['M200']> 1e9],color = 'grey',linewidth=0,alpha=0.6)
    ax13.fill_between(upper_fid_Mhalo, np.interp( upper_fid_Mhalo, lower_fid_Mhalo, lower_fid_Mstar ), upper_fid_Mstar, facecolor='y',linewidth=0,zorder=0,alpha=0.3)#Nadler outer  facecolor="#33669A"
    ax13.fill_between(upper_fid_inner_Mhalo, np.interp( upper_fid_inner_Mhalo, lower_fid_inner_Mhalo, lower_fid_inner_Mstar ),upper_fid_inner_Mstar, facecolor='y',linewidth=0,zorder=0,alpha=0.3) #Nadler inner  facecolor="#33669A"
    ax13.fill_between(Jethwa_upper_outer_Mhalo, np.interp( Jethwa_upper_outer_Mhalo, Jethwa_lower_outer_Mhalo, Jethwa_lower_outer_Mstar ),Jethwa_upper_outer_Mstar, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3)#Jethwa outer
    ax13.fill_between(Jethwa_upper_inner_Mhalo,np.interp( Jethwa_upper_inner_Mhalo, Jethwa_lower_inner_Mhalo, Jethwa_lower_inner_Mstar ), Jethwa_upper_inner_Mstar, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3) #Jethwa inner
    read, = ax13.plot( [], [],  color="grey", linestyle="solid", linewidth=4) #, label="Read+ 2017" )
    nadler, = ax13.plot( [], [],  color="#33669A", linestyle="solid", linewidth=8) #, label="Nadler+ 2020" )
    jethwa, = ax13.plot( [], [],  color="#66CC66", linestyle="solid", linewidth=8) #, label="Jethwa+ 2018" )    
    #ax13.scatter(objs_rom['Mpeak'],objs_rom[mstar_key],c = objs_rom['massiveDist'],cmap = cmx_r, norm = cNorm,alpha = 0.4, s = markersize*0.5)
    cen_plt = ax13.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Central'],c = objs_pd_comb['max_tindex'][objs_pd_comb['type']=='Central'], cmap = cmx_r, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax13.scatter(objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb[mstar_key][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = objs_pd_comb['max_tindex'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], cmap = cmx_r, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]*1.5, marker = "*", linewidths = edgewidth)     
    sat_plt = ax13.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Satellite'],c = objs_pd_comb['max_tindex'][objs_pd_comb['type']=='Satellite'], cmap = cmx_r, norm = cNorm,edgecolor = 'k', marker = "*", s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    #ax13.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax13.set_xscale('log')
    ax13.set_yscale('log')
    ax13.set_ylabel(mstar_axes_label)
    ax13.set_xlabel(r'M$_{\mathrm{vir, peak}}$/M$_\odot$')
    ax13.axis([1e8, 2e11, 2e2, 5e9]) #ax13.axis([1e8, 2e13, 2e2, 5e9])
    legend1 = ax13.legend([cen_plt,sat_plt],['Central','Satellite'],scatterpoints = 1,facecolor = 'white',loc = 0,framealpha = 0,frameon = False)
    legend2 = ax13.legend([read,nadler,jethwa],['Read+ 2017','Nadler+ 2020','Jethwa+ 2018'],facecolor = 'white',loc = 4,framealpha = 0,frameon = False)
    ax13.add_artist(legend1)
    ax13.add_artist(legend2)      
    cb = mpl.colorbar.ColorbarBase(ax13sub, cmap=cmx_r, norm=cNorm)
    cb.set_label(r"Tidal Index")
    fig13.tight_layout()
    fig13.show() 
    fig13.savefig(outfile_base + '_SMHM_rtindex_Mpeak.png',dpi = dpi)
    
    plt.close('all')
    #fig13.clear()
    fig13 = plt.figure(13,figsize=(plt_width,plt_width*aspect_ratio))
    fig13.set_size_inches(plt_width,plt_width*aspect_ratio)
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax13 = fig13.add_subplot(gs[0])
    ax13sub = fig13.add_subplot(gs[1])
    cNorm  = colors.LogNorm(vmin=1e4, vmax = 1e8)
    ax13.fill_between(read_abunmatch['M200'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_low'][read_abunmatch['M200']> 1e9]/read_abunmatch['M200'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_high'][read_abunmatch['M200']> 1e9]/read_abunmatch['M200'][read_abunmatch['M200']> 1e9],color = 'grey',linewidth=0,alpha = 0.6)
    ax13.fill_between(upper_fid_Mhalo, np.interp( upper_fid_Mhalo, lower_fid_Mhalo, lower_fid_Mstar )/upper_fid_Mhalo, upper_fid_Mstar/upper_fid_Mhalo, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3)#Nadler outer
    ax13.fill_between(upper_fid_inner_Mhalo, np.interp( upper_fid_inner_Mhalo, lower_fid_inner_Mhalo, lower_fid_inner_Mstar)/upper_fid_inner_Mhalo,upper_fid_inner_Mstar/upper_fid_inner_Mhalo, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3) #Nadler inner
    ax13.fill_between(Jethwa_upper_outer_Mhalo, np.interp( Jethwa_upper_outer_Mhalo, Jethwa_lower_outer_Mhalo, Jethwa_lower_outer_Mstar )/Jethwa_upper_outer_Mhalo,Jethwa_upper_outer_Mstar/Jethwa_upper_outer_Mhalo, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3)#Jethwa outer
    ax13.fill_between(Jethwa_upper_inner_Mhalo,np.interp( Jethwa_upper_inner_Mhalo, Jethwa_lower_inner_Mhalo, Jethwa_lower_inner_Mstar )/Jethwa_upper_inner_Mhalo, Jethwa_upper_inner_Mstar/Jethwa_upper_inner_Mhalo, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3) #Jethwa inner
    read, = ax13.plot( [], [],  color="grey", linestyle="solid", linewidth=4) #, label="Read+ 2017" )
    nadler, = ax13.plot( [], [],  color="#33669A", linestyle="solid", linewidth=8) #, label="Nadler+ 2020" )
    jethwa, = ax13.plot( [], [],  color="#66CC66", linestyle="solid", linewidth=8) #, label="Jethwa+ 2018" )        
    ax13.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Central']/objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],c = (objs_pd_comb['max_tindex'][objs_pd_comb['type']=='Central']), cmap = cmx_r, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax13.scatter(objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb[mstar_key][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]/objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = (objs_pd_comb['max_tindex'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]), cmap = cmx_r, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]*1.5, marker = "*", linewidths = edgewidth)
    ax13.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Satellite']/objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],c = (objs_pd_comb['max_tindex'][objs_pd_comb['type']=='Satellite']), cmap = cmx_r, norm = cNorm,edgecolor = 'k', marker = "*", s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    #ax13.scatter(readdata['M200']*1e10,readdata['M_*']*1e7/(readdata['M200']*1e10),marker = '+',c = 'k')
    ax13.set_xscale('log')
    ax13.set_yscale('log')
    ax13.set_ylabel(mstar_ratio_axes_label)
    ax13.set_xlabel(r'M$_{\mathrm{vir, peak}}$/M$_\odot$')
    ax13.axis([1e8, 2e11, 1e-6, 0.05])
    legend1 = ax13.legend([cen_plt,sat_plt],['Central','Satellite'],scatterpoints = 1,facecolor = 'white',loc = 4,framealpha = 0,frameon = False)
    legend2 = ax13.legend([read,nadler,jethwa],['Read+ 2017','Nadler+ 2020','Jethwa+ 2018'],facecolor = 'white',loc = 0,framealpha = 0,frameon = False)
    ax13.add_artist(legend1)
    ax13.add_artist(legend2)      
    cb = mpl.colorbar.ColorbarBase(ax13sub, cmap=cmx_r, norm=cNorm)
    cb.set_label(r"Tidal Index")
    fig13.tight_layout()
    fig13.show() 
    fig13.savefig(outfile_base + '_SMHMr_rtindex_Mpeak.png',dpi = dpi)

    plt.close('all')
    #fig13.clear()
    #fig2 = plt.figure(2, figsize=(plt_width*2,plt_width*aspect_ratio))
    #gs = gridspec.GridSpec(1,3,width_ratios=[15,15,1])
    #ax2a = fig13.add_subplot(gs[0])
    #ax2b = fig13.add_subplot(gs[1])
    #ax2sub = fig13.add_subplot(gs[2])
    fig13 = plt.figure(figsize=(plt_width,plt_width*aspect_ratio*2))
    gs = fig13.add_gridspec(2, hspace=0)
    axs13 = gs.subplots(sharex = True) #, constrained_layout=True)
    axs13 = axs13.flatten()
    ax13a = axs13[0]
    ax13b = axs13[1]
    cNorm  = colors.LogNorm(vmin=1e2, vmax = 1e8)
    ax13a.fill_between(read_abunmatch['M200'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_low'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_high'][read_abunmatch['M200']> 1e9],color = 'grey',linewidth=0,alpha = 0.6)
    ax13a.fill_between(upper_fid_Mhalo, np.interp( upper_fid_Mhalo, lower_fid_Mhalo, lower_fid_Mstar ), upper_fid_Mstar, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3)#Nadler outer
    ax13a.fill_between(upper_fid_inner_Mhalo, np.interp( upper_fid_inner_Mhalo, lower_fid_inner_Mhalo, lower_fid_inner_Mstar ),upper_fid_inner_Mstar, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3) #Nadler inner
    ax13a.fill_between(Jethwa_upper_outer_Mhalo, np.interp( Jethwa_upper_outer_Mhalo, Jethwa_lower_outer_Mhalo, Jethwa_lower_outer_Mstar ),Jethwa_upper_outer_Mstar, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3)#Jethwa outer
    ax13a.fill_between(Jethwa_upper_inner_Mhalo,np.interp( Jethwa_upper_inner_Mhalo, Jethwa_lower_inner_Mhalo, Jethwa_lower_inner_Mstar ), Jethwa_upper_inner_Mstar, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3) #Jethwa inner
    read, = ax13b.plot( [], [],  color="grey", linestyle="solid", linewidth=4) #, label="Read+ 2017" )
    nadler, = ax13b.plot( [], [],  color="#33669A", linestyle="solid", linewidth=8) #, label="Nadler+ 2020" )
    jethwa, = ax13b.plot( [], [],  color="#66CC66", linestyle="solid", linewidth=8) #, label="Jethwa+ 2018" )     
    cen_plt = ax13a.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Central'],c = (objs_pd_comb['max_tindex'][objs_pd_comb['type']=='Central']), cmap = cmx_r, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax13a.scatter(objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb[mstar_key][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = (objs_pd_comb['max_tindex'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]), cmap = cmx_r, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]*1.5, marker = "*", linewidths = edgewidth)    
    sat_plt =  ax13a.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Satellite'],c = (objs_pd_comb['max_tindex'][objs_pd_comb['type']=='Satellite']), cmap = cmx_r, norm = cNorm,edgecolor = 'k', marker = "*", s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    #ax13.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax13a.set_xscale('log')
    ax13a.set_yscale('log')
    ax13a.set_ylabel(mstar_axes_label)
    #ax13a.set_xlabel(r'M$_{\mathrm{vir, peak}}$/M$_\odot$')
    ax13a.axis([1e8, 2e11, 2e2, 5e9])
    #legend1 = ax13a.legend([cen_plt,sat_plt],['Central','Satellite'],scatterpoints = 1,facecolor = 'white',loc = 0,framealpha = 0,frameon = False)
    legend2 = ax13a.legend([read,nadler,jethwa],['Read+ 2017','Nadler+ 2020','Jethwa+ 2018'],facecolor = 'white',loc = 4,framealpha = 0,frameon = False)
    #ax13.add_artist(legend1)
    ax13a.add_artist(legend2)
    ax13a.label_outer()

    ax13b.fill_between(read_abunmatch['M200'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_low'][read_abunmatch['M200']> 1e9]/read_abunmatch['M200'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_high'][read_abunmatch['M200']> 1e9]/read_abunmatch['M200'][read_abunmatch['M200']> 1e9],color = 'grey',linewidth=0,alpha = 0.6)
    ax13b.fill_between(upper_fid_Mhalo, np.interp( upper_fid_Mhalo, lower_fid_Mhalo, lower_fid_Mstar )/upper_fid_Mhalo, upper_fid_Mstar/upper_fid_Mhalo, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3)#Nadler outer
    ax13b.fill_between(upper_fid_inner_Mhalo, np.interp( upper_fid_inner_Mhalo, lower_fid_inner_Mhalo, lower_fid_inner_Mstar)/upper_fid_inner_Mhalo,upper_fid_inner_Mstar/upper_fid_inner_Mhalo, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3) #Nadler inner
    ax13b.fill_between(Jethwa_upper_outer_Mhalo, np.interp( Jethwa_upper_outer_Mhalo, Jethwa_lower_outer_Mhalo, Jethwa_lower_outer_Mstar )/Jethwa_upper_outer_Mhalo,Jethwa_upper_outer_Mstar/Jethwa_upper_outer_Mhalo, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3)#Jethwa outer
    ax13b.fill_between(Jethwa_upper_inner_Mhalo,np.interp( Jethwa_upper_inner_Mhalo, Jethwa_lower_inner_Mhalo, Jethwa_lower_inner_Mstar )/Jethwa_upper_inner_Mhalo, Jethwa_upper_inner_Mstar/Jethwa_upper_inner_Mhalo, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3) #Jethwa inner
    read, = ax13b.plot( [], [],  color="grey", linestyle="solid", linewidth=4) #, label="Read+ 2017" )
    nadler, = ax13b.plot( [], [],  color="#33669A", linestyle="solid", linewidth=8) #, label="Nadler+ 2020" )
    jethwa, = ax13b.plot( [], [],  color="#66CC66", linestyle="solid", linewidth=8) #, label="Jethwa+ 2018" )        
    ax13b.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Central']/objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],c = (objs_pd_comb['max_tindex'][objs_pd_comb['type']=='Central']), cmap = cmx_r, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax13b.scatter(objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb[mstar_key][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]/objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = (objs_pd_comb['max_tindex'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]), cmap = cmx_r, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]*1.5, marker = "*", linewidths = edgewidth)    
    ax13b.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Satellite']/objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],c = (objs_pd_comb['max_tindex'][objs_pd_comb['type']=='Satellite']), cmap = cmx_r, norm = cNorm,edgecolor = 'k', marker = "*", s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    #ax13.scatter(readdata['M200']*1e10,readdata['M_*']*1e7/(readdata['M200']*1e10),marker = '+',c = 'k')
    ax13b.set_xscale('log')
    ax13b.set_yscale('log')
    ax13b.set_ylabel(mstar_ratio_axes_label)
    ax13b.set_xlabel(r'M$_{\mathrm{vir, peak}}$/M$_\odot$')
    ax13b.axis([1e8, 2e11, 1e-6, 0.05])
    legend1 = ax13b.legend([cen_plt, sat_plt],['Central','Satellite/Backsplash'],scatterpoints = 1,facecolor = 'white',loc = 0,framealpha = 0,frameon = False)
    #legend2 = ax13b.legend([read,nadler,jethwa],['Read+ 2017','Nadler+ 2020','Jethwa+ 2018'],facecolor = 'white',loc = 0,framealpha = 0,frameon = False)
    ax13b.add_artist(legend1)
    #ax13b.add_artist(legend2)
    ax13b.label_outer()
    fig13.tight_layout()
    
    sm = plt.cm.ScalarMappable(cmap=cmx_r, norm=cNorm)
    cbar = plt.colorbar(sm, location="right",ax=[ax13a,ax13b],aspect = 40) #pad=0.04)
    #cbar = mpl.colorbar.ColorbarBase(ax13sub, cmap=cmx_r, norm=cNorm, location="right",ax=[ax13a,ax13b]) #pad=0.04)
    cbar.set_label(r'Tidal Index') #, rotation=0)
    #fig13.tight_layout()
    #sm = plt.cm.ScalarMappable(cmap=cmx_r, norm=cNorm)
    #sm.set_array([])
    #cbar = plt.colorbar(sm, location="top", ax=axs13.ravel().tolist(), fraction = 0.1) #pad=0.04)
    #cbar.set_label(r'D$_{\mathrm{massive}}$ [kpc]') #, rotation=0)
    #cbar.ax.get_xaxis().labelpad = 15
    fig13.show()
    fig13.savefig(outfile_base + '_SMHM2_rtindex_Mpeak.png',dpi = dpi)

    ###################### Evolution with Time #################
    plt.close('all')
    #fig14.clear()
    #fig2 = plt.figure(2, figsize=(plt_width*2,plt_width*aspect_ratio))
    #gs = gridspec.GridSpec(1,3,width_ratios=[15,15,1])
    #ax2a = fig14.add_subplot(gs[0])
    #ax2b = fig14.add_subplot(gs[1])
    #ax2sub = fig14.add_subplot(gs[2])
    fig14, ax14 = plt.subplots(figsize=(plt_width,plt_width*aspect_ratio))
    #axs14 = axs14.flatten()
    #ax14 = axs14[0]
    cmx = plt.get_cmap("viridis")
    cNorm  = colors.Normalize(vmin=0, vmax = 13.8)
    sm = plt.cm.ScalarMappable(cmap=cmx, norm=cNorm)
    dmax = 11; dmin = 8; nbins = 8
    xaxis, data = bin_plt(np.log10(np.array(objs_pd_comb.loc[objs_pd_comb['vmass_z4']> 0, 'vmass_z4'].tolist())), np.array(objs_pd_comb.loc[objs_pd_comb['vmass_z4']> 0, 'smass_z4'].tolist()), xmax = dmax, xmin = dmin, nbins = nbins)
    ax14.fill_between(10**xaxis[data[:,4] >= 0], data[data[:,4] >= 0,0], data[data[:,4] >= 0,4], color=sm.to_rgba(time_z4), alpha=0.1, linewidth=0)
    ax14.fill_between(10**xaxis[data[:,3] >= 0], data[data[:,3] >= 0,1], data[data[:,3] >= 0,3], color=sm.to_rgba(time_z4), alpha=0.1, linewidth=0)
    pltz4, = ax14.plot(10**xaxis[data[:,2] >= 0], data[data[:,2] >= 0,2], color=sm.to_rgba(time_z4), label = "z = 4")
    xaxis, data = bin_plt(np.log10(np.array(objs_pd_comb.loc[objs_pd_comb['vmass_z2']> 0, 'vmass_z2'].tolist())), np.array(objs_pd_comb.loc[objs_pd_comb['vmass_z2']> 0, 'smass_z2'].tolist()), xmax = dmax, xmin = dmin, nbins = nbins)
    ax14.fill_between(10**xaxis[data[:,4] >= 0], data[data[:,4] >= 0,0], data[data[:,4] >= 0,4], color=sm.to_rgba(time_z2), alpha=0.1, linewidth=0)
    ax14.fill_between(10**xaxis[data[:,3] >= 0], data[data[:,3] >= 0,1], data[data[:,3] >= 0,3], color=sm.to_rgba(time_z2), alpha=0.1, linewidth=0)
    pltz2, = ax14.plot(10**xaxis[data[:,2] >= 0], data[data[:,2] >= 0,2], color=sm.to_rgba(time_z2), label = "z = 2")
    xaxis, data = bin_plt(np.log10(np.array(objs_pd_comb.loc[objs_pd_comb['vmass_z1']> 0, 'vmass_z1'].tolist())), np.array(objs_pd_comb.loc[objs_pd_comb['vmass_z1']> 0, 'smass_z1'].tolist()), xmax = dmax, xmin = dmin, nbins = nbins)
    ax14.fill_between(10**xaxis[data[:,4] >= 0], data[data[:,4] >= 0,0], data[data[:,4] >= 0,4], color=sm.to_rgba(time_z1), alpha=0.1, linewidth=0)
    ax14.fill_between(10**xaxis[data[:,3] >= 0], data[data[:,3] >= 0,1], data[data[:,3] >= 0,3], color=sm.to_rgba(time_z1), alpha=0.1, linewidth=0)
    pltz1, = ax14.plot(10**xaxis[data[:,2] >= 0], data[data[:,2] >= 0,2], color=sm.to_rgba(time_z1), label = "z = 1")
    xaxis, data = bin_plt(np.log10(np.array(objs_pd_comb.loc[objs_pd_comb['vmass_z0_5']> 0, 'vmass_z0_5'].tolist())), np.array(objs_pd_comb.loc[objs_pd_comb['vmass_z0_5']> 0, 'smass_z0_5'].tolist()), xmax = dmax, xmin = dmin, nbins = nbins)
    ax14.fill_between(10**xaxis[data[:,4] >= 0], data[data[:,4] >= 0,0], data[data[:,4] >= 0,4], color=sm.to_rgba(time_z0_5), alpha=0.1, linewidth=0)
    ax14.fill_between(10**xaxis[data[:,3] >= 0], data[data[:,3] >= 0,1], data[data[:,3] >= 0,3], color=sm.to_rgba(time_z0_5), alpha=0.1, linewidth=0)
    pltz05, = ax14.plot(10**xaxis[data[:,2] >= 0], data[data[:,2] >= 0,2], color=sm.to_rgba(time_z0_5), label = "z = 0.5")
    xaxis, data = bin_plt(np.log10(np.array(objs_pd_comb.loc[objs_pd_comb['Mpeak_snap'] == 4096, 'mvir'].tolist())), np.array(objs_pd_comb.loc[objs_pd_comb['Mpeak_snap']==4096, 'Mstar_z0'].tolist()), xmax = dmax, xmin = dmin, nbins = nbins)
    ax14.fill_between(10**xaxis[data[:,4] >= 0], data[data[:,4] >= 0,0], data[data[:,4] >= 0,4], color=sm.to_rgba(13.8), alpha=0.1, linewidth=0)
    ax14.fill_between(10**xaxis[data[:,3] >= 0], data[data[:,3] >= 0,1], data[data[:,3] >= 0,3], color=sm.to_rgba(13.8), alpha=0.1, linewidth=0)
    pltz0, = ax14.plot(10**xaxis[data[:,2] >= 0], data[data[:,2] >= 0,2], color=sm.to_rgba(13.8), label = "z = 0")
    
    cen_plt = ax14.scatter( objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Central'],c = (objs_pd_comb['Mpeak_snap'][objs_pd_comb['type']=='Central'])*13.8/4096, cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax14.scatter(           objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb[mstar_key][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = (objs_pd_comb['Mpeak_snap'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')])*13.8/4096, cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]*1.5, marker = "*", linewidths = edgewidth)
    sat_plt =  ax14.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Satellite'],c = (objs_pd_comb['Mpeak_snap'][objs_pd_comb['type']=='Satellite'])*13.8/4096, cmap = cmx, norm = cNorm,edgecolor = 'k', marker = "*", s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    ax14.set_xscale('log')
    ax14.set_yscale('log')
    ax14.set_ylabel(mstar_axes_label)
    ax14.set_xlabel(r'M$_{\mathrm{vir, peak}}$/M$_\odot$')
    ax14.axis([1e8, 2e11, 2e2, 5e9])
    legend1 = ax14.legend([cen_plt, sat_plt],['Central','Satellite/Backsplash'],scatterpoints = 1,facecolor = 'white',loc = 'lower right',framealpha = 0,frameon = False)
    ax14.add_artist(legend1)
    legend2 = ax14.legend(handles = [pltz4, pltz2, pltz1, pltz05, pltz0], loc = 0, frameon = False)
    ax14.add_artist(legend2)
    #ax14.label_outer()
    fig14.tight_layout()
    
    sm = plt.cm.ScalarMappable(cmap=cmx, norm=cNorm)
    cbar = plt.colorbar(sm, location="right",ax=[ax14],aspect = 40) #pad=0.04)
    #cbar = mpl.colorbar.ColorbarBase(ax14sub, cmap=cmx, norm=cNorm, location="right",ax=[ax14a,ax14b]) #pad=0.04)
    cbar.set_label(r'Time of Peak M$_{\mathrm{vir}}$ [Gyr]') #, rotation=0)
    #fig14.tight_layout()
    #sm = plt.cm.ScalarMappable(cmap=cmx, norm=cNorm)
    #sm.set_array([])
    #cbar = plt.colorbar(sm, location="top", ax=axs14.ravel().tolist(), fraction = 0.1) #pad=0.04)
    #cbar.set_label(r'D$_{\mathrm{massive}}$ [kpc]') #, rotation=0)
    #cbar.ax.get_xaxis().labelpad = 15
    fig14.show()
    fig14.savefig(outfile_base + '_SMHM_tmax_Mpeak.png',dpi = dpi)

    plt.close('all')
    #fig14.clear()
    #fig2 = plt.figure(2, figsize=(plt_width*2,plt_width*aspect_ratio))
    #gs = gridspec.GridSpec(1,3,width_ratios=[15,15,1])
    #ax2a = fig14.add_subplot(gs[0])
    #ax2b = fig14.add_subplot(gs[1])
    #ax2sub = fig14.add_subplot(gs[2])
    fig14 = plt.figure(figsize=(plt_width,plt_width*aspect_ratio*2))
    gs = fig14.add_gridspec(2, hspace=0)
    axs14 = gs.subplots(sharex = True) #, constrained_layout=True)
    axs14 = axs14.flatten()
    ax14a = axs14[0]
    ax14b = axs14[1]
    cmx = plt.get_cmap("viridis")
    cNorm  = colors.Normalize(vmin=0, vmax = 13.8)
    sm = plt.cm.ScalarMappable(cmap=cmx, norm=cNorm)
    """
    ax14a.fill_between(read_abunmatch['M200'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_low'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_high'][read_abunmatch['M200']> 1e9],color = 'grey',linewidth=0,alpha = 0.6)
    ax14a.fill_between(upper_fid_Mhalo, np.interp( upper_fid_Mhalo, lower_fid_Mhalo, lower_fid_Mstar ), upper_fid_Mstar, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3)#Nadler outer
    ax14a.fill_between(upper_fid_inner_Mhalo, np.interp( upper_fid_inner_Mhalo, lower_fid_inner_Mhalo, lower_fid_inner_Mstar ),upper_fid_inner_Mstar, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3) #Nadler inner
    ax14a.fill_between(Jethwa_upper_outer_Mhalo, np.interp( Jethwa_upper_outer_Mhalo, Jethwa_lower_outer_Mhalo, Jethwa_lower_outer_Mstar ),Jethwa_upper_outer_Mstar, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3)#Jethwa outer
    ax14a.fill_between(Jethwa_upper_inner_Mhalo,np.interp( Jethwa_upper_inner_Mhalo, Jethwa_lower_inner_Mhalo, Jethwa_lower_inner_Mstar ), Jethwa_upper_inner_Mstar, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3) #Jethwa inner
    read, = ax14b.plot( [], [],  color="grey", linestyle="solid", linewidth=4) #, label="Read+ 2017" )
    nadler, = ax14b.plot( [], [],  color="#33669A", linestyle="solid", linewidth=8) #, label="Nadler+ 2020" )
    jethwa, = ax14b.plot( [], [],  color="#66CC66", linestyle="solid", linewidth=8) #, label="Jethwa+ 2018" )   
    """
    """
    ax14a.scatter(           objs_pd_comb['vmass_z2'][objs_pd_comb['type']=='Central'],objs_pd_comb['smass_z2'][objs_pd_comb['type']=='Central'],c = time_z2+np.zeros(np.sum(objs_pd_comb['type']=='Central')), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth, alpha = 0.5)
    ax14a.scatter(           objs_pd_comb['vmass_z1'][objs_pd_comb['type']=='Central'],objs_pd_comb['smass_z1'][objs_pd_comb['type']=='Central'],c = time_z1+np.zeros(np.sum(objs_pd_comb['type']=='Central')), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth, alpha = 0.5)
    ax14a.scatter(           objs_pd_comb['vmass_z0_5'][objs_pd_comb['type']=='Central'],objs_pd_comb['smass_z0_5'][objs_pd_comb['type']=='Central'],c = time_z0_5+np.zeros(np.sum(objs_pd_comb['type']=='Central')), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth, alpha = 0.5)
    ax14a.scatter(           objs_pd_comb['vmass_z2'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb['smass_z2'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = time_z2+np.zeros(np.sum((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]*1.5, marker = "*", linewidths = edgewidth, alpha = 0.5)
    ax14a.scatter(           objs_pd_comb['vmass_z1'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb['smass_z1'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = time_z1+np.zeros(np.sum((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]*1.5, marker = "*", linewidths = edgewidth, alpha = 0.5)
    ax14a.scatter(           objs_pd_comb['vmass_z0_5'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb['smass_z0_5'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = time_z0_5+np.zeros(np.sum((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]*1.5, marker = "*", linewidths = edgewidth, alpha = 0.5)    
    ax14a.scatter(           objs_pd_comb['vmass_z2'][objs_pd_comb['type']=='Satellite'],objs_pd_comb['smass_z2'][objs_pd_comb['type']=='Satellite'],c = time_z2+np.zeros(np.sum(objs_pd_comb['type']=='Satellite')), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite'], linewidths = edgewidth, alpha = 0.5)
    ax14a.scatter(           objs_pd_comb['vmass_z1'][objs_pd_comb['type']=='Satellite'],objs_pd_comb['smass_z1'][objs_pd_comb['type']=='Satellite'],c = time_z1+np.zeros(np.sum(objs_pd_comb['type']=='Satellite')), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite'], linewidths = edgewidth, alpha = 0.5)
    ax14a.scatter(           objs_pd_comb['vmass_z0_5'][objs_pd_comb['type']=='Satellite'],objs_pd_comb['smass_z0_5'][objs_pd_comb['type']=='Satellite'],c = time_z0_5+np.zeros(np.sum(objs_pd_comb['type']=='Satellite')), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite'], linewidths = edgewidth, alpha = 0.5)
    """
    dmax = 11; dmin = 8
    xaxis, data = bin_plt(np.log10(np.array(objs_pd_comb.loc[objs_pd_comb['vmass_z4']> 0, 'vmass_z4'].tolist())), np.array(objs_pd_comb.loc[objs_pd_comb['vmass_z4']> 0, 'smass_z4'].tolist()), xmax = dmax, xmin = dmin, nbins = nbins)
    ax14a.fill_between(10**xaxis[data[:,4] >= 0], data[data[:,4] >= 0,0], data[data[:,4] >= 0,4], color=sm.to_rgba(time_z4), alpha=0.1, linewidth=0)
    ax14a.fill_between(10**xaxis[data[:,3] >= 0], data[data[:,3] >= 0,1], data[data[:,3] >= 0,3], color=sm.to_rgba(time_z4), alpha=0.1, linewidth=0)
    pltz4, = ax14a.plot(10**xaxis[data[:,2] >= 0], data[data[:,2] >= 0,2], color=sm.to_rgba(time_z4), label = "z = 4")
    xaxis, data = bin_plt(np.log10(np.array(objs_pd_comb.loc[objs_pd_comb['vmass_z2']> 0, 'vmass_z2'].tolist())), np.array(objs_pd_comb.loc[objs_pd_comb['vmass_z2']> 0, 'smass_z2'].tolist()), xmax = dmax, xmin = dmin, nbins = nbins)
    ax14a.fill_between(10**xaxis[data[:,4] >= 0], data[data[:,4] >= 0,0], data[data[:,4] >= 0,4], color=sm.to_rgba(time_z2), alpha=0.1, linewidth=0)
    ax14a.fill_between(10**xaxis[data[:,3] >= 0], data[data[:,3] >= 0,1], data[data[:,3] >= 0,3], color=sm.to_rgba(time_z2), alpha=0.1, linewidth=0)
    pltz2, = ax14a.plot(10**xaxis[data[:,2] >= 0], data[data[:,2] >= 0,2], color=sm.to_rgba(time_z2), label = "z = 2")
    xaxis, data = bin_plt(np.log10(np.array(objs_pd_comb.loc[objs_pd_comb['vmass_z1']> 0, 'vmass_z1'].tolist())), np.array(objs_pd_comb.loc[objs_pd_comb['vmass_z1']> 0, 'smass_z1'].tolist()), xmax = dmax, xmin = dmin, nbins = nbins)
    ax14a.fill_between(10**xaxis[data[:,4] >= 0], data[data[:,4] >= 0,0], data[data[:,4] >= 0,4], color=sm.to_rgba(time_z1), alpha=0.1, linewidth=0)
    ax14a.fill_between(10**xaxis[data[:,3] >= 0], data[data[:,3] >= 0,1], data[data[:,3] >= 0,3], color=sm.to_rgba(time_z1), alpha=0.1, linewidth=0)
    pltz1, ax14a.plot(10**xaxis[data[:,2] >= 0], data[data[:,2] >= 0,2], color=sm.to_rgba(time_z1), label = "z = 1")
    xaxis, data = bin_plt(np.log10(np.array(objs_pd_comb.loc[objs_pd_comb['vmass_z0_5']> 0, 'vmass_z0_5'].tolist())), np.array(objs_pd_comb.loc[objs_pd_comb['vmass_z0_5']> 0, 'smass_z0_5'].tolist()), xmax = dmax, xmin = dmin, nbins = nbins)
    ax14a.fill_between(10**xaxis[data[:,4] >= 0], data[data[:,4] >= 0,0], data[data[:,4] >= 0,4], color=sm.to_rgba(time_z0_5), alpha=0.1, linewidth=0)
    ax14a.fill_between(10**xaxis[data[:,3] >= 0], data[data[:,3] >= 0,1], data[data[:,3] >= 0,3], color=sm.to_rgba(time_z0_5), alpha=0.1, linewidth=0)
    pltz05, = ax14a.plot(10**xaxis[data[:,2] >= 0], data[data[:,2] >= 0,2], color=sm.to_rgba(time_z0_5), label = "z = 0.5")
    xaxis, data = bin_plt(np.log10(np.array(objs_pd_comb.loc[objs_pd_comb['Mpeak_snap'] == 4096, 'mvir'].tolist())), np.array(objs_pd_comb.loc[objs_pd_comb['Mpeak_snap']==4096, 'Mstar_z0'].tolist()), xmax = dmax, xmin = dmin, nbins = nbins)
    ax14a.fill_between(10**xaxis[data[:,4] >= 0], data[data[:,4] >= 0,0], data[data[:,4] >= 0,4], color=sm.to_rgba(13.8), alpha=0.1, linewidth=0)
    ax14a.fill_between(10**xaxis[data[:,3] >= 0], data[data[:,3] >= 0,1], data[data[:,3] >= 0,3], color=sm.to_rgba(13.8), alpha=0.1, linewidth=0)
    pltz0, = ax14a.plot(10**xaxis[data[:,2] >= 0], data[data[:,2] >= 0,2], color=sm.to_rgba(13.8), label = "z = 0")
    cen_plt = ax14a.scatter( objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Central'],c = (objs_pd_comb['Mpeak_snap'][objs_pd_comb['type']=='Central'])*13.8/4096, cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax14a.scatter(           objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb[mstar_key][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = (objs_pd_comb['Mpeak_snap'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')])*13.8/4096, cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]*1.5, marker = "*", linewidths = edgewidth)
    sat_plt =  ax14a.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Satellite'],c = (objs_pd_comb['Mpeak_snap'][objs_pd_comb['type']=='Satellite'])*13.8/4096, cmap = cmx, norm = cNorm,edgecolor = 'k', marker = "*", s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    #ax14.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax14a.set_xscale('log')
    ax14a.set_yscale('log')
    ax14a.set_ylabel(mstar_axes_label)
    #ax14a.set_xlabel(r'M$_{\mathrm{vir, peak}}$/M$_\odot$')
    ax14a.axis([1e8, 2e11, 2e2, 5e9])
#    legend2 = ax14a.legend(handles = [pltz4, pltz2, pltz1, pltz05, pltz0], loc = 0, frameon = False)
    legend2 = ax14a.legend(handles = [pltz4, pltz2, pltz1, pltz05, pltz0], loc = 0, frameon = False)
    ax14a.add_artist(legend2)

    """
    ax14b.fill_between(read_abunmatch['M200'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_low'][read_abunmatch['M200']> 1e9]/read_abunmatch['M200'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_high'][read_abunmatch['M200']> 1e9]/read_abunmatch['M200'][read_abunmatch['M200']> 1e9],color = 'grey',linewidth=0,alpha = 0.6)
    ax14b.fill_between(upper_fid_Mhalo, np.interp( upper_fid_Mhalo, lower_fid_Mhalo, lower_fid_Mstar )/upper_fid_Mhalo, upper_fid_Mstar/upper_fid_Mhalo, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3)#Nadler outer
    ax14b.fill_between(upper_fid_inner_Mhalo, np.interp( upper_fid_inner_Mhalo, lower_fid_inner_Mhalo, lower_fid_inner_Mstar)/upper_fid_inner_Mhalo,upper_fid_inner_Mstar/upper_fid_inner_Mhalo, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3) #Nadler inner
    ax14b.fill_between(Jethwa_upper_outer_Mhalo, np.interp( Jethwa_upper_outer_Mhalo, Jethwa_lower_outer_Mhalo, Jethwa_lower_outer_Mstar )/Jethwa_upper_outer_Mhalo,Jethwa_upper_outer_Mstar/Jethwa_upper_outer_Mhalo, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3)#Jethwa outer
    ax14b.fill_between(Jethwa_upper_inner_Mhalo,np.interp( Jethwa_upper_inner_Mhalo, Jethwa_lower_inner_Mhalo, Jethwa_lower_inner_Mstar )/Jethwa_upper_inner_Mhalo, Jethwa_upper_inner_Mstar/Jethwa_upper_inner_Mhalo, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3) #Jethwa inner
    read, = ax14b.plot( [], [],  color="grey", linestyle="solid", linewidth=4) #, label="Read+ 2017" )
    nadler, = ax14b.plot( [], [],  color="#33669A", linestyle="solid", linewidth=8) #, label="Nadler+ 2020" )
    jethwa, = ax14b.plot( [], [],  color="#66CC66", linestyle="solid", linewidth=8) #, label="Jethwa+ 2018" )  
    """
    """
    ax14b.scatter(           objs_pd_comb['vmass_z2'][objs_pd_comb['type']=='Central'],objs_pd_comb['smass_z2'][objs_pd_comb['type']=='Central']/objs_pd_comb['vmass_z2'][objs_pd_comb['type']=='Central'],c = time_z2+np.zeros(np.sum(objs_pd_comb['type']=='Central')), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth, alpha = 0.5)
    ax14b.scatter(           objs_pd_comb['vmass_z1'][objs_pd_comb['type']=='Central'],objs_pd_comb['smass_z1'][objs_pd_comb['type']=='Central']/objs_pd_comb['vmass_z1'][objs_pd_comb['type']=='Central'],c = time_z1+np.zeros(np.sum(objs_pd_comb['type']=='Central')), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth, alpha = 0.5)
    ax14b.scatter(           objs_pd_comb['vmass_z0_5'][objs_pd_comb['type']=='Central'],objs_pd_comb['smass_z0_5'][objs_pd_comb['type']=='Central']/objs_pd_comb['vmass_z0_5'][objs_pd_comb['type']=='Central'],c = time_z0_5+np.zeros(np.sum(objs_pd_comb['type']=='Central')), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth, alpha = 0.5)
    ax14b.scatter(           objs_pd_comb['vmass_z2'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb['smass_z2'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]/objs_pd_comb['vmass_z2'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = time_z2+np.zeros(np.sum((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]*1.5, marker = "*", linewidths = edgewidth, alpha = 0.5)
    ax14b.scatter(           objs_pd_comb['vmass_z1'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb['smass_z1'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]/objs_pd_comb['vmass_z1'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = time_z1+np.zeros(np.sum((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]*1.5, marker = "*", linewidths = edgewidth, alpha = 0.5)
    ax14b.scatter(           objs_pd_comb['vmass_z0_5'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb['smass_z0_5'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]/objs_pd_comb['vmass_z0_5'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = time_z0_5+np.zeros(np.sum((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback'))), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]*1.5, marker = "*", linewidths = edgewidth, alpha = 0.5)    
    ax14b.scatter(           objs_pd_comb['vmass_z2'][objs_pd_comb['type']=='Satellite'],objs_pd_comb['smass_z2'][objs_pd_comb['type']=='Satellite']/objs_pd_comb['vmass_z2'][objs_pd_comb['type']=='Satellite'],c = time_z2+np.zeros(np.sum(objs_pd_comb['type']=='Satellite')), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite'], linewidths = edgewidth, alpha = 0.5)
    ax14b.scatter(           objs_pd_comb['vmass_z1'][objs_pd_comb['type']=='Satellite'],objs_pd_comb['smass_z1'][objs_pd_comb['type']=='Satellite']/objs_pd_comb['vmass_z1'][objs_pd_comb['type']=='Satellite'],c = time_z1+np.zeros(np.sum(objs_pd_comb['type']=='Satellite')), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite'], linewidths = edgewidth, alpha = 0.5)
    ax14b.scatter(           objs_pd_comb['vmass_z0_5'][objs_pd_comb['type']=='Satellite'],objs_pd_comb['smass_z0_5'][objs_pd_comb['type']=='Satellite']/objs_pd_comb['vmass_z0_5'][objs_pd_comb['type']=='Satellite'],c = time_z0_5+np.zeros(np.sum(objs_pd_comb['type']=='Satellite')), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite'], linewidths = edgewidth, alpha = 0.5)
    """
    dmax = 11; dmin = 8
    xaxis, data = bin_plt(np.log10(np.array(objs_pd_comb.loc[objs_pd_comb['vmass_z4']> 0, 'vmass_z4'].tolist())), np.array(objs_pd_comb.loc[objs_pd_comb['vmass_z4']> 0, 'smass_z4'].tolist())/np.array(objs_pd_comb.loc[objs_pd_comb['vmass_z4']> 0, 'vmass_z4'].tolist()), xmax = dmax, xmin = dmin, nbins = nbins)
    ax14b.fill_between(10**xaxis[data[:,4] >= 0], data[data[:,4] >= 0,0], data[data[:,4] >= 0,4], color=sm.to_rgba(time_z4), alpha=0.1, linewidth=0)
    ax14b.fill_between(10**xaxis[data[:,3] >= 0], data[data[:,3] >= 0,1], data[data[:,3] >= 0,3], color=sm.to_rgba(time_z4), alpha=0.1, linewidth=0)
    ax14b.plot(10**xaxis[data[:,2] >= 0], data[data[:,2] >= 0,2], color=sm.to_rgba(time_z4))    
    xaxis, data = bin_plt(np.log10(np.array(objs_pd_comb.loc[objs_pd_comb['vmass_z2']> 0, 'vmass_z2'].tolist())), np.array(objs_pd_comb.loc[objs_pd_comb['vmass_z2']> 0, 'smass_z2'].tolist())/np.array(objs_pd_comb.loc[objs_pd_comb['vmass_z2']> 0, 'vmass_z2'].tolist()), xmax = dmax, xmin = dmin, nbins = nbins)
    ax14b.fill_between(10**xaxis[data[:,4] >= 0], data[data[:,4] >= 0,0], data[data[:,4] >= 0,4], color=sm.to_rgba(time_z2), alpha=0.1, linewidth=0)
    ax14b.fill_between(10**xaxis[data[:,3] >= 0], data[data[:,3] >= 0,1], data[data[:,3] >= 0,3], color=sm.to_rgba(time_z2), alpha=0.1, linewidth=0)
    ax14b.plot(10**xaxis[data[:,2] >= 0], data[data[:,2] >= 0,2], color=sm.to_rgba(time_z2))
    xaxis, data = bin_plt(np.log10(np.array(objs_pd_comb.loc[objs_pd_comb['vmass_z1']> 0, 'vmass_z1'].tolist())), np.array(objs_pd_comb.loc[objs_pd_comb['vmass_z1']> 0, 'smass_z1'].tolist())/np.array(objs_pd_comb.loc[objs_pd_comb['vmass_z1']> 0, 'vmass_z1'].tolist()), xmax = dmax, xmin = dmin, nbins = nbins)
    ax14b.fill_between(10**xaxis[data[:,4] >= 0], data[data[:,4] >= 0,0], data[data[:,4] >= 0,4], color=sm.to_rgba(time_z1), alpha=0.1, linewidth=0)
    ax14b.fill_between(10**xaxis[data[:,3] >= 0], data[data[:,3] >= 0,1], data[data[:,3] >= 0,3], color=sm.to_rgba(time_z1), alpha=0.1, linewidth=0)
    ax14b.plot(10**xaxis[data[:,2] >= 0], data[data[:,2] >= 0,2], color=sm.to_rgba(time_z1))
    xaxis, data = bin_plt(np.log10(np.array(objs_pd_comb.loc[objs_pd_comb['vmass_z0_5']> 0, 'vmass_z0_5'].tolist())), np.array(objs_pd_comb.loc[objs_pd_comb['vmass_z0_5']> 0, 'smass_z0_5'].tolist())/np.array(objs_pd_comb.loc[objs_pd_comb['vmass_z0_5']> 0, 'vmass_z0_5'].tolist()), xmax = dmax, xmin = dmin, nbins = nbins)
    ax14b.fill_between(10**xaxis[data[:,4] >= 0], data[data[:,4] >= 0,0], data[data[:,4] >= 0,4], color=sm.to_rgba(time_z0_5), alpha=0.1, linewidth=0)
    ax14b.fill_between(10**xaxis[data[:,3] >= 0], data[data[:,3] >= 0,1], data[data[:,3] >= 0,3], color=sm.to_rgba(time_z0_5), alpha=0.1, linewidth=0)
    ax14b.plot(10**xaxis[data[:,2] >= 0], data[data[:,2] >= 0,2], color=sm.to_rgba(time_z0_5))
    xaxis, data = bin_plt(np.log10(np.array(objs_pd_comb.loc[objs_pd_comb['Mpeak_snap'] == 4096, 'mvir'].tolist())), np.array(objs_pd_comb.loc[objs_pd_comb['Mpeak_snap']==4096, 'Mstar_z0'].tolist())/np.array(objs_pd_comb.loc[objs_pd_comb['Mpeak_snap'] == 4096, 'mvir'].tolist()), xmax = dmax, xmin = dmin, nbins = nbins)
    ax14b.fill_between(10**xaxis[data[:,4] >= 0], data[data[:,4] >= 0,0], data[data[:,4] >= 0,4], color=sm.to_rgba(13.8), alpha=0.1, linewidth=0)
    ax14b.fill_between(10**xaxis[data[:,3] >= 0], data[data[:,3] >= 0,1], data[data[:,3] >= 0,3], color=sm.to_rgba(13.8), alpha=0.1, linewidth=0)
    pltz0, = ax14b.plot(10**xaxis[data[:,2] >= 0], data[data[:,2] >= 0,2], color=sm.to_rgba(13.8), label = "z = 0")    
    ax14b.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Central']/objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],c = (objs_pd_comb['Mpeak_snap'][objs_pd_comb['type']=='Central'])*13.8/4096, cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax14b.scatter(objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb[mstar_key][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]/objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = (objs_pd_comb['Mpeak_snap'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')])*13.8/4096, cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]*1.5, marker = "*", linewidths = edgewidth)    
    ax14b.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Satellite']/objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],c = (objs_pd_comb['Mpeak_snap'][objs_pd_comb['type']=='Satellite'])*13.8/4096, cmap = cmx, norm = cNorm,edgecolor = 'k', marker = "*", s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    #ax14.scatter(readdata['M200']*1e10,readdata['M_*']*1e7/(readdata['M200']*1e10),marker = '+',c = 'k')
    ax14b.set_xscale('log')
    ax14b.set_yscale('log')
    ax14b.set_ylabel(mstar_ratio_axes_label)
    ax14b.set_xlabel(r'M$_{\mathrm{vir, peak}}$/M$_\odot$')
    ax14b.axis([1e8, 2e11, 1e-6, 0.05])
    legend1 = ax14b.legend([cen_plt, sat_plt],['Central','Satellite/Backsplash'],scatterpoints = 1,facecolor = 'white',loc = 0,framealpha = 0,frameon = False)
    #legend2 = ax14b.legend([read,nadler,jethwa],['Read+ 2017','Nadler+ 2020','Jethwa+ 2018'],facecolor = 'white',loc = 0,framealpha = 0,frameon = False)
    ax14b.add_artist(legend1)
    #ax14b.add_artist(legend2)
    ax14b.label_outer()
    fig14.tight_layout()
    
    sm = plt.cm.ScalarMappable(cmap=cmx, norm=cNorm)
    cbar = plt.colorbar(sm, location="right",ax=[ax14a,ax14b],aspect = 40) #pad=0.04)
    #cbar = mpl.colorbar.ColorbarBase(ax14sub, cmap=cmx, norm=cNorm, location="right",ax=[ax14a,ax14b]) #pad=0.04)
    cbar.set_label(r'Time of Peak M$_{\mathrm{vir}}$ [Gyr]') #, rotation=0)
    #fig14.tight_layout()
    #sm = plt.cm.ScalarMappable(cmap=cmx, norm=cNorm)
    #sm.set_array([])
    #cbar = plt.colorbar(sm, location="top", ax=axs14.ravel().tolist(), fraction = 0.1) #pad=0.04)
    #cbar.set_label(r'D$_{\mathrm{massive}}$ [kpc]') #, rotation=0)
    #cbar.ax.get_xaxis().labelpad = 15
    fig14.show()
    fig14.savefig(outfile_base + '_SMHM2_tmax_Mpeak.png',dpi = dpi)

    
    ########################## Fraction of Stellar and Virial mass loss ############################33
    plt.close('all')
    fig15 = plt.figure(15,figsize=(plt_width,plt_width*aspect_ratio))
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax15 = fig15.add_subplot(gs[0])
    #ax15sub = fig15.add_subplot(gs[1])
    cNorm  = colors.LogNorm(vmin=10, vmax = 2e3) 
    ax15.scatter(objs_pd_comb['mass'][objs_pd_comb['type']=='Central']/objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],objs_pd_comb['Mstar_z0'][objs_pd_comb['type']=='Central']/objs_pd_comb['Mstar_Mpeak'][objs_pd_comb['type']=='Central'],c = (objs_pd_comb['min_dist'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central']*1, linewidths = edgewidth, alpha = 0.5)
    ax15.scatter(objs_pd_comb['mass'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]/objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb['Mstar_z0'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]/objs_pd_comb['Mstar_Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = (objs_pd_comb['min_dist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]*1, linewidths = edgewidth, alpha = 0.5)
    ax15.scatter(objs_pd_comb['mass'][objs_pd_comb['type']=='Satellite']/objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],objs_pd_comb['Mstar_z0'][objs_pd_comb['type']=='Satellite']/objs_pd_comb['Mstar_Mpeak'][objs_pd_comb['type']=='Satellite'],c = (objs_pd_comb['min_dist'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k', marker = '*', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth, alpha = 0.5)
    #ax15.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax15.set_xscale('log')
    ax15.set_yscale('log')
    ax15.set_ylabel(r'M$_*$/M$_{\mathrm{*, peak}}$')
    ax15.set_xlabel(r'M$_{\mathrm{vir}}$/M$_{\mathrm{vir, peak}}$')
    fig15.tight_layout()
    ax15.legend([cen_plt,sat_plt],['Central','Satellite'],scatterpoints = 1,facecolor = 'white',loc = 0,framealpha = 0,frameon = False)
    sm = plt.cm.ScalarMappable(cmap=cmx, norm=cNorm)
    cbar = plt.colorbar(sm, ax=[ax15],aspect = 20) #pad=0.04)
    cbar.set_label(r'Min(D$_{\mathrm{massive}}$) [kpc]') #, rotation=0)
    #cb = mpl.colorbar.ColorbarBase(ax15sub, cmap=cmx, norm=cNorm)
    #cb.set_label(r"D$_{\mathrm{massive}}$ [kpc]")
    fig15.show()
    fig15.set_size_inches(plt_width,plt_width*aspect_ratio)
    fig15.show()
    fig15.savefig(outfile_base + '_fracMassloss_rClosest.png',dpi = dpi)

    plt.close('all')
    fig15 = plt.figure(15,figsize=(plt_width,plt_width*aspect_ratio))
    gs = fig15.add_gridspec(1, 2,  width_ratios=(8, 1.5),wspace=0.05) #, hspace=0.05)
    #                      left=0.1, right=0.9, bottom=0.1, top=0.9)
    ax15 = fig15.add_subplot(gs[0,0])
    #axs15 = gs.subplots(sharey = True) #, constrained_layout=True) 
    #axs15 = axs15.flatten()
    #ax15 = axs15[0]
    #ax15_hist = axs15[1]
    
    ax15.scatter(objs_pd_comb['min_dist'][objs_pd_comb['type']=='Central'],np.log10((objs_pd_comb['Mhalo_z0'][objs_pd_comb['type']=='Central']/objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central']).tolist()),facecolor = 'k', edgecolor = 'k', s = ((objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'])**0.33/10*ms_scale).tolist(), linewidths = edgewidth, alpha = 0.3)
    f_bary = ax15.scatter(objs_pd_comb['min_dist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],np.log10((objs_pd_comb['Mhalo_z0'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]/objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]).tolist()), facecolor = 'k', edgecolor = 'k', marker = '*', s =  ((objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')])**0.33/10*ms_scale*1.5).tolist(), linewidths = edgewidth, alpha = 0.3)
    ax15.scatter(objs_pd_comb['min_dist'][objs_pd_comb['type']=='Satellite'],np.log10((objs_pd_comb['Mhalo_z0'][objs_pd_comb['type']=='Satellite']/objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite']).tolist()),facecolor = 'k', edgecolor = 'k', marker = '*', s = ((objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'])**0.33/10*ms_scale*1.5).tolist(), linewidths = edgewidth, alpha = 0.3)

    cent = ax15.scatter(objs_pd_comb['min_dist'][objs_pd_comb['type']=='Central'],np.log10((objs_pd_comb['Mhalo_z0'][objs_pd_comb['type']=='Central']/objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central']).tolist()),facecolor = 'none', edgecolor = 'k', s = ((objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'])**0.33/10*ms_scale).tolist(), linewidths = edgewidth)
    ax15.scatter(objs_pd_comb['min_dist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],np.log10((objs_pd_comb['Mhalo_z0'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]/objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]).tolist()), facecolor = 'none', edgecolor = 'k', marker = '*', s =  ((objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')])**0.33/10*ms_scale*1.5).tolist(), linewidths = edgewidth)
    sat = ax15.scatter(objs_pd_comb['min_dist'][objs_pd_comb['type']=='Satellite'],np.log10((objs_pd_comb['Mhalo_z0'][objs_pd_comb['type']=='Satellite']/objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite']).tolist()),facecolor = 'none', edgecolor = 'k', marker = '*', s = ((objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'])**0.33/10*ms_scale*1.5).tolist(), linewidths = edgewidth)
    
    ax15.scatter(objs_pd_comb['min_dist'][objs_pd_comb['type']=='Central'],np.log10((objs_pd_comb['Mstar_z0'][objs_pd_comb['type']=='Central']/objs_pd_comb['Mstar_Mpeak'][objs_pd_comb['type']=='Central']).tolist()),facecolor = 'none', edgecolor = 'r', s = ((objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'])**0.33/10*ms_scale).tolist(), linewidths = edgewidth)
    f_star = ax15.scatter(objs_pd_comb['min_dist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],np.log10((objs_pd_comb['Mstar_z0'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]/objs_pd_comb['Mstar_Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]).tolist()), facecolor = 'none', edgecolor = 'r', marker = '*', s =  ((objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')])**0.33/10*ms_scale*1.5).tolist(), linewidths = edgewidth)
    ax15.scatter(objs_pd_comb['min_dist'][objs_pd_comb['type']=='Satellite'],np.log10((objs_pd_comb['Mstar_z0'][objs_pd_comb['type']=='Satellite']/objs_pd_comb['Mstar_Mpeak'][objs_pd_comb['type']=='Satellite']).tolist()),facecolor = 'none', edgecolor = 'r', marker = '*', s = ((objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'])**0.33/10*ms_scale*1.5).tolist(), linewidths = edgewidth)

    ax15.errorbar(objs_pd_comb['min_dist'][objs_pd_comb['Mstar_Mpeak'] == 0],objs_pd_comb['min_dist'][objs_pd_comb['Mstar_Mpeak'] == 0]*0+0.8,yerr = 0.2,fmt='none',ecolor = 'r',lolims = True, capsize = 3)
    
    ax15.set_xscale('log')
    ax15.set_ylabel(r'Log(M$_{z = 0}$/M$_{\mathrm{peak}}$)')
    ax15.set_xlabel(r"Minimum distance to massive galaxy [kpc]")    
    ax15.legend([f_star,f_bary],[r'Stellar Mass',r'Virial Mass'],loc = 4)
    #ax15.legend([cent,sat],[r'Central',r'Satellite/Backsplash'],loc = 3)
    
    binwidth = 0.3
    x_star = objs_pd_comb['Mstar_z0']/objs_pd_comb['Mstar_Mpeak']
    x_vir = objs_pd_comb['Mhalo_z0']/objs_pd_comb['Mpeak']
    min_x = -3.75 # np.log10(np.max(x_star,x_vir))
    max_x = 1.15 # np.log10(np.max(x_star,x_vir))
    ax15.axis([10,4e3 ,min_x, max_x]) # 1.3e3
    #fig15.tight_layout()

    ax15_hist = fig15.add_subplot(gs[0,1], sharey=ax15)
    #xymax = np.max(np.abs(x_star,x_vir))
    #lim = (int(xym/binwidth) + 1) * binwidth
    bins = np.arange(min_x, max_x + binwidth, binwidth)
    ax15_hist.hist(np.log10(x_vir), bins = bins, facecolor = 'none', edgecolor = 'k', orientation='horizontal')
    ax15_hist.hist(np.log10(x_star), bins = bins, facecolor = 'none', edgecolor = 'r', orientation='horizontal')
    ax15_hist.axis('off')
    #fig15.tight_layout()
    #ax15.legend([cen_plt,sat_plt],['Central','Satellite'],scatterpoints = 1,facecolor = 'white',loc = 0,framealpha = 0,frameon = False)
    fig15.show()
    fig15.set_size_inches(plt_width,plt_width*aspect_ratio)
    fig15.show()    
    fig15.savefig(outfile_base + '_fracMassloss_rMassGal.png',dpi = dpi,bbox_inches='tight')  
    
    plt.close('all')
    fig15 = plt.figure(15,figsize=(plt_width,plt_width*aspect_ratio))
    gs = fig15.add_gridspec(1, 2,  width_ratios=(8, 1.5),wspace=0.05) #, hspace=0.05)
    #                      left=0.1, right=0.9, bottom=0.1, top=0.9)
    ax15 = fig15.add_subplot(gs[0,0])
    #axs15 = gs.subplots(sharey = True) #, constrained_layout=True) 
    #axs15 = axs15.flatten()
    #ax15 = axs15[0]
    #ax15_hist = axs15[1]
    
    ax15.scatter(objs_pd_comb['min_dist'][objs_pd_comb['type']=='Central'],np.log10((objs_pd_comb['Mhalo_z0'][objs_pd_comb['type']=='Central']/objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central']).tolist()),facecolor = 'k', edgecolor = 'k', s = ((objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'])**0.33/10*ms_scale).tolist(), linewidths = edgewidth, alpha = 0.3)
    ax15.scatter(objs_pd_comb['min_dist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],np.log10((objs_pd_comb['Mhalo_z0'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]/objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]).tolist()), facecolor = 'k', edgecolor = 'k', s =  ((objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')])**0.33/10*ms_scale*1.5).tolist(), linewidths = edgewidth, alpha = 0.3)
    ax15.scatter(objs_pd_comb['min_dist'][objs_pd_comb['type']=='Satellite'],np.log10((objs_pd_comb['Mhalo_z0'][objs_pd_comb['type']=='Satellite']/objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite']).tolist()),facecolor = 'k', edgecolor = 'k', marker = '*', s = ((objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'])**0.33/10*ms_scale*1.5).tolist(), linewidths = edgewidth, alpha = 0.3)

    cent = ax15.scatter(objs_pd_comb['min_dist'][objs_pd_comb['type']=='Central'],np.log10((objs_pd_comb['Mhalo_z0'][objs_pd_comb['type']=='Central']/objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central']).tolist()),facecolor = 'none', edgecolor = 'k', s = ((objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'])**0.33/10*ms_scale).tolist(), linewidths = edgewidth)
    f_bary = ax15.scatter(objs_pd_comb['min_dist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],np.log10((objs_pd_comb['Mhalo_z0'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]/objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]).tolist()), facecolor = 'none', edgecolor = 'k', s =  ((objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')])**0.33/10*ms_scale*1.5).tolist(), linewidths = edgewidth)
    sat = ax15.scatter(objs_pd_comb['min_dist'][objs_pd_comb['type']=='Satellite'],np.log10((objs_pd_comb['Mhalo_z0'][objs_pd_comb['type']=='Satellite']/objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite']).tolist()),facecolor = 'none', edgecolor = 'k', marker = '*', s = ((objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'])**0.33/10*ms_scale*1.5).tolist(), linewidths = edgewidth)
    
    ax15.scatter(objs_pd_comb['min_dist'][objs_pd_comb['type']=='Central'],np.log10((objs_pd_comb['Mstar_z0'][objs_pd_comb['type']=='Central']/objs_pd_comb['Mstar_Mpeak'][objs_pd_comb['type']=='Central']).tolist()),facecolor = 'none', edgecolor = 'r', s = ((objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'])**0.33/10*ms_scale).tolist(), linewidths = edgewidth)
    f_star = ax15.scatter(objs_pd_comb['min_dist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],np.log10((objs_pd_comb['Mstar_z0'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]/objs_pd_comb['Mstar_Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]).tolist()), facecolor = 'none', edgecolor = 'r', s =  ((objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')])**0.33/10*ms_scale*1.5).tolist(), linewidths = edgewidth)
    ax15.scatter(objs_pd_comb['min_dist'][objs_pd_comb['type']=='Satellite'],np.log10((objs_pd_comb['Mstar_z0'][objs_pd_comb['type']=='Satellite']/objs_pd_comb['Mstar_Mpeak'][objs_pd_comb['type']=='Satellite']).tolist()),facecolor = 'none', edgecolor = 'r', marker = '*', s = ((objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'])**0.33/10*ms_scale*1.5).tolist(), linewidths = edgewidth)
    ax15.set_xscale('log')
    ax15.set_ylabel(r'Log(M$_{z = 0}$/M$_{\mathrm{peak}}$)')
    ax15.set_xlabel(r"Minimum distance to massive galaxy [kpc]")    
    ax15.legend([fstar,fbary],[r'Stellar Mass',r'Virial Mass'],loc = 4)
    #ax15.legend([cent,sat],[r'Central',r'Satellite/Backsplash'],loc = 3)
    
    binwidth = 0.25
    x_star = objs_pd_comb['Mstar_z0']/objs_pd_comb['Mstar_Mpeak']
    x_vir = objs_pd_comb['Mhalo_z0']/objs_pd_comb['Mpeak']
    min_x = -4 # np.log10(np.max(x_star,x_vir))
    max_x = 2 # np.log10(np.max(x_star,x_vir))
    ax15.axis([10,1.5e3,min_x, max_x])
    #fig15.tight_layout()

    ax15_hist = fig15.add_subplot(gs[0,1], sharey=ax15)
    #xymax = np.max(np.abs(x_star,x_vir))
    #lim = (int(xym/binwidth) + 1) * binwidth
    bins = np.arange(min_x, max_x + binwidth, binwidth)
    ax15_hist.hist(np.log10(x_vir), bins = bins, facecolor = 'none', edgecolor = 'k', orientation='horizontal')
    ax15_hist.hist(np.log10(x_star), bins = bins, facecolor = 'none', edgecolor = 'r', orientation='horizontal')
    ax15_hist.axis('off')
    #fig15.tight_layout()
    #ax15.legend([cen_plt,sat_plt],['Central','Satellite'],scatterpoints = 1,facecolor = 'white',loc = 0,framealpha = 0,frameon = False)
    fig15.show()
    fig15.set_size_inches(plt_width,plt_width*aspect_ratio)
    fig15.show()    
    fig15.savefig(outfile_base + '_fracMassloss_rClosest.png',dpi = dpi,bbox_inches='tight')   
        
    ###### tau 90 vs distance
    plt.close('all')
    plt.clf()
    fig16 = plt.figure(1,figsize = (plt_width,plt_width*aspect_ratio*1.2))
    fig16.set_size_inches(plt_width,plt_width*aspect_ratio*1.2)
    fig16.clear()
    gs = gridspec.GridSpec(ncols = 1,nrows = 2, figure=fig16, height_ratios=[1,15])
    ax16 = fig16.add_subplot(gs[1])
    ax16sub = fig16.add_subplot(gs[0])

    cmx = plt.get_cmap("cool_r")
    dcmap = discrete_cmap(256, cmx)
    cNorm  = colors.Normalize(vmin=0, vmax = 14)
    #ax16.scatter(objs_rom['massiveDist'],objs_rom['mstar'],c = objs_rom['tau90'],cmap = cmx, norm = cNorm,alpha = 0.2, s = markersize*0.5)
    quench_cond = np.array(list(zip(*objs_pd_comb['SFR']))[0])/objs_pd_comb['M_star'] < 1e-11
    sf_cond =  np.array(list(zip(*objs_pd_comb['SFR']))[0])/objs_pd_comb['M_star'] >= 1e-11
    sf = ax16.scatter(objs_pd_comb[sf_cond * (objs_pd_comb['type']=='Central').tolist()]['massiveDist'],objs_pd_comb[mstar_key][sf_cond * (objs_pd_comb['type']=='Central').tolist()],s = ((objs_pd_comb['Mpeak'][sf_cond * (objs_pd_comb['type']=='Central').tolist()])**0.33/15*ms_scale).tolist(), c = objs_pd_comb['tau90_vir'][sf_cond * (objs_pd_comb['type']=='Central').tolist()], cmap = cmx, norm = cNorm,edgecolor = 'k', linewidths = edgewidth)
    ax16.scatter(objs_pd_comb[sf_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')).tolist()]['massiveDist'],objs_pd_comb[mstar_key][sf_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')).tolist()],s = ((objs_pd_comb['Mpeak'][sf_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')).tolist()])**0.33/15*ms_scale*1.5).tolist(), c = objs_pd_comb['tau90_vir'][sf_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')).tolist()], cmap = cmx, norm = cNorm,edgecolor = 'k', linewidths = edgewidth, marker = "*")
    ax16.scatter(objs_pd_comb[sf_cond * (objs_pd_comb['type']=='Satellite').tolist()]['massiveDist'],objs_pd_comb[mstar_key][sf_cond * (objs_pd_comb['type']=='Satellite').tolist()],s = ((objs_pd_comb['Mpeak'][sf_cond * (objs_pd_comb['type']=='Satellite').tolist()])**0.33/15*ms_scale*1.5).tolist(), c = objs_pd_comb['tau90_vir'][sf_cond * (objs_pd_comb['type']=='Satellite').tolist()], cmap = cmx, norm = cNorm,edgecolor = 'k', linewidths = edgewidth, marker = "*")
    
    q = ax16.scatter(objs_pd_comb[quench_cond * (objs_pd_comb['type']=='Central').tolist()]['massiveDist'],objs_pd_comb[mstar_key][quench_cond * (objs_pd_comb['type']=='Central').tolist()],s = ((objs_pd_comb['Mpeak'][quench_cond * (objs_pd_comb['type']=='Central').tolist()])**0.33/15*ms_scale).tolist(), edgecolor = dcmap(cNorm(objs_pd_comb['tau90_vir'][quench_cond * (objs_pd_comb['type']=='Central').tolist()])), cmap = cmx, norm = cNorm,facecolor = 'none', linewidths = edgewidth*2)
    ax16.scatter(objs_pd_comb[quench_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')).tolist()]['massiveDist'],objs_pd_comb[mstar_key][quench_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')).tolist()],s = ((objs_pd_comb['Mpeak'][quench_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')).tolist()])**0.33/15*ms_scale*1.5).tolist(), edgecolor = dcmap(cNorm(objs_pd_comb['tau90_vir'][quench_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')).tolist()])), cmap = cmx, norm = cNorm,facecolor = 'none', linewidths = edgewidth*2, marker = "*")
    ax16.scatter(objs_pd_comb[quench_cond * (objs_pd_comb['type']=='Satellite').tolist()]['massiveDist'],objs_pd_comb[mstar_key][quench_cond * (objs_pd_comb['type']=='Satellite').tolist()],s = ((objs_pd_comb['Mpeak'][quench_cond * (objs_pd_comb['type']=='Satellite').tolist()])**0.33/15*ms_scale*1.5).tolist(), edgecolor = dcmap(cNorm(objs_pd_comb['tau90_vir'][quench_cond * (objs_pd_comb['type']=='Satellite').tolist()])), cmap = cmx, norm = cNorm,facecolor = 'none', linewidths = 2*edgewidth, marker = "*")

    #plt.scatter(objs_pd_e['h1dist'],objs_pd_e['M_star'])
    lgnd = ax16.legend([q,sf],['Quenched','Star forming'],scatterpoints = 1,facecolor = 'white',loc = 3,framealpha = 0,frameon = False)
    lgnd.legendHandles[0]._sizes = [markersize]
    lgnd.legendHandles[1]._sizes = [markersize]
    ax16.set_xscale('log')
    ax16.set_yscale('log')
    ax16.axis([40, 8e3, 1e2, 1e10])
    ax16.set_ylabel(mstar_axes_label)
    ax16.set_xlabel(r'D$_{\mathrm{massive}}$ [kpc]')
    sm = plt.cm.ScalarMappable(cmap=cmx, norm=cNorm)
    #ax16sub = fig16.add_subplot([0.15,0.87,0.35,0.03])
    #cb = plt.colorbar(sm, ax=ax16sub, orientation='horizontal',aspect = 20)    
    cb = mpl.colorbar.ColorbarBase(ax16sub, cmap=cmx, norm=cNorm, orientation='horizontal')
    cb.set_label(r"$\tau_{90}$ (Gyr)")
    #fig16.subplots_adjust(hspace = 0.3)
    fig16.tight_layout()
    fig16.subplots_adjust(hspace = 0.4)
    fig16.show()
    fig16.savefig(outfile_base + '_rMassGal_smass_t90.png',dpi = dpi)    

    ###### tau 90 vs distance
    plt.close('all')
    plt.clf()
    fig16 = plt.figure(1,figsize = (plt_width,plt_width*aspect_ratio*1.2))
    fig16.set_size_inches(plt_width,plt_width*aspect_ratio*1.2)
    fig16.clear()
    gs = gridspec.GridSpec(ncols = 1,nrows = 2, figure=fig16, height_ratios=[1,15])
    ax16 = fig16.add_subplot(gs[1])
    ax16sub = fig16.add_subplot(gs[0])

    cmx = plt.get_cmap("cool_r")
    dcmap = discrete_cmap(256, cmx)
    cNorm  = colors.Normalize(vmin=0, vmax = 14)
    #ax16.scatter(objs_rom['min_dist'],objs_rom['mstar'],c = objs_rom['tau90'],cmap = cmx, norm = cNorm,alpha = 0.2, s = markersize*0.5)
    quench_cond = np.array(list(zip(*objs_pd_comb['SFR']))[0])/objs_pd_comb['M_star'] < 1e-11
    sf_cond =  np.array(list(zip(*objs_pd_comb['SFR']))[0])/objs_pd_comb['M_star'] >= 1e-11
    sf = ax16.scatter(objs_pd_comb[sf_cond * (objs_pd_comb['type']=='Central').tolist()]['min_dist'],objs_pd_comb[mstar_key][sf_cond * (objs_pd_comb['type']=='Central').tolist()],s = ((objs_pd_comb['Mpeak'][sf_cond * (objs_pd_comb['type']=='Central').tolist()])**0.33/15*ms_scale).tolist(), c = objs_pd_comb['tau90_vir'][sf_cond * (objs_pd_comb['type']=='Central').tolist()], cmap = cmx, norm = cNorm,edgecolor = 'k', linewidths = edgewidth)
    ax16.scatter(objs_pd_comb[sf_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')).tolist()]['min_dist'],objs_pd_comb[mstar_key][sf_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')).tolist()],s = ((objs_pd_comb['Mpeak'][sf_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')).tolist()])**0.33/15*ms_scale*1.5).tolist(), c = objs_pd_comb['tau90_vir'][sf_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')).tolist()], cmap = cmx, norm = cNorm,edgecolor = 'k', linewidths = edgewidth, marker = "*")
    ax16.scatter(objs_pd_comb[sf_cond * (objs_pd_comb['type']=='Satellite').tolist()]['min_dist'],objs_pd_comb[mstar_key][sf_cond * (objs_pd_comb['type']=='Satellite').tolist()],s = ((objs_pd_comb['Mpeak'][sf_cond * (objs_pd_comb['type']=='Satellite').tolist()])**0.33/15*ms_scale*1.5).tolist(), c = objs_pd_comb['tau90_vir'][sf_cond * (objs_pd_comb['type']=='Satellite').tolist()], cmap = cmx, norm = cNorm,edgecolor = 'k', linewidths = edgewidth, marker = "*")
    
    q = ax16.scatter(objs_pd_comb[quench_cond * (objs_pd_comb['type']=='Central').tolist()]['min_dist'],objs_pd_comb[mstar_key][quench_cond * (objs_pd_comb['type']=='Central').tolist()],s = ((objs_pd_comb['Mpeak'][quench_cond * (objs_pd_comb['type']=='Central').tolist()])**0.33/15*ms_scale).tolist(), edgecolor = dcmap(cNorm(objs_pd_comb['tau90_vir'][quench_cond * (objs_pd_comb['type']=='Central').tolist()])), cmap = cmx, norm = cNorm,facecolor = 'none', linewidths = edgewidth*2)
    ax16.scatter(objs_pd_comb[quench_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')).tolist()]['min_dist'],objs_pd_comb[mstar_key][quench_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')).tolist()],s = ((objs_pd_comb['Mpeak'][quench_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')).tolist()])**0.33/15*ms_scale*1.5).tolist(), edgecolor = dcmap(cNorm(objs_pd_comb['tau90_vir'][quench_cond * ((objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')).tolist()])), cmap = cmx, norm = cNorm,facecolor = 'none', linewidths = edgewidth*2, marker = "*")
    ax16.scatter(objs_pd_comb[quench_cond * (objs_pd_comb['type']=='Satellite').tolist()]['min_dist'],objs_pd_comb[mstar_key][quench_cond * (objs_pd_comb['type']=='Satellite').tolist()],s = ((objs_pd_comb['Mpeak'][quench_cond * (objs_pd_comb['type']=='Satellite').tolist()])**0.33/15*ms_scale*1.5).tolist(), edgecolor = dcmap(cNorm(objs_pd_comb['tau90_vir'][quench_cond * (objs_pd_comb['type']=='Satellite').tolist()])), cmap = cmx, norm = cNorm,facecolor = 'none', linewidths = 2*edgewidth, marker = "*")

    #plt.scatter(objs_pd_e['h1dist'],objs_pd_e['M_star'])
    lgnd = ax16.legend([q,sf],['Quenched','Star forming'],scatterpoints = 1,facecolor = 'white',loc = 3,framealpha = 0,frameon = False)
    lgnd.legendHandles[0]._sizes = [markersize]
    lgnd.legendHandles[1]._sizes = [markersize]
    ax16.set_xscale('log')
    ax16.set_yscale('log')
    ax16.axis([10, 2e3, 1e2, 1e10])
    ax16.set_ylabel(mstar_axes_label)
    ax16.set_xlabel(r'Min(D$_{\mathrm{massive}}$) [kpc]')
    sm = plt.cm.ScalarMappable(cmap=cmx, norm=cNorm)
    #ax16sub = fig16.add_subplot([0.15,0.87,0.35,0.03])
    #cb = plt.colorbar(sm, ax=ax16sub, orientation='horizontal',aspect = 20)    
    cb = mpl.colorbar.ColorbarBase(ax16sub, cmap=cmx, norm=cNorm, orientation='horizontal')
    cb.set_label(r"$\tau_{90}$ (Gyr)")
    #fig16.subplots_adjust(hspace = 0.3)
    fig16.tight_layout()
    fig16.subplots_adjust(hspace = 0.4)
    fig16.show()
    fig16.savefig(outfile_base + '_rClosest_smass_t90.png',dpi = dpi)

    plt.close('all')
    #fig17.clear()
    #fig2 = plt.figure(2, figsize=(plt_width*2,plt_width*aspect_ratio))
    #gs = gridspec.GridSpec(1,3,width_ratios=[15,15,1])
    #ax2a = fig17.add_subplot(gs[0])
    #ax2b = fig17.add_subplot(gs[1])
    #ax2sub = fig17.add_subplot(gs[2])
    fig17 = plt.figure(figsize=(plt_width,plt_width*aspect_ratio*2))
    gs = fig17.add_gridspec(2, hspace=0)
    axs17 = gs.subplots(sharex = True) #, constrained_layout=True)
    axs17 = axs17.flatten()
    ax17a = axs17[0]
    ax17b = axs17[1]
    #cNorm  = colors.LogNorm(vmin=1e2, vmax = 1e8)
    cNorm  = colors.Normalize(5, vmax = 20)
    cmx = plt.get_cmap("viridis") 
    ax17a.fill_between(read_abunmatch['M200'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_low'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_high'][read_abunmatch['M200']> 1e9],color = 'grey',linewidth=0,alpha = 0.6)
    ax17a.fill_between(upper_fid_Mhalo, np.interp( upper_fid_Mhalo, lower_fid_Mhalo, lower_fid_Mstar ), upper_fid_Mstar, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3)#Nadler outer
    ax17a.fill_between(upper_fid_inner_Mhalo, np.interp( upper_fid_inner_Mhalo, lower_fid_inner_Mhalo, lower_fid_inner_Mstar ),upper_fid_inner_Mstar, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3) #Nadler inner
    ax17a.fill_between(Jethwa_upper_outer_Mhalo, np.interp( Jethwa_upper_outer_Mhalo, Jethwa_lower_outer_Mhalo, Jethwa_lower_outer_Mstar ),Jethwa_upper_outer_Mstar, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3)#Jethwa outer
    ax17a.fill_between(Jethwa_upper_inner_Mhalo,np.interp( Jethwa_upper_inner_Mhalo, Jethwa_lower_inner_Mhalo, Jethwa_lower_inner_Mstar ), Jethwa_upper_inner_Mstar, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3) #Jethwa inner
    read, = ax17b.plot( [], [],  color="grey", linestyle="solid", linewidth=4) #, label="Read+ 2017" )
    nadler, = ax17b.plot( [], [],  color="#33669A", linestyle="solid", linewidth=8) #, label="Nadler+ 2020" )
    jethwa, = ax17b.plot( [], [],  color="#66CC66", linestyle="solid", linewidth=8) #, label="Jethwa+ 2018" )     
    cen_plt = ax17a.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Central'],c = (objs_pd_comb['concentration'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax17a.scatter(objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb[mstar_key][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = (objs_pd_comb['concentration'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]*1.5, marker = "*", linewidths = edgewidth)    
    sat_plt =  ax17a.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Satellite'],c = (objs_pd_comb['concentration'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k', marker = "*", s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    #ax17.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax17a.set_xscale('log')
    ax17a.set_yscale('log')
    ax17a.set_ylabel(mstar_axes_label)
    #ax17a.set_xlabel(r'M$_{\mathrm{vir, peak}}$/M$_\odot$')
    ax17a.axis([1e8, 2e11, 2e2, 5e9])
    #legend1 = ax17a.legend([cen_plt,sat_plt],['Central','Satellite'],scatterpoints = 1,facecolor = 'white',loc = 0,framealpha = 0,frameon = False)
    legend2 = ax17a.legend([read,nadler,jethwa],['Read+ 2017','Nadler+ 2020','Jethwa+ 2018'],facecolor = 'white',loc = 4,framealpha = 0,frameon = False)
    #ax17.add_artist(legend1)
    ax17a.add_artist(legend2)
    ax17a.label_outer()

    ax17b.fill_between(read_abunmatch['M200'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_low'][read_abunmatch['M200']> 1e9]/read_abunmatch['M200'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_high'][read_abunmatch['M200']> 1e9]/read_abunmatch['M200'][read_abunmatch['M200']> 1e9],color = 'grey',linewidth=0,alpha = 0.6)
    ax17b.fill_between(upper_fid_Mhalo, np.interp( upper_fid_Mhalo, lower_fid_Mhalo, lower_fid_Mstar )/upper_fid_Mhalo, upper_fid_Mstar/upper_fid_Mhalo, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3)#Nadler outer
    ax17b.fill_between(upper_fid_inner_Mhalo, np.interp( upper_fid_inner_Mhalo, lower_fid_inner_Mhalo, lower_fid_inner_Mstar)/upper_fid_inner_Mhalo,upper_fid_inner_Mstar/upper_fid_inner_Mhalo, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3) #Nadler inner
    ax17b.fill_between(Jethwa_upper_outer_Mhalo, np.interp( Jethwa_upper_outer_Mhalo, Jethwa_lower_outer_Mhalo, Jethwa_lower_outer_Mstar )/Jethwa_upper_outer_Mhalo,Jethwa_upper_outer_Mstar/Jethwa_upper_outer_Mhalo, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3)#Jethwa outer
    ax17b.fill_between(Jethwa_upper_inner_Mhalo,np.interp( Jethwa_upper_inner_Mhalo, Jethwa_lower_inner_Mhalo, Jethwa_lower_inner_Mstar )/Jethwa_upper_inner_Mhalo, Jethwa_upper_inner_Mstar/Jethwa_upper_inner_Mhalo, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3) #Jethwa inner
    read, = ax17b.plot( [], [],  color="grey", linestyle="solid", linewidth=4) #, label="Read+ 2017" )
    nadler, = ax17b.plot( [], [],  color="#33669A", linestyle="solid", linewidth=8) #, label="Nadler+ 2020" )
    jethwa, = ax17b.plot( [], [],  color="#66CC66", linestyle="solid", linewidth=8) #, label="Jethwa+ 2018" )        
    ax17b.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Central']/objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],c = (objs_pd_comb['concentration'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax17b.scatter(objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb[mstar_key][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]/objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = (objs_pd_comb['concentration'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]*1.5, marker = "*", linewidths = edgewidth)    
    ax17b.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Satellite']/objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],c = (objs_pd_comb['concentration'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k', marker = "*", s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    #ax17.scatter(readdata['M200']*1e10,readdata['M_*']*1e7/(readdata['M200']*1e10),marker = '+',c = 'k')
    ax17b.set_xscale('log')
    ax17b.set_yscale('log')
    ax17b.set_ylabel(mstar_ratio_axes_label)
    ax17b.set_xlabel(r'M$_{\mathrm{vir, peak}}$/M$_\odot$')
    ax17b.axis([1e8, 2e11, 1e-6, 0.05])
    legend1 = ax17b.legend([cen_plt, sat_plt],['Central','Satellite/Backsplash'],scatterpoints = 1,facecolor = 'white',loc = 0,framealpha = 0,frameon = False)
    #legend2 = ax17b.legend([read,nadler,jethwa],['Read+ 2017','Nadler+ 2020','Jethwa+ 2018'],facecolor = 'white',loc = 0,framealpha = 0,frameon = False)
    ax17b.add_artist(legend1)
    #ax17b.add_artist(legend2)
    ax17b.label_outer()
    fig17.tight_layout()
    
    sm = plt.cm.ScalarMappable(cmap=cmx, norm=cNorm)
    cbar = plt.colorbar(sm, location="right",ax=[ax17a,ax17b],aspect = 40) #pad=0.04)
    #cbar = mpl.colorbar.ColorbarBase(ax17sub, cmap=cmx, norm=cNorm, location="right",ax=[ax17a,ax17b]) #pad=0.04)
    cbar.set_label(r'Concentration') #, rotation=0)
    #fig17.tight_layout()
    #sm = plt.cm.ScalarMappable(cmap=cmx, norm=cNorm)
    #sm.set_array([])
    #cbar = plt.colorbar(sm, location="top", ax=axs17.ravel().tolist(), fraction = 0.1) #pad=0.04)
    #cbar.set_label(r'D$_{\mathrm{massive}}$ [kpc]') #, rotation=0)
    #cbar.ax.get_xaxis().labelpad = 15
    fig17.show()
    fig17.savefig(outfile_base + '_SMHM2_concent_Mpeak.png',dpi = dpi)

    plt.close('all') #339924903
    fig18 = plt.figure(2,figsize=(plt_width,plt_width*aspect_ratio))
    fig18.set_size_inches(plt_width,plt_width*aspect_ratio)
    fig18.clear()
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax18 = fig18.add_subplot(gs[0])
    ax18sub = fig18.add_subplot(gs[1])
    #cNorm  = colors.LogNorm(vmin=10**1.5, vmax = 10**4)
    cmx = plt.get_cmap("viridis")
    cNorm  = colors.LogNorm(vmin=1e4, vmax = 1e8)
    cen_plt = ax18.scatter(objs_pd_comb['concentration'][objs_pd_comb['type']=='Central'],13.8 - (objs_pd_comb['Mpeak_snap'][objs_pd_comb['type']=='Central'])*13.8/4096,c = (objs_pd_comb['max_tindex'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax18.scatter(objs_pd_comb['concentration'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],13.8 - (objs_pd_comb['Mpeak_snap'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')])*13.8/4096,c = (objs_pd_comb['max_tindex'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]*1.5, marker = "*", linewidths = edgewidth)    
    sat_plt = ax18.scatter(objs_pd_comb['concentration'][objs_pd_comb['type']=='Satellite'],13.8 - (objs_pd_comb['Mpeak_snap'][objs_pd_comb['type']=='Satellite'])*13.8/4096,c = (objs_pd_comb['max_tindex'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, marker = "*", linewidths = edgewidth)
    """
    tau90_plt = ax18.scatter(objs_pd_comb['tau90_vir'][objs_pd_comb['type']=='Central'],objs_pd_comb['concentration'][objs_pd_comb['type']=='Central'],c = (objs_pd_comb['max_tindex'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth, marker = '+')
    ax18.scatter(objs_pd_comb['tau90_vir'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb['concentration'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = (objs_pd_comb['max_tindex'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], linewidths = edgewidth, marker = '+')    
    ax18.scatter(objs_pd_comb['tau90_vir'][objs_pd_comb['type']=='Satellite'],objs_pd_comb['concentration'][objs_pd_comb['type']=='Satellite'],c = (objs_pd_comb['max_tindex'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k',s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite'], linewidths = edgewidth, marker = '+')    
    """
    ax18.set_ylabel(r'Lookback Time of Peak M$_{\mathrm{vir}}$ [Gyr]')
    ax18.set_xlabel('Concentration')
    #ax18.set_yscale('log')
    #ax18.set_xscale('log')
    #ax18.axis([-0.5, 14, -0.5, 30])
    legend1 = ax18.legend([cen_plt,sat_plt],['Central','Satellite'],scatterpoints = 1,facecolor = 'white',loc = 4,framealpha = 0,frameon = False)
    ax18.add_artist(legend1)
    cb = mpl.colorbar.ColorbarBase(ax18sub, cmap=cmx, norm=cNorm)
    cb.set_label(r"Tidal Index")
    fig18.tight_layout()
    fig18.show()    
    fig18.savefig(outfile_base + '_concent_tmax.png',dpi = dpi)


    plt.close('all') #339924903
    fig18 = plt.figure(2,figsize=(plt_width,plt_width*aspect_ratio))
    fig18.set_size_inches(plt_width,plt_width*aspect_ratio)
    fig18.clear()
    gs = gridspec.GridSpec(2,1,height_ratios=[1,8],hspace = 0.05)

    ax18_hist = fig18.add_subplot(gs[0], sharey=ax18)
    #xymax = np.max(np.abs(x_star,x_vir))                                                                                                           
    #lim = (int(xym/binwidth) + 1) * binwidth                                                                                                       
    bins = np.arange(1, 3.6 + 0.2, 0.2)
    values, mid_bin = np.histogram(objs_pd_comb['min_dist'][objs_pd_comb['Mpeak_snap']==4096], bins = 10**bins)
    ax18_hist.hist(np.log10((objs_pd_comb['min_dist'][objs_pd_comb['type']=='Central']).tolist()), bins = bins, facecolor = 'none', edgecolor = 'k', linestyle='--')
    ax18_hist.hist(np.log10((objs_pd_comb['min_dist'][objs_pd_comb['Mpeak_snap']==4096]).tolist()), bins = bins, facecolor = 'none', edgecolor = 'k')
    ax18_hist.axis([1, 3.6, 0, np.max(np.histogram(np.log10((objs_pd_comb['min_dist'][objs_pd_comb['type']=='Central']).tolist()), bins = bins)[0])])
    ax18_hist.axis('off')
    
    ax18 = fig18.add_subplot(gs[1])
    #ax18sub = fig18.add_subplot(gs[1])
    cmx = plt.get_cmap("winter_r")
    cNorm  = colors.LogNorm(vmin=1e8, vmax = 1e11)
    cen_plt = ax18.scatter(np.log10(objs_pd_comb['min_dist'][objs_pd_comb['type']=='Central'].tolist()),13.8 - (objs_pd_comb['Mpeak_snap'][objs_pd_comb['type']=='Central'])*13.8/4096,c = (objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax18.scatter(np.log10(objs_pd_comb['min_dist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')].tolist()),13.8 - (objs_pd_comb['Mpeak_snap'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')])*13.8/4096,c = (objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]*1.5, marker = "*", linewidths = edgewidth)    
    sat_plt = ax18.scatter(np.log10(objs_pd_comb['min_dist'][objs_pd_comb['type']=='Satellite'].tolist()),13.8 - (objs_pd_comb['Mpeak_snap'][objs_pd_comb['type']=='Satellite'])*13.8/4096,c = (objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, marker = "*", linewidths = edgewidth)
    xaxis, data = bin_plt(np.log10(np.array(objs_pd_comb.loc[objs_pd_comb['type']=='Central', 'min_dist'].tolist())), 13.8 - np.array(objs_pd_comb.loc[objs_pd_comb['type']=='Central', 'Mpeak_snap'].tolist())*13.8/4096, xmax = 3.6, xmin = 1, nbins = 8)
    ax18.fill_between(xaxis[data[:,4] >= 0], data[data[:,4] >= 0,0], data[data[:,4] >= 0,4], color='k', alpha=0.1, linewidth=0)
    ax18.fill_between(xaxis[data[:,3] >= 0], data[data[:,3] >= 0,1], data[data[:,3] >= 0,3], color='k', alpha=0.1, linewidth=0)
    ax18.plot(xaxis[data[:,2] >= 0], data[data[:,2] >= 0,2], color='k')
    ax18.axis([1,3.6,-0.5, 13])
    ax18.set_ylabel(r'Lookback Time of Peak M$_{\mathrm{vir}}$ [Gyr]')
    ax18.set_xlabel('Min(D$_{\mathrm{massive}}$) [kpc]')
    #ax18.set_yscale('log')
    #ax18.set_xscale('log')
    #ax18.axis([-0.5, 14, -0.5, 30])
    legend1 = ax18.legend([cen_plt,sat_plt],['Central','Satellite'],scatterpoints = 1,facecolor = 'white',loc = 'lower left',frameon = False)
    ax18.add_artist(legend1)
    
    sm = plt.cm.ScalarMappable(cmap=cmx, norm=cNorm)
    cbar = plt.colorbar(sm, location="right",ax=[ax18_hist, ax18],aspect = 40)
    #cb = mpl.colorbar.ColorbarBase(ax18sub, cmap=cmx, norm=cNorm)
    cbar.set_label(r"M$_{\mathrm{vir, peak}}$ [M$_\odot]$")
    #fig18.tight_layout()
    fig18.show()    
    fig18.savefig(outfile_base + '_rClosest_tmax.png',dpi = dpi)

    #### vmax
    cmx = plt.get_cmap("viridis") 
    fig19 = plt.figure(19,figsize=(plt_width,plt_width*aspect_ratio))
    fig19.set_size_inches(plt_width,plt_width*aspect_ratio)
    fig19.clear()
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax19 = fig19.add_subplot(gs[0])
    ax19sub = fig19.add_subplot(gs[1])
    #cNorm  = colors.LogNorm(vmin=10**1.5, vmax = 10**4)
    cNorm  = colors.LogNorm(vmin=10, vmax = 2e3)
    cen_plt = ax19.scatter(objs_pd_comb['vmax_mmax'][objs_pd_comb['type']=='Central'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Central'],c = objs_pd_comb['min_dist'][objs_pd_comb['type']=='Central'], cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    sat_plt = ax19.scatter(objs_pd_comb['vmax_mmax'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb[mstar_key][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = objs_pd_comb['min_dist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]*1.5,marker = '*', linewidths = edgewidth)    
    ax19.scatter(objs_pd_comb['vmax_mmax'][objs_pd_comb['type']=='Satellite'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Satellite'],c = objs_pd_comb['min_dist'][objs_pd_comb['type']=='Satellite'], cmap = cmx, norm = cNorm,edgecolor = 'k',marker = '*', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    ax19.set_xscale('log')
    ax19.set_yscale('log') 
    ax19.set_ylabel(r'M$_*$/M$_\odot$')
    ax19.set_xlabel(r'V$_{\mathrm{max}}$ [km/s]')
    #ax19.axis([2e6, 2e11, 2e2, 5e9])
    legend = ax19.legend([cen_plt,sat_plt],['Central','Satellite'],scatterpoints = 1,facecolor = 'white',loc = 4,framealpha = 0,frameon = False)
    ax19.add_artist(legend)
    cb = mpl.colorbar.ColorbarBase(ax19sub, cmap=cmx, norm=cNorm)
    cb.set_label(r"Min(D$_{\mathrm{massive}}$) [kpc]")
    fig19.tight_layout()
    fig19.show()    
    fig19.savefig(outfile_base + '_SMVmax_rClosestDist.png',dpi = dpi)
       
    fig19 = plt.figure(19,figsize=(plt_width,plt_width*aspect_ratio))
    fig19.set_size_inches(plt_width,plt_width*aspect_ratio)
    fig19.clear()
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax19 = fig19.add_subplot(gs[0])
    ax19sub = fig19.add_subplot(gs[1])
    cmx = plt.get_cmap("viridis")
    cNorm  = colors.Normalize(vmin=0, vmax = 13.8)
    sm = plt.cm.ScalarMappable(cmap=cmx, norm=cNorm) 
    cen_plt = ax19.scatter(objs_pd_comb['vmax_mmax'][objs_pd_comb['type']=='Central'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Central'],c = objs_pd_comb['Mpeak_snap'][objs_pd_comb['type']=='Central']*13.8/4096, cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    sat_plt = ax19.scatter(objs_pd_comb['vmax_mmax'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb[mstar_key][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = objs_pd_comb['Mpeak_snap'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]*13.8/4096, cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]*1.5,marker = '*', linewidths = edgewidth)    
    ax19.scatter(objs_pd_comb['vmax_mmax'][objs_pd_comb['type']=='Satellite'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Satellite'],c = objs_pd_comb['Mpeak_snap'][objs_pd_comb['type']=='Satellite']*13.8/4096, cmap = cmx, norm = cNorm,edgecolor = 'k',marker = '*', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    ax19.set_xscale('log')
    ax19.set_yscale('log') 
    ax19.set_ylabel(r'M$_*$/M$_\odot$')
    ax19.set_xlabel(r'V$_{\mathrm{max}}$ [km/s]')
    #ax19.axis([2e6, 2e11, 2e2, 5e9])
    legend = ax19.legend([cen_plt,sat_plt],['Central','Satellite'],scatterpoints = 1,facecolor = 'white',loc = 4,framealpha = 0,frameon = False)
    ax19.add_artist(legend)
    cb = mpl.colorbar.ColorbarBase(ax19sub, cmap=cmx, norm=cNorm)
    cb.set_label(r'Time of Peak M$_{\mathrm{vir}}$ [Gyr]') #, rotation=0)
    fig19.tight_layout()
    fig19.show()    
    fig19.savefig(outfile_base + '_SMVmax_tmax.png',dpi = dpi)    

    cmx = plt.get_cmap("viridis") 
    fig19 = plt.figure(19,figsize=(plt_width,plt_width*aspect_ratio))
    fig19.set_size_inches(plt_width,plt_width*aspect_ratio)
    fig19.clear()
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax19 = fig19.add_subplot(gs[0])
    ax19sub = fig19.add_subplot(gs[1])
    cNorm  = colors.LogNorm(vmin=10**1.5, vmax = 10**4)
    #cNorm  = colors.LogNorm(vmin=10, vmax = 2e3)
    ax19.scatter(objs_pd_comb['vmax_max'][objs_pd_comb['type']=='Central'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Central'],c = objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central'], cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax19.scatter(objs_pd_comb['vmax_max'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb[mstar_key][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = objs_pd_comb['massiveDist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]*1.5,marker = '*', linewidths = edgewidth)    
    ax19.scatter(objs_pd_comb['vmax_max'][objs_pd_comb['type']=='Satellite'],objs_pd_comb[mstar_key][objs_pd_comb['type']=='Satellite'],c = objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite'], cmap = cmx, norm = cNorm,edgecolor = 'k',marker = '*', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    ax19.set_xscale('log')
    ax19.set_yscale('log') 
    ax19.set_ylabel(r'M$_*$/M$_\odot$')
    ax19.set_xlabel(r'V$_{\mathrm{max}}$ [km/s]')
    #ax19.axis([2e6, 2e11, 2e2, 5e9])
    legend = ax19.legend([cen_plt,sat_plt],['Central','Satellite'],scatterpoints = 1,facecolor = 'white',loc = 4,framealpha = 0,frameon = False)
    ax19.add_artist(legend)
    cb = mpl.colorbar.ColorbarBase(ax19sub, cmap=cmx, norm=cNorm)
    cb.set_label(r"D$_{\mathrm{massive}}$ [kpc]")
    fig19.tight_layout()
    fig19.show()    
    fig19.savefig(outfile_base + '_SMVmax_rMassGal.png',dpi = dpi)

    cmx = plt.get_cmap("viridis") 
    fig19 = plt.figure(19,figsize=(plt_width,plt_width*aspect_ratio))
    fig19.set_size_inches(plt_width,plt_width*aspect_ratio)
    fig19.clear()
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax19 = fig19.add_subplot(gs[0])
    ax19sub = fig19.add_subplot(gs[1])
    cNorm  = colors.LogNorm(vmin=10**1.5, vmax = 10**4)
    #cNorm  = colors.LogNorm(vmin=10, vmax = 2e3)
    ax19.scatter(objs_pd_comb['vmax_max'][objs_pd_comb['type']=='Central'],(objs_pd_comb['Mstar_z0'][objs_pd_comb['type']=='Central'] + objs_pd_comb['mHI'][objs_pd_comb['type']=='Central']),c = objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central'], cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax19.scatter(objs_pd_comb['vmax_max'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],(objs_pd_comb['mstar'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')] + objs_pd_comb['mHI'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]),c = objs_pd_comb['massiveDist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]*1.5,marker = '*', linewidths = edgewidth)    
    ax19.scatter(objs_pd_comb['vmax_max'][objs_pd_comb['type']=='Satellite'],(objs_pd_comb['Mstar_z0'][objs_pd_comb['type']=='Satellite'] + objs_pd_comb['mHI'][objs_pd_comb['type']=='Satellite']),c = objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite'], cmap = cmx, norm = cNorm,edgecolor = 'k',marker = '*', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    ax19.set_xscale('log')
    ax19.set_yscale('log') 
    ax19.set_ylabel(r'(M$_*$ + M$_{HI}$/M$_\odot$')
    ax19.set_xlabel(r'V$_{\mathrm{max}}$ [km/s]')
    #ax19.axis([2e6, 2e11, 2e2, 5e9])
    legend = ax19.legend([cen_plt,sat_plt],['Central','Satellite'],scatterpoints = 1,facecolor = 'white',loc = 4,framealpha = 0,frameon = False)
    ax19.add_artist(legend)
    cb = mpl.colorbar.ColorbarBase(ax19sub, cmap=cmx, norm=cNorm)
    cb.set_label(r"D$_{\mathrm{massive}}$ [kpc]")
    fig19.tight_layout()
    fig19.show()    
    fig19.savefig(outfile_base + '_BMVmax_rMassGal.png',dpi = dpi)  
    
    cmx = plt.get_cmap("viridis") 
    fig20 = plt.figure(20,figsize=(plt_width,plt_width*aspect_ratio))
    fig20.set_size_inches(plt_width,plt_width*aspect_ratio)
    fig20.clear()
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax20 = fig20.add_subplot(gs[0])
    ax20sub = fig20.add_subplot(gs[1])    
    cmx = plt.get_cmap("viridis")
    cNorm  = colors.LogNorm(vmin=10, vmax = 2e3)
    sm = plt.cm.ScalarMappable(cmap=cmx, norm=cNorm) 
    cen_plt = ax20.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],objs_pd_comb['bmass_z0'][objs_pd_comb['type']=='Central']/objs_pd_comb['bmass_peak'][objs_pd_comb['type']=='Central'],c = objs_pd_comb['min_dist'][objs_pd_comb['type']=='Central'], cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    sat_plt = ax20.scatter(objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb['bmass_z0'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]/objs_pd_comb['bmass_peak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = objs_pd_comb['min_dist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]*1.5,marker = '*', linewidths = edgewidth)    
    ax20.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],objs_pd_comb['bmass_z0'][objs_pd_comb['type']=='Satellite']/objs_pd_comb['bmass_peak'][objs_pd_comb['type']=='Satellite'],c = objs_pd_comb['min_dist'][objs_pd_comb['type']=='Satellite'], cmap = cmx, norm = cNorm,edgecolor = 'k',marker = '*', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    """
    cen_plt = ax20.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],objs_pd_comb['mstar'][objs_pd_comb['type']=='Central']/objs_pd_comb['Mstar_Mpeak'][objs_pd_comb['type']=='Central'],c = objs_pd_comb['min_dist'][objs_pd_comb['type']=='Central'], cmap = cmx, norm = cNorm,facecolor = 'none', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    sat_plt = ax20.scatter(objs_pd_comb['Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],objs_pd_comb['mstar'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]/objs_pd_comb['Mstar_Mpeak'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')],c = objs_pd_comb['min_dist'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')], cmap = cmx, norm = cNorm,facecolor = 'none', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') | (objs_pd_comb['type']=='Splashback')]*1.5,marker = '*', linewidths = edgewidth)    
    ax20.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],objs_pd_comb['mstar'][objs_pd_comb['type']=='Satellite']/objs_pd_comb['Mstar_Mpeak'][objs_pd_comb['type']=='Satellite'],c = objs_pd_comb['min_dist'][objs_pd_comb['type']=='Satellite'], cmap = cmx, norm = cNorm,facecolor = 'none',marker = '*', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    """    
    ax20.set_xscale('log')
    ax20.set_yscale('log') 
    ax20.set_ylabel(r'M$_{bary}$/M$_{bary, peak}$')
    ax20.set_xlabel(r'M$_{peak}$/M$_\odot$')
    #ax20.axis([2e6, 2e11, 2e2, 5e9])
    legend = ax20.legend([cen_plt,sat_plt],['Central','Satellite'],scatterpoints = 1,facecolor = 'white',loc = 4,framealpha = 0,frameon = False)
    ax20.add_artist(legend)
    cb = mpl.colorbar.ColorbarBase(ax20sub, cmap=cmx, norm=cNorm)
    cb.set_label(r"Min(D$_{\mathrm{massive}}$) [kpc]") #, rotation=0)
    fig20.tight_layout()
    fig20.show()
    fig20.savefig(outfile_base + '_z0vsPeakBary_rClosestDist.png',dpi = dpi) 
    
if __name__ == '__main__':
    if (socket.gethostname() == "quirm.math.grinnell.edu"):
        prefix = '/home/christenc/Data/Sims/'
        prefix_outfile = '/home/christenc/Figures/marvel/'
    else:
        prefix = '/home/christensen/Storage2/UW/MolecH/Cosmo/'
        prefix_outfile = '/home/christensen/Plots/marvel/'

    tfile_base_cm = 'cptmarvel.cosmo25cmb.4096g5HbwK1BH'
    tfile_cm = prefix + 'cptmarvel.cosmo25cmb/cptmarvel.cosmo25cmb.4096g5HbwK1BH/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096/ahf_200/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096' #'cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096'

    tfile_r = prefix + 'rogue.cosmo25cmb/rogue.cosmo25cmb.4096g5HbwK1BH/rogue.cosmo25cmb.4096g5HbwK1BH.004096/ahf_200/rogue.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_base_r = 'rogue.cosmo25cmb.4096g5HbwK1BH'

    tfile_e = prefix + 'elektra.cosmo25cmb/elektra.cosmo25cmb.4096g5HbwK1BH/elektra.cosmo25cmb.4096g5HbwK1BH.004096/ahf_200/elektra.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_base_e = 'elektra.cosmo25cmb.4096g5HbwK1BH'

    tfile_s = prefix + 'storm.cosmo25cmb/storm.cosmo25cmb.4096g5HbwK1BH/storm.cosmo25cmb.4096g5HbwK1BH.004096/ahf_200/storm.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_base_s = 'storm.cosmo25cmb.4096g5HbwK1BH'
    
    tfile_base_1 = 'h148.cosmo50PLK.3072g3HbwK1BH'
    tfile_1 = prefix + 'h148.cosmo50PLK.3072g/h148.cosmo50PLK.3072g3HbwK1BH/snapshots_200bkgdens/h148.cosmo50PLK.3072g3HbwK1BH.004096' #
    tfile_1 = prefix + 'h148.cosmo50PLK.3072g/h148.cosmo50PLK.3072g3HbwK1BH/h148.cosmo50PLK.3072g3HbwK1BH.004096/ahf_200/h148.cosmo50PLK.3072g3HbwK1BH.004096'
    
    tfile_base_1hr = 'h148.cosmo50PLK.6144g3HbwK1BH/'
    tfile_1hr = prefix + 'h148.cosmo50PLK.6144g/h148.cosmo50PLK.6144g3HbwK1BH/h148.cosmo50PLK.6144g3HbwK1BH.004096/ahf_200/h148.cosmo50PLK.6144g3HbwK1BH.004096'
    
    tfile_base_2 = 'h229.cosmo50PLK.3072gst5HbwK1BH'
    tfile_2 = prefix + 'h229.cosmo50PLK.3072g/h229.cosmo50PLK.3072gst5HbwK1BH/snapshots_200bkgdens/h229.cosmo50PLK.3072gst5HbwK1BH.004096' #
    tfile_2 = prefix + 'h229.cosmo50PLK.3072g/h229.cosmo50PLK.3072gst5HbwK1BH/h229.cosmo50PLK.3072gst5HbwK1BH.004096/ahf_200/h229.cosmo50PLK.3072gst5HbwK1BH.004096'    
    
    tfile_base_3 = 'h242.cosmo50PLK.3072gst5HbwK1BH'
    tfile_3 = prefix + 'h242.cosmo50PLK.3072g/h242.cosmo50PLK.3072gst5HbwK1BH/snapshots_200bkgdens/h242.cosmo50PLK.3072gst5HbwK1BH.004096' #
    tfile_3 = prefix + 'h242.cosmo50PLK.3072g/h242.cosmo50PLK.3072gst5HbwK1BH/h242.cosmo50PLK.3072gst5HbwK1BH.004096/ahf_200/h242.cosmo50PLK.3072gst5HbwK1BH.004096'

    tfile_base_4 = 'h329.cosmo50PLK.3072gst5HbwK1BH'
    tfile_4 = prefix + 'h329.cosmo50PLK.3072g/h329.cosmo50PLK.3072gst5HbwK1BH/snapshots_200bkgdens/h329.cosmo50PLK.3072gst5HbwK1BH.004096' #
    tfile_4 = prefix + 'h329.cosmo50PLK.3072g/h329.cosmo50PLK.3072gst5HbwK1BH/h329.cosmo50PLK.3072gst5HbwK1BH.004096/ahf_200/h329.cosmo50PLK.3072gst5HbwK1BH.004096'
    
    tfile_base_4hr = 'h329.cosmo50PLK.6144g5HbwK1BH'
    tfile_4hr = prefix + 'h329.cosmo50PLK.6144g/h329.cosmo50PLK.6144g5HbwK1BH/h329.cosmo50PLK.6144g5HbwK1BH.004096/ahf_200/h329.cosmo50PLK.6144g5HbwK1BH.004096' #    
    
    outfile_base = prefix_outfile + 'marvelJL'
    #SMHMv__distance([tfile_1,tfile_2,tfile_3,tfile_4],outfile_base,[tfile_base_1,tfile_base_2,tfile_base_3,tfile_base_4])
    tfiles = [tfile_1, tfile_2, tfile_3, tfile_4]
    tfile_base = [tfile_base_1, tfile_base_2, tfile_base_3, tfile_base_4]

    tfiles = [tfile_cm, tfile_e, tfile_r, tfile_s, tfile_1, tfile_2, tfile_3, tfile_4]
    tfile_base = [tfile_base_cm, tfile_base_e, tfile_base_r, tfile_base_s, tfile_base_1, tfile_base_2, tfile_base_3, tfile_base_4]
    
    tfiles = [tfile_cm, tfile_e, tfile_r, tfile_s, tfile_1hr, tfile_2, tfile_3, tfile_4hr]
    tfile_base = [tfile_base_cm, tfile_base_e, tfile_base_r, tfile_base_s, tfile_base_1hr, tfile_base_2, tfile_base_3, tfile_base_4hr]
    
    objs_pd_comb = SMHM_v_distance_data(tfiles, outfile_base, tfile_base)
    print(objs_pd_comb[(objs_pd_comb['Mstar_z0_photo'] > 5e5) & (objs_pd_comb['Mstar_z0_photo'] < 5e7) & (objs_pd_comb['massiveDist'] > 800) ]['haloid']
    #SMHM_v_distance_plts(tfiles,outfile_base,tfile_base,objs_pd_comb)

    
