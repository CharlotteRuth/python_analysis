#Charlotte Christensen

#8/13/19
#Plot the SMHM relation for the marvel and Justice League runs, coloring points according to distance/tau_90

#SMHM for environment paper

#%run /home/christenc/Code/python/python_analysis/SMHM_v_distance
import matplotlib as mpl
mpl.use('tkagg') #Also can try mpl.use('Agg') #for using mpl over ssh    
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


# Read in a pickled data file for python 3
def pickle_read(file):

    objs = []
    f=open(file, 'rb')
    while 1:
        try:
            u = pickle._Unpickler(f)
            u.encoding = 'latin1'
            p = u.load()
            objs.append(p)
            #objs.append(pickle.load(f))
        except EOFError:
            break
        
    f.close()

    return pd.DataFrame(objs)


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
            possible_match = (float(objs_pd_sim[objs_pd_sim['haloid'] == halo]['M_star']) < fdmdata[fdmdata['simname']==sim]['Mstar_z0']*(1 + smass_tol)) & (float(objs_pd_sim[objs_pd_sim['haloid'] == halo]['M_star']) > fdmdata[fdmdata['simname']==sim]['Mstar_z0']*(1 - smass_tol))& (float(objs_pd_sim[objs_pd_sim['haloid'] == halo]['mass']) < fdmdata[fdmdata['simname']==sim]['Mhalo_z0']*(1 + vmass_tol)) & (float(objs_pd_sim[objs_pd_sim['haloid'] == halo]['mass']) > fdmdata[fdmdata['simname']==sim]['Mhalo_z0']*(1 - vmass_tol))
            if sum(possible_match) == 0:
                print(sim, halo, 'XXX', float(objs_pd[(objs_pd['sim'] == sim) & (objs_pd['haloid'] == halo)]['M_star']), float(objs_pd[(objs_pd['sim'] == sim) & (objs_pd['haloid'] == halo)]['mass']), 'No Match')
            else:
                #print(fdmdata[fdmdata['simname']==sim][possible_match]['halogrp_z0'])
                arg_best_match = np.argmin(np.abs(fdmdata[fdmdata['simname']==sim][possible_match]['Mstar_z0'] - float(objs_pd_sim[objs_pd_sim['haloid'] == halo]['M_star'])))
                index_best_match = (fdmdata[fdmdata['simname']==sim][possible_match]).index[arg_best_match]
                objs_pd.loc[(objs_pd['sim'] == sim) & (objs_pd['haloid'] == halo),'m200_haloid'] = fdmdata.loc[index_best_match]['halogrp_z0']
                print(sim, halo, fdmdata.loc[index_best_match]['halogrp_z0'], float(objs_pd[(objs_pd['sim'] == sim) & (objs_pd['haloid'] == halo)]['M_star']), fdmdata.loc[index_best_match]['Mstar_z0'], float(objs_pd[(objs_pd['sim'] == sim) & (objs_pd['haloid'] == halo)]['mass']), fdmdata.loc[index_best_match]['Mhalo_z0'])
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
        if s=='h148' or s=='h229' or s=='h242' or s=='h329': # or s=='h148_6144' or s=='h329_6144': # if sat simulation, find distance to halo 1
            if s != sprev:
                tfile_it = tfile_it + 1
                sprev = s
            h1dist = data['h1dist'].tolist()[i]*0.6776942783267969
            massiveDist.append(h1dist)
       
            h1rvir = data['Rvir'][(data.sim==s) & (data.haloid==1)].tolist()[0]*0.6776942783267969
            #hostrvirs.append(h1rvir)
           
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

def SMHM_v_distance(tfiles,outfile_base,tfile_base,*halo_nums):
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

#From McConnachie data
    MW_M31_dist = [58,184,75,40,110,218,110,104,133,180,162,174,45,1350,1066,64,40,218,161,107,681,2030,1862,149,252,520,2266,2436,1945,2187,803,155,422,179,50,23,142,187,2078,1301,2860,1930,452,474,415,86,28,1435,1430,61,882,1367,2583,2288,2387,43,836]
    MW_M31_MHI_LV = [0,0,0,0,0,0,0,0,0,0,0,0,0,0.5588855,3.6373396,0,0,0,0,0,0,0.1292644,0.2144484,0.0085202,0.595621,1.3706083,0.3962263,0.9891736,0.624481,4.7199784,0.205775,0,2.1045441,0,0.3095693,0,0,0.0016497,0.8856712,3.2901259,0.2195973,0.583369,1.2762723,0.8449909,0.1581908,0.0978189,0,2.6453962,1.7043947,1.0156822,0,0.1163385,1.8810399,2.3728621,55.782561,0,2.2978932,]
    Calc_V_Mag = [-11.66,-12.37,-9.97,-8.12,-9.14,-12.61,-7.63,-6.9,-6.4,-6.7,-8.43,-9.4,-8.7,-10.48,-10.32,-3.71,-5.5,-8.59,-4.92,-9.11,-11.19,-11.24,-11.52,-13.44,-15,-14.38,-13.98,-15.84,-15.55,-9.5,-14.51,-5.84,-8,-5.25,-18.12,-16.45,-14.65,-14.75,-18.46,-15.53,-14.08,-18.56,-15.21,-12.3,-9.89,-11.07,-1.5,-13.85,-13.88,-16.83,-9.54,-12.47,-13.16,-12.39,-13.16,-2.7,-13.75]

    MW_M31_dist = np.array(MW_M31_dist)
    MW_M31_MHI_LV = np.array(MW_M31_MHI_LV)
    Calc_V_Mag = np.array(Calc_V_Mag)
    
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

    #Read data from Justin Read's abundance matching, provided by Ferah
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
    Jethwa_upper_inner_Mhalo, Jethwa_upper_inner_Mstar = np.transpose(np.genfromtxt(dataprefix+"webplotdigitizer/Jethwa inner upper.csv", dtype="float", usecols=(0,1), skip_header=0))
    Jethwa_lower_inner_Mhalo, Jethwa_lower_inner_Mstar = np.transpose(np.genfromtxt(dataprefix+"webplotdigitizer/Jethwa inner lower.csv", dtype="float", usecols=(0,1), skip_header=0))

    Jethwa_upper_outer_Mhalo, Jethwa_upper_outer_Mstar = np.transpose(np.genfromtxt(dataprefix+"webplotdigitizer/Jethwa outer upper.csv", dtype="float", usecols=(0,1), skip_header=0))
    Jethwa_lower_outer_Mhalo, Jethwa_lower_outer_Mstar = np.transpose(np.genfromtxt(dataprefix+"webplotdigitizer/Jethwa outer lower.csv", dtype="float", usecols=(0,1), skip_header=0))

    Read_upper_dashed_Mhalo, Read_upper_dashed_Mstar = np.transpose(np.genfromtxt(dataprefix+"webplotdigitizer/Read upper dashed.csv", dtype="float", usecols=(0,1), skip_header=0))
    Read_lower_dashed_Mhalo, Read_lower_dashed_Mstar = np.transpose(np.genfromtxt(dataprefix+"webplotdigitizer/Read lower dashed.csv", dtype="float", usecols=(0,1), skip_header=0))
    Read_upper_solid_Mhalo, Read_upper_solid_Mstar = np.transpose(np.genfromtxt(dataprefix+"webplotdigitizer/Read upper solid.csv", dtype="float", usecols=(0,1), skip_header=0))
    Read_lower_solid_Mhalo, Read_lower_solid_Mstar = np.transpose(np.genfromtxt(dataprefix+"webplotdigitizer/Read lower solid.csv", dtype="float", usecols=(0,1), skip_header=0))
    
    #f = open(dataprefix+'mstar_vs_mhalo_4Charlotte.txt', 'r')
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

    #read dataprefix+'/assembly_histories.npy'


    #read dataprefix+'/reduced_time_series_data.npy'
    #time_series = []
    #f=open(dataprefix+'/reduced_time_series_data.npy', 'rb')
    #while 1:
    #    try:
    #        time_series.append(pickle.load(f))
    #    except EOFError:
    #        break        
    #f.close()   

    
    objs_pd = None 
    for tfile, base in zip(tfiles, tfile_base):
        objs_dat = []
        print(tfile)
        '''
        f=open(tfile + '.MAP.data', 'rb')
        while 1:
            try:
                objs_dat.append(pickle.load(f))
            except EOFError:
                break        
        f.close()
        '''
        objs_dat = pd.read_csv(tfile + '.MAP.data.csv')
        if len(objs_dat) == 1:
            temp = pd.DataFrame(objs_dat[0])
        else:
            temp = pd.DataFrame(objs_dat)
        simname = base.split('.')[0]
        if (base.split('.')[2])[0] == '6':
            simname = simname+'_6144'
            temp['M_star'] = temp['mstar']
            temp['mass'] = temp['mvir']
            
        temp['sim'] = [simname]*len(temp)
        if not 'massiveDist' in temp:
            temp = distance_to_nearest_host(temp,[tfile])
            #temp.to_pickle(tfile + '.MAP.data')
            temp.to_csv(tfile + '.MAP.data.csv')

        #temp.to_csv(tfile + '.MAP.data.csv', index=False)

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
 
    ind = 0
    tau90 = np.empty(len(objs_pd))            
    for index, row in objs_pd.iterrows():
        row['sfh'] = row['sfh'].replace('  ',' ')
        row['sfh'] = row['sfh'].replace('   ',' ')
        row['sfh'] = row['sfh'].replace('    ',' ')
        row['sfh'] = row['sfh'].replace('     ',' ')
        row['sfh'] = row['sfh'].replace('      ',' ')
        row['sfh'] = row['sfh'].replace('       ',' ')
        row['sfh'] = row['sfh'].replace('        ',' ')
        row['sfh'] = row['sfh'].replace('         ',' ')
        row['sfh'] = row['sfh'].replace('          ',' ')
        row['sfh'] = row['sfh'].replace('           ',' ')
        row['sfh'] = row['sfh'].replace('            ',' ')        
        row['sfh'] = row['sfh'].replace('             ',' ')
        sfh_str = ((row['sfh'])[2:-1].replace('\n','')).split(' ')
        sfh = np.array([float(s) for s in sfh_str if s != ''])
        #sfh = np.array([float(x) for x in sfh_str])
        row['sfhbins'] = row['sfhbins'].replace('  ',' ')
        row['sfhbins'] = row['sfhbins'].replace('   ',' ')
        row['sfhbins'] = row['sfhbins'].replace('    ',' ')
        row['sfhbins'] = row['sfhbins'].replace('     ',' ')
        sfhbins_str = ((row['sfhbins'])[2:-1].replace('\n','')).split(' ')
        sfhbins = np.array([float(s) for s in sfhbins_str if s != ''])
        #sfhbins = np.array([float(x) for x in sfhbins_str])
        #sfh = row['sfh']
        #sfhbins = row['sfhbins']
        
        if len(sfhbins) != len(sfh):
            xarr = sfhbins[1:] - (sfhbins[1] - sfhbins[0])
        else:
            xarr = sfhbins[:]
        print(min(xarr),max(xarr),len(xarr))
        yarr = np.cumsum(sfh)/max(np.cumsum(sfh))
        if (yarr[0] >= 0.9):
            tau90[ind] = xarr[0]
        else:
            interp = interp1d(yarr, xarr) #, kind='cubic')
            if np.isnan(interp(0.9)):
                tau90[ind] = 0
            else:
                tau90[ind] = float(interp(0.9))
        ind = ind + 1            

    objs_pd['tau90'] = tau90
    
    halo_label = {'halo_label': ""}
    objs_pd = objs_pd.join(pd.DataFrame(columns=halo_label))
    objs_pd.loc[~objs_pd['m200_haloid'].isnull(),'halo_label'] = objs_pd[~objs_pd['m200_haloid'].isnull()]['sim']+objs_pd[~objs_pd['m200_haloid'].isnull()]['m200_haloid'].astype(str)
    objs_pd = objs_pd.set_index('halo_label')

    fdmdata = fdmdata.join(pd.DataFrame(columns=halo_label))
    fdmdata['halo_label'] = fdmdata['simname']+fdmdata['halogrp_z0'].astype(str)
    fdmdata = fdmdata.set_index('halo_label') 

    objs_pd_comb = pd.concat([objs_pd,fdmdata], join="inner", axis=1)

    #Set pointsizes so that lower res simulations are smaller
    objs_pd_comb['p_size'] = np.ones(len(objs_pd))*markersize
    mask = objs_pd_comb[(objs_pd_comb['simname']=='h148') | (objs_pd_comb['simname']=='h229') | (objs_pd_comb['simname']=='h242') | (objs_pd_comb['simname']=='h329')].index
    objs_pd_comb.loc[mask,'p_size'] = markersize*0.5

    #Calculate the metals in the gas
    objs_pd_comb['mZgas'] = objs_pd_comb['mZISM'] + objs_pd_comb['mZCool'] + objs_pd_comb['mZwarm'] + objs_pd_comb['mZHot']
    
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
    read, = ax1.plot( [], [], 'g', color="grey", linestyle="solid", linewidth=4) #, label="Read+ 2017" )
    nadler, = ax1.plot( [], [], 'g', color="#33669A", linestyle="solid", linewidth=8) #, label="Nadler+ 2020" )
    jethwa, = ax1.plot( [], [], 'g', color="#66CC66", linestyle="solid", linewidth=8) #, label="Jethwa+ 2018" )    
    #ax1.plot( Read_upper_dashed_Mhalo, Read_upper_dashed_Mstar, 'g', linestyle='--', color="#FF9966", linewidth=2 )
    #ax1.plot( Read_lower_dashed_Mhalo, Read_lower_dashed_Mstar, 'g', linestyle='--', color="#FF9966", linewidth=2 )
    #ax1.plot( Read_upper_solid_Mhalo, Read_upper_solid_Mstar, 'g', linestyle='solid', color="#FF9966", linewidth=4 )
    #ax1.plot( Read_lower_solid_Mhalo, Read_lower_solid_Mstar, 'g', linestyle='solid', color="#FF9966", linewidth=4 )
    cen_v_tau90 = ax1.scatter(objs_pd_comb['mass'][objs_pd_comb['type']=='Central'],objs_pd_comb['Mstar_Mpeak'][objs_pd_comb['type']=='Central'],c = objs_pd_comb['tau90'][objs_pd_comb['type']=='Central'], cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    cen_v_tau90 = ax1.scatter(objs_pd_comb['mass'][objs_pd_comb['type']=='Backsplash'],objs_pd_comb['Mstar_Mpeak'][objs_pd_comb['type']=='Backsplash'],c = objs_pd_comb['tau90'][objs_pd_comb['type']=='Backsplash'], cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Backsplash'], linewidths = edgewidth)
    sat_v_tau90 = ax1.scatter(objs_pd_comb['mass'][objs_pd_comb['type']=='Satellite'],objs_pd_comb['Mstar_Mpeak'][objs_pd_comb['type']=='Satellite'],c = objs_pd_comb['tau90'][objs_pd_comb['type']=='Satellite'], cmap = cmx, norm = cNorm, edgecolor = 'k',marker = "*", s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    #ax1.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_ylabel(r'M$_*$/M$_\odot$')
    ax1.set_xlabel(r'M$_{\mathrm{vir}}$/M$_\odot$')
    ax1.axis([2e6, 2e11, 2e2, 5e9])
    legend1 = ax1.legend([cen_v_tau90,sat_v_tau90],['Central','Satellite'],scatterpoints = 1,facecolor = 'white',loc = 4,framealpha = 0,frameon = False)
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
    read, = ax1.plot( [], [], 'g', color="grey", linestyle="solid", linewidth=4) #, label="Read+ 2017" )
    nadler, = ax1.plot( [], [], 'g', color="#33669A", linestyle="solid", linewidth=8) #, label="Nadler+ 2020" )
    jethwa, = ax1.plot( [], [], 'g', color="#66CC66", linestyle="solid", linewidth=8) #, label="Jethwa+ 2018" ) 
    cen_v_tau90 = ax1.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],objs_pd_comb['Mstar_z0'][objs_pd_comb['type']=='Central'],c = objs_pd_comb['tau90'][objs_pd_comb['type']=='Central'], cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    sat_v_tau90 = ax1.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],objs_pd_comb['Mstar_z0'][objs_pd_comb['type']=='Satellite'],c = objs_pd_comb['tau90'][objs_pd_comb['type']=='Satellite'], cmap = cmx, norm = cNorm,edgecolor = 'k',marker="*", s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    #ax1.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_ylabel(r'M$_{*, z = 0}$/M$_\odot$')
    ax1.set_xlabel(r'M$_{\mathrm{vir, peak}}$/M$_\odot$')
    ax1.axis([1e8, 2e11, 2e2, 5e9])
    legend1 = ax1.legend([cen_v_tau90,sat_v_tau90],['Central','Satellite'],scatterpoints = 1,facecolor = 'white',loc = 0,framealpha = 0,frameon = False)
    legend2 = ax1.legend([read,nadler,jethwa],['Read+ 2017','Nadler+ 2020','Jethwa+ 2018'],facecolor = 'white',loc = 4,framealpha = 0,frameon = False)
    ax1.add_artist(legend1)
    ax1.add_artist(legend2)    
    cb = mpl.colorbar.ColorbarBase(ax1sub, cmap=cmx, norm=cNorm)
    cb.set_label(r"$\tau_{90}$ (Gyr)")
    fig1.tight_layout()
    fig1.show()    
    fig1.savefig(outfile_base + '_SMHM_t90_Mpeak.png',dpi = dpi)
    
    #SMHM colored by distance to massive galaxy
    #plt.clf()
    fig2 = plt.figure(2,figsize=(plt_width,plt_width*aspect_ratio))
    fig2.set_size_inches(plt_width,plt_width*aspect_ratio)
    fig2.clear()
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax2 = fig2.add_subplot(gs[0])
    ax2sub = fig2.add_subplot(gs[1])
    cNorm  = colors.Normalize(vmin=1.5, vmax = 4)
    ax2.fill_between(read_abunmatch['M200'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_low'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_high'][read_abunmatch['M200']> 1e9],color = 'grey',linewidth=0,alpha=0.6)
    ax2.fill_between(upper_fid_Mhalo, np.interp( upper_fid_Mhalo, lower_fid_Mhalo, lower_fid_Mstar ), upper_fid_Mstar, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3)#Nadler outer
    ax2.fill_between(upper_fid_inner_Mhalo, np.interp( upper_fid_inner_Mhalo, lower_fid_inner_Mhalo, lower_fid_inner_Mstar ),upper_fid_inner_Mstar, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3) #Nadler inner
    ax2.fill_between(Jethwa_upper_outer_Mhalo, np.interp( Jethwa_upper_outer_Mhalo, Jethwa_lower_outer_Mhalo, Jethwa_lower_outer_Mstar ),Jethwa_upper_outer_Mstar, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3)#Jethwa outer
    ax2.fill_between(Jethwa_upper_inner_Mhalo,np.interp( Jethwa_upper_inner_Mhalo, Jethwa_lower_inner_Mhalo, Jethwa_lower_inner_Mstar ), Jethwa_upper_inner_Mstar, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3) #Jethwa inner
    read, = ax1.plot( [], [], 'g', color="grey", linestyle="solid", linewidth=4) #, label="Read+ 2017" )
    nadler, = ax1.plot( [], [], 'g', color="#33669A", linestyle="solid", linewidth=8) #, label="Nadler+ 2020" )
    jethwa, = ax1.plot( [], [], 'g', color="#66CC66", linestyle="solid", linewidth=8) #, label="Jethwa+ 2018" )     
    ax2.scatter(objs_pd_comb['mass'][objs_pd_comb['type']=='Central'],objs_pd_comb['Mstar_Mpeak'][objs_pd_comb['type']=='Central'],c = np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax2.scatter(objs_pd_comb['mass'][objs_pd_comb['type']=='Backsplash'],objs_pd_comb['Mstar_Mpeak'][objs_pd_comb['type']=='Backsplash'],c = np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Backsplash']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Backsplash'], linewidths = edgewidth)    
    ax2.scatter(objs_pd_comb['mass'][objs_pd_comb['type']=='Satellite'],objs_pd_comb['Mstar_Mpeak'][objs_pd_comb['type']=='Satellite'],c = np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k',marker = '*', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    #ax2.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_ylabel(r'M$_*$/M$_\odot$')
    ax2.set_xlabel(r'M$_{\mathrm{vir}}$/M$_\odot$')
    ax2.axis([2e6, 2e11, 2e2, 5e9])
    legend1 = ax2.legend([cen_v_tau90,sat_v_tau90],['Central','Satellite'],scatterpoints = 1,facecolor = 'white',loc = 4,framealpha = 0,frameon = False)
    legend2 = ax2.legend([read,nadler,jethwa],['Read+ 2017','Nadler+ 2020','Jethwa+ 2018'],facecolor = 'white',loc = 0,framealpha = 0,frameon = False)
    ax2.add_artist(legend1)
    ax2.add_artist(legend2)    
    cb = mpl.colorbar.ColorbarBase(ax2sub, cmap=cmx, norm=cNorm)
    cb.set_label(r"Log(D$_{\mathrm{massive}}$/1 kpc)")
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
    cNorm  = colors.Normalize(vmin=1.5, vmax = 4)
    ax2.fill_between(read_abunmatch['M200'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_low'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_high'][read_abunmatch['M200']> 1e9],color = 'grey',linewidth=0,alpha=0.6)
    ax2.fill_between(upper_fid_Mhalo, np.interp( upper_fid_Mhalo, lower_fid_Mhalo, lower_fid_Mstar ), upper_fid_Mstar, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3)#Nadler outer
    ax2.fill_between(upper_fid_inner_Mhalo, np.interp( upper_fid_inner_Mhalo, lower_fid_inner_Mhalo, lower_fid_inner_Mstar ),upper_fid_inner_Mstar, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3) #Nadler inner
    ax2.fill_between(Jethwa_upper_outer_Mhalo, np.interp( Jethwa_upper_outer_Mhalo, Jethwa_lower_outer_Mhalo, Jethwa_lower_outer_Mstar ),Jethwa_upper_outer_Mstar, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3)#Jethwa outer
    ax2.fill_between(Jethwa_upper_inner_Mhalo,np.interp( Jethwa_upper_inner_Mhalo, Jethwa_lower_inner_Mhalo, Jethwa_lower_inner_Mstar ), Jethwa_upper_inner_Mstar, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3) #Jethwa inner
    read, = ax1.plot( [], [], 'g', color="grey", linestyle="solid", linewidth=4) #, label="Read+ 2017" )
    nadler, = ax1.plot( [], [], 'g', color="#33669A", linestyle="solid", linewidth=8) #, label="Nadler+ 2020" )
    jethwa, = ax1.plot( [], [], 'g', color="#66CC66", linestyle="solid", linewidth=8) #, label="Jethwa+ 2018" )    
    ax2.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],objs_pd_comb['Mstar_z0'][objs_pd_comb['type']=='Central'],c = np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax2.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Backsplash'],objs_pd_comb['Mstar_z0'][objs_pd_comb['type']=='Backsplash'],c = np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Backsplash']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Backsplash'], linewidths = edgewidth)     
    ax2.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],objs_pd_comb['Mstar_z0'][objs_pd_comb['type']=='Satellite'],c = np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k', marker = "*", s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    #ax2.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_ylabel(r'M$_{*, z = 0}$/M$_\odot$')
    ax2.set_xlabel(r'M$_{\mathrm{vir, peak}}$/M$_\odot$')
    ax2.axis([1e8, 2e11, 2e2, 5e9])
    legend1 = ax2.legend([cen_v_tau90,sat_v_tau90],['Central','Satellite'],scatterpoints = 1,facecolor = 'white',loc = 0,framealpha = 0,frameon = False)
    legend2 = ax2.legend([read,nadler,jethwa],['Read+ 2017','Nadler+ 2020','Jethwa+ 2018'],facecolor = 'white',loc = 4,framealpha = 0,frameon = False)
    ax2.add_artist(legend1)
    ax2.add_artist(legend2)      
    cb = mpl.colorbar.ColorbarBase(ax2sub, cmap=cmx, norm=cNorm)
    cb.set_label(r"Log(D$_{\mathrm{massive}}$/1 kpc)")
    fig2.tight_layout()
    fig2.show() 
    fig2.savefig(outfile_base + '_SMHM_rMassGal_Mpeak.png',dpi = dpi)

    plt.close('all')
    #fig2.clear()
    fig2 = plt.figure(2,figsize=(plt_width,plt_width*aspect_ratio))
    fig2.set_size_inches(plt_width,plt_width*aspect_ratio)
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax2 = fig2.add_subplot(gs[0])
    ax2sub = fig2.add_subplot(gs[1])
    cNorm  = colors.Normalize(vmin=1.5, vmax = 4)
    ax2.fill_between(read_abunmatch['M200'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_low'][read_abunmatch['M200']> 1e9]/read_abunmatch['M200'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_high'][read_abunmatch['M200']> 1e9]/read_abunmatch['M200'][read_abunmatch['M200']> 1e9],color = 'grey',linewidth=0,alpha = 0.6)
    ax2.fill_between(upper_fid_Mhalo, np.interp( upper_fid_Mhalo, lower_fid_Mhalo, lower_fid_Mstar )/upper_fid_Mhalo, upper_fid_Mstar/upper_fid_Mhalo, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3)#Nadler outer
    ax2.fill_between(upper_fid_inner_Mhalo, np.interp( upper_fid_inner_Mhalo, lower_fid_inner_Mhalo, lower_fid_inner_Mstar)/upper_fid_inner_Mhalo,upper_fid_inner_Mstar/upper_fid_inner_Mhalo, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3) #Nadler inner
    ax2.fill_between(Jethwa_upper_outer_Mhalo, np.interp( Jethwa_upper_outer_Mhalo, Jethwa_lower_outer_Mhalo, Jethwa_lower_outer_Mstar )/Jethwa_upper_outer_Mhalo,Jethwa_upper_outer_Mstar/Jethwa_upper_outer_Mhalo, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3)#Jethwa outer
    ax2.fill_between(Jethwa_upper_inner_Mhalo,np.interp( Jethwa_upper_inner_Mhalo, Jethwa_lower_inner_Mhalo, Jethwa_lower_inner_Mstar )/Jethwa_upper_inner_Mhalo, Jethwa_upper_inner_Mstar/Jethwa_upper_inner_Mhalo, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3) #Jethwa inner
    read, = ax2.plot( [], [], 'g', color="grey", linestyle="solid", linewidth=4) #, label="Read+ 2017" )
    nadler, = ax2.plot( [], [], 'g', color="#33669A", linestyle="solid", linewidth=8) #, label="Nadler+ 2020" )
    jethwa, = ax2.plot( [], [], 'g', color="#66CC66", linestyle="solid", linewidth=8) #, label="Jethwa+ 2018" )        
    ax2.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],objs_pd_comb['Mstar_Mpeak'][objs_pd_comb['type']=='Central']/objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],c = np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax2.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Backsplash'],objs_pd_comb['Mstar_Mpeak'][objs_pd_comb['type']=='Backsplash']/objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Backsplash'],c = np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Backsplash']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Backsplash'], linewidths = edgewidth)
    ax2.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],objs_pd_comb['Mstar_Mpeak'][objs_pd_comb['type']=='Satellite']/objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],c = np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k', marker = "*", s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    #ax2.scatter(readdata['M200']*1e10,readdata['M_*']*1e7/(readdata['M200']*1e10),marker = '+',c = 'k')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_ylabel(r'M$_{*, z = 0}$/M$_{\mathrm{vir, peak}}$')
    ax2.set_xlabel(r'M$_{\mathrm{vir, peak}}$/M$_\odot$')
    ax2.axis([1e8, 2e11, 1e-6, 0.05])
    legend1 = ax2.legend([cen_v_tau90,sat_v_tau90],['Central','Satellite'],scatterpoints = 1,facecolor = 'white',loc = 4,framealpha = 0,frameon = False)
    legend2 = ax2.legend([read,nadler,jethwa],['Read+ 2017','Nadler+ 2020','Jethwa+ 2018'],facecolor = 'white',loc = 0,framealpha = 0,frameon = False)
    ax2.add_artist(legend1)
    ax2.add_artist(legend2)      
    cb = mpl.colorbar.ColorbarBase(ax2sub, cmap=cmx, norm=cNorm)
    cb.set_label(r"Log(D$_{\mathrm{massive}}$/1 kpc)")
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
    cNorm  = colors.Normalize(vmin=1.5, vmax = 4)
    ax2a.fill_between(read_abunmatch['M200'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_low'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_high'][read_abunmatch['M200']> 1e9],color = 'grey',linewidth=0,alpha = 0.6)
    ax2a.fill_between(upper_fid_Mhalo, np.interp( upper_fid_Mhalo, lower_fid_Mhalo, lower_fid_Mstar ), upper_fid_Mstar, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3)#Nadler outer
    ax2a.fill_between(upper_fid_inner_Mhalo, np.interp( upper_fid_inner_Mhalo, lower_fid_inner_Mhalo, lower_fid_inner_Mstar ),upper_fid_inner_Mstar, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3) #Nadler inner
    ax2a.fill_between(Jethwa_upper_outer_Mhalo, np.interp( Jethwa_upper_outer_Mhalo, Jethwa_lower_outer_Mhalo, Jethwa_lower_outer_Mstar ),Jethwa_upper_outer_Mstar, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3)#Jethwa outer
    ax2a.fill_between(Jethwa_upper_inner_Mhalo,np.interp( Jethwa_upper_inner_Mhalo, Jethwa_lower_inner_Mhalo, Jethwa_lower_inner_Mstar ), Jethwa_upper_inner_Mstar, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3) #Jethwa inner
    read, = ax1.plot( [], [], 'g', color="grey", linestyle="solid", linewidth=4) #, label="Read+ 2017" )
    nadler, = ax1.plot( [], [], 'g', color="#33669A", linestyle="solid", linewidth=8) #, label="Nadler+ 2020" )
    jethwa, = ax1.plot( [], [], 'g', color="#66CC66", linestyle="solid", linewidth=8) #, label="Jethwa+ 2018" )     
    ax2a.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],objs_pd_comb['Mstar_z0'][objs_pd_comb['type']=='Central'],c = np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax2a.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Backsplash'],objs_pd_comb['Mstar_z0'][objs_pd_comb['type']=='Backsplash'],c = np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Backsplash']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Backsplash'], linewidths = edgewidth)    
    ax2a.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],objs_pd_comb['Mstar_z0'][objs_pd_comb['type']=='Satellite'],c = np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k', marker = "*", s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    #ax2.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax2a.set_xscale('log')
    ax2a.set_yscale('log')
    ax2a.set_ylabel(r'M$_{*, z = 0}$/M$_\odot$')
    #ax2a.set_xlabel(r'M$_{\mathrm{vir, peak}}$/M$_\odot$')
    ax2a.axis([1e8, 2e11, 2e2, 5e9])
    #legend1 = ax2a.legend([cen_v_tau90,sat_v_tau90],['Central','Satellite'],scatterpoints = 1,facecolor = 'white',loc = 0,framealpha = 0,frameon = False)
    legend2 = ax2a.legend([read,nadler,jethwa],['Read+ 2017','Nadler+ 2020','Jethwa+ 2018'],facecolor = 'white',loc = 4,framealpha = 0,frameon = False)
    #ax2.add_artist(legend1)
    ax2a.add_artist(legend2)
    ax2a.label_outer()

    cNorm  = colors.Normalize(vmin=1.5, vmax = 4)
    ax2b.fill_between(read_abunmatch['M200'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_low'][read_abunmatch['M200']> 1e9]/read_abunmatch['M200'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_high'][read_abunmatch['M200']> 1e9]/read_abunmatch['M200'][read_abunmatch['M200']> 1e9],color = 'grey',linewidth=0,alpha = 0.6)
    ax2b.fill_between(upper_fid_Mhalo, np.interp( upper_fid_Mhalo, lower_fid_Mhalo, lower_fid_Mstar )/upper_fid_Mhalo, upper_fid_Mstar/upper_fid_Mhalo, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3)#Nadler outer
    ax2b.fill_between(upper_fid_inner_Mhalo, np.interp( upper_fid_inner_Mhalo, lower_fid_inner_Mhalo, lower_fid_inner_Mstar)/upper_fid_inner_Mhalo,upper_fid_inner_Mstar/upper_fid_inner_Mhalo, facecolor="#33669A",linewidth=0,zorder=0,alpha=0.3) #Nadler inner
    ax2b.fill_between(Jethwa_upper_outer_Mhalo, np.interp( Jethwa_upper_outer_Mhalo, Jethwa_lower_outer_Mhalo, Jethwa_lower_outer_Mstar )/Jethwa_upper_outer_Mhalo,Jethwa_upper_outer_Mstar/Jethwa_upper_outer_Mhalo, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3)#Jethwa outer
    ax2b.fill_between(Jethwa_upper_inner_Mhalo,np.interp( Jethwa_upper_inner_Mhalo, Jethwa_lower_inner_Mhalo, Jethwa_lower_inner_Mstar )/Jethwa_upper_inner_Mhalo, Jethwa_upper_inner_Mstar/Jethwa_upper_inner_Mhalo, facecolor="#66CC66",linewidth=0,zorder=0,alpha=0.3) #Jethwa inner
    read, = ax1.plot( [], [], 'g', color="grey", linestyle="solid", linewidth=4) #, label="Read+ 2017" )
    nadler, = ax1.plot( [], [], 'g', color="#33669A", linestyle="solid", linewidth=8) #, label="Nadler+ 2020" )
    jethwa, = ax1.plot( [], [], 'g', color="#66CC66", linestyle="solid", linewidth=8) #, label="Jethwa+ 2018" )        
    ax2b.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],objs_pd_comb['Mstar_Mpeak'][objs_pd_comb['type']=='Central']/objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],c = np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax2b.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Backsplash'],objs_pd_comb['Mstar_Mpeak'][objs_pd_comb['type']=='Backsplash']/objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Backsplash'],c = np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Backsplash']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Backsplash'], linewidths = edgewidth)    
    ax2b.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],objs_pd_comb['Mstar_Mpeak'][objs_pd_comb['type']=='Satellite']/objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],c = np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k', marker = "*", s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    #ax2.scatter(readdata['M200']*1e10,readdata['M_*']*1e7/(readdata['M200']*1e10),marker = '+',c = 'k')
    ax2b.set_xscale('log')
    ax2b.set_yscale('log')
    ax2b.set_ylabel(r'M$_{*, z = 0}$/M$_{\mathrm{vir, peak}}$')
    ax2b.set_xlabel(r'M$_{\mathrm{vir, peak}}$/M$_\odot$')
    ax2b.axis([1e8, 2e11, 1e-6, 0.05])
    legend1 = ax2b.legend([cen_v_tau90,sat_v_tau90],['Central','Satellite'],scatterpoints = 1,facecolor = 'white',loc = 4,framealpha = 0,frameon = False)
    #legend2 = ax2b.legend([read,nadler,jethwa],['Read+ 2017','Nadler+ 2020','Jethwa+ 2018'],facecolor = 'white',loc = 0,framealpha = 0,frameon = False)
    ax2b.add_artist(legend1)
    #ax2b.add_artist(legend2)
    ax2b.label_outer()
    fig2.tight_layout()
    
    sm = plt.cm.ScalarMappable(cmap=cmx, norm=cNorm)
    cbar = plt.colorbar(sm, location="right",ax=[ax2a,ax2b],aspect = 40) #pad=0.04)
    #cbar = mpl.colorbar.ColorbarBase(ax2sub, cmap=cmx, norm=cNorm, location="right",ax=[ax2a,ax2b]) #pad=0.04)
    cbar.set_label(r'Log(D$_{\mathrm{massive}}$/1 kpc)') #, rotation=0)
    #fig2.tight_layout()

    #sm = plt.cm.ScalarMappable(cmap=cmx, norm=cNorm)
    #sm.set_array([])
    #cbar = plt.colorbar(sm, location="top", ax=axs2.ravel().tolist(), fraction = 0.1) #pad=0.04)
    #cbar.set_label(r'Log(D$_{\mathrm{massive}}$/1 kpc)') #, rotation=0)
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
    read, = ax3.plot( [], [], 'g', color="grey", linestyle="solid", linewidth=4) #, label="Read+ 2017" )
    nadler, = ax3.plot( [], [], 'g', color="#33669A", linestyle="solid", linewidth=8) #, label="Nadler+ 2020" )
    jethwa, = ax3.plot( [], [], 'g', color="#66CC66", linestyle="solid", linewidth=8) #, label="Jethwa+ 2018" )
    ax3.scatter(objs_pd_comb['mass'][objs_pd_comb['type']=='Central'],objs_pd_comb['Mstar_Mpeak'][objs_pd_comb['type']=='Central'] ,edgecolor = 'k',facecolor = 'none', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)    
    ax3.scatter(objs_pd_comb['mass'][objs_pd_comb['type']=='Central'],objs_pd_comb['Mstar_Mpeak'][objs_pd_comb['type']=='Central'],c = np.log10(np.array(objs_pd_comb['mHI'][objs_pd_comb['type']=='Central'].tolist())), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax3.scatter(objs_pd_comb['mass'][objs_pd_comb['type']=='Backsplash'],objs_pd_comb['Mstar_Mpeak'][objs_pd_comb['type']=='Backsplash'] ,edgecolor = 'k',facecolor = 'none', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Backsplash'], linewidths = edgewidth)    
    ax3.scatter(objs_pd_comb['mass'][objs_pd_comb['type']=='Backsplash'],objs_pd_comb['Mstar_Mpeak'][objs_pd_comb['type']=='Backsplash'],c = np.log10(np.array(objs_pd_comb['mHI'][objs_pd_comb['type']=='Backsplash'].tolist())), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Backsplash'], linewidths = edgewidth)      
    ax3.scatter(objs_pd_comb['mass'][objs_pd_comb['type']=='Satellite'],objs_pd_comb['Mstar_Mpeak'][objs_pd_comb['type']=='Satellite'] ,edgecolor = 'k',facecolor = 'none',marker = '*', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)    
    ax3.scatter(objs_pd_comb['mass'][objs_pd_comb['type']=='Satellite'],objs_pd_comb['Mstar_Mpeak'][objs_pd_comb['type']=='Satellite'],c = np.log10(np.array(objs_pd_comb['mHI'][objs_pd_comb['type']=='Satellite'].tolist())), cmap = cmx, norm = cNorm,edgecolor = 'k',marker = '*', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)  
    #ax3.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax3.set_xscale('log')
    ax3.set_yscale('log')
    ax3.set_ylabel(r'M$_*$/M$_\odot$')
    ax3.set_xlabel(r'M$_{\mathrm{vir}}$/M$_\odot$')
    ax3.axis([2e6, 2e11, 2e2, 5e9])
    legend1 = ax3.legend([cen_v_tau90,sat_v_tau90],['Central','Satellite'],scatterpoints = 1,facecolor = 'white',loc = 4,framealpha = 0,frameon = False)
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
    read, = ax3.plot( [], [], 'g', color="grey", linestyle="solid", linewidth=4) #, label="Read+ 2017" )
    nadler, = ax3.plot( [], [], 'g', color="#33669A", linestyle="solid", linewidth=8) #, label="Nadler+ 2020" )
    jethwa, = ax3.plot( [], [], 'g', color="#66CC66", linestyle="solid", linewidth=8) #, label="Jethwa+ 2018" )
    ax3.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],objs_pd_comb['Mstar_Mpeak'][objs_pd_comb['type']=='Central'] ,edgecolor = 'k',facecolor = 'none', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)    
    ax3.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],objs_pd_comb['Mstar_Mpeak'][objs_pd_comb['type']=='Central'],c = np.log10(np.array(objs_pd_comb['mHI'][objs_pd_comb['type']=='Central'].tolist())), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax3.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Backsplash'],objs_pd_comb['Mstar_Mpeak'][objs_pd_comb['type']=='Backsplash'] ,edgecolor = 'k',facecolor = 'none', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Backsplash'], linewidths = edgewidth)    
    ax3.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Backsplash'],objs_pd_comb['Mstar_Mpeak'][objs_pd_comb['type']=='Backsplash'],c = np.log10(np.array(objs_pd_comb['mHI'][objs_pd_comb['type']=='Backsplash'].tolist())), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Backsplash'], linewidths = edgewidth)
    ax3.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],objs_pd_comb['Mstar_Mpeak'][objs_pd_comb['type']=='Satellite'] ,edgecolor = 'k',facecolor = 'none', marker = '*', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)    
    ax3.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],objs_pd_comb['Mstar_Mpeak'][objs_pd_comb['type']=='Satellite'],c = np.log10(np.array(objs_pd_comb['mHI'][objs_pd_comb['type']=='Satellite'].tolist())), cmap = cmx, norm = cNorm,edgecolor = 'k', marker = '*', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    #ax3.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax3.set_xscale('log')
    ax3.set_yscale('log')
    ax3.set_ylabel(r'M$_*$/M$_\odot$')
    ax3.set_xlabel(r'M$_{\mathrm{vir, peak}}$/M$_\odot$')
    ax3.axis([1e8, 2e11, 2e2, 5e9])
    legend1 = ax3.legend([cen_v_tau90,sat_v_tau90],['Central','Satellite'],scatterpoints = 1,facecolor = 'white',loc = 0,framealpha = 0,frameon = False)
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
    read, = ax4.plot( [], [], 'g', color="grey", linestyle="solid", linewidth=4) #, label="Read+ 2017" )
    nadler, = ax4.plot( [], [], 'g', color="#33669A", linestyle="solid", linewidth=8) #, label="Nadler+ 2020" )
    jethwa, = ax4.plot( [], [], 'g', color="#66CC66", linestyle="solid", linewidth=8) #, label="Jethwa+ 2018" )    
    ax4.scatter(objs_pd_comb['mass'][objs_pd_comb['type']=='Central'],objs_pd_comb['Mstar_Mpeak'][objs_pd_comb['type']=='Central'],c = objs_pd_comb['mHI'][objs_pd_comb['type']=='Central']/(objs_pd_comb['Mstar_Mpeak'][objs_pd_comb['type']=='Central'] + objs_pd_comb['mHI'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax4.scatter(objs_pd_comb['mass'][objs_pd_comb['type']=='Backsplash'],objs_pd_comb['Mstar_Mpeak'][objs_pd_comb['type']=='Backsplash'],c = objs_pd_comb['mHI'][objs_pd_comb['type']=='Backsplash']/(objs_pd_comb['Mstar_Mpeak'][objs_pd_comb['type']=='Backsplash'] + objs_pd_comb['mHI'][objs_pd_comb['type']=='Backsplash']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Backsplash'], linewidths = edgewidth)
    ax4.scatter(objs_pd_comb['mass'][objs_pd_comb['type']=='Satellite'],objs_pd_comb['Mstar_Mpeak'][objs_pd_comb['type']=='Satellite'],c = objs_pd_comb['mHI'][objs_pd_comb['type']=='Satellite']/(objs_pd_comb['Mstar_Mpeak'][objs_pd_comb['type']=='Satellite'] + objs_pd_comb['mHI'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k', marker = '*', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    #ax4.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax4.set_xscale('log')
    ax4.set_yscale('log')
    ax4.set_ylabel(r'M$_*$/M$_\odot$')
    ax4.set_xlabel(r'M$_{\mathrm{vir}}$/M$_\odot$')
    ax4.axis([2e6, 2e11, 2e2, 5e9])
    legend1 = ax4.legend([cen_v_tau90,sat_v_tau90],['Central','Satellite'],scatterpoints = 1,facecolor = 'white',loc = 4,framealpha = 0,frameon = False)
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
    read, = ax4.plot( [], [], 'g', color="grey", linestyle="solid", linewidth=4) #, label="Read+ 2017" )
    nadler, = ax4.plot( [], [], 'g', color="#33669A", linestyle="solid", linewidth=8) #, label="Nadler+ 2020" )
    jethwa, = ax4.plot( [], [], 'g', color="#66CC66", linestyle="solid", linewidth=8) #, label="Jethwa+ 2018" ) 
    ax4.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],objs_pd_comb['Mstar_Mpeak'][objs_pd_comb['type']=='Central'],c = objs_pd_comb['mHI'][objs_pd_comb['type']=='Central']/(objs_pd_comb['Mstar_Mpeak'][objs_pd_comb['type']=='Central'] + objs_pd_comb['mHI'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax4.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Backsplash'],objs_pd_comb['Mstar_Mpeak'][objs_pd_comb['type']=='Backsplash'],c = objs_pd_comb['mHI'][objs_pd_comb['type']=='Backsplash']/(objs_pd_comb['Mstar_Mpeak'][objs_pd_comb['type']=='Backsplash'] + objs_pd_comb['mHI'][objs_pd_comb['type']=='Backsplash']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Backsplash'], linewidths = edgewidth)    
    ax4.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],objs_pd_comb['Mstar_Mpeak'][objs_pd_comb['type']=='Satellite'],c = objs_pd_comb['mHI'][objs_pd_comb['type']=='Satellite']/(objs_pd_comb['Mstar_Mpeak'][objs_pd_comb['type']=='Satellite'] + objs_pd_comb['mHI'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k',marker = '*', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    #ax4.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax4.set_xscale('log')
    ax4.set_yscale('log')
    ax4.set_ylabel(r'M$_*$/M$_\odot$')
    ax4.set_xlabel(r'M$_{\mathrm{vir, peak}}$/M$_\odot$')
    ax4.axis([1e8, 2e11, 2e2, 5e9])
    legend1 = ax4.legend([cen_v_tau90,sat_v_tau90],['Central','Satellite'],scatterpoints = 1,facecolor = 'white',loc = 0,framealpha = 0,frameon = False)
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
    cNorm  = colors.Normalize(vmin=3, vmax = 11.5)
    ax5.scatter(objs_pd_comb['mass'][objs_pd_comb['type']=='Central'],(objs_pd_comb['Mstar_Mpeak'][objs_pd_comb['type']=='Central'] + objs_pd_comb['mHI'][objs_pd_comb['type']=='Central']),c = np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central']*1, linewidths = edgewidth)
    ax5.scatter(objs_pd_comb['mass'][objs_pd_comb['type']=='Backsplash'],(objs_pd_comb['Mstar_Mpeak'][objs_pd_comb['type']=='Backsplash'] + objs_pd_comb['mHI'][objs_pd_comb['type']=='Backsplash']),c = np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Backsplash']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Backsplash']*1, linewidths = edgewidth)
    ax5.scatter(objs_pd_comb['mass'][objs_pd_comb['type']=='Satellite'],(objs_pd_comb['Mstar_Mpeak'][objs_pd_comb['type']=='Satellite'] + objs_pd_comb['mHI'][objs_pd_comb['type']=='Satellite']),c = np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k', marker = '*', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    #ax5.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax5.set_xscale('log')
    ax5.set_yscale('log')
    ax5.set_ylabel(r'M$_*$ + M$_{\mathrm{HI}}$ [M$_\odot$]')
    ax5.set_xlabel(r'M$_{\mathrm{vir}}$ [M$_\odot$]')
    ax5.axis([2e6, 2e11, 1e3, 1e10])
    fig5.tight_layout()
    ax5.legend([cen_v_tau90,sat_v_tau90],['Central','Satellite'],scatterpoints = 1,facecolor = 'white',loc = 0,framealpha = 0,frameon = False)
    sm = plt.cm.ScalarMappable(cmap=cmx, norm=cNorm)
    cbar = plt.colorbar(sm, ax=[ax5],aspect = 20) #pad=0.04)
    cbar.set_label(r'Log(D$_{\mathrm{massive}}$/1 kpc)') #, rotation=0)
    #cb = mpl.colorbar.ColorbarBase(ax5sub, cmap=cmx, norm=cNorm)
    #cb.set_label(r"Log(D$_{\mathrm{massive}}$/1 kpc)")
    fig5.show()
    fig5.set_size_inches(plt_width,plt_width*aspect_ratio)
    fig5.show()
    fig5.savefig(outfile_base + '_BMHM_rMassGal.png',dpi = dpi)

    plt.close('all')
    fig5 = plt.figure(5,figsize=(plt_width,plt_width*aspect_ratio))
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax5 = fig5.add_subplot(gs[0])
    #ax5sub = fig5.add_subplot(gs[1])
    cNorm  = colors.Normalize(vmin=1.5, vmax = 4)    
    ax5.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],(objs_pd_comb['Mstar_Mpeak'][objs_pd_comb['type']=='Central'] + objs_pd_comb['mHI'][objs_pd_comb['type']=='Central']),c = np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax5.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Backsplash'],(objs_pd_comb['Mstar_Mpeak'][objs_pd_comb['type']=='Backsplash'] + objs_pd_comb['mHI'][objs_pd_comb['type']=='Backsplash']),c = np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Backsplash']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Backsplash'], linewidths = edgewidth)    
    ax5.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],(objs_pd_comb['Mstar_Mpeak'][objs_pd_comb['type']=='Satellite'] + objs_pd_comb['mHI'][objs_pd_comb['type']=='Satellite']),c = np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k', marker="*", s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    #ax5.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax5.set_xscale('log')
    ax5.set_yscale('log')
    ax5.set_ylabel(r'M$_*$ + M$_{\mathrm{HI}}$ [M$_\odot$]')
    ax5.set_xlabel(r'M$_{\mathrm{vir, peak}}$ [M$_\odot$]')
    ax5.axis([1e8, 2e11, 1e3, 1e10])
    fig5.tight_layout()
    ax5.legend([cen_v_tau90,sat_v_tau90],['Central','Satellite'],scatterpoints = 1,facecolor = 'white',loc = 0,framealpha = 0,frameon = False)
    sm = plt.cm.ScalarMappable(cmap=cmx, norm=cNorm)
    cbar = plt.colorbar(sm, ax=[ax5],aspect = 20) #pad=0.04)
    cbar.set_label(r'Log(D$_{\mathrm{massive}}$/1 kpc)') #, rotation=0)    
    #cb = mpl.colorbar.ColorbarBase(ax5sub, cmap=cmx, norm=cNorm)
    #cb.set_label(r"Log(D$_{\mathrm{massive}}$/1 kpc)")
    fig5.set_size_inches(plt_width,plt_width*aspect_ratio)
    fig5.show()
    fig5.savefig(outfile_base + '_BMHM_rMassGal_Mpeak.png',dpi = dpi)
    
    # Baryon fraction vs distance ###############################
    plt.close('all')
    fig6 = plt.figure(6,figsize=(plt_width,plt_width*aspect_ratio))
    fig6.clear()
    gs = gridspec.GridSpec(1,1)
    ax6 = fig6.add_subplot(gs[0])
    #Baryonic mass fraction of central galaxies
    ax6.scatter(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central'],(objs_pd_comb['M_star'][objs_pd_comb['type']=='Central'] + objs_pd_comb['mHI'][objs_pd_comb['type']=='Central'])/objs_pd_comb['mass'][objs_pd_comb['type']=='Central'],edgecolor = 'k', facecolor = 'none',          s = (objs_pd_comb['Rvir'][objs_pd_comb['type']=='Central']*2*ms_scale).tolist(), linewidths = edgewidth)
    ax6.scatter(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central'],(objs_pd_comb['M_star'][objs_pd_comb['type']=='Central'] + objs_pd_comb['mHI'][objs_pd_comb['type']=='Central'])/objs_pd_comb['mass'][objs_pd_comb['type']=='Central'],edgecolor = 'k', facecolor = 'k',alpha = 0.3, s = (objs_pd_comb['Rvir'][objs_pd_comb['type']=='Central']*2*ms_scale).tolist(), linewidths = edgewidth)
    ax6.scatter(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Backsplash'],(objs_pd_comb['M_star'][objs_pd_comb['type']=='Backsplash'] + objs_pd_comb['mHI'][objs_pd_comb['type']=='Backsplash'])/objs_pd_comb['mass'][objs_pd_comb['type']=='Backsplash'],edgecolor = 'k', facecolor = 'none',          s = (objs_pd_comb['Rvir'][objs_pd_comb['type']=='Backsplash']*2*ms_scale).tolist(), linewidths = edgewidth)
    ax6.scatter(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Backsplash'],(objs_pd_comb['M_star'][objs_pd_comb['type']=='Backsplash'] + objs_pd_comb['mHI'][objs_pd_comb['type']=='Backsplash'])/objs_pd_comb['mass'][objs_pd_comb['type']=='Backsplash'],edgecolor = 'k', facecolor = 'k',alpha = 0.3, s = (objs_pd_comb['Rvir'][objs_pd_comb['type']=='Backsplash']*2*ms_scale).tolist(), linewidths = edgewidth)
    #Baryonic mass fraction of satellites
    ax6.scatter(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite'],        (objs_pd_comb['M_star'][objs_pd_comb['type']=='Satellite'] + objs_pd_comb['mHI'][objs_pd_comb['type']=='Satellite'])/objs_pd_comb['mass'][objs_pd_comb['type']=='Satellite'],facecolor = 'none',edgecolor = 'k',             marker = '*', s = (objs_pd_comb['Rvir'][objs_pd_comb['type']=='Satellite']*2*ms_scale*1.5).tolist(), linewidths = edgewidth)
    fbary = ax6.scatter(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite'],(objs_pd_comb['M_star'][objs_pd_comb['type']=='Satellite'] + objs_pd_comb['mHI'][objs_pd_comb['type']=='Satellite'])/objs_pd_comb['mass'][objs_pd_comb['type']=='Satellite'],facecolor = 'k',   edgecolor = 'k',alpha = 0.3, marker = '*', s = (objs_pd_comb['Rvir'][objs_pd_comb['type']=='Satellite']*2*ms_scale*1.5).tolist(), linewidths = edgewidth)
    #Stellar mass fraction of central galaxies
    ax6.scatter(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central'],objs_pd_comb['M_star'][objs_pd_comb['type']=='Central']/objs_pd_comb['mass'][objs_pd_comb['type']=='Central'],    edgecolor = 'red',facecolor = 'none',               s = (objs_pd_comb['Rvir'][objs_pd_comb['type']=='Central']*2*ms_scale).tolist(),     linewidths = edgewidth)
    ax6.scatter(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Backsplash'],objs_pd_comb['M_star'][objs_pd_comb['type']=='Backsplash']/objs_pd_comb['mass'][objs_pd_comb['type']=='Backsplash'],    edgecolor = 'red',facecolor = 'none',               s = (objs_pd_comb['Rvir'][objs_pd_comb['type']=='Backsplash']*2*ms_scale).tolist(),     linewidths = edgewidth)    
    #Stellar mass fraction of satellites
    fstar = ax6.scatter(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite'],objs_pd_comb['M_star'][objs_pd_comb['type']=='Satellite']/objs_pd_comb['mass'][objs_pd_comb['type']=='Satellite'],edgecolor = 'red',facecolor = 'none', marker = '*', s = (objs_pd_comb['Rvir'][objs_pd_comb['type']=='Satellite']*2*ms_scale*1.5).tolist(), linewidths = edgewidth)
    #ax6.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    #ax6.scatter(objs_pd['massiveDist'],objs_pd['mHI']/objs_pd['mass'],facecolor = 'none',edgecolor = 'k')
    ax6.set_xscale('log')
    ax6.set_yscale('log')
    ax6.set_ylabel(r'M/M$_{\mathrm{vir}}$')
    ax6.set_xlabel(r"Distance to massive galaxy (kpc)")
    ax6.axis([17, 7e3, 5e-6, 0.2])
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
    #Baryonic mass fraction of central galaxies
    ax6.scatter(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central'],(objs_pd_comb['M_star'][objs_pd_comb['type']=='Central'] + objs_pd_comb['mHI'][objs_pd_comb['type']=='Central'])/objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],edgecolor = 'k', facecolor = 'none',          s = ((objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'])**0.33/10*ms_scale).tolist(), linewidths = edgewidth)
    ax6.scatter(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central'],(objs_pd_comb['M_star'][objs_pd_comb['type']=='Central'] + objs_pd_comb['mHI'][objs_pd_comb['type']=='Central'])/objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],edgecolor = 'k', facecolor = 'k',alpha = 0.3, s = ((objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'])**0.33/10*ms_scale).tolist(), linewidths = edgewidth)
    ax6.scatter(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Backsplash'],(objs_pd_comb['M_star'][objs_pd_comb['type']=='Backsplash'] + objs_pd_comb['mHI'][objs_pd_comb['type']=='Backsplash'])/objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Backsplash'],edgecolor = 'k', facecolor = 'none',          s = ((objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Backsplash'])**0.33/10*ms_scale).tolist(), linewidths = edgewidth)
    ax6.scatter(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Backsplash'],(objs_pd_comb['M_star'][objs_pd_comb['type']=='Backsplash'] + objs_pd_comb['mHI'][objs_pd_comb['type']=='Backsplash'])/objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Backsplash'],edgecolor = 'k', facecolor = 'k',alpha = 0.3, s = ((objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Backsplash'])**0.33/10*ms_scale).tolist(), linewidths = edgewidth)
    #Baryonic mass fraction of satellites
    ax6.scatter(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite'],        (objs_pd_comb['M_star'][objs_pd_comb['type']=='Satellite'] + objs_pd_comb['mHI'][objs_pd_comb['type']=='Satellite'])/objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],facecolor = 'none',edgecolor = 'k',             marker = '*', s = ((objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'])**0.33/10*ms_scale*1.5).tolist(), linewidths = edgewidth)
    fbary = ax6.scatter(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite'],(objs_pd_comb['M_star'][objs_pd_comb['type']=='Satellite'] + objs_pd_comb['mHI'][objs_pd_comb['type']=='Satellite'])/objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],facecolor = 'k',   edgecolor = 'k',alpha = 0.3, marker = '*', s = ((objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'])**0.33/10*ms_scale*1.5).tolist(), linewidths = edgewidth)
    #Stellar mass fraction of central galaxies
    ax6.scatter(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central'],          objs_pd_comb['M_star'][objs_pd_comb['type']=='Central']/objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],    edgecolor = 'red',facecolor = 'none',               s = ((objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'])**0.33/10*ms_scale).tolist(),     linewidths = edgewidth)
    ax6.scatter(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Backsplash'],          objs_pd_comb['M_star'][objs_pd_comb['type']=='Backsplash']/objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Backsplash'],    edgecolor = 'red',facecolor = 'none',               s = ((objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Backsplash'])**0.33/10*ms_scale).tolist(),     linewidths = edgewidth)
    #Stellar mass fraction of satellites
    fstar = ax6.scatter(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite'],objs_pd_comb['M_star'][objs_pd_comb['type']=='Satellite']/objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],edgecolor = 'red',facecolor = 'none', marker = '*', s = ((objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'])**0.33/10*ms_scale*1.5).tolist(), linewidths = edgewidth)
    #ax6.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax6.set_xscale('log')
    ax6.set_yscale('log')
    ax6.set_ylabel(r'M/M$_{\mathrm{vir, peak}}$')
    ax6.set_xlabel(r"Distance to massive galaxy (kpc)")
    ax6.axis([17, 7e3, 1e-6, 0.08])
    ax6.legend([fstar,fbary],[r'M$_*$/M$_{\mathrm{vir, peak}}$',r'(M$_*$ + M$_{\mathrm{HI}})$/M$_{\mathrm{vir, peak}}$'],loc = 3)
    fig6.tight_layout()
    fig6.show()
    fig6.set_size_inches(plt_width,plt_width*aspect_ratio)
    fig6.show()
    fig6.savefig(outfile_base + '_fbary_rMassGal_Mpeak.png',dpi = dpi)

    #Color vs virial or stellar mass
    plt.close('all')
    fig7 = plt.figure(7,figsize=(plt_width,plt_width*aspect_ratio))
    fig7.clear()
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax7 = fig7.add_subplot(gs[0])
    ax7sub = fig7.add_subplot(gs[1]) 
    cNorm  = colors.Normalize(vmin=1.5, vmax = 4)    
    ax7.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],(objs_pd_comb['B-V'][objs_pd_comb['type']=='Central']),c = np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax7.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Backsplash'],(objs_pd_comb['B-V'][objs_pd_comb['type']=='Backsplash']),c = np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Backsplash']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Backsplash'], linewidths = edgewidth)
    ax7.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],(objs_pd_comb['B-V'][objs_pd_comb['type']=='Satellite']),c = np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k', marker = '*', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    ax7.set_xscale('log')
    ax7.set_ylabel(r'B-V')
    ax7.set_xlabel(r'M$_{\mathrm{vir}}$ [M$_\odot$]')
    ax7.axis([1e8, 1e11, 0, 0.8])
    cb = mpl.colorbar.ColorbarBase(ax7sub, cmap=cmx, norm=cNorm)
    cb.set_label(r"Log(D$_{\mathrm{massive}}$/1 kpc)")
    fig7.set_size_inches(plt_width,plt_width*aspect_ratio)
    fig7.tight_layout()
    fig7.show()
    fig7.savefig(outfile_base + '_BV_Mpeak.png',dpi = dpi)

    plt.close('all')
    fig7 = plt.figure(7,figsize=(plt_width,plt_width*aspect_ratio))
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax7 = fig7.add_subplot(gs[0])
    ax7sub = fig7.add_subplot(gs[1])
    cNorm  = colors.Normalize(vmin=1.5, vmax = 4)
    ax7.scatter(objs_pd_comb['mass'][objs_pd_comb['type']=='Central'],(objs_pd_comb['B-V'][objs_pd_comb['type']=='Central']),c = np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax7.scatter(objs_pd_comb['mass'][objs_pd_comb['type']=='Backsplash'],(objs_pd_comb['B-V'][objs_pd_comb['type']=='Backsplash']),c = np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Backsplash']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Backsplash'], linewidths = edgewidth)
    ax7.scatter(objs_pd_comb['mass'][objs_pd_comb['type']=='Satellite'],(objs_pd_comb['B-V'][objs_pd_comb['type']=='Satellite']),c = np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k', marker = '*', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    ax7.set_xscale('log')
    ax7.set_ylabel(r'B-V')
    ax7.set_xlabel(r'M$_{\mathrm{vir}}$ [M$_\odot$]')
    ax7.axis([1e8, 1e11, 0, 0.8])
    cb = mpl.colorbar.ColorbarBase(ax7sub, cmap=cmx, norm=cNorm)
    cb.set_label(r"Log(D$_{\mathrm{massive}}$/1 kpc)")
    fig7.set_size_inches(plt_width,plt_width*aspect_ratio)
    fig7.tight_layout()
    fig7.show()
    fig7.savefig(outfile_base + '_BV_Mvir.png',dpi = dpi)

    plt.close('all')
    fig7 = plt.figure(7,figsize=(plt_width,plt_width*aspect_ratio))
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax7 = fig7.add_subplot(gs[0])
    ax7sub = fig7.add_subplot(gs[1]) 
    cNorm  = colors.Normalize(vmin=1.5, vmax = 4)
    ax7.scatter(objs_pd_comb['Mstar_z0'][objs_pd_comb['type']=='Central'],(objs_pd_comb['B-V'][objs_pd_comb['type']=='Central']),c = np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax7.scatter(objs_pd_comb['Mstar_z0'][objs_pd_comb['type']=='Backsplash'],(objs_pd_comb['B-V'][objs_pd_comb['type']=='Backsplash']),c = np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Backsplash']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Backsplash'], linewidths = edgewidth)
    ax7.scatter(objs_pd_comb['Mstar_z0'][objs_pd_comb['type']=='Satellite'],(objs_pd_comb['B-V'][objs_pd_comb['type']=='Satellite']),c = np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k', marker = '*', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    ax7.set_xscale('log')
    ax7.set_ylabel(r'B-V')
    ax7.set_xlabel(r'M$_{\mathrm{star}}$ [M$_\odot$]')
    ax7.axis([1e3, 5e9, 0.05, 0.8])   
    cb = mpl.colorbar.ColorbarBase(ax7sub, cmap=cmx, norm=cNorm)
    cb.set_label(r"Log(D$_{\mathrm{massive}}$/1 kpc)")
    fig7.set_size_inches(plt_width,plt_width*aspect_ratio)
    fig7.tight_layout()
    fig7.show()
    fig7.savefig(outfile_base + '_BV_Mstar.png',dpi = dpi)

    plt.close('all')
    fig7 = plt.figure(7,figsize=(plt_width,plt_width*aspect_ratio))
    gs = gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax7 = fig7.add_subplot(gs[0])
    ax7sub = fig7.add_subplot(gs[1]) 
    cNorm  = colors.Normalize(vmin=1.5, vmax = 4)
    ax7.scatter(objs_pd_comb['Mstar_Mpeak'][objs_pd_comb['type']=='Central'],(objs_pd_comb['B-V'][objs_pd_comb['type']=='Central']),c = np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax7.scatter(objs_pd_comb['Mstar_Mpeak'][objs_pd_comb['type']=='Backsplash'],(objs_pd_comb['B-V'][objs_pd_comb['type']=='Backsplash']),c = np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Backsplash']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Backsplash'], linewidths = edgewidth)
    ax7.scatter(objs_pd_comb['Mstar_Mpeak'][objs_pd_comb['type']=='Satellite'],(objs_pd_comb['B-V'][objs_pd_comb['type']=='Satellite']),c = np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k', marker = '*', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    #ax5.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax7.set_xscale('log')
    ax7.set_ylabel(r'B-V')
    ax7.set_xlabel(r'M$_{\mathrm{star, peak}}$ [M$_\odot$]')
    ax7.axis([1e3, 5e9, 0.05, 0.8])
    cb = mpl.colorbar.ColorbarBase(ax7sub, cmap=cmx, norm=cNorm)
    cb.set_label(r"Log(D$_{\mathrm{massive}}$/1 kpc)")
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
    cNorm  = colors.Normalize(vmin=1.5, vmax = 4)
    ax7a.scatter(objs_pd_comb['Mstar_z0'][objs_pd_comb['type']=='Central'],(objs_pd_comb['B-V'][objs_pd_comb['type']=='Central']),c = np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax7a.scatter(objs_pd_comb['Mstar_z0'][objs_pd_comb['type']=='Backsplash'],(objs_pd_comb['B-V'][objs_pd_comb['type']=='Backsplash']),c = np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Backsplash']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Backsplash'], linewidths = edgewidth)
    ax7a.scatter(objs_pd_comb['Mstar_z0'][objs_pd_comb['type']=='Satellite'],(objs_pd_comb['B-V'][objs_pd_comb['type']=='Satellite']),c = np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k', marker = '*', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    #ax5.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax7a.set_ylabel(r'B-V')
    ax7a.set_xscale('log')
    ax7a.set_xlabel(r'M$_{\mathrm{star}}$ [M$_\odot$]')
    ax7a.axis([1e3, 5e9, 0.05, 0.85])    
    ax7b.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Central'],(objs_pd_comb['B-V'][objs_pd_comb['type']=='Central']),c = np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax7b.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Backsplash'],(objs_pd_comb['B-V'][objs_pd_comb['type']=='Backsplash']),c = np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Backsplash']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Backsplash'], linewidths = edgewidth)
    ax7b.scatter(objs_pd_comb['Mpeak'][objs_pd_comb['type']=='Satellite'],(objs_pd_comb['B-V'][objs_pd_comb['type']=='Satellite']),c = np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k', marker = '*', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    ax7b.set_xscale('log')
    ax7b.set_xlabel(r'M$_{\mathrm{vir, peak}}$ [M$_\odot$]')
    ax7b.axis([1e8, 1e11, 0.05, 0.85])
    fig7.tight_layout()
    sm = plt.cm.ScalarMappable(cmap=cmx, norm=cNorm)
    cbar = plt.colorbar(sm, location="right",ax=[ax7a,ax7b],aspect = 20) #pad=0.04)
    #cbar = mpl.colorbar.ColorbarBase(ax2sub, cmap=cmx, norm=cNorm, location="right",ax=[ax2a,ax2b]) #pad=0.04)
    cbar.set_label(r'Log(D$_{\mathrm{massive}}$/1 kpc)') #, rotation=0)
    fig7.show()
    fig7.set_size_inches(plt_width*2,plt_width*aspect_ratio)
    fig7.show()
    fig7.savefig(outfile_base + '_BV_Mvir_Mstar.png',dpi = dpi)

# Observational comparisons    
    lum = {'Lum': ""}
    objs_pd_comb = objs_pd_comb.join(pd.DataFrame(columns=lum))
    objs_pd_comb['Lum'] = 10**((objs_pd_comb['V_mag'] - 4.81)/-2.5)
    
    plt.close('all')
    fig8 = plt.figure(8,figsize=(plt_width,plt_width*aspect_ratio*2))
    gs = fig8.add_gridspec(2,hspace=0)
    axs8 = gs.subplots(sharex = True) #, constrained_layout=True) 
    axs8 = axs8.flatten()
    ax8a = axs8[0]
    ax8b = axs8[1]
    cmxr = plt.get_cmap("viridis_r")
    cNorm  = colors.Normalize(vmin=-20, vmax = -8)
    sm = plt.cm.ScalarMappable(cmap=cmxr, norm=cNorm)
    ax8a.scatter(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central'],(objs_pd_comb['B-V'][objs_pd_comb['type']=='Central']),c = objs_pd_comb['V_mag'][objs_pd_comb['type']=='Central'], cmap = cmxr, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax8a.scatter(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Backsplash'],(objs_pd_comb['B-V'][objs_pd_comb['type']=='Backsplash']),c = objs_pd_comb['V_mag'][objs_pd_comb['type']=='Backsplash'], cmap = cmxr, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Backsplash'], linewidths = edgewidth)
    ax8a.scatter(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite'],(objs_pd_comb['B-V'][objs_pd_comb['type']=='Satellite']),c = objs_pd_comb['V_mag'][objs_pd_comb['type']=='Satellite'], cmap = cmxr, norm = cNorm,edgecolor = 'k', marker = '*', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    #ax5.scatter(readdata['M200']*1e10,readdata['M_*']*1e7,marker = '+',c = 'k')
    ax8a.set_ylabel(r'B-V')
    ax8a.set_xscale('log')
    ax8a.axis([17, 7e3, 0.09, 0.81])
    #ax8a.set_xlabel(r'Distance to massive galaxy (kpc)')
    nogas = ~np.isfinite(np.array([objs_pd_comb['Lum']/objs_pd_comb['mHI']]))
    ax8b.scatter(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central'],(objs_pd_comb['mHI'][objs_pd_comb['type']=='Central']/objs_pd_comb['Lum'][objs_pd_comb['type']=='Central']),c = objs_pd_comb['V_mag'][objs_pd_comb['type']=='Central'], cmap = cmxr, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    colorVal = sm.to_rgba(objs_pd_comb['V_mag'][(objs_pd_comb['type']=='Central') & nogas])
    for (dist, c, ms) in zip(objs_pd_comb['massiveDist'][(objs_pd_comb['type']=='Central') & nogas], colorVal, objs_pd_comb['p_size'][(objs_pd_comb['type']=='Central') & nogas]):
        ax8b.errorbar(dist, 1e-4,yerr =  8e-5,ecolor = c, c=c, marker = "o",uplims = True,zorder = 1, capsize = 3)
        #ax8b.scatter(dist, 1e-4,c=c, marker = "o",edgecolor = 'k')
    ax8b.scatter(objs_pd_comb['massiveDist'][(objs_pd_comb['type']=='Central') & nogas],(objs_pd_comb['mHI'][(objs_pd_comb['type']=='Central') & nogas])*0 + 1e-4,c = objs_pd_comb['V_mag'][(objs_pd_comb['type']=='Central') & nogas],cmap = cmxr, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Central') & nogas], linewidths = edgewidth,zorder = 2)
    
    ax8b.scatter(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Backsplash'],(objs_pd_comb['mHI'][objs_pd_comb['type']=='Backsplash']/objs_pd_comb['Lum'][objs_pd_comb['type']=='Backsplash']),c = objs_pd_comb['V_mag'][objs_pd_comb['type']=='Backsplash'], cmap = cmxr, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Backsplash'], linewidths = edgewidth)
    colorVal = sm.to_rgba(objs_pd_comb['V_mag'][(objs_pd_comb['type']=='Backsplash') & nogas])
    for (dist, c) in zip(objs_pd_comb['massiveDist'][(objs_pd_comb['type']=='Backsplash') & nogas], colorVal):
        ax8b.errorbar(dist, 1e-4,yerr =  8e-5,ecolor = c, c=c, marker = "o", uplims = True,zorder = 1, capsize = 3)
    ax8b.scatter(objs_pd_comb['massiveDist'][(objs_pd_comb['type']=='Backsplash') & nogas],(objs_pd_comb['mHI'][(objs_pd_comb['type']=='Backsplash') & nogas])*0 + 1e-4,c = objs_pd_comb['V_mag'][(objs_pd_comb['type']=='Backsplash') & nogas], cmap = cmxr, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Backsplash') & nogas], linewidths = edgewidth,zorder = 2)
    
    ax8b.scatter(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite'],(objs_pd_comb['mHI'][objs_pd_comb['type']=='Satellite']/objs_pd_comb['Lum'][objs_pd_comb['type']=='Satellite']),c = objs_pd_comb['V_mag'][objs_pd_comb['type']=='Satellite'], cmap = cmxr, norm = cNorm,edgecolor = 'k', marker = '*', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    ax8b.scatter(objs_pd_comb['massiveDist'][(objs_pd_comb['type']=='Satellite') & nogas],(objs_pd_comb['mHI'][(objs_pd_comb['type']=='Satellite') & nogas])*0 + 1e-4,c = objs_pd_comb['V_mag'][(objs_pd_comb['type']=='Satellite') & nogas], cmap = cmxr, norm = cNorm,edgecolor = 'k', marker = '*', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Satellite') & nogas]*1.5, linewidths = edgewidth)
    colorVal = sm.to_rgba(objs_pd_comb['V_mag'][(objs_pd_comb['type']=='Satellite') & nogas])
    for (dist, c) in zip(objs_pd_comb['massiveDist'][(objs_pd_comb['type']=='Satellite') & nogas], colorVal):
        ax8b.errorbar(dist, 1e-4,yerr =  8e-5,ecolor = c, c=c, uplims = True, marker = '*',zorder = 1, capsize = 3)        
    ax8b.scatter(objs_pd_comb['massiveDist'][(objs_pd_comb['type']=='Satellite') & nogas],(objs_pd_comb['mHI'][(objs_pd_comb['type']=='Satellite') & nogas])*0 + 1e-4,c = objs_pd_comb['V_mag'][(objs_pd_comb['type']=='Satellite') & nogas], cmap = cmxr, norm = cNorm,edgecolor = 'k', marker = '*', s = objs_pd_comb['p_size'][(objs_pd_comb['type']=='Satellite') & nogas]*1.5, linewidths = edgewidth,zorder = 2)

    ax8b.scatter(MW_M31_dist, MW_M31_MHI_LV, marker='x', c=Calc_V_Mag, vmin=-20, vmax = -8, cmap=cmxr, s = markersize)
    #ax8b.scatter(MW_M31_dist, MW_M31_MHI_LV, marker=',', c='k') #, s = markersize*0.5)
    
    ax8b.set_xscale('log')
    ax8b.set_yscale('log')
    ax8b.set_xlabel(r'Log(D$_{\mathrm{massive}}$/1 kpc)')
    ax8b.set_ylabel(r'M$_{HI}$/L$_V$ [M$_{\odot}$/L$_{\odot}$]')
    ax8b.axis([17, 7e3, 1e-5, 30])
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
    cNorm  = colors.Normalize(vmin=1.5, vmax = 4)
    ax9a.scatter(objs_pd_comb['Mstar_z0'][objs_pd_comb['type']=='Central'],(objs_pd_comb['mZgas'][objs_pd_comb['type']=='Central']),c = np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax9a.scatter(objs_pd_comb['Mstar_z0'][objs_pd_comb['type']=='Backsplash'],(objs_pd_comb['mZgas'][objs_pd_comb['type']=='Backsplash']),c = np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Backsplash']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Backsplash'], linewidths = edgewidth)
    ax9a.scatter(objs_pd_comb['Mstar_z0'][objs_pd_comb['type']=='Satellite'],(objs_pd_comb['mZgas'][objs_pd_comb['type']=='Satellite']),c = np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k', marker = '*', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    #ax5.scatter(readdata['M200']*1e10,readdata['M_*']*1e9,marker = '+',c = 'k')
    ax9a.set_ylabel(r'Z$_{\mathrm{gas}}$ [M$_\odot$]')
    ax9a.set_xscale('log')
    ax9a.set_yscale('log')
    ax9a.axis([1e2, 5e9, 1, 1e8])    
    ax9b.scatter(objs_pd_comb['Mstar_z0'][objs_pd_comb['type']=='Central'],(objs_pd_comb['mZgas'][objs_pd_comb['type']=='Central']/objs_pd_comb['mgas'][objs_pd_comb['type']=='Central']),c = np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Central']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Central'], linewidths = edgewidth)
    ax9b.scatter(objs_pd_comb['Mstar_z0'][objs_pd_comb['type']=='Backsplash'],(objs_pd_comb['mZgas'][objs_pd_comb['type']=='Backsplash'])/(objs_pd_comb['mgas'][objs_pd_comb['type']=='Backsplash']),c = np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Backsplash']), cmap = cmx, norm = cNorm,edgecolor = 'k', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Backsplash'], linewidths = edgewidth)
    ax9b.scatter(objs_pd_comb['Mstar_z0'][objs_pd_comb['type']=='Satellite'],(objs_pd_comb['mZgas'][objs_pd_comb['type']=='Satellite'])/(objs_pd_comb['mgas'][objs_pd_comb['type']=='Satellite']),c = np.log10(objs_pd_comb['massiveDist'][objs_pd_comb['type']=='Satellite']), cmap = cmx, norm = cNorm,edgecolor = 'k', marker = '*', s = objs_pd_comb['p_size'][objs_pd_comb['type']=='Satellite']*1.5, linewidths = edgewidth)
    ax9b.set_xscale('log')
    ax9b.set_yscale('log')
    ax9b.set_xlabel(r'M$_{\mathrm{star}}$ [M$_\odot$]')
    #ax9b.axis([1e8, 1e11, 0.05, 0.85])
    fig9.tight_layout()
    sm = plt.cm.ScalarMappable(cmap=cmx, norm=cNorm)
    cbar = plt.colorbar(sm, location="right",ax=[ax9a,ax9b],aspect = 20) #pad=0.04)
    #cbar = mpl.colorbar.ColorbarBase(ax2sub, cmap=cmx, norm=cNorm, location="right",ax=[ax2a,ax2b]) #pad=0.04)
    cbar.set_label(r'Log(D$_{\mathrm{massive}}$/1 kpc)') #, rotation=0)
    fig9.show()
    fig9.set_size_inches(plt_width,plt_width*aspect_ratio*2)
    fig9.show()
    fig9.savefig(outfile_base + '_BV_Mvir_Mstar.png',dpi = dpi)    
        
if __name__ == '__main__':
    if (socket.gethostname() == "quirm.math.grinnell.edu"):
        prefix = '/home/christenc/Data/Sims/'
        prefix_outfile = '/home/christenc/Figures/marvel/'
    else:
        prefix = '/home/christensen/Storage2/UW/MolecH/Cosmo/'
        prefix_outfile = '/home/christensen/Plots/marvel/'

    tfile_base_cm = 'cptmarvel.cosmo25cmb.4096g5HbwK1BH'
    tfile_cm = prefix + 'cptmarvel.cosmo25cmb/cptmarvel.cosmo25cmb.4096g5HbwK1BH/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096' #'cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096'

    tfile_r = prefix + 'rogue.cosmo25cmb/rogue.cosmo25cmb.4096g5HbwK1BH/rogue.cosmo25cmb.4096g5HbwK1BH.004096/rogue.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_base_r = 'rogue.cosmo25cmb.4096g5HbwK1BH'

    tfile_e = prefix + 'elektra.cosmo25cmb/elektra.cosmo25cmb.4096g5HbwK1BH/elektra.cosmo25cmb.4096g5HbwK1BH.004096/elektra.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_base_e = 'elektra.cosmo25cmb.4096g5HbwK1BH'

    tfile_s = prefix + 'storm.cosmo25cmb/storm.cosmo25cmb.4096g5HbwK1BH/storm.cosmo25cmb.4096g5HbwK1BH.004096/storm.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_base_s = 'storm.cosmo25cmb.4096g5HbwK1BH'
    
    tfile_base_1 = 'h148.cosmo50PLK.3072g3HbwK1BH'
    tfile_1 = prefix + 'h148.cosmo50PLK.3072g/h148.cosmo50PLK.3072g3HbwK1BH/snapshots_200bkgdens/h148.cosmo50PLK.3072g3HbwK1BH.004096' #
    tfile_1 = prefix + 'h148.cosmo50PLK.3072g/h148.cosmo50PLK.3072g3HbwK1BH/h148.cosmo50PLK.3072g3HbwK1BH.004096/h148.cosmo50PLK.3072g3HbwK1BH.004096'
    
    tfile_base_1hr = 'h148.cosmo50PLK.6144g3HbwK1BH/'
    tfile_1hr = prefix + 'h148.cosmo50PLK.6144g/h148.cosmo50PLK.6144g3HbwK1BH/h148.cosmo50PLK.6144g3HbwK1BH.004096/ahf_200/h148.cosmo50PLK.6144g3HbwK1BH.004096'
    
    tfile_base_2 = 'h229.cosmo50PLK.3072gst5HbwK1BH'
    tfile_2 = prefix + 'h229.cosmo50PLK.3072g/h229.cosmo50PLK.3072gst5HbwK1BH/snapshots_200bkgdens/h229.cosmo50PLK.3072gst5HbwK1BH.004096' #
    tfile_2 = prefix + 'h229.cosmo50PLK.3072g/h229.cosmo50PLK.3072gst5HbwK1BH/h229.cosmo50PLK.3072gst5HbwK1BH.004096/h229.cosmo50PLK.3072gst5HbwK1BH.004096'    
    
    tfile_base_3 = 'h242.cosmo50PLK.3072gst5HbwK1BH'
    tfile_3 = prefix + 'h242.cosmo50PLK.3072g/h242.cosmo50PLK.3072gst5HbwK1BH/snapshots_200bkgdens/h242.cosmo50PLK.3072gst5HbwK1BH.004096' #
    tfile_3 = prefix + 'h242.cosmo50PLK.3072g/h242.cosmo50PLK.3072gst5HbwK1BH/h242.cosmo50PLK.3072gst5HbwK1BH.004096/h242.cosmo50PLK.3072gst5HbwK1BH.004096'

    tfile_base_4 = 'h329.cosmo50PLK.3072gst5HbwK1BH'
    tfile_4 = prefix + 'h329.cosmo50PLK.3072g/h329.cosmo50PLK.3072gst5HbwK1BH/snapshots_200bkgdens/h329.cosmo50PLK.3072gst5HbwK1BH.004096' #
    tfile_4 = prefix + 'h329.cosmo50PLK.3072g/h329.cosmo50PLK.3072gst5HbwK1BH/h329.cosmo50PLK.3072gst5HbwK1BH.004096/h329.cosmo50PLK.3072gst5HbwK1BH.004096'
    
    tfile_base_4hr = 'h329.cosmo50PLK.6144g5HbwK1BH'
    tfile_4hr = prefix + 'h329.cosmo50PLK.6144g/h329.cosmo50PLK.6144g5HbwK1BH/h329.cosmo50PLK.6144g5HbwK1BH.004096/ahf_200/h329.cosmo50PLK.6144g5HbwK1BH.004096' #    
    
    outfile_base = prefix_outfile + 'marvelJL'
    #SMHMv__distance([tfile_1,tfile_2,tfile_3,tfile_4],outfile_base,[tfile_base_1,tfile_base_2,tfile_base_3,tfile_base_4])
    tfiles = [tfile_cm, tfile_e, tfile_r, tfile_s, tfile_1, tfile_1hr, tfile_2, tfile_3, tfile_4, tfile_4hr]
    tfile_base = [tfile_base_cm, tfile_base_e, tfile_base_r, tfile_base_s, tfile_base_1, tfile_base_1hr, tfile_base_2, tfile_base_3, tfile_base_4, tfile_base_4hr]
    SMHM_v_distance([tfile_cm, tfile_e, tfile_r, tfile_s, tfile_1, tfile_1hr, tfile_2, tfile_3, tfile_4, tfile_4hr],outfile_base,[tfile_base_cm, tfile_base_e, tfile_base_r, tfile_base_s, tfile_base_1, tfile_base_1hr, tfile_base_2, tfile_base_3, tfile_base_4, tfile_base_4hr])

