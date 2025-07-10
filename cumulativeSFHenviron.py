#Charlotte Christensen

# 8/2/19
#For the given simulations, calculate the star formation histories at different distances from a massive galaxy

#%run /home/christenc/Code/python/python_analysis/cumulativeSFHenviron
import matplotlib as mpl
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
import tangos
import sys, os, glob, pickle
sys.path.append(os.path.abspath("/home/christenc/Code/python/python_analysis/"))
from modules.user_tools import task

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
                #print(sim, halo, fdmdata.loc[index_best_match]['halogrp_z0'], float(objs_pd[(objs_pd['sim'] == sim) & (objs_pd['haloid'] == halo)]['M_star']), fdmdata.loc[index_best_match]['Mstar_z0'], float(objs_pd[(objs_pd['sim'] == sim) & (objs_pd['haloid'] == halo)]['mass']), fdmdata.loc[index_best_match]['Mhalo_z0'])
                fdmdata.loc[index_best_match,'Mstar_z0'] = 0 #set stellar mass to zero so it isn't used again

    return objs_pd

#Calculates the distance to the nearest massive galaxy
def distance_to_nearest_host(data, tfiles):
    massiveDist = []
    massiveDist_ID = []
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
                haloid = []
                for AHFhalo in h_dummy:
                    properties = AHFhalo.properties            
                    if (properties['mass'] > min_massiveHalo):
#                        print('Halo id: ',properties['halo_id'])
                        loc.append(np.array([properties['Xc']/properties['h'], properties['Yc']/properties['h'], properties['Zc']/properties['h']]))
                        rvir.append(properties['Rvir'])
                        haloid.append(properties['halo_id'])
                loc = np.array(loc)
                rvir = np.array(rvir)
                #for index, halo in data.iterrows():

            properties = h_dummy[int(data['haloid'].tolist()[i])].properties
            massiveDist.append(min(((properties['Xc']/properties['h'] - loc[:,0])**2 + (properties['Yc']/properties['h'] - loc[:,1])**2 + (properties['Zc']/properties['h'] - loc[:,2])**2)**(0.5)))
            #minind = np.where((((properties['Xc'] - loc[:,0])**2 + (properties['Yc'] - loc[:,1])**2 + (properties['Zc'] - loc[:,2])**2)**(0.5)) == massiveDist[-1])
            minind = np.argmin(((properties['Xc']/properties['h'] - loc[:,0])**2 + (properties['Yc']/properties['h'] - loc[:,1])**2 + (properties['Zc']/properties['h'] - loc[:,2])**2)**(0.5))
            massiveDist_ID.append(haloid[minind])
            #hostrvirs.append(rvir[minind]*0.6776942783267969)

    data['massiveDist'] = massiveDist
    data['massiveDist_ID'] = massiveDist_ID
    return data

def cumulativeSFHenviron(tfiles,outfile_base,tfile_base,*halo_nums):
    presentation = False
    if presentation:
        outfile_base = outfile_base + '_pres'
        plt.style.use(['default','/home/christenc/.config/matplotlib/presentation.mplstyle'])
        plt_width = 16 #8 #inches
        aspect_ratio = 1.0/4.0
        legendsize = 16
        dpi = 100
        lw = mpl.rcParams['lines.linewidth'] - 1        
    else:
        plt.style.use(['default','/home/christenc/.config/matplotlib/article.mplstyle'])
        plt_width = 8 #inches
        aspect_ratio = 1.0/4.0
        legendsize = 5
        dpi = 300
        lw = mpl.rcParams['lines.linewidth']

    if (socket.gethostname() == "ozma.grinnell.edu"):
        dataprefix = '/home/christensen/Code/Datafiles/' 
    else:
        dataprefix = '/home/christenc/Code/Datafiles/'
        
    #min_mass = 1e9 #Minimum halo mass for analysis in solar masses
    min_nstar =  100 #100 #Minimum number of stars for
    min_npart = 100
    min_noncontamFrac = 0.9
    min_massiveHalo = 10**11.5

    use_m200 = 1
    if use_m200:
        ext = '.m200.dist.'
    else:
        ext = '.MAP.'
        
    objs_pd = None 
    for tfile, base in zip(tfiles, tfile_base):
        objs_dat = []
        print(tfile)
        f=open(tfile + ext + 'data', 'rb')
        while 1:
            try:
                objs_dat.append(pickle.load(f))
            except EOFError:
                break        
        f.close()
        
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
#        if not 'massiveDist' in temp:
        if 1:
            temp = distance_to_nearest_host(temp,[tfile])
            temp.to_pickle(tfile + ext + 'data')
            temp.to_csv(tfile + ext + 'data.csv', index=False)            
        
        s = pynbody.load(tfile)
        h_dummy = s.halos(dummy = True)
        loc = []
        fMhires_ahf = []
        nstar_ahf = []
        haloid_ahf = []
        npart_ahf = []
        for AHFhalo in h_dummy:
            properties = AHFhalo.properties
            fMhires_ahf.append(AHFhalo.properties['fMhires'])
            nstar_ahf.append(AHFhalo.properties['n_star'])
            haloid_ahf.append(AHFhalo.properties['halo_id'])
            npart_ahf.append(AHFhalo.properties['npart'])
            #if (properties['mass'] > min_massiveHalo):
            #    loc.append(np.array([properties['Xc']/properties['h'], properties['Yc']/properties['h'], properties['Zc']/properties['h']]))
        #valid_halos_ahf = (np.array(haloid_ahf))[(np.array(fMhires_ahf) >= min_noncontamFrac) & (np.array(nstar_ahf) >= min_nstar)] # Minimum star particles
        valid_halos_ahf = (np.array(haloid_ahf))[(np.array(fMhires_ahf) >= min_noncontamFrac) & (np.array(npart_ahf) >= min_npart)]  #  Minimum total particles
        """
        loc = np.array(loc)
        massiveDist = [] #Distance to nearest massive (M_vir > min_massiveHalo) galaxy
        for halo in objs_dat:
            massiveDist.append(min(((halo['Xc'] - loc[:,0])**2 + (halo['Yc'] - loc[:,1])**2 + (halo['Zc'] - loc[:,2])**2)**(0.5)))
        """
        
        sim = ''
        if tfile.split('/')[-1] == 'cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096':
            sim = 'Cpt. Marvel'
        if tfile.split('/')[-1] == 'rogue.cosmo25cmb.4096g5HbwK1BH.004096':
            sim = 'Rogue'
        if tfile.split('/')[-1] == 'elektra.cosmo25cmb.4096g5HbwK1BH.004096':
            sim = 'Elektra'
        if tfile.split('/')[-1] == 'storm.cosmo25cmb.4096g5HbwK1BH.004096':
            sim = 'Storm'
        if tfile.split('/')[-1] == 'h329.cosmo50PLK.3072gst5HbwK1BH.004096':
            sim = 'Elena'
        if tfile.split('/')[-1] == 'h242.cosmo50PLK.3072gst5HbwK1BH.004096':
            sim = 'Sonia'
        if tfile.split('/')[-1] == 'h229.cosmo50PLK.3072gst5HbwK1BH.004096':
            sim = 'Ruth'
        if tfile.split('/')[-1] == 'h148.cosmo50PLK.3072g3HbwK1BH.004096':
            sim = 'Sandra'
        if tfile.split('/')[-1] == 'h329.cosmo50PLK.6144g5HbwK1BH.004096':
            sim = 'Elena Mint'
        if tfile.split('/')[-1] == 'h148.cosmo50PLK.6144g3HbwK1BH.004096':
            sim = 'Sandra Mint'

        #temp = pd.DataFrame(objs_dat)

        #temp['sim'] = [sim]*len(objs_dat)
        #temp['massiveDist'] = massiveDist
        """
        valid_halos = np.array(temp[(temp['n_star'] >= min_nstar) &(temp['fMhires']>=min_noncontamFrac )]['haloid'])
        intersec = np.intersect1d(valid_halos,valid_halos_ahf)
        if len(valid_halos) is not len(intersec):
            print(tfile)
            print('Number of valid halos in *.data: ',len(valid_halos))
            print('Extra halos in *.dat')
            print(np.setdiff1d(valid_halos,intersec))
            print('Number of valid halos in AHF: ',len(valid_halos_ahf))
            print('Extra halos in AHF')
            print(np.setdiff1d(valid_halos_ahf,intersec))
        """
        if objs_pd is None: 
            objs_pd = temp
        else:
            objs_pd = objs_pd.append(temp, ignore_index = True)

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

    assem_history = task(
        np.load,
        start_text="Loading assembly history data",
        end_text="Loaded assembly history data",
        fail_text="Failed to load assembly history data",
        exit_on_fail=True
    )(dataprefix + "/assembly_histories.npy", encoding="latin1", allow_pickle=True).item()
    assem_history = pd.DataFrame.from_dict(assem_history)
    
    #Match halos between my and Ferah's data
    fdmdata_mod = fdmdata.copy()
    objs_pd = match_halos(objs_pd, fdmdata_mod)
    #remove entries (rows) from objs_pd that have no match in Ferah's data
    index_rm = objs_pd[(objs_pd['m200_haloid']).isnull()].index
    print('index rm:', index_rm)
    objs_pd = objs_pd.drop(index_rm)
    
    halo_label = {'halo_label': ""}
    objs_pd = objs_pd.join(pd.DataFrame(columns=halo_label))
    objs_pd.loc[~objs_pd['m200_haloid'].isnull(),'halo_label'] = objs_pd[~objs_pd['m200_haloid'].isnull()]['sim']+objs_pd[~objs_pd['m200_haloid'].isnull()]['m200_haloid'].astype(str)
    objs_pd = objs_pd.set_index('halo_label')

    fdmdata = fdmdata.join(pd.DataFrame(columns=halo_label))
    fdmdata['halo_label'] = fdmdata['simname']+fdmdata['halogrp_z0'].astype(str)
    fdmdata = fdmdata.set_index('halo_label') 

    objs_pd_comb = pd.concat([objs_pd,fdmdata], join="inner", axis=1)

    #### SET plot parameters
    maxtime = 13.7

    use_min_dist = 1
    
    mass_color_scale = 1
    use_virial_mass = 1 #1
    use_maxv_mass = 1
    if mass_color_scale:
        if use_virial_mass:
            min_vmass = 7e7 #1e8
            max_vmass = 1e11 #3e12
            outfile_base_ext = outfile_base + '_vmass'
        elif use_maxv_mass:
            min_vmass = 1e8 #1e8
            max_vmass = 1e11 #3e12
            outfile_base_ext = outfile_base + '_vmmax'            
        else: #use stellar mass instead
            min_vmass = 1e3 #min_nstar*5e3 #1e8
            max_vmass = 2.5e9 #3e12
            outfile_base_ext = outfile_base + '_smass'
        outfile_base_ext = outfile_base_ext + '_masscolor'
        cmx = plt.get_cmap("winter_r") 
            #values = range(0,20)
            #cNorm  = colors.Normalize(vmin=0, vmax = values[-1]+1)
            #scalarMap = cm.ScalarMappable(norm=cNorm, cmap=cmx)
            
        #cNorm  = colors.Normalize(vmin=np.log10(min_vmass), vmax = math.ceil(np.log10(max_vmass)))
        cNorm  = colors.Normalize(vmin=np.log10(min_vmass), vmax = np.log10(max_vmass))
        scalarMap = cm.ScalarMappable(norm=cNorm, cmap=cmx)     
    else:
        cmx = plt.get_cmap("tab20")
        values = range(0,20)
        cNorm  = colors.Normalize(vmin=0, vmax = values[-1]+1)
        scalarMap = cm.ScalarMappable(norm=cNorm, cmap=cmx)
    
    objs_pd_comb_bk = objs_pd_comb.copy()
    objs_pd_comb = objs_pd_comb_bk[objs_pd_comb_bk['Mstar_z0'] < 1e6].copy()
    objs_pd_comb[objs_pd_comb[dist_key]> 320]['mISM']/objs_pd_comb[objs_pd_comb[dist_key]> 320][dist_key]**2*1e6
    
# Plot histories
    plt.clf()
    plt.close('all')
    fig1 = plt.figure(1,figsize=(plt_width,plt_width*aspect_ratio))
    fig1.set_size_inches(plt_width,plt_width*aspect_ratio)
    if mass_color_scale:
        gs =  gridspec.GridSpec(1,2,width_ratios=[45,1])
    else:
        gs =  gridspec.GridSpec(1,1)
    
    gs00 = gridspec.GridSpecFromSubplotSpec(1, 6, subplot_spec=gs[0],wspace=0.0)
    ax1 = fig1.add_subplot(gs00[0])
    ax2 = fig1.add_subplot(gs00[1],sharey=ax1)
    ax3 = fig1.add_subplot(gs00[2],sharey=ax1)
    ax4 = fig1.add_subplot(gs00[3],sharey=ax1)
    ax5 = fig1.add_subplot(gs00[4],sharey=ax1)
    if mass_color_scale:
        ax1sub = fig1.add_subplot(gs[1])
    axes = [ax1,ax2,ax3,ax4,ax5]

    fig2 = plt.figure(2,figsize=(plt_width,plt_width*aspect_ratio))
    fig2.set_size_inches(plt_width,plt_width*aspect_ratio)
    if mass_color_scale:
        gsa =  gridspec.GridSpec(1,2,width_ratios=[45,1])
    else:
        gsa =  gridspec.GridSpec(1,1)
    
    gs00a = gridspec.GridSpecFromSubplotSpec(1, 6, subplot_spec=gsa[0],wspace=0.0)
    ax1a = fig2.add_subplot(gs00a[0])
    ax2a = fig2.add_subplot(gs00a[1],sharey=ax1)
    ax3a = fig2.add_subplot(gs00a[2],sharey=ax1)
    ax4a = fig2.add_subplot(gs00a[3],sharey=ax1)
    ax5a = fig2.add_subplot(gs00a[4],sharey=ax1)
    if mass_color_scale:
        ax1asub = fig2.add_subplot(gsa[1])
    axesa = [ax1a,ax2a,ax3a,ax4a,ax5a]
    
    legendlines = []

    #maxradii_arr = [150,300,750,7000]
    #minradii_arr = [0.1,maxradii_arr[0],maxradii_arr[1],maxradii_arr[2]]
    #titles = ['r < 150 kpc','150 kpc < r < 300 kpc','300 kpc < r < 750 kpc','750 kpc < r']
    
    if use_min_dist:
        maxradii_arr = (10**np.array([1.0, 1.5, 2, 2.5, 3, 4]))[1:] # everything is actually below 10^3.5 kpc
        minradii_arr = (10**np.array([1.0, 1.5, 2, 2.5, 3, 4]))[:-1]
        titles = ['$r_{\mathrm{min}} < 32$ kpc', '32 $< r_{\mathrm{min}} < 100$ kpc','100 $< r_{\mathrm{min}} < 320$ kpc','320 $< r_{\mathrm{min}} < 1000$ kpc','$r_{\mathrm{min}} >$ 1000 kpc']        
        outfilename = outfile_base_ext + '_rClosest_sfhcum.png'
        dist_key = 'min_dist'
    else:
        maxradii_arr = (10**np.array([1.0, 2, 2.5, 3, 3.5, 4]))[1:] # everything is actually above 10^1.5 kpc
        minradii_arr = (10**np.array([1.0, 2, 2.5, 3, 3.5, 4]))[:-1]
        titles = ['r < 100 kpc','100 kpc < r < 320 kpc','320 kpc < r < 1 Mpc','1 Mpc < r < 3.2 Mpc','3.2 Mpc < r']
        outfilename = outfile_base_ext + '_rMassGal_sfhcum.png'
        dist_key = 'massiveDist'
        
    for minradii, maxradii, ax, axa, title in zip(minradii_arr, maxradii_arr, axes, axesa, titles):
        i = 0
        labels = []
        print(maxradii, len(objs_pd_comb[(objs_pd_comb[dist_key] < maxradii) & (objs_pd_comb[dist_key] >= minradii)]))
        print(max(objs_pd_comb[(objs_pd_comb[dist_key] < maxradii) & (objs_pd_comb[dist_key] >= minradii)][dist_key]))
        print(min(objs_pd_comb[(objs_pd_comb[dist_key] < maxradii) & (objs_pd_comb[dist_key] >= minradii)][dist_key]))
        for index, halo in objs_pd_comb[(objs_pd_comb[dist_key] < maxradii) & (objs_pd_comb[dist_key] >= minradii)].iterrows(): #len(h)-1):
            #if (halo['n_star'] < min_nstar):
            #    continue
            #if (halo['fMhires'] < min_noncontamFrac):
            #    continue
            #if (halo['n_star'] > 10*(halo['n_particles'] - halo['n_star'] - halo['n_gas'])):
            #    continue            

            quenched = 0
            lw_plot = lw #*2
            linestyle = '-'
                          
            #if (halo['SFR']/halo['M_star'] > 1e-11):
            if (halo['type'] == 'Central'):
                quenched = 0
                lw_plot = lw #*2
                linestyle = '-'
            else:
                quenched = 1
                lw_plot = lw
                linestyle = '--'
            
            
            #            if (halo['hostVirialR'] > 0):
            #                linestyle = '--'
            #            else:
            #                linestyle = '-'
                
            labels.append(str(halo['haloid']) + ' ' + halo['sim'])
            dtime = halo['sfhbins'][1] - halo['sfhbins'][0] #Gyr
            sfhcum = np.cumsum(halo['sfh']*dtime)
            
            if mass_color_scale:
                if use_virial_mass:
                    colorVal = scalarMap.to_rgba(np.log10(halo['mass']))
                elif use_maxv_mass:
                    colorVal = scalarMap.to_rgba(np.log10(halo['Mpeak']))                    
                else:
                    colorVal = scalarMap.to_rgba(np.log10(halo['M_star']))
                legendline = ax.plot((halo['sfhbins'][:-1] + halo['sfhbins'][1:])/2,sfhcum/max(sfhcum),color = colorVal,linestyle = linestyle,linewidth = lw_plot)
            else:
                colorVal = scalarMap.to_rgba((i % 20) + 1)
                legendline = ax.plot((halo['sfhbins'][:-1] + halo['sfhbins'][1:])/2,sfhcum/max(sfhcum),color = colorVal,linestyle = linestyle,linewidth = lw_plot) #,c = i) #, cmap = 'tab20')

            ind5 = np.argmin(np.abs((halo['sfhbins'][:-1] + halo['sfhbins'][1:])/2 - 5))

            #if (((halo['sfhbins'][:-1] + halo['sfhbins'][1:])/2)[-1] > 5) and (sfhcum[ind5]/max(sfhcum) < 0.1):
            if True:
                print("Late Starting SF: ",halo['sim'],halo['haloid'],halo[dist_key],np.log10(halo['Mstar_z0']),np.log10(halo['Mhalo_z0']))
                if (halo['sim'] == 'cptmarvel') or (halo['sim'] == 'elektra') or (halo['sim'] == 'rogue') or (halo['sim'] == 'storm'):
                    tangos_db = '/home/christenc/Storage/tangos_db/Marvel_r200.db'
                    tangos_db = '/home/christenc/Storage/tangos_db/Marvel_r200_N100min.db'
                    sim_key = halo['sim']
                elif (halo['sim'] == 'h148') or (halo['sim'] == 'h229') or (halo['sim'] == 'h242') or (halo['sim'] == 'h329'):
                    tangos_db = '/home/christenc/Storage/tangos_db/JL_r200_N100min.db'
                    sim_key = halo['sim']
                elif (halo['sim'] == 'h148_6144') or (halo['sim'] == 'h329_6144'):
                    tangos_db = '/home/christenc/Storage/tangos_db/JLmint_r200_N100.db'
                    sim_key = halo['sim'][0:4]+'mint'
                else:
                    print("Database not available")
                tangos.init_db(tangos_db)
                halot = tangos.get_halo("snapshots_200crit_" + sim_key + "/%4096/halo_" + halo['haloid'])
                halot1 = tangos.get_halo("snapshots_200crit_" + sim_key + "/%4096/halo_" + str(halo['massiveDist_ID']))
                Xc, Yc, Zc, t = halot.calculate_for_progenitors("Xc","Yc","Zc","t()")
                Xc1, Yc1, Zc1, t1 = halot1.calculate_for_progenitors("Xc","Yc","Zc", "t()")
                if len(Xc1) >= len(Xc):
                    dist_h1 = np.sqrt((Xc - Xc1[-1*len(Xc):])**2 + (Yc - Yc1[-1*len(Xc):])**2 +(Zc - Zc1[-1*len(Xc):])**2)
                    timet = t
                else:
                    dist_h1 = np.sqrt((Xc[-1*len(Xc1):] - Xc1)**2 + (Yc[-1*len(Xc1):] - Yc1)**2 +(Zc[-1*len(Xc1):] - Zc1)**2)
                    timet = t1                    
                
                sim_name = halo['sim']
                haloid_name = halo['haloid']
                cond = ((assem_history['Volume'] == sim_name) & (assem_history['halo grp @ z=0'] == int(haloid_name)))
                if np.sum(cond) == 0:
                    print("Missing halo: ",sim_name,int(haloid_name))
                    continue
                xarr = np.array(((assem_history['Time'][cond]).tolist())[0])
                assem_history_halo = np.array(((assem_history['Assembly History'][cond]).tolist())[0])
                star_history_halo = np.array(((assem_history['Stellar Mass History'][cond]).tolist())[0])
                gas_history_halo = np.array(((assem_history['Gas Mass History'][cond]).tolist())[0])
                yarr = (assem_history_halo)/np.max(assem_history_halo)
                axa.plot(xarr, yarr,color = colorVal)
                axa.plot(xarr, star_history_halo/np.max(star_history_halo), linestyle="--",color = colorVal)
                axa.plot(xarr, gas_history_halo/np.max(gas_history_halo), linestyle="-.",color = colorVal)
                axa.plot(timet, dist_h1/1e3,linestyle = ':', color = colorVal)
            ax.label_outer()
            legendlines.append(legendline)
            fig1.show()
            ax.set_ylim([0,1])
            ax.set_title(title,fontsize= 'x-small')
            i = i + 1
        if not mass_color_scale:
            ax.legend(labels,loc = 1,fontsize = 4) #'xx-small')
            
    ax3.set_xlabel('Time [Gyr]')
    ax1.set_ylabel(r'$M_*( t_{form} < t)/M_*$')
    #ax1.axis([0.1, 20, 0, 1])
    #ax1.set_xscale('log')
    #ax1.set_title(tfile_base)
    if mass_color_scale:
        cb = mpl.colorbar.ColorbarBase(ax1sub, cmap=cmx, norm=cNorm)
        if use_virial_mass:
            cb.set_label('Log(M$_{\mathrm{vir}}$/1 M$_\odot$)')
        elif use_maxv_mass:
            cb.set_label('Log(Peak M$_{\mathrm{vir}}$/1 M$_\odot$)')            
        else:
            cb.set_label('Log(M$_{\mathrm{*}}$/1 M$_\odot$)')

    plt.tight_layout()
    fig1.tight_layout()
    fig1.show()
    fig1.savefig(outfilename,dpi = dpi)
    fig1.clear()

    
if __name__ == '__main__':
    if (socket.gethostname() == "quirm.math.grinnell.edu"):
        prefix = '/home/christenc/Data/Sims/'
        prefix_outfile = '/home/christenc/Figures/marvel/'
    else:
        prefix = '/home/christensen/Storage2/UW/MolecH/Cosmo/'
        prefix_outfile = '/home/christensen/Plots/marvel/'

    tfile_base_cm = 'cptmarvel.cosmo25cmb.4096g5HbwK1BH'
    tfile_cm = prefix + 'cptmarvel.cosmo25cmb/cptmarvel.cosmo25cmb.4096g5HbwK1BH/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096' #'cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_cm = prefix + 'cptmarvel.cosmo25cmb/cptmarvel.cosmo25cmb.4096g5HbwK1BH/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096/ahf_200/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096'

    tfile_r = prefix + 'rogue.cosmo25cmb/rogue.cosmo25cmb.4096g5HbwK1BH/rogue.cosmo25cmb.4096g5HbwK1BH.004096/rogue.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_r = prefix + 'rogue.cosmo25cmb/rogue.cosmo25cmb.4096g5HbwK1BH/rogue.cosmo25cmb.4096g5HbwK1BH.004096/ahf_200/rogue.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_base_r = 'rogue.cosmo25cmb.4096g5HbwK1BH'

    tfile_e = prefix + 'elektra.cosmo25cmb/elektra.cosmo25cmb.4096g5HbwK1BH/elektra.cosmo25cmb.4096g5HbwK1BH.004096/elektra.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_e = prefix + 'elektra.cosmo25cmb/elektra.cosmo25cmb.4096g5HbwK1BH/elektra.cosmo25cmb.4096g5HbwK1BH.004096/ahf_200/elektra.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_base_e = 'elektra.cosmo25cmb.4096g5HbwK1BH'

    tfile_s = prefix + 'storm.cosmo25cmb/storm.cosmo25cmb.4096g5HbwK1BH/storm.cosmo25cmb.4096g5HbwK1BH.004096/storm.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_s = prefix + 'storm.cosmo25cmb/storm.cosmo25cmb.4096g5HbwK1BH/storm.cosmo25cmb.4096g5HbwK1BH.004096/ahf_200/storm.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_base_s = 'storm.cosmo25cmb.4096g5HbwK1BH'
    
    tfile_base_1 = 'h148.cosmo50PLK.3072g3HbwK1BH'
    tfile_1 = prefix + 'h148.cosmo50PLK.3072g/h148.cosmo50PLK.3072g3HbwK1BH/snapshots_200bkgdens/h148.cosmo50PLK.3072g3HbwK1BH.004096' #
    tfile_1 = prefix + 'h148.cosmo50PLK.3072g/h148.cosmo50PLK.3072g3HbwK1BH/h148.cosmo50PLK.3072g3HbwK1BH.004096/ahf_200/h148.cosmo50PLK.3072g3HbwK1BH.004096'
    #tfile_1 = prefix + 'h148.cosmo50PLK.3072g/h148.cosmo50PLK.3072g3HbwK1BH/h148.cosmo50PLK.3072g3HbwK1BH.004096/h148.cosmo50PLK.3072g3HbwK1BH.004096'    

    tfile_base_1hr = 'h148.cosmo50PLK.6144g3HbwK1BH/'
    tfile_1hr = prefix + 'h148.cosmo50PLK.6144g/h148.cosmo50PLK.6144g3HbwK1BH/h148.cosmo50PLK.6144g3HbwK1BH.004096/ahf_200/h148.cosmo50PLK.6144g3HbwK1BH.004096'
    
    tfile_base_2 = 'h229.cosmo50PLK.3072gst5HbwK1BH'
    tfile_2 = prefix + 'h229.cosmo50PLK.3072g/h229.cosmo50PLK.3072gst5HbwK1BH/snapshots_200bkgdens/h229.cosmo50PLK.3072gst5HbwK1BH.004096' #
    tfile_2 = prefix + 'h229.cosmo50PLK.3072g/h229.cosmo50PLK.3072gst5HbwK1BH/h229.cosmo50PLK.3072gst5HbwK1BH.004096/ahf_200/h229.cosmo50PLK.3072gst5HbwK1BH.004096' 
    #tfile_2 = prefix + 'h229.cosmo50PLK.3072g/h229.cosmo50PLK.3072gst5HbwK1BH/h229.cosmo50PLK.3072gst5HbwK1BH.004096/h229.cosmo50PLK.3072gst5HbwK1BH.004096'    
    
    tfile_base_3 = 'h242.cosmo50PLK.3072gst5HbwK1BH'
    tfile_3 = prefix + 'h242.cosmo50PLK.3072g/h242.cosmo50PLK.3072gst5HbwK1BH/snapshots_200bkgdens/h242.cosmo50PLK.3072gst5HbwK1BH.004096' #
    tfile_3 = prefix + 'h242.cosmo50PLK.3072g/h242.cosmo50PLK.3072gst5HbwK1BH/h242.cosmo50PLK.3072gst5HbwK1BH.004096/ahf_200/h242.cosmo50PLK.3072gst5HbwK1BH.004096'
    #tfile_3 = prefix + 'h242.cosmo50PLK.3072g/h242.cosmo50PLK.3072gst5HbwK1BH/h242.cosmo50PLK.3072gst5HbwK1BH.004096/h242.cosmo50PLK.3072gst5HbwK1BH.004096'

    tfile_base_4 = 'h329.cosmo50PLK.3072gst5HbwK1BH'
    tfile_4 = prefix + 'h329.cosmo50PLK.3072g/h329.cosmo50PLK.3072gst5HbwK1BH/snapshots_200bkgdens/h329.cosmo50PLK.3072gst5HbwK1BH.004096' #
    tfile_4 = prefix + 'h329.cosmo50PLK.3072g/h329.cosmo50PLK.3072gst5HbwK1BH/h329.cosmo50PLK.3072gst5HbwK1BH.004096/ahf_200/h329.cosmo50PLK.3072gst5HbwK1BH.004096'
    #tfile_4 = prefix + 'h329.cosmo50PLK.3072g/h329.cosmo50PLK.3072gst5HbwK1BH/h329.cosmo50PLK.3072gst5HbwK1BH.004096/h329.cosmo50PLK.3072gst5HbwK1BH.004096'

    tfile_base_4hr = 'h329.cosmo50PLK.6144g5HbwK1BH'
    tfile_4hr = prefix + 'h329.cosmo50PLK.6144g/h329.cosmo50PLK.6144g5HbwK1BH/h329.cosmo50PLK.6144g5HbwK1BH.004096/ahf_200/h329.cosmo50PLK.6144g5HbwK1BH.004096' #  
    
    outfile_base = prefix_outfile + 'marvelJL'
    tfiles = [tfile_1,tfile_2,tfile_3,tfile_4]
    tfile_base = [tfile_base_1,tfile_base_2,tfile_base_3,tfile_base_4]
    tfiles = [tfile_cm, tfile_e, tfile_r, tfile_s, tfile_1, tfile_2, tfile_3, tfile_4]
    tfile_base = [tfile_base_cm, tfile_base_e, tfile_base_r, tfile_base_s, tfile_base_1, tfile_base_2, tfile_base_3, tfile_base_4]

    tfiles = [tfile_cm, tfile_e, tfile_r, tfile_s, tfile_1hr, tfile_2, tfile_3, tfile_4hr]
    tfile_base = [tfile_base_cm, tfile_base_e, tfile_base_r, tfile_base_s, tfile_base_1hr, tfile_base_2, tfile_base_3, tfile_base_4hr]    
    
    cumulativeSFHenviron(tfiles,outfile_base,tfile_base)
    #cumulativeSFHenviron([tfile_cm, tfile_e, tfile_r, tfile_s, tfile_1,tfile_2,tfile_3,tfile_4],outfile_base,[tfile_base_cm, tfile_base_e, tfile_base_r, tfile_base_s, tfile_base_1,tfile_base_2,tfile_base_3,tfile_base_4])

        

