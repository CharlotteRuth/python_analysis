# Charlotte Christensen
# 12/6/21
# Plot the SMHM relation, mimicing abundance matching
# To do this, Tangos is used to link halos from dark matter to SPH sims
# The dark matter masses used are the 1) virial masses from the SPH,
#     2) the virial masses from the corresponding halos in the DMO sim
#     3) the sorted virial masses from the DMO sim
#
# Before running, use the tangos to do the following:
# tangos add snapshots_200crit_$id
# tangos import-properties Mvir Rvir Satellite --for snapshots_200crit_$id
# tangos write contamination_fraction --for snapshots_200crit_$id
# tangos write star_mass_profile --with-prerequisites --include-only="contamination_fraction<0.01" --include-only="NDM()>1000" --for snapshots_200crit_$id
# tangos add snapshots_200crit_$iddm
# tangos import-properties Mvir Rvir Satellite --for snapshots_200crit_$iddm
# tangos crosslink snapshots_200crit_$id snapshots_200crit_$iddm

import tangos
import pylab as p
import pandas as pd
import pynbody
import socket, os
import numpy as np
import matplotlib.pyplot as plt

def match_arr(x, y):
    # Match the halos from Ferah with the tangos database
    # Find indicies of elements in y that match x

    index = np.argsort(x)
    sorted_x = x[index] # sort the x array
    sorted_index = np.searchsorted(sorted_x, y)
    yindex = np.take(index, sorted_index, mode="clip")
    mask = x[yindex] != y
    result_y = np.ma.array(yindex, mask=mask)
    # use like y[~result_y.mask]. This will return the elements in y that match x

    xprime = y[~result_y.mask]
    yprime = x
    index = np.argsort(xprime)
    sorted_x = xprime[index] # sort the x array
    sorted_index = np.searchsorted(sorted_x, yprime)
    yindex = np.take(index, sorted_index, mode="clip")
    mask = xprime[yindex] != yprime
    result_x = np.ma.array(yindex, mask=mask)
    # use like x[~result_x.mask]. This will return the elements in x that match y

    # Will x[~result_x.mask] always equal y[~result_y.mask]?
    if np.sum(x[~result_x.mask] !=  y[~result_y.mask]):
        print("Warning: matched arrays have different orders")
    
    return result_x, result_y
    
#---------------------------------------------------------------

if (socket.gethostname() == "ozma.grinnell.edu"):
    dataprefix = '/home/christensen/Code/Datafiles/' 
else:
    dataprefix = '/home/christenc/Code/Datafiles/'

if (socket.gethostname() == "quirm.math.grinnell.edu"):
    prefix = '/home/christenc/Data/Sims/'
    prefix_outfile = '/home/christenc/Figures/marvel/'
else:
    prefix = '/home/christensen/Storage2/UW/MolecH/Cosmo/'
    prefix_outfile = '/home/christensen/Plots/marvel/'


hubble = 0.6776942783267969

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

# Initiate Plots
plt.close('all')
fig1 = plt.figure(figsize=(plt_width*2,plt_width*aspect_ratio*2))
ax1 = fig1.add_subplot()
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_ylabel(r'M$_*$/M$_\odot$')
ax1.set_xlabel(r'M$_{\mathrm{vir}}$/M$_\odot$')
ax1.axis([1e8,3e12,1e3,4e11])
ax1.fill_between(read_abunmatch['M200'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_low'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_high'][read_abunmatch['M200']> 1e9],color = 'grey', alpha = 0.5, zorder=0)

fig4 = plt.figure(figsize=(plt_width*2,plt_width*2))
ax4 = fig4.add_subplot()
ax4.set_xscale('log')
ax4.set_yscale('log')
ax4.set_ylabel(r'M$_*$/M$_\odot$')
ax4.set_xlabel(r'M$_{\mathrm{vir}}$/M$_\odot$')
ax4.axis([1e8,5e11,1e4,5e9])
ax4.fill_between(read_abunmatch['M200'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_low'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_high'][read_abunmatch['M200']> 1e9],color = 'grey', alpha = 0.5, zorder=0)

fig2 = plt.figure(figsize=(plt_width*2,plt_width*aspect_ratio*2))
ax2 = fig2.add_subplot()
ax2.plot([4e7,4e12],[1,1],color = 'k')
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_ylabel(r'M$_{vir,~measure}$/M$_{vir}$')
ax2.set_xlabel(r'M$_{\mathrm{vir}}$/M$_\odot$')

fig3 = plt.figure(figsize=(plt_width*2,plt_width*aspect_ratio*2))
ax3 = fig3.add_subplot()
ax3.plot([1e3,5e11],[1,1],color = 'k')
ax3.set_xscale('log')
ax3.set_yscale('log')
ax3.set_ylabel(r'M$_{vir,~measure}$/M$_{vir}$')
ax3.set_xlabel(r'M$_{*}$/M$_\odot$')

fig5 = plt.figure(figsize=(plt_width*2,plt_width*2))
ax5 = fig5.add_subplot()
ax5.set_xscale('log')
ax5.set_yscale('log')
ax5.set_ylabel(r'M$_*$/M$_\odot$')
ax5.set_xlabel(r'M$_{\mathrm{vir}}$/M$_\odot$')
ax5.axis([1e8,5e11,1e4,5e9])
ax5.fill_between(read_abunmatch['M200'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_low'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_high'][read_abunmatch['M200']> 1e9],color = 'grey', alpha = 0.5, zorder=0)

fig6 = plt.figure(figsize=(plt_width*2,plt_width*2))
ax6 = fig6.add_subplot()
ax6.set_xscale('log')
ax6.set_yscale('log')
ax6.set_ylabel(r'M$_*$/M$_\odot$')
ax6.set_xlabel(r'M$_{\mathrm{vir}}$/M$_\odot$')
ax6.axis([1e8,5e11,1e4,5e9])
ax6.fill_between(read_abunmatch['M200'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_low'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_high'][read_abunmatch['M200']> 1e9],color = 'grey', alpha = 0.5, zorder=0)


tangos.core.init_db("/home/christenc/Storage/tangos_db/JL_r200.db")
outfile_base = prefix_outfile + 'JL'
#export TANGOS_DB_CONNECTION="/home/christenc/Storage/tangos_db/JL_r200.db"
#os.environ["TANGOS_DB_CONNECTION"] = "/home/christenc/Storage/tangos_db/JL_r200.db"
sim = "h148"
sims = ["h148","h229","h242","h329"]

for sim in sims:
    timestep = tangos.get_timestep("snapshots_200crit_" + sim + "/%4096")
    timestep_dm = tangos.get_timestep("snapshots_200crit_" + sim + "dm/%4096")
    
    # Only if first run tangos write star_mass_profile --with-prerequisites --include-only="contamination_fraction<0.01" --include-only="NDM()>1000" --for snapshots_200crit_$id
    # Mvir, Mstar, sat, grp = timestep.calculate_all("Mvir","star_mass_profile[-1]","Satellite","Grp")
    # Mvir_match, Mstar_match, sat_match, grp_match, Mvir_match_dm = timestep.calculate_all("Mvir","star_mass_profile[-1]","Satellite","Grp","match('snapshots_200crit_" + sim + "dm').Mvir")
    # otherwise

    # SPH: virial mass, stellar mass, sat(y/n), grp
    Mvir, Mstar, sat, grp = timestep.calculate_all("Mvir","StarMass","Satellite","Grp") # add number of dark matter particles

    #halo = tangos.get_halo("snapshots_200crit_" + sim + "/%4096/halo_1")
    #halo_dm = tangos.get_halo("snapshots_200crit_" + sim + "dm/%4096/halo_1")
    #halo.keys()    

    #Mvir_dm, fMhires = timestep_dm.calculate_all("Mvir","fMhires")
    #Mvir_dm, fMhires, grpdm = timestep_dm.calculate_all("Mvir","Contam","Grp") # for some reason, this ignores all halos for which Contam is zero
    
    # DMO: virial mass, grp
    Mvir_dm,grpdm = timestep_dm.calculate_all("Mvir","Grp")
    

    # This section of the code is in case the tangos linking isn't working.
    tangos_link = 0
    if not tangos_link:
        filename = timestep.filename
        filename = '/'.join(filename.split('/')[3:])
        filename_dm = timestep_dm.filename
        filename_dm = '/'.join(filename_dm.split('/')[3:])    
        basename = '.'.join((filename.split('/')[-1]).split('.')[0:2])
        longbase = '.'.join((filename.split('/')[-1]).split('.')[0:3])
    
        simname = pynbody.load('/home/christenc/Data/Sims/' + basename + '.3072g/' + longbase + '/' + filename)
        halos = simname.halos()
        sim_dm = pynbody.load('/home/christenc/Data/Sims/' + basename + '.3072g/' + longbase + '/' + filename_dm)
        halos_dm = sim_dm.halos()
        sim_dm.dark['iord'] = sim_dm.dark['iord'] - (sim_dm.dark['iord'][0] - simname.dark['iord'][0])

        try:
            sim_dm['amiga.grp']
        except:
            halos_dm.make_grp()
            grpoutfile = sim_dm.filename + '.amiga.grp'
            formatarr = "%" + str(len(str(halos_dm._nhalos)) + 1) + "d"
            grp = np.append(np.array([len(sim_dm['grp'])]),sim_dm['grp']) #problems here
            np.savetxt(grpoutfile, grp, fmt = formatarr)    

        # since the grp file contains halos numbes larger than tangos tracks, force eliminate them
        # note, this happens because Alyon's code to create amiga files eliminate halos from low res region
        sim_dm['amiga.grp'][sim_dm['amiga.grp'] > max(grpdm)] = 0
        nhalos_dm =  max(sim_dm['amiga.grp'])
        all_halos_dm = np.histogram(sim_dm.dark['amiga.grp'],bins = nhalos_dm,range = (0,nhalos_dm))
        grp_match_dm = grp*0
        Mvir_match_dm = grp*0
        it = 0
        for grp_halo in grp[Mstar> 1e3]:
            halo = halos[grp_halo]
            ind = np.searchsorted(sim_dm.dark['iord'],halo.dark['iord'])
            # ind is an array such that sim_dm.dark['iord'][ind] == halo.dark['iord']
            matchhist = np.histogram((sim_dm.dark['amiga.grp'])[ind],bins = nhalos_dm,range = (0,nhalos_dm))
            min_dm_part = (Mvir[grp == grp_halo]/1e5)/10 #Match must have at least 10% the dark matter particle mass
            matchhist[0][matchhist[0]<min_dm_part] = 0
            ratio = matchhist[0]/all_halos_dm[0]
            arg_match = argmax(ratio[(ratio<1.5) & (ratio > 0)])
            grp_match_dm[(np.where(Mstar>1e3))[0][it]] = (matchhist[1][:-1])[(ratio<1.5) & (ratio > 0)][arg_match]
            print((matchhist[0])[(ratio<1.5) & (ratio > 0)],(matchhist[1][:-1])[(ratio<1.5) & (ratio > 0)],ratio[(ratio<1.5) & (ratio > 0)])
            print(grp_halo,grp_match_dm[it],(matchhist[0][:])[(ratio<1.5) & (ratio > 0)][arg_match])
            if (np.sum(ratio>1)>1):
                print("Problem: more matching particles than there are particles/ratio>1")
                print((matchhist[0])[ratio>1.5],(matchhist[1][:-1])[ratio>1.5])
            it = it + 1
        # Set the virial masses of the matched halo, leaving as zero any unmatched halos
        Mvir_match_dm[grp_match_dm> 0] = Mvir_dm[np.searchsorted(grpdm,grp_match_dm[grp_match_dm> 0])]
        Mvir_match = Mvir
        Mstar_match = Mstar
        sat_match = sat
        grp_match = grp
    else:    
    # matching SPH:  virial mass, stellar mass, sat(y/n), grp; matching DMO: virial mass
        Mvir_match, Mstar_match, sat_match, grp_match, Mvir_match_dm, grp_match_dm  = timestep.calculate_all("Mvir","StarMass","Satellite","Grp","match('snapshots_200crit_" + sim + "dm').Mvir","match('snapshots_200crit_" + sim + "dm').Grp")

        #Mvir_match_dm2, grp_match_dm2, Mvir_match2, Mstar_match2, sat_match2, grp_match2  = timestep_dm.calculate_all("Mvir","Grp","match('snapshots_200crit_" + sim + "').Mvir","match('snapshots_200crit_" + sim + "').StarMass","match('snapshots_200crit_" + sim + "').Satellite","match('snapshots_200crit_" + sim + "').Grp")
        #Mvir_match_dm, grp_match_dm, Mvir_match, Mstar_match, sat_match, grp_match = timestep_dm.calculate_all("Mvir","ID","match('snapshots_200crit_" + sim + "').Mvir","match('snapshots_200crit_" + sim + "').StarMass","match('snapshots_200crit_" + sim + "').Satellite","match('snapshots_200crit_" + sim + "').Grp")

        
    fMhires = np.empty(len(Mvir_dm)) + 1.0
    it = 0
    for grpdm_halo in grpdm:
        #halo = tangos.get_halo("snapshots_200crit_" + sim + "dm/%4096/halo_" + str(grpdm_halo))
        halo = timestep_dm.halos[it]
        fMhires[it] = halo['Contam']
        if it > 1291: break
        #Mvir_dm[i] < 1e5: break
        it = it + 1

    Mvir_dm_res = Mvir_dm[fMhires < 0.1]

    # Read in data of peak masses
    dm_peak_data = pd.read_csv(dataprefix + sim + '.cosmo50PLK.3072gst_Mpeak.csv')
    dm_peak_data = dm_peak_data.set_index('halogrp_z0')   

    # Match Ferah's information to my tangos DB
    match_ind_fdm, match_ind_tangos = match_arr(np.array(fdmdata[fdmdata['simname'] == sim]['halogrp_z0']),np.array(grp))
    Mvir = Mvir[~match_ind_tangos.mask]
    Mstar = Mstar[~match_ind_tangos.mask]
    sat = sat[~match_ind_tangos.mask]
    grp = grp[~match_ind_tangos.mask]
    fdm_sim = fdmdata[fdmdata['simname'] == sim][~match_ind_fdm.mask]

    # Select only those most massive of the dm halos
    sort_arg = np.flip(np.argsort(Mstar_match))
    Mvir_match_dm = Mvir_match_dm[sort_arg]
    grp_match_dm = grp_match_dm[sort_arg]
    Mvir_match = Mvir_match[sort_arg]
    Mstar_match = Mstar_match[sort_arg]
    sat_match = sat_match[sort_arg]
    grp_match = grp_match[sort_arg]

    Mvir_match_dm, unique_ind = np.unique(Mvir_match_dm, return_index=True)
    grp_match_dm = grp_match_dm[unique_ind]
    Mvir_match = Mvir_match[unique_ind]
    Mstar_match = Mstar_match[unique_ind]
    sat_match = sat_match[unique_ind]
    grp_match = grp_match[unique_ind]
    
    # Select those halos from the DMO simulation that match Ferah's halos
    match_ind_DMO, temp = match_arr(grp_match,grp)
    Mvir_match_dm = Mvir_match_dm[~match_ind_DMO.mask]
    Mvir_match = Mvir_match[~match_ind_DMO.mask]
    Mstar_match = Mstar_match[~match_ind_DMO.mask]
    sat_match = sat_match[~match_ind_DMO.mask]
    grp_match_dm = grp_match_dm[~match_ind_DMO.mask]
    
    # Find the peak halo masses for the halos in the DMO simulations that match Ferah's halos
    # Find the intersection between the halos in the tangos directory and those for which I have peak masses
    grp_match_intersect, comm1, comm2 = np.intersect1d(grp_match_dm,dm_peak_data.index,return_indices=True)
    Mvir_match_dm_peak = dm_peak_data.loc[grp_match_intersect,'Mpeak']
    # Mpeak comes from AHF, rather than amiga.stat, divide by hubble parameter

    # Sort the DM simulation by mass so it can be paired to sorted stellar masses
    Mvir_dm_sort = (np.flip(np.sort(Mvir_dm_res)))[0:min([len(Mstar),len(Mvir_dm_res)])]
    Mstar_sort = (np.flip(np.sort(Mstar)))[0:min([len(Mstar),len(Mvir_dm_res)])]
    Mvir_sort = (np.flip(np.sort(Mvir)))[0:min([len(Mvir),len(Mvir_dm_res)])]
    sort_arg = np.argsort(Mstar)
    sat_sort = (np.flip(sat[sort_arg]))[0:min([len(Mstar),len(Mvir_dm_res)])]

    # SMHM
    direct = ax1.scatter(Mvir[sat == 0], Mstar[sat == 0],marker = 'o',color = 'navy', zorder = 1, s = markersize)
    ax1.scatter(Mvir[sat != 0], Mstar[sat != 0],marker = 'D',color = 'navy', zorder = 1, s = markersize)
    peak = ax1.scatter(fdm_sim[fdm_sim['type']!='Satellite']['Mpeak'],fdm_sim[fdm_sim['type']!='Satellite']['Mstar_z0'],marker = 'o',facecolor = 'none',edgecolor = 'blueviolet',linewidths = 2*edgewidth, zorder = 2, s = markersize)
    ax1.scatter(fdm_sim[fdm_sim['type']=='Satellite']['Mpeak'],fdm_sim[fdm_sim['type']=='Satellite']['Mstar_z0'],marker = 'D',facecolor = 'none',edgecolor = 'blueviolet',linewidths = 2*edgewidth, zorder = 2, s = markersize)
    dmo = ax1.scatter(Mvir_match_dm_peak[sat_match[comm1] == 0], Mstar_match[comm1][[sat_match[comm1] == 0]],marker = 'o',facecolor = 'purple',edgecolor = 'purple',linewidths = edgewidth, zorder = 3, s = markersize, alpha = 0.5)
    ax1.scatter(Mvir_match_dm_peak[sat_match[comm1] != 0], Mstar_match[comm1][[sat_match[comm1] != 0]],marker = 'D',facecolor = 'purple',edgecolor = 'purple',linewidths = edgewidth, zorder = 3, s = markersize, alpha = 0.5)
    #dmo = ax1.scatter(Mvir_match_dm[sat_match == 0], Mstar_match[sat_match == 0],marker = 'o',facecolor = 'blueviolet',edgecolor = 'purple',linewidths = edgewidth, zorder = 3, s = markersize, alpha = 0.5)
    #ax1.scatter(Mvir_match_dm[sat_match != 0], Mstar_match[sat_match != 0],marker = 'D',facecolor = 'blueviolet',edgecolor = 'purple',linewidths = edgewidth, zorder = 3, s = markersize, alpha = 0.5)
    sort = ax1.scatter(Mvir_dm_sort[sat_sort == 0],Mstar_sort[sat_sort == 0],marker = 'o',facecolor = 'royalblue',edgecolor = 'royalblue',linewidths = edgewidth, zorder = 3, s = markersize, alpha = 0.5)
    ax1.scatter(Mvir_dm_sort[sat_sort != 0],Mstar_sort[sat_sort != 0],marker = 'D',facecolor = 'royalblue',edgecolor = 'royalblue',linewidths = edgewidth, zorder = 3, s = markersize, alpha = 0.5)
    legend = ax1.legend([direct,peak,dmo,sort],['M$_{vir}$','Peak M$_{vir}$','DMO M$_{vir}$','Sorted DMO M$_{vir}$'], loc=2)

     # SMHM
    direct = ax4.scatter(Mvir[sat == 0], Mstar[sat == 0],marker = 'o',color = 'navy', zorder = 1, s = markersize)
    ax4.scatter(Mvir[sat != 0], Mstar[sat != 0],marker = 'D',color = 'navy', zorder = 1, s = markersize)
    peak = ax4.scatter(fdm_sim[fdm_sim['type']!='Satellite']['Mpeak'],fdm_sim[fdm_sim['type']!='Satellite']['Mstar_z0'],marker = 'o',facecolor = 'none',edgecolor = 'blueviolet',linewidths = 2*edgewidth, zorder = 2, s = markersize)
    ax4.scatter(fdm_sim[fdm_sim['type']=='Satellite']['Mpeak'],fdm_sim[fdm_sim['type']=='Satellite']['Mstar_z0'],marker = 'D',facecolor = 'none',edgecolor = 'blueviolet',linewidths = 2*edgewidth, zorder = 2, s = markersize)
    legend = ax4.legend([direct,peak,dmo,sort],['M$_{vir}$','Peak M$_{vir}$','DMO M$_{vir}$','Sorted DMO M$_{vir}$'], loc=2)
    dmo = ax4.scatter(Mvir_match_dm_peak[sat_match[comm1] == 0], Mstar_match[comm1][sat_match[comm1] == 0],marker = 'o',facecolor = 'purple',edgecolor = 'purple',linewidths = edgewidth, zorder = 3, s = markersize, alpha = 0.5)
    ax4.scatter(Mvir_match_dm_peak[sat_match[comm1] != 0], Mstar_match[comm1][sat_match[comm1] != 0],marker = 'D',facecolor = 'purple',edgecolor = 'purple',linewidths = edgewidth, zorder = 3, s = markersize, alpha = 0.5)    
    #dmo = ax4.scatter(Mvir_match_dm[sat_match == 0], Mstar_match[sat_match == 0],marker = 'o',facecolor = 'purple',edgecolor = 'purple',linewidths = edgewidth, zorder = 3, s = markersize, alpha = 0.5)
    #ax4.scatter(Mvir_match_dm[sat_match != 0], Mstar_match[sat_match != 0],marker = 'D',facecolor = 'purple',edgecolor = 'purple',linewidths = edgewidth, zorder = 3, s = markersize, alpha = 0.5)   
    sort = ax4.scatter(Mvir_dm_sort[sat_sort == 0],Mstar_sort[sat_sort == 0],marker = 'o',facecolor = 'royalblue',edgecolor = 'royalblue',linewidths = edgewidth, zorder = 3, s = markersize, alpha = 0.5)
    ax4.scatter(Mvir_dm_sort[sat_sort != 0],Mstar_sort[sat_sort != 0],marker = 'D',facecolor = 'royalblue',edgecolor = 'royalblue',linewidths = edgewidth, zorder = 3, s = markersize, alpha = 0.5)   

    # Ratio of DM masses vs DM mass
    ax2.scatter(fdm_sim[(fdm_sim['simname'] == sim) & (fdm_sim['type']!='Satellite')]['Mhalo_z0'],fdm_sim[(fdm_sim['simname'] == sim) & (fdm_sim['type']!='Satellite')]['Mpeak']/fdm_sim[(fdm_sim['simname'] == sim) & (fdm_sim['type']!='Satellite')]['Mhalo_z0'],marker = 'o',facecolor = 'none',edgecolor = 'blueviolet',linewidths = edgewidth*2, zorder = 2, s = markersize)
    ax2.scatter(fdm_sim[(fdm_sim['simname'] == sim) & (fdm_sim['type']=='Satellite')]['Mhalo_z0'],fdm_sim[(fdm_sim['simname'] == sim) & (fdm_sim['type']=='Satellite')]['Mpeak']/fdm_sim[(fdm_sim['simname'] == sim) & (fdm_sim['type']=='Satellite')]['Mhalo_z0'],marker = 'D',facecolor = 'none',edgecolor = 'blueviolet',linewidths = edgewidth*2, zorder = 2, s = markersize)    
    ax2.scatter(Mvir_match_dm_peak[sat_match[comm1] == 0],Mvir_match_dm_peak[sat_match[comm1] == 0]/Mvir_match[comm1][sat_match[comm1] == 0],marker = 'o',facecolor = 'purple',edgecolor = 'purple',linewidths = edgewidth, zorder = 3, s = markersize, alpha = 0.5)
    ax2.scatter(Mvir_match_dm_peak[sat_match[comm1] != 0],Mvir_match_dm_peak[sat_match[comm1] != 0]/Mvir_match[comm1][sat_match[comm1] != 0],marker = 'D',facecolor = 'purple',edgecolor = 'purple',linewidths = edgewidth, zorder = 3, s = markersize, alpha = 0.5)
    #ax2.scatter(Mvir_match[sat_match == 0],Mvir_match_dm[sat_match == 0]/Mvir_match[sat_match == 0],marker = 'o',facecolor = 'purple',edgecolor = 'purple',linewidths = edgewidth, zorder = 3, s = markersize, alpha = 0.5)
    #ax2.scatter(Mvir_match[sat_match != 0],Mvir_match_dm[sat_match != 0]/Mvir_match[sat_match != 0],marker = 'D',facecolor = 'purple',edgecolor = 'purple',linewidths = edgewidth, zorder = 3, s = markersize, alpha = 0.5)   
    ax2.scatter(Mvir_sort[sat_sort == 0],Mvir_dm_sort[sat_sort == 0]/Mvir_sort[sat_sort == 0],marker = 'o',facecolor = 'royalblue',edgecolor = 'royalblue',linewidths = edgewidth, zorder = 3, s = markersize, alpha = 0.5)
    ax2.scatter(Mvir_sort[sat_sort != 0],Mvir_dm_sort[sat_sort != 0]/Mvir_sort[sat_sort != 0],marker = 'D',facecolor = 'royalblue',edgecolor = 'royalblue',linewidths = edgewidth, zorder = 3, s = markersize, alpha = 0.5)    
    legend = ax2.legend([peak,dmo,sort],['Peak M$_{vir}$','DMO M$_{vir}$','Sorted DMO M$_{vir}$'])
    #ax.axis([3e8,2e10,0.1,40])

    # Ratio of DM masses vs stellar mass
    ax3.scatter(fdm_sim[(fdm_sim['simname'] == sim) & (fdm_sim['type']!='Satellite')]['Mstar_z0'],fdm_sim[(fdm_sim['simname'] == sim) & (fdm_sim['type']!='Satellite')]['Mpeak']/fdm_sim[(fdm_sim['simname'] == sim) & (fdm_sim['type']!='Satellite')]['Mhalo_z0'],marker = 'o',facecolor = 'none',edgecolor = 'blueviolet',linewidths = edgewidth*2, zorder = 2, s = markersize)
    ax3.scatter(fdm_sim[(fdm_sim['simname'] == sim) & (fdm_sim['type']=='Satellite')]['Mstar_z0'],fdm_sim[(fdm_sim['simname'] == sim) & (fdm_sim['type']=='Satellite')]['Mpeak']/fdm_sim[(fdm_sim['simname'] == sim) & (fdm_sim['type']=='Satellite')]['Mhalo_z0'],marker = 'D',facecolor = 'none',edgecolor = 'blueviolet',linewidths = edgewidth*2, zorder = 2, s = markersize)    
    ax3.scatter(Mstar_match[comm1][sat_match[comm1] == 0],Mvir_match_dm_peak[sat_match[comm1] == 0]/Mvir_match[comm1][sat_match[comm1] == 0],marker = 'o',facecolor = 'purple',edgecolor = 'purple',linewidths = edgewidth, zorder = 3, s = markersize, alpha = 0.5)
    ax3.scatter(Mstar_match[comm1][sat_match[comm1] != 0],Mvir_match_dm_peak[sat_match[comm1] != 0]/Mvir_match[comm1][sat_match[comm1] != 0],marker = 'D',facecolor = 'purple',edgecolor = 'purple',linewidths = edgewidth, zorder = 3, s = markersize, alpha = 0.5)
    #ax3.scatter(Mstar_match[sat_match == 0],Mvir_match_dm[sat_match == 0]/Mvir_match[sat_match == 0],marker = 'o',facecolor = 'purple',edgecolor = 'purple',linewidths = edgewidth, zorder = 3, s = markersize, alpha = 0.5)
    #ax3.scatter(Mstar_match[sat_match != 0],Mvir_match_dm[sat_match != 0]/Mvir_match[sat_match != 0],marker = 'D',facecolor = 'purple',edgecolor = 'purple',linewidths = edgewidth, zorder = 3, s = markersize, alpha = 0.5)    
    ax3.scatter(Mstar_sort[sat_sort == 0],Mvir_dm_sort[sat_sort == 0]/Mvir_sort[sat_sort == 0],marker = 'o',facecolor = 'royalblue',edgecolor = 'royalblue',linewidths = edgewidth, zorder = 3, s = markersize, alpha = 0.5)
    ax3.scatter(Mstar_sort[sat_sort != 0],Mvir_dm_sort[sat_sort != 0]/Mvir_sort[sat_sort != 0],marker = 'D',facecolor = 'royalblue',edgecolor = 'royalblue',linewidths = edgewidth, zorder = 3, s = markersize, alpha = 0.5)
    legend = ax3.legend([peak,dmo,sort],['Peak M$_{vir}$','DMO M$_{vir}$','Sorted DMO M$_{vir}$'])

    # SMHM for sorted values
    sort = ax5.scatter(Mvir_dm_sort[sat_sort == 0],Mstar_sort[sat_sort == 0],marker = 'o',facecolor = 'royalblue',edgecolor = 'royalblue',linewidths = edgewidth, zorder = 3, s = markersize, alpha = 0.5)
    ax5.scatter(Mvir_dm_sort[sat_sort != 0],Mstar_sort[sat_sort != 0],marker = 'D',facecolor = 'royalblue',edgecolor = 'royalblue',linewidths = edgewidth, zorder = 3, s = markersize, alpha = 0.5)
    ax5.scatter(fdm_sim[fdm_sim['type']=='Satellite']['Mpeak'],fdm_sim[fdm_sim['type']=='Satellite']['Mstar_z0'],marker = 'D',facecolor = 'none',edgecolor = 'blueviolet',linewidths = 2*edgewidth, zorder = 2, s = markersize)
    legend = ax5.legend([direct,peak,dmo,sort],['M$_{vir}$','Peak M$_{vir}$','DMO M$_{vir}$','Sorted DMO M$_{vir}$'], loc=2)

    #SMHM for DMO matched
    dmo = ax6.scatter(Mvir_match_dm_peak[sat_match[comm1] == 0], Mstar_match[comm1][sat_match[comm1] == 0],marker = 'o',facecolor = 'purple',edgecolor = 'purple',linewidths = edgewidth, zorder = 3, s = markersize, alpha = 0.5)
    ax6.scatter(Mvir_match_dm_peak[sat_match[comm1] != 0], Mstar_match[comm1][sat_match[comm1] != 0],marker = 'D',facecolor = 'purple',edgecolor = 'purple',linewidths = edgewidth, zorder = 3, s = markersize, alpha = 0.5)
    ax6.scatter(fdm_sim[fdm_sim['type']=='Satellite']['Mpeak'],fdm_sim[fdm_sim['type']=='Satellite']['Mstar_z0'],marker = 'D',facecolor = 'none',edgecolor = 'blueviolet',linewidths = 2*edgewidth, zorder = 2, s = markersize)
    legend = ax6.legend([direct,peak,dmo,sort],['M$_{vir}$','Peak M$_{vir}$','DMO M$_{vir}$','Sorted DMO M$_{vir}$'], loc=2)
    
fig1.savefig(outfile_base + '_SMHM_multiMvir.png',dpi = dpi)
fig4.savefig(outfile_base + '_SMHM_multiMvir_read.png',dpi = dpi)
fig2.savefig(outfile_base + '_multiMvir_Mvir.png',dpi = dpi)
fig3.savefig(outfile_base + '_multiMvir_Mstar.png',dpi = dpi)

outfile_base = prefix_outfile + 'marvel'

fig1 = plt.figure(figsize=(plt_width*2,plt_width*aspect_ratio*2))
fig1.clear()
ax1 = fig1.add_subplot()
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_ylabel(r'M$_*$/M$_\odot$')
ax1.set_xlabel(r'M$_{\mathrm{vir}}$/M$_\odot$')
ax1.axis([1e8,3e12,1e3,4e11])
ax1.fill_between(read_abunmatch['M200'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_low'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_high'][read_abunmatch['M200']> 1e9],color = 'grey', alpha = 0.5, zorder=0)

fig4 = plt.figure(figsize=(plt_width*2,plt_width*2))
fig4.clear()
ax4 = fig4.add_subplot()
ax4.set_xscale('log')
ax4.set_yscale('log')
ax4.set_ylabel(r'M$_*$/M$_\odot$')
ax4.set_xlabel(r'M$_{\mathrm{vir}}$/M$_\odot$')
ax4.axis([1e8,5e11,1e4,5e9])
ax4.fill_between(read_abunmatch['M200'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_low'][read_abunmatch['M200']> 1e9],read_abunmatch['Mstar_high'][read_abunmatch['M200']> 1e9],color = 'grey', alpha = 0.5, zorder=0)

fig2 = plt.figure(figsize=(plt_width*2,plt_width*aspect_ratio*2))
fig2.clear()
ax2 = fig2.add_subplot()
ax2.plot([4e7,4e12],[1,1],color = 'k')
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_ylabel(r'M$_{vir,~measure}$/M$_{vir}$')
ax2.set_xlabel(r'M$_{\mathrm{vir}}$/M$_\odot$')

fig3 = plt.figure(figsize=(plt_width*2,plt_width*aspect_ratio*2))
fig3.clear()
ax3 = fig3.add_subplot()
ax3.plot([1e3,5e11],[1,1],color = 'k')
ax3.set_xscale('log')
ax3.set_yscale('log')
ax3.set_ylabel(r'M$_{vir,~measure}$/M$_{vir}$')
ax3.set_xlabel(r'M$_{*}$/M$_\odot$')
#export TANGOS_DB_CONNECTION="/home/christenc/Storage/tangos_db/Marvel_r200.db"
#os.environ["TANGOS_DB_CONNECTION"] = "/home/christenc/Storage/tangos_db/Marvel_r200.db"
tangos.core.init_db("/home/christenc/Storage/tangos_db/Marvel_r200.db")
sims = ["cptmarvel","storm","elektra","rogue"]

for sim in sims:
    timestep = tangos.get_timestep("snapshots_200crit_" + sim + "/%4096")
    # Only if first run tangos write star_mass_profile --with-prerequisites --include-only="contamination_fraction<0.01" --include-only="NDM()>1000" --for snapshots_200crit_$id
    # Mvir, Mstar, sat, grp = timestep.calculate_all("Mvir","star_mass_profile[-1]","Satellite","Grp")
    # Mvir_match, Mstar_match, sat_match, grp_match, Mvir_match_dm = timestep.calculate_all("Mvir","star_mass_profile[-1]","Satellite","Grp","match('snapshots_200crit_" + sim + "dm').Mvir")
    # otherwise

    # SPH: virial mass, stellar mass, sat(y/n), grp
    Mvir, Mstar, sat, grp = timestep.calculate_all("Mvir","StarMass","Satellite","Grp")

    # matching SPH:  virial mass, stellar mass, sat(y/n), grp; matching DMO: virial mass
    Mvir_match, Mstar_match, sat_match, grp_match, Mvir_match_dm, grp_match_dm  = timestep.calculate_all("Mvir","StarMass","Satellite","Grp","match('snapshots_200crit_" + sim + "dm').Mvir","match('snapshots_200crit_" + sim + "dm').Grp")
    
    timestep_dm = tangos.get_timestep("snapshots_200crit_" + sim + "dm/%4096")
    #halo = tangos.get_halo("snapshots_200crit_" + sim + "dm/%4096/halo_1")
    #halo.keys()    
    #Mvir_match_dm, grp_match_dm, Mvir_match, Mstar_match, sat_match, grp_match = timestep_dm.calculate_all("Mvir","ID","match('snapshots_200crit_" + sim + "').Mvir","match('snapshots_200crit_" + sim + "').StarMass","match('snapshots_200crit_" + sim + "').Satellite","match('snapshots_200crit_" + sim + "').Grp")
    #Mvir_dm, fMhires = timestep_dm.calculate_all("Mvir","fMhires")
    Mvir_dm, fMhires = timestep_dm.calculate_all("Mvir","Contam")
    Mvir_dm = Mvir_dm[fMhires < 0.01]   

    # Match Ferah's information to my tangos DB
    match_ind_fdm, match_ind_tangos = match_arr(np.array(fdmdata[fdmdata['simname'] == sim]['halogrp_z0']),np.array(grp))
    Mvir = Mvir[~match_ind_tangos.mask]
    Mstar = Mstar[~match_ind_tangos.mask]
    sat = sat[~match_ind_tangos.mask]
    grp = grp[~match_ind_tangos.mask]
    fdm_sim = fdmdata[fdmdata['simname'] == sim][~match_ind_fdm.mask]

    # Select only those most massive of the dm halos
    sort_arg = np.flip(np.argsort(Mstar_match))
    Mvir_match_dm = Mvir_match_dm[sort_arg]
    grp_match_dm = grp_match_dm[sort_arg]
    Mvir_match = Mvir_match[sort_arg]
    Mstar_match = Mstar_match[sort_arg]
    sat_match = sat_match[sort_arg]
    grp_match = grp_match[sort_arg]

    Mvir_match_dm, unique_ind = np.unique(Mvir_match_dm, return_index=True)
    grp_match_dm = grp_match_dm[unique_ind]
    Mvir_match = Mvir_match[unique_ind]
    Mstar_match = Mstar_match[unique_ind]
    sat_match = sat_match[unique_ind]
    grp_match = grp_match[unique_ind]
    
    # Select those halos from the DMO simulation that match Ferah's halos
    match_ind_DMO, temp = match_arr(grp_match,grp)
    Mvir_match_dm = Mvir_match_dm[~match_ind_DMO.mask]
    Mvir_match = Mvir_match[~match_ind_DMO.mask]
    Mstar_match = Mstar_match[~match_ind_DMO.mask]
    sat_match = sat_match[~match_ind_DMO.mask]
    grp_match_dm = grp_match_dm[~match_ind_DMO.mask]
    
    # Find the peak halo masses for the halos in the DMO simulations that match Ferah's halos
    """
    Mvir_match_dm_peak = np.empty(len(grp_match_dm))
    i = 0
    for halo in grp_match_dm:
        tangos_halo = tangos.get_halo("snapshots_200crit_" + sim + "dm/%.004096/halo_" + str(halo))
        if not tangos_halo is None:
            time_prog, ids_progs, mvir_prog = tangos_halo.calculate_for_progenitors("t()", "halo_number()", "Mvir")
            #time, ids, mvir = tangos_halo.calculate_for_progenitors("t()", "halo_number()", "Mvir")
            if len(mvir_prog) > 0:
                print(halo,time_prog[np.argmax(mvir_prog)],mvir_prog[0]/max(mvir_prog))
                Mvir_match_dm_peak[i] = max(mvir_prog)
        i = i + 1
    """
        
    
    # Find the intersection between the halos in the tangos directory and those for which I have peak masses
    #grp_match_intersect, comm1, comm2 = np.intersect1d(grp_match_dm,dm_peak_data.index,return_indices=True)
    #Mvir_match_dm_peak = dm_peak_data.loc[grp_match_intersect,'Mpeak']/hubble #Divide by hubble parameter because the Mpeak comes from AHF, rather than amiga.stat    


    
    # Sort the DM simulation by mass so it can be paired to sorted stellar masses
    Mvir_dm_sort = (np.flip(np.sort(Mvir_dm)))[0:min([len(Mstar),len(Mvir_dm)])]
    Mstar_sort = (np.flip(np.sort(Mstar)))[0:min([len(Mstar),len(Mvir_dm)])]
    Mvir_sort = (np.flip(np.sort(Mvir)))[0:min([len(Mvir),len(Mvir_dm)])]
    sort_arg = np.argsort(Mstar)
    sat_sort = (np.flip(sat[sort_arg]))[0:min([len(Mstar),len(Mvir_dm)])]

    # SMHM
    direct = ax1.scatter(Mvir[sat == 0], Mstar[sat == 0],marker = 'o',color = 'maroon', zorder = 1, s = markersize)
    ax1.scatter(Mvir[sat != 0], Mstar[sat != 0],marker = 'D',color = 'maroon', zorder = 1, s = markersize)
    peak = ax1.scatter(fdm_sim[fdm_sim['type']!='Satellite']['Mpeak'],fdm_sim[fdm_sim['type']!='Satellite']['Mstar_z0'],marker = 'o',facecolor = 'none',edgecolor = 'darkorange',linewidths = 2*edgewidth, zorder = 2, s = markersize)
    ax1.scatter(fdm_sim[fdm_sim['type']=='Satellite']['Mpeak'],fdm_sim[fdm_sim['type']=='Satellite']['Mstar_z0'],marker = 'D',facecolor = 'none',edgecolor = 'darkorange',linewidths = 2*edgewidth, zorder = 2, s = markersize)
    #dmo = ax1.scatter(Mvir_match_dm_peak[sat_match[comm1] == 0], Mstar_match[comm1][[sat_match[comm1] == 0]],marker = 'o',facecolor = 'blueviolet',edgecolor = 'blueviolet',linewidths = edgewidth, zorder = 3, s = markersize, alpha = 0.5)
    #ax1.scatter(Mvir_match_dm_peak[sat_match[comm1] != 0], Mstar_match[comm1][[sat_match[comm1] != 0]],marker = 'D',facecolor = 'blueviolet',edgecolor = 'blueviolet',linewidths = edgewidth, zorder = 3, s = markersize, alpha = 0.5)
    dmo = ax1.scatter(Mvir_match_dm[sat_match == 0], Mstar_match[sat_match == 0],marker = 'o',facecolor = 'coral',edgecolor = 'coral',linewidths = edgewidth, zorder = 3, s = markersize, alpha = 0.5)
    ax1.scatter(Mvir_match_dm[sat_match != 0], Mstar_match[sat_match != 0],marker = 'D',facecolor = 'coral',edgecolor = 'coral',linewidths = edgewidth, zorder = 3, s = markersize, alpha = 0.5)
    sort = ax1.scatter(Mvir_dm_sort[sat_sort == 0],Mstar_sort[sat_sort == 0],marker = 'o',facecolor = 'gold',edgecolor = 'gold',linewidths = edgewidth, zorder = 3, s = markersize, alpha = 0.5)
    ax1.scatter(Mvir_dm_sort[sat_sort != 0],Mstar_sort[sat_sort != 0],marker = 'D',facecolor = 'gold',edgecolor = 'gold',linewidths = edgewidth, zorder = 3, s = markersize, alpha = 0.5)
    legend = ax1.legend([direct,peak,dmo,sort],['M$_{vir}$','Peak M$_{vir}$','DMO M$_{vir}$','Sorted DMO M$_{vir}$'], loc=2)    

     # SMHM
    direct = ax4.scatter(Mvir[sat == 0], Mstar[sat == 0],marker = 'o',color = 'maroon', zorder = 1, s = markersize)
    ax4.scatter(Mvir[sat != 0], Mstar[sat != 0],marker = 'D',color = 'maroon', zorder = 1, s = markersize)
    peak = ax4.scatter(fdm_sim[fdm_sim['type']!='Satellite']['Mpeak'],fdm_sim[fdm_sim['type']!='Satellite']['Mstar_z0'],marker = 'o',facecolor = 'none',edgecolor = 'darkorange',linewidths = 2*edgewidth, zorder = 2, s = markersize)
    ax4.scatter(fdm_sim[fdm_sim['type']=='Satellite']['Mpeak'],fdm_sim[fdm_sim['type']=='Satellite']['Mstar_z0'],marker = 'D',facecolor = 'none',edgecolor = 'darkorange',linewidths = 2*edgewidth, zorder = 2, s = markersize)
    legend = ax4.legend([direct,peak,dmo,sort],['M$_{vir}$','Peak M$_{vir}$','DMO M$_{vir}$','Sorted DMO M$_{vir}$'], loc=2)
    #dmo = ax4.scatter(Mvir_match_dm_peak[sat_match[comm1] == 0], Mstar_match[comm1][sat_match[comm1] == 0],marker = 'o',facecolor = 'm',edgecolor = 'm',linewidths = edgewidth, zorder = 3, s = markersize, alpha = 0.5)
    #ax4.scatter(Mvir_match_dm_peak[sat_match[comm1] != 0], Mstar_match[comm1][sat_match[comm1] != 0],marker = 'D',facecolor = 'm',edgecolor = 'm',linewidths = edgewidth, zorder = 3, s = markersize, alpha = 0.5)    
    dmo = ax4.scatter(Mvir_match_dm[sat_match == 0], Mstar_match[sat_match == 0],marker = 'o',facecolor = 'coral',edgecolor = 'coral',linewidths = edgewidth, zorder = 3, s = markersize, alpha = 0.5)
    ax4.scatter(Mvir_match_dm[sat_match != 0], Mstar_match[sat_match != 0],marker = 'D',facecolor = 'coral',edgecolor = 'coral',linewidths = edgewidth, zorder = 3, s = markersize, alpha = 0.5)   
    sort = ax4.scatter(Mvir_dm_sort[sat_sort == 0],Mstar_sort[sat_sort == 0],marker = 'o',facecolor = 'gold',edgecolor = 'gold',linewidths = edgewidth, zorder = 3, s = markersize, alpha = 0.5)
    ax4.scatter(Mvir_dm_sort[sat_sort != 0],Mstar_sort[sat_sort != 0],marker = 'D',facecolor = 'gold',edgecolor = 'gold',linewidths = edgewidth, zorder = 3, s = markersize, alpha = 0.5)   

    # Ratio of DM masses vs DM mass
    ax2.scatter(fdm_sim[(fdm_sim['simname'] == sim) & (fdm_sim['type']!='Satellite')]['Mhalo_z0'],fdm_sim[(fdm_sim['simname'] == sim) & (fdm_sim['type']!='Satellite')]['Mpeak']/fdm_sim[(fdm_sim['simname'] == sim) & (fdm_sim['type']!='Satellite')]['Mhalo_z0'],marker = 'o',facecolor = 'none',edgecolor = 'darkorange',linewidths = edgewidth*2, zorder = 2, s = markersize)
    ax2.scatter(fdm_sim[(fdm_sim['simname'] == sim) & (fdm_sim['type']=='Satellite')]['Mhalo_z0'],fdm_sim[(fdm_sim['simname'] == sim) & (fdm_sim['type']=='Satellite')]['Mpeak']/fdm_sim[(fdm_sim['simname'] == sim) & (fdm_sim['type']=='Satellite')]['Mhalo_z0'],marker = 'D',facecolor = 'none',edgecolor = 'darkorange',linewidths = edgewidth*2, zorder = 2, s = markersize)    
    #ax2.scatter(Mvir_match_dm_peak[sat_match[comm1] == 0],Mvir_match_dm_peak[sat_match[comm1] == 0]/Mvir_match[comm1][sat_match[comm1] == 0],marker = 'o',facecolor = 'coral',edgecolor = 'coral',linewidths = edgewidth, zorder = 3, s = markersize, alpha = 0.5)
    #ax2.scatter(Mvir_match_dm_peak[sat_match[comm1] != 0],Mvir_match_dm_peak[sat_match[comm1] != 0]/Mvir_match[comm1][sat_match[comm1] != 0],marker = 'D',facecolor = 'coral',edgecolor = 'coral',linewidths = edgewidth, zorder = 3, s = markersize, alpha = 0.5)
    ax2.scatter(Mvir_match[sat_match == 0],Mvir_match_dm[sat_match == 0]/Mvir_match[sat_match == 0],marker = 'o',facecolor = 'coral',edgecolor = 'coral',linewidths = edgewidth, zorder = 3, s = markersize, alpha = 0.5)
    ax2.scatter(Mvir_match[sat_match != 0],Mvir_match_dm[sat_match != 0]/Mvir_match[sat_match != 0],marker = 'D',facecolor = 'coral',edgecolor = 'coral',linewidths = edgewidth, zorder = 3, s = markersize, alpha = 0.5)   
    ax2.scatter(Mvir_sort[sat_sort == 0],Mvir_dm_sort[sat_sort == 0]/Mvir_sort[sat_sort == 0],marker = 'o',facecolor = 'gold',edgecolor = 'gold',linewidths = edgewidth, zorder = 3, s = markersize, alpha = 0.5)
    ax2.scatter(Mvir_sort[sat_sort != 0],Mvir_dm_sort[sat_sort != 0]/Mvir_sort[sat_sort != 0],marker = 'D',facecolor = 'gold',edgecolor = 'gold',linewidths = edgewidth, zorder = 3, s = markersize, alpha = 0.5)    
    legend = ax2.legend([peak,dmo,sort],['Peak M$_{vir}$','DMO M$_{vir}$','Sorted DMO M$_{vir}$'])
    #ax.axis([3e8,2e10,0.1,40])

    # Ratio of DM masses vs stellar mass
    ax3.scatter(fdm_sim[(fdm_sim['simname'] == sim) & (fdm_sim['type']!='Satellite')]['Mstar_z0'],fdm_sim[(fdm_sim['simname'] == sim) & (fdm_sim['type']!='Satellite')]['Mpeak']/fdm_sim[(fdm_sim['simname'] == sim) & (fdm_sim['type']!='Satellite')]['Mhalo_z0'],marker = 'o',facecolor = 'none',edgecolor = 'darkorange',linewidths = edgewidth*2, zorder = 2, s = markersize)
    ax3.scatter(fdm_sim[(fdm_sim['simname'] == sim) & (fdm_sim['type']=='Satellite')]['Mstar_z0'],fdm_sim[(fdm_sim['simname'] == sim) & (fdm_sim['type']=='Satellite')]['Mpeak']/fdm_sim[(fdm_sim['simname'] == sim) & (fdm_sim['type']=='Satellite')]['Mhalo_z0'],marker = 'D',facecolor = 'none',edgecolor = 'darkorange',linewidths = edgewidth*2, zorder = 2, s = markersize)    
    #ax3.scatter(Mstar_match[comm1][sat_match[comm1] == 0],Mvir_match_dm_peak[sat_match[comm1] == 0]/Mvir_match[comm1][sat_match[comm1] == 0],marker = 'o',facecolor = 'm',edgecolor = 'm',linewidths = edgewidth, zorder = 3, s = markersize, alpha = 0.5)
    #ax3.scatter(Mstar_match[comm1][sat_match[comm1] != 0],Mvir_match_dm_peak[sat_match[comm1] != 0]/Mvir_match[comm1][sat_match[comm1] != 0],marker = 'D',facecolor = 'm',edgecolor = 'm',linewidths = edgewidth, zorder = 3, s = markersize, alpha = 0.5)
    ax3.scatter(Mstar_match[sat_match == 0],Mvir_match_dm[sat_match == 0]/Mvir_match[sat_match == 0],marker = 'o',facecolor = 'coral',edgecolor = 'coral',linewidths = edgewidth, zorder = 3, s = markersize, alpha = 0.5)
    ax3.scatter(Mstar_match[sat_match != 0],Mvir_match_dm[sat_match != 0]/Mvir_match[sat_match != 0],marker = 'D',facecolor = 'coral',edgecolor = 'coral',linewidths = edgewidth, zorder = 3, s = markersize, alpha = 0.5)    
    ax3.scatter(Mstar_sort[sat_sort == 0],Mvir_dm_sort[sat_sort == 0]/Mvir_sort[sat_sort == 0],marker = 'o',facecolor = 'gold',edgecolor = 'gold',linewidths = edgewidth, zorder = 3, s = markersize, alpha = 0.5)
    ax3.scatter(Mstar_sort[sat_sort != 0],Mvir_dm_sort[sat_sort != 0]/Mvir_sort[sat_sort != 0],marker = 'D',facecolor = 'gold',edgecolor = 'gold',linewidths = edgewidth, zorder = 3, s = markersize, alpha = 0.5)
    legend = ax3.legend([peak,dmo,sort],['Peak M$_{vir}$','DMO M$_{vir}$','Sorted DMO M$_{vir}$'])

    # SMHM for sorted values
    peak = ax5.scatter(fdm_sim[fdm_sim['type']!='Satellite']['Mpeak'],fdm_sim[fdm_sim['type']!='Satellite']['Mstar_z0'],marker = 'o',facecolor = 'none',edgecolor = 'darkorange',linewidths = 2*edgewidth, zorder = 2, s = markersize)
    ax5.scatter(fdm_sim[fdm_sim['type']=='Satellite']['Mpeak'],fdm_sim[fdm_sim['type']=='Satellite']['Mstar_z0'],marker = 'D',facecolor = 'none',edgecolor = 'darkorange',linewidths = 2*edgewidth, zorder = 2, s = markersize)
    sort = ax5.scatter(Mvir_dm_sort[sat_sort == 0],Mstar_sort[sat_sort == 0],marker = 'o',facecolor = 'gold',edgecolor = 'gold',linewidths = edgewidth, zorder = 3, s = markersize, alpha = 0.5)
    ax5.scatter(Mvir_dm_sort[sat_sort != 0],Mstar_sort[sat_sort != 0],marker = 'D',facecolor = 'gold',edgecolor = 'gold',linewidths = edgewidth, zorder = 3, s = markersize, alpha = 0.5)     

    #SMHM for DMO matched
    peak = ax6.scatter(fdm_sim[fdm_sim['type']!='Satellite']['Mpeak'],fdm_sim[fdm_sim['type']!='Satellite']['Mstar_z0'],marker = 'o',facecolor = 'none',edgecolor = 'darkorange',linewidths = 2*edgewidth, zorder = 2, s = markersize)
    ax6.scatter(fdm_sim[fdm_sim['type']=='Satellite']['Mpeak'],fdm_sim[fdm_sim['type']=='Satellite']['Mstar_z0'],marker = 'D',facecolor = 'none',edgecolor = 'darkorange',linewidths = 2*edgewidth, zorder = 2, s = markersize)
    dmo = ax6.scatter(Mvir_match_dm[sat_match == 0], Mstar_match[sat_match == 0],marker = 'o',facecolor = 'coral',edgecolor = 'coral',linewidths = edgewidth, zorder = 3, s = markersize, alpha = 0.5)
    ax6.scatter(Mvir_match_dm[sat_match != 0], Mstar_match[sat_match != 0],marker = 'D',facecolor = 'coral',edgecolor = 'coral',linewidths = edgewidth, zorder = 3, s = markersize, alpha = 0.5)   


    
fig1.savefig(outfile_base + '_SMHM_multiMvir.png',dpi = dpi)
fig4.savefig(outfile_base + '_SMHM_multiMvir_read.png',dpi = dpi)
fig2.savefig(outfile_base + '_multiMvir_Mvir.png',dpi = dpi)
fig3.savefig(outfile_base + '_multiMvir_Mstar.png',dpi = dpi)

fig5.savefig(prefix_outfile + 'JL_Marvel_SMHM_abundsort.png',dpi = dpi)
fig6.savefig(prefix_outfile + 'JL_Marvel_SMHM_abundmatch.png',dpi = dpi)