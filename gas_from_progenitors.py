# Charlotte Christensen
# 3/2/22
# This program analyzes the gas that was accreted within satellies/progenitors (i.e., clumpy accretion)
# It uses information from prog_gas_acc.py


#mpl.use('tkagg') #Also can try mpl.use('Agg') #for using mpl over ssh    
import matplotlib as mpl
import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
import pynbody
import numpy as np
import pandas as pd
import socket, sys, os, glob, pickle
import matplotlib.colors as colors
sys.path.append(os.path.abspath("/home/christenc/Code/python/python_analysis/"))

def match_arr(x, y):
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
    
    return result_x, result_y # elements of x; elements of y
    

datapath = "/home/christenc/Code/students/DebPathak/mergers/"
simpath = "/home/christenc/Data/Sims/"
outfile_base =  '/home/christenc/Figures/marvel/'
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

finalstep = 4096
    
gas_data = pd.read_csv("/home/christenc/Code/Datafiles/gas_accr_NearMintJL.csv")

gas_sims = ["h" + str(uniqID)[:3] for uniqID in gas_data['uniqID']]

#tstep_disrupt = np.array([int(str(uniqID)[3:7]) for uniqID in gas_data['uniqID']])

#gas_data['tstep_disrupt'] = tstep_disrupt
# sim_unit = pynbody.load(filenames['h329']  + ".00" + str(finalstep))
# sim_unit['mass'].units # check the mass units
gas_data['metal_mass'] = gas_data['mass']*(2.09*gas_data['OxMassFrac'] + 1.06*gas_data['FeMassFrac'])*1.59e+16

# number of gas particles accreted from surviving satellites:
# np.sum(gas_data['tstep_disrupt'] == finalstep)

for filename in filenames.keys():
    print(filenames[filename]  + ".00" + str(finalstep))
    sim = pynbody.load(filenames[filename]  + ".00" + str(finalstep))
    sim.physical_units()
    h = sim.halos()
    halo = h.load_copy(1)
    #halo.gas['metal_mass'] = halo.gas['mass'].in_units('Msol')*(2.09*halo.gas['OxMassFrac'] + 1.06*halo.gas['FeMassFrac'])
    halo['metal_mass'] = halo['mass'].in_units('Msol')*(2.09*halo['OxMassFrac'] + 1.06*halo['FeMassFrac'])
    halo['iord_match'] = halo['iord']
    halo.star['iord_match'] = halo.star['igasorder']

    rho_cut = '0.1 m_p cm**-3' # m_p cm**-3
    T_cut = '2e4 K' # K
    f_hot = pynbody.filt.HighPass('temp', T_cut)
    f_rarified = pynbody.filt.LowPass('rho', rho_cut)
    cgm = halo.gas[f_hot | f_rarified]
    
    # What mass in z = 0, halo 1 gas came from (surviving) satellites?
        # What fraction of CGM/disk gas in halo 1 came from (surviving) satellites
    # What mass in z = 0, halo 1 stars came from (surviving) satellites
    accr_gas_data = gas_data[(gas_data['tstep_disrupt'] == finalstep) & (gas_data['Simulation'] == filename)]
    # Match the gas accreted from satellites and the stars formed from that gas
    result_halo1_sat, result_sat = match_arr(halo['iord_match'], np.array(accr_gas_data['iord'].tolist()))

    print("\nMass accreted from satellites: {0:4.2e} Msol".format(np.sum(halo[~result_halo1_sat.mask]['mass'].in_units('Msol'))))
    print("Percent of main halo from satellite\n\tTotal: {0:4.2f}%; CGM: {1:4.2f}%; Disk Gas: {2:4.2f}%; Stars: {3:4.2f}%".format(
        np.sum(halo[~result_halo1_sat.mask]['mass'].in_units('Msol'))/np.sum(halo['mass'].in_units('Msol'))*100,
        np.sum(halo[~result_halo1_sat.mask].gas[f_hot | f_rarified]['mass'].in_units('Msol'))/np.sum(halo.gas[f_hot | f_rarified]['mass'].in_units('Msol'))*100,
        np.sum(halo[~result_halo1_sat.mask].gas[~f_hot & ~f_rarified]['mass'].in_units('Msol'))/np.sum(halo.gas[~f_hot & ~f_rarified]['mass'].in_units('Msol'))*100,
        np.sum(halo[~result_halo1_sat.mask].stars['mass'].in_units('Msol'))/np.sum(halo.stars['mass'].in_units('Msol'))*100))
    print("Distribution of material accreted from satellites\n\tCGM: {0:4.2f}%; Disk Gas: {1:4.2f}%; Stars: {2:4.2f}%".format(
        np.sum(halo[~result_halo1_sat.mask].gas[f_hot | f_rarified]['mass'].in_units('Msol'))/np.sum(halo[~result_halo1_sat.mask]['mass'].in_units('Msol'))*100,
        np.sum(halo[~result_halo1_sat.mask].gas[~f_hot & ~f_rarified]['mass'].in_units('Msol'))/np.sum(halo[~result_halo1_sat.mask]['mass'].in_units('Msol'))*100,
        np.sum(halo[~result_halo1_sat.mask].stars['mass'].in_units('Msol'))/np.sum(halo[~result_halo1_sat.mask]['mass'].in_units('Msol'))*100))
    print("\nMetal Mass accreted from satellites: {0:4.2e} Msol".format(np.sum(accr_gas_data[~result_sat.mask]['metal_mass'])))
    print("Percent of main halo metals from satellite\n\tTotal: {0:4.2f}%".format(
        np.sum(accr_gas_data[~result_sat.mask]['metal_mass'])/np.sum(halo['metal_mass'].in_units('Msol'))*100))

    
    prog_gas_data = gas_data[(gas_data['tstep_disrupt'] != finalstep) & (gas_data['Simulation'] == filename)]
    result_halo1_prog, result_prog = match_arr(halo['iord_match'], np.array(prog_gas_data['iord'].tolist()))
    print("\nMass accreted from progenitors: {0:4.2e} Msol".format(np.sum(halo[~result_halo1_prog.mask]['mass'].in_units('Msol'))))
    print("Percent of main halo from progenitors\n\tTotal: {0:4.2f}%; CGM: {1:4.2f}%; Disk Gas: {2:4.2f}%; Stars: {3:4.2f}%".format(
        np.sum(halo[~result_halo1_prog.mask]['mass'].in_units('Msol'))/np.sum(halo['mass'].in_units('Msol'))*100,
        np.sum(halo[~result_halo1_prog.mask].gas[f_hot | f_rarified]['mass'].in_units('Msol'))/np.sum(halo.gas[f_hot | f_rarified]['mass'].in_units('Msol'))*100,
        np.sum(halo[~result_halo1_prog.mask].gas[~f_hot & ~f_rarified]['mass'].in_units('Msol'))/np.sum(halo.gas[~f_hot & ~f_rarified]['mass'].in_units('Msol'))*100,
        np.sum(halo[~result_halo1_prog.mask].stars['mass'].in_units('Msol'))/np.sum(halo.stars['mass'].in_units('Msol'))*100))
    print("Distribution of material accreted from progenitors\n\tCGM: {0:4.2f}%; Disk Gas: {1:4.2f}%; Stars: {2:4.2f}%".format(
        np.sum(halo[~result_halo1_prog.mask].gas[f_hot | f_rarified]['mass'].in_units('Msol'))/np.sum(halo[~result_halo1_prog.mask]['mass'].in_units('Msol'))*100,
        np.sum(halo[~result_halo1_prog.mask].gas[~f_hot & ~f_rarified]['mass'].in_units('Msol'))/np.sum(halo[~result_halo1_prog.mask]['mass'].in_units('Msol'))*100,
        np.sum(halo[~result_halo1_prog.mask].stars['mass'].in_units('Msol'))/np.sum(halo[~result_halo1_prog.mask]['mass'].in_units('Msol'))*100))
    print("\nMetal Mass accreted from progenitors: {0:4.2e} Msol".format(np.sum(prog_gas_data[~result_prog.mask]['metal_mass'])))
    print("Percent of main halo metals from progenitors\n\tTotal: {0:4.2f}%".format(
        np.sum(prog_gas_data[~result_prog.mask]['metal_mass'])/np.sum(halo['metal_mass'].in_units('Msol'))*100))

    x_bins = np.linspace(-7, 3, 50)
    y_bins = np.linspace(1, 8, 50)
    fig, ax = plt.subplots()
    ax.hist2d(np.log10(halo.gas['rho'].in_units('m_p cm**-3')), np.log10(halo.gas['temp']), weights = halo.gas['mass'], bins = [x_bins, y_bins], norm = colors.LogNorm())
    
    H, xedges, yedges = np.histogram2d(np.log10(halo[~result_halo1_sat.mask].gas['rho'].in_units('m_p cm**-3')), np.log10(halo[~result_halo1_sat.mask].gas['temp']), bins=(x_bins,y_bins), weights = halo[~result_halo1_sat.mask].gas['mass'])
    #x_bin_sizes = (xedges[1:] - xedges[:-1]).reshape((1,nbins_x))
    #y_bin_sizes = (yedges[1:] - yedges[:-1]).reshape((nbins_y,1))
    X, Y = 0.5*(xedges[1:]+xedges[:-1]), 0.5*(yedges[1:]+yedges[:-1])
    contour = ax.contour(X, Y, np.log10(H.T))
    fig.savefig(outfile_base + 'histrhoT' + filename + '.png', dpi = dpi)

    fig1, ax1 = plt.subplots()
    ax1.hist(np.log10(halo.gas['temp']), weights = halo.gas['mass'],bins = 20, range = (1,7))
    ax1.hist( np.log10(halo[~result_halo1_sat.mask].gas['temp']), weights = halo[~result_halo1_sat.mask].gas['mass'],alpha = 0.5, bins = 20, range = (1,7))
    ax1.set_xlabel('Log(T/1 K)')
    fig1.savefig(outfile_base + 'histT' + filename + '.png', dpi = dpi)

    x_bins = np.linspace(-7, 3, 50)
    y_bins = np.linspace(1, 8, 50)
    fig, ax = plt.subplots()
    ax.hist2d(np.log10(halo.gas['rho'].in_units('m_p cm**-3')), np.log10(halo.gas['temp']), weights = halo.gas['mass'], bins = [x_bins, y_bins], norm = colors.LogNorm())
    
    H, xedges, yedges = np.histogram2d(np.log10(halo[~result_halo1_sat.mask].gas['rho'].in_units('m_p cm**-3')), np.log10(halo[~result_halo1_sat.mask].gas['temp']), bins=(x_bins,y_bins), weights = halo[~result_halo1_sat.mask].gas['metal_mass'])
    #x_bin_sizes = (xedges[1:] - xedges[:-1]).reshape((1,nbins_x))
    #y_bin_sizes = (yedges[1:] - yedges[:-1]).reshape((nbins_y,1))
    X, Y = 0.5*(xedges[1:]+xedges[:-1]), 0.5*(yedges[1:]+yedges[:-1])
    contour = ax.contour(X, Y, np.log10(H.T))
    fig.savefig(outfile_base + 'histrhoTZ' + filename + '.png', dpi = dpi)

    fig1, ax1 = plt.subplots()
    ax1.hist(np.log10(halo.gas['temp']), weights = halo.gas['metal_mass'],bins = 20, range = (1,7))
    ax1.hist( np.log10(halo[~result_halo1_sat.mask].gas['temp']), weights = halo[~result_halo1_sat.mask].gas['metal_mass'],alpha = 0.5, bins = 20, range = (1,7))
    ax1.set_xlabel('Log(T/1 K)')
    fig1.savefig(outfile_base + 'histTZ' + filename + '.png', dpi = dpi)
    x_bins = np.linspace(-7, 3, 50)
    y_bins = np.linspace(1, 8, 50)
    fig, ax = plt.subplots()
    ax.hist2d(np.log10(halo.gas['rho'].in_units('m_p cm**-3')), np.log10(halo.gas['temp']), weights = halo.gas['metal_mass'], bins = [x_bins, y_bins], norm = colors.LogNorm())
    
