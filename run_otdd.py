# Calculate the Earth Mover/Wasserstin/Optimal Transport distance for pair of accreted satellites
# October 14th, 2021
# Charlotte Christensen

# Note to self: compare results with L1 statistical distance and the Wasserstein 2 distance

import os
import pandas as pd
import pynbody
import matplotlib.gridspec as gridspec
from matplotlib import pyplot as plt
import matplotlib.colors as colors
import ot # https://pythonot.github.io/
import ot.plot
import re
import socket
import numpy as np

if (socket.gethostname() == "quirm.math.grinnell.edu"):
    prefix = '/home/christenc/Data/Sims/'
    prefix_outfile = '/home/christenc/Figures/marvel/'
else:
    prefix = '/home/christensen/Storage2/UW/MolecH/Cosmo/'
    prefix_outfile = '/home/christensen/Plots/marvel/'

outfile_base = prefix_outfile + 'marvelJL'
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

# Path name to simulation
filepath = '/home/christenc/Data/Sims/h148.cosmo50PLK.3072g/h148.cosmo50PLK.3072g3HbwK1BH/'
filename = 'h148.cosmo50PLK.3072g3HbwK1BH'

# Read in list stars labeled by progenitor
infallstars = pd.read_csv(filepath + 'h148haloIDinfo_revised.csv')
infallstars['unique_halo_ID'] =  (infallstars['unique_halo_ID'].astype('int')).astype('string')

infall_step = {'infall_step': ""}
infallstars = infallstars.join(pd.DataFrame(columns=infall_step))
infall_id = {'infall_id': ""}
infallstars = infallstars.join(pd.DataFrame(columns=infall_id))

for haloid in infallstars['unique_halo_ID'].unique():
    # Print unique haloid, the infall step, the haloid, and the number of stars in it
   print(haloid,haloid[3:7],', ',haloid[7:],np.sum(infallstars['unique_halo_ID'] == haloid))
   infallstars.loc[infallstars['unique_halo_ID'] == haloid,'infall_step'] = haloid[3:7]
   infallstars.loc[infallstars['unique_halo_ID'] == haloid,'infall_id'] = haloid[7:]

# Select the timestep and the haloids at that timestep
unique_halo_ID_1 = '14807012'
unique_halo_ID_2 = '148070125'
simname = 'Sandra'
infallstep = '0701'
haloid_1 = '2'
haloid_2 =  '25'

'''
unique_halo_ID_1 = '148102421'
unique_halo_ID_2 = '148102429'
simname = 'Sandra'
infallstep = '1024'
haloid_1 = '21'
haloid_2 = '29'
'''

# Search simulation directory to find all available steps
dirs = os.listdir(filepath)
steps = []
pattern = re.compile(filename + "\.00....")
for dir in dirs:
    if pattern.fullmatch(dir):
        print((re.split(r'\.00',dir))[1])
        steps.append((re.split(r'\.00',dir))[1])
steps.sort()
start_step_ind = steps.index(infallstep) - 1

# Set up arrays to store distance information
dists = []
dists_com = []
times = []

# Setting up a grid of figures
plt.clf()
ncol = 5
nrows = int(len(steps[start_step_ind - 1:])/ncol)
if (mod(len(steps[start_step_ind - 1:]),ncol) > 0):
    nrows = nrows + 1
gs = gridspec.GridSpec(nrows, ncol)
plt_width = plt_width*2
fig1 = plt.figure(1,figsize = (plt_width,plt_width*nrows/ncols))
fig1.set_size_inches(plt_width,plt_width*nrows/ncols, forward=True)
fig2 = plt.figure(2,figsize = (plt_width,plt_width*nrows/ncols))
fig3 = plt.figure(3,figsize = (plt_width,plt_width*nrows/ncols))
fig4 = plt.figure(4,figsize = (plt_width,plt_width*nrows/ncols))

ct = 0
for step in steps[start_step_ind - 1:]:
    # Load that timestep using pynbody
    print(filepath + filename + '.00' + step + '/' + filename + '.00' + step)
    sim = pynbody.load(filepath + filename + '.00' + step + '/' + filename + '.00' + step)
    h  = sim.halos(dummy=True))

    mainhalo = 1 # Need to use main progenitor branch to actually find the halo id for the main halo at a given snapshot

 ##### Label pynbody star objects with the appropriate unique halo ID
    stars = sim.stars[pynbody.filt.HighPass('tform','0 Gyr')] # Exclude BHs
    '''
    # This selection only works for the snapshot steps[start_step_ind - 1]
    halo_1 = h[int(haloid_1)]
    halo_2 = h[int(haloid_2)]
    halo_1.physical_units()
    halo_2.physical_units()
    '''
    
    # current_infallstars = infallstars[:len(stars.star)] #Only cuts off the later stars, presumably formed after the snapshot. If stars are missing earlier on so I'm not sure that this works. The method below should be more robust
    
    # Create an array for the stars corresponding to the unique halo ids
    x = np.array(infallstars['iord']) # Array of iords for the stars in the infall array
    y = np.array(stars.star['iord']) # Array of iords for stars in the current snapshot
    
    index = np.argsort(x)
    sorted_x = x[index] # Sorted list of iords for intall array
    sorted_index = np.searchsorted(sorted_x, y)
    yindex = np.take(index, sorted_index, mode="clip")
    mask = x[yindex] != y
    result = np.ma.array(yindex, mask=mask)
    current_infallstars = infallstars.loc[result[~result.mask].data] # Array of all iords
    
    stars.star['unique_Halo_ID'] = list(map(int, (current_infallstars['unique_halo_ID'][:len(stars.star)]).tolist()))

    # Use the unique halo ids to select the stars in each of the halos of interest
    indicies = np.where(stars.star['unique_Halo_ID'] == int(unique_halo_ID_1))
    halo_1 = stars.star[indicies]
    halo_1.physical_units()
    indicies = np.where(stars.star['unique_Halo_ID'] == int(unique_halo_ID_2))
    halo_2 = stars.star[indicies]
    halo_2.physical_units()
    
    # Make a plot of the stars in the halos at that timestep colored by Fe/H
    cm = plt.cm.get_cmap('gist_ncar')
    cNorm  = colors.Normalize(vmin=-6, vmax = 1)
    
    FeH_1 = halo_1.star['feh']
    idx_sort_1 = np.argsort(FeH_1)[::1]
    FeH_2 = halo_2.star['feh']
    idx_sort_2 = np.argsort(FeH_2)[::1]

    ax1 = fig1.add_subplot(gs[ct])
    ax1.scatter(halo_1.star['x'][idx_sort_1],halo_1.star['y'][idx_sort_1], s=0.5, c=FeH_1[idx_sort_1], cmap=cm, norm=cNorm)
    ax1.scatter(halo_2.star['x'][idx_sort_2],halo_2.star['y'][idx_sort_2], s=0.5, c=FeH_1[idx_sort_2], cmap=cm, norm=cNorm)
    ax1.text(0.1,0.8,"{:5.2f}".format(sim.properties['time'].in_units('Gyr')),transform=ax1.transAxes)
    fig1.show()

    # Calculate the center of mass distance
    halo_1_com = np.array([np.sum(halo_1.star['x']*halo_1.star['mass']),np.sum(halo_1.star['y']*halo_1.star['mass']),np.sum(halo_1.star['z']*halo_1.star['mass'])])/float(np.sum(halo_1.star['mass']))
    halo_2_com = np.array([np.sum(halo_2.star['x']*halo_2.star['mass']),np.sum(halo_2.star['y']*halo_2.star['mass']),np.sum(halo_2.star['z']*halo_2.star['mass'])])/float(np.sum(halo_2.star['mass']))
    dists_com.append(np.sqrt(np.sum((halo_1_com - halo_2_com)*(halo_1_com - halo_2_com)))/h[1].properties['Rvir'])
    print(f'COM D = {dists_com[-1]}')
        
    # Draw randomly from the distribution
    n = 100

    ind_1 = np.sort(np.random.choice(np.arange(0,len(halo_1.star)),n,p=np.array(halo_1.star['mass']/np.sum(halo_1.star['mass']))))
    #ind_1a = np.sort(np.random.choice(np.arange(0,len(halo_1.star)),n,p=np.array(halo_1.star['mass']/np.sum(halo_1.star['mass']))))     # Draw randomly from the same distribution to check method
    ind_2 = np.sort(np.random.choice(np.arange(0,len(halo_2.star)),n,p=np.array(halo_2.star['mass']/np.sum(halo_2.star['mass']))))

    data1 = (np.array([np.array(halo_1.star[ind_1]['x']),np.array(halo_1.star[ind_1]['y']),np.array(halo_1.star[ind_1]['z'])]).T)/h[1].properties['Rvir'] # Scale distances by the virial radius of the host
    #data1a =  np.array([np.array(halo_1.star[ind_1a]['x']),np.array(halo_1.star[ind_1a]['y']),np.array(halo_1.star[ind_1a]['z'])]).T
    data2 = (np.array([np.array(halo_2.star[ind_2]['x']),np.array(halo_2.star[ind_2]['y']),np.array(halo_2.star[ind_2]['z'])]).T)/h[1].properties['Rvir'] 

    # Scale velocities by h[1].properties['Vmax']
    
    # Calculate cost matrix
    M = ot.dist(data1,data2)
    #M = ot.dist(data1,data1a)

    ax2 = fig2.add_subplot(gs[ct])
    ax2.imshow(M,interpolation='nearest')
    ax2.axes.xaxis.set_ticks([])
    ax2.axes.yaxis.set_ticks([])   
    ax2.text(0.1,0.8,step,transform=ax2.transAxes)
    fig2.show()

    # Use cost matrix to calculate optimal transport distance
    a, b = np.ones((n,)) / n, np.ones((n,)) / n  # uniform distribution on samples
    G0 = ot.emd(a, b, M)

    # Plot pairs for optimal transport distance
    ax3 = fig3.add_subplot(gs[ct])
    ax3.imshow(G0, interpolation='nearest')
    #plt.title('OT matrix G0')
    ax3.axes.xaxis.set_ticks([])
    ax3.axes.yaxis.set_ticks([])
    ax3.text(0.1,0.8,step,transform=ax3.transAxes)
    fig3.show()

    # Plot visualization for optimal transport distance
    ax4 = fig4.add_subplot(gs[ct])
    #ot.plot.plot2D_samples_mat(data1, data1a, G0, c=[.5, .5, 1])
    ax4.plot(data1[:, 0], data1[:, 1], '+b', label='Source samples')
    #plt.plot(data1a[:, 0], data1a[:, 1], '+r', label='Source samples')
    ax4.plot(data2[:, 0], data2[:, 1], 'xr', label='Target samples')
    ot.plot.plot2D_samples_mat(data1, data2, G0, c=[.5, .5, 1])
    #plt.legend(loc=0)
    #plt.title('OT matrix with samples')
    ax4.text(0.1,0.8,step,transform=ax4.transAxes)
    fig4.show()

    # Add costs for each pair to calculate total optimal transport distance
    dists.append(np.sum(G0*M))
    times.append(sim.properties['time'].in_units('Gyr'))

    print(f'OTDD = {dists[-1]}')

    ct = ct+1 #Iterate figure

fig1.tight_layout()
fig1.show()
fig1.savefig(outfile_base + 'image_' + unique_halo_ID_1 + '_' + unique_halo_ID_2 + '.png',dpi = dpi)

fig3.tight_layout()
fig3.show()
fig3.savefig(outfile_base + 'OTDDgrid_' + unique_halo_ID_1 + '_' + unique_halo_ID_2 + '.png',dpi = dpi)

fig2.tight_layout()
fig2.show()
fig2.savefig(outfile_base + 'costgrid_' + unique_halo_ID_1 + '_' + unique_halo_ID_2 + '.png',dpi = dpi)

fig4.tight_layout()
fig4.show()
fig4.savefig(outfile_base + 'OTDDimage_' + unique_halo_ID_1 + '_' + unique_halo_ID_2 + '.png',dpi = dpi)

# Plot distances over time
fig5 = plt.figure(5,figsize = (plt_width,plt_width*aspect_ratio))
ax5 = fig5.add_subplot()
otdd_line = ax5.plot(times,dists,color='k', linewidth=edgewidth*3, label="Distance")
comd_line = ax5.plot(times,dists_com,color='k', linewidth=edgewidth*3, linestyle = 'dashed', label="Time [Gyr]")
ax5.set_ylabel(r'Distance')
ax5.set_xlabel(r'Time [Gyr]')
ax5.set_title(simname + ' ' + infallstep + ': ' + haloid_1 + ',  ' + haloid_2)
ax5.set_yscale('log')
legend5 = ax5.legend(['OTD Distance','COM Distance'],frameon = False) #[otdd_line, comd_line]) #,['OTD Distance','COM Distance']) #,loc = 4,framealpha = 0,frameon = False)
fig5.savefig(outfile_base + 'OTDD_COM_vtime_' + unique_halo_ID_1 + '_' + unique_halo_ID_2 + '.png',dpi = dpi)

