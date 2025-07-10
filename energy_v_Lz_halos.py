# Calculate the Earth Mover/Wasserstin/Optimal Transport distance for pair of accreted satellites
# October 14th, 2021
# Charlotte Christensen

# Note to self: compare results with L1 statistical distance and the Wasserstein 2 distance

import os
import pandas as pd
import pynbody
import tangos
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import ot # https://pythonot.github.io/
import ot.plot
import re
import socket
import matplotlib.colors as colors
from matplotlib import path
import numpy as np
from numpy.random import default_rng
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap
try:
    from astropy.convolution import Gaussian2DKernel, convolve
    astro_smooth = True
except ImportError as IE:
    astro_smooth = False

    
if (socket.gethostname() == "quirm.math.grinnell.edu"):
    prefix = '/home/christenc/Data/Sims/'
    prefix_outfile = '/home/christenc/Figures/marvel/'
else:
    prefix = '/home/christensen/Storage2/UW/MolecH/Cosmo/'
    prefix_outfile = '/home/christensen/Plots/marvel/'

outfile_base = prefix_outfile
presentation = False
if presentation:
    outfile_base = outfile_base + 'pres'
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
"""
filepath = '/home/christenc/Data/Sims/h148.cosmo50PLK.3072g/h148.cosmo50PLK.3072g3HbwK1BH/'
filename = 'h148.cosmo50PLK.3072g3HbwK1BH'
simname = 'Sandra'
simshort = 'h148'
"""

filepath = '/home/christenc/Data/Sims/h229.cosmo50PLK.3072g/h229.cosmo50PLK.3072gst5HbwK1BH/'
filename = 'h229.cosmo50PLK.3072gst5HbwK1BH'
simname = 'Ruth'
simshort = 'h229'


outfile_base = outfile_base + simname

# Use Tangos to get the main progenitor branch
tangos.init_db(filepath + simshort + '.db')

tangos_halo = tangos.get_halo('snapshots/' + filename + '.004096/halo_1')
timesteps = tangos.get_simulation("snapshots").timesteps
snap_nums = [re.findall(r'.00+[\d]+', str(snap))[0][3:] for snap in timesteps]
#time, ids, tss = tangos_halo.calculate_for_progenitors("t()", "halo_number()", "timestep_id()")
time, ids = tangos_halo.calculate_for_progenitors("t()", "halo_number()")
snap_nums = snap_nums[len(snap_nums) - len(time):] # Truncate the steps where the halo is not traced. Could be done more robustly

# Read in energy offset
e_offset = pd.read_csv(filepath + simshort + '_e_lz_offsets.csv')

# Read in list stars labeled by progenitor
infallstars = pd.read_csv(filepath + simshort + 'haloIDinfo_revised.csv')
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
timestep = '4096'

# Search simulation directory to find all available steps

dirs = os.listdir(filepath)
steps = []
pattern = re.compile(filename + "\.00....")
for dir in dirs:
    if pattern.fullmatch(dir):
        print((re.split(r'\.00',dir))[1])
        steps.append((re.split(r'\.00',dir))[1])
steps.sort()
'''
start_step_ind = steps.index(infallstep) - 1
'''
#steps = [timestep]

# Setting up a grid of figures
plt.ion()
plt.clf()
'''
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
'''

ct = 0
for step in reversed(steps):
    try:
        tangos_ind = snap_nums.index(step)
    except:
        print(step + ' not in tangos db')
        continue

    # Load that timestep using pynbody
    print(filepath + filename + '.00' + step + '/' + filename + '.00' + step)
    sim = pynbody.load(filepath + filename + '.00' + step + '/' + filename + '.00' + step)
    h  = sim.halos(dummy=True)

    mainhalo = (np.flip(ids))[tangos_ind] # Need to use main progenitor branch to actually find the halo id for the main halo at a given snapshot

 ##### Label pynbody star objects with the appropriate unique halo ID
    stars = sim.stars[pynbody.filt.HighPass('tform','0 Gyr')] # Exclude BHs
    pynbody.analysis.halo.center(h[mainhalo])
    if ct == 0:
        if (len(h[mainhalo].gas) > 0):
            cen = h[mainhalo].gas[pynbody.filt.Sphere("5 kpc")]
        else:
            cen = h[mainhalo][pynbody.filt.Sphere("5 kpc")]
        #trans = pynbody.analysis.calc_sideon_matrix(ang_mom_vec(cen))
        angmom = (cen['mass'].reshape((len(cen), 1)) *
            np.cross(cen['pos'], cen['vel'])).sum(axis=0).view(np.ndarray)
        vec_in = np.asarray(angmom)
        vec_in = vec_in / np.sum(vec_in ** 2).sum() ** 0.5
        vec_p1 = np.cross([1, 0, 0], vec_in)
        vec_p1 = vec_p1 / np.sum(vec_p1 ** 2).sum() ** 0.5
        vec_p2 = np.cross(vec_in, vec_p1)
        matr = np.concatenate((vec_p2, vec_in, vec_p1)).reshape((3, 3))
        
        #pynbody.analysis.angmom.faceon(h[mainhalo], cen=(0,0,0))
    sim = pynbody.transformation.transform(sim, matr)
        
    # current_infallstars = infallstars[:len(stars.star)] #Only cuts off the later stars, presumably formed after the snapshot. If stars are missing earlier on so I'm not sure that this works. The method below should be more robust
    
    # Create an array for the stars corresponding to the unique halo ids
    x = np.array(infallstars['iord']) # Array of iords for the stars in the infall array
    y = np.array(stars.star['iord']) # Array of iords for stars in the current snapshot
    
    index = np.argsort(x)
    sorted_x = x[index] # Sorted list of iords for infall array
    sorted_index = np.searchsorted(sorted_x, y)
    yindex = np.take(index, sorted_index, mode="clip")
    mask = x[yindex] != y
    result = np.ma.array(yindex, mask=mask)
    current_infallstars = infallstars.loc[result[~result.mask].data] # Array of all iords
    
    stars.star['unique_Halo_ID'] = list(map(int, (current_infallstars['unique_halo_ID'][:len(stars.star)]).tolist()))

    if ct == 0:
        unique_Halo_ID_list, unique_Halo_ID_list_ind, unique_Halo_ID_list_ct = np.unique(stars.star['unique_Halo_ID'], return_index=True, return_counts=True)
        infallsteps = np.array([int(str(a)[3:7]) for a in (unique_Halo_ID_list.tolist())])
        halo_select_cond = (unique_Halo_ID_list_ct >= 100) & (unique_Halo_ID_list_ct <= 1e5) & (infallsteps > 2./14.*4096) # Select halos within a certain number of stars and that accreted after 2 Gyr
        plt_order = np.argsort(infallsteps[halo_select_cond])
        #plt_order = (np.flip(numpy.argsort(unique_Halo_ID_list_ct)))[1:] # Reverse plot by number of stars at infall, ignoring the first one, since it presumably wasn't cut by mass

        cmapbase = plt.cm.turbo #plt.cm.get_cmap('turbo')
        nhalos = len(plt_order)
        #rng = default_rng()
        color_order = np.arange(nhalos)
        #rng.shuffle(color_order)    
    
    #cm = plt.cm.get_cmap('gist_ncar')

    #FeH = stars.star['feh']
    #idx_sort = np.argsort(FeH)[::1]
    #phasehalo1 = plt.scatter(stars.star['jz'][idx_sort]/(10**3), stars.star['te'][idx_sort]/(10**6), s=.5, vmin=-2, vmax=0.8, c=FeH[idx_sort], cmap=cm)

    # Plot Energy vs Lz
    fig1 = plt.figure(1) #,figsize = (plt_width,plt_width*aspect_ratio))
    ax1 = fig1.add_subplot(111)

    # Plot Polar coordinates   
    fig2 = plt.figure(2) #,figsize = (plt_width,plt_width*aspect_ratio))
    ax2 = fig2.add_subplot(111) #, projection='aitoff')
    #ax2.grid(True)
    
    # Phase space
    fig3 = plt.figure(3) #,figsize = (plt_width,plt_width*aspect_ratio))
    ax3 = fig3.add_subplot(111)
    
    halo_ct = 0
    ax1.clear()
    ax2.clear()
    ax3.clear()    
    for halo_ID_ind in plt_order[:]:
        # Use the unique halo ids to select the stars in each of the halos of interest
        halo_ID =  (unique_Halo_ID_list[halo_select_cond])[halo_ID_ind]
        #halo_ID = (unique_Halo_ID_list[halo_select_cond])[plt_order[halo_ct]]
        print('Unique Halo ID: ',halo_ID)
        indicies = np.where(stars.star['unique_Halo_ID'] == halo_ID)
        halostars = stars.star[indicies]
        halostars.physical_units()

        # Define a color map
        N = 256
        basecolor_it = np.int(np.mod(np.int(np.floor(cmapbase.N/nhalos)*color_order[halo_ct] + np.floor(cmapbase.N/nhalos/2)) + N/2,N))
        basecolor = cmapbase(basecolor_it) #np.array([0,0,1]) #blue
        print(basecolor_it, basecolor)
        my_cmap = np.empty((N,4))
        my_cmap[:,0] = basecolor[0]
        my_cmap[:,1] = basecolor[1]
        my_cmap[:,2] = basecolor[2]
        my_cmap[:,3] = np.linspace(0.05, 1, N)
        newcmp = ListedColormap(my_cmap)

        jz_range = [-60,60]
        e_range = [-0.25,0.05]
        x = np.array(halostars.star['jz']/(10**3))
        try:
            e_off = str((e_offset['00' + step])[0])
        except:
            e_off = str((e_offset['00' + np.flip(steps)[ct + 1]])[0])
        try:
            e_off = float(e_off)
        except:
            disallowed_characters = "[]"
            for character in disallowed_characters:
                e_off = e_off.replace(character, "")
            e_off = float(e_off)

        y = np.array(halostars.star['te']/(10**6)) - (e_off/10**6)
        H, xedges, yedges = np.histogram2d(x,y, bins=(100,80), range = [jz_range, e_range])
        print('Lz vs E density: ',np.log10(len(halostars)),np.log10(np.max(H)))
        xmesh, ymesh = np.meshgrid(xedges[:-1] + (xedges[1] - xedges[0]) / 2, yedges[:-1] + (yedges[1] - yedges[0]) / 2)

        # Smooth the contours (if astropy is installed)
        if astro_smooth:
            kernel = Gaussian2DKernel(stddev=1.)
            H=convolve(H,kernel)
        
        #fig,ax = plt.subplots(1, figsize=(7,6))
        levels = np.array([0.2,1.2,2.2,3.2,4.2,5.2])
        levels = np.array([0.5,1,1.5,2,2.5,3,3.5,4])
        levels = np.array([0.001, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
        levels = np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
        #clevels = ax1.contourf(xmesh,ymesh,np.log10(H.T),cmap=newcmp,levels = levels,extend=max)
        #ax1.contour(xmesh,ymesh,np.log10(H.T),cmap=newcmp,levels = levels,zorder = 1)
        if (np.max(H) > 10**levels[0]):
            clevels = ax1.contourf(xmesh,ymesh,H.T,cmap=newcmp,levels = 10**levels,extend='max')
            ax1.contour(xmesh, ymesh, H.T, levels=10**levels, zorder=1, colors=[basecolor], alpha=0.2)    #cmap=newcmp,                     
            # Identify points within contours
            inside = np.full_like(x,False,dtype=bool)
            for collect in clevels.collections:
                p = collect.get_paths()
                for level in p:
                    inside |= level.contains_points((np.array([x,y])).T)
            
            ax1.plot(x[~inside], y[~inside], marker='.', color=basecolor, alpha=0.2, lw=0, linestyle="")
        else:
            ax1.plot(x, y, marker='.', color=basecolor, alpha=0.2, markeredgewidth=0, lw=0, linestyle="")
        #ax1.show(block=False)

        halostars['theta'] = np.arcsin(halostars['y']/halostars['r'])
        halostars['phi_test'] = np.arctan2(halostars['z'],halostars['x'])

        phi_range = [-np.pi,np.pi]
        theta_range =[-np.pi/2,np.pi/2] 
        x = np.array(halostars.star['phi_test'])
        y = np.array(halostars.star['theta'])
        #x = np.array(halostars.star['x'])
        #y = np.array(halostars.star['y'])
        H, xedges, yedges = np.histogram2d(x,y, bins=(100,80), range = [phi_range, theta_range])
        print('Stellar density in skymap: ',np.log10(len(halostars)),np.log10(np.max(H)))
        xmesh, ymesh = np.meshgrid(xedges[:-1] + (xedges[1] - xedges[0]) / 2, yedges[:-1] + (yedges[1] - yedges[0]) / 2)

        # Smooth the contours (if astropy is installed)
        if astro_smooth:
            kernel = Gaussian2DKernel(stddev=1.)
            H=convolve(H,kernel)
        
        #fig,ax = plt.subplots(1, figsize=(7,6))
        levels = np.array([0.5,1,1.5,2,2.5,3,3.5,4])
        levels = np.array([0.0001, 0.5, 1.0, 1.5, 2.0])
        levels = np.array([1, 1.5, 2, 2.5, 3])
        #clevels = ax2.contourf(xmesh,ymesh,np.log10(H.T),cmap=newcmp,levels = levels)#,zorder=90)
        #ax2.contour(xmesh,ymesh,np.log10(H.T),cmap=newcmp,levels = levels,zorder = 1, extend=max)
        if (np.max(H) > 10**levels[0]):
            clevels = ax2.contourf(xmesh,ymesh,H.T,cmap=newcmp,levels = 10**levels, extend='max')#,zorder=90)
            ax2.contour(xmesh, ymesh, H.T, levels=10**levels, zorder=1, colors=[basecolor], alpha=0.2)        #cmap=newcmp,                  
        # Identify points within contours
            inside = np.full_like(x,False,dtype=bool)
            for collect in clevels.collections:
                p = collect.get_paths()
                for level in p:
                    inside |= level.contains_points((np.array([x,y])).T)
            
            ax2.plot(x[~inside],y[~inside],marker = '.', color = basecolor, alpha = 0.1, lw=0, linestyle="")
        else:
            ax2.plot(x,y,marker = '.', color = basecolor, alpha = 0.1, lw=0, linestyle="")

        # Radial phase space
        r_range = [0,400]
        vr_range = [-750,750]
        x = np.array(halostars.star['r'])
        y = np.array(halostars.star['vr'])
        H, xedges, yedges = np.histogram2d(x,y, bins=(100,80), range = [r_range,vr_range])
        print('Radial phase space: ',np.log10(len(halostars)),np.log10(np.max(H)))
        xmesh, ymesh = np.meshgrid(xedges[:-1] + (xedges[1] - xedges[0]) / 2, yedges[:-1] + (yedges[1] - yedges[0]) / 2)

        # Smooth the contours (if astropy is installed)
        if astro_smooth:
            kernel = Gaussian2DKernel(stddev=1.)
            H=convolve(H,kernel)
        
        #fig,ax = plt.subplots(1, figsize=(7,6))
        levels = np.array([0.5, 1, 1.5, 2, 2.5, 3, 3.5, 5])
        levels = np.array([0.001, 0.5, 1, 1.5, 2, 2.5, 3])
        levels = np.array([0.5, 1, 1.5, 2, 2.5, 3])
        #clevels = ax3.contourf(xmesh,ymesh,np.log10(H.T),cmap=newcmp,levels = levels, extend=max)
        #ax3.contour(xmesh,ymesh,np.log10(H.T),cmap=newcmp,levels = levels,zorder = 1)
        if (np.max(H) > 10**levels[0]):
            clevels = ax3.contourf(xmesh,ymesh,H.T,cmap=newcmp,levels = 10**levels, extend='max')
            ax3.contour(xmesh, ymesh, H.T, levels = 10**levels, zorder = 1, colors=[basecolor], alpha=0.2) #cmap=newcmp, 
                         
            # Identify points within contours
            inside = np.full_like(x,False,dtype=bool)
            for collect in clevels.collections:
                p = collect.get_paths()
                for level in p:
                    inside |= level.contains_points((np.array([x,y])).T)
            
            ax3.plot(x[~inside],y[~inside],marker = '.', color = basecolor, alpha = 0.2, lw=0, linestyle="")
        else:
            ax3.plot(x,y,marker = '.', color = basecolor, alpha = 0.2, lw=0, linestyle="")

        halo_ct = halo_ct+1
    
    ct = ct+1 #Iterate figure

    ax1.set_xlabel(r'L$_z$')
    ax1.set_ylabel(r'Energy')
    ax1.set_title(step)
    ax1.axis([jz_range[0], jz_range[1], e_range[0], e_range[1]])
    fig1.tight_layout()
    fig1.show()
    fig1.savefig(outfile_base + 'ELz_' + step + '.png',dpi = dpi)

    ax2.set_title(step)
    fig2.tight_layout()
    fig2.show()
    fig2.savefig(outfile_base + 'skymapr_' + step + '.png',dpi = dpi)

    ax3.set_xlabel(r'r [kpc]')
    ax3.set_ylabel(r'V$_r$ [km s$^{-1}$]')
    ax3.set_title(step)
    ax3.axis([r_range[0], r_range[1], vr_range[0], vr_range[1]])
    fig3.tight_layout()
    fig3.show()
    fig3.savefig(outfile_base + 'phaseSpace_' + step + '.png',dpi = dpi)

