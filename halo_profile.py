#Charlotte Christensen

#2/13/18
#Generate profiles of the halo

#Run with
#%run /home/christensen/Code/python/python_analysis/halo_profile.py
#or, on quirm,
#%run /home/christenc/Code/python/python_analysis/halo_profile.py
#ipython --pylab


import numpy as np
import matplotlib.pylab as plt
import matplotlib.mlab as mlab
import matplotlib.colors as colors
import matplotlib.cm as cm
import pynbody
import pynbody.plot as pp
import pynbody.plot.sph
import pynbody.filt as filt
import pynbody.units as units
import pynbody.analysis.profile as profile
import sys, os, glob, pickle
import socket

def halo_profile(tfile,halo_nums):
    s = pynbody.load(tfile)
    s.physical_units()
    h = s.halos()
    f_bar = 0.16510

    i = 0
    cmx = plt.get_cmap("tab20")
    values = range(0,20)
    cNorm  = colors.Normalize(vmin=0, vmax = values[-1]+1)
    scalarMap = cm.ScalarMappable(norm=cNorm, cmap=cmx)
    #f, axs = plt.subplots(1,3,figsize=(20,6))
    fig = plt.figure(4,figsize=(20.,6.))
    #fig.subplots_adjust(left=0.08, bottom=0.05, right=0.9, top=0.97, wspace=0.12, hspace=0.13)
    axs = [fig.add_subplot(1,3,1), fig.add_subplot(1,3,2), fig.add_subplot(1,3,3)]
        
    for halo_num in halo_nums:
        print(halo_num)
        halo = h[int(halo_num)]
        #pynbody.analysis.halo.center(halo,vel= False)
        pynbody.analysis.angmom.faceon(h[int(halo_num)])
        
        rvir = pynbody.array.SimArray(np.sqrt(np.max(halo['x'].in_units('kpc')**2 + halo['y'].in_units('kpc')**2 + halo['z'].in_units('kpc')**2)),'kpc')
        centre = (0,0,0)
        sphere = halo[pynbody.filt.Sphere(rvir, centre)]

        fig1 = plt.figure(1)
        ax1 = fig1.add_subplot(111)
        pynbody.plot.image(s.g, width=rvir*3, cmap='Blues',subplot=ax1, show_cbar=False);
        circle_rvir = plt.Circle((0, 0), rvir, color='k',fill=False)
        ax1.add_artist(circle_rvir)
        fig1.show
        fig1.savefig(tfile+'.'+halo_num+'.imageDens.png')
        
        fig2 = plt.figure(2)
        ax2 = fig2.add_subplot(111)
        pynbody.plot.image(s.g,qty="temp",width=rvir*3,cmap="YlOrRd", denoise=True,approximate_fast=False,subplot=ax2, show_cbar=False);
        #ax2.add_artist(circle_rvir)

        fig3 = plt.figure(3)
        ax3 = fig3.add_subplot(111)
        pynbody.plot.image(s.g,qty="temp",width=rvir*3,cmap="YlOrRd", denoise=True,approximate_fast=False,vmin=1e3,vmax=1e7,subplot=ax3, show_cbar=False);
        fig3.savefig(tfile+'.'+halo_num+'.imageT.png')

        den_range = [-8,2.0]
        temp_range = [2,9]
        fig5 = plt.figure(5)
        ax5 = fig5.add_subplot(111)
        pynbody.plot.hist2d(
            np.log10(sphere.gas['rho'].in_units('m_p cm**-3')), np.log10(sphere.gas['temp']), scalemin=1e1, scalemax=2e4, logscale = True,
            weights=sphere.gas['metals']*sphere.gas['mass'].in_units('Msol'),
            xrange=den_range, 
            ret_im=True,clear=False,subplot=ax5)
        ax5.set_xlabel('n [cm$^{-3}$]')
        ax5.set_ylabel('Temperature [K]')
        fig5.savefig(tfile+'.'+halo_num+'.phase_Z.png')

        colorVal = scalarMap.to_rgba(i % 20)
        p = profile.Profile(s.g,ndim=3,min='.01 kpc', max = str(2.5*rvir)+' kpc')
        p_dm = profile.Profile(s.d,ndim=3,min='.01 kpc', max = str(2.5*rvir)+' kpc')
        p_s = profile.Profile(s.s,ndim=3,min='.01 kpc', max = str(2.5*rvir)+' kpc')

        axs[0].plot(p['rbins']/rvir,(p['mass']+p_s['mass'])/(np.sum(p['mass']) + np.sum(p_s['mass'])),color = colorVal)
        #axs[0].plot(p['rbins']/rvir,(p_s['mass'])/(np.sum(p['mass']) + np.sum(p_s['mass'])), 'k',linestyle = ':')
        #axs[0].plot(p['rbins']/rvir,(p['mass'])/(np.sum(p['mass']) + np.sum(p_s['mass'])), 'k',linestyle = '--')
        axs[0].semilogy()
        axs[0].set_xlabel('R [kpc]')
        #axs[0].set_ylabel(r'$\Sigma_{\star}$ [M$_{\odot}$ kpc$^{-2}$]')

        axs[1].plot(p['rbins']/rvir,(p['mass_enc'] + p_s['mass_enc'])/np.max(p['mass_enc']+p_s['mass_enc']),color = colorVal)
        #axs[1].plot(p['rbins']/rvir,p_s['mass_enc']/np.max(p['mass_enc']+p_s['mass_enc']), 'k',linestyle = ':')
        #axs[1].plot(p['rbins']/rvir,p['mass_enc']/np.max(p['mass_enc']+p_s['mass_enc']), 'k',linestyle='--')
        #axs[1].semilogy()
        axs[1].set_xlabel('R [kpc]')
        #axs[1].set_ylabel(r'$\rho_{DM}$ [M$_{\odot}$ kpc$^{-3}$]')

        axs[2].plot(p['rbins']/rvir,(p['mass_enc']+p_s['mass_enc'])/p_dm['mass_enc'],color = colorVal)
        #axs[2].plot(p['rbins']/rvir,(p_s['mass_enc'])/p_dm['mass_enc'], 'k',linestyle = ':')
        #axs[2].plot(p['rbins']/rvir,(p['mass_enc'])/p_dm['mass_enc'], 'k', linestyle='--')
        axs[2].semilogy()
        axs[2].set_xlabel('R [kpc]')
        
        i += 1

    fig.savefig(tfile+'.halo_prof.png')
    fig.clear()

    

if __name__ == '__main__':
    if (socket.gethostname() == "quirm"):
        prefix = '/home/christenc/Data/Sims/'
        prefix_outfile = '/home/christenc/Figures/marvel/'
    else:
        prefix = '/home/christensen/Storage2/UW/MolecH/Cosmo/'
        prefix_outfile = '/home/christensen/Plots/marvel/'

    halo_nums = ['1','2','4','5','6','7','10','11','13','14','27']
    tfile_base = 'cptmarvel.cosmo25cmb.4096g5HbwK1BH'
    tfile = prefix + 'cptmarvel.cosmo25cmb/cptmarvel.cosmo25cmb.4096g5HbwK1BH/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096' #'cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096'
    outfile_base = prefix_outfile + tfile_base
    halo_profile(tfile,halo_nums)

    halo_nums = ['1','3','7','8','10','11','12','16','17','18','30','32','34','36','61','123']
    tfile = prefix + 'rogue.cosmo25cmb/rogue.cosmo25cmb.4096g5HbwK1BH/rogue.cosmo25cmb.4096g5HbwK1BH.004096/rogue.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_base = 'rogue.cosmo25cmb.4096g5HbwK1BH'
    outfile_base = prefix_outfile + tfile_base
    halo_profile(tfile,halo_nums)

    halo_nums = ['1','2','3','4','5','8','9','10','11','12','17','18','37','75']
    tfile = prefix + 'elektra.cosmo25cmb/elektra.cosmo25cmb.4096g5HbwK1BH/elektra.cosmo25cmb.4096g5HbwK1BH.004096/elektra.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_base = 'elektra.cosmo25cmb.4096g5HbwK1BH'
    outfile_base = prefix_outfile + tfile_base
    halo_profile(tfile,halo_nums)

    halo_nums = ['1','2','3','4','5','6','7','8','10','11','13','14','15','16','17','23','24','28','34','35','43','49','50','60','109','124','125','192','208'] 
    tfile = prefix + 'storm.cosmo25cmb/storm.cosmo25cmb.4096g5HbwK1BH/storm.cosmo25cmb.4096g5HbwK1BH.004096/storm.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_base = 'storm.cosmo25cmb.4096g5HbwK1BH'
    outfile_base = prefix_outfile + tfile_base
    halo_profile(tfile,halo_nums)
    
    tfile_base = 'h148.cosmo50PLK.3072g3HbwK1BH'
    tfile = prefix + 'h148.cosmo50PLK.3072g/h148.cosmo50PLK.3072g3HbwK1BH/h148.cosmo50PLK.3072g3HbwK1BH.004096/h148.cosmo50PLK.3072g3HbwK1BH.004096'    
    outfile_base = prefix_outfile + tfile_base
    #halo_profile(tfile,outfile_base,tfile_base,halo_nums)
     
    tfile_base = 'h329.cosmo50PLK.3072gst5HbwK1BH.004096'
    tfile = prefix + 'h329.cosmo50PLK.3072g/h329.cosmo50PLK.3072gst5HbwK1BH/h329.cosmo50PLK.3072gst5HbwK1BH.004096/h329.cosmo50PLK.3072gst5HbwK1BH.004096'    
    outfile_base = prefix_outfile + tfile_base
    #halo_profile(tfile,outfile_base,tfile_base,halo_nums)
    
    tfile_base = 'h229.cosmo50PLK.3072gst5HbwK1BH.004096'
    tfile = prefix + 'h229.cosmo50PLK.3072g/h229.cosmo50PLK.3072gst5HbwK1BH/h229.cosmo50PLK.3072gst5HbwK1BH.004096/h229.cosmo50PLK.3072gst5HbwK1BH.004096'    
    outfile_base = prefix_outfile + tfile_base
    #halo_profile(tfile,outfile_base,tfile_base,halo_nums)
    
    tfile_base = 'h242.cosmo50PLK.3072gst5HbwK1BH.004096'
    tfile = prefix + 'h242.cosmo50PLK.3072g/h242.cosmo50PLK.3072gst5HbwK1BH/h242.cosmo50PLK.3072gst5HbwK1BH.004096/h242.cosmo50PLK.3072gst5HbwK1BH.004096'    
    outfile_base = prefix_outfile + tfile_base
    #halo_profile(tfile,outfile_base,tfile_base,halo_nums)

        
