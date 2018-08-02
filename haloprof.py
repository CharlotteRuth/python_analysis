#6/21/2018
#Charlotte Christensen
#Plot the halo metallicity, temperature, and density of individual galaxies

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.colors as colors
import matplotlib.cm as cm
import matplotlib as mpl
import pynbody, sys
from pynbody.analysis import profile, angmom, halo
from pynbody import units, config
import pynbody.filt as f
import socket
import matplotlib.gridspec as gridspec

from pynbody.analysis import luminosity as lum
import os, glob, pickle
import os.path
import pandas as pd

#Run with
#%run /home/christensen/Code/python/python_analysis/haloprof.py
#or, on quirm,
#%run /home/christenc/Code/python/python_analysis/haloprof.py
#ipython --pylab

def pickle_read(file):

    objs = []
    f=open(file, 'rb')
    while 1:
        try:
            objs.append(pickle.load(f))
        except EOFError:
            break
        
    f.close()

    return pd.DataFrame(objs)


def haloprof_plot(tfile,halo_nums,halotype = 'all',normalize = False):        
    vmax = 12
    min_vmass = 1e8
    
    hfb = pynbody.load(tfile)
    h = hfb.halos()
    if normalize:
        h_dummy = hfb.halos(dummy = True)
    halo_data = pickle_read(tfile + ".data")

    for halo_num in halo_nums:        
        cmx = plt.get_cmap("viridis_r") 
        cNorm  = colors.Normalize(vmin=np.log10(min_vmass), vmax = vmax)
        scalarMap = cm.ScalarMappable(norm=cNorm, cmap=cmx)
        if (os.path.isfile(tfile + "_" + halo_num + "_halo" + ".data") == True):
            p_gas = pickle_read(tfile + "_" + halo_num + "_halo" + ".data")
            if normalize:
                Tvir = (0.59)*1.6726e-27*(h_dummy[int(halo_num)].properties['Vmax']*1000)**2/2/1.38064852e-23 #Fo completely ionized primordial gas
                dens_norm = p_gas['dens'][0][1][-1]
                T_norm = Tvir
                Z_norm = p_gas['metals'][0][1][0]
            else:
                dens_norm = 1
                T_norm = 1
                Z_norm = 1
            colorVal = scalarMap.to_rgba(np.log10(h[int(halo_num)].properties['mass']))
            if (halotype == 'all'):
                axs[0].plot(p_gas['rrvir'][0],p_gas['dens'][0][:,1]/dens_norm, color = colorVal)
                axs[1].plot(p_gas['rrvir'][0],p_gas['temp'][0][:,1]/T_norm, color = colorVal)
                axs[2].plot(p_gas['rrvir'][0],p_gas['metals'][0][:,1]/Z_norm, color = colorVal)
                plt.show()
            else:
                if (len(halo_data[halo_data['haloid'] == halo_num]) > 0) :
                    if (float(halo_data[halo_data['haloid'] == halo_num]['SFR']) > 1e-4) and (halotype == 'sf'):
                        axs[0].plot(p_gas['rrvir'][0],p_gas['dens'][0][:,1]/dens_norm, color = colorVal)
                        axs[1].plot(p_gas['rrvir'][0],p_gas['temp'][0][:,1]/T_norm, color = colorVal)
                        axs[2].plot(p_gas['rrvir'][0],p_gas['metals'][0][:,1]/Z_norm, color = colorVal)
                        plt.show()
                    if (float(halo_data[halo_data['haloid'] == halo_num]['SFR']) < 1e-4) and (halotype == 'quenched'):
                        axs[0].plot(p_gas['rrvir'][0],p_gas['dens'][0][:,1]/dens_norm, color = colorVal)
                        axs[1].plot(p_gas['rrvir'][0],p_gas['temp'][0][:,1]/T_norm, color = colorVal)
                        axs[2].plot(p_gas['rrvir'][0],p_gas['metals'][0][:,1]/Z_norm, color = colorVal)
                        plt.show()

                        
def haloprof_cumplot(tfile,halo_nums,halotype = 'all',normalize = False):
    vmax = 12
    min_vmass = 1e8
    f_bar = 0.16510
    Zyield = 0.02788242
    
    hfb = pynbody.load(tfile)
    h = hfb.halos()
    halo_data = pickle_read(tfile + ".data")

    for halo_num in halo_nums:        
        cmx = plt.get_cmap("viridis_r") 
        cNorm  = colors.Normalize(vmin=np.log10(min_vmass), vmax = vmax)
        scalarMap = cm.ScalarMappable(norm=cNorm, cmap=cmx)
        if (os.path.isfile(tfile + "_" + halo_num + "_halo" + ".data") == True):
            p_gas = pickle_read(tfile + "_" + halo_num + "_halo" + ".data")
            mass_norm = h[int(halo_num)].properties['mass']*f_bar
            Z_norm = np.sum(h[int(halo_num)].star['massform'].in_units('Msol'))*Zyield
            colorVal = scalarMap.to_rgba(np.log10(h[int(halo_num)].properties['mass']))
            if (halotype == 'all'):
                axs[0].plot(p_gas['rrvir'][0],p_gas['mass_enc'][0]/mass_norm, color = colorVal)
                axs[1].plot(p_gas['rrvir'][0],p_gas['metals_enc'][0]/Z_norm, color = colorVal)
                plt.show()
            else:
                if (len(halo_data[halo_data['haloid'] == halo_num]) > 0) :
                    if (float(halo_data[halo_data['haloid'] == halo_num]['SFR']) > 1e-4) and (halotype == 'sf'):
                        axs[0].plot(p_gas['rrvir'][0],p_gas['mass_enc'][0]/mass_norm, color = colorVal)
                        axs[1].plot(p_gas['rrvir'][0],p_gas['metals_enc'][0]/Z_norm, color = colorVal)
                        plt.show()
                    if (float(halo_data[halo_data['haloid'] == halo_num]['SFR']) < 1e-4) and (halotype == 'quenched'):
                        axs[0].plot(p_gas['rrvir'][0],p_gas['mass_enc'][0]/mass_norm, color = colorVal)
                        axs[1].plot(p_gas['rrvir'][0],p_gas['metals_enc'][0]/Z_norm, color = colorVal)
                        plt.show()    

                        
def haloprof(tfile,halo_nums):

    zsolar = 0.0130215
    hfb = pynbody.load(tfile)
    h = hfb.halos()
    plt_width = 3.5
    max_radii_scale = 2.0
    nbins = 100
    
    for halo_num in halo_nums:
        print(tfile+"_" + halo_num)
        halo_num_int = int(halo_num)
        hfb.physical_units()        
        pynbody.analysis.angmom.faceon(h[halo_num_int])
        hfbrvir = np.max(h[halo_num_int].dark['r'])
        halo = hfb[pynbody.filt.Sphere(max_radii_scale*hfbrvir, (0,0,0))]
        
        p_gas = pynbody.analysis.profile.QuantileProfile(halo.g,min=0,max=max_radii_scale*hfbrvir,ndim=3,nbins = nbins, weights = halo.g['mass'])
        p_gas_single = pynbody.analysis.profile.Profile(halo.g,min=0,max=max_radii_scale*hfbrvir,ndim=3,nbins = nbins)
        p_star_single = pynbody.analysis.profile.Profile(halo.s,min=0,max=max_radii_scale*hfbrvir,ndim=3,nbins = nbins)
        halo.g['radii'] = np.sqrt(halo.g['x']**2 + halo.g['y']**2 + halo.g['z']**2)
        halo.s['radii'] = np.sqrt(halo.s['x']**2 + halo.s['y']**2 + halo.s['z']**2)
        metals_enc_gas, loc = np.histogram(halo.g['radii'], bins  = nbins, weights = halo.g['metals']*halo.g['mass'].in_units('Msol'), range = [0,max_radii_scale*hfbrvir])
        metals_enc_stars, loc = np.histogram(halo.s['radii'], bins  = nbins, weights = halo.s['metals']*halo.s['mass'].in_units('Msol'), range = [0,max_radii_scale*hfbrvir])
        metals_enc = np.cumsum(metals_enc_gas) + np.cumsum(metals_enc_stars)

        if (len(halo.g) != 0):
            f, axs = plt.subplots(1,3,figsize=(plt_width*2.,plt_width*6./7.))
            axs[0].plot(p_gas['rbins']/hfbrvir,p_gas['rho'][:,1], 'k')
            axs[0].fill_between(p_gas['rbins']/hfbrvir, p_gas['rho'][:,0], p_gas['rho'][:,2], color = 'Grey', alpha=0.5)
            axs[0].semilogy()
            axs[0].set_xlabel(r'R [kpc]/R$_{vir}$')
            axs[0].set_ylabel(r'$\rho_{gas}$ [M$_{\odot}$ kpc$^{-3}$]')
            axs[0].set_xlim([0,1])
            
            axs[1].plot(p_gas['rbins']/hfbrvir,p_gas['temp'][:,1], 'k')
            axs[1].fill_between(p_gas['rbins']/hfbrvir, p_gas['temp'][:,0], p_gas['temp'][:,2], color = 'Grey', alpha=0.5)
            axs[1].semilogy()
            axs[1].set_xlabel(r'R [kpc]/R$_{vir}$')
            axs[1].set_ylabel(r'Temperature [K]')
            axs[1].set_xlim([0,1])

            axs[2].plot(p_gas['rbins']/hfbrvir,p_gas['metals'][:,1]/zsolar, 'k')
            axs[2].fill_between(p_gas['rbins']/hfbrvir, p_gas['metals'][:,0]/zsolar, p_gas['metals'][:,2]/zsolar, color = 'Grey', alpha=0.5)
            axs[2].semilogy()
            axs[2].set_xlabel(r'R [kpc]/R$_{vir}$')
            axs[2].set_ylabel(r'$Z/Z_{\odot}$')
            axs[2].set_xlim([0,1])
            plt.tight_layout()
            Z_grad = (np.polyfit(p_gas['rbins'][p_gas['rbins'] < hfbrvir]/hfbrvir,np.log10(p_gas['metals'][p_gas['rbins'] < hfbrvir,1]/zsolar),1))[1] #Note that I scale by R_vir here but I should use R_e
            plt.savefig(tfile + "_" + halo_num + "_haloprof.png")
            f.clear()
            f.clf()
            plt.close('all')
        
            pickle_file=open(tfile+"_" + halo_num + "_halo" + ".data","wb")
            pickle.dump({'rrvir':p_gas['rbins']/hfbrvir,
                    'dens':p_gas['rho'],     
                    'temp':p_gas['temp'],
                    'metals':p_gas['metals']/zsolar,
                    'mass_enc':p_gas_single['mass_enc'] + p_star_single['mass_enc'],
                    'metals_enc':metals_enc,
                    'Z_grad':Z_grad
                },pickle_file, pickle.HIGHEST_PROTOCOL)
            pickle_file.close()



if __name__ == '__main__':    
    if (socket.gethostname() == "quirm"):
        prefix = '/home/christenc/Data/Sims/'
        prefix_outfile = '/home/christenc/Figures/marvel/'
    else:
        prefix = '/home/christensen/Storage2/UW/MolecH/Cosmo/'
        prefix_outfile = '/home/christensen/Plots/marvel/'

    halo_nums_cm = ['1','2','4','5','6','7','10','11','13','14','27']
    tfile_cm = prefix + 'cptmarvel.cosmo25cmb/cptmarvel.cosmo25cmb.4096g5HbwK1BH/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_name_cm = 'cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096'
    
    halo_nums_r = ['1','3','7','8','10','11','12','16','17','18','30','32','34','36','61','123']
    tfile_r = prefix + 'rogue.cosmo25cmb/rogue.cosmo25cmb.4096g5HbwK1BH/rogue.cosmo25cmb.4096g5HbwK1BH.004096/rogue.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_name_r = 'rogue.cosmo25cmb.4096g5HbwK1BH.004096'

    tfile_e = prefix + 'elektra.cosmo25cmb/elektra.cosmo25cmb.4096g5HbwK1BH/elektra.cosmo25cmb.4096g5HbwK1BH.004096/elektra.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_name_e = 'elektra.cosmo25cmb.4096g5HbwK1BH.004096'
    halo_nums_e = ['1','2','3','4','5','8','9','10','11','12','17','18','37','75']

    tfile_s = prefix + 'storm.cosmo25cmb/storm.cosmo25cmb.4096g5HbwK1BH/storm.cosmo25cmb.4096g5HbwK1BH.004096/storm.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_name_s = 'storm.cosmo25cmb.4096g5HbwK1BH'
    halo_nums_s = ['1','2','3','4','5','6','7','8','10','11','13','14','15','16','17','23','24','28','34','35','43','49','50','60','109','124','125','192','208']   

    #haloprof(tfile_cm,halo_nums_cm)
    #haloprof(tfile_r,halo_nums_r)
    #haloprof(tfile_e,halo_nums_e)    
    #haloprof(tfile_s,halo_nums_s)
    
#-------------------------------------------
    outfile_base = prefix_outfile  + 'marvel_halos'
    presentation = False
    if presentation:
        outfile_base = outfile_base + '_pres'
        plt.style.use(['default','/home/christenc/.config/matplotlib/presentation.mplstyle'])
        plt_width = 8 #inches
        aspect_ratio = 3.0/4.0
        legendsize = 16
        dpi = 100
    else:
        plt.style.use(['default','/home/christenc/.config/matplotlib/article.mplstyle'])
        plt_width = 3.5 #inches
        aspect_ratio = 3.0/4.0
        legendsize = 5
        dpi = 300

    #All halos, cumulative mass and metals
    f, axs = plt.subplots(1,2,figsize=(plt_width*2.,plt_width*6./7.))
    axs[0].semilogy()
    axs[0].set_xlabel(r'R/R$_{vir}$')
    axs[0].set_ylabel(r'$M/(f_{bary} \times M_{vir})$')
    axs[1].semilogy()
    axs[1].set_xlabel(r'R/R$_{vir}$')
    axs[1].set_ylabel(r'$M_Z/(y\times M_*)$')
    haloprof_cumplot(tfile_cm,halo_nums_cm)
    haloprof_cumplot(tfile_e,halo_nums_e)
    haloprof_cumplot(tfile_s,halo_nums_s)
    haloprof_cumplot(tfile_r,halo_nums_r)
    axs[0].set_ylim([1e-4,3])
    axs[1].set_ylim([1e-3,100])
    plt.tight_layout(w_pad=1.4)
    f.subplots_adjust(right=0.8)
    b_ax = f.add_axes([0.85, f.subplotpars.bottom, 0.02, (f.subplotpars.top - f.subplotpars.bottom)])
    vmax = 12
    min_vmass = 1e8    
    cmx = plt.get_cmap("viridis_r") 
    cNorm  = colors.Normalize(vmin=np.log10(min_vmass), vmax = vmax)
    cb = mpl.colorbar.ColorbarBase(b_ax, cmap=cmx, norm=cNorm)
    cb.set_label('Log Virial Mass [M$_\odot$]')
    plt.savefig(outfile_base + '_baryZprof.png')
    plt.close()
        
    #All halos, not normalized
    f, axs = plt.subplots(1,3,figsize=(plt_width*2.,plt_width*6./7.))
    axs[0].semilogy()
    axs[0].set_xlabel(r'R/R$_{vir}$')
    axs[0].set_ylabel(r'$\rho_{gas}$ [M$_{\odot}$ kpc$^{-3}$]')
    axs[1].semilogy()
    axs[1].set_xlabel(r'R/R$_{vir}$')
    axs[1].set_ylabel(r'Temperature [K]')
    axs[2].semilogy()
    axs[2].set_xlabel(r'R/R$_{vir}$')
    axs[2].set_ylabel(r'$Z/Z_{\odot}$')
    haloprof_plot(tfile_cm,halo_nums_cm)
    haloprof_plot(tfile_e,halo_nums_e)
    haloprof_plot(tfile_s,halo_nums_s)
    haloprof_plot(tfile_r,halo_nums_r)
    axs[0].set_ylim([10,1e6])
    axs[1].set_ylim([1e4,1e5])
    axs[2].set_ylim([1e-4,2])
    axs[0].set_xlim([0,1])
    axs[1].set_xlim([0,1])
    axs[2].set_xlim([0,1])    
    plt.tight_layout(w_pad=1.4)
    f.subplots_adjust(right=0.8)
    b_ax = f.add_axes([0.85, f.subplotpars.bottom, 0.02, (f.subplotpars.top - f.subplotpars.bottom)])
    vmax = 12
    min_vmass = 1e8    
    cmx = plt.get_cmap("viridis_r") 
    cNorm  = colors.Normalize(vmin=np.log10(min_vmass), vmax = vmax)
    cb = mpl.colorbar.ColorbarBase(b_ax, cmap=cmx, norm=cNorm)
    cb.set_label('Log Virial Mass [M$_\odot$]')
    plt.savefig(outfile_base + '.png')
    plt.close()

    #Star forming halos, not normalized
    f, axs = plt.subplots(1,3,figsize=(plt_width*2.,plt_width*6./7.))
    axs[0].semilogy()
    axs[0].set_xlabel(r'R/R$_{vir}$')
    axs[0].set_ylabel(r'$\rho_{gas}$ [M$_{\odot}$ kpc$^{-3}$]')
    axs[1].semilogy()
    axs[1].set_xlabel(r'R/R$_{vir}$')
    axs[1].set_ylabel(r'Temperature [K]')
    axs[2].semilogy()
    axs[2].set_xlabel(r'R/R$_{vir}$')
    axs[2].set_ylabel(r'$Z/Z_{\odot}$')
    haloprof_plot(tfile_cm,halo_nums_cm,halotype = 'sf')
    haloprof_plot(tfile_e,halo_nums_e,halotype = 'sf')
    haloprof_plot(tfile_s,halo_nums_s,halotype = 'sf')
    haloprof_plot(tfile_r,halo_nums_r,halotype = 'sf')
    axs[0].set_ylim([100,1e6])
    axs[1].set_ylim([1e4,1e5])
    axs[2].set_ylim([1e-4,2])
    axs[0].set_xlim([0,1])
    axs[1].set_xlim([0,1])
    axs[2].set_xlim([0,1])     
    plt.tight_layout(w_pad=1.4)
    f.subplots_adjust(right=0.8)
    b_ax = f.add_axes([0.85, f.subplotpars.bottom, 0.02, (f.subplotpars.top - f.subplotpars.bottom)])
    vmax = 12
    min_vmass = 1e8    
    cmx = plt.get_cmap("viridis_r") 
    cNorm  = colors.Normalize(vmin=np.log10(min_vmass), vmax = vmax)
    cb = mpl.colorbar.ColorbarBase(b_ax, cmap=cmx, norm=cNorm)
    cb.set_label('Log Virial Mass [M$_\odot$]')
    plt.savefig(outfile_base + '.sfgal.png')
    plt.close()

    #Quenched halos, not normalized
    f, axs = plt.subplots(1,3,figsize=(plt_width*2.,plt_width*6./7.))
    axs[0].semilogy()
    axs[0].set_xlabel(r'R/R$_{vir}$')
    axs[0].set_ylabel(r'$\rho_{gas}$ [M$_{\odot}$ kpc$^{-3}$]')
    axs[1].semilogy()
    axs[1].set_xlabel(r'R/R$_{vir}$')
    axs[1].set_ylabel(r'Temperature [K]')
    axs[2].semilogy()
    axs[2].set_xlabel(r'R/R$_{vir}$')
    axs[2].set_ylabel(r'$Z/Z_{\odot}$')
    haloprof_plot(tfile_cm,halo_nums_cm,halotype = 'quenched')
    haloprof_plot(tfile_e,halo_nums_e,halotype = 'quenched')
    haloprof_plot(tfile_s,halo_nums_s,halotype = 'quenched')
    haloprof_plot(tfile_r,halo_nums_r,halotype = 'quenched')
    axs[0].set_ylim([10,1e6])
    axs[1].set_ylim([1e4,1e5])
    axs[2].set_ylim([1e-4,0.2])
    axs[0].set_xlim([0,1])
    axs[1].set_xlim([0,1])
    axs[2].set_xlim([0,1])     
    plt.tight_layout(w_pad=1.4)
    f.subplots_adjust(right=0.8)
    b_ax = f.add_axes([0.85, f.subplotpars.bottom, 0.02, (f.subplotpars.top - f.subplotpars.bottom)])
    vmax = 12
    min_vmass = 1e8    
    cmx = plt.get_cmap("viridis_r") 
    cNorm  = colors.Normalize(vmin=np.log10(min_vmass), vmax = vmax)
    cb = mpl.colorbar.ColorbarBase(b_ax, cmap=cmx, norm=cNorm)
    cb.set_label('Log Virial Mass [M$_\odot$]')
    plt.savefig(outfile_base + '.quenched.png')
    plt.close()    

    #All halos, normalized
    f, axs = plt.subplots(1,3,figsize=(plt_width*2.,plt_width*6./7.))
    axs[0].semilogy()
    axs[0].set_xlabel(r'R/R$_{vir}$')
    axs[0].set_ylabel(r'$\rho_{gas}/\rho_{vir}$')
    axs[1].semilogy()
    axs[1].set_xlabel(r'R/R$_{vir}$')
    axs[1].set_ylabel(r'T/T$_{vir}$')
    axs[2].semilogy()
    axs[2].set_xlabel(r'R/R$_{vir}$')
    axs[2].set_ylabel(r'$Z/Z_{cen}$')
    haloprof_plot(tfile_cm,halo_nums_cm, normalize = True)
    haloprof_plot(tfile_e,halo_nums_e, normalize = True)
    haloprof_plot(tfile_s,halo_nums_s, normalize = True)
    haloprof_plot(tfile_r,halo_nums_r, normalize = True)
    axs[0].set_ylim([1e-6,1e2])
    axs[1].set_ylim([1e-2,1e2])
    axs[2].set_ylim([1e-4,2])
    axs[0].set_xlim([0,1])
    axs[1].set_xlim([0,1])
    axs[2].set_xlim([0,1])     
    plt.tight_layout(w_pad=1.4)
    plt.subplots_adjust(right=0.8)
    b_ax = f.add_axes([0.85, f.subplotpars.bottom, 0.02, (f.subplotpars.top - f.subplotpars.bottom)]) 
    vmax = 12
    min_vmass = 1e8
    cmx = plt.get_cmap("viridis_r") 
    cNorm  = colors.Normalize(vmin=np.log10(min_vmass), vmax = vmax)
    cb = mpl.colorbar.ColorbarBase(b_ax, cmap=cmx, norm=cNorm)
    cb.set_label('Log Virial Mass [M$_\odot$]')
    plt.savefig(outfile_base + '_norm.png')
    plt.close()

    #All halos, star forming
    f, axs = plt.subplots(1,3,figsize=(plt_width*2.,plt_width*6./7.))
    axs[0].semilogy()
    axs[0].set_xlabel(r'R/R$_{vir}$')
    axs[0].set_ylabel(r'$\rho_{gas}/\rho_{vir}$')
    axs[1].semilogy()
    axs[1].set_xlabel(r'R/R$_{vir}$')
    axs[1].set_ylabel(r'T/T$_{vir}$')
    axs[2].semilogy()
    axs[2].set_xlabel(r'R/R$_{vir}$')
    axs[2].set_ylabel(r'$Z/Z_{cen}$')
    haloprof_plot(tfile_cm,halo_nums_cm, normalize = True, halotype = 'sf')
    haloprof_plot(tfile_e,halo_nums_e, normalize = True, halotype = 'sf')
    haloprof_plot(tfile_s,halo_nums_s, normalize = True, halotype = 'sf')
    haloprof_plot(tfile_r,halo_nums_r, normalize = True, halotype = 'sf')
    axs[0].set_ylim([1e-6,1e2])
    axs[1].set_ylim([1e-2,1e2])
    axs[2].set_ylim([1e-4,2])
    axs[0].set_xlim([0,1])
    axs[1].set_xlim([0,1])
    axs[2].set_xlim([0,1])     
    plt.tight_layout(w_pad=1.4)
    plt.subplots_adjust(right=0.8)
    b_ax = f.add_axes([0.85, f.subplotpars.bottom, 0.02, (f.subplotpars.top - f.subplotpars.bottom)]) 
    vmax = 12
    min_vmass = 1e8
    cmx = plt.get_cmap("viridis_r") 
    cNorm  = colors.Normalize(vmin=np.log10(min_vmass), vmax = vmax)
    cb = mpl.colorbar.ColorbarBase(b_ax, cmap=cmx, norm=cNorm)
    cb.set_label('Log Virial Mass [M$_\odot$]')
    plt.savefig(outfile_base + '_norm.sfgal.png')
    plt.close()
    
    #All halos, quenched
    f, axs = plt.subplots(1,3,figsize=(plt_width*2.,plt_width*6./7.))
    axs[0].semilogy()
    axs[0].set_xlabel(r'R/R$_{vir}$')
    axs[0].set_ylabel(r'$\rho_{gas}/\rho_{vir}$')
    axs[1].semilogy()
    axs[1].set_xlabel(r'R/R$_{vir}$')
    axs[1].set_ylabel(r'T/T$_{vir}$')
    axs[2].semilogy()
    axs[2].set_xlabel(r'R/R$_{vir}$')
    axs[2].set_ylabel(r'$Z/Z_{cen}$')
    haloprof_plot(tfile_cm,halo_nums_cm, normalize = True, halotype = 'quenched')
    haloprof_plot(tfile_e,halo_nums_e, normalize = True, halotype = 'quenched')
    haloprof_plot(tfile_s,halo_nums_s, normalize = True, halotype = 'quenched')
    haloprof_plot(tfile_r,halo_nums_r, normalize = True, halotype = 'quenched')
    axs[0].set_ylim([1e-6,1e2])
    axs[1].set_ylim([1e-2,1e2])
    axs[2].set_ylim([1e-4,2])
    axs[0].set_xlim([0,1])
    axs[1].set_xlim([0,1])
    axs[2].set_xlim([0,1])     
    plt.tight_layout(w_pad=1.4)
    plt.subplots_adjust(right=0.8)
    b_ax = f.add_axes([0.85, f.subplotpars.bottom, 0.02, (f.subplotpars.top - f.subplotpars.bottom)]) 
    vmax = 12
    min_vmass = 1e8
    cmx = plt.get_cmap("viridis_r") 
    cNorm  = colors.Normalize(vmin=np.log10(min_vmass), vmax = vmax)
    cb = mpl.colorbar.ColorbarBase(b_ax, cmap=cmx, norm=cNorm)
    cb.set_label('Log Virial Mass [M$_\odot$]')
    plt.savefig(outfile_base + '_norm.quenched.png')
    plt.close()
    
