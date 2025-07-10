import numpy as np
import matplotlib.pyplot as plt

import matplotlib
import pynbody, sys
from pynbody.analysis import profile, angmom, halo
from pynbody import units, config
import pynbody.filt as filt
import socket
import pandas as pd

from pynbody.analysis import luminosity as lum
import os, glob

# Run with the conda_env_py3 environment. threading for rogue breaks down with conda_env_py38, for an unknown reason


#Run with
#%run /home/christensen/Code/python/python_analysis/plotColumnDen.py
#or, on quirm,
#%run /home/christenc/Code/python/python_analysis/plotColumnDen.py
#ipython --pylab

def make_rs(im, width): # width should be width of box in kpc
    """ Create radial grid. Each pixel in im will be assigned a radial distance """
    xsize, ysize = np.shape(im)
    x = np.arange(-xsize/2, xsize/2)/float(xsize)*width
    y = np.arange(-ysize/2, ysize/2)/float(xsize)*width
    xs, ys = np.meshgrid(x,y)
    return np.sqrt(xs**2 + ys**2)

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

def createSpectra(tfile,halo_num,bordoloidata,burchettdata):
    presentation = True
    if presentation:
        plt.style.use(['default','/home/christenc/.config/matplotlib/presentation.mplstyle'])
        plt_width = 8.5 #inches
        aspect_ratio = 3.0/4.0
        legendsize = 16
        dpi = 100
    else:
        plt.style.use(['default','/home/christenc/.config/matplotlib/article.mplstyle'])
        plt_width = 3.5 #inches
        aspect_ratio = 3.0/4.0
        legendsize = 5
        dpi = 300

    pynbody.config['number_of_threads'] = 1
        
#tfile = sys.argv[1]
    Carbon   = 0.00213
    Oxygen   = 0.00541
    Silicon  = 0.00067
    Iron     = 0.00117

    hfb_all = pynbody.load(tfile)
    h = hfb_all.halos()

    halo_num = int(halo_num)
       
    hfbsmass = np.sum(h[halo_num].stars['mass'].in_units('Msol'))
    hfbvmass = np.sum(h[halo_num]['mass'].in_units('Msol'))
    hfb_all.physical_units()
    pynbody.analysis.angmom.edgeon(h[halo_num], move_all = True)
            
    hfbrvir = np.max(h[halo_num].gas['r'])
    hfb_all['velmag'] = np.sqrt(hfb_all['vx']**2 + hfb_all['vy']**2 + hfb_all['vz']**2)
    #velf = filt.LowPass('velmag',h[halo_num].properties['v_esc'])
    velf = filt.BandPass('vz',-1*h[halo_num].properties['v_esc'],h[halo_num].properties['v_esc'])
    hfb = hfb_all[velf]
    spafx = filt.BandPass('x',-1*np.round(hfbrvir/10)*10*2.5, np.round(hfbrvir/10)*10*2.5)
    spafy = filt.BandPass('y',-1*np.round(hfbrvir/10)*10*2.5, np.round(hfbrvir/10)*10*2.5)
    spafz = filt.BandPass('z',-1*np.round(hfbrvir/10)*10*2.5, np.round(hfbrvir/10)*10*2.5)
    hfb = hfb[spafx & spafy]# & spafz]
    hfb.properties['boxsize'] = pynbody.units.Unit("1000 Mpc")
    notdiskf = filt.Not(filt.Disc('8 kpc','2 kpc'))

    hfb.gas['hden'] = hfb.gas['rho']*hfb.gas['hydrogen']
    hfb.gas['hden'].units = hfb.gas['rho'].units    
    
    vres = 10 # km/s
    vmax = np.ceil(h[halo_num].properties['v_esc']/vres)*vres
    #numvbins = np.ceil(h[halo_num].properties['v_esc']/vres)*2
    vbin_start = np.arange(-1*vmax, vmax, vres)
    spect_stack = []
    for vbin_min in vbin_start:
        velf = filt.BandPass('vz',vbin_min,vbin_min + vres)
        hfb_bin = hfb[velf & notdiskf]
        print('Gas between ' + str(vbin_min) + ' and ' + str(vbin_min + vres) + ': ', len(hfb_bin.gas))
        plt.close()
        hfbhim = pynbody.plot.image(hfb_bin.gas, qty='hden', width=np.round(hfbrvir/10)*10*2.5, units='m_p cm**-2',log=True) #, num_threads = 1)
        if spect_stack == []:
            spect_stack = hfbhim
        else:
            spect_stack = np.dstack((spect_stack,np.array(hfbhim)))
        #wait = input('continue')

    plt.clf()
    fig1 = plt.figure(1)
    ax = fig1.add_subplot(1,1,1)
    ax.plot(vbin_start + vres/2,spect_stack[250,250,:])
    
    plt.close()
    hfbhim = pynbody.plot.image(hfb.gas[notdiskf], qty='hden', width=np.round(hfbrvir/10)*10*2.5, units='m_p cm**-2',log=True) #, num_threads = 1)
    
    hfb.gas['hiif'] = pynbody.analysis.ionfrac.calculate(hfb.gas,ion='hi')
    hfb.gas['hiden'] = hfb.gas['rho']*hfb.gas['hydrogen']*hfb.gas['hiif']
    hfb.gas['hiden'].units = hfb.gas['rho'].units
    plt.close()
    #hfbhiim = pynbody.plot.image(hfb.gas[notdiskf],qty='hiden', clear=False,
    #    units='m_p cm^-2', width=np.round(hfbrvir/10)*10*2.5, show_cbar=False, 
    #    vmin=1e12,vmax=1e15, log=True, threaded=False)
    hfbhiim = pynbody.plot.image(hfb.gas[notdiskf], qty='hiden', width=np.round(hfbrvir/10)*10*2.5, units='m_p cm**-2',log=True) #, threaded=False)

    hfb.gas['cden'] = hfb.gas['rho']*hfb.gas['OxMassFrac']*Carbon/Oxygen
    hfb.gas['cden'].units = hfb.gas['rho'].units
    plt.close()
    hfbcim = pynbody.plot.image(hfb.gas[notdiskf], qty='cden', width=np.round(hfbrvir/10)*10*2.5, units='m_p cm**-2',log=True) #, threaded=False)    
    
    hfb.gas['civif'] = pynbody.analysis.ionfrac.calculate(hfb.gas,ion='civ')
    hfb.gas['civden'] = hfb.gas['rho']*hfb.gas['OxMassFrac']*Carbon/Oxygen*hfb.gas['civif']
    hfb.gas['civden'].units = hfb.gas['rho'].units
    plt.close()
    #hfbcivim = pynbody.plot.image(hfb.gas[notdiskf],qty='civden', clear=False,
    #    units='12 m_p cm^-2', width=np.round(hfbrvir/10)*10*2.5, show_cbar=False, 
    #    vmin=1e12,vmax=1e15, log=True)
    hfbcivim = pynbody.plot.image(hfb.gas[notdiskf],qty='civden',width=np.round(hfbrvir/10)*10*2.5, units='m_p cm**-2',log=True) #, threaded=False)

    hfb.gas['siivif'] = pynbody.analysis.ionfrac.calculate(hfb.gas,ion='siiv')
    hfb.gas['siivden'] = hfb.gas['rho']*hfb.gas['OxMassFrac']*Silicon/Oxygen*hfb.gas['siivif']
    hfb.gas['siivden'].units = hfb.gas['rho'].units
    plt.close()
    #hfbcivim = pynbody.plot.image(hfb.gas[notdiskf],qty='civden', clear=False,
    #    units='12 m_p cm^-2', width=np.round(hfbrvir/10)*10*2.5, show_cbar=False, 
    #    vmin=1e12,vmax=1e15, log=True)
    hfbsiivim = pynbody.plot.image(hfb.gas[notdiskf],qty='siivden',width=np.round(hfbrvir/10)*10*2.5, units='m_p cm**-2',log=True) #, threaded=False)

    hfb.gas['oviif'] = pynbody.analysis.ionfrac.calculate(hfb.gas,ion='ovi')
    hfb.gas['oviden'] = hfb.gas['rho']*hfb.gas['OxMassFrac']*hfb.gas['oviif']
    hfb.gas['oviden'].units = hfb.gas['rho'].units
    plt.close()
    #hfbcivim = pynbody.plot.image(hfb.gas[notdiskf],qty='civden', clear=False,
    #    units='12 m_p cm^-2', width=np.round(hfbrvir/10)*10*2.5, show_cbar=False, 
    #    vmin=1e12,vmax=1e15, log=True)
    hfboviim = pynbody.plot.image(hfb.gas[notdiskf],qty='oviden',width=np.round(hfbrvir/10)*10*2.5, units='m_p cm**-2',log=True) #, threaded=False)
    
    # ifs = np.load("/home/christenc/Code/python/pynbody/pynbody/analysis/ionfracs.npz"
    # ovi?
    
    rs = make_rs(hfbcivim, np.round(hfbrvir/10)*10*2.5)/hfbrvir
    #rad_axes = np.log10(rs.flatten() + 1e-6)
    rad_axes = (rs.flatten() + 1e-6)

    radius_range = [0,2.5/2]
    H_range = [18, 20]
    C_range = [13, 17]
    HI_range = [12.5,17]
    CIV_range = [11.8,15]
    SIIV_range = [11, 15]
    OVI_range = [12.5, 14.5]

    # Plot Hydrogen vs radius
    plt.close()    
    plt.figure(1,figsize=(plt_width,plt_width*aspect_ratio))
    h_rad = pynbody.plot.hist2d(rad_axes, np.log10(hfbhim.flatten()), cmap=plt.cm.gist_yarg)
    xaxis, data = bin_plt(rad_axes, np.log10(hfbhim.flatten()), xmin = 0, xmax = radius_range[1], nbins = 20)
    plot(xaxis, data[:,2], c='c')
    plot(xaxis, data[:,1], c='c', linestyle = "--", linewidth=1)
    plot(xaxis, data[:,3], c='c', linestyle = "--", linewidth=1)
    plot(xaxis, data[:,0], c='c', linestyle = ":", linewidth=1)
    plot(xaxis, data[:,4], c='c', linestyle = ":", linewidth=1)        

    plt.plot([0,0],H_range,color = 'k')
    plt.axis([radius_range[0],radius_range[1],H_range[0],H_range[1]])
    plt.xlabel(r'$\rho/R_{vir}$')
    plt.ylabel(r'log N [H] [cm$^{-2}]$')
    plt.text(0.6,19.75,r'Log($M_*/M_\odot$) = %.1f'%(np.log10(hfbsmass)),size = legendsize,color = 'k')
    plt.tight_layout()
    plt.savefig(tfile+'.'+str(halo_num)+'.Hprof.png',dpi = dpi)

    # Plot Carbon vs radius
    plt.close()    
    plt.figure(1,figsize=(plt_width,plt_width*aspect_ratio))
    c_rad = pynbody.plot.hist2d(rad_axes, np.log10(hfbcim.flatten()), cmap=plt.cm.gist_yarg)
    xaxis, data = bin_plt(rad_axes, np.log10(hfbcim.flatten()), xmin = 0, xmax = radius_range[1], nbins = 20)
    plot(xaxis, data[:,2], c='c')
    plot(xaxis, data[:,1], c='c', linestyle = "--", linewidth=1)
    plot(xaxis, data[:,3], c='c', linestyle = "--", linewidth=1)
    plot(xaxis, data[:,0], c='c', linestyle = ":", linewidth=1)
    plot(xaxis, data[:,4], c='c', linestyle = ":", linewidth=1)        

    plt.plot([0,0],C_range,color = 'k')
    plt.axis([radius_range[0],radius_range[1],C_range[0],C_range[1]])
    plt.xlabel(r'$\rho/R_{vir}$')
    plt.ylabel(r'log N [C] [cm$^{-2}]$')
    plt.text(0.6,16.5,r'Log($M_*/M_\odot$) = %.1f'%(np.log10(hfbsmass)),size = legendsize,color = 'k')
    plt.tight_layout()
    plt.savefig(tfile+'.'+str(halo_num)+'.Cprof.png',dpi = dpi)    
    
    # Plot HI column density vs radius
    plt.close()
    plt.figure(1,figsize=(plt_width,plt_width*aspect_ratio))
    hi_rad = pynbody.plot.hist2d(rad_axes, np.log10(hfbhiim.flatten()), cmap=plt.cm.gist_yarg)
    xaxis, data = bin_plt(rad_axes, np.log10(hfbhiim.flatten()), xmin = 0, xmax = radius_range[1], nbins = 20)
    plot(xaxis, data[:,2], c='c')
    plot(xaxis, data[:,1], c='c', linestyle = "--", linewidth=1)
    plot(xaxis, data[:,3], c='c', linestyle = "--", linewidth=1)
    plot(xaxis, data[:,0], c='c', linestyle = ":", linewidth=1)
    plot(xaxis, data[:,4], c='c', linestyle = ":", linewidth=1)    
    #bordoloi_plt = plt.errorbar(bordoloidata[bordoloidata['logN_limit'] == 'no']['R']/bordoloidata[bordoloidata['logN_limit'] == 'no']['Rvir'],bordoloidata[bordoloidata['logN_limit'] == 'no']['logN'],yerr=bordoloidata[bordoloidata['logN_limit'] == 'no']['logN_err'],color = 'k',fmt = 'o')
    #plt.errorbar(bordoloidata[bordoloidata['logN_limit'] == 'upper']['R']/bordoloidata[bordoloidata['logN_limit'] == 'upper']['Rvir'],bordoloidata[bordoloidata['logN_limit'] == 'upper']['logN'],xerr=None,yerr=.05,fmt='_',color = 'k',ecolor='k',uplims=True) #,capsize=3,mew=0)
    #plt.errorbar(bordoloidata[bordoloidata['logN_limit'] == 'lower']['R']/bordoloidata[bordoloidata['logN_limit'] == 'lower']['Rvir'],bordoloidata[bordoloidata['logN_limit'] == 'lower']['logN'],xerr=None,yerr=.05,fmt='_',color = 'k',ecolor='k',lolims=True) #,capsize=3,mew=0)
    burchett_plt = plt.errorbar(burchettdata[burchettdata['flag_HI'] == 1]['Rrvir'],burchettdata[burchettdata['flag_HI'] == 1]['logN_HI'],yerr=burchettdata[burchettdata['flag_HI'] == 1]['logN_HI_err'],color = 'grey',fmt= 'o')
    plt.errorbar(burchettdata[burchettdata['flag_HI'] == 3]['Rrvir'],burchettdata[burchettdata['flag_HI'] == 3]['logN_HI'],xerr=None,yerr=.05,fmt='_',color = 'grey',ecolor='grey',uplims=True)
    plt.errorbar(burchettdata[burchettdata['flag_HI'] == 2]['Rrvir'],burchettdata[burchettdata['flag_HI'] == 2]['logN_HI'],xerr=None,yerr=.05,fmt='_',color = 'grey',ecolor='grey',lolims=True)    
    
    deltaLogMstar = 0.2
    #eqMstar_bor = np.all([np.array(bordoloidata['logMstar'] <= np.log10(hfbsmass) + deltaLogMstar),np.array(bordoloidata['logMstar'] >= np.log10(hfbsmass) - deltaLogMstar)],axis=0)
    eqMstar_bur = np.all([np.array(burchettdata['logMstar'] <= np.log10(hfbsmass) + deltaLogMstar),np.array(burchettdata['logMstar'] >= np.log10(hfbsmass) - deltaLogMstar)],axis=0)    
    #plt.errorbar(bordoloidata[eqMstar_bor][bordoloidata[eqMstar_bor]['logN_limit'] == 'no']['R']/bordoloidata[eqMstar_bor][bordoloidata[eqMstar_bor]['logN_limit'] == 'no']['Rvir'],bordoloidata[eqMstar_bor][bordoloidata[eqMstar_bor]['logN_limit'] == 'no']['logN'],yerr=bordoloidata[eqMstar_bor][bordoloidata[eqMstar_bor]['logN_limit'] == 'no']['logN_err'],color = 'r',fmt = 'o')
    #plt.errorbar(bordoloidata[eqMstar_bor][bordoloidata[eqMstar_bor]['logN_limit'] == 'upper']['R']/bordoloidata[eqMstar_bor][bordoloidata[eqMstar_bor]['logN_limit'] == 'upper']['Rvir'],bordoloidata[eqMstar_bor][bordoloidata[eqMstar_bor]['logN_limit'] == 'upper']['logN'],xerr=None,yerr=.05,fmt='_',color = 'r',ecolor='r',uplims=True) #,capsize=3,mew=0)
    #plt.errorbar(bordoloidata[eqMstar_bor][bordoloidata[eqMstar_bor]['logN_limit'] == 'lower']['R']/bordoloidata[eqMstar_bor][bordoloidata[eqMstar_bor]['logN_limit'] == 'lower']['Rvir'],bordoloidata[eqMstar_bor][bordoloidata[eqMstar_bor]['logN_limit'] == 'lower']['logN'],xerr=None,yerr=.05,fmt='_',color = 'r',ecolor='r',lolims=True) #,capsize=3,mew=0)
    plt.errorbar(burchettdata[eqMstar_bur][burchettdata[eqMstar_bur]['flag_HI'] == 1]['Rrvir'],burchettdata[eqMstar_bur][burchettdata[eqMstar_bur]['flag_HI'] == 1]['logN_HI'],yerr=burchettdata[eqMstar_bur][burchettdata[eqMstar_bur]['flag_HI'] == 1]['logN_HI_err'],color = 'magenta',fmt= 'o')
    plt.errorbar(burchettdata[eqMstar_bur][burchettdata[eqMstar_bur]['flag_HI'] == 3]['Rrvir'],burchettdata[eqMstar_bur][burchettdata[eqMstar_bur]['flag_HI'] == 3]['logN_HI'],xerr=None,yerr=.05,fmt='_',color = 'magenta',ecolor='magenta',uplims=True)
    plt.errorbar(burchettdata[eqMstar_bur][burchettdata[eqMstar_bur]['flag_HI'] == 2]['Rrvir'],burchettdata[eqMstar_bur][burchettdata[eqMstar_bur]['flag_HI'] == 2]['logN_HI'],xerr=None,yerr=.05,fmt='_',color = 'magenta',ecolor='magenta',lolims=True)    
    
    plt.plot([0,0],HI_range,color = 'k')
    plt.axis([radius_range[0],radius_range[1],HI_range[0],HI_range[1]])
    plt.xlabel(r'$\rho/R_{vir}$')
    plt.ylabel(r'log N [HI] [cm$^{-2}]$')
    plt.text(0.8,16.5,r'Log($M_*/M_\odot$) = %.1f'%(np.log10(hfbsmass)),size = legendsize,color = 'k')
    plt.tight_layout()
    plt.savefig(tfile+'.'+str(halo_num)+'.HIprof.png',dpi = dpi)

    
    plt.close()    
    plt.figure(1,figsize=(plt_width,plt_width*aspect_ratio))
    civ_rad = pynbody.plot.hist2d(rad_axes, np.log10(hfbcivim.flatten()), cmap=plt.cm.gist_yarg)
    xaxis, data = bin_plt(rad_axes, np.log10(hfbcivim.flatten()), xmin = 0, xmax = radius_range[1], nbins = 20)
    plot(xaxis, data[:,2], c='c')
    plot(xaxis, data[:,1], c='c', linestyle = "--", linewidth=1)
    plot(xaxis, data[:,3], c='c', linestyle = "--", linewidth=1)
    plot(xaxis, data[:,0], c='c', linestyle = ":", linewidth=1)
    plot(xaxis, data[:,4], c='c', linestyle = ":", linewidth=1)        
    bordoloi_plt = plt.errorbar(bordoloidata[bordoloidata['logN_limit'] == 'no']['R']/bordoloidata[bordoloidata['logN_limit'] == 'no']['Rvir'],bordoloidata[bordoloidata['logN_limit'] == 'no']['logN'],yerr=bordoloidata[bordoloidata['logN_limit'] == 'no']['logN_err'],color = 'k',fmt = 'o')
    plt.errorbar(bordoloidata[bordoloidata['logN_limit'] == 'upper']['R']/bordoloidata[bordoloidata['logN_limit'] == 'upper']['Rvir'],bordoloidata[bordoloidata['logN_limit'] == 'upper']['logN'],xerr=None,yerr=.05,fmt='_',color = 'k',ecolor='k',uplims=True) #,capsize=3,mew=0)
    plt.errorbar(bordoloidata[bordoloidata['logN_limit'] == 'lower']['R']/bordoloidata[bordoloidata['logN_limit'] == 'lower']['Rvir'],bordoloidata[bordoloidata['logN_limit'] == 'lower']['logN'],xerr=None,yerr=.05,fmt='_',color = 'k',ecolor='k',lolims=True) #,capsize=3,mew=0)
    burchett_plt = plt.errorbar(burchettdata[burchettdata['flag_CIV'] == 1]['Rrvir'],burchettdata[burchettdata['flag_CIV'] == 1]['logN_CIV'],yerr=burchettdata[burchettdata['flag_CIV'] == 1]['logN_CIV_err'],color = 'grey',fmt= 'o')
    plt.errorbar(burchettdata[burchettdata['flag_CIV'] == 3]['Rrvir'],burchettdata[burchettdata['flag_CIV'] == 3]['logN_CIV'],xerr=None,yerr=.05,fmt='_',color = 'grey',ecolor='grey',uplims=True)
    plt.errorbar(burchettdata[burchettdata['flag_CIV'] == 2]['Rrvir'],burchettdata[burchettdata['flag_CIV'] == 2]['logN_CIV'],xerr=None,yerr=.05,fmt='_',color = 'grey',ecolor='grey',lolims=True)    
    
    eqMstar_bor = np.all([np.array(bordoloidata['logMstar'] <= np.log10(hfbsmass) + deltaLogMstar),np.array(bordoloidata['logMstar'] >= np.log10(hfbsmass) - deltaLogMstar)],axis=0)
    eqMstar_bur = np.all([np.array(burchettdata['logMstar'] <= np.log10(hfbsmass) + deltaLogMstar),np.array(burchettdata['logMstar'] >= np.log10(hfbsmass) - deltaLogMstar)],axis=0)    
    plt.errorbar(bordoloidata[eqMstar_bor][bordoloidata[eqMstar_bor]['logN_limit'] == 'no']['R']/bordoloidata[eqMstar_bor][bordoloidata[eqMstar_bor]['logN_limit'] == 'no']['Rvir'],bordoloidata[eqMstar_bor][bordoloidata[eqMstar_bor]['logN_limit'] == 'no']['logN'],yerr=bordoloidata[eqMstar_bor][bordoloidata[eqMstar_bor]['logN_limit'] == 'no']['logN_err'],color = 'r',fmt = 'o')
    plt.errorbar(bordoloidata[eqMstar_bor][bordoloidata[eqMstar_bor]['logN_limit'] == 'upper']['R']/bordoloidata[eqMstar_bor][bordoloidata[eqMstar_bor]['logN_limit'] == 'upper']['Rvir'],bordoloidata[eqMstar_bor][bordoloidata[eqMstar_bor]['logN_limit'] == 'upper']['logN'],xerr=None,yerr=.05,fmt='_',color = 'r',ecolor='r',uplims=True) #,capsize=3,mew=0)
    plt.errorbar(bordoloidata[eqMstar_bor][bordoloidata[eqMstar_bor]['logN_limit'] == 'lower']['R']/bordoloidata[eqMstar_bor][bordoloidata[eqMstar_bor]['logN_limit'] == 'lower']['Rvir'],bordoloidata[eqMstar_bor][bordoloidata[eqMstar_bor]['logN_limit'] == 'lower']['logN'],xerr=None,yerr=.05,fmt='_',color = 'r',ecolor='r',lolims=True) #,capsize=3,mew=0)
    plt.errorbar(burchettdata[eqMstar_bur][burchettdata[eqMstar_bur]['flag_CIV'] == 1]['Rrvir'],burchettdata[eqMstar_bur][burchettdata[eqMstar_bur]['flag_CIV'] == 1]['logN_CIV'],yerr=burchettdata[eqMstar_bur][burchettdata[eqMstar_bur]['flag_CIV'] == 1]['logN_CIV_err'],color = 'magenta',fmt= 'o')
    plt.errorbar(burchettdata[eqMstar_bur][burchettdata[eqMstar_bur]['flag_CIV'] == 3]['Rrvir'],burchettdata[eqMstar_bur][burchettdata[eqMstar_bur]['flag_CIV'] == 3]['logN_CIV'],xerr=None,yerr=.05,fmt='_',color = 'magenta',ecolor='magenta',uplims=True)
    plt.errorbar(burchettdata[eqMstar_bur][burchettdata[eqMstar_bur]['flag_CIV'] == 2]['Rrvir'],burchettdata[eqMstar_bur][burchettdata[eqMstar_bur]['flag_CIV'] == 2]['logN_CIV'],xerr=None,yerr=.05,fmt='_',color = 'magenta',ecolor='magenta',lolims=True)    
    
    plt.plot([0,0],CIV_range,color = 'k')
    plt.axis([radius_range[0],radius_range[1],CIV_range[0],CIV_range[1]])
    plt.xlabel(r'$\rho/R_{vir}$')
    plt.ylabel(r'log N [CIV] [cm$^{-2}]$')
    plt.text(0.6,14.5,r'Log($M_*/M_\odot$) = %.1f'%(np.log10(hfbsmass)),size = legendsize,color = 'k')
    plt.tight_layout()
    plt.savefig(tfile+'.'+str(halo_num)+'.CIVprof.png',dpi = dpi)

    ## SI IV
    plt.close()    
    plt.figure(1,figsize=(plt_width,plt_width*aspect_ratio))
    siiv_rad = pynbody.plot.hist2d(rad_axes, np.log10(hfbsiivim.flatten()), cmap=plt.cm.gist_yarg)
    xaxis, data = bin_plt(rad_axes, np.log10(hfbsiivim.flatten()), xmin = 0, xmax = radius_range[1], nbins = 20)
    plot(xaxis, data[:,2], c='c')
    plot(xaxis, data[:,1], c='c', linestyle = "--", linewidth=1)
    plot(xaxis, data[:,3], c='c', linestyle = "--", linewidth=1)
    plot(xaxis, data[:,0], c='c', linestyle = ":", linewidth=1)
    plot(xaxis, data[:,4], c='c', linestyle = ":", linewidth=1)        

    plt.plot([0,0],SIIV_range,color = 'k')
    plt.axis([radius_range[0],radius_range[1],SIIV_range[0],SIIV_range[1]])
    plt.xlabel(r'$\rho/R_{vir}$')
    plt.ylabel(r'log N [SIIV] [cm$^{-2}]$')
    plt.text(0.6,14.5,r'Log($M_*/M_\odot$) = %.1f'%(np.log10(hfbsmass)),size = legendsize,color = 'k')
    plt.tight_layout()
    plt.savefig(tfile+'.'+str(halo_num)+'.SIIVprof.png',dpi = dpi)

    ## O VI
    plt.close()    
    plt.figure(1,figsize=(plt_width,plt_width*aspect_ratio))
    ovi_rad = pynbody.plot.hist2d(rad_axes, np.log10(hfboviim.flatten()), cmap=plt.cm.gist_yarg)
    xaxis, data = bin_plt(rad_axes, np.log10(hfboviim.flatten()), xmin = 0, xmax = radius_range[1], nbins = 20)
    plot(xaxis, data[:,2], c='c')
    plot(xaxis, data[:,1], c='c', linestyle = "--", linewidth=1)
    plot(xaxis, data[:,3], c='c', linestyle = "--", linewidth=1)
    plot(xaxis, data[:,0], c='c', linestyle = ":", linewidth=1)
    plot(xaxis, data[:,4], c='c', linestyle = ":", linewidth=1)        

    plt.plot([0,0],OVI_range,color = 'k')
    plt.axis([radius_range[0],radius_range[1],OVI_range[0],OVI_range[1]])
    plt.xlabel(r'$\rho/R_{vir}$')
    plt.ylabel(r'log N [OVI] [cm$^{-2}]$')
    plt.text(0.6,14.25,r'Log($M_*/M_\odot$) = %.1f'%(np.log10(hfbsmass)),size = legendsize,color = 'k')
    plt.tight_layout()
    plt.savefig(tfile+'.'+str(halo_num)+'.OVIprof.png',dpi = dpi)

    ind = np.argsort(rad_axes)
    data = np.vstack((rad_axes[ind], np.log10(hfbhim.flatten())[ind], np.log10(hfbcim.flatten())[ind], np.log10(hfbhiim.flatten())[ind], np.log10(hfbcivim.flatten())[ind], np.log10(hfbsiivim.flatten())[ind], np.log10(hfboviim.flatten())[ind]))
    data = data.T

    header = "{} {}: M_vir = {:5.4e} [Msol]; R_vir = {:.3f} [kpc]; M_* = {:5.4e} [Msol]\n{:12s}{:12s}{:12s}{:12s}{:12s}{:12s}{:12s}\n{:12s}{:12s}{:12s}{:12s}{:12s}{:12s}{:12s}".format(tfile.split("/")[6],halo_num,hfbvmass,hfbrvir,hfbsmass,"Dist", "H", "C", "HI", "CIV", "SIIV", "OVI", "[Rvir^{-1}]", "[cm^{-2}]", "[cm^{-2}]", "[cm^{-2}]", "[cm^{-2}]", "[cm^{-2}]", "[cm^{-2}]")
    np.savetxt(tfile + "." + str(halo_num) + ".columnDen.txt", data, fmt='%11.4e', header = header)
    
    return civ_rad

if __name__ == '__main__':    
    if (socket.gethostname() == "quirm.math.grinnell.edu"):
        prefix = '/home/christenc/Data/Sims/'
        outprefix = '/home/christenc/Figures/marvel/'
        dataprefix = '/home/christenc/Code/Datafiles/'           
    else:
        prefix = '/home/christensen/Storage2/UW/MolecH/Cosmo/'
        outprefix = '/home/christensen/Storage2/UW/MolecH/Cosmo/'
        dataprefix = '/home/christensen/Code/Datafiles/'          

    f = open(dataprefix+'Bordoloi2014.txt', 'r')
    bordoloidata = []
    for line in f:
        line = line.strip()
        columns = line.split()
        if len(columns) == 19 and columns[0] != 'QSOName':
            source = {}
            source['qsoname'] = columns[0]
            source['galaxy'] = columns[1]
            source['alpha'] = columns[2]
            source['delta'] = columns[3]
            source['zsys'] = float(columns[4])
            source['L'] = float(columns[5])
            source['logMstar'] = float(columns[6])
            source['R'] = float(columns[7])
            source['Rvir'] = float(columns[8])
            source['logsSFR'] = float(columns[9])
            source['logsSFr_limit'] = columns[10]
            source['logN'] = float(columns[11])
            source['logN_err'] = float(columns[12])
            source['logN_limit'] = columns[13]
            source['Wr'] = float(columns[14])
            source['Wr_err'] = float(columns[15])
            source['Wr_limit'] = columns[16]
            source['phi'] = float(columns[17])
            source['CIR3sigma'] = float(columns[18])
            bordoloidata.append(source)
    f.close()
    bordoloidata = pd.DataFrame(bordoloidata)

    f = open(dataprefix+'CIV_Burchett/Burchett2016_CIV_HI_virselect.dat', 'r')
    burchettdata = []
    for line in f:
        line = line.strip()
        columns = line.split()
        if len(columns) == 24 and columns[0] != 'field':
            source = {}
            source['field'] = columns[0]
            source['ra_gal'] = float(columns[1])
            source['dec_gal'] = float(columns[2])
            source['zgal'] = float(columns[3])
            source['NSAidx'] = columns[4]
            source['logMstar'] = float(columns[5])
            source['Rrvir'] = float(columns[6])
            source['R'] = float(columns[7])
            source['logN_CIV'] = float(columns[8])
            source['logN_CIV_err'] = float(columns[9])
            source['Wr_CIV'] = float(columns[10])
            source['Wr_CIV_err'] = float(columns[11])
            source['flag_CIV'] = float(columns[12])
            source['flag_HI'] = float(columns[13])
            source['Wr_HI'] = float(columns[14])
            source['Wr_HI_err'] = float(columns[15])
            source['logN_HI'] = float(columns[16])
            source['logN_HI_err'] = float(columns[17])            
            source['ra_qso'] = float(columns[18])
            source['dec_qso'] = float(columns[19])
            source['SFR'] = float(columns[20])
            source['SFR_err'] = float(columns[21])
            source['logMHalo'] = float(columns[22])
            source['logMHalo_single'] = float(columns[23])            
            burchettdata.append(source)
    f.close()
    burchettdata = pd.DataFrame(burchettdata)
        
    
    halo_nums = ['1', '2', '3', '5', '6', '7', '10', '11', '13']
    tfile = prefix + 'cptmarvel.cosmo25cmb/cptmarvel.cosmo25cmb.4096g5HbwK1BH/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096/ahf_200/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_name = 'cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096'
    for halo_num in halo_nums:
        print(tfile,halo_num)
        civ_rad_0 = createSpectra(tfile,halo_num,bordoloidata,burchettdata)
    
    tfile = prefix + 'rogue.cosmo25cmb/rogue.cosmo25cmb.4096g5HbwK1BH/rogue.cosmo25cmb.4096g5HbwK1BH.004096/ahf_200/rogue.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_name = 'rogue.cosmo25cmb.4096g5HbwK1BH.004096'
    #halo_nums = ['1','3', '8', '10'
    halo_nums = ['7','11','12','15','16','17','28','31']
    for halo_num in halo_nums:
        print(tfile,halo_num)        
        civ_rad_0 = createSpectra(tfile,halo_num,bordoloidata,burchettdata)
    #civ_rad_1 = createSpectra(tfile,halo_nums[0],bordoloidata,burchettdata)
    #civ_rad_2 = createSpectra(tfile,halo_nums[1],bordoloidata,burchettdata)

    tfile = prefix + 'elektra.cosmo25cmb/elektra.cosmo25cmb.4096g5HbwK1BH/elektra.cosmo25cmb.4096g5HbwK1BH.004096/ahf_200/elektra.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_name = 'elektra.cosmo25cmb.4096g5HbwK1BH.004096'
    halo_nums = ['1','2','3','4','5','10','12']
    for halo_num in halo_nums:
        print(tfile,halo_num)        
        civ_rad_0 = createSpectra(tfile,halo_num,bordoloidata,burchettdata)
    #civ_rad_3 = createSpectra(tfile,halo_nums[0],bordoloidata,burchettdata)
    #civ_rad_4 = createSpectra(tfile,halo_nums[1],bordoloidata,burchettdata)    

    tfile = prefix + 'storm.cosmo25cmb/storm.cosmo25cmb.4096g5HbwK1BH/storm.cosmo25cmb.4096g5HbwK1BH.004096/ahf_200/storm.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_name = 'storm.cosmo25cmb.4096g5HbwK1BH'
    halo_nums = ['1','2','3','4','5','6','7','8','10','11','12','14','15','23','44']
    for halo_num in halo_nums:
        print(tfile,halo_num)        
        civ_rad_0 = createSpectra(tfile,halo_num,bordoloidata,burchettdata)
    #civ_rad_5 = createSpectra(tfile,halo_nums[0],bordoloidata,burchettdata)
    #civ_rad_6 = createSpectra(tfile,halo_nums[1],bordoloidata,burchettdata)
