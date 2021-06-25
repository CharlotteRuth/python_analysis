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

#Run with
#%run /home/christensen/Code/python/python_analysis/plotColumnDen.py
#or, on quirm,
#%run /home/christenc/Code/python/python_analysis/plotColumnDen.py
#ipython --pylab

def make_rs(im,width): #width should be width of box in kpc
    xsize, ysize = np.shape(im)
    x = np.arange(-xsize/2, xsize/2)/float(xsize)*width
    y = np.arange(-ysize/2, ysize/2)/float(xsize)*width
    xs, ys = np.meshgrid(x,y)
    return np.sqrt(xs**2 + ys**2)

def plotColumnDen(tfile,halo_num,bordoloidata,burchettdata):
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
    
#tfile = sys.argv[1]
    Carbon   = 0.00213
    Oxygen   = 0.00541
    Silicon  = 0.00067
    Iron     = 0.00117

    hfb_all = pynbody.load(tfile)
    h = hfb_all.halos()

    halo_num = int(halo_num)
       
    hfbsmass = np.sum(h[halo_num].stars['mass'].in_units('Msol'))
    hfb_all.physical_units()
    pynbody.analysis.angmom.faceon(h[halo_num])
            
    hfbrvir = np.max(h[halo_num].gas['r'])
    hfb_all['velmag'] = np.sqrt(hfb_all['vx']**2 + hfb_all['vy']**2 + hfb_all['vz']**2)
    #velf = filt.LowPass('velmag',h[halo_num].properties['v_esc'])
    velf = filt.BandPass('vz',-1*h[halo_num].properties['v_esc'],h[halo_num].properties['v_esc'])
    hfb = hfb_all[velf]
    notdiskf = filt.Not(filt.Disc('8 kpc','2 kpc'))
        
    hfb.gas['hiif'] = pynbody.analysis.ionfrac.calculate(hfb.gas,ion='hi')
    hfb.gas['hiden'] = hfb.gas['rho']*hfb.gas['hydrogen']*hfb.gas['hiif']
    hfb.gas['hiden'].units = hfb.gas['rho'].units
    plt.close()
    hfbhiim = pynbody.plot.image(hfb.gas[notdiskf],qty='hiden', clear=False,
        units='m_p cm^-2', width=np.round(hfbrvir/10)*10*2.5, show_cbar=False, 
        vmin=1e12,vmax=1e15, log=True, threaded=False)    
    hfb.gas['civif'] = pynbody.analysis.ionfrac.calculate(hfb.gas,ion='civ')
    hfb.gas['civden'] = hfb.gas['rho']*hfb.gas['OxMassFrac']*Carbon/Oxygen*hfb.gas['civif']
    hfb.gas['civden'].units = hfb.gas['rho'].units
    plt.close()
    hfbcivim = pynbody.plot.image(hfb.gas[notdiskf],qty='civden', clear=False,
        units='12 m_p cm^-2', width=np.round(hfbrvir/10)*10*2.5, show_cbar=False, 
        vmin=1e12,vmax=1e15, log=True)

    rs = make_rs(hfbcivim,np.round(hfbrvir/10)*10*2.5)/hfbrvir
    #rad_axes = np.log10(rs.flatten() + 1e-6)
    rad_axes = (rs.flatten() + 1e-6)

    radius_range = [0,2.5/2]    
    HI_range = [12.5,17]
    CIV_range = [11.8,15]
    plt.close()
    plt.figure(1,figsize=(plt_width,plt_width*aspect_ratio))
    hi_rad = pynbody.plot.hist2d(rad_axes,np.log10(hfbhiim.flatten()), cmap=plt.cm.gist_yarg)
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
    civ_rad = pynbody.plot.hist2d(rad_axes,np.log10(hfbcivim.flatten()), cmap=plt.cm.gist_yarg)
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
        
    
    halo_nums = ['1']
    tfile = prefix + 'cptmarvel.cosmo25cmb/cptmarvel.cosmo25cmb.4096g5HbwK1BH/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_name = 'cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096'
    civ_rad_0 = plotColumnDen(tfile,halo_nums[0],bordoloidata,burchettdata)
    
    tfile = prefix + 'rogue.cosmo25cmb/rogue.cosmo25cmb.4096g5HbwK1BH/rogue.cosmo25cmb.4096g5HbwK1BH.004096/rogue.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_name = 'rogue.cosmo25cmb.4096g5HbwK1BH.004096'
    halo_nums = ['1','3']
    civ_rad_1 = plotColumnDen(tfile,halo_nums[0],bordoloidata,burchettdata)
    civ_rad_2 = plotColumnDen(tfile,halo_nums[1],bordoloidata,burchettdata)

    tfile = prefix + 'elektra.cosmo25cmb/elektra.cosmo25cmb.4096g5HbwK1BH/elektra.cosmo25cmb.4096g5HbwK1BH.004096/elektra.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_name = 'elektra.cosmo25cmb.4096g5HbwK1BH.004096'
    halo_nums = ['1','2']
    civ_rad_3 = plotColumnDen(tfile,halo_nums[0],bordoloidata,burchettdata)
    civ_rad_4 = plotColumnDen(tfile,halo_nums[1],bordoloidata,burchettdata)    

    tfile = prefix + 'storm.cosmo25cmb/storm.cosmo25cmb.4096g5HbwK1BH/storm.cosmo25cmb.4096g5HbwK1BH.004096/storm.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_name = 'storm.cosmo25cmb.4096g5HbwK1BH'
    halo_nums = ['1','2']
    civ_rad_5 = plotColumnDen(tfile,halo_nums[0],bordoloidata,burchettdata)
    civ_rad_6 = plotColumnDen(tfile,halo_nums[1],bordoloidata,burchettdata)
