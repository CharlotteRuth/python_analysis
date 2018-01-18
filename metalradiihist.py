#Charlotte Christensen
#6/22/16
#Plot the cumulative radial histrograms of mass and metals that were once part of the disk. This program relies on output from metalradii.py

#Run with
#%run /home/christensen/Code/python/python_analysis/metalradiihist.py


import matplotlib as mpl
import matplotlib.pylab as plt
import matplotlib.mlab as mlab
import matplotlib.colors as colors
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import numpy as np
import math
import resource
import multiprocessing
from scipy import interpolate
from mpl_toolkits.mplot3d import Axes3D
import pynbody
import pynbody.plot as pp
import pynbody.snapshot.tipsy
import pyfits
import sys, os, glob, pynbody.bridge
import time
#import matplotlib_ref_density

if __name__ == '__main__':
    zsolar = 0.0130215
    prefix = '/home/christensen/Storage2/UW/MolecH/Cosmo/'
    finalstep = '00512'

    dir799 = prefix + 'h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g14HBWK/'
    file799 = 'h799.cosmo25cmb.3072g14HBWK'
    key799 = 'h799'
    dir516 = prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/'
    file516 = 'h516.cosmo25cmb.3072g14HBWK'
    key516 = 'h516'
    dir986 = prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/'
    file986 = 'h986.cosmo50cmb.3072g14HBWK'
    key986 = 'h986'
    dir603 = prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/'
    file603 = 'h603.cosmo50cmb.3072g14HBWK'
    key603 = 'h603'
    dir258 = prefix + 'h258.cosmo50cmb.3072g/h258.cosmo50cmb.3072g14HMbwK/'
    file258 = 'h258.cosmo50cmb.3072g14HMbwK'
    key258 = 'h258'
    dir285 = prefix + 'h285.cosmo50cmb.3072g/h285.cosmo50cmb.3072g14HMbwK/'
    file285 = 'h285.cosmo50cmb.3072g14HMbwK'
    key285 = 'h285'
    dir239 = prefix + 'h239.cosmo50cmb.3072g/h239.cosmo50cmb.3072g14HMbwK/'
    file239 = 'h239.cosmo50cmb.3072g14HMbwK'
    key239 = 'h239'
    dirs    = np.array([ dir799, dir799, dir799, dir516, dir516, dir986, dir986, dir986, dir986, dir986, dir986, dir603, dir603, dir603, dir258, dir258, dir285, dir285, dir285, dir239])
    files   = np.array([file799,file799,file799,file516,file516,file986,file986,file986,file986,file986,file986,file603,file603,file603,file258,file258,file285,file285,file285,file239])
    haloid =  np.array([ '1'   , '4'   , '6'   , '1'   , '2'   , '1'   , '2'   , '3'   , '8'   , '15'  , '16'  , '1'   , '2'    , '3'  , '1'   , '4'   , '1'   , '4'   , '9'   , '1'   ])
    masssort =np.array([10,      2,      9,      1,      8,      15,     18,     4,      0,      13,     17,     3,      7,      6,      12,     5,      11,     14,     16,     19])
    rvir = np.array([74.0411,48.9178,42.3562,86.9863,62.8630,146.7808,100.5479,86.6712,56.5890,42.4110,37.9315,180.7260,121.0822,63.5753,237.1644,57.4247,247.5753,83.4658,58.8356,250.5206])
    rvir = rvir[masssort]
    dirs = dirs[masssort]
    files = files[masssort]
    haloid = haloid[masssort]
    #['h986.cosmo50cmb.3072g14HBWK' 'h799.cosmo25cmb.3072g14HBWK' 'h986.cosmo50cmb.3072g14HBWK' 'h799.cosmo25cmb.3072g14HBWK' 'h986.cosmo50cmb.3072g14HBWK' 'h258.cosmo50cmb.3072g14HMbwK' 'h285.cosmo50cmb.3072g14HMbwK' 'h516.cosmo25cmb.3072g14HBWK' 'h799.cosmo25cmb.3072g14HBWK' 'h603.cosmo50cmb.3072g14HBWK' 'h285.cosmo50cmb.3072g14HMbwK' 'h516.cosmo25cmb.3072g14HBWK' 'h986.cosmo50cmb.3072g14HBWK' 'h986.cosmo50cmb.3072g14HBWK' 'h603.cosmo50cmb.3072g14HBWK' 'h986.cosmo50cmb.3072g14HBWK' 'h603.cosmo50cmb.3072g14HBWK' 'h258.cosmo50cmb.3072g14HMbwK' 'h285.cosmo50cmb.3072g14HMbwK' 'h239.cosmo50cmb.3072g14HMbwK']                   
    #['16' '6' '15' '4' '8' '4' '9' '2' '1' '3' '4' '1' '3' '2' '2' '1' '1' '1' '1' '1']                                                                                                                  
    #[0,    1,  2,   3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]                                                                                                                   

    #create a color map that will assign a color to a line representing each simulation                                       
    #cols = []
    #for x in np.linspace(0,1, 254):
    #    rcol = (0.472-0.567*x+4.05*x**2)/(1.+8.72*x-19.17*x**2+14.1*x**3)
    #    gcol = 0.108932-1.22635*x+27.284*x**2-98.577*x**3+163.3*x**4-131.395*x**5+40.634*x**6
    #    bcol = 1./(1.97+3.54*x-68.5*x**2+243*x**3-297*x**4+125*x**5)
    #    cols.append((rcol, gcol, bcol))
    #cols.append((1,1,1))
    #cm_rainbow = matplotlib.colors.LinearSegmentedColormap.from_list("PaulT_rainbow", cols)
                                                                           
    cmx = plt.get_cmap("viridis_r") #plt.get_cmap(cm_rainbow) 
    # values = range(0,len(dirs))
    cNorm  = colors.Normalize(vmin=9.5, vmax = 12) #12) #0, vmax=values[-1]+1)
    scalarMap = cm.ScalarMappable(norm=cNorm, cmap=cmx)
    #    print scalarMap.get_clim()   

    #    histmet, bin_edges = np.histogram(disk_parts.gas['r']/hrvir, bins = 500, range=(0,max(disk_parts.gas['r']/hrvir)), weights=disk_parts.gas['mass']*disk_parts.gas['metals'])
    #    histcummet = np.cumsum(histmet)
    #    histox, bin_edges = np.histogram(disk_parts.gas['r']/hrvir, bins = 500, range=(0,max(disk_parts.gas['r']/hrvir)), weights=disk_parts.gas['mass']*disk_parts.gas['OxMassFrac'])
    #    histcumox = np.cumsum(histox)
    #    histfe, bin_edges = np.histogram(disk_parts.gas['r']/hrvir, bins = 500, range=(0,max(disk_parts.gas['r']/hrvir)), weights=disk_parts.gas['mass']*disk_parts.gas['FeMassFrac'])
    #    histcumfe = np.cumsum(histfe)

    fig1 = plt.figure(1,figsize=(14,4))
    fig2 = plt.figure(2)
    fig3 = plt.figure(3)
    fig4 = plt.figure(4)
    fig5 = plt.figure(5,figsize=(14,4))
    gs = gridspec.GridSpec(1,3,width_ratios=[15,15,1])
    gs1 =  gridspec.GridSpec(1,2,width_ratios=[15,1])
    ax1 = fig1.add_subplot(gs[0])
    ax1b = fig1.add_subplot(gs[1])
    ax1sub = fig1.add_subplot(gs[2])
    ax2 = fig2.add_subplot(gs1[0])
    ax2sub = fig2.add_subplot(gs1[1])
    ax3 = fig3.add_subplot(111)
    ax4 = fig4.add_subplot(111)
    ax5 = fig5.add_subplot(gs[0])
    ax5b = fig5.add_subplot(gs[1])
    ax5sub = fig5.add_subplot(gs[2])

    first = ax1.plot([1,1],[0,1],linestyle = '--',color = 'k') #,cmap = cm_rainbow)
    ax1b.plot([1,1],[0,1],linestyle = '--',color = 'k')
    ax2.plot([1,1],[0,1e3*zsolar],linestyle = '--',color = 'k')
    ax2.plot([1e-1,1e3],[1,1],linestyle = '--',color = 'k')
    ax3.plot([1,1],[0,1],linestyle = '--',color = 'k')
    ax4.plot([1,1],[0,1],linestyle = '--',color = 'k')
#   first = ax5.plot([1,1],[0,1],linestyle = '--',color = 'k') #,cmap = cm_rainbow)
#   ax5b.plot([1,1],[0,1],linestyle = '--',color = 'k')
    mtot_list = []
    mtot_list = [3.2e9, 4.4e9, 4.4e9, 6.8e9, 1.1e10,1.1e10,1.2e10,1.5e10,2.4e10,2.9e10,3.4e10,3.8e10,3.8e10,5.9e10,1.0e11,1.9e11,3.4e11,7.7e11,8.8e11,9.1e11]

    for i in reversed(range(0,14)): #len(dirs))):
        #tfile = dirs[i] + files[i] + '.' + '00512' + '/' + files[i] + '.' + '00512'
        #s = pynbody.load(tfile)
        #hs = s.halos()       
        #h = hs[int(haloid[i])]
        #mtot = h['mass'].sum().in_units('Msol') #h.properties['mass']
        #mtot_list.append(mtot)
        #'h986.16' 'h799.6' 'h986.15' 'h799.4' 'h986.8' 'h258.4' 'h285.9' 'h516.2''h799.1' 'h603.3'
        #'h285.4' 'h516.1''h986.3' 'h986.2' 'h603.2' 'h986.1''h603.1' 'h258.1''h285.1' 'h239.1'

        outfilebase = dirs[i] + files[i] + '.grp' + haloid[i]
        data = np.loadtxt(outfilebase + '_rhistz.txt')
        data_log =  np.loadtxt(outfilebase + '_rhistz_log.txt')
        print(outfilebase)
        colorVal = scalarMap.to_rgba(np.log10(mtot_list[i])) #values[i])
        #linecolor = (alog10(mtot_t) - 9.5)/2.5
        masshist = np.subtract(data_log[:,1], np.append([0],data_log[0:-1,1]))
        metalhist = np.subtract(data_log[:,2], np.append([0],data_log[0:-1,2]))
        ax1.plot(data[:,0],data[:,1],color = colorVal,linestyle = '-',lw = 2)
        ax1b.plot(data[:,0],data[:,2],color = colorVal,linestyle = '-',lw = 2)
        ax2.plot(10**data_log[:,0],metalhist/masshist/zsolar,color = colorVal,linestyle = '-',lw = 2)
        if max(metalhist/masshist/zsolar) > 1:
            print("Error: ",i,outfilebase)
        ax3.plot(data[:,0],data[:,3],color = colorVal,linestyle = '-',lw = 2)
        ax4.plot(data[:,0],data[:,4],color = colorVal,linestyle = '-',lw = 2)
        ax5.plot(data[:,0]*rvir[i],data[:,1],color = colorVal,linestyle = '-',lw = 2)
        ax5b.plot(data[:,0]*rvir[i],data[:,2],color = colorVal,linestyle = '-',lw = 2)
        print(data[14,2]) #1 virial radii
#        print(data[240,2]) #17.5 virial radii

    outfilebase = '/home/christensen/Plots/outflowz_images/hgal'
    ax1.set_xlabel('Radius/R$_{vir}$')
    ax1.set_ylabel('Cumulative histogram of gas')
    ax1.axis([0.1, 30, 0, 1])
    ax1.set_xscale('log')
  
    ax1b.set_xlabel('Radius/R$_{vir}$')
    ax1b.set_ylabel('Cumulative histogram of metals')
    ax1b.axis([0.1, 30, 0, 1])
    ax1b.set_xscale('log')

    cb = mpl.colorbar.ColorbarBase(ax1sub, cmap=cmx, norm=cNorm)
    cb.set_label('Log Virial Mass [M$_\odot$]')
    fig1.show()
    fig1.savefig(outfilebase + '_rhistmass_z.png')

    ax2.set_xlabel('Radius/R$_{vir}$')
    ax2.set_ylabel('Metallicity [Z/Z$_{\odot}$]')
    ax2.axis([0.1, 30, 1e-3, 1e2])
#    ax2.axis([0.1, 30, 1e-1/zsolar, 3e1/zsolar])
    ax2.set_xscale('log')
    ax2.set_yscale('log')

    cb = mpl.colorbar.ColorbarBase(ax2sub, cmap=cmx, norm=cNorm)
    cb.set_label('Log Virial Mass [M$_\odot$]')
    fig2.show()
    fig2.savefig(outfilebase + '_rhistmass_ratio.png')

    ax3.set_xlabel('Radius/R$_{vir}$')
    ax3.set_ylabel('Cumulative histogram of oxygen')
    ax3.axis([0.1, 20, 0, 1])
    ax3.set_xscale('log')
    fig3.show()
    fig3.savefig(outfilebase + '_rhistox.png')

    ax4.set_xlabel('Radius/R$_{vir}$')
    ax4.set_ylabel('Cumulative histogram of iron')
    ax4.axis([0.1, 20, 0, 1])
    ax4.set_xscale('log')
    fig4.show()
    fig4.savefig(outfilebase + '_rhistfe.png')

    ax5.set_xlabel('Radius [kpc]')
    ax5.set_ylabel('Cumulative histogram of gas')
    ax5.axis([1, 3000, 0, 1])
    ax5.set_xscale('log')
  
    ax5b.set_xlabel('Radius [kpc]')
    ax5b.set_ylabel('Cumulative histogram of metals')
    ax5b.axis([1, 3000, 0, 1])
    ax5b.set_xscale('log')

    cb = mpl.colorbar.ColorbarBase(ax5sub, cmap=cmx, norm=cNorm)
    cb.set_label('Log Virial Mass [M$_\odot$]')
    fig5.show()
    fig5.savefig(outfilebase + '_rhistmass_z.png')

