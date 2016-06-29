#Charlotte Christensen
#6/22/16
#Plot the cumulative radial histrograms of mass and metals that were once part of the disk. This program relies on output from metalradii.py

import matplotlib as mpl
import matplotlib.pylab as plt
import matplotlib.mlab as mlab
import matplotlib.colors as colors
import matplotlib.cm as cm
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

if __name__ == '__main__':
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
    dirs = dirs[masssort]
    files = files[masssort]
    haloid = haloid[masssort]
    #['h986.cosmo50cmb.3072g14HBWK' 'h799.cosmo25cmb.3072g14HBWK' 'h986.cosmo50cmb.3072g14HBWK' 'h799.cosmo25cmb.3072g14HBWK' 'h986.cosmo50cmb.3072g14HBWK' 'h258.cosmo50cmb.3072g14HMbwK' 'h285.cosmo50cmb.3072g14HMbwK' 'h516.cosmo25cmb.3072g14HBWK' 'h799.cosmo25cmb.3072g14HBWK' 'h603.cosmo50cmb.3072g14HBWK' 'h285.cosmo50cmb.3072g14HMbwK' 'h516.cosmo25cmb.3072g14HBWK' 'h986.cosmo50cmb.3072g14HBWK' 'h986.cosmo50cmb.3072g14HBWK' 'h603.cosmo50cmb.3072g14HBWK' 'h986.cosmo50cmb.3072g14HBWK' 'h603.cosmo50cmb.3072g14HBWK' 'h258.cosmo50cmb.3072g14HMbwK' 'h285.cosmo50cmb.3072g14HMbwK' 'h239.cosmo50cmb.3072g14HMbwK']                   
    #['16' '6' '15' '4' '8' '4' '9' '2' '1' '3' '4' '1' '3' '2' '2' '1' '1' '1' '1' '1']                                                                                                                  
    #[0,    1,  2,   3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]                                                                                                                   

    #create a color map that will assign a color to a line representing each simulation                                                                                                                   
    jet = cmx = plt.get_cmap('jet')
    values = range(len(dirs))
    cNorm  = colors.Normalize(vmin=0, vmax=values[-1])
    scalarMap = cm.ScalarMappable(norm=cNorm, cmap=jet)
    #print scalarMap.get_clim()   

#    histmet, bin_edges = np.histogram(disk_parts.gas['r']/hrvir, bins = 500, range=(0,max(disk_parts.gas['r']/hrvir)), weights=disk_parts.gas['mass']*disk_parts.gas['metals'])
#    histcummet = np.cumsum(histmet)
#    histox, bin_edges = np.histogram(disk_parts.gas['r']/hrvir, bins = 500, range=(0,max(disk_parts.gas['r']/hrvir)), weights=disk_parts.gas['mass']*disk_parts.gas['OxMassFrac'])
#    histcumox = np.cumsum(histox)
#    histfe, bin_edges = np.histogram(disk_parts.gas['r']/hrvir, bins = 500, range=(0,max(disk_parts.gas['r']/hrvir)), weights=disk_parts.gas['mass']*disk_parts.gas['FeMassFrac'])
#    histcumfe = np.cumsum(histfe)

    fig1 = plt.figure(1)
    fig2 = plt.figure(2)
    fig3 = plt.figure(3)
    fig4 = plt.figure(4)
    ax1 = fig1.add_subplot(111)
    ax2 = fig2.add_subplot(111)
    ax3 = fig3.add_subplot(111)
    ax4 = fig4.add_subplot(111)

    for i in range(0,len(dirs)):
        outfilebase = dirs[i] + files[i] + '.grp' + haloid[i]
        data = np.loadtxt(outfilebase + '_rhistz.txt')
        colorVal = scalarMap.to_rgba(values[i])

        ax1.plot(data[:,0],data[:,1],color = colorVal)
        ax2.plot(data[:,0],data[:,2],color = colorVal)
        ax3.plot(data[:,0],data[:,3],color = colorVal)
        ax4.plot(data[:,0],data[:,4],color = colorVal)

    outfilebase = '/home/christensen/Plots/outflowz_images/hgal'
    ax1.set_xlabel('Radius/R_vir')
    ax1.set_ylabel('Cumulative histogram of gas')
    ax1.axis([0, 20, 0, 1])
    fig1.savefig(outfilebase + '_rhistmass.png')
  
    ax2.set_xlabel('Radius/R_vir')
    ax2.set_ylabel('Cumulative histogram of metals')
    ax2.axis([0, 20, 0, 1])
    fig2.savefig(outfilebase + '_rhistz.png')

    ax3.set_xlabel('Radius/R_vir')
    ax3.set_ylabel('Cumulative histogram of oxygen')
    ax3.axis([0, 20, 0, 1])
    fig3.savefig(outfilebase + '_rhistox.png')

    ax4.set_xlabel('Radius/R_vir')
    ax4.set_ylabel('Cumulative histogram of iron')
    ax4.axis([0, 20, 0, 1])
    fig4.savefig(outfilebase + '_rhistfe.png')


