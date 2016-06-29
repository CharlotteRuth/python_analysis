#Charlotte Christensen                                                                                                                                                                                            
# 6/22/16   
#Calculate radial distribution of mass and metals that were once part of the disk

import matplotlib as mpl
mpl.use('Agg') #This command ensures that the plots are not displayed, allowing it to work in parallel. Because of it, you can not run python with the --pylab option
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

def loadfiles(task):
#Load an individual simulation
    tbegin = time.clock()

    dirs,files,grp,finalstep = task
    tfile = dirs + files + '.' + finalstep + '/' + files + '.' + finalstep
    print(tfile + '.' + grp)
    s = pynbody.load(tfile)
    s.physical_units()
    print("simulation loaded")
    hs = s.halos()       
    h = hs[int(grp)]
    print("halos loaded")
    pynbody.analysis.halo.center(h,mode='com') #not as accurate as hyb but saves memory and time
    
    pynbody.analysis.angmom.sideon(h)
    print('halo aligned')
    hrvir = np.max(h.gas['r'])      
    #Read in iords of all particle that have been in the disk
    disk_iord_file = pyfits.open(dirs + 'grp' + grp + '.reaccrdisk_iord.fits')
    disk_iord = disk_iord_file[0].data
    indicies = np.in1d(s['iord'],disk_iord)
    disk_parts = s[np.nonzero(indicies)]
    #Make images
    outfilebase = dirs + files + '.grp' + grp
    print('Disk particles selected, now write to ' + outfilebase)
    outflowImages(s,disk_parts,hrvir,outfilebase)
    pltZvsR(disk_parts,hrvir,outfilebase)
    print 'Memory usage: %s (kb)' % resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    tend = time.clock()
    print 'Time elapsed: %f ' % (float(tend) - float(tbegin))


def outflowImages(s,disk_parts,hrvir,outfilebase):
#Plots all material that has been in the disk of the galaxy
    width = 20*hrvir
    fig = plt.figure(0)
    ax = fig.add_subplot(1,1,1)
    sim = pp.sph.image(s.gas,
                    units='m_p cm^-2', width=width, resolution=500, clear=False,
                    title='Outflow',show_cbar=True,vmin=1e16,vmax=1e22,
                    subplot = ax, cmap='ocean_r');
    sim_outflow = pp.sph.image(disk_parts.gas,
                            units='m_p cm^-2', width=width, resolution = 500, clear=False,
                            title='Outflow',show_cbar=True,vmin=1e16,vmax=1e22,
                            subplot = ax, cmap='ocean_r',noplot  = 1);
    circ=plt.Circle((0,0), radius = hrvir, color='b', fill=False)#,ec='black',fc='none')
    dx = width/500
    xarr = np.linspace(-1*width/2 + dx,1*width/2 - dx,num=500)
    yarr = xarr
    plt.contour(xarr,yarr,np.log10(sim_outflow),levels=[16,17,18,19,20,21,22],colors='k',)
    ax.add_patch(circ)
    fig.savefig(outfilebase + '_outflow2.png')
    plt.show()
    plt.close()
    
#Plots all material that has been in the disk of the galaxy
    fig = plt.figure(1)
    ax = fig.add_subplot(1,1,1)
    sim = pp.sph.image(s.gas,
                    units='m_p cm^-2', width=width, resolution=500, clear=False,
                    title='Outflow',show_cbar=True,vmin=1e16,vmax=1e22,
                    subplot = ax, cmap='ocean_r',noplot  = 1);
    sim_outflow = pp.sph.image(disk_parts.gas,
                            units='m_p cm^-2', width=width, resolution = 500, clear=False,
                            title='Outflow',show_cbar=True,vmin=1e16,vmax=1e22,
                            subplot = ax, cmap='ocean_r');
    circ=plt.Circle((0,0), radius = hrvir, color='b', fill=False)#,ec='black',fc='none')
    plt.contour(xarr,yarr,np.log10(sim),levels=[16,17,18,19,20,21,22],colors='k',)
    ax.add_patch(circ)
    fig.savefig(outfilebase + '_outflow.png')
    plt.show()
    plt.close()

#Plots all material that has been in the disk of the galaxy
    fig = plt.figure(2)
    ax = fig.add_subplot(1,1,1)
    disk_parts.gas['oxden'] = disk_parts.gas['rho']*disk_parts.gas['OxMassFrac']
    disk_parts.gas['oxden'].units = disk_parts.gas['rho'].units
    sim = pp.sph.image(disk_parts.gas,qty='oxden',
                    units='16 m_p cm^-2', width=width,  clear=False,
                    title='Oxygen in Outflow',show_cbar=True,vmin=1e12,vmax=1e18,
                    subplot = ax, cmap='ocean_r');
    circle=plt.Circle((0,0), radius = hrvir,ec='black',fc='none')
    ax.add_artist(circle)
    fig.savefig(outfilebase + '_ox.png')
    plt.show()
    plt.close()

#Plot oxygen fraction
#Plots all material that has been in the disk of the galaxy
    fig = plt.figure(3)
    ax = fig.add_subplot(1,1,1)
    sim = pp.sph.image(disk_parts.gas,qty='OxMassFrac',
                    width=width/2,  clear=False,
                    title='Oxygen Fraction in Outflow',show_cbar=True,vmin=1e-4,vmax=1e-1,
                    subplot = ax, cmap='ocean_r');
    circle=plt.Circle((0,0), radius = hrvir,ec='black',fc='none')
    ax.add_artist(circle)
    fig.savefig(outfilebase + '_oxfrac.png')
    plt.show()
    plt.close()

def pltZvsR(disk_parts,hrvir,outfilebase,color = 'k'):
#Histogram of gas vs radius
#disk_parts.gas['radius'] = sqrt(disk_parts.gas['x'].in_units('kpc')**2 + disk_parts.gas['y'].in_units('kpc')**2 + disk_parts.gas['z'].in_units('kpc')**2)
#
    hist, bin_edges = np.histogram(disk_parts.gas['r']/hrvir, bins = 500, range=(0,max(disk_parts.gas['r']/hrvir)), weights=disk_parts.gas['mass'])
    histcum = np.cumsum(hist)
#plt.bar(bin_edges[:-1], hist, width = min(diff(bin_edges)))
#plt.semilogy((bin_edges[:-1] + bin_edges[1:])/2,histcum/max(histcum))
    fig1 = plt.figure(4)
    ax1 = fig1.add_subplot(111)
    ax1.plot((bin_edges[:-1] + bin_edges[1:])/2,histcum/max(histcum),color = color)
    ax1.set_xlabel('Radius/R_vir')
    ax1.set_ylabel('Cumulative histogram of gas')
    ax1.axis([0, 20, 0, 1])
    fig1.savefig(outfilebase + '_rhistmass.png')
  
    histmet, bin_edges = np.histogram(disk_parts.gas['r']/hrvir, bins = 500, range=(0,max(disk_parts.gas['r']/hrvir)), weights=disk_parts.gas['mass']*disk_parts.gas['metals'])
    histcummet = np.cumsum(histmet)
    fig2 = plt.figure(5)
    ax2 = fig2.add_subplot(111)
    ax2.plot((bin_edges[:-1] + bin_edges[1:])/2,histcummet/max(histcummet),color = color)
    ax2.set_xlabel('Radius/R_vir')
    ax2.set_ylabel('Cumulative histogram of metals')
    ax2.axis([0, 20, 0, 1])
    fig2.savefig(outfilebase + '_rhistz.png')

    histox, bin_edges = np.histogram(disk_parts.gas['r']/hrvir, bins = 500, range=(0,max(disk_parts.gas['r']/hrvir)), weights=disk_parts.gas['mass']*disk_parts.gas['OxMassFrac'])
    histcumox = np.cumsum(histox)
#    plt.plot((bin_edges[:-1] + bin_edges[1:])/2,histcumox/max(histcumox))

    histfe, bin_edges = np.histogram(disk_parts.gas['r']/hrvir, bins = 500, range=(0,max(disk_parts.gas['r']/hrvir)), weights=disk_parts.gas['mass']*disk_parts.gas['FeMassFrac'])
    histcumfe = np.cumsum(histfe)
#    plt.plot((bin_edges[:-1] + bin_edges[1:])/2,histcumfe/max(histcumfe))

    np.savetxt(outfilebase + '_rhistz.txt', np.c_[(bin_edges[:-1] + bin_edges[1:])/2,histcum/max(histcum),histcummet/max(histcummet),histcumox/max(histcumox),histcumfe/max(histcumfe)])


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
    
    tasks = []
    #for i in range(0,len(dirs)): 
    for i in range(0,len(dirs)):   
        tasks.append((dirs[i],files[i],haloid[i],finalstep))

    pool = multiprocessing.Pool(processes=1,maxtasksperchild=1)
    results = pool.map_async(loadfiles,tasks)
    pool.close()
    pool.join()

"""
    for i in range(0,len(dirs)-1):
        
        loadfiles(tasks[i])
        #colorVal = scalarMap.to_rgba(values[i])
        #pltZvsR(disk_parts,hrvir,'/home/christensen/Plots/outflowz_images/hgal',color = colorVal)
        s = 0
        hs = 0
        h = 0
        disk_parts = 0
        print 'Memory usage: %s (kb)' % resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
""" 
       


