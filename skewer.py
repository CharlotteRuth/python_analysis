#Charlotte Christensen
#5/23/18
#Select particles along a line of sight and make a velocity profile


#Run with
#%run /home/christensen/Code/python/python_analysis/skewer.py
#or, on quirm,
#%run /home/christenc/Code/python/python_analysis/skewer.py
#ipython --pylab

import numpy as np
import pynbody
import pynbody.filt as filt
import pynbody.units as units
from pynbody.filt import *
from pynbody import tipsy
import pynbody.plot.sph as sph
from array import array
import matplotlib.pylab as plt
import matplotlib.mlab as mlab
from pylab import *
import matplotlib.pyplot  as pyplot
from scipy.optimize import curve_fit
from numpy import sqrt, pi, exp, linspace, random
from scipy.integrate import trapz as trapz
import sys, os, glob, pickle
import socket
import pandas as pd
import math

def plt_skewer(path,tfile,haloid):
    sphere = pynbody.load(path + tfile + "." + str(haloid) + ".tipsy.smooth.std")
    sphere.physical_units()

    files = glob.glob(path + "CGM." + haloid + "/specaim*kpc")
    x0_arr = np.zeros(len(files))
    z0_arr = np.zeros(len(files))
    ind = 0
    for f in files:
        f = (f.split("/"))[-1]
        filebase = (f.split("specaim."))[1]
        binzfile = open(path+"CGM." + haloid + "/binzfile." + filebase,"r")
        temp = binzfile.readline()
        line = binzfile.readline()
        binzfile.close()
        x0_arr[ind] = int(float(line.split()[12])*sphere.properties['boxsize'].in_units('kpc'))
        z0_arr[ind] = int(float(line.split()[13])*sphere.properties['boxsize'].in_units('kpc'))
        ind += 1

    hiif = pynbody.analysis.ionfrac.calculate(sphere.gas,ion='hi')
    sphere.gas['hiden'] = sphere.gas['rho']*sphere.gas['hydrogen']*hiif
    sphere.gas['hiden'].units = sphere.gas['rho'].units
    spherehiim = pynbody.plot.image(sphere.gas,qty='hiden',
                    units='m_p cm^-2',  clear=False, width = 200)
#                            title='log$_{10}$ N$_{HI}$ cm$^{-2}$',
#                   vmin=12,vmax=20);
    savefig(path + tfile + '.'+haloid + '.HImap.png')
    
    spherevel = pynbody.plot.image(sphere.gas,qty='vz',cmap = "RdYlBu",
                    units='km s^-1',  clear=True, width = 200,log=False)
    savefig(path + tfile + '.'+haloid + '.velmap_faceon.png')
    pynbody.analysis.angmom.sideon(sphere, cen=(0,0,0))
    spherevel = pynbody.plot.sph.velocity_image(sphere.gas,qty='vz',vector_color = 'black',cmap = "RdYlBu",units='km s^-1',  clear=True, width = 200,log=False)
    savefig(path + tfile + '.'+haloid + '.velmap_edgeon.png')
    plt.clear()
    
    width = 400
    resolution = width
    vmin = -150
    vmax = 150

    #x0 = -10 #-52.0
    #x0_arr = [0,-10,-20,-30,-40,-50,-60,-70,-80,-90,-100,-110,-120]

    ind = 0
    for x0 in x0_arr:
        z0 = z0_arr[ind] #30.0
        sphere.gas['z'] = sphere.gas['z'] - z0
        distance = np.sqrt((sphere.gas['x'].in_units('kpc') - x0)**2 + sphere.gas['z'].in_units('kpc')**2)
        sphere.gas['reldist'] = sphere.gas['smooth'].in_units('kpc') - distance
        
        skewer = sphere.gas[pynbody.filt.HighPass('reldist', '0 kpc')]
    
#pynbody.plot.image(skewer,clear=True, width = 1000,log=False)
        spherevel = pynbody.plot.image(skewer,qty='vy',cmap = "RdYlBu", units='km s^-1', clear=True, width = width,resolution = resolution,log=False,vmin = vmin,vmax = vmax,noplot = True)
#savefig(path + tfile + '.'+haloid + 'x%.1f_kpc_y%.1f_kpc.velmap.png' % (x0,z0))

        fhi = np.copy(pynbody.analysis.ionfrac.calculate(skewer,ion='hi'))      
        normalize = np.copy(skewer['mass'])*fhi*len(skewer['mass'])/np.sum(np.copy(skewer['mass'])*fhi)
        skewer['mass'] = skewer['mass']*normalize
        spherevelHI = pynbody.plot.image(skewer,qty='vy',cmap = "RdYlBu",
                            units='km s^-1',  clear=True, width = width,resolution = resolution,log=False,vmin = vmin, vmax = vmax,noplot=True)
        skewer['mass'] = skewer['mass']/normalize

        civ = np.copy(pynbody.analysis.ionfrac.calculate(skewer,ion='civ'))
        normalize = np.copy(skewer['mass'])*civ*len(skewer['mass'])/np.sum(np.copy(skewer['mass'])*civ)
        skewer['mass'] = skewer['mass']*normalize
        spherevelciv = pynbody.plot.image(skewer,qty='vy',cmap = "RdYlBu",
                            units='km s^-1',  clear=True, width = width,resolution = resolution,log=False,vmin = vmin, vmax = vmax,noplot = True)
        skewer['mass'] = skewer['mass']/normalize

        ovi = np.copy(pynbody.analysis.ionfrac.calculate(skewer,ion='ovi'))
        normalize = np.copy(skewer['mass'])*ovi*len(skewer['mass'])/np.sum(np.copy(skewer['mass'])*ovi)
        skewer['mass'] = skewer['mass']*normalize
        spherevelovi = pynbody.plot.image(skewer,qty='vy',cmap = "RdYlBu",
                            units='km s^-1',  clear=True, width = width,resolution = resolution,log=False,vmin = vmin, vmax = vmax,noplot = True)
        skewer['mass'] = skewer['mass']/normalize

        fig_hist = plt.figure(1)
        vy,v_hist = np.histogram(skewer['vy'].in_units('km s^-1'),weights = skewer['mass'],range = (vmin,vmax),bins = 100)
        ax1 = fig_hist.add_subplot(411)
        ax1.set_title('x = %.1f; y = %.1f' % (x0,z0))
        ax1.hist(vy,v_hist)
        ax1.set_yscale("log")
        ax1.text(-120,20,"Total gas")
        vy,v_hi_hist = np.histogram(skewer['vy'].in_units('km s^-1'),weights = skewer['mass']*fhi,range = (vmin,vmax),bins = 100)        
        ax2 = fig_hist.add_subplot(412)
        ax2.hist(vy,v_hi_hist)
        ax2.set_yscale("log")
        ax2.text(-120,20,"HI")
        vy,v_civ_hist = np.histogram(skewer['vy'].in_units('km s^-1'),weights = skewer['mass']*civ,range = (vmin,vmax),bins = 100)        
        ax3 = fig_hist.add_subplot(413)
        ax3.hist(vy,v_civ_hist)
        ax3.set_yscale("log")
        ax3.text(-120,20,"CIV")
        vy,v_ovi_hist = np.histogram(skewer['vy'].in_units('km s^-1'),weights = skewer['mass']*ovi,range = (vmin,vmax),bins = 100)        
        ax4 = fig_hist.add_subplot(414)
        ax4.hist(vy,v_ovi_hist)
        ax4.set_yscale("log")
        ax4.text(-120,20,"OVI")
        ax4.set_xlabel("Velocity [km/s]")
        plt.show()
        savefig(path + tfile + '.'+haloid + 'x%.1f_kpc_y%.1f_kpc.hist4.png' % (x0,z0))
        plt.close(fig_hist)
        
        dx = 20 #shape(spherevel)[0]/10*2.0
        xmin = int((shape(spherevel)[0]/2 - dx/2) + x0)
        xmax = int((shape(spherevel)[0]/2 + dx/2) + x0)
        fig = plt.figure(figsize = (6,4))
        sps = [fig.add_subplot(1,4,1),fig.add_subplot(1,4,2),fig.add_subplot(1,4,3),fig.add_subplot(1,4,4)]
        ims = sps[0].imshow(spherevel[:,xmin:xmax],extent = (dx/2*-1,dx/2,width/2*-1,width/2),cmap = "RdYlBu",vmin = vmin, vmax = vmax)
        sps[1].imshow(spherevelHI[:,xmin:xmax],extent = (dx/2*-1,dx/2,width/2*-1,width/2),cmap = "RdYlBu",vmin = vmin, vmax = vmax)
        sps[2].imshow(spherevelciv[:,xmin:xmax],extent = (dx/2*-1,dx/2,width/2*-1,width/2),cmap = "RdYlBu",vmin = vmin, vmax = vmax)
        sps[3].imshow(spherevelovi[:,xmin:xmax],extent = (dx/2*-1,dx/2,width/2*-1,width/2),cmap = "RdYlBu",vmin = vmin, vmax = vmax)
        fig.colorbar(ims).set_label("v_y [km s^-1]")
        plt.show()
        plt.close(fig)
        savefig(path + tfile + '.'+haloid + 'x%.1f_kpc_y%.1f_kpc.velmap4.png' % (x0,z0))
        sphere.gas['z'] = sphere.gas['z'] + z0
        ind += 1
        
    fig.clear()


if __name__ == '__main__':    
    if (socket.gethostname() == "quirm"):
        prefix = '/home/christenc/Data/Sims/'
    else:
        prefix = '/home/christensen/Storage2/UW/MolecH/Cosmo/'
    
    path = prefix + "rogue.cosmo25cmb/rogue.cosmo25cmb.4096g5HbwK1BH/rogue.cosmo25cmb.4096g5HbwK1BH.004096/"
    tfile = "rogue.cosmo25cmb.4096g5HbwK1BH.004096"
    haloid = "1"
    
    plt_skewer(path,tfile,haloid)
