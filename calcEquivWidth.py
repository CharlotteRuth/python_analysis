#Charlotte Christensen
#2/27/18
#Iterate through the binzfiles produced by specexbin and determine the equivalent widths for the ions


#Run with
#%run /home/christensen/Code/python/python_analysis/calcEquivWidth.py
#or, on quirm,
#%run /home/christenc/Code/python/python_analysis/calcEquivWidth.py
#ipython --pylab

import numpy as np
import pynbody
import pynbody.filt as filt
import pynbody.units as units
from pynbody.filt import *
from pynbody import tipsy
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

def tick_function(x_arr, lambda0):
    vel = (x_arr - lambda0)/lambda0*3.0e5
    return ["%.1f" % v for v in vel]


#For each ion, calculate the equivalent width
def calc_eq_width(lambda0, tau, redshift, lambda_sp, flux_sp, title):
    intensity = exp(-1 * tau)
    lambda_arr = redshift*lambda0 + lambda0
    area = trapz(1 - intensity, dx=(abs(lambda_arr[1] - lambda_arr[0])))
    eq_width = area/ (max(intensity))
    return eq_width

    """
    new_tick_locations = np.array([-1000.0,-800.0,-600.0,-400.0,-200.0,0,200.0,400.0,600.0,800.0,1000.0])/3e5*lambda0 + lambda0
    fig1 = plt.figure(1)
    ax1 = fig1.add_subplot(111)
    ax2 = ax1.twiny()
    ax1.plot(lambda_arr,intensity)
    ax2.set_xlim(ax1.get_xlim())
    ax2.set_xticks(new_tick_locations)
    ax2.set_xticklabels(tick_function(new_tick_locations,lambda0))
    ax1.plot(lambda_sp,flux_sp)
    ax1.plot([lambda0, lambda0],ax1.get_ylim(),color = 'k',linestyle = 'dashed')
    ax1.text(min(lambda_sp),0.6,title)
    ax1.set_xlabel("Lambda [Angstrom]")
    ax2.set_xlabel("Velocity [km/s]")
    ax1.set_ylabel("Normalized Flux")
    ax1.text(min(lambda_sp),0.5,"%.2f Angstrom" % eq_width)
    plt.show()
    wait = raw_input("PRESS ENTER TO CONTINUE.")
    fig1.clear()
    """

def analyze_los(path,filebase,tfile):
    #For each line of sight file, list impact parameter, angle, and equivalent widths of ions

    #Read in simulation file to get the length units
    s = pynbody.load(tfile)

    #Read in the binzfile to determine the location of the line of sight
    binzfile = open(path+"/binzfile." + filebase,"r")
    temp = binzfile.readline()
    line = binzfile.readline()
    binzfile.close()
    x_pos = float(line.split()[12])*s.properties['boxsize'].in_units('kpc a')
    y_pos = float(line.split()[13])*s.properties['boxsize'].in_units('kpc a')
    b_param = np.sqrt(x_pos**2 + y_pos**2)
    angle = math.atan2(y_pos,x_pos)

    #Read in the specaim file to determine the optical depth
    r0 = np.loadtxt(path+"/specaim." + filebase)
    redshift = r0[:,0]
    ion9 = 1
    if ion9:
        HI_tau = r0[:,7] #optical depth
        HeI_tau = r0[:,11]
        CIII_tau = r0[:,15]
        CIV_tau = r0[:,19]
        OIV_tau = r0[:,23]
        OVI_tau = r0[:,27]
        NeIII_tau = r0[:,31]
        MgII_tau = r0[:,35]
        SiIV_tau = r0[:,39]

    s0_HI = np.loadtxt(path+"/specaim." + filebase + "_H1216.raw")
    lambda_sp_HI = s0_HI[:,0]
    vel_sp_HI = s0_HI[:,1]
    flux_sp_HI = s0_HI[:,2]
    noice_sp_HI = s0_HI[:,3]
    rho_sp_HI = s0_HI[:,4]
    temp_sp_HI = s0_HI[:,5]
    met_sp_HI = s0_HI[:,6]

    s0_civ = np.loadtxt(path+"/specaim." + filebase + "_CIV1548.raw")
    lambda_sp_civ = s0_civ[:,0]
    vel_sp_civ = s0_civ[:,1]
    flux_sp_civ = s0_civ[:,2]
    noice_sp_civ = s0_civ[:,3]
    rho_sp_civ = s0_civ[:,4]
    temp_sp_civ = s0_civ[:,5]
    met_sp_civ = s0_civ[:,6]

    s0_ovi = np.loadtxt(path+"/specaim." + filebase + "_OVI1032.raw")
    lambda_sp_ovi = s0_ovi[:,0]
    vel_sp_ovi = s0_ovi[:,1]
    flux_sp_ovi = s0_ovi[:,2]
    noice_sp_ovi = s0_ovi[:,3]
    rho_sp_ovi = s0_ovi[:,4]
    temp_sp_ovi = s0_ovi[:,5]
    met_sp_ovi = s0_ovi[:,6]

    s0_siiv = np.loadtxt(path+"/specaim." + filebase + "_SiIV1394.raw")
    lambda_sp_siiv = s0_siiv[:,0]
    vel_sp_siiv = s0_siiv[:,1]
    flux_sp_siiv = s0_siiv[:,2]
    noice_sp_siiv = s0_siiv[:,3]
    rho_sp_siiv = s0_siiv[:,4]
    temp_sp_siiv = s0_siiv[:,5]
    met_sp_siiv = s0_siiv[:,6]    

    fig1 = plt.figure(1)
    ax = fig1.add_subplot(111)
    ax.set_ylabel("Normalized Flux")
    ax.set_xlabel("Lambda [Angstrom]")
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
    new_tick_locations = np.array([-1000.0,-800.0,-600.0,-400.0,-200.0,0,200.0,400.0,600.0,800.0,1000.0])/3e5
    
    #HI
    lambda0_HI = 1215.6701
    title_HI = "%.1f kpc, %.2f rad" % (b_param, angle)
    eq_width_hi = calc_eq_width(lambda0_HI,HI_tau,redshift, lambda_sp_HI, flux_sp_HI,title_HI)
    intensity = exp(-1 * HI_tau)
    lambda_arr = redshift*lambda0_HI + lambda0_HI

    ax1 = fig1.add_subplot(411)
    ax2 = ax1.twiny()
    ax1.plot(lambda_arr,intensity)
    ax2.set_xlim(ax1.get_xlim())
    ax2.set_xticks(new_tick_locations*lambda0_HI + lambda0_HI)
    ax2.set_xticklabels(tick_function(new_tick_locations*lambda0_HI + lambda0_HI,lambda0_HI))
    ax1.plot(lambda_sp_HI,flux_sp_HI)
    ax1.plot([lambda0_HI, lambda0_HI],ax1.get_ylim(),color='k',linestyle='dashed')
    ax1.text(min(lambda_sp_HI),0.6,title_HI)
    #ax1.set_xlabel("Lambda [Angstrom]")
    ax2.set_xlabel("Velocity [km/s]")
    ax1.text(min(lambda_sp_HI),0.3,"HI 1215: %.2f Angstrom" % eq_width_hi)

    #SiIV
    lambda0_siiv = 1393.755
    title_siiv = "SiIV1394: %.1f kpc, %.2f rad" % (b_param, angle)
    eq_width_siiv = calc_eq_width(lambda0_siiv,SiIV_tau,redshift, lambda_sp_siiv, flux_sp_siiv,title_siiv)
    intensity = exp(-1 * SiIV_tau)
    lambda_arr = redshift*lambda0_siiv + lambda0_siiv
    
    ax1 = fig1.add_subplot(412)
    #ax2 = ax1.twiny()
    ax1.plot(lambda_arr,intensity)
    #ax2.set_xlim(ax1.get_xlim())
    #ax2.set_xticks(new_tick_locations*lambda0_siiv + lambda0_siiv)
    #ax2.set_xticklabels(tick_function(new_tick_locations*lambda0_siiv + lambda0_siiv,lambda0_siiv))
    ax1.plot(lambda_sp_siiv,flux_sp_siiv)
    ax1.plot([lambda0_siiv,lambda0_siiv],ax1.get_ylim(),color='k',linestyle='dashed')
    #ax1.text(min(lambda_sp_siiv),0.6,title_siiv)
    #ax1.set_xlabel("Lambda [Angstrom]")
    #ax2.set_xlabel("Velocity [km/s]")
    #ax1.set_ylabel("Normalized Flux")
    ax1.text(min(lambda_sp_siiv),0.5,"SiIV 1394: %.2f Angstrom" % eq_width_siiv)
    
    #CIV
    lambda0_civ = 1548.195
    title_civ = "CIV1548: %.1f kpc, %.2f rad" % (b_param, angle)
    eq_width_civ = calc_eq_width(lambda0_civ,CIV_tau,redshift, lambda_sp_civ, flux_sp_civ,title_civ)
    intensity = exp(-1 * CIV_tau)
    lambda_arr = redshift*lambda0_civ + lambda0_civ
    
    ax1 = fig1.add_subplot(413)
    #ax2 = ax1.twiny()
    ax1.plot(lambda_arr,intensity)
    #ax2.set_xlim(ax1.get_xlim())
    #ax2.set_xticks(new_tick_locations*lambda0_civ + lambda0_civ)
    #ax2.set_xticklabels(tick_function(new_tick_locations*lambda0_civ + lambda0_civ,lambda0_civ))
    ax1.plot(lambda_sp_civ,flux_sp_civ)
    ax1.plot([lambda0_civ, lambda0_civ],ax1.get_ylim(),color='k',linestyle='dashed')
    #ax1.text(min(lambda_sp_civ),0.6,title_civ)
    #ax1.set_xlabel("Lambda [Angstrom]")
    #ax2.set_xlabel("Velocity [km/s]")
    #ax1.set_ylabel("Normalized Flux")
    ax1.text(min(lambda_sp_civ),0.5,"CIV 1548: %.2f Angstrom" % eq_width_civ)
    
    #OVI
    lambda0_ovi = 1031.927
    title_ovi = "OVI1032: %.1f kpc, %.2f rad" % (b_param, angle)
    eq_width_ovi = calc_eq_width(lambda0_ovi,OVI_tau,redshift, lambda_sp_ovi, flux_sp_ovi,title_ovi)
    intensity = exp(-1 * OVI_tau)
    lambda_arr = redshift*lambda0_ovi + lambda0_ovi
    
    ax1 = fig1.add_subplot(414)
    #ax2 = ax1.twiny()
    ax1.plot(lambda_arr,intensity)
    #ax2.set_xlim(ax1.get_xlim())
    #ax2.set_xticks(new_tick_locations*lambda0_ovi + lambda0_ovi)
    #ax2.set_xticklabels(tick_function(new_tick_locations*lambda0_ovi + lambda0_ovi,lambda0_ovi))
    ax1.plot(lambda_sp_ovi,flux_sp_ovi)
    ax1.plot([lambda0_ovi, lambda0_ovi],ax1.get_ylim(),color='k',linestyle='dashed')
    #ax1.text(min(lambda_sp_ovi),0.6,title_ovi)
    ax1.set_xlabel("Lambda [Angstrom]")
    #ax2.set_xlabel("Velocity [km/s]")
    #ax1.set_ylabel("Normalized Flux")
    ax1.text(min(lambda_sp_ovi),0.5,"OVI 1032: %.2f Angstrom" % eq_width_ovi)

    plt.show()
    plt.savefig(path +  '/spect.' + filebase + '.png')
    #wait = raw_input("PRESS ENTER TO CONTINUE.")
    fig1.clear()
    
    return {'b': b_param,
                     'angle': angle,
                     'eq_width_hi': eq_width_hi,
                     'eq_width_civ': eq_width_civ,                     
                     'eq_width_ovi': eq_width_ovi,
                     'eq_width_siiv': eq_width_siiv,
                     }


if __name__ == '__main__':    
    if (socket.gethostname() == "quirm"):
        prefix = '/home/christenc/Data/Sims/'
    else:
        prefix = '/home/christensen/Storage2/UW/MolecH/Cosmo/'
 
#Find all line of sight files in the directory
    path = prefix + "elektra.cosmo25cmb/elektra.cosmo25cmb.4096g5HbwK1BH/elektra.cosmo25cmb.4096g5HbwK1BH.004096/"
    tfile = "elektra.cosmo25cmb.4096g5HbwK1BH.004096"
    halo = "1"
    halo = "2"
    
    path = prefix + "rogue.cosmo25cmb/rogue.cosmo25cmb.4096g5HbwK1BH/rogue.cosmo25cmb.4096g5HbwK1BH.004096/"
    tfile = "rogue.cosmo25cmb.4096g5HbwK1BH.004096"
    halo = "1"
    halo = "3"
    
    path = prefix + "storm.cosmo25cmb/storm.cosmo25cmb.4096g5HbwK1BH/storm.cosmo25cmb.4096g5HbwK1BH.004096/"
    tfile = "storm.cosmo25cmb.4096g5HbwK1BH.004096"
    halo = "1"
    halo = "2"
    
    #path = prefix + "cptmarvel.cosmo25cmb/cptmarvel.cosmo25cmb.4096g5HbwK1BH/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096/"
    #tfile = "cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096"    
    #halo = "1"
        
    pickle_file=open(path + tfile+"_" + halo + "_CGM" + ".data","wb")
    files = glob.glob(path + "CGM." + halo + "/specaim*kpc")
    for f in files:
        f = (f.split("/"))[-1]
        filebase = (f.split("specaim."))[1]
        print(filebase)
        los_data = analyze_los(path + "CGM." + halo + "/",filebase,path + tfile)
        pickle.dump(los_data,pickle_file, pickle.HIGHEST_PROTOCOL)
    pickle_file.close()

    pickle_file=open(path + tfile+"_" + halo + "_CGMvoight" + ".data","wb")
    files = glob.glob(path + "CGM." + halo + "/specaim*H1216.vpm")
    for f in files:
        f = (f.split("/"))[-1]
        filebase = (f.split("_"))[0] + '_' + (f.split("_"))[1]
        x_pos = float(((filebase.split(".x")[1]).split("kpc"))[0])
        y_pos = float(((filebase.split("_y")[1]).split("kpc"))[0])
        b_param = np.sqrt(x_pos**2 + y_pos**2)
        angle = math.atan2(y_pos,x_pos)
        
        fh = open(path + "CGM." + halo + '/' + filebase + '_H1216.vpm','r')
        line = fh.readlines()
        if len(line) != 0 :
            y = [value for value in line[0].split()]
            N_hi = float(y[1])*10**13
            eq_width_hi = float(y[7])
        fh.close()

        fh = open(path + "CGM." + halo + '/' + filebase + '_CIV1548.vpm','r')
        line = fh.readlines()
        if len(line) != 0 :
            y = [value for value in line[0].split()]
            N_civ = float(y[1])*10**13
            eq_width_civ = float(y[7])
        fh.close()

        fh = open(path + "CGM." + halo + '/' + filebase + '_OVI1032.vpm','r')
        line = fh.readlines()
        if len(line) != 0 :
            y = [value for value in line[0].split()]
            N_ovi = float(y[1])*10**13
            eq_width_ovi = float(y[7])
        fh.close()

        fh = open(path + "CGM." + halo + '/' + filebase + '_SiIV1394.vpm','r')
        line = fh.readlines()
        if len(line) != 0 :
            y = [value for value in line[0].split()]
            N_siiv = float(y[1])*10**13
            eq_width_siiv = float(y[7])
        fh.close()

        #print('eq_width_hi: ',eq_width_hi)
        #print('eq_width_civ: ',eq_width_civ)
        #print('eq_width_ovi: ',eq_width_ovi)
        #print('eq_width_siiv: ',eq_width_siiv)
        
        pickle.dump({'b': b_param,
                     'angle': angle,
                     'eq_width_hi': eq_width_hi,
                     'N_hi': N_hi,                    
                     'eq_width_civ': eq_width_civ,
                     'N_civ': N_civ,                     
                     'eq_width_ovi': eq_width_ovi,
                     'N_ovi': N_ovi,                     
                     'eq_width_siiv': eq_width_siiv,
                     'N_siiv': N_siiv,                     
                     } ,pickle_file, pickle.HIGHEST_PROTOCOL)
    pickle_file.close()     
