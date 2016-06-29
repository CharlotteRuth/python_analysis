#ipython --pylab

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pynbody, sys
from pynbody.analysis import profile, angmom, halo
from pynbody import units, config
import pynbody.filt as f

from pynbody.analysis import luminosity as lum
import os, glob

XSOLfe=0.117E-2
XSOLO=0.96E-2
XSOLH=0.706

#halo_num = '6'
#tfile = '/home/christensen/Storage2/UW/MolecH/Cosmo/h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g14HBWK/steps/h799.cosmo25cmb.3072g14HBWK_2.00512.dir/h799.cosmo25cmb.3072g14HBWK.00512.halo.' + halo_num
#tfile_name = 'h799.cosmo25cmb.3072g14HBWK.00512.halo.' + halo_num

#halo_num = '1'
#tfile = '/home/christensen/Storage2/UW/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK_2.00512.dir/h516.cosmo25cmb.3072g14HBWK.00512.halo.' + halo_num
#tfile_name = 'h516.cosmo25cmb.3072g14HBWK.00512.halo.' + halo_num

#1,2,3,8,15,16
halo_num = '1'
tfile = '/home/christensen/Storage2/UW/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/h986.cosmo50cmb.3072g14HBWK_2.00512.dir/h986.cosmo50cmb.3072g14HBWK.00512.halo'#.' + halo_num
tfile_name = 'h986.cosmo50cmb.3072g14HBWK'#.halo.' + halo_num

#1,2,3
halo_num = '1'
tfile = '/home/christensen/Storage2/UW/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.00512.dir/h603.cosmo50cmb.3072g14HBWK.00512'#.halo.' + halo_num
tfile_name = 'h603.cosmo50cmb.3072g14HBWK.00512'#.halo.' + halo_num


#h1 =  pynbody.load(tfile)
#h1.physical_units()
s = pynbody.load(filename)
h = s.halos()
h.physical_units()
h7 = h.load_copy(7)
halo = h7

#plt.scatter(halo.star['age'].in_units('Gyr'),np.log10(halo.star['OxMassFrac'])- np.log10(XSOLO/XSOLH))
#plt.axis([0, 14, -4, 0])
#xlabel('Age [Gyr]')
#plt.show()
#plt.scatter(halo.star['age'].in_units('Gyr'),np.log10(halo.star['FeMassFrac'])- np.log10(XSOLfe/XSOLH))
#plt.axis([0, 14, -4, 0])
#xlabel('Age [Gyr]')
#plt.show()

nbins = 50
xmin = -1e-4
xmax = 14

dt = 1.0
agebin = arange(xmin,xmax,dt) + dt/2.0
meanOx = empty(size(agebin))
stdevOx = empty(size(agebin))
meanFe = empty(size(agebin))
stdevFe = empty(size(agebin))

stars = halo.star
stars[f.LowPass('OxMassFrac',1e-7)].star['OxMassFrac'] = 1e-7
stars[f.LowPass('FeMassFrac',1e-8)].star['FeMassFrac'] = 1e-8

for i in range(size(agebin)):
    meanOx[i] = np.mean(np.log10(stars[f.BandPass('age',agebin[i] - dt/2,agebin[i] + dt/2)].star['OxMassFrac']) - np.log10(XSOLO/XSOLH))
    stdevOx[i] = np.std(np.log10(stars[f.BandPass('age',agebin[i] - dt/2,agebin[i] + dt/2)].star['OxMassFrac']) - np.log10(XSOLO/XSOLH))
    meanFe[i] = np.mean(np.log10(stars[f.BandPass('age',agebin[i] - dt/2,agebin[i] + dt/2)].star['FeMassFrac']) - np.log10(XSOLfe/XSOLH))
    stdevFe[i] = np.std(np.log10(stars[f.BandPass('age',agebin[i] - dt/2,agebin[i] + dt/2)].star['FeMassFrac']) - np.log10(XSOLfe/XSOLH))

ymin = -1.85
ymax = 0 #-0.85
ymin = -2.2

aspectratio = (xmax - xmin)/(ymax - ymin)*0.7
binnedOxMassFrac, xedges, yedges = np.histogram2d(np.log10(halo.star['OxMassFrac'])- np.log10(XSOLO/XSOLH), halo.star['age'].in_units('Gyr'), bins = [nbins,nbins], range = [[ymin,ymax],[xmin,xmax]])
x = linspace(xmin, xmax, num = nbins + 1)
y = linspace(ymin, ymax, num = nbins + 1)
xcenter = (x[0:-1]+x[1:])/2.0
ycenter = (y[0:-1]+y[1:])/2.0
plt.imshow(log10(binnedOxMassFrac), cmap="Blues", extent = [xmin, xmax, ymin, ymax], interpolation = 'nearest', origin = 'lower',aspect = aspectratio)
#plot(agebin,meanOx,'ko')
errorbar(agebin,meanOx, xerr=dt/2, yerr=stdevOx, fmt ='ko')
#plt.contour(xcenter, ycenter, log10(binnedOxMassFrac))
axis([xmin,xmax,ymin,ymax])
xlabel('Age [Gyr]')
ylabel('[Ox/H]')
title(tfile_name)
plt.show()
plt.savefig(tfile+'.OxVsAge_log.png')

plt.close()

ymin = -1.65
ymax = 0 #-0.65
ymin = -2

aspectratio = (xmax - xmin)/(ymax - ymin)*0.7
binnedFeMassFrac, xedges, yedges = np.histogram2d(np.log10(halo.star['FeMassFrac'])- np.log10(XSOLfe/XSOLH), halo.star['age'].in_units('Gyr'), bins = [nbins,nbins], range = [[ymin,ymax],[xmin,xmax]])
x = linspace(xmin, xmax, num = nbins + 1)
y = linspace(ymin, ymax, num = nbins + 1)
xcenter = (x[0:-1]+x[1:])/2.0
ycenter = (y[0:-1]+y[1:])/2.0
plt.imshow(log10(binnedFeMassFrac), cmap="Blues", extent = [xmin, xmax, ymin, ymax], interpolation = 'nearest', origin = 'lower',aspect = aspectratio)
errorbar(agebin,meanFe, xerr=dt/2, yerr=stdevFe, fmt ='ko')
#plt.contour(xcenter, ycenter, log10(binnedOxMassFrac))
xlabel('Age [Gyr]')
ylabel('[Fe/H]')
title(tfile_name)
axis([xmin,xmax,ymin,ymax])
plt.show()
plt.savefig(tfile+'.FeVsAge_log.png')

plt.close()
