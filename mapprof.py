import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pynbody, sys
from pynbody.analysis import profile, angmom, halo
from pynbody import units, config
import pynbody.filt as f

from pynbody.analysis import luminosity as lum
import os, glob

def make_rs(im):
    xsize, ysize = np.shape(im)
    x = np.arange(-xsize/2, xsize/2)
    y = np.arange(-ysize/2, ysize/2)
    xs, ys = np.meshgrid(x,y)
    # 2.0 because using 2 kpc pixels
    return 2.0*np.sqrt(xs**2 + ys**2)

tfile = sys.argv[1]

fig = plt.figure(figsize=(8.,12.))
fig.subplots_adjust(left=0.08, bottom=0.05, right=0.9, top=0.97, wspace=0.12, 
                    hspace=0.13)

sps = [fig.add_subplot(3,2,1), fig.add_subplot(3,2,2), fig.add_subplot(3,2,3),
       fig.add_subplot(3,2,4), fig.add_subplot(3,2,5), fig.add_subplot(3,2,6)]

hfb = pynbody.load(tfile)
h = hfb.halos()
hfbsmass = np.sum(h[1].stars['mass'].in_units('Msol'))
#hfblstar = 10.0**(0.4*(-21 - pynbody.analysis.luminosity.halo_mag(h[1])))
hfb.physical_units()
pynbody.analysis.angmom.faceon(h[1])
hfbrvir = np.max(h[1].gas['r'])
notdiskf = f.Not(f.Disc('40 kpc','3 kpc'))

# Upper left:  HI map
#hfb.gas['hiden'] = hfb.gas['rho']*hfb.gas['HI']
hiif = pynbody.analysis.ionfrac.calculate(hfb.gas,ion='hi')
hfb.gas['hiden'] = hfb.gas['rho']*hfb.gas['hydrogen']*hiif
hfb.gas['hiden'].units = hfb.gas['rho'].units
hfbhiim = pynbody.plot.image(hfb.gas[notdiskf],qty='hiden', clear=False,
                   units='m_p cm^-2', width=400, show_cbar=False, 
                   subplot=sps[0],vmin=12,vmax=20)
sps[0].set_xlabel('x [kpc]')
sps[0].set_ylabel('y [kpc]')

# Upper right:  OVI map
oviif = pynbody.analysis.ionfrac.calculate(hfb.gas)
hfb.gas['oviden'] = hfb.gas['rho']*hfb.gas['OxMassFrac']*oviif
hfb.gas['oviden'].units = hfb.gas['rho'].units
hfboviim = pynbody.plot.image(hfb.gas[notdiskf],qty='oviden', clear=False,
                   units='16 m_p cm^-2', width=1000,
                   subplot=sps[1],vmin=12,vmax=17)
sps[1].set_xlabel('x [kpc]')



obsrho=np.array([265, 126, 255, 181, 156, 124, 155, 130, 140, 87, 294, 64, 
                 83, 209, 221, 32])
obsnhi=np.array([13.8,15.14,14.64,13.46, 14.04, 13.64, 14.4, 15.67,
                 15.73, 14.8, 14.83, 15.35, 15.29, 15.01, 14.3, 15.11])
obsnovi = np.array([13.74, 14.9, 13.71, 13.6, 13.35, 13.77, 14.08, 14.3, 
                    14.16, 14.7, 14.38, 13.75, 14.16, 13.93, 13.61, 13.1])

# Lower left: HI radial profile
#import ipdb; ipdb.set_trace()
rs = make_rs(hfboviim)

pynbody.plot.hist2d(np.log10(rs.flatten() + 1e-6),hfbhiim.flatten(),xrange=[0.3,3],
                    yrange=[12,19], subplot=sps[2], cmap=plt.cm.gist_yarg)
#pynbody.plot.hist2d(rs.flatten(),hfbhiim.flatten(),xrange=[0.3,3],
#                    yrange=[12,19], subplot=sps[2], cmap=plt.cm.gist_yarg)
sps[2].plot(np.log10(obsrho),obsnhi,'o', label='z=0 0.1 L* < L < L* Prochaska et al (2011)')
sps[2].set_xlabel(r'log$_{10}$ ($\rho$ [kpc])')
sps[2].plot(np.log10([hfbrvir,hfbrvir]),[12,19],label=r'$r_{vir}$')
sps[2].legend(loc=0,prop=matplotlib.font_manager.FontProperties(size='small'))
sps[2].set_ylim(12,19)
sps[2].set_xlim(0.3,3)
sps[2].set_xlabel(r'log$_{10}$($\rho$ [kpc])')

# Lower right: OVI radial profile
pynbody.plot.hist2d(np.log10(rs.flatten() + 1e-6),hfboviim.flatten(),xrange=[0.3,3],
                    yrange=[10,17], subplot=sps[3], cmap=plt.cm.gist_yarg)
#pynbody.plot.hist2d(rs.flatten(),hfboviim.flatten(),xrange=[0.3,3],
#                    yrange=[10,17], subplot=sps[3], cmap=plt.cm.gist_yarg)
sps[3].plot(np.log10(obsrho),obsnovi,'o', label='0.1 L* < L < L* Prochaska et al (2011)')
#hfboxmean, hfboxmin, hfboxmax, hfboxbins = pynbody.plot.sph.image_radial_profile(hfboxim)
#sps[3].plot(np.log10(hfboxbins*2), np.log10(hfboxmean), '-', label='total Ox mean')
sps[3].plot(np.log10([hfbrvir,hfbrvir]),[10,17],label=r'$r_{vir}$')
sps[3].set_ylim(10,17)
sps[3].set_xlim(0.3,3)
sps[3].set_xlabel(r'log$_{10}$($\rho$ [kpc])')

plt.setp(sps[1].get_yticklabels(), visible=False)
sps[2].set_ylabel('log$_{10}$ (N[cm$^{-2}$])')
sps[0].set_title('HI')
sps[1].set_title('OVI')
#sps[3].text(0.7,0.9,'L = %4.2f L*'%(hfblstar),transform=sps[3].transAxes)
sps[3].text(0.05,0.9,'z = %4.2f'%(hfb.properties['z']),transform=sps[3].transAxes)

# Total Mass weighted phase diagram
totphaseim = pynbody.plot.hist2d(
    np.log10(h[1].gas['rho'].in_units('m_p cm**-3')), 
    np.log10(h[1].gas['temp']), ret_im=True, scalemin=5,scalemax=8.5, 
    weights=h[1].gas['mass'].in_units('Msol'),
    subplot=sps[4], clear=False, xrange=[-6,2.5],
    yrange=[3,7])
sps[4].set_xlabel('n [cm$^{-3}$]')
sps[4].set_ylabel('Temperature [K]')

hiax = fig.add_axes([0.38,0.07,0.02,0.25])
cb1 = fig.colorbar(totphaseim,cax=hiax)
cb1.set_label('Mass [M$_\odot$])')

# OVI mass weighted phase diagram
hfb.gas['oviif'] = pynbody.analysis.ionfrac.calculate(hfb.gas)
logrho = np.log10(h[1].gas['rho'].in_units('m_p cm**-3'))
logtemp = np.log10(h[1].gas['temp'])
oviphaseim = pynbody.plot.hist2d(
    logrho, logtemp, ret_im=True, scalemin=-3,scalemax=4, 
    weights=h[1].gas['OxMassFrac']*h[1].gas['oviif']*h[1].gas['mass'].in_units('Msol'), 
    subplot=sps[5], clear=False, xrange=[-6,2.5], yrange=[3,7])
sps[5].set_xlabel('n [cm$^{-3}$]')

oviax = fig.add_axes([0.90,0.07,0.02,0.25])
cb2 = fig.colorbar(oviphaseim,cax=oviax)
cb2.set_label('Mass O$_{VI}$ [M$_\odot$])')
#np.savez(tfile+'.phaseims.npz',oviphaseim, totphaseim)

plt.savefig(tfile+'.mapprof.png')

dum = plt.figure()
dumsp = fig.add_subplot(1,1,1)

# Not plotted, but used later:  Total Oxygen Map
hfb.gas['oxden'] = hfb.gas['rho']*hfb.gas['OxMassFrac']
hfb.gas['oxden'].units = hfb.gas['rho'].units
hfboxim=pynbody.plot.sideon_image(hfb.gas,qty='oxden',
                                  units='16 m_p cm^-2', width=1000, 
                                  center=False, subplot=dumsp, 
                                  title='z = %4.2f'%hfb.properties['z'],
                                  qtytitle='log$_{10}$ N$_{O}$',noplot=True)


#np.savez(tfile+'.ims.npz',hfboviim, hfbhiim, hfbsmass,hfblstar, hfbrvir, hfboxim)
np.savez(tfile+'.ims.npz',hfboviim, hfbhiim, hfbsmass,hfbrvir, hfboxim)

