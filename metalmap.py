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

tfile = '/home/christensen/Storage2/UW/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/h986.cosmo50cmb.3072g14HBWK_2.00512.dir/h986.cosmo50cmb.3072g14HBWK.00512.halo.1'
#tfile = '/home/christensen/Storage2/UW/MolecH/Cosmo/h239.cosmo50cmb.3072g/h239.cosmo50cmb.3072g14HMbwK/steps/h239.cosmo50cmb.3072g14HMbwK.00512.dir/h239.cosmo50cmb.3072g14HMbwK.00512.halo.1'
#tfile = sys.argv[1]
h1 =  pynbody.load(tfile)
h1.physical_units()
pynbody.analysis.angmom.sideon(h1)
hrvir = np.max(h1.gas['r'])
xymax = hrvir/2**(0.5)

#-------------------------------------------------------------------
fig = plt.figure(figsize=(7.,2.5))
fig.subplots_adjust(left=0.08, bottom=0.05, right=0.9, top=0.97, wspace=0.12, 
                    hspace=0.13)
gs = matplotlib.gridspec.GridSpec(1, 4, width_ratios=[4,4,4,1.4]) 
sps = [fig.add_subplot(gs[0]), fig.add_subplot(gs[1]), fig.add_subplot(gs[2]), fig.add_subplot(gs[3])] 
ax1 = (sps[0])
#fig,axes = plt.subplots(sps[0],sharey=True)
#fig.tight_layout()

#HI map
hiif = pynbody.analysis.ionfrac.calculate(h1.gas,ion='hi')
h1.gas['hiden'] = h1.gas['rho']*h1.gas['hydrogen']*hiif
h1.gas['hiden'].units = h1.gas['rho'].units
h1hiim = pynbody.plot.image(h1.gas,qty='hiden',
                   units='m_p cm^-2', width=xymax,  clear=False,
                            title='log$_{10}$ N$_{HI}$ cm$^{-2}$',
                   subplot=sps[0],vmin=12,vmax=20);
for label in ax1.yaxis.get_ticklabels():
    label.set_fontsize(8)

for label in ax1.xaxis.get_ticklabels():
    label.set_fontsize(8)

ax1.set_title(ax1.get_title(),fontsize = 12)

#OVI map
ax2 = subplot(sps[1],sharey = ax1)
ax2.set_fontsize(12)
oviif = pynbody.analysis.ionfrac.calculate(h1.gas)
h1.gas['oviden'] = h1.gas['rho']*h1.gas['OxMassFrac']*oviif
h1.gas['oviden'].units = h1.gas['rho'].units
h1oviim = pynbody.plot.image(h1.gas,qty='oviden',
                   units='16 m_p cm^-2', width=xymax, 
                             title='log$_{10}$ N$_{OVI}$ cm$^{-2}$',show_cbar=False,
                   subplot=sps[1],vmin=12,vmax=20);
setp(ax2.get_yticklabels(), visible=False)
setp(ax2.set_ylabel(''))
for label in ax2.xaxis.get_ticklabels():
    label.set_fontsize(8)

ax2.set_title(ax2.get_title(),fontsize = 12)

#Si IV map
ax3 = subplot(sps[2],sharey = ax1)
siivden = pynbody.analysis.ionfrac.calculate(h1.gas,ion='siiv')
h1.gas['siivden'] = h1.gas['rho']*(0.00067/0.00117)*h1.gas['FeMassFrac']*siivden #assuming solar O/Mg abundance ratio
h1.gas['siivden'].units = h1.gas['rho'].units
h1mgiiim = pynbody.plot.image(h1.gas,qty='siivden',
                   units='24 m_p cm^-2', width=xymax,  clear=False,
                             title='log$_{10}$ N$_{SiIV}$ cm$^{-2}$',show_cbar=False,
                   subplot=sps[2],vmin=12,vmax=20);
setp(ax3.get_yticklabels(), visible=False)
setp(ax3.set_ylabel(''))
for label in ax3.xaxis.get_ticklabels():
    label.set_fontsize(8)

sps[3].axis('off')
ax3.set_title(ax3.get_title(),fontsize = 12)
show()
savefig(tfile + '.ion3map.png')

#Also possible, CIV, HI, OI, OII, NV, SiIV, MGII

#-----------------------------------------------------------------------
fig = plt.figure(figsize=(8.,12.))
fig.subplots_adjust(left=0.08, bottom=0.05, right=0.9, top=0.97, wspace=0.12, 
                    hspace=0.13)

sps = [fig.add_subplot(3,2,1), fig.add_subplot(3,2,2), fig.add_subplot(3,2,3),
       fig.add_subplot(3,2,4), fig.add_subplot(3,2,5), fig.add_subplot(3,2,6)] #, fig.add_subplot(3,2,7), fig.add_subplot(3,2,8)]

#h1g = pynbody.plot.image(h1.g, width=xymax, cmap='Blues');

#Total Hydrogen Map
h1.gas['h'] = h1.gas['rho']*h1.gas['hydrogen']
h1.gas['h'].units = h1.gas['rho'].units
hfboxim=pynbody.plot.image(h1.gas,qty='h',
                           units='m_p cm^-2', width=xymax, clear=False,
                           qtytitle='log$_{10}$ N$_{H}$ cm$^{-2}$',
                           title = 'log$_{10}$ N$_{H}$ cm$^{-2}$',
                           subplot=sps[0],vmin=12,vmax=20);

#Total Oxygen Map
h1.gas['oxden'] = h1.gas['rho']*h1.gas['OxMassFrac']
h1.gas['oxden'].units = h1.gas['rho'].units
h1oxim=pynbody.plot.image(h1.gas,qty='oxden',
                           units='16 m_p cm^-2', width=xymax, clear=False,
                           qtytitle='log$_{10}$ N$_{O}$ cm$^{-2}$', show_cbar=False,
                          title = 'log$_{10}$ N$_{O}$ cm$^{-2}$',
                           subplot=sps[2],vmin=12,vmax=20);

#Total Iron Map
h1.gas['feden'] = h1.gas['rho']*h1.gas['FeMassFrac']
h1.gas['feden'].units = h1.gas['rho'].units
h1feim=pynbody.plot.image(h1.gas,qty='feden',
                           units='56 m_p cm^-2', width=xymax, clear=False, 
                           qtytitle='log$_{10}$ N$_{Fe}$ cm$^{-2}$',show_cbar=False,
                          title = 'log$_{10}$ N$_{Fe}$ cm$^{-2}$',
                           subplot=sps[4],vmin=12,vmax=20);

#HI map
hiif = pynbody.analysis.ionfrac.calculate(h1.gas,ion='hi')
h1.gas['hiden'] = h1.gas['rho']*h1.gas['hydrogen']*hiif
h1.gas['hiden'].units = h1.gas['rho'].units
h1hiim = pynbody.plot.image(h1.gas,qty='hiden',
                   units='m_p cm^-2', width=xymax,  clear=False,
                            title='log$_{10}$ N$_{HI}$ cm$^{-2}$',show_cbar=False,
                   subplot=sps[1],vmin=12,vmax=20);

#OVI map
oviif = pynbody.analysis.ionfrac.calculate(h1.gas)
h1.gas['oviden'] = h1.gas['rho']*h1.gas['OxMassFrac']*oviif
h1.gas['oviden'].units = h1.gas['rho'].units
h1oviim = pynbody.plot.image(h1.gas,qty='oviden',
                   units='16 m_p cm^-2', width=xymax, 
                             title='log$_{10}$ N$_{OVI}$ cm$^{-2}$',show_cbar=False,
                   subplot=sps[3],vmin=12,vmax=20);

#MgII map
mgiiif = pynbody.analysis.ionfrac.calculate(h1.gas,ion='mgii')
h1.gas['mgiiden'] = h1.gas['rho']*0.08128*h1.gas['OxMassFrac']*mgiiif #assuming solar O/Mg abundance ratio
h1.gas['mgiiden'].units = h1.gas['rho'].units
h1mgiiim = pynbody.plot.image(h1.gas,qty='mgiiden',
                   units='24 m_p cm^-2', width=xymax,  clear=False,
                             title='log$_{10}$ N$_{MgII}$ cm$^{-2}$',show_cbar=False,
                   subplot=sps[5],vmin=12,vmax=20);

#MgIV map
#mgivif = pynbody.analysis.ionfrac.calculate(h1.gas,ion='mgiv')
#h1.gas['mgivden'] = h1.gas['rho']*0.08128*h1.gas['OxMassFrac']*mgivif #assuming solar O/Mg abundance ratio
#h1.gas['mgivden'].units = h1.gas['rho'].units
#h1mgivim = pynbody.plot.image(h1.gas,qty='mgivden',
#                   units='24 m_p cm^-2', width=xymax, 
#                   vmin=12,vmax=20);

plt.savefig(tfile + '.ionmap.png')


fig = plt.figure(figsize=(8.,8.))
#fig.subplots_adjust(left=0.08, bottom=0.05, right=0.9, top=0.97, wspace=0.12, 
#                    hspace=0.13)

sps = [fig.add_subplot(2,2,1), fig.add_subplot(2,2,2), fig.add_subplot(2,2,3)]

#Total Mass Map
h1mim=pynbody.plot.image(h1.gas,qty='rho',
                           units='Msol pc^-2', width=xymax,  clear=False,
                           qtytitle='log$_{10}$ M Msol pc$^{-2}$',show_cbar=False,
                         title ='log$_{10}$ M Msol pc$^{-2}$',
                         subplot=sps[0],vmin = -5,vmax=1);

#Total Metal Map
h1.gas['zden'] = h1.gas['rho']*h1.gas['FeMassFrac']*1.06 + h1.gas['rho']*h1.gas['OxMassFrac']*2.09
h1.gas['zden'].units = h1.gas['rho'].units
h1zim=pynbody.plot.image(h1.gas,qty='zden',
                           units='Msol pc^-2', width=xymax,  clear=False,
                           qtytitle='log$_{10}$ M$_Z$ Msol pc$^{-2}$',show_cbar=False,
                         title='log$_{10}$ M$_Z$ Msol pc$^{-2}$',
                         subplot=sps[1],vmin = -5,vmax=1);

#Metallicity Fraction
h1.gas['zfrac'] = h1.gas['zden']/h1.gas['rho']/0.0130215
h1zmim=pynbody.plot.image(h1.gas,qty='zfrac',
                          av_z ='rho',
                           width=xymax, clear=False,
                           qtytitle='log$_{10}$ Z/Z$_{solar}$',
                         subplot=sps[2],vmin = -2,vmax=0);

plt.savefig(tfile+'.zmap.png')
