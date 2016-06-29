import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pynbody, sys
import pynbody.plot.sph as sph
from pynbody.analysis import profile, angmom, halo
from pynbody import units, config
import pynbody.filt as f

from pynbody.analysis import luminosity as lum
import os, glob

#tfile = '/home/christensen/Storage2/UW/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/h986.cosmo50cmb.3072g14HBWK_2.00512.dir/h986.cosmo50cmb.3072g14HBWK.00512.halo.1'
tfile = '/home/christensen/Storage2/UW/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/h986.cosmo50cmb.3072g14HBWK_2.00512.dir/h986.cosmo50cmb.3072g14HBWK.00512'
#tfile = '/home/christensen/Storage2/UW/MolecH/Cosmo/h239.cosmo50cmb.3072g/h239.cosmo50cmb.3072g14HMbwK/steps/h239.cosmo50cmb.3072g14HMbwK.00512.dir/h239.cosmo50cmb.3072g14HMbwK.00512.halo.1'

#h1 =  pynbody.load(tfile)
#h1.physical_units()
s =  pynbody.load(tfile)
s.physical_units()
h = s.halos()
h1 = h[1]
pynbody.analysis.halo.center(h1,mode='hyb')
pynbody.analysis.angmom.sideon(h1)
hrvir = np.max(h1.gas['r'])
xymax = hrvir/2**(0.5)

#-------------------------- Read data on outflows -------------------------
outflow1_file = open('/home/christensen/Storage2/UW/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/disk_.1.txt')
outflow1 = outflow1_file.readlines()
outflow1_file.close()
outflow1 = map(long, outflow1)
outflow1 = np.array(outflow1)
indicies = np.in1d(s['iord'],outflow1)
outflow1_part = s[np.nonzero(indicies)]

outflow2_file = open('/home/christensen/Storage2/UW/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/disk_1.txt')
outflow2 = outflow2_file.readlines()
outflow2_file.close()
outflow2 = map(long, outflow2)
outflow2 = np.array(outflow2)
indicies = np.in1d(s['iord'],outflow2)
outflow2_part = s[np.nonzero(indicies)]

outflow3_file = open('/home/christensen/Storage2/UW/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/disk_2.txt')
outflow3 = outflow3_file.readlines()
outflow3_file.close()
outflow3 = map(long, outflow3)
outflow3 = np.array(outflow3)
indicies = np.in1d(s['iord'],outflow3)
outflow3_part = s[np.nonzero(indicies)]

#sph.image(outflow1_part.gas,units='m_p cm^-2', width=5*xymax)

#------------------------------------------------------
width = 3*xymax

fig = plt.figure(figsize=(7.,2.5))
fig.subplots_adjust(left=0.08, bottom=0.05, right=0.9, top=0.97, wspace=0.12, 
                    hspace=0.13)
gs = matplotlib.gridspec.GridSpec(1, 4, width_ratios=[4,4,4,1.4]) 
sps = [fig.add_subplot(gs[0]), fig.add_subplot(gs[1]), fig.add_subplot(gs[2]), fig.add_subplot(gs[3])] 
ax1 = (sps[0])
sim = sph.image(outflow1_part.gas,
                   units='m_p cm^-2', width=width,  clear=False,
                            title='0.1 Gyr ago',show_cbar=True,vmin=1e16,vmax=1e22,
                   subplot=sps[0],cmap='ocean_r');
for label in ax1.yaxis.get_ticklabels():
    label.set_fontsize(8)

for label in ax1.xaxis.get_ticklabels():
    label.set_fontsize(8)

ax1.set_title(ax1.get_title(),fontsize = 12)
circle1=plt.Circle((0,0),xymax,ec='black',fc='none')
ax1.add_artist(circle1)

#Ejected 1 Myr ago
ax2 = subplot(sps[1],sharey = ax1)
#ax2.set_fontsize(12)
sim = sph.image(outflow2_part.gas,
                   units='m_p cm^-2', width=width,  clear=False,
                            title='1 Gyr ago',show_cbar=False,
                   subplot=sps[1],vmin=1e16,vmax=1e22,cmap='ocean_r');
setp(ax2.get_yticklabels(), visible=False)
setp(ax2.set_ylabel(''))
for label in ax2.xaxis.get_ticklabels():
    label.set_fontsize(8)

ax2.set_title(ax2.get_title(),fontsize = 12)
circle2=plt.Circle((0,0),xymax,ec='black',fc='none')
ax2.add_artist(circle2)

#Ejected 2 Myr ago
ax3 = subplot(sps[2],sharey = ax1)
sim = sph.image(outflow3_part.gas,
                   units='m_p cm^-2', width=width,  clear=False,
                            title='2 Gyr ago',show_cbar=False,
                   subplot=sps[2],vmin=1e16,vmax=1e22,cmap='ocean_r');

setp(ax3.get_yticklabels(), visible=False)
setp(ax3.set_ylabel(''))
for label in ax3.xaxis.get_ticklabels():
    label.set_fontsize(8)

sps[3].axis('off')
ax3.set_title(ax3.get_title(),fontsize = 12)
circle3=plt.Circle((0,0),xymax,ec='black',fc='none')
ax3.add_artist(circle3)
show()
savefig(tfile + '.outflow_hist.png')

#------------------------------------------------
levels = np.array((1e16,1e17,1e18,1e19,1e20,1e21,1e22))
width = 3*xymax
circle=plt.Circle((0,0),xymax,ec='black',fc='none')

fig = plt.figure(figsize=(7.,2.5))
fig.subplots_adjust(left=0.08, bottom=0.05, right=0.9, top=0.97, wspace=0.12, 
                    hspace=0.13)
#gs = matplotlib.gridspec.GridSpec(1, 4, width_ratios=[4,4,4,1.4]) 
gs = matplotlib.gridspec.GridSpec(1, 3, width_ratios=[4,4,4])
sps = [fig.add_subplot(gs[0]), fig.add_subplot(gs[1]), fig.add_subplot(gs[2])]#, fig.add_subplot(gs[3])] 

#Ejected 0.1 Myr ago
ax1 = (sps[0])
#h1hiim = sph.image(s.gas,qty='hiden',
#                   units='m_p cm^-2', width=width,  clear=False,
#                            title='N$_{HI}$ cm$^{-2}$',
#                   vmin=1e12,vmax=1e20);
sim = sph.image(s.gas,
                   units='m_p cm^-2', width=width,  clear=False,
                            title='N$_{gas}$',show_cbar=False,vmin=1e16,vmax=1e22,cmap='Greys'
                   ,subplot=sps[0]);

im = sph.image(outflow1_part.gas,units='m_p cm^-2',noplot = True)
res = im.shape
x,y = np.meshgrid(np.linspace(-width.in_units('cm')/2,width.in_units('cm')/2,res[0]),np.linspace(-width.in_units('cm')/2,width.in_units('cm')/2,res[0]))
show()
plt.contour(x,y,im,levels=levels,colors = 'r')
#sph.contour(outflow1_part.gas,units='m_p cm^-2', width=width, levels=levels,subplot=sps[0])
for label in ax1.yaxis.get_ticklabels():
    label.set_fontsize(8)

ax1.set_title(ax1.get_title(),fontsize = 12)
ax1.add_artist(circle)

#Ejected 1 Myr ago
ax2 = subplot(sps[1],sharey = ax1)
ax2.set_fontsize(12)
#h1hiim = sph.image(s.gas,qty='hiden',
#                   units='m_p cm^-2', width=width,  clear=False,
#                            title='N$_{HI}$ cm$^{-2}$',
#                   vmin=1e12,vmax=1e20);
sim = sph.image(s.gas,
                   units='m_p cm^-2', width=width,  clear=False,
                            title='N$_{gas}$',show_cbar=False,
                   subplot=sps[1],vmin=1e16,vmax=1e22,cmap='Greys');

sph.contour(outflow2_part.gas,units='m_p cm^-2', width=width, levels=levels,subplot=sps[1])
setp(ax2.get_yticklabels(), visible=False)
setp(ax2.set_ylabel(''))
for label in ax2.xaxis.get_ticklabels():
    label.set_fontsize(8)

ax2.set_title(ax2.get_title(),fontsize = 12)
ax2.add_artist(circle)

#Ejected 2 Myr ago
ax3 = subplot(sps[2],sharey = ax1)
#h1hiim = sph.image(s.gas,qty='hiden',
#                   units='m_p cm^-2', width=width,  clear=False,
#                            title='N$_{HI}$ cm$^{-2}$',
#                   vmin=1e12,vmax=1e20);
sim = sph.image(s.gas,
                   units='m_p cm^-2', width=width,  clear=False,
                            title='N$_{gas}$',show_cbar=False,
                   subplot=sps[1],vmin=1e16,vmax=1e22,cmap='Greys');

sph.contour(outflow3_part.gas,units='m_p cm^-2', width=width, levels=levels,subplot=sps[2])
setp(ax3.get_yticklabels(), visible=False)
setp(ax3.set_ylabel(''))
for label in ax3.xaxis.get_ticklabels():
    label.set_fontsize(8)

sps[3].axis('off')
ax3.set_title(ax3.get_title(),fontsize = 12)
ax3.add_artist(circle)
show()

#---------------------------------------------------------

sph.velocity_image(outflow1_part.gas,vector_color="cyan",qty="temp", width=xymax,cmap="YlOrRd", approximate_fast=False,fill_nan=True,fill_val = 0.0)

sph.image(s.gas,width=500)


#HI map
hiif = pynbody.analysis.ionfrac.calculate(s.gas,ion='hi')
s.gas['hiden'] = s.gas['rho']*s.gas['hydrogen']*hiif
s.gas['hiden'].units = s.gas['rho'].units
#shiim = sph.image(s.gas,qty='hiden',
#                   units='m_p cm^-2', width=5*xymax,  clear=False,
#                            title='N$_{HI}$ cm$^{-2}$',
#                   vmin=1e12,vmax=1e20);
sim = sph.image(s.gas,
                   units='m_p cm^-2', width=5*xymax,  clear=False,
                            title='N$_{gas}$ cm$^{-2}$',
                   vmin=1e18,vmax=1e22);

sph.velocity_image(s.gas,vector_color="cyan",qty="hiden",units='m_p cm^-2',width=5*xymax, approximate_fast=False,fill_nan=True,fill_val = 0.0)
