import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pynbody, sys
from pynbody.analysis import profile, angmom, halo
from pynbody import units, config
import pynbody.filt as f

from pynbody.analysis import luminosity as lum
import os, glob

#def make_rs(im):
#    xsize, ysize = np.shape(im)
#    x = np.arange(-xsize/2, xsize/2)
#    y = np.arange(-ysize/2, ysize/2)
#    xs, ys = np.meshgrid(x,y)
#    return 2.0*np.sqrt(xs**2 + ys**2)

tfile = sys.argv[1]

fig = plt.figure(figsize=(8.,8.))
#fig.subplots_adjust(left=0.08, bottom=0.05, right=0.9, top=0.97, wspace=0.12, 
#                    hspace=0.13)

#sps = [fig.add_subplot(3,2,1), fig.add_subplot(3,2,2), fig.add_subplot(3,2,3),
#       fig.add_subplot(3,2,4), fig.add_subplot(3,2,5), fig.add_subplot(3,2,6)]

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
                   vmin=12,vmax=20)
hfbhiim.set_xlabel('x [kpc]')
hfbhiim.set_ylabel('y [kpc]')

plt.savefig(tfile+'.mapprofHI.png')

