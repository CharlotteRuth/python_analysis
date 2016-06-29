#run from ipython command line:
#run /home/christensen/Code/python/phasecompare.py

import pynbody, pickle
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()
fig.subplots_adjust(hspace=0,wspace=0)
gs = matplotlib.gridspec.GridSpec(2,2,width_ratios=[4,1],height_ratios=[1,4])
ax = [plt.subplot(gs[0]),plt.subplot(gs[2]),
      plt.subplot(gs[3])]

nmin, nmax, Tmin, Tmax =1,3.5,1.5,4
prefix = '/home/christensen/Storage1/UW/MolecH/Cosmo/'

names = [prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/h986.cosmo50cmb.3072g14HBWK.00512/h986.cosmo50cmb.3072g14HBWK.00512']

cs = ['r']
labels = ['Fiducial']
for i,n in enumerate(names):
    s = pynbody.load(n)
    h = s.halos()
    pynbody.analysis.halo.center(h[1])
    in2kpc = pynbody.filt.LowPass('r','2 kpc')

    x = np.log10(h[1][in2kpc].gas['rho'].in_units('m_p cm^-3'))
    y = np.log10(h[1][in2kpc].gas['temp'])
    mass = h[1][in2kpc].gas['mass'].in_units('Msol')
    xs = np.linspace(nmin,nmax,500)
    ys = np.linspace(Tmin,Tmax,500)
    good = np.where((x > nmin) & (x<nmax) & (y>Tmin) & (y<Tmax))
    arr = pynbody.plot.util.fast_kde(x[good],y[good],extents=[nmin,nmax,Tmin,Tmax],
                                     gridsize=(500,500),
                                     weights=mass[good])
    levels = np.logspace(3,6,6)
    contp = ax[1].contour(xs,ys,arr,levels=levels,colors=[cs[i]]*10)
    #plt.clabel(contp,fmt='%.1e',inline_space=1)
 
    ax[0].hist(x[good],weights=mass[good],range=[nmin,nmax],bins=100,
               histtype='step',color=cs[i],log=True)
    
    ax[2].hist(y[good],weights=mass[good],range=[Tmin,Tmax],bins=100,log=True,
               histtype='step',orientation='horizontal',color=cs[i])


ax[1].set_xlabel('log$_{10}$(n [cm$^{-3}$])')
ax[1].set_ylabel(r'log$_{10}$(Temperature [K])')
ax[1].set_xlim(nmin,nmax)
ax[0].set_xlim(nmin,nmax)
ax[1].set_ylim(Tmin,Tmax-0.1)
ax[2].set_ylim(Tmin,Tmax-0.1)
ax[0].set_ylabel('Mass [M$_\odot$]')
ax[2].set_xlabel('Mass [M$_\odot$]')
plt.setp(ax[0].get_xticklabels(),visible=False)
plt.setp(ax[2].get_yticklabels(),visible=False)
ax[2].set_xlim(9e4,2e8)
#ax[2].xaxis.set_major_locator(matplotlib.ticker.FixedLocator([1e9,2e9,3e9]))
#plt.colorbar(im).set_label('log_{10}(Mass)')
#plt.legend(loc=0)


plt.savefig('2kpcphasecompare.eps')
