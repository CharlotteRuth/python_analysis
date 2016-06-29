import numpy as np, sys
import matplotlib.pyplot as plt
import matplotlib
import pynbody, sys, math
from pynbody.analysis import profile, angmom, halo
from pynbody import units, config
import pynbody.filt as f
import diskfitting

from pynbody.analysis import luminosity as lum
import os, glob, pickle

rmin = '6 kpc'
rmax = '10 kpc'
zmin = '0 kpc'
zmax = '3 kpc'

annulus = f.Disc(rmax,zmax) & ~f.Disc(rmin,zmax) & ~f.Disc(rmax,zmin)
disk = pynbody.filt.Disc('30 kpc','3 kpc')
bins = 50

s = pynbody.load(sys.argv[1])
h = s.halos()
pynbody.analysis.angmom.faceon(h[1])

diskstars = h[1].star[disk]
annstars = h[1].star[annulus]
n_left = len(annstars)
starsperbin = n_left / bins
print "stars per bin: "+str(starsperbin)

hzs = []
rexps = []
hzerrs = []
hrerrs = []
ages = []
ofes = []
fehs = []
sigmavzs = []
mass = []
vprofzbins = []
vprofmasses = []
vprofn = []
rprofrbins = []
rprofmasses = []
rprofn = []
rprofdr = []
rprofbs = []

tfile = '01024/g1536.01024'
n_done=0
i=0
while n_left>0:
    print i
    n_block = min(n_left,starsperbin)
    n_left -= n_block
    the_particles = annstars[n_done:n_done+n_block]
    #thesestars = diskstars[n_done:n_done+n_block]
    #the_particles = thesestars[annulus]

    fitnum = len(the_particles)

    if fitnum > 100:
        print "Fitting %d particles"%(fitnum)

        hr,hz = diskfitting.fit(the_particles['rxy'].in_units('kpc'),
                                np.abs(the_particles['pos'][:,2].in_units('kpc')))
        print "Powell hr: %g, hz: %g"%(hr,hz)
    if np.isfinite(hr):
        age = np.mean(thesestars['age'].in_units('Gyr'))
        ages.append({'mean':age,
                     'sigma':diskfitting.quantilesigma(thesestars['age'].in_units('Gyr')),
                     'stddev':np.std(thesestars['age'].in_units('Gyr'))})
        ofe = np.mean(thesestars['ofe'])
        ofes.append({'mean':ofe,'stddev':np.std(thesestars['ofe']),
                     'sigma':diskfitting.quantilesigma(thesestars['ofe'])})
        feh = np.mean(thesestars['feh'])
        fehs.append({'mean':feh,'stddev':np.std(thesestars['feh']),
                     'sigma':diskfitting.quantilesigma(thesestars['feh'])})

        mass.append(np.sum(thesestars['mass'].in_units('Msol')))
        
        hr, hz, hrerr, hzerr = diskfitting.mcerrors(
            the_particles['rxy'].in_units('kpc'),
            np.abs(the_particles['z'].in_units('kpc')),
            [hr,hz], plot=tfile+'.paramwalk%02d.png'%i,
            title="OFe: %2.2f, FeH: %2.2f, age: %2.1f Gyr, N:%d"%(ofe,feh,age,fitnum))

        print "emcee hr: %g, hz: %g"%(hr,hz)
        hzs.append(hz)
        rexps.append(hr)

        print "hrerr: %g, hzerr: %g"%(hrerr,hzerr)
        hzerrs.append(hzerr)
        hrerrs.append(hrerr)

        vertstarprof = profile.VerticalProfile(thesestars, rmin, rmax, zmax,ndim=2,
                                               nbins=20)
        vprofzbins.append(vertstarprof['rbins'].in_units('kpc'))
        vprofmasses.append(vertstarprof['mass'].in_units('Msol'))
        vprofn.append(vertstarprof['n'])

        radstarprof = profile.Profile(thesestars,min=rmin,max=rmax,nbins=20)
        rprofrbins.append(radstarprof['rbins'].in_units('kpc'))
        rprofmasses.append(radstarprof['mass'].in_units('Msol'))
        rprofdr.append(radstarprof['dr'].in_units('kpc'))
        rprofbs.append(radstarprof._binsize.in_units('kpc^2'))
        rprofn.append(radstarprof['n'])

        sigmavzs.append(diskfitting.velocity_dispersion(thesestars))
        

        diskfitting.plot_two_profiles(vertstarprof,radstarprof,
                                      rmin=rmin,rmax=rmax,zmin=zmin,zmax=zmax,
                                      hz=hz,hr=hr,hzerr=hzerr,hrerr=hrerr,
                                      units='Msol pc^-2',
                                      outfig=tfile+'.profs%02d'%i,
          title="OFe: %2.2f, FeH: %2.2f, age: %2.1f Gyr, N:%d"%(ofe,feh,age,fitnum))
        i+=1

    n_done+=n_block
    
    
pickle.dump({'hz':np.array(hzs), 'rexp': np.array(rexps), 'mass':np.array(mass),
             'hzerr':np.array(hzerrs), 'rexperr': np.array(hrerrs),
             'sigmavz': np.array(sigmavzs), 'age':ages, 'ofe':ofes,'feh':fehs,
             'rprofrbins':rprofrbins,'vprofzbins':vprofzbins,
             'rprofmasses':rprofmasses,'vprofmasses':vprofmasses,
             'vprofn':vprofn,'rprofn':rprofn,'rprofdr':rprofdr,
             'rprofbs':rprofbs},
            open(tfile+'.z0agedecomp.dat','w'))
