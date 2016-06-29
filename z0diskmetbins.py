import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pynbody, sys, math
from pynbody.analysis import profile, angmom, halo
from pynbody import units, config
import pynbody.filt as f
import diskfitting, emcee

from pynbody.analysis import luminosity as lum
import os, glob, pickle


tfile = '01024/g1536.01024'
s = pynbody.load(tfile)
h = s.halos()
pynbody.analysis.angmom.faceon(h[1])

sform=pynbody.load('starlogform/new.starlog.tipsy')
s.s['rxyform'] = sform.s['rxy']

hzs = []
rexps = []
hzerrs = []
hrerrs = []
age = []
ofes = []
fehs = []
sigmavzs = []
mass = []
rforms = []
vprofzbins = []
vprofmasses = []
vprofn = []
rprofrbins = []
rprofmasses = []
rprofn = []
rprofdr = []
rprofbs = []

minofe = 0 #np.min(diskstars['ofe'])
maxofe = 0.35 #np.max([diskstars['ofe'],0.5])
ofestep = 0.05
minfeh = -3 #np.min([diskstars['feh'],-3])
maxfeh = 0.5 #np.max([diskstars['feh'],0.5])
fehstep = 0.1

i=0
rmin = '6 kpc'
rmax = '10 kpc'
zmin = '0 kpc'
zmax = '3 kpc'

disk = pynbody.filt.SolarNeighborhood(r1=units.Unit(rmin),height=units.Unit(zmax))
diskstars = h[1].star[disk]

for ofe in np.arange(minofe,maxofe,ofestep):
    for feh in np.arange(minfeh,maxfeh,fehstep):
        print "ofe: %g, feh: %g"%(ofe,feh)
        
        metfilt = (f.BandPass('ofe',ofe,ofe+ofestep) & 
                   f.BandPass('feh',feh,feh+fehstep))
        thesestars = diskstars[metfilt]

        annulus = f.Disc(rmax,zmax) & ~f.Disc(rmin,zmax) & ~f.Disc(rmax,zmin)
        
        the_particles = thesestars[annulus]
        fitnum = len(the_particles)
            
        if fitnum > 100:
            
            print "Fitting %d particles"%(fitnum)
            
            hr,hz = diskfitting.fit(the_particles['rxy'].in_units('kpc'),
                                    np.abs(the_particles['z'].in_units('kpc')))
            print "Powell hr: %g, hz: %g"%(hr,hz)

            #if np.isfinite(hr):
            print "Using emcee to find errors"
            hr, hz, hrerr, hzerr = diskfitting.mcerrors(
                the_particles['rxy'].in_units('kpc'),
                np.abs(the_particles['z'].in_units('kpc')),
                [hr,hz], plot=tfile+'.mbparamwalk%02d.png'%i,
                title="OFe: %2.2f, FeH: %2.2f, N:%d"%(ofe,feh,fitnum))
            
            #if hrerr > 2:
            #    if isinstance(rmin, str): rmini=units.Unit(rmin)
            #    if isinstance(rmax, str): rmaxi=units.Unit(rmax)
            #    if isinstance(zmin, str): zmini=units.Unit(zmin)
            #    if isinstance(zmax, str): zmaxi=units.Unit(zmax)
                
            #    annulus = f.Disc(rmax,zmax) & ~f.Disc(rmin,zmax) & ~f.Disc(rmax,zmin)
            #    nwalkers, ndim = 6, 2
            #    sampler = emcee.EnsembleSampler(nwalkers, ndim, 
            #                                    diskfitting.twoexp_likelihood,
            #                 args=(thesestars.s[annulus]['rxy'].in_units('kpc'),
            #                       thesestars.s[annulus]['z'].in_units('kpc'),
            #                       rmini,rmaxi,zmini,zmaxi))
            #    p0 = [(np.random.rand(ndim)-[0.5,0.5])*0.05*[hr,hz]+[hr,hz] for idum in xrange(nwalkers)]
            #    sampler.run_mcmc(p0,1000)
            #    hrfc = sampler.flatchain[:,0]
            #    hr = sorted(hrfc)[int(np.floor(0.01*len(hrfc)))]
            
            print "hr: %g, hz: %g"%(hr,hz)
            hzs.append(hz)
            rexps.append(hr)
            
            print "hrerr: %g, hzerr: %g"%(hrerr,hzerr)
            hzerrs.append(hzerr)
            hrerrs.append(hrerr)
            
            vertstarprof = profile.VerticalProfile(thesestars, rmin, rmax, zmax,
                                                   ndim=2, nbins=20)
            vprofzbins.append(vertstarprof['rbins'].in_units('kpc'))
            vprofmasses.append(vertstarprof['mass'].in_units('Msol'))
            vprofn.append(vertstarprof['n'])
            
            radstarprof = profile.Profile(thesestars,min=rmin,max=rmax,nbins=20)
            rprofrbins.append(radstarprof['rbins'].in_units('kpc'))
            rprofmasses.append(radstarprof['mass'].in_units('Msol'))
            rprofdr.append(radstarprof['dr'].in_units('kpc'))
            rprofbs.append(radstarprof._binsize.in_units('kpc^2'))
            rprofn.append(radstarprof['n'])
            
            diskfitting.plot_two_profiles(vertstarprof,radstarprof,
                                          rmin=rmin,rmax=rmax,zmin=zmin,
                                          zmax=zmax,
                                          hz=hz,hr=hr,hzerr=hzerr,hrerr=hrerr,
                                          outfig=tfile+'.metbinprofs%02d'%i,
                                          units='Msol pc^-2',
                                          title="OFe: %2.2f, FeH: %2.2f, N:%d"%(ofe,feh,fitnum))
            
            sigmavzs.append(diskfitting.velocity_dispersion(thesestars))
            
            age.append({'mean':np.mean(thesestars['age'].in_units('Gyr')),
                        'sigma':diskfitting.quantilesigma(thesestars['age'].in_units('Gyr')),
                        'stddev':np.std(thesestars['age'].in_units('Gyr'))})
            rforms.append({'mean':np.mean(thesestars['rxyform'].in_units('kpc')),
                           'sigma':diskfitting.quantilesigma(thesestars['rxyform'].in_units('kpc')),
                           'stddev':np.std(thesestars['rxyform'].in_units('kpc'))})
            ofes.append(ofe)
            fehs.append(feh)
            
            mass.append(np.sum(thesestars['mass'].in_units('Msol')))
            i = i+1
            
        
pickle.dump({'hz':np.array(hzs), 'rexp': np.array(rexps), 'mass':np.array(mass),
             'hzerr':np.array(hzerrs), 'rexperr': np.array(hrerrs),
             'sigmavz': np.array(sigmavzs), 'age':age, 'ofe':ofes,'feh':fehs,
             'ofestep':ofestep,'fehstep':fehstep,'rform':rforms,
             'rprofrbins':rprofrbins,'vprofzbins':vprofzbins,
             'rprofmasses':rprofmasses,'vprofmasses':vprofmasses,
             'vprofn':vprofn,'rprofn':rprofn,'rprofdr':rprofdr,
             'rprofbs':rprofbs},
            open(tfile+'.z0metbinagedecomp.dat','w'))
