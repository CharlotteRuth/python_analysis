import pickle, matplotlib, matplotlib.pyplot as plt
import numpy as np
import monoAbundanceMW as mam

colormap = matplotlib.cm.jet

def _squeeze(o,omin,omax):
    return (o-omin)/(omax-omin)

d = pickle.load(open('01024/g1536.01024.z0agedecomp.dat'))
hz = d['hz']*1000
hzerr = d['hzerr']*1000
rexp = d['rexp'][(d['hzerr'] < 99)]
rexperr = d['rexperr'][(d['hzerr'] < 99)]
mass = d['mass'][(d['hzerr'] < 99)]

#Bovy stuff
bovyhz,bovyhzerr=zip(*[mam.hz(u[0],u[1],err=True) for u in zip(mam.fehs(),mam.afes())])
bovyrexp,bovyrexperr=zip(*[mam.hr(u[0],u[1],err=True) for u in zip(mam.fehs(),mam.afes())])
bovymass=np.array([mam.abundanceDist(u[0],u[1]) for u in zip(mam.fehs(),mam.afes())])

maxhr = 5
maxhz = 2000
bovyplotrexp = np.copy(bovyrexp)
bighr = bovyrexp > maxhr
bovyplotrexp[bighr] = maxhr-0.3

ax = plt.subplot(1,1,1)

ages = [i['mean'] for i in d['age']]
sizes = 10*(d['mass']/4e8)
print sizes
scat = plt.scatter(d['rexp'],hz,c=ages,s=sizes,edgecolors='none')
plt.xlim(1,6.5)
plt.ylim(0,2000)

ax.set_xlabel('$r_{exp}$ [kpc]')
ax.set_ylabel('$h_z$ [pc]')
cbar = plt.colorbar(scat)
cbar.set_label('Age [Gyr]')

for i in range(len(hz)):
    if (d['rexperr'][i] < 1):
        plt.errorbar(d['rexp'][i],hz[i],
                 xerr=d['rexperr'][i], yerr=hzerr[i],
                 color=colormap(_squeeze(ages[i],np.amin(ages),np.amax(ages))),
                 elinewidth=1.0,capsize=3,zorder=0)

plt.savefig('hzrexp.png')

plt.clf()

ax = plt.subplot(1,1,1)

ofes = [i['mean'] for i in d['ofe']]
scat = plt.scatter(d['rexp'],hz,c=ofes)
plt.xlim(0,5)
plt.ylim(0,2000)

ax.set_xlabel('$r_{exp}$ [kpc]')
ax.set_ylabel('$h_z$ [pc]')
cbar = plt.colorbar(scat)
cbar.set_label('[O/Fe]')

plt.savefig('hzrexpofe.png')

plt.clf()

ax = plt.subplot(1,1,1)

ages = [i['mean'] for i in d['age']]
fehs = [i['mean'] for i in d['feh']]
ofes = [i['mean'] for i in d['ofe']]
scat = plt.scatter(fehs,ofes,c=ages)

ax.set_xlabel('[Fe/H]')
ax.set_ylabel('[O/Fe]')
cbar = plt.colorbar(scat)
cbar.set_label('Age [Gyr]')

plt.savefig('ofefeh.png')

plt.clf()

ax = plt.subplot(1,1,1)

scat = plt.scatter(fehs,ofes,c=hz)

ax.set_xlabel('[Fe/H]')
ax.set_ylabel('[O/Fe]')
cbar = plt.colorbar(scat)
cbar.set_label('$h_z$ [kpc]')

plt.savefig('ofefehhz.png')

plt.clf()

ax = plt.subplot(1,1,1)

scat = plt.scatter(fehs,ofes,c=d['rexp'])

ax.set_xlabel('[Fe/H]')
ax.set_ylabel('[O/Fe]')
cbar = plt.colorbar(scat)
cbar.set_label('$r_{exp}$ [kpc]')

plt.savefig('ofefehrexp.png')

plt.clf()

import math
from matplotlib.ticker import FormatStrFormatter
rmin=6
rmax=10
zmin=0 
zmax=3
sigR0s = []
for i,hr in enumerate(rexp):
    hzexp = 2*hz[i] / 1000
    norm_int = (4.0*math.pi*hr*hzexp*(-np.exp(-rmax/hr)*(hr + rmax) + 
                                    np.exp(-rmin/hr)*(hr + rmin))*
                (math.tanh(zmax/hzexp)-math.tanh(zmin/hzexp)))
    rho0 = mass[i] / norm_int
    sigmaR0 = rho0*np.exp(-8 / hr)*2*(math.tanh(zmax/hzexp)-math.tanh(zmin/hzexp)) / 1e6
    sigR0s.append(sigmaR0)
    if mass[i] > 3*np.min(mass):
        #import pdb; pdb.set_trace()
        plt.errorbar(hz[i],sigmaR0,
                     xerr=hzerr[i],
                     color=colormap(_squeeze(ages[i],np.amin(ages),np.amax(ages))),
                     elinewidth=1.0,capsize=3,zorder=0)

plt.semilogy()
scat = plt.scatter(hz,np.array(sigR0s),c=ages,s=sizes,edgecolors='none')
cbar = plt.colorbar(scat)
cbar.set_label('Age [Gyr]')

plt.hist(hz,weights=np.array(sigR0s),bins=15,histtype='step',color='k',
         linewidth=2,range=[165,1800],label='MaGICC simulation')
plt.hist(bovyhz,range=[165,1200],weights=bovymass*0.13/0.07,color='#444444',
         histtype='step',bins=10,lw=1.,zorder=10,label='Bovy+ (2012b)')
plt.xlim(0,maxhz)
plt.ylim(0.4,17)
ax = plt.gca()
ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))
plt.xlabel('$h_z$ [pc]')
plt.ylabel('$\Sigma_{R_0}(h_z)$ [M$_\odot$ pc$^{-2}$]')
plt.legend(loc=0)

plt.savefig('agesurfdenshz.eps')
