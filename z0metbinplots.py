import pickle, matplotlib, matplotlib.pyplot as plt
import numpy as np
import diskfitting
import monoAbundanceMW as mam

colormap = matplotlib.cm.jet

d = pickle.load(open('01024/g1536.01024.z0metbinagedecomp.dat'))
hz = d['hz'][(d['hzerr'] < 99)]*1000
hzerr = d['hzerr'][(d['hzerr'] < 99)]*1000
rexp = d['rexp'][(d['hzerr'] < 99)]
rexperr = d['rexperr'][(d['hzerr'] < 99)]
mass = d['mass'][(d['hzerr'] < 99)]
maxhr = 5
maxhz = 1990

sizes = 5*mass/np.min(mass)

ages = np.array([i['mean'] for i in d['age']])
ages = ages[(d['hzerr'] < 99)]

plotrexp = np.copy(rexp)
bighr = rexp > maxhr
plotrexp[bighr] = maxhr-0.3

#Bovy stuff
bovyhz,bovyhzerr=zip(*[mam.hz(u[0],u[1],err=True) for u in zip(mam.fehs(),mam.afes())])
bovyrexp,bovyrexperr=zip(*[mam.hr(u[0],u[1],err=True) for u in zip(mam.fehs(),mam.afes())])
bovymass=np.array([mam.abundanceDist(u[0],u[1]) for u in zip(mam.fehs(),mam.afes())])

bovyplotrexp = np.copy(bovyrexp)
bighr = bovyrexp > maxhr
bovyplotrexp[bighr] = maxhr-0.3

def _squeeze(o,omin,omax):
    return (o-omin)/(omax-omin)

def dotplot(plotrexp,hz,rexp=rexp,rexperr=rexperr,hzerr=hzerr, mass=mass,
            plotc=ages,sizes=sizes,alpha=1.0,axes=False,
            cbar_label=False,cmap=colormap):
    if not axes: axes = plt.gca()
    scat = axes.scatter(plotrexp,hz,s=sizes,c=plotc,cmap=cmap,alpha=alpha,
                       edgecolors='none')
    axes.set_xlim(0.5,maxhr)
    axes.set_ylim(0,maxhz)
    
    axes.set_xlabel('$r_{exp}$ [kpc]')
    axes.set_ylabel('$h_z$ [pc]')
    if cbar_label:
        cbar = plt.colorbar(scat,ax=axes)
        cbar.set_label(cbar_label)

    for i in range(len(hz)):
        if type(plotc) is str: col=plotc
        else: col=cmap(_squeeze(plotc[i],np.amin(plotc),np.amax(plotc)))
        if rexp[i] < maxhr:
            if mass[i] > np.min(mass):
                axes.errorbar(rexp[i],hz[i], xerr=rexperr[i], yerr=hzerr[i],
                             color=col, elinewidth=1.0,capsize=3,zorder=0)
        else:
            axes.errorbar(maxhr-0.2,hz[i],xerr=0.1, xuplims=True,
                         color=col,elinewidth=1.0,capsize=3)


f,ax=plt.subplots(2,1)
f.subplots_adjust(left=0.1,right=0.98,top=0.97,hspace=0.)
dotplot(bovyplotrexp,bovyhz,rexperr=bovyrexperr,hzerr=bovyhzerr,
        rexp=bovyrexp,mass=bovymass,axes=ax[0],
        plotc=mam.afes(),sizes=bovymass*200,cbar_label=r'[$\alpha$/Fe]')
ax[0].text(0.6,0.9,'Milky Way (Bovy+ 2012)',transform=ax[0].transAxes)
plt.setp(ax[0].get_xticklabels(),visible=False)
ofes = np.array(d['ofe'])[(d['hzerr'] < 99)]
dotplot(plotrexp,hz,plotc=ofes,axes=ax[1],cbar_label='[O/Fe]')
ax[1].text(0.6,0.9,'MaGICC Simulation',transform=ax[1].transAxes)

plt.savefig('mbhzrexpofe.eps')

plt.clf()

dotplot(plotrexp,hz,plotc=ages)

plt.savefig('mbhzrexpage.png')

plt.clf()

fehs = np.array(d['feh'])[(d['hzerr'] < 99)]
dotplot(plotrexp,hz,plotc=fehs,cbar_label='[Fe/H]')

plt.savefig('mbhzrexpfeh.png')

plt.clf()

diskfitting.block_histogram(qtytitle='r$_{exp}$ [kpc]')
plt.savefig('rexpblockofefeh.png')
plt.clf()

diskfitting.block_histogram(qty='hz',qtytitle='$h_z$ [kpc]',vmin=0.2,vmax=2)
plt.savefig('hzblockofefeh.png')
plt.clf()

diskfitting.block_histogram(qty='meanage',qtytitle='Age [Gyr]',vmin=0,vmax=14)
plt.savefig('ageblockofefeh.png')
plt.clf()

diskfitting.block_histogram(qty='agedisp',qtytitle='$1\sigma$ Age [Gyr] (blocks)',vmin=0,vmax=3)
d = pickle.load(open('01024/g1536.01024.z0agedecomp.dat'))
aages = [i['mean'] for i in d['age']]
afehs = [i['mean'] for i in d['feh']]
aofes = [i['mean'] for i in d['ofe']]
scat = plt.scatter(afehs,aofes,c=aages,edgecolors='none')
cbar = plt.colorbar(scat)
cbar.set_label('Age [Gyr] (points)')
plt.savefig('agedispblockofefeh.png')
plt.clf()

diskfitting.block_histogram(qty='meanrform',qtytitle='Rform [kpc]',vmin=0,vmax=12)
plt.savefig('rformblockofefeh.pdf')
plt.clf()

diskfitting.block_histogram(qty='rformdisp',qtytitle='$1\sigma$ Rform [kpc]',vmin=0,vmax=5)
plt.savefig('rformdispblockofefeh.pdf')
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
                     color=colormap(_squeeze(ofes[i],np.amin(ofes),np.amax(ofes))),
                     elinewidth=1.0,capsize=3,zorder=0)

plt.semilogy()
scat = plt.scatter(hz,np.array(sigR0s),c=ofes,s=sizes,edgecolors='none')
cbar = plt.colorbar(scat)
cbar.set_label('[O/Fe]')

plt.hist(hz,weights=np.array(sigR0s),bins=8,histtype='step',color='k',
         linewidth=2,range=[0,1800])
plt.hist(bovyhz,range=[165,1200],weights=bovymass*0.13/0.07,color='#444444',
         histtype='step',bins=12,lw=1.,zorder=10)
plt.xlim(0,maxhz)
plt.ylim(0.02,12)
ax = plt.gca()
ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))
plt.xlabel('$h_z$ [pc]')
plt.ylabel('$\Sigma_{R_0}(h_z)$ [M$_\odot$ pc$^{-2}$]')

plt.savefig('surfdenshz.png')
