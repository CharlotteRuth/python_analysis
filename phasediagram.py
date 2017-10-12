#Charlotte Christensen
#5/20/17
# This code will create a phase diagram of a simulation that includes multiphase particles from superbubble feedback.
# I will be using it to determine the mass fraction in the hot phase for gas particles that are very dense but unable to form stars.

#run from ipython command line:
#run /home/christensen/Code/python/phasediagram.py

import pynbody
import pynbody.plot as pp
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

#Load the simulation and adjust units
#sim = '/home/christensen/Storage1/UW/MolecH/12M/Disk_Collapse_1e6_sb_H2/old/Disk_Collapse_1e6.00100.scalez1.SB.metal.NR.000001'
#sim = '/home/christensen/Storage1/UW/MolecH/12M/Disk_Collapse_1e6_sb_H2/Disk_Collapse_1e6_sb_H2.000001'
#sim = '/home/christensen/Storage1/UW/MolecH/12M/Disk_Collapse_1e5/scalez1.SB.H2.NR.noJeans/Disk_Collapse_1e5.00100.scalez1.SB.H2.NR.000005'
sim = '/home/christensen/Storage1/UW/MolecH/12M/Disk_Collapse_1e5/scalez1.SB.H2.NR.noJeans/Disk_Collapse_1e5.00100.scalez1.SB.H2.NR.000050'
sim = '/home/christensen/Storage1/UW/MolecH/12M/Disk_Collapse_1e5/scalez1.SB.H2.NR/Disk_Collapse_1e5.00100.scalez1.SB.H2.NR.000047'
sim = '/home/christensen/Storage1/UW/MolecH/12M/Disk_Collapse_1e5/scalez1.SB.H2.NR.LowT/Disk_Collapse_1e5.00100.scalez1.SB.H2.NR.000049'
sim = '/home/christensen/Storage1/UW/MolecH/12M/Disk_Collapse_1e5/scalez1.SB.H2.NR.noJeans/Disk_Collapse_1e5.00100.scalez1.SB.H2.NR.NJ.000047'
sim = '/home/christensen/Storage1/UW/MolecH/12M/Disk_Collapse_1e5/scalez1.SB.H2.NR.noJeans.LowT/Disk_Collapse_1e5.00100.scalez1.SB.H2.NR.NJ.lowT.000052'
sim = '/home/christensen/Storage1/UW/MolecH/12M/Disk_Collapse_1e6.00100.scalez1.SB.H2.NR.lowT/Disk_Collapse_1e6.00100.scalez1.SB.H2.NR.000005'
sim = '/home/christensen/Storage1/UW/MolecH/12M/Disk_Collapse_1e6.00100.scalez1.SB.H2.NR.noJeans.lowT/Disk_Collapse_1e6.00100.scalez1.SB.H2.NR.000005'
sim = '/home/christensen/Storage1/UW/MolecH/12M/Disk_Collapse_1e5/Disk_Collapse_1e5.00100.scalez1.SB.H2.codetest/Disk_Collapse_1e5.00100.scalez1.SB.H2.000005'
sim = '/home/christensen/Storage1/UW/MolecH/12M/Disk_Collapse_1e5/scalez1.SB.H2.NR/Disk_Collapse_1e5.00100.scalez1.SB.H2.NR.000006'
sim = '/home/christensen/Storage1/UW/MolecH/12M/Disk_Collapse_1e5/Disk_Collapse_1e5.00100.scalez1.SB.H2.NR.bk/Disk_Collapse_1e5.00100.scalez1.SB.H2.NR.000005'
sim = '/home/christensen/Storage1/UW/MolecH/12M/Disk_Collapse_1e5/Disk_Collapse_1e5.00100.scalez1.SB.H2.NR.codetest/Disk_Collapse_1e5.00100.scalez1.SB.H2.NR.000018'
sim = '/home/christensen/Storage1/UW/MolecH/12M/Disk_Collapse_1e5/Disk_Collapse_1e5.00100.scalez1.SB.metal.NR/Disk_Collapse_1e5.00100.scalez1.SB.metalNR.000002'
sim = '/home/christensen/Storage1/UW/MolecH/12M/Disk_Collapse_1e5/Disk_Collapse_1e5.00100.scalez1.SB.metal.NR/Disk_Collapse_1e5.00100.scalez1.SB.metalNR.000027'

s = pynbody.load(sim)
print(max(s.g['H2']))
print(min(s.g['temp']))
index = np.where(s.g['massHot'] > 0)
multiPhaseGas = s.g[index]

#gammam1 = 1.66 # For a monotomic gas
#PoverRho = gammam1*(p->uHotPred()*frac+p->uPred()*(1-frac));
#fDensity = p->fDensity*PoverRho/(gammam1*p->uPred());

#Create a phase diagram
nmin, nmax, Tmin, Tmax =-9,9,0.5,8
plt.figure(1)
pp.rho_T(s,rho_range=[nmin,nmax],t_range = [Tmin,Tmax],rho_units='m_p cm^-3')
print(max(s.gas['rho'].in_units('m_p cm^-3')))
plt.savefig(sim+'.phase.png')

#Create a histogram of gas densities
x = np.log10(s.gas['rho'].in_units('m_p cm^-3'))
mass = s.gas['mass'].in_units('Msol')
plt.figure(2)
hist(x,weights=mass,range=[nmin,nmax],bins=100,histtype='step',log=True)
plt.xlabel(r'log$_{10}(\rho$/m$_p$cm$^{-3}$)') #cm$^{-3}$
plt.ylabel('Mass of gas')
plt.savefig(sim+'.rhohist.png')

#Select that dense gas and determine what the hot gas fraction is
densecut_units = '3000000 m_p cm^-3'
densecut = np.log10(3000000)

densegas = s.g[pynbody.filt.HighPass('rho',densecut_units)]
#densegas = s.g[pynbody.filt.HighPass('rho','1000000 m_p cm^-3')]
plt.figure(3)
#hist(np.log10(densegas['rho'].in_units('m_p cm^-3')),weights=densegas['mass'].in_units('Msol'),bins=100,histtype='step',log=True)
hist(np.log10(densegas['rho'].in_units('m_p cm^-3')),bins=100,histtype='step',range=[densecut,nmax])
plt.xlabel(r'log$_{10}(\rho$/m$_p$cm$^{-3}$)') #cm$^{-3}$
plt.ylabel('Number of Particles') 
plt.savefig(sim+'.densrhohist.png')

plt.figure(4)
hist(densegas['massHot'],range=[0,1],bins=100,histtype='step',log=True)
#hist(s.g['massHot'],range=[0,1],bins=100,histtype='step',log=True)
plt.xlabel('Mass fraction in hot phase') #cm$^{-3}$
plt.ylabel('Number of Particles') 
plt.savefig(sim+'.hotphasehist.png')
print(densegas['massHot'])

plt.figure(5)
plot(np.log10(s.g['rho'].in_units('m_p cm^-3')),s.g['massHot'],'.')
ylim(-.01,.3)
xlim(6,7)
plt.xlabel(r'log$_{10}(\rho$/m$_p$cm$^{-3}$)') #cm$^{-3}$
plt.ylabel('Mass fraction in hot phase') 
plt.savefig(sim+'.hotphase_v_rho.png')

plt.figure(6)
plot(np.log10(s.g['rho'].in_units('m_p cm^-3')),np.log10(s.g['massHot']),'.')
plt.xlabel(r'log$_{10}(\rho$/m$_p$cm$^{-3}$)') #cm$^{-3}$
plt.ylabel(r'log$_{10}$(Mass fraction in hot phase)') 
plot([6.5,6.5],[-4,0])
plt.savefig(sim+'.hotphase_v_rho_zoom.png')
print(densegas['H2'])

plt.figure(7)
plot(np.log10(s.g['rho'].in_units('m_p cm^-3')),np.log10(2*s.g['H2']/(2*s.g['H2'] + 2*s.g['HI'])),'.')
plot(np.log10(multiPhaseGas['rho'].in_units('m_p cm^-3')),np.log10(2*multiPhaseGas['H2']/(2*multiPhaseGas['H2'] + 2*multiPhaseGas['HI'])),'.')
plt.xlabel(r'log$_{10}(\rho$/m$_p$cm$^{-3}$)') #cm$^{-3}$
plt.ylabel(r'log$_{10}$H2/(HI + H2)') 
plt.savefig(sim+'.H2fractionlog.png')

plt.figure(8)
plot(np.log10(s.g['rho'].in_units('m_p cm^-3')),2*s.g['H2']/(2*s.g['H2'] + 2*s.g['HI']),'.')
plot(np.log10(multiPhaseGas['rho'].in_units('m_p cm^-3')),2*multiPhaseGas['H2']/(2*multiPhaseGas['H2'] + 2*multiPhaseGas['HI']),'.')
plt.xlabel(r'log$_{10}(\rho$/m$_p$cm$^{-3}$)') #cm$^{-3}$
plt.ylabel('H2/(HI + H2)') 
plt.savefig(sim+'.H2fraction.png')
print(densegas['H2'])
