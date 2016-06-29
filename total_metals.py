import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pynbody, sys
from pynbody.analysis import profile, angmom, halo
from pynbody import units, config
import pynbody.filt as f

from pynbody.analysis import luminosity as lum
import os, glob

XSOLfe=0.117E-2
XSOLO=0.96E-2
XSOLH=0.706
ZSUN = 0.02
grav = 6.6725985e-11 # m^3 kg^-1 s^-2
s_per_yr = 3.15569e7

#filename = '/home/christensen/Storage1/UW/MolecH/Cosmo/h937.cosmo25cmb.4096g/h937.cosmo25cmb.4096g1MbwK1C52/h937.cosmo25cmb.4096g1MbwK1C52.004096/h937.cosmo25cmb.4096g1MbwK1C52.004096'
#haloind = 7

#filename = '/home/christensen/Storage2/UW/MolecH/Cosmo/h239.cosmo50cmb.3072g/h239.cosmo50cmb.3072g14HMbwK/steps/h239.cosmo50cmb.3072g14HMbwK.00512.dir/h239.cosmo50cmb.3072g14HMbwK.00512'
#haloind = 11

#Use amiga, not AHF here
filename = '/home/christensen/Storage2/UW/MolecH/Cosmo/h258.cosmo50cmb.3072g/h258.cosmo50cmb.3072g14HMbwK/steps/h258.cosmo50cmb.3072g14HMbwK.00512.dir/h258.cosmo50cmb.3072g14HMbwK.00512'
haloind = 15

#Use amiga, not AHF here
filename = '/home/christensen/Storage2/UW/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/h986.cosmo50cmb.3072g14HBWK_2.00512.dir/h986.cosmo50cmb.3072g14HBWK.00512'
haloind = 22

tempcut = 1e4
denscut = 100 #amu/cc
cstar = 0.1
deltat = 1e6*s_per_yr #yr
s = pynbody.load(filename)
h = s.halos()
#halo = h.load_copy(haloind) 
#OR
halo = h[haloind]

print('Stellar Mass: %e' % float(sum(halo.stars['mass'].in_units('Msol'))))
lum = pynbody.analysis.luminosity.halo_lum(halo,band='v')
print('V-Band Lum: %e' % float(lum))
print('V-Band Mag: %e' % float(pynbody.analysis.luminosity.halo_mag(halo,band='v')))
print('B-Band Mag: %e' % float(pynbody.analysis.luminosity.halo_mag(halo,band='b')))
print('HI Mass: %e' % float(sum(halo.gas['mass'].in_units('Msol')*halo.gas['HI'])))
print('f_gas: %f' % float(sum(halo.gas['mass'].in_units('Msol')*halo.gas['HI'])/(sum(halo.gas['mass'].in_units('Msol')*halo.gas['HI']) + sum(halo.stars['mass'].in_units('Msol')))))
print('Ox Mass: %e' % float(sum(halo.stars['mass'].in_units('Msol')*halo.stars['OxMassFrac'])))

#Align halo
pynbody.analysis.angmom.faceon(halo) 

#Image of the galaxy
halo.gas['hiden'] = halo.gas['rho']*halo.gas['HI']
plt.figure(0)
pynbody.plot.image(halo.gas,qty='hiden',units='m_p cm^-2')
plt.close

#Star formation history
sfh = pynbody.plot.stars.sfh(halo.s)
sfhcum = np.cumsum(sfh[0])
plt.figure(0)
plt.plot(pynbody.analysis.cosmology.age(s, unit='Gyr') - sfh[1][1:len(sfhcum)+1],sfhcum/max(sfhcum))
plt.xlabel('Age [Gyr]'); plt.ylabel('Cumulative SFH')
plt.axis([14, 0, 0, 1])
plt.close

#Calculating the metallicity
halo.gas['SFeff'] = 1.0 - exp(-1.0*cstar*sqrt(4.0*pi*grav*units.Unit("m**3 kg**-1 s**-2")*halo.gas['rho'].in_units('kg m**-3'))*deltat*units.Unit('s'))
lowdense = pynbody.filt.LowPass('rho','100 m_p cm**-3 a**-3')
halo.gas[lowdense]['SFeff'] = 0
hot = pynbody.filt.HighPass('temp','10000 K')
halo.gas[hot]['SFeff'] = 0

allhe_mass = halo.gas['mass'].in_units('Msol')*halo.gas['HeI'] + halo.gas['mass'].in_units('Msol')*halo.gas['HeII']
allhe_plus_h_mass = (1-halo.gas['metals'])*halo.gas['mass'].in_units('Msol')
allh_mass = allhe_plus_h_mass - (4*allhe_mass)
allnh = (allh_mass/1.6726e-24)*(2.0e33)
alln_ox = (halo.gas['OxMassFrac']*halo.gas['mass'].in_units('Msol'))*2.0e33/2.66e-23
all_ox_abund = alln_ox/allnh
print("12 + logOx (total): %f" % float(12+log10(mean(all_ox_abund))))

#halo.gas['rho'].in_units('m_p cm**-3 a**-3')
dense = pynbody.filt.HighPass('rho','0.1 m_p cm**-3 a**-3')
densegas = halo.gas[dense]
cold = pynbody.filt.LowPass('temp','12000 K')
colddensegas = densegas[cold]
allhe_mass = colddensegas['mass'].in_units('Msol')*colddensegas['HeI'] + colddensegas['mass'].in_units('Msol')*colddensegas['HeII']
allhe_plus_h_mass = (1-colddensegas['metals'])*colddensegas['mass'].in_units('Msol')
allh_mass = allhe_plus_h_mass - (4*allhe_mass)
allnh = (allh_mass/1.6726e-24)*(2.0e33)
alln_ox = (colddensegas['OxMassFrac']*colddensegas['mass'].in_units('Msol'))*2.0e33/2.66e-23
all_ox_abund = alln_ox/allnh
print("12 + logOx (cold gas): %f" % float(12+log10(mean(all_ox_abund))))

p_age = profile.Profile(halo.s, calc_x = lambda x: x.s['age'].in_units('Gyr'), max = '14 Gyr')
plt.plot(p_age['rbins'], log10(p_age['metals']/ZSUN), 'k', label = 'mean z',hold=False)
plt.xlabel('Age [Gyr]'); plt.ylabel('Log(Z/Z_sun)')
plt.axis([14, 0, -3, -1])
