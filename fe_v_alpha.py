
# on emu
# cd /home/christenc/REPOSITORY/e12Gals/h148.cosmo50PLK.6144g3HbwK1BH/snapshots_200crit_h148mint


import pynbody
s = pynbody.load('h148.cosmo50PLK.6144g3HbwK1BH.004096')
h = s.halos()
halo = h.load_copy(2)

stars = halo.star
stars[pynbody.filt.LowPass('OxMassFrac',1e-7)].star['OxMassFrac'] = 1e-7
stars[pynbody.filt.LowPass('FeMassFrac',1e-8)].star['FeMassFrac'] = 1e-8

XSOLFe = 0.125E-2         # 1.31e-3
# Looks very wrong ([O/Fe] ~ 0.2-0.3 higher than solar),
# probably because SN ejecta are calculated with
# Woosley + Weaver (1995) based on Anders + Grevesse (1989)
# XSOLO=0.59E-2           # 5.8e-2
XSOLO = 0.84E-2
XSOLH = 0.706             # 0.74

# Redefined to use the set minimum values of FeMassFrac and OxMassFrac
stars['feh'] = np.log10(stars['FeMassFrac'] / stars['hydrogen']) - np.log10(XSOLFe / XSOLH) 
stars['oxh'] = np.log10(stars['OxMassFrac'] / stars['hydrogen']) - np.log10(XSOLO / XSOLH)
stars['ofe'] = np.log10(stars['OxMassFrac'] / stars['FeMassFrac']) - np.log10(XSOLO / XSOLFe)

# Calculating the mean of the log (folloing Kirby+ 2013, figure 1
meanOxH = np.sum(stars['mass']*stars['oxh'])/np.sum(stars['mass'])

meanFeH = np.sum(stars['mass']*stars['feh'])/np.sum(stars['mass'])

meanOxFe = np.sum(stars['mass']*stars['ofe'])/np.sum(stars['mass'])
