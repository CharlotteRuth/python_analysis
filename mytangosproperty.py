# export TANGOS_PROPERTY_MODULES=mytangosproperty

from tangos.properties import PropertyCalculation
from tangos.properties.pynbody import PynbodyPropertyCalculation
import numpy as np
import pynbody

"""
class ExampleHaloProperty(PropertyCalculation):
    names = "myproperty"

    def calculate(self, particle_data, existing_properties):
        return 42.0
"""

"""
class ExampleHaloProperty(PropertyCalculation):
    names = 'my_x'

    def calculate(self, particle_data, existing_properties):
        return existing_properties['Xc']
"""
    

class ExampleHaloProperty(PynbodyPropertyCalculation):
    names = "meanFeH", "meanOxH", "meanOxFe"
    
    def calculate(self, particle_data, existing_properties):
        XSOLFe = 0.125E-2         # 1.31e-3
        # Looks very wrong ([O/Fe] ~ 0.2-0.3 higher than solar),
        # probably because SN ejecta are calculated with
        # Woosley + Weaver (1995) based on Anders + Grevesse (1989)
        # XSOLO=0.59E-2           # 5.8e-2
        XSOLO = 0.84E-2
        XSOLH = 0.706             # 0.74

        stars = particle_data.stars
        if len(stars) == 0:
            return -10, -10, -10
        stars[pynbody.filt.LowPass('OxMassFrac',1e-7)].star['OxMassFrac'] = 1e-7
        stars[pynbody.filt.LowPass('FeMassFrac',1e-8)].star['FeMassFrac'] = 1e-8

        #stars['metals'] =2.09*stars['OxMassFrac'] + 1.06*stars['FeMassFrac'] 
        #stars['hydrogen'] = 0.764 - (3.1 * self['metals'])
        hydrogen = stars.star['hydrogen']
        hydrogen = np.array(hydrogen)
        hydrogen[hydrogen < 1e-2] = 0.01
        #stars[pynbody.filt.LowPass('hydrogen',1e-2)].star['hydrogen'] = 0.01 # Not sure why some super-metallicity stars

        
        stars['feh'] = np.log10(stars['FeMassFrac'] / hydrogen) - np.log10(XSOLFe / XSOLH) 
        stars['oxh'] = np.log10(stars['OxMassFrac'] / hydrogen) - np.log10(XSOLO / XSOLH)
        # Following Applebaum+ 2021
        stars[pynbody.filt.LowPass('oxh',-4)]['oxh'] = -4
        stars[pynbody.filt.LowPass('feh',-4)]['feh'] = -4
	#stars[pynbody.filt.LowPass('feh',-4)]['feh']= -4
        stars['ofe'] = np.log10(stars['OxMassFrac'] / stars['FeMassFrac']) - np.log10(XSOLO / XSOLFe)

        # Calculating the mean of the log (folloing Kirby+ 2013, figure 1
        meanFeH = np.sum(stars['mass']*stars['feh'])/np.sum(stars['mass'])
        meanOxH = np.sum(stars['mass']*stars['oxh'])/np.sum(stars['mass'])
        meanOxFe = np.sum(stars['mass']*stars['ofe'])/np.sum(stars['mass'])
        
        return meanFeH, meanOxH, meanOxFe
    
    # names = "FeH", "OxH", "FeOx"
#    names = "velocity_dispersion"

#    def calculate(self, particle_data, existing_properties):
#        return np.std(particle_data['vel'])

    
