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


class ExampleHaloProperty(PynbodyPropertyCalculation):
    names = 'my_x'

    def calculate(self, particle_data, existing_properties):
        return existing_properties['Xc']

    
"""
class ExampleHaloProperty(PynbodyPropertyCalculation):
    names = "velocity_dispersion"

    def calculate(self, particle_data, existing_properties):
        return np.std(particle_data['vel'])
    
"""
