
# coding: utf-8

# In[8]:

#Importing
import numpy as np
import pynbody
import pynbody.filt as filt
import pynbody.units as units
get_ipython().magic(u'matplotlib inline')
from pynbody.filt import *
from pynbody import tipsy
from array import array
import matplotlib.pylab as plt
import matplotlib.mlab as mlab
from pylab import *
import matplotlib.pyplot  as pyplot
from __future__ import division
from scipy.optimize import curve_fit
from numpy import sqrt, pi, exp, linspace, random
from scipy.integrate import trapz as trapz

radius = [0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1]

def make_intensity (redshift, *args):
    init_eq_width = [] # this will hold the inital eq width values which will then be averaged for each radius
    final_eq_width = [] # this will be the final array that holds the average eq widths at each radius 
    SD_array = [] # this will hold the standward deviations of each of the eq widths at a specific radius
    n = 1 # this is a counter so see which argument in *args we are in 
    for i in args:
        intensity = exp(-1 * i)
        area = trapz(i, dx=(abs(redshift[1] - redshift[0])))
        eq_width = area/ (max(i))
        init_eq_width.append(eq_width)

        if n == 1:
            avg_eq = mean(init_eq_width)
            sd = np.std(init_eq_width)
            final_eq_width.append(avg_eq)
            SD_array.append(sd)
            init_eq_width = []
        
        if n >= 2 and n <= 5:
            if n == 5:
                avg_eq = mean(init_eq_width)
                sd = np.std(init_eq_width)
                final_eq_width.append(avg_eq)
                SD_array.append(sd)
                init_eq_width = []
                
        if n >= 6 and n <= 9:
            if n == 9:
                avg_eq = mean(init_eq_width)
                sd = np.std(init_eq_width)
                final_eq_width.append(avg_eq)
                SD_array.append(sd)
                init_eq_width = []
        
        if n >= 10 and n <= 15:
            if n == 15:
                avg_eq = mean(init_eq_width)
                sd = np.std(init_eq_width)
                final_eq_width.append(avg_eq)
                SD_array.append(sd)
                init_eq_width = []
            
        if n >= 16 and n <= 24:
            if n == 24:
                avg_eq = mean(init_eq_width)
                sd = np.std(init_eq_width)
                final_eq_width.append(avg_eq)
                SD_array.append(sd)
                init_eq_width = []
            
        if n >= 25 and n <= 35:
            if n == 35:
                avg_eq = mean(init_eq_width)
                sd = np.std(init_eq_width)
                final_eq_width.append(avg_eq)
                SD_array.append(sd)
                init_eq_width = []
            
        if n >= 36 and n <= 47:
            if n == 47:
                avg_eq = mean(init_eq_width)
                sd = np.std(init_eq_width)
                final_eq_width.append(avg_eq)
                SD_array.append(sd)
                init_eq_width = []
            
        if n >= 48 and n <= 60:
            if n == 60:
                avg_eq = mean(init_eq_width)
                sd = np.std(init_eq_width)
                final_eq_width.append(avg_eq)
                SD_array.append(sd)
                init_eq_width = []
            
        if n >= 61 and n <= 73:
            if n == 73:
                avg_eq = mean(init_eq_width)
                sd = np.std(init_eq_width)
                final_eq_width.append(avg_eq)
                SD_array.append(sd)
                init_eq_width = []
            
        if n >= 74 and n <= 91:
            if n == 91:
                avg_eq = mean(init_eq_width)
                sd = np.std(init_eq_width)
                final_eq_width.append(avg_eq)
                SD_array.append(sd)
                init_eq_width = []
            
        if n >= 92 and n <= 111:
            if n == 111:
                avg_eq = mean(init_eq_width)
                sd = np.std(init_eq_width)
                final_eq_width.append(avg_eq)
                SD_array.append(sd)
                init_eq_width = []
        n+=1       
    print ("                                   Average Equivalent Widths")
    print(final_eq_width)
    print ("                                   Standard Deviation")
    print(SD_array)
    
        
def equivalent_width(x, width, stand_dev):
    plot(x, width)
    errorbar(x, width, yerr = stand_dev, fmt = 'bo')
    title("Equivalent Width")
    xlabel("Radius (in fraction of Radius)")
    ylabel("Equivalent Width")


# In[9]:

r0 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS1/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.0.0")
redshift = r0[:,0]
HI_r0 = r0[:,7]
HeI_r0 = r0[:,11]
CIV_r0 = r0[:,27]
OI_r0 = r0[:,47]
OVI_r0 = r0[:,67]
MgII_r0 = r0[:,87]
SiII_r0 = r0[:,103]
FeII_r0 = r0[:,127]

r1_1 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS1/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.1.1")
redshift = r1_1[:,0]
HI_r1_1 = r1_1[:,7]
HeI_r1_1 = r1_1[:,11]
CIV_r1_1 = r1_1[:,27]
OI_r1_1 = r1_1[:,47]
OVI_r1_1 = r1_1[:,67]
MgII_r1_1 = r1_1[:,87]
SiII_r1_1 = r1_1[:,103]
FeII_r1_1 = r1_1[:,127]

r1_2 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS1/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.2.1")
redshift = r1_2[:,0]
HI_r1_2 = r1_2[:,7]
HeI_r1_2 = r1_2[:,11]
CIV_r1_2 = r1_2[:,27]
OI_r1_2 = r1_2[:,47]
OVI_r1_2 = r1_2[:,67]
MgII_r1_2 = r1_2[:,87]
SiII_r1_2 = r1_2[:,103]
FeII_r1_2 = r1_2[:,127]

r1_3 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS1/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.3.1")
redshift = r1_3[:,0]
HI_r1_3 = r1_3[:,7]
HeI_r1_3 = r1_3[:,11]
CIV_r1_3 = r1_3[:,27]
OI_r1_3 = r1_3[:,47]
OVI_r1_3 = r1_3[:,67]
MgII_r1_3 = r1_3[:,87]
SiII_r1_3 = r1_3[:,103]
FeII_r1_3 = r1_3[:,127]

r1_4 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS1/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.4.1")
redshift = r1_4[:,0]
HI_r1_4 = r1_4[:,7]
HeI_r1_4 = r1_4[:,11]
CIV_r1_4 = r1_4[:,27]
OI_r1_4 = r1_4[:,47]
OVI_r1_4 = r1_4[:,67]
MgII_r1_4 = r1_4[:,87]
SiII_r1_4 = r1_4[:,103]
FeII_r1_4 = r1_4[:,127]

r2_2 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS2/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.2.2")
redshift = r2_2[:,0]
HI_r2_2 = r2_2[:,7]
HeI_r2_2 = r2_2[:,11]
CIV_r2_2 = r2_2[:,27]
OI_r2_2 = r2_2[:,47]
OVI_r2_2 = r2_2[:,67]
MgII_r2_2 = r2_2[:,87]
SiII_r2_2 = r2_2[:,103]
FeII_r2_2 = r2_2[:,127]

r2_4 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS2/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.4.2")
redshift = r2_4[:,0]
HI_r2_4 = r2_4[:,7]
HeI_r2_4 = r2_4[:,11]
CIV_r2_4 = r2_4[:,27]
OI_r2_4 = r2_4[:,47]
OVI_r2_4 = r2_4[:,67]
MgII_r2_4 = r2_4[:,87]
SiII_r2_4 = r2_4[:,103]
FeII_r2_4 = r2_4[:,127]

r2_6 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS2/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.6.2")
redshift = r2_6[:,0]
HI_r2_6 = r2_6[:,7]
HeI_r2_6 = r2_6[:,11]
CIV_r2_6 = r2_6[:,27]
OI_r2_6 = r2_6[:,47]
OVI_r2_6 = r2_6[:,67]
MgII_r2_6 = r2_6[:,87]
SiII_r2_6 = r2_6[:,103]
FeII_r2_6 = r2_6[:,127]

r2_8 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS2/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.8.2")
redshift = r2_8[:,0]
HI_r2_8 = r2_8[:,7]
HeI_r2_8 = r2_8[:,11]
CIV_r2_8 = r2_8[:,27]
OI_r2_8 = r2_8[:,47]
OVI_r2_8 = r2_8[:,67]
MgII_r2_8 = r2_8[:,87]
SiII_r2_8 = r2_8[:,103]
FeII_r2_8 = r2_8[:,127]

r3_1 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS3/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.1.3")
redshift = r3_1[:,0]
HI_r3_1 = r3_1[:,7]
HeI_r3_1 = r3_1[:,11]
CIV_r3_1 = r3_1[:,27]
OI_r3_1 = r3_1[:,47]
OVI_r3_1 = r3_1[:,67]
MgII_r3_1 = r3_1[:,87]
SiII_r3_1 = r3_1[:,103]
FeII_r3_1 = r3_1[:,127]

r3_2 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS3/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.2.3")
redshift = r3_2[:,0]
HI_r3_2 = r3_2[:,7]
HeI_r3_2 = r3_2[:,11]
CIV_r3_2 = r3_2[:,27]
OI_r3_2 = r3_2[:,47]
OVI_r3_2 = r3_2[:,67]
MgII_r3_2 = r3_2[:,87]
SiII_r3_2 = r3_2[:,103]
FeII_r3_2 = r3_2[:,127]

r3_3 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS3/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.3.3")
redshift = r3_3[:,0]
HI_r3_3 = r3_3[:,7]
HeI_r3_3 = r3_3[:,11]
CIV_r3_3 = r3_3[:,27]
OI_r3_3 = r3_3[:,47]
OVI_r3_3 = r3_3[:,67]
MgII_r3_3 = r3_3[:,87]
SiII_r3_3 = r3_3[:,103]
FeII_r3_3 = r3_3[:,127]

r3_6 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS3/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.6.3")
redshift = r3_6[:,0]
HI_r3_6 = r3_6[:,7]
HeI_r3_6 = r3_6[:,11]
CIV_r3_6 = r3_6[:,27]
OI_r3_6 = r3_6[:,47]
OVI_r3_6 = r3_6[:,67]
MgII_r3_6 = r3_6[:,87]
SiII_r3_6 = r3_6[:,103]
FeII_r3_6 = r3_6[:,127]

r3_9 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS3/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.9.3")
redshift = r3_9[:,0]
HI_r3_9 = r3_9[:,7]
HeI_r3_9 = r3_9[:,11]
CIV_r3_9 = r3_9[:,27]
OI_r3_9 = r3_9[:,47]
OVI_r3_9 = r3_9[:,67]
MgII_r3_9 = r3_9[:,87]
SiII_r3_9 = r3_9[:,103]
FeII_r3_9 = r3_9[:,127]

r3_12 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS3/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.12.3")
redshift = r3_12[:,0]
HI_r3_12 = r3_12[:,7]
HeI_r3_12 = r3_12[:,11]
CIV_r3_12 = r3_12[:,27]
OI_r3_12 = r3_12[:,47]
OVI_r3_12 = r3_12[:,67]
MgII_r3_12 = r3_12[:,87]
SiII_r3_12 = r3_12[:,103]
FeII_r3_12 = r3_12[:,127]

r4_1 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS4/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.1.4")
redshift = r4_1[:,0]
HI_r4_1 = r4_1[:,7]
HeI_r4_1 = r4_1[:,11]
CIV_r4_1 = r4_1[:,27]
OI_r4_1 = r4_1[:,47]
OVI_r4_1 = r4_1[:,67]
MgII_r4_1 = r4_1[:,87]
SiII_r4_1 = r4_1[:,103]
FeII_r4_1 = r4_1[:,127]

r4_2 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS4/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.2.4")
redshift = r4_2[:,0]
HI_r4_2 = r4_2[:,7]
HeI_r4_2 = r4_2[:,11]
CIV_r4_2 = r4_2[:,27]
OI_r4_2 = r4_2[:,47]
OVI_r4_2 = r4_2[:,67]
MgII_r4_2 = r4_2[:,87]
SiII_r4_2 = r4_2[:,103]
FeII_r4_2 = r4_2[:,127]

r4_4 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS4/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.4.4")
redshift = r4_4[:,0]
HI_r4_4 = r4_4[:,7]
HeI_r4_4 = r4_4[:,11]
CIV_r4_4 = r4_4[:,27]
OI_r4_4 = r4_4[:,47]
OVI_r4_4 = r4_4[:,67]
MgII_r4_4 = r4_4[:,87]
SiII_r4_4 = r4_4[:,103]
FeII_r4_4 = r4_4[:,127]

r4_6 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS4/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.6.4")
redshift = r4_6[:,0]
HI_r4_6 = r4_6[:,7]
HeI_r4_6 = r4_6[:,11]
CIV_r4_6 = r4_6[:,27]
OI_r4_6 = r4_6[:,47]
OVI_r4_6 = r4_6[:,67]
MgII_r4_6 = r4_6[:,87]
SiII_r4_6 = r4_6[:,103]
FeII_r4_6 = r4_6[:,127]

r4_8 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS4/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.8.4")
redshift = r4_8[:,0]
HI_r4_8 = r4_8[:,7]
HeI_r4_8 = r4_8[:,11]
CIV_r4_8 = r4_8[:,27]
OI_r4_8 = r4_8[:,47]
OVI_r4_8 = r4_8[:,67]
MgII_r4_8 = r4_8[:,87]
SiII_r4_8 = r4_8[:,103]
FeII_r4_8 = r4_8[:,127]

r4_10 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS4/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.10.4")
redshift = r4_10[:,0]
HI_r4_10 = r4_10[:,7]
HeI_r4_10 = r4_10[:,11]
CIV_r4_10 = r4_10[:,27]
OI_r4_10 = r4_10[:,47]
OVI_r4_10 = r4_10[:,67]
MgII_r4_10 = r4_10[:,87]
SiII_r4_10 = r4_10[:,103]
FeII_r4_10 = r4_10[:,127]

r4_12 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS4/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.12.4")
redshift = r4_12[:,0]
HI_r4_12 = r4_12[:,7]
HeI_r4_12 = r4_12[:,11]
CIV_r4_12 = r4_12[:,27]
OI_r4_12 = r4_12[:,47]
OVI_r4_12 = r4_12[:,67]
MgII_r4_12 = r4_12[:,87]
SiII_r4_12 = r4_12[:,103]
FeII_r4_12 = r4_12[:,127]

r4_14 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS4/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.14.4")
redshift = r4_14[:,0]
HI_r4_14 = r4_14[:,7]
HeI_r4_14 = r4_14[:,11]
CIV_r4_14 = r4_14[:,27]
OI_r4_14 = r4_14[:,47]
OVI_r4_14 = r4_14[:,67]
MgII_r4_14 = r4_14[:,87]
SiII_r4_14 = r4_14[:,103]
FeII_r4_14 = r4_14[:,127]

r4_16 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS4/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.16.4")
redshift = r4_16[:,0]
HI_r4_16 = r4_16[:,7]
HeI_r4_16 = r4_16[:,11]
CIV_r4_16 = r4_16[:,27]
OI_r4_16 = r4_16[:,47]
OVI_r4_16 = r4_16[:,67]
MgII_r4_16 = r4_16[:,87]
SiII_r4_16 = r4_16[:,103]
FeII_r4_16 = r4_16[:,127]

r5_1 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS5/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.1.5")
redshift = r5_1[:,0]
HI_r5_1 = r5_1[:,7]
HeI_r5_1 = r5_1[:,11]
CIV_r5_1 = r5_1[:,27]
OI_r5_1 = r5_1[:,47]
OVI_r5_1 = r5_1[:,67]
MgII_r5_1 = r5_1[:,87]
SiII_r5_1 = r5_1[:,103]
FeII_r5_1 = r5_1[:,127]

r5_3 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS5/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.3.5")
redshift = r5_3[:,0]
HI_r5_3 = r5_3[:,7]
HeI_r5_3 = r5_3[:,11]
CIV_r5_3 = r5_3[:,27]
OI_r5_3 = r5_3[:,47]
OVI_r5_3 = r5_3[:,67]
MgII_r5_3 = r5_3[:,87]
SiII_r5_3 = r5_3[:,103]
FeII_r5_3 = r5_3[:,127]

r5_4 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS5/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.4.5")
redshift = r5_4[:,0]
HI_r5_4 = r5_4[:,7]
HeI_r5_3 = r5_4[:,11]
CIV_r5_4 = r5_4[:,27]
OI_r5_4 = r5_4[:,47]
OVI_r5_4 = r5_4[:,67]
MgII_r5_4 = r5_4[:,87]
SiII_r5_4 = r5_4[:,103]
FeII_r5_4 = r5_4[:,127]

r5_7 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS5/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.7.5")
redshift = r5_7[:,0]
HI_r5_7 = r5_7[:,7]
HeI_r5_7 = r5_7[:,11]
CIV_r5_7 = r5_7[:,27]
OI_r5_7 = r5_7[:,47]
OVI_r5_7 = r5_7[:,67]
MgII_r5_7 = r5_7[:,87]
SiII_r5_7 = r5_7[:,103]
FeII_r5_7 = r5_7[:,127]

r5_8 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS5/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.8.5")
redshift = r5_8[:,0]
HI_r5_8 = r5_8[:,7]
HeI_r5_8 = r5_8[:,11]
CIV_r5_8 = r5_8[:,27]
OI_r5_8 = r5_8[:,47]
OVI_r5_8 = r5_8[:,67]
MgII_r5_8 = r5_8[:,87]
SiII_r5_8 = r5_8[:,103]
FeII_r5_8 = r5_8[:,127]

r5_11 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS5/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.11.5")
redshift = r5_11[:,0]
HI_r5_11 = r5_11[:,7]
HeI_r5_11 = r5_11[:,11]
CIV_r5_11 = r5_11[:,27]
OI_r5_11 = r5_11[:,47]
OVI_r5_11 = r5_11[:,67]
MgII_r5_11 = r5_11[:,87]
SiII_r5_11 = r5_11[:,103]
FeII_r5_11 = r5_11[:,127]

r5_12 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS5/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.12.5")
redshift = r5_12[:,0]
HI_r5_12 = r5_12[:,7]
HeI_r5_12 = r5_12[:,11]
CIV_r5_12 = r5_12[:,27]
OI_r5_12 = r5_12[:,47]
OVI_r5_12 = r5_12[:,67]
MgII_r5_12 = r5_12[:,87]
SiII_r5_12 = r5_12[:,103]
FeII_r5_12 = r5_12[:,127]

r5_15 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS5/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.15.5")
redshift = r5_15[:,0]
HI_r5_15 = r5_15[:,7]
HeI_r5_15 = r5_15[:,11]
CIV_r5_15 = r5_15[:,27]
OI_r5_15 = r5_15[:,47]
OVI_r5_15 = r5_15[:,67]
MgII_r5_15 = r5_15[:,87]
SiII_r5_15 = r5_15[:,103]
FeII_r5_15 = r5_15[:,127]

r5_16 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS5/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.16.5")
redshift = r5_16[:,0]
HI_r5_16 = r5_16[:,7]
HeI_r5_16 = r5_16[:,11]
CIV_r5_16 = r5_16[:,27]
OI_r5_16 = r5_16[:,47]
OVI_r5_16 = r5_16[:,67]
MgII_r5_16 = r5_16[:,87]
SiII_r5_16 = r5_16[:,103]
FeII_r5_16 = r5_16[:,127]

r5_19 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS5/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.19.5")
redshift = r5_19[:,0]
HI_r5_19 = r5_19[:,7]
HeI_r5_19 = r5_19[:,11]
CIV_r5_19 = r5_19[:,27]
OI_r5_19 = r5_19[:,47]
OVI_r5_19 = r5_19[:,67]
MgII_r5_19 = r5_19[:,87]
SiII_r5_19 = r5_19[:,103]
FeII_r5_19 = r5_19[:,127]

r5_20 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS5/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.20.5")
redshift = r5_20[:,0]
HI_r5_20 = r5_20[:,7]
HeI_r5_20 = r5_20[:,11]
CIV_r5_20 = r5_20[:,27]
OI_r5_20 = r5_20[:,47]
OVI_r5_20 = r5_20[:,67]
MgII_r5_20 = r5_20[:,87]
SiII_r5_20 = r5_20[:,103]
FeII_r5_20 = r5_20[:,127]

r6_1 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS6/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.1.6")
redshift = r6_1[:,0]
HI_r6_1 = r6_1[:,7]
HeI_r6_1 = r6_1[:,11]
CIV_r6_1 = r6_1[:,27]
OI_r6_1 = r6_1[:,47]
OVI_r6_1 = r6_1[:,67]
MgII_r6_1 = r6_1[:,87]
SiII_r6_1 = r6_1[:,103]
FeII_r6_1 = r6_1[:,127]

r6_2 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS6/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.2.6")
redshift = r6_2[:,0]
HI_r6_2 = r6_2[:,7]
HeI_r6_2 = r6_2[:,11]
CIV_r6_2 = r6_2[:,27]
OI_r6_2 = r6_2[:,47]
OVI_r6_2 = r6_2[:,67]
MgII_r6_2 = r6_2[:,87]
SiII_r6_2 = r6_2[:,103]
FeII_r6_2 = r6_2[:,127]

r6_5 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS6/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.5.6")
redshift = r6_5[:,0]
HI_r6_5 = r6_5[:,7]
HeI_r6_5 = r6_5[:,11]
CIV_r6_5 = r6_5[:,27]
OI_r6_5 = r6_5[:,47]
OVI_r6_5 = r6_5[:,67]
MgII_r6_5 = r6_5[:,87]
SiII_r6_5 = r6_5[:,103]
FeII_r6_5 = r6_5[:,127]

r6_7 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS6/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.7.6")
redshift = r6_7[:,0]
HI_r6_7 = r6_7[:,7]
HeI_r6_7 = r6_7[:,11]
CIV_r6_7 = r6_7[:,27]
OI_r6_7 = r6_7[:,47]
OVI_r6_7 = r6_7[:,67]
MgII_r6_7 = r6_7[:,87]
SiII_r6_7 = r6_7[:,103]
FeII_r6_7 = r6_7[:,127]

r6_9 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS6/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.9.6")
redshift = r6_9[:,0]
HI_r6_9 = r6_9[:,7]
HeI_r6_9 = r6_9[:,11]
CIV_r6_9 = r6_9[:,27]
OI_r6_9 = r6_9[:,47]
OVI_r6_9 = r6_9[:,67]
MgII_r6_9 = r6_9[:,87]
SiII_r6_9 = r6_9[:,103]
FeII_r6_9 = r6_9[:,127]

r6_11 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS6/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.11.6")
redshift = r6_11[:,0]
HI_r6_11 = r6_11[:,7]
HeI_r6_11 = r6_11[:,11]
CIV_r6_11 = r6_11[:,27]
OI_r6_11 = r6_11[:,47]
OVI_r6_11 = r6_11[:,67]
MgII_r6_11 = r6_11[:,87]
SiII_r6_11 = r6_11[:,103]
FeII_r6_11 = r6_11[:,127]

r6_13 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS6/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.13.6")
redshift = r6_13[:,0]
HI_r6_13 = r6_13[:,7]
HeI_r6_13 = r6_13[:,11]
CIV_r6_13 = r6_13[:,27]
OI_r6_13 = r6_13[:,47]
OVI_r6_13 = r6_13[:,67]
MgII_r6_13 = r6_13[:,87]
SiII_r6_13 = r6_13[:,103]
FeII_r6_13 = r6_13[:,127]

r6_15 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS6/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.15.6")
redshift = r6_15[:,0]
HI_r6_15 = r6_15[:,7]
HeI_r6_15 = r6_15[:,11]
CIV_r6_15 = r6_15[:,27]
OI_r6_15 = r6_15[:,47]
OVI_r6_15 = r6_15[:,67]
MgII_r6_15 = r6_15[:,87]
SiII_r6_15 = r6_15[:,103]
FeII_r6_15 = r6_15[:,127]

r6_17 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS6/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.17.6")
redshift = r6_17[:,0]
HI_r6_17 = r6_17[:,7]
HeI_r6_17 = r6_17[:,11]
CIV_r6_17 = r6_17[:,27]
OI_r6_17 = r6_17[:,47]
OVI_r6_17 = r6_17[:,67]
MgII_r6_17 = r6_17[:,87]
SiII_r6_17 = r6_17[:,103]
FeII_r6_17 = r6_17[:,127]

r6_19 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS6/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.19.6")
redshift = r6_19[:,0]
HI_r6_19 = r6_19[:,7]
HeI_r6_19 = r6_19[:,11]
CIV_r6_19 = r6_19[:,27]
OI_r6_19 = r6_19[:,47]
OVI_r6_19 = r6_19[:,67]
MgII_r6_19 = r6_19[:,87]
SiII_r6_19 = r6_19[:,103]
FeII_r6_19 = r6_19[:,127]

r6_21 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS6/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.21.6")
redshift = r6_21[:,0]
HI_r6_21 = r6_21[:,7]
HeI_r6_21 = r6_21[:,11]
CIV_r6_21 = r6_21[:,27]
OI_r6_21 = r6_21[:,47]
OVI_r6_21 = r6_21[:,67]
MgII_r6_21 = r6_21[:,87]
SiII_r6_21 = r6_21[:,103]
FeII_r6_21 = r6_21[:,127]

r6_22 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS6/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.22.6")
redshift = r6_22[:,0]
HI_r6_22 = r6_22[:,7]
HeI_r6_22 = r6_22[:,11]
CIV_r6_22 = r6_22[:,27]
OI_r6_22 = r6_22[:,47]
OVI_r6_22 = r6_22[:,67]
MgII_r6_22 = r6_22[:,87]
SiII_r6_22 = r6_22[:,103]
FeII_r6_22 = r6_22[:,127]

r7_1 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS7/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.1.7")
redshift = r7_1[:,0]
HI_r7_1 = r7_1[:,7]
HeI_r7_1 = r7_1[:,11]
CIV_r7_1 = r7_1[:,27]
OI_r7_1 = r7_1[:,47]
OVI_r7_1 = r7_1[:,67]
MgII_r7_1 = r7_1[:,87]
SiII_r7_1 = r7_1[:,103]
FeII_r7_1 = r7_1[:,127]

r7_3 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS7/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.3.7")
redshift = r7_3[:,0]
HI_r7_3 = r7_3[:,7]
HeI_r7_3 = r7_3[:,11]
CIV_r7_3 = r7_3[:,27]
OI_r7_3 = r7_3[:,47]
OVI_r7_3 = r7_3[:,67]
MgII_r7_3 = r7_3[:,87]
SiII_r7_3 = r7_3[:,103]
FeII_r7_3 = r7_3[:,127]

r7_4 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS7/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.4.7")
redshift = r7_4[:,0]
HI_r7_4 = r7_4[:,7]
HeI_r7_4 = r7_4[:,11]
CIV_r7_4 = r7_4[:,27]
OI_r7_4 = r7_4[:,47]
OVI_r7_4 = r7_4[:,67]
MgII_r7_4 = r7_4[:,87]
SiII_r7_4 = r7_4[:,103]
FeII_r7_4 = r7_4[:,127]

r7_6 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS7/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.6.7")
redshift = r7_6[:,0]
HI_r7_6 = r7_6[:,7]
HeI_r7_6 = r7_6[:,11]
CIV_r7_6 = r7_6[:,27]
OI_r7_6 = r7_6[:,47]
OVI_r7_6 = r7_6[:,67]
MgII_r7_6 = r7_6[:,87]
SiII_r7_6 = r7_6[:,103]
FeII_r7_6 = r7_6[:,127]

r7_8 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS7/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.8.7")
redshift = r7_8[:,0]
HI_r7_8 = r7_8[:,7]
HeI_r7_8 = r7_8[:,11]
CIV_r7_8 = r7_8[:,27]
OI_r7_8 = r7_8[:,47]
OVI_r7_8 = r7_8[:,67]
MgII_r7_8 = r7_8[:,87]
SiII_r7_8 = r7_8[:,103]
FeII_r7_8 = r7_8[:,127]

r7_10 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS7/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.10.7")
redshift = r7_10[:,0]
HI_r7_10 = r7_10[:,7]
HeI_r7_10 = r7_10[:,11]
CIV_r7_10 = r7_10[:,27]
OI_r7_10 = r7_10[:,47]
OVI_r7_10 = r7_10[:,67]
MgII_r7_10 = r7_10[:,87]
SiII_r7_10 = r7_10[:,103]
FeII_r7_10 = r7_10[:,127]

r7_12 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS7/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.12.7")
redshift = r7_12[:,0]
HI_r7_12 = r7_12[:,7]
HeI_r7_12 = r7_12[:,11]
CIV_r7_12 = r7_12[:,27]
OI_r7_12 = r7_12[:,47]
OVI_r7_12 = r7_12[:,67]
MgII_r7_12 = r7_12[:,87]
SiII_r7_12 = r7_12[:,103]
FeII_r7_12 = r7_12[:,127]

r7_15 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS7/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.15.7")
redshift = r7_15[:,0]
HI_r7_15 = r7_15[:,7]
HeI_r7_15 = r7_15[:,11]
CIV_r7_15 = r7_15[:,27]
OI_r7_15 = r7_15[:,47]
OVI_r7_15 = r7_15[:,67]
MgII_r7_15 = r7_15[:,87]
SiII_r7_15 = r7_15[:,103]
FeII_r7_15 = r7_15[:,127]

r7_18 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS7/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.18.7")
redshift = r7_18[:,0]
HI_r7_18 = r7_18[:,7]
HeI_r7_18 = r7_18[:,11]
CIV_r7_18 = r7_18[:,27]
OI_r7_18 = r7_18[:,47]
OVI_r7_18 = r7_18[:,67]
MgII_r7_18 = r7_18[:,87]
SiII_r7_18 = r7_18[:,103]
FeII_r7_18 = r7_18[:,127]

r7_20 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS7/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.20.7")
redshift = r7_20[:,0]
HI_r7_20 = r7_20[:,7]
HeI_r7_20 = r7_20[:,11]
CIV_r7_20 = r7_20[:,27]
OI_r7_20 = r7_20[:,47]
OVI_r7_20 = r7_20[:,67]
MgII_r7_20 = r7_20[:,87]
SiII_r7_20 = r7_20[:,103]
FeII_r7_20 = r7_20[:,127]

r7_23 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS7/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.23.7")
redshift = r7_23[:,0]
HI_r7_23 = r7_23[:,7]
HeI_r7_23 = r7_23[:,11]
CIV_r7_23 = r7_23[:,27]
OI_r7_23 = r7_23[:,47]
OVI_r7_23 = r7_23[:,67]
MgII_r7_23 = r7_23[:,87]
SiII_r7_23 = r7_23[:,103]
FeII_r7_23 = r7_23[:,127]

r7_26 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS7/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.26.7")
redshift = r7_26[:,0]
HI_r7_26 = r7_26[:,7]
HeI_r7_26 = r7_26[:,11]
CIV_r7_26 = r7_26[:,27]
OI_r7_26 = r7_26[:,47]
OVI_r7_26 = r7_26[:,67]
MgII_r7_26 = r7_26[:,87]
SiII_r7_26 = r7_26[:,103]
FeII_r7_26 = r7_26[:,127]

r7_28 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS7/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.28.7")
redshift = r7_28[:,0]
HI_r7_28 = r7_28[:,7]
HeI_r7_28 = r7_28[:,11]
CIV_r7_28 = r7_28[:,27]
OI_r7_28 = r7_28[:,47]
OVI_r7_28 = r7_28[:,67]
MgII_r7_28 = r7_28[:,87]
SiII_r7_28 = r7_28[:,103]
FeII_r7_28 = r7_28[:,127]

r8_1 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS8/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.1.8")
redshift = r8_1[:,0]
HI_r8_1 = r8_1[:,7]
HeI_r8_1 = r8_1[:,11]
CIV_r8_1 = r8_1[:,27]
OI_r8_1 = r8_1[:,47]
OVI_r8_1 = r8_1[:,67]
MgII_r8_1 = r8_1[:,87]
SiII_r8_1 = r8_1[:,103]
FeII_r8_1 = r8_1[:,127]

r8_4 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS8/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.4.8")
redshift = r8_4[:,0]
HI_r8_4 = r8_4[:,7]
HeI_r8_4 = r8_4[:,11]
CIV_r8_4 = r8_4[:,27]
OI_r8_4 = r8_4[:,47]
OVI_r8_4 = r8_4[:,67]
MgII_r8_4 = r8_4[:,87]
SiII_r8_4 = r8_4[:,103]
FeII_r8_4 = r8_4[:,127]

r8_6 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS8/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.6.8")
redshift = r8_6[:,0]
HI_r8_6 = r8_6[:,7]
HeI_r8_6 = r8_6[:,11]
CIV_r8_6 = r8_6[:,27]
OI_r8_6 = r8_6[:,47]
OVI_r8_6 = r8_6[:,67]
MgII_r8_6 = r8_6[:,87]
SiII_r8_6 = r8_6[:,103]
FeII_r8_6 = r8_6[:,127]

r8_8 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS8/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.8.8")
redshift = r8_8[:,0]
HI_r8_8 = r8_8[:,7]
HeI_r8_8 = r8_8[:,11]
CIV_r8_8 = r8_8[:,27]
OI_r8_8 = r8_8[:,47]
OVI_r8_8 = r8_8[:,67]
MgII_r8_8 = r8_8[:,87]
SiII_r8_8 = r8_8[:,103]
FeII_r8_8 = r8_8[:,127]

r8_10 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS8/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.10.8")
redshift = r8_10[:,0]
HI_r8_10 = r8_10[:,7]
HeI_r8_10 = r8_10[:,11]
CIV_r8_10 = r8_10[:,27]
OI_r8_10 = r8_10[:,47]
OVI_r8_10 = r8_10[:,67]
MgII_r8_10 = r8_10[:,87]
SiII_r8_10 = r8_10[:,103]
FeII_r8_10 = r8_10[:,127]

r8_12 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS8/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.12.8")
redshift = r8_12[:,0]
HI_r8_12 = r8_12[:,7]
HeI_r8_12 = r8_12[:,11]
CIV_r8_12 = r8_12[:,27]
OI_r8_12 = r8_12[:,47]
OVI_r8_12 = r8_12[:,67]
MgII_r8_12 = r8_12[:,87]
SiII_r8_12 = r8_12[:,103]
FeII_r8_12 = r8_12[:,127]

r8_14 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS8/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.14.8")
redshift = r8_14[:,0]
HI_r8_14 = r8_14[:,7]
HeI_r8_14 = r8_14[:,11]
CIV_r8_14 = r8_14[:,27]
OI_r8_14 = r8_14[:,47]
OVI_r8_14 = r8_14[:,67]
MgII_r8_14 = r8_14[:,87]
SiII_r8_14 = r8_14[:,103]
FeII_r8_14 = r8_14[:,127]

r8_16 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS8/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.16.8")
redshift = r8_16[:,0]
HI_r8_16 = r8_16[:,7]
HeI_r8_16 = r8_16[:,11]
CIV_r8_16 = r8_16[:,27]
OI_r8_16 = r8_16[:,47]
OVI_r8_16 = r8_16[:,67]
MgII_r8_16 = r8_16[:,87]
SiII_r8_16 = r8_16[:,103]
FeII_r8_16 = r8_16[:,127]

r8_18 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS8/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.18.8")
redshift = r8_18[:,0]
HI_r8_18 = r8_18[:,7]
HeI_r8_18 = r8_18[:,11]
CIV_r8_18 = r8_18[:,27]
OI_r8_18 = r8_18[:,47]
OVI_r8_18 = r8_18[:,67]
MgII_r8_18 = r8_18[:,87]
SiII_r8_18 = r8_18[:,103]
FeII_r8_18 = r8_18[:,127]

r8_20 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS8/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.20.8")
redshift = r8_20[:,0]
HI_r8_20 = r8_20[:,7]
HeI_r8_20 = r8_20[:,11]
CIV_r8_20 = r8_20[:,27]
OI_r8_20 = r8_20[:,47]
OVI_r8_20 = r8_20[:,67]
MgII_r8_20 = r8_20[:,87]
SiII_r8_20 = r8_20[:,103]
FeII_r8_20 = r8_20[:,127]

r8_22 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS8/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.22.8")
redshift = r8_22[:,0]
HI_r8_22 = r8_22[:,7]
HeI_r8_22 = r8_22[:,11]
CIV_r8_22 = r8_22[:,27]
OI_r8_22 = r8_22[:,47]
OVI_r8_22 = r8_22[:,67]
MgII_r8_22 = r8_22[:,87]
SiII_r8_22 = r8_22[:,103]
FeII_r8_22 = r8_22[:,127]

r8_24 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS8/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.24.8")
redshift = r8_24[:,0]
HI_r8_24 = r8_24[:,7]
HeI_r8_24 = r8_24[:,11]
CIV_r8_24 = r8_24[:,27]
OI_r8_24 = r8_24[:,47]
OVI_r8_24 = r8_24[:,67]
MgII_r8_24 = r8_24[:,87]
SiII_r8_24 = r8_24[:,103]
FeII_r8_24 = r8_24[:,127]

r8_26 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS8/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.26.8")
redshift = r8_26[:,0]
HI_r8_26 = r8_26[:,7]
HeI_r8_26 = r8_26[:,11]
CIV_r8_26 = r8_26[:,27]
OI_r8_26 = r8_26[:,47]
OVI_r8_26 = r8_26[:,67]
MgII_r8_26 = r8_26[:,87]
SiII_r8_26 = r8_26[:,103]
FeII_r8_26 = r8_26[:,127]

r8_28 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS8/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.28.8")
redshift = r8_28[:,0]
HI_r8_28 = r8_28[:,7]
HeI_r8_28 = r8_28[:,11]
CIV_r8_28 = r8_28[:,27]
OI_r8_28 = r8_28[:,47]
OVI_r8_28 = r8_28[:,67]
MgII_r8_28 = r8_28[:,87]
SiII_r8_28 = r8_28[:,103]
FeII_r8_28 = r8_28[:,127]

r8_30 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS8/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.30.8")
redshift = r8_30[:,0]
HI_r8_30 = r8_30[:,7]
HeI_r8_30 = r8_30[:,11]
CIV_r8_30 = r8_30[:,27]
OI_r8_30 = r8_30[:,47]
OVI_r8_30 = r8_30[:,67]
MgII_r8_30 = r8_30[:,87]
SiII_r8_30 = r8_30[:,103]
FeII_r8_30 = r8_30[:,127]

r8_32 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS8/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.32.8")
redshift = r8_32[:,0]
HI_r8_32 = r8_32[:,7]
HeI_r8_32 = r8_32[:,11]
CIV_r8_32 = r8_32[:,27]
OI_r8_32 = r8_32[:,47]
OVI_r8_32 = r8_32[:,67]
MgII_r8_32 = r8_32[:,87]
SiII_r8_32 = r8_32[:,103]
FeII_r8_32 = r8_32[:,127]

r9_1 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS9/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.1.9")
redshift = r9_1[:,0]
HI_r9_1 = r9_1[:,7]
HeI_r9_1 = r9_1[:,11]
CIV_r9_1 = r9_1[:,27]
OI_r9_1 = r9_1[:,47]
OVI_r9_1 = r9_1[:,67]
MgII_r9_1 = r9_1[:,87]
SiII_r9_1 = r9_1[:,103]
FeII_r9_1 = r9_1[:,127]

r9_4 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS9/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.4.9")
redshift = r9_4[:,0]
HI_r9_4 = r9_4[:,7]
HeI_r9_4 = r9_4[:,11]
CIV_r9_4 = r9_4[:,27]
OI_r9_4 = r9_4[:,47]
OVI_r9_4 = r9_4[:,67]
MgII_r9_4 = r9_4[:,87]
SiII_r9_4 = r9_4[:,103]
FeII_r9_4 = r9_4[:,127]

r9_6 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS9/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.6.9")
redshift = r9_6[:,0]
HI_r9_6 = r9_6[:,7]
HeI_r9_6 = r9_6[:,11]
CIV_r9_6 = r9_6[:,27]
OI_r9_6 = r9_6[:,47]
OVI_r9_6 = r9_6[:,67]
MgII_r9_6 = r9_6[:,87]
SiII_r9_6 = r9_6[:,103]
FeII_r9_6 = r9_6[:,127]

r9_8 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS9/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.8.9")
redshift = r9_8[:,0]
HI_r9_8 = r9_8[:,7]
HeI_r9_8 = r9_8[:,11]
CIV_r9_8 = r9_8[:,27]
OI_r9_8 = r9_8[:,47]
OVI_r9_8 = r9_8[:,67]
MgII_r9_8 = r9_8[:,87]
SiII_r9_8 = r9_8[:,103]
FeII_r9_8 = r9_8[:,127]

r9_10 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS9/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.10.9")
redshift = r9_10[:,0]
HI_r9_10 = r9_10[:,7]
HeI_r9_10 = r9_10[:,11]
CIV_r9_10 = r9_10[:,27]
OI_r9_10 = r9_10[:,47]
OVI_r9_10 = r9_10[:,67]
MgII_r9_10 = r9_10[:,87]
SiII_r9_10 = r9_10[:,103]
FeII_r9_10 = r9_10[:,127]

r9_12 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS9/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.12.9")
redshift = r9_12[:,0]
HI_r9_12 = r9_12[:,7]
HeI_r9_12 = r9_12[:,11]
CIV_r9_12 = r9_12[:,27]
OI_r9_12 = r9_12[:,47]
OVI_r9_12 = r9_12[:,67]
MgII_r9_12 = r9_12[:,87]
SiII_r9_12 = r9_12[:,103]
FeII_r9_12 = r9_12[:,127]

r9_14 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS9/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.14.9")
redshift = r9_14[:,0]
HI_r9_14 = r9_14[:,7]
HeI_r9_14 = r9_14[:,11]
CIV_r9_14 = r9_14[:,27]
OI_r9_14 = r9_14[:,47]
OVI_r9_14 = r9_14[:,67]
MgII_r9_14 = r9_14[:,87]
SiII_r9_14 = r9_14[:,103]
FeII_r9_14 = r9_14[:,127]

r9_16 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS9/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.16.9")
redshift = r9_16[:,0]
HI_r9_16 = r9_16[:,7]
HeI_r9_16 = r9_16[:,11]
CIV_r9_16 = r9_16[:,27]
OI_r9_16 = r9_16[:,47]
OVI_r9_16 = r9_16[:,67]
MgII_r9_16 = r9_16[:,87]
SiII_r9_16 = r9_16[:,103]
FeII_r9_16 = r9_16[:,127]

r9_18 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS9/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.18.9")
redshift = r9_18[:,0]
HI_r9_18 = r9_18[:,7]
HeI_r9_18 = r9_18[:,11]
CIV_r9_18 = r9_18[:,27]
OI_r9_18 = r9_18[:,47]
OVI_r9_18 = r9_18[:,67]
MgII_r9_18 = r9_18[:,87]
SiII_r9_18 = r9_18[:,103]
FeII_r9_18 = r9_18[:,127]

r9_20 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS9/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.20.9")
redshift = r9_20[:,0]
HI_r9_20 = r9_20[:,7]
HeI_r9_20 = r9_20[:,11]
CIV_r9_20 = r9_20[:,27]
OI_r9_20 = r9_20[:,47]
OVI_r9_20 = r9_20[:,67]
MgII_r9_20 = r9_20[:,87]
SiII_r9_20 = r9_20[:,103]
FeII_r9_20 = r9_20[:,127]

r9_22 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS9/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.22.9")
redshift = r9_22[:,0]
HI_r9_22 = r9_22[:,7]
HeI_r9_22 = r9_22[:,11]
CIV_r9_22 = r9_22[:,27]
OI_r9_22 = r9_22[:,47]
OVI_r9_22 = r9_22[:,67]
MgII_r9_22 = r9_22[:,87]
SiII_r9_22 = r9_22[:,103]
FeII_r9_22 = r9_22[:,127]

r9_24 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS9/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.24.9")
redshift = r9_24[:,0]
HI_r9_24 = r9_24[:,7]
HeI_r9_24 = r9_24[:,11]
CIV_r9_24 = r9_24[:,27]
OI_r9_24 = r9_24[:,47]
OVI_r9_24 = r9_24[:,67]
MgII_r9_24 = r9_24[:,87]
SiII_r9_24 = r9_24[:,103]
FeII_r9_24 = r9_24[:,127]

r9_26 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS9/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.26.9")
redshift = r9_26[:,0]
HI_r9_26 = r9_26[:,7]
HeI_r9_26 = r9_26[:,11]
CIV_r9_26 = r9_26[:,27]
OI_r9_26 = r9_26[:,47]
OVI_r9_26 = r9_26[:,67]
MgII_r9_26 = r9_26[:,87]
SiII_r9_26 = r9_26[:,103]
FeII_r9_26 = r9_26[:,127]

r9_28 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS9/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.28.9")
redshift = r9_28[:,0]
HI_r9_28 = r9_28[:,7]
HeI_r9_28 = r9_28[:,11]
CIV_r9_28 = r9_28[:,27]
OI_r9_28 = r9_28[:,47]
OVI_r9_28 = r9_28[:,67]
MgII_r9_28 = r9_28[:,87]
SiII_r9_28 = r9_28[:,103]
FeII_r9_28 = r9_28[:,127]

r9_30 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS9/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.30.9")
redshift = r9_30[:,0]
HI_r9_30 = r9_30[:,7]
HeI_r9_30 = r9_30[:,11]
CIV_r9_30 = r9_30[:,27]
OI_r9_30 = r9_30[:,47]
OVI_r9_30 = r9_30[:,67]
MgII_r9_30 = r9_30[:,87]
SiII_r9_30 = r9_30[:,103]
FeII_r9_30 = r9_30[:,127]

r9_32 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS9/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.32.9")
redshift = r9_32[:,0]
HI_r9_32 = r9_32[:,7]
HeI_r9_32 = r9_32[:,11]
CIV_r9_32 = r9_32[:,27]
OI_r9_32 = r9_32[:,47]
OVI_r9_32 = r9_32[:,67]
MgII_r9_32 = r9_32[:,87]
SiII_r9_32 = r9_32[:,103]
FeII_r9_32 = r9_32[:,127]

r9_34 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS9/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.34.9")
redshift = r9_34[:,0]
HI_r9_34 = r9_34[:,7]
HeI_r9_34 = r9_34[:,11]
CIV_r9_34 = r9_34[:,27]
OI_r9_34 = r9_34[:,47]
OVI_r9_34 = r9_34[:,67]
MgII_r9_34 = r9_34[:,87]
SiII_r9_34 = r9_34[:,103]
FeII_r9_34 = r9_34[:,127]

r9_36 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS9/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.36.9")
redshift = r9_36[:,0]
HI_r9_36 = r9_36[:,7]
HeI_r9_36 = r9_36[:,11]
CIV_r9_36 = r9_36[:,27]
OI_r9_36 = r9_36[:,47]
OVI_r9_36 = r9_36[:,67]
MgII_r9_36 = r9_36[:,87]
SiII_r9_36 = r9_36[:,103]
FeII_r9_36 = r9_36[:,127]

r10_1 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS10/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.1.10")
redshift = r10_1[:,0]
HI_r10_1 = r10_1[:,7]
HeI_r10_1 = r10_1[:,11]
CIV_r10_1 = r10_1[:,27]
OI_r10_1 = r10_1[:,47]
OVI_r10_1 = r10_1[:,67]
MgII_r10_1 = r10_1[:,87]
SiII_r10_1 = r10_1[:,103]
FeII_r10_1 = r10_1[:,127]

r10_4 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS10/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.4.10")
redshift = r10_4[:,0]
HI_r10_4 = r10_4[:,7]
HeI_r10_4 = r10_4[:,11]
CIV_r10_4 = r10_4[:,27]
OI_r10_4 = r10_4[:,47]
OVI_r10_4 = r10_4[:,67]
MgII_r10_4 = r10_4[:,87]
SiII_r10_4 = r10_4[:,103]
FeII_r10_4 = r10_4[:,127]

r10_6 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS10/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.6.10")
redshift = r10_6[:,0]
HI_r10_6 = r10_6[:,7]
HeI_r10_6 = r10_6[:,11]
CIV_r10_6 = r10_6[:,27]
OI_r10_6 = r10_6[:,47]
OVI_r10_6 = r10_6[:,67]
MgII_r10_6 = r10_6[:,87]
SiII_r10_6 = r10_6[:,103]
FeII_r10_6 = r10_6[:,127]

r10_8 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS10/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.8.10")
redshift = r10_8[:,0]
HI_r10_8 = r10_8[:,7]
HeI_r10_8 = r10_8[:,11]
CIV_r10_8 = r10_8[:,27]
OI_r10_8 = r10_8[:,47]
OVI_r10_8 = r10_8[:,67]
MgII_r10_8 = r10_8[:,87]
SiII_r10_8 = r10_8[:,103]
FeII_r10_8 = r10_8[:,127]

r10_10 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS10/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.10.10")
redshift = r10_10[:,0]
HI_r10_10 = r10_10[:,7]
HeI_r10_10 = r10_10[:,11]
CIV_r10_10 = r10_10[:,27]
OI_r10_10 = r10_10[:,47]
OVI_r10_10 = r10_10[:,67]
MgII_r10_10 = r10_10[:,87]
SiII_r10_10 = r10_10[:,103]
FeII_r10_10 = r10_10[:,127]

r10_12 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS10/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.12.10")
redshift = r10_12[:,0]
HI_r10_12 = r10_12[:,7]
HeI_r10_12 = r10_12[:,11]
CIV_r10_12 = r10_12[:,27]
OI_r10_12 = r10_12[:,47]
OVI_r10_12 = r10_12[:,67]
MgII_r10_12 = r10_12[:,87]
SiII_r10_12 = r10_12[:,103]
FeII_r10_12 = r10_12[:,127]

r10_14 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS10/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.14.10")
redshift = r10_14[:,0]
HI_r10_14 = r10_14[:,7]
HeI_r10_14 = r10_14[:,11]
CIV_r10_14 = r10_14[:,27]
OI_r10_14 = r10_14[:,47]
OVI_r10_14 = r10_14[:,67]
MgII_r10_14 = r10_14[:,87]
SiII_r10_14 = r10_14[:,103]
FeII_r10_14 = r10_14[:,127]

r10_16 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS10/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.16.10")
redshift = r10_16[:,0]
HI_r10_16 = r10_16[:,7]
HeI_r10_16 = r10_16[:,11]
CIV_r10_16 = r10_16[:,27]
OI_r10_16 = r10_16[:,47]
OVI_r10_16 = r10_16[:,67]
MgII_r10_16 = r10_16[:,87]
SiII_r10_16 = r10_16[:,103]
FeII_r10_16 = r10_16[:,127]

r10_18 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS10/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.18.10")
redshift = r10_18[:,0]
HI_r10_18 = r10_18[:,7]
HeI_r10_18 = r10_18[:,11]
CIV_r10_18 = r10_18[:,27]
OI_r10_18 = r10_18[:,47]
OVI_r10_18 = r10_18[:,67]
MgII_r10_18 = r10_18[:,87]
SiII_r10_18 = r10_18[:,103]
FeII_r10_18 = r10_18[:,127]

r10_20 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS10/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.20.10")
redshift = r10_20[:,0]
HI_r10_20 = r10_20[:,7]
HeI_r10_20 = r10_20[:,11]
CIV_r10_20 = r10_20[:,27]
OI_r10_20 = r10_20[:,47]
OVI_r10_20 = r10_20[:,67]
MgII_r10_20 = r10_20[:,87]
SiII_r10_20 = r10_20[:,103]
FeII_r10_20 = r10_20[:,127]

r10_22 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS10/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.22.10")
redshift = r10_22[:,0]
HI_r10_22= r10_22[:,7]
HeI_r10_22 = r10_22[:,11]
CIV_r10_22 = r10_22[:,27]
OI_r10_22 = r10_22[:,47]
OVI_r10_22 = r10_22[:,67]
MgII_r10_22 = r10_22[:,87]
SiII_r10_22 = r10_22[:,103]
FeII_r10_22 = r10_22[:,127]

r10_24 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS10/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.24.10")
redshift = r10_24[:,0]
HI_r10_24 = r10_24[:,7]
HeI_r10_24 = r10_24[:,11]
CIV_r10_24 = r10_24[:,27]
OI_r10_24 = r10_24[:,47]
OVI_r10_24 = r10_24[:,67]
MgII_r10_24 = r10_24[:,87]
SiII_r10_24 = r10_24[:,103]
FeII_r10_24 = r10_24[:,127]

r10_26 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS10/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.26.10")
redshift = r10_26[:,0]
HI_r10_26 = r10_26[:,7]
HeI_r10_26 = r10_26[:,11]
CIV_r10_26 = r10_26[:,27]
OI_r10_26 = r10_26[:,47]
OVI_r10_26 = r10_26[:,67]
MgII_r10_26 = r10_26[:,87]
SiII_r10_26 = r10_26[:,103]
FeII_r10_26 = r10_26[:,127]

r10_28 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS10/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.28.10")
redshift = r10_28[:,0]
HI_r10_28 = r10_28[:,7]
HeI_r10_28 = r10_28[:,11]
CIV_r10_28 = r10_28[:,27]
OI_r10_28 = r10_28[:,47]
OVI_r10_28 = r10_28[:,67]
MgII_r10_28 = r10_28[:,87]
SiII_r10_28 = r10_28[:,103]
FeII_r10_28 = r10_28[:,127]

r10_30 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS10/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.30.10")
redshift = r10_30[:,0]
HI_r10_30 = r10_30[:,7]
HeI_r10_30 = r10_30[:,11]
CIV_r10_30 = r10_30[:,27]
OI_r10_30 = r10_30[:,47]
OVI_r10_30 = r10_30[:,67]
MgII_r10_30 = r10_30[:,87]
SiII_r10_30 = r10_30[:,103]
FeII_r10_30 = r10_30[:,127]

r10_32 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS10/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.32.10")
redshift = r10_32[:,0]
HI_r10_32 = r10_32[:,7]
HeI_r10_32 = r10_32[:,11]
CIV_r10_32 = r10_32[:,27]
OI_r10_32 = r10_32[:,47]
OVI_r10_32 = r10_32[:,67]
MgII_r10_32 = r10_32[:,87]
SiII_r10_32 = r10_32[:,103]
FeII_r10_32 = r10_32[:,127]

r10_34 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS10/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.34.10")
redshift = r10_34[:,0]
HI_r10_34 = r10_34[:,7]
HeI_r10_34 = r10_34[:,11]
CIV_r10_34 = r10_34[:,27]
OI_r10_34 = r10_34[:,47]
OVI_r10_34 = r10_34[:,67]
MgII_r10_34 = r10_34[:,87]
SiII_r10_34 = r10_34[:,103]
FeII_r10_34 = r10_34[:,127]

r10_36 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS10/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.36.10")
redshift = r10_36[:,0]
HI_r10_36 = r10_36[:,7]
HeI_r10_36 = r10_36[:,11]
CIV_r10_36 = r10_36[:,27]
OI_r10_36 = r10_36[:,47]
OVI_r10_36 = r10_36[:,67]
MgII_r10_36 = r10_36[:,87]
SiII_r10_36 = r10_36[:,103]
FeII_r10_36 = r10_36[:,127]

r10_38 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS10/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.38.10")
redshift = r10_38[:,0]
HI_r10_38 = r10_38[:,7]
HeI_r10_38 = r10_38[:,11]
CIV_r10_38 = r10_38[:,27]
OI_r10_38 = r10_38[:,47]
OVI_r10_38 = r10_38[:,67]
MgII_r10_38 = r10_38[:,87]
SiII_r10_38 = r10_38[:,103]
FeII_r10_38 = r10_38[:,127]

r10_40 = np.loadtxt("/home/sheehank/Map_research/MAP_Research/SimulationData/Cptmarvel/Halo1/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.dir/LOS10/specaim.cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.LOSFile.40.10")
redshift = r10_40[:,0]
HI_r10_40 = r10_40[:,7]
HeI_r10_40 = r10_40[:,11]
CIV_r10_40 = r10_40[:,27]
OI_r10_40 = r10_40[:,47]
OVI_r10_40 = r10_40[:,67]
MgII_r10_40 = r10_40[:,87]
SiII_r10_40 = r10_40[:,103]
FeII_r10_40 = r10_40[:,127]


# **Hydrogen I**

# In[3]:

make_intensity (redshift, HI_r0, HI_r1_1, HI_r1_2, HI_r1_3, HI_r1_4, HI_r2_2, HI_r2_4, HI_r2_6, HI_r2_8, HI_r3_1, HI_r3_2, HI_r3_3, HI_r3_6, HI_r3_9, HI_r3_12, HI_r4_1, HI_r4_2, HI_r4_4, HI_r4_6, HI_r4_8, HI_r4_10, HI_r4_12, HI_r4_14, HI_r4_16, HI_r5_1, HI_r5_3, HI_r5_4, HI_r5_7, HI_r5_8, HI_r5_11, HI_r5_12, HI_r5_15, HI_r5_16, HI_r5_19, HI_r5_20, HI_r6_1, HI_r6_2, HI_r6_5, HI_r6_7, HI_r6_9, HI_r6_11, HI_r6_13, HI_r6_15, HI_r6_17, HI_r6_19, HI_r6_21, HI_r6_22, HI_r7_1, HI_r7_3, HI_r7_4, HI_r7_6, HI_r7_8, HI_r7_10, HI_r7_12, HI_r7_15, HI_r7_18, HI_r7_20, HI_r7_23, HI_r7_26, HI_r7_28, HI_r8_1, HI_r8_4, HI_r8_6, HI_r8_8, HI_r8_10, HI_r8_12, HI_r8_14, HI_r8_16, HI_r8_18, HI_r8_20, HI_r8_22, HI_r8_24, HI_r8_26, HI_r8_28, HI_r8_30, HI_r8_32, HI_r9_1, HI_r9_4, HI_r9_6, HI_r9_8, HI_r9_10, HI_r9_12, HI_r9_14, HI_r9_16, HI_r9_18, HI_r9_20, HI_r9_22, HI_r9_24, HI_r9_26, HI_r9_28, HI_r9_30, HI_r9_32, HI_r9_34, HI_r9_36, HI_r10_1, HI_r10_4, HI_r10_6, HI_r10_8, HI_r10_10, HI_r10_12, HI_r10_14, HI_r10_16, HI_r10_18, HI_r10_20, HI_r10_22, HI_r10_24, HI_r10_26, HI_r10_28, HI_r10_30, HI_r10_32, HI_r10_34, HI_r10_36, HI_r10_38, HI_r10_40)


# In[4]:

HI_width = [5.1428539724938146e-05, 0.00011121764482738497, 0.00013454049892814567, 0.00014601444537996316, 0.00013899783091470995, 0.00013194587810028068, 0.00012863481121124727, 0.00012239304715220921, 0.00011843440294819343, 0.0001149884103039194, 0.00011142796507729387]
HI_SD = [0.0, 2.5673014343785385e-06, 6.5214459892204986e-06, 9.5135400985214572e-06, 8.3129325309790935e-06, 5.0670774909211386e-06, 7.6205377824422522e-06, 4.691481179307774e-06, 1.8199468097906216e-06, 3.03708078706868e-06, 3.4155598037487991e-06]


# In[5]:

plot(radius, HI_width)
errorbar(radius, HI_width, yerr = HI_SD, fmt = 'bo')
title("Hydrogen I Equivalent Width\n Z = 0")
xlabel("Radius (R/Rvir)")
ylabel("Equivalent Width")


# **Carcon IV**

# In[20]:

make_intensity(redshift, OVI_r0, OVI_r1_1, OVI_r1_2, OVI_r1_3, OVI_r1_4, OVI_r2_2, OVI_r2_4, OVI_r2_6, OVI_r2_8, OVI_r3_1, OVI_r3_2, OVI_r3_3, OVI_r3_6, OVI_r3_9, OVI_r3_12, OVI_r4_1, OVI_r4_2, OVI_r4_4, OVI_r4_6, OVI_r4_8, OVI_r4_10, OVI_r4_12, OVI_r4_14, OVI_r4_16, OVI_r5_1, OVI_r5_3, OVI_r5_4, OVI_r5_7, OVI_r5_8, OVI_r5_11, OVI_r5_12, OVI_r5_15, OVI_r5_16, OVI_r5_19, OVI_r5_20, OVI_r6_1, OVI_r6_2, OVI_r6_5, OVI_r6_7, OVI_r6_9, OVI_r6_11, OVI_r6_13, OVI_r6_15, OVI_r6_17, OVI_r6_19, OVI_r6_21, OVI_r6_22, OVI_r7_1, OVI_r7_3, OVI_r7_4, OVI_r7_6, OVI_r7_8, OVI_r7_10, OVI_r7_12, OVI_r7_15, OVI_r7_18, OVI_r7_20, OVI_r7_23, OVI_r7_26, OVI_r7_28, OVI_r8_1, OVI_r8_4, OVI_r8_6, OVI_r8_8, OVI_r8_10, OVI_r8_12, OVI_r8_14, OVI_r8_16, OVI_r8_18, OVI_r8_20, OVI_r8_22, OVI_r8_24, OVI_r8_26, OVI_r8_28, OVI_r8_30, OVI_r8_32, OVI_r9_1, OVI_r9_4, OVI_r9_6, OVI_r9_8, OVI_r9_10, OVI_r9_12, OVI_r9_14, OVI_r9_16, OVI_r9_18, OVI_r9_20, OVI_r9_22, OVI_r9_24, OVI_r9_26, OVI_r9_28, OVI_r9_30, OVI_r9_32, OVI_r9_34, OVI_r9_36, OVI_r10_1, OVI_r10_4, OVI_r10_6, OVI_r10_8, OVI_r10_10, OVI_r10_12, OVI_r10_14, OVI_r10_16, OVI_r10_18, OVI_r10_20, OVI_r10_22, OVI_r10_24, OVI_r10_26, OVI_r10_28, OVI_r10_30, OVI_r10_32, OVI_r10_34, OVI_r10_36, OVI_r10_38, OVI_r10_40)


# In[21]:

CIV_width = [9.414924810876614e-05, 7.8696641285401226e-05, 7.791847379015569e-05, 7.2006564252225606e-05, 6.706523541428752e-05, 6.0851059211055571e-05, 4.5281864496212018e-05, 4.1921515759684886e-05, 3.4038495250205607e-05, 3.3388707397500223e-05, 3.1429558301990854e-05]
CIV_SD = [0.0, 2.4743216188748063e-05, 1.6982418742215769e-05, 1.9608000156097438e-05, 1.4919990777198485e-05, 1.3527306543761415e-05, 7.1612811777481039e-06, 1.0282464029715337e-05, 5.8506862234880293e-06, 4.855767316517154e-06, 3.987779881685937e-06]


# In[23]:

plot(radius, CIV_width)
errorbar(radius, CIV_width, yerr = CIV_SD, fmt = 'bo')
title("Carbon IV Equivalent Width\n Z = 0")
xlabel("Radius (R/Rvir)")
ylabel("Equivalent Width")


# **Oxygen VI**

# In[10]:

make_intensity(redshift, OVI_r0, OVI_r1_1, OVI_r1_2, OVI_r1_3, OVI_r1_4, OVI_r2_2, OVI_r2_4, OVI_r2_6, OVI_r2_8, OVI_r3_1, OVI_r3_2, OVI_r3_3, OVI_r3_6, OVI_r3_9, OVI_r3_12, OVI_r4_1, OVI_r4_2, OVI_r4_4, OVI_r4_6, OVI_r4_8, OVI_r4_10, OVI_r4_12, OVI_r4_14, OVI_r4_16, OVI_r5_1, OVI_r5_3, OVI_r5_4, OVI_r5_7, OVI_r5_8, OVI_r5_11, OVI_r5_12, OVI_r5_15, OVI_r5_16, OVI_r5_19, OVI_r5_20, OVI_r6_1, OVI_r6_2, OVI_r6_5, OVI_r6_7, OVI_r6_9, OVI_r6_11, OVI_r6_13, OVI_r6_15, OVI_r6_17, OVI_r6_19, OVI_r6_21, OVI_r6_22, OVI_r7_1, OVI_r7_3, OVI_r7_4, OVI_r7_6, OVI_r7_8, OVI_r7_10, OVI_r7_12, OVI_r7_15, OVI_r7_18, OVI_r7_20, OVI_r7_23, OVI_r7_26, OVI_r7_28, OVI_r8_1, OVI_r8_4, OVI_r8_6, OVI_r8_8, OVI_r8_10, OVI_r8_12, OVI_r8_14, OVI_r8_16, OVI_r8_18, OVI_r8_20, OVI_r8_22, OVI_r8_24, OVI_r8_26, OVI_r8_28, OVI_r8_30, OVI_r8_32, OVI_r9_1, OVI_r9_4, OVI_r9_6, OVI_r9_8, OVI_r9_10, OVI_r9_12, OVI_r9_14, OVI_r9_16, OVI_r9_18, OVI_r9_20, OVI_r9_22, OVI_r9_24, OVI_r9_26, OVI_r9_28, OVI_r9_30, OVI_r9_32, OVI_r9_34, OVI_r9_36, OVI_r10_1, OVI_r10_4, OVI_r10_6, OVI_r10_8, OVI_r10_10, OVI_r10_12, OVI_r10_14, OVI_r10_16, OVI_r10_18, OVI_r10_20, OVI_r10_22, OVI_r10_24, OVI_r10_26, OVI_r10_28, OVI_r10_30, OVI_r10_32, OVI_r10_34, OVI_r10_36, OVI_r10_38, OVI_r10_40)


# In[11]:

OVI_width = [9.414924810876614e-05, 7.8696641285401226e-05, 7.791847379015569e-05, 7.2006564252225606e-05, 6.706523541428752e-05, 6.0851059211055571e-05, 4.5281864496212018e-05, 4.1921515759684886e-05, 3.4038495250205607e-05, 3.3388707397500223e-05, 3.1429558301990854e-05]
OVI_SD = [0.0, 2.4743216188748063e-05, 1.6982418742215769e-05, 1.9608000156097438e-05, 1.4919990777198485e-05, 1.3527306543761415e-05, 7.1612811777481039e-06, 1.0282464029715337e-05, 5.8506862234880293e-06, 4.855767316517154e-06, 3.987779881685937e-06]


# In[12]:

plot(radius, OVI_width)
errorbar(radius, OVI_width, yerr = OVI_SD, fmt = 'bo')
title("Oxygen IV Equivalent Width\n Z = 0")
xlabel("Radius (R/Rvir)")
ylabel("Equivalent Width")


# In[ ]:



