import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors as colors
import socket
import pandas as pd
import sys, os, glob, pickle


#Begum+, 2008 Table 1
f = open(dataprefix+'FiggsTable1.txt','r')
begumdata_T1 = []
for line in f:
    line = line.strip()
    columns = line.split()
    if len(columns) == 15 and columns[0] != 'Galaxy':
        source = {}
        source['gal'] = columns[0]
        source['RA_1'] = float(columns[1]) # hour
        source['RA_2'] = float(columns[2]) # min
        source['RA_3'] = float(columns[3]) # s
        source['dec_1'] = float(columns[4]) # degree
        source['dec_2'] = float(columns[5]) # arcmin
        source['dec_3'] = float(columns[6]) # arcsec
        source['M_B'] = float(columns[7]) # mag
        source['D_H0'] = float(columns[8]) # arcmin
        source['B_V'] = float(columns[9]) # mag == -1 if unknown
        source['d'] = float(columns[10]) # Mpc
        source['D_estm'] = columns[11]
        source['group'] = columns[12]
        source['i_opt'] = float(columns[13]) # degree
        source['ref'] = columns[14]
        begumdata_T1.append(source)
f.close()
begumdata_T1 = pd.DataFrame(begumdata_T1)

#Begum+, 2008 Table 3
f = open(dataprefix+'FiggsTable3.txt','r')
#Conversion between Mag_B and stellar mass from McGaugh+ 2008
M_B_sun = 5.44 #Mann & von Braun, 2015 PASP 127,102
MLa = -0.942
MLb = 1.69
begumdata = []
for line in f:
    line = line.strip()
    columns = line.split()
    if len(columns) == 17 and columns[0] != 'Galaxy':
        source = {}
        source['gal'] = columns[0]
        source['FI_GMRT'] = float(columns[1]) #Jy km s^-1
        source['FI_GMRT_err_sign'] = columns[2]
        source['FI_GMRT_err'] = float(columns[3]) #Jy km s^-1
        source['V_sys'] = float(columns[4]) #km s^-1
        source['DeltaV_50'] = float(columns[5]) #km s^-1
        source['D_HI'] = float(columns[6]) #arcmin
        source['M_HI'] = float(columns[7]) #10^6 Msun
        source['M_HI_L_B'] = float(columns[8])
        source['D_HI_D_H0'] = float(columns[9])
        source['FI_GMRT_FI_SD'] = float(columns[10])#-1 means data was not available
        source['FI_GMRT_FI_SD_err_sign'] = columns[11]
        source['FI_GMRT_FI_SD_err'] = float(columns[12])
        source['i_HI'] = float(columns[13]) #-1 means data was not available
        source['i_HI_err_sign'] = columns[14]
        source['i_HI_err'] = float(columns[15])
        source['ref'] = columns[16]
        if any(begumdata_T1['gal'] == columns[0]):
            source['M_B'] = float(begumdata_T1['M_B'][begumdata_T1['gal'] == columns[0]]) # mag
            source['B_V'] = float(begumdata_T1['B_V'][begumdata_T1['gal'] == columns[0]]) # mag
            if source['B_V'] != -1:
                source['L_B'] = 10**((source['M_B'] - M_B_sun)/-2.5)
                MLfit = 10**(MLa + MLb*source['B_V'])
                source['smass'] = MLfit*source['L_B']
            else:
                source['L_B'] = -1
                source['smass'] = -1                
        else:
            source['M_B'] = -1
            source['B_V'] = -1
            source['L_B'] = -1
            source['smass'] = -1
        begumdata.append(source)
f.close()
begumdata = pd.DataFrame(begumdata)

#Plot HI gas fraction (M_HI/M_*) versus stellar mass from the FIGGS dwarf galaxy sample
plt.figure(1)
begum_plt = plt.scatter(begumdata['smass'],begumdata['M_HI']*1e6/begumdata['smass'],s = markersize,c = "grey", marker = "D")

