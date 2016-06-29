#main program to call various functions

import pynbody
import matplotlib.pylab as plt
import matplotlib.mlab as mlab
import matplotlib.cm as cm
import numpy as np
import math
from scipy import interpolate
from mpl_toolkits.mplot3d import Axes3D
import v_functions as w #most of the custom functions used are in this file
import sys, os, glob, pynbody.bridge
import pynbody.snapshot.tipsy

#Establishes galaxy to be used for the rest of simulation. Re-run the code to
#use a different galaxy

sim_number = raw_input('Galaxy Number: ')

pathbase = '/home/christensen/Storage1/UW/MolecH/Cosmo/'

path = {'239':pathbase + 'h239.cosmo50cmb.3072g/',
               '258':pathbase + 'h258.cosmo50cmb.3072g/',
               '285':pathbase + 'h285.cosmo50cmb.3072g/',
               '516':pathbase + 'h516.cosmo25cmb.3072g/',
               '516m':pathbase + 'h516.cosmo25cmb.3072g/',
               '603':pathbase + 'h603.cosmo50cmb.3072g/',
               '799':pathbase + 'h799.cosmo25cmb.3072g/',
               '986':pathbase + 'h986.cosmo50cmb.3072g/'}


simulations = {'239':path[sim_number] + 'h239.cosmo50cmb.3072g14HMbwK/h239.cosmo50cmb.3072g14HMbwK.00512/h239.cosmo50cmb.3072g14HMbwK.00512',
               '258':path[sim_number] + 'h258.cosmo50cmb.3072g14HMbwK/h258.cosmo50cmb.3072g14HMbwK.00512/h258.cosmo50cmb.3072g14HMbwK.00512',
               '285':path[sim_number] + 'h285.cosmo50cmb.3072g14HMbwK/h285.cosmo50cmb.3072g14HMbwK.00512/h285.cosmo50cmb.3072g14HMbwK.00512',
               '516':path[sim_number] + 'h516.cosmo25cmb.3072g14HBWK_2/h516.cosmo25cmb.3072g14HBWK.00512/h516.cosmo25cmb.3072g14HBWK.00512',
               '516m':path[sim_number] + 'h516.cosmo25cmb.3072g1MBWK/h516.cosmo25cmb.3072g1MBWK.00512/h516.cosmo25cmb.3072g1MBWK.00512',
               '603':path[sim_number] + 'h603.cosmo50cmb.3072g14HBWK/h603.cosmo50cmb.3072g14HBWK.00512/h603.cosmo50cmb.3072g14HBWK.00512',
               '799':path[sim_number] + 'h799.cosmo25cmb.3072g14HBWK/h799.cosmo25cmb.3072g14HBWK.00512/h799.cosmo25cmb.3072g14HBWK.00512',
               '986':path[sim_number] + 'h986.cosmo50cmb.3072g14HBWK/h986.cosmo50cmb.3072g14HBWK.00512/h986.cosmo50cmb.3072g14HBWK.00512'}

starlogs = {'239':path[sim_number] + 'h239.cosmo50cmb.3072g14HMbwK/h239.cosmo50cmb.3072g14HMbwK.starlog',
            '258':path[sim_number] + 'h258.cosmo50cmb.3072g14HMbwK/h258.cosmo50cmb.3072g14HMbwK.starlog',
            '285':path[sim_number] + 'h285.cosmo50cmb.3072g14HMbwK/h285.cosmo50cmb.3072g14HMbwK.starlog.fits',
            '516':path[sim_number] + 'h516.cosmo25cmb.3072g14HBWK_2/h516.cosmo25cmb.3072g14HBWK.starlog.fits',
            '516m':path[sim_number] + 'h516.cosmo25cmb.3072g1MBWK/h516.cosmo25cmb.3072g1MBWK.starlog',
            '603':path[sim_number] + 'h603.cosmo50cmb.3072g14HBWK/h603.cosmo50cmb.3072g14HBWK.starlog.fits',
            '799':path[sim_number] + 'h799.cosmo25cmb.3072g14HBWK/h799.cosmo25cmb.3072g14HBWK_2merge.starlog.fits',
            '986':path[sim_number] + 'h986.cosmo50cmb.3072g14HBWK/h986.cosmo50cmb.3072g14HBWK.starlog.fits'}
    
#saved alignment data
alignments = {'239':path[sim_number] + 'h239.cosmo50cmb.3072g14HMbwK/grp1.alignment.txt',
              '258':path[sim_number] + 'h258.cosmo50cmb.3072g14HMbwK/grp1.alignment.txt',
              '285':path[sim_number] + 'h285.cosmo50cmb.3072g14HMbwK/grp1.alignment.txt',
              '516':path[sim_number] + 'h516.cosmo25cmb.3072g14HBWK/grp1.alignment.txt',
              '516m':path[sim_number] + 'h516.cosmo25cmb.3072g1MBWK/grp1.alignment.txt',
              '603':path[sim_number] + 'h603.cosmo50cmb.3072g14HBWK/grp1.alignment.txt',
              '799':path[sim_number] + 'h799.cosmo25cmb.3072g14HBWK/grp1.alignment.vdisp.txt',
              '986':path[sim_number] + 'h986.cosmo50cmb.3072g14HBWK/grp1.alignment.txt'}

#where to save alignment data
save_data = {'239':path[sim_number] + 'h239.cosmo50cmb.3072g14HMbwK/',
              '258':path[sim_number] + 'h258.cosmo50cmb.3072g14HMbwK/',
              '285':path[sim_number] + 'h285.cosmo50cmb.3072g14HMbwK/',
              '516':path[sim_number] + 'h516.cosmo25cmb.3072g14HBWK/',
              '516m':path[sim_number] + 'h516.cosmo25cmb.3072g1MBWK/',
              '603':path[sim_number] + 'h603.cosmo50cmb.3072g14HBWK/',
              '799':path[sim_number] + 'h799.cosmo25cmb.3072g14HBWK/',
              '986':path[sim_number] + 'h986.cosmo50cmb.3072g14HBWK/'}

#information on mergers
mergersx = {'239':[13.732138263148142-0.141262*3.97e+01,13.732138263148142-0.141262*3.97e+01], 
           '258':[13.732138263148142-0.132968*3.97e+01,13.732138263148142-0.132968*3.97e+01], 
           '285':[13.732138263148142-0.0721415*3.97e+01,13.732138263148142-0.0721415*3.97e+01],
           '516':[13.731484849641138-0.155087*3.88e+01,13.731484849641138-0.155087*3.88e+01],
           '799':[13.730648257851353-0.0859657*3.88e+01,13.730648257851353-0.0859657*3.88e+01],
           '986':[13.732138263148142-0.0859657*3.97e+01,13.732138263148142-0.0859657*3.97e+01,13.732138263148142-0.168911*3.97e+01,13.732138263148142-0.168911*3.97e+01]}

mergersy = {'239':[1,200], 
           '258':[1,200], 
           '285':[1,200],
           '516':[1,200],
           '799':[1,200],
           '986':[1,200,200,1]}

mergersy2 = {'239':[1,140], 
             '258':[1,140], 
             '285':[1,140],
             '516':[1,30],
             '799':[1,25],
             '986':[1,50,50,1]}


s = pynbody.load(simulations[sim_number])
h = s.halos()

if (sim_number == '258') or (sim_number == '239'):
    x = os.path.abspath(save_data[sim_number])
#    x = os.path.dirname(x)
    l = glob.glob(os.path.join(x,"*.starlog"))
    if (len(l)) :
        for filename in l :
            sl = pynbody.tipsy.StarLog(filename)
else:
    sl = w.fits_starlog(starlogs[sim_number], sim_number)

#selects appropriate stars and aligns galaxy
bin_size = .25
age_min = 0
age_max = 14

h1 = h[1].s[(h[1].s['tform'].in_units('Gyr') >= 2.1891415550120792)]

al_cx,al_cy,al_cz,al_vcx,al_vcy,al_vcz, al_am = w.center_values(h1, 'mass')

t1, t2 = w.align_sim(al_am[0], al_am[1], al_am[2])

h1['x'] = h1['x']-al_cx
h1['y'] = h1['y']-al_cy
h1['z'] = h1['z']-al_cz

h1['vx'] = h1['vx']-al_vcx
h1['vy'] = h1['vy']-al_vcy
h1['vz'] = h1['vz']-al_vcz

h1['pos'] = np.dot(np.dot(h1['pos'],t1),t2)
h1['vel'] = np.dot(np.dot(h1['vel'],t1),t2)

#access starlog
c = np.in1d(sl['iord'],h1['iord'])
h1sl = sl[c]
age = h1['tform'].in_units('Gyr').max() - h1['tform'].in_units('Gyr')

posx_list,posy_list,posz_list,indices = w.align_spline(h1sl, alignments[sim_number])

satellites = {'239':[0.002,0.006,0.002], 
              '258':[0.008,0.008,0.008], 
              '285':[0.0025,0.0025,0.0025],
              '516':[0.00015,0.00025,0.0002],
              '516m':[1,1,1],
              '603':[1,1,1],
              '799':[0.0001,0.0001,0.0001],
              '986':[0.0004,0.0003,0.0004]}

satellites2 = {'239':[0.003], 
              '258':[0.002], 
              '285':[0.0025],
              '516':[0.00025],
              '516m':[1],
              '603':[1],
              '799':[0.0001],
              '986':[0.0008]}

index1 = np.sqrt(np.absolute(h1sl['x']-posx_list[indices])**2 + np.absolute(h1sl['y']-posy_list[indices])**2 + np.absolute(h1sl['z']-posz_list[indices])**2) <= satellites2[sim_number]

h1sl = h1sl[index1]

#uncomment this to align stars
h1sl['pos'],h1sl['vel'] = w.align_stars3(h1sl,.1, alignments[sim_number])

h1['vel'].convert_units('km s**-1')
h1sl['vel'].convert_units('km s**-1 aform')

#uncomment this to save new alignment data

np.save(save_data[sim_number]+'vel', h1sl['vel'])
np.save(save_data[sim_number]+'pos', h1sl['pos'])

h1sl['vel'] = np.load(save_data[sim_number]+'vel.npy')
h1sl['pos'] = np.load(save_data[sim_number]+'pos.npy')

h1slbk = h1sl
age2bk = h1slbk['tform'].in_units('Gyr').max() - h1slbk['tform'].in_units('Gyr')
h1slr = (h1slbk['x'].in_units('kpc aform')**2 + h1slbk['y'].in_units('kpc aform')**2 + h1slbk['z'].in_units('kpc aform')**2)**0.5*h1slbk['x'].sim['aform']
h1slr2 = (h1slbk['x'].in_units('kpc aform')**2 + h1slbk['y'].in_units('kpc aform')**2)**0.5*h1slbk['x'].sim['aform']
h1slv2 = (h1slbk['vx'].in_units('km s**-1 aform')**2 + h1slbk['vy'].in_units('km s**-1 aform')**2)**0.5*h1slbk['x'].sim['aform']
#plot(age2bk,h1slr,'.b')
plot(age2bk,h1slr2,'.r')
plot(age2bk,h1slv2,'.r')
plot(age2bk,np.absolute(h1slbk['z'].in_units('kpc aform')*h1slbk['x'].sim['aform']),'b.')
plot(age2bk,h1slbk['vz'].in_units('km s**-1 aform')*h1slbk['x'].sim['aform'],'b.')

index2 = np.absolute(h1sl['z'].in_units('kpc aform')*h1sl['x'].sim['aform']) <= 0.2
h1sl = h1sl[index2]
age2 = h1sl['tform'].in_units('Gyr').max() - h1sl['tform'].in_units('Gyr')
plot(age2,np.absolute(h1sl['z'].in_units('kpc aform')*h1sl['x'].sim['aform']),'.g')
plot(age2,h1sl['vz'].in_units('km s**-1 aform')*h1sl['x'].sim['aform'],'.g')
plot(age2,sqrt(h1sl['x'].in_units('kpc aform')**2 + h1sl['y'].in_units('kpc aform')**2)*h1sl['x'].sim['aform'],'.y')

#this section calculates velocity dispersion for various objects using functions

#ISM dispersion
sigma, years = w.ISM(sim_number, s)

#z=0
v_disp, v_age = w.dispersion_data(bin_size, age_min, age_max, h1, age)

#at formation
v_disp2, v_age2 = w.dispersion_data(bin_size, age_min, age_max, h1sl, age2)

#stars formed ex situ
hsat = h1[(np.in1d(h1['iord'],np.setdiff1d(h1['iord'], h1sl['iord'])))]
age3 = h1['tform'].in_units('Gyr').max() - hsat['tform'].in_units('Gyr')
v_disp3, v_age3 = w.dispersion_data(bin_size, age_min, age_max, hsat, age3)

#sigma2, years2 = w.ISM2(sim_number, h[1].g, v_disp2, v_age2)
sigma2, years2 = w.gas_part_grav_sof(sim_number, v_disp2, v_age2)

age_points = []
disp_points = []
disp2_points = []

#major mergers
x = mergersx[sim_number]
y = mergersy[sim_number]
y2 = mergersy2[sim_number]

obs_age = [0,2,5,11.5]
obs_disp = [5.25,10,14.5,17.5]
xerr = [.5,2,1.67,2.25]
yerr = [4,3,3,5.5]

#uncomment this section to plot galaxies
#w.plot_tform(h1sl, 0, 14, 'x','y', 'kpc aform')
#w.plot_tform(h1, 0, 14, 'x', 'y', 'kpc')
#w.plot_tform(hsat, 0, 14, 'x', 'z', 'kpc')
#w.plot_tform(h1sl, 0, 14, 'x','z','kpc aform')
#w.plot_tform(h1, 0, 14, 'x','z', 'kpc')
#w.plot_tform_rect(h1, 0, 14, 'x','z', 'kpc')
#w.plot_tform_rect(h1sl, 0, 14, 'x','z','kpc aform')
#w.velocity_hist(7,8, h1, h1sl, 'vz')

#velocity dispersion plot
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(v_age, v_disp, 'k', label = ('z = 0'), linewidth = 2)
ax.plot(v_age2, v_disp2, 'k', label = ('At Formation'), linewidth = 1)
#ax.plot(v_age3, v_disp3, 'k--', label = ('Stars formed Ex Situ, z = 0'), linewidth = 1)
#ax.plot(years, sigma, 'b', label = ('ISM'))
ax.plot(years2, sigma2, 'green', label = ('ISM + Dynamical Heating'))
ax.plot(x, y, 'r', label = ('Major Merger'))
#ax.set_yscale('log')
ax.grid(b=True, which='major')
ax.grid(b=True, which='minor')
ax.set_xlim(0,13)
ax.set_ylim(0,30)
plt.xlabel('Age [Gyr]')
plt.ylabel('Velocity Dispersion (km/s)')
plt.legend(loc=2, prop={'size':10})

w.data_table('dispersion_data_table'+sim_number+'.txt', sim_number, v_age, v_disp2, v_disp)

#linear plot
#fig = plt.figure()
#ax2 = fig.add_subplot(1,1,1)
#ax2.plot(v_age, v_disp, 'k', label = ('z = 0'), linewidth = 2)
#ax2.plot(v_age2, v_disp2, 'k', label = ('At Formation'), linewidth = 1)
#ax2.plot(v_age3, v_disp3, 'k--', label = ('Stars formed Ex Situ, z = 0'), linewidth = 1)
#ax2.plot(years, sigma, 'b', label = ('ISM'))
#ax2.plot(years2, sigma2, 'green', label = ('ISM + Dynamical Heating'))
#ax2.errorbar(obs_age, obs_disp,xerr=xerr, yerr=yerr, fmt='o', label=('WLM observation'))
#ax2.plot(x, y, 'r', label = ('Major Merger'))
#ax2.grid(b=True, which='major')
#ax2.grid(b=True, which='minor')
#ax2.set_xlim(0,13)
#ax2.set_ylim(0,30)
#plt.xlabel('Age [Gyr]')
#plt.ylabel('Velocity Dispersion (km/s)')
#plt.legend(loc=2, prop={'size':10})

#plot heating vs cooling (I don't think this ever actually worked out)
## fig = plt.figure()
## ax3 = fig.add_subplot(1,1,1)
## rat = [10.5900831496/4.09380998495,11.7682659385/7.05834632684,22.7632034449/13.3361145074,44.0133264522/42.0437522664,39.5531516632/36.5477203063,43.7119681057/41.6717635497]
## massrat = [2.4e10,3.8e10,1.9e11,6.8e11,7.7e11,8.8e11]
## ax3.plot(massrat, rat, 'k')
## ax.set_xscale('log')
## ax3.grid(b=True, which='major')
## ax3.grid(b=True, which='minor')
## plt.ylabel('Sigma_heat/Sigma_cool')
## plt.xlabel('log 10 Solar Masses')

## print 'Sigma heating: '
## print np.mean(sigma2[~np.isnan(sigma2)])
## print 'Sigma cooling: '
## v_disp2 = np.array(v_disp2)
## print np.mean(v_disp2[~np.isnan(v_disp2)])
pause = raw_input('Continue? ')
