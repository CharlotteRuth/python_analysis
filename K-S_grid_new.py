#First draft by Lindsey Byrne

import sys
import math
import pylab
import numpy
import pynbody
import matplotlib.pylab as plt

def single_sfr(s): #returns sfr for the past 1 Myr in Msol/yr
    if len(s.s)==0: #if there are no stars in the (sub)sim
        return 0 #then sfr should be 0
    else:
        s.s['mass'].convert_units('Msol') #make sure mass is in Msol
        s.s['age'].convert_units('yr') #make sure time is in yr
        filtr = pynbody.filt.LowPass('age', 50*10**6)
        return s.s[filtr]['mass'].sum()/(50*10**6) #Msol/yr

def cubefilt(sim, x1, y1, x2, y2):
    return (sim[(sim["x"] > x1) & (sim["x"] < x2) & (sim["y"] > y1) & (sim["y"] < y2)])

array = numpy.array

kennicutt = numpy.loadtxt('/home/byrnelin/MAP/kennicutt.txt')
HI = 10**kennicutt[:,0]
H2 = 10**kennicutt[:,2]
gas_sd = numpy.log10(HI+2*H2)
allgas = kennicutt[:,2]
SFR = kennicutt[:,3]

bigiel = numpy.loadtxt('/home/byrnelin/MAP/bigiel.txt')
HI_b = 10**bigiel[:,0]
H2_b = 10**bigiel[:,2]
gas_sd_b = numpy.log10(HI_b+2*H2_b)
SFR_b = bigiel[:,4]

#plt.plot(gas_sd, SFR, 'd', color='black')
plt.plot(gas_sd_b, SFR_b, ',', color='black')
print(len(sys.argv))
for x in numpy.arange(1, len(sys.argv)):
    print(sys.argv[x])
#    s = pynbody.load(sys.argv[x])
#    s.physical_units()
    sim = pynbody.load(sys.argv[x])
    h = sim.halos()
    s = h[1]
    s.physical_units()
    pynbody.analysis.halo.center(s)
    pynbody.analysis.angmom.faceon(s)
    stars_formed = s.s[s.s['age'].in_units('Myr')<50]
    r_max = stars_formed.s['r'].in_units('kpc').max()
 
    print 'r_max is ', r_max
    r_max_pc = s.s['r'].in_units('pc').max()
    z_max = (((s.g['pos'][:,2])**2).max())**0.5 + 0.01 #z-axis confinement, 0.01=safety_extension
    grid_intvl = .75 #in kpc
    num_boxes = int(math.ceil(r_max/grid_intvl)) #of boxes along side of grid b/t center and edge
    print 'number of boxes is ', (2*num_boxes)*(2*num_boxes) #total number of boxes
    area = grid_intvl*grid_intvl
    area_pc = (grid_intvl*1000)*(grid_intvl*1000)
    hor1 = -num_boxes*grid_intvl #first x coordinate
    ver1 = -num_boxes*grid_intvl #first y coordinate
    Gas_SFR = array([])
    xcoords = numpy.arange(hor1, -hor1, grid_intvl)
    def kernel(hcoord):
#        cube1 = pynbody.filt.Cuboid(hcoord,[ver1],[-z_max],[hcoord+grid_intvl],[ver1+grid_intvl],[z_max])
#        c = s[cube1]
        c = cubefilt(s,hcoord,ver1,hcoord+grid_intvl,ver1+grid_intvl)
        SigmaSFR_now = single_sfr(c)/area
        SigmaGas_now = 0.0
        if len(c.g) > 0:
            SigmaGas_now = (c.g['mass']*(c.g['HI']+2*c.g['H2'])).sum()/area_pc
        return array([SigmaGas_now, SigmaSFR_now])
    for i in xrange(0, 2*num_boxes):
        Gas_SFR_new = array(map(kernel, xcoords))
        Gas_SFR = numpy.append(Gas_SFR, Gas_SFR_new)
        print ver1
        ver1 += grid_intvl
    SigmaGas = numpy.log10(Gas_SFR[0::2])
    SigmaSFR = numpy.log10(Gas_SFR[1::2])
#    print SigmaGas
#    print SigmaSFR
#    print Gas_SFR
#    print len(Gas_SFR)
#    plt.plot(SigmaGas, SigmaSFR,'x', color='red',lw=2)
    plt.plot(SigmaGas, SigmaSFR,'x')
#plt.plot(gas_sd, SFR, 'd', color='black')
#plt.plot(gas_sd_b, SFR_b, ',', color='black')
plt.xlabel('$\Sigma_{HI+H2}$ [M$_\odot$/pc^2]')
plt.ylabel('$\Sigma_{SFR}$ [M$_\odot$/(yr*kpc^2)]')
plt.ylim(-5,0)
plt.xlim(-3,4)
plt.title('ks_cosmo_shield_750pc')
plt.savefig("ks_cosmo_shield_750pc_2.png")
