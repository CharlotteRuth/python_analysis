import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib

ion = 'hi'
iion = 0

ion = 'civ'
iion = 3

ion = 'ovi'
iion = 5

nhlow = -9.0
nhpts = 240
deltanh = 0.05
tlow = 2.5
tpts = 140
deltat = 0.05
nion = 9


filename = "/home/christenc/Code/specexbin/ionfiles/lt00HM12_i9"
data = []
f = open(filename,'r')
for line in f.readlines():
    y = [float(value) for value in line.split()]
    data.append( y )
f.close()
data = np.array(data)

assert (shape(data)[0] == nhpts*tpts)
assert (shape(data)[1] == nion)
data2 = data.reshape(nhpts,tpts,nion) #density, temp, ions

dens_vals_cgmspec = np.linspace(nhlow,nhlow + deltanh*(nhpts - 1),num = nhpts)
temp_vals_cgmspec = np.linspace(tlow,tlow + deltat*(tpts - 1),num = tpts)

ifs = np.load("/home/christenc/Code/python/pynbody/lib/python2.7/site-packages/pynbody-0.43-py2.7-linux-x86_64.egg/pynbody/analysis/ionfracs.npz") 
redshift_vals = ifs['redshiftvals'].view(np.ndarray)
temp_vals = ifs['tempvals'].view(np.ndarray)
dens_vals = ifs['denvals'].view(np.ndarray)
vals = ifs[ion + 'if'].view(np.ndarray) #redshift, temp, density

dens = 0
temp = 7.5
nhind_cgmspec = (np.abs(dens_vals_cgmspec - dens)).argmin()
nhind_pynbody = (np.abs(dens_vals - dens)).argmin()
tind_cgmspec = (np.abs(temp_vals_cgmspec - temp)).argmin()
tind_pynbody = (np.abs(temp_vals - temp)).argmin()
print(vals[0,tind_pynbody,nhind_pynbody])
print(data2[nhind_cgmspec,tind_cgmspec,0])

plt.close()
plt.figure(1)
plt.imshow(vals[0,:,:],vmin = -9, vmax = 0,origin = 'lower',extent = (dens_vals[0],dens_vals[-1],temp_vals[0],temp_vals[-1]))
plt.show()
#plt.close()

plt.figure(2)
plt.imshow(transpose(data2[:,:,iion]),vmin = -9, vmax = 0,origin = 'lower',extent = (dens_vals_cgmspec[0],dens_vals_cgmspec[-1],temp_vals_cgmspec[0],temp_vals_cgmspec[-1]))
plt.axis([dens_vals[0],dens_vals[-1],temp_vals[0],temp_vals[-1]])
plt.show()
