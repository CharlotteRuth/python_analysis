import os
import numpy as np
import pynbody
import pynbody.filt as filt
import pynbody.units as units
from pynbody.filt import *
from pynbody import tipsy      #Importing tipsy

##### Writing the line of sight files
# **Make 10 Line of Sight folders and make sure that there is a copy of the .tab file in each**

def writeLOS(los_name,kpcunit, b_arr, nangle):
    x = np.zeros(len(b_arr)*nangle)
    y = np.zeros(len(b_arr)*nangle)
    z = np.zeros(len(b_arr)*nangle)
    axis = np.zeros(len(b_arr)*nangle) + 2
    angles = np.linspace(0,2*np.pi,nangle,endpoint=False)
    
    for i in range(len(b_arr)):
        x[i*nangle: nangle + i*nangle] = np.cos(angles)*b_arr[i]/kpcunit
        y[i*nangle: nangle + i*nangle] = np.sin(angles)*b_arr[i]/kpcunit

    losfile = open(los_name,"w")
    for i in range(len(x)):
        x_str = "{0:0.0f}kpc".format(x[i]*kpcunit)
        y_str = "{0:0.0f}kpc".format(y[i]*kpcunit)
        if x[i] >= 0:
            x_str = '+' + x_str
        if y[i] >= 0:
            y_str = '+' + y_str
        losfile.write("{0:6.5f} {1:6.5f} {2:6.5f} {3:3d} {4:>16}\n".format(x[i],y[i],z[i],int(axis[i]),'x' + x_str + '_y' + y_str))
    
    losfile.close()      


#Define the line of sight parameters
nangle = 12
b_arr = np.linspace(10,120,num = 12) #in kpc

#DEFINE SIMULATION
path = '/home/christenc/Data/Sims/elektra.cosmo25cmb/elektra.cosmo25cmb.4096g5HbwK1BH/elektra.cosmo25cmb.4096g5HbwK1BH.004096/'
tfile = 'elektra.cosmo25cmb.4096g5HbwK1BH.004096'
haloid = 1

os.chdir(path)
#os.mkdir(path + 'CGM')

##### Writing the simulation file
s = pynbody.load(path + tfile)
#s.physical_units()

h = s.halos()     #Loading halos 
halo = h[haloid]     
pynbody.analysis.halo.center(halo,mode='hyb')   #Centering halo
pynbody.analysis.angmom.faceon(halo, cen=(0,0,0))
Rvir = halo.properties['Rvir']/halo.properties['h']*halo.properties['a']
Rvir_sim = np.max(halo['r'].in_units('kpc'))
if not (Rvir_sim < 1.1*Rvir and Rvir_sim > 0.9*Rvir):
    print('Problem with definition of Rvir',Rvir,Rvir_sim)
else:
    print(Rvir)

#Define a filter for a sphere of certain radius (in simulation units)
#sphere = halo[filt.Sphere(Rvir*5 + ' kpc', cen=(0,0,0))] 
sphere_write = s[pynbody.filt.Sphere(str(Rvir*5) + ' kpc', (0,0,0))]

#Defining the filtered simulation as a .TIPSY file 
sphere_write.write(fmt=tipsy.TipsySnap, filename=path + tfile + "." + str(haloid) + ".tipsy")


print("In a terminal, run:")
print("> cd " + path)
#print("> /home/arredond/MAP_Research/Code/smooth/smooth -s 32g -o " + tfile + ".tipsy hsmooth < " + tfile + ".tipsy")
print("> /home/christenc/Code/smooth/smooth -s 32g -o " + tfile + "." + str(haloid) + ".tipsy hsmooth < " + tfile + "." + str(haloid) + ".tipsy")
input("When done, press Enter to continue...")


#Loading the tipsy file 
sphere = pynbody.load(path + tfile + "." + str(haloid) + '.tipsy')

np.copyto(sphere.gas['eps'],sphere.gas['hsm'])
#sphere.gas['eps'] = sphere.gas['hsm']
print(min(sphere.gas['eps']))
print(max(sphere.gas['eps']))

#Writing the tipsy file as a standard file ##MAKE SURE TO NAME THIS .SMOOTH.STD
sphere.write(fmt=tipsy.TipsySnap, filename=path + tfile + "." + str(haloid) + ".tipsy.smooth.std")

#Loading the standard file
sphere2 = pynbody.load(path + tfile + "." + str(haloid) + ".tipsy.smooth.std")

print(min(sphere2.gas['eps']))
print(max(sphere2.gas['eps']))

print("In a terminal, run:")
print("> /home/christenc/Code/tipsy_tools/std2ascii < " + tfile + "." + str(haloid) + ".tipsy.smooth.std > " + tfile + "." + str(haloid) + ".tipsy.smooth.asc")
print("> /home/christenc/Code/tipsy_tools/ascii2bin < " + tfile + "." + str(haloid) + ".tipsy.smooth.asc > " + tfile + "." + str(haloid) + ".tipsy.smooth.bin")
print("> ln -s " + path + tfile + "." + str(haloid) + ".tipsy.smooth.bin CGM/" + tfile + "." + str(haloid) + ".tipsy.smooth.bin")
input("When done, press Enter to continue...")

##### Writing the .aux file
#Define mass fraction of Oxygen, Iron, Carbon, and Silicon within the filter
#Create an array of each mass fraction
Carbon   = 0.00213
Oxygen   = 0.00541
Silicon  = 0.00067
Iron     = 0.00117

oxlist = np.array(sphere_write.gas['OxMassFrac'])
felist = np.array(sphere_write.gas['FeMassFrac'])
clist = np.array(sphere_write.gas['OxMassFrac'])*(Carbon/Oxygen) #sphere_write.gas['c2'])
silist = np.array(sphere_write.gas['FeMassFrac'])*(Silicon/Iron) #sphere_write.gas['si2'])
# The star formation rate is only used in some versions of specexbin
#deltat = max(sphere_write.star['tform'].in_units('yr'))/4096
#sphere_write.gas['tdyn'] = (4.0*np.pi*sphere_write.gas['rho'].in_units('g m^-3')*6.67408*10**(-11))**(-0.5)/3.154e7
#sphere_write.gas['sfr'] = 1.0 - np.exp(-1.0*0.1*deltat/sphere_write.gas['tdyn']*2*sphere_write.gas['H2']/(2*sphere_write.gas['H2']+sphere_write.gas['HI']))
#All gas particles with temperatures more than 10^3 K or densities less than 0.1 amu/cc should never form stars. 
#Therefore, set sphere_write.gas['sfr'] for those particles to zero
#sfrlist = np.array(sphere_write.gas['sfr'])

#Define an auxiliary array containing the mass fraction array of each particle
aux_array = np.array([clist, oxlist, silist, felist]) #, sfrlist])
aux_array_1d = aux_array.flatten(order = 'F')
write_array = array(aux_array_1d,dtype='f')
f = open(path + 'CGM/' + tfile + "." + str(haloid) + ".tipsy.smooth.aux", 'w+b')
write_array.tofile(f)
f.close()


##### Writing the .tab file
#Defining the first entry in our final array ##THIS IS NAMED .TIPSY.SMOOTH
sim_name = path + 'CGM/' + tfile + "." + str(haloid) + ".tipsy.smooth"
#Defining the second entry in our final array 
Redshift = sphere2.properties['z']
#Defining the minimum redshift that we observe for our galaxy 
minRedshift = (min(sphere2['vz'].in_units('km s^-1'))/(300000.0)) + Redshift
#Defining the maximum redshift that we observe for our galaxy 
maxRedshift = (max(sphere2['vz'].in_units('km s^-1'))/(300000.0)) + Redshift
#Listing infile information
output = "{} {} {} {} {}\n"
#Saving the output information as .tab file 
with open(path + 'CGM/' + tfile + "." + str(haloid) + ".tipsy.smooth.tab", 'w') as f:
    f.write(output.format(sim_name,str(Redshift), str(minRedshift), str(maxRedshift), tfile))

los_name = sim_name+".LOSFile"
writeLOS(los_name,s.properties['boxsize'].in_units('kpc a'), b_arr, nangle)

boxsize = s.properties['boxsize'].in_units('Mpc a') #in_units('Mpc') #in_units('Mpc a')
print("Now run the following in your destination folder in terminal:\n")
print("/home/arredond/MAP_Research/specexbin/contspecexbin_v8 {} {} {} {} 1 2\n".format(sim_name,los_name,Redshift,str(boxsize*sphere_write.properties['h'])))
