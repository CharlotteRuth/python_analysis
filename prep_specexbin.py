#Charlotte Christensen
#2/27/18
#Prep a tipsy file to run specexbin on it

#Run with
#%run /home/christensen/Code/python/python_analysis/prep_specexbin.py
#or, on quirm,
#%run /home/christenc/Code/python/python_analysis/prep_specexbin.py
#ipython --pylab

import os
import numpy as np
import pynbody
import pynbody.filt as filt
import pynbody.units as units
from pynbody.filt import *
from pynbody import tipsy      #Importing tipsy
import socket

##### Writing the line of sight files
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

def prep_specexbin_files(path,tfile,haloid):
    os.chdir(path)
    if not os.path.isdir(path + 'CGM.' + str(haloid)):
        os.mkdir(path + 'CGM.' + str(haloid))
        
    #Define the line of sight parameters
    nangle = 12
    b_arr = np.linspace(10,120,num = 12) #in kpc

    ##### Writing the simulation file
    s = pynbody.load(path + tfile)
    #s.physical_units()

    h = s.halos()     #Loading halos 
    halo = h[haloid]     
    pynbody.analysis.halo.center(halo,mode='hyb')   #Centering halo
    pynbody.analysis.angmom.faceon(halo, cen=(0,0,0))
    s.rotate_x(45)
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
    print("cd " + path)
    #print("> /home/arredond/MAP_Research/Code/smooth/smooth -s 32g -o " + tfile + ".tipsy hsmooth < " + tfile + ".tipsy")
    print("/home/christenc/Code/smooth/smooth -s 32g -o " + tfile + "." + str(haloid) + ".tipsy hsmooth < " + tfile + "." + str(haloid) + ".tipsy")
    wait = raw_input("When done, press Enter to continue...")


    #Loading the tipsy file 
    sphere = pynbody.load(path + tfile + "." + str(haloid) + '.tipsy')
    
    np.copyto(sphere.gas['eps'],sphere.gas['hsm'])
    #sphere.gas['eps'] = sphere.gas['hsm']
    print("Minimum softening: " + str(min(sphere.gas['eps'])))
    print("Maximum softening: " + str(max(sphere.gas['eps'])))
    
    #Writing the tipsy file as a standard file ##MAKE SURE TO NAME THIS .SMOOTH.STD
    sphere.write(fmt=tipsy.TipsySnap, filename=path + tfile + "." + str(haloid) + ".tipsy.smooth.std")
    
    #Loading the standard file
    sphere2 = pynbody.load(path + tfile + "." + str(haloid) + ".tipsy.smooth.std")
    
    print("Minimum softening (should equal above): " + str(min(sphere2.gas['eps'])))
    print("Maximum softening (should equal above): " + str(max(sphere2.gas['eps'])))
    
    print("In a terminal, run:")
    print("rm CGM." + str(haloid) + "/" + tfile + "." + str(haloid) + ".tipsy.smooth.bin CGM." + str(haloid) + "/" + tfile + "." + str(haloid) + ".tipsy.smooth "  + tfile + "." + str(haloid) + ".tipsy.smooth.bin " + tfile + "." + str(haloid) + ".tipsy.smooth.asc")
    print("/home/christenc/Code/tipsy_tools/std2ascii < " + tfile + "." + str(haloid) + ".tipsy.smooth.std > " + tfile + "." + str(haloid) + ".tipsy.smooth.asc")
    print("/home/christenc/Code/tipsy_tools/ascii2bin < " + tfile + "." + str(haloid) + ".tipsy.smooth.asc > " + tfile + "." + str(haloid) + ".tipsy.smooth.bin")
    print("ln -s " + path + tfile + "." + str(haloid) + ".tipsy.smooth.bin CGM." + str(haloid) + "/" + tfile + "." + str(haloid) + ".tipsy.smooth.bin")
    print("ln -s " + path + tfile + "." + str(haloid) + ".tipsy.smooth.bin CGM." + str(haloid) + "/" + tfile + "." + str(haloid) + ".tipsy.smooth")
    print("cd CGM." + str(haloid))
    wait2 = raw_input("When done, press Enter to continue...")
    
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
    write_array = np.array(aux_array_1d,dtype='f')
    f = open(path + 'CGM.' + str(haloid)+'/' + tfile + "." + str(haloid) + ".tipsy.smooth.aux", 'w+b')
    write_array.tofile(f)
    f.close()
    
    
    ##### Writing the .tab file
    #Defining the first entry in our final array ##THIS IS NAMED .TIPSY.SMOOTH
    sim_name = path + 'CGM.'+ str(haloid)+'/' + tfile + "." + str(haloid) + ".tipsy.smooth"
    #Defining the second entry in our final array 
    #Redshift = sphere2.properties['z']
    Redshift = halo.properties['z']
    #Defining the minimum redshift that we observe for our galaxy 
    minRedshift = (min(halo['vz'].in_units('km s^-1'))/(300000.0)) + Redshift
    #minRedshift = (min(sphere2['vz'].in_units('km s^-1'))/(300000.0)) + Redshift
    #Defining the maximum redshift that we observe for our galaxy 
    #maxRedshift = (max(sphere2['vz'].in_units('km s^-1'))/(300000.0)) + Redshift
    maxRedshift = (max(halo['vz'].in_units('km s^-1'))/(300000.0)) + Redshift
        #Listing infile information
    output = "{} {} {} {} {}\n"
    #Saving the output information as .tab file 
    with open(path + 'CGM.'+str(haloid)+'/' + tfile + "." + str(haloid) + ".tipsy.smooth.tab", 'w') as f:
        f.write(output.format(sim_name,str(Redshift), str(minRedshift), str(maxRedshift), tfile))

    f.close()
    los_name = sim_name+".LOSFile"
    writeLOS(los_name,s.properties['boxsize'].in_units('kpc a'), b_arr, nangle)

    boxsize = s.properties['boxsize'].in_units('Mpc a') #in_units('Mpc') #in_units('Mpc a')
    print("Now run the following in your destination folder in terminal:\n")
    print("/home/christenc/Code/specexbin/contspecexbin_v8 {} {} {} {} 1 2\n".format(sim_name,los_name,Redshift,str(boxsize*halo.properties['h'])))
    f = open(path + 'CGM.'+str(haloid)+'/' + tfile + "." + str(haloid) + '.CGM.sh',"w")
    f.write("/home/christenc/Code/specexbin/contspecexbin_v8 {} {} {} {} 1 2\n".format(sim_name,los_name,Redshift,str(boxsize*halo.properties['h'])))
    f.close()

if __name__ == '__main__':
    if (socket.gethostname() == "quirm"):
        prefix = '/home/christenc/Data/Sims/'
    else:
        prefix = '/home/christensen/Storage2/UW/MolecH/Cosmo/'
        
    #DEFINE SIMULATION
    path = prefix + 'cptmarvel.cosmo25cmb/cptmarvel.cosmo25cmb.4096g5HbwK1BH/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096/'
    tfile = 'cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096'
    haloid = 1
    #prep_specexbin_files(path,tfile,haloid)
    
    path = prefix + 'elektra.cosmo25cmb/elektra.cosmo25cmb.4096g5HbwK1BH/elektra.cosmo25cmb.4096g5HbwK1BH.004096/'
    tfile = 'elektra.cosmo25cmb.4096g5HbwK1BH.004096'
    haloid = 1
    #prep_specexbin_files(path,tfile,haloid)
    haloid = 2
    #prep_specexbin_files(path,tfile,haloid)
    
    path = prefix + "rogue.cosmo25cmb/rogue.cosmo25cmb.4096g5HbwK1BH/rogue.cosmo25cmb.4096g5HbwK1BH.004096/"
    tfile = "rogue.cosmo25cmb.4096g5HbwK1BH.004096"
    haloid = 1
    #prep_specexbin_files(path,tfile,haloid)
    haloid = 3
    #prep_specexbin_files(path,tfile,haloid)
    
    path = prefix + 'storm.cosmo25cmb/storm.cosmo25cmb.4096g5HbwK1BH/storm.cosmo25cmb.4096g5HbwK1BH.004096/'
    tfile = 'storm.cosmo25cmb.4096g5HbwK1BH.004096'
    haloid = 1
    #prep_specexbin_files(path,tfile,haloid)
    haloid = 2    
    prep_specexbin_files(path,tfile,haloid)
