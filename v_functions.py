#These are all of the functions that are potentially called by v_sim2.py

import pynbody
import matplotlib.pylab as plt
import matplotlib.mlab as mlab
import matplotlib.cm as cm
import numpy as np
import math
from scipy import interpolate
from mpl_toolkits.mplot3d import Axes3D

#********************************************************************************
#aligns simulation
def align_sim(amvx, amvy, amvz):
    angmom = [amvx, amvy, amvz]

    phi = math.atan(angmom[2]/angmom[0])
    t1 = -np.array([(-math.sin(phi),0,math.cos(phi)),
                   (0,1,0),
                   (math.cos(phi),0,math.sin(phi))])
    angmom2 = np.dot(angmom,t1)

    theta = math.atan(angmom2[2]/angmom2[1])
    t2 = np.array([(1,0,0),
                   (0,-math.sin(theta),math.cos(theta)),
                   (0,math.cos(theta),math.sin(theta))])
    angmom3 = np.dot(angmom2,t2)
    return t1, t2

#aligns simulation
def align_sim_list(angmom):

    phi = math.atan(angmom[2]/angmom[0])
    t1 = -np.array([(-math.sin(phi),0,math.cos(phi)),
                   (0,1,0),
                   (math.cos(phi),0,math.sin(phi))])
    angmom2 = np.dot(angmom,t1)

    theta = math.atan(angmom2[2]/angmom2[1])
    t2 = np.array([(1,0,0),
                   (0,-math.sin(theta),math.cos(theta)),
                   (0,math.cos(theta),math.sin(theta))])
    angmom3 = np.dot(angmom2,t2)
    return t1, t2
    

#returns velocity dispersion
def dispersion_data(bin_size, age_min, age_max, sim, age):
    v_disp = []
    v_age = []
    for i in range (0, int((age_max-age_min)/bin_size)):
        age_bin = sim[(age >= age_min+(bin_size*i)) &
                      (age < age_min+(bin_size*(i+1)))]
        #v_age.append(i*bin_size)
        v_age.append(float(age_max - np.sum(age_bin['tform'].in_units('Gyr'))/len(age_bin)))
        v_disp.append(float(np.std(age_bin['vz'])))
    return v_disp, v_age

#I made this as a more intuitive version of python's spline function
def spline(time_array, var, time):
    rep = interpolate.splrep(time_array, var, s=0)
    return interpolate.splev(time, rep, der=0)

def transform(type, type1, type2, type3, sim, t):
    h1slx = np.sum(sim[type]*t[::3][indices], axis=1)
    h1sly = np.sum(sim[type]*t[1::3][indices], axis=1)
    h1slz = np.sum(sim[type]*t[2::3][indices], axis=1)
    sim[type1] = h1slx
    sim[type2] = h1sly
    sim[type3] = h1slz

#plots a histogram of velocities
def velocity_hist(age_min, age_max, sim1, sim2, vel_type):
    age_bin1 = sim1[(sim1['tform'].in_units('Gyr').max()-sim1['tform'].in_units('Gyr') >= age_min) & (sim1['tform'].in_units('Gyr').max()-sim1['tform'].in_units('Gyr') < age_max)]
    #age_bin2 = sim2[(sim2['tform'].in_units('Gyr').max()-sim2['tform'].in_units('Gyr') >= age_min) & (sim2['tform'].in_units('Gyr').max()-sim2['tform'].in_units('Gyr') < age_max)]
    
    d1 = np.std(age_bin1['vz'].in_units('km s**-1'))
    #d2 = np.std(age_bin2['vz'].in_units('km s**-1 aform'))

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.hist(age_bin1[vel_type].in_units('km s**-1'), bins = 50, color = 'white', label = ('z = 0: %f'%d1))
    #ax.hist(age_bin2[vel_type].in_units('km s**-1 aform'), bins = 50,  color = 'grey', label = ('At time of formation: %f'%d2), alpha = 0.5)
    plt.xlabel('Velocity ['+vel_type+'] (km/s)')
    plt.ylabel('Number of stars')
    plt.title('Stars in range %d to %d Gyrs'%(age_min, age_max))
    #plt.suptitle(sim_number)
    plt.legend(loc = 2, prop={'size':10})

#plots a galaxy with stars in an age range, colors based on age
#sim: which simulation should be used
#age_mi: minimum age of star particles
#age_ma: maximum age of star particles
#d1: dimension 1 ('x','y','z')
#d2: dimension 2 ('x','y','z')
#units: units of distance ('kpc' or 'kpc aform' (for starlogs))
def plot_tform(sim, age_mi, age_ma, d1, d2, units):
    age_bin = sim[(sim['tform'].in_units('Gyr') >= age_mi) &
                  (sim['tform'].in_units('Gyr') < age_ma)]
    fig = plt.figure()
    fig.add_subplot(1,1,1,adjustable='box', aspect=0.3)
    plt.subplot('111', axisbg = (0,0,.05))
    #plt.scatter(age_bin[d1].in_units(units),age_bin[d2].in_units(units), c=abs(age_bin['tform'].in_units('Gyr')),marker = '.',cmap=cm.get_cmap('drew'), edgecolors= 'none', alpha = 0.5)
    plt.scatter(age_bin[d1].in_units(units),age_bin[d2].in_units(units), c=abs(age_bin['tform'].in_units('Gyr')),marker = '.', edgecolors= 'none', alpha = 0.5)
    cbar = plt.colorbar()
    cbar.set_label('Gyr')
    plt.axis((-50,50,-50,50))
    plt.xlabel('%s [kpc]'%d1)
    plt.ylabel('%s [kpc]'%d2)

#same as plot_tform, but it makes a smaller rectangle. For side on views.
def plot_tform_rect(sim, age_mi, age_ma, d1, d2, units):
    age_bin = sim[(sim['tform'].in_units('Gyr') >= age_mi) &
                  (sim['tform'].in_units('Gyr') < age_ma)]
    fig = plt.figure()
    ax = fig.add_axes([.125,.1,.625,0.25])
    ax.set_axis_bgcolor((0,0,0.05))
    ax.scatter(age_bin[d1].in_units(units),age_bin[d2].in_units(units), c=abs(age_bin['tform'].in_units('Gyr')),marker = '.', edgecolors= 'none', alpha = 0.5)
    plt.axis((-50,50,-20,20))
    plt.xlabel('%s [kpc]'%d1)
    plt.ylabel('%s [kpc]'%d2)

#I think this is the same as plot_tform, but colors based on velocity
def plot_mag(sim, age_mi, age_ma, d1, d2, units):
    age_bin = sim[(sim['tform'].in_units('Gyr') >= age_mi) &
                  (sim['tform'].in_units('Gyr') < age_ma)]
    plt.figure()
    color = np.sqrt(np.sum((age_bin['vel'].in_units(units))**2, axis = 1))
    plt.scatter(age_bin[d1].in_units('kpc aform'),age_bin[d2].in_units('kpc aform'), c=color,marker = '.',cmap=cm.get_cmap('copper'), edgecolors= 'none', alpha = 0.5)
    cbar = plt.colorbar()
    cbar.set_label(units)
    plt.xlabel('%s [kpc]'%d1)
    plt.ylabel('%s [kpc]'%d2)

#same as plot_tform, but a three dimesional plot.
def plot_tform_3d(sim, age_mi, age_ma):
    age_bin = sim[(sim['tform'].in_units('Gyr') >= age_mi) &
                  (sim['tform'].in_units('Gyr') < age_ma)]
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(age_bin['x'].in_units('kpc aform'),age_bin['y'].in_units('kpc aform'),age_bin['z'].in_units('kpc aform'), c=age_bin['tform'].in_units('Gyr'),marker = '.', edgecolors= 'none', alpha = 0.5)
    #cbar = ax.colorbar()
    #cbar.set_label('Time of Formation [Gyr]')
    ax.axis([-6,6,-6,6])
    ax.set_zlim3d(-6,6)
    ax.set_xlabel('x [kpc]')
    ax.set_ylabel('y [kpc]')
    ax.set_zlabel('z [kpc]')

#returns the center of angluar momentum for a group of particles
def center_values(sim,mass_type):
    cx = np.sum(sim['x'])/len(sim['x'])
    cy = np.sum(sim['y'])/len(sim['y'])
    cz = np.sum(sim['z'])/len(sim['z'])
    cpos = sim['pos'] - [cx,cy,cz]
    vcx = np.sum(sim['vx'])/len(sim['vx'])
    vcy = np.sum(sim['vy'])/len(sim['vy'])
    vcz = np.sum(sim['vz'])/len(sim['vz'])
    cvelx = sim['vx'] - vcx
    cvely = sim['vy'] - vcy
    cvelz = sim['vz'] - vcz
    momentum = np.concatenate((np.reshape(sim[mass_type]*cvelx,(len(cvelx),1)),np.reshape(sim[mass_type]*cvely,(len(cvely),1)),np.reshape(sim[mass_type]*cvelz,(len(cvelz),1))), axis=1)
    angmom = np.sum(np.cross(cpos,momentum),axis=0)
    angmom_norm = angmom/np.sqrt(np.sum(angmom**2))
    return float(cx), float(cy), float(cz), float(vcx), float(vcy), float(vcz), angmom_norm

#aligns a simulation using spline function
def align_spline(sim, alignment_name):
    al_time = []
    al_cx = []
    al_cy = []
    al_cz = []
    al_vcx = []
    al_vcy = []
    al_vcz = []
    al_amvx = []
    al_amvy = []
    al_amvz = []
    al_a = []
    with open(alignment_name) as data:
        for line in data:
            f = line.split()
            al_time.append(float(f[2]))
            a = (1/(1+float(f[3])))
            al_a.append(a)
            al_cx.append(float(f[4]))
            al_cy.append(float(f[5]))
            al_cz.append(float(f[6]))
            al_vcx.append(a*float(f[7]))
            al_vcy.append(a*float(f[8]))
            al_vcz.append(a*float(f[9]))
            al_amvx.append(float(f[10]))
            al_amvy.append(float(f[11]))
            al_amvz.append(float(f[12]))  
            
    u, indices = np.unique(sim['tform'], return_inverse = True)
    posx_list = []
    posy_list = []
    posz_list = []
    for u_time in u:
        
        cx = spline(al_time, al_cx, u_time)
        cy = spline(al_time, al_cy, u_time)
        cz = spline(al_time, al_cz, u_time)
        posx_list.append(cx)
        posy_list.append(cy)
        posz_list.append(cz)
                
    posx_list = np.array(posx_list)
    posy_list = np.array(posy_list)
    posz_list = np.array(posz_list)
    return posx_list,posy_list,posz_list,indices

def unique_rows(data):
    uniq, indices = np.unique(data.view(data.dtype.descr * data.shape[1]), return_inverse = True )
    return uniq.view(data.dtype).reshape(-1, data.shape[1]), indices

def unique_rows(data):
    ncols = data.shape[1]
    dtype = data.dtype.descr * ncols 
    struct = data.view(dtype)

def align_stars3(sim, bin_size, filename):
    age = 13.7345996799999508 - sim['tform'].in_units('Gyr')
    age_max = age.max()
    age_min = age.min()
    cx = []
    cy = []
    cz = []
    vcx = []
    vcy = []
    vcz = []
    amx = []
    amy = []
    amz = []
    time = []
    a = []

    with open(filename) as data:
        for line in data:
            f = line.split()
            time.append(float(f[2]))
            a.append(1/(1+float(f[3])))
            cx.append(float(f[4]))
            cy.append(float(f[5]))
            cz.append(float(f[6]))
            vcx.append(float(f[7]))
            vcy.append(float(f[8]))
            vcz.append(float(f[9]))
            amx.append(float(f[10]))
            amy.append(float(f[11]))
            amz.append(float(f[12]))


    x = sim['tform']

    vcx = np.array(vcx)
    vcy = np.array(vcy)
    vcz = np.array(vcz)
    a = np.array(a)

    a2 = spline(time, a, x)

    vcx2 = spline(time,vcx*a,x)
    vcy2 = spline(time,vcy*a,x)
    vcz2 = spline(time,vcz*a,x)

    cx2 = spline(time,cx,x)
    cy2 = spline(time,cy,x)
    cz2 = spline(time,cz,x)  

    amx2 = spline(time,amx,x)
    amy2 = spline(time,amy,x)
    amz2 = spline(time,amz,x)

    sim['x'] = sim['x']-cx2
    sim['y'] = sim['y']-cy2
    sim['z'] = sim['z']-cz2

    sim['vx'] = (sim['vx'])/a2-vcx2
    sim['vy'] = (sim['vy'])/a2-vcy2
    sim['vz'] = (sim['vz'])/a2-vcz2

#    u_list, indices = unique_rows(np.concatenate((np.reshape(amx2, (len(amx2), 1)),np.reshape(amy2, (len(amy2), 1)),np.reshape(amz2, (len(amz2), 1))), axis=1))
    u_list = np.concatenate((np.reshape(amx2, (len(amx2), 1)),np.reshape(amy2, (len(amy2), 1)),np.reshape(amz2, (len(amz2), 1))), axis=1)

    t1_list = []
    t2_list = []
    for u in u_list:
        #print u
        t1, t2 = align_sim(u[0], u[1], u[2])
        t1_list.append(t1)
        t2_list.append(t2)

#    t1_list = np.array(t1_list)[indices]
#    t2_list = np.array(t2_list)[indices]
        
    sim['pos'] = np.sum(t2_list*np.sum(t1_list*sim['pos'].reshape(len(sim['pos']),1,3), axis = 2).reshape(len(sim['pos']),1,3), axis = 2)
    sim['vel'] = np.sum(t2_list*np.sum(t1_list*sim['vel'].reshape(len(sim['vel']),1,3), axis = 2).reshape(len(sim['vel']),1,3), axis = 2)

    return sim['pos'],sim['vel'] #,cx2,cy2,cz2,vcx2,vcy2,vcz2,amx2,amy2,amz2, sim['tform']

def align_stars2(sim, bin_size):
    age = 13.7345996799999508 - sim['tform'].in_units('Gyr')
    age_max = age.max()
    age_min = age.min()
    cx = []
    cy = []
    cz = []
    vcx = []
    vcy = []
    vcz = []
    amx = []
    amy = []
    amz = []
    time = []
    a = []
    t1_list = []
    t2_list = []
    for i in range (0, int((age_max-age_min)/bin_size)):
        age_bin = sim[(age >= age_min+(bin_size*i)) &
                      (age < age_min+(bin_size*(i+1)))]
        print len(age_bin)#,age_bin['tform'].in_units('Gyr').min()
        al_cx,al_cy,al_cz,al_vcx,al_vcy,al_vcz, al_am = center_values(age_bin, 'massform')
        
        t1, t2 = align_sim(al_am[0], al_am[1], al_am[2])
        cx.append(al_cx)
        cy.append(al_cy)
        cz.append(al_cz)
        vcx.append(al_vcx)
        vcy.append(al_vcy)
        vcz.append(al_vcz)
        amx.append(al_am[0])
        amy.append(al_am[1])
        amz.append(al_am[2])
        t1_list.append(t1)
        t2_list.append(t2)
        time.append(np.sum(age_bin['tform'].in_units('Gyr'))/len(age_bin))

    x = sim['tform'].in_units('Gyr')

    vcx2 = spline(time[::-1],vcx[::-1],x)
    vcy2 = spline(time[::-1],vcy[::-1],x)
    vcz2 = spline(time[::-1],vcz[::-1],x)

    cx2 = spline(time[::-1],cx[::-1],x)
    cy2 = spline(time[::-1],cy[::-1],x)
    cz2 = spline(time[::-1],cz[::-1],x)  

    amx2 = spline(time[::-1],amx[::-1],x)
    amy2 = spline(time[::-1],amy[::-1],x)
    amz2 = spline(time[::-1],amz[::-1],x)

    sim['x'] = sim['x']-cx2
    sim['y'] = sim['y']-cy2
    sim['z'] = sim['z']-cz2

    sim['vx'] = sim['vx']-vcx2
    sim['vy'] = sim['vy']-vcy2
    sim['vz'] = sim['vz']-vcz2

    u_list, indices = unique_rows(np.concatenate((np.reshape(amx2, (len(amx2), 1)),np.reshape(amy2, (len(amy2), 1)),np.reshape(amz2, (len(amz2), 1))), axis=1))


    t1_list = []
    t2_list = []
    for u in u_list:
        #print u
        t1, t2 = align_sim(u[0], u[1], u[2])
        t1_list.append(t1)
        t2_list.append(t2)

    t1_list = np.array(t1_list)[indices]
    t2_list = np.array(t2_list)[indices]

    sim['pos'] = np.sum(t2_list*np.sum(t1_list*sim['pos'].reshape(len(sim['pos']),1,3), axis = 2).reshape(len(sim['pos']),1,3), axis = 2)
    sim['vel'] = np.sum(t2_list*np.sum(t1_list*sim['vel'].reshape(len(sim['vel']),1,3), axis = 2).reshape(len(sim['vel']),1,3), axis = 2)

    return sim['pos'],sim['vel'],cx2,cy2,cz2,vcx2,vcy2,vcz2,amx2,amy2,amz2, sim['tform']
   

def align_stars(sim, bin_size):
    age = 13.7345996799999508 - sim['tform'].in_units('Gyr')
    age_max = age.max()
    age_min = age.min()
    cx = []
    cy = []
    cz = []
    vcx = []
    vcy = []
    vcz = []
    amx = []
    amy = []
    amz = []
    time = []
    for i in range (0, int((age_max-age_min)/bin_size)):
        age_bin = sim[(age >= age_min+(bin_size*i)) &
                      (age < age_min+(bin_size*(i+1)))]
                   
    ## t = np.arange(0, len(sim), bin_size)
    ## if len(sim)%bin_size != 0:
    ##     t = np.append(t,t.max()+(len(sim)%bin_size)-1)
    ## for i in range (0, len(t)-1):
    ##     if t[i+1] == max(t): 
    ##         age_bin = sim[t[i]:t[i+1]+1]
    ##     else:
    ##         age_bin = sim[t[i]:t[i+1]]
    
## bin_size = 10
## u, indices = np.unique(sim['tform'].in_units('Gyr'), return_inverse = True)
## ulist = np.arange(0, len(u), bin_size)
## if len(u)%bin_size != 0:
##    ulist = np.append(ulist,ulist.max()+(len(u)%bin_size)-1)
## for i in range(0, len(ulist)-1):
##    if u[ulist[i+1]] == max(u):
##        age_bin = sim[(sim['tform'].in_units('Gyr') >= u[ulist[i]])&(sim['tform'].in_units('Gyr') <= u[ulist[i+1]])]
##    else:
##        age_bin = sim[(sim['tform'].in_units('Gyr') >= u[ulist[i]])&(sim['tform'].in_units('Gyr') < u[ulist[i+1]])]

        print len(age_bin),age_bin['tform'].in_units('Gyr').min()
        al_cx,al_cy,al_cz,al_vcx,al_vcy,al_vcz, al_am = center_values(age_bin, 'massform')
        
        t1, t2 = align_sim(al_am[0], al_am[1], al_am[2])
        
        age_bin['x'] = age_bin['x']-al_cx
        age_bin['y'] = age_bin['y']-al_cy
        age_bin['z'] = age_bin['z']-al_cz
        
        age_bin['vx'] = age_bin['vx']-al_vcx
        age_bin['vy'] = age_bin['vy']-al_vcy
        age_bin['vz'] = age_bin['vz']-al_vcz
        
        age_bin['pos'] = np.dot(np.dot(age_bin['pos'],t1),t2)
        age_bin['vel'] = np.dot(np.dot(age_bin['vel'],t1),t2)
        
        k = np.in1d(sim['iord'],age_bin['iord'])
        sim['pos'][k] = age_bin['pos']
        sim['vel'][k] = age_bin['vel']
        cx.append(al_cx)
        cy.append(al_cy)
        cz.append(al_cz)
        vcx.append(al_vcx)
        vcy.append(al_vcy)
        vcz.append(al_vcz)
        amx.append(al_am[0])
        amy.append(al_am[1])
        amz.append(al_am[2])
        time.append(np.sum(age_bin['tform'].in_units('Gyr'))/len(age_bin))

    return sim['pos'], sim['vel'],cx,cy,cz,vcx,vcy,vcz,amx,amy,amz, time

def units(value, unit):
    unit_dict ={
        'kpc':2.5e+04,
        'km s**-1':6.30e+02,
        'Gyr':3.88e+01
        }
    return value*unit_dict[unit]

#opens fits starlogs
def fits_starlog(filename, sim_number):
#    from astropy.io import fits
#    sl = fits.open(filename)
    import pyfits
    sl = pyfits.open(filename)
    hsl = sl[1].data

    units = {'239':['5.00e+04 kpc aform','1.26e+03 km s**-1 aform','3.97e+01 Gyr','1.85e+16 Msol'], 
              '258':['5.00e+04 kpc aform','1.26e+03 km s**-1 aform','3.97e+01 Gyr','1.85e+16 Msol'], 
              '285':['5.00e+04 kpc aform','1.26e+03 km s**-1 aform','3.97e+01 Gyr','1.85e+16 Msol'],
              '516':['2.5e+04 kpc aform','6.30e+02 km s**-1 aform','3.88e+01 Gyr','2.31e+15 Msol'],
              '516m':['2.5e+04 kpc aform','6.30e+02 km s**-1 aform','3.88e+01 Gyr','2.31e+15 Msol'],
              '799':['2.5e+04 kpc aform','6.30e+02 km s**-1 aform','3.88e+01 Gyr','2.31e+15 Msol'],
              '986':['5.00e+04 kpc aform','1.26e+03 km s**-1 aform','3.97e+01 Gyr','1.85e+16 Msol']}

    h1sl = pynbody.new(n_particles=len(hsl))
    h1sl['pos'] = pynbody.array.SimArray(np.concatenate((pynbody.array.SimArray(hsl['x'], units[sim_number][0]).reshape(len(hsl),1),(pynbody.array.SimArray(hsl['y'], units[sim_number][0]).reshape(len(hsl),1)),(pynbody.array.SimArray(hsl['z'], units[sim_number][0]).reshape(len(hsl),1))), axis = 1), units[sim_number][0])
    h1sl['vel'] = pynbody.array.SimArray(np.concatenate((pynbody.array.SimArray(hsl['vx'], units[sim_number][1]).reshape(len(hsl),1),(pynbody.array.SimArray(hsl['vy'], units[sim_number][1]).reshape(len(hsl),1)),(pynbody.array.SimArray(hsl['vz'], units[sim_number][1]).reshape(len(hsl),1))), axis = 1),units[sim_number][1])
    h1sl['tform'] = pynbody.array.SimArray(hsl['timeform'],units[sim_number][2])
    h1sl['iord'] = pynbody.array.SimArray(hsl['iorderstar'])
    h1sl['massform'] = pynbody.array.SimArray(hsl['massform'], units[sim_number][3])
    
    sl.close
    return h1sl

def snap_starlog(file_name):
    import sys, os, glob, pynbody.bridge
    import pynbody.snapshot.tipsy
    x = os.path.abspath(s._filename)
    x = os.path.dirname(x)
    l = glob.glob(os.path.join(x,"*.starlog"))
    if (len(l)) :
        for filename in l :
            sl = pynbody.tipsy.StarLog(filename)
    return sl

#I think this returns the velocity dispersion of the ISM for a simulation
#sim_number: simulation number
#s is the simulation
def ISM(sim_number, s):
    pathbase = '/home/christensen/Storage1/UW/MolecH/Cosmo/'

    path = {'239':pathbase + 'h239.cosmo50cmb.3072g/',
               '258':pathbase + 'h258.cosmo50cmb.3072g/',
               '285':pathbase + 'h285.cosmo50cmb.3072g/',
               '516':pathbase + 'h516.cosmo25cmb.3072g/',
               '516m':pathbase + 'h516.cosmo25cmb.3072g/',
               '603':pathbase + 'h603.cosmo50cmb.3072g/',
               '799':pathbase + 'h799.cosmo25cmb.3072g/',
               '986':pathbase + 'h986.cosmo50cmb.3072g/'}
    
    prof = {'239':path[sim_number] + 'h239.cosmo50cmb.3072g14HMbwK/h239.cosmo50cmb.3072g14HMbwK.1_prof.txt',
               '258':path[sim_number] + 'h258.cosmo50cmb.3072g14HMbwK/h258.cosmo50cmb.3072g14HMbwK.1_prof.txt',
               '285':path[sim_number] + 'h285.cosmo50cmb.3072g14HMbwK/h285.cosmo50cmb.3072g14HMbwK.1_prof.txt',
               '516':path[sim_number] + 'h516.cosmo25cmb.3072g14HBWK/h516.cosmo25cmb.3072g14HBWK.1_prof.txt',
               '603':path[sim_number] + 'h603.cosmo50cmb.3072g14HBWK/h603.cosmo50cmb.3072g14HBWK.1_prof.txt',
               '799':path[sim_number] + 'h799.cosmo25cmb.3072g14HBWK/h799.cosmo25cmb.3072g14HBWK.1_prof.txt',
               '986':path[sim_number] + 'h986.cosmo50cmb.3072g14HBWK/h986.cosmo50cmb.3072g14HBWK.1_prof.txt'}

    redshift = []
    V_circ = []
    firstline = True

    with open(prof[sim_number]) as data:
        for line in data:
            f = line.split()
            if firstline == False:
                V_circ.append(float(f[3]))
                redshift.append(float(f[2]))
            else:
                firstline = False

    #calculates sigma                
    V_rot_prime = []
    sigma = []
    years = []
    for i in range(0, len(redshift)):
        z = redshift[i]
        m_log = 10.5
        t_dep = 1.5*math.pow(1+z,-1)
        a = -1.73 + 1.26/(1+math.exp((10.49-m_log)/-.25))
        b = 1.85 + 1.57/(1+math.exp((10.35-m_log)/.19))
        sfr = math.pow(10, a)*math.pow(1+z,b)
        f_gas = 1/(1 + math.pow(t_dep*sfr,-1))
        V_rot_prime = math.sqrt(math.pow(V_circ[i],2)/((1+math.pow(f_gas,2))*(3-1)))
        sigma.append((1/math.sqrt(2))*V_rot_prime*f_gas)

    u_age = pynbody.analysis.cosmology.age(s)
    for r in redshift:
        years.append(u_age - pynbody.analysis.cosmology.age(s, z = r))

    return sigma, years

#finds velocity dispersion due to scattering.
def ISM2(sim_number, si, v_disp2, v_age2):
    sim = si[(si['temp']<1e4) & (si['rho'].in_units('g cm**-3') > 1.66e-24)]
    G = pynbody.array.SimArray(6.674e-8)
    G.units = 'cm**3 g**-1 s**-2'
    rho = sim['rho'].in_units('g cm**-3')
    k_b = pynbody.array.SimArray(1.380658e-16)
    k_b.units = 'g cm**2 s**-2 K**-1'
    T = sim['temp'].in_units('K')
    HII =  max(sim['HI']+sim['H2']) - sim['HI'] - sim['H2']
    HeIII =  (1-HII) - sim['HeI'] - sim['HeII']
    mu = pynbody.array.SimArray((1.67e-24*sim['HI']) + ((2*1.67e-24)*sim['H2']) + ((.5*1.67e-24)*(HII)) + ((4*1.67e-24)*sim['HeI']) + ((2*1.67e-24)*sim['HeII']) + (((4/3)*1.67e-24)*(HeIII)))
    mu.units = 'g'

    R_j = .5*np.sqrt((15*k_b*T)/(4*math.pi*G*mu*rho))

    mass_j = ((4*math.pi)/3)*rho*(R_j**3)

    sigma_not = np.array(v_disp2) #initial birth velocity of stars
    t = v_age2
    t2 = np.array(t)*10e9
    #n_c = pynbody.array.SimArray(3.078278150751079e-63)#pynbody.array.SimArray(1.5391390753755398e-63)#
    #n_c = pynbody.array.SimArray(min(sim['rho'].in_units('g cm**-3')/sim['mass'].in_units('g')))
    #n_c.units = 'cm**-3'
    M_c = pynbody.array.SimArray(min(mass_j))
    #pynbody.array.SimArray(1.9491415999999998e+39)#pynbody.array.SimArray(9.745708e+38) #pynbody.array.SimArray(max(mass_j))#characteristic mass
    M_c.units = 'g'

    n_c = pynbody.array.SimArray(3e-24/(M_c))
    n_c.units = 'cm**-3'
    ln_C = 3 #not sure what this is
    F_B_w = 0.98 #constant?
    gamma = (3*(math.pi**(3/2))*(G*G)*(M_c*M_c)*n_c*ln_C*F_B_w)/4
    
    sigma = (sigma_not**3 + (3/2)*gamma.in_units('km**3 s**-3 yr**-1')*t2)**(1/3.0)

    #fig = plt.figure()
    #ax = fig.add_subplot(1,1,1)
    #ax.scatter(T, mass_j, marker='.', edgecolors='none', alpha=0.5)
    #ax.set_yscale('log')
    #ax.set_xscale('log')
    #ax.set_xlim(1e-30,1e-21)
    #ax.set_ylim(1e12,1e21)
    #plt.ylabel('Jeans mass and Jeans length relation [g**2 cm**-3]')
    #plt.xlabel('Temp [g/cm]')
    #plt.title('%s' %sim_number)


    return sigma, np.array(t)

#uses the softening length for find scattering data
def gas_part_grav_sof(sim_number, v_disp2, v_age2):
    particle_mass = {'239':2.7e4,
                     '258':2.7e4,
                     '285':2.7e4,
                     '516':3.3e3,
                     '799':3.3e3,
                     '986':2.7e4}
    softening_length = {'239':170,
                        '258':170,
                        '285':170,
                        '516':87,
                        '799':87,
                        '986':170}
    G = pynbody.array.SimArray(6.674e-8)
    G.units = 'cm**3 g**-1 s**-2'

    R_j = pynbody.array.SimArray(softening_length[sim_number])
    R_j.units = 'pc'

    mass_j = pynbody.array.SimArray(particle_mass[sim_number])
    mass_j.units = 'Msol'
    
    sigma_not = np.array(v_disp2) #initial birth velocity of stars
    t = v_age2
    t2 = np.array(t)*10e9

    n_c = pynbody.array.SimArray(3e-24/(mass_j.in_units('g')))
    n_c.units = 'cm**-3'
    print 'n_c'
    print n_c
    ln_C = 3 #not sure what this is
    F_B_w = 0.98 #constant?
    gamma = (3*(math.pi**(3/2))*(G*G)*(mass_j.in_units('g')*mass_j.in_units('g'))*n_c*ln_C*F_B_w)/4
    print 'gamma'
    print gamma.in_units('km**3 s**-3 yr**-1')
    
    sigma = (sigma_not**3 + (3/2)*gamma.in_units('km**3 s**-3 yr**-1')*t2)**(1/3.0)
    return sigma, np.array(t)

def data_table(filename, sim_number, v1, v2, v3):
    f = open(filename, 'w')
    f.truncate()

    f.write('h'+sim_number+'\n')
    f.write('Time_in_Gyrs\t Sigma_{ISM}\t Sigma_{star,z=0}\n')
    print len(v1)
    print len(v2)
    print len(v3)
    for n in range(0, len(v1)):
        f.write(str(v1[n])+'\t'+str(v2[n])+'\t'+str(v3[n])+'\n')
    f.close()
#************************************************************************
