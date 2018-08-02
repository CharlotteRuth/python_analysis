import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pynbody, sys
from pynbody.analysis import profile, angmom, halo
from pynbody import units, config
import pynbody.filt as f
import socket

from pynbody.analysis import luminosity as lum
import os, glob

#Run with
#%run /home/christensen/Code/python/python_analysis/mapprof.py
#or, on quirm,
#%run /home/christenc/Code/python/python_analysis/mapprof.py
#ipython --pylab

def make_rs(im,width): #width should be width of box in kpc
    xsize, ysize = np.shape(im)
    x = np.arange(-xsize/2, xsize/2)/float(xsize)*width
    y = np.arange(-ysize/2, ysize/2)/float(xsize)*width
    xs, ys = np.meshgrid(x,y)
    # 2.0 because using 2 kpc pixels
    return 2.0*np.sqrt(xs**2 + ys**2)

def mapprof(tfile,halo_nums):
#tfile = sys.argv[1]
    Carbon   = 0.00213
    Oxygen   = 0.00541
    Silicon  = 0.00067
    Iron     = 0.00117

    hfb = pynbody.load(tfile)
    h = hfb.halos()

    for halo_num in halo_nums:
        halo_num = int(halo_num)
       
        hfbsmass = np.sum(h[halo_num].stars['mass'].in_units('Msol'))
        #hfblstar = 10.0**(0.4*(-21 - pynbody.analysis.luminosity.halo_mag(h[halo_num])))
        hfb.physical_units()
        pynbody.analysis.angmom.faceon(h[halo_num])
        #if len(h[halo_num].gas) == 0:
        #    continue
            
        hfbrvir = np.max(h[halo_num].gas['r'])
        notdiskf = f.Not(f.Disc('8 kpc','2 kpc'))
        
        hiif = pynbody.analysis.ionfrac.calculate(hfb.gas,ion='hi')
        hfb.gas['hiden'] = hfb.gas['rho']*hfb.gas['hydrogen']*hiif
        hfb.gas['hiden'].units = hfb.gas['rho'].units
        hfb.gas['oviif'] = pynbody.analysis.ionfrac.calculate(hfb.gas)
        hfb.gas['oviden'] = hfb.gas['rho']*hfb.gas['OxMassFrac']*hfb.gas['oviif']
        hfb.gas['oviden'].units = hfb.gas['rho'].units
        hfb.gas['civif'] = pynbody.analysis.ionfrac.calculate(hfb.gas,ion='civ')
        hfb.gas['civden'] = hfb.gas['rho']*hfb.gas['OxMassFrac']*Carbon/Oxygen*hfb.gas['civif']
        hfb.gas['civden'].units = hfb.gas['rho'].units

        fig = plt.figure(figsize=(12.,12.))
        fig.subplots_adjust(left=0.08, bottom=0.05, right=0.9, top=0.97, wspace=0.12, 
                            hspace=0.13)
    
        sps = [fig.add_subplot(3,3,1), fig.add_subplot(3,3,2), fig.add_subplot(3,3,3),
            fig.add_subplot(3,3,4), fig.add_subplot(3,3,5), fig.add_subplot(3,3,6),
            fig.add_subplot(3,3,7), fig.add_subplot(3,3,8), fig.add_subplot(3,3,9)] 

        width = round(hfbrvir/10)*10*2.5
        # Upper left:  HI map
        hfbhiim = pynbody.plot.image(hfb.gas[notdiskf],qty='hiden', clear=False,
                        units='m_p cm^-2', width=width, show_cbar=False, 
                        subplot=sps[0],vmin=1e13,vmax=1e19, log=True)
        sps[0].set_xlabel('x [kpc]')
        sps[0].set_ylabel('y [kpc]')
        circle_rvir = plt.Circle((0, 0), hfbrvir, color='w',fill=False)
        sps[0].add_artist(circle_rvir)
        
        # Upper center:  OVI map
        hfboviim = pynbody.plot.image(hfb.gas[notdiskf],qty='oviden', clear=False,
                        units='16 m_p cm^-2', width=width, show_cbar=False, 
                        subplot=sps[1],vmin=1e12,vmax=1e15, log=True)
        sps[1].set_xlabel('x [kpc]')
        #sps[1].add_artist(circle_rvir)

        # Upper right: CIV map
        hfbcivim = pynbody.plot.image(hfb.gas[notdiskf],qty='civden', clear=False,
                        units='12 m_p cm^-2', width=width, show_cbar=False, 
                        subplot=sps[2],vmin=1e12,vmax=1e15, log=True)
        sps[2].set_xlabel('x [kpc]')
        #sps[2].add_artist(circle_rvir)
                
        obsrho=np.array([265, 126, 255, 181, 156, 124, 155, 130, 140, 87, 294, 64, 
                        83, 209, 221, 32])
        obsnhi=np.array([13.8,15.14,14.64,13.46, 14.04, 13.64, 14.4, 15.67,
                        15.73, 14.8, 14.83, 15.35, 15.29, 15.01, 14.3, 15.11])
        obsnovi = np.array([13.74, 14.9, 13.71, 13.6, 13.35, 13.77, 14.08, 14.3, 
                            14.16, 14.7, 14.38, 13.75, 14.16, 13.93, 13.61, 13.1])
        
        # Middle left: HI radial profile
        rs = make_rs(hfboviim,width)
        #rad_axes = np.log10(rs.flatten() + 1e-6)
        rad_axes = rs.flatten() + 1e-6
        
        #radius_range = [-1.0,0.2]#[0.3,3]
        radius_range = [0,2.5/2]
        HI_range = [13,19]
        pynbody.plot.hist2d(rad_axes,np.log10(hfbhiim.flatten()),
                            subplot=sps[3], cmap=plt.cm.gist_yarg)
        #pynbody.plot.hist2d(rs.flatten(),hfbhiim.flatten(),xrange=[0.3,3],
        #                    yrange=[12,19], subplot=sps[2], cmap=plt.cm.gist_yarg)
        #sps[3].plot(np.log10(obsrho),obsnhi,'o', label='z=0 0.1 L* < L < L* Prochaska et al (2011)')
        #sps[3].plot(np.log10([hfbrvir,hfbrvir]),[9,19],label=r'$r_{vir}$',color = 'k')
        sps[3].plot([0,0],HI_range,color = 'k')
        #sps[3].legend(loc=0,prop=matplotlib.font_manager.FontProperties(size='small'))
        sps[3].set_ylim(HI_range)
        sps[3].set_xlim(radius_range)
        #sps[3].set_xlabel(r'log$_{10}$ ($\rho$ [kpc])')
        #sps[3].set_xlabel(r'log$_{10}$ ($\rho/R_{vir}$)')
        sps[3].set_xlabel(r'$\rho/R_{vir}$')
        
        # Middle center: OVI radial profile
        OVI_range = [12,15]
        pynbody.plot.hist2d(rad_axes,np.log10(hfboviim.flatten()),
                            subplot=sps[4], cmap=plt.cm.gist_yarg)
        #pynbody.plot.hist2d(rs.flatten(),hfboviim.flatten(),xrange=[0.3,3],
        #                    yrange=[10,17], subplot=sps[3], cmap=plt.cm.gist_yarg)
        #sps[4].plot(np.log10(obsrho),obsnovi,'o', label='0.1 L* < L < L* Prochaska et al (2011)')
        #hfboxmean, hfboxmin, hfboxmax, hfboxbins = pynbody.plot.sph.image_radial_profile(hfboxim)
        #sps[4].plot(np.log10(hfboxbins*2), np.log10(hfboxmean), '-', label='total Ox mean')
        #sps[4].plot(np.log10([hfbrvir,hfbrvir]),OVI_range,label=r'$r_{vir}$')
        sps[4].plot([0,0],OVI_range,color = 'k')
        sps[4].set_ylim(OVI_range)
        sps[4].set_xlim(radius_range)
        #sps[4].set_xlabel(r'log$_{10}$($\rho$ [kpc])')
        #sps[4].set_xlabel(r'log$_{10}$ ($\rho/R_{vir}$)')
        sps[4].set_xlabel(r'$\rho/R_{vir}$')

        # Middle center: CIV radial profile
        CIV_range = [12,15]
        pynbody.plot.hist2d(rad_axes,np.log10(hfbcivim.flatten()),
                            subplot=sps[5], cmap=plt.cm.gist_yarg)
        #pynbody.plot.hist2d(rs.flatten(),hfboviim.flatten(),xrange=[0.3,3],
        #                    yrange=[10,17], subplot=sps[3], cmap=plt.cm.gist_yarg)
        #sps[5].plot(np.log10(obsrho),obsnovi,'o', label='0.1 L* < L < L* Prochaska et al (2011)')
        #hfboxmean, hfboxmin, hfboxmax, hfboxbins = pynbody.plot.sph.image_radial_profile(hfboxim)
        #sps[4].plot(np.log10(hfboxbins*2), np.log10(hfboxmean), '-', label='total Ox mean')
        #sps[5].plot(np.log10([hfbrvir,hfbrvir]),[10,15],label=r'$r_{vir}$')
        sps[5].plot([0,0],CIV_range,color = 'k')
        sps[5].set_ylim(CIV_range)
        sps[5].set_xlim(radius_range)
        #sps[5].set_xlabel(r'log$_{10}$($\rho$ [kpc])')
        #sps[5].set_xlabel(r'log$_{10}$ ($\rho/R_{vir}$)')
        sps[5].set_xlabel(r'$\rho/R_{vir}$')

             
        plt.setp(sps[1].get_yticklabels(), visible=False)
        sps[3].set_ylabel('log$_{10}$ (N[cm$^{-2}$])')
        sps[0].set_title('HI')
        sps[1].set_title('OVI')
        sps[2].set_title('CIV')
        #sps[3].text(0.7,0.9,'L = %4.2f L*'%(hfblstar),transform=sps[3].transAxes)
        sps[2].text(0.05,0.9,'z = %4.2f'%(hfb.properties['z']),transform=sps[2].transAxes,color = 'w')
        sps[0].text(0.05,0.9,'Halo = %i'%(halo_num),transform=sps[0].transAxes,color = 'w')
        
        den_range = [-8,2.0]
        temp_range = [2,9]
        # Total Mass weighted phase diagram
        centre = (0,0,0)
        radius = hfbrvir*1.25
        sphere = hfb[pynbody.filt.Sphere(radius, centre)]
        currenttime = hfb.properties['time'].in_units('yr')
        #coolon = sphere[pynbody.filt.HighPass('coolontime', str(currenttime)+' yr')]
        logrho = np.log10(sphere.gas['rho'].in_units('m_p cm**-3'))
        logtemp = np.log10(sphere.gas['temp'])
        #logrho = np.log10(h[halo_num].gas['rho'].in_units('m_p cm**-3'))
        #logtemp = np.log10(h[halo_num].gas['temp'])
        #logrho_coolon = np.log10(coolon.gas['rho'].in_units('m_p cm**-3'))
        #logtemp_coolon = np.log10(coolon.gas['temp'])
        
        totphaseim = pynbody.plot.hist2d(
            logrho, logtemp, scalemin=1e3, scalemax=1e7, logscale = True,
            weights=sphere.gas['mass'].in_units('Msol'),
            xrange=den_range, 
            ret_im=True,clear=False,subplot=sps[6])
        #totphaseim = pynbody.plot.hist2d(
        #    logrho_coolon, logtemp_coolon, scalemin=1e3, scalemax=1e7, logscale = True,
        #    weights=sphere.gas['mass'].in_units('Msol'),
        #    xrange=den_range, yrange=temp_range,
        #    ret_im=True,clear=False,subplot=sps[6])        
        sps[6].set_ylim(temp_range)
        sps[6].set_xlim(den_range)
        sps[6].set_xlabel('n [cm$^{-3}$]')
        sps[6].set_ylabel('Temperature [K]')
        
        #hiax = fig.add_axes([0.38,0.07,0.02,0.25])
        #cb1 = fig.colorbar(totphaseim,cax=hiax)
        #cb1.set_label('Mass [M$_\odot$])')

        # OVI mass weighted phase diagram

        oviphaseim = pynbody.plot.hist2d(
            logrho, logtemp, scalemin=1e-7, scalemax=1e3, logscale = True,
            weights=sphere.gas['OxMassFrac']*sphere.gas['oviif']*sphere.gas['mass'].in_units('Msol'), 
            xrange=den_range,
            ret_im=True, subplot=sps[7], clear=False)
        #oviphaseim = pynbody.plot.hist2d(
        #    logrho_coolon, logtemp_coolon, scalemin=1e-7, scalemax=1e3, logscale = True,
        #    weights=sphere.gas['OxMassFrac']*sphere.gas['oviif']*sphere.gas['mass'].in_units('Msol'), 
        #    xrange=den_range, yrange=temp_range,
        #    ret_im=True, subplot=sps[7], clear=False)
        sps[7].set_xlabel('n [cm$^{-3}$]')
        sps[7].set_ylim(temp_range)
        sps[7].set_xlim(den_range)

        # CIV mass weighted phase diagram
        civphaseim = pynbody.plot.hist2d(
            logrho, logtemp, scalemin=1e-7, scalemax=1e3, logscale = True,
            weights=sphere.gas['OxMassFrac']*sphere.gas['civif']*sphere.gas['mass'].in_units('Msol'), 
            xrange=den_range,
            ret_im=True, subplot=sps[8], clear=False)
        #civphaseim = pynbody.plot.hist2d(
        #    logrho_coolon, logtemp_coolon, scalemin=1e-7, scalemax=1e3, logscale = True,
        #    weights=sphere.gas['OxMassFrac']*sphere.gas['civif']*sphere.gas['mass'].in_units('Msol'), 
        #    xrange=den_range, yrange=temp_range,
        #    ret_im=True, subplot=sps[8], clear=False)        
        sps[8].set_xlabel('n [cm$^{-3}$]')
        sps[8].set_ylim(temp_range)
        sps[8].set_xlim(den_range)
                
        #oviax = fig.add_axes([0.90,0.07,0.02,0.25])
        #cb2 = fig.colorbar(oviphaseim,cax=oviax)
        #cb2.set_label('Mass O$_{VI}$ [M$_\odot$])')
        #np.savez(tfile+'.phaseims.npz',oviphaseim, totphaseim)
        
        plt.savefig(tfile+'.'+str(halo_num)+'.mapprof.png')
        fig.clf()
        
        #dum = plt.figure()
        #dumsp = fig.add_subplot(1,1,1)
        
        ## Not plotted, but used later:  Total Oxygen Map
        #hfb.gas['oxden'] = hfb.gas['rho']*hfb.gas['OxMassFrac']
        #hfb.gas['oxden'].units = hfb.gas['rho'].units
        #hfboxim=pynbody.plot.sideon_image(hfb.gas,qty='oxden',
        #                                units='16 m_p cm^-2', width=1000, 
        #                                center=False, subplot=dumsp, 
        #                                title='z = %4.2f'%hfb.properties['z'],
        #                                qtytitle='log$_{10}$ N$_{O}$',noplot=True)
        
        
        #np.savez(tfile+'.ims.npz',hfboviim, hfbhiim, hfbsmass,hfblstar, hfbrvir, hfboxim)
        #np.savez(tfile+'.'+string(halo_num)+'.ims.npz',hfboviim, hfbhiim, hfbsmass,hfbrvir, hfboxim)


        

if __name__ == '__main__':    
    if (socket.gethostname() == "quirm"):
        prefix = '/home/christenc/Data/Sims/'
    else:
        prefix = '/home/christensen/Storage2/UW/MolecH/Cosmo/'

    halo_nums = ['1','2','4','5','6','7','10','11','13','14','27']
    tfile = prefix + 'cptmarvel.cosmo25cmb/cptmarvel.cosmo25cmb.4096g5HbwK1BH/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_name = 'cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096'
    #mapprof(tfile,halo_nums)
    
    halo_nums = ['1','3','7','8','10','11','12','16','17','18','30','32','34','36','61','123']
    tfile = prefix + 'rogue.cosmo25cmb/rogue.cosmo25cmb.4096g5HbwK1BH/rogue.cosmo25cmb.4096g5HbwK1BH.004096/rogue.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_name = 'rogue.cosmo25cmb.4096g5HbwK1BH.004096'
    mapprof(tfile,halo_nums)

    tfile = prefix + 'elektra.cosmo25cmb/elektra.cosmo25cmb.4096g5HbwK1BH/elektra.cosmo25cmb.4096g5HbwK1BH.004096/elektra.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_name = 'elektra.cosmo25cmb.4096g5HbwK1BH.004096'
    halo_nums = ['1','2','3','4','5','8','9','10','11','12','17','18','37','75']
    mapprof(tfile,halo_nums)

    tfile = prefix + 'storm.cosmo25cmb/storm.cosmo25cmb.4096g5HbwK1BH/storm.cosmo25cmb.4096g5HbwK1BH.004096/storm.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_name = 'storm.cosmo25cmb.4096g5HbwK1BH'
    halo_nums = ['1','2','3','4','5','6','7','8','10','11','13','14','15','16','17','23','24','28','34','35','43','49','50','60','109','124','125','192','208']   
    mapprof(tfile,halo_nums)
