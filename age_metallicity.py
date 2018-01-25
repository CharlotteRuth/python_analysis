# Charlotte Christensen
# 1/25/18
# Plot the stellar metallicity vs age relationship for a galaxy with [O/H] and [Fe/H]


#Run with
#%run /home/christensen/Code/python/python_analysis/age_metallicity.py
#or, on quirm,
#%run /home/christenc/Code/python/python_analysis/age_metallicity.py
#ipython --pylab



import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pynbody, sys
from pynbody.analysis import profile, angmom, halo
from pynbody import units, config
import pynbody.filt as f
from pynbody.analysis import luminosity as lum
import os, glob
import socket

XSOLfe=0.117E-2
XSOLO=0.96E-2
XSOLH=0.706


def pltZvsAge(tfile,halo_nums):
#h1 =  pynbody.load(tfile)
#h1.physical_units()
    s = pynbody.load(tfile)
    s.physical_units()
    h = s.halos()

    for halo_num in halo_nums:
        halo = h.load_copy(int(halo_num))

        stars = halo.star
#stars[f.LowPass('OxMassFrac',1e-7)].star['OxMassFrac'] = 1e-7
#stars[f.LowPass('FeMassFrac',1e-8)].star['FeMassFrac'] = 1e-8

#oxh: http://pynbody.github.io/pynbody/_modules/pynbody/snapshot/tipsy.html#oxh
#Sets metallicities of zero to the minimum metallicity of the halo and then calculates [O/H] with Asplund et al 09 solar values
    #plt.scatter(halo.star['age'].in_units('Gyr'),halo.star['oxh']) 
    #plt.axis([0, 14, -4, 0])
    #xlabel('Age [Gyr]')
    #plt.show()
#feh: http://pynbody.github.io/pynbody/_modules/pynbody/snapshot/tipsy.html#feh
#Sets metallicities of zero to the minimum metallicity of the halo and then calculates [Fe/H] with Asplund et al 09 solar values
    #plt.scatter(halo.star['age'].in_units('Gyr'),halo.star['feh'])
    #plt.axis([0, 14, -4, 0])
    #xlabel('Age [Gyr]')
    #plt.show()
    #plt.close()

        nbins = 50
        xmin = -1e-4
        xmax = 14

        dt = 1.0
        agebin = np.arange(xmin,xmax,dt) + dt/2.0
        meanOx = np.empty(np.size(agebin))
        stdevOx = np.empty(np.size(agebin))
        meanFe = np.empty(np.size(agebin))
        stdevFe = np.empty(np.size(agebin))

        for i in range(np.size(agebin)):
            minage = str(agebin[i] - dt/2) + ' Gyr'
            maxage = str(agebin[i] + dt/2) + ' Gyr'
            meanOx[i] = np.log10(np.mean(10**(stars[f.BandPass('age',minage,maxage)].star['oxh'])))
            stdevOx[i] = np.std(stars[f.BandPass('age',minage,maxage)].star['oxh'])
            meanFe[i] = np.log10(np.mean(10**(stars[f.BandPass('age',minage,maxage)].star['feh'])))
            stdevFe[i] = np.std(stars[f.BandPass('age',minage,maxage)].star['feh'])
            
        ymax = max(meanOx + stdevOx) #max(halo.star['oxh']) #0 #-0.85
        ymin = min(meanOx - stdevOx) #min(halo.star['oxh']) #-2.2 #-1.85
        
        aspectratio = (xmax - xmin)/(ymax - ymin)*0.7
        binnedOxMassFrac, xedges, yedges = np.histogram2d(halo.star['oxh'], halo.star['age'].in_units('Gyr'), bins = [nbins,nbins], range = [[ymin,ymax],[xmin,xmax]])
        x = np.linspace(xmin, xmax, num = nbins + 1)
        y = np.linspace(ymin, ymax, num = nbins + 1)
        xcenter = (x[0:-1]+x[1:])/2.0
        ycenter = (y[0:-1]+y[1:])/2.0
        plt.imshow(np.log10(binnedOxMassFrac), cmap="Blues", extent = [xmin, xmax, ymin, ymax], interpolation = 'nearest', origin = 'lower',aspect = aspectratio)
        #plot(agebin,meanOx,'ko')
        plt.errorbar(agebin,meanOx, xerr=dt/2, yerr=stdevOx, fmt ='ko')
        #plt.contour(xcenter, ycenter, np.log10(binnedOxMassFrac))
        plt.axis([xmin,xmax,ymin,ymax])
        plt.xlabel('Age [Gyr]')
        plt.ylabel('[O/H]')
        plt.title(tfile_name + ', ' + halo_num)
        plt.text(1,(ymax - ymin)/5 + ymin,'O/H = ' + "{:.2f}".format(np.log10(np.mean(10**halo.star['oxh']))))
        plt.show()
        plt.savefig(tfile + '.' + halo_num + '.OxVsAge_log.png')
    
        plt.close()

        ymax = max(meanFe + stdevFe) #max(halo.star['feh']) #0 #-0.65 
        ymin = min(meanFe - stdevFe) #min(halo.star['feh']) #-2 ymin = -1.65

        aspectratio = (xmax - xmin)/(ymax - ymin)*0.7
        binnedFeMassFrac, xedges, yedges = np.histogram2d(halo.star['feh'], halo.star['age'].in_units('Gyr'), bins = [nbins,nbins], range = [[ymin,ymax],[xmin,xmax]])
        x = np.linspace(xmin, xmax, num = nbins + 1)
        y = np.linspace(ymin, ymax, num = nbins + 1)
        xcenter = (x[0:-1]+x[1:])/2.0
        ycenter = (y[0:-1]+y[1:])/2.0
        plt.imshow(np.log10(binnedFeMassFrac), cmap="Blues", extent = [xmin, xmax, ymin, ymax], interpolation = 'nearest', origin = 'lower',aspect = aspectratio)
        plt.errorbar(agebin,meanFe, xerr=dt/2, yerr=stdevFe, fmt ='ko')
        #plt.contour(xcenter, ycenter, np.log10(binnedOxMassFrac))
        plt.xlabel('Age [Gyr]')
        plt.ylabel('[Fe/H]')
        plt.title(tfile_name + ', ' + halo_num)
        plt.text(1,(ymax - ymin)/5 + ymin,'Fe/H = ' + "{:.2f}".format(np.log10(np.mean(10**halo.star['feh']))))
        plt.axis([xmin,xmax,ymin,ymax])
        plt.show()
        plt.savefig(tfile + '.' + halo_num + '.FeVsAge_log.png')
    
        plt.close()


if __name__ == '__main__':    
    if (socket.gethostname() == "quirm"):
        prefix = '/home/christenc/Data/Sims/'
    else:
        prefix = '/home/christensen/Storage2/UW/MolecH/Cosmo/'


    #halo_num = '6'
    #tfile = '/home/christensen/Storage2/UW/MolecH/Cosmo/h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g14HBWK/steps/h799.cosmo25cmb.3072g14HBWK_2.00512.dir/h799.cosmo25cmb.3072g14HBWK.00512.halo.' + halo_num
    #tfile_name = 'h799.cosmo25cmb.3072g14HBWK.00512.halo.' + halo_num

    #halo_num = '1'
    #tfile = '/home/christensen/Storage2/UW/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK_2.00512.dir/h516.cosmo25cmb.3072g14HBWK.00512.halo.' + halo_num
    #tfile_name = 'h516.cosmo25cmb.3072g14HBWK.00512.halo.' + halo_num

    #1,2,3,8,15,16
    #halo_num = '1'
    #tfile = '/home/christensen/Storage2/UW/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/h986.cosmo50cmb.3072g14HBWK_2.00512.dir/h986.cosmo50cmb.3072g14HBWK.00512.halo'#.' + halo_num
    #tfile_name = 'h986.cosmo50cmb.3072g14HBWK'#.halo.' + halo_num

    #1,2,3
    #halo_num = '1'
    #tfile = '/home/christensen/Storage2/UW/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.00512.dir/h603.cosmo50cmb.3072g14HBWK.00512'#.halo.' + halo_num
    #tfile_name = 'h603.cosmo50cmb.3072g14HBWK.00512'#.halo.' + halo_num

    halo_nums = ['1','2','4','5','6','7','10','11','13','14','27']
    tfile = prefix + 'cptmarvel.cosmo25cmb/cptmarvel.cosmo25cmb.4096g5HbwK1BH/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_name = 'cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096'

#    tfile = prefix + 'elektra.cosmo25cmb/elektra.cosmo25cmb.4096g5HbwK1BH/elektra.cosmo25cmb.4096g5HbwK1BH.004096/elektra.cosmo25cmb.4096g5HbwK1BH.004096'
#    tfile_name = 'elektra.cosmo25cmb.4096g5HbwK1BH.004096'

    halo_nums = ['1','3','7','8','10','11','12','16','17','18','30','32','34','36','61','123']
    tfile = prefix + 'rogue.cosmo25cmb/rogue.cosmo25cmb.4096g5HbwK1BH/rogue.cosmo25cmb.4096g5HbwK1BH.004096/rogue.cosmo25cmb.4096g5HbwK1BH.004096'
    tfile_name = 'rogue.cosmo25cmb.4096g5HbwK1BH.004096'
    
    pltZvsAge(tfile,halo_nums)

        
