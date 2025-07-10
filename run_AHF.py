#Run AHF on all the files within a directory (use in the snapshots/snapshots_200bkgdens directory)

# mv snapshots_200bkgdens snapshots_200bkgdens_suspect
# mkdir snapshots_200bkgdens
# cd snapshots_200bkgdens
# base=h229.cosmo50PLK.3072gst5HbwK1BH
# ln -s ../$base.00????/ahf_200/* .
# steps=`ls -d ../$base.00???? | cut -d '.' -f6`
# for i in $steps
# do
#     ln -s ../$base.$i/$base.$i ../$base.$i/$base.$i.[cEFHilmMO]* .
# done

import pynbody
import numpy as np
import sys, os, glob

filename_base = "h148.cosmo50PLK.6144g3HbwK1BH"
filename_base = "h329.cosmo50PLK.3072gst5HbwK1BH"
#filename_base = "h242.cosmo50PLK.3072gst5HbwK1BH"
#filename_base = "h229.cosmo50PLK.3072gst5HbwK1BH"
filename_base = "h148.cosmo50PLK.3072g3HbwK1BH"
filename_base = "h148.cosmo50PLK.3072gst"
filename_base = "h229.cosmo50PLK.3072gst"
filename_base = "h242.cosmo50PLK.3072gst"
filename_base = "h329.cosmo50PLK.3072gst"

snapshots = glob.glob(filename_base + ".00????")
path = os.getcwd()
path = '/'.join(path.split('/')[:-1])
for tfile in snapshots:
    print(path + '/' + tfile + '/ahf_200/' + tfile + '.*.AHF_halos')
    if len(glob.glob(path + '/' + tfile + '/ahf_200/' + tfile + '.*.AHF_halos')) != 0:
        continue
    s = pynbody.load(tfile)
    #if os.path.exists('AHF.in'):
    #os.system("cp ../AHF.in AHF.in")
    mineps = min(s.dark['eps'].in_units('kpc a')) #Following Alyson here on comoving vs physical units for softening
    print(s.properties['boxsize'].in_units('kpc a'),s.properties['a'],mineps,(s.properties['boxsize'].in_units('kpc a')*s.properties['a']/mineps),np.log10(s.properties['boxsize'].in_units('kpc a')*s.properties['a']/mineps))
    LgridMax = 2.0**round(np.log10(s.properties['boxsize'].in_units('kpc a')*s.properties['a']/mineps)/np.log10(2.0)) #2.0*round(np.log10(boxsize*h.time/soft/np.log10(2.0)))

    f = open(r"AHF.in", "w")
    f.write("[AHF] \n")
    f.write("ic_filename = " + tfile + "\n")
    f.write("ic_filetype = 90 \n")
    f.write("outfile_prefix = " + tfile + "\n")
    f.write("LgridDomain = 512\n") #256\n")
    f.write("LgridMax = " + str(int(LgridMax)) + "\n")
    f.write("NperDomCell = 8\nNperRefCell = 8\nVescTune = 1\nNminPerHalo = 64\nRhoVir = 0\nDvir = 200\nMaxGatherRad = 1 Mpc/h\n")
    #f.write("NperDomCell = 5\nNperRefCell = 5\nVescTune = 1.5\nNminPerHalo = 50\nRhoVir = 0\nDvir = 200\nMaxGatherRad = 10.0\n")    
    f.write("[TIPSY] \nTIPSY_OMEGA0 = 0.308600 \nTIPSY_LAMBDA0 = 0.691400 \nTIPSY_BOXSIZE = 3.388471e+01 \nTIPSY_VUNIT = 1.170700e+03 \nTIPSY_MUNIT = 1.079489e+16 \nTIPSY_EUNIT = 2.500000e-02\n") #According to Alyson, eunit is always 0.025; see email for discussion
    f.close()
    
    #h = s.halos()

    os.system("mv AHF.in " + tfile + ".AHF.in")


