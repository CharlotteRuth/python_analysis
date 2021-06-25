#Run AHF on all the files within a directory

import pynbody
import sys, os, glob

filename_base = "h329.cosmo50PLK.3072gst5HbwK1BH"
filename_base = "h242.cosmo50PLK.3072gst5HbwK1BH"
filename_base = "h229.cosmo50PLK.3072gst5HbwK1BH"
filename_base = "h148.cosmo50PLK.3072g3HbwK1BH"
snapshots = glob.glob(filename_base + ".00????")

for tfile in snapshots:
    if os.path.exists(tfile + '.*.AHF_fpos'):
        continue
    print(tfile)
    s = pynbody.load(tfile)
    h = s.halos()
    if os.path.exists('AHF.in'):
        os.system("mv AHF.in " + tfile + ".AHF.in")

    
