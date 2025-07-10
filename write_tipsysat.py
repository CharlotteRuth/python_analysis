# Writes the .amiga.stat etc files from the AHF output

import pynbody
import glob
import numpy as np

#filebase = './h148.cosmo50PLK.3072g3HbwK1BH'
#filebase = './h242.cosmo50PLK.3072gst5HbwK1BH'
#filebase = './h229.cosmo50PLK.3072gst5HbwK1BH'
#filebase = './h329.cosmo50PLK.3072gst5HbwK1BH'

filebase = './h329.cosmo50PLK.3072gst'
filebase = './h229.cosmo50PLK.3072gst'
filebase = './h242.cosmo50PLK.3072gst'
filebase = './h148.cosmo50PLK.3072gst'

filebase = './h148.cosmo50PLK.6144g3HbwK1BH'

#filebase = 'cptmarvel.cosmo25cmb.4096'

files = glob.glob(filebase + '.00????.*AHF_halos')

ct = 0
for file_i in files:
    file_i = '.'.join((file_i.split('.'))[0:5])
    print(ct)
    ct = ct + 1
    if len(glob.glob(file_i + '*.amiga.stat')) != 0:
        continue
    print(file_i)
    sim = pynbody.load(file_i)
    halos = sim.halos(dummy=True, dosort = True)
    halos.writestat(sim, halos, sim.filename+ ".amiga.stat")

    snapshot = sim.ancestor                           
    try:
        snapshot['grp']
    except:
        halos.make_grp()
    grpoutfile = sim.filename + '.amiga.grp'
    formatarr = "%" + str(len(str(halos._nhalos)) + 1) + "d"
    grp = np.append(np.array([len(snapshot['grp'])]),snapshot['grp'])
    np.savetxt(grpoutfile, snapshot['grp'], fmt = formatarr)

    halos.writetipsy(sim, halos, sim.filename + ".amiga.gtp") 


#edit the amiga.stat file to remove the "?" from Satellite?
    
# Because writegrp (and, therfore, writehalos) doesn't work
"""
    snapshot = sim.ancestor
    try:
        snapshot['grp']
    except:
        halos.make_grp()
    grpoutfile = sim.filename + '.amiga.grp'
    formatarr = "%" + str(len(str(halos._nhalos)) + 1) + "d"
    grp = np.append(np.array([len(snapshot['grp'])]),snapshot['grp'])
    np.savetxt(grpoutfile, snapshot['grp'], fmt = formatarr)

    halos.writetipsy(sim, halos, sim.filename + ".amiga.gtp")
"""
