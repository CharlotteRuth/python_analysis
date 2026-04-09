import pynbody
import glob
import numpy as np
import os

import sys
sys.path.append("/home/christenc/Code/python_analysis/")
import write_tipsyarr

path = "/home/christenc/REPOSITORY/e12Gals/h148.cosmo50PLK.6144g3HbwK1BH/snapshots_200crit_h148mint/"
filename = 'h148.cosmo50PLK.6144g3HbwK1BH'
step_addmass = ['000896','003360','003456','003744']
step_addigord = ['000777','003360','003456','003744']

finalstep = 'h148.cosmo50PLK.6144g3HbwK1BH.004096'

sim_final = pynbody.load(path + finalstep)

os.chdir(path)
for step in step_addmass:
    sim = pynbody.load(path + filename + '.' + step)
    sim.gas['massform'] = np.zeros(len(sim.gas))
    sim.dark['massform'] = np.zeros(len(sim.dark))
    sim.star['massform'] = (sim_final.star['massform'])[0:len(sim.star)]
    write_tipsyarr.write_tipsyarr(sim, sim['massform'], ext="massform", typefmt="e")
    del sim

del sim_final['massform']

for step in step_addigord:
    sim = pynbody.load(path + filename + '.' + step)
    sim.gas['igasorder'] = np.zeros(len(sim.gas))
    sim.dark['igasorder'] = np.zeros(len(sim.dark))
    sim.star['igasorder'] = (sim_final.star['igasorder'])[0:len(sim.star)]
    write_tipsyarr.write_tipsyarr(sim, sim['igasorder'], ext="igasorder")
    del sim

del sim_final['igasorder']
