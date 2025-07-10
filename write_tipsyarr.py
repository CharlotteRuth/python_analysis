import pynbody
import glob
import numpy as np
import os

"""
sim.gas['iord'] = np.empty(len(sim.gas))
sim.dark['iord'] = sim3.dark['iord']
sim.star['iord'] = (sim3.star['iord'])[0:len(sim.star)]

write_tipsyarr(sim, sim[iord], ext = "iord")
"""

# typefmt="d" for integers, "f" for float, "b" for binary, and "s" string
def write_tipsyarr(sim, arr, ext, typefmt="d"):
    #halos = sim.halos(dummy=True, dosort=True)
    formatarr = "%" + str(len(str(np.max(arr))) + 1) + typefmt
    arr_write = np.append(np.array([len(arr)]),arr)
    np.savetxt(sim.filename + '.' + ext, arr_write, fmt = formatarr) 
