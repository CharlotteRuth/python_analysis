

import numpy as np
import pynbody
import halo_trace as ht
import os
import pandas as pd
import sys
sys.path.append('/home/christenc/Code/halo_trace/')

sim_base = "/home/christenc/REPOSITORY/e12Gals/h148.cosmo50PLK.6144g3HbwK1BH/snapshots_200crit_h148mint/"
trace = ht.tracing.trace_halos(sim_base=sim_base, min_ntot=100)
