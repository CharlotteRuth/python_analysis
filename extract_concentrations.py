# %% Initialization Cell

import time

print("Loading libraries...") #end and flush removed
t1 = time.time()

import sys, os, gc
import numpy as np

sys.path.append(os.path.abspath("."))
from modules.user_tools import task

t2 = time.time()
print("Libraries loaded in "+str( round(t2-t1, 3) )+"s.")

#clear any unused variables in RAM
gc.collect()

print("Working directory:", os.getcwd())

data = task(
    np.load,
    start_text="Loading data",
    end_text="Loaded data",
    fail_text="Failed to load data",
    exit_on_fail=True
)("reduced_time_series_data10.npy", encoding="latin1", allow_pickle=True).item()

new_data = { volume: { halo_grp: {} for halo_grp in data[volume].keys() } for volume in data.keys() }

Marvel_volumes = ["cptmarvel","elektra","rogue","storm"]

volumes = data.keys()
for volume in Marvel_volumes:
    halo_grps = data[volume].keys()
    for halo_grp in halo_grps:
        halo = data[volume][halo_grp]

        labels = ["concentration","scale radius"] #,"DM core slope"]
        for label in labels:
            try:
                item = data[volume][halo_grp][label]["data"][0]

                if label=="concentration" and isinstance(item,dict):
                    if 'r_s' in item.keys():
                        print(f"{volume} {halo_grp} '{label}' has the output for 'scale radius' for some reason.")
                    else:
                        print(f"{volume} {halo_grp} '{label}' has the output for something else.")
                    continue
                elif label=="scale radius" and isinstance(item,float):
                    if isinstance(item,float):
                        print(f"{volume} {halo_grp} '{label}' has the output for 'concentration' for some reason.")
                    else:
                        print(f"{volume} {halo_grp} '{label}' has the output for something else.")
                    continue

                # if label == "scale radius":
                #     print("r_s:",item["r_s"])


                new_data[volume][halo_grp][label] = item
            except:
                print(f"{volume} halo {halo_grp} is missing '{label}'.")

task(
    np.save,
    start_text="Saving data",
    end_text="Saved data",
    fail_text="Failed to save data",
    exit_on_fail=True
)("plots/calculate_beta/concentrations_and_einasto.npy", new_data)
