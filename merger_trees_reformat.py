import pandas as pd
import numpy as np
import sys
import os 
import tqdm

import pynbody

pynbody.config['halo-class-priority'] =  [pynbody.halo.ahf.AHFCatalogue,
                                          pynbody.halo.GrpCatalogue,
                                          pynbody.halo.AmigaGrpCatalogue,
                                          pynbody.halo.legacy.RockstarIntermediateCatalogue,
                                          pynbody.halo.rockstar.RockstarCatalogue,
                                          pynbody.halo.subfind.SubfindCatalogue, pynbody.halo.hop.HOPCatalogue]


def get_shortname_base(sim):
    if sim=='h148':
        shortname='h148.cosmo50PLK.3072g'
        base='h148.cosmo50PLK.3072g3HbwK1BH'
    if sim=='h229':
        shortname='h229.cosmo50PLK.3072g'
        base='h229.cosmo50PLK.3072gst5HbwK1BH'
    if sim=='h242':
        shortname='h242.cosmo50PLK.3072g'
        base='h242.cosmo50PLK.3072gst5HbwK1BH'
    if sim=='h329':
        shortname='h329.cosmo50PLK.3072g'
        base='h329.cosmo50PLK.3072gst5HbwK1BH'
    return shortname, base



new_keys = ['h148-12','h148-27','h148-34','h148-38','h148-55','h148-65','h148-251','h148-249','h148-282',
    'h229-14','h229-18','h229-20','h229-22','h229-49','h242-21','h242-38','h242-69','h329-29','h329-117']

new_keys = np.append(new_keys, ['h148-1','h229-1','h242-1','h329-1'])

haloids_h148 = np.char.replace(new_keys[np.char.startswith(new_keys, 'h148')], 'h148-', '').astype(int)
haloids_h148 = dict(zip(haloids_h148, [[]]*len(haloids_h148)))

haloids_h229 = np.char.replace(new_keys[np.char.startswith(new_keys, 'h229')], 'h229-', '').astype(int)
haloids_h229 = dict(zip(haloids_h229, [[]]*len(haloids_h229)))

haloids_h242 = np.char.replace(new_keys[np.char.startswith(new_keys, 'h242')], 'h242-', '').astype(int)
haloids_h242 = dict(zip(haloids_h242, [[]]*len(haloids_h242)))

haloids_h329 = np.char.replace(new_keys[np.char.startswith(new_keys, 'h329')], 'h329-', '').astype(int)
haloids_h329 = dict(zip(haloids_h329, [[]]*len(haloids_h329)))

for key in tqdm.tqdm(new_keys):
    sim = key[:4]
    z0haloid = int(key[5:])

    shortname, base = get_shortname_base(sim)

    os.environ['TANGOS_DB_CONNECTION'] = '/home/akinshol/Tangos-2022/databases/JL.db'
    os.environ['TANGOS_SIMULATION_FOLDER'] = f'/home/christenc/Data/Sims/{shortname}/{base}/'
    import tangos

    print(f'Computing merger trees for {sim}-{z0haloid}')


    tangos_halo = tangos.get_halo(f'snapshots_200crit_{sim}/{base}.004096/halo_{z0haloid}')
    time, haloid = tangos_halo.calculate_for_progenitors("t()","halo_number()")

    print(haloid)
    
    exec(f'haloids_{sim}[z0haloid] = haloid')
    
import glob
for sim in ['h148','h229','h242','h329']:
    shortname, base = get_shortname_base(sim)
    filepaths = np.array(glob.glob(f'/home/christenc/Data/Sims/{shortname}/{base}/snapshots_200crit_{sim}/{base}.00????'))
    s = np.argsort([-int(f[-4:]) for f in filepaths])
    filepaths = filepaths[s]
    exec(f'filepaths_{sim} = filepaths')




output = dict(
    filepaths = dict(
        h148 = filepaths_h148, 
        h229 = filepaths_h229, 
        h242 = filepaths_h242, 
        h329 = filepaths_h329
    ),
    haloids = dict(
        h148 = haloids_h148, 
        h229 = haloids_h229, 
        h242 = haloids_h242,
        h329 = haloids_h329
    )
)

import pickle
with open('filepaths_haloids_2022.pickle','wb') as f:
    pickle.dump(output, f)
    



#assert 1==2
#tangos_halo = tangos.get_halo("snapshots_200crit_" + simshort + '/' + base + f'.004096/halo_{haloid}')
#time, a, id_main, mvir_main, rvir_main = tangos_halo.calculate_for_progenitors("t()", "a()", "halo_number()",
#"Mvir", "Rvir")




